! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module par_partition_class
  ! Serial modules
  use types
  use memor
  use fem_partition_class
  use fem_import_class
  use fem_element_import_class
  use partition_import
  use maps_class
  use stdio
  ! solve_pardiso_mkl_class is required because 
  ! free_only_struct, free_clean constants are declared
  ! whithin it
  !use solve_pardiso_mkl_class 
  use par_io

  ! Trilinos-Interfaces modules
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use psb_penv_mod
  use par_context_class

# include "debug.i90"
  implicit none
  private

  ! ***IMPORTANT NOTE***: I am assuming that the
  ! constructor of par_partition_class is responsible
  ! for creating/destroying all the members/objects it 
  ! contains

  ! Par_partition will be pointed by all objects
  ! derived from par_mesh (graph + matrix + vector). This
  ! would allow LA algorithms (such as BDDC) to access to 
  ! the data contained in f_partition (objects, objects on the
  ! on the interfaces etc.). Besides, this design minimizes the number of 
  ! changes in the codes as the relationship among derived 
  ! data types is kept. The main disadvantage of this approach
  ! is that f_import, maps_state, row_map, and col_map members
  ! of par_partition are not conceptually related to a mesh, but actually
  ! with linear algebra objects. Can we live with this ?

  integer(ip), parameter :: max_ndofs       = 10
  integer(ip), parameter :: map_non_created =  0
  integer(ip), parameter :: map_created     =  1

  integer(ip), parameter :: overlapped  = 0

  ! (Recursively) distributed partition
  type par_partition
    type (fem_partition)      :: f_part
    type (fem_import)         :: f_import               ! Created from f_part
    type (fem_element_import) :: f_el_import            ! Created from f_part

    type (par_context), pointer  :: p_context => NULL() ! Fine process
    type (par_context), pointer  :: g_context => NULL() ! Fine_to coarse comm
    type (par_context), pointer  :: c_context => NULL() ! Coarse process

    type (par_context), pointer  :: q_context => NULL() ! Available (unused) processes 
    type (par_context), pointer  :: d_context => NULL() ! Available (coarse unused) process
    type (par_context), pointer  :: b_context => NULL() ! Intercommunicator betwen p_context and q_context (bcast and recursive call)
    type (par_context), pointer  :: w_context => NULL() ! World communicator (all process involved).

    ! Recursive reference
    !integer(ip)  :: this_level = 1
    !integer(ip)  :: num_levels = 2

    !integer(ip), allocatable :: lpart(num_levels)

    !type (par_partition), pointer  :: c_part => NULL()

     ! Multilevel basic data
     integer(ip) ::                &
     !   this_level = 1,            &
        num_levels = 2            
     ! If we impose that the last level is always serial these arrays are of size num_levels-1
     ! and we recover the current situation for num_levels=2.
     integer(ip), allocatable ::   &
        id_parts(:),               &    ! Part identifier (num_levels)
        num_parts(:)                    ! Number of parts (num_levels)

    ! Trilinos data, could be avoided using ENABLE_TRILINOS?
    integer(ip)          :: maps_state(max_ndofs)  ! State of maps (non-created vs. created) 
    !type (epetra_map)    :: row_map(max_ndofs)     ! Generated from f_part%nmap+ndofs
    !type (epetra_map)    :: col_map(max_ndofs)     ! Generated from f_part%nmap+ndofs

    ! Object which handles nearest neighbours
    ! communication in order to get data associated
    ! to external vertices in a vertex-based decomposition
    !type (epetra_import) :: importer(max_ndofs)    ! Created from row_map + col_map

  end type par_partition

  interface par_partition_bcast
     module procedure par_partition_bcast_lg, par_partition_bcast_ip, par_partition_bcast_vec
  end interface par_partition_bcast

  interface par_partition_create
     module procedure par_partition_create_standard_execution_model,   &
                      par_partition_create_multilevel_execution_model, &
                      par_partition_create_old
  end interface par_partition_create
  
  interface par_partition_free
     module procedure par_partition_free_one_shot, par_partition_free_progressively
  end interface par_partition_free
  

  ! Types
  public :: par_partition, map_non_created, map_created, max_ndofs

  ! Functions
  public :: par_partition_create, par_partition_free, par_partition_print, par_partition_bcast
  !public :: par_partition_map, par_partition_out

contains

  subroutine par_partition_bcast_lg(p_part,conv)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
    ! Parameters
    type(par_partition), intent(in)    :: p_part
    logical            , intent(inout) :: conv

    ! Locals
    integer :: mpi_comm_b, info


    assert ( associated(p_part%p_context) ) 
    assert ( p_part%p_context%created .eqv. .true.)

    if ( associated(p_part%w_context) ) then
       ! If w_context is associated, then b_context should be also associated
       ! b_context is an intercomm among p_context & q_context
       ! Therefore the semantics of the mpi_bcast subroutine slightly changes
       ! P_0 in p_context is responsible for bcasting conv to all the processes
       ! in q_context

       ! Get MPI communicator associated to icontxt_b (in
       ! the current implementation of our wrappers
       ! to the MPI library icontxt and mpi_comm are actually 
       ! the same)
       call psb_get_mpicomm (p_part%b_context%icontxt, mpi_comm_b)
       
       if (p_part%p_context%iam >=0) then
          if ( p_part%p_context%iam == psb_root_ ) then
             call mpi_bcast(conv,1,MPI_LOGICAL,MPI_ROOT,mpi_comm_b,info)
             check( info == mpi_success )
          else
             call mpi_bcast(conv,1,MPI_LOGICAL,MPI_PROC_NULL,mpi_comm_b,info)
             check( info == mpi_success )
          end if
       else if (p_part%q_context%iam >=0) then
          call mpi_bcast(conv,1,MPI_LOGICAL,psb_root_,mpi_comm_b,info)
          check( info == mpi_success )
       end if
       
    else

       ! Only if global context has been created (i.e.,
       ! there are coarse MPI Tasks) it is required to
       ! broadcast conv
       if ( associated(p_part%g_context) ) then
          assert ( p_part%g_context%created .eqv. .true.)
          
          if (p_part%g_context%iam == 0) then
             call psb_snd(p_part%g_context%icontxt, conv, p_part%g_context%np-1)
             !else if (p_part%g_context%iam == p_part%g_context%np-p_part%c_context%np) then
          else if (p_part%g_context%iam == p_part%g_context%np-1) then
             call psb_rcv(p_part%g_context%icontxt, conv, 0)
          end if
          
          ! IMPORTANT NOTE: task id == 0 in the global context
          !                 MUST BE A FINE TASK (if intracommunicators are used)
          ! call psb_bcast(p_part%g_context%icontxt,conv,root=0)
       end if
    end if

  end subroutine par_partition_bcast_lg

  subroutine par_partition_bcast_ip(p_part,conv)
    implicit none
    type(par_partition), intent(in)    :: p_part
    integer(ip)        , intent(inout) :: conv

    assert ( associated(p_part%p_context) ) 
    assert ( p_part%p_context%created .eqv. .true.)

    ! Only if global context has been created (i.e.,
    ! there are coarse MPI Tasks) it is required to
    ! broadcast conv
    if ( associated(p_part%g_context) ) then
       assert ( p_part%g_context%created .eqv. .true.)

       if (p_part%g_context%iam == 0) then
          call psb_snd(p_part%g_context%icontxt, conv, p_part%g_context%np-1)
       !else if (p_part%g_context%iam == p_part%g_context%np-p_part%c_context%np) then
       else if (p_part%g_context%iam == p_part%g_context%np-1) then
          call psb_rcv(p_part%g_context%icontxt, conv, 0)
       end if

       ! IMPORTANT NOTE: task id == 0 in the global context
       !                 MUST BE A FINE TASK (if intracommunicators are used)
       ! call psb_bcast(p_part%g_context%icontxt,conv,root=0)
    end if

  end subroutine par_partition_bcast_ip

  subroutine par_partition_bcast_vec(p_part,ivec)
    implicit none
    type(par_partition), intent(in)    :: p_part
    integer(ip)        , intent(inout) :: ivec(:)

    assert ( associated(p_part%p_context) ) 
    assert ( p_part%p_context%created .eqv. .true.)

    ! Only if global context has been created (i.e.,
    ! there are coarse MPI Tasks) it is required to
    ! broadcast conv
    if ( associated(p_part%g_context) ) then
       assert ( p_part%g_context%created .eqv. .true.)

       if (p_part%g_context%iam == 0) then
          call psb_snd(p_part%g_context%icontxt, ivec, p_part%g_context%np-1)
       !else if (p_part%g_context%iam == p_part%g_context%np-p_part%c_context%np) then
       else if (p_part%g_context%iam == p_part%g_context%np-1) then
          call psb_rcv(p_part%g_context%icontxt, ivec, 0)
       end if

       ! IMPORTANT NOTE: task id == 0 in the global context
       !                 MUST BE A FINE TASK (if intracommunicators are used)
       ! call psb_bcast(p_part%g_context%icontxt,ivec,root=0)
    end if

  end subroutine par_partition_bcast_vec


  !=============================================================================

!   recursive subroutine par_partition_out(p_part)
!     implicit none
!     type(par_partition), intent(in) :: p_part
!     write(*,*) 'LEVEL: ', p_part%this_level, 'p_context: ',p_part%p_context%iam, p_part%p_context%np
!     write(*,*) 'LEVEL: ', p_part%this_level, 'c_context: ',p_part%c_context%iam, p_part%c_context%np
!     write(*,*) 'LEVEL: ', p_part%this_level, 'g_context: ',p_part%g_context%iam, p_part%g_context%np
!     write(*,*) 'LEVEL: ', p_part%this_level, 'b_context: ',p_part%b_context%iam, p_part%b_context%np
!     if(p_part%this_level<p_part%num_levels-1) call par_partition_out(p_part%c_part)
!   end subroutine par_partition_out

!   subroutine par_partition_map(num_levels,num_procs,map_procs,w_context,p_part)
! #ifdef MPI_MOD
!     use mpi
! #endif
!     implicit none 
! #ifdef MPI_H
!     include 'mpif.h'
! #endif
!     integer(ip)        , intent(in)            :: num_levels
!     integer(ip)        , intent(inout)         :: num_procs(num_levels)
!     integer(ip)        , intent(in)            :: map_procs
!     type(par_context)  , intent(in)            :: w_context
!     type(par_partition), intent(inout), target :: p_part
!     ! Map the partition onto a set of MPI ranks, creating appropriate communicators
!     integer(ip)           :: num_procs_total, my_level, i, ierror
!     integer, allocatable  :: comms(:), fc_comms(:), inter_comms(:)
!     integer               :: my_rank, my_color
!     type(par_partition), pointer :: i_part

!     ! At least one level domain decomposition
!     assert ( num_levels > 0 ) 
!     call memalloc (num_levels, comms, __FILE__,__LINE__ )
!     call memalloc (num_levels, fc_comms, __FILE__,__LINE__ )
!     call memalloc (num_levels, inter_comms, __FILE__,__LINE__ )


!     ! Create communicators according to the desired map, to be coded here.
!     ! We use key=0 in the calls to mpi_comm_split (then the rank order is inherited from parent comm)
!     if(map_procs == overlapped) then

!        my_level=0
!        num_procs_total = 0
!        do i=1,num_levels
!           num_procs_total = num_procs_total + num_procs(i)
!           if(my_level==0 .and. w_context%iam < num_procs_total) my_level = i
!        end do

!        assert ( num_procs_total == w_context%np )
!        i_part => p_part
!        comms(1) = w_context%icontxt
!        do i=1,num_levels-1
!           write(*,*) my_level,i
!           ! First create, for each level, fine and higher level intracomms.
!           ! Also create intercomms between them to implement broadcasts.
!           if(my_level==i) then ! my_color=1, I'm a fine task
!              call mpi_comm_split(comms(i), 1, 0 , comms(i+1), ierror)
!              assert(ierror==MPI_SUCCESS)
!              call mpi_intercomm_create(comms(i+1), 0, comms(i), num_procs(i), i, inter_comms(i+1), ierror)
!              assert(ierror==MPI_SUCCESS)
!           else if(my_level>i) then ! my_color=2, I'm a coarse or higher task
!              call mpi_comm_split(comms(i), 2, 0, comms(i+1), ierror)
!              assert(ierror==MPI_SUCCESS)
!              call mpi_intercomm_create(comms(i+1), 0, comms(i),            0, i, inter_comms(i+1), ierror)
!              assert(ierror==MPI_SUCCESS)
!           end if

!           !call mpi_barrier(w_context%icontxt,ierror)
!           !assert(ierror==MPI_SUCCESS)

!           ! Now create fine to coarse intracomms.
!           ! Here is where the process assignment to levels really occurs.
!           if(my_level==i) then
!              call mpi_comm_rank(comms(i),my_rank,ierror) 
!              assert(ierror==MPI_SUCCESS)
!              my_color = my_rank / num_procs(i+1)
!              call mpi_comm_split(comms(i), my_color, 0 , fc_comms(i+1), ierror)
!           else if(my_level==i+1) then
!              call mpi_comm_rank(comms(i+1),my_color,ierror) 
!              assert(ierror==MPI_SUCCESS)
!              call mpi_comm_split(comms(i), my_color, 0 , fc_comms(i+1), ierror)
!           else if(my_level>i+1) then
!              call mpi_comm_split(comms(i), mpi_undefined, 0 , fc_comms(i+1), ierror)
!           end if

!           !call mpi_barrier(w_context%icontxt,ierror)
!           !assert(ierror==MPI_SUCCESS)
!        end do

!        do i=1,num_levels-1

!           ! Store communicators in par_partition
!           i_part%num_levels = num_levels
!           i_part%this_level = i
!           allocate(i_part%p_context,stat=ierror)
!           assert(ierror==0)
!           allocate(i_part%g_context,stat=ierror)
!           assert(ierror==0)
!           allocate(i_part%c_context,stat=ierror)
!           assert(ierror==0)
!           allocate(i_part%b_context,stat=ierror)
!           assert(ierror==0)
!           if(my_level==i) then                                 ! Fine tasks
!              i_part%p_context%icontxt = comms(i+1)
!              i_part%g_context%icontxt = fc_comms(i+1)
!              i_part%b_context%icontxt = inter_comms(i+1)
!              i_part%c_context%icontxt = mpi_comm_null
!           else if(my_level==i+1) then                          ! Coarse tasks
!              i_part%p_context%icontxt = mpi_comm_null
!              i_part%g_context%icontxt = fc_comms(i+1)
!              i_part%b_context%icontxt = inter_comms(i+1)
!              ! if(i<=num_levels-2) then
!              !    i_part%c_context%icontxt = comms(i+2)
!              ! else
!              !    i_part%c_context%icontxt = comms(i+1)
!              ! end if
!              i_part%c_context%icontxt = mpi_comm_null
!           else if(my_level>i+1) then                           ! Higher tasks
!              i_part%p_context%icontxt = mpi_comm_null
!              i_part%g_context%icontxt = mpi_comm_null
!              i_part%b_context%icontxt = inter_comms(i+1)
!              i_part%c_context%icontxt = mpi_comm_null
!           else                                                 ! Lower tasks
!              i_part%p_context%icontxt = mpi_comm_null
!              i_part%g_context%icontxt = mpi_comm_null
!              i_part%b_context%icontxt = mpi_comm_null
!              i_part%c_context%icontxt = mpi_comm_null
!            end if
!            i_part%p_context%created = .true.
!            i_part%p_context%handler = inhouse
!            call psb_info ( i_part%p_context%icontxt, i_part%p_context%iam, i_part%p_context%np )
!            call psb_info ( i_part%c_context%icontxt, i_part%c_context%iam, i_part%c_context%np )
!            call psb_info ( i_part%g_context%icontxt, i_part%g_context%iam, i_part%g_context%np )
!            call psb_info ( i_part%b_context%icontxt, i_part%b_context%iam, i_part%b_context%np )

!            if(i<=num_levels-2) then
!               allocate(i_part%c_part)
!               i_part => i_part%c_part
!            end if
!        end do
!     end if

!   end subroutine par_partition_map

  !--------------------------------------------------------------------------------------------
  subroutine par_partition_create_standard_execution_model ( p_context, p_part )
    implicit none 
    ! Parameters
    type(par_context)  , target , intent(in)  :: p_context
    type(par_partition)         , intent(out) :: p_part

    assert(p_context%created .eqv. .true.)
    assert(p_context%handler == inhouse .or. p_context%handler == trilinos )
    p_part%p_context => p_context
    nullify(p_part%w_context)
    nullify(p_part%q_context)
    nullify(p_part%b_context)
    nullify(p_part%c_context)
    nullify(p_part%d_context)

  end subroutine par_partition_create_standard_execution_model

  subroutine par_partition_create_multilevel_execution_model ( w_context,  & ! Intracomm including p & q
                                                               p_context,  & ! Fine-tasks intracomm
                                                               q_context,  & ! Rest-of-tasks intracomm
                                                               b_context,  & ! Intercomm p<=>q
                                                               num_levels, &
                                                               id_parts,   &
                                                               num_parts,  &
                                                               p_part)
    implicit none 
    ! Parameters
    type(par_context), target, intent(in)  :: w_context
    type(par_context), target, intent(in)  :: p_context
    type(par_context), target, intent(in)  :: q_context
    type(par_context), target, intent(in)  :: b_context
    integer(ip)              , intent(in)  :: num_levels
    integer(ip)              , intent(in)  :: id_parts(:)
    integer(ip)              , intent(in)  :: num_parts(:)
    type(par_partition)      , intent(out) :: p_part

    assert(p_context%created .eqv. .true.)
    assert(p_context%handler == inhouse .or. p_context%handler == trilinos )
    p_part%p_context => p_context

    assert(w_context%created .eqv. .true.)
    assert(w_context%handler == inhouse .or. w_context%handler == trilinos )
    p_part%w_context => w_context

    assert(q_context%created .eqv. .true.)
    assert(q_context%handler == inhouse .or. q_context%handler == trilinos )
    p_part%q_context => q_context

    assert(b_context%created .eqv. .true.)
    assert(b_context%handler == inhouse .or. b_context%handler == trilinos )
    p_part%b_context => b_context

    assert(size(id_parts) == num_levels)
    assert(size(num_parts) == num_levels)
    p_part%num_levels = num_levels
    call memalloc(p_part%num_levels, p_part%id_parts,__FILE__,__LINE__ )
    call memalloc(p_part%num_levels, p_part%num_parts,__FILE__,__LINE__ )
    p_part%id_parts = id_parts
    p_part%num_parts = num_parts

  end subroutine par_partition_create_multilevel_execution_model


  subroutine par_partition_read ( dir_path, prefix, p_part )
    implicit none 
    ! Parameters
    character *(*)     , intent(in)    :: dir_path
    character *(*)     , intent(in)    :: prefix
    type(par_partition), intent(inout) :: p_part

    ! Locals
    integer         :: iam, num_procs, lunio
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert ( associated(p_part%w_context) )
    assert ( associated(p_part%p_context) )
    assert( p_part%p_context%created .eqv. .true.)
    assert( p_part%p_context%handler == inhouse .or. p_part%p_context%handler == trilinos )

    if(p_part%p_context%iam<0) return

    call par_context_info ( p_part%p_context, iam, num_procs )

    call fem_partition_compose_name ( prefix, name )

    call par_filename( p_part%p_context, name)

    !lunio = io_open (trim(dir_path) // '/' // trim(name), 'read', status='old')

    ! Read fem_partition data from path_file file
    lunio = io_open (trim(dir_path) // '/' // trim(name))
    ! call fem_partition_read ( trim(dir_path) // '/' // trim(name), p_part%f_part ) !old partition
    !write(*,*) 'par_partition_create:'
    call fem_partition_read ( lunio, p_part%f_part )
    !write(*,*) 'par_partition_create: ',p_part%f_part%nparts, num_procs
    call io_close (lunio)

    assert ((p_part%p_context%handler==inhouse.and.p_part%f_part%ptype==element_based).or.(p_part%p_context%handler == trilinos.and.p_part%f_part%ptype == vertex_based))
    !write(*,*) p_part%f_part%nparts, num_procs
    assert ( p_part%f_part%nparts == num_procs )

    if ( p_part%f_part%pinfo == interfaces ) then
       call partition_to_import ( p_part%f_part, p_part%f_import )
    end if

    ! call fem_import_print ( 6, p_part%f_import ) ! DBG:

    ! write(*,*) inhouse, element_based, trilinos, vertex_based, max_ndofs ! DBG: 

    if ( p_part%p_context%handler == trilinos ) then

       do id_map=1, max_ndofs
          p_part%maps_state(id_map) = map_non_created
       end do

    end if

  end subroutine par_partition_read

  subroutine par_partition_create_old ( dir_path, prefix, p_context, p_part, g_context, c_context )
    implicit none 
    ! Parameters
    character *(*)              , intent(in)  :: dir_path
    character *(*)              , intent(in)  :: prefix
    type(par_context)  , target , intent(in)  :: p_context
    type(par_partition)         , intent(out) :: p_part
    type(par_context)  , target , intent(in), optional  :: g_context
    type(par_context)  , target , intent(in), optional  :: c_context

    ! Locals
    integer         :: iam, num_procs, lunio
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id

    assert(p_context%created .eqv. .true.)
    assert(p_context%handler == inhouse .or. p_context%handler == trilinos )

    p_part%p_context => p_context
    if(present(g_context)) then
       p_part%g_context => g_context
    else
       nullify(p_part%g_context)
    end if

    if(present(c_context)) then
       p_part%c_context => c_context
    else
       nullify(p_part%c_context)
    end if

    if(p_context%iam<0) return

    call par_context_info ( p_context, iam, num_procs )

    call fem_partition_compose_name ( prefix, name )

    call par_filename( p_context, name)

    !lunio = io_open (trim(dir_path) // '/' // trim(name), 'read', status='old')

    ! Read fem_partition data from path_file file
    lunio = io_open (trim(dir_path) // '/' // trim(name))
    ! call fem_partition_read ( trim(dir_path) // '/' // trim(name), p_part%f_part ) !old partition
    write(*,*) 'par_partition_create:'
    call fem_partition_read ( lunio, p_part%f_part )
    !write(*,*) 'par_partition_create: ',p_part%f_part%nparts, num_procs
    call io_close (lunio)

    assert ((p_context%handler==inhouse.and.p_part%f_part%ptype==element_based).or.(p_context%handler == trilinos.and.p_part%f_part%ptype == vertex_based))
    !write(*,*) p_part%f_part%nparts, num_procs
    assert ( p_part%f_part%nparts == num_procs )

    if ( p_part%f_part%pinfo == interfaces ) then
       call partition_to_import ( p_part%f_part, p_part%f_import )
    end if

    ! call fem_import_print ( 6, p_part%f_import ) ! DBG:

    ! write(*,*) inhouse, element_based, trilinos, vertex_based, max_ndofs ! DBG: 

    if ( p_context%handler == trilinos ) then

       do id_map=1, max_ndofs
          p_part%maps_state(id_map) = map_non_created
       end do

    end if

  end subroutine par_partition_create_old

  !=============================================================================
  subroutine par_partition_free_one_shot ( p_part )
    implicit none 
    ! Parameters
    type(par_partition), intent(inout) :: p_part
    integer(ip)                        :: id_map
    
    call par_partition_free_progressively ( p_part, free_only_struct )
    call par_partition_free_progressively ( p_part, free_clean )
    
  end subroutine par_partition_free_one_shot

  !=============================================================================
  subroutine par_partition_free_progressively ( p_part, mode )
    implicit none 
    ! Parameters
    type(par_partition), intent(inout) :: p_part
    integer(ip)        , intent(in)    :: mode
    integer(ip)                        :: id_map

    ! Parallel context required
    assert ( associated(p_part%p_context) )
    assert ( p_part%p_context%created .eqv. .true.)
    assert ( mode == free_clean .or. mode == free_only_struct )

    if ( mode == free_clean ) then
       if ( allocated(p_part%id_parts) ) then
          call memfree(p_part%id_parts,__FILE__,__LINE__ )
       end if

       if ( allocated(p_part%num_parts) ) then
          call memfree(p_part%num_parts,__FILE__,__LINE__ )
       end if

       if(p_part%p_context%iam<0) then 
          nullify ( p_part%w_context ) 
          nullify ( p_part%p_context )
          nullify ( p_part%q_context ) 
          nullify ( p_part%c_context ) 
          nullify ( p_part%d_context ) 
          nullify ( p_part%g_context )
          nullify ( p_part%b_context )
          return          
       end  if

       if ( p_part%p_context%handler == trilinos ) then
          do id_map=1, max_ndofs
             if (p_part%maps_state(id_map) == map_created) then
                ! Destruct col-map
                !call epetra_map_destruct   ( p_part%col_map(id_map) )

                ! Destruct row-map
                !call epetra_map_destruct   ( p_part%row_map(id_map) )

                ! Destruct import
                !call epetra_import_destruct ( p_part%importer(id_map) ) 
             end if
             p_part%maps_state(id_map) = map_non_created
          end do
       end if
       nullify ( p_part%w_context ) 
       nullify ( p_part%p_context )
       nullify ( p_part%q_context ) 
       nullify ( p_part%c_context ) 
       nullify ( p_part%d_context ) 
       nullify ( p_part%g_context )
       nullify ( p_part%b_context )

    else if (mode == free_only_struct) then
       if(p_part%p_context%iam<0) then 
          return
       end if

       ! Free fem_partition data
       call fem_partition_free ( p_part%f_part )
       
       if ( p_part%f_part%pinfo == interfaces ) then
          ! Free fem_import data
          call fem_import_free ( p_part%f_import )
       end if
    end if

  end subroutine par_partition_free_progressively

  !=============================================================================
  subroutine par_partition_print(lunou, p_part)
    implicit none
    type(par_partition)  ,  intent(in) :: p_part
    integer(ip)          ,  intent(in) :: lunou
    integer(ip)                        :: id_map

    ! p_part%p_context is required within this subroutine
    assert ( associated(p_part%p_context) )
    assert ( p_part%p_context%created .eqv. .true.)
    if(p_part%p_context%iam<0) return

    write(lunou,'(a)') '*** begin par_partition data structure ***'
    call fem_partition_print (lunou, p_part%f_part)
    if ( p_part%p_context%handler == trilinos ) then 
       do id_map=1, max_ndofs
          if (p_part%maps_state(id_map) == map_created) then
             !call epetra_map_print ( p_part%col_map(id_map) )
             !call epetra_map_print ( p_part%row_map(id_map) )
             !call epetra_import_print ( p_part%importer(id_map) ) 
          end if
       end do
    end if
    write(lunou,'(a)') '*** end par_partition data structure ***'
  end subroutine par_partition_print

end module par_partition_class
