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
module par_environment_names
  ! Serial modules
  use types_names
  use memor_names
  use stdio_names
  use FPL

  ! Abstract modules
  use environment_names

  ! Parallel modules
  use psb_penv_mod_names
  use par_context_names
  
#ifdef MPI_MOD
  use mpi
#endif
  implicit none
#ifdef MPI_H
  include 'mpif.h'
#endif
  
# include "debug.i90"
  private

  ! This type manages the assigment of tasks to different levels as well as
  ! the dutties assigned to each of them. The array parts_mapping gives the
  ! part assigned to this task and the parts of coarser levels it belongs to. 
  ! Different levels in the hierarchy are managed recursively.
  type, extends(environment_t) ::  par_environment_t
     !private 
     logical                          :: has_been_created = .false.  ! Has the parallel environment been created?
     type (par_context_t)             :: l1_context                  ! 1st lev MPI tasks context
     type (par_context_t)             :: lgt1_context                ! > 1st lev MPI tasks context
     type (par_context_t)             :: l1_lgt1_context             ! Intercommunicator among l1 and lgt1 context
     type (par_context_t)             :: l1_to_l2_context            ! Subcommunicators for l1 to/from l2 data transfers
     
     ! Number of levels in the multilevel hierarchy of tasks
     integer(ip)                      :: task = 0 ! A unique global ID for the current task, independent of the execution context
     integer(ip)                      :: num_levels = 0
     integer(ip), allocatable         :: parts_mapping(:), num_parts_per_level(:)
     
     type(par_environment_t), pointer :: next_level 
   contains
     procedure :: read                                => par_environment_read_file
     procedure :: write                               => par_environment_write_file
     procedure :: create                              => par_environment_create
     procedure :: free                                => par_environment_free
     procedure :: print                               => par_environment_print
     procedure :: created                             => par_environment_created
     procedure :: get_next_level                      => par_environment_get_next_level
     procedure :: get_l1_rank                         => par_environment_get_l1_rank
     procedure :: get_l1_size                         => par_environment_get_l1_size
     procedure :: get_l1_to_l2_rank                   => par_environment_get_l1_to_l2_rank
     procedure :: get_l1_to_l2_size                   => par_environment_get_l1_to_l2_size
     procedure :: get_lgt1_rank                       => par_environment_get_lgt1_rank
     procedure :: get_l1_context                      => par_environment_get_l1_context
     procedure :: get_lgt1_context                    => par_environment_get_lgt1_context
     procedure :: am_i_lgt1_task                      => par_environment_am_i_lgt1_task
     procedure :: am_i_l1_to_l2_task                  => par_environment_am_i_l1_to_l2_task
     procedure :: am_i_l1_to_l2_root                  => par_environment_am_i_l1_to_l2_root
     
     procedure :: get_l2_part_id_l1_task_is_mapped_to => par_environment_get_l2_part_id_l1_task_is_mapped_to
     
     procedure, private :: par_environment_l1_neighbours_exchange_rp
     procedure, private :: par_environment_l1_neighbours_exchange_ip
     procedure, private :: par_environment_l1_neighbours_exchange_igp
     procedure, private :: par_environment_l1_neighbours_exchange_single_ip
     procedure, private :: par_environment_l1_neighbours_exchange_wo_pack_unpack_ieep
     generic   :: l1_neighbours_exchange      => par_environment_l1_neighbours_exchange_rp, &
                                                 par_environment_l1_neighbours_exchange_ip,&
                                                 par_environment_l1_neighbours_exchange_igp,&
                                                 par_environment_l1_neighbours_exchange_single_ip, &
                                                 par_environment_l1_neighbours_exchange_wo_pack_unpack_ieep
     
     procedure, private :: par_environment_l1_scatter_scalar_ip
     generic   :: l1_scatter => par_environment_l1_scatter_scalar_ip                                                 
     
     procedure, private :: par_environment_l1_gather_scalar_ip
     generic   :: l1_gather => par_environment_l1_gather_scalar_ip 
     
     procedure, private :: par_environment_l2_from_l1_gather_ip
     procedure, private :: par_environment_l2_from_l1_gather_igp
     procedure, private :: par_environment_l2_from_l1_gather_ip_1D_array
     procedure, private :: par_environment_l2_from_l1_gatherv_ip_1D_array
     procedure, private :: par_environment_l2_from_l1_gatherv_igp_1D_array
     procedure, private :: par_environment_l2_from_l1_gatherv_rp_1D_array
     procedure, private :: par_environment_l2_from_l1_gatherv_rp_2D_array
     generic   :: l2_from_l1_gather => par_environment_l2_from_l1_gather_ip, &
                                       par_environment_l2_from_l1_gather_igp, &
                                       par_environment_l2_from_l1_gather_ip_1D_array, &
                                       par_environment_l2_from_l1_gatherv_ip_1D_array, &
                                       par_environment_l2_from_l1_gatherv_igp_1D_array, &
                                       par_environment_l2_from_l1_gatherv_rp_1D_array, &
                                       par_environment_l2_from_l1_gatherv_rp_2D_array
                                       
     procedure, private :: par_environment_l2_to_l1_scatterv_rp_1D_array                                  
     generic   :: l2_to_l1_scatter => par_environment_l2_to_l1_scatterv_rp_1D_array
                                       
                                       
     procedure, private :: par_environment_l1_to_l2_transfer_ip
     procedure, private :: par_environment_l1_to_l2_transfer_ip_1D_array
     generic  :: l1_to_l2_transfer => par_environment_l1_to_l2_transfer_ip, &
                                      par_environment_l1_to_l2_transfer_ip_1D_array             
     
     procedure, private :: par_environment_l1_bcast_scalar_ip
     generic   :: l1_bcast => par_environment_l1_bcast_scalar_ip 
     
                                                 
     ! Deferred TBPs inherited from class(environment_t)
     procedure :: info                        => par_environment_info
     procedure :: am_i_l1_task                => par_environment_am_i_l1_task
     procedure :: l1_lgt1_bcast               => par_environment_l1_lgt1_bcast
     procedure :: l1_barrier                  => par_environment_l1_barrier
     procedure :: l1_sum_scalar_rp            => par_environment_l1_sum_scalar_rp
     procedure :: l1_sum_vector_rp            => par_environment_l1_sum_vector_rp
     procedure :: l1_max_scalar_rp            => par_environment_l1_max_scalar_rp
     procedure :: l1_max_vector_rp            => par_environment_l1_max_vector_rp
  end type par_environment_t

  ! Types
  public :: par_environment_t, par_environment_write_files

contains

  !=============================================================================
  subroutine par_environment_compose_name ( prefix, name ) 
    implicit none
    character (len=*), intent(in)    :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.env'
  end subroutine par_environment_compose_name

  !=============================================================================
  subroutine par_environment_write_files ( parameter_list, envs )
    implicit none
    ! Parameters
    type(ParameterList_t)  , intent(in) :: parameter_list
    type(par_environment_t), intent(in)  :: envs(:)

    ! Locals
    integer(ip)          :: nenvs
    integer(ip)          :: istat
    logical              :: is_present
    character(len=256)   :: dir_path
    character(len=256)   :: prefix
    character(len=:), allocatable :: name, rename
    integer(ip)          :: lunio
    integer(ip)          :: i

    nenvs = size(envs)

    ! Mandatory parameters
    is_present = .true.
    is_present =  is_present.and. parameter_list%isPresent(key = dir_path_out_key)
    is_present =  is_present.and. parameter_list%isPresent(key = prefix_key)
    assert(is_present)
     
    istat = 0
    istat = istat + parameter_list%get(key = dir_path_out_key, value = dir_path)
    istat = istat + parameter_list%get(key = prefix_key  , value = prefix)
    check(istat==0)

    call par_environment_compose_name ( prefix, name )
    
    do i=1,nenvs
       rename=name
       call numbered_filename_compose(i,nenvs,rename)
       lunio = io_open (trim(dir_path) // '/' // trim(rename))
       call envs(i)%write (lunio)
       call io_close (lunio)
    end do

    ! name, and rename should be automatically deallocated by the compiler when they
    ! go out of scope. Should we deallocate them explicitly for safety reasons?
  end subroutine  par_environment_write_files
  
  !=============================================================================
  subroutine par_environment_read_file ( this,lunio)
    implicit none 
    class(par_environment_t), intent(inout) :: this
    integer(ip)             , intent(in)    :: lunio

    read ( lunio, '(10i10)' ) this%task
    read ( lunio, '(10i10)' ) this%num_levels
    call memalloc ( this%num_levels, this%num_parts_per_level,__FILE__,__LINE__  )
    call memalloc ( this%num_levels, this%parts_mapping,__FILE__,__LINE__  )
    read ( lunio, '(10i10)' ) this%num_parts_per_level
    read ( lunio, '(10i10)' ) this%parts_mapping
    
  end subroutine par_environment_read_file

  !=============================================================================
  subroutine par_environment_write_file ( this,lunio)
    implicit none 
    class(par_environment_t), intent(in) :: this
    integer(ip)             , intent(in) :: lunio
    write ( lunio, '(10i10)' ) this%task
    write ( lunio, '(10i10)' ) this%num_levels
    write ( lunio, '(10i10)' ) this%num_parts_per_level
    write ( lunio, '(10i10)' ) this%parts_mapping    
  end subroutine par_environment_write_file

  !=============================================================================
  recursive subroutine par_environment_create ( this, world_context, num_levels, num_parts_per_level, parts_mapping)
    implicit none 
    ! Parameters
    class(par_environment_t)   , intent(inout) :: this
    type(par_context_t)        , intent(in)    :: world_context
    integer(ip)                , intent(in)    :: num_levels
    integer(ip)                , intent(in)    :: num_parts_per_level(num_levels)
    integer(ip)                , intent(in)    :: parts_mapping(num_levels)
    integer                                    :: my_color
    integer(ip)                                :: istat
    
    assert ( num_levels >= 1 )
    assert ( world_context%get_rank() >= 0 )
    
    call this%free()
    
    this%num_levels = num_levels
    call memalloc(this%num_levels, this%parts_mapping,__FILE__,__LINE__ )
    call memalloc(this%num_levels, this%num_parts_per_level,__FILE__,__LINE__ )
    this%parts_mapping = parts_mapping
    this%num_parts_per_level = num_parts_per_level
    
    ! Create this%l1_context and this%lgt1_context by splitting world_context
    call world_context%split ( world_context%get_rank() < this%num_parts_per_level(1), this%l1_context, this%lgt1_context )

    ! Create l1_to_l2_context, where inter-level data transfers actually occur
    if ( this%num_levels > 1 ) then
      if(this%l1_context%get_rank() >= 0) then
         my_color = this%parts_mapping(2)
      else if( this%lgt1_context%get_rank() < this%num_parts_per_level(2)  ) then
         my_color = this%lgt1_context%get_rank()+1
      else
         my_color = mpi_undefined
      end if
      call world_context%split ( my_color, this%l1_to_l2_context )
    else
      call this%l1_to_l2_context%nullify()
    end if
    
    ! Create l1_lgt1_context as an intercommunicator among l1_context <=> lgt1_context 
    if ( this%num_levels > 1 ) then
      call this%l1_lgt1_context%create ( world_context, this%l1_context, this%lgt1_context )
    else
      call this%l1_lgt1_context%nullify()
    end if

    if ( this%num_levels > 1 .and. this%lgt1_context%get_rank() >= 0 ) then
      allocate(this%next_level, stat=istat)
      check(istat == 0)
      call this%next_level%create( this%lgt1_context, &
                                   this%num_levels-1, &
                                   this%num_parts_per_level(2:), &
                                   this%parts_mapping(2:) )
    else
      nullify(this%next_level)
    end if
    
    this%has_been_created = .true. 
  end subroutine par_environment_create

  !=============================================================================
  recursive subroutine par_environment_free ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(inout) :: this
    integer(ip)                             :: istat

    if (this%has_been_created) then
       if ( this%num_levels > 1 .and. this%lgt1_context%get_rank() >= 0 ) then
         call this%next_level%free()
         deallocate ( this%next_level, stat = istat )
         assert ( istat == 0 )
       end if       
       this%num_levels = 0
       call memfree(this%parts_mapping , __FILE__, __LINE__ )
       call memfree(this%num_parts_per_level, __FILE__, __LINE__ )
       call this%l1_context%free(finalize=.false.)
       call this%lgt1_context%free(finalize=.false.)
       call this%l1_lgt1_context%free(finalize=.false.)
       call this%l1_to_l2_context%free(finalize=.false.)
       this%has_been_created = .false.
    end if
  end subroutine par_environment_free
  
  !=============================================================================
  recursive subroutine par_environment_print ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer(ip)                          :: istat

    if (this%has_been_created) then
      write(*,*) 'LEVELS: ', this%num_levels, 'l1_context      : ',this%l1_context%get_rank(), this%l1_context%get_size()
      write(*,*) 'LEVELS: ', this%num_levels, 'lgt1_context    : ',this%lgt1_context%get_rank(), this%lgt1_context%get_size()
      write(*,*) 'LEVELS: ', this%num_levels, 'l1_lgt1_context : ',this%l1_lgt1_context%get_rank(), this%l1_lgt1_context%get_size()
      write(*,*) 'LEVELS: ', this%num_levels, 'l1_to_l2_context: ',this%l1_to_l2_context%get_rank(), this%l1_to_l2_context%get_size()
   
      if ( this%num_levels > 1 .and. this%lgt1_context%get_rank() >= 0 ) then
        call this%next_level%print()
      end if
    end if
  end subroutine par_environment_print

  !=============================================================================
  function par_environment_created ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    logical                              :: par_environment_created
    par_environment_created =  this%has_been_created 
  end function par_environment_created
  
  !=============================================================================
  function par_environment_get_next_level( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target, intent(in) :: this
    type(par_environment_t)         , pointer    :: par_environment_get_next_level
    
    par_environment_get_next_level => this%next_level
  end function par_environment_get_next_level
  
  !=============================================================================
  function par_environment_get_l1_rank ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer                              :: par_environment_get_l1_rank
    par_environment_get_l1_rank = this%l1_context%get_rank()
  end function par_environment_get_l1_rank
  
  !=============================================================================
  function par_environment_get_l1_size ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer                              :: par_environment_get_l1_size
    par_environment_get_l1_size = this%l1_context%get_size()
  end function par_environment_get_l1_size
  
    !=============================================================================
  function par_environment_get_l1_to_l2_rank ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer                              :: par_environment_get_l1_to_l2_rank
    par_environment_get_l1_to_l2_rank = this%l1_to_l2_context%get_rank()
  end function par_environment_get_l1_to_l2_rank
  
  !=============================================================================
  pure function par_environment_get_l1_to_l2_size ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer                              :: par_environment_get_l1_to_l2_size
    par_environment_get_l1_to_l2_size = this%l1_to_l2_context%get_size()
  end function par_environment_get_l1_to_l2_size
  
  !=============================================================================
  function par_environment_get_lgt1_rank ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in) :: this
    integer                              :: par_environment_get_lgt1_rank
    par_environment_get_lgt1_rank = this%lgt1_context%get_rank()
  end function par_environment_get_lgt1_rank
  
  !=============================================================================
  function par_environment_get_l1_context ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target,  intent(in) :: this
    type(par_context_t)     , pointer             :: par_environment_get_l1_context
    par_environment_get_l1_context => this%l1_context
  end function par_environment_get_l1_context
  
  !=============================================================================
  function par_environment_get_lgt1_context ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target,  intent(in) :: this
    type(par_context_t)     , pointer             :: par_environment_get_lgt1_context
    par_environment_get_lgt1_context => this%lgt1_context
  end function par_environment_get_lgt1_context
  
  !=============================================================================
  function par_environment_am_i_lgt1_task(this) 
    implicit none
    class(par_environment_t) ,intent(in)  :: this
    logical                               :: par_environment_am_i_lgt1_task
    par_environment_am_i_lgt1_task = (this%lgt1_context%get_rank() >= 0)
  end function par_environment_am_i_lgt1_task
  
  function par_environment_am_i_l1_to_l2_task(this)
    implicit none
    class(par_environment_t), intent(in) :: this
    logical                              :: par_environment_am_i_l1_to_l2_task 
    par_environment_am_i_l1_to_l2_task = (this%l1_to_l2_context%get_rank() >= 0)
  end function par_environment_am_i_l1_to_l2_task
  
  function par_environment_am_i_l1_to_l2_root(this)
    implicit none
    class(par_environment_t), intent(in) :: this
    logical                              :: par_environment_am_i_l1_to_l2_root
    par_environment_am_i_l1_to_l2_root = (this%l1_to_l2_context%get_rank() == this%l1_to_l2_context%get_size()-1)
  end function par_environment_am_i_l1_to_l2_root
  
  !=============================================================================
  function par_environment_get_l2_part_id_l1_task_is_mapped_to (this)
    implicit none
    class(par_environment_t), intent(in) :: this
    integer(ip) :: par_environment_get_l2_part_id_l1_task_is_mapped_to 
    par_environment_get_l2_part_id_l1_task_is_mapped_to = this%parts_mapping(2) 
  end function par_environment_get_l2_part_id_l1_task_is_mapped_to 

  subroutine par_environment_l1_barrier(this) 
    implicit none
    ! Dummy arguments
    class(par_environment_t),intent(in)  :: this
    ! Local variables
    integer :: mpi_comm_p, ierr
    assert ( this%am_i_l1_task() )
    call psb_get_mpicomm (this%l1_context%get_icontxt(), mpi_comm_p)
    call mpi_barrier ( mpi_comm_p, ierr)
    check ( ierr == mpi_success )
  end subroutine par_environment_l1_barrier

  subroutine par_environment_info(this,me,np) 
    implicit none
    class(par_environment_t),intent(in)  :: this
    integer(ip)           ,intent(out) :: me
    integer(ip)           ,intent(out) :: np
    me = this%l1_context%get_rank()
    np = this%l1_context%get_size()
  end subroutine par_environment_info
  
  function par_environment_am_i_l1_task(this) 
    implicit none
    class(par_environment_t) ,intent(in)  :: this
    logical                             :: par_environment_am_i_l1_task 
    par_environment_am_i_l1_task = (this%l1_context%get_rank() >= 0)
  end function par_environment_am_i_l1_task

  subroutine par_environment_l1_lgt1_bcast(this,condition)
    implicit none 
    ! Parameters
    class(par_environment_t), intent(in)    :: this
    logical                 , intent(inout) :: condition
    ! Locals
    integer :: mpi_comm_b, info
    assert ( this%created() )
    if ( this%num_levels > 1 ) then
       ! b_context is an intercomm among p_context & q_context
       ! Therefore the semantics of the mpi_bcast subroutine slightly changes
       ! P_0 in p_context is responsible for bcasting condition to all the processes
       ! in q_context

       ! Get MPI communicator associated to icontxt_b (in
       ! the current implementation of our wrappers
       ! to the MPI library icontxt and mpi_comm are actually 
       ! the same)
       call psb_get_mpicomm (this%l1_lgt1_context%get_icontxt(), mpi_comm_b)
       
       if (this%l1_context%get_rank() >=0) then
          if ( this%l1_context%get_rank() == psb_root_ ) then
             call mpi_bcast(condition,1,MPI_LOGICAL,MPI_ROOT,mpi_comm_b,info)
             check( info == mpi_success )
          else
             call mpi_bcast(condition,1,MPI_LOGICAL,MPI_PROC_NULL,mpi_comm_b,info)
             check( info == mpi_success )
          end if
       else if (this%lgt1_context%get_rank() >=0) then
          call mpi_bcast(condition,1,MPI_LOGICAL,psb_root_,mpi_comm_b,info)
          check( info == mpi_success )
       end if
    end if
  end subroutine par_environment_l1_lgt1_bcast
  
  subroutine par_environment_l1_sum_scalar_rp (this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call psb_sum(this%l1_context%get_icontxt(), alpha)
  end subroutine par_environment_l1_sum_scalar_rp 
     
 subroutine par_environment_l1_sum_vector_rp(this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call psb_sum(this%l1_context%get_icontxt(), alpha)
 end subroutine par_environment_l1_sum_vector_rp
  
 subroutine par_environment_l1_max_scalar_rp (this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call psb_max(this%l1_context%get_icontxt(), alpha)
  end subroutine par_environment_l1_max_scalar_rp 
     
 subroutine par_environment_l1_max_vector_rp(this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call psb_max(this%l1_context%get_icontxt(), alpha)
 end subroutine par_environment_l1_max_vector_rp
 
 ! When packing   (gathering) ,    buffer <- alpha * x
 ! When unpacking (scattering),    x <- beta*x + buffer
 subroutine par_environment_l1_neighbours_exchange_rp ( this, & 
                                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                          alpha, beta, x)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in) :: this
     
   ! Control info to receive
   integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
   integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

   ! Control info to send
   integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
   integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

   ! Floating point data
   real(rp), intent(in)    :: alpha, beta
   real(rp), intent(inout) :: x(:)
     
   ! Communication related locals 
   integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
   integer :: the_mpi_comm,  iret
   integer :: p2pstat(mpi_status_size)
   integer :: icontxt

   ! Request handlers for non-blocking receives
   integer, allocatable :: rcvhd(:)

   ! Request handlers for non-blocking receives
   integer, allocatable :: sndhd(:)

   real(rp), allocatable :: sndbuf(:) 
   real(rp), allocatable :: rcvbuf(:)
   
   assert (this%am_i_l1_task())
   
   icontxt   = this%l1_context%get_icontxt()
   my_pid    = this%l1_context%get_rank()
   num_procs = this%l1_context%get_size()

   call psb_get_mpicomm (icontxt, the_mpi_comm)

   call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
   call memalloc (num_snd, sndhd, __FILE__,__LINE__)

   call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
   call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)

   ! Pack send buffers
   call pack_rp ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, alpha, x, sndbuf )

   ! First post all the non blocking receives   
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
              &  psb_mpi_real, proc_to_comm, &
              &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Secondly post all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be sent
      sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then 
         call mpi_isend(sndbuf(snd_ptrs(i)), sizmsg, &
              & psb_mpi_real, proc_to_comm,    &
              & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Wait on all non-blocking receives
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_wait(rcvhd(i), p2pstat, iret)

      else if ( list_rcv(i)-1 == my_pid ) then
         if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
            write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                 & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
         end if

         rcvbuf( rcv_ptrs(i):rcv_ptrs(i)+sizmsg-1) = &
              sndbuf( snd_ptrs(i): snd_ptrs(i)+sizmsg-1 )
      end if
   end do

   ! Finally wait on all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then
         call mpi_wait(sndhd(i), p2pstat, iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Unpack recv buffers
   call unpack_rp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, beta, rcvbuf, x )

   call memfree (rcvhd,__FILE__,__LINE__) 
   call memfree (sndhd,__FILE__,__LINE__)

   call memfree (sndbuf,__FILE__,__LINE__)
   call memfree (rcvbuf,__FILE__,__LINE__)
   
 end subroutine par_environment_l1_neighbours_exchange_rp
 
 ! When packing   (gathering) ,    buffer <- alpha * x
 ! When unpacking (scattering),    x <- beta*x + buffer
 subroutine par_environment_l1_neighbours_exchange_ip ( this, & 
                                                        num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                                        num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                        x, &
                                                        chunk_size)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in)    :: this
   ! Control info to receive
   integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
   integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
   ! Control info to send
   integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
   integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
   ! Raw data to be exchanged
   integer(ip)             , intent(inout) :: x(:)
   integer(ip)   , optional, intent(in)    :: chunk_size
     
   ! Communication related locals 
   integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
   integer :: the_mpi_comm,  iret
   integer :: p2pstat(mpi_status_size)
   integer :: icontxt

   ! Request handlers for non-blocking receives
   integer, allocatable :: rcvhd(:)

   ! Request handlers for non-blocking receives
   integer, allocatable :: sndhd(:)

   integer(ip), allocatable :: sndbuf(:) 
   integer(ip), allocatable :: rcvbuf(:)
   
   integer(ip) :: chunk_size_
   
   assert( this%am_i_l1_task() )
   
   if ( present(chunk_size) ) then
     chunk_size_ = chunk_size
   else
     chunk_size_ = 1
   end if
   
   
   icontxt   = this%l1_context%get_icontxt()
   my_pid    = this%l1_context%get_rank()
   num_procs = this%l1_context%get_size()


   call psb_get_mpicomm (icontxt, the_mpi_comm)

   call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
   call memalloc (num_snd, sndhd, __FILE__,__LINE__)

   call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
   call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

   ! Pack send buffers
   call pack_ip ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

   ! First post all the non blocking receives   
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
              &  psb_mpi_integer, proc_to_comm, &
              &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Secondly post all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be sent
      sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then 
         call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
              & psb_mpi_integer, proc_to_comm,    &
              & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Wait on all non-blocking receives
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_wait(rcvhd(i), p2pstat, iret)
         check ( iret == mpi_success )
      else if ( list_rcv(i)-1 == my_pid ) then
         if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
            write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                 & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
            check(.false.)
         end if
         rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
              sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
      end if
   end do

   ! Finally wait on all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_
      
      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then
         call mpi_wait(sndhd(i), p2pstat, iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Unpack recv buffers
   call unpack_ip (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, x )

   call memfree (rcvhd,__FILE__,__LINE__) 
   call memfree (sndhd,__FILE__,__LINE__)

   call memfree (sndbuf,__FILE__,__LINE__)
   call memfree (rcvbuf,__FILE__,__LINE__)
 end subroutine par_environment_l1_neighbours_exchange_ip

 subroutine par_environment_l1_neighbours_exchange_igp ( this, & 
                                                         num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                                         num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                         x, &
                                                         chunk_size)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in)    :: this
   ! Control info to receive
   integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
   integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
   ! Control info to send
   integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
   integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
   ! Raw data to be exchanged
   integer(igp)            , intent(inout) :: x(:)
   integer(ip)   , optional, intent(in)    :: chunk_size
   
     
   ! Communication related locals 
   integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
   integer :: the_mpi_comm,  iret
   integer :: p2pstat(mpi_status_size)
   integer :: icontxt

   ! Request handlers for non-blocking receives
   integer, allocatable :: rcvhd(:)

   ! Request handlers for non-blocking receives
   integer, allocatable :: sndhd(:)

   integer(igp), allocatable :: sndbuf(:) 
   integer(igp), allocatable :: rcvbuf(:)
   
   integer(ip) :: chunk_size_
   
   assert( this%am_i_l1_task() )
   
   
   if ( present(chunk_size) ) then
     chunk_size_ = chunk_size
   else
     chunk_size_ = 1
   end if
  
   
   icontxt   = this%l1_context%get_icontxt()
   my_pid    = this%l1_context%get_rank()
   num_procs = this%l1_context%get_size()


   call psb_get_mpicomm (icontxt, the_mpi_comm)

   call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
   call memalloc (num_snd, sndhd, __FILE__,__LINE__)

   call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1))*chunk_size_, sndbuf, __FILE__,__LINE__)
   call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1))*chunk_size_, rcvbuf, __FILE__,__LINE__)

   ! Pack send buffers
   call pack_igp ( snd_ptrs(num_snd+1)-snd_ptrs(1), chunk_size_, pack_idx, x, sndbuf )

   ! First post all the non blocking receives   
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_irecv(  rcvbuf((rcv_ptrs(i)-1)*chunk_size_+1), sizmsg,        &
              &  psb_mpi_long_integer, proc_to_comm, &
              &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Secondly post all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be sent
      sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then 
         call mpi_isend(sndbuf((snd_ptrs(i)-1)*chunk_size_+1), sizmsg, &
              & psb_mpi_long_integer, proc_to_comm,    &
              & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Wait on all non-blocking receives
   do i=1, num_rcv
      proc_to_comm = list_rcv(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (rcv_ptrs(i+1)-rcv_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
         call mpi_wait(rcvhd(i), p2pstat, iret)
         check ( iret == mpi_success )
      else if ( list_rcv(i)-1 == my_pid ) then
         if ( sizmsg /= (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_ ) then 
            write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                 & sizmsg, snd_ptrs(i+1)-snd_ptrs(i) 
            check(.false.)     
         end if
         rcvbuf( (rcv_ptrs(i)-1)*chunk_size_+1:(rcv_ptrs(i)-1)*chunk_size_+sizmsg) = &
              sndbuf( (snd_ptrs(i)-1)*chunk_size_+1:(snd_ptrs(i)-1)*chunk_size_+sizmsg )
      end if
   end do

   ! Finally wait on all non-blocking sends
   do i=1, num_snd
      proc_to_comm = list_snd(i)

      ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
      call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

      ! Message size to be received
      sizmsg = (snd_ptrs(i+1)-snd_ptrs(i))*chunk_size_

      if ( (sizmsg > 0) .and. (list_snd(i)-1 /= my_pid) ) then
         call mpi_wait(sndhd(i), p2pstat, iret)
         check ( iret == mpi_success )
      end if
   end do

   ! Unpack recv buffers
   call unpack_igp (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), chunk_size_, unpack_idx, rcvbuf, x )

   call memfree (rcvhd,__FILE__,__LINE__) 
   call memfree (sndhd,__FILE__,__LINE__)

   call memfree (sndbuf,__FILE__,__LINE__)
   call memfree (rcvbuf,__FILE__,__LINE__)
 end subroutine par_environment_l1_neighbours_exchange_igp

 subroutine par_environment_l1_neighbours_exchange_single_ip ( this, & 
                                                               num_neighbours, &
                                                               list_neighbours, &
                                                               input_data,&
                                                               output_data)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in) :: this
        
   integer                 , intent(in)    :: num_neighbours
   integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
   integer(ip)             , intent(in)    :: input_data
   integer(ip)             , intent(inout) :: output_data(num_neighbours)
   
   integer(ip), allocatable :: ptrs(:)        ! How much data does the part send/recv to/from each neighbour?
   integer(ip), allocatable :: unpack_idx(:)  ! Where the data received from each neighbour is copied/added 
                                              ! on the local vectors of the part ?
   integer(ip), allocatable :: pack_idx(:)    ! Where is located the data to be sent to 
                                              ! each neighbour on the local vectors of the part ?
   
   integer(ip), allocatable :: buffer(:)  
   integer(ip)              :: i 
   
   assert  (this%am_i_l1_task())
   call memalloc ( num_neighbours+1, ptrs, __FILE__, __LINE__ )
   ptrs(1)=1
   do i=2, num_neighbours+1
      ptrs(i)=ptrs(i-1)+1
   end do

   call memalloc ( ptrs(num_neighbours+1)-1, pack_idx, __FILE__, __LINE__ )
   pack_idx = 1 

   call memalloc ( ptrs(num_neighbours+1)-1, unpack_idx, __FILE__, __LINE__ )
   do i=1, ptrs(num_neighbours+1)-1
      unpack_idx(i) = i + 1
   end do

   call memalloc ( num_neighbours+1, buffer, __FILE__, __LINE__ )
   buffer(1) = input_data

   call this%l1_neighbours_exchange ( num_neighbours,    &
                                      list_neighbours,   &
                                      ptrs,              &
                                      unpack_idx,        &  
                                      num_neighbours,    &
                                      list_neighbours,   &
                                      ptrs,              &
                                      pack_idx,          &
                                      buffer )
   
   output_data = buffer(2:)
   
   call memfree (buffer    , __FILE__, __LINE__ )
   call memfree (pack_idx  , __FILE__, __LINE__ )
   call memfree (unpack_idx, __FILE__, __LINE__ )
   call memfree (ptrs      , __FILE__, __LINE__ )
  end subroutine par_environment_l1_neighbours_exchange_single_ip
 
 subroutine pack_rp ( n, pack_idx, alpha, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)   :: n
     integer (ip), intent(in)   :: pack_idx(n)
     real    (rp), intent(in)   :: alpha
     real    (rp), intent(in)   :: x(*)
     real    (rp), intent(inout):: y(*)

     ! Locals
     integer(ip) :: i

     if (alpha == 0.0_rp) then 
        ! do nothing
     else if (alpha == 1.0_rp) then 
        do i=1,n
           y(i) = x(pack_idx(i))
        end do
     else if (alpha == -1.0_rp) then 
        do i=1,n
           y(i) = x(pack_idx(i))
        end do
     else  
        do i=1,n
           y(i) = alpha*x(pack_idx(i))
        end do
     end if

   end subroutine pack_rp

   subroutine unpack_rp ( n, unpack_idx, beta, x, y )
     implicit none

     ! Parameters
     integer(ip), intent(in)    :: n
     integer(ip), intent(in)    :: unpack_idx(n)
     real(rp)   , intent(in)    :: beta
     real(rp)   , intent(in)    :: x(*)
     real(rp)   , intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i

     if (beta == 0.0_rp) then
        do i=1,n
           y(unpack_idx(i)) = x(i)
        end do
     else if (beta == 1.0_rp) then
        do i=1,n
           y(unpack_idx(i)) = y(unpack_idx(i)) + x(i)
        end do
     else
        do i=1,n
           y(unpack_idx(i)) = beta*y(unpack_idx(i)) + x(i)
        end do
     end if
   end subroutine unpack_rp
   
  subroutine pack_ip ( n, chunk_size, pack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)     :: n
     integer (ip), intent(in)     :: chunk_size
     integer (ip), intent(in)     :: pack_idx(n)
     integer (ip), intent(in)    :: x(*)
     integer (ip), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, j, startx, endx
     integer(ip) :: current
     current=1
     do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
         y(current) = x(j)
         current = current + 1
       end do
     end do
   end subroutine pack_ip
   
   subroutine unpack_ip ( n, chunk_size, unpack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)     :: n
     integer (ip), intent(in)     :: chunk_size
     integer (ip), intent(in)     :: unpack_idx(n)
     integer (ip), intent(in)    :: x(*)
     integer (ip), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, j, starty, endy, current
     current = 1
     do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
         y(j) = x(current)
         current = current + 1
       end do
     end do
   end subroutine unpack_ip
    
  subroutine pack_igp ( n, chunk_size, pack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)     :: n
     integer (ip), intent(in)     :: chunk_size
     integer (ip), intent(in)     :: pack_idx(n)
     integer (igp), intent(in)    :: x(*)
     integer (igp), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, j, startx, endx
     integer(ip) :: current
     current=1
     do i=1,n
       startx = (pack_idx(i)-1)*chunk_size + 1
       endx   = startx + chunk_size - 1
       do j=startx, endx
         y(current) = x(j)
         current = current + 1
       end do
     end do
   end subroutine pack_igp
   
   subroutine unpack_igp ( n, chunk_size, unpack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)     :: n
     integer (ip), intent(in)     :: chunk_size
     integer (ip), intent(in)     :: unpack_idx(n)
     integer (igp), intent(in)    :: x(*)
     integer (igp), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, j, starty, endy, current
     current = 1
     do i=1,n
       starty = (unpack_idx(i)-1)*chunk_size + 1
       endy   = starty + chunk_size - 1
       do j=starty, endy
         y(j) = x(current)
         current = current + 1
       end do
     end do
   end subroutine unpack_igp
   
   subroutine par_environment_l1_neighbours_exchange_wo_pack_unpack_ieep ( this, &
                                                                           number_neighbours, &
                                                                           neighbour_ids, &
                                                                           snd_ptrs, &
                                                                           snd_buf, & 
                                                                           rcv_ptrs, &
                                                                           rcv_buf )
    use psb_const_mod_names
    use psb_penv_mod_names
#ifdef MPI_MOD
    use mpi
#endif
     implicit none
#ifdef MPI_H
     include 'mpif.h'
#endif
     ! Parameters
     class(par_environment_t)  , intent(in)    :: this 
     integer(ip)               , intent(in)    :: number_neighbours
     integer(ip)               , intent(in)    :: neighbour_ids(number_neighbours)
     integer(ip)               , intent(in)    :: snd_ptrs(number_neighbours+1)
     integer(ieep)             , intent(in)    :: snd_buf(snd_ptrs(number_neighbours+1)-1)   
     integer(ip)               , intent(in)    :: rcv_ptrs(number_neighbours+1)
     integer(ieep)             , intent(out)   :: rcv_buf(rcv_ptrs(number_neighbours+1)-1)
     
     ! Communication related locals 
     integer(ip) :: icontxt 
     integer     :: my_pid, proc_to_comm, sizmsg
     integer     :: the_mpi_comm,  iret
     integer     :: p2pstat(mpi_status_size)
     
     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: rcvhd

     ! Request handlers for non-blocking receives
     integer, allocatable, dimension(:) :: sndhd
     
     type(par_context_t), pointer :: l1_context
     integer(ip) :: i
     
     assert ( this%am_i_l1_task() )
     
     l1_context => this%get_l1_context()
     icontxt    = l1_context%get_icontxt()
     my_pid     = l1_context%get_rank()
     
     call memalloc (number_neighbours, rcvhd, __FILE__,__LINE__)
     call memalloc (number_neighbours, sndhd, __FILE__,__LINE__)
     
     
     call psb_get_mpicomm (icontxt, the_mpi_comm)
       
     ! First post all the non blocking receives   
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_irecv(  rcv_buf(rcv_ptrs(i)), sizmsg, &
                        &  psb_mpi_integer1, proc_to_comm, &
                        &  psb_double_swap_tag, the_mpi_comm, rcvhd(i), iret)
          check ( iret == mpi_success )
       end if
     end do

     ! Secondly post all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be sent
        sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
    
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then 
             call mpi_isend(snd_buf(snd_ptrs(i)), sizmsg, &
                     & psb_mpi_integer1, proc_to_comm, &
                     & psb_double_swap_tag, the_mpi_comm, sndhd(i), iret)
             check ( iret == mpi_success )
        end if
     end do

     ! Wait on all non-blocking receives
     do i=1, number_neighbours
       proc_to_comm = neighbour_ids(i)
         
       ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
       call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
       ! Message size to be received
       sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
      
       if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(rcvhd(i), p2pstat, iret)
          check (iret == mpi_success)
       else if ( neighbour_ids(i)-1 == my_pid ) then
          if ( sizmsg /= snd_ptrs(i+1)-snd_ptrs(i) ) then 
             write(0,*) 'Fatal error in single_exchange: mismatch on self sendf', & 
                     & sizmsg, snd_ptrs(i+1)-snd_ptrs(i)
          end if
          rcv_buf( rcv_ptrs(i):rcv_ptrs(i+1)-1 ) = &
                 snd_buf( snd_ptrs(i):snd_ptrs(i+1)-1 )
       end if
     end do

     ! Finally wait on all non-blocking sends
     do i=1, number_neighbours
        proc_to_comm = neighbour_ids(i)
          
        ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
        call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
          
        ! Message size to be received
        sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
      
        if ( (sizmsg > 0) .and. (neighbour_ids(i)-1 /= my_pid) ) then
          call mpi_wait(sndhd(i), p2pstat, iret)
          check ( iret == mpi_success )
        end if
     end do
     
     call memfree (rcvhd ,__FILE__,__LINE__) 
     call memfree (sndhd ,__FILE__,__LINE__)
   end subroutine par_environment_l1_neighbours_exchange_wo_pack_unpack_ieep 
   
   
   subroutine par_environment_l1_gather_scalar_ip ( this, root, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer                 , intent(in)   :: root
      integer(ip)             , intent(in)   :: input_data
      integer(ip)             , intent(out)  :: output_data(this%l1_context%get_size())
      integer(ip) :: icontxt
      integer     :: the_mpi_comm, iret
      assert ( this%am_i_l1_task() )
      icontxt = this%l1_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gather( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l1_gather_scalar_ip
   
   subroutine par_environment_l1_scatter_scalar_ip ( this, root, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer                 , intent(in)   :: root
      integer(ip)             , intent(in)   :: input_data(this%l1_context%get_size())
      integer(ip)             , intent(out)  :: output_data
      integer(ip) :: icontxt
      integer     :: the_mpi_comm, iret
      assert( this%am_i_l1_task() )
      icontxt = this%l1_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_scatter( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l1_scatter_scalar_ip
   
   subroutine par_environment_l1_bcast_scalar_ip ( this, root, data )
      implicit none
      class(par_environment_t), intent(in)    :: this
      integer                 , intent(in)    :: root
      integer(ip)             , intent(inout) :: data
      integer(ip) :: icontxt
      assert ( this%am_i_l1_task() )
      icontxt = this%l1_context%get_icontxt()
      call psb_bcast ( icontxt, data, root=root)
   end subroutine par_environment_l1_bcast_scalar_ip
   
   subroutine par_environment_l1_to_l2_transfer_ip ( this, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)      :: this
      integer(ip)             , intent(in)      :: input_data
      integer(ip)             , intent(inout)   :: output_data
      integer(ip) :: icontxt
      integer(ip) :: recv_rank
      integer(ip) :: send_rank
      assert ( this%am_i_l1_to_l2_task() )
      send_rank = 0
      recv_rank = this%l1_to_l2_context%get_size()-1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      if ( this%get_l1_to_l2_rank() == send_rank ) then
        call psb_snd(icontxt, input_data, recv_rank)
      else if (this%get_l1_to_l2_rank() == recv_rank) then
        call psb_rcv(icontxt, output_data, send_rank)
      end if
   end subroutine par_environment_l1_to_l2_transfer_ip
   
   subroutine par_environment_l1_to_l2_transfer_ip_1D_array ( this, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)      :: this
      integer(ip)             , intent(in)      :: input_data(:)
      integer(ip)             , intent(inout)   :: output_data(:)
      integer(ip) :: icontxt
      integer(ip) :: recv_rank
      integer(ip) :: send_rank
      assert ( this%am_i_l1_to_l2_task() )
      send_rank = 0
      recv_rank = this%l1_to_l2_context%get_size()-1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      if ( this%get_l1_to_l2_rank() == send_rank ) then
        call psb_snd(icontxt, input_data, recv_rank)
      else if (this%get_l1_to_l2_rank() == recv_rank) then
        call psb_rcv(icontxt, output_data, send_rank)
      end if
   end subroutine par_environment_l1_to_l2_transfer_ip_1D_array
  
   
   subroutine par_environment_l2_from_l1_gather_ip ( this, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: input_data
      integer(ip)             , intent(out)  :: output_data(this%l1_to_l2_context%get_size())
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      assert ( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gather( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gather_ip

   subroutine par_environment_l2_from_l1_gather_igp ( this, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(igp)            , intent(in)   :: input_data
      integer(igp)            , intent(out)  :: output_data(this%l1_to_l2_context%get_size())
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      assert ( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gather( input_data, 1, psb_mpi_long_integer, output_data, 1, psb_mpi_long_integer, root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gather_igp
   
   subroutine par_environment_l2_from_l1_gather_ip_1D_array ( this, input_data_size, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: input_data_size
      integer(ip)             , intent(in)   :: input_data(input_data_size)
      integer(ip)             , intent(out)  :: output_data(*)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gather( input_data,  input_data_size, psb_mpi_integer, &
                       output_data, input_data_size, psb_mpi_integer, &
                       root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gather_ip_1D_array
   
   subroutine par_environment_l2_from_l1_gatherv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: input_data_size
      integer(ip)             , intent(in)   :: input_data(input_data_size)
      integer(ip)             , intent(in)   :: recv_counts(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: displs(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(out)  :: output_data(*)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gatherv( input_data, input_data_size, psb_mpi_integer, &
           output_data, recv_counts, displs, psb_mpi_integer, &
           root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gatherv_ip_1D_array
   
   subroutine par_environment_l2_from_l1_gatherv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: input_data_size
      integer(igp)            , intent(in)   :: input_data(input_data_size)
      integer(ip)             , intent(in)   :: recv_counts(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: displs(this%l1_to_l2_context%get_size())
      integer(igp)            , intent(out)  :: output_data(*)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gatherv( input_data, input_data_size, psb_mpi_long_integer, &
                        output_data, recv_counts, displs, psb_mpi_long_integer, &
                        root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gatherv_igp_1D_array

   subroutine par_environment_l2_from_l1_gatherv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: input_data_size
      real(rp)                , intent(in)   :: input_data(input_data_size)
      integer(ip)             , intent(in)   :: recv_counts(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: displs(this%l1_to_l2_context%get_size())
      real(rp)                , intent(out)  :: output_data(*)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gatherv( input_data , input_data_size, psb_mpi_real, &
                        output_data, recv_counts, displs, psb_mpi_real, &
                        root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gatherv_rp_1D_array
   
   subroutine par_environment_l2_from_l1_gatherv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      real(rp)                , intent(in)   :: input_data(:,:)
      integer(ip)             , intent(in)   :: recv_counts(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: displs(this%l1_to_l2_context%get_size())
      real(rp)                , intent(out)  :: output_data(*)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_gatherv( input_data , size(input_data,1)*size(input_data,2), psb_mpi_real, &
                        output_data, recv_counts, displs, psb_mpi_real, &
                        root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_from_l1_gatherv_rp_2D_array
   
   
     subroutine par_environment_l2_to_l1_scatterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      real(rp)                , intent(in)   :: input_data(*)
      integer(ip)             , intent(in)   :: send_counts(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: displs(this%l1_to_l2_context%get_size())
      integer(ip)             , intent(in)   :: output_data_size
      real(rp)                , intent(out)  :: output_data(output_data_size)
      
      integer(ip)            :: icontxt
      integer                :: the_mpi_comm, iret, root
      
      assert( this%am_i_l1_to_l2_task() )
      root    = this%l1_to_l2_context%get_size() - 1 
      icontxt = this%l1_to_l2_context%get_icontxt()
      call psb_get_mpicomm (icontxt, the_mpi_comm)
      call mpi_scatterv( input_data, send_counts, displs, psb_mpi_real, &
                         output_data, output_data_size, psb_mpi_real, &
                         root, the_mpi_comm, iret)
      check( iret == mpi_success )
   end subroutine par_environment_l2_to_l1_scatterv_rp_1D_array
   
end module par_environment_names
