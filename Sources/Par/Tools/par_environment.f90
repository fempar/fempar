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

  type, extends(environment_t) ::  par_environment_t
     private 
     logical                          :: has_been_created = .false.  ! Has the parallel environment been created?
     type (par_context_t)             :: l1_context                  ! 1st lev MPI tasks context
     type (par_context_t)             :: lgt1_context                ! > 1st lev MPI tasks context
     type (par_context_t)             :: l1_lgt1_context             ! Intercommunicator among l1 and lgt1 context
     type (par_context_t)             :: l1_to_l2_context            ! Subcommunicators for l1 to/from l2 data transfers
     
     ! Number of levels in the multilevel hierarchy of MPI tasks
     integer(ip)                      :: num_levels = 0
     integer(ip), allocatable         :: parts_mapping(:), num_parts_per_level(:)
     
     type(par_environment_t), pointer :: next_level 
   contains
     procedure :: create                      => par_environment_create
     procedure :: free                        => par_environment_free
     procedure :: print                       => par_environment_print
     procedure :: created                     => par_environment_created
     procedure :: get_l1_rank                 => par_environment_get_l1_rank 
     procedure :: get_lgt1_rank               => par_environment_get_lgt1_rank
     procedure :: get_l1_context              => par_environment_get_l1_context
     procedure :: get_lgt1_context            => par_environment_get_lgt1_context
     procedure :: am_i_lgt1_task              => par_environment_am_i_lgt1_task
     
     procedure, private :: par_environment_l1_neighbours_exchange_real
     procedure, private :: par_environment_l1_neighbours_exchange_integer
     procedure, private :: par_environment_l1_neighbours_exchange_single_integer     
     generic   :: l1_neighbours_exchange      => par_environment_l1_neighbours_exchange_real, &
                                                 par_environment_l1_neighbours_exchange_integer,&
                                                 par_environment_l1_neighbours_exchange_single_integer
     
     procedure, private :: par_environment_l1_scatter_scalar_integer
     generic   :: l1_scatter => par_environment_l1_scatter_scalar_integer                                                 
     
     procedure, private :: par_environment_l1_gather_scalar_integer
     generic   :: l1_gather => par_environment_l1_gather_scalar_integer 
     
     procedure, private :: par_environment_l1_bcast_scalar_integer
     generic   :: l1_bcast => par_environment_l1_bcast_scalar_integer 
                                                 
     ! Deferred TBPs inherited from class(environment_t)
     procedure :: info                        => par_environment_info
     procedure :: am_i_l1_task                => par_environment_am_i_l1_task
     procedure :: l1_lgt1_bcast               => par_environment_l1_lgt1_bcast
     procedure :: l1_barrier                  => par_environment_l1_barrier
     procedure :: l1_sum_real_scalar          => par_environment_l1_sum_real_scalar
     procedure :: l1_sum_real_vector          => par_environment_l1_sum_real_vector
  end type par_environment_t

  ! Types
  public :: par_environment_t

contains

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
  function par_environment_get_l1_rank ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target,  intent(in) :: this
    integer                                       :: par_environment_get_l1_rank
    par_environment_get_l1_rank = this%l1_context%get_rank()
  end function par_environment_get_l1_rank
  
  !=============================================================================
  function par_environment_get_lgt1_rank ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target,  intent(in) :: this
    integer                                       :: par_environment_get_lgt1_rank
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
    assert ( this%created() )
    par_environment_am_i_lgt1_task = (this%lgt1_context%get_rank() >= 0)
  end function par_environment_am_i_lgt1_task

  subroutine par_environment_l1_barrier(this) 
    implicit none
    ! Dummy arguments
    class(par_environment_t),intent(in)  :: this

    ! Local variables
    integer :: mpi_comm_p, ierr

    ! Parallel environment MUST BE already created
    assert ( this%created() )

    if ( this%am_i_l1_task() ) then
      call psb_get_mpicomm (this%l1_context%get_icontxt(), mpi_comm_p)
      call mpi_barrier ( mpi_comm_p, ierr)
      check ( ierr == mpi_success )
    end if
  end subroutine par_environment_l1_barrier

  subroutine par_environment_info(this,me,np) 
    implicit none
    class(par_environment_t),intent(in)  :: this
    integer(ip)           ,intent(out) :: me
    integer(ip)           ,intent(out) :: np
    assert ( this%created() )
    me = this%l1_context%get_rank()
    np = this%l1_context%get_size()
  end subroutine par_environment_info
  
  function par_environment_am_i_l1_task(this) 
    implicit none
    class(par_environment_t) ,intent(in)  :: this
    logical                             :: par_environment_am_i_l1_task 
    assert ( this%created() )
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
  
  subroutine par_environment_l1_sum_real_scalar (this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    if ( this%am_i_l1_task() ) then
      call psb_sum(this%l1_context%get_icontxt(), alpha)
    end if
  end subroutine par_environment_l1_sum_real_scalar 
     
 subroutine par_environment_l1_sum_real_vector(this,alpha)
    implicit none
    class(par_environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    if ( this%am_i_l1_task() ) then
      call psb_sum(this%l1_context%get_icontxt(), alpha)
    end if
 end subroutine par_environment_l1_sum_real_vector
 
 ! When packing   (gathering) ,    buffer <- alpha * x
 ! When unpacking (scattering),    x <- beta*x + buffer
 subroutine par_environment_l1_neighbours_exchange_real ( this, & 
                                                          num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                                          num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                          alpha, beta, x)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in) :: this
     
   ! Control info to receive
   integer                 , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
   integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

   ! Control info to send
   integer                 , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
   integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

   ! Floating point data
   real(rp), intent(in)    :: alpha, beta
   real(rp), intent(inout) :: x(:)
     
   ! Communication related locals 
   integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
   integer :: mpi_comm,  iret
   integer :: p2pstat(mpi_status_size)
   integer :: icontxt

   ! Arrays required by mpi_all_to_all
   integer, allocatable, dimension(:) :: sndidx, rcvidx, &
                                      &  sndsiz, rcvsiz

   ! Request handlers for non-blocking receives
   integer, allocatable :: rcvhd(:)

   ! Request handlers for non-blocking receives
   integer, allocatable :: sndhd(:)

   real(rp), allocatable :: sndbuf(:) 
   real(rp), allocatable :: rcvbuf(:)
   
   if ( this%am_i_l1_task() ) then
      icontxt   = this%l1_context%get_icontxt()
      my_pid    = this%l1_context%get_rank()
      num_procs = this%l1_context%get_size()

      ! Get MPI communicator associated to icontxt (in
      ! the current implementation of our wrappers
      ! to the MPI library icontxt and mpi_comm are actually 
      ! the same)
      call psb_get_mpicomm (icontxt, mpi_comm)

      call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
      call memalloc (num_snd, sndhd, __FILE__,__LINE__)

      call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
      call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)

      ! Pack send buffers
      call pack_real ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, alpha, x, sndbuf )

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
                 &  psb_double_swap_tag, mpi_comm, rcvhd(i), iret)
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
                 & psb_double_swap_tag, mpi_comm, sndhd(i), iret)
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
      call unpack_real (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, beta, rcvbuf, x )

      call memfree (rcvhd,__FILE__,__LINE__) 
      call memfree (sndhd,__FILE__,__LINE__)

      call memfree (sndbuf,__FILE__,__LINE__)
      call memfree (rcvbuf,__FILE__,__LINE__)
   end if
 end subroutine par_environment_l1_neighbours_exchange_real
 
 ! When packing   (gathering) ,    buffer <- alpha * x
 ! When unpacking (scattering),    x <- beta*x + buffer
 subroutine par_environment_l1_neighbours_exchange_integer ( this, & 
                                                             num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
                                                             num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                             x)
   use psb_const_mod_names
   use psb_penv_mod_names
   implicit none
   class(par_environment_t), intent(in) :: this
     
   ! Control info to receive
   integer                 , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
   integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

   ! Control info to send
   integer                 , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
   integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

   ! Floating point data
   integer(ip), intent(inout) :: x(:)
     
   ! Communication related locals 
   integer :: my_pid, num_procs, i, proc_to_comm, sizmsg
   integer :: mpi_comm,  iret
   integer :: p2pstat(mpi_status_size)
   integer :: icontxt

   ! Arrays required by mpi_all_to_all
   integer, allocatable, dimension(:) :: sndidx, rcvidx, &
                                      &  sndsiz, rcvsiz

   ! Request handlers for non-blocking receives
   integer, allocatable :: rcvhd(:)

   ! Request handlers for non-blocking receives
   integer, allocatable :: sndhd(:)

   integer(ip), allocatable :: sndbuf(:) 
   integer(ip), allocatable :: rcvbuf(:)
   
   if ( this%am_i_l1_task() ) then
      icontxt   = this%l1_context%get_icontxt()
      my_pid    = this%l1_context%get_rank()
      num_procs = this%l1_context%get_size()

      ! Get MPI communicator associated to icontxt (in
      ! the current implementation of our wrappers
      ! to the MPI library icontxt and mpi_comm are actually 
      ! the same)
      call psb_get_mpicomm (icontxt, mpi_comm)

      call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
      call memalloc (num_snd, sndhd, __FILE__,__LINE__)

      call memalloc ((snd_ptrs(num_snd+1)-snd_ptrs(1)), sndbuf, __FILE__,__LINE__)
      call memalloc ((rcv_ptrs(num_rcv+1)-rcv_ptrs(1)), rcvbuf, __FILE__,__LINE__)

      ! Pack send buffers
      call pack_integer ( snd_ptrs(num_snd+1)-snd_ptrs(1), pack_idx, x, sndbuf )
      
      ! First post all the non blocking receives   
      do i=1, num_rcv
         proc_to_comm = list_rcv(i)

         ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
         call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

         ! Message size to be received
         sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

         if ( (sizmsg > 0) .and. (list_rcv(i)-1 /= my_pid) ) then
            call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
                 &  psb_mpi_integer, proc_to_comm, &
                 &  psb_double_swap_tag, mpi_comm, rcvhd(i), iret)
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
                 & psb_mpi_integer, proc_to_comm,    &
                 & psb_double_swap_tag, mpi_comm, sndhd(i), iret)
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
            check ( iret == mpi_success )
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
      call unpack_integer (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), unpack_idx, rcvbuf, x )

      call memfree (rcvhd,__FILE__,__LINE__) 
      call memfree (sndhd,__FILE__,__LINE__)

      call memfree (sndbuf,__FILE__,__LINE__)
      call memfree (rcvbuf,__FILE__,__LINE__)
   end if
 end subroutine par_environment_l1_neighbours_exchange_integer
 
 ! When packing   (gathering) ,    buffer <- alpha * x
 ! When unpacking (scattering),    x <- beta*x + buffer
 subroutine par_environment_l1_neighbours_exchange_single_integer ( this, & 
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
   
   if ( this%am_i_l1_task() ) then
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
     
   end if   
     
   
  end subroutine par_environment_l1_neighbours_exchange_single_integer
 
 subroutine pack_real ( neq, pack_idx, alpha, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)   :: neq
     integer (ip), intent(in)   :: pack_idx(:)
     real    (rp), intent(in)   :: alpha
     real    (rp), intent(in)   :: x(*)
     real    (rp), intent(inout):: y(*)

     ! Locals
     integer(ip) :: i, idof, base_l, base_r

     if (alpha == 0.0_rp) then 
        ! do nothing
     else if (alpha == 1.0_rp) then 
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else if (alpha == -1.0_rp) then 
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = x( base_r )
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else  
        do i=1,neq
           base_l = i
           base_r = pack_idx(i)
           y(base_l) = alpha*x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     end if

   end subroutine pack_real

   subroutine unpack_real ( neq, unpack_idx, beta, x, y )
     implicit none

     ! Parameters
     integer(ip), intent(in)    :: neq
     integer(ip), intent(in)    :: unpack_idx(*)
     real(rp)   , intent(in)    :: beta
     real(rp), intent(in)    :: x(*)
     real(rp), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, base_l, base_r

     if (beta == 0.0_rp) then
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = x(base_r)
           ! write (*,*) 'UNPACK', base_l, base_r, y(base_l), x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else if (beta == 1.0_rp) then
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = y(base_l) + x(base_r)
           ! write (*,*) 'UNPACK', base_l, base_r, y(base_l), x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     else
        do i=1,neq
           base_r = i
           base_l = unpack_idx(i)
           y(base_l) = beta*y(base_l) + x(base_r)
           base_l = base_l + 1
           base_r = base_r + 1
        end do
     end if
   end subroutine unpack_real
   
    subroutine pack_integer ( neq, pack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)    :: neq
     integer (ip), intent(in)    :: pack_idx(:)
     integer (ip), intent(in)    :: x(*)
     integer (ip), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, base_l, base_r

     do i=1,neq
       base_l = i
       base_r = pack_idx(i)
       y(base_l) = x(base_r)
       base_l = base_l + 1
       base_r = base_r + 1
     end do
   end subroutine pack_integer
   
   subroutine unpack_integer ( neq, unpack_idx, x, y )
     implicit none

     ! Parameters
     integer (ip), intent(in)    :: neq
     integer (ip), intent(in)    :: unpack_idx(:)
     integer (ip), intent(in)    :: x(*)
     integer (ip), intent(inout) :: y(*)

     ! Locals
     integer(ip) :: i, base_l, base_r
     
     do i=1,neq
       base_r = i
       base_l = unpack_idx(i)
       y(base_l) = x(base_r)
       base_l = base_l + 1
       base_r = base_r + 1
     end do
     
   end subroutine unpack_integer
   
   subroutine par_environment_l1_gather_scalar_integer ( this, root, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: root
      integer(ip)             , intent(in)   :: input_data
      integer(ip)             , intent(out)  :: output_data(*)
      
      integer(ip) :: icontxt
      integer     :: mpi_comm, iret
      
      if ( this%am_i_l1_task() ) then
        icontxt = this%l1_context%get_icontxt()
        call psb_get_mpicomm (icontxt, mpi_comm)
        call mpi_gather( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, mpi_comm, iret)
        check( iret == mpi_success )
      end if
   end subroutine par_environment_l1_gather_scalar_integer
   
   subroutine par_environment_l1_scatter_scalar_integer ( this, root, input_data, output_data )
      implicit none
      class(par_environment_t), intent(in)   :: this
      integer(ip)             , intent(in)   :: root
      integer(ip)             , intent(in)   :: input_data(*)
      integer(ip)             , intent(out)  :: output_data
      
      integer(ip) :: icontxt
      integer     :: mpi_comm, iret
      
      if ( this%am_i_l1_task() ) then
        icontxt = this%l1_context%get_icontxt()
        call psb_get_mpicomm (icontxt, mpi_comm)
        call mpi_scatter( input_data, 1, psb_mpi_integer, output_data, 1, psb_mpi_integer, root, mpi_comm, iret)
        check( iret == mpi_success )
      end if
   end subroutine par_environment_l1_scatter_scalar_integer
   
   subroutine par_environment_l1_bcast_scalar_integer ( this, root, data )
      implicit none
      class(par_environment_t), intent(in)    :: this
      integer(ip)             , intent(in)    :: root
      integer(ip)             , intent(inout) :: data
      
      integer(ip) :: icontxt
      integer     :: mpi_comm, iret
      
      if ( this%am_i_l1_task() ) then
        icontxt = this%l1_context%get_icontxt()
        call psb_bcast ( icontxt, data, root=root)
      end if
   end subroutine par_environment_l1_bcast_scalar_integer
   
   
end module par_environment_names
