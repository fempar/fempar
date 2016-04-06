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
     procedure :: get_l1_context              => par_environment_get_l1_context
     procedure :: get_lgt1_context            => par_environment_get_l1_context
     procedure :: am_i_lgt1_task              => par_environment_am_i_lgt1_task

     ! Deferred TBPs inherited from class(environment_t)
     procedure :: info                        => par_environment_info
     procedure :: am_i_l1_task                => par_environment_am_i_l1_task
     procedure :: bcast                       => par_environment_bcast
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
  function par_environment_get_l1_context ( this )
    implicit none 
    ! Parameters
    class(par_environment_t), target,  intent(in) :: this
    type(par_context_t)     , pointer             :: par_environment_get_l1_context
    par_environment_get_l1_context => this%l1_context
  end function par_environment_get_l1_context
  
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
      check ( ierr == 0 )
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

  subroutine par_environment_bcast(this,condition)
#ifdef MPI_MOD
    use mpi
#endif
    implicit none 
#ifdef MPI_H
    include 'mpif.h'
#endif
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
  end subroutine par_environment_bcast
  
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

end module par_environment_names
