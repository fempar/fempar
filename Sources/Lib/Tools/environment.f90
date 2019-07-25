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
module environment_names
  use types_names
  use memor_names
  use stdio_names
  use FPL
  use timer_names
  ! Parallel modules
  use execution_context_names
  use mpi_context_names
  use mpi_omp_context_names
  use serial_context_names
  use environment_parameters_names
  implicit none

# include "debug.i90"
  private

  ! This type manages the assigment of tasks to different levels as well as
  ! the duties assigned to each of them. The array task_aggregation_among_levels gives the
  ! part assigned to this task and the parts of coarser levels it belongs to. 
  ! Different levels in the hierarchy are managed recursively.
  type ::  environment_t
     private 
     integer(ip)                       :: num_levels
     integer(ip)         , allocatable :: num_tasks_x_level(:)
     integer(ip)         , allocatable :: task_aggregation_among_levels(:)
     class(execution_context_t), allocatable :: world_context          ! All tasks context (=l1+lgt1 tasks)
     class(execution_context_t), allocatable :: l1_context             ! 1st lev tasks context
     class(execution_context_t), allocatable :: lgt1_context           ! > 1st lev tasks context
     class(execution_context_t), allocatable :: l1_to_l2_context       ! Mixed l1 and l2 tasks context 
     type(environment_t), pointer :: next_level                        ! Pointer to the next level environment
   contains     
     procedure          :: create                                      => environment_create
     procedure, private :: process_parameters                          => environment_process_parameters
     procedure, private :: setup_world_context                         => environment_setup_world_context
     procedure, private :: setup_world_l1_and_lgt1_contexts_all_levels => environment_setup_world_l1_and_lgt1_contexts_all_levels
     procedure          :: setup_l1_to_l2_context_all_levels           => environment_setup_l1_to_l2_context_all_levels
     procedure, private :: l1_to_l2_context_all_levels_free            => environment_l1_to_l2_context_all_levels_free

     procedure :: free                           => environment_free
     procedure :: print                          => environment_print
     procedure :: report_times                   => environment_report_times
     procedure :: create_l1_timer                => environment_create_l1_timer

     ! Getters
     procedure :: get_num_levels                   => environment_get_num_levels
     procedure :: get_num_tasks                    => environment_get_num_tasks
     procedure :: get_num_tasks_x_level            => environment_get_num_tasks_x_level
     procedure :: get_next_level                   => environment_get_next_level
     procedure :: get_w_context                    => environment_get_w_context
     procedure :: get_l1_context                   => environment_get_l1_context
     procedure :: get_lgt1_context                 => environment_get_lgt1_context
     procedure :: get_l1_rank                      => environment_get_l1_rank
     procedure :: get_l1_size                      => environment_get_l1_size
     procedure :: get_l1_to_l2_rank                => environment_get_l1_to_l2_rank
     procedure :: get_l1_to_l2_size                => environment_get_l1_to_l2_size
     procedure :: get_lgt1_rank                    => environment_get_lgt1_rank
     procedure :: get_l2_rank_l1_rank_is_mapped_to => environment_get_l2_rank_l1_rank_is_mapped_to
     
     ! Who am I?
     procedure :: am_i_lgt1_task                   => environment_am_i_lgt1_task
     procedure :: am_i_l1_to_l2_task               => environment_am_i_l1_to_l2_task
     procedure :: am_i_l1_to_l2_root               => environment_am_i_l1_to_l2_root
     procedure :: am_i_l1_root                     => environment_am_i_l1_root

     procedure, private :: environment_l1_neighbours_exchange_rp
     procedure, private :: environment_l1_neighbours_exchange_wo_alpha_beta_rp
     procedure, private :: environment_l1_neighbours_exchange_wo_alpha_beta_variable_rp
     procedure, private :: environment_l1_neighbours_exchange_ip
     procedure, private :: environment_l1_neighbours_exchange_igp
     procedure, private :: environment_l1_neighbours_exchange_single_ip
     procedure, private :: environment_l1_neighbours_exchange_wo_pack_unpack_ieep
     procedure, private :: environment_l1_neighbours_exchange_wo_unpack_ip
     procedure, private :: environment_l1_neighbours_exchange_variable_igp
     procedure, private :: environment_l1_neighbours_exchange_variable_ip
     generic   :: l1_neighbours_exchange      => environment_l1_neighbours_exchange_rp, &
                                                 environment_l1_neighbours_exchange_wo_alpha_beta_rp, &
                                                 environment_l1_neighbours_exchange_wo_alpha_beta_variable_rp, &
                                                 environment_l1_neighbours_exchange_ip,&
                                                 environment_l1_neighbours_exchange_igp,&
                                                 environment_l1_neighbours_exchange_single_ip, &
                                                 environment_l1_neighbours_exchange_wo_pack_unpack_ieep, &
                                                 environment_l1_neighbours_exchange_wo_unpack_ip, &
                                                 environment_l1_neighbours_exchange_variable_igp, &
                                                 environment_l1_neighbours_exchange_variable_ip

     procedure, private :: environment_l1_scatter_scalar_ip
     procedure, private :: environment_l1_scatter_scalar_igp
     generic   :: l1_scatter => environment_l1_scatter_scalar_ip, &
                                environment_l1_scatter_scalar_igp

     procedure, private :: environment_l1_gather_scalar_ip
     procedure, private :: environment_l1_gather_scalar_igp
     generic   :: l1_gather => environment_l1_gather_scalar_ip, &
                               environment_l1_gather_scalar_igp

     procedure, private :: environment_l1_bcast_scalar_ip
     procedure, private :: environment_l1_bcast_scalar_igp
     generic   :: l1_bcast => environment_l1_bcast_scalar_ip, &
                              environment_l1_bcast_scalar_igp 

     procedure, private :: environment_l2_from_l1_gather_ip
     procedure, private :: environment_l2_from_l1_gather_igp
     procedure, private :: environment_l2_from_l1_gather_ip_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_ip_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_igp_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_rp_1D_array
     procedure, private :: environment_l2_from_l1_gatherv_rp_2D_array
     generic   :: l2_from_l1_gather => environment_l2_from_l1_gather_ip, &
          environment_l2_from_l1_gather_igp, &
          environment_l2_from_l1_gather_ip_1D_array, &
          environment_l2_from_l1_gatherv_ip_1D_array, &
          environment_l2_from_l1_gatherv_igp_1D_array, &
          environment_l2_from_l1_gatherv_rp_1D_array, &
          environment_l2_from_l1_gatherv_rp_2D_array

     procedure, private :: environment_l2_to_l1_scatter_ip
     procedure, private :: environment_l2_to_l1_scatterv_ip_1D_array                                  
     procedure, private :: environment_l2_to_l1_scatterv_rp_1D_array                                  
     generic   :: l2_to_l1_scatter => environment_l2_to_l1_scatter_ip,           &
                                      environment_l2_to_l1_scatterv_ip_1D_array, &
                                      environment_l2_to_l1_scatterv_rp_1D_array


     procedure, private :: environment_l1_to_l2_transfer_ip
     procedure, private :: environment_l1_to_l2_transfer_ip_1D_array
     procedure, private :: environment_l1_to_l2_transfer_rp
     procedure, private :: environment_l1_to_l2_transfer_rp_1D_array
     procedure, private :: environment_l1_to_l2_transfer_logical
     generic  :: l1_to_l2_transfer => environment_l1_to_l2_transfer_ip, &
          environment_l1_to_l2_transfer_ip_1D_array, &
          environment_l1_to_l2_transfer_rp, &
          environment_l1_to_l2_transfer_rp_1D_array, &
          environment_l1_to_l2_transfer_logical

     ! Deferred TBPs inherited from class(environment_t)
     !procedure :: info                        => environment_info
     procedure :: am_i_l1_task                => environment_am_i_l1_task
     procedure :: l1_lgt1_bcast               => environment_l1_lgt1_bcast
     procedure :: l1_barrier                  => environment_l1_barrier
     procedure :: l1_sum_scalar_rp            => environment_l1_sum_scalar_rp
     procedure :: l1_sum_vector_rp            => environment_l1_sum_vector_rp
     procedure :: l1_max_scalar_rp            => environment_l1_max_scalar_rp
     procedure :: l1_max_vector_rp            => environment_l1_max_vector_rp
     procedure :: l1_max_scalar_ip            => environment_l1_max_scalar_ip
     procedure :: l1_sum_scalar_igp           => environment_l1_sum_scalar_igp
     procedure :: l1_sum_vector_igp           => environment_l1_sum_vector_igp
     generic  :: l1_sum                       => l1_sum_scalar_rp, l1_sum_vector_rp, l1_sum_scalar_igp, l1_sum_vector_igp
     generic  :: l1_max                       => l1_max_scalar_rp, l1_max_vector_rp, l1_max_scalar_ip
     
     procedure :: w_sum_scalar_ip            => environment_w_sum_scalar_ip
     procedure :: w_max_scalar_ip            => environment_w_max_scalar_ip
     procedure :: w_sum_scalar_rp            => environment_w_sum_scalar_rp
     procedure :: w_max_scalar_rp            => environment_w_max_scalar_rp
     procedure :: w_sum_vector_rp            => environment_w_sum_vector_rp
     procedure :: w_max_vector_rp            => environment_w_max_vector_rp
     procedure :: w_sum_scalar_igp           => environment_w_sum_scalar_igp
     procedure :: w_max_scalar_igp           => environment_w_max_scalar_igp
     procedure :: w_sum_vector_igp           => environment_w_sum_vector_igp
     procedure :: w_max_vector_igp           => environment_w_max_vector_igp
     generic  :: w_sum                       => w_sum_scalar_ip, w_sum_scalar_rp, w_sum_vector_rp, w_sum_scalar_igp, w_sum_vector_igp
     generic  :: w_max                       => w_max_scalar_ip, w_max_scalar_rp, w_max_vector_rp, w_max_scalar_igp, w_max_vector_igp    

  end type environment_t

  ! Types
  public :: environment_t, environment_compose_name

contains

  !=============================================================================
  subroutine environment_compose_name ( prefix, name ) 
    implicit none
    character (len=*), intent(in)    :: prefix 
    character (len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.env'
  end subroutine environment_compose_name

  subroutine environment_create( this, world_context, parameters)
    implicit none 
    class(environment_t)      , intent(inout) :: this
    class(execution_context_t), intent(in)    :: world_context
    type(ParameterList_t)     , intent(in)    :: parameters 
    integer(ip) :: istat
    call this%free()
    ! Extract num_levels and num_tasks_x_level from parameters
    call this%process_parameters(world_context,parameters)
    call this%setup_world_l1_and_lgt1_contexts_all_levels(world_context)
  end subroutine environment_create
  
  subroutine environment_setup_world_context(this, world_context)
    implicit none
    class(environment_t)      , intent(inout) :: this
    class(execution_context_t), intent(in)    :: world_context
    integer(ip) :: istat
    allocate(this%world_context,mold=world_context,stat=istat); check(istat==0)
    this%world_context = world_context
  end subroutine  environment_setup_world_context
  
  subroutine environment_process_parameters(this, world_context, parameters)
    implicit none
    class(environment_t)      , intent(inout) :: this
    class(execution_context_t), intent(in)    :: world_context
    type(ParameterList_t)     , intent(in)    :: parameters
    
    integer(ip), allocatable :: array_size(:)
    integer(ip) :: istat
    character(len=:), allocatable :: error_message
    
    ! Mandatory
    assert(parameters%isAssignable(environment_num_levels_key, this%num_levels))
    istat = parameters%get(key = environment_num_levels_key, value = this%num_levels)
    assert(istat==0)

    ! Mandatory (array)
    assert(parameters%isPresent(key = environment_num_tasks_x_level_key))
    istat = parameters%GetShape(key = environment_num_tasks_x_level_key, shape = array_size); check(istat==0)
    
    error_message = 'The size of the ' // & 
                    environment_num_tasks_x_level_key // & 
                    ' (' // ch(array_size(1)) // ') ' // &
                    'array must match the one provided to ' // &
                    environment_num_levels_key // & 
                    ' (' // ch(this%num_levels) // ') '
    massert(array_size(1)==this%num_levels, error_message)
    
    call memalloc(array_size(1), this%num_tasks_x_level,__FILE__,__LINE__)
    istat = parameters%get(key = environment_num_tasks_x_level_key, value = this%num_tasks_x_level); check(istat==0)
    
    error_message = 'At present, only a single task & 
                   supported for the last level provided to ' // & 
                   environment_num_tasks_x_level_key
    massert ( this%num_tasks_x_level(this%num_levels) == 1,  error_message)   
    
    error_message = 'The sum of all tasks provided to ' // & 
                   environment_num_tasks_x_level_key // & 
                   ' (' // ch(sum(this%num_tasks_x_level)) // ') ' // & 
                   'must match the total number of tasks that ENVIRONMENT is handling' // & 
                   ' (' // ch(world_context%get_num_tasks()) // ') '
    massert ( world_context%get_num_tasks() == sum(this%num_tasks_x_level), error_message)
  end subroutine environment_process_parameters
  
  !=============================================================================
  recursive subroutine environment_setup_world_l1_and_lgt1_contexts_all_levels(this, world_context)
    implicit none 
    ! Parameters
    class(environment_t)      , intent(inout) :: this
    class(execution_context_t), intent(in)    :: world_context
    integer                             :: my_color
    integer(ip)                         :: istat

    assert ( this%num_levels >= 1 )
    assert ( world_context%get_current_task() >= 0 )
    
    call this%setup_world_context(world_context)

    ! Create this%l1_context and this%lgt1_context by splitting world_context
    call this%world_context%split_by_condition ( this%world_context%get_current_task() < this%num_tasks_x_level(1), &
                                                 this%l1_context, &
                                                 this%lgt1_context )

    if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
       allocate(this%next_level, stat=istat);check(istat == 0)
       this%next_level%num_levels = this%num_levels-1
       call memalloc(this%next_level%num_levels, this%next_level%num_tasks_x_level,__FILE__,__LINE__ )
       this%next_level%num_tasks_x_level(:) = this%num_tasks_x_level(2:)
       call this%next_level%setup_world_l1_and_lgt1_contexts_all_levels(this%lgt1_context)
    else
       nullify(this%next_level)
    end if
  end subroutine environment_setup_world_l1_and_lgt1_contexts_all_levels
  
  ! Assuming that we have 7 MPI tasks in total,  4 2 1 MPI tasks per level, and the following:
  ! * L1 MPI tasks 0,1 (0,1 in WORLD COMMUNICATOR) are aggregated into L2 MPI task 0 (4 in WORLD COMMUNICATOR)
  ! * L1 MPI tasks 2,3 (3,2 in WORLD COMMUNICATOR) are aggregated into L2 MPI task 1 (5 in WORLD COMMUNICATOR)
  ! * L2 MPI tasks 0,1 (4,5 in WORLD COMMUNICATOR) are aggregated into L3 MPI task 0 (6 in WORLD_COMMUNICATOR)
  ! Then, the task_aggregation_among_levels(NUM_LEVELS) array looks like at each MPI task as follows (on the root of the recursion):
  ! ("x" means that the corresponding array entry is not actually consumed by environment_setup_l1_to_l2_context_all_levels)
  ! * MPI task 0 in WORLD_COMMUNICATOR: (x,1,x)
  ! * MPI task 1 in WORLD_COMMUNICATOR: (x,1,x)
  ! * MPI task 2 in WORLD_COMMUNICATOR: (x,2,x)
  ! * MPI task 3 in WORLD_COMMUNICATOR: (x,2,x)
  ! * MPI task 4 in WORLD_COMMUNICATOR: (x,x,1)
  ! * MPI task 5 in WORLD_COMMUNICATOR: (x,x,1)
  ! * MPI task 6 in WORLD_COMMUNICATOR: (x,x,x)
  recursive subroutine environment_setup_l1_to_l2_context_all_levels( this, task_aggregation_among_levels ) 
    implicit none
    ! Parameters
    class(environment_t)      , intent(inout) :: this
    integer(ip)               , intent(in)    :: task_aggregation_among_levels(:)
    integer(ip) :: istat, my_color
    
    call this%l1_to_l2_context_all_levels_free()
    
    assert ( this%num_levels == size(task_aggregation_among_levels) )
        
    call memalloc ( this%num_levels, this%task_aggregation_among_levels, __FILE__, __LINE__ )
    this%task_aggregation_among_levels(:) = task_aggregation_among_levels
    
    ! Create l1_to_l2_context, where inter-level data transfers actually occur
    if ( this%num_levels > 1 ) then
       if(this%l1_context%get_current_task() >= 0) then
          my_color = this%task_aggregation_among_levels(2)
       else if( this%lgt1_context%get_current_task() < this%num_tasks_x_level(2)  ) then
          my_color = this%lgt1_context%get_current_task()+1
       else
          my_color = undefined_color
       end if
       call this%world_context%split_by_color ( my_color, this%l1_to_l2_context )
    else
       allocate(this%l1_to_l2_context,mold=this%world_context,stat=istat);check(istat==0)
       call this%l1_to_l2_context%nullify()
    end if
    
    if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
       assert (associated(this%next_level) )
       call this%next_level%setup_l1_to_l2_context_all_levels(task_aggregation_among_levels(2:))
    end if
  end subroutine environment_setup_l1_to_l2_context_all_levels

  !=============================================================================
  recursive subroutine environment_l1_to_l2_context_all_levels_free ( this ) 
    implicit none 
    ! Parameters
    class(environment_t), intent(inout) :: this
    integer(ip) :: istat
    if ( allocated(this%task_aggregation_among_levels) ) & 
      call memfree ( this%task_aggregation_among_levels, __FILE__, __LINE__ )
    
    ! Free l1_to_l2_context
    if ( allocated(this%l1_to_l2_context) ) then
       call this%l1_to_l2_context%free(finalize=.false.)
       deallocate(this%l1_to_l2_context,stat=istat);check(istat==0)
    end if 
   
    if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
       assert (associated(this%next_level) )
       call this%next_level%l1_to_l2_context_all_levels_free()
    end if  
  end subroutine environment_l1_to_l2_context_all_levels_free 
  
  !=============================================================================
  recursive subroutine environment_free ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(inout) :: this
    integer(ip) :: istat
    if ( allocated( this%world_context ) ) then
       if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
         call this%next_level%free()
         deallocate ( this%next_level, stat = istat ); assert ( istat == 0 )
       end if
       if ( allocated(this%l1_context) ) then
         call this%l1_context%free(finalize=.false.)
         deallocate ( this%l1_context, stat = istat ); assert ( istat == 0 )
       end if  
       if ( allocated(this%lgt1_context) ) then
         call this%lgt1_context%free(finalize=.false.)
         deallocate ( this%lgt1_context, stat = istat ); assert ( istat == 0 )
       end if 
       if ( allocated(this%l1_to_l2_context) ) then
         call this%l1_to_l2_context%free(finalize=.false.)
         deallocate ( this%l1_to_l2_context, stat = istat ); assert ( istat == 0 )
       end if 
       if ( allocated(this%world_context) ) then
         call this%world_context%free(finalize=.false.)
         deallocate ( this%world_context, stat = istat ); assert ( istat == 0 )
       end if
       this%num_levels = 0
       if (allocated(this%task_aggregation_among_levels)) call memfree(this%task_aggregation_among_levels , __FILE__, __LINE__ )
       if (allocated(this%num_tasks_x_level)) call memfree(this%num_tasks_x_level, __FILE__, __LINE__ )
    end if
  end subroutine environment_free

  subroutine environment_report_times ( this, show_header, luout )
    implicit none 
    class(environment_t), intent(inout) :: this
    logical, intent(in), optional      :: show_header 
    integer(ip), intent(in), optional  :: luout
    call this%l1_context%report_times(show_header, luout)
  end subroutine environment_report_times
 
  !=============================================================================
  recursive subroutine environment_print ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer(ip) :: istat
    write(*,*) 'LEVELS: ', this%num_levels, 'l1_context      : ',this%l1_context%get_current_task(), this%l1_context%get_num_tasks()
    write(*,*) 'LEVELS: ', this%num_levels, 'lgt1_context    : ',this%lgt1_context%get_current_task(), this%lgt1_context%get_num_tasks()
    !write(*,*) 'LEVELS: ', this%num_levels, 'l1_lgt1_context : ',this%l1_lgt1_context%get_rank(), this%l1_lgt1_context%get_size()
    write(*,*) 'LEVELS: ', this%num_levels, 'l1_to_l2_context: ',this%l1_to_l2_context%get_current_task(), this%l1_to_l2_context%get_num_tasks()
    if ( this%num_levels > 1 .and. this%lgt1_context%get_current_task() >= 0 ) then
      call this%next_level%print()
    end if
  end subroutine environment_print

  !=============================================================================
  function environment_get_num_levels (this)
    implicit none 
    class(environment_t), intent(in) :: this
    integer(ip) :: environment_get_num_levels
    environment_get_num_levels = this%num_levels 
  end function environment_get_num_levels
  
  !=============================================================================
  function environment_get_num_tasks (this)
    implicit none 
    class(environment_t), intent(in) :: this
    integer(ip) :: ilevel, environment_get_num_tasks
    assert(this%num_levels>0)
    environment_get_num_tasks = 0
    do ilevel=1,this%num_levels
       environment_get_num_tasks = environment_get_num_tasks + this%num_tasks_x_level(ilevel)
    end do
  end function environment_get_num_tasks
  
  !=============================================================================
  function environment_get_num_tasks_x_level (this)
    implicit none 
    class(environment_t), target, intent(in) :: this
    integer(ip), pointer :: environment_get_num_tasks_x_level(:)
    environment_get_num_tasks_x_level => this%num_tasks_x_level 
  end function environment_get_num_tasks_x_level

  !=============================================================================
  function environment_get_next_level( this )
    implicit none 
    ! Parameters
    class(environment_t), target, intent(in) :: this
    type(environment_t)         , pointer    :: environment_get_next_level

    environment_get_next_level => this%next_level
  end function environment_get_next_level
  
  !=============================================================================
  function environment_get_l1_context ( this ) result(l1_context)
    implicit none 
    ! Parameters
    class(environment_t),       target, intent(in) :: this
    class(execution_context_t), pointer            :: l1_context
    l1_context => this%l1_context
  end function environment_get_l1_context
  
  !=============================================================================
  function environment_get_lgt1_context ( this ) result(lgt1_context)
    implicit none 
    ! Parameters
    class(environment_t),       target, intent(in) :: this
    class(execution_context_t), pointer            :: lgt1_context
    lgt1_context => this%lgt1_context
  end function environment_get_lgt1_context

  !=============================================================================
  function environment_get_l1_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_rank
    environment_get_l1_rank = this%l1_context%get_current_task()
  end function environment_get_l1_rank

  !=============================================================================
  function environment_get_l1_size ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_size
    environment_get_l1_size = this%l1_context%get_num_tasks()
  end function environment_get_l1_size

  !=============================================================================
  function environment_get_l1_to_l2_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_to_l2_rank
    environment_get_l1_to_l2_rank = this%l1_to_l2_context%get_current_task()
  end function environment_get_l1_to_l2_rank

  !=============================================================================
  pure function environment_get_l1_to_l2_size ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_l1_to_l2_size
    environment_get_l1_to_l2_size = this%l1_to_l2_context%get_num_tasks()
  end function environment_get_l1_to_l2_size

  !=============================================================================
  function environment_get_lgt1_rank ( this )
    implicit none 
    ! Parameters
    class(environment_t), intent(in) :: this
    integer                              :: environment_get_lgt1_rank
    environment_get_lgt1_rank = this%lgt1_context%get_current_task()
  end function environment_get_lgt1_rank

  !=============================================================================
  function environment_am_i_lgt1_task(this) 
    implicit none
    class(environment_t) ,intent(in)  :: this
    logical                               :: environment_am_i_lgt1_task
    environment_am_i_lgt1_task = (this%lgt1_context%get_current_task() >= 0)
  end function environment_am_i_lgt1_task

  !=============================================================================
  function environment_am_i_l1_to_l2_task(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_to_l2_task 
    environment_am_i_l1_to_l2_task = (this%l1_to_l2_context%get_current_task() >= 0)
  end function environment_am_i_l1_to_l2_task

  !=============================================================================
  function environment_am_i_l1_to_l2_root(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_to_l2_root
    environment_am_i_l1_to_l2_root = (this%l1_to_l2_context%get_current_task() == this%l1_to_l2_context%get_num_tasks()-1)
  end function environment_am_i_l1_to_l2_root

  !=============================================================================
  function environment_am_i_l1_root(this)
    implicit none
    class(environment_t), intent(in) :: this
    logical                              :: environment_am_i_l1_root
    !environment_am_i_l1_root = (this%l1_context%get_current_task() == 0)
    environment_am_i_l1_root = this%l1_context%am_i_root()
  end function environment_am_i_l1_root

  !=============================================================================
  function environment_get_l2_rank_l1_rank_is_mapped_to (this)
    implicit none
    class(environment_t), intent(in) :: this
    integer(ip) :: environment_get_l2_rank_l1_rank_is_mapped_to 
    environment_get_l2_rank_l1_rank_is_mapped_to = this%task_aggregation_among_levels(2) 
  end function environment_get_l2_rank_l1_rank_is_mapped_to

  !=============================================================================
  subroutine environment_l1_barrier(this) 
    implicit none
    class(environment_t),intent(in)  :: this
    call this%l1_context%barrier()
  end subroutine environment_l1_barrier

  !=============================================================================
  subroutine environment_info(this,me,np) 
    implicit none
    class(environment_t),intent(in)  :: this
    integer(ip)           ,intent(out) :: me
    integer(ip)           ,intent(out) :: np
    me = this%l1_context%get_current_task()
    np = this%l1_context%get_num_tasks()
  end subroutine environment_info

  !=============================================================================
  function environment_am_i_l1_task(this) 
    implicit none
    class(environment_t) ,intent(in)  :: this
    logical                             :: environment_am_i_l1_task 
    environment_am_i_l1_task = (this%l1_context%get_current_task() >= 0)
  end function environment_am_i_l1_task

  !=============================================================================
  subroutine environment_l1_lgt1_bcast(this,condition)
    implicit none 
    ! Parameters
    class(environment_t), intent(in)    :: this
    logical                 , intent(inout) :: condition
    call this%world_context%bcast_subcontext(this%l1_context,this%lgt1_context,condition)
  end subroutine environment_l1_lgt1_bcast

  !=============================================================================
  subroutine environment_create_l1_timer(this, message, mode, timer)
    class(environment_t)      , intent(in)    :: this
    character(len=*)          , intent(in)    :: message
    character(len=*), optional, intent(in)    :: mode
    type(timer_t)             , intent(inout) :: timer
    call timer%free()
    call timer%create(this%l1_context, message, mode)
  end subroutine environment_create_l1_timer  
  
  !=============================================================================
  subroutine environment_l1_sum_scalar_rp (this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call this%l1_context%sum_scalar_rp(alpha)
  end subroutine environment_l1_sum_scalar_rp
  
  !=============================================================================
  subroutine environment_w_sum_scalar_ip(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: alpha
    call this%world_context%sum_scalar_ip(alpha)
  end subroutine environment_w_sum_scalar_ip
  
  !=============================================================================
  subroutine environment_w_max_scalar_ip(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: alpha
    call this%world_context%max_scalar_ip(alpha)
  end subroutine environment_w_max_scalar_ip
  
  !=============================================================================
  subroutine environment_w_sum_scalar_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    call this%world_context%sum_scalar_rp(alpha)
  end subroutine environment_w_sum_scalar_rp

  !=============================================================================
  subroutine environment_w_max_scalar_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha
    call this%world_context%max_scalar_rp(alpha)
  end subroutine environment_w_max_scalar_rp

  !=============================================================================
  subroutine environment_l1_sum_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call this%l1_context%sum_vector_rp(alpha)
  end subroutine environment_l1_sum_vector_rp
  
  !=============================================================================
  subroutine environment_w_sum_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    call this%world_context%sum_vector_rp(alpha)
  end subroutine environment_w_sum_vector_rp
  
  !=============================================================================
  subroutine environment_w_max_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)             , intent(inout) :: alpha(:) 
    call this%world_context%max_vector_rp(alpha)
  end subroutine environment_w_max_vector_rp

  !=============================================================================
  subroutine environment_l1_max_scalar_rp (this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha
    assert ( this%am_i_l1_task() )
    call this%l1_context%max_scalar_rp(alpha)
  end subroutine environment_l1_max_scalar_rp

  !=============================================================================
  subroutine environment_l1_max_vector_rp(this,alpha)
    implicit none
    class(environment_t) , intent(in)    :: this
    real(rp)                 , intent(inout) :: alpha(:) 
    assert (this%am_i_l1_task())
    call this%l1_context%max_vector_rp(alpha)
  end subroutine environment_l1_max_vector_rp
  
  !=============================================================================
  subroutine environment_l1_max_scalar_ip (this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(ip)          , intent(inout) :: n
    assert ( this%am_i_l1_task() )
    call this%l1_context%max_scalar_ip(n)
  end subroutine environment_l1_max_scalar_ip
  
  !=============================================================================
  subroutine environment_l1_sum_scalar_igp (this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n
    assert ( this%am_i_l1_task() )
    call this%l1_context%sum_scalar_igp(n)
  end subroutine environment_l1_sum_scalar_igp
  
  !=============================================================================
  subroutine environment_w_sum_scalar_igp (this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n
    call this%world_context%sum_scalar_igp(n)
  end subroutine environment_w_sum_scalar_igp
  
  !=============================================================================
  subroutine environment_w_max_scalar_igp (this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n
    call this%world_context%max_scalar_igp(n)
  end subroutine environment_w_max_scalar_igp
  
  !=============================================================================
  subroutine environment_l1_sum_vector_igp (this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n(:)
    assert ( this%am_i_l1_task() )
    call this%l1_context%sum_vector_igp(n)
  end subroutine environment_l1_sum_vector_igp
  
  !=============================================================================
  subroutine environment_w_sum_vector_igp(this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n(:)
    call this%world_context%sum_vector_igp(n)
  end subroutine environment_w_sum_vector_igp
  
  !=============================================================================
  subroutine environment_w_max_vector_igp(this,n)
    implicit none
    class(environment_t) , intent(in)    :: this
    integer(igp)         , intent(inout) :: n(:)
    call this%world_context%max_vector_igp(n)
  end subroutine environment_w_max_vector_igp  

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_rp ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       alpha, beta, x, y)
    implicit none
    class(environment_t), intent(inout) :: this

    ! Control info to receive
    integer(ip)             , intent(in) :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in) :: unpack_idx (rcv_ptrs(num_rcv+1)-1)

    ! Control info to send
    integer(ip)             , intent(in) :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in) :: pack_idx (snd_ptrs(num_snd+1)-1)

    ! Floating point data
    real(rp), intent(in)    :: alpha, beta
    real(rp), intent(in)    :: x(:)
    real(rp), intent(inout) :: y(:)

    assert (this%am_i_l1_task())
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                      num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                      alpha, beta, x, y)

  end subroutine environment_l1_neighbours_exchange_rp
  
  !=============================================================================
  subroutine environment_l1_neighbours_exchange_wo_alpha_beta_rp ( this, & 
       &                                                           num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                           num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                           x, y, chunk_size)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    real(rp)                , intent(in)    :: x(:)
    real(rp)                , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size  
  
    assert( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x,y,chunk_size)
  end subroutine environment_l1_neighbours_exchange_wo_alpha_beta_rp
  
  !=============================================================================
  subroutine environment_l1_neighbours_exchange_wo_alpha_beta_variable_rp ( this, & 
       &                                                           num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                                           num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                                           x, y, ptr_chunk_size_snd, ptr_chunk_size_rcv)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    real(rp)                , intent(in)    :: x(:)
    real(rp)                , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size_snd(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size_rcv(:)
  
    assert( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x, y, ptr_chunk_size_snd, ptr_chunk_size_rcv )
  end subroutine environment_l1_neighbours_exchange_wo_alpha_beta_variable_rp

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_ip ( this, & 
       &                                             num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       &                                             num_snd, list_snd, snd_ptrs, pack_idx,   &
       &                                              x, y, chunk_size)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size

    assert( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x,y,chunk_size)

  end subroutine environment_l1_neighbours_exchange_ip

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_igp ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       x, y, chunk_size, mask)
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
    integer(igp)  , optional, intent(in)    :: mask

    assert( this%am_i_l1_task() )

    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x, y, chunk_size, mask)

  end subroutine environment_l1_neighbours_exchange_igp

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_single_ip ( this, & 
       num_neighbours, &
       list_neighbours, &
       input_data,&
       output_data)
    implicit none
    class(environment_t), intent(in)    :: this

    integer                 , intent(in)    :: num_neighbours
    integer(ip)             , intent(in)    :: list_neighbours (num_neighbours)
    integer(ip)             , intent(in)    :: input_data
    integer(ip)             , intent(inout) :: output_data(num_neighbours)

    assert  (this%am_i_l1_task())
    call this%l1_context%neighbours_exchange ( num_neighbours, &
         &                                     list_neighbours, &
         &                                     input_data,&
         &                                     output_data)

  end subroutine environment_l1_neighbours_exchange_single_ip

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_wo_pack_unpack_ieep ( this, &
       &                                                              num_neighbours, &
       &                                                              neighbour_ids, &
       &                                                              snd_ptrs, &
       &                                                              snd_buf, & 
       &                                                              rcv_ptrs, &
       &                                                              rcv_buf )
    ! Parameters
    class(environment_t)  , intent(in)    :: this 
    integer(ip)               , intent(in)    :: num_neighbours
    integer(ip)               , intent(in)    :: neighbour_ids(num_neighbours)
    integer(ip)               , intent(in)    :: snd_ptrs(num_neighbours+1)
    integer(ieep)             , intent(in)    :: snd_buf(snd_ptrs(num_neighbours+1)-1)   
    integer(ip)               , intent(in)    :: rcv_ptrs(num_neighbours+1)
    integer(ieep)             , intent(out)   :: rcv_buf(rcv_ptrs(num_neighbours+1)-1)

    assert ( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_neighbours, &
         &                                     neighbour_ids, &
         &                                     snd_ptrs, &
         &                                     snd_buf, & 
         &                                     rcv_ptrs, &
         &                                     rcv_buf )

  end subroutine environment_l1_neighbours_exchange_wo_pack_unpack_ieep

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_wo_unpack_ip ( this, &
                                                               num_rcv, list_rcv, rcv_ptrs, rcv_buf, &
                                                               num_snd, list_snd, snd_ptrs, pack_idx,   &
                                                               x, chunk_size)
    implicit none
    class(environment_t)    , intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(out)   :: rcv_buf(:)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)   , optional, intent(in)    :: chunk_size
    assert ( this%am_i_l1_task() )
    call this%l1_context%neighbours_exchange ( num_rcv, &
         &                                     list_rcv, &
         &                                     rcv_ptrs, &
         &                                     rcv_buf, & 
         &                                     num_snd, &
         &                                     list_snd,& 
         &                                     snd_ptrs,& 
         &                                     pack_idx,& 
         &                                     x, chunk_size )
  end subroutine environment_l1_neighbours_exchange_wo_unpack_ip

  !=============================================================================
  subroutine environment_l1_neighbours_exchange_variable_igp ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       x, y, ptr_chunk_size, mask )
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(igp)            , intent(in)    :: x(:)
    integer(igp)            , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(igp)  , optional, intent(in)    :: mask

    assert( this%am_i_l1_task() )

    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x, y, ptr_chunk_size, mask)
  end subroutine environment_l1_neighbours_exchange_variable_igp
  
  !=============================================================================
  subroutine environment_l1_neighbours_exchange_variable_ip ( this, & 
       num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
       num_snd, list_snd, snd_ptrs, pack_idx,   &
       x, y, ptr_chunk_size, mask )
    implicit none
    class(environment_t), intent(in)    :: this
    ! Control info to receive
    integer(ip)             , intent(in)    :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip)             , intent(in)    :: unpack_idx (rcv_ptrs(num_rcv+1)-1)
    ! Control info to send
    integer(ip)             , intent(in)    :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip)             , intent(in)    :: pack_idx (snd_ptrs(num_snd+1)-1)
    ! Raw data to be exchanged
    integer(ip)             , intent(in)    :: x(:)
    integer(ip)             , intent(inout) :: y(:)
    integer(ip)             , intent(in)    :: ptr_chunk_size(:)
    integer(ip)   , optional, intent(in)    :: mask

    assert( this%am_i_l1_task() )

    call this%l1_context%neighbours_exchange ( num_rcv, list_rcv, rcv_ptrs, unpack_idx, & 
         &                                     num_snd, list_snd, snd_ptrs, pack_idx,   &
         &                                     x, y, ptr_chunk_size, mask)
  end subroutine environment_l1_neighbours_exchange_variable_ip
  
  
  !=============================================================================
  subroutine environment_l1_gather_scalar_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:) ! ( this%l1_context%get_num_tasks())
    assert ( this%am_i_l1_task() )
    call this%l1_context%gather ( input_data, output_data )
  end subroutine environment_l1_gather_scalar_ip
  
  !=============================================================================
  subroutine environment_l1_gather_scalar_igp ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(igp)             , intent(in)   :: input_data
    integer(igp)             , intent(out)  :: output_data(:) ! ( this%l1_context%get_num_tasks())
    assert ( this%am_i_l1_task() )
    call this%l1_context%gather ( input_data, output_data )
  end subroutine environment_l1_gather_scalar_igp
  
  !=============================================================================
  subroutine environment_l1_scatter_scalar_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data(:) ! ( this%l1_context%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data
    assert( this%am_i_l1_task() )
    call this%l1_context%scatter ( input_data, output_data )
  end subroutine environment_l1_scatter_scalar_ip
  
  !=============================================================================
  subroutine environment_l1_scatter_scalar_igp ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(igp)        , intent(in)   :: input_data(:) ! ( this%l1_context%get_num_tasks())
    integer(igp)        , intent(out)  :: output_data
    assert( this%am_i_l1_task() )
    call this%l1_context%scatter ( input_data, output_data )
  end subroutine environment_l1_scatter_scalar_igp
  
  subroutine environment_l1_bcast_scalar_ip ( this, data )
    implicit none
    class(environment_t), intent(in)    :: this
    integer(ip)             , intent(inout) :: data
    assert ( this%am_i_l1_task() )
    call this%l1_context%bcast ( data )
  end subroutine environment_l1_bcast_scalar_ip
  
  subroutine environment_l1_bcast_scalar_igp ( this, data )
    implicit none
    class(environment_t), intent(in)    :: this
    integer(igp)             , intent(inout) :: data
    assert ( this%am_i_l1_task() )
    call this%l1_context%bcast ( data )
  end subroutine environment_l1_bcast_scalar_igp

  !=============================================================================
  subroutine environment_l1_to_l2_transfer_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data
    integer(ip)             , intent(inout)   :: output_data
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_ip

  !=============================================================================
  subroutine environment_l1_to_l2_transfer_ip_1D_array ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    integer(ip)             , intent(in)      :: input_data(:)
    integer(ip)             , intent(inout)   :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_ip_1D_array

  !=============================================================================
  subroutine environment_l1_to_l2_transfer_rp ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    real(rp)                , intent(in)      :: input_data
    real(rp)                , intent(inout)   :: output_data
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_rp

  !=============================================================================
  subroutine environment_l1_to_l2_transfer_rp_1D_array ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this  
    real(rp)                , intent(in)      :: input_data(:)
    real(rp)                , intent(inout)   :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_rp_1D_array
  
  !=============================================================================
  subroutine environment_l1_to_l2_transfer_logical ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)      :: this
    logical             , intent(in)      :: input_data
    logical             , intent(inout)   :: output_data
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%root_send_master_rcv(input_data, output_data )
  end subroutine environment_l1_to_l2_transfer_logical
  
  !=============================================================================
  subroutine environment_l2_from_l1_gather_ip ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data
    integer(ip)             , intent(out)  :: output_data(:) ! ( this%l1_to_l2_context%get_num_tasks())
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data,output_data )
  end subroutine environment_l2_from_l1_gather_ip

  !=============================================================================
  subroutine environment_l2_from_l1_gather_igp ( this, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(igp)            , intent(in)   :: input_data
    integer(igp)            , intent(out)  :: output_data(:) ! ( this%l1_to_l2_context%get_num_tasks())
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data,output_data )
  end subroutine environment_l2_from_l1_gather_igp

  !=============================================================================
  subroutine environment_l2_from_l1_gather_ip_1D_array ( this, input_data_size, input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data,output_data )
  end subroutine environment_l2_from_l1_gather_ip_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_ip_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(ip)             , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_ip_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_igp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    integer(igp)            , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(igp)            , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_igp_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_rp_1D_array ( this, input_data_size, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)             , intent(in)   :: input_data_size
    real(rp)                , intent(in)   :: input_data(input_data_size)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data_size, input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_rp_1D_array

  !=============================================================================
  subroutine environment_l2_from_l1_gatherv_rp_2D_array ( this, input_data, recv_counts, displs, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:,:)
    integer(ip)             , intent(in)   :: recv_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    real(rp)                , intent(out)  :: output_data(:)
    assert ( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%gather_to_master(input_data, recv_counts, displs, output_data )
  end subroutine environment_l2_from_l1_gatherv_rp_2D_array
  
  !=============================================================================
  subroutine environment_l2_to_l1_scatter_ip ( this,  input_data, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)            , intent(in)   :: input_data(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)            , intent(out)  :: output_data 
    assert( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%scatter_from_master (input_data, output_data )
  end subroutine environment_l2_to_l1_scatter_ip
  
  !=============================================================================
  subroutine environment_l2_to_l1_scatterv_ip_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    integer(ip)              , intent(in)   :: input_data(:)
    integer(ip)             , intent(in)   :: send_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: output_data_size
    integer(ip)             , intent(out)  :: output_data(output_data_size)
    assert( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%scatter_from_master (input_data, send_counts, displs, output_data_size, output_data )
  end subroutine environment_l2_to_l1_scatterv_ip_1D_array
  
  !=============================================================================
  subroutine environment_l2_to_l1_scatterv_rp_1D_array ( this, input_data, send_counts, displs, output_data_size, output_data )
    implicit none
    class(environment_t), intent(in)   :: this
    real(rp)                , intent(in)   :: input_data(:)
    integer(ip)             , intent(in)   :: send_counts(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: displs(:) ! ( this%l1_to_l2_context%get_num_tasks())
    integer(ip)             , intent(in)   :: output_data_size
    real(rp)                , intent(out)  :: output_data(output_data_size)
    assert( this%am_i_l1_to_l2_task() )
    call this%l1_to_l2_context%scatter_from_master (input_data, send_counts, displs, output_data_size, output_data )
  end subroutine environment_l2_to_l1_scatterv_rp_1D_array
  
 !=============================================================================
  function environment_get_w_context ( this ) result(w_context)
    implicit none 
    ! Parameters
    class(environment_t),       target, intent(in) :: this
    class(execution_context_t), pointer            :: w_context
    w_context => this%world_context
  end function environment_get_w_context
  
end module environment_names
