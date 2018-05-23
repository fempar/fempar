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
! ==============================================================================================
module hts_theta_method_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
  
  type :: theta_method_t
     private
					type(environment_t), pointer :: environment
     real(rp)                     :: theta
     real(rp)                     :: initial_time     
     real(rp)                     :: current_time
     real(rp)                     :: final_time     
     real(rp)                     :: time_step
     real(rp)                     :: max_time_step 
     real(rp)                     :: min_time_step 
     integer(ip)                  :: save_solution_every_n_steps 
     integer(ip)                  :: current_step      
     real(rp)                     :: next_time_to_be_printed  
     integer(ip)                  :: num_steps
   contains
     procedure, non_overridable :: create                    => theta_method_create
     procedure, non_overridable :: update_solutions          => theta_method_update_solutions
     procedure, non_overridable :: move_time_forward         => theta_method_move_time_forward
     procedure, non_overridable :: move_time_backwards       => theta_method_move_time_backwards
     procedure, non_overridable :: print                     => theta_method_print
     procedure, non_overridable :: print_this_step           => theta_method_print_this_step 
     procedure, non_overridable :: update_time_to_be_printed => theta_method_update_time_to_be_printed   
     procedure, non_overridable :: get_theta                 => theta_method_get_theta
     procedure, non_overridable :: get_initial_time          => theta_method_get_initial_time
     procedure, non_overridable :: get_current_time          => theta_method_get_current_time
     procedure, non_overridable :: get_final_time            => theta_method_get_final_time
     procedure, non_overridable :: get_time_step             => theta_method_get_time_step
     procedure, non_overridable :: finished                  => theta_method_finished
  end type theta_method_t
  
  public :: theta_method_t

contains
  !===============================================================================================
  subroutine theta_method_create(this, environment, theta, initial_time, final_time, num_time_steps, max_time_step, min_time_step, save_every_n_steps )
    implicit none
    class(theta_method_t), intent(inout) :: this
				type(environment_t)  , target, intent(in)     :: environment
    real(rp)             , intent(in)    :: theta
    real(rp)             , intent(in)    :: initial_time
    real(rp)             , intent(in)    :: final_time
    integer(ip)          , intent(in)    :: num_time_steps
    real(rp)             , intent(in)    :: max_time_step 
    real(rp)             , intent(in)    :: min_time_step 
    integer(ip)          , intent(in)    :: save_every_n_steps 
    
				this%environment                 => environment
    this%theta                       = theta
    this%initial_time                = initial_time
    this%final_time                  = final_time
    this%num_steps                   = num_time_steps
    this%time_step                   = ( this%final_time - this%initial_time ) / real(this%num_steps,rp)
    this%current_time                = this%time_step
    this%current_step                = 1
    this%max_time_step               = max_time_step 
    this%min_time_step               = min_time_step 
    this%save_solution_every_n_steps = save_every_n_steps  
    this%next_time_to_be_printed     = final_time/real(save_every_n_steps) 
    
  end subroutine theta_method_create

  !===============================================================================================
  subroutine theta_method_update_solutions(this, current_solution, previous_solution )
    implicit none
    class(theta_method_t), intent(inout) :: this
    type(fe_function_t)  , intent(inout) :: current_solution 
    type(fe_function_t)  , intent(inout) :: previous_solution
    class(vector_t), pointer :: dof_values_current_solution
    class(vector_t), pointer :: dof_values_previous_solution

    
    assert ( .not. this%finished() )
    assert( (this%theta <= 1.0_rp) .and. (this%theta > 0.0_rp) )
    
    dof_values_current_solution => current_solution%get_free_dof_values()
    dof_values_previous_solution => previous_solution%get_free_dof_values()
    
    ! u^{n+1} = ( 1/theta ) * ( u^{n+theta} - (1-theta) * u^{n} ) (BCs not touched)
    dof_values_current_solution  = (1.0_rp/this%theta) * (dof_values_current_solution - (1.0_rp-this%theta)* dof_values_previous_solution)
    
    ! u^{n} <- u^{n+1} (including BCs)
    call previous_solution%copy(current_solution)
  end subroutine theta_method_update_solutions
 
  !===============================================================================================
  subroutine theta_method_move_time_forward(this, num_iterations, ideal_num_iterations )
    implicit none
    class(theta_method_t), intent(inout) :: this
    integer(ip)  , optional , intent(in)    :: num_iterations
    integer(ip)  , optional , intent(in)    :: ideal_num_iterations 
    
    ! Update time step with nonlinear convergence history 
	if ( present(num_iterations) .and. present(ideal_num_iterations) ) then 
    this%time_step = min( real(ideal_num_iterations,rp)/real(num_iterations,rp), 5.0_rp )*this%time_step
    this%time_step = max(this%time_step, this%min_time_step) 
    this%time_step = min(this%time_step, this%max_time_step) 
	end if 
    
    ! Update theta-method scheme 
    this%current_time = this%current_time + this%time_step
    this%current_step = this%current_step + 1
  end subroutine theta_method_move_time_forward
  
    !===============================================================================================
  subroutine theta_method_move_time_backwards(this, current_solution, previous_solution)
    implicit none
    class(theta_method_t), intent(inout) :: this
    type(fe_function_t)  , intent(inout) :: current_solution 
    type(fe_function_t)  , intent(inout) :: previous_solution
    real(rp)                             :: updated_time_step 
    class(vector_t), pointer :: dof_values_current_solution
    class(vector_t), pointer :: dof_values_previous_solution
   
    write(*,*) 'Nonlinear Convergence Did not succeed or aborted (large Residual), reducing time step and cycling'
    
    ! Update time step with nonlinear convergence history 
    updated_time_step = this%time_step/2.0_rp 
    
    ! Theta method backwards with new time step 
    this%current_time = this%current_time - this%time_step + updated_time_step 
    
    ! Assign current time step size the updated one 
    this%time_step = updated_time_step
    
    ! Back to the previous converged solution  
    dof_values_previous_solution => previous_solution%get_free_dof_values()
    dof_values_current_solution =>  current_solution%get_free_dof_values()
    dof_values_current_solution = dof_values_previous_solution 
    
  end subroutine theta_method_move_time_backwards
  
  !===============================================================================================
  subroutine theta_method_print(this,luout)
    implicit none
    class(theta_method_t), intent(in) :: this
    integer(ip)          , intent(in) :: luout
				if (this%environment%get_l1_rank() == 0) then
    write(luout,*) '========================================================================'
    write(luout,'(a10,i6,a20, e10.3,a12,e10.3)') 'Time step ', this%current_step, ': Solving for t=', this%current_time, 'with dt', this%time_step
				end if 
  end subroutine theta_method_print
  
  !=============================================================================================== 
    function theta_method_print_this_step(this)
    implicit none
    class(theta_method_t), intent(inout) :: this
    logical :: theta_method_print_this_step 

    theta_method_print_this_step = ( this%current_time .ge. this%next_time_to_be_printed )
    
  end function theta_method_print_this_step 
  
    !=============================================================================================== 
   subroutine theta_method_update_time_to_be_printed(this)
    implicit none
    class(theta_method_t), intent(inout) :: this
    logical :: theta_method_print_this_step 
     
    this%next_time_to_be_printed = this%next_time_to_be_printed + &
                                   this%final_time/real(this%save_solution_every_n_steps,rp) 
                                   
  end subroutine theta_method_update_time_to_be_printed  
  
  !===============================================================================================
  function theta_method_get_theta(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    real(rp)                          :: theta_method_get_theta
    theta_method_get_theta = this%theta
  end function theta_method_get_theta
  
  !===============================================================================================
  function theta_method_get_current_time(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    real(rp)                          :: theta_method_get_current_time
    theta_method_get_current_time = this%current_time
  end function theta_method_get_current_time
  
  !===============================================================================================
  function theta_method_get_initial_time(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    real(rp)                          :: theta_method_get_initial_time
    theta_method_get_initial_time = this%initial_time
  end function theta_method_get_initial_time
  
  !===============================================================================================
  function theta_method_get_final_time(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    real(rp)                          :: theta_method_get_final_time
    theta_method_get_final_time = this%final_time
  end function theta_method_get_final_time
  
  !===============================================================================================
  function theta_method_get_time_step(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    real(rp)                          :: theta_method_get_time_step
    theta_method_get_time_step = this%time_step
  end function theta_method_get_time_step
    
  !===============================================================================================
  function theta_method_finished(this)
    implicit none
    class(theta_method_t), intent(in) :: this
    logical                           :: theta_method_finished
    theta_method_finished = (this%current_time > this%final_time)
  end function theta_method_finished
  
end module hts_theta_method_names




