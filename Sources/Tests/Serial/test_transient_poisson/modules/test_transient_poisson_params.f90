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
!* title: Test Parameters
module test_transient_poisson_params_names
  use fempar_names
# include "debug.i90"
  implicit none
  private

  character(len=*), parameter :: print_nonlinear_iteration_key = 'print_nonlinear_iteration'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key            = 'write_solution'        
  character(len=*), parameter :: initial_time_key              = 'initial_time'
  character(len=*), parameter :: final_time_key                = 'final_time'
  character(len=*), parameter :: time_step_key                 = 'time_step'
  character(len=*), parameter :: num_time_steps_key            = 'num_time_steps'
  character(len=*), parameter :: time_integration_scheme_key   = 'time_integration_scheme'
  character(len=*), parameter :: is_test_key                   = 'is_test'

  type:: test_transient_poisson_params_t  
   private
     logical :: print_nonlinear_iteration = .false.
   contains
     procedure, non_overridable             :: process_parameters
     procedure, non_overridable             :: get_parameter_list
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_triangulation_type
     procedure, non_overridable             :: get_initial_time
     procedure, non_overridable             :: get_final_time
     procedure, non_overridable             :: get_time_step
     procedure, non_overridable             :: get_time_integration_scheme
     procedure, non_overridable             :: get_is_test
     procedure, non_overridable             :: get_print_nonlinear_iteration
  end type test_transient_poisson_params_t  

  ! Types
  public :: test_transient_poisson_params_t

contains

  !==================================================================================================
  subroutine test_transient_poisson_define_user_parameters()
    implicit none
    ! IO parameters
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order') 
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution') 

    call parameter_handler%add(initial_time_key, '--initial-time', 0.0_rp, 'Initial time: t0', switch_ab='-t0') 
    call parameter_handler%add(final_time_key, '--final-time', 1.0_rp, 'Final time: tf', switch_ab='-tf') 
    call parameter_handler%add(time_step_key, '--time-step', 1.0_rp, 'Time step size: dt', switch_ab='-dt') 
    call parameter_handler%add(num_time_steps_key, '--num-time-steps', 1, 'Maximum number of time steps: nt', switch_ab='-nt') 
    call parameter_handler%add(time_integration_scheme_key, '--time-integration-scheme', 'backward_euler', 'Time disctetization scheme of the DIRK solver.', switch_ab='-rk-scheme') 

    call parameter_handler%add(is_test_key, '--is-test', .false., 'Test convergence order of the runge kutta scheme', switch_ab='-test') 
    call parameter_handler%add(print_nonlinear_iteration_key, '--print-nonlinear-iteration', .false., 'Print nonlinear iteration output', switch_ab='-nl-it') 
  end subroutine test_transient_poisson_define_user_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(test_transient_poisson_define_user_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                      :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    integer(ip)                                         :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                             :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                       :: get_triangulation_type
    call parameter_handler%GetAsString(key = static_triang_generate_from_key, string = get_triangulation_type)
  end function get_triangulation_type 
  
  !==================================================================================================
  function get_initial_time(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                            :: get_initial_time
    call parameter_handler%Get(key = initial_time_key, Value = get_initial_time)
  end function get_initial_time
  
    !==================================================================================================
  function get_final_time(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                            :: get_final_time
    call parameter_handler%Get(key = final_time_key, Value = get_final_time)
  end function get_final_time
  
    !==================================================================================================
  function get_time_step(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                            :: get_time_step
    call parameter_handler%Get(key = time_step_key, Value = get_time_step)
  end function get_time_step
  
    !==================================================================================================
  function get_time_integration_scheme(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                       :: get_time_integration_scheme
    call parameter_handler%GetAsString(key = time_integration_scheme_key, string = get_time_integration_scheme)
  end function get_time_integration_scheme
    
    !==================================================================================================
  function get_is_test(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                             :: get_is_test
    call parameter_handler%Get(key = is_test_key, Value = get_is_test)
  end function get_is_test
    
    !==================================================================================================
  function get_print_nonlinear_iteration(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                   :: get_print_nonlinear_iteration
    get_print_nonlinear_iteration = this%print_nonlinear_iteration
  end function get_print_nonlinear_iteration
  
end module test_transient_poisson_params_names
