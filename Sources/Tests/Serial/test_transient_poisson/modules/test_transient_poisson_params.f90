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
  character(len=*), parameter :: fe_formulation_key            = 'fe_formulation'
  character(len=*), parameter :: laplacian_type_key            = 'laplacian_type'
  character(len=*), parameter :: reference_fe_geo_order_key    = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key            = 'write_solution'        
  character(len=*), parameter :: use_void_fes_key              = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key         = 'use_void_fes_case'
  character(len=*), parameter :: initial_time_key              = 'initial_time'
  character(len=*), parameter :: final_time_key                = 'final_time'
  character(len=*), parameter :: time_step_key                 = 'time_step'
  character(len=*), parameter :: num_time_steps_key            = 'num_time_steps'
  character(len=*), parameter :: time_integration_scheme_key   = 'time_integration_scheme'
  character(len=*), parameter :: is_test_key                   = 'is_test'

  type, extends(fempar_parameter_handler_t) :: test_transient_poisson_params_t  
   private
     logical :: print_nonlinear_iteration = .false.
   contains
     procedure                              :: define_parameters   => test_transient_poisson_define_parameters
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_fe_formulation
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_laplacian_type
     procedure, non_overridable             :: get_triangulation_type
     procedure, non_overridable             :: get_num_dims
     procedure, non_overridable             :: get_use_void_fes
     procedure, non_overridable             :: get_use_void_fes_case
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
  subroutine test_transient_poisson_define_parameters(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(inout) :: this

    ! IO parameters
    call this%add(fe_formulation_key, '--fe-formulation', 'cG', 'cG or dG FE formulation for Poisson problem', switch_ab='-f')
    call this%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1,'Order of the triangulation reference fe', switch_ab='-gorder')
    call this%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order') 
    call this%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution') 
    call this%add(laplacian_type_key, '--laplacian-type', 'scalar', 'Scalar or Vector-Valued Laplacian PDE (scalar,vector)?', switch_ab='-lt') 

    call this%add(use_void_fes_key, '--use-void-fes', .false., 'Use a hybrid FE space formed by full and void FEs', switch_ab='-use-voids') 
    call this%add(use_void_fes_case_key, '--use-void-fes-case', 'popcorn', 'Select where to put void fes using one of the predefined patterns.', switch_ab='-use-voids-case (popcorn,half,quarter)') 

    call this%add(initial_time_key, '--initial-time', 0.0_rp, 'Initial time: t0', switch_ab='-t0') 
    call this%add(final_time_key, '--final-time', 1.0_rp, 'Final time: tf', switch_ab='-tf') 
    call this%add(time_step_key, '--time-step', 1.0_rp, 'Time step size: dt', switch_ab='-dt') 
    call this%add(num_time_steps_key, '--num-time-steps', 1, 'Maximum number of time steps: nt', switch_ab='-nt') 
    call this%add(time_integration_scheme_key, '--time-integration-scheme', 'backward_euler', 'Time disctetization scheme of the DIRK solver.', switch_ab='-rk-scheme') 

    call this%add(is_test_key, '--is-test', .false., 'Test convergence order of the runge kutta scheme', switch_ab='-test') 
    call this%add(print_nonlinear_iteration_key, '--print-nonlinear-iteration', .false., 'Print nonlinear iteration output', switch_ab='-nl-it') 
  end subroutine test_transient_poisson_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  function get_dir_path_out(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path_out
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_out_key, 'string'))
    error = list%GetAsString(key = dir_path_out_key, string = get_dir_path_out)
    assert(error==0)
  end function get_dir_path_out

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix

    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_geo_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_geo_order_key, get_reference_fe_geo_order))
    error = list%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
    assert(error==0)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    logical                                       :: is_present
    logical                                       :: same_data_type
    integer(ip), allocatable                      :: shape(:)
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_use_void_fes(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                       :: get_use_void_fes
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_void_fes_key, get_use_void_fes))
    error = list%Get(key = use_void_fes_key, Value = get_use_void_fes)
    assert(error==0)
  end function get_use_void_fes

  !==================================================================================================
  function get_use_void_fes_case(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_use_void_fes_case
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_void_fes_case_key, 'string'))
    error = list%GetAsString(key = use_void_fes_case_key, string = get_use_void_fes_case)
    assert(error==0)
  end function get_use_void_fes_case
  
  !==================================================================================================
  function get_initial_time(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                      :: get_initial_time
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(initial_time_key, get_initial_time))
    error = list%Get(key = initial_time_key, Value = get_initial_time)
    assert(error==0)
  end function get_initial_time
  
    !==================================================================================================
  function get_final_time(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                      :: get_final_time
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(final_time_key, get_final_time))
    error = list%Get(key = final_time_key, Value = get_final_time)
    assert(error==0)
  end function get_final_time
  
    !==================================================================================================
  function get_time_step(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                      :: get_time_step
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(time_step_key, get_time_step))
    error = list%Get(key = time_step_key, Value = get_time_step)
    assert(error==0)
  end function get_time_step
  
    !==================================================================================================
  function get_time_integration_scheme(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_time_integration_scheme
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(time_integration_scheme_key, 'string'))
    error = list%GetAsString(key = time_integration_scheme_key, string = get_time_integration_scheme)
    assert(error==0)
  end function get_time_integration_scheme
    
    !==================================================================================================
  function get_is_test(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                       :: get_is_test
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(is_test_key, get_is_test))
    error = list%Get(key = is_test_key, Value = get_is_test)
    assert(error==0)
  end function get_is_test

  !==================================================================================================
  function get_fe_formulation(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                       :: get_fe_formulation
    type(ParameterList_t), pointer                      :: list
    integer(ip)                                         :: error
    list  => this%get_values()
    assert(list%isAssignable(fe_formulation_key, 'string'))
    error = list%GetAsString(key = fe_formulation_key, string = get_fe_formulation)
    assert(error==0)
  end function get_fe_formulation
  
  
  !==================================================================================================
  function get_laplacian_type(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_laplacian_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(laplacian_type_key, 'string'))
    error = list%GetAsString(key = laplacian_type_key, string = get_laplacian_type)
    assert(error==0)
  end function get_laplacian_type 

  !==================================================================================================
  function get_num_dims(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    integer(ip) :: get_num_dims
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_num_dims_key, get_num_dims))
    error = list%Get(key = struct_hex_triang_num_dims_key, value = get_num_dims)
    assert(error==0)
  end function get_num_dims
  
    !==================================================================================================
  function get_print_nonlinear_iteration(this)
    implicit none
    class(test_transient_poisson_params_t) , intent(in) :: this
    logical                                   :: get_print_nonlinear_iteration
    get_print_nonlinear_iteration = this%print_nonlinear_iteration
  end function get_print_nonlinear_iteration
  
end module test_transient_poisson_params_names
