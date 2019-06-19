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
module test_poisson_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private

  character(len=*), parameter :: fe_formulation_key            = 'fe_formulation'
  character(len=*), parameter :: laplacian_type_key            = 'laplacian_type'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key            = 'write_solution'
  character(len=*), parameter :: write_matrix_key              = 'write_matrix'
  character(len=*), parameter :: write_error_norms_key         = 'write_error_norms'
  character(len=*), parameter :: write_aggr_info_key           = 'write_aggr_info'
  character(len=*), parameter :: max_level_key                 = 'max_level'  
  character(len=*), parameter :: check_solution_key            = 'check_solution'
  character(len=*), parameter :: is_dirichlet_key              = 'is_dirichlet'
  character(len=*), parameter :: is_beta_constant_key          = 'is_beta_constant'
  character(len=*), parameter :: use_constraints_key           = 'use_constraints'
  character(len=*), parameter :: levelset_type_key             = 'levelset_type'
  character(len=*), parameter :: levelset_tol_key              = 'levelset_tol'
  character(len=*), parameter :: domain_limits_key             = 'domain_limits'
  character(len=*), parameter :: only_setup_key                = 'only_setup'
  character(len=*), parameter :: strong_dirichlet_key          = 'strong_dirichlet'
  character(len=*), parameter :: refinement_pattern_key        = 'refinement_pattern'
  character(len=*), parameter :: is_in_fe_space_key            = 'is_in_fe_space'
  character(len=*), parameter :: dir_path_out_key              = 'dir_path_out'
  character(len=*), parameter :: out_prefix_key                = 'out_prefix'

  type :: test_poisson_params_t  
   contains
     procedure, non_overridable             :: process_parameters
     procedure, non_overridable             :: get_parameter_list
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_output_handler_prefix
     procedure, non_overridable             :: get_output_handler_dir_path
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_write_matrix
     procedure, non_overridable             :: get_write_error_norms
     procedure, non_overridable             :: get_write_aggr_info
     procedure, non_overridable             :: get_laplacian_type
     procedure, non_overridable             :: get_max_level
     procedure, non_overridable             :: is_in_fe_space
     procedure, non_overridable             :: are_checks_active
     procedure, non_overridable             :: get_unfitted_boundary_is_dirichlet
     procedure, non_overridable             :: get_is_constant_nitches_beta
     procedure, non_overridable             :: get_use_constraints
     procedure, non_overridable             :: get_levelset_function_type
     procedure, non_overridable             :: get_levelset_tolerance
     procedure, non_overridable             :: get_domain_limits
     procedure, non_overridable             :: get_only_setup
     procedure, non_overridable             :: is_strong_dirichlet_on_fitted_boundary
     procedure, non_overridable             :: get_refinement_pattern
     procedure, non_overridable             :: get_struct_hex_mesh_generator_num_dims_key
  end type test_poisson_params_t  

  ! Types
  public :: test_poisson_params_t

contains
  
  !==================================================================================================
  subroutine test_poisson_define_user_parameters()
    implicit none
    ! IO parameters
    call parameter_handler%add(dir_path_out_key, '--dir-path-out', '.', 'Output path',  switch_ab='-dpo') 
    call parameter_handler%add(out_prefix_key, '--prefix', 'output', 'Filename prefix for output results',  switch_ab='-p') 
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe',  switch_ab='-order') 
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution') 
    call parameter_handler%add(write_matrix_key, '--write-matrix', .false., 'Write matrix in matrix market format', switch_ab='-wmatrix') 
    call parameter_handler%add(write_error_norms_key, '--write-error-norms', .false., 'Write error norms in csv format', switch_ab='-werrornorms') 
    call parameter_handler%add(write_aggr_info_key, '--write-aggr-info', .false., 'Write info about the aggregates in csv format', switch_ab='-waggrinfo') 
    call parameter_handler%add(laplacian_type_key, '--laplacian-type', 'scalar', 'Scalar or Vector-Valued Laplacian PDE? (scalar,vector)', switch_ab='-lt') 
    call parameter_handler%add(max_level_key, '--max_level', 3, 'Maximum h-refinement level allowed', switch_ab='-maxl') 
    call parameter_handler%add(check_solution_key, '--check_solution', .true., 'Check or not the solution', switch_ab='-check') 
    call parameter_handler%add(is_dirichlet_key, '--is_dirichlet', .true., 'True if the unfitted boundary is dirichlet', switch_ab='-is_diri') 
    call parameter_handler%add(is_beta_constant_key, '--is_beta_constant', .false., 'True if the Nitsches beta is constant', switch_ab='-is_bconst') 
    call parameter_handler%add(use_constraints_key, '--use_constraints', .true., 'Use or not the constraints provided by the cut cell aggregation', switch_ab='-uconstraints') 
    call parameter_handler%add(levelset_type_key, '--levelset-type', 'sphere', 'Name of the levelset function', switch_ab='-lstype') 
    call parameter_handler%add(levelset_tol_key, '--levelset-tol', 1.0e-6_rp,'Tolerance for the levelset function', switch_ab='-lstol') 
    call parameter_handler%add(domain_limits_key, '--domain-limits', [ 0.0_rp, 1.0_rp], 'Info about the domain limits', switch_ab='-dom') 
    call parameter_handler%add(only_setup_key, '--only-setup', .false., &
                'True if compute only the setup of the problem, i.e., skip discrete integration and linear solver', &
                switch_ab='-osetup') 
    call parameter_handler%add(strong_dirichlet_key, '--strong_dirichlet', .true., &
                'True if strong dirichlet conditions are imposed on the body-fitted boundary', &
                switch_ab='-sdiri') 
    call parameter_handler%add(refinement_pattern_key, '--refinement_pattern', 'uniform', 'name of the refinement pattern to use', switch_ab='-rpattern') 
    call parameter_handler%add(is_in_fe_space_key, '--solution_in_fe_space', .true., 'Is the solution in fe space', switch_ab='-in_space') 

  end subroutine test_poisson_define_user_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(test_poisson_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(test_poisson_define_user_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    type(ParameterList_t), pointer            :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list
  
  ! GETTERS *****************************************************************************************

  !==================================================================================================
  function get_dir_path(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_dir_path
    get_dir_path = parameter_handler%get_dir_path()
  end function get_dir_path

   !==================================================================================================
  function get_prefix(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_prefix
    call parameter_handler%getasstring(out_prefix_key, get_prefix)
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_dir_path_out
    call parameter_handler%getasstring(dir_path_out_key, get_dir_path_out)
  end function get_dir_path_out
  
  !==================================================================================================
  function get_output_handler_prefix(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_output_handler_prefix
    call parameter_handler%getasstring(output_handler_prefix_key, get_output_handler_prefix)
  end function get_output_handler_prefix

  !==================================================================================================
  function get_output_handler_dir_path(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_output_handler_dir_path
    call parameter_handler%getasstring(output_handler_dir_path_key, get_output_handler_dir_path)
  end function get_output_handler_dir_path

  !==================================================================================================
  function is_in_fe_space(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: is_in_fe_space
    call parameter_handler%Get(key = is_in_fe_space_key, Value = is_in_fe_space)
  end function is_in_fe_space
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    integer(ip)                               :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution

  !==================================================================================================
  function get_write_matrix(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_write_matrix
    call parameter_handler%Get(key = write_matrix_key, Value = get_write_matrix)
  end function get_write_matrix

  !==================================================================================================
  function get_write_error_norms(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_write_error_norms
    call parameter_handler%Get(key = write_error_norms_key, Value = get_write_error_norms)
  end function get_write_error_norms

  !==================================================================================================
  function get_write_aggr_info(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_write_aggr_info
    call parameter_handler%Get(key = write_aggr_info_key, Value = get_write_aggr_info)
  end function get_write_aggr_info
  
  !==================================================================================================
  function get_laplacian_type(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_laplacian_type
    call parameter_handler%GetAsString(key = laplacian_type_key, string = get_laplacian_type)
  end function get_laplacian_type 

  !==================================================================================================
  function get_max_level(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    integer(ip)                               :: get_max_level
    call parameter_handler%Get(key = max_level_key, Value = get_max_level)
  end function get_max_level

  !==================================================================================================
  function are_checks_active(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: are_checks_active
    call parameter_handler%Get(key = check_solution_key, Value = are_checks_active)
  end function are_checks_active

  !==================================================================================================
  function get_unfitted_boundary_is_dirichlet(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_unfitted_boundary_is_dirichlet
    call parameter_handler%Get(key = is_dirichlet_key, Value = get_unfitted_boundary_is_dirichlet)
  end function get_unfitted_boundary_is_dirichlet

  !==================================================================================================
  function get_is_constant_nitches_beta(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_is_constant_nitches_beta
    call parameter_handler%Get(key = is_beta_constant_key, Value = get_is_constant_nitches_beta)
  end function get_is_constant_nitches_beta

  !==================================================================================================
  function get_use_constraints(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_use_constraints
    call parameter_handler%Get(key = use_constraints_key, Value = get_use_constraints)
  end function get_use_constraints

  !==================================================================================================
  function get_levelset_function_type(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_levelset_function_type
    call parameter_handler%GetAsString(key = levelset_type_key, string = get_levelset_function_type)
  end function get_levelset_function_type 

  !==================================================================================================
  function get_levelset_tolerance(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    real(rp)                                  :: get_levelset_tolerance
    call parameter_handler%Get(key = levelset_tol_key, Value = get_levelset_tolerance)
  end function get_levelset_tolerance

  !==================================================================================================
  function get_domain_limits(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    real(rp)                                  :: get_domain_limits(2)
    call parameter_handler%Get(key = domain_limits_key, Value = get_domain_limits)
  end function get_domain_limits 

  !==================================================================================================
  function get_only_setup(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: get_only_setup
    call parameter_handler%Get(key = only_setup_key, Value = get_only_setup)
  end function get_only_setup

  !==================================================================================================
  function is_strong_dirichlet_on_fitted_boundary(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    logical                                   :: is_strong_dirichlet_on_fitted_boundary
    call parameter_handler%Get(key = strong_dirichlet_key, Value = is_strong_dirichlet_on_fitted_boundary)
  end function is_strong_dirichlet_on_fitted_boundary

  !==================================================================================================
  function get_refinement_pattern(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable             :: get_refinement_pattern
    call parameter_handler%GetAsString(key = refinement_pattern_key, string = get_refinement_pattern)
  end function get_refinement_pattern 

  !==================================================================================================
  function get_struct_hex_mesh_generator_num_dims_key(this)
    implicit none
    class(test_poisson_params_t) , intent(in) :: this
    integer(ip)                               :: get_struct_hex_mesh_generator_num_dims_key
    call parameter_handler%Get(key = struct_hex_mesh_generator_num_dims_key, value = get_struct_hex_mesh_generator_num_dims_key)
  end function get_struct_hex_mesh_generator_num_dims_key

end module test_poisson_params_names
