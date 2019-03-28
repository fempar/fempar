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
module stokes_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private

  character(len=*), parameter :: fe_formulation_key            = 'fe_formulation'
  character(len=*), parameter :: laplacian_type_key            = 'laplacian_type'
  character(len=*), parameter :: reference_fe_geo_order_key    = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key            = 'write_solution'
  character(len=*), parameter :: write_matrix_key              = 'write_matrix'
  character(len=*), parameter :: write_error_norms_key         = 'write_error_norms'
  character(len=*), parameter :: write_aggr_info_key           = 'write_aggr_info'
  character(len=*), parameter :: max_level_key                 = 'max_level'
  character(len=*), parameter :: case_id_key                   = 'case_id'        
  character(len=*), parameter :: bc_case_id_key                = 'bc_case_id'        
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
  character(len=*), parameter :: lin_solver_type_key           = 'lin_solver_type'
  character(len=*), parameter :: use_levelset_complement_key   = 'use_levelset_compoment'


  type, extends(fempar_parameter_handler_t) :: stokes_params_t  
   contains
     procedure                              :: define_parameters   => stokes_define_parameters
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_write_matrix
     procedure, non_overridable             :: get_write_error_norms
     procedure, non_overridable             :: get_write_aggr_info
     procedure, non_overridable             :: get_laplacian_type
     procedure, non_overridable             :: get_triangulation_type
     procedure, non_overridable             :: get_num_dims
     procedure, non_overridable             :: get_max_level
     procedure, non_overridable             :: get_case_id
     procedure, non_overridable             :: get_bc_case_id
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
     procedure, non_overridable             :: get_lin_solver_type
     procedure, non_overridable             :: get_use_levelset_complement
  end type stokes_params_t  

  ! Types
  public :: stokes_params_t

contains
  
  !==================================================================================================
  subroutine stokes_define_parameters(this)
    implicit none
    class(stokes_params_t) , intent(inout) :: this

    ! IO parameters
    call this%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1, 'Order of the triangulation reference fe', switch_ab='-gorder')
    call this%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe',  switch_ab='-order') 
    call this%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution') 
    call this%add(write_matrix_key, '--write-matrix', .false., 'Write matrix in matrix market format', switch_ab='-wmatrix') 
    call this%add(write_error_norms_key, '--write-error-norms', .false., 'Write error norms in csv format', switch_ab='-werrornorms') 
    call this%add(write_aggr_info_key, '--write-aggr-info', .false., 'Write info about the aggregates in csv format', switch_ab='-waggrinfo') 
    call this%add(laplacian_type_key, '--laplacian-type', 'scalar', 'Scalar or Vector-Valued Laplacian PDE? (scalar,vector)', switch_ab='-lt') 
    call this%add(max_level_key, '--max_level', 3, 'Maximum h-refinement level allowed', switch_ab='-maxl') 
    call this%add(case_id_key, '--case_id', 1, 'Id of the functions used in the run', switch_ab='-cid') 
    call this%add(bc_case_id_key, '--bc_case_id', 1, 'Id of the boundary setup used in the run', switch_ab='-bcid') 
    call this%add(check_solution_key, '--check_solution', .true., 'Check or not the solution', switch_ab='-check') 
    call this%add(is_dirichlet_key, '--is_dirichlet', .true., 'True if the unfitted boundary is dirichlet', switch_ab='-is_diri') 
    call this%add(is_beta_constant_key, '--is_beta_constant', .false., 'True if the Nitsches beta is constant', switch_ab='-is_bconst') 
    call this%add(use_constraints_key, '--use_constraints', .true., 'Use or not the constraints provided by the cut cell aggregation', switch_ab='-uconstraints') 
    call this%add(levelset_type_key, '--levelset-type', 'sphere', 'Name of the levelset function', switch_ab='-lstype') 
    call this%add(levelset_tol_key, '--levelset-tol', 1.0e-6_rp,'Tolerance for the levelset function', switch_ab='-lstol') 
    call this%add(domain_limits_key, '--domain-limits', [ 0.0_rp, 1.0_rp], 'Info about the domain limits', switch_ab='-dom') 
    call this%add(only_setup_key, '--only-setup', .false., &
                'True if compute only the setup of the problem, i.e., skip discrete integration and linear solver', &
                switch_ab='-osetup') 
    call this%add(strong_dirichlet_key, '--strong_dirichlet', .true., &
                'True if strong dirichlet conditions are imposed on the body-fitted boundary', &
                switch_ab='-sdiri') 
    call this%add(refinement_pattern_key, '--refinement_pattern', 'uniform', 'name of the refinement pattern to use', switch_ab='-rpattern') 
    call this%add(lin_solver_type_key, '--lin_solver_type', 'pardiso', 'name of the linear solver to use', switch_ab='-lsolver') 
    call this%add(use_levelset_complement_key, '--use_levelset_complement', .false., 'if true then we use the complement of the levelset', switch_ab='-ulscomp') 
  end subroutine stokes_define_parameters

  ! GETTERS *****************************************************************************************
  
  !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_reference_fe_geo_order
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_geo_order_key, get_reference_fe_geo_order))
    error = list%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
    assert(error==0)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_reference_fe_order
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_write_solution
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_write_matrix(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_write_matrix
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(write_matrix_key, get_write_matrix))
    error = list%Get(key = write_matrix_key, Value = get_write_matrix)
    assert(error==0)
  end function get_write_matrix

  !==================================================================================================
  function get_write_error_norms(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_write_error_norms
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(write_error_norms_key, get_write_error_norms))
    error = list%Get(key = write_error_norms_key, Value = get_write_error_norms)
    assert(error==0)
  end function get_write_error_norms

  !==================================================================================================
  function get_write_aggr_info(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_write_aggr_info
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(write_aggr_info_key, get_write_aggr_info))
    error = list%Get(key = write_aggr_info_key, Value = get_write_aggr_info)
    assert(error==0)
  end function get_write_aggr_info
  
  !==================================================================================================
  function get_laplacian_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable       :: get_laplacian_type
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(laplacian_type_key, get_laplacian_type))
    error = list%GetAsString(key = laplacian_type_key, string = get_laplacian_type)
    assert(error==0)
  end function get_laplacian_type 
  
  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_triangulation_type
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
  !==================================================================================================
  function get_num_dims(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_num_dims
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_num_dims_key, get_num_dims))
    error = list%Get(key = struct_hex_triang_num_dims_key, Value = get_num_dims)
    assert(error==0)
  end function get_num_dims

  !==================================================================================================
  function get_max_level(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_max_level
   type(ParameterList_t), pointer       :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(max_level_key, get_max_level))
    error = list%Get(key = max_level_key, Value = get_max_level)
    assert(error==0)
  end function get_max_level

  !==================================================================================================
  function get_case_id(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_case_id
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(case_id_key, get_case_id))
    error = list%Get(key = case_id_key, Value = get_case_id)
    assert(error==0)
  end function get_case_id

  !==================================================================================================
  function get_bc_case_id(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip)                         :: get_bc_case_id
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(bc_case_id_key, get_bc_case_id))
    error = list%Get(key = bc_case_id_key, Value = get_bc_case_id)
    assert(error==0)
  end function get_bc_case_id

  !==================================================================================================
  function are_checks_active(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: are_checks_active
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(check_solution_key, are_checks_active))
    error = list%Get(key = check_solution_key, Value = are_checks_active)
    assert(error==0)
  end function are_checks_active

  !==================================================================================================
  function get_unfitted_boundary_is_dirichlet(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_unfitted_boundary_is_dirichlet
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(is_dirichlet_key, get_unfitted_boundary_is_dirichlet))
    error = list%Get(key = is_dirichlet_key, Value = get_unfitted_boundary_is_dirichlet)
    assert(error==0)
  end function get_unfitted_boundary_is_dirichlet

  !==================================================================================================
  function get_is_constant_nitches_beta(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_is_constant_nitches_beta
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(is_beta_constant_key, get_is_constant_nitches_beta))
    error = list%Get(key = is_beta_constant_key, Value = get_is_constant_nitches_beta)
    assert(error==0)
  end function get_is_constant_nitches_beta

  !==================================================================================================
  function get_use_constraints(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_use_constraints
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(use_constraints_key, get_use_constraints))
    error = list%Get(key = use_constraints_key, Value = get_use_constraints)
    assert(error==0)
  end function get_use_constraints

  !==================================================================================================
  function get_levelset_function_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable       :: get_levelset_function_type
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(levelset_type_key, get_levelset_function_type))
    error = list%GetAsString(key = levelset_type_key, string = get_levelset_function_type)
    assert(error==0)
  end function get_levelset_function_type 

  !==================================================================================================
  function get_levelset_tolerance(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    real(rp)                            :: get_levelset_tolerance
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(levelset_tol_key, get_levelset_tolerance))
    error = list%Get(key = levelset_tol_key, Value = get_levelset_tolerance)
    assert(error==0)
  end function get_levelset_tolerance

  !==================================================================================================
  function get_domain_limits(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    real(rp)                            :: get_domain_limits(2)
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(domain_limits_key, get_domain_limits))
    error = list%Get(key = domain_limits_key, Value = get_domain_limits)
    assert(error==0)
  end function get_domain_limits 

  !==================================================================================================
  function get_only_setup(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: get_only_setup
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(only_setup_key, get_only_setup))
    error = list%Get(key = only_setup_key, Value = get_only_setup)
    assert(error==0)
  end function get_only_setup

  !==================================================================================================
  function is_strong_dirichlet_on_fitted_boundary(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical                             :: is_strong_dirichlet_on_fitted_boundary
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(strong_dirichlet_key, is_strong_dirichlet_on_fitted_boundary))
    error = list%Get(key = strong_dirichlet_key, Value = is_strong_dirichlet_on_fitted_boundary)
    assert(error==0)
  end function is_strong_dirichlet_on_fitted_boundary

  !==================================================================================================
  function get_refinement_pattern(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable       :: get_refinement_pattern
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(refinement_pattern_key, get_refinement_pattern))
    error = list%GetAsString(key = refinement_pattern_key, string = get_refinement_pattern)
    assert(error==0)
  end function get_refinement_pattern 

  !==================================================================================================
  function get_lin_solver_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable       :: get_lin_solver_type
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(lin_solver_type_key, get_lin_solver_type))
    error = list%GetAsString(key = lin_solver_type_key, string = get_lin_solver_type)
    assert(error==0)
  end function get_lin_solver_type 

  !==================================================================================================
  function get_use_levelset_complement(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_use_levelset_complement
    type(ParameterList_t), pointer      :: list
    integer(ip)                         :: error
    list  => this%get_values()
    assert(list%isAssignable(use_levelset_complement_key, get_use_levelset_complement))
    error = list%Get(key = use_levelset_complement_key, Value = get_use_levelset_complement)
    assert(error==0)
  end function get_use_levelset_complement

end module stokes_params_names
