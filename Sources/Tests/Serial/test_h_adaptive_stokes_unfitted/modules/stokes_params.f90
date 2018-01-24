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

  type stokes_params_t  
     private 
     ! IO parameters
     character(len=:), allocatable :: default_dir_path 
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     character(len=:), allocatable :: default_fe_formulation
     character(len=:), allocatable :: default_reference_fe_geo_order
     character(len=:), allocatable :: default_reference_fe_order
     character(len=:), allocatable :: default_write_solution
     character(len=:), allocatable :: default_write_matrix
     character(len=:), allocatable :: default_write_error_norms
     character(len=:), allocatable :: default_write_aggr_info
     character(len=:), allocatable :: default_laplacian_type
     character(len=:), allocatable :: default_triangulation_type
     character(len=:), allocatable :: default_num_dims
     character(len=:), allocatable :: default_nx
     character(len=:), allocatable :: default_ny
     character(len=:), allocatable :: default_nz
     character(len=:), allocatable :: default_is_periodic_in_x
     character(len=:), allocatable :: default_is_periodic_in_y
     character(len=:), allocatable :: default_is_periodic_in_z
     character(len=:), allocatable :: default_max_level
     character(len=:), allocatable :: default_case_id
     character(len=:), allocatable :: default_bc_case_id
     character(len=:), allocatable :: default_check_solution
     character(len=:), allocatable :: default_unfitted_boundary_is_dirichlet
     character(len=:), allocatable :: default_is_constant_nitches_beta
     character(len=:), allocatable :: default_use_constraints
     character(len=:), allocatable :: default_levelset_function_type
     character(len=:), allocatable :: default_levelset_tolerance
     character(len=:), allocatable :: default_domain_limits
     character(len=:), allocatable :: default_only_setup
     character(len=:), allocatable :: default_strong_dirichlet_on_fitted_boundary
     character(len=:), allocatable :: default_refinement_pattern
     character(len=:), allocatable :: default_lin_solver_type
     character(len=:), allocatable :: default_use_levelset_complement

     type(Command_Line_Interface):: cli 

     ! IO parameters
     character(len=str_cla_len)    :: dir_path
     character(len=str_cla_len)    :: prefix
     character(len=str_cla_len)    :: dir_path_out
     character(len=str_cla_len)    :: fe_formulation
     integer(ip)                   :: reference_fe_geo_order
     integer(ip)                   :: reference_fe_order
     logical                       :: write_solution
     logical                       :: write_matrix
     logical                       :: write_error_norms
     logical                       :: write_aggr_info
     character(len=str_cla_len)    :: laplacian_type

     character(len=str_cla_len)    :: triangulation_type
     integer(ip) :: num_dims     
     integer(ip) :: num_of_cells_x_dir(0:SPACE_DIM-1)
     integer(ip) :: is_dir_periodic(0:SPACE_DIM-1)
     integer(ip) :: max_level
     integer(ip) :: case_id
     integer(ip) :: bc_case_id
     logical :: check_sol
     logical :: unfitted_boundary_is_dirichlet
     logical :: is_constant_nitches_beta
     logical :: use_constraints
     character(len=str_cla_len)    :: levelset_function_type
     real(rp) :: levelset_tolerance
     real(rp)    :: domain_limits(2)
     logical :: only_setup
     logical :: strong_dirichlet_on_fitted_boundary
     character(len=str_cla_len)    :: refinement_pattern
     character(len=str_cla_len)    :: lin_solver_type
     logical :: use_levelset_complement

   contains
     procedure, non_overridable             :: create       => stokes_create
     procedure, non_overridable, private    :: set_default  => stokes_set_default
     procedure, non_overridable, private    :: add_to_cli   => stokes_add_to_cli
     procedure, non_overridable             :: parse        => stokes_parse 
     procedure, non_overridable             :: free         => stokes_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
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

  subroutine stokes_create(this)
    implicit none
    class(stokes_params_t), intent(inout) :: this
    
    call this%free()
    
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'test_poisson',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR test to solve the 2D Poisson PDE with known analytical solution. &
                                  Boundary set ID 1 MUST BE ASSIGNED to the whole boundary.', &
         &        examples    = ['test_poisson -h  ', 'test_poisson -h  ' ])
    
    call this%set_default()
    call this%add_to_cli()
  end subroutine stokes_create
  
  subroutine stokes_set_default(this)
    implicit none
    class(stokes_params_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = 'data/'
    this%default_prefix         = 'square'
    this%default_dir_path_out   = 'output/'
    this%default_reference_fe_geo_order = '1'
    this%default_reference_fe_order = '1'
    this%default_write_solution = '.false.'
    this%default_write_matrix   = '.false.'
    this%default_write_error_norms = '.false.'
    this%default_write_aggr_info   = '.false.'
    this%default_laplacian_type = 'scalar'
    
    this%default_triangulation_type = 'unstructured'
    this%default_num_dims = '2'
    this%default_nx = '1'
    this%default_ny = '1'
    this%default_nz = '1'
    this%default_is_periodic_in_x = '0'
    this%default_is_periodic_in_y = '0'
    this%default_is_periodic_in_z = '0'
    this%default_max_level = '3'
    this%default_case_id = '1'
    this%default_bc_case_id = '1'
    this%default_check_solution = '.true.'
    this%default_unfitted_boundary_is_dirichlet = '.true.'
    this%default_is_constant_nitches_beta       = '.false.'
    this%default_use_constraints = '.true.'
    this%default_levelset_function_type = 'sphere'
    this%default_levelset_tolerance = '1.0e-6'
    this%default_domain_limits = '0.0 1.0'
    this%default_only_setup = '.false.'
    this%default_strong_dirichlet_on_fitted_boundary = '.true.'
    this%default_refinement_pattern = 'uniform'
    this%default_lin_solver_type = 'pardiso'
    this%default_use_levelset_complement = '.false.'
    
  end subroutine stokes_set_default
  
  !==================================================================================================
  subroutine stokes_add_to_cli(this)
    implicit none
    class(stokes_params_t) , intent(inout) :: this

    ! Locals
    integer(ip) :: error

    ! IO parameters
    call this%cli%add(switch='--dir-path',switch_ab='-d',                              &
         &            help='Directory of the source files',required=.false., act='store',                &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Name of the GiD files',  &
         &            required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--dir-path-out',switch_ab='-o',help='Output Directory',&
         &            required=.false.,act='store',def=trim(this%default_dir_path_out),error=error)
    check(error==0)  
    call this%cli%add(switch='--reference-fe-geo-order',switch_ab='-gorder',help='Order of the triangulation reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_geo_order),error=error)
    check(error==0)  
    call this%cli%add(switch='--reference-fe-order',switch_ab='-order',help='Order of the fe space reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_order),error=error) 
    check(error==0) 
    call this%cli%add(switch='--write-solution',switch_ab='-wsolution',help='Write solution in VTK format',&
         &            required=.false.,act='store',def=trim(this%default_write_solution),error=error) 
    check(error==0)
    call this%cli%add(switch='--write-matrix',switch_ab='-wmatrix',help='Write matrix in matrix market format',&
         &            required=.false.,act='store',def=trim(this%default_write_matrix),error=error) 
    check(error==0)
    call this%cli%add(switch='--write-error-norms',switch_ab='-werrornorms',help='Write error norms in csv format',&
         &            required=.false.,act='store',def=trim(this%default_write_error_norms),error=error) 
    check(error==0)
    call this%cli%add(switch='--write-aggr-info',switch_ab='-waggrinfo',help='Write info about the aggregates in csv format',&
         &            required=.false.,act='store',def=trim(this%default_write_aggr_info),error=error) 
    check(error==0)
    call this%cli%add(switch='--laplacian-type',switch_ab='-lt',help='Scalar or Vector-Valued Laplacian PDE?',&
         &            required=.false.,act='store',def=trim(this%default_laplacian_type),choices='scalar,vector',error=error) 
    check(error==0)
    call this%cli%add(switch='--triangulation-type',switch_ab='-tt',help='Structured or unstructured (GiD) triangulation?',&
         &            required=.false.,act='store',def=trim(this%default_triangulation_type),choices='structured,unstructured',error=error) 
    check(error==0) 
    call this%cli%add(switch='--num_of_dims',switch_ab='-dim',help='Number of space dimensions',&
         &            required=.false.,act='store',def=trim(this%default_num_dims),error=error) 
    check(error==0) 
    call this%cli%add(switch='--num_of_cells_in_x',switch_ab='-nx',help='Number of cells in x',&
         &            required=.false.,act='store',def=trim(this%default_nx),error=error) 
    check(error==0) 
    call this%cli%add(switch='--num_of_cells_in_y',switch_ab='-ny',help='Number of cells in y',&
         &            required=.false.,act='store',def=trim(this%default_ny),error=error) 
    check(error==0) 
    call this%cli%add(switch='--num_of_cells_in_z',switch_ab='-nz',help='Number of cells in z',&
         &            required=.false.,act='store',def=trim(this%default_nz),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_x',switch_ab='-px',help='Is the mesh periodic in x',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_x),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_y',switch_ab='-py',help='Is the mesh periodic in y',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_y),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_z',switch_ab='-pz',help='Is the mesh periodic in z',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_z),error=error) 
    check(error==0) 
    call this%cli%add(switch='--max_level',switch_ab='-maxl',help='Maximum h-refinement level allowed',&
         &            required=.false.,act='store',def=trim(this%default_max_level),error=error) 
    check(error==0) 
    call this%cli%add(switch='--case_id',switch_ab='-cid',help='Id of the functions used in the run',&
         &            required=.false.,act='store',def=trim(this%default_case_id),error=error) 
    check(error==0) 
    call this%cli%add(switch='--bc_case_id',switch_ab='-bcid',help='Id of the boundary setup used in the run',&
         &            required=.false.,act='store',def=trim(this%default_bc_case_id),error=error) 
    check(error==0) 
    call this%cli%add(switch='--check_solution',switch_ab='-check',help='Check or not the solution',&
         &            required=.false.,act='store',def=trim(this%default_check_solution),error=error) 
    check(error==0) 
    call this%cli%add(switch='--is_dirichlet',switch_ab='-is_diri',help='True if the unfitted boundary is dirichlet',&
         &            required=.false.,act='store',def=trim(this%default_unfitted_boundary_is_dirichlet),error=error) 
    check(error==0) 
    call this%cli%add(switch='--is_beta_constant',switch_ab='-is_bconst',help='True if the Nitsches beta is constant',&
         &            required=.false.,act='store',def=trim(this%default_is_constant_nitches_beta),error=error) 
    check(error==0) 
    call this%cli%add(switch='--use_constraints',switch_ab='-uconstraints',help='Use or not the constraints provided by the cut cell aggregation',&
         &            required=.false.,act='store',def=trim(this%default_use_constraints),error=error) 
    check(error==0) 
    call this%cli%add(switch='--levelset-type',switch_ab='-lstype',help='Name of the levelset function',&
         &            required=.false.,act='store',def=trim(this%default_levelset_function_type),error=error) 
    check(error==0) 
    call this%cli%add(switch='--levelset-tol',switch_ab='-lstol',help='Tolerance for the levelset function',&
         &            required=.false.,act='store',def=trim(this%default_levelset_tolerance),error=error) 
    check(error==0) 
    call this%cli%add(switch='--domain-limits',switch_ab='-dom',help='Info about the domain limits',&
         &            required=.false.,act='store',def=trim(this%default_domain_limits),error=error,nargs='2') 
    check(error==0) 
    call this%cli%add(switch='--only-setup',switch_ab='-osetup',help='True if compute only the setup of the problem, i.e., skip discrete integration and linear solver',&
         &            required=.false.,act='store',def=trim(this%default_only_setup),error=error) 
    check(error==0) 
    call this%cli%add(switch='--strong_dirichlet',switch_ab='-sdiri',help='True if strong dirichlet conditions are imposed on the body-fitted boundary',&
         &            required=.false.,act='store',def=trim(this%default_strong_dirichlet_on_fitted_boundary),error=error) 
    call this%cli%add(switch='--refinement_pattern',switch_ab='-rpattern',help='name of the refinement pattern to use',&
         &            required=.false.,act='store',def=trim(this%default_refinement_pattern),error=error) 
    call this%cli%add(switch='--lin_solver_type',switch_ab='-lsolver',help='name of the linear solver to use',&
         &            required=.false.,act='store',def=trim(this%default_lin_solver_type),error=error) 
    check(error==0) 
    call this%cli%add(switch='--use_levelset_complement',switch_ab='-ulscomp',help='if true then we use the complement of the levelset',&
         &            required=.false.,act='store',def=trim(this%default_use_levelset_complement),error=error) 
    check(error==0) 
  end subroutine stokes_add_to_cli
  
  subroutine stokes_parse(this,parameter_list)
    implicit none
    class(stokes_params_t), intent(inout) :: this
    type(ParameterList_t)       , intent(inout) :: parameter_list
    integer(ip) :: istat
    
    call this%cli%parse(error=istat); check(istat==0)
    
    ! IO parameters
    call this%cli%get(switch='-d',val=this%dir_path    ,error=istat); check(istat==0)
    call this%cli%get(switch='-p',val=this%prefix       ,error=istat); check(istat==0)
    call this%cli%get(switch='-o',val=this%dir_path_out,error=istat); check(istat==0)
    call this%cli%get(switch='-gorder',val=this%reference_fe_geo_order,error=istat); check(istat==0)
    call this%cli%get(switch='-order',val=this%reference_fe_order,error=istat); check(istat==0)
    call this%cli%get(switch='-wsolution',val=this%write_solution,error=istat); check(istat==0)
    call this%cli%get(switch='-wmatrix',val=this%write_matrix,error=istat); check(istat==0)
    call this%cli%get(switch='-werrornorms',val=this%write_error_norms,error=istat); check(istat==0)
    call this%cli%get(switch='-waggrinfo',val=this%write_aggr_info,error=istat); check(istat==0)
    call this%cli%get(switch='-lt',val=this%laplacian_type,error=istat); check(istat==0)
    call this%cli%get(switch='-tt',val=this%triangulation_type,error=istat); check(istat==0)
    call this%cli%get(switch='-dim',val=this%num_dims,error=istat); check(istat==0)
    call this%cli%get(switch='-nx',val=this%num_of_cells_x_dir(0),error=istat); check(istat==0)
    call this%cli%get(switch='-ny',val=this%num_of_cells_x_dir(1),error=istat); check(istat==0)
    call this%cli%get(switch='-nz',val=this%num_of_cells_x_dir(2),error=istat); check(istat==0)
    call this%cli%get(switch='-px',val=this%is_dir_periodic(0),error=istat); check(istat==0)
    call this%cli%get(switch='-py',val=this%is_dir_periodic(1),error=istat); check(istat==0)
    call this%cli%get(switch='-pz',val=this%is_dir_periodic(2),error=istat); check(istat==0)
    call this%cli%get(switch='-maxl',val=this%max_level,error=istat); check(istat==0)
    call this%cli%get(switch='-cid',val=this%case_id,error=istat); check(istat==0)
    call this%cli%get(switch='-bcid',val=this%bc_case_id,error=istat); check(istat==0)
    call this%cli%get(switch='-check',val=this%check_sol,error=istat); check(istat==0)
    call this%cli%get(switch='-is_diri',val=this%unfitted_boundary_is_dirichlet,error=istat); check(istat==0)
    call this%cli%get(switch='-is_bconst',val=this%is_constant_nitches_beta,error=istat); check(istat==0)
    call this%cli%get(switch='-uconstraints',val=this%use_constraints,error=istat); check(istat==0)
    call this%cli%get(switch='-lstype',val=this%levelset_function_type,error=istat); check(istat==0)
    call this%cli%get(switch='-lstol',val=this%levelset_tolerance,error=istat); check(istat==0)
    call this%cli%get(switch='-dom',val=this%domain_limits,error=istat); check(istat==0)
    call this%cli%get(switch='-osetup',val=this%only_setup,error=istat); check(istat==0)
    call this%cli%get(switch='-sdiri',val=this%strong_dirichlet_on_fitted_boundary,error=istat); check(istat==0)
    call this%cli%get(switch='-rpattern',val=this%refinement_pattern,error=istat); check(istat==0)
    call this%cli%get(switch='-lsolver',val=this%lin_solver_type,error=istat); check(istat==0)
    call this%cli%get(switch='-ulscomp',val=this%use_levelset_complement,error=istat); check(istat==0)

    call parameter_list%init()
    istat = 0
    istat = istat + parameter_list%set(key = dir_path_key, value = this%dir_path)
    istat = istat + parameter_list%set(key = prefix_key  , value = this%prefix)
    istat = istat + parameter_list%set(key = geometry_interpolation_order_key  , value = this%reference_fe_geo_order)
    check(istat==0)

    if(trim(this%triangulation_type)=='unstructured') then
       istat = parameter_list%set(key = triangulation_generate_key, value = triangulation_generate_from_mesh)
    else if(trim(this%triangulation_type)=='structured') then
       istat = parameter_list%set(key = triangulation_generate_key         , value = triangulation_generate_structured)
       istat = istat + parameter_list%set(key = num_dims_key   , value = this%num_dims)
       istat = istat + parameter_list%set(key = num_cells_x_dir_key, value = this%num_of_cells_x_dir)
       istat = istat + parameter_list%set(key = is_dir_periodic_key        , value = this%is_dir_periodic)
    end if
    check(istat==0)
  end subroutine stokes_parse  

  subroutine stokes_free(this)
    implicit none
    class(stokes_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_dir_path_out)) deallocate(this%default_dir_path_out)
    if(allocated(this%default_reference_fe_geo_order)) deallocate(this%default_reference_fe_geo_order)
    if(allocated(this%default_reference_fe_order)) deallocate(this%default_reference_fe_order)
    if(allocated(this%default_write_solution)) deallocate(this%default_write_solution)
    if(allocated(this%default_write_matrix)) deallocate(this%default_write_matrix)
    if(allocated(this%default_write_error_norms)) deallocate(this%default_write_error_norms)
    if(allocated(this%default_write_aggr_info)) deallocate(this%default_write_aggr_info)
    if(allocated(this%default_laplacian_type)) deallocate(this%default_laplacian_type)
    if(allocated(this%default_max_level)) deallocate(this%default_max_level)
    if(allocated(this%default_case_id)) deallocate(this%default_case_id)
    if(allocated(this%default_bc_case_id)) deallocate(this%default_bc_case_id)
    if(allocated(this%default_check_solution)) deallocate(this%default_check_solution)
    if(allocated(this%default_unfitted_boundary_is_dirichlet)) deallocate(this%default_unfitted_boundary_is_dirichlet)
    if(allocated(this%default_is_constant_nitches_beta)) deallocate(this%default_is_constant_nitches_beta)
    if(allocated(this%default_use_constraints)) deallocate(this%default_use_constraints)
    if(allocated(this%default_levelset_function_type)) deallocate(this%default_levelset_function_type)
    if(allocated(this%default_levelset_tolerance)) deallocate(this%default_levelset_tolerance)
    if(allocated(this%default_domain_limits)) deallocate(this%default_domain_limits)
    if(allocated(this%default_only_setup)) deallocate(this%default_only_setup)
    if(allocated(this%default_strong_dirichlet_on_fitted_boundary)) deallocate(this%default_strong_dirichlet_on_fitted_boundary)
    if(allocated(this%default_refinement_pattern)) deallocate(this%default_refinement_pattern)
    if(allocated(this%default_lin_solver_type)) deallocate(this%default_lin_solver_type)
    if(allocated(this%default_use_levelset_complement)) deallocate(this%default_use_levelset_complement)
    call this%cli%free()
  end subroutine stokes_free

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_dir_path
    get_dir_path = trim(this%dir_path)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_prefix
    get_prefix = trim(this%prefix)
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_dir_path_out
    get_dir_path_out = trim(this%dir_path_out)
  end function get_dir_path_out
  
  !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_geo_order
    get_reference_fe_geo_order = this%reference_fe_geo_order
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_order
    get_reference_fe_order = this%reference_fe_order
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_write_solution
    get_write_solution = this%write_solution
  end function get_write_solution

  !==================================================================================================
  function get_write_matrix(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_write_matrix
    get_write_matrix = this%write_matrix
  end function get_write_matrix

  !==================================================================================================
  function get_write_error_norms(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_write_error_norms
    get_write_error_norms = this%write_error_norms
  end function get_write_error_norms

  !==================================================================================================
  function get_write_aggr_info(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_write_aggr_info
    get_write_aggr_info = this%write_aggr_info
  end function get_write_aggr_info
  
  !==================================================================================================
  function get_laplacian_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_laplacian_type
    get_laplacian_type = trim(this%laplacian_type)
  end function get_laplacian_type 
  
  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_triangulation_type
    get_triangulation_type = trim(this%triangulation_type)
  end function get_triangulation_type 
  
  !==================================================================================================
  function get_num_dims(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_num_dims
    get_num_dims = this%num_dims
  end function get_num_dims

  !==================================================================================================
  function get_max_level(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_max_level
    get_max_level = this%max_level
  end function get_max_level

  !==================================================================================================
  function get_case_id(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_case_id
    get_case_id = this%case_id
  end function get_case_id

  !==================================================================================================
  function get_bc_case_id(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    integer(ip) :: get_bc_case_id
    get_bc_case_id = this%bc_case_id
  end function get_bc_case_id

  !==================================================================================================
  function are_checks_active(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: are_checks_active
    are_checks_active = this%check_sol
  end function are_checks_active

  !==================================================================================================
  function get_unfitted_boundary_is_dirichlet(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_unfitted_boundary_is_dirichlet
    get_unfitted_boundary_is_dirichlet = this%unfitted_boundary_is_dirichlet
  end function get_unfitted_boundary_is_dirichlet

  !==================================================================================================
  function get_is_constant_nitches_beta(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_is_constant_nitches_beta
    get_is_constant_nitches_beta = this%is_constant_nitches_beta
  end function get_is_constant_nitches_beta

  !==================================================================================================
  function get_use_constraints(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_use_constraints
    get_use_constraints = this%use_constraints
  end function get_use_constraints

  !==================================================================================================
  function get_levelset_function_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_levelset_function_type
    get_levelset_function_type = trim(this%levelset_function_type)
  end function get_levelset_function_type 

  !==================================================================================================
  function get_levelset_tolerance(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    real(rp) :: get_levelset_tolerance
    get_levelset_tolerance = this%levelset_tolerance
  end function get_levelset_tolerance

  !==================================================================================================
  function get_domain_limits(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    real(rp) :: get_domain_limits(2)
    get_domain_limits = this%domain_limits
  end function get_domain_limits 

  !==================================================================================================
  function get_only_setup(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_only_setup
    get_only_setup = this%only_setup
  end function get_only_setup

  !==================================================================================================
  function is_strong_dirichlet_on_fitted_boundary(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: is_strong_dirichlet_on_fitted_boundary
    is_strong_dirichlet_on_fitted_boundary = this%strong_dirichlet_on_fitted_boundary
  end function is_strong_dirichlet_on_fitted_boundary

  !==================================================================================================
  function get_refinement_pattern(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_refinement_pattern
    get_refinement_pattern = trim(this%refinement_pattern)
  end function get_refinement_pattern 

  !==================================================================================================
  function get_lin_solver_type(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_lin_solver_type
    get_lin_solver_type = trim(this%lin_solver_type)
  end function get_lin_solver_type 

  !==================================================================================================
  function get_use_levelset_complement(this)
    implicit none
    class(stokes_params_t) , intent(in) :: this
    logical :: get_use_levelset_complement
    get_use_levelset_complement = this%use_levelset_complement
  end function get_use_levelset_complement

end module stokes_params_names
