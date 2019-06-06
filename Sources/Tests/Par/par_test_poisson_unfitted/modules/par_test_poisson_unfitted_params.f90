module par_test_poisson_unfitted_params_names
  use fempar_names
  use level_set_functions_gallery_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: coarse_fe_handler_type_key = 'coarse_fe_handler_type'    
  character(len=*), parameter :: level_set_function_type_key= 'level_set_function_type'    
  character(len=*), parameter :: use_preconditioner_key     = 'use_preconditioner'    
  character(len=*), parameter :: unfitted_boundary_type_key = 'unfitted_boundary_type'    
  character(len=*), parameter :: nitsche_beta_factor_key    = 'nitsche_beta_factor'    
  character(len=*), parameter :: levelset_tolerance_key     = 'levelset_tolerance'    
  character(len=*), parameter :: num_runs_key            = 'num_runs'    
  character(len=*), parameter :: is_in_fe_space_key         = 'is_in_fe_space'    
  character(len=*), parameter :: are_checks_active_key      = 'are_checks_active'    

  character(len=*), public, parameter :: unfitted_coarse_fe_handler_value = 'unfitted'    
  character(len=*), public, parameter :: standard_coarse_fe_handler_value = 'standard'    
  character(len=*), public, parameter :: stiffness_coarse_fe_handler_value = 'stiffness'

  type :: par_test_poisson_unfitted_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_coarse_fe_handler_type
       procedure, non_overridable             :: get_level_set_function_type
       procedure, non_overridable             :: get_use_preconditioner
       procedure, non_overridable             :: get_unfitted_boundary_type
       procedure, non_overridable             :: get_nitsche_beta_factor
       procedure, non_overridable             :: get_levelset_tolerance
       procedure, non_overridable             :: get_num_runs
       procedure, non_overridable             :: print
       procedure, non_overridable, private    :: print_character_switch
       procedure, non_overridable, private    :: print_integer_switch
       procedure, non_overridable, private    :: print_real_switch
       procedure, non_overridable, private    :: print_logical_switch
       procedure, non_overridable             :: is_in_fe_space
       procedure, non_overridable             :: are_checks_active
  end type par_test_poisson_unfitted_params_t

  ! Types
  public :: par_test_poisson_unfitted_params_t

contains

  !==================================================================================================
  subroutine par_test_poisson_unfitted_params_define_parameters()
    implicit none
    ! common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! specific
    call parameter_handler%add(coarse_fe_handler_type_key, '--coarse-fe-handler', unfitted_coarse_fe_handler_value, 'Type of coarse fe handler. `'//&
      standard_coarse_fe_handler_value//'`, `'//&
      unfitted_coarse_fe_handler_value//'` or `'//&
      stiffness_coarse_fe_handler_value//'` ?', switch_ab='-chandler')
    call parameter_handler%add(level_set_function_type_key, '--level-set-function', level_set_sphere_str, 'Type of levelset to be used.'//&
      ' The possible values are the public character constants defined in the `level_set_functions_gallery_names` module', switch_ab='-levelset')
    call parameter_handler%add(use_preconditioner_key, '--use-preconditioner', .true., 'Use (T) or not (F) a preconditioner', switch_ab='-precond')
    call parameter_handler%add(unfitted_boundary_type_key, '--unfitted-boundary', 'dirichlet', 'Use (dirichlet) or not (neumann) boundary conditions on the unfitted boundary', switch_ab='-uboundary')
    call parameter_handler%add(nitsche_beta_factor_key, '--nitsche-beta', 2.0_rp, 'Set the value of the factor to compute nitches beta', switch_ab='-beta')
    call parameter_handler%add(levelset_tolerance_key, '--level-set-tol',  1.0e-2_rp, 'Set the tolerance of the levelset', switch_ab='-levelsettol')
    call parameter_handler%add(num_runs_key, '--number-runs', 1, 'Number of times the simulation is repeated (useful when measuring times)', switch_ab='-nruns')
    call parameter_handler%add(is_in_fe_space_key, '--solution_in_fe_space', .true., 'Is the solution in fe space', switch_ab='-in_space')
    call parameter_handler%add(are_checks_active_key, '--check_solution', .true., 'Check or not the solution', '-check')
  end subroutine par_test_poisson_unfitted_params_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_poisson_unfitted_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                      :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable                     :: get_dir_path
    call parameter_handler%GetAsString(key = dir_path_key, string = get_dir_path)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable                     :: get_prefix
    call parameter_handler%GetAsString(key = prefix_key, string = get_prefix)
  end function get_prefix
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                            :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                                :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution

  !==================================================================================================
  function get_coarse_fe_handler_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable                     :: get_coarse_fe_handler_type
    call parameter_handler%GetAsString(key = coarse_fe_handler_type_key, string = get_coarse_fe_handler_type)
  end function get_coarse_fe_handler_type

  !==================================================================================================
  function get_level_set_function_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable                     :: get_level_set_function_type
    call parameter_handler%GetAsString(key = level_set_function_type_key, string = get_level_set_function_type)
  end function get_level_set_function_type

  !==================================================================================================
  function get_use_preconditioner(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                                :: get_use_preconditioner
    call parameter_handler%Get(key = use_preconditioner_key, Value = get_use_preconditioner)
  end function get_use_preconditioner

  !==================================================================================================
  function get_unfitted_boundary_type(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=:),      allocatable                     :: get_unfitted_boundary_type
    call parameter_handler%GetAsString(key = unfitted_boundary_type_key, string = get_unfitted_boundary_type)
  end function get_unfitted_boundary_type

  !==================================================================================================
  function get_nitsche_beta_factor(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    real(rp)                                               :: get_nitsche_beta_factor
    call parameter_handler%Get(key = nitsche_beta_factor_key, Value = get_nitsche_beta_factor)
  end function get_nitsche_beta_factor

  !==================================================================================================
  function get_levelset_tolerance(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    real(rp)                                               :: get_levelset_tolerance
    call parameter_handler%Get(key = levelset_tolerance_key, Value = get_levelset_tolerance)
  end function get_levelset_tolerance

  !==================================================================================================
  function get_num_runs(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    integer(ip)                                            :: get_num_runs
    call parameter_handler%Get(key = num_runs_key, Value = get_num_runs)
  end function get_num_runs

  !==================================================================================================
  subroutine print(this, environment)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    class(environment_t) :: environment

    if (environment%get_l1_rank()==0) then
      call this%print_character_switch(coarse_fe_handler_type_key )
      call this%print_character_switch(level_set_function_type_key)
      call this%print_logical_switch  (use_preconditioner_key     )
      call this%print_character_switch(unfitted_boundary_type_key )
      call this%print_real_switch     (nitsche_beta_factor_key    )
      call this%print_real_switch     (levelset_tolerance_key     )
      call this%print_integer_switch  (num_runs_key            )
    end if

  end subroutine print

  !==================================================================================================
  subroutine print_character_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    character(len=:), allocatable                          :: val
    call parameter_handler%GetAsString(key = switch_key, string = val     )
    write(*,'(a30,a20)') switch_key, val
  end subroutine print_character_switch

  !==================================================================================================
  subroutine print_integer_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    integer(ip)                                            :: val
    call parameter_handler%Get(key = switch_key, value  = val     )
    write(*,'(a30,i20)') switch_key, val
  end subroutine print_integer_switch

  !==================================================================================================
  subroutine print_real_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    real(rp)                                               :: val
    call parameter_handler%Get(key = switch_key, value  = val     )
    write(*,'(a30,e20.5)') switch_key, val
  end subroutine print_real_switch

  !==================================================================================================
  subroutine print_logical_switch(this,switch_key)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    character(len=*),                           intent(in) :: switch_key
    logical                                                :: val
    call parameter_handler%Get(key = switch_key, value  = val)
    write(*,'(a30,l20)') switch_key, val
  end subroutine print_logical_switch

  !==================================================================================================
  function is_in_fe_space(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                                :: is_in_fe_space
    call parameter_handler%Get(key = is_in_fe_space_key, Value = is_in_fe_space)
  end function is_in_fe_space

  !==================================================================================================
  function are_checks_active(this)
    implicit none
    class(par_test_poisson_unfitted_params_t) , intent(in) :: this
    logical                                                :: are_checks_active
    call parameter_handler%Get(key = are_checks_active_key, Value = are_checks_active)
  end function are_checks_active

end module par_test_poisson_unfitted_params_names
