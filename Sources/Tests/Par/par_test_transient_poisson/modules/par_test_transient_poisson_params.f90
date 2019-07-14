module par_test_transient_poisson_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key      = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key          = 'write_solution'        
  character(len=*), parameter :: use_void_fes_key            = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key       = 'use_void_fes_case'
  character(len=*), parameter :: initial_time_key            = 'initial_time'
  character(len=*), parameter :: final_time_key              = 'final_time'
  character(len=*), parameter :: time_step_key               = 'time_step'
  character(len=*), parameter :: time_integration_scheme_key = 'time_integration_scheme'
  character(len=*), parameter :: is_test_key                 = 'is_test'
  
  type :: par_test_transient_poisson_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_use_void_fes
       procedure, non_overridable             :: get_use_void_fes_case
       procedure, non_overridable             :: get_initial_time
       procedure, non_overridable             :: get_final_time
       procedure, non_overridable             :: get_time_step
       procedure, non_overridable             :: get_time_integration_scheme
       procedure, non_overridable             :: get_is_test
  end type par_test_transient_poisson_params_t

  ! Types
  public :: par_test_transient_poisson_params_t

contains

  !==================================================================================================
  subroutine par_test_transient_poisson_params_define_parameters()
    implicit none

    ! Common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! Specific
    call parameter_handler%add(use_void_fes_key, '--use-void-fes', .false., 'Use a hybrid FE space formed by full and void FEs', switch_ab='-use-voids')
    call parameter_handler%add(use_void_fes_case_key, '--use-void-fes-case', 'popcorn', 'Select where to put void fes using one of the predefined patterns. Possible values: `popcorn`, `half`, `quarter` ', switch_ab='-use-voids-case')
    call parameter_handler%add(initial_time_key, '--initial-time', 0.0_rp, 'Initial time: t0', switch_ab='-t0')
    call parameter_handler%add(final_time_key, '--final-time', 1.0_rp, 'Final time: tf', switch_ab='-tf')
    call parameter_handler%add(time_step_key, '--time-step', 1.0_rp, 'Time step size: dt', switch_ab='-dt')
    call parameter_handler%add(time_integration_scheme_key, '--time-integration-schem', 'backward_euler', 'Time disctetization scheme of the DIRK solver.', switch_ab='-rk-scheme')
    call parameter_handler%add(is_test_key, '--is-test', .false., 'Test the convergence order', switch_ab='-test')
  end subroutine par_test_transient_poisson_params_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_transient_poisson_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                      :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    integer(ip)                                             :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    logical                                                 :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution

  !==================================================================================================
  function get_use_void_fes(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    logical                                                 :: get_use_void_fes
    call parameter_handler%Get(key = use_void_fes_key, Value = get_use_void_fes)
  end function get_use_void_fes

  !==================================================================================================
  function get_use_void_fes_case(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                           :: get_use_void_fes_case
    call parameter_handler%GetAsString(key = use_void_fes_case_key, string = get_use_void_fes_case)
  end function get_use_void_fes_case
  
  !==================================================================================================
  function get_initial_time(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                                :: get_initial_time
    call parameter_handler%Get(key = initial_time_key, Value = get_initial_time)
  end function get_initial_time
  
    !==================================================================================================
  function get_final_time(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                                :: get_final_time
    call parameter_handler%Get(key = final_time_key, Value = get_final_time)
  end function get_final_time
  
    !==================================================================================================
  function get_time_step(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    real(rp)                                                :: get_time_step
    call parameter_handler%Get(key = time_step_key, Value = get_time_step)
  end function get_time_step
  
    !==================================================================================================
  function get_time_integration_scheme(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                           :: get_time_integration_scheme
    call parameter_handler%GetAsString(key = time_integration_scheme_key, string = get_time_integration_scheme)
  end function get_time_integration_scheme
    
    !==================================================================================================
  function get_is_test(this)
    implicit none
    class(par_test_transient_poisson_params_t) , intent(in) :: this
    logical                                                 :: get_is_test
    call parameter_handler%Get(key = is_test_key, Value = get_is_test)
  end function get_is_test
 
end module par_test_transient_poisson_params_names
