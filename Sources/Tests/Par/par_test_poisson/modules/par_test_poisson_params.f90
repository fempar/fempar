module par_test_poisson_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'
  character(len=*), parameter :: use_void_fes_key           = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key      = 'use_void_fes_case'
  character(len=*), parameter :: preconditioner_type_key    = 'preconditioner_type'

  type :: par_test_poisson_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_use_void_fes
       procedure, non_overridable             :: get_use_void_fes_case
       procedure, non_overridable             :: get_preconditioner_type
       !procedure, non_overridable             :: get_num_dims
  end type par_test_poisson_params_t

  ! Types
  public :: par_test_poisson_params_t
  
contains

  !==================================================================================================
  subroutine par_test_poisson_params_define_parameters()
    implicit none
    ! Common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! Overwritten value
    call parameter_handler%update(static_triang_generate_from_key, static_triang_generate_from_mesh_data_files)

    ! Specific
    call parameter_handler%add(use_void_fes_key, '--use-void-fes', .false., 'Use a hybrid FE space formed by full and void FEs', '-use-voids')
    call parameter_handler%add(use_void_fes_case_key, '--use-void-fes-case', 'popcorn', 'Select where to put void fes using one of the predefined patterns. Possible values: `popcorn`, `half`, `quarter` ', switch_ab='-use-voids-case')
    call parameter_handler%add(preconditioner_type_key, '--preconditioner-type', 'mlbddc', 'Select preconditioner type. Possible values: `identity`, `jacobi`, `mlbddc`.', switch_ab='-prec-type')

  end subroutine par_test_poisson_params_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_poisson_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                      :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    logical                                       :: is_present
    logical                                       :: same_data_type
    integer(ip), allocatable                      :: shape(:)
    list  => parameter_handler%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_triangulation_type
    call parameter_handler%GetAsString(key = static_triang_generate_from_key, string = get_triangulation_type)
  end function get_triangulation_type 

  !==================================================================================================
  function get_use_void_fes(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    logical                                       :: get_use_void_fes
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(use_void_fes_key, get_use_void_fes))
    error = list%Get(key = use_void_fes_key, Value = get_use_void_fes)
    assert(error==0)
  end function get_use_void_fes

  !==================================================================================================
  function get_use_void_fes_case(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_use_void_fes_case
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(use_void_fes_case_key, 'string'))
    error = list%GetAsString(key = use_void_fes_case_key, string = get_use_void_fes_case)
    assert(error==0)
  end function get_use_void_fes_case
  
  !==================================================================================================
  function get_preconditioner_type(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_preconditioner_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
    assert(list%isAssignable(preconditioner_type_key, 'string'))
    error = list%GetAsString(key = preconditioner_type_key, string = get_preconditioner_type)
    assert(error==0)
  end function get_preconditioner_type

end module par_test_poisson_params_names
