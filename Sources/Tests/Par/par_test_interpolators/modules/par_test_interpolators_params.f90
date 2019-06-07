module par_test_interpolators_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key          = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key              = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key          = 'triangulation_type'    
  
  type :: par_test_interpolators_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_triangulation_type
  end type par_test_interpolators_params_t

  ! Types
  public :: par_test_interpolators_params_t

contains

  !==================================================================================================
  subroutine par_test_interpolators_params_define_parameters()
    implicit none
    ! Common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')
  end subroutine par_test_interpolators_params_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_interpolators_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                      :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    character(len=:),      allocatable                  :: get_dir_path
    call parameter_handler%GetAsString(key = dir_path_key, string = get_dir_path)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    character(len=:),      allocatable                  :: get_prefix
    call parameter_handler%GetAsString(key = prefix_key, string = get_prefix)
  end function get_prefix
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    integer(ip)                                         :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    character(len=:), allocatable                       :: get_triangulation_type
    call parameter_handler%GetAsString(key = static_triang_generate_from_key, string = get_triangulation_type)
  end function get_triangulation_type 

end module par_test_interpolators_params_names
