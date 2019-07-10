module par_test_interpolators_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key          = 'reference_fe_order'    
  
  type :: par_test_interpolators_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_reference_fe_order
  end type par_test_interpolators_params_t

  ! Types
  public :: par_test_interpolators_params_t

contains

  !==================================================================================================
  subroutine par_test_interpolators_params_define_parameters()
    implicit none
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
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
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    integer(ip)                                         :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order

end module par_test_interpolators_params_names
