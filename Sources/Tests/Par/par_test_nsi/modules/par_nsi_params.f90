module par_nsi_params_names
  use fempar_names
  use nsi_discrete_integration_names
  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key        = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key            = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key        = 'triangulation_type'    
  character(len=*), parameter :: discrete_integration_type_key = 'discrete_integration_type' 
  character(len=*), parameter :: viscosity_key                 = 'viscosity'  

  type :: par_nsi_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_orders
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_viscosity
  end type par_nsi_params_t

  ! Types
  public :: par_nsi_params_t

contains

  !==================================================================================================
  subroutine par_nsi_params_define_parameters()
    implicit none
    integer(ip)                            :: default_reference_fe_orders(2)
    default_reference_fe_orders = [2, 1]
    ! Common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', default_reference_fe_orders, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! Specific
    call parameter_handler%add(discrete_integration_type_key, '--discrete-integration-type', discrete_integration_type_galerkin, &
                  'Discrete formulation (only Galerkin)', switch_ab='-di' )
    call parameter_handler%add(viscosity_key, '--viscosity', 1.0_rp, 'Viscosity', switch_ab='-visco' )
  end subroutine par_nsi_params_define_parameters

  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_nsi_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_nsi_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    type(ParameterList_t), pointer       :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    character(len=:),      allocatable   :: get_dir_path
    call parameter_handler%GetAsString(key = dir_path_key, string = get_dir_path)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    character(len=:),      allocatable   :: get_prefix
    call parameter_handler%GetAsString(key = prefix_key, string = get_prefix)
  end function get_prefix
  
  !==================================================================================================
  function get_reference_fe_orders(this, ifield)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    integer(ip)             , intent(in) :: ifield
    integer(ip)                          :: reference_fe_orders(2), get_reference_fe_orders
    call parameter_handler%Get(key = reference_fe_order_key, Value = reference_fe_orders)
    get_reference_fe_orders = reference_fe_orders(ifield)
  end function get_reference_fe_orders
  
  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    logical                              :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_triangulation_type
    call parameter_handler%GetAsString(key = static_triang_generate_from_key, string = get_triangulation_type)
  end function get_triangulation_type 
  
  !==================================================================================================
  function get_viscosity(this)
    implicit none
    class(par_nsi_params_t) , intent(in) :: this
    real(rp)                             :: get_viscosity
    call parameter_handler%Get(key = viscosity_key, value = get_viscosity)
  end function get_viscosity

end module par_nsi_params_names
