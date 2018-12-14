module par_test_interpolators_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key      = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key          = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key              = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key          = 'triangulation_type'    
  character(len=*), parameter :: variable_degree_within_field_key  = 'variable_degree_within_field'
  
  type, extends(parameter_handler_t) :: par_test_interpolators_params_t
     private
     contains
       procedure :: define_parameters  => par_test_interpolators_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_is_variable_degree_within_field
  end type par_test_interpolators_params_t

  ! Types
  public :: par_test_interpolators_params_t

contains

  !==================================================================================================
  subroutine par_test_interpolators_params_define_parameters(this)
    implicit none
    class(par_test_interpolators_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key            , value = '.')      ; check(error==0)
    error = list%set(key = prefix_key              , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key        , value = '.')      ; check(error==0)
    error = list%set(key = struct_hex_triang_num_dims_key            , value =  2)       ; check(error==0)
    !error = list%set(key = struct_hex_triang_domain_limits_key        , value =  [0,1,0,1,0,1])       ; check(error==0)      
    error = list%set(key = struct_hex_triang_num_cells_dir     , value =  [12,12,12])          ; check(error==0)
    error = list%set(key = struct_hex_triang_is_dir_periodic_key     , value =  [0,0,0])             ; check(error==0)
    error = list%set(key = struct_hex_triang_num_levels_key          , value =  3)                   ; check(error==0)
    error = list%set(key = struct_hex_triang_num_parts_x_dir_key     , value =  [4,4,0,2,2,0,1,1,0]) ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triang_generate_key        , value =  triangulation_generate_from_mesh) ; check(error==0)
    error = list%set(key = variable_degree_within_field_key , value =  .false.)                      ; check(error==0)
 
    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                 ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                   ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')             ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_dims_key                  , value = '--dim')                      ; check(error==0)
 !error = switches%set(key = struct_hex_triang_domain_limits_key    , value = '--domain_limits')            ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_cells_dir           , value = '--number_of_cells')          ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_levels_key                , value = '--number_of_levels')         ; check(error==0)
    error = switches%set(key = struct_hex_triang_num_parts_x_dir_key   , value = '--number_of_parts_per_dir')  ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')   ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order'    )   ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution'        )   ; check(error==0)
    error = switches%set(key = triang_generate_key    , value = '--trinagulation-type'    )   ; check(error==0)
    error = switches%set(key = variable_degree_within_field_key, value = '--variable_degree'    )    ; check(error==0)
                                                             
    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = struct_hex_triang_num_dims_key               , value = '-dm')       ; check(error==0)
 !error = switches_ab%set(key = struct_hex_triang_domain_limits_key , value = '-dl')       ; check(error==0)
    error = switches_ab%set(key = struct_hex_triang_num_cells_dir        , value = '-n')        ; check(error==0) 
    error = switches_ab%set(key = struct_hex_triang_num_levels_key             , value = '-l')        ; check(error==0)
    error = switches_ab%set(key = struct_hex_triang_num_parts_x_dir_key, value = '-np')       ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key     , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triang_generate_key , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = variable_degree_within_field_key      , value = '-variable_degree')       ; check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')            ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                    ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                         ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_dims_key                   , value = 'Number of space dimensions')               ; check(error==0)
 !error = helpers%set(key = struct_hex_triang_domain_limits_key     , value = 'Domain limits of the mesh')                ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_cells_dir            , value = 'Number of cells per dir')                  ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_levels_key                 , value = 'Number of levels')                         ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_parts_x_dir_key            , value = 'Number of parts per dir and per level')    ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')  ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')       ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')             ; check(error==0)
    error = helpers%set(key = variable_degree_within_field_key , value = 'Different order Reference_fes are considered for the same field'); check(error==0)
    
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triang_generate_key     , value = msg)  ; check(error==0)
    
    error = required%set(key = dir_path_key                  , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                    , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key              , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_dims_key                  , value = .false.) ; check(error==0)
 !error = required%set(key = struct_hex_triang_domain_limits_key    , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_cells_dir           , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_levels_key                , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_parts_x_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key       , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key           , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key               , value = .false.) ; check(error==0)
    error = required%set(key = triang_generate_key       , value = .false.) ; check(error==0)
    error = required%set(key = variable_degree_within_field_key , value = .false.) ; check(error==0)

  end subroutine par_test_interpolators_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
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
    class(par_test_interpolators_params_t) , intent(in) :: this
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
    class(par_test_interpolators_params_t) , intent(in) :: this
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
    class(par_test_interpolators_params_t) , intent(in) :: this
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
    class(par_test_interpolators_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_is_variable_degree_within_field(this)
    implicit none
    class(par_test_interpolators_params_t) , intent(in) :: this
    logical                                       :: get_is_variable_degree_within_field
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    logical                                       :: is_present
    logical                                       :: same_data_type
    integer(ip), allocatable                      :: shape(:)
    list  => this%get_values()
    assert(list%isAssignable(variable_degree_within_field_key, get_is_variable_degree_within_field))
    error = list%Get(key = variable_degree_within_field_key, Value = get_is_variable_degree_within_field)
    assert(error==0)
  end function get_is_variable_degree_within_field

end module par_test_interpolators_params_names
