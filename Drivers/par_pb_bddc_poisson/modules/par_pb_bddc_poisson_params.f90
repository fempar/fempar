module par_pb_bddc_poisson_params_names
  use fempar_names
  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'    
  character(len=*), parameter :: jump_key                   = 'jump'    
  character(len=*), parameter :: inclusion_key              = 'inclusion'    

  type, extends(parameter_generator_t) :: par_pb_bddc_poisson_params_t
     private
     contains
       procedure                              :: set_default  => par_pb_bddc_poisson_params_set_default
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_jump
       procedure, non_overridable             :: get_inclusion
       !procedure, non_overridable             :: get_num_dimensions
  end type par_pb_bddc_poisson_params_t

  ! Types
  public :: par_pb_bddc_poisson_params_t

contains

  !==================================================================================================
  subroutine par_pb_bddc_poisson_params_set_default(this)
    implicit none
    class(par_pb_bddc_poisson_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_parameters()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key      , value = '.') ; check(error==0)
    error = list%set(key = prefix_key        , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key  , value = '.') ; check(error==0)
    error = list%set(key = number_of_dimensions_key          , value =  2)                   ; check(error==0)
    error = list%set(key = number_of_cells_per_dir_key       , value =  [12,12,12])          ; check(error==0)
    error = list%set(key = is_dir_periodic_key               , value =  [0,0,0])             ; check(error==0)
    error = list%set(key = number_of_levels_key              , value =  3)                   ; check(error==0)
    error = list%set(key = number_of_parts_per_dir_key       , value =  [4,4,0,2,2,0,1,1,0]) ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triangulation_generate_key        , value =  triangulation_generate_from_mesh) ; check(error==0)
    error = list%set(key = execution_context_key             , value =  mpi_context)                      ; check(error==0)
    error = list%set(key = jump_key                          , value =  1)  ; check(error==0)
    error = list%set(key = inclusion_key                     , value =  1)  ; check(error==0)

    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                  ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')            ; check(error==0)
    error = switches%set(key = number_of_cells_per_dir_key   , value = '--number_of_cells')         ; check(error==0)
    error = switches%set(key = number_of_levels_key          , value = '--number_of_levels')        ; check(error==0)
    error = switches%set(key = number_of_parts_per_dir_key   , value = '--number_of_parts_per_dir') ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')  ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order')      ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution')          ; check(error==0)
    error = switches%set(key = triangulation_generate_key    , value = '--triangulation-type')      ; check(error==0)
    error = switches%set(key = execution_context_key         , value = '--execution_context')       ; check(error==0)
    error = switches%set(key = jump_key                      , value = '--jump')                    ; check(error==0)
    error = switches%set(key = inclusion_key                 , value = '--inclusion')               ; check(error==0)
                                                             
    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = number_of_cells_per_dir_key, value = '-n')        ; check(error==0) 
    error = switches_ab%set(key = number_of_levels_key       , value = '-l')        ; check(error==0)
    error = switches_ab%set(key = number_of_parts_per_dir_key, value = '-np')       ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key     , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triangulation_generate_key , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = execution_context_key      , value = '-exe')      ; check(error==0)
    error = switches_ab%set(key = jump_key                   , value = '-j')        ; check(error==0)
    error = switches_Ab%set(key = inclusion_key              , value = '-i')        ; check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')               ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                       ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                            ; check(error==0)
    error = helpers%set(key = number_of_cells_per_dir_key    , value = 'Number of cells per dir')                     ; check(error==0)
    error = helpers%set(key = number_of_levels_key           , value = 'Number of levels')                            ; check(error==0)
    error = helpers%set(key = number_of_parts_per_dir_key    , value = 'Number of parts per dir and per level')       ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')     ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')          ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')                ; check(error==0)
    error = helpers%set(key = jump_key                       , value = 'Jump of physical parameter in the inclusion') ; check(error==0)
    error = helpers%set(key = inclusion_key                  , value = 'Inclusion type')                              ; check(error==0)

    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triangulation_generate_key     , value = msg)  ; check(error==0)
    
    msg = 'serial (*) or mpi (*) context?'
    write(msg(9:9),'(i1)') serial_context
    write(msg(20:20),'(i1)') mpi_context
    error = helpers%set(key = execution_context_key     , value = msg)  ; check(error==0)

    
    error = required%set(key = dir_path_key                  , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                    , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key              , value = .false.) ; check(error==0)
    error = required%set(key = number_of_cells_per_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = number_of_levels_key          , value = .false.) ; check(error==0)
    error = required%set(key = number_of_parts_per_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key    , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key        , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key            , value = .false.) ; check(error==0)
    error = required%set(key = triangulation_generate_key    , value = .false.) ; check(error==0)
    error = required%set(key = execution_context_key         , value = .false.) ; check(error==0)
    error = required%set(key = jump_key                      , value = .false.) ; check(error==0)
    error = required%set(key = inclusion_key                 , value = .false.) ; check(error==0)

  end subroutine par_pb_bddc_poisson_params_set_default

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, dir_path_key, get_dir_path))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    check(error==0)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, prefix_key, get_prefix))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    check(error==0)
  end function get_prefix

    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_geo_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, reference_fe_geo_order_key, get_reference_fe_geo_order))
    error = list%Get(key = reference_fe_geo_order_key, Value = get_reference_fe_geo_order)
    check(error==0)
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    check(error==0)
  end function get_reference_fe_order
  
  !==========================================================================================par_pb_bddc_poisson_params_t========
  function get_write_solution(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    check(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    check(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_jump(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_jump
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, jump_key, get_jump))
    error = list%Get(key = jump_key, Value = get_jump)
    check(error==0)
  end function get_jump

  !==================================================================================================
  function get_inclusion(this)
    implicit none
    class(par_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_inclusion
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    assert(parameter_consistency(list, inclusion_key, get_inclusion))
    error = list%Get(key = inclusion_key, Value = get_inclusion)
    check(error==0)
  end function get_inclusion

end module par_pb_bddc_poisson_params_names
