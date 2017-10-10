module par_pb_bddc_linear_elasticity_params_names
  use fempar_names
  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'  
  character(len=*), parameter :: write_matrices_key         = 'write_matrices'  
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'    
  character(len=*), parameter :: jump_key                   = 'jump'    
  character(len=*), parameter :: inclusion_key              = 'inclusion'  
  character(len=*), parameter :: coarse_fe_handler_type_key = 'coarse_fe_handler_type_key' 
  character(len=*), parameter :: standard_bddc              = 'standard_bddc' 
  character(len=*), parameter :: pb_bddc                    = 'pb_bddc'
  character(len=*), parameter :: heterogeneous              = 'heterogeneous'
  character(len=*), parameter :: homogeneous                = 'homogeneous'
  character(len=*), parameter :: discrete_integration_type_key  = 'discrete_integration_type_key'
  character(len=*), parameter :: nchannel_x_direction_key = 'nchannel_x_direction' 
  character(len=*), parameter :: nparts_with_channels_key   = 'nparts_with_channels' 

  type, extends(parameter_handler_t) :: par_pb_bddc_linear_elasticity_params_t
     private
     contains
       procedure                              :: define_parameters  => par_pb_bddc_linear_elasticity_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_dir_path_out
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_write_matrices
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_jump
       procedure, non_overridable             :: get_inclusion
       procedure, non_overridable             :: get_discrete_integration_type
       procedure, non_overridable             :: get_coarse_fe_handler_type
       procedure, non_overridable             :: get_nchannel_x_direction
       procedure, non_overridable             :: get_nparts_with_channels
       procedure, non_overridable             :: get_nparts
       !procedure, non_overridable             :: get_num_dims
  end type par_pb_bddc_linear_elasticity_params_t

  ! Types
  public :: par_pb_bddc_linear_elasticity_params_t, standard_bddc, pb_bddc

contains

  !==================================================================================================
  subroutine par_pb_bddc_linear_elasticity_params_define_parameters(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key                      , value = '.') ; check(error==0)
    error = list%set(key = prefix_key                        , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key                  , value = '.') ; check(error==0)
    error = list%set(key = num_dims_key                      , value =  3)                   ; check(error==0)
    error = list%set(key = num_cells_x_dir_key               , value =  [12,12,12])          ; check(error==0)
    error = list%set(key = is_dir_periodic_key               , value =  [0,0,0])             ; check(error==0)
    error = list%set(key = num_levels_key                    , value =  3)                   ; check(error==0)
    error = list%set(key = num_parts_x_dir_key               , value =  [4,4,0,2,2,0,1,1,0]) ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = write_matrices_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triangulation_generate_key        , value =  triangulation_generate_from_mesh) ; check(error==0)
    error = list%set(key = execution_context_key             , value =  mpi_context)                      ; check(error==0)
    error = list%set(key = jump_key                          , value =  1)  ; check(error==0)
    error = list%set(key = inclusion_key                     , value =  1)  ; check(error==0)
    error = list%set(key = coarse_space_use_vertices_key     , value =  .true.)                      ; check(error==0)
    error = list%set(key = coarse_space_use_edges_key        , value =  .false.)                      ; check(error==0)
    error = list%set(key = coarse_space_use_faces_key        , value =  .true.)                      ; check(error==0)
    error = list%set(key = coarse_fe_handler_type_key        , value =  pb_bddc)                     ; check(error==0)
    error = list%set(key = discrete_integration_type_key     , value =  heterogeneous)               ; check(error==0)
    error = list%set(key = nchannel_x_direction_key          , value = [1,1,1])                      ; check(error==0)
    error = list%set(key = nparts_with_channels_key          , value = [1,1,1])                      ; check(error==0)


    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                  ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')            ; check(error==0)
    error = switches%set(key = num_dims_key                  , value = '--dim')                     ; check(error==0)
    error = switches%set(key = num_cells_x_dir_key           , value = '--num_cells')         ; check(error==0)
    error = switches%set(key = num_levels_key                , value = '--num_levels')        ; check(error==0)
    error = switches%set(key = num_parts_x_dir_key           , value = '--num_parts_x_dir') ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')  ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order')      ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution')          ; check(error==0)
    error = switches%set(key = write_matrices_key            , value = '--write-matrices')          ; check(error==0)
    error = switches%set(key = triangulation_generate_key    , value = '--triangulation-type')      ; check(error==0)
    error = switches%set(key = execution_context_key         , value = '--execution_context')       ; check(error==0)
    error = switches%set(key = jump_key                      , value = '--jump')                    ; check(error==0)
    error = switches%set(key = inclusion_key                 , value = '--inclusion')               ; check(error==0)
    error = switches%set(key = coarse_space_use_vertices_key , value = '--coarse-space-use-vertices'); check(error==0)
    error = switches%set(key = coarse_space_use_edges_key    , value = '--coarse-space-use-edges' )  ; check(error==0)
    error = switches%set(key = coarse_space_use_faces_key    , value = '--coarse-space-use-faces' )  ; check(error==0)
    error = switches%set(key = coarse_fe_handler_type_key    , value = '--coarse-fe-handler')        ; check(error==0)
    error = switches%set(key = discrete_integration_type_key , value = '--discrete_integration')     ; check(error==0)
    error = switches%set(key = nchannel_x_direction_key      , value = '--nchannel_x_direction')     ; check(error==0)
    error = switches%set(key = nparts_with_channels_key      , value = '--nparts_with_channels')     ; check(error==0)


                                                             
    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = num_dims_key               , value = '-dm')       ; check(error==0)
    error = switches_ab%set(key = num_cells_x_dir_key        , value = '-n')        ; check(error==0) 
    error = switches_ab%set(key = num_levels_key             , value = '-l')        ; check(error==0)
    error = switches_ab%set(key = num_parts_x_dir_key        , value = '-np')       ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key     , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = write_matrices_key         , value = '-wmatrices'); check(error==0)
    error = switches_ab%set(key = triangulation_generate_key , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = execution_context_key      , value = '-exe')      ; check(error==0)
    error = switches_ab%set(key = jump_key                   , value = '-j')        ; check(error==0)
    error = switches_ab%set(key = inclusion_key              , value = '-i')        ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_vertices_key , value = '-use-vertices'); check(error==0)
    error = switches_ab%set(key = coarse_space_use_edges_key    , value = '-use-edges' )  ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_faces_key    , value = '-use-faces' )  ; check(error==0)
    error = switches_ab%set(key = discrete_integration_type_key , value = '-integration') ; check(error==0)
    error = switches_ab%set(key = coarse_fe_handler_type_key    , value = '-coarse-handler')        ; check(error==0)
    error = switches_ab%set(key = nchannel_x_direction_key      , value = '-nc')        ; check(error==0)
    error = switches_ab%set(key = nparts_with_channels_key      , value = '-npwc')      ; check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')               ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                       ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                            ; check(error==0)
    error = helpers%set(key = num_dims_key                   , value = 'Number of space dimensions')               ; check(error==0)
    error = helpers%set(key = num_cells_x_dir_key            , value = 'Number of cells per dir')                     ; check(error==0)
    error = helpers%set(key = num_levels_key                 , value = 'Number of levels')                            ; check(error==0)
    error = helpers%set(key = num_parts_x_dir_key            , value = 'Number of parts per dir and per level')       ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')     ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')          ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')                ; check(error==0)
    error = helpers%set(key = write_matrices_key             , value = 'Write local-to-subdomain sparse matrices  in matrix market format')                ; check(error==0)
    error = helpers%set(key = jump_key                       , value = 'Jump of physical parameter in the inclusion') ; check(error==0)
    error = helpers%set(key = inclusion_key                  , value = 'Inclusion type')                              ; check(error==0)
    error = helpers%set(key = coarse_space_use_vertices_key  , value  = 'Include vertex coarse DoFs in coarse FE space'); check(error==0)
    error = helpers%set(key = coarse_space_use_edges_key     , value  = 'Include edge coarse DoFs in coarse FE space' )  ; check(error==0)
    error = helpers%set(key = coarse_space_use_faces_key     , value  = 'Include face coarse DoFs in coarse FE space' )  ; check(error==0)
    error = helpers%set(key = discrete_integration_type_key  , value  = 'Which type of discrete integration to use? homogeneous or heterogeneous')        ; check(error==0)
    error = helpers%set(key = coarse_fe_handler_type_key     , value  = 'Which coarse fe handler to use?')        ; check(error==0)
    error = helpers%set(key = nchannel_x_direction_key       , value  = 'Number of channels per direction')       ; check(error==0)
    error = helpers%set(key = nparts_with_channels_key       , value  = 'Number of parts per with channels')      ; check(error==0)
    
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
    error = required%set(key = num_cells_x_dir_key           , value = .false.) ; check(error==0)
    error = required%set(key = num_dims_key                  , value = .false.) ; check(error==0)
    error = required%set(key = num_levels_key                , value = .false.) ; check(error==0)
    error = required%set(key = num_parts_x_dir_key           , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key    , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key        , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key            , value = .false.) ; check(error==0)
    error = required%set(key = write_matrices_key            , value = .false.) ; check(error==0)
    error = required%set(key = triangulation_generate_key    , value = .false.) ; check(error==0)
    error = required%set(key = execution_context_key         , value = .false.) ; check(error==0)
    error = required%set(key = jump_key                      , value = .false.) ; check(error==0)
    error = required%set(key = inclusion_key                 , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_vertices_key , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_edges_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_faces_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_fe_handler_type_key    , value = .false.) ; check(error==0)
    error = required%set(key = discrete_integration_type_key , value = .false.) ; check(error==0)
    error = required%set(key = nchannel_x_direction_key      , value = .false.) ; check(error==0)
    error = required%set(key = nparts_with_channels_key      , value = .false.) ; check(error==0)


  end subroutine par_pb_bddc_linear_elasticity_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path 
  
  ! GETTERS *****************************************************************************************
  function get_dir_path_out(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path_out
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_out_key, 'string'))
    error = list%GetAsString(key = dir_path_out_key, string = get_dir_path_out)
    assert(error==0)
  end function get_dir_path_out

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==========================================================================================par_pb_bddc_linear_elasticity_params_t========
  function get_write_solution(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    check(error==0)
  end function get_write_solution
  
  !==========================================================================================par_pb_bddc_linear_elasticity_params_t========
  function get_write_matrices(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    logical                                       :: get_write_matrices
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(write_matrices_key, get_write_matrices))
    error = list%Get(key = write_matrices_key, Value = get_write_matrices)
    check(error==0)
  end function get_write_matrices

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_jump(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_jump
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(jump_key, get_jump))
    error = list%Get(key = jump_key, Value = get_jump)
    assert(error==0)
  end function get_jump

  !==================================================================================================
  function get_inclusion(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_inclusion
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(inclusion_key, get_inclusion))
    error = list%Get(key = inclusion_key, Value = get_inclusion)
    assert(error==0)
  end function get_inclusion

  !==================================================================================================  
  function get_coarse_fe_handler_type(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_coarse_fe_handler_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(coarse_fe_handler_type_key, get_coarse_fe_handler_type))
    error = list%GetAsString(key = coarse_fe_handler_type_key, string = get_coarse_fe_handler_type)
    assert(error==0)
  end function get_coarse_fe_handler_type 

  !==================================================================================================  
  function get_discrete_integration_type(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_discrete_integration_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(discrete_integration_type_key, get_discrete_integration_type))
    error = list%GetAsString(key = discrete_integration_type_key, string = get_discrete_integration_type)
    assert(error==0)
  end function get_discrete_integration_type 
  
  !==================================================================================================
  function get_nchannel_x_direction(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_nchannel_x_direction(3)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(nchannel_x_direction_key, get_nchannel_x_direction))
    error = list%Get(key = nchannel_x_direction_key, Value = get_nchannel_x_direction)
    assert(error==0)
  end function get_nchannel_x_direction

  !==================================================================================================
  function get_nparts_with_channels(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_nparts_with_channels(3)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(nparts_with_channels_key, get_nparts_with_channels))
    error = list%Get(key = nparts_with_channels_key, Value = get_nparts_with_channels)
    assert(error==0)
  end function get_nparts_with_channels

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: num_levels
    integer(ip)                                   :: get_nparts(3)
    integer(ip), allocatable :: num_parts_x_dir(:) ! 0:SPACE_DIM-1)
    integer(ip), allocatable :: array_size(:)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(num_levels_key, num_levels))
    error = list%Get(key = num_levels_key, Value = num_levels)
    assert(error==0)       
    error = list%GetShape(key = num_parts_x_dir_key   , shape = array_size); 
    check(error==0)
    assert(array_size(1) >= num_levels*SPACE_DIM)
    call memalloc(array_size(1), num_parts_x_dir)
    error = list%get(key = num_parts_x_dir_key , value = num_parts_x_dir) 
    check(error==0)
    get_nparts=num_parts_x_dir(1:3)
    if (allocated(array_size)) deallocate(array_size) 
    call memfree(num_parts_x_dir)

  end function get_nparts
  
 
end module par_pb_bddc_linear_elasticity_params_names
