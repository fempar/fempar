module par_pb_bddc_maxwell_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key           = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key               = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key                   = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key               = 'triangulation_type'    
  character(len=*), parameter :: bddc_edge_continuity_algorithm_key   = 'bddc_edge_continuity_algorithm'
  character(len=*), parameter :: bddc_weighting_function_case_key     = 'bddc_weighting_function_case'
  character(len=*), parameter :: mass_coeff_white_key                 = 'mass_coeff_white'
  character(len=*), parameter :: curl_curl_coeff_white_key            = 'curl_curl_coeff_white '
  character(len=*), parameter :: mass_coeff_black_key                 = 'mass_coeff_black'
  character(len=*), parameter :: curl_curl_coeff_black_key            = 'curl_curl_coeff_black'
  character(len=*), parameter :: materials_distribution_case_key      = 'materials_distribution_case'
  character(len=*), parameter :: materials_coefficient_case_key       = 'materials_coefficient_case'
  character(len=*), parameter :: channels_ratio_key                   = 'channels_ratio' 
  character(len=*), parameter :: rpb_bddc_threshold_key               = 'rpb_bddc_threshold'
  character(len=*), parameter :: boundary_mass_trick_key              = 'boundary_mass_trick'
  character(len=*), parameter :: num_peaks_curl_curl_coeff_key        = 'num_peaks_curl_curl_coeff' 
  character(len=*), parameter :: num_peaks_mass_coeff_key             = 'num_peaks_mass_coeff' 
  
  character(len=*), parameter :: homogeneous       = 'homogeneous'
  character(len=*), parameter :: checkerboard      = 'checkerboard'       
  character(len=*), parameter :: channels          = 'channels'  
  character(len=*), parameter :: heterogeneous     = 'heterogeneous' 
  
  character(len=*), parameter :: unit_coefficients = 'unit_coefficients'
  character(len=*), parameter :: constant          = 'constant' 
  character(len=*), parameter :: sinusoidal        = 'sinusoidal'
  
  type, extends(parameter_handler_t) :: par_pb_bddc_maxwell_params_t
     private
     contains
       procedure :: define_parameters  => par_test_maxwell_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_mass_coeff_white 
       procedure, non_overridable             :: get_curl_curl_coeff_white  
       procedure, non_overridable             :: get_mass_coeff_black 
       procedure, non_overridable             :: get_curl_curl_coeff_black 
       procedure, non_overridable             :: get_materials_distribution_case 
       procedure, non_overridable             :: get_materials_coefficient_case
       procedure, non_overridable             :: get_channels_ratio 
       procedure, non_overridable             :: get_rpb_bddc_threshold 
       procedure, non_overridable             :: get_num_peaks_curl_curl_coeff
       procedure, non_overridable             :: get_num_peaks_mass_coeff 
       procedure, non_overridable             :: get_boundary_mass_trick 
       procedure, non_overridable             :: get_nparts 
  end type par_pb_bddc_maxwell_params_t

  ! Types
  public :: par_pb_bddc_maxwell_params_t
  
  public :: checkerboard, channels, homogeneous, unit_coefficients, heterogeneous, radial
  public :: constant, sinusoidal 

contains

  !==================================================================================================
  subroutine par_test_maxwell_params_define_parameters(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key                           , value = '.'); check(error==0)
    error = list%set(key = prefix_key                             , value = 'square'); check(error==0)
    error = list%set(key = dir_path_out_key                       , value = '.'); check(error==0)
    error = list%set(key = struct_hex_triang_num_dims_key         , value =  3); check(error==0)     
    error = list%set(key = struct_hex_triang_num_cells_dir        , value =  [8,8,8]); check(error==0)
    error = list%set(key = struct_hex_triang_is_dir_periodic_key  , value =  [0,0,0]); check(error==0)
    error = list%set(key = struct_hex_triang_num_levels_key       , value =  2); check(error==0)
    error = list%set(key = struct_hex_triang_num_parts_x_dir_key  , value =  [2,2,2,1,1,1,0,0,0]); check(error==0)
    error = list%set(key = reference_fe_geo_order_key             , value =  1); check(error==0)
    error = list%set(key = reference_fe_order_key                 , value =  1); check(error==0)
    error = list%set(key = write_solution_key                     , value =  .false.); check(error==0)
    error = list%set(key = triang_generate_key                    , value =  triangulation_generate_structured); check(error==0)
    error = list%set(key = coarse_space_use_vertices_key          , value =  .true.); check(error==0)
    error = list%set(key = coarse_space_use_edges_key             , value =  .true.); check(error==0)
    error = list%set(key = coarse_space_use_faces_key             , value =  .false.); check(error==0)
    error = list%set(key = bddc_edge_continuity_algorithm_key     , value =  tangential_average_and_first_order_moment ) ; check(error==0)
    error = list%set(key = bddc_weighting_function_case_key       , value =  cardinality ) ; check(error==0)
    error = list%set(key = mass_coeff_white_key                   , value =  1.0 ); check(error==0)
    error = list%set(key = curl_curl_coeff_white_key              , value =  1.0 ); check(error==0)
    error = list%set(key = mass_coeff_black_key                   , value =  1.0 ); check(error==0)
    error = list%set(key = curl_curl_coeff_black_key              , value =  1.0 ); check(error==0)
    error = list%set(key = materials_distribution_case_key        , value = checkerboard); check(error==0) 
    error = list%set(key = materials_coefficient_case_key         , value = constant); check(error==0) 
    error = list%set(key = channels_ratio_key                     , value =  0.1 ); check(error==0)
    error = list%set(key = num_peaks_curl_curl_coeff_key          , value =  3)   ; check(error==0)
    error = list%set(key = num_peaks_mass_coeff_key               , value =  3)   ; check(error==0)
    error = list%set(key = rpb_bddc_threshold_key                 , value = 10.0 ); check(error==0)
    error = list%set(key = boundary_mass_trick_key                , value =  .false.); check(error==0)
    
    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                       , value = '--dir-path'); check(error==0)
    error = switches%set(key = prefix_key                         , value = '--prefix'); check(error==0)
    error = switches%set(key = dir_path_out_key                   , value = '--dir-path-out'); check(error==0)
    error = switches%set(key = struct_hex_triang_num_dims_key     , value = '--dim'); check(error==0)
    error = switches%set(key = struct_hex_triang_num_cells_dir    , value = '--number_of_cells'); check(error==0)
    error = switches%set(key = struct_hex_triang_num_levels_key   , value = '--number_of_levels'); check(error==0)
    error = switches%set(key = struct_hex_triang_num_parts_x_dir_key , value = '--number_of_parts_per_dir')  ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key         , value = '--reference-fe-geo-order')   ; check(error==0)
    error = switches%set(key = reference_fe_order_key             , value = '--reference-fe-order'    )   ; check(error==0)
    error = switches%set(key = write_solution_key                 , value = '--write-solution'        )   ; check(error==0)
    error = switches%set(key = triang_generate_key                , value = '--triangulation-type'    )   ; check(error==0)
    error = switches%set(key = coarse_space_use_vertices_key      , value = '--coarse-space-use-vertices'); check(error==0)
    error = switches%set(key = coarse_space_use_edges_key         , value = '--coarse-space-use-edges' )  ; check(error==0)
    error = switches%set(key = coarse_space_use_faces_key         , value = '--coarse-space-use-faces' )  ; check(error==0)
    error = switches%set(key = bddc_edge_continuity_algorithm_key , value = '--BDDC_edge_continuity_algorithm' ) ; check(error==0)
    error = switches%set(key = bddc_weighting_function_case_key   , value = '--BDDC_weighting_function_case' ) ; check(error==0)
    error = switches%set(key = mass_coeff_white_key               , value = '--mass_coeff_white' )  ; check(error==0)
    error = switches%set(key = curl_curl_coeff_white_key          , value = '--curl_curl_coeff_white ' )  ; check(error==0)
    error = switches%set(key = mass_coeff_black_key               , value = '--mass_coeff_black' )  ; check(error==0)
    error = switches%set(key = curl_curl_coeff_black_key          , value = '--curl_curl_coeff_black' )  ; check(error==0)
    error = switches%set(key = materials_distribution_case_key    , value = '--materials_distribution_case' )  ; check(error==0)
    error = switches%set(key = materials_coefficient_case_key     , value = '--materials_coefficient_case' )  ; check(error==0)
    error = switches%set(key = channels_ratio_key                 , value = '--channels_ratio' )  ; check(error==0)
    error = switches%set(key = num_peaks_curl_curl_coeff_key      , value = '--num_peaks_curl_curl_coeff' )  ; check(error==0)
    error = switches%set(key = num_peaks_mass_coeff_key           , value = '--num_peaks_mass_coeff' )  ; check(error==0)
    error = switches%set(key = rpb_bddc_threshold_key             , value = '--rpb_bddc_threshold' )  ; check(error==0)
    error = switches%set(key = boundary_mass_trick_key            , value = '--boundary_mass_trick' )  ; check(error==0)
                                                             
    error = switches_ab%set(key = dir_path_key                    , value = '-d'); check(error==0) 
    error = switches_ab%set(key = prefix_key                      , value = '-p'); check(error==0) 
    error = switches_ab%set(key = dir_path_out_key                , value = '-o'); check(error==0) 
    error = switches_ab%set(key = struct_hex_triang_num_dims_key  , value = '-dm') ; check(error==0)
    error = switches_ab%set(key = struct_hex_triang_num_cells_dir , value = '-n') ; check(error==0) 
    error = switches_ab%set(key = struct_hex_triang_num_levels_key, value = '-l') ; check(error==0)
    error = switches_ab%set(key = struct_hex_triang_num_parts_x_dir_key, value = '-np') ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key      , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key          , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key              , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triang_generate_key             , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_vertices_key   , value = '-use-vertices'); check(error==0)
    error = switches_ab%set(key = coarse_space_use_edges_key      , value = '-use-edges' )  ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_faces_key      , value = '-use-faces' )  ; check(error==0)
    error = switches_ab%set(key = bddc_edge_continuity_algorithm_key, value = '-edge_cont' )  ; check(error==0)
    error = switches_ab%set(key = bddc_weighting_function_case_key, value = '-bddc_weights' )  ; check(error==0)
    error = switches_ab%set(key = mass_coeff_white_key            , value = '-mass_coeff_white' )  ; check(error==0)
    error = switches_ab%set(key = curl_curl_coeff_white_key       , value = '-curl_curl_coeff_white ' )  ; check(error==0)
    error = switches_ab%set(key = mass_coeff_black_key            , value = '-mass_coeff_black' )  ; check(error==0)
    error = switches_ab%set(key = curl_curl_coeff_black_key       , value = '-curl_curl_coeff_black' )  ; check(error==0)
    error = switches_ab%set(key = materials_distribution_case_key , value = '-materials_case' )  ; check(error==0)
    error = switches_ab%set(key = materials_coefficient_case_key  , value = '-coefficient_case' )  ; check(error==0)
    error = switches_ab%set(key = channels_ratio_key              , value = '-channels_ratio' )  ; check(error==0)
    error = switches_ab%set(key = num_peaks_curl_curl_coeff_key   , value = '-num_peaks_curl_curl_coeff' )  ; check(error==0)
    error = switches_ab%set(key = num_peaks_mass_coeff_key        , value = '-num_peaks_mass_coeff' )  ; check(error==0)
    error = switches_ab%set(key = rpb_bddc_threshold_key          , value = '-rpb_bddc_threshold' )  ; check(error==0)
    error = switches_ab%set(key = boundary_mass_trick_key         , value = '-bmass_trick' )  ; check(error==0)

    error = helpers%set(key = dir_path_key                        , value = 'Directory of the source files') ; check(error==0)
    error = helpers%set(key = prefix_key                          , value = 'Name of the GiD files'); check(error==0)
    error = helpers%set(key = dir_path_out_key                    , value = 'Output Directory') ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_dims_key      , value = 'Number of space dimensions'); check(error==0)
    error = helpers%set(key = struct_hex_triang_num_cells_dir     , value = 'Number of cells per dir') ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_levels_key    , value = 'Number of levels') ; check(error==0)
    error = helpers%set(key = struct_hex_triang_num_parts_x_dir_key, value = 'Number of parts per dir and per level') ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key          , value = 'Order of the triangulation reference fe'); check(error==0)
    error = helpers%set(key = reference_fe_order_key              , value = 'Order of the fe space reference fe') ; check(error==0)
    error = helpers%set(key = write_solution_key                  , value = 'Write solution in VTK format') ; check(error==0)
    error = helpers%set(key = coarse_space_use_vertices_key       , value  = 'Include vertex coarse DoFs in coarse FE space'); check(error==0)
    error = helpers%set(key = coarse_space_use_edges_key          , value  = 'Include edge coarse DoFs in coarse FE space' ); check(error==0)
    error = helpers%set(key = coarse_space_use_faces_key          , value  = 'Include face coarse DoFs in coarse FE space' ); check(error==0)
    
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triang_generate_key     , value = msg)  ; check(error==0)
 
    msg = 'Specify BDDC space continuity: tangential_average, tangential_average_and_first_order_moment, all_dofs_in_coarse_edges' 
    error = helpers%set(key = bddc_edge_continuity_algorithm_key  , value = msg)  ; check(error==0)
    msg = 'Define BDDC weighting function from: cardinality (inverse of the cardinality of each dof), curl_curl_coeff, mass_coeff, stiffness (diagonal entries of the operator).'
    error = helpers%set(key = bddc_weighting_function_case_key, value = msg  ); check(error==0) 
                        
    error = helpers%set(key = mass_coeff_white_key                 , value  = 'mass_coeff_white value' ) ; check(error==0)
    error = helpers%set(key = curl_curl_coeff_white_key            , value  = 'curl_curl_coeff_white  value' )  ; check(error==0)
    error = helpers%set(key = mass_coeff_black_key                 , value  = 'mass_coeff_black value' ) ; check(error==0)
    error = helpers%set(key = curl_curl_coeff_black_key            , value  = 'curl_curl_coeff_black value' )  ; check(error==0)
    error = helpers%set(key = materials_distribution_case_key      , value  = 'Materials distribution case: choose between: checkerboard, channels, radial, heterogeneous' )  ; check(error==0)
    error = helpers%set(key = materials_coefficient_case_key       , value  = 'Materials coefficient case: choose between: constant, sinusoidal' )  ; check(error==0)
    error = helpers%set(key = channels_ratio_key                   , value  = 'Ratio channel/non-channel of the cross section for every direction)' ) ; check(error==0)
    error = helpers%set(key = num_peaks_curl_curl_coeff_key        , value  = 'Number of peaks for the sinusoidal function describing the curl_curl_coeff' ) ; check(error==0)
    error = helpers%set(key = num_peaks_mass_coeff_key             , value  = 'Number of peaks for the sinusoidal function describing the mass_coeff' ) ; check(error==0)
    error = helpers%set(key = rpb_bddc_threshold_key               , value  = 'Threshold for the relaxed PB-BDDC subparts partition' ) ; check(error==0)
    error = helpers%set(key = boundary_mass_trick_key              , value  = 'Is the boundary mass trick active?' ); check(error==0)
    
    error = required%set(key = dir_path_key                        , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                          , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key                    , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_dims_key      , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_cells_dir     , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_levels_key    , value = .false.) ; check(error==0)
    error = required%set(key = struct_hex_triang_num_parts_x_dir_key, value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key          , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key              , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key                  , value = .false.) ; check(error==0)
    error = required%set(key = triang_generate_key                 , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_vertices_key       , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_edges_key          , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_faces_key          , value = .false.) ; check(error==0)
    error = required%set(key = bddc_edge_continuity_algorithm_key  , value = .false.) ; check(error==0)
    error = required%set(key = bddc_weighting_function_case_key    , value = .false.) ; check(error==0)
    error = required%set(key = mass_coeff_white_key                , value = .false.) ; check(error==0)
    error = required%set(key = curl_curl_coeff_white_key           , value = .false.) ; check(error==0)
    error = required%set(key = mass_coeff_black_key                , value = .false.) ; check(error==0)
    error = required%set(key = curl_curl_coeff_black_key           , value = .false.) ; check(error==0)
    error = required%set(key = materials_distribution_case_key     , value = .false.) ; check(error==0)
    error = required%set(key = materials_coefficient_case_key      , value = .false.) ; check(error==0)
    error = required%set(key = channels_ratio_key                  , value = .false.) ; check(error==0)
    error = required%set(key = num_peaks_curl_curl_coeff_key       , value = .false. ) ; check(error==0) 
    error = required%set(key = num_peaks_mass_coeff_key            , value = .false. ) ; check(error==0)
    error = required%set(key = rpb_bddc_threshold_key              , value = .false.) ; check(error==0)
    error = required%set(key = boundary_mass_trick_key             , value = .false.) ; check(error==0)

  end subroutine par_test_maxwell_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
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
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix
  
  ! =======================================================================================
   function get_nparts(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t)           :: this
    integer(ip)                                   :: num_levels
    integer(ip)                                   :: get_nparts(3)
    integer(ip), allocatable :: num_parts_x_dir(:) 
    integer(ip), allocatable :: array_size(:)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_num_levels_key , num_levels))
    error = list%Get(key = struct_hex_triang_num_levels_key , Value = num_levels)
    assert(error==0)       
    error = list%GetShape(key = struct_hex_triang_num_parts_x_dir_key   , shape = array_size); 
    check(error==0)
    assert(array_size(1) >= num_levels*SPACE_DIM)
    call memalloc(array_size(1), num_parts_x_dir)
    error = list%get(key = struct_hex_triang_num_parts_x_dir_key , value = num_parts_x_dir) 
    check(error==0)
    get_nparts=num_parts_x_dir(1:3)
    if (allocated(array_size)) deallocate(array_size) 
    call memfree(num_parts_x_dir)

  end function get_nparts
  
    !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
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
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
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
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
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
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_mass_coeff_white(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_mass_coeff_white
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(mass_coeff_white_key, get_mass_coeff_white))
    error = list%Get(key = mass_coeff_white_key, Value = get_mass_coeff_white)
    assert(error==0)
  end function get_mass_coeff_white
  
      !==================================================================================================
  function get_mass_coeff_black(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_mass_coeff_black
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(mass_coeff_black_key, get_mass_coeff_black))
    error = list%Get(key = mass_coeff_black_key, Value = get_mass_coeff_black)
    assert(error==0)
  end function get_mass_coeff_black
  
     !==================================================================================================
  function get_curl_curl_coeff_white (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_curl_curl_coeff_white 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(curl_curl_coeff_white_key, get_curl_curl_coeff_white ))
    error = list%Get(key = curl_curl_coeff_white_key, Value = get_curl_curl_coeff_white )
    assert(error==0)
  end function get_curl_curl_coeff_white 
  
       !==================================================================================================
  function get_curl_curl_coeff_black (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_curl_curl_coeff_black 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(curl_curl_coeff_black_key, get_curl_curl_coeff_black ))
    error = list%Get(key = curl_curl_coeff_black_key, Value = get_curl_curl_coeff_black)
    assert(error==0)
  end function get_curl_curl_coeff_black 
  
      !==================================================================================================
  function get_materials_distribution_case(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                    :: get_materials_distribution_case
    type(ParameterList_t), pointer                   :: list
    integer(ip)                                      :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(materials_distribution_case_key, dummy_string))
    error = list%GetAsString(key = materials_distribution_case_key, string = get_materials_distribution_case)
    assert(error==0)
  end function get_materials_distribution_case
  
        !==================================================================================================
  function get_materials_coefficient_case(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                    :: get_materials_coefficient_case
    type(ParameterList_t), pointer                   :: list
    integer(ip)                                      :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(materials_coefficient_case_key, dummy_string))
    error = list%GetAsString(key = materials_coefficient_case_key, string = get_materials_coefficient_case)
    assert(error==0)
  end function get_materials_coefficient_case

       !==================================================================================================
  function get_channels_ratio (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_channels_ratio 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(channels_ratio_key, get_channels_ratio ))
    error = list%Get(key = channels_ratio_key, Value = get_channels_ratio )
    assert(error==0)
  end function get_channels_ratio
  
      !==================================================================================================
  function get_rpb_bddc_threshold(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_rpb_bddc_threshold
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(rpb_bddc_threshold_key, get_rpb_bddc_threshold))
    error = list%Get(key = rpb_bddc_threshold_key, Value = get_rpb_bddc_threshold)
    assert(error==0)
  end function get_rpb_bddc_threshold
  
    !==================================================================================================
  function get_boundary_mass_trick(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t), intent(in) :: this
    logical                                       :: get_boundary_mass_trick
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(boundary_mass_trick_key, get_boundary_mass_trick))
    error = list%Get(key = boundary_mass_trick_key, Value = get_boundary_mass_trick)
    assert(error==0)
  end function get_boundary_mass_trick
  
    !==================================================================================================
  function get_num_peaks_curl_curl_coeff(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_num_peaks_curl_curl_coeff
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(num_peaks_curl_curl_coeff_key, get_num_peaks_curl_curl_coeff))
    error = list%Get(key = num_peaks_curl_curl_coeff_key, Value = get_num_peaks_curl_curl_coeff)
    assert(error==0)
  end function get_num_peaks_curl_curl_coeff
  
      !==================================================================================================
  function get_num_peaks_mass_coeff(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_num_peaks_mass_coeff
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(num_peaks_mass_coeff_key, get_num_peaks_mass_coeff))
    error = list%Get(key = num_peaks_mass_coeff_key, Value = get_num_peaks_mass_coeff)
    assert(error==0)
  end function get_num_peaks_mass_coeff
  
end module par_pb_bddc_maxwell_params_names
