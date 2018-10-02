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
  character(len=*), parameter :: permeability_white_key               = 'permeability_white'
  character(len=*), parameter :: resistivity_white_key                = 'resistivity_white '
  character(len=*), parameter :: permeability_black_key               = 'permeability_black'
  character(len=*), parameter :: resistivity_black_key                = 'resistivity_black'
  character(len=*), parameter :: materials_distribution_case_key      = 'materials_distribution_case'
  character(len=*), parameter :: channels_ratio_key                   = 'channels_ratio' 
  character(len=*), parameter :: rpb_bddc_threshold_key               = 'rpb_bddc_threshold'
  character(len=*), parameter :: boundary_mass_trick_key              = 'boundary_mass_trick'
  
  character(len=*), parameter :: homogeneous     = 'homogeneous'
  character(len=*), parameter :: checkerboard    = 'checkerboard'       
  character(len=*), parameter :: channels        = 'channels'  
  character(len=*), parameter :: heterogeneous   = 'heterogeneous' 
  
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
       procedure, non_overridable             :: get_permeability_white 
       procedure, non_overridable             :: get_resistivity_white  
       procedure, non_overridable             :: get_permeability_black 
       procedure, non_overridable             :: get_resistivity_black 
       procedure, non_overridable             :: get_materials_distribution_case 
       procedure, non_overridable             :: get_channels_ratio 
       procedure, non_overridable             :: get_rpb_bddc_threshold 
       procedure, non_overridable             :: get_boundary_mass_trick 
       procedure, non_overridable             :: get_nparts 
  end type par_pb_bddc_maxwell_params_t

  ! Types
  public :: par_pb_bddc_maxwell_params_t
  
  public :: checkerboard, channels, homogeneous, heterogeneous 

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

    error = list%set(key = dir_path_key            , value = '.')      ; check(error==0)
    error = list%set(key = prefix_key              , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key        , value = '.')      ; check(error==0)
    error = list%set(key = num_dims_key            , value =  2)       ; check(error==0)
	   !error = list%set(key = hex_mesh_domain_limits_key        , value =  [0,1,0,1,0,1])       ; check(error==0)      
    error = list%set(key = num_cells_x_dir_key     , value =  [12,12,12])          ; check(error==0)
    error = list%set(key = is_dir_periodic_key     , value =  [0,0,0])             ; check(error==0)
    error = list%set(key = num_levels_key          , value =  3)                   ; check(error==0)
    error = list%set(key = num_parts_x_dir_key     , value =  [4,4,0,2,2,0,1,1,0]) ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triangulation_generate_key        , value =  triangulation_generate_from_mesh) ; check(error==0)
    error = list%set(key = coarse_space_use_vertices_key     , value =  .true.)                                    ; check(error==0)
    error = list%set(key = coarse_space_use_edges_key        , value =  .true.)                                    ; check(error==0)
    error = list%set(key = coarse_space_use_faces_key        , value =  .false.)                                   ; check(error==0)
	   error = list%set(key = bddc_edge_continuity_algorithm_key, value =  tangential_average_and_first_order_moment ) ; check(error==0)
    error = list%set(key = bddc_weighting_function_case_key, value =  cardinality ) ; check(error==0)
    error = list%set(key = permeability_white_key   , value =  1.0 ); check(error==0)
    error = list%set(key = resistivity_white_key    , value =  1.0 ); check(error==0)
    error = list%set(key = permeability_black_key   , value =  1.0 ); check(error==0)
    error = list%set(key = resistivity_black_key    , value =  1.0 ); check(error==0)
    error = list%set(key = materials_distribution_case_key, value = checkerboard); check(error==0) 
    error = list%set(key = channels_ratio_key   , value =  0.1 ); check(error==0)
    error = list%set(key = rpb_bddc_threshold_key   , value = 10.0 ); check(error==0)
    error = list%set(key = boundary_mass_trick_key, value =  .false.); check(error==0)
    
    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                 ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                   ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')             ; check(error==0)
    error = switches%set(key = num_dims_key                  , value = '--dim')                      ; check(error==0)
	!error = switches%set(key = hex_mesh_domain_limits_key    , value = '--domain_limits')            ; check(error==0)
    error = switches%set(key = num_cells_x_dir_key           , value = '--number_of_cells')          ; check(error==0)
    error = switches%set(key = num_levels_key                , value = '--number_of_levels')         ; check(error==0)
    error = switches%set(key = num_parts_x_dir_key   , value = '--number_of_parts_per_dir')  ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')   ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order'    )   ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution'        )   ; check(error==0)
    error = switches%set(key = triangulation_generate_key    , value = '--trinagulation-type'    )   ; check(error==0)
    error = switches%set(key = coarse_space_use_vertices_key , value = '--coarse-space-use-vertices'); check(error==0)
    error = switches%set(key = coarse_space_use_edges_key    , value = '--coarse-space-use-edges' )  ; check(error==0)
    error = switches%set(key = coarse_space_use_faces_key    , value = '--coarse-space-use-faces' )  ; check(error==0)
	   error = switches%set(key = bddc_edge_continuity_algorithm_key , value = '--BDDC_edge_continuity_algorithm' ) ; check(error==0)
    error = switches%set(key = bddc_weighting_function_case_key , value = '--BDDC_weighting_function_case' ) ; check(error==0)
    error = switches%set(key = permeability_white_key  , value = '--permeability_white' )  ; check(error==0)
    error = switches%set(key = resistivity_white_key   , value = '--resistivity_white ' )  ; check(error==0)
    error = switches%set(key = permeability_black_key  , value = '--permeability_black' )  ; check(error==0)
    error = switches%set(key = resistivity_black_key   , value = '--resistivity_black' )  ; check(error==0)
    error = switches%set(key = materials_distribution_case_key   , value = '--materials_distribution_case' )  ; check(error==0)
    error = switches%set(key = channels_ratio_key  , value = '--channels_ratio' )  ; check(error==0)
    error = switches%set(key = rpb_bddc_threshold_key  , value = '--rpb_bddc_threshold' )  ; check(error==0)
    error = switches%set(key = boundary_mass_trick_key  , value = '--boundary_mass_trick' )  ; check(error==0)
                                                             
    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = num_dims_key               , value = '-dm')       ; check(error==0)
	   !error = switches_ab%set(key = hex_mesh_domain_limits_key , value = '-dl')       ; check(error==0)
    error = switches_ab%set(key = num_cells_x_dir_key        , value = '-n')        ; check(error==0) 
    error = switches_ab%set(key = num_levels_key             , value = '-l')        ; check(error==0)
    error = switches_ab%set(key = num_parts_x_dir_key, value = '-np')       ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key     , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triangulation_generate_key , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_vertices_key , value = '-use-vertices'); check(error==0)
    error = switches_ab%set(key = coarse_space_use_edges_key    , value = '-use-edges' )  ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_faces_key    , value = '-use-faces' )  ; check(error==0)
	   error = switches_ab%set(key = bddc_edge_continuity_algorithm_key , value = '-edge_cont' )  ; check(error==0)
    error = switches_ab%set(key = bddc_weighting_function_case_key , value = '-bddc_weights' )  ; check(error==0)
    error = switches_ab%set(key = permeability_white_key   , value = '-permeability_white' )  ; check(error==0)
    error = switches_ab%set(key = resistivity_white_key    , value = '-resistivity_white ' )  ; check(error==0)
    error = switches_ab%set(key = permeability_black_key   , value = '-permeability_black' )  ; check(error==0)
    error = switches_ab%set(key = resistivity_black_key    , value = '-resistivity_black' )  ; check(error==0)
    error = switches_ab%set(key = materials_distribution_case_key, value = '-materials_case' )  ; check(error==0)
    error = switches_ab%set(key = channels_ratio_key    , value = '-channels_ratio' )  ; check(error==0)
    error = switches_ab%set(key = rpb_bddc_threshold_key    , value = '-rpb_bddc_threshold' )  ; check(error==0)
    error = switches_ab%set(key = boundary_mass_trick_key    , value = '-bmass_trick' )  ; check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')            ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                    ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                         ; check(error==0)
    error = helpers%set(key = num_dims_key                   , value = 'Number of space dimensions')               ; check(error==0)
	   !error = helpers%set(key = hex_mesh_domain_limits_key     , value = 'Domain limits of the mesh')                ; check(error==0)
    error = helpers%set(key = num_cells_x_dir_key            , value = 'Number of cells per dir')                  ; check(error==0)
    error = helpers%set(key = num_levels_key                 , value = 'Number of levels')                         ; check(error==0)
    error = helpers%set(key = num_parts_x_dir_key            , value = 'Number of parts per dir and per level')    ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')  ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')       ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')             ; check(error==0)
    error = helpers%set(key = coarse_space_use_vertices_key , value  = 'Include vertex coarse DoFs in coarse FE space'); check(error==0)
    error = helpers%set(key = coarse_space_use_edges_key    , value  = 'Include edge coarse DoFs in coarse FE space' )  ; check(error==0)
    error = helpers%set(key = coarse_space_use_faces_key    , value  = 'Include face coarse DoFs in coarse FE space' )  ; check(error==0)
    
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triangulation_generate_key     , value = msg)  ; check(error==0)
	
   	msg = 'Specify BDDC space continuity: tangential_average, tangential_average_and_first_order_moment, all_dofs_in_coarse_edges' 
    error = helpers%set(key = bddc_edge_continuity_algorithm_key  , value = msg)  ; check(error==0)
    msg = 'Define BDDC weighting function from: cardinality (inverse of the cardinality of each dof), resistivity, permeability, stiffness (diagonal entries of the operator).'
    error = helpers%set(key = bddc_weighting_function_case_key, value = msg  ); check(error==0) 
                        
    error = helpers%set(key = permeability_white_key   , value  = 'permeability_white value' ) ; check(error==0)
    error = helpers%set(key = resistivity_white_key    , value  = 'resistivity_white  value' )  ; check(error==0)
    error = helpers%set(key = permeability_black_key   , value  = 'permeability_black value' ) ; check(error==0)
    error = helpers%set(key = resistivity_black_key    , value  = 'resistivity_black value' )  ; check(error==0)
    error = helpers%set(key = materials_distribution_case_key, value  = 'Materials distribution case: choose between: checkerboard, channels' )  ; check(error==0)
    error = helpers%set(key = channels_ratio_key   , value  = 'Ratio channel/non-channel of the cross section for every direction)' ) ; check(error==0)
    error = helpers%set(key = rpb_bddc_threshold_key   , value  = 'Threshold for the relaxed PB-BDDC subparts partition' ) ; check(error==0)
    error = helpers%set(key = boundary_mass_trick_key   , value  = 'Is the boundary mass trick active?' ); check(error==0)
    
    error = required%set(key = dir_path_key                  , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                    , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key              , value = .false.) ; check(error==0)
    error = required%set(key = num_dims_key                  , value = .false.) ; check(error==0)
	   !error = required%set(key = hex_mesh_domain_limits_key    , value = .false.) ; check(error==0)
    error = required%set(key = num_cells_x_dir_key           , value = .false.) ; check(error==0)
    error = required%set(key = num_levels_key                , value = .false.) ; check(error==0)
    error = required%set(key = num_parts_x_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key    , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key        , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key            , value = .false.) ; check(error==0)
    error = required%set(key = triangulation_generate_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_vertices_key , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_edges_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_faces_key    , value = .false.) ; check(error==0)
	   error = required%set(key = bddc_edge_continuity_algorithm_key , value = .false.) ; check(error==0)
    error = required%set(key = bddc_weighting_function_case_key , value = .false.) ; check(error==0)
    error = required%set(key = permeability_white_key   , value = .false.) ; check(error==0)
    error = required%set(key = resistivity_white_key    , value = .false.) ; check(error==0)
    error = required%set(key = permeability_black_key   , value = .false.) ; check(error==0)
    error = required%set(key = resistivity_black_key    , value = .false.) ; check(error==0)
    error = required%set(key = materials_distribution_case_key, value = .false.) ; check(error==0)
    error = required%set(key = channels_ratio_key    , value = .false.) ; check(error==0)
    error = required%set(key = rpb_bddc_threshold_key    , value = .false.) ; check(error==0)
    error = required%set(key = boundary_mass_trick_key    , value = .false.) ; check(error==0)

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
    assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_permeability_white(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_permeability_white
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(permeability_white_key, get_permeability_white))
    error = list%Get(key = permeability_white_key, Value = get_permeability_white)
    assert(error==0)
  end function get_permeability_white
  
      !==================================================================================================
  function get_permeability_black(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_permeability_black
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(permeability_black_key, get_permeability_black))
    error = list%Get(key = permeability_black_key, Value = get_permeability_black)
    assert(error==0)
  end function get_permeability_black
  
     !==================================================================================================
  function get_resistivity_white (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_resistivity_white 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(resistivity_white_key, get_resistivity_white ))
    error = list%Get(key = resistivity_white_key, Value = get_resistivity_white )
    assert(error==0)
  end function get_resistivity_white 
  
       !==================================================================================================
  function get_resistivity_black (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_resistivity_black 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(resistivity_black_key, get_resistivity_black ))
    error = list%Get(key = resistivity_black_key, Value = get_resistivity_black)
    assert(error==0)
  end function get_resistivity_black 
  
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

end module par_pb_bddc_maxwell_params_names
