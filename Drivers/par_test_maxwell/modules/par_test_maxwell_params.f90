module par_test_maxwell_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key           = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key               = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key                   = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key               = 'triangulation_type'    
  character(len=*), parameter :: bddc_edge_continuity_algorithm_key   = 'bddc_edge_continuity_algorithm'
  character(len=*), parameter :: permeability_key                     = 'permeability'
  character(len=*), parameter :: resistivity_key                      = 'resistivity'
  
  type, extends(parameter_handler_t) :: par_test_maxwell_params_t
     private
     contains
       procedure :: define_parameters  => par_test_maxwell_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_permeability 
       procedure, non_overridable             :: get_resistivity 
       !procedure, non_overridable             :: get_num_dimensions
  end type par_test_maxwell_params_t

  ! Types
  public :: par_test_maxwell_params_t

contains

  !==================================================================================================
  subroutine par_test_maxwell_params_define_parameters(this)
    implicit none
    class(par_test_maxwell_params_t), intent(inout) :: this
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
    error = list%set(key = execution_context_key             , value =  mpi_context)                      ; check(error==0)
    error = list%set(key = coarse_space_use_vertices_key     , value =  .true.)                                    ; check(error==0)
    error = list%set(key = coarse_space_use_edges_key        , value =  .true.)                                    ; check(error==0)
    error = list%set(key = coarse_space_use_faces_key        , value =  .false.)                                   ; check(error==0)
	   error = list%set(key = bddc_edge_continuity_algorithm_key, value =  tangential_average_and_first_order_moment) ; check(error==0)
    error = list%set(key = permeability_key   , value =  1.0 ); check(error==0)
    error = list%set(key = resistivity_key    , value =  1.0 ); check(error==0)
 
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
    error = switches%set(key = execution_context_key         , value = '--execution_context'    )    ; check(error==0)
    error = switches%set(key = coarse_space_use_vertices_key , value = '--coarse-space-use-vertices'); check(error==0)
    error = switches%set(key = coarse_space_use_edges_key    , value = '--coarse-space-use-edges' )  ; check(error==0)
    error = switches%set(key = coarse_space_use_faces_key    , value = '--coarse-space-use-faces' )  ; check(error==0)
	   error = switches%set(key = bddc_edge_continuity_algorithm_key , value = '--BDDC_edge_continuity_algorithm' ) ; check(error==0)
    error = switches%set(key = permeability_key  , value = '--permeability' )  ; check(error==0)
    error = switches%set(key = resistivity_key   , value = '--resistivity' )  ; check(error==0)
                                                             
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
    error = switches_ab%set(key = execution_context_key      , value = '-exe')       ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_vertices_key , value = '-use-vertices'); check(error==0)
    error = switches_ab%set(key = coarse_space_use_edges_key    , value = '-use-edges' )  ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_faces_key    , value = '-use-faces' )  ; check(error==0)
	   error = switches_ab%set(key = bddc_edge_continuity_algorithm_key , value = '-edge_cont' )  ; check(error==0)
    error = switches_ab%set(key = permeability_key   , value = '-permeability' )  ; check(error==0)
    error = switches_ab%set(key = resistivity_key    , value = '-resistivity' )  ; check(error==0)

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
    
    msg = 'serial (*) or mpi (*) context?'
    write(msg(9:9),'(i1)') serial_context
    write(msg(20:20),'(i1)') mpi_context
    error = helpers%set(key = execution_context_key     , value = msg)  ; check(error==0)
	
	msg = 'Specify BDDC space continuity: Tangent component on coarse edges (*), tangent component + first order moment (*) or one-to-one over all fine edges (*) '
    write(msg(67:67),'(i1)') tangential_average 
    write(msg(111:111),'(i1)') tangential_average_and_first_order_moment 
	write(msg(149:149), '(i1)') all_dofs_in_coarse_edges  
    error = helpers%set(key = bddc_edge_continuity_algorithm_key  , value = msg)  ; check(error==0)
    error = helpers%set(key = permeability_key   , value  = 'Permeability value' ) ; check(error==0)
    error = helpers%set(key = resistivity_key    , value  = 'Resistivity value' )  ; check(error==0)

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
    error = required%set(key = execution_context_key         , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_vertices_key , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_edges_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_faces_key    , value = .false.) ; check(error==0)
	   error = required%set(key = bddc_edge_continuity_algorithm_key , value = .false.) ; check(error==0)
    error = required%set(key = permeability_key   , value = .false.) ; check(error==0)
    error = required%set(key = resistivity_key    , value = .false.) ; check(error==0)

  end subroutine par_test_maxwell_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_maxwell_params_t) , intent(in) :: this
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
    class(par_test_maxwell_params_t) , intent(in) :: this
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
    class(par_test_maxwell_params_t) , intent(in) :: this
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
    class(par_test_maxwell_params_t) , intent(in) :: this
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
    class(par_test_maxwell_params_t) , intent(in) :: this
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
    class(par_test_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_permeability(this)
    implicit none
    class(par_test_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_permeability
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(permeability_key, get_permeability))
    error = list%Get(key = permeability_key, Value = get_permeability)
    assert(error==0)
  end function get_permeability
  
     !==================================================================================================
  function get_resistivity(this)
    implicit none
    class(par_test_maxwell_params_t) , intent(in) :: this
    real(rp)                                      :: get_resistivity
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(resistivity_key, get_resistivity))
    error = list%Get(key = resistivity_key, Value = get_resistivity)
    assert(error==0)
  end function get_resistivity


end module par_test_maxwell_params_names
