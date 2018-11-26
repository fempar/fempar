module par_test_h_adaptive_poisson_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private
  
  character(len=*), parameter :: even_cells     = 'even_cells'       
  character(len=*), parameter :: inner_region   = 'inner_region' 
  character(len=*), parameter :: uniform        = 'uniform' 
  
  character(len=*), parameter :: reference_fe_geo_order_key     = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key         = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key             = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key         = 'triangulation_type'
  character(len=*), parameter :: use_void_fes_key               = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key          = 'use_void_fes_case'
  character(len=*), parameter :: coupling_criteria_key          = 'coupling_criteria'
  
  ! Meshing parameters 
  character(len=*), parameter :: refinement_pattern_case_key   = 'refinement_pattern_case'
  character(len=*), parameter :: domain_limits_key             = 'domain_limits'
  character(len=*), parameter :: inner_region_size_key         = 'inner_region_size '
  character(len=*), parameter :: num_refinements_key           = 'num_refinements'
  character(len=*), parameter :: min_num_refinements_key       = 'min_num_refinements'
  

  type, extends(parameter_handler_t) :: par_test_h_adaptive_poisson_params_t
     private
     contains
       procedure :: define_parameters  => par_test_h_adaptive_poisson_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_use_void_fes
       procedure, non_overridable             :: get_use_void_fes_case
       procedure, non_overridable             :: get_refinement_pattern_case 
       procedure, non_overridable             :: get_domain_limits
       procedure, non_overridable             :: get_inner_region_size 
       procedure, non_overridable             :: get_num_refinements 
       procedure, non_overridable             :: get_min_num_refinements
       procedure, non_overridable             :: get_subparts_coupling_criteria 
       !procedure, non_overridable             :: get_num_dims
  end type par_test_h_adaptive_poisson_params_t

  ! Parameters 
  public :: even_cells, inner_region, uniform   
  
  ! Types
  public :: par_test_h_adaptive_poisson_params_t

contains

  !==================================================================================================
  subroutine par_test_h_adaptive_poisson_params_define_parameters(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    error = list%set(key = dir_path_key            , value = '.') ; check(error==0)
    error = list%set(key = prefix_key              , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key        , value = '.') ; check(error==0)
    error = list%set(key = num_dims_key          , value =  2)                   ; check(error==0)
    error = list%set(key = num_cells_x_dir_key       , value =  [12,12,12])          ; check(error==0)
    error = list%set(key = is_dir_periodic_key               , value =  [0,0,0])             ; check(error==0)
    error = list%set(key = num_levels_key              , value =  3)                   ; check(error==0)
    error = list%set(key = num_parts_x_dir_key       , value =  [4,4,0,2,2,0,1,1,0]) ; check(error==0)
    error = list%set(key = reference_fe_geo_order_key        , value =  1)                   ; check(error==0)
    error = list%set(key = reference_fe_order_key            , value =  1)                   ; check(error==0)
    error = list%set(key = write_solution_key                , value =  .false.)             ; check(error==0)
    error = list%set(key = triangulation_generate_key        , value =  triangulation_generate_from_mesh) ; check(error==0)
    error = list%set(key = coarse_space_use_vertices_key     , value =  .true.)                      ; check(error==0)
    error = list%set(key = coarse_space_use_edges_key        , value =  .true.)                      ; check(error==0)
    error = list%set(key = coarse_space_use_faces_key        , value =  .true.)                      ; check(error==0)
    error = list%set(key = use_void_fes_key                  , value =  .false.)                     ; check(error==0)
    error = list%set(key = use_void_fes_case_key             , value =  'popcorn')                   ; check(error==0)
    error = list%set(key = refinement_pattern_case_key       , value = inner_region  )          ; check(error==0)
    error = list%set(key = domain_limits_key                 , value = [0.0,1.0,0.0,1.0,0.0,1.0]) ; check(error==0)
    error = list%set(key = inner_region_size_key , value = [0.1,0.1,0.1]) ; check(error==0)
    error = list%set(key = num_refinements_key               , value = 3) ; check(error==0)
    error = list%set(key = min_num_refinements_key           , value = 1) ; check(error==0)
    error = list%set(key = coupling_criteria_key    , value = loose_coupling) ; check(error==0) 

    ! Only some of them are controlled from cli
    error = switches%set(key = dir_path_key                  , value = '--dir-path')                 ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')                   ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')             ; check(error==0)
    error = switches%set(key = num_dims_key      , value = '--dim')                      ; check(error==0)
    error = switches%set(key = num_cells_x_dir_key   , value = '--num_cells')          ; check(error==0)
    error = switches%set(key = num_levels_key          , value = '--num_levels')         ; check(error==0)
    error = switches%set(key = num_parts_x_dir_key   , value = '--num_parts_x_dir')  ; check(error==0)
    error = switches%set(key = reference_fe_geo_order_key    , value = '--reference-fe-geo-order')   ; check(error==0)
    error = switches%set(key = reference_fe_order_key        , value = '--reference-fe-order'    )   ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write-solution'        )   ; check(error==0)
    error = switches%set(key = triangulation_generate_key    , value = '--trinagulation-type'    )   ; check(error==0)
    error = switches%set(key = coarse_space_use_vertices_key , value = '--coarse-space-use-vertices'); check(error==0)
    error = switches%set(key = coarse_space_use_edges_key    , value = '--coarse-space-use-edges' )  ; check(error==0)
    error = switches%set(key = coarse_space_use_faces_key    , value = '--coarse-space-use-faces' )  ; check(error==0)
    error = switches%set(key = use_void_fes_key              , value = '--use-void-fes' )            ; check(error==0)
    error = switches%set(key = use_void_fes_case_key         , value = '--use-void-fes-case' )       ; check(error==0)
    error = switches%set(key = refinement_pattern_case_key   , value = '--refinement_pattern_case' )       ; check(error==0)
    error = switches%set(key = domain_limits_key      , value = '--domain_limits')     ; check(error==0)
    error = switches%set(key = inner_region_size_key   , value = '--inner_region_size') ; check(error==0)
    error = switches%set(key = num_refinements_key    , value = '--num_refinements') ; check(error==0)
    error = switches%set(key = min_num_refinements_key, value = '--min_num_refinements') ; check(error==0)
    error = switches%set(key = coupling_criteria_key, value = '--subparts_coupling_criteria') ; check(error==0)

    error = switches_ab%set(key = dir_path_key               , value = '-d')        ; check(error==0) 
    error = switches_ab%set(key = prefix_key                 , value = '-p')        ; check(error==0) 
    error = switches_ab%set(key = dir_path_out_key           , value = '-o')        ; check(error==0) 
    error = switches_ab%set(key = num_dims_key   , value = '-dm')      ; check(error==0)
    error = switches_ab%set(key = num_cells_x_dir_key, value = '-n')        ; check(error==0) 
    error = switches_ab%set(key = num_levels_key       , value = '-l')        ; check(error==0)
    error = switches_ab%set(key = num_parts_x_dir_key, value = '-np')       ; check(error==0)
    error = switches_ab%set(key = reference_fe_geo_order_key , value = '-gorder')   ; check(error==0)
    error = switches_ab%set(key = reference_fe_order_key     , value = '-order')    ; check(error==0)
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution'); check(error==0)
    error = switches_ab%set(key = triangulation_generate_key , value = '-tt')       ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_vertices_key , value = '-use-vertices'); check(error==0)
    error = switches_ab%set(key = coarse_space_use_edges_key    , value = '-use-edges' )  ; check(error==0)
    error = switches_ab%set(key = coarse_space_use_faces_key    , value = '-use-faces' )  ; check(error==0)
    error = switches_ab%set(key = use_void_fes_key              , value = '-use-voids' )  ; check(error==0)
    error = switches_ab%set(key = use_void_fes_case_key         , value = '-use-voids-case' ); check(error==0)
    error = switches_ab%set(key = refinement_pattern_case_key   , value = '-refinement-pattern-case' ); check(error==0)
    error = switches_ab%set(key = domain_limits_key          , value = '-dl')       ; check(error==0)
    error = switches_ab%set(key = inner_region_size_key       , value = '-ir_size')    ; check(error==0)
    error = switches_ab%set(key = num_refinements_key        , value = '-num_refs')    ; check(error==0)
    error = switches_ab%set(key = min_num_refinements_key    , value = '-min_num_refs')    ; check(error==0)
    error = switches_ab%set(key = coupling_criteria_key    , value = '-subparts_coupling')    ; check(error==0)

    error = helpers%set(key = dir_path_key                   , value = 'Directory of the source files')            ; check(error==0)
    error = helpers%set(key = prefix_key                     , value = 'Name of the GiD files')                    ; check(error==0)
    error = helpers%set(key = dir_path_out_key               , value = 'Output Directory')                         ; check(error==0)
    error = helpers%set(key = num_dims_key       , value = 'Number of space dimensions')               ; check(error==0)
    error = helpers%set(key = num_cells_x_dir_key    , value = 'Number of cells per dir')                  ; check(error==0)
    error = helpers%set(key = num_levels_key           , value = 'Number of levels')                         ; check(error==0)
    error = helpers%set(key = num_parts_x_dir_key    , value = 'Number of parts per dir and per level')    ; check(error==0)
    error = helpers%set(key = reference_fe_geo_order_key     , value = 'Order of the triangulation reference fe')  ; check(error==0)
    error = helpers%set(key = reference_fe_order_key         , value = 'Order of the fe space reference fe')       ; check(error==0)
    error = helpers%set(key = write_solution_key             , value = 'Write solution in VTK format')             ; check(error==0)
    error = helpers%set(key = coarse_space_use_vertices_key , value  = 'Include vertex coarse DoFs in coarse FE space'); check(error==0)
    error = helpers%set(key = coarse_space_use_edges_key    , value  = 'Include edge coarse DoFs in coarse FE space' )  ; check(error==0)
    error = helpers%set(key = coarse_space_use_faces_key    , value  = 'Include face coarse DoFs in coarse FE space' )  ; check(error==0)
    error = helpers%set(key = use_void_fes_key              , value  = 'Use a hybrid FE space formed by full and void FEs' )  ; check(error==0)
    error = helpers%set(key = use_void_fes_case_key         , value  = 'Select where to put void fes using one of the predefined patterns. Possible values: `popcorn`, `half`, `quarter` ' ); check(error==0)
    error = helpers%set(key = refinement_pattern_case_key   , value  = 'Select refinement pattern. Possible values: even_cells, centered_refinement' ); check(error==0)
    error = helpers%set(key = domain_limits_key     , value = 'Domain limits of the mesh')                ; check(error==0)
    error = helpers%set(key = inner_region_size_key  , value = 'Concentric with the domain refined area length) ') ; check(error==0)
    error = helpers%set(key = num_refinements_key     , value = 'Number of adaptive mesh refinements from a plain cell') ; check(error==0)
    error = helpers%set(key = min_num_refinements_key , value = 'Minimum number of adaptive mesh refinements for any cell') ; check(error==0)
    error = helpers%set(key = coupling_criteria_key, value = 'Criteria to decide whether two subparts are connected or not and identify disconnected parts accordingly') ; check(error==0)
 
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triangulation_generate_key     , value = msg)  ; check(error==0)
    

    error = required%set(key = dir_path_key                  , value = .false.) ; check(error==0)
    error = required%set(key = prefix_key                    , value = .false.) ; check(error==0)
    error = required%set(key = dir_path_out_key              , value = .false.) ; check(error==0)
    error = required%set(key = num_dims_key      , value = .false.) ; check(error==0)
    error = required%set(key = num_cells_x_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = num_levels_key          , value = .false.) ; check(error==0)
    error = required%set(key = num_parts_x_dir_key   , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_geo_order_key    , value = .false.) ; check(error==0)
    error = required%set(key = reference_fe_order_key        , value = .false.) ; check(error==0)
    error = required%set(key = write_solution_key            , value = .false.) ; check(error==0)
    error = required%set(key = triangulation_generate_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_vertices_key , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_edges_key    , value = .false.) ; check(error==0)
    error = required%set(key = coarse_space_use_faces_key    , value = .false.) ; check(error==0)
    error = required%set(key = use_void_fes_key              , value = .false.) ; check(error==0)
    error = required%set(key = use_void_fes_case_key         , value = .false.) ; check(error==0)
    error = required%set(key = refinement_pattern_case_key   , value = .false.) ; check(error==0)
    error = required%set(key = domain_limits_key,                   value = .false.) ; check(error==0)
    error = required%set(key = inner_region_size_key  , value = .false.)  ; check(error==0)
    error = required%set(key = num_refinements_key,                 value = .false.)  ; check(error==0)
    error = required%set(key = min_num_refinements_key,             value = .false.)  ; check(error==0)
    error = required%set(key = coupling_criteria_key,               value = .false.)  ; check(error==0)

  end subroutine par_test_h_adaptive_poisson_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_use_void_fes(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    logical                                       :: get_use_void_fes
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_void_fes_key, get_use_void_fes))
    error = list%Get(key = use_void_fes_key, Value = get_use_void_fes)
    assert(error==0)
  end function get_use_void_fes

  !==================================================================================================
  function get_use_void_fes_case(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_use_void_fes_case
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_void_fes_case_key, 'string'))
    error = list%GetAsString(key = use_void_fes_case_key, string = get_use_void_fes_case)
    assert(error==0)
  end function get_use_void_fes_case
  
    !==================================================================================================
  function get_refinement_pattern_case(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                            :: get_refinement_pattern_case
    type(ParameterList_t), pointer                           :: list
    integer(ip)                                              :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(refinement_pattern_case_key, dummy_string))
    error = list%GetAsString(key = refinement_pattern_case_key, string = get_refinement_pattern_case)
    assert(error==0)
  end function get_refinement_pattern_case
    
  !==================================================================================================
  function get_domain_limits(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    real(rp)                                  :: get_domain_limits(6)
    type(ParameterList_t), pointer            :: list
    integer(ip)                               :: error
    list  => this%get_values()
    assert(list%isAssignable(domain_limits_key, get_domain_limits))
    error = list%Get(key = domain_limits_key, Value = get_domain_limits)
    assert(error==0)
  end function get_domain_limits

  !==================================================================================================
  function get_inner_region_size(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    real(rp)                                  :: get_inner_region_size(0:SPACE_DIM-1)
    type(ParameterList_t), pointer            :: list
    integer(ip)                               :: error
    list  => this%get_values()
    assert(list%isAssignable(inner_region_size_key , get_inner_region_size ))
    error = list%Get(key = inner_region_size_key , Value = get_inner_region_size )
    assert(error==0)
  end function get_inner_region_size

  !==================================================================================================
  function get_num_refinements(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_num_refinements
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(num_refinements_key, get_num_refinements))
    error = list%Get(key = num_refinements_key, Value = get_num_refinements)
    assert(error==0)
  end function get_num_refinements

  !==================================================================================================
  function get_min_num_refinements(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_min_num_refinements
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(min_num_refinements_key, get_min_num_refinements))
    error = list%Get(key = min_num_refinements_key, Value = get_min_num_refinements)
    assert(error==0)
  end function get_min_num_refinements
  
  !==================================================================================================
  function get_subparts_coupling_criteria(this)
    implicit none
    class(par_test_h_adaptive_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                            :: get_subparts_coupling_criteria
    type(ParameterList_t), pointer                           :: list
    integer(ip)                                              :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(coupling_criteria_key, dummy_string))
    error = list%GetAsString(key = coupling_criteria_key, string = get_subparts_coupling_criteria)
    assert(error==0)
  end function get_subparts_coupling_criteria

end module par_test_h_adaptive_poisson_params_names
