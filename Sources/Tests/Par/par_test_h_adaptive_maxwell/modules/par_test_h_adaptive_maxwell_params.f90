module par_test_h_adaptive_maxwell_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: even_cells     = 'even_cells'       
  character(len=*), parameter :: inner_region   = 'inner_region'   
  character(len=*), parameter :: inner_sphere   = 'inner_sphere'
  character(len=*), parameter :: uniform        = 'uniform'
  character(len=*), parameter :: error_based    = 'error_based'
  
  character(len=*), parameter :: reference_fe_geo_order_key = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'
  character(len=*), parameter :: use_void_fes_key           = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key      = 'use_void_fes_case'
  character(len=*), parameter :: coupling_criteria_key      = 'coupling_criteria'
    ! Meshing parameters 
  character(len=*), parameter :: refinement_pattern_case_key   = 'refinement_pattern_case'
  character(len=*), parameter :: domain_limits_key             = 'domain_limits'
  character(len=*), parameter :: inner_region_size_key         = 'inner_region_size'
  character(len=*), parameter :: refinement_radius_key         = 'refinement_radius'
  character(len=*), parameter :: num_refinements_key           = 'num_refinements'
  character(len=*), parameter :: min_num_refinements_key       = 'min_num_refinements'
  ! Solution 
  character(len=*), parameter :: analytical_function_case_key  = 'analytical_function_case'

  type, extends(parameter_handler_t) :: par_test_h_adaptive_maxwell_params_t
     private
     contains
       procedure :: define_parameters  => par_test_h_adaptive_maxwell_params_define_parameters
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
       procedure, non_overridable             :: get_refinement_radius
       procedure, non_overridable             :: get_num_refinements 
       procedure, non_overridable             :: get_min_num_refinements
       procedure, non_overridable             :: get_subparts_coupling_criteria
       procedure, non_overridable             :: get_analytical_function_case
  end type par_test_h_adaptive_maxwell_params_t

  ! Parameters 
  public :: even_cells, inner_region, inner_sphere, uniform, error_based
  
  ! Types
  public :: par_test_h_adaptive_maxwell_params_t

contains

  !==================================================================================================
  subroutine par_test_h_adaptive_maxwell_params_define_parameters(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t), intent(inout) :: this
    character(len=:), allocatable                              :: msg

    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh

    ! Common
    call this%add(dir_path_key,'--dir-path', '.', 'Directory of the source files', switch_ab='-d')
    call this%add(prefix_key, '--prefix', 'square', 'Name of the GiD files', switch_ab='-p')
    call this%add(dir_path_out_key, '--dir-path-out', '.', 'Output Directory', switch_ab='-o')
    call this%add(struct_hex_triang_num_dims_key, '--dim', 2, 'Number of space dimensions', switch_ab='-dm')
    call this%add(struct_hex_triang_num_cells_dir, '--num_cells', [12,12,12], 'Number of cells per dir', switch_ab='-n')
    call this%add(struct_hex_triang_is_dir_periodic_key, '--STRUCT_HEX_TRIANG_IS_DIR_PERIODIC', [0,0,0], 'Is the mesh periodic for every dimension')           
    call this%add(struct_hex_triang_num_levels_key, '--num_levels', 3, 'Number of levels', switch_ab='-l')
    call this%add(struct_hex_triang_num_parts_x_dir_key, '--num_parts_x_dir', [4,4,0,2,2,0,1,1,0], 'Number of parts per dir and per level', switch_ab='-np')
    call this%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1, 'Order of the triangulation reference fe', switch_ab='-gorder')
    call this%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call this%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')
    call this%add(triang_generate_key, '--trinagulation-type', triangulation_generate_from_mesh, msg, switch_ab='-tt')
    call this%add(coarse_space_use_vertices_key, '--coarse-space-use-vertices', .true., 'Include vertex coarse DoFs in coarse FE space', switch_ab='-use-vertices')
    call this%add(coarse_space_use_edges_key, '--coarse-space-use-edges', .true., 'Include edge coarse DoFs in coarse FE space', switch_ab='-use-edges')
    call this%add(coarse_space_use_faces_key, '--coarse-space-use-faces', .true., 'Include face coarse DoFs in coarse FE space', switch_ab='-use-faces')

    ! Specific
    call this%add(use_void_fes_key, '--use-void-fes', .true., 'Use a hybrid FE space formed by full and void FEs', switch_ab='-use-voids')
    call this%add(use_void_fes_case_key, '--use-void-fes-case', 'popcorn', &
                 'Select where to put void fes using one of the predefined patterns. Possible values: `popcorn`, `half`, `quarter`', &
                 switch_ab='-use-voids-case')
    call this%add(refinement_pattern_case_key, '--refinement_pattern_case', inner_region, &
                'Select refinement pattern. Possible values: even_cells, inner_region, inner_sphere, uniform, error_based', &
                switch_ab='-refinement-pattern-case' )
    call this%add(domain_limits_key, '--domain_limits', [-1.0,1.0,-1.0,1.0,-1.0,1.0], 'Domain limits of the mesh', switch_ab='-dl')
    call this%add(inner_region_size_key, '--inner_region_size', [0.1,0.1,0.1], 'Concentric with the domain refined area length)', switch_ab='-ir_size')
    call this%add(num_refinements_key, '--num_refinements', 2, 'Number of adaptive mesh refinements from a plain cell', switch_ab='-num_refs')
    call this%add(min_num_refinements_key, '--min_num_refinements', 2, 'Minimum number of adaptive mesh refinements for any cell', switch_ab='-min_num_refs')
    call this%add(coupling_criteria_key, '--subparts_coupling_criteria', loose_coupling, &
                  'Criteria to decide whether two subparts are connected or not and identify disconnected parts accordingly', &
                  switch_ab='-subparts_coupling')    



    call this%add(refinement_radius_key, '--refinement_radius', 0.1, 'Concentric with the domain refined area Radius)', switch_ab='-R_ref')
    call this%add(analytical_function_case_key, '--analytical_function_case', 'in_fe_space', &
                  'Select analytical solution case. Possible values: in_fe_space, fichera_2D, fichera_3D', &
                  switch_ab='-function-case')

  end subroutine par_test_h_adaptive_maxwell_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_use_void_fes(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
    real(rp)                                  :: get_inner_region_size(0:SPACE_DIM-1)
    type(ParameterList_t), pointer            :: list
    integer(ip)                               :: error
    list  => this%get_values()
    assert(list%isAssignable(inner_region_size_key , get_inner_region_size ))
    error = list%Get(key = inner_region_size_key , Value = get_inner_region_size )
    assert(error==0)
  end function get_inner_region_size
  
    !==================================================================================================
  function get_refinement_radius(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
    real(rp)                                  :: get_refinement_radius
    type(ParameterList_t), pointer            :: list
    integer(ip)                               :: error
    list  => this%get_values()
    assert(list%isAssignable(refinement_radius_key , get_refinement_radius ))
    error = list%Get(key = refinement_radius_key , Value = get_refinement_radius )
    assert(error==0)
  end function get_refinement_radius

  !==================================================================================================
  function get_num_refinements(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
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
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                            :: get_subparts_coupling_criteria
    type(ParameterList_t), pointer                           :: list
    integer(ip)                                              :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(coupling_criteria_key, dummy_string))
    error = list%GetAsString(key = coupling_criteria_key, string = get_subparts_coupling_criteria)
    assert(error==0)
  end function get_subparts_coupling_criteria
  
        !==================================================================================================
  function get_analytical_function_case(this)
    implicit none
    class(par_test_h_adaptive_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                            :: get_analytical_function_case
    type(ParameterList_t), pointer                           :: list
    integer(ip)                                              :: error
    character(1) :: dummy_string
    list  => this%get_values()
    assert(list%isAssignable(analytical_function_case_key, dummy_string))
    error = list%GetAsString(key = analytical_function_case_key, string = get_analytical_function_case)
    assert(error==0)
  end function get_analytical_function_case

end module par_test_h_adaptive_maxwell_params_names
