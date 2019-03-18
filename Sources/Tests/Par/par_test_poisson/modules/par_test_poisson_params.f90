module par_test_poisson_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_geo_order_key = 'reference_fe_geo_order'
  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key     = 'triangulation_type'
  character(len=*), parameter :: use_void_fes_key           = 'use_void_fes'
  character(len=*), parameter :: use_void_fes_case_key      = 'use_void_fes_case'
  character(len=*), parameter :: preconditioner_type_key    = 'preconditioner_type'

  type, extends(parameter_handler_t) :: par_test_poisson_params_t
     private
     contains
       procedure :: define_parameters  => par_test_poisson_params_define_parameters
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_reference_fe_geo_order
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_use_void_fes
       procedure, non_overridable             :: get_use_void_fes_case
       procedure, non_overridable             :: get_preconditioner_type
       !procedure, non_overridable             :: get_num_dims
  end type par_test_poisson_params_t

  ! Types
  public :: par_test_poisson_params_t
  
contains

  !==================================================================================================
  subroutine par_test_poisson_params_define_parameters(this)
    implicit none
    class(par_test_poisson_params_t), intent(inout) :: this
    character(len=:), allocatable                   :: msg

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
    call this%add(triang_generate_key, '--triangulation-type', triangulation_generate_from_mesh, msg, switch_ab='-tt')
    call this%add(coarse_space_use_vertices_key, '--coarse-space-use-vertices', .true., 'Include vertex coarse DoFs in coarse FE space', switch_ab='-use-vertices')
    call this%add(coarse_space_use_edges_key, '--coarse-space-use-edges', .true., 'Include edge coarse DoFs in coarse FE space', switch_ab='-use-edges')
    call this%add(coarse_space_use_faces_key, '--coarse-space-use-faces', .true., 'Include face coarse DoFs in coarse FE space', switch_ab='-use-faces')

    ! Specific
    call this%add(use_void_fes_key, '--use-void-fes', .false., 'Use a hybrid FE space formed by full and void FEs', '-use-voids')
    call this%add(use_void_fes_case_key, '--use-void-fes-case', 'popcorn', 'Select where to put void fes using one of the predefined patterns. Possible values: `popcorn`, `half`, `quarter` ', switch_ab='-use-voids-case')
    call this%add(preconditioner_type_key, '--preconditioner-type', 'mlbddc', 'Select preconditioner type. Possible values: `identity`, `jacobi`, `mlbddc`.', switch_ab='-prec-type')

  end subroutine par_test_poisson_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
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
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_use_void_fes_case
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(use_void_fes_case_key, 'string'))
    error = list%GetAsString(key = use_void_fes_case_key, string = get_use_void_fes_case)
    assert(error==0)
  end function get_use_void_fes_case
  
  !==================================================================================================
  function get_preconditioner_type(this)
    implicit none
    class(par_test_poisson_params_t) , intent(in) :: this
    character(len=:), allocatable                 :: get_preconditioner_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(preconditioner_type_key, 'string'))
    error = list%GetAsString(key = preconditioner_type_key, string = get_preconditioner_type)
    assert(error==0)
  end function get_preconditioner_type

end module par_test_poisson_params_names
