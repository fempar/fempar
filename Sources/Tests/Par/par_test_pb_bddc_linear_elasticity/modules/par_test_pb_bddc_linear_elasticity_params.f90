module par_test_pb_bddc_linear_elasticity_params_names
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
  character(len=*), parameter :: nchannel_x_direction_key   = 'nchannel_x_direction' 
  character(len=*), parameter :: nparts_with_channels_key   = 'nparts_with_channels'
  character(len=*), parameter :: is_a_beam_key              = 'is_a_beam'
  character(len=*), parameter :: size_sub_object_key        = 'size_sub_object_key'
  

  type, extends(parameter_handler_t) :: par_test_pb_bddc_linear_elasticity_params_t
     private
     contains
       procedure                              :: define_parameters  => par_test_pb_bddc_linear_elasticity_params_define_parameters
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
       procedure, non_overridable             :: get_hex_mesh_domain_limits
       procedure, non_overridable             :: get_is_a_beam
       procedure, non_overridable             :: get_nparts
       procedure, non_overridable             :: get_num_cells_x_dir
       procedure, non_overridable             :: get_size_sub_object
       !procedure, non_overridable             :: get_num_dims
  end type par_test_pb_bddc_linear_elasticity_params_t

  ! Types
  public :: par_test_pb_bddc_linear_elasticity_params_t, standard_bddc, pb_bddc

contains

  !==================================================================================================
  subroutine par_test_pb_bddc_linear_elasticity_params_define_parameters(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t), intent(inout) :: this
    character(len=:), allocatable                                     :: msg

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
    call this%add(write_matrices_key, '--write-matrices', .false., 'Write local-to-subdomain sparse matrices  in matrix market format', switch_ab='-wmatrices') 
    call this%add(jump_key, '--jump', 1, 'Jump of physical parameter in the inclusion', switch_ab='-j')
    call this%add(inclusion_key, '--inclusion', 1, 'Inclusion type', switch_ab='-i')
    call this%add(coarse_fe_handler_type_key, '--coarse-fe-handler', pb_bddc, 'Which coarse fe handler to use?', switch_ab='-coarse-handler')
    call this%add(nchannel_x_direction_key, '--nchannel_x_direction', [1,1,1], 'Number of channels per direction', switch_ab='-nc')
    call this%add(nparts_with_channels_key, '--nparts_with_channels', [1,1,1], 'Number of parts per with channels', switch_ab='-npwc')
    call this%add(size_sub_object_key, '--size_sub_object', 1, 'Size of subobject in number of mesh size (must be an integer)', switch_ab='-sso')

    call this%add(discrete_integration_type_key, '--discrete_integration', heterogeneous, &
                  'Which type of discrete integration to use? homogeneous or heterogeneous', &
                  switch_ab='-integration')
    call this%add(struct_hex_triang_domain_limits_key, '--hex_mesh_domain_limits', [0.0_rp,2.0_rp,0.0_rp,0.5_rp,0.0_rp,0.5_rp], &
                  'Limits of the domain', &
                  switch_ab='-domain_limits')
    call this%add(is_a_beam_key, '--is_a_beam', .true., 'Is the studied object a beam', switch_ab='-is_a_beam')
    call this%add(size_sub_object_key, '--size_sub_object', 1, 'Size of subobject in number of mesh size (must be an integer)', switch_ab='-sso')

  end subroutine par_test_pb_bddc_linear_elasticity_params_define_parameters

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(reference_fe_order_key, get_reference_fe_order))
    error = list%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
    assert(error==0)
  end function get_reference_fe_order
  
  !==========================================================================================par_test_pb_bddc_linear_elasticity_params_t========
  function get_write_solution(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    check(error==0)
  end function get_write_solution
  
  !==========================================================================================par_test_pb_bddc_linear_elasticity_params_t========
  function get_write_matrices(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(triang_generate_key, get_triangulation_type))
    error = list%Get(key = triang_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  function get_jump(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
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
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_nparts_with_channels(3)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(nparts_with_channels_key, get_nparts_with_channels))
    error = list%Get(key = nparts_with_channels_key, Value = get_nparts_with_channels)
    assert(error==0)
  end function get_nparts_with_channels
!==================================================================================================
  function get_hex_mesh_domain_limits(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    real(rp)                                   :: get_hex_mesh_domain_limits(6)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_domain_limits_key, get_hex_mesh_domain_limits))
    error = list%Get(key = struct_hex_triang_domain_limits_key, Value = get_hex_mesh_domain_limits)
    assert(error==0)
  end function get_hex_mesh_domain_limits

  !==================================================================================================
  function get_is_a_beam(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    logical                                       :: get_is_a_beam
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(is_a_beam_key, get_is_a_beam))
    error = list%Get(key = is_a_beam_key, Value = get_is_a_beam)
    check(error==0)
  end function get_is_a_beam

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: num_levels
    integer(ip)                                   :: get_nparts(3)
    integer(ip), allocatable :: num_parts_x_dir(:) ! 0:SPACE_DIM-1)
    integer(ip), allocatable :: array_size(:)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_num_levels_key, num_levels))
    error = list%Get(key = struct_hex_triang_num_levels_key, Value = num_levels)
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
  function get_num_cells_x_dir(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_num_cells_x_dir(3)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(struct_hex_triang_num_cells_dir, get_num_cells_x_dir))
    error = list%Get(key = struct_hex_triang_num_cells_dir, Value = get_num_cells_x_dir)
    assert(error==0)
  end function get_num_cells_x_dir
  !==================================================================================================
  function get_size_sub_object(this)
    implicit none
    class(par_test_pb_bddc_linear_elasticity_params_t) , intent(in) :: this
    integer(ip)                                   :: get_size_sub_object
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(size_sub_object_key, get_size_sub_object))
    error = list%Get(key = size_sub_object_key, Value = get_size_sub_object)
    assert(error==0)
  end function get_size_sub_object
  
  !==================================================================================================

 
end module par_test_pb_bddc_linear_elasticity_params_names
