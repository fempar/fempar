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
  
  character(len=*), parameter :: unit_constant     = 'unit_constant'
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
  
  public :: checkerboard, channels, homogeneous, heterogeneous
  public :: unit_constant, constant, sinusoidal

contains

  !==================================================================================================
  subroutine par_test_maxwell_params_define_parameters(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t), intent(inout) :: this
    character(len=:), allocatable                      :: msg

    
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh

    ! Common
    call this%add(dir_path_key,'--dir-path', '.', 'Directory of the source files', switch_ab='-d')
    call this%add(prefix_key, '--prefix', 'square', 'Name of the GiD files', switch_ab='-p')
    call this%add(dir_path_out_key, '--dir-path-out', '.', 'Output Directory', switch_ab='-o')
    call this%add(struct_hex_triang_num_dims_key, '--dim', 3, 'Number of space dimensions', switch_ab='-dm')
    call this%add(struct_hex_triang_num_cells_dir, '--num_cells', [8,8,8], 'Number of cells per dir', switch_ab='-n')
    call this%add(struct_hex_triang_is_dir_periodic_key, '--STRUCT_HEX_TRIANG_IS_DIR_PERIODIC', [0,0,0], 'Is the mesh periodic for every dimension')           
    call this%add(struct_hex_triang_num_levels_key, '--num_levels', 2, 'Number of levels', switch_ab='-l')
    call this%add(struct_hex_triang_num_parts_x_dir_key, '--num_parts_x_dir', [2,2,2,1,1,1,0,0,0], 'Number of parts per dir and per level', switch_ab='-np')
    call this%add(reference_fe_geo_order_key, '--reference-fe-geo-order', 1, 'Order of the triangulation reference fe', switch_ab='-gorder')
    call this%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call this%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')
    call this%add(triang_generate_key, '--triangulation-type', triangulation_generate_structured, msg, switch_ab='-tt')
    call this%add(coarse_space_use_vertices_key, '--coarse-space-use-vertices', .true., 'Include vertex coarse DoFs in coarse FE space', switch_ab='-use-vertices')
    call this%add(coarse_space_use_edges_key, '--coarse-space-use-edges', .true., 'Include edge coarse DoFs in coarse FE space', switch_ab='-use-edges')
    call this%add(coarse_space_use_faces_key, '--coarse-space-use-faces', .false., 'Include face coarse DoFs in coarse FE space', switch_ab='-use-faces')

    ! Specific
    call this%add(bddc_edge_continuity_algorithm_key, '--BDDC_edge_continuity_algorithm', tangential_average_and_first_order_moment, &
                  'Specify BDDC space continuity: tangential_average, tangential_average_and_first_order_moment, all_dofs_in_coarse_edges', &
                  switch_ab='-edge_cont')
    call this%add(bddc_scaling_function_case_key, '--BDDC_scaling_function_case', cardinality, &
                  'Specify BDDC space continuity: tangential_average, tangential_average_and_first_order_moment, all_dofs_in_coarse_edges', &
                  switch_ab='-bddc_weights')
    call this%add(mass_coeff_white_key, '--mass_coeff_white', 1.0_rp, 'mass_coeff_white value', switch_ab='-mass_coeff_white')
    call this%add(curl_curl_coeff_white_key, '--curl_curl_coeff_white', 1.0_rp, 'curl_curl_coeff_white  value', switch_ab='-curl_curl_coeff_white')
    call this%add(mass_coeff_black_key, '--mass_coeff_black', 1.0_rp, 'mass_coeff_black value', switch_ab='-mass_coeff_black')
    call this%add(curl_curl_coeff_black_key, '--curl_curl_coeff_black', 1.0_rp, 'curl_curl_coeff_black value', switch_ab='-curl_curl_coeff_black')
    call this%add(materials_distribution_case_key, '--materials_distribution_case', homogeneous, &
                  'Materials distribution case: choose between: checkerboard, channels, radial, heterogeneous', &
                  switch_ab='-materials_case')
    call this%add(materials_coefficient_case_key, '--materials_coefficient_case', unit_constant, &
                  'Materials coefficient case: choose between: constant, sinusoidal', &
                  switch_ab='-coefficient_case')
    call this%add(channels_ratio_key, '--channels_ratio', 0.1_rp, &
                  'Ratio channel/non-channel of the cross section for every direction)', &
                  switch_ab='-channels_ratio')
    call this%add(num_peaks_curl_curl_coeff_key, '--num_peaks_curl_curl_coeff', 3, &
                  'Number of peaks for the sinusoidal function describing the curl_curl_coeff', &
                  switch_ab='-num_peaks_curl_curl_coeff')
    call this%add(num_peaks_mass_coeff_key, '--num_peaks_mass_coeff', 3, &
                  'Number of peaks for the sinusoidal function describing the mass_coeff', &
                  switch_ab='-num_peaks_mass_coeff')
    call this%add(rpb_bddc_threshold_key, '--rpb_bddc_threshold', 10.0_rp, &
                  'Threshold for the relaxed PB-BDDC subparts partition', &
                  switch_ab='-rpb_bddc_threshold')
    call this%add(boundary_mass_trick_key, '--boundary_mass_trick', .false., 'Is the boundary mass trick active?', switch_ab='-bmass_trick')

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
