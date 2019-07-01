module par_pb_bddc_maxwell_params_names
  use fempar_names

  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key               = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key                   = 'write_solution'        
  character(len=*), parameter :: triangulation_type_key               = 'triangulation_type'    
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
  
  type :: par_pb_bddc_maxwell_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
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
       procedure, non_overridable             :: get_values
       procedure, non_overridable             :: free
  end type par_pb_bddc_maxwell_params_t

  ! Types
  public :: par_pb_bddc_maxwell_params_t
  
  public :: checkerboard, channels, homogeneous, heterogeneous
  public :: unit_constant, constant, sinusoidal

contains

  !==================================================================================================
  subroutine par_test_maxwell_params_define_user_parameters()
    implicit none
    ! Common
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! Specific
    call parameter_handler%add(bddc_edge_continuity_algorithm_key, bddc_edge_continuity_algorithm_cla_name, tangential_average_and_first_order_moment, &
                  'Specify BDDC space continuity: tangential_average, tangential_average_and_first_order_moment, all_dofs_in_coarse_edges', &
                  switch_ab='-edge_cont')
    call parameter_handler%add(mass_coeff_white_key, '--mass_coeff_white', 1.0_rp, 'mass_coeff_white value', switch_ab='-mass_coeff_white')
    call parameter_handler%add(curl_curl_coeff_white_key, '--curl_curl_coeff_white', 1.0_rp, 'curl_curl_coeff_white  value', switch_ab='-curl_curl_coeff_white')
    call parameter_handler%add(mass_coeff_black_key, '--mass_coeff_black', 1.0_rp, 'mass_coeff_black value', switch_ab='-mass_coeff_black')
    call parameter_handler%add(curl_curl_coeff_black_key, '--curl_curl_coeff_black', 1.0_rp, 'curl_curl_coeff_black value', switch_ab='-curl_curl_coeff_black')
    call parameter_handler%add(materials_distribution_case_key, '--materials_distribution_case', homogeneous, &
                  'Materials distribution case: choose between: checkerboard, channels, radial, heterogeneous', &
                  switch_ab='-materials_case')
    call parameter_handler%add(materials_coefficient_case_key, '--materials_coefficient_case', unit_constant, &
                  'Materials coefficient case: choose between: constant, sinusoidal', &
                  switch_ab='-coefficient_case')
    call parameter_handler%add(channels_ratio_key, '--channels_ratio', 0.1_rp, &
                  'Ratio channel/non-channel of the cross section for every direction)', &
                  switch_ab='-channels_ratio')
    call parameter_handler%add(num_peaks_curl_curl_coeff_key, '--num_peaks_curl_curl_coeff', 3, &
                  'Number of peaks for the sinusoidal function describing the curl_curl_coeff', &
                  switch_ab='-num_peaks_curl_curl_coeff')
    call parameter_handler%add(num_peaks_mass_coeff_key, '--num_peaks_mass_coeff', 3, &
                  'Number of peaks for the sinusoidal function describing the mass_coeff', &
                  switch_ab='-num_peaks_mass_coeff')
    call parameter_handler%add(rpb_bddc_threshold_key, '--rpb_bddc_threshold', 10.0_rp, &
                  'Threshold for the relaxed PB-BDDC subparts partition', &
                  switch_ab='-rpb_bddc_threshold')
    call parameter_handler%add(boundary_mass_trick_key, '--boundary_mass_trick', .false., 'Is the boundary mass trick active?', switch_ab='-bmass_trick')

  end subroutine par_test_maxwell_params_define_user_parameters


  !==================================================================================================

  subroutine process_parameters(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_maxwell_params_define_user_parameters)
  end subroutine process_parameters

  !==================================================================================================

  function get_parameter_list(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                   :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list


  subroutine free(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t), intent(inout) :: this
    call parameter_handler%free()
  end subroutine


  function get_values(this) result(values)
    implicit none
    class(par_pb_bddc_maxwell_params_t), target, intent(in) :: this
    type(ParameterList_t),      pointer            :: values
    values => parameter_handler%get_values()
  end function get_values



  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
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
    list  => parameter_handler%get_values()
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
    call parameter_handler%Get(key = struct_hex_mesh_generator_num_levels_key , Value = num_levels) 
    call parameter_handler%GetAsArray(key = struct_hex_mesh_generator_num_parts_x_dim_key, Value = num_parts_x_dir)
    get_nparts=num_parts_x_dir(1:3)
  end function get_nparts
   
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                   :: get_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => parameter_handler%get_values()
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
    list  => parameter_handler%get_values()
    assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key = write_solution_key, Value = get_write_solution)
    assert(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable :: get_triangulation_type
    call parameter_handler%GetAsString(key = static_triang_generate_from_key, string = get_triangulation_type)
  end function get_triangulation_type 
  
    !==================================================================================================
  function get_mass_coeff_white(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_mass_coeff_white
    call parameter_handler%Get(key = mass_coeff_white_key, Value = get_mass_coeff_white)
  end function get_mass_coeff_white
  
      !==================================================================================================
  function get_mass_coeff_black(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_mass_coeff_black
    call parameter_handler%Get(key = mass_coeff_black_key, Value = get_mass_coeff_black)
  end function get_mass_coeff_black
  
     !==================================================================================================
  function get_curl_curl_coeff_white (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_curl_curl_coeff_white 
    call parameter_handler%Get(key = curl_curl_coeff_white_key, Value = get_curl_curl_coeff_white )
  end function get_curl_curl_coeff_white 
  
       !==================================================================================================
  function get_curl_curl_coeff_black (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_curl_curl_coeff_black 
    call parameter_handler%Get(key = curl_curl_coeff_black_key, Value = get_curl_curl_coeff_black)
  end function get_curl_curl_coeff_black 
  
      !==================================================================================================
  function get_materials_distribution_case(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                    :: get_materials_distribution_case
    call parameter_handler%GetAsString(key = materials_distribution_case_key, string = get_materials_distribution_case)
  end function get_materials_distribution_case
  
        !==================================================================================================
  function get_materials_coefficient_case(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    character(len=:), allocatable                    :: get_materials_coefficient_case
    call parameter_handler%GetAsString(key = materials_coefficient_case_key, string = get_materials_coefficient_case)
  end function get_materials_coefficient_case

       !==================================================================================================
  function get_channels_ratio (this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_channels_ratio 
    call parameter_handler%Get(key = channels_ratio_key, Value = get_channels_ratio )
  end function get_channels_ratio
  
      !==================================================================================================
  function get_rpb_bddc_threshold(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    real(rp)                                         :: get_rpb_bddc_threshold
    call parameter_handler%Get(key = rpb_bddc_threshold_key, Value = get_rpb_bddc_threshold)
  end function get_rpb_bddc_threshold
  
    !==================================================================================================
  function get_boundary_mass_trick(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t), intent(in) :: this
    logical                                         :: get_boundary_mass_trick
    call parameter_handler%Get(key = boundary_mass_trick_key, Value = get_boundary_mass_trick)
  end function get_boundary_mass_trick
  
    !==================================================================================================
  function get_num_peaks_curl_curl_coeff(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                      :: get_num_peaks_curl_curl_coeff
    call parameter_handler%Get(key = num_peaks_curl_curl_coeff_key, Value = get_num_peaks_curl_curl_coeff)
  end function get_num_peaks_curl_curl_coeff
  
      !==================================================================================================
  function get_num_peaks_mass_coeff(this)
    implicit none
    class(par_pb_bddc_maxwell_params_t) , intent(in) :: this
    integer(ip)                                      :: get_num_peaks_mass_coeff
    call parameter_handler%Get(key = num_peaks_mass_coeff_key, Value = get_num_peaks_mass_coeff)
  end function get_num_peaks_mass_coeff
  
end module par_pb_bddc_maxwell_params_names
