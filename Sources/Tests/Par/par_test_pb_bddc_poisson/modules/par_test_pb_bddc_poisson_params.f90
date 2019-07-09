module pb_bddc_poisson_params_names
  use fempar_names
  implicit none
#include "debug.i90" 
  private

  character(len=*), parameter :: reference_fe_order_key     = 'reference_fe_order'    
  character(len=*), parameter :: write_solution_key         = 'write_solution'  
  character(len=*), parameter :: write_matrices_key         = 'write_matrices'  
  character(len=*), parameter :: jump_key                   = 'jump'    
  character(len=*), parameter :: inclusion_key              = 'inclusion'  
  character(len=*), parameter :: coarse_fe_handler_type_key = 'coarse_fe_handler_type_key' 
  character(len=*), parameter :: standard_bddc              = 'standard_bddc' 
  character(len=*), parameter :: pb_bddc                    = 'pb_bddc' 
  character(len=*), parameter :: nchannel_x_direction_key   = 'nchannel_x_direction' 
  character(len=*), parameter :: nparts_with_channels_key   = 'nparts_with_channels' 
  character(len=*), parameter :: size_sub_object_key        = 'size_sub_object_key' 
  character(len=*), parameter :: dir_path_out_key           = 'dir_path_out' 

  type :: par_test_pb_bddc_poisson_params_t
     private
     contains
       procedure, non_overridable             :: process_parameters
       procedure, non_overridable             :: get_parameter_list
       procedure, non_overridable             :: get_dir_path_out
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_output_handler_dir_path
       procedure, non_overridable             :: get_output_handler_prefix
       procedure, non_overridable             :: get_reference_fe_order
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_write_matrices
       procedure, non_overridable             :: get_jump
       procedure, non_overridable             :: get_inclusion
       procedure, non_overridable             :: get_coarse_fe_handler_type
       procedure, non_overridable             :: get_nchannel_x_direction
       procedure, non_overridable             :: get_nparts_with_channels
       procedure, non_overridable             :: get_nparts
       procedure, non_overridable             :: get_num_cells_x_dir
       procedure, non_overridable             :: get_size_sub_object
       procedure, non_overridable             :: get_hex_mesh_domain_limits
  end type par_test_pb_bddc_poisson_params_t

  ! Types
  public :: par_test_pb_bddc_poisson_params_t, standard_bddc, pb_bddc

contains

  !==================================================================================================
  subroutine par_test_pb_bddc_poisson_params_define_parameters()
    implicit none
    ! Common
    call parameter_handler%add(dir_path_out_key, '--dir_path_out', '.', 'Default output directory', switch_ab='-dpo')
    call parameter_handler%add(reference_fe_order_key, '--reference-fe-order', 1, 'Order of the fe space reference fe', switch_ab='-order')
    call parameter_handler%add(write_solution_key, '--write-solution', .false., 'Write solution in VTK format', switch_ab='-wsolution')

    ! Specific
    call parameter_handler%add(write_matrices_key, '--write-matrices', .false., 'Write local-to-subdomain sparse matrices  in matrix market format', switch_ab='-wmatrices') 
    call parameter_handler%add(jump_key, '--jump', 1, 'Jump of physical parameter in the inclusion', switch_ab='-j')
    call parameter_handler%add(inclusion_key, '--inclusion', 1, 'Inclusion type', switch_ab='-i')
    call parameter_handler%add(coarse_fe_handler_type_key, '--coarse-fe-handler', pb_bddc, 'Which coarse fe handler to use?', switch_ab='-coarse-handler')
    call parameter_handler%add(nchannel_x_direction_key, '--nchannel_x_direction', [1,1,1], 'Number of channels per direction', switch_ab='-nc')
    call parameter_handler%add(nparts_with_channels_key, '--nparts_with_channels', [1,1,1], 'Number of parts per with channels', switch_ab='-npwc')
    call parameter_handler%add(size_sub_object_key, '--size_sub_object', 1, 'Size of subobject in number of mesh size (must be an integer)', switch_ab='-sso')
        
  end subroutine par_test_pb_bddc_poisson_params_define_parameters

  !==================================================================================================
  subroutine process_parameters(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in)  :: this
    call parameter_handler%process_parameters(par_test_pb_bddc_poisson_params_define_parameters)
  end subroutine process_parameters

  !==================================================================================================
  function get_parameter_list(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    type(ParameterList_t), pointer                        :: get_parameter_list
    get_parameter_list  => parameter_handler%get_values()
  end function get_parameter_list

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable                    :: get_dir_path_out
    call parameter_handler%GetAsString(key = dir_path_out_key, string = get_dir_path_out)
  end function get_dir_path_out
  
  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable                    :: get_prefix
    call parameter_handler%GetAsString(key = prefix_key, string = get_prefix)
  end function get_prefix
 
  !==================================================================================================
  function get_output_handler_dir_path(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable                    :: get_output_handler_dir_path
    call parameter_handler%GetAsString(key = output_handler_dir_path_key, string = get_output_handler_dir_path)
  end function get_output_handler_dir_path

  !==================================================================================================
  function get_output_handler_prefix(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable                    :: get_output_handler_prefix
    call parameter_handler%GetAsString(key = output_handler_prefix_key, string = get_output_handler_prefix)
  end function get_output_handler_prefix
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_reference_fe_order
    call parameter_handler%Get(key = reference_fe_order_key, Value = get_reference_fe_order)
  end function get_reference_fe_order
  
  !==========================================================================================par_test_pb_bddc_poisson_params_t========
  function get_write_solution(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    logical                                               :: get_write_solution
    call parameter_handler%Get(key = write_solution_key, Value = get_write_solution)
  end function get_write_solution
  
  !==========================================================================================par_test_pb_bddc_poisson_params_t========
  function get_write_matrices(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    logical                                               :: get_write_matrices
    call parameter_handler%Get(key = write_matrices_key, Value = get_write_matrices)
  end function get_write_matrices

  !==================================================================================================
  function get_jump(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_jump
    call parameter_handler%Get(key = jump_key, Value = get_jump)
  end function get_jump

  !==================================================================================================
  function get_inclusion(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_inclusion
    call parameter_handler%Get(key = inclusion_key, Value = get_inclusion)
  end function get_inclusion
  
  !==================================================================================================
  function get_coarse_fe_handler_type(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    character(len=:),      allocatable                    :: get_coarse_fe_handler_type
    call parameter_handler%GetAsString(key = coarse_fe_handler_type_key, string = get_coarse_fe_handler_type)
  end function get_coarse_fe_handler_type 
  
  !==================================================================================================
  function get_nchannel_x_direction(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_nchannel_x_direction(3)
    call parameter_handler%Get(key = nchannel_x_direction_key, Value = get_nchannel_x_direction)
  end function get_nchannel_x_direction

  !==================================================================================================
  function get_nparts_with_channels(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_nparts_with_channels(3)
    call parameter_handler%Get(key = nparts_with_channels_key, Value = get_nparts_with_channels)
  end function get_nparts_with_channels

  !==================================================================================================
  function get_nparts(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                   :: get_nparts(3)
    integer(ip), allocatable                      :: num_parts_x_dir(:)
    integer(ip), allocatable                      :: array_size(:)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    call parameter_handler%GetAsArray(key = struct_hex_mesh_generator_num_parts_x_dim_x_level_key, Value = num_parts_x_dir)
    get_nparts=num_parts_x_dir(1:3)
  end function get_nparts
  
  !==================================================================================================
  function get_hex_mesh_domain_limits(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    real(rp)                                              :: get_hex_mesh_domain_limits(6)
    call parameter_handler%Get(key = struct_hex_mesh_generator_domain_limits_key, Value = get_hex_mesh_domain_limits)
  end function get_hex_mesh_domain_limits
  
  !==================================================================================================
    function get_num_cells_x_dir(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: num_levels
    integer(ip)                                           :: get_num_cells_x_dir(3)
    integer(ip), allocatable                              :: num_cells_x_dir(:)
    integer(ip)                                           :: error
    call parameter_handler%GetAsArray(key = struct_hex_mesh_generator_num_cells_x_dim_key, Value = num_cells_x_dir)
    get_num_cells_x_dir = 0
    get_num_cells_x_dir(1:size(num_cells_x_dir))=num_cells_x_dir(:)
  end function get_num_cells_x_dir
  !==================================================================================================
  function get_size_sub_object(this)
    implicit none
    class(par_test_pb_bddc_poisson_params_t) , intent(in) :: this
    integer(ip)                                           :: get_size_sub_object
    call parameter_handler%Get(key = size_sub_object_key, Value = get_size_sub_object)
  end function get_size_sub_object
  
  !==================================================================================================

end module pb_bddc_poisson_params_names
