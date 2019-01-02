! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fempar_parameter_handler_names
  use types_names
  use iterative_linear_solver_parameters_names
  use direct_solver_parameters_names
  use environment_names
  use mesh_distribution_names
  use metis_interface_names
  use triangulation_names
  use uniform_hex_mesh_generator_names
  use fe_space_names
  use parameter_handler_names
  use reference_fe_names
  use FPL
# include "debug.i90"
  implicit none
  private

  character(*), parameter :: FLAP_HELP_MESSAGE_TABULATOR = "    "
  character(*), parameter :: BRK_LINE = NEW_LINE(FLAP_HELP_MESSAGE_TABULATOR) // FLAP_HELP_MESSAGE_TABULATOR
 
  type, extends(parameter_handler_t) :: fempar_parameter_handler_t 
     private 
     procedure(define_user_parameters), pointer :: define_user_parameters => NULL()
   contains
     procedure                           :: process_parameters                      => fph_process_parameters
     procedure, private, non_overridable :: define_fempar_parameters                => fph_define_fempar_parameters
     procedure                           :: define_parameters                       => fph_define_parameters
     procedure                           :: free                                    => fph_free
     procedure                           :: get_dir_path => fph_get_dir_path
     procedure                           :: get_dir_path_out => fph_get_dir_path_out
     procedure                           :: get_prefix => fph_get_prefix
  end type fempar_parameter_handler_t
  
  interface
      subroutine define_user_parameters(this)
        import :: fempar_parameter_handler_t 
        class(fempar_parameter_handler_t ), intent(inout)    :: this
      end subroutine define_user_parameters
  end interface
  public :: fempar_parameter_handler_t
contains
  subroutine fph_process_parameters(this,define_user_parameters_procedure,parse_cla)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    procedure(define_user_parameters), optional :: define_user_parameters_procedure
    logical,            optional, intent(in)    :: parse_cla         !< Parse command line arguments
    logical :: parse_cla_
    parse_cla_ = .true.
    if ( present(parse_cla) ) parse_cla_ = parse_cla
    call this%free()
    if ( parse_cla_ ) then
      call this%init_cli()
    end if    
    call this%initialize_lists()
    call this%define_fempar_parameters()
    call this%define_parameters()
    if (present(define_user_parameters_procedure)) then
      this%define_user_parameters => define_user_parameters_procedure
      call this%define_user_parameters()
    end if
    
#ifdef DEBUG
    call this%assert_lists_consistency()
#endif
    if ( parse_cla_ ) then 
      call this%add_to_cli()
      call this%parse()
    end if
  end subroutine fph_process_parameters
  
  subroutine fph_define_parameters(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    !call this%define_fempar_parameters()
    !if (associated(this%define_user_parameters)) then
    !  call this%define_user_parameters()
    !end if
  end subroutine fph_define_parameters
  
  subroutine fph_free(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    call parameter_handler_free(this) ! this%parameter_handler_t%free()
    nullify(this%define_user_parameters)
  end subroutine fph_free  
  
  !subroutine fph_set_define_user_parameters_procedure(this, define_user_parameters_procedure)
  !  implicit none
  !  class(fempar_parameter_handler_t), intent(inout) :: this
  !  procedure(define_user_parameters) :: define_user_parameters_procedure
  !  this%define_user_parameters => define_user_parameters_procedure
  !  !call this%define_user_parameters()
  !end subroutine fph_set_define_user_parameters_procedure  
  
  subroutine fph_define_fempar_parameters(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    type(parameterlist_t), pointer :: values, switches, switches_ab, helpers, required 
    integer(ip) :: error, i
    character(len=:), allocatable :: help_string
#ifdef UMFPACK
    integer(ip) :: umfpack_control(UMFPACK_INFO)
#endif 

    values      => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()

    error = 0


    ! General
    ! We should use some kind of style to define all these keys. We should include the acceptable values and explain 
    ! what do they mean. In some keys we use magic numbers 0, 1, etc instead of descriptive labels. 
    ! We could also eliminate the switch_ab, since it is not very practical, and when we will have many keys, it will
    ! be hard to get short abbreviations without conflict.
    help_string = "Directory of the source files" // BRK_LINE // & 
                  "To read from GiD mesh" // BRK_LINE // &
                  "To read from GiD mesh" 
                  
    error = error + helpers%set(key = dir_path_key     , Value= help_string)
    error = error + switches%set(key = dir_path_key    , Value= '--DIR_PATH')
    error = error + values%set(key = dir_path_key      , Value= '.')

    error = error + helpers%set(key = prefix_key     , Value= 'Name of the GiD files')
    error = error + switches%set(key = prefix_key    , Value= '--PREFIX')
    error = error + values%set(key = prefix_key      , Value= 'A')

    error = error + helpers%set(key = dir_path_out_key     , Value= 'Output Directory')
    error = error + switches%set(key = dir_path_out_key    , Value= '--DIR_PATH_OUT')
    error = error + values%set(key = dir_path_out_key      , Value= '.')

   
    ! Iterative linear solver keys
    error = error + helpers%set(key = ils_type_key     , value = 'Iterative linear solver type')
    error = error + switches%set(key = ils_type_key    , value = '--ILS_TYPE')
    error = error + values%set(key = ils_type_key      , value = rgmres_name )

    error = error + helpers%set(key = ils_rtol_key     , value = 'Relative tolerance for the stopping criteria')
    error = error + switches%set(key = ils_rtol_key    , value = '--ILS_RTOL')
    error = error + values%set(key = ils_rtol_key      , value = default_rtol )

    error = error + helpers%set(key = ils_atol_key     , value = 'Absolute tolerance for the stopping criteria')
    error = error + switches%set(key = ils_atol_key    , value = '--ILS_ATOL')
    error = error + values%set(key = ils_atol_key      , value = default_atol )

    error = error + helpers%set(key = ils_stopping_criterium_key     , value = 'Stopping criterium type')
    error = error + switches%set(key = ils_stopping_criterium_key    , value = '--ILS_STOPPING_CRITERIUM')
    error = error + values%set(key = ils_stopping_criterium_key      , value = default_rgmres_stopping_criteria )

    error = error + helpers%set(key = ils_output_frequency_key     , value = 'Frequency for output printing') 
    error = error + switches%set(key = ils_output_frequency_key    , value = '--ILS_OUTPUT_FREQUENCY') 
    error = error + values%set(key = ils_output_frequency_key      , value = default_output_frequency ) 

    error = error + helpers%set(key = ils_max_num_iterations_key     , value = 'Maximum number of iterations')
    error = error + switches%set(key = ils_max_num_iterations_key    , value = '--ILS_MAX_NUM_ITERATIONS')
    error = error + values%set(key = ils_max_num_iterations_key      , value = default_max_num_iterations )

    error = error + helpers%set(key = ils_track_convergence_history_key     , value = 'Track convergence history')
    error = error + switches%set(key = ils_track_convergence_history_key    , value = '--ILS_TRACK_CONVERGENCE_HISTORY')
    error = error + values%set(key = ils_track_convergence_history_key      , value = default_track_convergence_history )

    error = error + helpers%set(key = ils_max_dim_krylov_basis_key     , value = 'Maximum dimension of Krylov basis') 
    error = error + switches%set(key = ils_max_dim_krylov_basis_key    , value = '--ILS_MAX_DIM_KRYLOV_BASIS') 
    error = error + values%set(key = ils_max_dim_krylov_basis_key      , value = default_dkrymax ) 

    error = error + helpers%set(key = ils_orthonorm_strategy_key     , value = 'Orthonormalization strategy (for GMRES)')
    error = error + switches%set(key = ils_orthonorm_strategy_key    , value = '--ILS_ORTHONORM_STRATEGY')
    error = error + values%set(key = ils_orthonorm_strategy_key      , value = default_orthonorm_strat)

    error = error + helpers%set(key = ils_relaxation_key     , value = 'Relaxation value (for Richardson)')
    error = error + switches%set(key = ils_relaxation_key    , value = '--ILS_RELAXATION')
    error = error + values%set(key = ils_relaxation_key      , value = default_richardson_relaxation )

    !error = error + helpers%set(key = ils_luout_key     , value = 'Write unit for solver report')
    !error = error + switches%set(key = ils_luout_key    , value = '--sol_out')
    !error = error + values%set(key = ils_luout_key      , value = default_luout )

    ! Sparse direct solver keys
    error = error + helpers%set(key = dls_type_key     , value = 'Direct solver type')
    error = error + switches%set(key = dls_type_key    , value = '--DLS_TYPE   ')
    error = error + values%set(key = dls_type_key      , value = pardiso_mkl )

    ! PARDISO MKL default values
    error = error + helpers%set(key = pardiso_mkl_message_level     , value = 'PARDISO message level')
    error = error + switches%set(key = pardiso_mkl_message_level    , value = '--PARDISO_messg_lev')
    error = error + values%set(key = pardiso_mkl_message_level      , value = pardiso_mkl_default_message_level )

    error = error + helpers%set(key = pardiso_mkl_iparm     , value = 'PARDISO parameters')     
    error = error + switches%set(key = pardiso_mkl_iparm    , value = '--PARDISO_params')     
    error = error + values%set(key = pardiso_mkl_iparm      , value = [ (pardiso_mkl_default_iparm, i=1,64) ] )     

    error = error + helpers%set(key = pardiso_mkl_matrix_type     , value = 'PARDISO matrix type')
    error = error + switches%set(key = pardiso_mkl_matrix_type    , value = '--PARDISO_mat_type')
    error = error + values%set(key = pardiso_mkl_matrix_type      , value = pardiso_mkl_default_matrix_type )
    
    ! UMFPACK keys
#ifdef UMFPACK
    call umfpack_di_defaults(umfpack_control)
    error = error + helpers%set(key = umfpack_control_params     , value = 'UMFPACK control parameters')
    error = error + switches%set(key = umfpack_control_params    , value = '--UMFPACK_cont_params')
    error = error + values%set(key = umfpack_control_params      , value = umfpack_control )
#endif 

    ! Finite element space 

    ! New keys to create FE space with ParameterList
    ! Reference FE types ( lagrangian, edge, face, B-spline...)

    ! Number of fields in the FE space
    error = error + helpers%set(key = fes_num_fields_key, value = 'Finite element space number of fields')
    error = error + switches%set(key = fes_num_fields_key, value = '--FES_NUM_FIELDS')
    error = error + values%set(key = fes_num_fields_key, value = 1)
    
    ! Number of reference FEs
    error = error + helpers%set(key = fes_num_ref_fes_key, value = 'Finite element space number of fields')
    error = error + switches%set(key = fes_num_ref_fes_key, value = '--FES_NUM_REF_FES')
    error = error + values%set(key = fes_num_ref_fes_key, value = 1)

    ! set_ids_to_reference_fes: Given a field ID and cell set ID, returns the desired reference FE ID
    error = error + helpers%set(key = fes_set_ids_ref_fes_key, value = 'Set IDs to reference FEs for every field')
    error = error + switches%set(key = fes_set_ids_ref_fes_key, value = '--FES_SET_IDS_REF_FES')
    error = error + values%set(key = fes_set_ids_ref_fes_key, value = [ 1 ])

    ! Reference FE IDs    
    error = error + helpers%set(key = fes_ref_fe_types_key, value = 'Reference finite element types')
    error = error + switches%set(key = fes_ref_fe_types_key, value = '--FES_REF_FE_TYPES')
    error = error + values%set(key = fes_ref_fe_types_key, value = fe_type_lagrangian)

    ! Reference FE order (0,1,2,...) fe_space_orders_key
    error = error + helpers%set(key = fes_ref_fe_orders_key, value = 'Reference finite element orders')
    error = error + switches%set(key = fes_ref_fe_orders_key, value = '--FES_REF_FE_ORDERS')
    error = error + values%set(key = fes_ref_fe_orders_key, value = [ 1 ])

    ! FE space conformities ( true = no face integration needed, false = face integration required )
    error = error + helpers%set(key = fes_ref_fe_conformities_key, value = 'Finite element space conformities')
    error = error + switches%set(key = fes_ref_fe_conformities_key, value = '--FES_REF_FE_CONFORMITIES')
    error = error + values%set(key = fes_ref_fe_conformities_key, value = [ .true. ])

    ! FE space continuities ( true = continuous FE space, false = otherwise )
    error = error + helpers%set(key = fes_ref_fe_continuities_key, value = 'Finite element space continuities')
    error = error + switches%set(key = fes_ref_fe_continuities_key, value = '--FES_REF_FE_CONTINUITIES')
    error = error + values%set(key = fes_ref_fe_continuities_key, value = [ .true. ])

    ! FE field types (scalar, vector, tensor)
    error = error + helpers%set(key = fes_field_types_key, value = 'Finite element space field types')
    error = error + switches%set(key = fes_field_types_key, value = '--FES_FIELD_TYPES')
    error = error + values%set(key = fes_field_types_key, value =  field_type_scalar )

    ! FE field blocks (scalar, vector, tensor)
    error = error + helpers%set(key = fes_field_blocks_key, value = 'Finite element space field blocks')
    error = error + switches%set(key = fes_field_blocks_key, value = '--FES_FIELD_BLOCKS')
    error = error + values%set(key = fes_field_blocks_key, value = [ 1 ])

    ! FE space construction type homogeneous/heterogeneous ( .true. = all cells same reference fe, .false. = otherwise ) 
    error = error + helpers%set(key = fes_same_ref_fes_all_cells_key, value = 'Finite element space fixed reference fe logical')
    error = error + switches%set(key = fes_same_ref_fes_all_cells_key, value = '--FES_SAME_REFS_ALL_CELLS')
    error = error + values%set(key = fes_same_ref_fes_all_cells_key, value = .true.)

    error = error + helpers%set(key = coarse_space_use_vertices_key     , value = 'Shape functions on vertices')
    error = error + switches%set(key = coarse_space_use_vertices_key    , value = '--FES_COARSE_SPACE_USE_VERTICES')
    error = error + values%set(key = coarse_space_use_vertices_key      , value = .true. )

    error = error + helpers%set(key = coarse_space_use_edges_key     , value = 'Shape functions on edges')
    error = error + switches%set(key = coarse_space_use_edges_key    , value = '--FES_COARSE_SPACE_USE_EDGES')
    error = error + values%set(key = coarse_space_use_edges_key      , value = .true. )

    error = error + helpers%set(key = coarse_space_use_faces_key     , value = 'Shape functions on faces')
    error = error + switches%set(key = coarse_space_use_faces_key    , value = '--FES_COARSE_SPACE_USE_FACES')
    error = error + values%set(key = coarse_space_use_faces_key      , value = .true. )

    ! Environment
    error = error + helpers%set(key = environment_type_key     , value = 'Type of environment')
    error = error + switches%set(key = environment_type_key    , value = '--ENV_TYPE')
    error = error + values%set(key = environment_type_key      , value = structured )

    ! Partitioner
    error = error + helpers%set(key = num_parts_key     , value = 'Number of parts to split mesh with a graph partitioner') 
    error = error + switches%set(key = num_parts_key    , value = '--PART_NUM_PARTS') 
    error = error + values%set(key = num_parts_key      , value = 1) 

    error = error + helpers%set(key = num_levels_distribution_key     , value = 'Number of levels of the parallel distribution') 
    error = error + switches%set(key = num_levels_distribution_key    , value = '--PART_NUM_LEVELS_DISTRIBUTION') 
    error = error + values%set(key = num_levels_distribution_key      , value = 1) 

    error = error + helpers%set(key = num_parts_x_level_key     , value = 'Number of parts per level') 
    error = error + switches%set(key = num_parts_x_level_key    , value = '--PART_NUM_PARTS_X_LEVEL') 
    error = error + values%set(key = num_parts_x_level_key      , value = [1]) 

    error = error + helpers%set(key = debug_key     , value = 'Debug key for partitioner') 
    error = error + switches%set(key = debug_key    , value = '--PART_DEBUG') 
    error = error + values%set(key = debug_key      , value = 0) 

    error = error + helpers%set(key = strategy_key     , value = 'Strategy key for partitioner') 
    error = error + switches%set(key = strategy_key    , value = '--PART_STRATEGY') 
    error = error + values%set(key = strategy_key      , value = part_kway ) 

    ! METIS keys
    error = error + helpers%set(key = metis_option_debug_key     , value = 'METIS debug key') 
    error = error + switches%set(key = metis_option_debug_key    , value = '--METIS_DEBUG') 
    error = error + values%set(key = metis_option_debug_key      , value = 2 ) 

    error = error + helpers%set(key = metis_option_ufactor_key     , value = 'METIS option ufactor') 
    error = error + switches%set(key = metis_option_ufactor_key    , value = '--METIS_OPTION_UFACTOR') 
    error = error + values%set(key = metis_option_ufactor_key      , value = 30 ) 

    error = error + helpers%set(key = metis_option_minconn_key     , value = 'METIS option minconn') 
    error = error + switches%set(key = metis_option_minconn_key    , value = '--METIS_OPTION_MINCONN') 
    error = error + values%set(key = metis_option_minconn_key      , value = 0 ) 

    error = error + helpers%set(key = metis_option_contig_key     , value = 'METIS option config') 
    error = error + switches%set(key = metis_option_contig_key    , value = '--METIS_OPTION_CONFIG')
    error = error + values%set(key = metis_option_contig_key      , value = 1 ) 

    error = error + helpers%set(key = metis_option_ctype_key     , value = 'METIS option ctype')
    error = error + switches%set(key = metis_option_ctype_key    , value = '--METIS_OPTION_CTYPE') 
    error = error + values%set(key = metis_option_ctype_key      , value = METIS_CTYPE_SHEM ) 

    error = error + helpers%set(key = metis_option_iptype_key     , value = 'METIS option iptype')
    error = error + switches%set(key = metis_option_iptype_key    , value = '--METIS_OPTION_IPTYPE')
    error = error + values%set(key = metis_option_iptype_key      , value = METIS_IPTYPE_EDGE )


    ! Triangulation keys
    error = error + helpers%set(key = triang_generate_key     , Value = 'Way to generate the triangulation')
    error = error + switches%set(key = triang_generate_key    , Value = '--TRIANG_GENERATE')
    error = error + values%set(key = triang_generate_key      , Value = triangulation_generate_structured )
    
    error = error + helpers%set(key = triang_geometric_interpolation_order_key     , value = 'Interpolation order for geometrical mapping' )
    error = error + switches%set(key = triang_geometric_interpolation_order_key    , value = '--TRIANG_GEOMETRIC_INTERPOLATION_ORDER' )
    error = error + values%set(key = triang_geometric_interpolation_order_key      , value = 1 )

    ! Uniform hexahedral mesh keys
    error = error + helpers%set(key = struct_hex_triang_num_dims_key     , value = 'Number of space dimensions')                   
    error = error + switches%set(key = struct_hex_triang_num_dims_key    , value = '--STRUCT_HEX_TRIANG_NUM_DIMS')             
    error = error + values%set(key = struct_hex_triang_num_dims_key      , value = 2)                   

    error = error + helpers%set(key = struct_hex_triang_num_cells_dir     , value = 'Number of cells per each dimension')  
    error = error + switches%set(key = struct_hex_triang_num_cells_dir    , value = '--STRUCT_HEX_TRIANG_NUM_CELLS_DIM')
    error = error + values%set(key = struct_hex_triang_num_cells_dir      , value = [10,10,10]) 

    error = error + helpers%set(key = struct_hex_triang_is_dir_periodic_key     , value = 'Is the mesh periodic for every dimension')           
    error = error + switches%set(key = struct_hex_triang_is_dir_periodic_key    , value = '--STRUCT_HEX_TRIANG_IS_DIR_PERIODIC')    
    error = error + values%set(key = struct_hex_triang_is_dir_periodic_key      , value = [0,0,0])           

    error = error + helpers%set(key = struct_hex_triang_num_levels_key     , value = 'Number of levels')                   
    error = error + switches%set(key = struct_hex_triang_num_levels_key    , value = '--STRUCT_HEX_TRIANG_NUM_LEVELS')          
    error = error + values%set(key = struct_hex_triang_num_levels_key      , value = 1)                   

    error = error + helpers%set(key = struct_hex_triang_num_parts_x_dir_key     , value = 'Number of parts per each dimension')
    error = error + switches%set(key = struct_hex_triang_num_parts_x_dir_key    , value = '--STRUCT_HEX_TRIANG_NUM_PARTS_X_DIR')
    error = error + values%set(key = struct_hex_triang_num_parts_x_dir_key      , value = [1,1,1])


    error = error + helpers%set(key = struct_hex_triang_domain_limits_key     , value = 'Domain interval per direction')
    error = error + switches%set(key = struct_hex_triang_domain_limits_key    , value = '--STRUCT_HEX_TRIANG_DOMAIN_LIMITS')
    error = error + values%set(key = struct_hex_triang_domain_limits_key      , value = [0.0,1.0,0.0,1.0,0.0,1.0] )

!    ! BDDC
!    ! I would say that all this should be handled by groups, one for the Dirichlet problem, 
!    ! one for the Neumann problem, one for the coarse problem, etc
!    !error = error + values%set(key = mlbddc_coarse_matrix_symmetric_storage , Value= .false. )
!    !error = error + values%set(key = mlbddc_coarse_matrix_is_symmetric      , Value= .false. )
!    !error = error + values%set(key = mlbddc_coarse_matrix_sign              , Value= SPARSE_MATRIX_SIGN_UNKNOWN )
!
!    ! Output Handler
!    !error = error + values%set(key = oh_staticgrid , Value= )
!    !error = error + values%set(key = xh5_Strategy  , Value= )
!    !error = error + values%set(key = xh5_Info      , Value= )
!
!    ! Sample
!    !error = error + values%set(key = , Value= )
    
    !if associated ( this%define_user_parameters) then
    !  call this%define_user_parameters()
    !end if 

  end subroutine fph_define_fempar_parameters

    ! GETTERS *****************************************************************************************
  function fph_get_dir_path(this)
    implicit none
    class(fempar_parameter_handler_t) , intent(in) :: this
    character(len=:),      allocatable            :: fph_get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = fph_get_dir_path)
    assert(error==0)
  end function fph_get_dir_path 
  
  ! GETTERS *****************************************************************************************
  function fph_get_dir_path_out(this)
    implicit none
    class(fempar_parameter_handler_t) , intent(in) :: this
    character(len=:),      allocatable            :: fph_get_dir_path_out
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_out_key, 'string'))
    error = list%GetAsString(key = dir_path_out_key, string = fph_get_dir_path_out)
    assert(error==0)
  end function fph_get_dir_path_out

  !==================================================================================================
  function fph_get_prefix(this)
    implicit none
    class(fempar_parameter_handler_t) , intent(in) :: this
    character(len=:),      allocatable            :: fph_get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_values()
    assert(list%isAssignable(prefix_key, 'string'))
    error = list%GetAsString(key = prefix_key, string = fph_get_prefix)
    assert(error==0)
  end function fph_get_prefix


end module fempar_parameter_handler_names 
