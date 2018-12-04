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
   contains
     procedure, private, non_overridable :: define_fempar_parameters => fph_define_fempar_parameters
     procedure                           :: define_parameters        => fph_define_parameters
     procedure                           :: define_user_parameters   => fph_define_user_parameters 
  end type fempar_parameter_handler_t

  public :: fempar_parameter_handler_t

contains
  subroutine fph_define_parameters(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    call this%define_fempar_parameters()
    call this%define_user_parameters()
  end subroutine fph_define_parameters
 
  subroutine fph_define_user_parameters(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this

  end subroutine fph_define_user_parameters

  
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
    error = error + switches%set(key = dir_path_key    , Value= '--dir_path')
    error = error + switches_ab%set(key = dir_path_key , Value= '-d')
    error = error + values%set(key = dir_path_key      , Value= '.')

    error = error + helpers%set(key = prefix_key     , Value= 'Name of the GiD files')
    error = error + switches%set(key = prefix_key    , Value= '--prefix')
    error = error + switches_ab%set(key = prefix_key , Value= '-p')
    error = error + values%set(key = prefix_key      , Value= '')

    error = error + helpers%set(key = dir_path_out_key     , Value= 'Output Directory')
    error = error + switches%set(key = dir_path_out_key    , Value= '--dir_path_out')
    error = error + switches_ab%set(key = dir_path_out_key , Value= '-o')
    error = error + values%set(key = dir_path_out_key      , Value= '.')

   
    ! Iterative linear solver keys
    error = error + helpers%set(key = ils_type     , value = 'Iterative solver type')
    error = error + switches%set(key = ils_type    , value = '--iter-sol')
    error = error + switches_ab%set(key = ils_type , value = '-is')
    error = error + values%set(key = ils_type      , value = rgmres_name )

    error = error + helpers%set(key = ils_rtol     , value = 'Relative tolerance for the stopping criteria')
    error = error + switches%set(key = ils_rtol    , value = '--rel_tol')
    error = error + switches_ab%set(key = ils_rtol , value = '-rt')
    error = error + values%set(key = ils_rtol      , value = default_rtol )

    error = error + helpers%set(key = ils_atol     , value = 'Absolute tolerance for the stopping criteria')
    error = error + switches%set(key = ils_atol    , value = '--abs_tol')
    error = error + switches_ab%set(key = ils_atol , value = '-at')
    error = error + values%set(key = ils_atol      , value = default_atol )

    error = error + helpers%set(key = ils_stopping_criteria     , value = 'Stopping criterium type')
    error = error + switches%set(key = ils_stopping_criteria    , value = '--stop_crit')
    error = error + switches_ab%set(key = ils_stopping_criteria , value = '-st')
    error = error + values%set(key = ils_stopping_criteria      , value = default_rgmres_stopping_criteria )

    error = error + helpers%set(key = ils_output_frequency     , value = 'Frequency for output printing') 
    error = error + switches%set(key = ils_output_frequency    , value = '--out_freq') 
    error = error + switches_ab%set(key = ils_output_frequency , value = '-of') 
    error = error + values%set(key = ils_output_frequency      , value = default_output_frequency ) 

    error = error + helpers%set(key = ils_max_num_iterations     , value = 'Maximum number of iterations')
    error = error + switches%set(key = ils_max_num_iterations    , value = '--max_num_iter')
    error = error + switches_ab%set(key = ils_max_num_iterations , value = '-mi')
    error = error + values%set(key = ils_max_num_iterations      , value = default_max_num_iterations )

    error = error + helpers%set(key = ils_track_convergence_history     , value = 'Track convergence history')
    error = error + switches%set(key = ils_track_convergence_history    , value = '--track_conv')
    error = error + switches_ab%set(key = ils_track_convergence_history , value = '-tc')
    error = error + values%set(key = ils_track_convergence_history      , value = default_track_convergence_history )

    error = error + helpers%set(key = ils_dkrymax     , value = 'Maximum dimension of Krylov basis') 
    error = error + switches%set(key = ils_dkrymax    , value = '--dim_kry') 
    error = error + switches_ab%set(key = ils_dkrymax , value = '-dk') 
    error = error + values%set(key = ils_dkrymax      , value = default_dkrymax ) 

    error = error + helpers%set(key = ils_orthonorm_strat     , value = 'Orthonormalization strategy (for GMRES)')
    error = error + switches%set(key = ils_orthonorm_strat    , value = '--ort_str')
    error = error + switches_ab%set(key = ils_orthonorm_strat , value = '-os')
    error = error + values%set(key = ils_orthonorm_strat      , value = default_orthonorm_strat)

    error = error + helpers%set(key = ils_relaxation     , value = 'Relaxation value (for Richardson)')
    error = error + switches%set(key = ils_relaxation    , value = '--rel_val')
    error = error + switches_ab%set(key = ils_relaxation , value = '-rv')
    error = error + values%set(key = ils_relaxation      , value = default_richardson_relaxation )

    error = error + helpers%set(key = ils_luout     , value = 'Write unit for solver report')
    error = error + switches%set(key = ils_luout    , value = '--sol_out')
    error = error + switches_ab%set(key = ils_luout , value = '-so')
    error = error + values%set(key = ils_luout      , value = default_luout )

    ! Sparse direct solver
    error = error + helpers%set(key = direct_solver_type     , value = 'Direct solver type')
    error = error + switches%set(key = direct_solver_type    , value = '--dir_sol')
    error = error + switches_ab%set(key = direct_solver_type , value = '-ds')
    error = error + values%set(key = direct_solver_type      , value = pardiso_mkl )

    ! PARDISO MKL default values
    error = error + helpers%set(key = pardiso_mkl_message_level     , value = 'PARDISO message level')
    error = error + switches%set(key = pardiso_mkl_message_level    , value = '--PARDISO_messg_lev')
    error = error + switches_ab%set(key = pardiso_mkl_message_level , value = '-pml')
    error = error + values%set(key = pardiso_mkl_message_level      , value = pardiso_mkl_default_message_level )

    error = error + helpers%set(key = pardiso_mkl_iparm     , value = 'PARDISO parameters')     
    error = error + switches%set(key = pardiso_mkl_iparm    , value = '--PARDISO_params')     
    error = error + switches_ab%set(key = pardiso_mkl_iparm , value = '-pp')     
    error = error + values%set(key = pardiso_mkl_iparm      , value = [ (pardiso_mkl_default_iparm, i=1,64) ] )     

    error = error + helpers%set(key = pardiso_mkl_matrix_type     , value = 'PARDISO matrix type')
    error = error + switches%set(key = pardiso_mkl_matrix_type    , value = '--PARDISO_mat_type')
    error = error + switches_ab%set(key = pardiso_mkl_matrix_type , value = '-pmt')
    error = error + values%set(key = pardiso_mkl_matrix_type      , value = pardiso_mkl_default_matrix_type )

#ifdef UMFPACK
    call umfpack_di_defaults(umfpack_control)
    error = error + helpers%set(key = umfpack_control_params     , value = 'UMFPACK control parameters')
    error = error + switches%set(key = umfpack_control_params    , value = '--UMFPACK_cont_params')
    error = error + switches_ab%set(key = umfpack_control_params , value = '-ucp')
    error = error + values%set(key = umfpack_control_params      , value = umfpack_control )
#endif 

    ! Finite element space 

    ! New keys to create FE space with ParameterList
    ! Reference FE types ( lagrangian, edge, face, B-spline...)

    ! FE field IDs (0, 1, 2...) to be consitent with mesh data
    error = error + helpers%set(key = fe_space_num_fields_key, value = 'Finite element space number of fields')
    error = error + switches%set(key = fe_space_num_fields_key, value = '--FE_SPACE_NUM_FIELD')
    error = error + values%set(key = fe_space_num_fields_key, value = 1)

    error = error + helpers%set(key = reference_fe_type_key, value = 'Reference finite element types')
    error = error + switches%set(key = reference_fe_type_key, value = '--REFERENCE_FE_TYPES')
    error = error + values%set(key = reference_fe_type_key, value = fe_type_lagrangian)

    ! Reference FE order (0,1,2,...)
    error = error + helpers%set(key = reference_fe_orders_key, value = 'Reference finite element orders')
    error = error + switches%set(key = reference_fe_orders_key, value = '--REFERENCE_FE_ORDERS')
    error = error + values%set(key = reference_fe_orders_key, value = [ 1 ])

    ! FE space conformities ( true = no face integration needed, false = face integration required )
    error = error + helpers%set(key = reference_fe_conformity_key, value = 'Finite element space conformities')
    error = error + switches%set(key = reference_fe_conformity_key, value = '--FE__SPACE_CONFORMITIES')
    error = error + values%set(key = reference_fe_conformity_key, value = [ .true. ])

    ! FE space continuities ( true = continuous FE space, false = otherwise )
    error = error + helpers%set(key = reference_fe_continuity_key, value = 'Finite element space continuities')
    error = error + switches%set(key = reference_fe_continuity_key, value = '--FE__SPACE_CONTINUITIES')
    error = error + values%set(key = reference_fe_continuity_key, value = [ .true. ])

    ! FE field types (scalar, vector, tensor)
    error = error + helpers%set(key = fe_space_field_types_key, value = 'Finite element space field types')
    error = error + switches%set(key = fe_space_field_types_key, value = '--FE_SPACE_FIELD_TYPES')
    error = error + values%set(key = fe_space_field_types_key, value =  field_type_scalar )

    ! FE field blocks (scalar, vector, tensor)
    error = error + helpers%set(key = fe_space_field_blocks_key, value = 'Finite element space field blocks')
    error = error + switches%set(key = fe_space_field_blocks_key, value = '--FE_SPACE_FIELD_BLOCKS')
    error = error + values%set(key = fe_space_field_blocks_key, value = [ 1 ])

    ! FE space construction type homogeneous/heterogeneous ( .true. = all cells same reference fe, .false. = otherwise ) 
    error = error + helpers%set(key = fe_space_same_reference_fe_all_cells_key, value = 'Finite element space fixed reference fe logical')
    error = error + switches%set(key = fe_space_same_reference_fe_all_cells_key, value = '--FE__SPACE_SAME_REFERENCE_FE_ALL_CELLS')
    error = error + values%set(key = fe_space_same_reference_fe_all_cells_key, value = .true.)

    ! FE field IDs (0, 1, 2...) to be consitent with mesh data
    error = error + helpers%set(key = fe_space_field_ids_key, value = 'Finite element space field ids')
    error = error + switches%set(key = fe_space_field_ids_key, value = '--FE_SPACE_FIELD_IDS')
    error = error + values%set(key = fe_space_field_ids_key, value = [ 1 ])


    error = error + helpers%set(key = coarse_space_use_vertices_key     , value = 'Shape functions on vertices')
    error = error + switches%set(key = coarse_space_use_vertices_key    , value = '--use_verts')
    error = error + switches_ab%set(key = coarse_space_use_vertices_key , value = '-uv')
    error = error + values%set(key = coarse_space_use_vertices_key      , value = .true. )

    error = error + helpers%set(key = coarse_space_use_edges_key     , value = 'Shape functions on edges')
    error = error + switches%set(key = coarse_space_use_edges_key    , value = '--use_edges')
    error = error + switches_ab%set(key = coarse_space_use_edges_key , value = '-ue')
    error = error + values%set(key = coarse_space_use_edges_key      , value = .true. )

    error = error + helpers%set(key = coarse_space_use_faces_key     , value = 'Shape functions on faces')
    error = error + switches%set(key = coarse_space_use_faces_key    , value = '--use_faces')
    error = error + switches_ab%set(key = coarse_space_use_faces_key , value = '-uf')
    error = error + values%set(key = coarse_space_use_faces_key      , value = .true. )

    ! Environment
    error = error + helpers%set(key = environment_type_key     , value = 'Type of environment')
    error = error + switches%set(key = environment_type_key    , value = '--env_type')
    error = error + switches_ab%set(key = environment_type_key , value = '-et')
    error = error + values%set(key = environment_type_key      , value = structured )

    ! Partitioner
    error = error + helpers%set(key = num_parts_key     , value = 'Number of parts to split mesh with a graph partitioner') 
    error = error + switches%set(key = num_parts_key    , value = '--num_parts_part') 
    error = error + switches_ab%set(key = num_parts_key , value = '-npp') 
    error = error + values%set(key = num_parts_key      , value = 1) 

    error = error + helpers%set(key = num_levels_distribution_key     , value = 'Number of levels of the parallel distribution') 
    error = error + switches%set(key = num_levels_distribution_key    , value = '--num_lev_dist') 
    error = error + switches_ab%set(key = num_levels_distribution_key , value = '-nld') 
    error = error + values%set(key = num_levels_distribution_key      , value = 1) 

    error = error + helpers%set(key = num_parts_x_level_key     , value = 'Number of parts per level') 
    error = error + switches%set(key = num_parts_x_level_key    , value = '--num_pars_lev') 
    error = error + switches_ab%set(key = num_parts_x_level_key , value = '-npl') 
    error = error + values%set(key = num_parts_x_level_key      , value = [1]) 

    error = error + helpers%set(key = debug_key     , value = 'Debug key for partitioner') 
    error = error + switches%set(key = debug_key    , value = '--debug_key_part') 
    error = error + switches_ab%set(key = debug_key , value = '-dkp') 
    error = error + values%set(key = debug_key      , value = 0) 

    error = error + helpers%set(key = strategy_key     , value = 'Strategy key for partitioner') 
    error = error + switches%set(key = strategy_key    , value = '--strat_key') 
    error = error + switches_ab%set(key = strategy_key , value = '-sk') 
    error = error + values%set(key = strategy_key      , value = part_kway ) 

    error = error + helpers%set(key = metis_option_debug_key     , value = 'METIS debug key') 
    error = error + switches%set(key = metis_option_debug_key    , value = '--METIS_debug_key') 
    error = error + switches_ab%set(key = metis_option_debug_key , value = '-mdk') 
    error = error + values%set(key = metis_option_debug_key      , value = 2 ) 

    error = error + helpers%set(key = metis_option_ufactor_key     , value = 'METIS option ufactor') 
    error = error + switches%set(key = metis_option_ufactor_key    , value = '--METIS_option_ufactor') 
    error = error + switches_ab%set(key = metis_option_ufactor_key , value = '-mou') 
    error = error + values%set(key = metis_option_ufactor_key      , value = 30 ) 

    error = error + helpers%set(key = metis_option_minconn_key     , value = 'METIS option minconn') 
    error = error + switches%set(key = metis_option_minconn_key    , value = '--METIS_option_minconn') 
    error = error + switches_ab%set(key = metis_option_minconn_key , value = '-mom') 
    error = error + values%set(key = metis_option_minconn_key      , value = 0 ) 

    error = error + helpers%set(key = metis_option_contig_key     , value = 'METIS option config') 
    error = error + switches%set(key = metis_option_contig_key    , value = '--METIS_option_config') 
    error = error + switches_ab%set(key = metis_option_contig_key , value = '-moc') 
    error = error + values%set(key = metis_option_contig_key      , value = 1 ) 

    error = error + helpers%set(key = metis_option_ctype_key     , value = 'METIS option ctype')
    error = error + switches%set(key = metis_option_ctype_key    , value = '--METIS_option_ctype') 
    error = error + switches_ab%set(key = metis_option_ctype_key , value = '-moct') 
    error = error + values%set(key = metis_option_ctype_key      , value = METIS_CTYPE_SHEM ) 

    error = error + helpers%set(key = metis_option_iptype_key     , value = 'METIS option iptype')
    error = error + switches%set(key = metis_option_iptype_key    , value = '--METIS_option_iptype')
    error = error + switches_ab%set(key = metis_option_iptype_key , value = '-moi')
    error = error + values%set(key = metis_option_iptype_key      , value = METIS_IPTYPE_EDGE )


    ! Triangulation
    error = error + helpers%set(key = triangulation_generate_key     , Value = 'Way to generate the triangulation')
    error = error + switches%set(key = triangulation_generate_key    , Value = '--triang_gen')
    error = error + switches_ab%set(key = triangulation_generate_key , Value = '-tg')
    error = error + values%set(key = triangulation_generate_key      , Value = triangulation_generate_structured )

    ! Uniform hexahedral mesh keys
    error = error + helpers%set(key = num_dims_key     , value = 'Number of space dimensions')                   
    error = error + switches%set(key = num_dims_key    , value = '--num_dims')                   
    error = error + switches_ab%set(key = num_dims_key , value = '-nd')                   
    error = error + values%set(key = num_dims_key      , value = 2)                   

    error = error + helpers%set(key = num_cells_x_dir_key     , value = 'Number of cells per each dimension')  
    error = error + switches%set(key = num_cells_x_dir_key    , value = '--num_cells_dim')
    error = error + switches_ab%set(key = num_cells_x_dir_key , value = '-nc')
    error = error + values%set(key = num_cells_x_dir_key      , value = [10,10,10]) 

    error = error + helpers%set(key = is_dir_periodic_key     , value = 'Is the mesh periodic for each dimension')           
    error = error + switches%set(key = is_dir_periodic_key    , value = '--is_periodic_dim')           
    error = error + switches_ab%set(key = is_dir_periodic_key , value = '-ip')           
    error = error + values%set(key = is_dir_periodic_key      , value = [0,0,0])           

    error = error + helpers%set(key = num_levels_key     , value = 'Number of levels')                   
    error = error + switches%set(key = num_levels_key    , value = '--num_levels')                   
    error = error + switches_ab%set(key = num_levels_key , value = '-nl')                   
    error = error + values%set(key = num_levels_key      , value = 1)                   

    error = error + helpers%set(key = num_parts_x_dir_key     , value = 'Number of parts per each dimension')
    error = error + switches%set(key = num_parts_x_dir_key    , value = '--num_parts_dim')
    error = error + switches_ab%set(key = num_parts_x_dir_key , value = '-np')
    error = error + values%set(key = num_parts_x_dir_key      , value = [1,1,1])

    error = error + helpers%set(key = geometric_interpolation_order_key     , value = 'Interpolation order for geometrical mapping' )
    error = error + switches%set(key = geometric_interpolation_order_key    , value = '--geo_int_order' )
    error = error + switches_ab%set(key = geometric_interpolation_order_key , value =  '-go')
    error = error + values%set(key = geometric_interpolation_order_key      , value = 1 )

    error = error + helpers%set(key = hex_mesh_domain_limits_key     , value = 'Domain interval per direction')
    error = error + switches%set(key = hex_mesh_domain_limits_key    , value = '--dom_int_dir')
    error = error + switches_ab%set(key = hex_mesh_domain_limits_key , value = '-di')
    error = error + values%set(key = hex_mesh_domain_limits_key      , value = [0.0,1.0,0.0,1.0,0.0,1.0] )

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

  end subroutine fph_define_fempar_parameters



end module fempar_parameter_handler_names 
