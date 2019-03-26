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
    use nonlinear_solver_names
    use environment_names
    use mesh_distribution_names
    use metis_interface_names
    use triangulation_names
    use uniform_hex_mesh_generator_names
    use p4est_triangulation_names
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
        procedure                           :: process_parameters       => fph_process_parameters
        procedure, private, non_overridable :: define_fempar_parameters => fph_define_fempar_parameters
        procedure                           :: define_parameters        => fph_define_parameters
        procedure                           :: free                     => fph_free
        procedure                           :: get_dir_path             => fph_get_dir_path
        procedure                           :: get_dir_path_out         => fph_get_dir_path_out
        procedure                           :: get_prefix               => fph_get_prefix
    end type fempar_parameter_handler_t
  
    interface
        subroutine define_user_parameters(this)
            import :: fempar_parameter_handler_t 
            class(fempar_parameter_handler_t ), intent(inout) :: this
        end subroutine define_user_parameters
    end interface

    public :: fempar_parameter_handler_t

contains

    subroutine fph_process_parameters(this,define_user_parameters_procedure,parse_cla)
    !------------------------------------------------------------------
    !< Initialize lists, publish parameters and fill values from CLI
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t),           intent(inout) :: this
        procedure(define_user_parameters), optional                :: define_user_parameters_procedure
        logical,                           optional, intent(in)    :: parse_cla
        logical :: parse_cla_
    !------------------------------------------------------------------
        call this%free()
        parse_cla_ = .true.; if ( present(parse_cla) ) parse_cla_ = parse_cla
        if ( parse_cla_ ) call this%init_cli()
        call this%initialize_lists()
        call this%define_fempar_parameters()
        call this%define_parameters()
        if (present(define_user_parameters_procedure)) then
            this%define_user_parameters => define_user_parameters_procedure
            call this%define_user_parameters()
            else
            nullify(this%define_user_parameters)
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
    !------------------------------------------------------------------
    !< Define parameters procedure
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t), intent(inout) :: this
    !------------------------------------------------------------------
        !call this%define_fempar_parameters()
        !if (associated(this%define_user_parameters)) then
        !    call this%define_user_parameters()
        !end if
    end subroutine fph_define_parameters


    subroutine fph_free(this)
    !------------------------------------------------------------------
    !< Free fempar_parameter_handler_t derived type
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t), intent(inout) :: this
    !------------------------------------------------------------------
        call parameter_handler_free(this) 
        nullify(this%define_user_parameters)
    end subroutine fph_free  
  
    subroutine fph_define_fempar_parameters(this)
    !------------------------------------------------------------------
    !< Publish FEMPAR parameters
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t), intent(inout) :: this
        character(len=:), allocatable                    :: help_string
        integer(ip)                                      :: i
#ifdef UMFPACK
        integer(ip)                                      :: umfpack_control(UMFPACK_INFO)
#endif 
    !------------------------------------------------------------------
        ! General
        ! We should use some kind of style to define all these keys. We should include the acceptable values and explain 
        ! what do they mean. In some keys we use magic numbers 0, 1, etc instead of descriptive labels. 
        ! We could also eliminate the switch_ab, since it is not very practical, and when we will have many keys, it will
        ! be hard to get short abbreviations without conflict.
        help_string = "Directory of the source files" // BRK_LINE // & 
                  "To read from GiD mesh" // BRK_LINE // &
                  "To read from GiD mesh" 
                  
        call this%add(dir_path_key, '--DIR_PATH', '.', help_string)
        call this%add(prefix_key, '--PREFIX', 'A', 'Name of the GiD files')
        call this%add(dir_path_out_key, '--DIR_PATH_OUT', '.', 'Output Directory')

        ! Iterative linear solver keys
        call this%add(ils_type_key, '--ILS_TYPE', rgmres_name,  'Iterative linear solver type')
        call this%add(ils_rtol_key, '--ILS_RTOL', default_rtol, 'Relative tolerance for the stopping criteria')
        call this%add(ils_atol_key, '--ILS_ATOL', default_atol, 'Absolute tolerance for the stopping criteria')

        call this%add(ils_stopping_criterium_key, '--ILS_STOPPING_CRITERIUM', default_rgmres_stopping_criteria, help='Stopping criterium type')
        call this%add(ils_output_frequency_key, '--ILS_OUTPUT_FREQUENCY', default_output_frequency, 'Frequency for output printing') 
        call this%add(ils_max_num_iterations_key, '--ILS_MAX_NUM_ITERATIONS', default_max_num_iterations, 'Maximum number of iterations')
        call this%add(ils_track_convergence_history_key,'--ILS_TRACK_CONVERGENCE_HISTORY', default_track_convergence_history,'Track convergence history')
        call this%add(ils_max_dim_krylov_basis_key, '--ILS_MAX_DIM_KRYLOV_BASIS', default_dkrymax, 'Maximum dimension of Krylov basis') 
        call this%add(ils_orthonorm_strategy_key, '--ILS_ORTHONORM_STRATEGY', default_orthonorm_strat, 'Orthonormalization strategy (for GMRES)')

        call this%add(ils_relaxation_key, '--ILS_RELAXATION', default_richardson_relaxation, 'Relaxation value (for Richardson)')

        !call this%add(ils_luout_key, '--sol_out', default_luout, 'Write unit for solver report')

        ! Sparse direct solver keys
        call this%add(dls_type_key, '--DLS_TYPE', pardiso_mkl, 'Direct solver type')

        ! PARDISO MKL default values
        call this%add(pardiso_mkl_message_level, '--PARDISO_messg_lev', pardiso_mkl_default_message_level, 'PARDISO message level')
        call this%add(pardiso_mkl_iparm, '--PARDISO_params', [ (pardiso_mkl_default_iparm, i=1,64) ], 'PARDISO parameters')     
        call this%add(pardiso_mkl_matrix_type, '--PARDISO_mat_type', pardiso_mkl_default_matrix_type, 'PARDISO matrix type')

        ! UMFPACK keys
#ifdef UMFPACK
        call umfpack_di_defaults(umfpack_control)
        call this%add(umfpack_control_params, '--UMFPACK_cont_params', umfpack_control, 'UMFPACK control parameters')
#endif 

        call this%add(nls_rtol_key, '--NLS_RTOL', default_nls_rtol, 'Relative tolerance for the nonlinear solvers stopping criteria')
        call this%add(nls_atol_key, '--NLS_ATOL', default_nls_atol, 'Absolute tolerance for the nonlinear solvers stopping criteria')

        call this%add(nls_stopping_criterium_key, '--NLS_STOPPING_CRITERIUM', default_nls_stopping_criterium, 'Nonlinear solvers stopping criterium type')
        call this%add(nls_max_num_iterations_key, '--NLS_MAX_NUM_ITERATIONS', default_nls_max_iter, 'Nonlinear solvers maximum number of iterations')
        call this%add(nls_print_iteration_output_key, '--NLS_PRINT_ITERATION_OUTPUT', default_nls_print_iteration_output, 'Print output per nonlinear solver iteration')

        ! Finite element space 

        ! New keys to create FE space with ParameterList
        ! Reference FE types ( lagrangian, edge, face, B-spline...)

        ! Number of fields in the FE space
        call this%add(fes_num_fields_key, '--FES_NUM_FIELDS', 1, 'Finite element space number of fields')

        ! Number of reference FEs
        call this%add(fes_num_ref_fes_key,  '--FES_NUM_REF_FES', 1, 'Finite element space number of fields')

        ! set_ids_to_reference_fes: Given a field ID and cell set ID, returns the desired reference FE ID
        call this%add(fes_set_ids_ref_fes_key, '--FES_SET_IDS_REF_FES', [ 1 ], 'Set IDs to reference FEs for every field')

        ! Reference FE IDs    
        call this%add(fes_ref_fe_types_key, '--FES_REF_FE_TYPES', fe_type_lagrangian, 'Reference finite element types')

        ! Reference FE order (0,1,2,...) fe_space_orders_key
        call this%add(fes_ref_fe_orders_key, '--FES_REF_FE_ORDERS', [ 1 ], 'Reference finite element orders')

        ! FE space conformities ( true = no face integration needed, false = face integration required )
        call this%add(fes_ref_fe_conformities_key, '--FES_REF_FE_CONFORMITIES', [ .true. ], 'Finite element space conformities')

        ! FE space continuities ( true = continuous FE space, false = otherwise )
        call this%add(fes_ref_fe_continuities_key, '--FES_REF_FE_CONTINUITIES', [ .true. ], 'Finite element space continuities')

        ! FE field types (scalar, vector, tensor)
        call this%add(fes_field_types_key, '--FES_FIELD_TYPES', field_type_scalar, 'Finite element space field types')

        ! FE field blocks (scalar, vector, tensor)
        call this%add(fes_field_blocks_key, '--FES_FIELD_BLOCKS', [ 1 ], 'Finite element space field blocks')

        ! FE space construction type homogeneous/heterogeneous ( .true. = all cells same reference fe, .false. = otherwise ) 
        call this%add(fes_same_ref_fes_all_cells_key, '--FES_SAME_REFS_ALL_CELLS', .true., 'Finite element space fixed reference fe logical')
        call this%add(coarse_space_use_vertices_key, '--FES_COARSE_SPACE_USE_VERTICES', .true., 'Shape functions on vertices')
        call this%add(coarse_space_use_edges_key, '--FES_COARSE_SPACE_USE_EDGES', .true., 'Shape functions on edges')
        call this%add(coarse_space_use_faces_key, '--FES_COARSE_SPACE_USE_FACES', .true., 'Shape functions on faces')

        ! Environment
        call this%add(environment_type_key, '--ENV_TYPE', structured, 'Type of environment')

        ! Partitioner
        call this%add(num_parts_key, '--PART_NUM_PARTS', 1, 'Number of parts to split mesh with a graph partitioner') 
        call this%add(num_levels_distribution_key, '--PART_NUM_LEVELS_DISTRIBUTION', 1, 'Number of levels of the parallel distribution') 
        call this%add(num_parts_x_level_key, '--PART_NUM_PARTS_X_LEVEL', [1], 'Number of parts per level') 
        call this%add(debug_key, '--PART_DEBUG', 0, 'Debug key for partitioner') 
        call this%add(strategy_key, '--PART_STRATEGY', part_kway, 'Strategy key for partitioner') 

        ! METIS keys
        call this%add(metis_option_debug_key, '--METIS_DEBUG', 2, 'METIS debug key') 
        call this%add(metis_option_ufactor_key, '--METIS_OPTION_UFACTOR', 30, 'METIS option ufactor') 
        call this%add(metis_option_minconn_key, '--METIS_OPTION_MINCONN', 0, 'METIS option minconn') 
        call this%add(metis_option_contig_key, '--METIS_OPTION_CONFIG', 1, 'METIS option config') 
        call this%add(metis_option_ctype_key, '--METIS_OPTION_CTYPE', METIS_CTYPE_SHEM, 'METIS option ctype')
        call this%add(metis_option_iptype_key, '--METIS_OPTION_IPTYPE', METIS_IPTYPE_EDGE, 'METIS option iptype')

        ! Triangulation keys
        call this%add(triang_generate_key, '--TRIANG_GENERATE', triangulation_generate_structured, 'Way to generate the triangulation')
        call this%add(triang_geometric_interpolation_order_key, '--TRIANG_GEOMETRIC_INTERPOLATION_ORDER', 1, 'Interpolation order for geometrical mapping' )

        ! Uniform hexahedral mesh keys
        call this%add(struct_hex_triang_num_dims_key, '--STRUCT_HEX_TRIANG_NUM_DIMS', 2, 'Number of space dimensions')                   
        call this%add(struct_hex_triang_num_cells_dir, '--STRUCT_HEX_TRIANG_NUM_CELLS_DIM', [10,10,10], 'Number of cells per each dimension')  
        call this%add(struct_hex_triang_is_dir_periodic_key, '--STRUCT_HEX_TRIANG_IS_DIR_PERIODIC', [0,0,0], 'Is the mesh periodic for every dimension')           
        call this%add(struct_hex_triang_num_levels_key, '--STRUCT_HEX_TRIANG_NUM_LEVELS', 1, 'Number of levels')
        call this%add(struct_hex_triang_num_parts_x_dir_key, '--STRUCT_HEX_TRIANG_NUM_PARTS_X_DIR', [1,1,1], 'Number of parts per each dimension')
        call this%add(struct_hex_triang_domain_limits_key, '--STRUCT_HEX_TRIANG_DOMAIN_LIMITS', [0.0,1.0,0.0,1.0,0.0,1.0], 'Domain interval per direction')

        call this%add(p4est_triang_log_level_key, '--P4EST_TRIANG_LOG_LEVEL', FEMPAR_SC_LP_DEFAULT, 'p4est library level of logging output')
        call this%add(p4est_triang_2_1_k_balance_key,'--P4EST_TRIANG_2_1_K_BALANCE', default_p4est_triang_2_1_k_balance,  'value of k for 2:1 k-balanced forest-of-octrees (use with care, at present, only k={0,1} supported/tested)')
        call this%add(p4est_triang_k_ghost_cells_key, '--P4EST_TRIANG_K_GHOST_CELLS', default_p4est_triang_k_ghost_cells, 'value of k for the k-ghost cells set of each processor (k=0 works for any FE space; k>0 should work depending on the FE space, although NOT tested, use with care)')

        ! BDDC
        call this%add(bddc_scaling_function_case_key, '--BDDC_SCALING_FUNCTION_CASE', cardinality, 'Scaling type for the BDDC weighting operator')

        !    ! I would say that all this should be handled by groups, one for the Dirichlet problem, 
        !    ! one for the Neumann problem, one for the coarse problem, etc
        !    !error = error + values%set(key=mlbddc_coarse_matrix_symmetric_storage , Value= .false. )
        !    !error = error + values%set(key=mlbddc_coarse_matrix_is_symmetric      , Value= .false. )
        !    !error = error + values%set(key=mlbddc_coarse_matrix_sign              , Value= SPARSE_MATRIX_SIGN_UNKNOWN )
        !
        !    ! Output Handler
        !    !error = error + values%set(key=oh_staticgrid , Value= )
        !    !error = error + values%set(key=xh5_Strategy  , Value= )
        !    !error = error + values%set(key=xh5_Info      , Value= )
        !
        !    ! Sample
        !    !error = error + values%set(key=, Value= )

        !if associated ( this%define_user_parameters) then
        !  call this%define_user_parameters()
        !end if 

    end subroutine fph_define_fempar_parameters


    function fph_get_dir_path(this)
    !------------------------------------------------------------------
    !< Get dir_path
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in) :: this
        character(len=:),      allocatable             :: fph_get_dir_path
        type(ParameterList_t), pointer                 :: list
        integer(ip)                                    :: error
    !------------------------------------------------------------------
        list  => this%get_values()
        assert(list%isAssignable(dir_path_key, 'string'))
        error = list%GetAsString(key=dir_path_key, string = fph_get_dir_path)
        assert(error==0)
    end function fph_get_dir_path 


    function fph_get_dir_path_out(this)
    !------------------------------------------------------------------
    !< Get dir_path_out
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in) :: this
        character(len=:),      allocatable             :: fph_get_dir_path_out
        type(ParameterList_t), pointer                 :: list
        integer(ip)                                    :: error
    !------------------------------------------------------------------
        list  => this%get_values()
        assert(list%isAssignable(dir_path_out_key, 'string'))
        error = list%GetAsString(key=dir_path_out_key, string = fph_get_dir_path_out)
        assert(error==0)
    end function fph_get_dir_path_out


    function fph_get_prefix(this)
    !------------------------------------------------------------------
    !< Get prefix
    !------------------------------------------------------------------
        implicit none
        class(fempar_parameter_handler_t) , intent(in) :: this
        character(len=:),      allocatable             :: fph_get_prefix
        type(ParameterList_t), pointer                 :: list
        integer(ip)                                    :: error
    !------------------------------------------------------------------
        list  => this%get_values()
        assert(list%isAssignable(prefix_key, 'string'))
        error = list%GetAsString(key=prefix_key, string = fph_get_prefix)
        assert(error==0)
    end function fph_get_prefix

end module fempar_parameter_handler_names 
