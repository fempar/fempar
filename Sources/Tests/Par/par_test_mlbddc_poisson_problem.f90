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


module command_line_parameters_names
  use serial_names
  use par_names
  use IR_Precision ! Integers and reals precision definition.
  use Data_Type_Command_Line_Interface, only: Type_Command_Line_Interface
  implicit none
#include "debug.i90" 
  private
  
  ! Types
  type par_test_mlbddc_poisson_problem_parameters_t
     ! Input problem location
     character(len=256) :: dir_path
     character(len=256) :: prefix
     
     ! Input problem number of dimensions
     integer(ip)                   :: ndime

     ! Number of parts in which the problem was split
     integer(ip)                   :: nparts

     ! Number of uniform refinement steps of input problem
     character(len=:), allocatable :: default_num_refinement_steps
     integer(ip)                   :: num_refinement_steps

     ! Graph Storage and Matrix properties
     character(len=:), allocatable :: default_symmetric_storage
     character(len=:), allocatable :: default_is_symmetric
     character(len=:), allocatable :: default_sign
     logical                       :: symmetric_storage
     logical                       :: is_symmetric
     integer(ip)                   :: sign

     ! BDDC parameters
     character(len=:), allocatable :: default_projection
     character(len=:), allocatable :: default_nn_sys_sol_strat
     character(len=:), allocatable :: default_pad_collectives
     character(len=:), allocatable :: default_schur_edge_lag_mult
     character(len=:), allocatable :: default_subd_elmat_calc

     integer(ip)                   :: projection
     integer(ip)                   :: nn_sys_sol_strat
     integer(ip)                   :: pad_collectives
     integer(ip)                   :: schur_edge_lag_mult
     integer(ip)                   :: subd_elmat_calc

     type(Type_Command_Line_Interface) :: cli

  contains
     procedure :: set_default_params => par_test_mlbddc_poisson_problem_parameters_set_default_params
     procedure :: set_cli            => par_test_mlbddc_poisson_problem_parameters_set_cli
     procedure :: parse              => par_test_mlbddc_poisson_problem_parameters_parse
  end type par_test_mlbddc_poisson_problem_parameters_t

  ! Types
  public :: par_test_mlbddc_poisson_problem_parameters_t

contains
  !==================================================================================================
  subroutine par_test_mlbddc_poisson_problem_parameters_set_default_params(params)
    implicit none
    class(par_test_mlbddc_poisson_problem_parameters_t), intent(inout) :: params

    params%dir_path = ''
    params%prefix = ''

    ! Number of uniform refinement steps of input problem
    params%default_num_refinement_steps = '3'

    ! Graph Storage and Matrix properties
    params%default_symmetric_storage   = '.true.'
    params%default_is_symmetric        = '.true.'
    params%default_sign                = 'positive_definite'
    params%default_projection          = 'galerkin'
    params%default_nn_sys_sol_strat    = 'corners_rest_part_solve_expl_schur'
    params%default_pad_collectives     = '.false.'
    params%default_schur_edge_lag_mult = 'reuse_from_phis'
    params%default_subd_elmat_calc     = 'phit_minus_c_i_t_lambda'

  end subroutine par_test_mlbddc_poisson_problem_parameters_set_default_params

  !==================================================================================================
  subroutine par_test_mlbddc_poisson_problem_parameters_set_cli(params)
    implicit none
    class(par_test_mlbddc_poisson_problem_parameters_t), intent(inout) :: params
    integer(I4P) :: error
    
    ! Initialize Command Line Interface
    call params%cli%init(progname    = 'par_test_mlbddc_poisson_problem', &
         &               version     = '',                               &
         &               authors     = '',                               &
         &               license     = '',                               &
         &               description = "FEMPAR parallel test driver to solve Poisson problem with analytical&
                                     & solution using Galerkin FEM, P1 FEs, and a 2-level PCG-BDDC solver. The&
                                     & input mesh might be uniformly refined if desired.", &
         &               examples    = ['par_test_mlbddc_poisson_problem -h'] )

    ! Set Command Line Arguments
    call params%cli%add(switch='--dir-path',switch_ab='-dir-path',help='Absolute or relative PATH to the partitioned&
                       & problem. Must end with /',required=.true., act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--prefix',switch_ab='-prefix',help='Prefix for all input files (mesh, conditions, etc.).& 
                       & E.g., if these files were generated from square.gid GiD project, then --prefix square.',& 
                       & required=.true., act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--ndime',switch_ab='-ndime',help='Number of dimensions of the problem at hand.',& 
                       & required=.true., act='store', choices='2,3', error=error)
    check(error==0)

    call params%cli%add(switch='--nparts',switch_ab='-nparts',help='Number of parts in which the problem was split.',& 
                       & required=.true., act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--num-refinement-steps',switch_ab='-num-refinement-steps',help='Number of uniform refinement steps.',& 
                       & required=.false., def=params%default_num_refinement_steps, act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--symmetric-storage',switch_ab='-symmetric-storage',help='Use symmetric storage for coefficient matrix.',& 
                       & required=.false., def=params%default_symmetric_storage, act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--is-symmetric',switch_ab='-is-symmetric',help='Is the coefficient matrix symmetric?',& 
                       & required=.false., def=params%default_is_symmetric, act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--sign',switch_ab='-sign',help='Sign of the matrix.',& 
                       & required=.false., def=params%default_sign, & 
                       & choices='positive_definite,indefinite,unknown', act='store', error=error)
    check(error==0)


    call params%cli%add(switch='--projection',switch_ab='-projection',help='Coarse-grid projection strategy.',& 
                       & required=.false., def=params%default_projection, & 
                       & choices='galerkin,petrov_galerkin', act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--nn-sys-sol-strat',switch_ab='-nn-sys-sol-strat',help='Strategy followed to solve the constrained Neumann problem.',& 
                       & required=.false., def=params%default_nn_sys_sol_strat, & 
                       & choices='corners_rest_part_solve_expl_schur,direct_solve_constrained_problem', act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--pad-collectives',switch_ab='-pad-collectives',help='Apply padding to inter-level collectives?',& 
                       & required=.false., def=params%default_pad_collectives, & 
                       & act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--schur-edge-lag-mult',switch_ab='-schur-edge-lag-mult',&
                        help='Are the data resulting from computation of the edge Schur complement&
                       & resulting from numerical set-up going to be re-used during application or computed from scratch?',& 
                       & required=.false., def=params%default_schur_edge_lag_mult, & 
                       & choices='reuse_from_phis,compute_from_scratch', act='store', error=error)
    check(error==0)

    call params%cli%add(switch='--subd-elmat-calc',switch_ab='-subd-elmat-calc',&
                        help='Controls the way the subdomain coarse-grid matrix is computed.',&
                       & required=.false., def=params%default_subd_elmat_calc, & 
                       & choices='phit_minus_c_i_t_lambda,phit_a_i_phi', act='store', error=error)
    check(error==0)
  end subroutine par_test_mlbddc_poisson_problem_parameters_set_cli

  !==================================================================================================
  subroutine par_test_mlbddc_poisson_problem_parameters_parse(params)
    implicit none
    class(par_test_mlbddc_poisson_problem_parameters_t), intent(inout) :: params
    integer(I4P)                                                       :: error
    character(len=256)                                                 :: aux_string
    logical                                                            :: aux_logical

    call params%cli%parse(error=error)
    check(error==0)

    call params%cli%get(switch='-dir-path',val=params%dir_path,error=error); check(error==0)
    call params%cli%get(switch='-prefix',val=params%prefix,error=error); check(error==0)
    call params%cli%get(switch='-ndime',val=params%ndime,error=error); check(error==0)
    call params%cli%get(switch='-nparts',val=params%nparts,error=error); check(error==0)
    call params%cli%get(switch='-num-refinement-steps',val=params%num_refinement_steps,error=error); check(error==0)
    call params%cli%get(switch='-symmetric-storage',val=params%symmetric_storage,error=error); check(error==0)
    call params%cli%get(switch='-is-symmetric',val=params%is_symmetric,error=error); check(error==0)
    aux_string(:)=' '
    call params%cli%get(switch='-sign',val=aux_string,error=error); check(error==0)
    if (trim(adjustl(aux_string)) .eq. 'positive_definite') then
       params%sign = positive_definite
    else if ( trim(adjustl(aux_string)) .eq. 'indefinite') then
       params%sign = indefinite
    else if ( trim(adjustl(aux_string)) .eq. 'unknown') then
       params%sign = unknown
    end if

    aux_string(:)=' '
    call params%cli%get(switch='-projection',val=aux_string,error=error); check(error==0)
    if (trim(adjustl(aux_string)) .eq. 'galerkin') then
       params%projection = galerkin
    else if ( trim(adjustl(aux_string)) .eq. 'petrov_galerkin') then
       params%projection = petrov_galerkin
    end if

    aux_string(:)=' '
    call params%cli%get(switch='-nn-sys-sol-strat',val=aux_string,error=error); check(error==0)
    if (trim(adjustl(aux_string)) .eq. 'corners_rest_part_solve_expl_schur') then
       params%nn_sys_sol_strat = corners_rest_part_solve_expl_schur
    else if ( trim(adjustl(aux_string)) .eq. 'direct_solve_constrained_problem') then
       params%nn_sys_sol_strat = direct_solve_constrained_problem
    end if

    call params%cli%get(switch='-pad-collectives',val=aux_logical,error=error); check(error==0)
    if (aux_logical) then
       params%pad_collectives = pad
    else 
       params%pad_collectives = nopad
    end if

    aux_string(:)=' '
    call params%cli%get(switch='-schur-edge-lag-mult',val=aux_string,error=error); check(error==0)
    if (trim(adjustl(aux_string)) .eq. 'reuse_from_phis') then
       params%schur_edge_lag_mult = reuse_from_phis 
    else if ( trim(adjustl(aux_string)) .eq. 'compute_from_scratch') then
       params%schur_edge_lag_mult = compute_from_scratch
    end if

    aux_string(:)=' '
    call params%cli%get(switch='-subd-elmat-calc',val=aux_string,error=error); check(error==0)
    if (trim(adjustl(aux_string)) .eq. 'phit_minus_c_i_t_lambda') then
       params%subd_elmat_calc = phit_minus_c_i_t_lambda
    else if ( trim(adjustl(aux_string)) .eq. 'phit_a_i_phi') then
       params%subd_elmat_calc = phit_a_i_phi
    end if

  end subroutine par_test_mlbddc_poisson_problem_parameters_parse

end module command_line_parameters_names



program par_test_mlbddc_poisson_problem
  !----------------------------------------------------------
  ! Parallel partitioner test
  !----------------------------------------------------------
  use serial_names
  use par_names
  use prob_names
  use command_line_parameters_names

  implicit none
#include "debug.i90" 

  ! Our data
  type(par_context_t)       :: w_context, p_context, q_context, b_context
  type(par_environment_t)   :: p_env
  type(par_mesh_t)          :: p_mesh
  type(par_triangulation_t) :: p_trian
  type(par_fe_space_t)      :: p_fe_space

  type(par_scalar_matrix_t), target                :: p_mat
  type(par_scalar_array_t), target                :: p_vec, p_unk
  class(vector_t) , pointer           :: x, y
  class(operator_t), pointer           :: A

  ! Preconditioner-related data structures
  type(par_preconditioner_dd_diagonal_t)           :: p_prec_dd_diag
  type(par_preconditioner_dd_mlevel_bddc_t), target :: p_mlevel_bddc
  type(par_preconditioner_dd_mlevel_bddc_params_t), target  :: p_mlevel_bddc_pars
  type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars
  integer(ip), allocatable :: kind_coarse_dofs(:)
  type(solver_control_t)            :: sctrl
  type(blocks_dof_distribution_t), pointer    :: blocks_dof_distribution
  type(dof_descriptor_t)            :: dof_descriptor
  type(par_conditions_t)            :: p_cond

  type(cdr_problem_t)               :: my_problem
  type(cdr_discrete_t)              :: my_discrete
  type(cdr_nonlinear_t), target     :: my_approximation
  type(par_scalar_t)                :: enorm
  type(error_norm_t)   , target     :: error_compute
  type(discrete_integration_pointer_t)  :: approximations(1)

  integer(ip)              :: num_levels
  integer(ip), allocatable :: id_parts(:), num_parts(:)

  ! Arguments
  integer(ip)                   :: lunio, lunou, istat
  character(len=:), allocatable :: name
  integer(ip)                   :: i, j, ierror, iblock
  type(par_timer_t)             :: par_uniform_refinement_timer, par_mesh_to_triangulation_timer, par_fe_space_create_timer

  integer(ip), allocatable :: order(:,:), material(:), problem(:)
  integer(ip), allocatable :: continuity(:,:)
  logical    , allocatable :: enable_face_integration(:,:)


  type(par_test_mlbddc_poisson_problem_parameters_t) :: test_params

  call meminit

  ! Start parallel execution
  call par_context_create (w_context)

#ifdef ENABLE_MKL
  call test_params%set_default_params()
  call test_params%set_cli()
  call test_params%parse()

  ! This test only works with two levels.
  num_levels = 2
  call memalloc(num_levels, id_parts , __FILE__, __LINE__)
  call memalloc(num_levels, num_parts, __FILE__, __LINE__)

  num_parts = (/test_params%nparts, 1/)
  id_parts = (/w_context%iam+1, 1/)


  ! Create p_context and q_context splitting w_context
  if(w_context%iam < num_parts(1)) then
     call par_context_create ( 1, p_context, q_context, w_context )
  else
     call par_context_create ( 2, q_context, p_context, w_context )
  end if
  assert ( (p_context%iam >=0 .and. q_context%iam < 0) .or. (p_context%iam < 0 .and. q_context%iam >= 0))
  
  ! Create b_context as an intercommunicator among p_context <=> q_context 
  call par_context_create ( w_context, p_context, q_context, b_context )

  ! Create parallel environment
  call par_environment_create( p_env,& 
                               w_context,& 
                               p_context,& 
                               q_context,&
                               b_context,&
                               num_levels,&
                               id_parts, & 
                               num_parts )

  call par_timer_create ( par_mesh_to_triangulation_timer, 'PAR_MESH_TO_TRIANGULATION', w_context%icontxt )
  call par_timer_create ( par_fe_space_create_timer, 'PAR_FE_SPACE_CREATE', w_context%icontxt )
  call par_timer_create ( par_uniform_refinement_timer, 'PAR_UNIFORM_REFINEMENT', w_context%icontxt )

  ! Read mesh
  call par_mesh_read ( test_params%dir_path, test_params%prefix, p_env, p_mesh )

  ! Read boundary conditions
  call par_conditions_read(test_params%dir_path, test_params%prefix, p_mesh%f_mesh%npoin, p_env, p_cond)

  do i=1, test_params%num_refinement_steps
     call par_timer_init (par_mesh_to_triangulation_timer)
     call par_timer_start (par_mesh_to_triangulation_timer)   
     call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)
     call par_timer_stop (par_mesh_to_triangulation_timer)   
     call par_timer_report (par_mesh_to_triangulation_timer)   

     call par_mesh_free(p_mesh)

     call par_timer_init (par_uniform_refinement_timer)
     call par_timer_start (par_uniform_refinement_timer) 
     call par_uniform_refinement ( p_trian, p_mesh, p_cond )
     call par_timer_stop (par_uniform_refinement_timer)  
     call par_timer_report (par_uniform_refinement_timer)

     call par_triangulation_free(p_trian)
  end do
  call par_timer_init (par_mesh_to_triangulation_timer)
  call par_timer_start (par_mesh_to_triangulation_timer)   
  call par_mesh_to_triangulation (p_mesh, p_trian, p_cond)
  call par_timer_stop (par_mesh_to_triangulation_timer)   
  call par_timer_report (par_mesh_to_triangulation_timer)   


  call dof_descriptor%create( 1, 1, 1 )

  call my_problem%create( p_trian%f_trian%num_dims )
  call my_discrete%create( my_problem )
  call my_approximation%create(my_problem,my_discrete)
  approximations(1)%p => my_approximation
  
  call dof_descriptor%set_problem( 1, my_discrete )
  ! ... for as many problems as we have

  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, continuity, __FILE__, __LINE__)
  continuity = 1
  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, enable_face_integration, __FILE__, __LINE__)
  enable_face_integration = .false.
  call memalloc( p_trian%f_trian%num_elems, dof_descriptor%nvars_global, order, __FILE__, __LINE__)
  order = 1
  call memalloc( p_trian%f_trian%num_elems, material, __FILE__, __LINE__)
  material = 1
  call memalloc( p_trian%f_trian%num_elems, problem, __FILE__, __LINE__)
  problem = 1

  call par_timer_start (par_fe_space_create_timer)

  call par_fe_space_create ( p_trian, dof_descriptor, p_fe_space, problem, &
                             p_cond, continuity, enable_face_integration, order, material, &
                             time_steps_to_store = 1, &
                             hierarchical_basis = .false., &
                             & static_condensation = .false., num_continuity = 1 )

  call par_timer_stop (par_fe_space_create_timer)
  call par_timer_report(par_fe_space_create_timer)

  call p_fe_space%set_analytical_code( spatial_code=(/1/), temporal_code=(/0/) ) 

  call par_create_distributed_dof_info ( p_fe_space )
  blocks_dof_distribution => p_fe_space%get_blocks_dof_distribution()
  
  call p_fe_space%make_coefficient_matrix ( test_params%symmetric_storage, & 
										    test_params%is_symmetric, &
											test_params%sign, &
											p_mat )

  call p_vec%create_and_allocate ( blocks_dof_distribution%get_block(1), p_env )
  p_vec%state = part_summed

  call p_unk%create_and_allocate ( blocks_dof_distribution%get_block(1), p_env )
  p_unk%state = full_summed

  call par_update_analytical_bcond( vars_of_unk=(/1/), ctime=0.0_rp, p_fe_space=p_fe_space)

  call par_volume_integral( approximations, p_fe_space, p_mat, p_vec )
 
  ! Define (recursive) parameters
  point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
  do i=1, num_levels-1
     point_to_p_mlevel_bddc_pars%ndime                = test_params%ndime
     point_to_p_mlevel_bddc_pars%unknowns             = all_unknowns
     point_to_p_mlevel_bddc_pars%pad_collectives      = test_params%pad_collectives
     point_to_p_mlevel_bddc_pars%projection           = test_params%projection
     point_to_p_mlevel_bddc_pars%schur_edge_lag_mult  = test_params%schur_edge_lag_mult
     point_to_p_mlevel_bddc_pars%subd_elmat_calc      = test_params%subd_elmat_calc             
     point_to_p_mlevel_bddc_pars%correction_mode      = additive_symmetric                  
     point_to_p_mlevel_bddc_pars%nn_sys_sol_strat     = test_params%nn_sys_sol_strat  

     if ( i < num_levels-1 ) then
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat = recursive_bddc
        point_to_p_mlevel_bddc_pars%ppars_harm%type      = pardiso_mkl_prec  
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type = pardiso_mkl_prec 
        
        if ( i == 1 ) then
           point_to_p_mlevel_bddc_pars%spars_coarse%method = direct
           point_to_p_mlevel_bddc_pars%spars_coarse%itmax  = 200
           point_to_p_mlevel_bddc_pars%spars_coarse%rtol   = 1.0e-08
           point_to_p_mlevel_bddc_pars%spars_coarse%trace  = 1
           point_to_p_mlevel_bddc_pars%correction_mode  = additive
        end if
        allocate(point_to_p_mlevel_bddc_pars%ppars_coarse_bddc, stat = ierror)
        check(ierror==0)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     else
        point_to_p_mlevel_bddc_pars%co_sys_sol_strat = serial_gather
        point_to_p_mlevel_bddc_pars%ppars_harm%type          =pardiso_mkl_prec  
        point_to_p_mlevel_bddc_pars%ppars_dirichlet%type     =pardiso_mkl_prec   
        point_to_p_mlevel_bddc_pars%ppars_coarse_serial%type =pardiso_mkl_prec  
        nullify ( point_to_p_mlevel_bddc_pars%ppars_coarse_bddc )
     end if
  end do


  ! Compute error norm
  call error_compute%create(my_problem,my_discrete)
  approximations(1)%p => error_compute
  call enorm%create(p_env)

  call memalloc ( test_params%ndime, kind_coarse_dofs, __FILE__, __LINE__ )
  kind_coarse_dofs(1) = corners
  kind_coarse_dofs(2) = corners_and_edges
  if ( test_params%ndime == 3 ) then
     kind_coarse_dofs(3) = corners_edges_and_faces
  end if

  sctrl%method=cg
  sctrl%trace=1
  sctrl%itmax=800
  sctrl%dkrymax=800
  sctrl%stopc=res_res
  sctrl%orto=icgs
  sctrl%rtol=1.0e-08

  do j=1,test_params%ndime
     point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
     do i=1, num_levels-1
        point_to_p_mlevel_bddc_pars%kind_coarse_dofs = kind_coarse_dofs(j)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     end do
     call p_unk%init(0.0_rp)
     ! Create multilevel bddc inverse 
     call par_preconditioner_dd_mlevel_bddc_create( p_mat, p_mlevel_bddc, p_mlevel_bddc_pars )
     ! Ass struct
     call par_preconditioner_dd_mlevel_bddc_ass_struct ( p_mat, p_mlevel_bddc )
     ! Fill val
     call par_preconditioner_dd_mlevel_bddc_fill_val ( p_mlevel_bddc )

     call abstract_solve(p_mat,p_mlevel_bddc,p_vec,p_unk,sctrl,p_env)

     call par_update_solution(p_unk,p_fe_space)

     call enorm%init()
     call par_volume_integral(approximations,p_fe_space,enorm)
     call enorm%reduce()
     if(w_context%iam==0) then
        check ( sqrt(enorm%get()) < 1.0e-06 )
     end if

     ! Free bddc inverse
     call par_preconditioner_dd_mlevel_bddc_free( p_mlevel_bddc, free_values)
     call par_preconditioner_dd_mlevel_bddc_free( p_mlevel_bddc, free_struct)
     call par_preconditioner_dd_mlevel_bddc_free( p_mlevel_bddc, free_clean)
  end do

  call memfree ( kind_coarse_dofs, __FILE__, __LINE__ )

  call p_mat%free()
  call p_vec%free()
  call p_unk%free()

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( enable_face_integration, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)

  call par_fe_space_free(p_fe_space) 
  call my_problem%free
  call my_discrete%free
  call my_approximation%free
  call error_compute%free
  call dof_descriptor_free (dof_descriptor)
  call par_triangulation_free(p_trian)
  call par_conditions_free (p_cond)
  call par_mesh_free (p_mesh)

  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)

  call par_environment_free (p_env)
  call par_context_free ( b_context, .false. )
  call par_context_free ( p_context, .false. )
  call par_context_free ( q_context, .false. )
#endif

  call par_context_free ( w_context )

  call memstatus

end program par_test_mlbddc_poisson_problem

