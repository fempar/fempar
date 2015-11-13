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
module par_command_line_parameters_names
  use types_names
  use Data_Type_Command_Line_Interface
# include "debug.i90"

  implicit none
  private

  type par_test_cdr_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out

     character(len=:), allocatable :: default_kfl_conv 
     character(len=:), allocatable :: default_kfl_tder 
     character(len=:), allocatable :: default_kfl_react 
     character(len=:), allocatable :: default_react
     character(len=:), allocatable :: default_diffu 
     character(len=:), allocatable :: default_space_solution_flag 
     character(len=:), allocatable :: default_tempo_solution_flag 

     character(len=:), allocatable :: default_kfl_stab
     character(len=:), allocatable :: default_ftime
     character(len=:), allocatable :: default_itime
     character(len=:), allocatable :: default_theta
     character(len=:), allocatable :: default_time_step

     character(len=:), allocatable :: default_continuity
     character(len=:), allocatable :: default_enable_face_integration
     character(len=:), allocatable :: default_order

  end type par_test_cdr_params_t

  type par_test_cdr_parallel_params_t
     ! Input problem number of dimensions
     integer(ip)                   :: ndime

     ! Number of parts in which the problem was split
     integer(ip)                   :: nparts

     ! Graph Storage and Matrix properties
     character(len=:), allocatable :: default_symmetric_storage
     character(len=:), allocatable :: default_is_symmetric
     character(len=:), allocatable :: default_sign
     logical                       :: symmetric_storage(1)
     logical                       :: is_symmetric(1)
     integer(ip)                   :: sign(1)

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

   contains
     procedure :: set_par_default_params => par_test_cdr_parallel_params_set_par_default_params
  end type par_test_cdr_parallel_params_t

  ! Types
  public :: par_test_cdr_params_t, par_test_cdr_parallel_params_t

  ! Functions
  public :: set_par_default_params,cli_add_par_params,set_par_default_params_analytical
  public :: set_par_default_params_transient

contains

  subroutine set_par_default_params(params)
    implicit none
    type(par_test_cdr_params_t), intent(inout) :: params

    ! IO parameters
    params%default_dir_path     = 'data/'
    params%default_prefix       = 'square_4x4'
    params%default_dir_path_out = 'output'

    ! Problem parameters
    params%default_kfl_conv            = '0' ! Enabling advection
    params%default_kfl_tder            = '0' ! Time derivative not computed 
    params%default_kfl_react           = '0' ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '3'
    params%default_tempo_solution_flag = '0'

    ! Solver parameter
    params%default_kfl_stab = '0'   ! Stabilization of convective term (0: Off, 2: OSS)

    ! Time integration
    params%default_itime           = '0'
    params%default_ftime           = '0'
    params%default_theta           = '1.0'
    params%default_time_step       = '0.0'

    ! FE Space parameters
    params%default_continuity              = '1'
    params%default_enable_face_integration = '.false.'
    params%default_order                   = '1'
  end subroutine set_par_default_params

 !==================================================================================================
  subroutine par_test_cdr_parallel_params_set_par_default_params(params)
    use serial_names
    implicit none
    class(par_test_cdr_parallel_params_t), intent(inout) :: params

    !Input problem number of dimensions
     params%ndime                      = 2

     ! Number of parts in which the problem was split
     params%nparts                     = 1
    
    ! Graph Storage and Matrix properties
    params%default_symmetric_storage   = '.true.'
    params%default_is_symmetric        = '.true.'
    params%default_sign                = 'positive_definite'
    params%symmetric_storage           = (/.false./)
    params%is_symmetric                = (/.false./)
    params%sign                        = (/positive_definite/) 

    params%default_projection          = 'galerkin'
    params%default_nn_sys_sol_strat    = 'corners_rest_part_solve_expl_schur'
    params%default_pad_collectives     = '.false.'
    params%default_schur_edge_lag_mult = 'reuse_from_phis'
    params%default_subd_elmat_calc     = 'phit_minus_c_i_t_lambda'

  end subroutine par_test_cdr_parallel_params_set_par_default_params

  !==================================================================================================
  subroutine cli_add_par_params(cli,params,par_params,group)
    implicit none
    type(Type_Command_Line_Interface)    , intent(inout) :: cli
    type(par_test_cdr_params_t)          , intent(in)    :: params
    character(*)                         , intent(in)    :: group
    class(par_test_cdr_parallel_params_t), intent(inout) :: par_params
    ! Locals
    integer(ip) :: error

    ! Set Command Line Arguments
    ! IO parameters
    call cli%add(group=trim(group),switch='--dir_path',switch_ab='-d',                              &
         &       help='Directory of the source files',required=.false., act='store',                &
         &       def=trim(params%default_dir_path),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--prefix',switch_ab='-pr',help='Name of the GiD files',  &
         &       required=.false.,act='store',def=trim(params%default_prefix),error=error) 
    check(error==0)
    call cli%add(group=trim(group),switch='--dir_path_out',switch_ab='-out',help='Output Directory',&
         &       required=.false.,act='store',def=trim(params%default_dir_path_out),error=error)
    check(error==0)

    ! Problem parameters
    call cli%add(group=trim(group),switch='--kfl_conv',switch_ab='-kconv',help='Convection flag',   &
         &       required=.false.,act='store',def=trim(params%default_kfl_conv),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_react',switch_ab='-kreac',help='Reaction flag',    &
         &       required=.false.,act='store',def=trim(params%default_kfl_react),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_tder',switch_ab='-ktd',                            &
         &       help='Temporal derivative computation flag',required=.false.,act='store',          &
         &       def=trim(params%default_kfl_tder),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--react',switch_ab='-reac',help='Constant Reaction Value'&
         &       ,required=.false.,act='store',def=trim(params%default_react),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--diffu',switch_ab='-diff',                              &
         &       help='Constant Diffusion Value',required=.false.,act='store',                      &
         &       def=trim(params%default_diffu),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--space_solution_flag',switch_ab='-ssol',                &
         &       help='Space analytical solution',required=.false.,act='store',nargs='*',           &
         &       def=trim(params%default_space_solution_flag),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--tempo_solution_flag',switch_ab='-tsol',                &
         &       help='Temporal analytical solution',required=.false.,act='store',nargs='*',        &
         &       def=trim(params%default_tempo_solution_flag),error=error)
    check(error==0)

    ! Solver parameters
    call cli%add(group=trim(group),switch='--kfl_stab',switch_ab='-kst',help='Stabilization flag',  &
         &       required=.false.,act='store',def=trim(params%default_kfl_stab),error=error)
    check(error==0)

    ! Time integration parameters
    call cli%add(group=trim(group),switch='--itime',switch_ab='-t0',help='Initial time',            &
         &       required=.false.,act='store',def=trim(params%default_itime),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--ftime',switch_ab='-tf',help='Final time',              &
         &       required=.false.,act='store',def=trim(params%default_ftime),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--tstep',switch_ab='-ts',help='Time step',               &
         &       required=.false.,act='store',def=trim(params%default_time_step),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--theta',switch_ab='-tht',help='Theta method',           &
         &       required=.false.,act='store',def=trim(params%default_theta),error=error)
    check(error==0)

    ! FE Space parameters 
    call cli%add(group=trim(group),switch='--continuity',switch_ab='-cg',                           &
         &       help='Flag for the continuity of the FE Space',required=.false.,act='store',       &
         &       def=trim(params%default_continuity),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--face_integ',switch_ab='-fi',                           &
         &       help='Allow face integration',required=.false.,act='store',                        &
         &       def=trim(params%default_enable_face_integration),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--order',switch_ab='-p',                                 &
         &       help='Initial order of the approximation',required=.false.,act='store',            &
         &       def=trim(params%default_order),error=error)
    check(error==0)

    ! Solver
    call cli%add(group=trim(group),switch='--projection',switch_ab='-projection',                   &
         &       help='Coarse-grid projection strategy.',required=.false.,                          &
         &       def=par_params%default_projection,choices='galerkin,petrov_galerkin',              &
         &              act='store', error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--nn-sys-sol-strat',switch_ab='-nn-sys-sol-strat',                         &
         &       help='Strategy followed to solve the constrained Neumann problem',required=.false.,&
         &       def=par_params%default_nn_sys_sol_strat,                                           & 
         &       choices='corners_rest_part_solve_expl_schur,direct_solve_constrained_problem',     &
         &       act='store', error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--pad-collectives',switch_ab='-pad-collectives',         &
         &       help='Apply padding to inter-level collectives?',required=.false.,                 &
         &       def=par_params%default_pad_collectives,act='store', error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--schur-edge-lag-mult',switch_ab='-schur-edge-lag-mult', &
         &       help='Are the data resulting from computation of the edge Schur complement         &
         &       resulting from numerical set-up going to be re-used during application or computed &
         &       from scratch?',required=.false., def=par_params%default_schur_edge_lag_mult,       &
         &       choices='reuse_from_phis,compute_from_scratch', act='store', error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--subd-elmat-calc',switch_ab='-subd-elmat-calc',                           &
         &       help='Controls the way the subdomain coarse-grid matrix is computed.',             &
         &       required=.false., def=par_params%default_subd_elmat_calc,                          &
         &       choices='phit_minus_c_i_t_lambda,phit_a_i_phi', act='store', error=error)
    check(error==0)

  end subroutine cli_add_par_params

  !==================================================================================================
  subroutine set_par_default_params_analytical(params)
    implicit none
    type(par_test_cdr_params_t), intent(inout) :: params

    ! Names
    params%default_kfl_conv            = '0'    ! Enabling advection
    params%default_kfl_tder            = '0'    ! Time derivative not computed 
    params%default_kfl_react           = '0'    ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '4'
    params%default_tempo_solution_flag = '0'
    
  end subroutine set_par_default_params_analytical

  !==================================================================================================
  subroutine set_par_default_params_transient(params)
    implicit none
    type(par_test_cdr_params_t), intent(inout) :: params

    ! Problem
    params%default_kfl_conv            = '11'   ! Enabling advection
    params%default_kfl_tder            = '0'    ! Time derivative not computed 
    params%default_kfl_react           = '0'    ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '3'
    params%default_tempo_solution_flag = '1'
    params%default_itime               = '0'
    params%default_ftime               = '1'

    ! Time integration method
    params%default_theta           = '1.0'
    params%default_time_step       = '0.1'

  end subroutine set_par_default_params_transient

end module par_command_line_parameters_names

!****************************************************************************************************
!****************************************************************************************************

program par_test_cdr
  use serial_names
  use par_names
  use prob_names
  use lib_vtk_io_interface_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
  use time_integration_names
  use theta_method_names
  use par_command_line_parameters_names

  implicit none
#include "debug.i90"
  ! Our data
  type(par_context_t)                   :: w_context, p_context, q_context, b_context
  type(par_environment_t)               :: p_env
  type(par_mesh_t)                      :: p_mesh
  type(par_triangulation_t)             :: p_trian
  type(par_conditions_t)                :: p_cond
  type(dof_descriptor_t)                :: dof_descriptor
  type(par_fe_space_t)                  :: p_fe_space

  type(cdr_problem_t)                   :: my_problem
  type(cdr_discrete_t)                  :: my_discrete
  type(theta_method_t)                  :: theta_integ
  type(cdr_transient_t)        , target :: cdr_matvec
  type(error_norm_t)           , target :: compute_error
  type(par_scalar_t)                    :: enorm
  type(vtk_t)                           :: fevtk
  integer(ip)                           :: num_approximations
  type(p_discrete_integration_t)        :: approximations(1)
  class(matrix_t)             , pointer :: matrix
  type(par_scalar_matrix_t)   , pointer :: p_my_matrix
  class(array_t)              , pointer :: array
  type(par_scalar_array_t)    , pointer :: p_my_array
  type(par_scalar_array_t)    , target  :: p_feunk
  type(fe_affine_operator_t)            :: fe_affine_operator
  type(fe_affine_operator_t)            :: fe_affine_operator_error

  type(solver_control_t)                :: sctrl
  integer(ip)                           :: num_levels
  integer(ip), allocatable              :: id_parts(:), num_parts(:)

  type(par_preconditioner_dd_mlevel_bddc_t)       , target  :: p_mlevel_bddc
  type(par_preconditioner_dd_mlevel_bddc_params_t), target  :: p_mlevel_bddc_pars
  type(par_preconditioner_dd_mlevel_bddc_params_t), pointer :: point_to_p_mlevel_bddc_pars
  type(par_test_cdr_parallel_params_t)                      :: test_params
  integer(ip)                                 , allocatable :: kind_coarse_dofs(:)
 ! type(preconditioner_t)        :: feprec
  !type(preconditioner_params_t) :: ppars

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip)              :: i, j, vars_prob(1) = 1, ierror, iblock
  integer(ip)              :: space_solution_flag(1), tempo_solution_flag(1)
  integer(ip), allocatable :: order(:,:), material(:), problem(:)

  integer(ip), allocatable :: continuity(:,:)
  logical    , allocatable :: enable_face_integration(:,:)

  logical                  :: diagonal_blocks_symmetric_storage(1)
  logical                  :: diagonal_blocks_symmetric(1)
  integer(ip)              :: diagonal_blocks_sign(1)

  integer(ip) :: lunio, istat

  type(Type_Command_Line_Interface):: cli 
  character(len=:), allocatable :: group

  call meminit

  ! Start parallel execution
  call par_context_create (w_context)

  ! Read IO parameters
  call read_flap_cli_par_test_cdr(cli,test_params)
  call cli%parse(error=istat)
  if(cli%run_command('analytical')) then
     group = 'analytical'
  elseif(cli%run_command('transient')) then
     group = 'transient'
  else
     group = 'analytical'
  end if
  
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
  assert ( (p_context%iam >=0 .and. q_context%iam < 0) .or.                                         &
       &   (p_context%iam < 0 .and. q_context%iam >= 0))
  
  ! Create b_context as an intercommunicator among p_context <=> q_context 
  call par_context_create ( w_context, p_context, q_context, b_context )

  ! Create parallel environment
  call par_environment_create( p_env,w_context,p_context,q_context,b_context,num_levels,            &
       &                       id_parts,num_parts )

  ! Read mesh
  call cli%get(group=trim(group),switch='-d'  ,val=dir_path    ,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-pr' ,val=prefix      ,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-out',val=dir_path_out,error=istat); check(istat==0)
  call par_mesh_read (dir_path, prefix, p_env, p_mesh)

  ! Read conditions 
  write(*,*) __FILE__,__LINE__,trim(dir_path), trim(prefix), p_mesh%f_mesh%npoin
  call par_conditions_read (dir_path, prefix, p_mesh%f_mesh%npoin, p_env, p_cond)

  ! Construc triangulation
  call par_mesh_to_triangulation ( p_mesh, p_trian, p_cond )

  ! Assign the DOFs
  vars_prob = 1
  !( dof_descriptor, nblocks, nprobs, nvars_global, vars_block, dof_coupl )
  call dof_descriptor%create( 1, 1, 1 )

  ! Create the physical problem
  call physical_cdr_problem_create(my_problem,p_trian,cli,space_solution_flag,tempo_solution_flag,  &
       &                           group)

  ! Create the discrete method associated the problem
  call discrete_cdr_problem_create(my_discrete,my_problem,cli,group)

  ! Define the method to solve the proble
  call approximation_cdr_problem_create(approximations,my_problem,my_discrete,theta_integ,          &
       &                                cdr_matvec,cli,p_trian%f_trian,order,continuity,            &
       &                                enable_face_integration,                                    &
       &                                material,problem,group,dof_descriptor)

  ! Create FE Space
  call p_fe_space%create( p_trian, dof_descriptor, problem, p_cond, continuity,                     &
       &                  enable_face_integration, order, material, time_steps_to_store = 3,        &
       &                  hierarchical_basis = .false., static_condensation = .false.,              &
       &                  num_continuity = 1 )
  call par_create_distributed_dof_info ( p_fe_space )

  ! Initialize VTK output
  call fevtk%initialize(p_trian%f_trian,p_fe_space%serial_fe_space,my_problem,p_env,dir_path_out,   &
       &                prefix,linear_order=.false.)
  
  ! Assign analytical solution
  call p_fe_space%set_analytical_code(space_solution_flag,tempo_solution_flag)

  !Initialize time parameters
  call cdr_matvec%time_integ%initialize()

  ! Initialize solution
  if (theta_integ%dtinv .ne. 0.0_rp) then
     call par_update_analytical_initial(vars_prob,theta_integ%itime,p_fe_space)
     istat = fevtk%write_VTK(t_step = theta_integ%real_time)
     call theta_integ%update()
     ! Store the solution in the previous time step
     call par_update_transient_solution(p_fe_space,vars_prob,my_discrete%get_current(),             &
          &                          my_discrete%get_prev_step(),theta_integ)
  end if

  ! Create vef2dof array
  call par_create_distributed_dof_info( p_fe_space )

  ! Create the operator
  call fe_affine_operator%create (test_params%symmetric_storage , test_params%is_symmetric,         &
       test_params%sign, p_fe_space%serial_fe_space, approximations)
  call fe_affine_operator%symbolic_setup()
   
  ! Create the matrix
  matrix => fe_affine_operator%get_matrix()
  select type(matrix)
     class is(par_scalar_matrix_t)
     p_my_matrix => matrix
     class default
     check(.false.)
  end select

  ! Create array 
  array => fe_affine_operator%get_array()
  select type(array)
     class is(par_scalar_array_t)
     p_my_array => array
     class default
     check(.false.)
  end select

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
           point_to_p_mlevel_bddc_pars%correction_mode  = additive_symmetric
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
  ! call preconditioner_create(p_my_matrix,feprec,ppars)
 ! call preconditioner_symbolic(p_my_matrix,feprec)

  ! Create the computation of the error
  call compute_error%create(my_problem,my_discrete)
  call enorm%create(p_env)
  
  call memalloc ( test_params%ndime, kind_coarse_dofs, __FILE__, __LINE__ )
  kind_coarse_dofs(1) = corners
  kind_coarse_dofs(2) = corners_and_edges
  if ( test_params%ndime == 3 ) then
     kind_coarse_dofs(3) = corners_edges_and_faces
  end if

  sctrl%method  = cg
  sctrl%trace   = 1
  sctrl%itmax   = 800
  sctrl%dkrymax = 800
  sctrl%stopc   = res_res
  sctrl%orto    = icgs
  sctrl%rtol    = 1.0e-08

  do while (.not. theta_integ%finished) 

     ! Print the time step
     call theta_integ%print(6)

     ! Update boundary conditions
     call par_update_strong_dirichlet_bcond( p_fe_space, p_cond )
     call par_update_analytical_bcond(vars_prob,theta_integ%ctime,p_fe_space)

     ! Initialize the matrix and vector
     call p_my_matrix%init(0.0_rp)
     call p_my_array%init(0.0_rp) 
     ! Create the unknown vector
     call p_feunk%clone(p_my_array) 

     ! Compute the matrix an vector of the problem
     approximations(1)%discrete_integration => cdr_matvec
     call fe_affine_operator%numerical_setup()

     ! Preconditioner setting
     point_to_p_mlevel_bddc_pars => p_mlevel_bddc_pars
     do i=1, num_levels-1
        point_to_p_mlevel_bddc_pars%kind_coarse_dofs = kind_coarse_dofs(j)
        point_to_p_mlevel_bddc_pars => point_to_p_mlevel_bddc_pars%ppars_coarse_bddc
     end do
     ! Create multilevel bddc inverse 
     call par_preconditioner_dd_mlevel_bddc_create( p_my_matrix, p_mlevel_bddc, p_mlevel_bddc_pars )
     ! Ass struct
     call par_preconditioner_dd_mlevel_bddc_ass_struct ( p_my_matrix, p_mlevel_bddc )
     ! Fill val
     call par_preconditioner_dd_mlevel_bddc_fill_val ( p_mlevel_bddc )

     ! Update the preconditioner
     !call preconditioner_numeric(feprec)
     !call preconditioner_log_info(feprec)
     
     ! Solve the matrix-vector problem
     call abstract_solve(p_my_matrix,p_mlevel_bddc,p_my_array,p_feunk,sctrl,p_env)

     ! Store the solution in FE space
     call par_update_solution(p_feunk, p_fe_space)

     ! Store the solution in the previous time step
     call par_update_transient_solution(p_fe_space,vars_prob,my_discrete%get_current(),             &
          &                             my_discrete%get_prev_step(),theta_integ)

     ! Compute Error Norm
     approximations(1)%discrete_integration => compute_error
     compute_error%ctime = theta_integ%real_time
     call enorm%create(p_env)
     call enorm%init(0.0_rp)
     call p_fe_space%volume_integral(approximations,enorm)
     call enorm%reduce()
     write(*,*) 'XXX Error wrt analytical solution XXX',  sqrt(enorm%get_value())

     ! Free bddc inverse
     call  p_mlevel_bddc%free()!par_preconditioner_dd_mlevel_bddc_free( p_mlevel_bddc )
     call solver_control_free_conv_his(sctrl)

     ! Print solution to VTK file
     istat = fevtk%write_VTK(t_step = theta_integ%real_time)
     istat = fevtk%write_PVTK(t_step = theta_integ%real_time)
 
     ! Update the time integration variables
     call theta_integ%update()
  end do

   istat = fevtk%write_PVD()

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( enable_face_integration, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)

  call memfree(id_parts , __FILE__, __LINE__)
  call memfree(num_parts, __FILE__, __LINE__)

  call memfree ( kind_coarse_dofs, __FILE__, __LINE__ )
  call p_feunk%free()
  call fe_affine_operator%free()
  call p_fe_space%free()
  call my_problem%free
  call my_discrete%free
  call cdr_matvec%free
  call compute_error%free
  call compute_error%free
  call dof_descriptor_free ( dof_descriptor )
  call fevtk%free
  call par_triangulation_free ( p_trian )
  call par_conditions_free ( p_cond )
  call par_mesh_free (p_mesh)
  call memstatus
contains
  !==================================================================================================
  subroutine read_flap_cli_par_test_cdr(cli,test_params)
    use Data_Type_Command_Line_Interface
    use par_command_line_parameters_names
    use serial_names
    implicit none
    type(Type_Command_Line_Interface)   , intent(out)   :: cli
    type(par_test_cdr_parallel_params_t), intent(inout) :: test_params
    
    ! Locals
    type(par_test_cdr_params_t) :: analytical_params,transient_params
    logical     :: authors_print
    integer(ip) :: istat

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'par_test_cdr',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         &    'Serial FEMPAR driver to solve transient CDR problems using Continuous-Galerkin .',   &
         &        examples    = ['par_test_cdr -h         ', 'par_test_cdr analytical -h ' ])

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical solution')
    call cli%add_group(group='transient',description='Solve a problem with an transient solution')

    ! Set Parallel parameters
    call test_params%set_par_default_params()

    ! Set Command Line Arguments for each group
    call set_par_default_params(analytical_params)
    call set_par_default_params_analytical(analytical_params)
    call cli_add_par_params(cli,analytical_params,test_params,'analytical') 

    call set_par_default_params(transient_params)
    call set_par_default_params_transient(transient_params)
    call cli_add_par_params(cli,transient_params,test_params,'transient')

    ! Solver
    call cli%get(switch='-projection'         ,val=test_params%projection         ,error=istat)
    check(istat==0)
    call cli%get(switch='-nn-sys-sol-strat'   ,val=test_params%nn_sys_sol_strat   ,error=istat) 
    check(istat==0)
    call cli%get(switch='-pad-collectives'    ,val=test_params%pad_collectives    ,error=istat) 
    check(istat==0)
    call cli%get(switch='-subd-elmat-calc'    ,val=test_params%subd_elmat_calc    ,error=istat)
    check(istat==0)
    call cli%get(switch='-schur-edge-lag-mult',val=test_params%schur_edge_lag_mult,error=istat)
    check(istat==0)

  end subroutine read_flap_cli_par_test_cdr

  !==================================================================================================
  subroutine physical_cdr_problem_create(my_problem,p_trian,cli,space_solution_flag,                &
       &                                 tempo_solution_flag,group)
    use serial_names
    use par_names
    use prob_names
    use Data_Type_Command_Line_Interface
    implicit none
    type(cdr_problem_t)              , intent(inout) :: my_problem
    type(par_triangulation_t)        , intent(in)    :: p_trian
    type(Type_Command_Line_Interface), intent(inout) :: cli
    integer(ip)                      , intent(inout) :: space_solution_flag(1)
    integer(ip)                      , intent(inout) :: tempo_solution_flag(1)
    character(len=:),allocatable     , intent(in)    :: group

    integer(ip) :: istat
 
    !Create the problem to solve 
    call my_problem%create( p_trian%f_trian%num_dims )
    call cli%get(group=trim(group),switch='-kconv',val=my_problem%kfl_conv ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-kreac',val=my_problem%kfl_react,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ktd'  ,val=my_problem%kfl_tder ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-reac' ,val=my_problem%reaction ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-diff' ,val=my_problem%diffusion ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ssol' ,val=space_solution_flag ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tsol' ,val=tempo_solution_flag ,error=istat)
    check(istat==0)
  end subroutine physical_cdr_problem_create

  !==================================================================================================
  subroutine discrete_cdr_problem_create(my_discrete,my_problem,cli,group)
    use serial_names
    use prob_names
    use Data_Type_Command_Line_Interface
    implicit none
    type(cdr_discrete_t)             , intent(inout) :: my_discrete
    type(cdr_problem_t)              , intent(inout) :: my_problem
    type(Type_Command_Line_Interface), intent(inout) :: cli
    character(len=:),allocatable     , intent(in)    :: group

    real(rp)    :: time_step
    integer(ip) :: istat

    call my_discrete%create( my_problem)
    call cli%get(group=trim(group),switch='-kst' ,val=my_discrete%kfl_stab,error=istat)
    check(istat==0)
  end subroutine discrete_cdr_problem_create

  !==================================================================================================
  subroutine approximation_cdr_problem_create(approximations, my_problem,my_discrete,theta_integ,   &
       &                                      cdr_matvec,cli,f_trian,order,continuity,              &
       &                                      enable_face_integration,material, problem,group,      &
       &                                      dof_descriptor)
    use serial_names
    use prob_names
    use theta_method_names
    use Data_Type_Command_Line_Interface
    implicit none
    type(p_discrete_integration_t)   , intent(inout)         :: approximations(1)
    type(cdr_problem_t)              , intent(in)            :: my_problem
    type(cdr_discrete_t)             , intent(in)            :: my_discrete
    type(theta_method_t)             , intent(in)   , target :: theta_integ
    type(cdr_transient_t)            , intent(inout), target :: cdr_matvec
    type(Type_Command_Line_Interface), intent(inout)         :: cli
    type(triangulation_t)            , intent(in)            :: f_trian
    integer(ip)         , allocatable, intent(inout)         :: order(:,:),continuity(:,:)
    logical             , allocatable, intent(inout)         :: enable_face_integration(:,:)
    integer(ip)         , allocatable, intent(inout)         :: material(:),problem(:)
    character(len=:),allocatable     , intent(in)            :: group
    type(dof_descriptor_t)           , intent(inout)         :: dof_descriptor

    integer(ip) :: continuity_flag,face_int_flag, order_flag, istat, num_approximations
    real(rp)    :: time_step

    ! Create the solver type
    call cdr_matvec%create(my_problem,my_discrete)
    cdr_matvec%time_integ => theta_integ
    num_approximations = 1
    approximations(1)%discrete_integration => cdr_matvec
    call dof_descriptor%set_problem( 1, my_discrete )

    ! Continuity set up
    call memalloc( f_trian%num_elems, dof_descriptor%nvars_global, continuity, __FILE__, __LINE__)
    call cli%get(group=trim(group),switch='-cg', val=continuity_flag, error=istat); check(istat==0)
    continuity = continuity_flag

    ! Face integration set up
    call memalloc( f_trian%num_elems, dof_descriptor%nvars_global, enable_face_integration,         &
         &        __FILE__, __LINE__)
    call cli%get(group=trim(group),switch='-fi', val=face_int_flag, error=istat); check(istat==0)
    enable_face_integration = face_int_flag ! (cG/ No face integration)

    ! Order set up
    call memalloc( f_trian%num_elems, dof_descriptor%nvars_global, order, __FILE__, __LINE__)
    call cli%get(group=trim(group),switch='-p', val=order_flag, error=istat); check(istat==0)
    order = order_flag

    ! Define the material
    call memalloc( f_trian%num_elems, material, __FILE__, __LINE__)
    material = 1

    ! Define the problem
    call memalloc( f_trian%num_elems, problem, __FILE__, __LINE__)
    problem = 1


    call cli%get(group=trim(group),switch='-t0'  ,val=cdr_matvec%time_integ%itime     ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tf'  ,val=cdr_matvec%time_integ%ftime     ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ts'  ,val=cdr_matvec%time_integ%time_step ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tht' ,val=cdr_matvec%time_integ%theta     ,error=istat)
    check(istat==0)

  end subroutine approximation_cdr_problem_create


end program par_test_cdr
