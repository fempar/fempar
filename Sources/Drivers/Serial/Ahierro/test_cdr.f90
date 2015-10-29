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
  use types_names
  use Data_Type_Command_Line_Interface
# include "debug.i90"

  implicit none
  private

  type test_cdr_params_t
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

     character(len=:), allocatable :: default_kfl_thet 
     character(len=:), allocatable :: default_kfl_stab
     character(len=:), allocatable :: default_kfl_proj 
     character(len=:), allocatable :: default_k1tau   
     character(len=:), allocatable :: default_k2tau   
     character(len=:), allocatable :: default_dtinv 
     character(len=:), allocatable :: default_ctime 
     character(len=:), allocatable :: default_tdimv 
     character(len=:), allocatable :: default_nvars

     character(len=:), allocatable :: default_continuity
     character(len=:), allocatable :: default_enable_face_integration
     character(len=:), allocatable :: default_order
     
  end type test_cdr_params_t

  ! Types
  public :: test_cdr_params_t

  ! Functions
  public :: set_default_params,cli_add_params,set_default_params_analytical

contains

  subroutine set_default_params(params)
    implicit none
    type(test_cdr_params_t), intent(inout) :: params

    ! IO parameters
    params%default_dir_path     = 'data'
    params%default_prefix       = 'square_4x4'
    params%default_dir_path_out = 'output'

    ! Problem parameters
    params%default_kfl_conv            = '1' ! Enabling advection
    params%default_kfl_tder            = '0' ! Time derivative not computed 
    params%default_kfl_react           = '0'! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '3'
    params%default_tempo_solution_flag = '0'

    ! Solver parameter
    params%default_kfl_thet = '0'   ! Theta-method time integration (BE=0, CN=0)
    params%default_kfl_stab = '0'   ! Stabilization of convective term (0: Off, 2: OSS)
    params%default_kfl_proj = '0'   ! Projections weighted with tau's (On=1, Off=0)
    params%default_k1tau    = '4.0' ! C1 constant on stabilization parameter
    params%default_k2tau    = '2.0' ! C2 constant on stabilization parameter
    params%default_dtinv    = '1.0' ! Inverse of time step
    params%default_ctime    = '0.0' ! Current time
    params%default_tdimv    = '2'   ! Number of temporal steps stored
    params%default_nvars    = '1'

    ! FE Space parameters
    params%default_continuity              = '1'
    params%default_enable_face_integration = '.false.'
    params%default_order                   = '1'
  end subroutine set_default_params
  !==================================================================================================

  subroutine cli_add_params(cli,params,group)
    implicit none
    type(Type_Command_Line_Interface)            , intent(inout) :: cli
    type(test_cdr_params_t)                      , intent(in)    :: params
    character(*)                                 , intent(in)    :: group
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
    call cli%add(group=trim(group),switch='--kfl_thet',switch_ab='-ktht',help='Theta method',       &
         &       required=.false.,act='store',def=trim(params%default_kfl_thet),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_stab',switch_ab='-kst',help='Stabilization flag',  &
         &       required=.false.,act='store',def=trim(params%default_kfl_stab),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_proj',switch_ab='-kpr',                            &
         &       help='Projections weighted with tau flag',required=.false.,act='store',            &
         &       def=trim(params%default_kfl_proj),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--k1tau',switch_ab='-k1t',help='Tau 1',                  &
         &       required=.false.,act='store',def=trim(params%default_k1tau),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--k2tau',switch_ab='-k2t',help='Tau 2',                  &
         &       required=.false.,act='store',def=trim(params%default_k2tau),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--dtinv',switch_ab='-tinv',help='Inverse of time step',  &
         &       required=.false.,act='store',def=trim(params%default_dtinv),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--ctime',switch_ab='-ct',help='Current time',            &
         &       required=.false.,act='store',def=trim(params%default_ctime),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--tdimv',switch_ab='-tss',                               &
         &       help='Number of temporal steps stored',required=.false.,act='store',               &
         &       def=trim(params%default_tdimv),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--nvars',switch_ab='-nv',                                &
         &       help='Number of variables of the problem',required=.false.,act='store',            &
         &       def=trim(params%default_nvars),error=error)
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
  end subroutine cli_add_params
  !==================================================================================================

  subroutine set_default_params_analytical(params)
    implicit none
    type(test_cdr_params_t), intent(inout) :: params

    ! Names
    params%default_kfl_conv            = '1'    ! Enabling advection
    params%default_kfl_tder            = '0'    ! Time derivative not computed 
    params%default_kfl_react           = '0'    ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '4'
    params%default_tempo_solution_flag = '0'

  end subroutine set_default_params_analytical

end module command_line_parameters_names

!****************************************************************************************************
!****************************************************************************************************

program test_cdr
  use serial_names
  use prob_names
  use lib_vtk_io_interface_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                          :: f_mesh
  type(triangulation_t)                 :: f_trian
  type(conditions_t)                    :: f_cond
  type(dof_descriptor_t)                :: dof_descriptor
  type(serial_fe_space_t)               :: fe_space
  type(cdr_problem_t)                   :: my_problem
  type(cdr_discrete_t)                  :: my_discrete
  type(cdr_nonlinear_t), target         :: cdr_matvec
  type(error_norm_t)       , target     :: compute_error
  type(serial_scalar_t)                 :: enorm
  type(vtk_t)                           :: fevtk
  integer(ip)                           :: num_approximations
  type(p_discrete_integration_t)        :: approximations(1)
  class(matrix_t)             , pointer :: matrix
  type(serial_scalar_matrix_t), pointer :: my_matrix
  class(array_t)              , pointer :: array
  type(serial_scalar_array_t) , pointer :: my_array
  type(serial_scalar_array_t) , target  :: feunk
  type(fe_affine_operator_t)            :: fe_affine_operator
  type(fe_affine_operator_t)            :: fe_affine_operator_error


  class(vector_t) , pointer :: x, y
  class(operator_t), pointer :: A

  type(preconditioner_t)        :: feprec
  type(preconditioner_params_t) :: ppars
  type(solver_control_t)        :: sctrl
  type(serial_environment_t)    :: senv

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

  ! Read IO parameters
  call read_flap_cli_test_cdr(cli)
  call cli%parse(error=istat)
  if(cli%run_command('analytical')) then
     group = 'analytical'
  else
     group = 'analytical'
  end if
  call cli%get(group=trim(group),switch='-d',val=dir_path,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-pr',val=prefix,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-out',val=dir_path_out,error=istat); check(istat==0)

  ! Read mesh
  call mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)

  ! Construc triangulation
  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  ! Assign the DOFs
  vars_prob = 1
  !( dof_descriptor, nblocks, nprobs, nvars_global, vars_block, dof_coupl )
  call dof_descriptor%create( 1, 1, 1 )

  ! Create the physical problem
  call physical_cdr_problem_create(my_problem,cli,space_solution_flag,tempo_solution_flag)

  ! Create the discrete method associated the problem
  call discrete_cdr_problem_create(my_discrete,cli)

  ! Define the method to solve the proble
  call approximation_cdr_problem_create(my_problem,my_discrete,cdr_matvec,cli, f_trian,order,     &
       &                                continuity,enable_face_integration,material,problem)

  ! Create FE Space
  call fe_space%create ( f_trian, dof_descriptor, problem, f_cond, continuity,                      &
       &                 enable_face_integration, order, material, time_steps_to_store = 1,         &
       &                 hierarchical_basis = .false., static_condensation = .false.,               &
       &                 num_continuity = 1 )
  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,my_problem,senv,dir_path_out,prefix, linear_order=.false.)

  ! Assign analytical solution
  call fe_space%set_analytical_code(space_solution_flag,tempo_solution_flag)

  ! Update boundary conditions
  call update_strong_dirichlet_bcond( fe_space, f_cond )
  call update_analytical_bcond((/1/),my_discrete%ctime,fe_space)

  call create_dof_info( fe_space )

  ! Create the operator
  diagonal_blocks_symmetric_storage = (/.true./)
  diagonal_blocks_symmetric         = (/.true./)
  diagonal_blocks_sign              = (/positive_definite/) 
  call fe_affine_operator%create (diagonal_blocks_symmetric_storage , diagonal_blocks_symmetric,    &
       diagonal_blocks_sign, fe_space, approximations)
  call fe_affine_operator%symbolic_setup()
  call fe_affine_operator%numerical_setup()

  ! Create the matrix
  matrix => fe_affine_operator%get_matrix()
  select type(matrix)
     class is(serial_scalar_matrix_t)
     my_matrix => matrix
     class default
     check(.false.)
  end select

  array => fe_affine_operator%get_array()
  select type(array)
     class is(serial_scalar_array_t)
     my_array => array
     class default
     check(.false.)
  end select
  call feunk%clone(my_array) 

  write(*,*) 'XXX error vs exact norm XXX', feunk%nrm2()

  ! Create preconditioners
  ppars%type = pardiso_mkl_prec
  call preconditioner_create(my_matrix,feprec,ppars)
  call preconditioner_symbolic(my_matrix,feprec)
  call preconditioner_numeric(feprec)
  call preconditioner_log_info(feprec)

  call abstract_solve(my_matrix,feprec,my_array,feunk,sctrl,senv)

  ! Store the solution in FE space
  call update_solution(feunk, fe_space)

  ! Print solution to VTK file
  istat = fevtk%write_VTK()

  ! Compute Error Norm
  call compute_error%create(my_problem,my_discrete)
  approximations(1)%discrete_integration => compute_error
  call enorm%create()
  call enorm%init(0.0_rp)
  call fe_space%volume_integral(approximations,enorm)
  call enorm%reduce()

  call solver_control_free_conv_his(sctrl)
  write(*,*) __FILE__,__LINE__,fe_space%finite_elements(1)%integ(1)%p%uint_phy%nlocs
  A => my_matrix
  x => my_array
  y => feunk
  call feunk%print(6)
  call my_array%print(6)
  call my_matrix%print(6)

  write(*,*) 'XXX error solver norm XXX',  sqrt(enorm%get_value()) !feunk%nrm2()

  call memfree( continuity, __FILE__, __LINE__)
  call memfree( enable_face_integration, __FILE__, __LINE__)
  call memfree( order, __FILE__, __LINE__)
  call memfree( material, __FILE__, __LINE__)
  call memfree( problem, __FILE__, __LINE__)

  call feunk%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call my_problem%free
  call my_discrete%free
  call cdr_matvec%free
  call compute_error%free
  call dof_descriptor_free ( dof_descriptor )
  call fevtk%free
  call triangulation_free ( f_trian )
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)
  call memstatus
contains
  !==================================================================================================

  subroutine read_flap_cli_test_cdr(cli)
    implicit none
    type(Type_Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(test_cdr_params_t) :: analytical_params
    logical     :: authors_print
    integer(ip) :: error

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'test_cdr',                                                         &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         &    'Parallel FEMPAR driver to solve transient CDR problems using Continuous-Galerkin .', &
         &        examples    = ['test_cdr -h         ', 'test_cdr analytical -h ' ])

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical solution')

    ! Set Command Line Arguments for each group
    call set_default_params(analytical_params)
    call set_default_params_analytical(analytical_params)
    call cli_add_params(cli,analytical_params,'analytical')

  end subroutine read_flap_cli_test_cdr
  !==================================================================================================

  subroutine physical_cdr_problem_create(my_problem,cli,space_solution_flag,tempo_solution_flag)
    implicit none
    type(cdr_problem_t)              , intent(inout) :: my_problem
    type(Type_Command_Line_Interface), intent(inout) :: cli
    integer(ip)                      , intent(inout) :: space_solution_flag(1)
    integer(ip)                      , intent(inout) :: tempo_solution_flag(1)

    !Create the problem to solve 
    call my_problem%create( f_trian%num_dims )
    call cli%get(group=trim(group),switch='-kconv',val=my_problem%kfl_conv ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-kreac',val=my_problem%kfl_react,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ktd'  ,val=my_problem%kfl_tder ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-reac' ,val=my_problem%react    ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-diff' ,val=my_problem%diffu    ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ssol' ,val=space_solution_flag ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tsol' ,val=tempo_solution_flag ,error=istat)
    check(istat==0)
  end subroutine physical_cdr_problem_create
  !==================================================================================================

  subroutine discrete_cdr_problem_create(my_discrete,cli)
    implicit none
    type(cdr_discrete_t)             , intent(inout) :: my_discrete
    type(Type_Command_Line_Interface), intent(inout) :: cli

    call my_discrete%create( my_problem)
    call cli%get(group=trim(group),switch='-ktht',val=my_discrete%kfl_thet,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-kst' ,val=my_discrete%kfl_stab,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-kpr' ,val=my_discrete%kfl_proj,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-k1t' ,val=my_discrete%k1tau   ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-k2t' ,val=my_discrete%k2tau   ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tinv',val=my_discrete%dtinv   ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-ct'  ,val=my_discrete%ctime   ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-tss' ,val=my_discrete%tdimv   ,error=istat)
    check(istat==0)
    call cli%get(group=trim(group),switch='-nv'  ,val=my_discrete%nvars   ,error=istat)
    check(istat==0)
  end subroutine discrete_cdr_problem_create
  !==================================================================================================

  subroutine approximation_cdr_problem_create(my_problem,my_discrete,cdr_matvec,cli,          &
       &                                      f_trian,order,continuity,enable_face_integration,     &
       &                                      material, problem)
    implicit none
    type(cdr_problem_t)              , intent(in)            :: my_problem
    type(cdr_discrete_t)             , intent(in)            :: my_discrete
    type(cdr_nonlinear_t)            , intent(inout), target :: cdr_matvec
    type(Type_Command_Line_Interface), intent(inout)         :: cli
    type(triangulation_t)            , intent(in)            :: f_trian
    integer(ip)         , allocatable, intent(inout)         :: order(:,:),continuity(:,:)
    logical             , allocatable, intent(inout)         :: enable_face_integration(:,:)
    integer(ip)         , allocatable, intent(inout)         :: material(:),problem(:)

    integer(ip) :: continuity_flag,face_int_flag, order_flag

    ! Create the solver type
    call cdr_matvec%create(my_problem,my_discrete)
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
  end subroutine approximation_cdr_problem_create


end program test_cdr
