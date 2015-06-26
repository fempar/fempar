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
program test_nsi_iss
  use fem_names
  use nsi_names
  use nsi_cg_iss_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"
  
  ! Types
  type(geom_data_t)                    :: gdata
  type(bound_data_t)                   :: bdata
  type(fem_fixed_info_t)               :: ginfo
  type(fem_triangulation_t)            :: f_trian
  type(fem_conditions_t)               :: f_cond
  type(dof_handler_t)                  :: dhand
  type(fe_space_t)                     :: fe_space  
  type(nsi_problem_t)                  :: myprob
  type(nsi_cg_iss_discrete_t) , target :: mydisc
  type(nsi_cg_iss_matvec_t)   , target :: matvec
  type(discrete_integration_pointer)   :: approx(1)
  type(fem_matrix_t)          , target :: femat
  type(fem_vector_t)          , target :: fevec,feunk
  type(fem_precond_t)         , target :: feprec
  type(fem_precond_params_t)           :: ppars
  type(solver_control_t)               :: sctrl
  type(serial_environment_t)           :: senv
  type(fem_vtk_t)                      :: fevtk
  class(base_operand_t)      , pointer :: x, b
  class(base_operator_t)     , pointer :: A, M
  type(fem_graph_t)          , pointer :: f_graph
  type(fem_block_graph_t)              :: f_blk_graph
  type(plain_vector_t)                 :: enorm
  type(nsi_cg_iss_error_t)    , target :: ecalc

  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: gtype(1) = (/ csr /)
  integer(ip) :: ibloc,jbloc,istat
  integer(ip) :: num_approximations = 1

  ! Allocatable
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: which_approx(:)

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez
  
  call meminit

  ! Read parameters from command-line
  call read_pars_cl_test_nsi(prefix,dir_path_out,nex,ney,nez)

  ! Generate geometry data
  call geom_data_create(gdata,nex,ney,nez)

  ! Generate boundary data
  call bound_data_create(gdata%ndime+1,gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp

  ! Generate element geometrical fixed info
  call finite_element_fixed_info_create(ginfo,Q_type_id,1,gdata%ndime,ginfo_state)

  ! Generate triangulation
  call gen_triangulation(1,gdata,bdata,ginfo,f_trian,f_cond,material)

  ! Create dof_handler
  call dhand%create(1,1,gdata%ndime+1)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call matvec%create(myprob,mydisc)
  call dhand%set_problem(1,mydisc)
  approx(1)%p       => matvec
  mydisc%dtinv      = 0.0_rp
  myprob%kfl_conv   = 1
  myprob%diffu      = 1.0_rp
  myprob%case_veloc = 1
  myprob%case_press = 1

  ! Allocate auxiliar elemental arrays
  call memalloc(f_trian%num_elems,dhand%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dhand%nvars_global,order,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  order(:,1:gdata%ndime) = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 
  
  ! Create fe_space
  call fe_space_create(f_trian,dhand,fe_space,problem,f_cond,continuity,order,material,which_approx, &
       &                time_steps_to_store=3, hierarchical_basis=.false.,             &
       &                static_condensation=.false.,num_continuity=1)

  ! Create plain vectors
  call fe_space_plain_vector_create((/2/),fe_space)
  call fe_space_plain_vector_point(2,fe_space)

  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,myprob,senv,dir_path_out,prefix,linear_order=.true.)

  ! Create dof info
  call create_dof_info(dhand,f_trian,fe_space,f_blk_graph,gtype)
  f_graph => f_blk_graph%get_block(1,1)

  ! Allocate matrices and vectors
  call fem_matrix_alloc(csr_mat,symm_false,f_graph,femat)
  call fem_vector_alloc(f_graph%nv,fevec)
  call fem_vector_alloc(f_graph%nv,feunk)
  call fevec%init(0.0_rp)

  ! Apply boundary conditions to unkno
  call fem_update_strong_dirichlet_bcond(fe_space,f_cond)
  call fem_update_analytical_bcond((/1:gdata%ndime/),myprob%case_veloc,0.0_rp,fe_space)
  call fem_update_analytical_bcond((/gdata%ndime+1/),myprob%case_press,0.0_rp,fe_space)

  ! Solver control parameters
  sctrl%method = direct
  sctrl%trace  = 100
  sctrl%track_conv_his = .true.

  ! Construct preconditioner
  ppars%type   = pardiso_mkl_prec
  call fem_precond_create(femat,feprec,ppars)
  call fem_precond_symbolic(femat,feprec)
  call fem_precond_log_info(feprec)

!!$  ! Define operators
!!$  A => femat
!!$  M => feprec
!!$  b => fevec
!!$  x => feunk

  ! Do nonlinear iterations
  call nonlinear_iteration(sctrl,1.0e-10_rp,10,senv,approx,fe_space,femat,feprec,fevec,feunk)
  
!!$  ! Integrate
!!$  call volume_integral(approx,fe_space,femat,fevec)
!!$
!!$  ! Solve
!!$  !call fem_vector_print(6,fevec)
!!$  call abstract_solve(femat,feprec,fevec,feunk,sctrl,senv)
!!$  !call abstract_solve(A,M,b,x,sctrl,senv)
!!$  call solver_control_log_conv_his(sctrl)
!!$  call solver_control_free_conv_his(sctrl)
!!$
!!$  ! Store solution to unkno
!!$  call fem_update_solution(feunk,fe_space)

  ! Free preconditioner
  call fem_precond_free(precond_free_struct,feprec)
  call fem_precond_free(precond_free_clean,feprec)

  ! Print solution to VTK file
  istat = fevtk%write_VTK()

  ! Compute error norm
  call enorm%create(2)
  call ecalc%create(myprob,mydisc)
  approx(1)%p => ecalc
  call volume_integral(approx,fe_space,enorm)
  write(*,*) 'Velocity error norm: ', sqrt(enorm%b(1))
  write(*,*) 'Pressure error norm: ', sqrt(enorm%b(2)) 

  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  call fevtk%free
  call f_blk_graph%free()
  call fem_vector_free(feunk)
  call fem_vector_free(fevec)
  call fem_matrix_free(femat) 
  call fe_space_free(fe_space) 
  call myprob%free
  call mydisc%free
  call matvec%free
  call ecalc%free
  call dof_handler_free(dhand)
  call fem_triangulation_free(f_trian)
  call fem_conditions_free(f_cond)
  call finite_element_fixed_info_free(ginfo)
  call bound_data_free(bdata)
  call enorm%free

  call memstatus

contains
  
  !==================================================================================================
  subroutine read_pars_cl_test_nsi(prefix,dir_path_out,nex,ney,nez)
    implicit none
    character*(*), intent(out) :: prefix, dir_path_out
    integer(ip)  , intent(out) :: nex,ney,nez
    character(len=256)         :: program_name
    character(len=256)         :: argument 
    integer                    :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==5) ) then
       write (6,*) 'Usage: ', trim(program_name), ' prefix dir_path_out nex ney nez'
       stop
    end if

    call getarg(1, argument)
    prefix = trim(argument)

    call getarg(2,argument)
    dir_path_out = trim(argument)

    call getarg(3, argument)
    read (argument,*) nex

    call getarg(4, argument)
    read (argument,*) ney

    call getarg(5, argument)
    read (argument,*) nez

  end subroutine read_pars_cl_test_nsi

  !==================================================================================================
  subroutine nonlinear_iteration( sctrl, nltol, maxit, env, approx, fe_space, A, M, b, x )
    implicit none
    type(solver_control_t)            , intent(inout) :: sctrl
    real(rp)                          , intent(in)    :: nltol
    integer(ip)                       , intent(in)    :: maxit    
    class(abstract_environment)       , intent(in)    :: env
    type(discrete_integration_pointer), intent(inout) :: approx(:)
    type(fe_space_t)                  , intent(inout) :: fe_space
    class(base_operator_t)            , intent(inout) :: A, M
    class(base_operand_t)             , intent(inout) :: x, b
    ! Locals
    integer(ip) :: iiter
    real(rp)    :: resnorm
    
    iiter = 0
    do while( iiter < maxit )

       ! Update counter
       iiter = iiter+1

       ! Initialize Matrix and vector
       ! ***************** Abstract procedure to initialize a base_operator ************************!
       select type (A)
       type is(fem_matrix_t)
          call fem_matrix_zero(A)
       class default
          check(.false.)
       end select
       !********************************************************************************************!
       call b%init(0.0_rp)

       ! Integrate system
       call volume_integral(approx,fe_space,A,b)

       ! Compute Numeric preconditioner
       ! ***************** Abstract procedure to compute precond numeric ***************************!
       select type (A)
       type is(fem_matrix_t)
          select type (M)
          type is(fem_precond_t)
             call fem_precond_numeric(A,M)
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
       !********************************************************************************************!

       ! Solve system
       call abstract_solve(A,M,b,x,sctrl,env)
       call solver_control_log_conv_his(sctrl)
       call solver_control_free_conv_his(sctrl)

       ! Free Numeric preconditioner
       ! ******************** Abstract procedure to free precond numeric ***************************!
       select type (M)
       type is(fem_precond_t)
          call fem_precond_free(precond_free_values,M)
          class default
          check(.false.)
       end select
       !********************************************************************************************!
       
       ! Store solution to unkno
       ! ***************** Abstract procedure to update from a base operant ************************!
       select type (x)
       type is(fem_vector_t)
          call fem_update_solution(x,fe_space)
       class default
          check(.false.)
       end select
       !********************************************************************************************!
       
       ! Store nonlinear iteration ( k+1 --> k )
       call fem_update_nonlinear(fe_space)

       ! Check convergence
       x = b - A*x
       resnorm = x%nrm2()
       write(*,*) 'Nonlinear error norm: ', resnorm
       if( resnorm < nltol) exit
       
    end do

  end subroutine nonlinear_iteration
  
end program test_nsi_iss
