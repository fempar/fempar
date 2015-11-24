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
module postprocess_field_nsi_names
  use serial_names
  use prob_names
# include "debug.i90"
  implicit none
  private

  type, extends(postprocess_field_t) :: postprocess_field_velocity_t
   contains
     procedure :: compute_field => nsi_compute_field_vel
  end type postprocess_field_velocity_t

  type, extends(postprocess_field_t) :: postprocess_field_pressure_t
   contains
     procedure :: compute_field => nsi_compute_field_pre
  end type postprocess_field_pressure_t

  public :: postprocess_field_velocity_t, postprocess_field_pressure_t
contains
  subroutine nsi_compute_field_vel(postprocess_field)
    !******** READ ME: *********************************************************!
    ! This subroutine nsi_compute_field has been created only for testing       !
    ! purposes and it is left here as an easy example of how to use the         !
    ! postprocess_fields_names module.                                          !
    !---------------------------------------------------------------------------!
    implicit none
    class(postprocess_field_velocity_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)        :: nelem, ndime, i

    ! Extract some parameter form fe_space
    nelem = postprocess_field%fe_space%g_trian%num_elems
    ndime = postprocess_field%fe_space%g_trian%num_dims

    do i=1,nelem
       postprocess_field%fe_postprocess_field(i)%nodal_properties =                              &
            &          postprocess_field%fe_space%finite_elements(i)%unkno(:,1:ndime,1)
    end do
    

  end subroutine nsi_compute_field_vel

  subroutine nsi_compute_field_pre(postprocess_field)
    !******** READ ME: *********************************************************!
    ! This subroutine nsi_compute_field has been created only for testing       !
    ! purposes and it is left here as an easy example of how to use the         !
    ! postprocess_fields_names module.                                          !
    !---------------------------------------------------------------------------!
    implicit none
    class(postprocess_field_pressure_t), intent(inout) :: postprocess_field
    ! Locals
    integer(ip)        :: nelem, ndime, i

    ! Extract some parameter form fe_space
    nelem = postprocess_field%fe_space%g_trian%num_elems
    ndime = postprocess_field%fe_space%g_trian%num_dims

    do i=1,nelem
       postprocess_field%fe_postprocess_field(i)%nodal_properties =                              &
            &          postprocess_field%fe_space%finite_elements(i)%unkno(:,1+ndime:1+ndime,1)
    end do
    
  end subroutine nsi_compute_field_pre
end module postprocess_field_nsi_names

program test_nsi_iss
  use serial_names
  use prob_names
  use postprocess_field_nsi_names
  implicit none
# include "debug.i90"
  
  ! Types
  type(uniform_mesh_descriptor_t)       :: gdata
  type(uniform_conditions_descriptor_t) :: bdata
  type(reference_element_t)             :: geo_reference_element
  type(triangulation_t)                 :: f_trian
  type(conditions_t)                    :: f_cond
  type(dof_descriptor_t)                :: dof_descriptor
  type(serial_fe_space_t)               :: fe_space  
  type(nsi_problem_t)                   :: myprob

  type(nsi_cg_iss_discrete_t) , target  :: mydisc
  type(nsi_cg_iss_matvec_t)   , target  :: cg_iss_matvec
  type(p_discrete_integration_t)        :: approx(1)
  class(matrix_t)             , pointer :: matrix
  type(serial_scalar_matrix_t), pointer :: femat
  class(array_t)              , pointer :: array
  type(serial_scalar_array_t) , pointer :: fevec
  type(serial_scalar_array_t)           :: feunk
  type(fe_affine_operator_t)            :: fe_affine_operator

  
  type(preconditioner_t)                :: feprec
  type(preconditioner_params_t)         :: ppars
  type(solver_control_t)                :: sctrl
  type(serial_environment_t)            :: senv
  type(vtk_t)                           :: fevtk
  type(serial_scalar_t)                 :: enorm_u, enorm_p
  type(error_norm_t)           , target :: error_compute
  type(postprocess_field_velocity_t)    :: postprocess_vel
  type(postprocess_field_pressure_t)    :: postprocess_pre
  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: ibloc,jbloc,istat,i

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2

  ! Allocatable
  integer(ip), allocatable :: continuity(:,:)
  logical    , allocatable :: enable_face_integration(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez
  
  call meminit

  ! Read parameters from command-line
  call read_pars_cl_test_nsi(prefix,dir_path_out,nex,ney,nez)

  ! Generate geometry data
  call uniform_mesh_descriptor_create(gdata,nex,ney,nez)

  ! Generate boundary data
  call uniform_conditions_descriptor_create(gdata%ndime+1,gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp

  ! Generate element geometrical fixed info
  call reference_element_create(geo_reference_element,Q_type_id,1,gdata%ndime)

  ! Generate triangulation
  call generate_uniform_triangulation(1,gdata,bdata,geo_reference_element,f_trian,f_cond,material)

  ! Create dof_descriptor
  call dof_descriptor%create(1,1,gdata%ndime+1)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call cg_iss_matvec%create(myprob,mydisc)
  call dof_descriptor%set_problem(1,mydisc)
  approx(1)%discrete_integration       => cg_iss_matvec
  mydisc%dtinv      = 0.0_rp
  myprob%kfl_conv   = 1
  myprob%diffu      = 1.0_rp

  ! Allocate auxiliar elemental arrays
  call memalloc(f_trian%num_elems,dof_descriptor%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dof_descriptor%nvars_global,enable_face_integration, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dof_descriptor%nvars_global,order,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,problem,__FILE__,__LINE__)
  continuity             = 1
  enable_face_integration = .false.
  order(:,1:gdata%ndime) = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  
  ! Create fe_space
  call fe_space%create(f_trian,dof_descriptor,problem,f_cond,continuity,enable_face_integration, &
       &                order,material,time_steps_to_store=3, hierarchical_basis=.false.,    &
       &                static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,myprob,senv,dir_path_out,prefix)!,linear_order=.true.)

  ! Create dof info
  call create_dof_info(fe_space)

  ! Assign analytical solution
  if(gdata%ndime==2) then
     call fe_space%set_analytical_code((/4,5,3/),(/0,0,0/))
  else
     write(*,*) 'analytical function not ready for 3D'
  end if
  
  ! Apply boundary conditions to unkno
  call update_strong_dirichlet_bcond(fe_space,f_cond)
  ! CORRECT
  call update_analytical_bcond(vars_of_unk=(/(i,i=1,gdata%ndime)/),ctime=0.0_rp,fe_space=fe_space)
  call update_analytical_bcond(vars_of_unk=(/gdata%ndime+1/),ctime=0.0_rp,fe_space=fe_space)
  
  call fe_affine_operator%create ( (/.false./), &
								   (/.false./), & 
								   (/indefinite/), &
								   fe_space, &
								   approx)
  
  matrix => fe_affine_operator%get_matrix()
  select type(matrix)
  class is(serial_scalar_matrix_t)
    femat => matrix
  class default
    check(.false.)
  end select 
  
  array => fe_affine_operator%get_array()
  select type(array)
  class is(serial_scalar_array_t)
    fevec => array
  class default
    check(.false.)
  end select 
  call feunk%clone(fevec) 

  ! Solver control parameters
  sctrl%method = direct
  sctrl%trace  = 100
  sctrl%track_conv_his = .true.

  ! Construct preconditioner
  ppars%type   = pardiso_mkl_prec
  call preconditioner_create(femat,feprec,ppars)
  call preconditioner_symbolic_setup(feprec)
  call preconditioner_log_info(feprec)

  ! Do nonlinear iterations
  !call nonlinear_iteration(sctrl,1.0e-10_rp,10,senv,approx,fe_space,femat,feprec,fevec,feunk)

  ! Free preconditioner
  call preconditioner_free_in_stages(feprec, free_symbolic_setup)
  call preconditioner_free_in_stages(feprec, free_clean)

  ! Compute postprocess field
  call postprocess_vel%create('velocity',gdata%ndime,fe_space,senv,                                 &
       &                      use_interpolation_order_from_variable,variable_identifier=1)
  call postprocess_pre%create('pressure',1,fe_space,senv,                                           &
       &                      use_interpolation_order_from_variable,variable_identifier=2)
  call postprocess_vel%compute_and_finalize_field()
  call postprocess_pre%compute_and_finalize_field() 

  ! Print solution to VTK file
  istat = fevtk%write_VTK_start()
  istat = fevtk%write_VTK_unknowns()
  istat = fevtk%write_VTK_field(postprocess_vel)
  istat = fevtk%write_VTK_field(postprocess_pre)
  istat = fevtk%write_VTK_end()

  call postprocess_vel%free
  call postprocess_pre%free

  ! Compute error norm
  call error_compute%create(myprob,mydisc)
  approx(1)%discrete_integration => error_compute
  error_compute%unknown_id = velocity
  call enorm_u%init(0.0_rp)
  call fe_space%volume_integral(approx,enorm_u)
  error_compute%unknown_id = pressure
  call enorm_p%init(0.0_rp)
  call fe_space%volume_integral(approx,enorm_p)
  write(*,*) 'Velocity error norm: ', sqrt(enorm_u%get_value())
  write(*,*) 'Pressure error norm: ', sqrt(enorm_p%get_value()) 

  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(enable_face_integration,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call fevtk%free
  call feunk%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call myprob%free
  call mydisc%free
  call error_compute%free
  call cg_iss_matvec%free
  call dof_descriptor_free(dof_descriptor)
  call triangulation_free(f_trian)
  call conditions_free(f_cond)
  call reference_element_free(geo_reference_element)
  call uniform_conditions_descriptor_free(bdata)

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

  !!==================================================================================================
  !subroutine nonlinear_iteration( sctrl, nltol, maxit, env, approx, fe_space, A, M, b, x )
  !  implicit none
  !  type(solver_control_t)              , intent(inout) :: sctrl
  !  real(rp)                            , intent(in)    :: nltol
  !  integer(ip)                         , intent(in)    :: maxit    
  !  class(abstract_environment_t)       , intent(in)    :: env
  !  type(p_discrete_integration_t), intent(inout) :: approx(:)
  !  type(serial_fe_space_t)                    , intent(inout) :: fe_space
  !  class(operator_t)              , intent(inout) :: A, M
  !  class(vector_t)               , intent(inout) :: x, b
  !  ! Locals
  !  integer(ip) :: iiter
  !  real(rp)    :: resnorm,ininorm
  !  
  !  iiter = 0
  !  do while( iiter < maxit )

  !     ! Update counter
  !     iiter = iiter+1

  !     ! Initialize Matrix and vector
  !     ! ***************** Abstract procedure to initialize a abstract_operator ************************!
  !     select type (A)
  !     type is(serial_scalar_matrix_t)
  !        call A%init(0.0_rp)
  !     class default
  !        check(.false.)
  !     end select
  !     !********************************************************************************************!
  !     call b%init(0.0_rp)

  !     ! Integrate system
  !     call volume_integral(approx,fe_space,A,b)

  !     ! Check convergence
  !     if(iiter==1) ininorm = b%nrm2()
  !     x = b - A*x
  !     resnorm = x%nrm2()
  !     if( resnorm < nltol*ininorm) then
  !        write(*,*) 'Nonlinear iterations: ', iiter
  !        write(*,*) 'Nonlinear error norm: ', resnorm
  !        exit
  !     end if

  !     ! Compute Numeric preconditioner
  !     ! ***************** Abstract procedure to compute precond numeric ***************************!
  !     select type (A)
  !     type is(serial_scalar_matrix_t)
  !        select type (M)
  !        type is(preconditioner_t)
  !           call preconditioner_numeric(M)
  !        class default
  !           check(.false.)
  !        end select
  !     class default
  !        check(.false.)
  !     end select
  !     !********************************************************************************************!

  !     ! Solve system
  !     call abstract_solve(A,M,b,x,sctrl,env)
  !     call solver_control_log_conv_his(sctrl)
  !     call solver_control_free_conv_his(sctrl)

  !     ! Free Numeric preconditioner
  !     ! ******************** Abstract procedure to free precond numeric ***************************!
  !     select type (M)
  !     type is(preconditioner_t)
  !        call preconditioner_free(preconditioner_free_values,M)
  !        class default
  !        check(.false.)
  !     end select
  !     !********************************************************************************************!
  !     
  !     ! Store solution to unkno
  !     ! ***************** Abstract procedure to update from a base operant ************************!
  !     select type (x)
  !     type is(serial_scalar_array_t)
  !        call update_solution(x,fe_space)
  !     class default
  !        check(.false.)
  !     end select
  !     !********************************************************************************************!
  !     
  !     ! CORRECT
  !     ! Store nonlinear iteration ( k+1 --> k )
  !     !call update_nonlinear(fe_space)
  !     
  !  end do

  !end subroutine nonlinear_iteration
  
end program test_nsi_iss
