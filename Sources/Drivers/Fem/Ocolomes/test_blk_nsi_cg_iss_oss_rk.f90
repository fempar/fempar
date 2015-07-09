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
module my_linear_algebra_names
  use serial_names
  implicit none
# include "debug.i90"
  private
  
  type my_linear_algebra_t 
     type(block_matrix_t)           :: mass_u_matrix
     type(block_matrix_t)           :: lapla_p_matrix
     type(block_matrix_t)           :: block_matrix_u
     type(block_matrix_t)           :: block_matrix_p
     type(block_vector_t)           :: block_vector_u
     type(block_vector_t)           :: block_vector_p
     type(block_vector_t)           :: block_unknown_u
     type(block_vector_t)           :: block_unknown_p
     type(block_operator_t)         :: block_operator_u
     type(block_operator_t)         :: block_operator_p
     type(block_operand_t)          :: block_operand_vec_u
     type(block_operand_t)          :: block_operand_vec_p
     type(block_operand_t)          :: block_operand_unk_u
     type(block_operand_t)          :: block_operand_unk_p
     type(block_preconditioner_l_t) :: block_preconditioner_u
     type(block_preconditioner_l_t) :: block_preconditioner_p
     type(preconditioner_t)         :: u_preconditioner
     type(preconditioner_t)         :: x_preconditioner
     type(preconditioner_t)         :: w_preconditioner
     type(preconditioner_t)         :: p_preconditioner
     type(preconditioner_params_t)  :: u_ppars
     type(preconditioner_params_t)  :: x_ppars
     type(preconditioner_params_t)  :: w_ppars
     type(preconditioner_params_t)  :: p_ppars     
   contains
     procedure :: build => build_linear_algebra
     procedure :: free  => free_linear_algebra
  end type my_linear_algebra_t

 ! Types
  public :: my_linear_algebra_t
  
contains

  !==================================================================================================
  subroutine build_linear_algebra(la,blk_graph)
    implicit none
    class(my_linear_algebra_t), intent(inout) :: la
    type(block_graph_t)       , intent(in)    :: blk_graph

    ! Allocate matrices and vectors
    call la%block_matrix_u%alloc(blk_graph)
    call la%block_matrix_u%set_block_to_zero(1,2)
    call la%block_matrix_u%set_block_to_zero(2,1)
    call la%block_matrix_u%set_block_to_zero(2,2)
    call la%block_matrix_u%set_block_to_zero(2,3)
    call la%block_matrix_u%set_block_to_zero(3,2)
    call la%block_vector_u%alloc(blk_graph)
    call la%block_unknown_u%alloc(blk_graph)
    call la%block_vector_u%init(0.0_rp)
    call la%block_matrix_p%alloc(blk_graph)
    call la%block_matrix_p%set_block_to_zero(1,3)
    call la%block_matrix_p%set_block_to_zero(2,2)
    call la%block_matrix_p%set_block_to_zero(2,3)
    call la%block_matrix_p%set_block_to_zero(3,1)
    call la%block_matrix_p%set_block_to_zero(3,2)
    call la%block_matrix_p%set_block_to_zero(3,3)
    call la%block_vector_p%alloc(blk_graph)
    call la%block_unknown_p%alloc(blk_graph)
    call la%block_vector_p%init(0.0_rp)

    ! Auxiliar matrices
    call la%mass_u_matrix%alloc(blk_graph)
    call la%mass_u_matrix%set_block_to_zero(1,2)
    call la%mass_u_matrix%set_block_to_zero(1,3)
    call la%mass_u_matrix%set_block_to_zero(2,1)
    call la%mass_u_matrix%set_block_to_zero(2,2)
    call la%mass_u_matrix%set_block_to_zero(2,3)
    call la%mass_u_matrix%set_block_to_zero(3,1)
    call la%mass_u_matrix%set_block_to_zero(3,2)
    call la%mass_u_matrix%set_block_to_zero(3,3)
    call la%lapla_p_matrix%alloc(blk_graph)
    call la%lapla_p_matrix%set_block_to_zero(1,1)
    call la%lapla_p_matrix%set_block_to_zero(1,2)
    call la%lapla_p_matrix%set_block_to_zero(1,3)
    call la%lapla_p_matrix%set_block_to_zero(2,1)
    call la%lapla_p_matrix%set_block_to_zero(2,3)
    call la%lapla_p_matrix%set_block_to_zero(3,1)
    call la%lapla_p_matrix%set_block_to_zero(3,2)
    call la%lapla_p_matrix%set_block_to_zero(3,3)
    
    ! Construct U-preconditioner (K^-1)
    la%u_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%block_matrix_u%get_block(1,1),la%u_preconditioner,la%u_ppars)
    call preconditioner_symbolic(la%block_matrix_u%get_block(1,1),la%u_preconditioner)
    call preconditioner_log_info(la%u_preconditioner)

    ! Construct X-preconditioner (Mx^-1)
    la%x_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%block_matrix_u%get_block(3,3),la%x_preconditioner,la%x_ppars)
    call preconditioner_symbolic(la%block_matrix_u%get_block(3,3),la%x_preconditioner)
    call preconditioner_log_info(la%x_preconditioner)

    ! Construct W-preconditioner (Mu^-1)
    la%w_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%mass_u_matrix%get_block(1,1),la%w_preconditioner,la%w_ppars)
    call preconditioner_symbolic(la%mass_u_matrix%get_block(1,1),la%w_preconditioner)
    call preconditioner_log_info(la%w_preconditioner)

    ! Construct P-preconditioner (Lp^-1)
    la%p_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%lapla_p_matrix%get_block(2,2),la%p_preconditioner,la%p_ppars)
    call preconditioner_symbolic(la%lapla_p_matrix%get_block(2,2),la%p_preconditioner)
    call preconditioner_log_info(la%p_preconditioner)
    
    ! Create U-X Block operator
    call la%block_operator_u%create(2,2)
    call la%block_operator_u%set_block(1,1,la%block_matrix_u%get_block(1,1))
    call la%block_operator_u%set_block(1,2,la%block_matrix_u%get_block(1,3))
    call la%block_operator_u%set_block(2,1,la%block_matrix_u%get_block(3,1))
    call la%block_operator_u%set_block(2,2,la%block_matrix_u%get_block(3,3))

    ! Create U-X Block operand
    call la%block_operand_vec_u%create(2)
    call la%block_operand_vec_u%set_block(1,la%block_vector_u%blocks(1))
    call la%block_operand_vec_u%set_block(2,la%block_vector_u%blocks(3))
    call la%block_operand_unk_u%create(2)
    call la%block_operand_unk_u%set_block(1,la%block_unknown_u%blocks(1))
    call la%block_operand_unk_u%set_block(2,la%block_unknown_u%blocks(3))

    ! Create U-X Block preconditioner
    call la%block_preconditioner_u%create(2)
    call la%block_preconditioner_u%set_block(1,1,la%u_preconditioner)
    call la%block_preconditioner_u%set_block(2,1,la%block_matrix_u%get_block(3,1))
    call la%block_preconditioner_u%set_block(2,2,la%x_preconditioner)
    
    ! Create W-P Block operator
    call la%block_operator_p%create(2,2)
    call la%block_operator_p%set_block(1,1,la%block_matrix_p%get_block(1,1))
    call la%block_operator_p%set_block(1,2,la%block_matrix_p%get_block(1,2))
    call la%block_operator_p%set_block(2,1,la%block_matrix_p%get_block(2,1))
    call la%block_operator_p%set_block(2,2,la%block_matrix_p%get_block(2,2))

    ! Create W-P Block operand
    call la%block_operand_vec_p%create(2)
    call la%block_operand_vec_p%set_block(1,la%block_vector_p%blocks(1))
    call la%block_operand_vec_p%set_block(2,la%block_vector_p%blocks(2))
    call la%block_operand_unk_p%create(2)
    call la%block_operand_unk_p%set_block(1,la%block_unknown_p%blocks(1))
    call la%block_operand_unk_p%set_block(2,la%block_unknown_p%blocks(2))

    ! Create W-P Block preconditioner
    call la%block_preconditioner_p%create(2)
    call la%block_preconditioner_p%set_block(1,1,la%w_preconditioner)
    call la%block_preconditioner_p%set_block(2,1,la%block_matrix_p%get_block(2,1))
    call la%block_preconditioner_p%set_block(2,2,la%p_preconditioner)
    
  end subroutine build_linear_algebra

  !==================================================================================================
  subroutine free_linear_algebra(la)
    implicit none
    class(my_linear_algebra_t), intent(inout) :: la

    ! Destroy W-P Block preconditioner
    call la%block_preconditioner_p%destroy()

    ! Destroy W-P Block operand
    call la%block_operand_vec_p%destroy()
    call la%block_operand_unk_p%destroy()

    ! Destroy W-P Block operand
    call la%block_operator_p%destroy()

    ! Destroy U-X Block preconditioner
    call la%block_preconditioner_u%destroy()

    ! Destroy U-X Block operand
    call la%block_operand_vec_u%destroy()
    call la%block_operand_unk_u%destroy()

    ! Destroy U-X Block operand
    call la%block_operator_u%destroy()

    ! Destroy P-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%p_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%p_preconditioner)

    ! Destroy W-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%w_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%w_preconditioner)

    ! Destroy X-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%x_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%x_preconditioner)

    ! Destroy U-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%u_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%u_preconditioner)

    ! Deallocate auxiliar matrices
    call la%mass_u_matrix%free()
    call la%lapla_p_matrix%free()
          
    ! Deallocate problem matrix
    call la%block_matrix_u%free()
    call la%block_matrix_p%free()

    ! Deallocate problem vectors
    call la%block_vector_u%free()
    call la%block_vector_p%free()
    call la%block_unknown_u%free()
    call la%block_unknown_p%free()

  end subroutine free_linear_algebra

end module my_linear_algebra_names

program test_blk_nsi_cg_iss_oss_rk
  use serial_names
  use my_linear_algebra_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use norm_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"

  ! Types
  type(uniform_mesh_descriptor_t)                     :: gdata
  type(uniform_conditions_descriptor_t)               :: bdata
  type(reference_element_t)                           :: geo_reference_element
  type(triangulation_t)                               :: f_trian
  type(conditions_t)                                  :: f_cond
  type(dof_descriptor_t)                              :: dof_descriptor
  type(fe_space_t)                                    :: fe_space  
  type(nsi_problem_t)                                 :: myprob
  type(nsi_cg_iss_oss_discrete_t)            , target :: mydisc
  type(nsi_cg_iss_oss_rk_momentum_t)         , target :: cg_iss_oss_rk_momentum
  type(nsi_cg_iss_oss_rk_pressure_t)         , target :: cg_iss_oss_rk_pressure
  type(nsi_cg_iss_oss_rk_momentum_update_t)  , target :: cg_iss_oss_rk_momentum_update
  type(nsi_cg_iss_oss_rk_projection_update_t), target :: cg_iss_oss_rk_projection_update
  type(discrete_integration_pointer_t)                :: approx(1)
  type(preconditioner_params_t)                       :: ppars
  type(solver_control_t)                              :: sctrl
  type(serial_environment_t)                          :: senv
  type(vtk_t)                                         :: fevtk
  class(base_operand_t)                     , pointer :: x, b
  class(base_operator_t)                    , pointer :: A, M
  type(block_graph_t)                                 :: blk_graph
  type(scalar_t)                                      :: enorm_u, enorm_p
  type(error_norm_t)                         , target :: error_compute
  type(nsi_cg_iss_oss_lapla_p_t)             , target :: lapla_p_integration
  type(nsi_cg_iss_oss_massu_t)               , target :: mass_u_integration
  type(my_linear_algebra_t)                           :: linear_algebra
  type(rungekutta_integrator_t)              , target :: rkinteg

  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat,i
  integer(ip) :: num_approximations = 1
  integer(ip) :: setterms(6,2),settable(3)

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2
  integer(ip), parameter :: linear=0, nonlinear=1
  integer(ip), parameter :: explicit=0, implicit=1

  ! Allocatable
  integer(ip), allocatable :: continuity(:,:)
  integer(ip), allocatable :: order(:,:)
  integer(ip), allocatable :: material(:)
  integer(ip), allocatable :: problem(:)
  integer(ip), allocatable :: which_approx(:)
  integer(ip), allocatable :: vars_block(:)
  integer(ip), allocatable :: dof_coupling(:,:)

  ! Arguments
  character(len=256) :: dir_path_out,prefix
  integer(ip)        :: nex,ney,nez,nstage,rk_order,rk_flag

  call meminit

  ! Read parameters from command-line
  call read_pars_cl_test_blk_nsi_cg_iss_oss_rk(prefix,dir_path_out,nex,ney,nez,nstage,rk_order,rk_flag)

  ! Generate geometry data
  call uniform_mesh_descriptor_create(gdata,nex,ney,nez)

  ! Generate boundary data
  call uniform_conditions_descriptor_create(2*gdata%ndime+1,2*gdata%ndime+1,gdata%ndime,bdata)
  bdata%poin%code(gdata%ndime+1,1:2**gdata%ndime-1) = 0
  bdata%line%code(gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+1,:) = 0
  bdata%poin%code(gdata%ndime+1,2**gdata%ndime) = 1
  bdata%poin%valu(gdata%ndime+1,2**gdata%ndime) = 0.0_rp
  bdata%poin%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%line%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%surf%valu(1:gdata%ndime,:) = 1.0_rp
  bdata%poin%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0
  bdata%line%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0
  bdata%surf%code(gdata%ndime+2:2*gdata%ndime+1,:) = 0

  ! Generate element geometrical fixed info
  call finite_element_fixed_info_create(geo_reference_element,Q_type_id,1,gdata%ndime)

  ! Generate triangulation
  call generate_uniform_triangulation(1,gdata,bdata,geo_reference_element,f_trian,f_cond,material)

  ! Define Runge-Kutta method
  settable      = (/nstage,rk_order,rk_flag/)
  setterms(1,:) = (/linear,implicit/)     ! Diffusion
  setterms(2,:) = (/nonlinear,implicit/)  ! Convection
  setterms(3,:) = (/linear,explicit/)     ! Pressure Gradient
  setterms(4,:) = (/nonlinear,implicit/)  ! OSS_vu
  setterms(5,:) = (/nonlinear,implicit/)  ! OSS_vx
  setterms(6,:) = (/linear,explicit/)     ! Force
  call rkinteg%create(setterms,settable)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call cg_iss_oss_rk_momentum%create(myprob,mydisc)
  call cg_iss_oss_rk_pressure%create(myprob,mydisc)
  cg_iss_oss_rk_momentum%rkinteg => rkinteg
  mydisc%dtinv    = 0.0_rp
  mydisc%kfl_proj = 1
  mydisc%kfl_lump = 1
  myprob%kfl_conv = 1
  myprob%diffu    = 1.0_rp

  ! Create dof_descriptor
  call dof_descriptor%create(3,1,mydisc%nvars,vars_block,dof_coupling)
  call dof_descriptor%set_problem(1,mydisc)

  ! Allocate auxiliar elemental arrays
  call memalloc(f_trian%num_elems,dof_descriptor%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dof_descriptor%nvars_global,order,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  order                  = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 
  
  ! Create fe_space
  call fe_space_create(f_trian,dof_descriptor,fe_space,problem,f_cond,continuity,order,material, &
       &               which_approx,time_steps_to_store=3, hierarchical_basis=.false.,           &
       &               static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,myprob,senv,dir_path_out,prefix,linear_order=.true.)

  ! Create dof info
  call create_dof_info(dof_descriptor,f_trian,fe_space,blk_graph,gtype)

  ! Assign analytical solution
  if(gdata%ndime==2) then
     call fe_space%set_analytical_code((/4,5,3,0,0/),(/0,0,0,0,0/))
  else
     write(*,*) 'analytical function not ready for 3D'
  end if

  ! Create linear algebra structures
  call linear_algebra%build(blk_graph)

  ! Compute auxiliar matrices
  call lapla_p_integration%create(myprob,mydisc)
  approx(1)%p => lapla_p_integration
  call volume_integral(approx,fe_space,linear_algebra%lapla_p_matrix)
  call mass_u_integration%create(myprob,mydisc)
  approx(1)%p => mass_u_integration
  call volume_integral(approx,fe_space,linear_algebra%mass_u_matrix)

  ! Do time steps
!!$  call time_steps_rk(sctrl,sctrl%rtol*1.0e2_rp,10,senv,fe_space,linear_algebra,                  &
!!$       &            cg_iss_oss_rk_momentum,cg_iss_oss_rk_pressure,cg_iss_oss_rk_momentum_update, &
!!$       &            cg_iss_oss_rk_projection_update)

  ! Print solution to VTK file
  istat = fevtk%write_VTK()

  ! Compute error norm
  call error_compute%create(myprob,mydisc)
  approx(1)%p => error_compute
  error_compute%unknown_id = velocity
  call enorm_u%init()
  call volume_integral(approx,fe_space,enorm_u)
  error_compute%unknown_id = pressure
  call enorm_p%init()
  call volume_integral(approx,fe_space,enorm_p)
  write(*,*) 'Velocity error norm: ', sqrt(enorm_u%get())
  write(*,*) 'Pressure error norm: ', sqrt(enorm_p%get()) 

  ! Free linear algebra structures
  call linear_algebra%free()

  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  call memfree(vars_block,__FILE__,__LINE__)
  call memfree(dof_coupling,__FILE__,__LINE__)
  call fevtk%free
  call blk_graph%free()
  call fe_space_free(fe_space) 
  call myprob%free
  call mydisc%free
  call error_compute%free
  call cg_iss_oss_rk_momentum%free
  call cg_iss_oss_rk_pressure%free
  call cg_iss_oss_rk_pressure%free
  call lapla_p_integration%free
  call mass_u_integration%free
  call rkinteg%free
  call dof_descriptor_free(dof_descriptor)
  call triangulation_free(f_trian)
  call conditions_free(f_cond)
  call finite_element_fixed_info_free(geo_reference_element)
  call uniform_conditions_descriptor_free(bdata)

  call memstatus

contains

  !==================================================================================================
  subroutine read_pars_cl_test_blk_nsi_cg_iss_oss_rk(prefix,dir_path_out,nex,ney,nez,nstage,order,flag)
    implicit none
    character*(*), intent(out) :: prefix, dir_path_out
    integer(ip)  , intent(out) :: nex,ney,nez,nstage,order,flag
    character(len=256)         :: program_name
    character(len=256)         :: argument 
    integer                    :: numargs,iargc

    numargs = iargc()
    call getarg(0, program_name)
    if (.not. (numargs==5) ) then
       write (6,*) 'Usage: ', trim(program_name), ' prefix dir_path_out nex ney nez nstage order flag'
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

    call getarg(6, argument)
    read (argument,*) nstage

    call getarg(7, argument)
    read (argument,*) order

    call getarg(8, argument)
    read (argument,*) flag

  end subroutine read_pars_cl_test_blk_nsi_cg_iss_oss_rk

  !==================================================================================================
  subroutine time_steps_rk(sctrl,nltol,maxit,sttol,maxst,env,fe_space,la,momentum,pressure, &
       &                   momentum_update,projection_update)
    implicit none
    type(solver_control_t)               , intent(inout) :: sctrl
    real(rp)                             , intent(in)    :: nltol,sttol
    integer(ip)                          , intent(in)    :: maxit,maxst
    class(abstract_environment_t)        , intent(in)    :: env
    type(fe_space_t)                     , intent(inout) :: fe_space
    type(my_linear_algebra_t)    , target, intent(inout) :: la
    class(discrete_integration_t), target, intent(inout) :: momentum
    class(discrete_integration_t), target, intent(inout) :: pressure
    class(discrete_integration_t), target, intent(inout) :: momentum_update 
    class(discrete_integration_t), target, intent(inout) :: projection_update
    ! Locals
    type(discrete_integration_pointer_t)   :: approx(1)
    type(rungekutta_integrator_t), pointer :: rkinteg
    integer(ip) :: istage,nstage,istep
    real(rp)    :: rtime,ctime,prevtime

    select type(momentum)
    class is(nsi_cg_iss_oss_rk_momentum_t)
       rkinteg => momentum%rkinteg
    class default
       write(0,'(a)') 'time_steps_rk: unsupported momentum class'
       check(1==0)
    end select

    ! Initialize time steps
    rtime = 1.0_rp
    istep = 0

!!$    ! Time steps loop
!!$    step: do while (momentum%discret%ctime<ftime.and.rtime>sttol.and.istep<=maxst)
!!$
!!$       ! Current time
!!$       prevtime = momentum%discret%ctime
!!$       momentum%discret%ctime = momentum%discret%ctime + 1.0_rp/momentum%discret%dtinv
!!$       write(*,*) '============================================================'
!!$       write(*,*) 'Time step: ',istep+1,',  Current time: ',momentum%discret%ctime
!!$
!!$       ! Loop over stages
!!$       stage: do istage=1,rkinteg%rktable(1)%p%stage
!!$
!!$          ! Set current time
!!$          ctime = prevtime + 1.0_rp/dtinv*rkinteg%rktable(1)%p%c(istage)
!!$          rkinteg%istge = istage
!!$          
!!$          ! Update boundary conditions
!!$          call update_strong_dirichlet_bcond(fe_space,f_cond)
!!$          call update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,fe_space)
!!$          
!!$          ! Momentum equation solution
!!$          !call nonlinear_iteration(sctrl,sctrl%rtol*1.0e2_rp,10,senv,approx,fe_space,linear_algebra)
!!$
!!$          ! Update boundary conditions (velocity derivative)
!!$
!!$          ! Pressure equation solution
!!$
!!$       end do stage
!!$
!!$       ! Update boundary conditions
!!$
!!$       ! Momentum update equation solution
!!$
!!$       ! Projection update equation solution
!!$       
!!$       ! Pressure equation solution
!!$
!!$       ! Check steady state
!!$
!!$    end do step

  end subroutine time_steps_rk
  
end program test_blk_nsi_cg_iss_oss_rk
