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
     type(block_matrix_t) :: lapla_p_matrix
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine build_linear_algebra

  !==================================================================================================
  subroutine free_linear_algebra(la)
    implicit none
    class(my_linear_algebra_t), intent(inout) :: la
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  type(uniform_mesh_descriptor_t)            :: gdata
  type(uniform_conditions_descriptor_t)      :: bdata
  type(reference_element_t)                  :: geo_reference_element
  type(triangulation_t)                      :: f_trian
  type(conditions_t)                         :: f_cond
  type(dof_descriptor_t)                     :: dof_descriptor
  type(fe_space_t)                           :: fe_space  
  type(nsi_problem_t)                        :: myprob
  type(nsi_cg_iss_oss_discrete_t)   , target :: mydisc
  type(nsi_cg_iss_oss_rk_momentum_t), target :: cg_iss_oss_rk_momentum
  type(nsi_cg_iss_oss_rk_pressure_t), target :: cg_iss_oss_rk_pressure
  type(discrete_integration_pointer_t)       :: approx(1)
  type(preconditioner_params_t)              :: ppars
  type(solver_control_t)                     :: sctrl
  type(serial_environment_t)                 :: senv
  type(vtk_t)                                :: fevtk
  class(base_operand_t)            , pointer :: x, b
  class(base_operator_t)           , pointer :: A, M
  type(block_graph_t)                        :: blk_graph
  type(scalar_t)                             :: enorm_u, enorm_p
  type(error_norm_t)                , target :: error_compute
  type(nsi_cg_iss_oss_lapla_p_t)    , target :: lapla_p_integration
  type(my_linear_algebra_t)                  :: linear_algebra
  type(rungekutta_integrator_t)     , target :: rkinteg

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
  call fe_space_create(f_trian,dof_descriptor,fe_space,problem,f_cond,continuity,order,material,which_approx, &
       &               time_steps_to_store=3, hierarchical_basis=.false.,                            &
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

  ! Compute auxiliar matrix
  call lapla_p_integration%create(myprob,mydisc)
  approx(1)%p => lapla_p_integration
  call volume_integral(approx,fe_space,linear_algebra%lapla_p_matrix)

  ! Do time steps
  call time_steps_rk(sctrl,sctrl%rtol*1.0e2_rp,10,senv,approx,fe_space,linear_algebra)

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
  subroutine time_steps_rk( sctrl, nltol, maxit, env, approx, fe_space, la )
    implicit none
    type(solver_control_t)              , intent(inout) :: sctrl
    real(rp)                            , intent(in)    :: nltol
    integer(ip)                         , intent(in)    :: maxit    
    class(abstract_environment_t)       , intent(in)    :: env
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    type(my_linear_algebra_t), target   , intent(inout) :: la

    !!!!!!!!!

  end subroutine time_steps_rk
  
end program test_blk_nsi_cg_iss_oss_rk
