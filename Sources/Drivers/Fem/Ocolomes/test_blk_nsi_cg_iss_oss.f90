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
program test_blk_nsi_cg_iss_oss
  use serial_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"

  ! Types
  type(geom_data_t)                       :: gdata
  type(bound_data_t)                      :: bdata
  type(reference_element_t)               :: geo_reference_element
  type(triangulation_t)                   :: f_trian
  type(conditions_t)                      :: f_cond
  type(dof_handler_t)                     :: dhand
  type(fe_space_t)                        :: fe_space  
  type(nsi_problem_t)                     :: myprob
  type(nsi_cg_iss_oss_discrete_t), target :: mydisc
  type(nsi_cg_iss_oss_matvec_t)  , target :: matvec
  type(discrete_integration_pointer)      :: approx(1)
  type(block_matrix_t)           , target :: blmat
  type(block_vector_t)           , target :: blvec,blunk
  type(block_precond_t)                   :: blprec
  type(precond_params_t)                  :: ppars
  type(solver_control_t)                  :: sctrl
  type(serial_environment_t)              :: senv
  type(vtk_t)                             :: fevtk
  class(base_operand_t)         , pointer :: x, b
  class(base_operator_t)        , pointer :: A, M
  type(graph_t)                 , pointer :: graph
  type(block_graph_t)                     :: blk_graph
  type(plain_vector_t)                    :: enorm
  type(nsi_cg_iss_oss_error_t)   , target :: ecalc

  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat
  integer(ip) :: num_approximations = 1

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
  integer(ip)        :: nex,ney,nez

  call meminit

    ! Generate geometry data
  call geom_data_create(gdata,nex,ney,nez)

  ! Generate boundary data
  call bound_data_create(2*gdata%ndime+1,2*gdata%ndime+1,gdata%ndime,bdata)
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
  call gen_triangulation(1,gdata,bdata,geo_reference_element,f_trian,f_cond,material)

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call matvec%create(myprob,mydisc)
  approx(1)%p       => matvec
  mydisc%dtinv      = 0.0_rp
  myprob%kfl_conv   = 1
  myprob%diffu      = 1.0_rp
  myprob%case_veloc = 1
  myprob%case_press = 1

  ! Create dof_handler
  call dhand%create(3,1,mydisc%nvars,vars_block,dof_coupling)
  call dhand%set_problem(1,mydisc)

  ! Allocate auxiliar elemental arrays
  call memalloc(f_trian%num_elems,dhand%nvars_global,continuity, __FILE__,__LINE__)
  call memalloc(f_trian%num_elems,dhand%nvars_global,order,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,problem,__FILE__,__LINE__)
  call memalloc(f_trian%num_elems,which_approx,__FILE__,__LINE__)
  continuity             = 1
  order                  = 2
  order(:,gdata%ndime+1) = 1
  problem                = 1
  which_approx           = 1 
  
  ! Create fe_space
  call fe_space_create(f_trian,dhand,fe_space,problem,f_cond,continuity,order,material,which_approx, &
       &               time_steps_to_store=3, hierarchical_basis=.false.,                            &
       &               static_condensation=.false.,num_continuity=1)

  ! Create plain vectors
  call fe_space_plain_vector_create((/2/),fe_space)
  call fe_space_plain_vector_point(2,fe_space)

  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,myprob,senv,dir_path_out,prefix,linear_order=.true.)

  ! Create dof info
  call create_dof_info(dhand,f_trian,fe_space,blk_graph,gtype)

  ! Allocate matrices and vectors
  call blmat%alloc(blk_graph)
  call blvec%alloc(blk_graph)
  call blunk%alloc(blk_graph)
  call blvec%init(0.0_rp)

  ! Apply boundary conditions to unkno
  call update_strong_dirichlet_bcond(fe_space,f_cond)
  call update_analytical_bcond((/1:gdata%ndime/),myprob%case_veloc,0.0_rp,fe_space)
  call update_analytical_bcond((/gdata%ndime+1/),myprob%case_press,0.0_rp,fe_space)

  ! Solver control parameters
  sctrl%method = direct
  sctrl%trace  = 100
  sctrl%track_conv_his = .true.
  
end program test_blk_nsi_cg_iss_oss
