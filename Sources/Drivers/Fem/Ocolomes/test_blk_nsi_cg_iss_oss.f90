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
  type(nsi_cg_iss_oss_massp_t)   , target :: mass_p_integration
  type(block_matrix_t)                    :: blmass_p
  type(block_operator_t)                  :: block_operator
  type(block_operand_t)                   :: block_operand
  type(block_preconditioner_l_t)          :: block_precond

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

  ! Read parameters from command-line
  call read_pars_cl_test_blk_nsi_cg_iss_oss(prefix,dir_path_out,nex,ney,nez)

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

  ! Auxiliar matrices
  call blmass_p%alloc(blk_graph)
  call blmass_p%set_block_to_zero(2,2)

  ! Compute auxiliar matrix
  call mass_p_compute%create(myprob,mydisc)
  call volume_integral(approx,fe_space,blmass_p)

  ! Apply boundary conditions to unkno
  call update_strong_dirichlet_bcond(fe_space,f_cond)
  call update_analytical_bcond((/1:gdata%ndime/),myprob%case_veloc,0.0_rp,fe_space)
  call update_analytical_bcond((/gdata%ndime+1/),myprob%case_press,0.0_rp,fe_space)

  ! Solver control parameters
  sctrl%method = direct
  sctrl%trace  = 100
  sctrl%track_conv_his = .true.

  ! Build block structures
  call build_block_structures(blmat,blvec,blmass_p%get_block(2,2),blk_operator,blk_operand,blk_precond)

  ! Point discrete integration
  approx(1)%p => matvec

  !!!!!!!!!!!!!!!!!1

contains

  !==================================================================================================
  subroutine read_pars_cl_test_blk_nsi_cg_iss_oss(prefix,dir_path_out,nex,ney,nez)
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

  end subroutine read_pars_cl_test_blk_nsi_cg_iss_oss

  !==================================================================================================
  subroutine build_block_structures(block_matrix,block_vector,mass_p_matrix,block_operator, &
       &                            block_operand,block_preconditioner)
    implicit none
    type(block_matrix_t)          , intent(in)  :: block_matrix
    type(block_vector_t)          , intent(in)  :: block_vector
    type(matrix_t)                , intent(in)  :: mass_p_matrix
    type(block_operator_t)        , intent(out) :: block_operator
    type(block_operand_t)         , intent(out) :: block_operand
    type(block_preconditioner_l_t), intent(out) :: block_preconditioner
    ! Locals
    type(preconditioner_t)        :: u_preconditioner
    type(preconditioner_t)        :: p_preconditioner
    type(preconditioner_t)        :: x_preconditioner
    type(preconditioner_params_t) :: u_ppars
    type(preconditioner_params_t) :: p_ppars
    type(preconditioner_params_t) :: x_ppars
    type(block_operator_t)        :: block_up_operator
    type(block_operator_t)        :: block_up_x_operator
    type(block_operator_t)        :: block_x_up_operator
    type(block_operand_t)         :: block_up_operand
    type(block_preconditioner_t)  :: block_up_preconditioner
    
    ! Construct U-preconditioner (K^-1)
    u_ppars%type = pardiso_mkl_prec
    call preconditioner_create(block_matrix%get_block(1,1),u_preconditioner,u_ppars)
    call preconditioner_symbolic(block_matrix%get_block(1,1),u_preconditioner)
    call preconditioner_log_info(u_preconditioner)

    ! Construct P-preconditioner (M^-1)
    p_ppars%type = pardiso_mkl_prec
    call preconditioner_create(mass_p_matrix,p_preconditioner,p_ppars)
    call preconditioner_symbolic(mass_p_matrix,p_preconditioner)
    call preconditioner_log_info(p_preconditioner)

    ! Construct X-preconditioner (Mx^-1)
    x_ppars%type = pardiso_mkl_prec
    call preconditioner_create(block_matrix%get_block(3,3),x_preconditioner,x_ppars)
    call preconditioner_symbolic(mass_p_matrix,x_preconditioner)
    call preconditioner_log_info(x_preconditioner)
    
    ! Create U-P Block operator
    call block_up_operator%create(2,2)
    call block_up_operator%set_block(1,1,block_matrix%get_block(1,1))
    call block_up_operator%set_block(1,2,block_matrix%get_block(1,2))
    call block_up_operator%set_block(2,1,block_matrix%get_block(2,1))
    call block_up_operator%set_block(2,2,block_matrix%get_block(2,2))

    ! Create U-P Block operand
    call block_up_operand%create(2)
    call block_up_operand%set_block(1,block_vector%blocks(1))
    call block_up_operand%set_block(2,block_vector%blocks(2))

    ! Create U-P Block preconditioner
    call block_up_preconditioner%create(2)
    call block_up_preconditioner%set_block(1,1,u_preconditioner)
    call block_up_preconditioner%set_block(2,1,block_matrix%get_block(2,1))
    call block_up_preconditioner%set_block(2,2,p_preconditioner)

    ! Create UP-X Block operator
    call block_up_x_operator%create(2,1)
    call block_up_x_operator%set_block(1,1,block_matrix%get_block(1,3)
    call block_up_x_operator%set_block_to_zero(2,1) 

    ! Create X-UP Block operator
    call block_x_up_operator%create(1,2)
    call block_x_up_operator%set_block(1,1,block_matrix%get_block(3,1)
    call block_x_up_operator%set_block_to_zero(1,2)

    ! Create Global Block operator
    call block_operator%create(2,2)
    call block_operator%set_block(1,1,block_up_operator)
    call block_operator%set_block(1,2,block_up_x_operator)
    call block_operator%set_block(2,1,block_x_up_operator)
    call block_operator%set_block(2,2,block_matrix%get_bloc(3,3))

    ! Create Global Block operand
    call block_operand%create(2)
    call block_operand%set_block(1,block_up_operand)
    call block_operand%set_block(2,block_vector%blocks(3))

    ! Create Global Block preconditioner
    call block_preconditioner%create(2,2)
    call block_preconditioner%set_block(1,1,block_up_preconditioner)
    call block_preconditioner%set_block(2,1,block_x_up_operator)
    call block_preconditioner%set_block(2,2,x_preconditioner)
    
  end subroutine build_block_structures

  !==================================================================================================
  subroutine compute_preconditioner(A,M)
    implicit none
    class(base_operator_t), intent(in)    :: A
    class(base_operator_t), intent(inout) :: M

    ! precond numeric (K^-1)
    select type (A)
    type is(block_operator_t)
       select type(A%get_block(1,1))
       type is(block_operator_t)
          select type(A%get_block(1,1)%get_block(1,1))
          type is(matrix_t)
             select type (M)
             type is(block_preconditioner_t)
                select type(M%get_block(1,1))
                type is(block_preconditioner_t)
                   select type(M%get_block(1,1)%get_block(1,1))
                   type is(preconditioner_t)
                      call preconditioner_numeric(A%get_block(1,1)%get_block(1,1), &
                           &                      M%get_block(1,1)%get_block(1,1))
                   class default
                      check(.false.)
                   end select
                class default
                   check(.false.)
                end select
             class default
                check(.false.)
             end select
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond numeric (Mp^-1)
    select type (A)
    type is(block_operator_t)
       select type(A%get_block(1,1))
       type is(block_operator_t)
          select type(A%get_block(1,1)%get_block(2,2))
          type is(matrix_t)
             select type (M)
             type is(block_preconditioner_t)
                select type(M%get_block(1,1))
                type is(block_preconditioner_t)
                   select type(M%get_block(1,1)%get_block(2,2))
                   type is(preconditioner_t)
                      call preconditioner_numeric(A%get_block(1,1)%get_block(2,2), &
                           &                      M%get_block(1,1)%get_block(2,2))
                   class default
                      check(.false.)
                   end select
                class default
                   check(.false.)
                end select
             class default
                check(.false.)
             end select
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond numeric (Mx^-1)
    select type (A)
    type is(block_operator_t)
       select type(A%get_block(2,2))
       type is(matrix_t)
          select type (M)
          type is(block_preconditioner_t)
             select type(M%get_block(2,2))
             type is(preconditioner_t)
                call preconditioner_numeric(A%get_block(2,2),M%get_block(2,2))
             class default
                check(.false.)
             end select
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select
          
  end subroutine compute_preconditioner
  
  !==================================================================================================
  subroutine free_preconditioner(M)
    implicit none
    class(base_operator_t), intent(inout) :: M

    ! precond free (K^-1)
    select type (M)
    type is(block_preconditioner_t)
       select type(M%get_block(1,1))
       type is(block_preconditioner_t)
          select type(M%get_block(1,1)%get_block(1,1))
          type is(preconditioner_t)
             call preconditioner_free(preconditioner_free_values,M%get_block(1,1)%get_block(1,1))
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond free (Mp^-1)
    select type (M)
    type is(block_preconditioner_t)
       select type(M%get_block(1,1))
       type is(block_preconditioner_t)
          select type(M%get_block(1,1)%get_block(2,2))
          type is(preconditioner_t)
             call preconditioner_free(preconditioner_free_values,M%get_block(1,1)%get_block(2,2))
          class default
             check(.false.)
          end select
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond free (Mx^-1)
    select type (M)
    type is(block_preconditioner_t)
       select type(M%get_block(2,2))
       type is(preconditioner_t)
          call preconditioner_free(preconditioner_free_values,M%get_block(2,2))
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select
          
  end subroutine free_preconditioner

  !==================================================================================================
  subroutine free_block_structures(block_operator,block_operand,block_preconditioner)
    implicit none
    type(block_operator_t)        , intent(inout) :: block_operator
    type(block_operand_t)         , intent(inout) :: block_operand
    type(block_preconditioner_l_t), intent(inout) :: block_preconditioner
    
    ! precond free (K^-1)
    select type(block_preconditioner%get_block(1,1))
    type is(block_preconditioner_t)
       select type(block_preconditioner%get_block(1,1)%get_block(1,1))
       type is(preconditioner_t)
          call preconditioner_free(preconditioner_free_struct, &
               &                   block_preconditioner%get_block(1,1)%get_block(1,1))
          call preconditioner_free(preconditioner_free_clean, &
               &                   block_preconditioner%get_block(1,1)%get_block(1,1))
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond free (Mp^-1)
    select type(block_preconditioner%get_block(1,1))
    type is(block_preconditioner_t)
       select type(block_preconditioner%get_block(1,1)%get_block(2,2))
       type is(preconditioner_t)
          call preconditioner_free(preconditioner_free_struct, &
               &                   block_preconditioner%get_block(1,1)%get_block(2,2))
          call preconditioner_free(preconditioner_free_clean, &
               &                   block_preconditioner%get_block(1,1)%get_block(2,2))
       class default
          check(.false.)
       end select
    class default
       check(.false.)
    end select

    ! precond free (Mx^-1)
    select type(block_preconditioner%get_block(2,2))
    type is(preconditioner_t)
       call preconditioner_free(preconditioner_free_struct,block_preconditioner%get_block(2,2))
       call preconditioner_free(preconditioner_free_clean,block_preconditioner%get_block(2,2))
    class default
       check(.false.)
    end select

    ! Free U-P Block operator
    select type(block_operator%get_block(1,1))
    type is(block_operator_t)
       call block_operator%get_block(1,1)%destroy()
    class default
       check(.false.)
    end select

    ! Free U-P Block operand
    call block_operand%blocks(1)%destroy()    

    ! Free UP-X Block operator
    select type(block_operator%get_block(1,2))
    type is(block_operator_t)
       call block_operator%get_block(1,2)%destroy()
    class default
       check(.false.)
    end select

    ! Free X-UP Block operator
    select type(block_operator%get_block(2,1))
    type is(block_operator_t)
       call block_operator%get_block(2,1)%destroy()
    class default
       check(.false.)
    end select

    ! Free Global operator
    call block_operator%destroy()

    ! Free Global operand
    call block_operand%destroy()
    
    ! Free Global preconditioner
    call block_preconditioner%destroy()
    
  end subroutine free_block_structures
  
  
end program test_blk_nsi_cg_iss_oss
