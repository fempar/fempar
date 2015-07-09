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
  
  type, extends(nonlinear_operator_t) :: my_linear_algebra_t 
     type(block_matrix_t)           :: block_matrix
     type(block_vector_t)           :: block_vector
     type(block_vector_t)           :: block_unknown
     type(block_matrix_t)           :: mass_p_matrix
     type(block_operator_t)         :: block_operator
     type(block_operand_t)          :: block_operand_vec
     type(block_operand_t)          :: block_operand_unk
     type(block_preconditioner_l_t) :: block_preconditioner
     type(preconditioner_t)         :: u_preconditioner
     type(preconditioner_t)         :: p_preconditioner
     type(preconditioner_t)         :: x_preconditioner
     type(preconditioner_params_t)  :: u_ppars
     type(preconditioner_params_t)  :: p_ppars
     type(preconditioner_params_t)  :: x_ppars
     type(block_operator_t)         :: block_up_operator
     type(block_operator_t)         :: block_up_x_operator
     type(block_operator_t)         :: block_x_up_operator
     type(block_operator_t)         :: block_x_operator
     type(block_operand_t)          :: block_up_operand_vec
     type(block_operand_t)          :: block_up_operand_unk
     type(block_operand_t)          :: block_x_operand_vec
     type(block_operand_t)          :: block_x_operand_unk
     type(block_preconditioner_l_t) :: block_up_preconditioner
     type(block_preconditioner_l_t) :: block_x_preconditioner
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
    class(my_linear_algebra_t), target, intent(inout) :: la
    type(block_graph_t)               , intent(in)    :: blk_graph

    ! Allocate matrices and vectors
    call la%block_matrix%alloc(blk_graph)
    call la%block_matrix%set_block_to_zero(2,2)
    call la%block_matrix%set_block_to_zero(2,3)
    call la%block_matrix%set_block_to_zero(3,2)
    call la%block_vector%alloc(blk_graph)
    call la%block_unknown%alloc(blk_graph)
    call la%block_vector%init(0.0_rp)

    ! Auxiliar matrices
    call la%mass_p_matrix%alloc(blk_graph)
    call la%mass_p_matrix%set_block_to_zero(1,1)
    call la%mass_p_matrix%set_block_to_zero(1,2)
    call la%mass_p_matrix%set_block_to_zero(1,3)
    call la%mass_p_matrix%set_block_to_zero(2,1)
    call la%mass_p_matrix%set_block_to_zero(2,3)
    call la%mass_p_matrix%set_block_to_zero(3,1)
    call la%mass_p_matrix%set_block_to_zero(3,2)
    call la%mass_p_matrix%set_block_to_zero(3,3)
    
    ! Construct U-preconditioner (K^-1)
    la%u_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%block_matrix%get_block(1,1),la%u_preconditioner,la%u_ppars)
    call preconditioner_symbolic(la%block_matrix%get_block(1,1),la%u_preconditioner)
    call preconditioner_log_info(la%u_preconditioner)

    ! Construct P-preconditioner (M^-1)
    la%p_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%mass_p_matrix%get_block(2,2),la%p_preconditioner,la%p_ppars)
    call preconditioner_symbolic(la%mass_p_matrix%get_block(2,2),la%p_preconditioner)
    call preconditioner_log_info(la%p_preconditioner)

    ! Construct X-preconditioner (Mx^-1)
    la%x_ppars%type = pardiso_mkl_prec
    call preconditioner_create(la%block_matrix%get_block(3,3),la%x_preconditioner,la%x_ppars)
    call preconditioner_symbolic(la%block_matrix%get_block(3,3),la%x_preconditioner)
    call preconditioner_log_info(la%x_preconditioner)
    
    ! Create U-P Block operator
    call la%block_up_operator%create(2,2)
    call la%block_up_operator%set_block(1,1,la%block_matrix%get_block(1,1))
    call la%block_up_operator%set_block(1,2,la%block_matrix%get_block(1,2))
    call la%block_up_operator%set_block(2,1,la%block_matrix%get_block(2,1))
    call la%block_up_operator%set_block_to_zero(2,2)

    ! Create U-P Block operand
    call la%block_up_operand_vec%create(2)
    call la%block_up_operand_vec%set_block(1,la%block_vector%blocks(1))
    call la%block_up_operand_vec%set_block(2,la%block_vector%blocks(2))
    call la%block_up_operand_unk%create(2)
    call la%block_up_operand_unk%set_block(1,la%block_unknown%blocks(1))
    call la%block_up_operand_unk%set_block(2,la%block_unknown%blocks(2))

    ! Create U-P Block preconditioner
    call la%block_up_preconditioner%create(2)
    call la%block_up_preconditioner%set_block(1,1,la%u_preconditioner)
    call la%block_up_preconditioner%set_block(2,1,la%block_matrix%get_block(2,1))
    call la%block_up_preconditioner%set_block(2,2,la%p_preconditioner)

    ! Create UP-X Block operator
    call la%block_up_x_operator%create(2,1)
    call la%block_up_x_operator%set_block(1,1,la%block_matrix%get_block(1,3))
    call la%block_up_x_operator%set_block_to_zero(2,1) 

    ! Create X-UP Block operator
    call la%block_x_up_operator%create(1,2)
    call la%block_x_up_operator%set_block(1,1,la%block_matrix%get_block(3,1))
    call la%block_x_up_operator%set_block_to_zero(1,2)

    ! Create X Block operator
    call la%block_x_operator%create(1,1)
    call la%block_x_operator%set_block(1,1,la%block_matrix%get_block(3,3))

    ! Create X Block operand
    call la%block_x_operand_vec%create(1)
    call la%block_x_operand_vec%set_block(1,la%block_vector%blocks(3))
    call la%block_x_operand_unk%create(1)
    call la%block_x_operand_unk%set_block(1,la%block_unknown%blocks(3))

    ! Create X Block preconditioner
    call la%block_x_preconditioner%create(1)
    call la%block_x_preconditioner%set_block(1,1,la%x_preconditioner)

    ! Create Global Block operator
    call la%block_operator%create(2,2)
    call la%block_operator%set_block(1,1,la%block_up_operator)
    call la%block_operator%set_block(1,2,la%block_up_x_operator)
    call la%block_operator%set_block(2,1,la%block_x_up_operator)
    call la%block_operator%set_block(2,2,la%block_x_operator)

    ! Create Global Block operand
    call la%block_operand_vec%create(2)
    call la%block_operand_vec%set_block(1,la%block_up_operand_vec)
    call la%block_operand_vec%set_block(2,la%block_x_operand_vec)
    call la%block_operand_unk%create(2)
    call la%block_operand_unk%set_block(1,la%block_up_operand_unk)
    call la%block_operand_unk%set_block(2,la%block_x_operand_unk)

    ! Create Global Block preconditioner
    call la%block_preconditioner%create(2)
    call la%block_preconditioner%set_block(1,1,la%block_up_preconditioner)
    call la%block_preconditioner%set_block(2,1,la%block_x_up_operator)
    call la%block_preconditioner%set_block(2,2,la%block_x_preconditioner)

    ! Assign operators
    la%A => la%block_operator
    la%M => la%block_preconditioner
    la%b => la%block_operand_vec
    la%x => la%block_operand_unk
    la%A_int => la%block_matrix
    la%b_int => la%block_vector
    la%x_sol => la%block_unknown
    
  end subroutine build_linear_algebra

  !==================================================================================================
  subroutine free_linear_algebra(la)
    implicit none
    class(my_linear_algebra_t), intent(inout) :: la

    ! Unassign operators
    la%A => null()
    la%M => null()
    la%b => null()
    la%x => null()
    la%A_int => null()
    la%b_int => null()
    la%x_sol => null()

    ! Destroy Global Block preconditioner
    call la%block_preconditioner%destroy()

    ! Destroy Global Block operand
    call la%block_operand_vec%destroy()
    call la%block_operand_unk%destroy()
    
    ! Destroy Global Block operator
    call la%block_operator%destroy()

    ! Destroy X Block preconditioner
    call la%block_x_preconditioner%destroy()

    ! Destroy X Block operand
    call la%block_x_operand_vec%destroy()
    call la%block_x_operand_unk%destroy()

    ! Destroy X Block operator
    call la%block_x_operator%destroy()

    ! Destroy X-UP Block operator
    call la%block_x_up_operator%destroy()

    ! Destroy UP-X Block operator
    call la%block_up_x_operator%destroy()

    ! Destroy U-P Block preconditioner
    call la%block_up_preconditioner%destroy()

    ! Destroy U-P Block operand
    call la%block_up_operand_vec%destroy()
    call la%block_up_operand_unk%destroy()

    ! Destroy U-P Block operator
    call la%block_up_operator%destroy()

    ! Destroy X-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%x_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%x_preconditioner)

    ! Destroy P-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%p_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%p_preconditioner)

    ! Destroy U-preconditioner
    call preconditioner_free(preconditioner_free_struct,la%u_preconditioner)
    call preconditioner_free(preconditioner_free_clean,la%u_preconditioner)

    ! Deallocate auxiliar matrices
    call la%mass_p_matrix%free()
          
    ! Deallocate problem matrix
    call la%block_matrix%free()

    ! Deallocate problem vectors
    call la%block_vector%free()
    call la%block_unknown%free()
    
  end subroutine free_linear_algebra

end module my_linear_algebra_names

program test_blk_nsi_cg_iss_oss
  use serial_names
  use my_linear_algebra_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use norm_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"

  ! Types
  type(uniform_mesh_descriptor_t)         :: gdata
  type(uniform_conditions_descriptor_t)   :: bdata
  type(reference_element_t)               :: geo_reference_element
  type(triangulation_t)                   :: f_trian
  type(conditions_t)                      :: f_cond
  type(dof_descriptor_t)                  :: dof_descriptor
  type(fe_space_t)                        :: fe_space  
  type(nsi_problem_t)                     :: myprob
  type(nsi_cg_iss_oss_discrete_t), target :: mydisc
  type(nsi_cg_iss_oss_matvec_t)  , target :: cg_iss_oss_matvec
  type(discrete_integration_pointer_t)    :: approx(1)
  type(preconditioner_params_t)           :: ppars
  type(solver_control_t)                  :: sctrl
  type(serial_environment_t)              :: senv
  type(vtk_t)                             :: fevtk
  class(base_operand_t)         , pointer :: x, b
  class(base_operator_t)        , pointer :: A, M
  type(block_graph_t)                     :: blk_graph
  type(scalar_t)                          :: enorm_u, enorm_p
  type(error_norm_t)             , target :: error_compute
  type(nsi_cg_iss_oss_massp_t)   , target :: mass_p_integration
  type(my_linear_algebra_t)               :: linear_algebra

  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat,i
  integer(ip) :: num_approximations = 1

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2

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

  ! Create problem
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call cg_iss_oss_matvec%create(myprob,mydisc)
  mydisc%dtinv      = 0.0_rp
  mydisc%kfl_proj   = 1
  mydisc%kfl_lump   = 1
  myprob%kfl_conv   = 0
  myprob%diffu      = 1.0_rp

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
  !call graph_print(6,blk_graph%get_block(3,3))
  !call fe_space_print(6,fe_space)

  ! Assign analytical solution
  if(gdata%ndime==2) then
     call fe_space%set_analytical_code((/4,5,3,0,0/),(/0,0,0,0,0/))
  else
     write(*,*) 'analytical function not ready for 3D'
  end if

  ! Create linear algebra structures
  call linear_algebra%build(blk_graph)
  linear_algebra%max_iter = 10
  linear_algebra%nltol    = 1.0e-12_rp

  ! Compute auxiliar matrix
  call mass_p_integration%create(myprob,mydisc)
  approx(1)%p => mass_p_integration
  call volume_integral(approx,fe_space,linear_algebra%mass_p_matrix)

  ! Apply boundary conditions to unkno
  call update_strong_dirichlet_bcond(fe_space,f_cond)
  call update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),0.0_rp,fe_space)

  ! Solver control parameters
  sctrl%method = rgmres
  sctrl%trace  = 100
  sctrl%rtol   = 1.0e-14_rp
  sctrl%track_conv_his = .false.

  ! Point discrete integration
  approx(1)%p => cg_iss_oss_matvec

  ! Do nonlinear iterations
  call nonlinear_iteration(sctrl,senv,approx,fe_space,linear_algebra)

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
  call cg_iss_oss_matvec%free
  call dof_descriptor_free(dof_descriptor)
  call triangulation_free(f_trian)
  call conditions_free(f_cond)
  call finite_element_fixed_info_free(geo_reference_element)
  call uniform_conditions_descriptor_free(bdata)

  call memstatus

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
  subroutine compute_preconditioner(la)
    implicit none
    type(my_linear_algebra_t), intent(inout) :: la

    ! precond numeric (K^-1)
    call preconditioner_numeric(la%block_matrix%get_block(1,1),la%u_preconditioner)

    ! precond numeric (Mp^-1)
    call preconditioner_numeric(la%mass_p_matrix%get_block(2,2),la%p_preconditioner)

    ! precond numeric (Mx^-1)
    call preconditioner_numeric(la%block_matrix%get_block(3,3),la%x_preconditioner)
          
  end subroutine compute_preconditioner
  
  !==================================================================================================
  subroutine free_preconditioner(la)
    implicit none
    type(my_linear_algebra_t), intent(inout) :: la

    ! precond free (K^-1)
    call preconditioner_free(preconditioner_free_values,la%u_preconditioner)

    ! precond free (Mp^-1)
    call preconditioner_free(preconditioner_free_values,la%p_preconditioner)

    ! precond free (Mx^-1)
    call preconditioner_free(preconditioner_free_values,la%x_preconditioner)
          
  end subroutine free_preconditioner

  !==================================================================================================
  subroutine nonlinear_iteration( sctrl, env, approx, fe_space, nlop )
    implicit none
    type(solver_control_t)              , intent(inout) :: sctrl
    class(abstract_environment_t)       , intent(in)    :: env
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    !class(nonlinear_operator_t), target , intent(inout) :: nlop
    type(my_linear_algebra_t), target , intent(inout) :: nlop
    ! Locals
    ! Locals
    integer(ip)           :: iiter
    real(rp)              :: resnorm,ininorm
    type(block_operand_t) :: y

    ! Checks
    check(associated(nlop%A))
    check(associated(nlop%M))
    check(associated(nlop%b))
    check(associated(nlop%x))
    check(associated(nlop%A_int))
    check(associated(nlop%b_int))
    check(associated(nlop%x_sol))
        
    iiter = 0
    do while( iiter < nlop%max_iter )

       ! Update counter
       iiter = iiter+1

       ! Initialize Matrix and vector
       ! ***************** Abstract procedure to initialize a base_operator ************************!
       ! call nlop%A_int%init()
       call block_matrix_zero(nlop%block_matrix)
       !********************************************************************************************!
       call nlop%b_int%init(0.0_rp)

       ! Integrate system
       call volume_integral(approx,fe_space,nlop%A_int,nlop%b_int)

       ! Check convergence
       if(iiter==1) ininorm = nlop%b%nrm2()   
       y = nlop%b - nlop%A*nlop%x
       resnorm = y%nrm2()
       if( resnorm < nlop%nltol*ininorm) then
          write(*,*) 'Nonlinear iterations: ', iiter
          write(*,*) 'Nonlinear error norm: ', resnorm
          exit
       end if

       ! Compute Numeric preconditioner
       ! ***************** Abstract procedure to compute precond numeric ***************************!
       ! call nlop%M%compute()
       call compute_preconditioner(nlop)
       !********************************************************************************************!

       ! Solve system
       call abstract_solve(nlop%A,nlop%M,nlop%b,nlop%x,sctrl,env)
       call solver_control_log_conv_his(sctrl)
       call solver_control_free_conv_his(sctrl)

       ! Free Numeric preconditioner
       ! ******************** Abstract procedure to free precond numeric ***************************!
       ! call nlop%M%free_values()
       call free_preconditioner(nlop)
       !********************************************************************************************!
       
       ! Store solution to unkno
       call update_solution(nlop%x_sol,fe_space)
       
       ! Store nonlinear iteration ( k+1 --> k )
       call update_nonlinear(fe_space)
       
    end do

    ! Deallocate
    call y%free()

  end subroutine nonlinear_iteration
  
end program test_blk_nsi_cg_iss_oss


