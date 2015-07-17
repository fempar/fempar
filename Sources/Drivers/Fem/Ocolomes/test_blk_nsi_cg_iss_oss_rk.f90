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
module my_nonlinear_operator_names
  use serial_names
  implicit none
# include "debug.i90"
  private
  
  type, extends(picard_nonlinear_operator_t) :: momentum_operator_t 
     type(block_matrix_t)              :: block_matrix_l
     type(block_matrix_t)              :: block_matrix_nl
     type(block_matrix_t)              :: block_matrix_m
     type(block_matrix_t)              :: block_u_matrix
     type(block_matrix_t)              :: mass_x_matrix
     type(block_vector_t)              :: block_vector_l
     type(block_vector_t)              :: block_vector_nl
     type(block_vector_t)              :: block_vector_m
     type(block_vector_t)              :: block_vector
     type(block_vector_t)              :: block_unknown
     type(block_operator_t)            :: Am, Al, Anl
     type(abs_operator_t)              :: A_uu
     type(block_operator_t)            :: A_momentum
     type(block_operand_t)             :: b_momentum, bm, bl, bnl, brhs
     type(block_operand_t)             :: x_momentum
     type(block_preconditioner_l_t)    :: M_momentum
     type(preconditioner_t)            :: u_preconditioner
     type(preconditioner_t)            :: x_preconditioner
     type(preconditioner_params_t)     :: u_ppars
     type(preconditioner_params_t)     :: x_ppars
     type(inverse_operator_t)          :: inv_A_uu
     type(solver_control_t)            :: sctrl_A_uu
     real(rp)                          :: dt, aii
   contains
     procedure :: build => build_momentum_operator
     procedure :: free  => free_momentum_operator
     procedure :: compute_vector => compute_vector_momentum_operator
  end type momentum_operator_t
  
  type, extends(picard_nonlinear_operator_t) :: pressure_operator_t 
     type(block_matrix_t)           :: mass_u_matrix
     type(block_matrix_t)           :: lapla_p_matrix
     type(block_matrix_t)           :: block_matrix
     type(block_vector_t)           :: block_vector
     type(block_vector_t)           :: block_unknown
     type(block_operator_t)         :: A_pressure
     type(block_operand_t)          :: b_pressure
     type(block_operand_t)          :: x_pressure
     type(block_preconditioner_l_t) :: M_pressure
     type(preconditioner_t)         :: w_preconditioner
     type(preconditioner_t)         :: p_preconditioner
     type(preconditioner_params_t)  :: w_ppars
     type(preconditioner_params_t)  :: p_ppars     
   contains
     procedure :: build => build_pressure_operator
     procedure :: free  => free_pressure_operator
     procedure :: compute_vector => compute_vector_pressure_operator
  end type pressure_operator_t
  
  type, extends(picard_nonlinear_operator_t) :: momentum_update_operator_t 
     type(block_matrix_t)           :: mass_u_matrix
     type(block_matrix_t)           :: block_matrix
     type(block_vector_t)           :: block_vector
     type(block_vector_t)           :: block_unknown
     type(block_operator_t)         :: A_update_u
     type(block_operand_t)          :: b_update_u
     type(block_operand_t)          :: x_update_u
     type(block_preconditioner_l_t) :: M_update_u
     type(preconditioner_t)         :: u_preconditioner
     type(preconditioner_params_t)  :: u_ppars
   contains
     procedure :: build => build_momentum_update_operator
     procedure :: free  => free_momentum_update_operator
     procedure :: compute_vector => compute_vector_momentum_update_operator
  end type momentum_update_operator_t
  
  type, extends(picard_nonlinear_operator_t) :: projection_update_operator_t 
     type(block_matrix_t)           :: mass_x_matrix
     type(block_matrix_t)           :: block_matrix
     type(block_vector_t)           :: block_vector
     type(block_vector_t)           :: block_unknown
     type(block_operator_t)         :: A_update_x
     type(block_operand_t)          :: b_update_x
     type(block_operand_t)          :: x_update_x
     type(block_preconditioner_l_t) :: M_update_x
     type(preconditioner_t)         :: x_preconditioner
     type(preconditioner_params_t)  :: x_ppars
   contains
     procedure :: build => build_projection_update_operator
     procedure :: free  => free_projection_update_operator
     procedure :: compute_vector => compute_vector_projection_update_operator
  end type projection_update_operator_t

 ! Types
  public :: momentum_operator_t, pressure_operator_t, momentum_update_operator_t, &
       &    projection_update_operator_t
  
contains

  !==================================================================================================
  subroutine build_momentum_operator(nlop,blk_graph,env)
    implicit none
    class(momentum_operator_t) , intent(inout) :: nlop
    type(block_graph_t)        , intent(in)    :: blk_graph
    class(serial_environment_t), intent(in)    :: env

    ! Modify default nonlinear parameters
    nlop%max_iter = 20

    ! Modify default inverse_operator solver parameters
    nlop%sctrl_A_uu%rtol = 1.0e-14_rp
    
    ! Initialize parameters
    nlop%dt  = 0.0_rp
    nlop%aii = 0.0_rp

    ! Allocate matrices and vectors
    ! Mass matrix
    call nlop%block_matrix_m%alloc(blk_graph)
    call nlop%block_matrix_m%set_block_to_zero(1,2)
    call nlop%block_matrix_m%set_block_to_zero(1,3)
    call nlop%block_matrix_m%set_block_to_zero(2,1)
    call nlop%block_matrix_m%set_block_to_zero(2,2)
    call nlop%block_matrix_m%set_block_to_zero(2,3)
    call nlop%block_matrix_m%set_block_to_zero(3,1)
    call nlop%block_matrix_m%set_block_to_zero(3,2)
    call nlop%block_matrix_m%set_block_to_zero(3,3)
    nlop%block_matrix_m%fill_values_stage = update_transient
    ! Linear matrix
    call nlop%block_matrix_l%alloc(blk_graph)
    call nlop%block_matrix_l%set_block_to_zero(1,2)
    call nlop%block_matrix_l%set_block_to_zero(1,3)
    call nlop%block_matrix_l%set_block_to_zero(2,1)
    call nlop%block_matrix_l%set_block_to_zero(2,2)
    call nlop%block_matrix_l%set_block_to_zero(2,3)
    call nlop%block_matrix_l%set_block_to_zero(3,1)
    call nlop%block_matrix_l%set_block_to_zero(3,2)
    call nlop%block_matrix_l%set_block_to_zero(3,3)
    nlop%block_matrix_l%fill_values_stage = update_transient
    ! Nonlinear matrix
    call nlop%block_matrix_nl%alloc(blk_graph)
    call nlop%block_matrix_nl%set_block_to_zero(1,2)
    call nlop%block_matrix_nl%set_block_to_zero(2,1)
    call nlop%block_matrix_nl%set_block_to_zero(2,2)
    call nlop%block_matrix_nl%set_block_to_zero(2,3)
    call nlop%block_matrix_nl%set_block_to_zero(3,2)
    nlop%block_matrix_nl%fill_values_stage = update_nonlinear
    ! Auxiliar matrix (Block-U preconditioner)
    call nlop%block_u_matrix%alloc(blk_graph)
    call nlop%block_u_matrix%set_block_to_zero(1,2)
    call nlop%block_u_matrix%set_block_to_zero(1,3)
    call nlop%block_u_matrix%set_block_to_zero(2,1)
    call nlop%block_u_matrix%set_block_to_zero(2,2)
    call nlop%block_u_matrix%set_block_to_zero(2,3)
    call nlop%block_u_matrix%set_block_to_zero(3,1)
    call nlop%block_u_matrix%set_block_to_zero(3,2)
    call nlop%block_u_matrix%set_block_to_zero(3,3)
    nlop%block_u_matrix%fill_values_stage = update_nonlinear
    ! Auxiliar matrix (Block-X preconditioner)
    call nlop%mass_x_matrix%alloc(blk_graph)
    call nlop%mass_x_matrix%set_block_to_zero(1,1)
    call nlop%mass_x_matrix%set_block_to_zero(1,2)
    call nlop%mass_x_matrix%set_block_to_zero(1,3)
    call nlop%mass_x_matrix%set_block_to_zero(2,1)
    call nlop%mass_x_matrix%set_block_to_zero(2,2)
    call nlop%mass_x_matrix%set_block_to_zero(2,3)
    call nlop%mass_x_matrix%set_block_to_zero(3,1)
    call nlop%mass_x_matrix%set_block_to_zero(3,2)
    nlop%block_matrix_m%fill_values_stage = update_transient
    ! Mass BCs vector
    call nlop%block_vector_m%alloc(blk_graph)
    ! Linear BCs vector
    call nlop%block_vector_l%alloc(blk_graph)
    ! Nonlinear BCs vector
    call nlop%block_vector_nl%alloc(blk_graph)
    ! RHS vector
    call nlop%block_vector%alloc(blk_graph)
    call nlop%block_vector%init(0.0_rp)
    ! Unknown vector
    call nlop%block_unknown%alloc(blk_graph)
    
    ! Construct U-preconditioner (K^-1)
    nlop%u_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%block_u_matrix%get_block(1,1),nlop%u_preconditioner,nlop%u_ppars)
    call preconditioner_symbolic(nlop%block_u_matrix%get_block(1,1),nlop%u_preconditioner)
    call preconditioner_log_info(nlop%u_preconditioner)

    ! Construct X-preconditioner (Mx^-1)
    nlop%x_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%mass_x_matrix%get_block(3,3),nlop%x_preconditioner,nlop%x_ppars)
    call preconditioner_symbolic(nlop%mass_x_matrix%get_block(3,3),nlop%x_preconditioner)
    call preconditioner_log_info(nlop%x_preconditioner)

    ! Create Mass operator
    call nlop%Am%create(2,2)
    call nlop%Am%set_block(1,1,nlop%block_matrix_m%get_block(1,1))
    call nlop%Am%set_block_to_zero(1,2)
    call nlop%Am%set_block_to_zero(2,1)
    call nlop%Am%set_block_to_zero(2,2)

    ! Create Linear operator
    call nlop%Al%create(2,2)
    call nlop%Al%set_block(1,1,nlop%block_matrix_l%get_block(1,1))
    call nlop%Al%set_block_to_zero(1,2)
    call nlop%Al%set_block_to_zero(2,1)
    call nlop%Al%set_block_to_zero(2,2)

    ! Create Nonlinear operator
    call nlop%Anl%create(2,2)
    call nlop%Anl%set_block(1,1,nlop%block_matrix_nl%get_block(1,1))
    call nlop%Anl%set_block(1,2,nlop%block_matrix_nl%get_block(1,3))
    call nlop%Anl%set_block(2,1,nlop%block_matrix_nl%get_block(3,1))
    call nlop%Anl%set_block(2,2,nlop%block_matrix_nl%get_block(3,3))

    ! Construct abstract operator for UU block
    nlop%A_uu =  nlop%dt * nlop%Am%get_block(1,1)  +  nlop%aii * nlop%Al%get_block(1,1) + &
         &       nlop%aii * nlop%Anl%get_block(1,1)

    ! Create inverse operator for Block UU
    call nlop%inv_A_uu%create(nlop%A_uu,nlop%u_preconditioner,nlop%sctrl_A_uu,env)
    
    ! Construct global operator
    call nlop%A_momentum%create(2,2)
    call nlop%A_momentum%set_block(1,1,nlop%A_uu)
    call nlop%A_momentum%set_block(1,2,nlop%aii*nlop%block_matrix_nl%get_block(1,3))
    call nlop%A_momentum%set_block(2,1,nlop%aii*nlop%block_matrix_nl%get_block(3,1))
    call nlop%A_momentum%set_block(2,2,nlop%aii*nlop%block_matrix_nl%get_block(3,3))  
    
    ! Create Mass BCs operand
    call nlop%bm%create(2)
    call nlop%bm%set_block(1,nlop%block_vector_m%blocks(1))
    call nlop%bm%set_block(2,nlop%block_vector_m%blocks(3))

    ! Create Linear BCs operand
    call nlop%bl%create(2)
    call nlop%bl%set_block(1,nlop%block_vector_l%blocks(1))
    call nlop%bl%set_block(2,nlop%block_vector_l%blocks(3))

    ! Create Nonlinear BCs operand
    call nlop%bnl%create(2)
    call nlop%bnl%set_block(1,nlop%block_vector_nl%blocks(1))
    call nlop%bnl%set_block(2,nlop%block_vector_nl%blocks(3))

    ! Create RHS operand
    call nlop%brhs%create(2)
    call nlop%brhs%set_block(1,nlop%block_vector%blocks(1))
    call nlop%brhs%set_block(2,nlop%block_vector%blocks(3))

    ! Construct solution operand
    call nlop%x_momentum%create(2)
    call nlop%x_momentum%set_block(1,nlop%block_unknown%blocks(1))
    call nlop%x_momentum%set_block(2,nlop%block_unknown%blocks(3))

    ! Construct global block preconditioner
    call nlop%M_momentum%create(2)
    call nlop%M_momentum%set_block(1,1,nlop%inv_A_uu)
    call nlop%M_momentum%set_block(2,1,nlop%aii*nlop%Anl%get_block(2,1))
    call nlop%M_momentum%set_block(2,2,nlop%x_preconditioner)

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_momentum, nlop%M_momentum, nlop%b_momentum, nlop%x_momentum, &
         &           nlop%block_unknown,                                                 &
         &           A_int_c = nlop%block_matrix_l, b_int_c = nlop%block_vector_l,     &
         &           A_int_t = nlop%block_matrix_m, b_int_t = nlop%block_vector_m,     &
         &           A_int_n = nlop%block_matrix_nl,b_int_n = nlop%block_vector_nl)
    
  end subroutine build_momentum_operator

  !==================================================================================================
  subroutine compute_vector_momentum_operator(nlop)
    implicit none
    class(momentum_operator_t), intent(inout) :: nlop
    
    ! Construct global operand
    nlop%b_momentum = nlop%dt * nlop%bm + nlop%aii * nlop%bl + nlop%aii * nlop%bnl + nlop%brhs
    
  end subroutine compute_vector_momentum_operator

  !==================================================================================================
  subroutine build_pressure_operator(nlop,blk_graph)
    implicit none
    class(pressure_operator_t), intent(inout) :: nlop
    type(block_graph_t)       , intent(in)    :: blk_graph

    ! Allocate matrices and vectors
    call nlop%block_matrix%alloc(blk_graph)
    call nlop%block_matrix%set_block_to_zero(1,3)
    call nlop%block_matrix%set_block_to_zero(2,2)
    call nlop%block_matrix%set_block_to_zero(2,3)
    call nlop%block_matrix%set_block_to_zero(3,1)
    call nlop%block_matrix%set_block_to_zero(3,2)
    call nlop%block_matrix%set_block_to_zero(3,3)
    nlop%block_matrix%fill_values_stage = update_transient
    call nlop%block_vector%alloc(blk_graph)
    call nlop%block_unknown%alloc(blk_graph)
    call nlop%block_vector%init(0.0_rp)

    ! Auxiliar matrices
    call nlop%mass_u_matrix%alloc(blk_graph)
    call nlop%mass_u_matrix%set_block_to_zero(1,2)
    call nlop%mass_u_matrix%set_block_to_zero(1,3)
    call nlop%mass_u_matrix%set_block_to_zero(2,1)
    call nlop%mass_u_matrix%set_block_to_zero(2,2)
    call nlop%mass_u_matrix%set_block_to_zero(2,3)
    call nlop%mass_u_matrix%set_block_to_zero(3,1)
    call nlop%mass_u_matrix%set_block_to_zero(3,2)
    call nlop%mass_u_matrix%set_block_to_zero(3,3)
    call nlop%lapla_p_matrix%alloc(blk_graph)
    call nlop%lapla_p_matrix%set_block_to_zero(1,1)
    call nlop%lapla_p_matrix%set_block_to_zero(1,2)
    call nlop%lapla_p_matrix%set_block_to_zero(1,3)
    call nlop%lapla_p_matrix%set_block_to_zero(2,1)
    call nlop%lapla_p_matrix%set_block_to_zero(2,3)
    call nlop%lapla_p_matrix%set_block_to_zero(3,1)
    call nlop%lapla_p_matrix%set_block_to_zero(3,2)
    call nlop%lapla_p_matrix%set_block_to_zero(3,3)

    ! Construct W-preconditioner (Mu^-1)
    nlop%w_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%mass_u_matrix%get_block(1,1),nlop%w_preconditioner,nlop%w_ppars)
    call preconditioner_symbolic(nlop%mass_u_matrix%get_block(1,1),nlop%w_preconditioner)
    call preconditioner_log_info(nlop%w_preconditioner)
    nlop%w_preconditioner%fill_values_stage = update_transient

    ! Construct P-preconditioner (Lp^-1)
    nlop%p_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%lapla_p_matrix%get_block(2,2),nlop%p_preconditioner,nlop%p_ppars)
    call preconditioner_symbolic(nlop%lapla_p_matrix%get_block(2,2),nlop%p_preconditioner)
    call preconditioner_log_info(nlop%p_preconditioner)
    nlop%p_preconditioner%fill_values_stage = update_transient
    
    ! Create W-P Block operator
    call nlop%A_pressure%create(2,2)
    call nlop%A_pressure%set_block(1,1,nlop%block_matrix%get_block(1,1))
    call nlop%A_pressure%set_block(1,2,nlop%block_matrix%get_block(1,2))
    call nlop%A_pressure%set_block(2,1,nlop%block_matrix%get_block(2,1))
    call nlop%A_pressure%set_block_to_zero(2,2)

    ! Create W-P Block operand
    call nlop%b_pressure%create(2)
    call nlop%b_pressure%set_block(1,nlop%block_vector%blocks(1))
    call nlop%b_pressure%set_block(2,nlop%block_vector%blocks(2))
    call nlop%x_pressure%create(2)
    call nlop%x_pressure%set_block(1,nlop%block_unknown%blocks(1))
    call nlop%x_pressure%set_block(2,nlop%block_unknown%blocks(2))

    ! Create W-P Block preconditioner
    call nlop%M_pressure%create(2)
    call nlop%M_pressure%set_block(1,1,nlop%w_preconditioner)
    call nlop%M_pressure%set_block(2,1,nlop%block_matrix%get_block(2,1))
    call nlop%M_pressure%set_block(2,2,nlop%p_preconditioner)
    nlop%M_pressure%fill_values_stage = update_transient

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_pressure, nlop%M_pressure, nlop%b_pressure, nlop%x_pressure, &
         &           nlop%block_unknown, A_int_t=nlop%block_matrix, b_int_t=nlop%block_vector)
    
  end subroutine build_pressure_operator

  !==================================================================================================
  subroutine compute_vector_pressure_operator(nlop)
    implicit none
    class(pressure_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_pressure_operator

  !==================================================================================================
  subroutine build_momentum_update_operator(nlop,blk_graph)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop
    type(block_graph_t)              , intent(in)    :: blk_graph

    ! Allocate integration matrices and vectors
    call nlop%block_matrix%alloc(blk_graph)
    call nlop%block_matrix%set_block_to_zero(1,2)
    call nlop%block_matrix%set_block_to_zero(1,3)
    call nlop%block_matrix%set_block_to_zero(2,1)
    call nlop%block_matrix%set_block_to_zero(2,2)
    call nlop%block_matrix%set_block_to_zero(2,3)
    call nlop%block_matrix%set_block_to_zero(3,1)
    call nlop%block_matrix%set_block_to_zero(3,2)
    call nlop%block_matrix%set_block_to_zero(3,3)
    nlop%block_matrix%fill_values_stage = update_transient
    call nlop%block_vector%alloc(blk_graph)
    call nlop%block_unknown%alloc(blk_graph)

    ! Auxiliar matrices
    call nlop%mass_u_matrix%alloc(blk_graph)
    call nlop%mass_u_matrix%set_block_to_zero(1,2)
    call nlop%mass_u_matrix%set_block_to_zero(1,3)
    call nlop%mass_u_matrix%set_block_to_zero(2,1)
    call nlop%mass_u_matrix%set_block_to_zero(2,2)
    call nlop%mass_u_matrix%set_block_to_zero(2,3)
    call nlop%mass_u_matrix%set_block_to_zero(3,1)
    call nlop%mass_u_matrix%set_block_to_zero(3,2)
    call nlop%mass_u_matrix%set_block_to_zero(3,3)
    
    ! Construct U-preconditioner (Mu^-1)
    nlop%u_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%mass_u_matrix%get_block(1,1),nlop%u_preconditioner,nlop%u_ppars)
    call preconditioner_symbolic(nlop%mass_u_matrix%get_block(1,1),nlop%u_preconditioner)
    call preconditioner_log_info(nlop%u_preconditioner)
    nlop%u_preconditioner%fill_values_stage = update_transient
    
    ! Create U Block operator
    call nlop%A_update_u%create(1,1)
    call nlop%A_update_u%set_block(1,1,nlop%block_matrix%get_block(1,1))

    ! Create U Block operand
    call nlop%b_update_u%create(1)
    call nlop%b_update_u%set_block(1,nlop%block_vector%blocks(1))
    call nlop%x_update_u%create(1)
    call nlop%x_update_u%set_block(1,nlop%block_unknown%blocks(1))

    ! Create U Block preconditioner
    call nlop%M_update_u%create(1)
    call nlop%M_update_u%set_block(1,1,nlop%u_preconditioner)
    nlop%M_update_u%fill_values_stage = update_transient

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_update_u, nlop%M_update_u, nlop%b_update_u, nlop%x_update_u, &
         &           nlop%block_unknown, A_int_t=nlop%block_matrix, b_int_t=nlop%block_vector)
    
  end subroutine build_momentum_update_operator

  !==================================================================================================
  subroutine compute_vector_momentum_update_operator(nlop)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_momentum_update_operator

 !==================================================================================================
  subroutine build_projection_update_operator(nlop,blk_graph)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop
    type(block_graph_t)                , intent(in)    :: blk_graph

    ! Allocate matrices and vectors
    call nlop%block_matrix%alloc(blk_graph)
    call nlop%block_matrix%set_block_to_zero(1,1)
    call nlop%block_matrix%set_block_to_zero(1,2)
    call nlop%block_matrix%set_block_to_zero(1,3)
    call nlop%block_matrix%set_block_to_zero(2,1)
    call nlop%block_matrix%set_block_to_zero(2,2)
    call nlop%block_matrix%set_block_to_zero(2,3)
    call nlop%block_matrix%set_block_to_zero(3,1)
    call nlop%block_matrix%set_block_to_zero(3,2)
    nlop%block_matrix%fill_values_stage = update_transient
    call nlop%block_vector%alloc(blk_graph)
    call nlop%block_unknown%alloc(blk_graph)

    ! Auxiliar matrices
    call nlop%mass_x_matrix%alloc(blk_graph)
    call nlop%mass_x_matrix%set_block_to_zero(1,1)
    call nlop%mass_x_matrix%set_block_to_zero(1,2)
    call nlop%mass_x_matrix%set_block_to_zero(1,3)
    call nlop%mass_x_matrix%set_block_to_zero(2,1)
    call nlop%mass_x_matrix%set_block_to_zero(2,2)
    call nlop%mass_x_matrix%set_block_to_zero(2,3)
    call nlop%mass_x_matrix%set_block_to_zero(3,1)
    call nlop%mass_x_matrix%set_block_to_zero(3,2)

    ! Construct X-preconditioner (Mx^-1)
    nlop%x_ppars%type = pardiso_mkl_prec
    call preconditioner_create(nlop%mass_x_matrix%get_block(3,3),nlop%x_preconditioner,nlop%x_ppars)
    call preconditioner_symbolic(nlop%mass_x_matrix%get_block(3,3),nlop%x_preconditioner)
    call preconditioner_log_info(nlop%x_preconditioner)
    nlop%x_preconditioner%fill_values_stage = update_transient
    
    ! Create U Block operator
    call nlop%A_update_x%create(1,1)
    call nlop%A_update_x%set_block(1,1,nlop%block_matrix%get_block(3,3))

    ! Create U Block operand
    call nlop%b_update_x%create(1)
    call nlop%b_update_x%set_block(1,nlop%block_vector%blocks(3))
    call nlop%x_update_x%create(1)
    call nlop%x_update_x%set_block(1,nlop%block_unknown%blocks(3))

    ! Create U Block preconditioner
    call nlop%M_update_x%create(1)
    call nlop%M_update_x%set_block(1,1,nlop%x_preconditioner)
    nlop%M_update_x%fill_values_stage = update_transient

    ! Fill picard nonlinear operator pointers
    call nlop%create(nlop%A_update_x, nlop%M_update_x, nlop%b_update_x, nlop%x_update_x, &
         &           nlop%block_unknown, A_int_t=nlop%block_matrix, b_int_t=nlop%block_vector)
    
  end subroutine build_projection_update_operator

  !==================================================================================================
  subroutine compute_vector_projection_update_operator(nlop)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop
    
    ! Dummy
    
  end subroutine compute_vector_projection_update_operator

  !==================================================================================================
  subroutine free_momentum_operator(nlop)
    implicit none
    class(momentum_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_momentum%destroy()

    ! Destroy global Block operand
    call nlop%b_momentum%destroy()
    call nlop%x_momentum%destroy()

    ! Destroy intermediate operands
    call nlop%bm%destroy()
    call nlop%bl%destroy()
    call nlop%bnl%destroy()
    call nlop%brhs%destroy()

    ! Destroy global Block operator
    call nlop%A_momentum%free()

    ! Destroy intermediate operators
    call nlop%Am%destroy()
    call nlop%Al%destroy()
    call nlop%Anl%destroy()

    ! Destroy X-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%x_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%x_preconditioner)

    ! Destroy U-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%u_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%u_preconditioner)

    ! Deallocate vectors
    call nlop%block_unknown%free()
    call nlop%block_vector%free()
    call nlop%block_vector_m%free()
    call nlop%block_vector_l%free()
    call nlop%block_vector_nl%free()
    
    ! Deallocate matrices
    call nlop%mass_x_matrix%free()
    call nlop%block_u_matrix%free()
    call nlop%block_matrix_m%free()
    call nlop%block_matrix_l%free()
    call nlop%block_matrix_nl%free()

  end subroutine free_momentum_operator

  !==================================================================================================
  subroutine free_pressure_operator(nlop)
    implicit none
    class(pressure_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_pressure%destroy()

    ! Destroy global Block operand
    call nlop%b_pressure%destroy()
    call nlop%x_pressure%destroy()

    ! Destroy global Block operator
    call nlop%A_pressure%destroy()

    ! Destroy P-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%p_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%p_preconditioner)

    ! Destroy W-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%w_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%w_preconditioner)

    ! Deallocate vectors
    call nlop%block_unknown%free()
    call nlop%block_vector%free()
    
    ! Deallocate matrices
    call nlop%mass_u_matrix%free()
    call nlop%lapla_p_matrix%free()
    call nlop%block_matrix%free()

  end subroutine free_pressure_operator

  !==================================================================================================
  subroutine free_momentum_update_operator(nlop)
    implicit none
    class(momentum_update_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_update_u%destroy()

    ! Destroy global Block operand
    call nlop%b_update_u%destroy()
    call nlop%x_update_u%destroy()

    ! Destroy global Block operator
    call nlop%A_update_u%destroy()

    ! Destroy U-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%u_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%u_preconditioner)

    ! Deallocate vectors
    call nlop%block_unknown%free()
    call nlop%block_vector%free()
    
    ! Deallocate matrices
    call nlop%mass_u_matrix%free()
    call nlop%block_matrix%free()

  end subroutine free_momentum_update_operator

  !==================================================================================================
  subroutine free_projection_update_operator(nlop)
    implicit none
    class(projection_update_operator_t), intent(inout) :: nlop

    ! Unassign picard nonlinear operator
    call nlop%unassign()

    ! Destroy global Block preconditioner
    call nlop%M_update_x%destroy()

    ! Destroy global Block operand
    call nlop%b_update_x%destroy()
    call nlop%x_update_x%destroy()

    ! Destroy global Block operator
    call nlop%A_update_x%destroy()

    ! Destroy X-preconditioner
    call preconditioner_free(preconditioner_free_struct,nlop%x_preconditioner)
    call preconditioner_free(preconditioner_free_clean,nlop%x_preconditioner)

    ! Deallocate vectors
    call nlop%block_unknown%free()
    call nlop%block_vector%free()
    
    ! Deallocate matrices
    call nlop%mass_x_matrix%free()
    call nlop%block_matrix%free()

  end subroutine free_projection_update_operator

end module my_nonlinear_operator_names

program test_blk_nsi_cg_iss_oss_rk
  use serial_names
  use my_nonlinear_operator_names
  use nsi_names
  use nsi_cg_iss_oss_names
  use norm_names
  use lib_vtk_io_interface_names
  implicit none
# include "debug.i90"

  ! Geometry & Integration Types
  type(uniform_mesh_descriptor_t)       :: gdata
  type(uniform_conditions_descriptor_t) :: bdata
  type(reference_element_t)             :: geo_reference_element
  type(triangulation_t)                 :: f_trian
  type(conditions_t)                    :: f_cond
  type(dof_descriptor_t)                :: dof_descriptor
  type(fe_space_t)                      :: fe_space  
  type(block_graph_t)                   :: blk_graph
  type(scalar_t)                        :: enorm_u, enorm_p

  ! Problem types
  type(nsi_problem_t)                                 :: myprob
  type(nsi_cg_iss_oss_discrete_t)                     :: mydisc
  type(nsi_cg_iss_oss_rk_momentum_t)         , target :: cg_iss_oss_rk_momentum
  type(nsi_cg_iss_oss_rk_momentum_rhs_t)     , target :: cg_iss_oss_rk_momentum_rhs
  type(nsi_cg_iss_oss_rk_pressure_t)         , target :: cg_iss_oss_rk_pressure
  type(nsi_cg_iss_oss_rk_momentum_update_t)  , target :: cg_iss_oss_rk_momentum_update
  type(nsi_cg_iss_oss_rk_projection_update_t), target :: cg_iss_oss_rk_projection_update
  type(nsi_cg_iss_oss_lapla_p_t)             , target :: lapla_p_integration
  type(nsi_cg_iss_oss_massu_t)               , target :: mass_u_integration
  type(nsi_cg_iss_oss_massx_t)               , target :: mass_x_integration
  type(error_norm_t)                         , target :: error_compute
  type(rungekutta_integrator_t)              , target :: rkinteg
  type(time_integration_t)                   , target :: tinteg
  type(discrete_integration_pointer_t)                :: approx(1)
  type(momentum_operator_t)                           :: momentum_operator
  type(pressure_operator_t)                           :: pressure_operator
  type(momentum_update_operator_t)                    :: momentum_update_operator
  type(projection_update_operator_t)                  :: projection_update_operator

  ! Solver types
  type(preconditioner_params_t) :: ppars
  type(solver_control_t)        :: sctrl
  type(serial_environment_t)    :: senv

  ! Postproces types
  type(vtk_t) :: fevtk

  ! Logicals
  logical :: ginfo_state

  ! Integers
  integer(ip) :: gtype(3) = (/ csr, csr, csr /)
  integer(ip) :: ibloc,jbloc,istat,i
  integer(ip) :: num_approximations = 1
  integer(ip) :: setterms(6,2),settable(3)

  ! Parameters
  integer(ip), parameter :: velocity=1, pressure=2
  integer(ip), parameter :: nonlinear=0, linear=1
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
  call reference_element_create(geo_reference_element,Q_type_id,1,gdata%ndime)

  ! Generate triangulation
  call generate_uniform_triangulation(1,gdata,bdata,geo_reference_element,f_trian,f_cond,material)

  ! Define Runge-Kutta method
  settable      = (/nstage,rk_order,rk_flag/)
  setterms(1,:) = (/update_transient,implicit/)  ! Diffusion
  setterms(2,:) = (/update_nonlinear,implicit/)  ! Convection
  setterms(3,:) = (/update_transient,explicit/)  ! Pressure Gradient
  setterms(4,:) = (/update_nonlinear,implicit/)  ! OSS_vu
  setterms(5,:) = (/update_nonlinear,implicit/)  ! OSS_vx
  setterms(6,:) = (/update_transient,explicit/)  ! Force
  call rkinteg%create(setterms,settable)

  ! Create problems
  call myprob%create(gdata%ndime)
  call mydisc%create(myprob)
  call mydisc%vars_block(myprob,vars_block)
  call mydisc%dof_coupling(myprob,dof_coupling)
  call cg_iss_oss_rk_momentum%create(myprob,mydisc)
  call cg_iss_oss_rk_momentum_rhs%create(myprob,mydisc)
  call cg_iss_oss_rk_pressure%create(myprob,mydisc)
  call cg_iss_oss_rk_momentum_update%create(myprob,mydisc)
  call cg_iss_oss_rk_projection_update%create(myprob,mydisc)
  call lapla_p_integration%create(myprob,mydisc)
  call mass_u_integration%create(myprob,mydisc)
  call mass_x_integration%create(myprob,mydisc)
  cg_iss_oss_rk_momentum%rkinteg => rkinteg
  cg_iss_oss_rk_momentum_rhs%rkinteg => rkinteg
  cg_iss_oss_rk_pressure%tinteg  => tinteg
  cg_iss_oss_rk_momentum_update%rkinteg => rkinteg
  rkinteg%dtinv   = 1.0_rp
  rkinteg%ftime   = 1.0_rp
  mydisc%kfl_proj = 0
  mydisc%kfl_lump = 0
  myprob%kfl_conv = 0
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
       &               which_approx,time_steps_to_store=3+rkinteg%rk_table_implicit%stage,       &
       &               hierarchical_basis=.false.,static_condensation=.false.,num_continuity=1)

  ! Initialize VTK output
  call fevtk%initialize(f_trian,fe_space,myprob,senv,dir_path_out,prefix,linear_order=.true.)

  ! Create dof info
  call create_dof_info(dof_descriptor,f_trian,fe_space,blk_graph,gtype)

  ! Assign analytical solution
  if(gdata%ndime==2) then
     call fe_space%set_analytical_code((/1,2,3,0,0/),(/1,1,0,0,0/))
  else
     write(*,*) 'analytical function not ready for 3D'
  end if

  ! Create picard nonlinear operators
  call momentum_operator%build(blk_graph,senv)
  call pressure_operator%build(blk_graph)
  call momentum_update_operator%build(blk_graph)
  call projection_update_operator%build(blk_graph)

  ! Do time steps
  call do_time_steps_rk_nsi(rkinteg,sctrl,1.0e-7_rp,100,senv,fe_space,momentum_operator,     &
       &                    cg_iss_oss_rk_momentum,pressure_operator,cg_iss_oss_rk_pressure, &
       &                    momentum_update_operator,cg_iss_oss_rk_momentum_update,          &
       &                    projection_update_operator,cg_iss_oss_rk_projection_update,      &
       &                    cg_iss_oss_rk_momentum_rhs,tinteg)

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

  ! Deallocate
  call memfree(continuity,__FILE__,__LINE__)
  call memfree(order,__FILE__,__LINE__)
  call memfree(material,__FILE__,__LINE__)
  call memfree(problem,__FILE__,__LINE__)
  call memfree(which_approx,__FILE__,__LINE__)
  call memfree(vars_block,__FILE__,__LINE__)
  call memfree(dof_coupling,__FILE__,__LINE__)
  call momentum_operator%free()
  call pressure_operator%free()
  call momentum_update_operator%free()
  call projection_update_operator%free()
  call fevtk%free
  call blk_graph%free()
  call fe_space_free(fe_space) 
  call myprob%free
  call mydisc%free
  call error_compute%free
  call cg_iss_oss_rk_momentum%free
  call cg_iss_oss_rk_momentum_rhs%free
  call cg_iss_oss_rk_pressure%free
  call cg_iss_oss_rk_momentum_update%free
  call cg_iss_oss_rk_projection_update%free
  call lapla_p_integration%free
  call mass_u_integration%free
  call mass_x_integration%free
  call rkinteg%free
  call dof_descriptor_free(dof_descriptor)
  call triangulation_free(f_trian)
  call conditions_free(f_cond)
  call reference_element_free(geo_reference_element)
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
    if (.not. (numargs==8) ) then
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
  subroutine do_time_steps_rk_nsi(rkinteg,sctrl,sttol,maxst,env,fe_space,momentum_operator,    &
       &                          momentum_integration,pressure_operator,pressure_integration, &
       &                          momentum_update_operator,momentum_update_integration,        &
       &                          projection_update_operator,projection_update_integration,    &
       &                          momentum_rhs_integration,tinteg)
    implicit none
    type(rungekutta_integrator_t)        , intent(inout) :: rkinteg
    type(solver_control_t)               , intent(inout) :: sctrl
    real(rp)                             , intent(in)    :: sttol
    integer(ip)                          , intent(in)    :: maxst
    class(abstract_environment_t)        , intent(in)    :: env
    type(fe_space_t)                     , intent(inout) :: fe_space
    class(momentum_operator_t)           , intent(inout) :: momentum_operator
    class(pressure_operator_t)           , intent(inout) :: pressure_operator
    class(momentum_update_operator_t)    , intent(inout) :: momentum_update_operator
    class(projection_update_operator_t)  , intent(inout) :: projection_update_operator
    class(discrete_integration_t), target, intent(inout) :: momentum_integration
    class(discrete_integration_t), target, intent(inout) :: pressure_integration
    class(discrete_integration_t), target, intent(inout) :: momentum_update_integration
    class(discrete_integration_t), target, intent(inout) :: projection_update_integration
    class(discrete_integration_t), target, intent(inout) :: momentum_rhs_integration
    class(time_integration_t)            , intent(inout) :: tinteg
    ! Locals
    type(discrete_integration_pointer_t) :: approx(1)
    integer(ip) :: istage,nstage,istep
    real(rp)    :: rtime,ctime,prevtime

    ! Initialize time steps
    rtime = 1.0_rp
    istep = 1

    ! Update steady boundary conditions
    call update_strong_dirichlet_bcond(fe_space,f_cond)
          
    ! Update analytical/time dependent boundary conditions
    call update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),rkinteg%ctime,fe_space)

    ! Compute constant operators
    ! Momentum
    approx(1)%p => momentum_integration
    momentum_integration%integration_stage = update_constant
    call momentum_operator%fill_constant(approx,fe_space) 
    ! Momentum update
    approx(1)%p => momentum_update_integration
    momentum_update_integration%integration_stage = update_constant
    call momentum_update_operator%fill_constant(approx,fe_space)
    !**************************** NSI specific tasks *****************************************!
    ! Pressure
    approx(1)%p => pressure_integration
    pressure_integration%integration_stage = update_constant
    call pressure_operator%fill_constant(approx,fe_space)
    ! Projection update
    approx(1)%p => projection_update_integration
    projection_update_integration%integration_stage = update_constant
    call projection_update_operator%fill_constant(approx,fe_space)
    ! Preconditioner matrices
    approx(1)%p => lapla_p_integration
    call volume_integral(approx,fe_space,pressure_operator%lapla_p_matrix)
    approx(1)%p => mass_u_integration
    call volume_integral(approx,fe_space,pressure_operator%mass_u_matrix)
    call volume_integral(approx,fe_space,momentum_update_operator%mass_u_matrix)
    approx(1)%p => mass_x_integration
    call volume_integral(approx,fe_space,projection_update_operator%mass_x_matrix)
    ! scal --> L*dt, M/dt
    !*****************************************************************************************!

    ! Time steps loop
    step: do while (rkinteg%ctime<rkinteg%ftime.and.rtime>sttol.and.istep<=maxst)

       ! Current time
       prevtime = rkinteg%ctime
       rkinteg%ctime = rkinteg%ctime + 1.0_rp/rkinteg%dtinv
       write(*,'(a)') '============================================================'
       write(*,'(a15,i8,a,a20,e15.8)') 'Time step:     ',istep,', ','Current time: ',rkinteg%ctime

       ! Loop over stages
       stage: do istage=1,1!rkinteg%rk_table(1)%p%stage

          ! Set current time
          ctime = prevtime + 1.0_rp/rkinteg%dtinv*rkinteg%rk_table(1)%p%c(istage)
          rkinteg%istage = istage
          write(*,'(a)') '------------------------------------------------------------'
          write(*,'(a20,i3,a,a21,e15.8)') 'Runge-Kutta stage:  ',istage,',','Current time: ',ctime

          ! Set momentum operator scalars
          momentum_operator%dt  = rkinteg%dtinv
          momentum_operator%aii = rkinteg%rk_table_implicit%A(istage,istage)
         
          ! Skip 1st stage
          if(istage==1.and.momentum_operator%aii==0) then
             write(*,*) 'First stage skipped for momentum equation.'
             ! (U_n --> U_1)
             call update_nonlinear_solution(fe_space,approx(1)%p%working_vars,3,3+1)
          else
          
             ! Update analytical/time dependent boundary conditions
             call update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,fe_space)

             ! Compute momentum transient operators
             approx(1)%p => momentum_integration
             momentum_integration%integration_stage = update_transient
             call momentum_operator%fill_transient(approx,fe_space) 

             ! Compute momentum operand (RHS)
             approx(1)%p => momentum_rhs_integration
             call volume_integral(approx,fe_space,momentum_operator%block_vector)

             ! Momentum equation solution
             approx(1)%p => momentum_integration
             call momentum_operator%apply(sctrl,senv,approx,fe_space)

             ! Store unkno to istage position (U --> U_i)
             call update_nonlinear_solution(fe_space,approx(1)%p%working_vars,1,3+istage)             

          end if

          !**************************** NSI specific tasks *****************************************!
          ! Update boundary conditions (velocity derivative)
          call update_analytical_bcond((/(i,i=1,gdata%ndime)/),ctime,fe_space,2)

          ! Compute pressure transient operators
          approx(1)%p => pressure_integration
          pressure_integration%integration_stage = update_transient
          call pressure_operator%fill_transient(approx,fe_space) 

          ! Pressure equation solution
          approx(1)%p => pressure_integration
          tinteg%ctime = ctime
          call pressure_operator%apply(sctrl,senv,approx,fe_space)

          ! Store unkno to istage position (P --> P_i)
          call update_nonlinear_solution(fe_space,approx(1)%p%working_vars,1,3+istage)        
          !*****************************************************************************************!        

       end do stage

!!$       ! Set current time
!!$       ctime = rkinteg%ctime
!!$       write(*,'(a)') '------------------------------------------------------------'
!!$       write(*,'(a24,a21,e15.8)') 'Runge-Kutta update      ','Current time: ',ctime
!!$
!!$       ! Update analytical/time dependent boundary conditions
!!$       call update_analytical_bcond((/(i,i=1,gdata%ndime+1)/),ctime,fe_space)
!!$
!!$       ! Compute momentum update transient operators
!!$       approx(1)%p => momentum_update_integration
!!$       momentum_update_integration%integration_stage = update_transient
!!$       call momentum_update_operator%fill_transient(approx,fe_space) 
!!$       
!!$       ! Momentum update equation solution
!!$       approx(1)%p => momentum_update_integration
!!$       call momentum_update_operator%apply(sctrl,senv,approx,fe_space)
!!$
!!$       !**************************** NSI specific tasks *****************************************!
!!$       ! Compute projection transient operators
!!$       approx(1)%p => projection_update_integration
!!$       projection_update_integration%integration_stage = update_transient
!!$       call projection_update_operator%fill_transient(approx,fe_space) 
!!$
!!$       ! Projection update equation solution
!!$       approx(1)%p => projection_update_integration
!!$       call projection_update_operator%apply(sctrl,senv,approx,fe_space)
!!$
!!$       ! Compute pressure transient operators
!!$       approx(1)%p => pressure_integration
!!$       pressure_integration%integration_stage = update_transient
!!$       call pressure_operator%fill_transient(approx,fe_space) 
!!$       
!!$       ! Pressure equation solution
!!$       approx(1)%p => pressure_integration
!!$       tinteg%ctime = ctime
!!$       call pressure_operator%apply(sctrl,senv,approx,fe_space)
!!$       !*****************************************************************************************!  

       ! Check steady state

    end do step

  end subroutine do_time_steps_rk_nsi
  
end program test_blk_nsi_cg_iss_oss_rk
