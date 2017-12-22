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
module time_stepping_names
  ! sbadia: to check whether all this needed
  use types_names
  use memor_names

  use vector_space_names
  use fe_space_names
  use operator_names
  use vector_names

  use assembler_names
  use sparse_assembler_names
  use block_sparse_assembler_names  
  use par_sparse_assembler_names

  use sparse_matrix_names, only: sparse_matrix_t
  use block_sparse_matrix_names
  use par_sparse_matrix_names

  use serial_scalar_array_names
  use serial_block_array_names
  use par_scalar_array_names

  use array_names
  use matrix_names
  use discrete_integration_names
  use environment_names
  use direct_solver_names
  use block_layout_names

  implicit none
# include "debug.i90"

  private

  character(:), parameter :: forward_euler  = "forward_euler"
  character(:), parameter :: backward_euler = "backward_euler"
  character(:), parameter :: crank_nicolson = "crank_nicolson"

  ! states to be defined
  integer(ip), parameter :: created             = 0
  integer(ip), parameter :: residual_computed   = 1 
  integer(ip), parameter :: tangent_computed    = 2 
  integer(ip), parameter :: assembler_computed  = 3

  ! sbadia: to make auto-documentation style
  !
  !This operator is the global RK operator as follows.
  !In RK methods, for a s-stage method, given u_0, we must compute
  ![ v_1, v_s ] st
  !v_1 + A ( t + c_1 dt , v_0 + a_11 v_1 ) = 0
  !...
  !v_i + A ( t + c_i dt, v_0 + \sum_{j=1}^s a_ij v_j ) = 0
  !...
  !to finally get
  !u_1 =  u_0 + \sum_{i=1}^s b_i v_i = 0.
  ! As a result, the RK problem can be stated in compact form as R(U)=0, where
  ! U = [v_1, ..., v_s]. The objective is to solve this problem and 
  ! extract u_1 = u_0 + \sum_{i=1}^s b_i v_i = 0.
  ! One must compute, given U, its application R(U). So, I want a apply method that, 
  ! given the stage i I want to get, it applies the i-stage R_i(U) application.
  ! However, in DIRK RK methods, the ones usually used, the coupling is weaker, and 
  ! should exploit it. R_i(U) only depends on v_j, for j \leq i.
  ! In this case, we want to solve at every stage
  ! R_i( v_0, ..., v_i ) = 0, where v_0, ..., v_i-1 are known, which we can denote as
  ! R*_i(v_i) = 0.
  ! In order to solve this potential nonlinear operator, we must provide its tangent too.
  ! We need the nonlinear operator: Given i, and v_1, ..., v_i-1, v_i+1,...,v_s, 
  ! R_ij(V) = R_i (v_1, ...,v_i-1,V,v_i+1,...,v_s).
  ! Next, I compute the apply and tangent of this operator, which is all I need to apply
  ! the whole operator.

  ! sbadia: I must provide a method to get the solution
  ! u_1 =  u_0 + \sum_{i=1}^s b_i v_i = 0.

  type:: time_stepping_operator_t ! , extends(operator_t) commented for the moment,
                                  ! no implicit RK implemented yet
     private
     ! sbadia: can we put here just the stage block
     class(time_stepping_stage_fe_operator_t)           :: fe_op
     type(time_stepping_scheme_t)                       :: scheme
     class(vector_t)                          , pointer :: initial_value => NULL()
     real(rp)                                           :: dt
     class(vector_t)                      , allocatable :: dofs_stages(:)

     ! sbadia: For the moment, we are not interested in full RK implementations,
     ! even though it would be an easy paper about preconditioning these schemes
     ! So, we don't really need to use the block assembler for the
     ! all-stages operator. We note that the matrix is not needed to be computed
     ! for every stage since it is always the same
     !class(assembler_t)                   , allocatable :: assembler
     !type(block_vector_t)                               :: dofs_stages_block_vector
   contains
     !procedure :: create_from_operators => time_stepping_operator_create_from_operators
     !procedure :: set_initial_data      => time_stepping_operator_set_initial_data
     !procedure :: set_time_step_size    => time_stepping_operator_set_time_step_size     
     !procedure, private :: allocate_dofs_stages      => time_stepping_operator_allocate_dofs_stages
     !!!procedure :: create             => time_stepping_operator_create
     ! sbadia: It must be defined since it is an operator, but for the moment
     !!! we do not want to use it. Dummy implementation...
     !!!procedure :: apply              => time_stepping_operator_apply
     !!! sbadia: to be implemented
     !!! procedure :: free               => time_stepping_operator_free
     !!!procedure, private :: apply_row                 => time_stepping_operator_apply_row
     !!!procedure, private :: compute_tangent_block     => time_stepping_operator_compute_tangent_block
     !!!procedure, private :: compute_mass_matrix       => time_stepping_operator_compute_mass_matrix
     !!!procedure, private :: set_evaluation_point_row  => time_stepping_operator_set_evaluation_point_row	 
  end type time_stepping_operator_t
  
  ! This operator represents the operator R_ij = R_i(x,x,v_j,x,x),
  ! where the x denotes that the unknown is fixed. It requires to add
  ! contributions from the nonlinear operator A and mass matrix M as
  ! alpha*A + beta*M
  type, extends(operator_t) :: time_stepping_operator_stage_t
     private
     type(time_stepping_operator_t), pointer :: ts_op 
     integer(ip) :: i
     integer(ip) :: j
   contains
     
  end type time_stepping_operator_stage_t
  
  type :: dirk_time_stepping_solver_t
     private
     type(time_stepping_operator_t) :: op
     type(time_stepping_operator_block_t) :: op_block
     type(nonlinear_solver_t) :: nl_solver
   contains
  end type dirk_time_stepping_solver_t

  type :: time_stepping_scheme_t
     private
     integer(ip) :: time_integrator
     integer(ip) :: order
     integer(ip) :: num_stages
     real(rp)    , allocatable :: a(:,:)
     real(rp)    , allocatable :: b(:)
     real(rp)    , allocatable :: c(:)
   contains
     !procedure :: create => time_stepping_scheme_create
  end type time_stepping_scheme_t
  
  
  public :: time_stepping_operator_t, dirk_time_stepping_solver_t

end module time_stepping_names
