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
  use sparse_matrix_parameters_names
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

  use fe_operator_names
  use nonlinear_solver_names
  
  use time_stepping_mass_discrete_integration_names
  implicit none
# include "debug.i90"

  private

  character(*), parameter :: forward_euler     = "forward_euler"
  character(*), parameter :: backward_euler    = "backward_euler"
  character(*), parameter :: trapezoidal_rule  = "trapezoidal_rule"
  character(*), parameter :: mid_point_implicit = "mid_point_implicit"
  character(*), parameter :: mid_point_explicit = "mid_point_explicit"
  character(*), parameter :: imex_rk_2_2_1_implicit = "imex_rk_2_2_1_implicit"
  character(*), parameter :: imex_rk_2_2_1_explicit = "imex_rk_2_2_1_explicit"
  character(*), parameter :: imex_rk_2_2_2_implicit = "imex_rk_2_2_2_implicit"
  character(*), parameter :: imex_rk_2_2_2_explicit = "imex_rk_2_2_2_explicit"
  character(*), parameter :: imex_rk_2_3_implicit = "imex_rk_2_3_implicit"
  character(*), parameter :: imex_rk_2_3_explicit = "imex_rk_2_3_explicit"
  character(*), parameter :: runge_kutta_3     = "runge_kutta_3"
  character(*), parameter :: runge_kutta_4     = "runge_kutta_4"
  character(*), parameter :: runge_kutta_4_3_8 = "runge_kutta_4_3_8"
  character(*), parameter :: imex_rk_4_3_implicit = "imex_rk_4_3_implicit"
  character(*), parameter :: imex_rk_4_3_explicit = "imex_rk_4_3_explicit" 

  ! states to be defined
  integer(ip), parameter :: created             = 0
  integer(ip), parameter :: residual_computed   = 1 
  integer(ip), parameter :: tangent_computed    = 2 
  integer(ip), parameter :: assembler_computed  = 3
  
  type :: butcher_tableau_t
    private
    character(:), allocatable :: time_integration_scheme
    integer(ip)               :: order      = 0
    integer(ip)               :: num_stages = 0
    real(rp)    , allocatable :: a(:,:)
    real(rp)    , allocatable :: b(:)
    real(rp)    , allocatable :: c(:)
  contains
    procedure          :: create          => butcher_tableau_create
    procedure, private :: allocate_arrays => butcher_tableau_allocate_arrays  
    procedure          :: free            => butcher_tableau_free
  end type butcher_tableau_t
  
  ! This operator represents the operator R_ij = R_i(x,x,v_j,x,x),
  ! where the x denotes that the unknown is fixed. It requires to add
  ! contributions from the nonlinear operator A and mass nonlinear operator 
  ! M as alpha*A + beta*M
  type, extends(fe_operator_t) :: time_stepping_stage_fe_operator_t
     private
     type(time_stepping_operator_t), pointer :: ts_op      => NULL()
     class(fe_operator_t), pointer :: fe_op   => NULL()
     type(fe_operator_t)           :: mass_op 
     type(mass_discrete_integration_t)       :: mass_integration
     integer(ip) :: i
     integer(ip) :: j
     ! Scratch data to support the efficient implementation of 
     ! time_stepping_stage_fe_operator_set_evaluation_point
     class(vector_t), allocatable    :: aux
   contains
     procedure :: create                => time_stepping_stage_fe_operator_create
     procedure :: create_from_operators => time_stepping_stage_fe_operator_create_from_operators
     procedure :: create_mass_operator  => time_stepping_stage_fe_operator_create_mass_operator
     procedure :: free                  => time_stepping_stage_fe_operator_free
     procedure :: set_row               => time_stepping_stage_fe_operator_set_row     
     procedure :: set_evaluation_point  => time_stepping_stage_fe_operator_set_evaluation_point
     procedure :: set_evaluation_time   => time_stepping_stage_fe_operator_set_evaluation_time
     procedure :: is_linear             => time_stepping_stage_fe_operator_is_linear
     procedure :: compute_internal_residual      => time_stepping_stage_fe_operator_compute_internal_residual
     procedure :: compute_tangent       => time_stepping_stage_fe_operator_compute_tangent
  end type time_stepping_stage_fe_operator_t

  ! sbadia: to make auto-documentation style
  !
  ! This operator is the global RK operator as follows.
  ! In RK methods, for a s-stage method, given u_0, we must compute
  ! [ v_1, ...,  v_s ] st
  ! v_1 + A ( t + c_1 dt , v_0 + a_11 v_1 ) = 0
  ! ...
  ! v_i + A ( t + c_i dt, v_0 + \sum_{j=1}^s a_ij v_j ) = 0
  ! ...
  ! to finally get
  ! u_1 =  u_0 + \sum_{i=1}^s b_i v_i.
  ! As a result, the RK problem can be stated in compact form as R(U)=0, where
  ! U = [v_1, ..., v_s]. The objective is to solve this problem and 
  ! extract u_1 = u_0 + \sum_{i=1}^s b_i v_i.
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

  type:: time_stepping_operator_t ! extends(operator_t) commented for the moment,
                                  ! no implicit RK implemented yet
     private
     type(butcher_tableau_t)                            :: scheme
     type(time_stepping_stage_fe_operator_t)            :: fe_op
     class(vector_t)                          , pointer :: initial_value => NULL()
     class(vector_t)                      , allocatable :: dofs_stages(:)
     real(rp)                                           :: dt, current_initial_time, final_time !time => initial_time (t0)
     
     ! sbadia: For the moment, we are not interested in full RK implementations,
     ! even though it would be an easy paper about preconditioning these schemes
     ! So, we don't really need to use the block assembler for the
     ! all-stages operator. We note that the matrix is not needed to be computed
     ! for every stage since it is always the same
   contains
     procedure :: create                          => time_stepping_operator_create
     procedure :: free                            => time_stepping_operator_free
     procedure :: update                          => time_stepping_operator_update
     procedure :: set_initial_data                => time_stepping_operator_set_initial_data
     procedure :: set_initial_time                => time_stepping_operator_set_initial_time 
     procedure :: set_final_time                  => time_stepping_operator_set_final_time 
     procedure :: set_time_step_size              => time_stepping_operator_set_time_step_size 
     procedure :: update_current_time             => time_stepping_operator_update_current_time 
     procedure :: get_time_step                   => time_stepping_operator_get_time_step_size
     procedure :: get_current_time                => time_stepping_operator_get_current_time
     procedure :: get_final_time                  => time_stepping_operator_get_final_time 
     procedure :: get_matrix                      => time_stepping_operator_get_matrix
     procedure :: get_order                       => time_stepping_operator_get_order 
     procedure :: has_finished                    => time_stepping_operator_has_finished
     procedure, private :: allocate_dofs_stages   => time_stepping_operator_allocate_dofs_stages
     procedure, private :: deallocate_dofs_stages => time_stepping_operator_deallocate_dofs_stages
     procedure, private :: get_stage_operator     => time_stepping_operator_get_stage_operator
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
  
  type :: dirk_solver_t
    private
    type(time_stepping_operator_t), pointer :: ts_op     => NULL()
    type(nonlinear_solver_t)      , pointer :: nl_solver => NULL()
  contains
    procedure :: create              => dirk_solver_create
    procedure :: apply               => dirk_solver_apply
    procedure :: advance_fe_function => dirk_solver_advance_fe_function
    procedure :: free                => dirk_solver_free
  end type dirk_solver_t
  
  public :: time_stepping_operator_t, dirk_solver_t
  
contains
  
#include "sbm_time_stepping_operator.i90"
#include "sbm_butcher_tableau.i90"
#include "sbm_time_stepping_stage_fe_operator.i90"
#include "sbm_dirk_solver.i90"
  
end module time_stepping_names
