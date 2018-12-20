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

module jacobi_preconditioner_names
 ! Tools
 use types_names
 use list_types_names
 use allocatable_array_names
 use FPL
 
 ! Integration related modules
 use triangulation_names
 use fe_space_names
 use fe_operator_names
 
 ! Linear Algebra related modules
 use operator_names
 use par_scalar_array_names
 use serial_scalar_array_names
 use vector_space_names
 use vector_names
 
 use matrix_names
 use base_sparse_matrix_names
 use sparse_matrix_parameters_names
 use sparse_matrix_names
 use par_sparse_matrix_names
 use direct_solver_names
 use direct_solver_parameters_names

#ifdef ENABLE_BLAS
 use blas77_interfaces_names
#endif

#ifdef ENABLE_LAPACK
 use lapack77_interfaces_names
#endif
  
 ! Parallel communication-related data structures
 use environment_names
 
 implicit none
# include "debug.i90"
 private
 
 integer(ip), parameter :: jacobi_preconditioner_STATE_START    = 0
 integer(ip), parameter :: jacobi_preconditioner_STATE_CREATED  = 1
 integer(ip), parameter :: jacobi_preconditioner_STATE_SYMBOLIC = 2 ! Symbolic data already computed
 integer(ip), parameter :: jacobi_preconditioner_STATE_NUMERIC  = 3 ! Numerical data already computed
 
  !-----------------------------------------------------------------
  ! State transition diagram for type(jacobi_preconditioner_t)
  !-----------------------------------------------------------------
  ! Input State         | Action                | Output State 
  !-----------------------------------------------------------------
  ! Start               | create                | Created
  ! Start               | free_clean            | Start
  ! Start               | free_symbolic         | Start
  ! Start               | free_numeric          | Start
  ! Start               | update_matrix         | Start ! it does nothing

 
  ! Created             | symbolic_setup        | Symbolic         ! perform symbolic_setup()
  ! Created             | numerical_setup       | Numeric          ! perform symbolic_setup()+numerical_setup()
  ! Created             | apply                 | Numeric          ! perform symbolic_setup()+numerical_setup()
  ! Created             | free_clean            | Start
  ! Created             | free_symbolic         | Create           ! it does nothing
  ! Created             | free_numeric          | Create           ! it does nothing
  ! Created             | update_matrix         | Create           ! it does nothing

  ! Symbolic            | symbolic_setup                        | Symbolic         ! it does nothing
  ! Symbolic            | numerical_setup                       | Numeric          ! perform numerical_setup() 
  ! Symbolic            | apply                                 | Numeric          ! perform numerical_setup()
  ! Symbolic            | free_clean                            | Start
  ! Symbolic            | free_symbolic                         | Created
  ! Symbolic            | free_numeric                          | Symbolic         ! it does nothing
  ! Symbolic            | update_matrix + same_nonzero_pattern  | Symbolic         ! it does nothing
  ! Symbolic            | update_matrix + !same_nonzero_pattern | Symbolic         ! free_symbolic()+symbolic_setup()
    
    
  ! Numeric             | symbolic_setup                        | Numeric          ! it does nothing
  ! Numeric             | numeric_setup                         | Numeric          ! it does nothing
  ! Numeric             | apply                                 | Numeric          ! it does nothing
  ! Numeric             | free_numeric                          | Symbolic
  ! Numeric             | free_symbolic                         | Created
  ! Numeric             | free_clean                            | Start
  ! Numeric             | update_matrix + same_nonzero_pattern  | Numeric          ! free_numerical_setup()+numerical_setup()
  ! Numeric             | update_matrix + !same_nonzero_pattern | Numeric          ! free_numerical_setup()+free_symbolic_setup()
                                                                                   ! symbolic_setup()+numeric_setup()
 type, extends(operator_t) :: jacobi_preconditioner_t
   private
   integer(ip)                                 :: state  = jacobi_preconditioner_STATE_START   
   class(environment_t), pointer               :: environment => NULL()

   ! Pointer to the fe_nonlinear_operator_t this jacobi_preconditioner_t instance has been created from
   type(fe_operator_t) , pointer               :: fe_nonlinear_operator => NULL()
   
   class(vector_t), allocatable                :: inverse_diagonal        
  
   ! This pointer is set-up during jacobi_preconditioner_t%create() and re-used in the rest of stages.
   ! Therefore, type(parameter_list_t) to which type(jacobi_preconditioner_t) points to MUST NOT BE
   ! freed before type(jacobi_preconditioner_t).
   type(parameterlist_t)         , pointer     :: jacobi_preconditioner_params   => NULL()
   
 contains
   procedure, non_overridable  :: jacobi_preconditioner_create_w_parameter_list
   procedure, non_overridable  :: jacobi_preconditioner_create_wo_parameter_list
   generic                     :: create => jacobi_preconditioner_create_w_parameter_list, &
                                            jacobi_preconditioner_create_wo_parameter_list 
   
   procedure,                  private :: create_vector_spaces                  => jacobi_preconditioner_create_vector_spaces
   procedure,                  private :: create_and_allocate_inverse_diagonal  => jacobi_preconditioner_create_and_allocate_inverse_diagonal
   procedure,                  private :: free_and_destroy_inverse_diagonal     => jacobi_preconditioner_free_and_destroy_inverse_diagonal
 
   ! State transition handling-related TBPs
   procedure, non_overridable, private :: set_state_start              => jacobi_preconditioner_set_state_start
   procedure, non_overridable, private :: set_state_created            => jacobi_preconditioner_set_state_created
   procedure, non_overridable, private :: set_state_symbolic           => jacobi_preconditioner_set_state_symbolic
   procedure, non_overridable, private :: set_state_numeric            => jacobi_preconditioner_set_state_numeric
   procedure, non_overridable, private :: state_is_start               => jacobi_preconditioner_state_is_start
   procedure, non_overridable, private :: state_is_created             => jacobi_preconditioner_state_is_created
   procedure, non_overridable, private :: state_is_symbolic            => jacobi_preconditioner_state_is_symbolic
   procedure, non_overridable, private :: state_is_numeric             => jacobi_preconditioner_state_is_numeric
 
   ! Symbolic setup-related TBPs
   procedure, non_overridable          :: symbolic_setup               => jacobi_preconditioner_symbolic_setup

   ! Numerical setup-related TBPs
   procedure, non_overridable          :: numerical_setup              => jacobi_preconditioner_numerical_setup

   ! Apply related TBPs
   procedure                           :: apply                        => jacobi_preconditioner_apply
   procedure                           :: apply_add                    => jacobi_preconditioner_apply_add
   
   ! Free-related TBPs
   procedure, non_overridable          :: free                         => jacobi_preconditioner_free
   procedure, non_overridable          :: free_clean                   => jacobi_preconditioner_free_clean
   procedure, non_overridable          :: free_symbolic_setup          => jacobi_preconditioner_free_symbolic_setup   
   procedure, non_overridable          :: free_numerical_setup         => jacobi_preconditioner_free_numerical_setup
 
   procedure, non_overridable, private :: am_i_l1_task                 => jacobi_preconditioner_am_i_l1_task
   procedure                           :: is_linear                    => jacobi_preconditioner_is_linear
   procedure, private                  :: get_par_environment          => jacobi_preconditioner_get_par_environment
   procedure, private                  :: set_par_environment          => jacobi_preconditioner_set_par_environment
   
   ! Miscellaneous 
   procedure, private                  :: get_par_sparse_matrix        => jacobi_preconditioner_get_par_sparse_matrix
   procedure, private                  :: get_fe_space                 => jacobi_preconditioner_get_fe_space
   procedure, private                  :: get_par_fe_space             => jacobi_preconditioner_get_par_fe_space
   procedure                 , private :: is_operator_associated       => jacobi_preconditioner_is_operator_associated
   procedure                 , private :: nullify_operator             => jacobi_preconditioner_nullify_operator 
     
   ! Update-matrix related TBPs
   procedure                           :: update_matrix                                   => jacobi_preconditioner_update_matrix   
end type jacobi_preconditioner_t

interface jacobi_preconditioner_t
  module procedure create_jacobi_preconditioner
end interface jacobi_preconditioner_t
 
contains

#include "sbm_jacobi_preconditioner.i90"

end module jacobi_preconditioner_names
