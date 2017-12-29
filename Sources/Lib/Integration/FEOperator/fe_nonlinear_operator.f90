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
module fe_nonlinear_operator_names
  use types_names
  use memor_names

  use vector_space_names
  use reference_fe_names
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

  integer(ip), parameter :: created             = 0
  integer(ip), parameter :: residual_computed   = 1 
  integer(ip), parameter :: tangent_computed    = 2 
  integer(ip), parameter :: assembler_computed  = 3

  integer(ip), parameter :: free_setup          = 0
  

  ! State transition diagram for type(fe_nonlinear_operator_t)
  ! -------------------------------------------------
  ! Input State | Action               | Output State 
  ! -------------------------------------------------
  ! start       | create               | created
  ! start       | free_setup          | start
  ! start       | free_setup  | start
  ! start       | free_clean           | start 
  ! start       | free                 | start 

  ! created     | setup       | setup
  ! created     | setup      | setup
  ! created     | get_tangent+         | setup
  !               get_translation+               "
  !               get_domain_vector_space+       "
  !               get_range_vector_space+        "
  !               abort_if_not_in_range+         "
  !               abort_if_not_in_domain         "
  ! created     | free_setup          | created
  ! created     | free_setup          | created
  ! created     | free_clean           | start
  ! created     | free                 | start

  ! setup    | setup       | setup
  ! setup    | setup      | setup
  ! setup    | free_setup          | setup
  ! setup    | free_setup          | created
  ! setup    | free_clean           | start
  ! setup    | free                 | start
  ! setup    | get_tangent+         | setup
  !                         get_translation+               "
  !                         get_domain_vector_space+       "
  !                         get_range_vector_space+        "
  !                         abort_if_not_in_range+         "
  !                         abort_if_not_in_domain         "

  ! setup    | setup       | setup
  ! setup    | setup      | setup
  ! setup    | free_setup          | setup
  ! setup    | free_setup  | created
  ! setup    | free                 | start
  ! setup    | free_clean           |  
  ! setup    | get_tangent+         | setup
  !                        get_translation+               "
  !                        get_domain_vector_space+       "
  !                        get_range_vector_space+        "
  !                        abort_if_not_in_range+         "
  !                        abort_if_not_in_domain         "


  type, extends(operator_t):: fe_nonlinear_operator_t
     private
     integer(ip)                       , pointer     :: state => NULL()
     character(:)                      , allocatable :: sparse_matrix_storage_format
     class(serial_fe_space_t)          , pointer     :: test_fe_space          => NULL() ! test_fe_space
     class(serial_fe_space_t)          , pointer     :: trial_fe_space         => NULL() ! To be used in the future
     class(discrete_integration_t)     , pointer     :: discrete_integration   => NULL()
     class(assembler_t)                , pointer     :: assembler => NULL()
     logical                           , allocatable :: diagonal_blocks_symmetric_storage(:)
     logical                           , allocatable :: diagonal_blocks_symmetric(:)
     integer(ip)                       , allocatable :: diagonal_blocks_sign(:)
   contains
     procedure :: create                                => fe_nonlinear_operator_create
     procedure :: free                                  => fe_nonlinear_operator_free
     procedure :: set_evaluation_point                  => fe_nonlinear_operator_set_evaluation_point  
     procedure :: reallocate_after_remesh               => fe_nonlinear_operator_reallocate_after_remesh

     procedure :: compute_tangent                       => fe_nonlinear_operator_compute_tangent
     procedure :: compute_residual                      => fe_nonlinear_operator_compute_residual  

     procedure :: apply                                 => fe_nonlinear_operator_apply
     procedure :: apply_add                             => fe_nonlinear_operator_apply_add
     procedure :: is_linear                             => fe_nonlinear_operator_is_linear

     procedure :: get_tangent                           => fe_nonlinear_operator_get_tangent 
     procedure :: get_matrix                            => fe_nonlinear_operator_get_matrix
     procedure :: get_translation                       => fe_nonlinear_operator_get_translation
     procedure :: get_fe_space                          => fe_nonlinear_operator_get_fe_space
     procedure :: get_trial_fe_space                    => fe_nonlinear_operator_get_trial_fe_space
     procedure :: get_discrete_integration              => fe_nonlinear_operator_get_discrete_integration
     procedure :: get_domain_vector_space               => fe_nonlinear_operator_get_domain_vector_space
     procedure :: get_range_vector_space                => fe_nonlinear_operator_get_range_vector_space
	    procedure :: get_state                             => fe_nonlinear_operator_get_state
	    procedure :: get_diagonal_blocks_symmetric_storage => fe_nonlinear_operator_get_diagonal_blocks_symmetric_storage
	    procedure :: get_diagonal_blocks_symmetric         => fe_nonlinear_operator_get_diagonal_blocks_symmetric
	    procedure :: get_diagonal_blocks_sign              => fe_nonlinear_operator_get_diagonal_blocks_sign
     procedure :: get_sparse_matrix_storage_format      => fe_nonlinear_operator_get_sparse_matrix_storage_format 
     procedure :: get_assembler                         => fe_nonlinear_operator_get_assembler        
	 
     procedure :: allocate_state                        => fe_nonlinear_operator_allocate_state
     procedure :: set_state                             => fe_nonlinear_operator_set_state
     procedure :: deallocate_state                      => fe_nonlinear_operator_deallocate_state
	    procedure :: abort_if_not_in_range                 => fe_nonlinear_operator_abort_if_not_in_range
     procedure :: abort_if_not_in_domain                => fe_nonlinear_operator_abort_if_not_in_domain

     procedure, private :: create_serial_assembler     => fe_nonlinear_operator_create_serial_assembler
     procedure, private :: create_par_assembler        => fe_nonlinear_operator_create_par_assembler
     procedure, private :: create_vector_spaces        => fe_nonlinear_operator_create_vector_spaces
  end type fe_nonlinear_operator_t
  
type, extends(fe_nonlinear_operator_t) :: fe_affine_operator_t
  private
contains
  procedure          :: compute                     => fe_affine_operator_compute
  procedure          :: apply                       => fe_affine_operator_apply
  procedure          :: apply_add                   => fe_affine_operator_apply_add
  procedure          :: is_linear                   => fe_affine_operator_is_linear
  procedure          :: get_tangent                 => fe_affine_operator_get_tangent
  procedure          :: get_translation             => fe_affine_operator_get_translation
  procedure          :: get_matrix                  => fe_affine_operator_get_matrix
end type fe_affine_operator_t

type, extends(fe_nonlinear_operator_t) :: scale_add_fe_operator_t
  ! This operator represents alpha*A+B, given the fe_nonlinear_operator_t's A and B
  ! and scalars alpha and beta
  private
  class(fe_nonlinear_operator_t), pointer :: op1 => NULL() ! A
  class(fe_nonlinear_operator_t), pointer :: op2 => NULL() ! B
  real(rp)                                :: alpha        ! alpha
contains
  procedure          :: create                      => scale_add_fe_operator_create
  procedure          :: create_from_operators       => scale_add_fe_operator_create_from_operators
  procedure          :: free                        => scale_add_fe_operator_free
  procedure          :: set_evaluation_point        => scale_add_fe_operator_set_evaluation_point
  procedure          :: compute_tangent             => scale_add_fe_operator_compute_tangent
  procedure          :: compute_residual            => scale_add_fe_operator_compute_residual  
  procedure          :: is_linear                   => scale_add_fe_operator_is_linear
  procedure          :: set_scalars                 => scale_add_fe_operator_set_scalars
end type scale_add_fe_operator_t


  ! Types
  public :: fe_nonlinear_operator_t, fe_affine_operator_t, scale_add_fe_operator_t

contains

  
#include "sbm_fe_nonlinear_operator.i90"
#include "sbm_fe_affine_operator.i90"
#include "sbm_scal_add_fe_operator.i90"  
end module
