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
module fe_operator_names
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

  ! NOTE: In the state transition diagram below we use the state "start" to refer
  !       to the instance state in which ".not. associated(this%state)" is .true.
  !       Let us not however, that "start" is not actually an state to which we 
  !       assign a parameter constant above.
  
  ! State transition diagram for type(fe_operator_t)
  ! -------------------------------------------------
  ! Input State | Action               | Output State 
  ! -------------------------------------------------
  ! start               | create               | created
  ! start               | free                 | start
  ! start               | is_linear            | start
 
  ! created             | compute_tangent      | tangent_computed
  ! created             | compute_residual     | residual_computed
  ! created             | set_evaluation_point | created           ! Re-initializes to zero LA data structures
  ! created             | reallocate_after_re* | created           ! Re-creates LA data structures
  ! created             | free                 | start
  ! created             | is_linear            | created
  ! created             | get*/abort*          | created
  
  ! tangent_computed    | compute_residual     | assembler_computed
  ! tangent_computed    | compute_tangent      | tangent_computed   ! Does nothing
  ! tangent_computed    | set_evaluation_point | created            ! Re-initializes to zero LA data structures
  ! tangent_computed    | reallocate_after_re* | created            ! Re-creates LA data structures
  ! tangent_computed    | free                 | start 
  ! tangent_computed    | is_linear            | tangent_computed
  ! tangent_computed    | get*/abort*          | tangent_computed
 
  ! residual_computed   | compute_residual     | residual_computed   ! Does nothing
  ! residual_computed   | compute_tangent      | assembler_computed  
  ! residual_computed   | set_evaluation_point | created             ! Re-initializes to zero LA data structures
  ! residual_computed   | reallocate_after_re* | created             ! Re-creates LA data structures
  ! residual_computed   | free                 | start 
  ! residual_computed   | is_linear            | residual_computed
  ! residual_computed   | get*/abort*          | residual_computed
  
  ! assembler_computed  | compute_residual     | assembler_computed   ! Does nothing
  ! assembler_computed  | compute_tangent      | assembler_computed   ! Does nothing
  ! assembler_computed  | set_evaluation_point | created              ! Re-initializes to zero LA data structures
  ! assembler_computed  | reallocate_after_re* | created              ! Re-creates LA data structures
  ! assembler_computed  | free                 | start 
  ! assembler_computed  | is_linear            | assembler_computed
  ! assembler_computed  | get*/abort*          | assembler_computed

  type, extends(operator_t):: fe_operator_t
     private
     integer(ip)                       , pointer     :: state                    => NULL()
     character(:)                      , allocatable :: sparse_matrix_storage_format
     class(serial_fe_space_t)          , pointer     :: test_fe_space            => NULL() ! test_fe_space
     class(serial_fe_space_t)          , pointer     :: trial_fe_space           => NULL() ! To be used in the future
     class(discrete_integration_t)     , pointer     :: discrete_integration     => NULL()
     class(assembler_t)                , pointer     :: assembler                => NULL()
     class(vector_t)                   , pointer     :: current_evaluation_point => NULL()

     logical                           , allocatable :: diagonal_blocks_symmetric_storage(:)
     logical                           , allocatable :: diagonal_blocks_symmetric(:)
     integer(ip)                       , allocatable :: diagonal_blocks_sign(:)
   contains
     procedure          :: create                                => fe_operator_create
     procedure          :: free                                  => fe_operator_free
     procedure          :: set_evaluation_point                  => fe_operator_set_evaluation_point  
     procedure          :: reallocate_after_remesh               => fe_operator_reallocate_after_remesh
                                                                
     procedure          :: compute_tangent                       => fe_operator_compute_tangent
     procedure          :: compute_residual                      => fe_operator_compute_residual
     procedure, private :: compute_internal_residual             => fe_operator_compute_internal_residual
     procedure          :: force_compute                         => fe_operator_force_compute
                                                                
     procedure          :: apply                                 => fe_operator_apply
     procedure          :: apply_add                             => fe_operator_apply_add
     procedure          :: is_linear                             => fe_operator_is_linear
                                                                
     procedure          :: get_tangent                           => fe_operator_get_tangent 
     procedure          :: get_matrix                            => fe_operator_get_matrix
     procedure          :: get_translation                       => fe_operator_get_translation
     procedure          :: get_fe_space                          => fe_operator_get_fe_space
     procedure          :: get_trial_fe_space                    => fe_operator_get_trial_fe_space
     procedure          :: get_discrete_integration              => fe_operator_get_discrete_integration
     procedure          :: get_domain_vector_space               => fe_operator_get_domain_vector_space
     procedure          :: get_range_vector_space                => fe_operator_get_range_vector_space
     procedure          :: get_environment                       => fe_operator_get_environment
     
     procedure          :: get_diagonal_blocks_symmetric_storage => fe_operator_get_diagonal_blocks_symmetric_storage
     procedure          :: get_diagonal_blocks_symmetric         => fe_operator_get_diagonal_blocks_symmetric
     procedure          :: get_diagonal_blocks_sign              => fe_operator_get_diagonal_blocks_sign
     procedure          :: get_sparse_matrix_storage_format      => fe_operator_get_sparse_matrix_storage_format
     procedure          :: get_state                             => fe_operator_get_state
     procedure          :: get_current_evaluation_point          => fe_operator_get_current_evaluation_point
     procedure          :: get_assembler                         => fe_operator_get_assembler
     
     procedure          :: allocate_state                        => fe_operator_allocate_state
     procedure          :: deallocate_state                      => fe_operator_deallocate_state
     procedure          :: set_state                             => fe_operator_set_state
     
     
     procedure          :: abort_if_not_in_range                 => fe_operator_abort_if_not_in_range
     procedure          :: abort_if_not_in_domain                => fe_operator_abort_if_not_in_domain
                                                               
     procedure, private :: create_serial_assembler               => fe_operator_create_serial_assembler
     procedure, private :: create_par_assembler                  => fe_operator_create_par_assembler
     procedure, private :: create_vector_spaces                  => fe_operator_create_vector_spaces
  end type fe_operator_t
  
type, extends(fe_operator_t):: fe_affine_operator_t
  private
 
contains
  procedure          :: compute                     => fe_affine_operator_compute
  procedure          :: compute_residual            => fe_affine_operator_compute_residual
  procedure          :: force_compute               => fe_affine_operator_force_compute
  procedure          :: set_evaluation_point        => fe_affine_operator_set_evaluation_point  
  procedure          :: is_linear                   => fe_affine_operator_is_linear
end type fe_affine_operator_t

  ! Types
  public :: fe_operator_t, fe_affine_operator_t

contains

  
#include "sbm_fe_operator.i90"
#include "sbm_fe_affine_operator.i90"
  
end module
