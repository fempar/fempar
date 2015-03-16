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
!=============================================================================
module fem_matrix_fem_precond_fem_vector_solver
  ! Serial modules
  use types
  use solver_base

  ! Specialize generic_mat data structure and associated methods
  use fem_matrix_names , only : &
       &  generic_mat => fem_matrix, generic_info => fem_matrix_info, &
       &  generic_matvec => fem_matvec

  ! Specialize generic_pre data structure and associated methods
  use fem_precond_names, only :  &
       &  generic_pre    => fem_precond,   generic_precond => fem_precond_apply, &
       &  generic_bcast  => fem_precond_bcast, generic_fine_task  => fem_precond_fine_task   

  ! Specialize generic_vec data structure and associated methods
  use fem_vector_names , only :  &
       &  generic_vec    => fem_vector,       generic_dot     => fem_vector_dot,    &    
       &  generic_copy   => fem_vector_copy,  generic_zero    => fem_vector_zero,   & 
       &  generic_scale  => fem_vector_scale, generic_mxpy    => fem_vector_mxpy,   & 
       &  generic_axpy   => fem_vector_axpy,  generic_aypx    => fem_vector_aypx,   & 
       &  generic_pxpy   => fem_vector_pxpy,  generic_pxmy    => fem_vector_pxmy,   &
       &  generic_nrm2   => fem_vector_nrm2,  generic_clone   => fem_vector_clone,  &
       &  generic_free   => fem_vector_free,  generic_comm    => fem_vector_comm
  
  ! Specialize generic_krylov_basis data structure and associated methods
  use fem_vector_krylov_basis_names, only :  &
       &  generic_krylov_basis              => fem_vector_krylov_basis              , &
       &  generic_krylov_basis_alloc        => fem_vector_krylov_basis_alloc        , & 
       &  generic_krylov_basis_free         => fem_vector_krylov_basis_free         , &
       &  generic_krylov_basis_extract_view => fem_vector_krylov_basis_extract_view , &
       &  generic_krylov_basis_multidot     => fem_vector_krylov_basis_multidot     , &
       &  generic_krylov_basis_multiaxpy    => fem_vector_krylov_basis_multiaxpy       

  implicit none
  private
  public :: generic_solve

# include "debug.i90"

contains

# include "solver.i90"

end module fem_matrix_fem_precond_fem_vector_solver

!=============================================================================
module fem_block_matrix_fem_block_precond_fem_block_vector_solver
  ! Serial modules
  use types
  use solver_base

  ! Specialize generic_mat data structure and associated methods
  use fem_block_matrix_names , only : &
       &  generic_mat => fem_block_matrix, generic_info => fem_block_matrix_info

  use fem_block_matrix_vector, only : generic_matvec => fem_block_matvec  

  ! Specialize generic_pre data structure and associated methods
  use fem_block_precond_names, only :  &
       &  generic_pre    => fem_block_precond,   generic_precond => fem_block_precond_apply, &
       &  generic_bcast  => fem_block_precond_bcast, generic_fine_task  => fem_block_precond_fine_task 
  
  ! Specialize generic_vec data structure and associated methods 
  use fem_block_vector_names , only :  &
       &  generic_vec   => fem_block_vector,       generic_dot   => fem_block_vector_dot,   &    
       &  generic_copy  => fem_block_vector_copy,  generic_zero  => fem_block_vector_zero,  & 
       &  generic_scale => fem_block_vector_scale, generic_mxpy  => fem_block_vector_mxpy,  & 
       &  generic_axpy  => fem_block_vector_axpy,  generic_aypx  => fem_block_vector_aypx,  & 
       &  generic_pxpy  => fem_block_vector_pxpy,  generic_pxmy  => fem_block_vector_pxmy,  &
       &  generic_nrm2  => fem_block_vector_nrm2,  generic_clone => fem_block_vector_clone, &
       &  generic_free  => fem_block_vector_free,  generic_comm  => fem_block_vector_comm
  
  ! Specialize generic_krylov_basis data structure and associated methods
  use fem_block_vector_krylov_basis_names, only :  &
       &  generic_krylov_basis              => fem_block_vector_krylov_basis              , &
       &  generic_krylov_basis_alloc        => fem_block_vector_krylov_basis_alloc        , & 
       &  generic_krylov_basis_free         => fem_block_vector_krylov_basis_free         , &
       &  generic_krylov_basis_extract_view => fem_block_vector_krylov_basis_extract_view , &
       &  generic_krylov_basis_multidot     => fem_block_vector_krylov_basis_multidot     , &
       &  generic_krylov_basis_multiaxpy    => fem_block_vector_krylov_basis_multiaxpy       

  implicit none
  private
  public :: generic_solve

# include "debug.i90"

contains

# include "solver.i90"

end module fem_block_matrix_fem_block_precond_fem_block_vector_solver

!=============================================================================
module solver
  use solver_base
  use fem_matrix_fem_precond_fem_vector_solver, only: & 
       &  fem_matrix_fem_precond_fem_vector_solve => generic_solve

  use fem_block_matrix_fem_block_precond_fem_block_vector_solver, only: & 
       &  fem_block_matrix_fem_block_precond_fem_block_vector_solve => generic_solve
  
  implicit none
  private
  
  interface solve
     module procedure fem_matrix_fem_precond_fem_vector_solve, &
          & fem_block_matrix_fem_block_precond_fem_block_vector_solve
  end interface solve
  
  public :: solve
  
end module solver
