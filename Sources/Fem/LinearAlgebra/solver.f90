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
module matrix_preconditioner_vector_solver_names
  ! Serial modules
  use types_names
  use solver_base_names

  ! Specialize generic_mat data structure and associated methods
  use matrix_names , only : &
       &  generic_mat => matrix_t, generic_info => matrix_info, &
       &  generic_matvec => matvec

  ! Specialize generic_pre data structure and associated methods
  use preconditioner_names, only :  &
       &  generic_pre    => preconditioner_t,   generic_preconditioner => preconditioner_apply, &
       &  generic_bcast  => preconditioner_bcast, generic_fine_task  => preconditioner_fine_task   

  ! Specialize generic_vec data structure and associated methods
  use vector_names , only :  &
       &  generic_vec    => vector_t,       generic_dot     => vector_dot,    &    
       &  generic_copy   => vector_copy,  generic_zero    => vector_zero,   & 
       &  generic_scale  => vector_scale, generic_mxpy    => vector_mxpy,   & 
       &  generic_axpy   => vector_axpy,  generic_aypx    => vector_aypx,   & 
       &  generic_pxpy   => vector_pxpy,  generic_pxmy    => vector_pxmy,   &
       &  generic_nrm2   => vector_nrm2,  generic_clone   => vector_clone,  &
       &  generic_free   => vector_free,  generic_comm    => vector_comm
  
  ! Specialize generic_krylov_basis data structure and associated methods
  use vector_krylov_basis_names, only :  &
       &  generic_krylov_basis              => vector_krylov_basis_t            , &
       &  generic_krylov_basis_alloc        => vector_krylov_basis_alloc        , & 
       &  generic_krylov_basis_free         => vector_krylov_basis_free         , &
       &  generic_krylov_basis_extract_view => vector_krylov_basis_extract_view , &
       &  generic_krylov_basis_multidot     => vector_krylov_basis_multidot     , &
       &  generic_krylov_basis_multiaxpy    => vector_krylov_basis_multiaxpy       

  implicit none
  private
  public :: generic_solve

# include "debug.i90"

contains

# include "solver.i90"

end module matrix_preconditioner_vector_solver_names

!=============================================================================
module block_matrix_block_preconditioner_block_vector_solver_names
  ! Serial modules
  use types_names
  use solver_base_names

  ! Specialize generic_mat data structure and associated methods
  use block_matrix_names , only : &
       &  generic_mat => block_matrix_t, generic_info => block_matrix_info

  use block_matrix_vector_names, only : generic_matvec => block_matvec  

  ! Specialize generic_pre data structure and associated methods
  use block_preconditioner_names, only :  &
       &  generic_pre    => block_preconditioner_t,   generic_preconditioner => block_preconditioner_apply, &
       &  generic_bcast  => block_preconditioner_bcast, generic_fine_task  => block_preconditioner_fine_task 
  
  ! Specialize generic_vec data structure and associated methods 
  use block_vector_names , only :  &
       &  generic_vec   => block_vector_t,       generic_dot   => block_vector_dot,   &    
       &  generic_copy  => block_vector_copy,  generic_zero  => block_vector_zero,  & 
       &  generic_scale => block_vector_scale, generic_mxpy  => block_vector_mxpy,  & 
       &  generic_axpy  => block_vector_axpy,  generic_aypx  => block_vector_aypx,  & 
       &  generic_pxpy  => block_vector_pxpy,  generic_pxmy  => block_vector_pxmy,  &
       &  generic_nrm2  => block_vector_nrm2,  generic_clone => block_vector_clone, &
       &  generic_free  => block_vector_free,  generic_comm  => block_vector_comm
  
  ! Specialize generic_krylov_basis data structure and associated methods
  use block_vector_krylov_basis_names, only :  &
       &  generic_krylov_basis              => block_vector_krylov_basis_t              , &
       &  generic_krylov_basis_alloc        => block_vector_krylov_basis_alloc        , & 
       &  generic_krylov_basis_free         => block_vector_krylov_basis_free         , &
       &  generic_krylov_basis_extract_view => block_vector_krylov_basis_extract_view , &
       &  generic_krylov_basis_multidot     => block_vector_krylov_basis_multidot     , &
       &  generic_krylov_basis_multiaxpy    => block_vector_krylov_basis_multiaxpy       

  implicit none
  private
  public :: generic_solve

# include "debug.i90"

contains

# include "solver.i90"

end module block_matrix_block_preconditioner_block_vector_solver_names

!=============================================================================
module solver_names
  use solver_base_names
  use matrix_preconditioner_vector_solver_names, only: & 
       &  matrix_preconditioner_vector_solve => generic_solve

  use block_matrix_block_preconditioner_block_vector_solver_names, only: & 
       &  block_matrix_block_preconditioner_block_vector_solve => generic_solve
  
  implicit none
  private
  
  interface solve
     module procedure matrix_preconditioner_vector_solve, &
          & block_matrix_block_preconditioner_block_vector_solve
  end interface solve
  
  public :: solve
  
end module solver_names
