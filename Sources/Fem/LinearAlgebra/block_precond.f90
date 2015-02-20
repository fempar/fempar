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
module fem_block_precond_class
  ! Serial modules
  use types
  use fem_block_matrix_class
  use fem_block_vector_class

# include "debug.i90"

  ! ** IMPORTANT NOTE **: This is
  ! just a dummy class to specialize
  ! Krylov methods on fem_block_matrices/vectors. 
  ! In the future, fem_block_precond_class 
  ! could become a useful class as long as 
  ! the set of serial solvers/preconditioners 
  ! becomes "rich enough"
  
  implicit none
  private

  type fem_block_precond
  end type fem_block_precond

  ! Types
  public :: fem_block_precond

  ! Functions
  public :: fem_block_precond_create, fem_block_precond_apply, fem_block_precond_bcast, fem_block_precond_fine_task
  contains

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine fem_block_precond_bcast(prec,conv)
    implicit none
    type(fem_block_precond) , intent(in)      :: prec
    logical                 , intent( inout ) :: conv
  end subroutine fem_block_precond_bcast

  ! Dummy method required to specialize Krylov subspace methods
  ! Needs to be filled with the abs operator machinery.
  function fem_block_precond_fine_task(prec)
    implicit none
    type(fem_block_precond) , intent(in) :: prec
    logical                        :: fem_block_precond_fine_task
    fem_block_precond_fine_task = .true. 
  end function fem_block_precond_fine_task

  !=============================================================================
  subroutine  fem_block_precond_create (f_b_prec)
    implicit none
    ! Parameters
    type(fem_block_precond),  intent(out) :: f_b_prec

    ! Compute preconditioner ... 
  end subroutine fem_block_precond_create

  !=============================================================================
  subroutine fem_block_precond_apply (f_b_mat, f_b_prec, x, y)
    implicit none
    ! Parameters
    type(fem_block_matrix) , intent(in)    :: f_b_mat
    type(fem_block_precond), intent(in)    :: f_b_prec
    type(fem_block_vector) , intent(in)    :: x
    type(fem_block_vector) , intent(inout) :: y

    call fem_block_vector_copy (x, y)
  end subroutine fem_block_precond_apply

end module fem_block_precond_class
