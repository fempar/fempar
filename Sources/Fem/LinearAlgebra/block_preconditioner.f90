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
module block_preconditioner_names
  ! Serial modules
use types_names
  use block_matrix_names
  use block_vector_names

# include "debug.i90"

  ! ** IMPORTANT NOTE **: This is
  ! just a dummy class to specialize
  ! Krylov methods on block_matrices/vectors. 
  ! In the future, block_preconditioner_names 
  ! could become a useful class as long as 
  ! the set of serial solvers/preconditioners 
  ! becomes "rich enough"
  
  implicit none
  private

  type block_preconditioner_t
  end type block_preconditioner_t

  ! Types
  public :: block_preconditioner_t

  ! Functions
  public :: block_preconditioner_create, block_preconditioner_apply, block_preconditioner_bcast, block_preconditioner_fine_task
  contains

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine block_preconditioner_bcast(prec,conv)
    implicit none
    type(block_preconditioner_t) , intent(in)      :: prec
    logical                 , intent( inout ) :: conv
  end subroutine block_preconditioner_bcast

  ! Dummy method required to specialize Krylov subspace methods
  ! Needs to be filled with the abs operator machinery.
  function block_preconditioner_fine_task(prec)
    implicit none
    type(block_preconditioner_t) , intent(in) :: prec
    logical                        :: block_preconditioner_fine_task
    block_preconditioner_fine_task = .true. 
  end function block_preconditioner_fine_task

  !=============================================================================
  subroutine  block_preconditioner_create (f_b_prec)
    implicit none
    ! Parameters
    type(block_preconditioner_t),  intent(out) :: f_b_prec

    ! Compute preconditioner ... 
  end subroutine block_preconditioner_create

  !=============================================================================
  subroutine block_preconditioner_apply (f_b_mat, f_b_prec, x, y)
    implicit none
    ! Parameters
    type(block_matrix_t) , intent(in)    :: f_b_mat
    type(block_preconditioner_t), intent(in)    :: f_b_prec
    type(block_vector_t) , intent(in)    :: x
    type(block_vector_t) , intent(inout) :: y

    call block_vector_copy (x, y)
  end subroutine block_preconditioner_apply

end module block_preconditioner_names
