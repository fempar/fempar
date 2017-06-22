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
module fempar_sm_linear_solver_names
  use fempar_names
  !use base_sparse_matrix_names
  implicit none
# include "debug.i90"
  private

  type, extends(iterative_linear_solver_t) :: linear_solver_t
     private
     type(fe_affine_operator_t), pointer :: A => NULL()
     type(mlbddc_t)            , pointer :: M => NULL()
   contains
     procedure :: setup_operators => linear_solver_setup_operators
     procedure :: update          => linear_solver_update
     procedure :: get_A
  end type linear_solver_t

  public :: linear_solver_t

contains

  subroutine linear_solver_setup_operators(this,A,M)
    implicit none
    class(linear_solver_t)            , intent(inout) :: this
    type(fe_affine_operator_t), target, intent(in) :: A
    type(mlbddc_t)            , target, intent(in) :: M
    this%A => A
    this%M => M
    call this%iterative_linear_solver_t%set_operators(A,M)
  end subroutine linear_solver_setup_operators

  subroutine linear_solver_update(this)
    implicit none
    class(linear_solver_t), intent(inout) :: this
    call this%A%numerical_setup()
    call this%M%free_numerical_setup()
    call this%M%numerical_setup()
  end subroutine linear_solver_update

  function get_A(this)
    implicit none
    class(linear_solver_t), target, intent(in) :: this
    type(fe_affine_operator_t), pointer        :: get_A
    get_A => this%A
  end function get_A

end module fempar_sm_linear_solver_names
