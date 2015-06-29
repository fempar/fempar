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
module matrix_preconditioner_solver_names
  ! Serial modules
  use types_names
  use solver_base_names
  use serial_environment_names
  use abstract_solver_names
  use matrix_names
  use vector_names
  use preconditioner_names
  implicit none
  private

  interface solve
     module procedure matrix_preconditioner_vector_solve, &
                      matrix_preconditioner_r1_solve, & 
                      matrix_preconditioner_r2_solve
  end interface solve

  public :: solve

# include "debug.i90"
contains

  subroutine matrix_preconditioner_vector_solve (A,M,b,x,pars)
    implicit none
    type(matrix_t)    ,intent(in)    :: A     ! Matrix
    type(preconditioner_t)   ,intent(in)    :: M     ! Preconditioner
    type(vector_t)    ,intent(in)    :: b     ! RHS
    type(vector_t)    ,intent(inout) :: x     ! Approximate solution
    type(solver_control_t),intent(inout) :: pars  ! Solver parameters

    ! Locals
    type(serial_environment_t) :: senv

    call abstract_solve (A, M, b, x, pars, senv)

  end subroutine matrix_preconditioner_vector_solve

  subroutine matrix_preconditioner_r1_solve (A,M,b,x,pars)
    implicit none
    type(matrix_t)            ,intent(in)    :: A          ! Matrix
    type(preconditioner_t)           ,intent(in)    :: M          ! Preconditioner
    real(rp)           , target ,intent(in)    :: b(A%gr%nv) ! RHS
    real(rp)           , target ,intent(inout) :: x(A%gr%nv) ! Approximate solution
    type(solver_control_t)        ,intent(inout) :: pars       ! Solver parameters

    ! Locals
    type(vector_t)         :: vector_b
    type(vector_t)         :: vector_x
    type(serial_environment_t) :: senv

    ! fill vector_b members
    vector_b%neq     =  A%gr%nv
    vector_b%mode    =  reference
    vector_b%b => b   

    ! fill vector_x members
    vector_x%neq     = A%gr%nv 
    vector_x%mode    = reference 
    vector_x%b       => x 

    call abstract_solve (A, M, vector_b, vector_x, pars, senv)

  end subroutine matrix_preconditioner_r1_solve

  subroutine matrix_preconditioner_r2_solve (A,M,b,ldb,x,ldx,pars)
    implicit none
    type(matrix_t)            ,intent(in)    :: A                 ! Matrix
    type(preconditioner_t)           ,intent(in)    :: M                 ! Preconditioner
    type(solver_control_t)        ,intent(inout) :: pars              ! Solver parameters
    integer(ip)                 ,intent(in)    :: ldb, ldx
    real(rp)           , target ,intent(in)    :: b(ldb, pars%nrhs) ! RHS
    real(rp)           , target ,intent(inout) :: x(ldx, pars%nrhs) ! Approximate solution

    ! Locals
    type(vector_t)         :: vector_b
    type(vector_t)         :: vector_x
    type(serial_environment_t) :: senv
    integer(ip)              :: k
    integer(ip)              :: tot_its

!!$    if (pars%method == direct) then
!!$       call preconditioner_apply (A, M, pars%nrhs, vector_b, ldb, vector_x, ldx)
!!$    else
    tot_its = 0 
    do k=1, pars%nrhs
       ! fill b members
       vector_b%neq     =  A%gr%nv
       vector_b%mode    =  reference
       vector_b%b       => b(:,k)
       
       ! fill vector_x members
       vector_x%neq     =  A%gr%nv
       vector_x%mode    =  reference 
       vector_x%b       => x(:,k)

       call abstract_solve (A, M, vector_b, vector_x, pars, senv)
       
       tot_its = tot_its + pars%it
    end do
    pars%it = tot_its
!!$    end if
  end subroutine matrix_preconditioner_r2_solve

end module matrix_preconditioner_solver_names


