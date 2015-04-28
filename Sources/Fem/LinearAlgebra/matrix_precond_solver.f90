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
module fem_matrix_fem_precond_solver
  ! Serial modules
  use types
  use solver_base
  use serial_environment_names
  use abstract_solver
  use fem_matrix_names
  use fem_vector_names
  use fem_precond_names
  implicit none
  private

  interface solve
     module procedure fem_matrix_fem_precond_fem_vector_solve, &
                      fem_matrix_fem_precond_r1_solve, & 
                      fem_matrix_fem_precond_r2_solve
  end interface solve

  public :: solve

# include "debug.i90"
contains

  subroutine fem_matrix_fem_precond_fem_vector_solve (A,M,b,x,pars)
    implicit none
    type(fem_matrix)    ,intent(in)    :: A     ! Matrix
    type(fem_precond)   ,intent(in)    :: M     ! Preconditioner
    type(fem_vector)    ,intent(in)    :: b     ! RHS
    type(fem_vector)    ,intent(inout) :: x     ! Approximate solution
    type(solver_control),intent(inout) :: pars  ! Solver parameters

    ! Locals
    type(serial_environment) :: senv

    call abstract_solve (A, M, b, x, pars, senv)

  end subroutine fem_matrix_fem_precond_fem_vector_solve

  subroutine fem_matrix_fem_precond_r1_solve (A,M,b,x,pars)
    implicit none
    type(fem_matrix)            ,intent(in)    :: A          ! Matrix
    type(fem_precond)           ,intent(in)    :: M          ! Preconditioner
    real(rp)           , target ,intent(in)    :: b(A%gr%nv) ! RHS
    real(rp)           , target ,intent(inout) :: x(A%gr%nv) ! Approximate solution
    type(solver_control)        ,intent(inout) :: pars       ! Solver parameters

    ! Locals
    type(fem_vector)         :: fem_b
    type(fem_vector)         :: fem_x
    type(serial_environment) :: senv

    ! fill b members
    fem_b%neq     =  A%gr%nv
    fem_b%mode    =  reference
    fem_b%b       => b   

    ! fill x members
    fem_x%neq     = A%gr%nv 
    fem_x%mode    = reference 
    fem_x%b       => x 

    call abstract_solve (A, M, fem_b, fem_x, pars, senv)

  end subroutine fem_matrix_fem_precond_r1_solve

  subroutine fem_matrix_fem_precond_r2_solve (A,M,b,ldb,x,ldx,pars)
    implicit none
    type(fem_matrix)            ,intent(in)    :: A                 ! Matrix
    type(fem_precond)           ,intent(in)    :: M                 ! Preconditioner
    type(solver_control)        ,intent(inout) :: pars              ! Solver parameters
    integer(ip)                 ,intent(in)    :: ldb, ldx
    real(rp)           , target ,intent(in)    :: b(ldb, pars%nrhs) ! RHS
    real(rp)           , target ,intent(inout) :: x(ldx, pars%nrhs) ! Approximate solution

    ! Locals
    type(fem_vector)         :: fem_b
    type(fem_vector)         :: fem_x
    type(serial_environment) :: senv
    integer(ip)              :: k
    integer(ip)              :: tot_its

!!$    if (pars%method == direct) then
!!$       call fem_precond_apply (A, M, pars%nrhs, b, ldb, x, ldx)
!!$    else
    tot_its = 0 
    do k=1, pars%nrhs
       ! fill b members
       fem_b%neq     =  A%gr%nv
       fem_b%mode    =  reference
       fem_b%b       => b(:,k)
       
       ! fill x members
       fem_x%neq     =  A%gr%nv
       fem_x%mode    =  reference 
       fem_x%b       => x(:,k)

       call abstract_solve (A, M, fem_b, fem_x, pars, senv)
       
       tot_its = tot_its + pars%it
    end do
    pars%it = tot_its
!!$    end if
  end subroutine fem_matrix_fem_precond_r2_solve

end module fem_matrix_fem_precond_solver


