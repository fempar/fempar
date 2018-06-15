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
! This module should encapsulate all what is needed to solve 
! the following problem: given $b_h \in V_h'$, compute the solution
! $u_h \in V_h$ of the problem $A_h(u_h) = b_h$ (if it does exists).
module fe_solver_names
  use fempar_names
  use solver_names
  use operator_names


# include "debug.i90"
  implicit none
  private

  type, extends(operator_t) :: fe_solver_t       
  private
  class(fe_operator_t),  pointer :: fe_nonlinear_operator => null()
  class(solver_t),                 pointer :: solver => null()

contains
  ! Deferred operator_t TBP
  procedure :: apply                          => fe_solver_apply
  procedure :: apply_add                      => fe_solver_apply_add
  procedure :: is_linear                      => fe_solver_is_linear
  procedure :: free                           => fe_solver_free

  ! Getters
  procedure :: get_fe_nonlinear_operator      => fe_solver_get_fe_nonlinear_operator
  procedure :: get_solver                     => fe_solver_get_solver

  ! Setters 
  procedure :: set_fe_nonlinear_operator      => fe_solver_set_fe_nonlinear_operator
  procedure :: set_solver                     => fe_solver_set_solver


end type fe_solver_t 

public :: fe_solver_t

contains


subroutine fe_solver_apply(this,x,y)
 implicit none
 class(fe_solver_t),       intent(inout) :: this
 class(vector_t),             intent(in) :: x
 class(vector_t),          intent(inout) :: y 
 call this%solver%apply(x,y)         
end subroutine fe_solver_apply


subroutine fe_solver_apply_add(this,x,y)
 implicit none
 class(fe_solver_t),       intent(inout) :: this
 class(vector_t),             intent(in) :: x
 class(vector_t),          intent(inout) :: y 
 call this%solver%apply_add(x,y)   
end subroutine fe_solver_apply_add

function fe_solver_is_linear(this) 
 implicit none
 class(fe_solver_t), intent(in)   :: this
 logical                          :: fe_solver_is_linear
 fe_solver_is_linear = .false.
end function fe_solver_is_linear

function fe_solver_get_fe_nonlinear_operator(this)
 implicit none
 class(fe_solver_t),        intent(inout) :: this
 class(fe_operator_t), pointer  :: fe_solver_get_fe_nonlinear_operator
 fe_solver_get_fe_nonlinear_operator => this%fe_nonlinear_operator
end function fe_solver_get_fe_nonlinear_operator

function fe_solver_get_solver(this)
 implicit none
 class(fe_solver_t),        intent(inout) :: this
 class(solver_t),            pointer      :: fe_solver_get_solver
 fe_solver_get_solver => this%solver
end function fe_solver_get_solver


subroutine fe_solver_set_fe_nonlinear_operator( this, fe_nonlinear_operator )
 implicit none
 class(fe_solver_t),                    intent(inout)    :: this
 class(fe_operator_t),  target,  intent(in)    :: fe_nonlinear_operator  
 this%fe_nonlinear_operator => fe_nonlinear_operator
end subroutine  fe_solver_set_fe_nonlinear_operator

subroutine fe_solver_set_solver ( this, solver )
 implicit none
 class(fe_solver_t),   intent(inout)    :: this
 class(solver_t),  target,  intent(in)  :: solver   
 this%solver => solver
end subroutine  fe_solver_set_solver

subroutine fe_solver_free(this)
 implicit none
 class(fe_solver_t), intent(inout)    :: this
 nullify(this%fe_nonlinear_operator)
 nullify(this%solver)
end subroutine fe_solver_free


end module fe_solver_names
