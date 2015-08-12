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
module par_nonlinear_operator_names
  use base_operand_names
  use abstract_solver_names
  use abstract_environment_names
  use problem_names
  use par_fe_space_names
  implicit none
# include "debug.i90"
  private

  type, abstract :: par_nonlinear_operator_t
     class(base_operand_t), pointer :: b => NULL()  ! Source operand
     class(base_operand_t), pointer :: x => NULL()  ! Result operand
   contains
     procedure (apply_interface), deferred :: apply
  end type par_nonlinear_operator_t

  ! Abstract interfaces
  abstract interface
     subroutine apply_interface(this,sctrl,env,approx,p_fe_space)
       import :: par_nonlinear_operator_t,solver_control_t,abstract_environment_t, &
            &    discrete_integration_pointer_t,par_fe_space_t
       implicit none
       class(par_nonlinear_operator_t), target , intent(inout) :: this
       type(solver_control_t)                  , intent(inout) :: sctrl
       class(abstract_environment_t)           , intent(in)    :: env
       type(discrete_integration_pointer_t)    , intent(inout) :: approx(:)
       type(par_fe_space_t)                    , intent(inout) :: p_fe_space
     end subroutine apply_interface
  end interface

  ! Types
  public :: par_nonlinear_operator_t

end module par_nonlinear_operator_names
