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
module nonlinear_operator_names
  use types_names
  use base_operator_names
  use base_operand_names
  implicit none
# include "debug.i90"
  private

  type, abstract :: nonlinear_operator_t
     integer(ip)                     :: max_iter = 1          ! Maximum nonlinear iterations
     real(rp)                        :: nltol    = 1.0e-12_rp ! Nonlinear tolerance
     class(base_operator_t), pointer :: A => NULL()           ! Matrix operator pointer
     class(base_operator_t), pointer :: M => NULL()           ! Preconditioner operator pointer
     class(base_operand_t) , pointer :: b => NULL()           ! RHS operand pointer
     class(base_operand_t) , pointer :: x => NULL()           ! Solution operand pointer
     class(base_operator_t), pointer :: A_int => NULL()       ! (par/block) matrix pointer
     class(base_operand_t) , pointer :: b_int => NULL()       ! (par/block) vector pointer
     class(base_operand_t) , pointer :: x_sol => NULL()       ! (par/block) vector pointer
  end type nonlinear_operator_t

  ! Types
  public :: nonlinear_operator_t

end module nonlinear_operator_names
