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
module base_integrable_operator_names
  use types_names
  use base_operator_names
  implicit none
# include "debug.i90"
  private
  
  ! Abstract integrable operator
  type, abstract, extends(base_operator_t) :: base_integrable_operator_t
   contains
     procedure :: fill_values => integrable_operator_fill_values
     procedure :: free_values => integrable_operator_free_values
     procedure (init_interface), deferred :: init
  end type base_integrable_operator_t

  ! Abstract interfaces
  abstract interface
     ! Initialize integrable operator
     subroutine init_interface(op)
       import :: base_integrable_operator_t
       implicit none
       class(base_integrable_operator_t), intent(inout) :: op
     end subroutine init_interface
  end interface

  ! Types
  public :: base_integrable_operator_t

contains
  
  !==================================================================================================
  subroutine integrable_operator_fill_values(op,stage)
    !-----------------------------------------------------------------------------------------------!
    !    Dummy subroutine to implement a deferred method for base_integrable_operator_t.            !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(base_integrable_operator_t), intent(inout) :: op
    integer(ip), optional            , intent(in)    :: stage
  end subroutine integrable_operator_fill_values
  
  !==================================================================================================
  subroutine integrable_operator_free_values(op)
    !-----------------------------------------------------------------------------------------------!
    !    Dummy subroutine to implement a deferred method for base_integrable_operator_t.            !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(base_integrable_operator_t), intent(inout) :: op
  end subroutine integrable_operator_free_values

end module base_integrable_operator_names
