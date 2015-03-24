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
module serial_operator_names
  use types
  use base_operator_names
  implicit none

  private

  ! Serial operator (overrides those methods that are 
  ! trivially implemented in the sequential case)
  type, abstract, extends(base_operator) :: serial_operator
   contains
     procedure  :: info => serial_operator_info
     procedure  :: am_i_fine_task => serial_operator_am_i_fine_task
     procedure  :: bcast => serial_operator_bcast
  end type serial_operator

  public :: serial_operator

contains
  subroutine serial_operator_info(op,me,np) 
    implicit none
    class(serial_operator), intent(in)  :: op
    integer(ip)           , intent(out) :: me
    integer(ip)           , intent(out) :: np
    me = 0
    np = 1
  end subroutine serial_operator_info

  function serial_operator_am_i_fine_task(op) 
    implicit none
    class(serial_operator), intent(in)  :: op
    logical                             :: serial_operator_am_i_fine_task
    serial_operator_am_i_fine_task = .true.
  end function serial_operator_am_i_fine_task
  
  subroutine serial_operator_bcast (op, condition)
    implicit none
    class(serial_operator), intent(in) :: op
    logical, intent(inout) :: condition
  end subroutine serial_operator_bcast

end module serial_operator_names
