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
module serial_scalar_names
  use types_names
  use integrable_names
  implicit none
  private

  type, extends(integrable_t) :: serial_scalar_t
     private
     real(rp) :: value 
  contains
    procedure :: create     => serial_scalar_create
	procedure :: init       => serial_scalar_init
	procedure :: sum        => serial_scalar_sum
    procedure :: reduce     => serial_scalar_reduce
	procedure :: get_value  => serial_scalar_get_value
	procedure :: free       => serial_scalar_free
  end type serial_scalar_t
  
  ! Types
  public :: serial_scalar_t

contains

  !==================================================================================================
  subroutine serial_scalar_create (this)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
  end subroutine serial_scalar_create

  !==================================================================================================
  subroutine serial_scalar_init (this,b)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
    real(rp)              ,  intent(in)   :: b
    this%value = b
  end subroutine serial_scalar_init

  !==================================================================================================
  subroutine serial_scalar_sum (this,b)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
    real(rp)              , intent(in)    :: b
    this%value = this%value + b
  end subroutine serial_scalar_sum
  
  !==================================================================================================
  subroutine serial_scalar_reduce (this)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
  end subroutine serial_scalar_reduce

  !==================================================================================================
  function serial_scalar_get_value (this)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
    real(rp)                              :: serial_scalar_get_value
    serial_scalar_get_value = this%value
  end function serial_scalar_get_value

  !==================================================================================================
  subroutine serial_scalar_free (this)
    implicit none
    class(serial_scalar_t), intent(inout) :: this
  end subroutine serial_scalar_free
  
end module serial_scalar_names
  
