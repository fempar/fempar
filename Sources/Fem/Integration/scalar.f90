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
module scalar_names
  use types_names
  use integrable_names

  type, extends(integrable_t) :: scalar_t
     private
     real(rp) :: a     ! Scalar
  contains
    procedure :: free => scalar_free
    procedure :: init => scalar_init
    procedure :: sum  => scalar_sum
    procedure :: get  => scalar_get
  end type scalar_t

contains

  !==================================================================================================
  subroutine scalar_free (this)
    !-----------------------------------------------------------------------------------------------!
    !   Dummy subroutine to specialize deferred procedure.                                          !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(scalar_t), intent(inout) :: this
  end subroutine scalar_free

  !==================================================================================================
  subroutine scalar_init (this)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine initialize the scalar.                                                      !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(scalar_t), intent(inout) :: this
    
    this%a = 0.0_rp

  end subroutine scalar_init

  !==================================================================================================
  subroutine scalar_sum (this,b)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine sum a real to the scalar.                                                   !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(scalar_t), intent(inout) :: this
    real(rp)       , intent(in)    :: b
    
    this%a = this%a + b

  end subroutine scalar_sum

  !==================================================================================================
  function scalar_get (this)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine returns the scalar value.                                                   !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(scalar_t), intent(inout) :: this
    real(rp)                       :: scalar_get
    
    scalar_get = this%a

  end function scalar_get

end module scalar_names
  
