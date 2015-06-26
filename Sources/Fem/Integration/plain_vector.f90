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
module plain_vector_names
  use types_names
  use memor_names
  use integrable_names
  implicit none
# include "debug.i90"

  type, extends(integrable_t) :: plain_vector_t
     integer(ip)           :: neq      ! Number of equations
     real(rp), allocatable :: b(:)     ! Vector
  contains
    procedure :: create => plain_vector_create
    procedure :: free   => plain_vector_free
  end type plain_vector_t

  ! Types
  public :: plain_vector_t

contains

  !==================================================================================================
  subroutine plain_vector_create(plain_vector,nd)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine creates a plain vector.                                                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(plain_vector_t), intent(out) :: plain_vector
    integer(ip)          , intent(in)  :: nd
    
    ! Fill plain_vector
    plain_vector%neq = nd
    call memalloc(nd,plain_vector%b,__FILE__,__LINE__)

  end subroutine plain_vector_create

  !==================================================================================================
  subroutine plain_vector_free (this)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates a plain vector.                                                 !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(plain_vector_t), intent(inout) :: this
    
    ! Deallocate plain_vector
    if(allocated(this%b)) call memfree(this%b,__FILE__,__LINE__)
    this%neq = 0

  end subroutine plain_vector_free
  
end module plain_vector_names
