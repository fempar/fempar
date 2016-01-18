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
module subsets_names
  use types_names
  use memor_names
  use stdio_names
  implicit none
# include "debug.i90"
  
  private

  type subsets_t
     integer(ip)                :: &
          nelem=0                          ! Number of elements
     integer(ip), allocatable   :: &
          list(:)                          ! Subset id for every element
  end type subsets_t

  ! Types
  public :: subsets_t

  ! Methods
  public :: subsets_create, subsets_free

contains

  !===============================================================================================
  subroutine subsets_create(nelem,subset)
    implicit none
    integer(ip)        , intent(in)  :: nelem
    type(subsets_t), intent(out) :: subset

    subset%nelem=nelem

    call memalloc (subset%nelem,subset%list, __FILE__,__LINE__)
    subset%list=0 

    return

  end subroutine subsets_create

  !===============================================================================================
  subroutine subsets_free(subset)
    implicit none
    type(subsets_t), intent(inout) :: subset

    subset%nelem= 0

    call memfree (subset%list,__FILE__,__LINE__)

    return

  end subroutine subsets_free

end module subsets_names
