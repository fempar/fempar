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
module bomap_names
  use types
  use memor
  implicit none
  private

  type bomap
     integer(ip)            :: ndime = 2
     integer(ip)            :: ndimb = 1
     integer(ip)            :: nlocs = 1      ! Evaluation (usually Gauss) points
     real(rp), allocatable  :: clocs(:,:)     ! Coordinates of evaluation points (ndime,nlocs)
     real(rp), allocatable  :: baloc(:,:,:)   ! Local base (ndime-1,ndime-1,nlocs)
     real(rp), allocatable  :: detjm(:)       ! Map Jacobian det     (nlocs)
  end type bomap

  ! Types
  public :: bomap

  ! Functions
  public :: bomap_create, bomap_free

contains
  !==============================================================================
  ! 
  !==============================================================================
  subroutine bomap_create(ndime,nlocs,map)
    implicit none
    integer(ip), intent(in)  :: ndime,nlocs
    type(bomap), intent(out) :: map
    integer(ip) :: iloc

    map%ndime=ndime
    map%ndimb=ndime-1
    map%nlocs=nlocs

    call memalloc(ndime,nlocs,map%clocs,__FILE__,__LINE__)
    call memalloc(ndime,ndime,nlocs,map%baloc,__FILE__,__LINE__)
    call memalloc(nlocs,map%detjm,__FILE__,__LINE__)

  end subroutine bomap_create

  !==============================================================================
  subroutine bomap_free(map)
    implicit none
    type(bomap), intent(inout) :: map

    call memfree(map%clocs,__FILE__,__LINE__)
    call memfree(map%baloc,__FILE__,__LINE__)
    call memfree(map%detjm,__FILE__,__LINE__)

  end subroutine bomap_free

end module bomap_names
