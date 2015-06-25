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
module femap_names
use types_names
use memor_names
  implicit none
  private

  type femap
     integer(ip)            :: khes = 0       ! Second derivatives of the map (1=yes 0=no)
     integer(ip)            :: ndime = 2
     integer(ip)            :: nlocs = 1      ! Evaluation (usually Gauss) points
     real(rp), allocatable  :: clocs(:,:)     ! Coordinates of evaluation points (ndime,nlocs)
     real(rp), allocatable  :: jacob(:,:,:)   ! Map Jacobian         (ndime,ndime,nlocs)
     real(rp), allocatable  :: jainv(:,:,:)   ! Map Jacobian inverse (ndime,ndime,nlocs)
     real(rp), allocatable  :: detjm(:)       ! Map Jacobian det     (nlocs)
     real(rp), allocatable  :: hleng(:,:)     ! Local base length    (ndime,nlocs)
     real(rp), allocatable  :: d2sdx(:,:,:,:) ! 2nd derivatives (ndime,ndime,ndime,nlocs)
  end type femap

  ! Types
  public :: femap

  ! Functions
  public :: femap_create, femap_free

contains
  !==============================================================================
  ! 
  !==============================================================================
  subroutine femap_create(khes,ndime,nlocs,map)
    implicit none
    integer(ip), intent(in)  :: khes,ndime,nlocs
    type(femap), intent(out) :: map
    integer(ip) :: iloc

    map%khes=khes
    map%ndime=ndime
    map%nlocs=nlocs

    call memalloc(ndime,nlocs,map%clocs,__FILE__,__LINE__)
    call memalloc(ndime,ndime,nlocs,map%jacob,__FILE__,__LINE__)
    call memalloc(ndime,ndime,nlocs,map%jainv,__FILE__,__LINE__)
    call memalloc(nlocs,map%detjm,__FILE__,__LINE__)
    call memalloc(ndime,nlocs,map%hleng,__FILE__,__LINE__)
    if(khes==1) call memalloc(ndime,ndime,ndime,nlocs,map%d2sdx,  __FILE__,__LINE__)

  end subroutine femap_create

  !==============================================================================
  subroutine femap_free(map)
    implicit none
    type(femap), intent(inout) :: map

    call memfree(map%clocs,__FILE__,__LINE__)
    call memfree(map%jacob,__FILE__,__LINE__)
    call memfree(map%jainv,__FILE__,__LINE__)
    call memfree(map%detjm,__FILE__,__LINE__)
    call memfree(map%hleng,__FILE__,__LINE__)
    if(map%khes==1) call memfree(map%d2sdx,__FILE__,__LINE__)

  end subroutine femap_free

end module femap_names
