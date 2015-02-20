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
module maps_class
  use types
  use memor
  implicit none
#include "debug.i90"
  private
  !================================================================================================
  ! As maps may have overlap there are 3 types of entries:
  !  i) internal (not replicated own entries ),
  !  b) boundary (own replicated entries) and 
  !  e) external (replicated not own entries)
  ! The number of local entries is
  ! nl=ni+nb+ne
  !================================================================================================
  type map
     integer(ip) :: nl=0 ! Number of local entries
     integer(ip) :: ng=0 ! Number of global entries
     integer(ip) :: ni=0 ! Internal entries
     integer(ip) :: nb=0 ! Boundary entries
     integer(ip) :: ne=0 ! External entries
     integer(ip), allocatable ::  l2g(:) ! Local to global list
  end type map

  type map_igp
     integer(ip)  :: nl=0 ! Number of local entries
     integer(igp) :: ng=0 ! Number of global entries
     integer(ip)  :: ni=0 ! Internal entries
     integer(ip)  :: nb=0 ! Boundary entries
     integer(ip)  :: ne=0 ! External entries
     integer(igp), allocatable ::  l2g(:) ! Local to global list
  end type map_igp

  interface map_alloc
    module procedure map_alloc_ip, map_alloc_igp
  end interface map_alloc
  interface map_free
    module procedure map_free_ip, map_free_igp
  end interface map_free
  interface map_write
    module procedure map_write_ip, map_write_igp
  end interface map_write
  interface map_read
    module procedure map_read_ip, map_read_igp
  end interface map_read

  ! Types
  public :: map, map_igp

  ! Functions
  public :: map_alloc, map_free, map_write, map_read !, map_one_to_zero_indexing

contains

  !================================================================================================
  subroutine map_alloc_ip(nl,ng,mp)
    implicit none
    integer(ip), intent(in)       :: nl,ng

    ! WARNING: the following parameter is set
    ! to inout because there may be other components
    ! of the structure (e.g., ne) that maye be set
    ! in advance to the call of the routine. I strongly
    ! recommend to leave it as inout
    type(map), intent(inout) :: mp

    assert(nl>=0)

    mp%nl=nl
    mp%ng=ng
    call memalloc (mp%nl, mp%l2g, __FILE__,__LINE__)
  end subroutine map_alloc_ip

  !================================================================================================
  subroutine map_free_ip(mp)
    implicit none
    type(map), intent(inout) :: mp

  ! assert(allocated(mp%l2g))
  if (mp%nl > 0) then 
   call memfree (mp%l2g,__FILE__,__LINE__)
  end if

  end subroutine map_free_ip

  !================================================================================================
  subroutine map_write_ip(lunio,mp)
    ! Parameters
    integer  , intent(in) :: lunio
    type(map), intent(in) :: mp

    write ( lunio, '(10i10)' ) mp%nl, mp%ng, mp%ni, mp%nb, mp%ne
    if(mp%nl>0) write ( lunio,'(10i10)') mp%l2g

  end subroutine map_write_ip

  !================================================================================================
  subroutine map_read_ip(lunio,mp)
    ! Parameters
    integer  , intent(in)    :: lunio
    type(map), intent(inout) :: mp

    read ( lunio, '(10i10)' ) mp%nl, mp%ng, mp%ni, mp%nb, mp%ne
    if(mp%nl>0) then
       call map_alloc(mp%nl, mp%ng, mp)
       read ( lunio,'(10i10)') mp%l2g
    else
       call map_alloc(0, 0, mp)
    end if
  end subroutine map_read_ip

  subroutine map_alloc_igp(nl,ng,mp)
    implicit none
    integer(ip) , intent(in) :: nl
    integer(igp), intent(in) :: ng

    ! WARNING: the following parameter is set
    ! to inout because there may be other components
    ! of the structure (e.g., ne) that maye be set
    ! in advance to the call of the routine. I strongly
    ! recommend to leave it as inout
    type(map_igp), intent(inout) :: mp

    assert(nl>=0)

    mp%nl=nl
    mp%ng=ng
    call memalloc (int(mp%nl,igp), mp%l2g, __FILE__,__LINE__)
  end subroutine map_alloc_igp

  subroutine map_free_igp(mp)
    implicit none
    type(map_igp), intent(inout) :: mp

    assert(allocated(mp%l2g))
    call memfree (mp%l2g,__FILE__,__LINE__)
  end subroutine map_free_igp

  !================================================================================================
  subroutine map_write_igp(lunio,mp)
    ! Parameters
    integer  , intent(in) :: lunio
    type(map_igp), intent(in) :: mp

    write ( lunio, '(10i10)' ) mp%nl, mp%ng, mp%ni, mp%nb, mp%ne
    if(mp%nl>0) write ( lunio,'(10i10)') mp%l2g
  end subroutine map_write_igp

  !================================================================================================
  subroutine map_read_igp(lunio,mp)
    ! Parameters
    integer  , intent(in)    :: lunio
    type(map_igp), intent(inout) :: mp

    read ( lunio, '(10i10)' ) mp%nl, mp%ng, mp%ni, mp%nb, mp%ne
    if(mp%nl>0) then
       call map_alloc_igp(mp%nl, mp%ng, mp)
       read ( lunio,'(10i10)') mp%l2g
    end if
  end subroutine map_read_igp

end module maps_class
