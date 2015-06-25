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
!***********************************************************************
! All allocatable routines
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc
!***********************************************************************
module array_ip1_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type array_ip1
     integer(ip)               :: nd1
     integer(ip), allocatable  :: a(:)
  end type array_ip1
  public :: array_ip1
# define var_type type(array_ip1)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module array_ip1_names
!=============================================================================
module array_ip2_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type array_ip2
     integer(ip)               :: nd1, nd2
     integer(ip), allocatable  :: a(:,:)
  end type array_ip2
  public :: array_ip2
# define var_type type(array_ip2)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module array_ip2_names
!=============================================================================
module array_rp1_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type array_rp1
     integer(ip)               :: nd1
     real(rp)    , allocatable :: a(:) ! Simple real 2D array
  end type array_rp1
  public :: array_rp1
# define var_type type(array_rp1)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module array_rp1_names
!=============================================================================
module array_rp2_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type array_rp2
     integer(ip)               :: nd1, nd2
     real(rp)    , allocatable :: a(:,:) ! Simple real 2D array
     contains
       procedure :: sum => sum_array_rp2_array_rp2
       generic   :: operator(+) => sum
  end type array_rp2
  public :: array_rp2
# define var_type type(array_rp2)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
  function sum_array_rp2_array_rp2(x,y) result(z)
    implicit none
    class(array_rp2), intent(in) :: x
    type(array_rp2), intent(in) :: y
    type(array_rp2) :: z
    !call x%GuardTemp()
    !call y%GuardTemp()
    !call z%SetTemp()
    assert(size(x%a,1)==size(y%a,1))
    assert(size(x%a,2)==size(y%a,2))
    call memalloc(size(x%a,1),size(x%a,2),z%a,__FILE__,__LINE__)
    z%a = x%a + y%a
    !call x%CleanTemp()
    !call y%CleanTemp()
  end function sum_array_rp2_array_rp2
end module array_rp2_names
!=============================================================================
module array_rp3_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type array_rp3
     integer(ip)               :: nd1, nd2,nd3
     real(rp)    , allocatable :: a(:,:,:) ! Simple real 2D array
  end type array_rp3
  public :: array_rp3
# define var_type type(array_rp3)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module array_rp3_names
!=============================================================================

module array_names
use types_names
use memor_names
  use array_ip1_names
  use array_ip2_names
  use array_rp1_names
  use array_rp2_names
  use array_rp3_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  type array_ip1_pointer
     type(array_ip1)          , pointer :: p => NULL()
  end type array_ip1_pointer

  type array_ip2_pointer
     type(array_ip2)          , pointer :: p => NULL()
  end type array_ip2_pointer

  type array_rp1_pointer
     type(array_rp1)          , pointer :: p => NULL()
  end type array_rp1_pointer

  type array_rp2_pointer
     type(array_rp2)          , pointer :: p => NULL()
  end type array_rp2_pointer

  type array_rp3_pointer
     type(array_rp3)          , pointer :: p => NULL()
  end type array_rp3_pointer

  ! Types
  public :: array_ip1_pointer, array_ip2_pointer, array_rp1_pointer, &
       &    array_rp2_pointer, array_rp3_pointer

  interface array_create
     module procedure array_ip1_create, array_ip2_create, array_rp1_create
     module procedure array_rp2_create, array_rp3_create
  end interface array_create

  interface array_free
     module procedure array_ip1_free, array_ip2_free
     module procedure array_rp1_free, array_rp2_free, array_rp3_free
  end interface array_free

  ! Functions
  public :: array_create, array_free, memalloc, memrealloc, memfree, memmovealloc
  ! public :: memalloc,  memrealloc,  memfree, memmovealloc
  public :: array_ip1, array_ip2, array_rp1, array_rp2, array_rp3

contains 

  !=============================================================================
  subroutine array_ip1_create(nd1,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1
    type(array_ip1), intent(out) :: array

    array%nd1 = nd1

    call memalloc(nd1,array%a,__FILE__,__LINE__)
    array%a = 0
  end subroutine array_ip1_create

  !=============================================================================
  subroutine array_ip2_create(nd1,nd2,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2
    type(array_ip2), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2

    call memalloc(nd1,nd2,array%a,__FILE__,__LINE__)
    array%a = 0
  end subroutine array_ip2_create

  !=============================================================================
  subroutine array_rp1_create(nd1,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1
    type(array_rp1), intent(out) :: array

    array%nd1 = nd1

    call memalloc(nd1,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine array_rp1_create

  !=============================================================================
  subroutine array_rp2_create(nd1,nd2,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2
    type(array_rp2), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2

    call memalloc(nd1,nd2,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine array_rp2_create

  !=============================================================================
  subroutine array_rp3_create(nd1,nd2,nd3,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2, nd3
    type(array_rp3), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2
    array%nd3 = nd3

    call memalloc(nd1,nd2,nd3,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine array_rp3_create

  !=============================================================================
  subroutine array_ip1_free(array)
    implicit none
    type(array_ip1), intent(inout) :: array
    array%nd1 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine array_ip1_free

  !=============================================================================
  subroutine array_ip2_free(array)
    implicit none
    type(array_ip2), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine array_ip2_free

  !=============================================================================
  subroutine array_rp1_free(array)
    implicit none
    type(array_rp1), intent(inout) :: array
    array%nd1 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine array_rp1_free

  !=============================================================================
  subroutine array_rp2_free(array)
    implicit none
    type(array_rp2), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine array_rp2_free

  !=============================================================================
  subroutine array_rp3_free(array)
    implicit none
    type(array_rp3), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    array%nd3 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine array_rp3_free

end module array_names
