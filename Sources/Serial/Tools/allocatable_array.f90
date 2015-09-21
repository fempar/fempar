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
module allocatable_array_ip1_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type allocatable_array_ip1_t
     integer(ip)               :: nd1
     integer(ip), allocatable  :: a(:)
  end type allocatable_array_ip1_t
  public :: allocatable_array_ip1_t
# define var_type type(allocatable_array_ip1_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module allocatable_array_ip1_names
!=============================================================================
module allocatable_array_ip2_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type allocatable_array_ip2_t
     integer(ip)               :: nd1, nd2
     integer(ip), allocatable  :: a(:,:)
  end type allocatable_array_ip2_t
  public :: allocatable_array_ip2_t
# define var_type type(allocatable_array_ip2_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module allocatable_array_ip2_names
!=============================================================================
module allocatable_array_rp1_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type allocatable_array_rp1_t
     integer(ip)               :: nd1
     real(rp)    , allocatable :: a(:) ! Simple real 2D array
  end type allocatable_array_rp1_t
  public :: allocatable_array_rp1_t
# define var_type type(allocatable_array_rp1_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module allocatable_array_rp1_names
!=============================================================================
module allocatable_array_rp2_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type allocatable_array_rp2_t
     integer(ip)               :: nd1, nd2
     real(rp)    , allocatable :: a(:,:) ! Simple real 2D array
     contains
       procedure :: sum => sum_allocatable_array_rp2_array_rp2
       generic   :: operator(+) => sum
  end type allocatable_array_rp2_t
  public :: allocatable_array_rp2_t
# define var_type type(allocatable_array_rp2_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
  function sum_allocatable_array_rp2_array_rp2(x,y) result(z)
    implicit none
    class(allocatable_array_rp2_t), intent(in) :: x
    type(allocatable_array_rp2_t), intent(in) :: y
    type(allocatable_array_rp2_t) :: z
    !call x%GuardTemp()
    !call y%GuardTemp()
    !call z%SetTemp()
    assert(size(x%a,1)==size(y%a,1))
    assert(size(x%a,2)==size(y%a,2))
    call memalloc(size(x%a,1),size(x%a,2),z%a,__FILE__,__LINE__)
    z%a = x%a + y%a
    !call x%CleanTemp()
    !call y%CleanTemp()
  end function sum_allocatable_array_rp2_array_rp2
end module allocatable_array_rp2_names
!=============================================================================
module allocatable_array_rp3_names
use types_names
use memor_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  type allocatable_array_rp3_t
     integer(ip)               :: nd1, nd2,nd3
     real(rp)    , allocatable :: a(:,:,:) ! Simple real 2D array
  end type allocatable_array_rp3_t
  public :: allocatable_array_rp3_t
# define var_type type(allocatable_array_rp3_t)
# define var_size 52
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module allocatable_array_rp3_names
!=============================================================================

module allocatable_array_names
use types_names
use memor_names
  use allocatable_array_ip1_names
  use allocatable_array_ip2_names
  use allocatable_array_rp1_names
  use allocatable_array_rp2_names
  use allocatable_array_rp3_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private

  type p_allocatable_array_ip1_t
     type(allocatable_array_ip1_t)          , pointer :: p => NULL()
  end type p_allocatable_array_ip1_t

  type p_allocatable_array_ip2_t
     type(allocatable_array_ip2_t)          , pointer :: p => NULL()
  end type p_allocatable_array_ip2_t

  type p_allocatable_array_rp1_t
     type(allocatable_array_rp1_t)          , pointer :: p => NULL()
  end type p_allocatable_array_rp1_t

  type p_allocatable_array_rp2_t
     type(allocatable_array_rp2_t)          , pointer :: p => NULL()
  end type p_allocatable_array_rp2_t

  type p_allocatable_array_rp3_t
     type(allocatable_array_rp3_t)          , pointer :: p => NULL()
  end type p_allocatable_array_rp3_t

  ! Types
  public :: p_allocatable_array_ip1_t, p_allocatable_array_ip2_t, p_allocatable_array_rp1_t, &
       &    p_allocatable_array_rp2_t, p_allocatable_array_rp3_t

  interface allocatable_array_create
     module procedure allocatable_array_ip1_create, allocatable_array_ip2_create, allocatable_array_rp1_create
     module procedure allocatable_array_rp2_create, allocatable_array_rp3_create
  end interface allocatable_array_create

  interface allocatable_array_free
     module procedure allocatable_array_ip1_free, allocatable_array_ip2_free
     module procedure allocatable_array_rp1_free, allocatable_array_rp2_free, allocatable_array_rp3_free
  end interface allocatable_array_free

  ! Functions
  public :: allocatable_array_create, allocatable_array_free, memalloc, memrealloc, memfree, memmovealloc
  public :: allocatable_array_ip1_t, allocatable_array_ip2_t, allocatable_array_rp1_t, allocatable_array_rp2_t, allocatable_array_rp3_t

contains 

  !=============================================================================
  subroutine allocatable_array_ip1_create(nd1,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1
    type(allocatable_array_ip1_t), intent(out) :: array

    array%nd1 = nd1

    call memalloc(nd1,array%a,__FILE__,__LINE__)
    array%a = 0
  end subroutine allocatable_array_ip1_create

  !=============================================================================
  subroutine allocatable_array_ip2_create(nd1,nd2,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2
    type(allocatable_array_ip2_t), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2

    call memalloc(nd1,nd2,array%a,__FILE__,__LINE__)
    array%a = 0
  end subroutine allocatable_array_ip2_create

  !=============================================================================
  subroutine allocatable_array_rp1_create(nd1,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1
    type(allocatable_array_rp1_t), intent(out) :: array

    array%nd1 = nd1

    call memalloc(nd1,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine allocatable_array_rp1_create

  !=============================================================================
  subroutine allocatable_array_rp2_create(nd1,nd2,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2
    type(allocatable_array_rp2_t), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2

    call memalloc(nd1,nd2,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine allocatable_array_rp2_create

  !=============================================================================
  subroutine allocatable_array_rp3_create(nd1,nd2,nd3,array)
    implicit none
    integer(ip)    , intent(in)  :: nd1, nd2, nd3
    type(allocatable_array_rp3_t), intent(out) :: array

    array%nd1 = nd1
    array%nd2 = nd2
    array%nd3 = nd3

    call memalloc(nd1,nd2,nd3,array%a,__FILE__,__LINE__)
    array%a = 0.0_rp
  end subroutine allocatable_array_rp3_create

  !=============================================================================
  subroutine allocatable_array_ip1_free(array)
    implicit none
    type(allocatable_array_ip1_t), intent(inout) :: array
    array%nd1 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine allocatable_array_ip1_free

  !=============================================================================
  subroutine allocatable_array_ip2_free(array)
    implicit none
    type(allocatable_array_ip2_t), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine allocatable_array_ip2_free

  !=============================================================================
  subroutine allocatable_array_rp1_free(array)
    implicit none
    type(allocatable_array_rp1_t), intent(inout) :: array
    array%nd1 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine allocatable_array_rp1_free

  !=============================================================================
  subroutine allocatable_array_rp2_free(array)
    implicit none
    type(allocatable_array_rp2_t), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine allocatable_array_rp2_free

  !=============================================================================
  subroutine allocatable_array_rp3_free(array)
    implicit none
    type(allocatable_array_rp3_t), intent(inout) :: array
    array%nd1 = 0
    array%nd2 = 0
    array%nd3 = 0
    call memfree(array%a,__FILE__,__LINE__)
  end subroutine allocatable_array_rp3_free

end module allocatable_array_names
