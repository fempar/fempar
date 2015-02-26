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
module array_class
  use types
  use memor
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

  type array_ip2
     integer(ip)               :: nd1, nd2
     integer(ip), allocatable  :: a(:,:)
  end type array_ip2

  type array_rp1
     integer(ip)               :: nd1
     real(rp)    , allocatable :: a(:) ! Simple real 2D array
  end type array_rp1

  type array_rp2
     integer(ip)               :: nd1, nd2
     real(rp)    , allocatable :: a(:,:) ! Simple real 2D array
  end type array_rp2

  type array_rp3
     integer(ip)               :: nd1, nd2,nd3
     real(rp)    , allocatable :: a(:,:,:) ! Simple real 2D array
  end type array_rp3

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
  public :: array_ip1, array_ip2, array_rp1, array_rp2, array_rp3
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

  ! Memalloc interfaces
  interface memalloc
     module procedure memalloc_array_rp11,memalloc_array_rp12,memalloc_array_rp13,memalloc_array_rp14
     module procedure memalloc_array_rp21,memalloc_array_rp22,memalloc_array_rp23,memalloc_array_rp24
     module procedure memalloc_array_rp31,memalloc_array_rp32,memalloc_array_rp33,memalloc_array_rp34
     module procedure memalloc_array_ip11,memalloc_array_ip12,memalloc_array_ip13,memalloc_array_ip14
     module procedure memalloc_array_ip21,memalloc_array_ip22,memalloc_array_ip23,memalloc_array_ip24
  end interface memalloc
  interface memrealloc
     module procedure memrealloc_array_rp11,memrealloc_array_rp12,memrealloc_array_rp13,            &
          &           memrealloc_array_rp14
     module procedure memrealloc_array_rp21,memrealloc_array_rp22,memrealloc_array_rp23,            &
          &           memrealloc_array_rp24
     module procedure memrealloc_array_rp31,memrealloc_array_rp32,memrealloc_array_rp33,            &
          &           memrealloc_array_rp34
     module procedure memrealloc_array_ip11,memrealloc_array_ip12,memrealloc_array_ip13,            &
          &           memrealloc_array_ip14
     module procedure memrealloc_array_ip21,memrealloc_array_ip22,memrealloc_array_ip23,            &
          &           memrealloc_array_ip24
  end interface memrealloc
  interface memfree
     module procedure memfree_array_rp11,memfree_array_rp12,memfree_array_rp13,memfree_array_rp14
     module procedure memfree_array_rp21,memfree_array_rp22,memfree_array_rp23,memfree_array_rp24
     module procedure memfree_array_rp31,memfree_array_rp32,memfree_array_rp33,memfree_array_rp34
     module procedure memfree_array_ip11,memfree_array_ip12,memfree_array_ip13,memfree_array_ip14
     module procedure memfree_array_ip21,memfree_array_ip22,memfree_array_ip23,memfree_array_ip24
  end interface memfree
  interface memmovealloc
     module procedure memmovealloc_array_rp11,memmovealloc_array_rp12,memmovealloc_array_rp13,    &
          &           memmovealloc_array_rp14
     module procedure memmovealloc_array_rp21,memmovealloc_array_rp22,memmovealloc_array_rp23,    &
          &           memmovealloc_array_rp24
     module procedure memmovealloc_array_rp31,memmovealloc_array_rp32,memmovealloc_array_rp33,    &
          &           memmovealloc_array_rp34
     module procedure memmovealloc_array_ip11,memmovealloc_array_ip12,memmovealloc_array_ip13,    &
          &           memmovealloc_array_ip14
     module procedure memmovealloc_array_ip21,memmovealloc_array_ip22,memmovealloc_array_ip23,    &
          &           memmovealloc_array_ip24
  end interface memmovealloc

  ! Functions
  public :: array_create, array_free, memalloc, memrealloc, memfree, memmovealloc

contains 

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(array_ip1)
  !***********************************************************************
  !***********************************************************************
# define var_type type(array_ip1)
# define var_size 52
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_array_ip11    
# define generic_memalloc_2      memalloc_array_ip12    
# define generic_memalloc_3      memalloc_array_ip13    
# define generic_memalloc_4      memalloc_array_ip14    
# define generic_memrealloc_1    memrealloc_array_ip11  
# define generic_memrealloc_2    memrealloc_array_ip12  
# define generic_memrealloc_3    memrealloc_array_ip13  
# define generic_memrealloc_4    memrealloc_array_ip14  
# define generic_memfree_1       memfree_array_ip11     
# define generic_memfree_2       memfree_array_ip12     
# define generic_memfree_3       memfree_array_ip13     
# define generic_memfree_4       memfree_array_ip14     
# define generic_memmovealloc_1  memmovealloc_array_ip11
# define generic_memmovealloc_2  memmovealloc_array_ip12
# define generic_memmovealloc_3  memmovealloc_array_ip13
# define generic_memmovealloc_4  memmovealloc_array_ip14
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(array_ip2)
  !***********************************************************************
  !***********************************************************************
# define var_type type(array_ip2)
# define var_size 52
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_array_ip21    
# define generic_memalloc_2      memalloc_array_ip22    
# define generic_memalloc_3      memalloc_array_ip23    
# define generic_memalloc_4      memalloc_array_ip24    
# define generic_memrealloc_1    memrealloc_array_ip21  
# define generic_memrealloc_2    memrealloc_array_ip22  
# define generic_memrealloc_3    memrealloc_array_ip23  
# define generic_memrealloc_4    memrealloc_array_ip24  
# define generic_memfree_1       memfree_array_ip21     
# define generic_memfree_2       memfree_array_ip22     
# define generic_memfree_3       memfree_array_ip23     
# define generic_memfree_4       memfree_array_ip24     
# define generic_memmovealloc_1  memmovealloc_array_ip21
# define generic_memmovealloc_2  memmovealloc_array_ip22
# define generic_memmovealloc_3  memmovealloc_array_ip23
# define generic_memmovealloc_4  memmovealloc_array_ip24
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(array_rp2)
  !***********************************************************************
  !***********************************************************************
# define var_type type(array_rp1)
# define var_size 52
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_array_rp11    
# define generic_memalloc_2      memalloc_array_rp12    
# define generic_memalloc_3      memalloc_array_rp13    
# define generic_memalloc_4      memalloc_array_rp14    
# define generic_memrealloc_1    memrealloc_array_rp11  
# define generic_memrealloc_2    memrealloc_array_rp12  
# define generic_memrealloc_3    memrealloc_array_rp13  
# define generic_memrealloc_4    memrealloc_array_rp14  
# define generic_memfree_1       memfree_array_rp11     
# define generic_memfree_2       memfree_array_rp12     
# define generic_memfree_3       memfree_array_rp13     
# define generic_memfree_4       memfree_array_rp14     
# define generic_memmovealloc_1  memmovealloc_array_rp11
# define generic_memmovealloc_2  memmovealloc_array_rp12
# define generic_memmovealloc_3  memmovealloc_array_rp13
# define generic_memmovealloc_4  memmovealloc_array_rp14
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(array_rp2)
  !***********************************************************************
  !***********************************************************************
# define var_type type(array_rp2)
# define var_size 52
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_array_rp21    
# define generic_memalloc_2      memalloc_array_rp22    
# define generic_memalloc_3      memalloc_array_rp23    
# define generic_memalloc_4      memalloc_array_rp24    
# define generic_memrealloc_1    memrealloc_array_rp21  
# define generic_memrealloc_2    memrealloc_array_rp22  
# define generic_memrealloc_3    memrealloc_array_rp23  
# define generic_memrealloc_4    memrealloc_array_rp24  
# define generic_memfree_1       memfree_array_rp21     
# define generic_memfree_2       memfree_array_rp22     
# define generic_memfree_3       memfree_array_rp23     
# define generic_memfree_4       memfree_array_rp24     
# define generic_memmovealloc_1  memmovealloc_array_rp21
# define generic_memmovealloc_2  memmovealloc_array_rp22
# define generic_memmovealloc_3  memmovealloc_array_rp23
# define generic_memmovealloc_4  memmovealloc_array_rp24
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! Specialization to allocate type(array_rp3)
  !***********************************************************************
  !***********************************************************************
# define var_type type(array_rp3)
# define var_size 52
# define var_attr allocatable,target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloc_array_rp31    
# define generic_memalloc_2      memalloc_array_rp32    
# define generic_memalloc_3      memalloc_array_rp33    
# define generic_memalloc_4      memalloc_array_rp34    
# define generic_memrealloc_1    memrealloc_array_rp31  
# define generic_memrealloc_2    memrealloc_array_rp32  
# define generic_memrealloc_3    memrealloc_array_rp33  
# define generic_memrealloc_4    memrealloc_array_rp34  
# define generic_memfree_1       memfree_array_rp31     
# define generic_memfree_2       memfree_array_rp32     
# define generic_memfree_3       memfree_array_rp33     
# define generic_memfree_4       memfree_array_rp34     
# define generic_memmovealloc_1  memmovealloc_array_rp31
# define generic_memmovealloc_2  memmovealloc_array_rp32
# define generic_memmovealloc_3  memmovealloc_array_rp33
# define generic_memmovealloc_4  memmovealloc_array_rp34
# include "memor.i90"

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

end module array_class
