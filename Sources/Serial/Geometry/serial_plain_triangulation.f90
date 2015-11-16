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
module serial_plain_triangulation_names
  use types_names
  use memor_names
  use fe_space_types_names
  use hash_table_names  
  use serial_triangulation_names
  implicit none
  private
# include "debug.i90"

  type, extends(element_iterator_t) :: plain_element_iterator_t
     private
     integer(ip)                                 :: ielem=0
     type(serial_plain_triangulation_t), pointer :: trian
   contains
     procedure :: next      => plain_element_iterator_next
     procedure :: has_next  => plain_element_iterator_has_next
  end type plain_element_iterator_t

  type, extends(serial_triangulation_t) :: serial_plain_triangulation_t
     private
     integer(ip) :: elem_array_len = -1  ! length that the elements array is allocated for. 
     type(element_topology_t), allocatable :: elems(:) ! array of elements in the mesh.
   contains
     procedure :: create_element_iterator
     procedure :: free_element_iterator
     procedure :: create => serial_plain_triangulation_create
     procedure :: free   => serial_plain_triangulation_free
  end type serial_plain_triangulation_t
  
  ! Types
  public :: serial_plain_triangulation_t

contains

  !=============================================================================
  subroutine create_element_iterator(this,iterator)
    implicit none
    class(serial_plain_triangulation_t), target     , intent(in)  :: this
    class(element_iterator_t)          , allocatable, intent(out) :: iterator
    integer(ip) :: istat
    allocate(plain_element_iterator_t :: iterator, stat=istat)
    check(istat==0)
    select type(iterator)
    class is(plain_element_iterator_t)
       iterator%trian => this
    end select
  end subroutine create_element_iterator

  !=============================================================================
  function plain_element_iterator_next (this) result(p)
    implicit none
    class(plain_element_iterator_t), intent(inout) :: this
    class(element_topology_t)      , pointer       :: p
    this%ielem = this%ielem + 1
    assert(this%ielem<=this%trian%elem_array_len)
    p => this%trian%elems(this%ielem)
  end function plain_element_iterator_next

  !=============================================================================
  function plain_element_iterator_has_next(this) result(res)
    implicit none
    class(plain_element_iterator_t), intent(inout) :: this
    logical :: res
    res = ( this%ielem < this%trian%num_elems )
    if(.not.res) this%ielem=0 ! Reset iterator for the next loop
  end function plain_element_iterator_has_next

  !=============================================================================
  subroutine free_element_iterator(this,iterator)
    implicit none
    class(serial_plain_triangulation_t)   , intent(in)    :: this
    class(element_iterator_t), allocatable, intent(inout) :: iterator
    ! No internal memory for this simple iterator.
    deallocate(iterator)
  end subroutine free_element_iterator

  !=============================================================================
  !=============================================================================
  !=============================================================================
  subroutine serial_plain_triangulation_create(trian,len)
    implicit none
    class(serial_plain_triangulation_t), intent(inout) :: trian
    integer(ip)                       , intent(in)    :: len
    integer(ip) :: istat
    trian%elem_array_len = len
    ! Allocate the element structure array 
    allocate(trian%elems(trian%elem_array_len), stat=istat)
    check(istat==0)
    !call trian%serial_triangulation_t%create(len)
  end subroutine serial_plain_triangulation_create

  !=============================================================================
  subroutine serial_plain_triangulation_free(trian)
    implicit none
    class(serial_plain_triangulation_t), intent(inout) :: trian
    integer(ip) :: istat

    !call trian%serial_triangulation_t%free()
    trian%elem_array_len = -1
    deallocate(trian%elems, stat=istat)
    check(istat==0)
  end subroutine serial_plain_triangulation_free

  !=============================================================================

end module serial_plain_triangulation_names
