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
module serial_triangulation_names
  use types_names
  use memor_names
  use element_id_names
  use migratory_element_names
  use fe_space_types_names
  use hash_table_names  
  use JP_triangulation_names
  implicit none
  private
# include "debug.i90"

  type, extends(JP_triangulation_t) :: serial_triangulation_t
     private
   contains
     procedure :: create => serial_triangulation_create
     procedure :: free   => serial_triangulation_free
     procedure :: create_element_iterator
     procedure :: free_element_iterator
     !procedure :: get_element_topology_pointer
  end type serial_triangulation_t

  ! Types
  public :: serial_triangulation_t

contains

  !=============================================================================
  subroutine create_element_iterator(this,iterator)
    implicit none
    class(serial_triangulation_t), target     , intent(in)  :: this
    class(migratory_element_iterator_t)    , allocatable, intent(out) :: iterator
    call this%element_set%create_iterator(iterator)
  end subroutine create_element_iterator
  !=============================================================================
  subroutine free_element_iterator(this,iterator)
    implicit none
    class(serial_triangulation_t)   , intent(in)    :: this
    class(migratory_element_iterator_t), allocatable, intent(inout) :: iterator
    call this%element_set%free_iterator(iterator)
  end subroutine free_element_iterator

  !=============================================================================
  !=============================================================================
  !=============================================================================
  subroutine serial_triangulation_create(trian,len)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip)                  , intent(in)    :: len
    type(JP_element_topology_t) :: mold

    allocate(plain_migratory_element_set_t :: trian%element_set)
    call trian%element_set%create(len,mold)

    ! Mother class function (not type bounded by standard restriction)
    call JP_triangulation_create(trian,len)

  end subroutine serial_triangulation_create

  !=============================================================================
  subroutine serial_triangulation_free(trian)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip) :: istat

    ! Mother class function (not type bounded by standard restriction)
    call JP_triangulation_free(trian)

    ! Deallocate the element structure array 
    call trian%element_set%free()
  end subroutine serial_triangulation_free

  !=============================================================================

end module serial_triangulation_names
