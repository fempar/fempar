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
  use JP_element_topology_names
  use JP_triangulation_names
  implicit none
  private
# include "debug.i90"

  type, extends(JP_triangulation_t) :: serial_triangulation_t
     private
   contains
     procedure :: create  => serial_triangulation_create
     procedure :: free    => serial_triangulation_free
     procedure :: to_dual => serial_triangulation_to_dual
     !procedure :: create_element_iterator
     !procedure :: free_element_iterator
  end type serial_triangulation_t

  ! Types
  public :: serial_triangulation_t

contains

  !=============================================================================
  subroutine serial_triangulation_create(trian,size)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip)                  , intent(in)    :: size
    ! Concrete types to select element and element_id in the set
    type(JP_element_topology_t) :: element_mold

    ! Allocate and create element_set
    allocate(plain_migratory_element_set_t :: trian%element_set)
    call trian%element_set%create(size,element_mold)

    ! Mother class function (not type bounded by standard restriction)
    !call JP_triangulation_create(trian)
    call trian%JP_triangulation_t%create(size)

  end subroutine serial_triangulation_create

  !=============================================================================
  subroutine serial_triangulation_free(trian)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian

    ! Mother class function (not type bounded by standard restriction)
    !call JP_triangulation_free(trian)
    call trian%JP_triangulation_t%free()

    ! Deallocate the element structure array 
    call trian%element_set%free()

  end subroutine serial_triangulation_free

  !=============================================================================
  subroutine serial_triangulation_to_dual(trian)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip) :: iobj,istat

    assert( trian%state == JP_triangulation_elements_filled )

    ! Allocate the vef structure array 
    allocate(vef_topology_t :: trian%vefs(trian%num_vefs), stat=istat)
    check(istat==0)
    do iobj=1, trian%num_vefs
       !call initialize_vef_topology(trian%vefs(iobj))
       call trian%vefs(iobj)%create()
    end do

    ! Invoke mother class function
    call trian%JP_triangulation_t%to_dual()

  end subroutine serial_triangulation_to_dual

end module serial_triangulation_names
