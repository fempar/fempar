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
module element_id_names
  use types_names
  use memor_names
#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
  private
# include "debug.i90"

  ! Abstract identifier of elements, used to acces information in element based structures.
  ! It is used to acces the hash_element data structure defined in has_element_set_*.i90
  type, abstract :: element_id_t
     private
   contains
     procedure(print_interface)     , deferred :: print
     procedure(assignment_interface), deferred :: assign
     procedure(set_index_interface) , deferred :: set_index
     generic :: assignment(=) => assign
  end type element_id_t

  abstract interface
     subroutine print_interface(this,unit)
       import :: ip, element_id_t
       implicit none
       class(element_id_t), intent(in) :: this
       integer(ip)        , intent(in) :: unit
     end subroutine print_interface
     subroutine assignment_interface(this,that)
       import :: element_id_t
       implicit none
       class(element_id_t), intent(inout) :: this
       class(element_id_t), intent(in)    :: that
     end subroutine assignment_interface
     subroutine set_index_interface(this,id)
       import :: element_id_t, ip
       implicit none
       class(element_id_t), intent(inout) :: this
       integer(ip)        , intent(in)    :: id
     end subroutine set_index_interface
  end interface

  ! Integer identifier
  type, extends(element_id_t) :: ip_element_id_t
     private
     integer(ip) :: ielem  = 0
   contains
     procedure  :: get_index       => ip_element_id_get_index
     procedure  :: set_index       => ip_element_id_set_index
     procedure  :: print           => ip_element_id_print
     procedure  :: assign          => ip_element_id_assign
  end type ip_element_id_t


  ! Unique identifier of elements in a forest of element trees
  ! I implemented an ordering given by (level, tree, morton)
  ! which permits to uniformly fill the hash table (under uniform refinement).
  type, extends(element_id_t) :: forest_element_id_t
     !private
     integer(ip)  :: level  = -1
     integer(ip)  :: tree   = -1
     integer(igp) :: morton = -1
   contains
     procedure  :: get_index       => forest_element_id_get_index
     procedure  :: get_morton      => forest_element_id_get_morton
     procedure  :: get_leveln      => forest_element_id_get_level
     procedure  :: set_index       => forest_element_id_set_index
     procedure  :: set_morton      => forest_element_id_set_morton
     procedure  :: set_leveln      => forest_element_id_set_level
     procedure  :: print           => forest_element_id_print
     procedure  :: assign          => forest_element_id_assign
  end type forest_element_id_t

  ! This constant is used to define a to_int function:
  integer(ip), parameter :: estimated_forest_size = 1e5

  public :: element_id_t, ip_element_id_t, forest_element_id_t

  public :: estimated_forest_size

contains

!# include "mem_body.i90"

  !=============================================================================
  function ip_element_id_get_index(this) result(id)
    implicit none
    class(ip_element_id_t), intent(in) :: this
    integer(igp) :: id
    id = this%ielem
  end function ip_element_id_get_index

  subroutine ip_element_id_set_index(this,id)
    implicit none
    class(ip_element_id_t), intent(inout) :: this
    integer(ip)           , intent(in)    :: id
    this%ielem = id 
  end subroutine ip_element_id_set_index

  subroutine ip_element_id_print(this,unit)
    implicit none
    class(ip_element_id_t), intent(in) :: this
    integer(ip), intent(in)            :: unit
    write(unit,'(a,3(2x,i10))') 'ip element ID: ',this%ielem
  end subroutine ip_element_id_print

  subroutine ip_element_id_assign(this,that)
    implicit none
    class(ip_element_id_t), intent(inout) :: this
    class(element_id_t)       , intent(in)    :: that
    logical :: flag
    select type(that)
    class is(ip_element_id_t)
       this%ielem=that%ielem
    class default
       write(*,*) 'Error calling ip_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end subroutine ip_element_id_assign

  !=============================================================================
  function forest_element_id_get_index(this) result(id)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(ip) :: id
    id = this%tree
  end function forest_element_id_get_index
  function forest_element_id_get_morton(this) result(id)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(igp) :: id
    id = this%morton
  end function forest_element_id_get_morton
  function forest_element_id_get_level(this) result(id)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(igp) :: id
    id = this%level
  end function forest_element_id_get_level

  subroutine forest_element_id_set_index(this,id)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: id
    this%tree = id
  end subroutine forest_element_id_set_index
  subroutine forest_element_id_set_morton(this,id)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    integer(igp)       , intent(in)    :: id
    this%morton = id
  end subroutine forest_element_id_set_morton
  subroutine forest_element_id_set_level(this,id)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    integer(igp)       , intent(in)    :: id
    this%level = id
  end subroutine forest_element_id_set_level

  subroutine forest_element_id_print(this,unit)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(ip), intent(in)                :: unit
    write(unit,'(a,3(2x,i10))') 'Forest element ID (tree, level, morton): ',this%tree, this%level, this%morton
  end subroutine forest_element_id_print

  subroutine forest_element_id_assign(this,that)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    class(element_id_t)       , intent(in)    :: that
    logical :: flag
    select type(that)
    class is(forest_element_id_t)
       this%tree=that%tree
       this%morton=that%morton
       this%level=that%level
    class default
       write(*,*) 'Error calling forest_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end subroutine forest_element_id_assign

end module element_id_names
