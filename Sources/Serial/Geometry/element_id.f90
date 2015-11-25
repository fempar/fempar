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
   contains
     procedure(comparison_interface), deferred :: is_smaller_than
     procedure(comparison_interface), deferred :: is_greater_than
     procedure(comparison_interface), deferred :: is_equal_to
     procedure(to_int_interface)    , deferred :: to_int
     procedure(print_interface)     , deferred :: print
     procedure(assignment_interface), deferred :: assign
     generic :: assignment(=) => assign
  end type element_id_t

  abstract interface
     function comparison_interface(this,that) result(flag)
       import :: element_id_t
       implicit none
       class(element_id_t), intent(in) :: this
       class(element_id_t), intent(in) :: that
       logical :: flag
     end function comparison_interface
     function to_int_interface(this) result(key)
       import :: element_id_t, ip
       implicit none
       class(element_id_t), intent(in) :: this
       integer(ip)                     :: key
     end function to_int_interface
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
  end interface

  ! Define memalloc and memfree for element_id_t allocatable
! # define var_attr allocatable, target
! # define point(a,b) call move_alloc(a,b)
! # define generic_status_test             allocated
! # define generic_memalloc_interface      memalloc
! # define generic_memrealloc_interface    memrealloc
! # define generic_memfree_interface       memfree
! # define generic_memmovealloc_interface  memmovealloc
! # define var_type type(element_id_t)
! # define var_size 12
! # define bound_kind ip
! # include "mem_header.i90"
!   public :: memalloc,  memrealloc,  memfree, memmovealloc

  ! Unique identifier of elements in a forest of element trees
  type, extends(element_id_t) :: forest_element_id_t
     integer(ip)  :: tree   = 0
     integer(ip)  :: level  = 0
     integer(igp) :: morton = 0
   contains
     procedure  :: get_tree
     procedure  :: get_morton
     procedure  :: set_tree
     procedure  :: set_morton
     procedure  :: is_smaller_than => forest_element_id_is_smaller_than
     procedure  :: is_greater_than => forest_element_id_is_greater_than
     procedure  :: is_equal_to     => forest_element_id_is_equal_to
     procedure  :: to_int          => forest_element_id_to_int
     procedure  :: print           => forest_element_id_print
     procedure  :: assign          => forest_element_id_assign
  end type forest_element_id_t

  ! This constant is used to define a to_int function:
  ! with the constraint estimated_forest_size*estimated_tree_size < 2**31
  integer(ip), parameter :: estimated_forest_size = 1e5
  integer(ip), parameter :: estimated_tree_size = 2e4

  public :: element_id_t, forest_element_id_t

contains

!# include "mem_body.i90"

  function get_tree(this) result(id)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(ip) :: id
    id = this%tree
  end function get_tree
  function get_morton(this) result(id)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(igp) :: id
    id = this%morton
  end function get_morton

  subroutine set_tree(this,id)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: id
    this%tree = id
  end subroutine set_tree
  subroutine set_morton(this,id)
    implicit none
    class(forest_element_id_t), intent(inout) :: this
    integer(igp)       , intent(in)    :: id
    this%morton = id
  end subroutine set_morton

  function forest_element_id_is_smaller_than(this,that) result(flag)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    class(element_id_t), intent(in)    :: that
    logical :: flag
    select type(that)
    class is(forest_element_id_t)
       flag = .false.
       if(this%tree<that%tree) then
          flag = .true.
       else if(this%tree==that%tree) then
          if(this%morton<that%morton) flag = .true.
       end if
    class default
       write(*,*) 'Error calling forest_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end function forest_element_id_is_smaller_than

  function forest_element_id_is_greater_than(this,that) result(flag)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    class(element_id_t), intent(in) :: that
    logical :: flag
    select type(that)
    class is(forest_element_id_t)
       flag = .false.
       if(this%tree>that%tree) then
          flag = .true.
       else if(this%tree==that%tree) then
          if(this%morton>that%morton) flag = .true.
       end if
    class default
       write(*,*) 'Error calling forest_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end function forest_element_id_is_greater_than

  function forest_element_id_is_equal_to(this,that) result(flag)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    class(element_id_t), intent(in) :: that
    logical :: flag
    select type(that)
    class is(forest_element_id_t)
       if(this%tree==that%tree.and.this%morton==that%morton) then
          flag = .true.
       else 
          flag = .false.
       end if
    class default
       write(*,*) 'Error calling forest_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end function forest_element_id_is_equal_to

  function forest_element_id_to_int(this) result(key)
    implicit none
    class(forest_element_id_t), intent(in) :: this
    integer(ip)                     :: key
    key =  1 + (this%tree-1) * estimated_forest_size + &
         &  abs (int( this%morton - int(this%morton/estimated_tree_size) * estimated_tree_size ) )
  end function forest_element_id_to_int

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
    class default
       write(*,*) 'Error calling forest_element_id_t comparison:'
       write(*,*) 'cannot compare against objects of another class'
       check(.false.)
    end select
  end subroutine forest_element_id_assign

end module element_id_names
