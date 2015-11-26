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
module migratory_element_names
  use types_names
  use hash_table_names
  use element_id_names
  implicit none
  private
# include "debug.i90"

  type, abstract :: migratory_element_t
     private
     class(element_id_t), allocatable :: id
   contains
     procedure (size_interface)  , deferred :: size
     procedure (pack_interface)  , deferred :: pack
     procedure (unpack_interface), deferred :: unpack
     procedure (free_interface)  , deferred :: free
     procedure (assignment_interface), deferred :: assign
     procedure :: get_id    => get_migratory_element_id
     procedure :: create_id => create_migratory_element_id
     generic   :: assignment(=) => assign
  end type migratory_element_t

  abstract interface
     subroutine size_interface(my,n)
       import :: migratory_element_t, ip
       implicit none
       class(migratory_element_t), intent(in)  :: my
       integer(ip)           , intent(out) :: n
     end subroutine size_interface

     subroutine pack_interface(my,n,buffer)
       import :: migratory_element_t, ip, ieep
       implicit none
       class(migratory_element_t), intent(in)  :: my
       integer(ip)             , intent(in)  :: n
       integer(ieep)           , intent(out) :: buffer(n)
     end subroutine pack_interface

     subroutine unpack_interface(my,n,buffer)
       import :: migratory_element_t, ip, ieep
       implicit none
       class(migratory_element_t), intent(inout) :: my
       integer(ip)             , intent(in)    :: n
       integer(ieep)           , intent(in)    :: buffer(n)
     end subroutine unpack_interface

     subroutine free_interface(my)
       import :: migratory_element_t
       implicit none
       class(migratory_element_t), intent(inout) :: my
     end subroutine free_interface

     subroutine assignment_interface(this,that)
       import :: migratory_element_t
       implicit none
       class(migratory_element_t), intent(inout) :: this
       class(migratory_element_t), intent(in)    :: that
     end subroutine assignment_interface

  end interface

  type, abstract :: migratory_element_pointer_t
     class(migratory_element_t), pointer :: p
  end type migratory_element_pointer_t

  public :: migratory_element_t, migratory_element_pointer_t

#define  template_element_t          migratory_element_t
#define  template_element_set_t      migratory_element_set_t
#define  template_element_iterator_t migratory_element_iterator_t
#include "element_set.i90"
  public :: migratory_element_set_t, migratory_element_iterator_t

#define  plain_template_element_set_t      plain_migratory_element_set_t
#define  plain_template_element_iterator_t plain_migratory_element_iterator_t
#include "plain_element_set_header.i90"
  public :: plain_migratory_element_set_t, plain_migratory_element_iterator_t

#define  hash_template_element_set_t      hash_migratory_element_set_t
#define  hash_template_element_iterator_t hash_migratory_element_iterator_t
#include "hash_element_set_header.i90"
  public :: hash_migratory_element_set_t, hash_migratory_element_iterator_t

contains

  function get_migratory_element_id(my) result(id)
    implicit none
    class(migratory_element_t), target, intent(in) :: my
    class(element_id_t)       , pointer    :: id
    id => my%id
  end function get_migratory_element_id
  subroutine create_migratory_element_id(my,id)
    implicit none
    class(migratory_element_t), intent(inout) :: my
    class(element_id_t)                       :: id
    if(allocated(my%id)) deallocate(my%id)
    allocate(my%id, mold=id)
  end subroutine create_migratory_element_id

#include "plain_element_set_body.i90"
#include "hash_element_set_body.i90"

end module migratory_element_names
