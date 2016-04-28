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
     !private
     !class(element_id_t), allocatable :: id
   contains
     !procedure (size_interface)  , deferred :: size
     !procedure (pack_interface)  , deferred :: pack
     !procedure (unpack_interface), deferred :: unpack
     procedure  :: size   => migratory_element_size
     procedure  :: pack   => migratory_element_pack
     procedure  :: unpack => migratory_element_unpack
     procedure (free_interface)  , deferred :: free
     procedure (assignment_interface), deferred :: assign
     !procedure :: get_id    => get_migratory_element_id
     !procedure :: create_id => create_migratory_element_id
     generic   :: assignment(=) => assign
  end type migratory_element_t

  abstract interface
     subroutine free_interface(this)
       import :: migratory_element_t
       implicit none
       class(migratory_element_t), intent(inout) :: this
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


contains
  !=============================================================================
  !=============================================================================
  subroutine migratory_element_size(this,n)
    implicit none
    class(migratory_element_t), intent(in)  :: this
    integer(ip)               , intent(out) :: n
    n=0
    write(*,*) 'Function size() of migratory_element_t is not implemented'
    write(*,*) 'It must be implemented in any derived class that invokes it'
    check(.false.)
  end subroutine migratory_element_size
  
  subroutine migratory_element_pack(this,n,buffer)
    implicit none
    class(migratory_element_t), intent(in)  :: this
    integer(ip)               , intent(in)  :: n
    integer(ieep)             , intent(out) :: buffer(n)
    buffer(1) = 0
    write(*,*) 'Function pack() of migratory_element_t is not implemented'
    write(*,*) 'It must be implemented in any derived class that invokes it'
    check(.false.)
  end subroutine migratory_element_pack
  
  subroutine migratory_element_unpack(this,n,buffer)
    implicit none
    class(migratory_element_t), intent(inout) :: this
    integer(ip)               , intent(in)    :: n
    integer(ieep)             , intent(in)    :: buffer(n)
    write(*,*) 'Function unpack() of migratory_element_t is not implemented'
    write(*,*) 'It must be implemented in any derived class that invokes it'
    check(.false.)
  end subroutine migratory_element_unpack
  
  !=============================================================================
  !=============================================================================
  ! function get_migratory_element_id(this) result(id)
  !   implicit none
  !   class(migratory_element_t), target, intent(in) :: this
  !   class(element_id_t)       , pointer    :: id
  !   id => this%id
  ! end function get_migratory_element_id
  ! subroutine create_migratory_element_id(this,id)
  !   implicit none
  !   class(migratory_element_t), intent(inout) :: this
  !   class(element_id_t)                       :: id
  !   if(allocated(this%id)) deallocate(this%id)
  !   allocate(this%id, mold=id)
  ! end subroutine create_migratory_element_id

end module migratory_element_names
