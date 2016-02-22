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
module SB_matrix_array_assembler_names
  use types_names
  use matrix_names
  use array_names
  use SB_assembler_names
  
  implicit none
# include "debug.i90"
  private

  type, abstract, extends(SB_assembler_t) :: SB_matrix_array_assembler_t
     private
     class(matrix_t), pointer :: matrix
     class(array_t) , pointer :: array
   contains
     procedure(matrix_array_assembler_compress_storage_interface), deferred :: compress_storage
     procedure :: set_matrix       => matrix_array_assembler_set_matrix
     procedure :: set_array        => matrix_array_assembler_set_array
     procedure :: allocate         => matrix_array_assembler_allocate
     procedure :: free_in_stages   => matrix_array_assembler_free_in_stages
     procedure :: get_matrix       => matrix_array_assembler_get_matrix
     procedure :: get_array        => matrix_array_assembler_get_array
  end type SB_matrix_array_assembler_t
		
		abstract interface
		   subroutine matrix_array_assembler_compress_storage_interface ( this, sparse_matrix_storage_format )
       import :: SB_matrix_array_assembler_t
       implicit none 
       class(SB_matrix_array_assembler_t), intent(inout) :: this 
       character(*)                      , intent(in)    :: sparse_matrix_storage_format 
     end subroutine matrix_array_assembler_compress_storage_interface
		end interface
		
	 
  ! Data types
  public :: SB_matrix_array_assembler_t
  
contains
  ! Sets the pointer to class(matrix_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine matrix_array_assembler_set_matrix(this,matrix)
    implicit none
    class(SB_matrix_array_assembler_t), intent(inout) :: this
    class(matrix_t), pointer, intent(in) :: matrix
    this%matrix => matrix
  end subroutine matrix_array_assembler_set_matrix
  
  ! Sets the pointer to class(array_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine matrix_array_assembler_set_array(this,array)
    implicit none
    class(SB_matrix_array_assembler_t), intent(inout) :: this
    class(array_t), pointer, intent(in) :: array
    this%array => array
  end subroutine matrix_array_assembler_set_array
  
  subroutine matrix_array_assembler_allocate(this)
    implicit none
    class(SB_matrix_array_assembler_t), intent(inout) :: this
    call this%array%allocate()  
  end subroutine matrix_array_assembler_allocate
  
  subroutine matrix_array_assembler_free_in_stages(this,action)
    implicit none
    class(SB_matrix_array_assembler_t), intent(inout) :: this
    integer(ip)                       , intent(in)    :: action
    call this%matrix%free_in_stages(action)
    call this%array%free_in_stages(action)
    if ( action == free_clean ) then
       deallocate(this%matrix)
       nullify(this%matrix)
       deallocate(this%array)
       nullify(this%array)
    end if
  end subroutine matrix_array_assembler_free_in_stages
  
    function matrix_array_assembler_get_matrix(this)
    implicit none
    class(SB_matrix_array_assembler_t), target, intent(in) :: this
    class(matrix_t), pointer :: matrix_array_assembler_get_matrix
    matrix_array_assembler_get_matrix => this%matrix 
  end function matrix_array_assembler_get_matrix
  
  function matrix_array_assembler_get_array(this)
    implicit none
    class(SB_matrix_array_assembler_t), target, intent(in) :: this
    class(array_t), pointer :: matrix_array_assembler_get_array
    matrix_array_assembler_get_array => this%array 
  end function matrix_array_assembler_get_array
  
end module SB_matrix_array_assembler_names
