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
module serial_block_matrix_array_assembler_names
  use types_names
  use finite_element_names
  use dof_descriptor_names
  use allocatable_array_names

  ! Abstract modules
  use matrix_array_assembler_names
  use matrix_names
  use array_names
  
  ! Concrete implementations
  use serial_scalar_matrix_array_assembler_names
  use serial_scalar_matrix_names
  use serial_block_matrix_names
  use serial_block_array_names
  
  implicit none
# include "debug.i90"
  private

  type, extends(matrix_array_assembler_t) :: serial_block_matrix_array_assembler_t
  contains
	procedure :: assembly => serial_block_matrix_array_assembler_assembly
  end type
	 
  ! Data types
  public :: serial_block_matrix_array_assembler_t
  
contains
  subroutine serial_block_matrix_array_assembler_assembly(this,dof_descriptor,finite_element) 
    implicit none
    class(serial_block_matrix_array_assembler_t), intent(inout) :: this
    type(dof_descriptor_t)                       , intent(in)    :: dof_descriptor
    type(finite_element_t)                       , intent(in)    :: finite_element
	
	class(matrix_t), pointer :: matrix
	class(array_t) , pointer :: array
	
	matrix => this%get_matrix()
	array  => this%get_array()
	
	select type(matrix)
    class is(serial_block_matrix_t)
	   call element_serial_block_matrix_assembly( dof_descriptor, finite_element, matrix )
	 class default
       check(.false.)
    end select  
	
    select type(array)
    class is(serial_block_array_t)
	   call element_serial_block_array_assembly( dof_descriptor, finite_element, array )
	 class default
       check(.false.)
    end select 
  end subroutine serial_block_matrix_array_assembler_assembly
  
  subroutine element_serial_block_matrix_assembly(  dof_descriptor, finite_element,  a ) 
    implicit none
	type(dof_descriptor_t)     , intent(in)    :: dof_descriptor
	type(finite_element_t)     , intent(in)    :: finite_element
    type(serial_block_matrix_t), intent(inout) :: a
	
    integer(ip) :: ivar, iblock, jblock
    type(serial_scalar_matrix_t), pointer :: f_matrix
    do iblock = 1, dof_descriptor%nblocks
       do jblock = 1, dof_descriptor%nblocks
          f_matrix => a%blocks(iblock,jblock)%serial_scalar_matrix
          if ( associated(f_matrix) ) then
             call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, f_matrix, iblock, jblock )
          end if 
       end do
    end do
  end subroutine element_serial_block_matrix_assembly
  
  subroutine element_serial_block_array_assembly(  dof_descriptor, finite_element, a ) 
    implicit none
	type(finite_element_t), intent(in) :: finite_element
    type(dof_descriptor_t), intent(in) :: dof_descriptor
    type(serial_block_array_t), intent(inout) :: a
    integer(ip) :: iblock
    do iblock = 1, dof_descriptor%nblocks
       call element_serial_scalar_array_assembly( dof_descriptor, & 
												  finite_element, & 
												  a%blocks(iblock), &
                                                  iblock )
    end do
  end subroutine element_serial_block_array_assembly
  
  subroutine face_element_serial_block_matrix_assembly(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_pointer_t), intent(in)    :: finite_element(2)
    type(serial_block_matrix_t)  , intent(inout) :: a

    integer(ip) :: iblock, jblock, i
    type(allocatable_array_ip1_t) :: start_aux
    type(serial_scalar_matrix_t), pointer :: f_matrix

    ! Auxiliar local start array to store start(2)
    call allocatable_array_create(dof_descriptor%problems(finite_element(2)%p%problem)%p%nvars+1,start_aux)
    start_aux%a = finite_element(2)%p%start%a

    finite_element(2)%p%start%a = finite_element(2)%p%start%a + finite_element(1)%p%start%a(dof_descriptor%problems(finite_element(1)%p%problem)%p%nvars+1) - 1

    do iblock = 1, dof_descriptor%nblocks
       do jblock = 1, dof_descriptor%nblocks
          f_matrix => a%blocks(iblock,jblock)%serial_scalar_matrix
          if ( associated(f_matrix) ) then
            do i = 1,2
               call element_serial_scalar_matrix_assembly( dof_descriptor, finite_element(i)%p, f_matrix, iblock, jblock )
            end do
            call face_element_serial_scalar_matrix_assembly( dof_descriptor, finite_element, fe_face, f_matrix, iblock, jblock )
          end if
       end do
    end do

    ! Restore start(2)
    finite_element(2)%p%start%a = start_aux%a

  end subroutine face_element_serial_block_matrix_assembly

  subroutine face_element_serial_block_array_assembly(  fe_face, finite_element, dof_descriptor, a ) 
    implicit none
    type(dof_descriptor_t), intent(in)    :: dof_descriptor
    type(fe_face_t)       , intent(in)    :: fe_face
    type(finite_element_t), intent(in)    :: finite_element
    type(serial_block_array_t)  , intent(inout) :: a

    integer(ip) :: iblock

    ! Note: This subroutine only has sense on interface / boundary faces, only
    ! related to ONE element.
    do iblock = 1, dof_descriptor%nblocks
       call face_element_serial_scalar_array_assembly( dof_descriptor, finite_element, fe_face, a%blocks(iblock), iblock )
    end do
  end subroutine face_element_serial_block_array_assembly

end module serial_block_matrix_array_assembler_names
