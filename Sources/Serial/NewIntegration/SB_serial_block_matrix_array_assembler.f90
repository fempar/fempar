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
module SB_serial_block_matrix_array_assembler_names
  use types_names
  use dof_descriptor_names
  use allocatable_array_names

  ! Abstract modules
  use SB_matrix_array_assembler_names
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

  type, extends(SB_matrix_array_assembler_t) :: SB_serial_block_matrix_array_assembler_t
contains
  procedure :: assembly => serial_block_matrix_array_assembler_assembly
end type

! Data types
public :: SB_serial_block_matrix_array_assembler_t

contains
subroutine serial_block_matrix_array_assembler_assembly( this, el2dof, elmat, elvec, &
                                                         number_fe_spaces, blocks, nodes ) 
 implicit none
 class(SB_serial_block_matrix_array_assembler_t), intent(inout) :: this
 real(rp), intent(in) :: elmat(:,:), elvec(:) 
 type(i1p_t), intent(in) :: el2dof(:)
 integer(ip), intent(in) :: blocks(:), number_fe_spaces, nodes(:)

 class(matrix_t), pointer :: matrix
 class(array_t) , pointer :: array

 matrix => this%get_matrix()
 array  => this%get_array()

 select type(matrix)
    class is(serial_block_matrix_t)
    call element_serial_block_matrix_assembly( matrix,  el2dof,  elmat, number_fe_spaces, blocks, nodes )
    class default
    check(.false.)
 end select

 select type(array)
    class is(serial_block_array_t)
    call element_serial_block_array_assembly( array,  el2dof,  elvec, number_fe_spaces, blocks, nodes )
    class default
    check(.false.)
 end select
end subroutine serial_block_matrix_array_assembler_assembly

subroutine element_serial_block_matrix_assembly( a, el2dof, elmat, number_fe_spaces, blocks, nodes )
  implicit none
  ! Parameters
  type(serial_block_matrix_t)         , intent(inout) :: a
  real(rp), intent(in) :: elmat(:,:) 
  type(i1p_t), intent(in) :: el2dof(:)
  integer(ip), intent(in) :: blocks(:), number_fe_spaces, nodes(:)

  integer(ip) :: c_i, ifem, iblock, inode, idof
  integer(ip) :: c_j, jfem, jblock, jnode, jdof, k

  c_i = 0
  c_j = 0
  do ifem = 1, number_fe_spaces
     iblock = blocks(ifem)
     do inode = 1, nodes(ifem)
        idof = el2dof(ifem)%l(inode)
        c_i = c_i + 1
        do jfem = 1,number_fe_spaces
           jblock = blocks(jfem)
           do jnode = 1, nodes(jfem)
              jdof = el2dof(jfem)%l(jnode)
              c_j = c_j + 1
              if (  (.not. a%blocks(iblock,jblock)%serial_scalar_matrix%graph%symmetric_storage) .and. jdof > 0 ) then
                 do k = a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof),a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof+1)-1
                    if ( a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ja(k) == jdof ) exit
                 end do
                 assert ( k < a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof+1) )
                 a%blocks(iblock,jblock)%serial_scalar_matrix%a(k) = a%blocks(iblock,jblock)%serial_scalar_matrix%a(k) + elmat(c_i,c_j)
              else if ( jdof >= idof ) then 
                 do k = a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof),a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof+1)-1
                    if ( a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ja(k) == jdof ) exit
                 end do
                 assert ( k < a%blocks(iblock,jblock)%serial_scalar_matrix%graph%ia(idof+1) )
                 a%blocks(iblock,jblock)%serial_scalar_matrix%a(k) = a%blocks(iblock,jblock)%serial_scalar_matrix%a(k) + elmat(c_i,c_j)
              end if
           end do
        end do
     end do
  end do

end subroutine element_serial_block_matrix_assembly

subroutine element_serial_block_array_assembly(  a, el2dof, elvec, number_fe_spaces, blocks, nodes )
  implicit none
  ! Parameters
  type(serial_block_array_t), intent(inout) :: a
  real(rp), intent(in) :: elvec(:)
  integer(ip), intent(in) :: blocks(:), number_fe_spaces, nodes(:)
  type(i1p_t), intent(in) :: el2dof(:)

  integer(ip) :: c_i, ifem, iblock, inode, idof

  c_i = 0
  do ifem = 1, number_fe_spaces
     iblock = blocks(ifem)
     do inode = 1, nodes(ifem)
        idof = el2dof(ifem)%l(inode) 
        c_i = c_i+1
        if ( idof  > 0 ) then
           a%blocks(iblock)%b(idof) =  a%blocks(iblock)%b(idof) + elvec(c_i)
        end if
     end do
  end do

end subroutine element_serial_block_array_assembly

end module SB_serial_block_matrix_array_assembler_names