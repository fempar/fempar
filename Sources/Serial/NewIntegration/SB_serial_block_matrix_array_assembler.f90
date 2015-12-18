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
subroutine serial_block_matrix_array_assembler_assembly( this, & 
                                                         number_fe_spaces, &
                                                         elem2dof, &
                                                         field_blocks, &
                                                         number_nodes, &
                                                         field_coupling, &
                                                         elmat, &
                                                         elvec )
 implicit none
 class(SB_serial_block_matrix_array_assembler_t), intent(inout) :: this
 integer(ip)                                    , intent(in)    :: number_fe_spaces
 integer(ip)                                    , intent(in)    :: number_nodes(number_fe_spaces)
 type(i1p_t)                                    , intent(in)    :: elem2dof(number_fe_spaces)
 integer(ip)                                    , intent(in)    :: field_blocks(number_fe_spaces)
 logical                                        , intent(in)    :: field_coupling(number_fe_spaces,number_fe_spaces)
 ! elmat MUST have as many rows/columns as \sum_{i=1}^{number_fe_spaces} number_nodes(i)
 real(rp)                                       , intent(in)    :: elmat(:,:) 
 ! elvec MUST have as many entries as \sum_{i=1}^{number_fe_spaces} number_nodes(i)
 real(rp)                                       , intent(in)    :: elvec(:)  

 class(matrix_t), pointer :: matrix
 class(array_t) , pointer :: array

 matrix => this%get_matrix()
 array  => this%get_array()

 select type(matrix)
    class is(serial_block_matrix_t)
    call element_serial_block_matrix_assembly( matrix,  elem2dof,  elmat, number_fe_spaces, field_blocks, number_nodes, field_coupling )
    class default
    check(.false.)
 end select

 select type(array)
    class is(serial_block_array_t)
    call element_serial_block_array_assembly( array,  elem2dof,  elvec, number_fe_spaces, field_blocks, number_nodes, field_coupling )
    class default
    check(.false.)
 end select
end subroutine serial_block_matrix_array_assembler_assembly

subroutine element_serial_block_matrix_assembly( a, el2dof, elmat, number_fe_spaces, &
     blocks, nodes, blocks_coupling )
  implicit none
  ! Parameters
  type(serial_block_matrix_t), intent(inout) :: a
  real(rp), intent(in) :: elmat(:,:) 
  type(i1p_t), intent(in) :: el2dof(:)
  integer(ip), intent(in) :: blocks(:), number_fe_spaces, nodes(:)
  logical, intent(in) :: blocks_coupling(:,:)

  integer(ip) :: c_i, ifem, iblock, inode, idof
  integer(ip) :: c_j, jfem, jblock, jnode, jdof, k
  type(serial_scalar_matrix_t), pointer :: mat

  c_i = 0
  do ifem = 1, number_fe_spaces
     iblock = blocks(ifem)
     do inode = 1, nodes(ifem)
        idof = el2dof(ifem)%p(inode)
        c_i = c_i + 1
        if ( idof > 0 ) then
           c_j = 0
           do jfem = 1,number_fe_spaces
              jblock = blocks(jfem)
              if ( blocks_coupling(iblock,jblock) ) then
                 mat => a%get_block(iblock,jblock)
                 do jnode = 1, nodes(jfem)
                    jdof = el2dof(jfem)%p(jnode)
                    c_j = c_j + 1
                    !write(*,*) 'ifem,iblock,inode,idof',ifem,iblock,inode,idof
                    !write(*,*) 'jfem,jblock,jnode,jdof',jfem,jblock,jnode,jdof
                    if (  (.not.mat%graph%symmetric_storage) .and. jdof > 0 ) then
                       do k = mat%graph%ia(idof),mat%graph%ia(idof+1)-1
                          if ( mat%graph%ja(k) == jdof ) exit
                       end do
                       assert ( k < mat%graph%ia(idof+1) )
                       mat%a(k) = mat%a(k) + elmat(c_i,c_j)
                    else if ( jdof >= idof ) then 
                       do k = mat%graph%ia(idof),mat%graph%ia(idof+1)-1
                          if ( mat%graph%ja(k) == jdof ) exit
                       end do
                       assert ( k < mat%graph%ia(idof+1) )
                       mat%a(k) = mat%a(k) + elmat(c_i,c_j)
                    end if
                 end do
              else
                 c_j = c_j + nodes(jfem)
              end if
           end do
        end if
     end do
  end do

end subroutine element_serial_block_matrix_assembly

subroutine element_serial_block_array_assembly( a, el2dof, elvec, number_fe_spaces, blocks, nodes, blocks_coupling )
  implicit none
  ! Parameters
  type(serial_block_array_t), intent(inout) :: a
  real(rp), intent(in) :: elvec(:)
  integer(ip), intent(in) :: blocks(:), number_fe_spaces, nodes(:)
  type(i1p_t), intent(in) :: el2dof(:)
  logical, intent(in) :: blocks_coupling(:,:)

  integer(ip) :: c_i, ifem, iblock, inode, idof

  c_i = 0
  do ifem = 1, number_fe_spaces
     iblock = blocks(ifem)
     do inode = 1, nodes(ifem)
        idof = el2dof(ifem)%p(inode) 
        c_i = c_i+1
        if ( idof  > 0 ) then
           a%blocks(iblock)%b(idof) =  a%blocks(iblock)%b(idof) + elvec(c_i)
        end if
     end do
  end do

end subroutine element_serial_block_array_assembly

end module SB_serial_block_matrix_array_assembler_names
