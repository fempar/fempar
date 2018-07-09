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
module sparse_assembler_names
  use types_names
  use allocatable_array_names

  ! Abstract modules
  use assembler_names
  use matrix_names
  use array_names

  ! Concrete implementations
  use sparse_matrix_names
  use serial_scalar_array_names

  implicit none
# include "debug.i90"
  private

  type, extends(assembler_t) :: sparse_assembler_t
  contains
    procedure :: assembly_array          => sparse_assembler_assembly_array
    procedure :: assembly_matrix         => sparse_assembler_assembly_matrix
    procedure :: allocate                => sparse_assembler_allocate
    procedure :: compress_storage_matrix => sparse_assembler_compress_storage_matrix 
  end type

! Data types
public :: sparse_assembler_t

contains

  subroutine sparse_assembler_assembly_array( this,           & 
                                              num_fields,  &
                                              field_blocks,   &
                                              num_dofs,    &
                                              fe_dofs,       &
                                              elvec )
    implicit none
    class(sparse_assembler_t), intent(inout) :: this
    integer(ip)                           , intent(in)    :: num_fields
    integer(ip)                           , intent(in)    :: field_blocks(num_fields)
    integer(ip)                           , intent(in)    :: num_dofs(num_fields)
    type(i1p_t)                           , intent(in)    :: fe_dofs(num_fields)
    real(rp)                              , intent(in)    :: elvec(:)

    class(array_t) , pointer :: array

    array  => this%get_array()
    select type(array)
      class is(serial_scalar_array_t)
      call element_serial_scalar_array_assembly( array,         &
                                                 num_fields, &
                                                 num_dofs,   &
                                                 fe_dofs,      &
                                                 elvec )
      class default
      check(.false.)
    end select
    
  end subroutine sparse_assembler_assembly_array
  
  subroutine sparse_assembler_assembly_matrix( this,            &
                                                            num_fields,   &
                                                            field_blocks,    &
                                                            field_coupling,  &
                                                            num_row_dofs, &
                                                            num_col_dofs, &
                                                            fe_dofs_row,   &
                                                            fe_dofs_col,   &
                                                            elmat )
    implicit none
    class(sparse_assembler_t), intent(inout) :: this
    integer(ip)                           , intent(in)    :: num_fields
    integer(ip)                           , intent(in)    :: field_blocks(num_fields)
    logical                               , intent(in)    :: field_coupling(num_fields,num_fields)
    integer(ip)                           , intent(in)    :: num_row_dofs(num_fields)
    integer(ip)                           , intent(in)    :: num_col_dofs(num_fields)
    type(i1p_t)                           , intent(in)    :: fe_dofs_row(num_fields)
    type(i1p_t)                           , intent(in)    :: fe_dofs_col(num_fields)
    real(rp)                              , intent(in)    :: elmat(:,:) 

    class(matrix_t), pointer :: matrix

    matrix => this%get_matrix()
    select type(matrix)
      class is(sparse_matrix_t)
      call element_sparse_matrix_assembly( matrix,          &
                                           num_fields,   &
                                           num_row_dofs, &
                                           num_col_dofs, &
                                           fe_dofs_row,   &
                                           fe_dofs_col,   &
                                           field_coupling,  &
                                           elmat )
      class default
      check(.false.)
    end select

  end subroutine sparse_assembler_assembly_matrix

  subroutine sparse_assembler_allocate( this )
    implicit none
    class(sparse_assembler_t), intent(inout) :: this
    class(array_t), pointer :: array
    array=>this%get_array()
    call array%allocate()
  end subroutine sparse_assembler_allocate

  subroutine sparse_assembler_compress_storage_matrix( this, & 
                                                             sparse_matrix_storage_format )
    implicit none
    class(sparse_assembler_t) , intent(inout) :: this
    character(*)                              , intent(in)    :: sparse_matrix_storage_format
    class(matrix_t), pointer :: matrix
    matrix=>this%get_matrix() 
    select type(matrix)
      class is(sparse_matrix_t)
      call matrix%convert(sparse_matrix_storage_format)
      class default
      check(.false.)
    end select
  end subroutine sparse_assembler_compress_storage_matrix 
  
  subroutine element_serial_scalar_array_assembly( array, num_fields, num_dofs, fe_dofs, elvec )
    implicit none
    ! Parameters
    type(serial_scalar_array_t), intent(inout) :: array
    integer(ip)                , intent(in)    :: num_fields
    integer(ip)                , intent(in)    :: num_dofs(num_fields)
    type(i1p_t)                , intent(in)    :: fe_dofs(num_fields)
    real(rp)                   , intent(in)    :: elvec(:) 
    
    integer(ip) :: inode, idof, ielvec, ife_space
    
    ielvec = 0
    do ife_space = 1, num_fields
      call array%add( num_dofs(ife_space), &
                      fe_dofs(ife_space)%p,   &
                      ielvec,                  &
                      elvec )
      ielvec = ielvec + num_dofs(ife_space)
    end do
    
  end subroutine element_serial_scalar_array_assembly

  subroutine element_sparse_matrix_assembly( matrix, num_fields, num_row_dofs,      &
                                             num_col_dofs, fe_dofs_row, fe_dofs_col,  &
                                             field_coupling, elmat )
    implicit none
    ! Parameters
    type(sparse_matrix_t), intent(inout) :: matrix
    integer(ip)          , intent(in)    :: num_fields
    integer(ip)          , intent(in)    :: num_row_dofs(num_fields)
    integer(ip)          , intent(in)    :: num_col_dofs(num_fields)
    type(i1p_t)          , intent(in)    :: fe_dofs_row(num_fields)
    type(i1p_t)          , intent(in)    :: fe_dofs_col(num_fields)
    logical              , intent(in)    :: field_coupling(num_fields,num_fields)
    real(rp)             , intent(in)    :: elmat(:,:) 

    integer(ip) :: ife_space, jfe_space
    integer(ip) :: idof, jdof 
    integer(ip) :: inode, jnode
    integer(ip) :: ielmat, jelmat

    ielmat=0
    do ife_space=1, num_fields
       jelmat=0
       do jfe_space=1, num_fields
          if ((field_coupling(ife_space,jfe_space))) then
             call matrix%insert(num_row_dofs(ife_space),num_col_dofs(jfe_space),                 &
                  &             fe_dofs_row(ife_space)%p,fe_dofs_col(jfe_space)%p,ielmat,jelmat,   &
                  &             elmat)
          end if
          jelmat=jelmat+num_col_dofs(jfe_space)
       end do
       ielmat=ielmat+num_row_dofs(ife_space)
    end do

  end subroutine element_sparse_matrix_assembly

end module sparse_assembler_names

