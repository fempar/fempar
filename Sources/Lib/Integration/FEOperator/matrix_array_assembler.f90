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
module matrix_array_assembler_names
  use types_names
  use matrix_names
  use array_names
  
  implicit none
# include "debug.i90"
  private

  type, abstract :: matrix_array_assembler_t
     private
     class(matrix_t), pointer :: matrix
     class(array_t) , pointer :: array
   contains
     procedure (assembly_array_interface)  , deferred :: assemble_array
     procedure (assembly_matrix_interface) , deferred :: assemble_matrix
     procedure (assembly_interface)        , deferred :: assembly
     procedure (face_assembly_interface)   , deferred :: face_assembly
     procedure (compress_storage_interface), deferred :: compress_storage
     procedure                                        :: set_matrix       => matrix_array_assembler_set_matrix
     procedure                                        :: set_array        => matrix_array_assembler_set_array
     procedure                                        :: allocate_array   => matrix_array_assembler_allocate_array
     procedure                                        :: allocate_matrix  => matrix_array_assembler_allocate_matrix
     procedure                                        :: init_array       => matrix_array_assembler_init_array
     procedure                                        :: init_matrix      => matrix_array_assembler_init_matrix
     procedure                                        :: free_in_stages   => matrix_array_assembler_free_in_stages
     procedure                                        :: get_matrix       => matrix_array_assembler_get_matrix
     procedure                                        :: get_array        => matrix_array_assembler_get_array
  end type matrix_array_assembler_t
		 
  abstract interface
     
     subroutine assembly_array_interface( this,           & 
                                          number_fields,  &
                                          field_blocks,   &
                                          field_coupling, &
                                          number_dofs,    &
                                          cell2dof,       &
                                          elvec )
       import :: matrix_array_assembler_t, rp, ip, i1p_t
       implicit none
       class(matrix_array_assembler_t) , intent(inout) :: this
       integer(ip)                     , intent(in)    :: number_fields
       integer(ip)                     , intent(in)    :: field_blocks(number_fields)
       logical                         , intent(in)    :: field_coupling(number_fields,number_fields)
       integer(ip)                     , intent(in)    :: number_dofs(number_fields)
       type(i1p_t)                     , intent(in)    :: cell2dof(number_fields)
       real(rp)                        , intent(in)    :: elvec(:)
     end subroutine assembly_array_interface
     
     subroutine assembly_matrix_interface( this,            &
                                           number_fields,   &
                                           field_blocks,    &
                                           field_coupling,  &
                                           number_row_dofs, &
                                           number_col_dofs, &
                                           cell2row_dofs,   &
                                           cell2col_dofs,   &
                                           elmat )
       import :: matrix_array_assembler_t, rp, ip, i1p_t
       implicit none
       class(matrix_array_assembler_t), intent(inout) :: this
       integer(ip)                    , intent(in)    :: number_fields
       integer(ip)                    , intent(in)    :: field_blocks(number_fields)
       logical                        , intent(in)    :: field_coupling(number_fields,number_fields)
       integer(ip)                    , intent(in)    :: number_row_dofs(number_fields)
       integer(ip)                    , intent(in)    :: number_col_dofs(number_fields)
       type(i1p_t)                    , intent(in)    :: cell2row_dofs(number_fields)
       type(i1p_t)                    , intent(in)    :: cell2col_dofs(number_fields)
       real(rp)                       , intent(in)    :: elmat(:,:) 
     end subroutine assembly_matrix_interface
     
     subroutine assembly_interface( this,             & 
                                    number_fields,    &
                                    number_dofs,     &
                                    elem2dof,         &
                                    field_blocks,     &
                                    field_coupling,   &
                                    elmat,            &
                                    elvec )
       import :: matrix_array_assembler_t, rp, ip, i1p_t
       implicit none
       class(matrix_array_assembler_t) , intent(inout) :: this
       integer(ip)           , intent(in)    :: number_fields
       integer(ip)           , intent(in)    :: number_dofs(number_fields)
       type(i1p_t)           , intent(in)    :: elem2dof(number_fields)
       integer(ip)           , intent(in)    :: field_blocks(number_fields)
       logical               , intent(in)    :: field_coupling(number_fields,number_fields)
       ! elmat MUST have as many rows/columns as \sum_{i=1}^{number_fields} number_dofs(i)
       real(rp)              , intent(in)    :: elmat(:,:) 
       ! elvec MUST have as many entries as \sum_{i=1}^{number_fields} number_dofs(i)
       real(rp)              , intent(in)    :: elvec(:)   
     end subroutine assembly_interface

     subroutine face_assembly_interface(this, &
                                        number_fields, &
                                        test_number_dofs, &
                                        trial_number_dofs, &
                                        test_elem2dof, &
                                        trial_elem2dof, &
                                        field_blocks, &
                                        field_coupling, &
                                        facemat, &
                                        facevec) 
       import :: matrix_array_assembler_t, rp, ip, i1p_t
       implicit none
       class(matrix_array_assembler_t), intent(inout) :: this
       integer(ip)                    , intent(in)    :: number_fields
       integer(ip)                    , intent(in)    :: test_number_dofs(number_fields)
       integer(ip)                    , intent(in)    :: trial_number_dofs(number_fields)
       type(i1p_t)                    , intent(in)    :: test_elem2dof(number_fields)
       type(i1p_t)                    , intent(in)    :: trial_elem2dof(number_fields)
       integer(ip)                    , intent(in)    :: field_blocks(number_fields)
       logical                        , intent(in)    :: field_coupling(number_fields,number_fields)
       real(rp)                       , intent(in)    :: facemat(:,:) 
       real(rp)                       , intent(in)    :: facevec(:)  
     end subroutine face_assembly_interface

     subroutine compress_storage_interface( this, & 
          &                                 sparse_matrix_storage_format )
       import :: matrix_array_assembler_t
       implicit none
       class(matrix_array_assembler_t) , intent(inout) :: this
       character(*)          , intent(in)    :: sparse_matrix_storage_format
     end subroutine compress_storage_interface
    end interface
    
  ! Data types
  public :: matrix_array_assembler_t
  
contains
  ! Sets the pointer to class(matrix_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine matrix_array_assembler_set_matrix(this,matrix)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    class(matrix_t), pointer, intent(in) :: matrix
    this%matrix => matrix
  end subroutine matrix_array_assembler_set_matrix
  
  ! Sets the pointer to class(array_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine matrix_array_assembler_set_array(this,array)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    class(array_t), pointer, intent(in) :: array
    this%array => array
  end subroutine matrix_array_assembler_set_array
  
  subroutine matrix_array_assembler_allocate_array(this)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    call this%array%allocate()  
  end subroutine matrix_array_assembler_allocate_array

  subroutine matrix_array_assembler_allocate_matrix(this)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    call this%matrix%allocate()
  end subroutine matrix_array_assembler_allocate_matrix

  subroutine matrix_array_assembler_init_array(this, value)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    real(rp),                        intent(in)    :: value
    call this%array%init(value)  
  end subroutine matrix_array_assembler_init_array

  subroutine matrix_array_assembler_init_matrix(this, value)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
    real(rp),                        intent(in)    :: value
    call this%matrix%init(value)  
  end subroutine matrix_array_assembler_init_matrix
  
  subroutine matrix_array_assembler_free_in_stages(this,action)
    implicit none
    class(matrix_array_assembler_t), intent(inout) :: this
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
    class(matrix_array_assembler_t), target, intent(in) :: this
    class(matrix_t), pointer :: matrix_array_assembler_get_matrix
    matrix_array_assembler_get_matrix => this%matrix 
  end function matrix_array_assembler_get_matrix
  
  function matrix_array_assembler_get_array(this)
    implicit none
    class(matrix_array_assembler_t), target, intent(in) :: this
    class(array_t), pointer :: matrix_array_assembler_get_array
    matrix_array_assembler_get_array => this%array 
  end function matrix_array_assembler_get_array
  
end module matrix_array_assembler_names
