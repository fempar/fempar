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
module assembler_names
  use types_names
  use matrix_names
  use array_names
  
  implicit none
# include "debug.i90"
  private

  type, abstract :: assembler_t
     private
     class(matrix_t), pointer :: matrix
     class(array_t) , pointer :: array
   contains
     procedure (assembly_array_interface)  , deferred :: assembly_array
     procedure (assembly_matrix_interface) , deferred :: assembly_matrix
     procedure (compress_storage_matrix_interface), deferred :: compress_storage_matrix 
     procedure                                        :: compress_storage_array => assembler_compress_storage_array 
     procedure                                        :: compress_storage       => assembler_compress_storage 
     procedure                                        :: set_matrix             => assembler_set_matrix
     procedure                                        :: set_array              => assembler_set_array
     procedure                                        :: allocate_array         => assembler_allocate_array
     procedure                                        :: allocate_matrix        => assembler_allocate_matrix
     procedure                                        :: init_array             => assembler_init_array
     procedure                                        :: init_matrix            => assembler_init_matrix
     procedure                                        :: free_in_stages         => assembler_free_in_stages
     procedure                                        :: free                   => assembler_free
     procedure                                        :: get_matrix             => assembler_get_matrix
     procedure                                        :: get_array              => assembler_get_array
  end type assembler_t
   
  abstract interface
     
     subroutine assembly_array_interface( this,           & 
                                          num_fields,     &
                                          field_blocks,   &
                                          num_dofs,       &
                                          fe_dofs,        &
                                          elvec )
       import :: assembler_t, rp, ip, i1p_t
       implicit none
       class(assembler_t) , intent(inout) :: this
       integer(ip)                     , intent(in)    :: num_fields
       integer(ip)                     , intent(in)    :: field_blocks(num_fields)
       integer(ip)                     , intent(in)    :: num_dofs(num_fields)
       type(i1p_t)                     , intent(in)    :: fe_dofs(num_fields)
       real(rp)                        , intent(in)    :: elvec(:)
     end subroutine assembly_array_interface
     
     subroutine assembly_matrix_interface( this,            &
                                           num_fields,   &
                                           field_blocks,    &
                                           field_coupling,  &
                                           num_row_dofs, &
                                           num_col_dofs, &
                                           fe_dofs_row,   &
                                           fe_dofs_col,   &
                                           elmat )
       import :: assembler_t, rp, ip, i1p_t
       implicit none
       class(assembler_t), intent(inout) :: this
       integer(ip)                    , intent(in)    :: num_fields
       integer(ip)                    , intent(in)    :: field_blocks(num_fields)
       logical                        , intent(in)    :: field_coupling(num_fields,num_fields)
       integer(ip)                    , intent(in)    :: num_row_dofs(num_fields)
       integer(ip)                    , intent(in)    :: num_col_dofs(num_fields)
       type(i1p_t)                    , intent(in)    :: fe_dofs_row(num_fields)
       type(i1p_t)                    , intent(in)    :: fe_dofs_col(num_fields)
       real(rp)                       , intent(in)    :: elmat(:,:) 
     end subroutine assembly_matrix_interface

     subroutine compress_storage_matrix_interface( this, & 
          &                                 sparse_matrix_storage_format )
       import :: assembler_t
       implicit none
       class(assembler_t) , intent(inout) :: this
       character(*)          , intent(in)    :: sparse_matrix_storage_format
     end subroutine compress_storage_matrix_interface
     
     subroutine compress_storage_array_interface( this )
       import :: assembler_t
       implicit none
       class(assembler_t) , intent(inout) :: this
     end subroutine compress_storage_array_interface
          
    end interface
    
  ! Data types
  public :: assembler_t
  public :: assembler_set_matrix, assembler_set_array
  
contains
  subroutine assembler_compress_storage_array(this)
    implicit none
    class(assembler_t), intent(inout) :: this
    ! Do-nothing unless being overriden 
  end subroutine assembler_compress_storage_array

  subroutine assembler_compress_storage(this, sparse_matrix_storage_format)
    implicit none
    class(assembler_t), intent(inout) :: this
    character(*)      , intent(in)    :: sparse_matrix_storage_format
    call this%compress_storage_matrix( sparse_matrix_storage_format ) 
    call this%compress_storage_array() 
  end subroutine assembler_compress_storage
  
  ! Sets the pointer to class(matrix_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine assembler_set_matrix(this,matrix)
    implicit none
    class(assembler_t), intent(inout) :: this
    class(matrix_t), pointer, intent(in) :: matrix
    this%matrix => matrix
  end subroutine assembler_set_matrix
  
  ! Sets the pointer to class(array_t) in such a way that this 
  ! can become reponsible to free it later on 
  subroutine assembler_set_array(this,array)
    implicit none
    class(assembler_t), intent(inout) :: this
    class(array_t), pointer, intent(in) :: array
    this%array => array
  end subroutine assembler_set_array
  
  subroutine assembler_allocate_array(this)
    implicit none
    class(assembler_t), intent(inout) :: this
    call this%array%allocate()  
  end subroutine assembler_allocate_array

  subroutine assembler_allocate_matrix(this)
    implicit none
    class(assembler_t), intent(inout) :: this
    call this%matrix%allocate()
  end subroutine assembler_allocate_matrix

  subroutine assembler_init_array(this, value)
    implicit none
    class(assembler_t), intent(inout) :: this
    real(rp),                        intent(in)    :: value
    call this%array%init(value)  
  end subroutine assembler_init_array

  subroutine assembler_init_matrix(this, value)
    implicit none
    class(assembler_t), intent(inout) :: this
    real(rp),                        intent(in)    :: value
    call this%matrix%init(value)  
  end subroutine assembler_init_matrix
  
  subroutine assembler_free_in_stages(this,action)
    implicit none
    class(assembler_t), intent(inout) :: this
    integer(ip)                       , intent(in)    :: action
    call this%matrix%free_in_stages(action)
    call this%array%free_in_stages(action)
    if ( action == free_clean ) then
       deallocate(this%matrix)
       nullify(this%matrix)   
       deallocate(this%array)
       nullify(this%array)  
    end if
  end subroutine assembler_free_in_stages
  
  subroutine assembler_free(this)
    implicit none
    class(assembler_t), intent(inout) :: this
    call this%free_in_stages(free_numerical_setup)
    call this%free_in_stages(free_symbolic_setup)
    call this%free_in_stages(free_clean)
  end subroutine assembler_free
  
    function assembler_get_matrix(this)
    implicit none
    class(assembler_t), target, intent(in) :: this
    class(matrix_t), pointer :: assembler_get_matrix
    assembler_get_matrix => this%matrix 
  end function assembler_get_matrix
  
  function assembler_get_array(this)
    implicit none
    class(assembler_t), target, intent(in) :: this
    class(array_t), pointer :: assembler_get_array
    assembler_get_array => this%array 
  end function assembler_get_array
  
end module assembler_names
