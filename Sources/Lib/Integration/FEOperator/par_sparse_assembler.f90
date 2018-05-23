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
module par_sparse_assembler_names
  use types_names
  use allocatable_array_names

  ! Abstract modules
  use assembler_names
  use matrix_names
  use array_names

  ! Concrete implementations
  use par_sparse_matrix_names
  use par_scalar_array_names

  implicit none
# include "debug.i90"
  private

  type, extends(assembler_t) :: par_sparse_assembler_t
    private
    type(par_sparse_matrix_t), pointer :: par_sparse_matrix => NULL()
    type(par_scalar_array_t) , pointer :: par_scalar_array  => NULL()
  contains
    procedure :: set_matrix              => par_sparse_assembler_set_matrix
    procedure :: set_array               => par_sparse_assembler_set_array
    procedure :: assembly_array          => par_sparse_assembler_assembly_array
    procedure :: assembly_matrix         => par_sparse_assembler_assembly_matrix
    procedure :: allocate                => par_sparse_assembler_allocate
    procedure :: compress_storage_matrix => par_sparse_assembler_compress_storage_matrix 
    procedure :: compress_storage_array  => par_sparse_assembler_compress_storage_array 
  end type

! Data types
public :: par_sparse_assembler_t

contains
  subroutine par_sparse_assembler_set_matrix ( this, matrix ) 
     implicit none
     class(par_sparse_assembler_t) , intent(inout) :: this
     class(matrix_t),       pointer, intent(in)    :: matrix
     call assembler_set_matrix(this,matrix)
     select type(matrix)
       class is(par_sparse_matrix_t)
         this%par_sparse_matrix => matrix
       class default
       check(.false.)
     end select
  end subroutine par_sparse_assembler_set_matrix 
  
  subroutine par_sparse_assembler_set_array ( this, array ) 
     implicit none
     class(par_sparse_assembler_t) , intent(inout) :: this
     class(array_t)  ,      pointer, intent(in)    :: array
     call assembler_set_array(this,array)
     select type(array)
       class is(par_scalar_array_t)
         this%par_scalar_array => array
       class default
       check(.false.)
     end select
  end subroutine par_sparse_assembler_set_array

  subroutine par_sparse_assembler_assembly_array( this,           & 
                                                     num_fields,  &
                                                     field_blocks,   &
                                                     field_coupling, &
                                                     num_dofs,    &
                                                     fe_dofs,       &
                                                     elvec )
    implicit none
    class(par_sparse_assembler_t), intent(inout) :: this
    integer(ip)                               , intent(in)    :: num_fields
    integer(ip)                               , intent(in)    :: field_blocks(num_fields)
    logical                                   , intent(in)    :: field_coupling(num_fields,num_fields)
    integer(ip)                               , intent(in)    :: num_dofs(num_fields)
    type(i1p_t)                               , intent(in)    :: fe_dofs(num_fields)
    real(rp)                                  , intent(in)    :: elvec(:)
    integer(ip) :: inode, idof, ielvec, ife_space
    assert ( associated(this%par_scalar_array) )
    ielvec = 0
    do ife_space = 1, num_fields
      call this%par_scalar_array%add( num_dofs(ife_space), &
                                      fe_dofs(ife_space)%p,  &
                                      ielvec,                 &
                                      elvec )
      ielvec = ielvec + num_dofs(ife_space)
    end do
  end subroutine par_sparse_assembler_assembly_array
      
  subroutine par_sparse_assembler_assembly_matrix( this,           &
                                                   num_fields,     &
                                                   field_blocks,   &
                                                   field_coupling, &
                                                   num_row_dofs,   &
                                                   num_col_dofs,   &
                                                   fe_dofs_row,    &
                                                   fe_dofs_col,    &
                                                   elmat )
    implicit none
    class(par_sparse_assembler_t), intent(inout) :: this
    integer(ip)                               , intent(in)    :: num_fields
    integer(ip)                               , intent(in)    :: field_blocks(num_fields)
    logical                                   , intent(in)    :: field_coupling(num_fields,num_fields)
    integer(ip)                               , intent(in)    :: num_row_dofs(num_fields)
    integer(ip)                               , intent(in)    :: num_col_dofs(num_fields)
    type(i1p_t)                               , intent(in)    :: fe_dofs_row(num_fields)
    type(i1p_t)                               , intent(in)    :: fe_dofs_col(num_fields)
    real(rp)                                  , intent(in)    :: elmat(:,:) 
    integer(ip) :: ife_space, jfe_space
    integer(ip) :: idof, jdof 
    integer(ip) :: inode, jnode
    integer(ip) :: ielmat, jelmat
    assert ( associated(this%par_sparse_matrix) )
    ielmat=0
    do ife_space=1, num_fields
       jelmat=0
       do jfe_space=1, num_fields
          if ((field_coupling(ife_space,jfe_space))) then
             call this%par_sparse_matrix%insert(num_row_dofs(ife_space),num_col_dofs(jfe_space),               &
                  &                             fe_dofs_row(ife_space)%p,fe_dofs_col(jfe_space)%p,ielmat,jelmat, &
                  &                             elmat)
          end if
          jelmat=jelmat+num_col_dofs(jfe_space)
       end do
       ielmat=ielmat+num_row_dofs(ife_space)
    end do    
  end subroutine par_sparse_assembler_assembly_matrix
    
  subroutine par_sparse_assembler_allocate( this )
    implicit none
    class(par_sparse_assembler_t), intent(inout) :: this
    assert ( associated(this%par_scalar_array) )
    call this%par_scalar_array%allocate()
  end subroutine par_sparse_assembler_allocate

  subroutine par_sparse_assembler_compress_storage_matrix( this, & 
                                                                 sparse_matrix_storage_format )
    implicit none
    class(par_sparse_assembler_t) , intent(inout) :: this
    character(*)                              , intent(in)    ::  sparse_matrix_storage_format
    assert ( associated(this%par_sparse_matrix) )
    call this%par_sparse_matrix%convert(sparse_matrix_storage_format)
  end subroutine par_sparse_assembler_compress_storage_matrix
  
  subroutine par_sparse_assembler_compress_storage_array( this )
    implicit none
    class(par_sparse_assembler_t) , intent(inout) :: this
    assert ( associated(this%par_scalar_array) )
    call this%par_scalar_array%comm() 
  end subroutine par_sparse_assembler_compress_storage_array

end module par_sparse_assembler_names

