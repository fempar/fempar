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
module matrix_names
  use types_names
  use operator_names

  implicit none
  private

  type, abstract, extends(operator_t) :: matrix_t
  contains
     procedure (allocate_interface)        , deferred :: allocate
     procedure (init_interface)            , deferred :: init
     procedure (scal_interface)            , deferred :: scal
     procedure (add_interface)             , deferred :: add
     procedure (copy_interface)            , deferred :: copy
     procedure (free_in_stages_interface)  , deferred :: free_in_stages
     procedure (create_iterator_interface) , deferred :: create_iterator
     procedure (extract_diagonal_interface), deferred :: extract_diagonal
     procedure (get_sign_interface)        , deferred :: get_sign
     
     procedure :: is_linear               => matrix_is_linear
     procedure :: reallocate_after_remesh => matrix_reallocate_after_remesh
     ! This subroutine is an instance of the Template Method pattern with
     ! free_in_stages being the primitive method. According to this pattern,
     ! template methods cannot be overrided by subclasses
     procedure :: free                    => matrix_free_template_method
  end type matrix_t

  type, abstract :: matrix_iterator_t
     contains
       procedure (matrix_iterator_free)        , deferred :: free
       procedure (matrix_iterator_next)        , deferred :: next
       procedure (matrix_iterator_has_finished), deferred :: has_finished
       procedure (matrix_iterator_get_row)     , deferred :: get_row
       procedure (matrix_iterator_get_column)  , deferred :: get_column
       procedure (matrix_iterator_get_entry)   , deferred :: get_entry
       procedure (matrix_iterator_set_entry)   , deferred :: set_entry
    end type matrix_iterator_t
  
  abstract interface
     ! Allocates the entries of the matrix once it has been created and symbolically set-up
     subroutine allocate_interface(this) 
       import :: matrix_t
       implicit none
       class(matrix_t)       , intent(inout) :: this
     end subroutine allocate_interface
     ! Init the entries of the matrix once it has been created and symbolically set-up
     subroutine init_interface(this, alpha) 
       import :: matrix_t
       import :: rp
       implicit none
       class(matrix_t)       , intent(inout) :: this
       real(rp)              , intent(in)    :: alpha
     end subroutine init_interface
     ! Scal the entries of the matrix up
     subroutine scal_interface(this, alpha) 
       import :: matrix_t
       import :: rp
       implicit none
       class(matrix_t)       , intent(inout) :: this
       real(rp)              , intent(in)    :: alpha
     end subroutine scal_interface
     ! Performs this = alpha*op1 + beta*op2
     subroutine add_interface(this, alpha, op1, beta, op2) 
       import :: matrix_t
       import :: rp
       implicit none
       class(matrix_t)       , intent(inout) :: this
       real(rp)              , intent(in)    :: alpha
       class(matrix_t)       , intent(in)    :: op1
       real(rp)              , intent(in)    :: beta
       class(matrix_t)       , intent(in)    :: op2
     end subroutine add_interface
     !Copy matrix
     subroutine copy_interface(this, op) 
       import :: matrix_t
       import :: rp
       implicit none
       class(matrix_t)       , intent(inout) :: this
       class(matrix_t)       , intent(in)    :: op
     end subroutine copy_interface
     
     ! Progressively free a matrix_t in three stages: action={free_numeric,free_symbolic,free_clean}
     subroutine free_in_stages_interface(this,action) 
       import :: matrix_t, ip
       implicit none
       class(matrix_t)       , intent(inout) :: this
       integer(ip)           , intent(in)    :: action
     end subroutine free_in_stages_interface

     subroutine create_iterator_interface(this, iblock, jblock, iterator)
       !-----------------------------------------------------------------
       !< Get a pointer to an iterator over the matrix entries
       !-----------------------------------------------------------------
       import :: matrix_t
       import :: matrix_iterator_t
       import :: ip
       class(matrix_t)                      , intent(in)    :: this
       integer(ip)                          , intent(in)    :: iblock 
       integer(ip)                          , intent(in)    :: jblock 
       class(matrix_iterator_t), allocatable, intent(inout) :: iterator
       !-----------------------------------------------------------------
     end subroutine create_iterator_interface
     
     subroutine extract_diagonal_interface(this, diagonal)
       !-----------------------------------------------------------------
       !< Extract the diagonal entries of a matrix in a rank-1 array
       !-----------------------------------------------------------------
       import :: matrix_t
       import :: rp
       class(matrix_t), intent(in)    :: this
       real(rp)       , intent(inout) :: diagonal(:)
       !-----------------------------------------------------------------
     end subroutine extract_diagonal_interface
     
     function get_sign_interface(this) result(sign)
       !-----------------------------------------------------------------
       !< Get the sign of the sparse matrix
       !-----------------------------------------------------------------
       import :: matrix_t
       import :: ip
       class(matrix_t), intent(in) :: this
       integer(ip)                 :: sign
       !-----------------------------------------------------------------
     end function get_sign_interface
  end interface
  
  !-----------------------------------------------------------------
  !< MATRIX_ITERATOR SUBROUTINES
  !-----------------------------------------------------------------
  abstract interface
      subroutine matrix_iterator_free(this)
       !-----------------------------------------------------------------
       !< Free the iterator
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       class(matrix_iterator_t), intent(inout) :: this
     end subroutine matrix_iterator_free

     subroutine matrix_iterator_next(this)
       !-----------------------------------------------------------------
       !< Set the pointer to the following entry of the matrix
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       class(matrix_iterator_t), intent(inout) :: this
     end subroutine matrix_iterator_next

     function matrix_iterator_has_finished(this)
       !-----------------------------------------------------------------
       !< Check if the pointer of the matrix has reached the end
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       class(matrix_iterator_t), intent(in) :: this
       logical :: matrix_iterator_has_finished
     end function matrix_iterator_has_finished

     function matrix_iterator_get_row(this)
       !-----------------------------------------------------------------
       !< Get the row index of the entry of the matrix
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       import :: ip
       class(matrix_iterator_t), intent(in) :: this
       integer(ip) :: matrix_iterator_get_row
     end function matrix_iterator_get_row

     function matrix_iterator_get_column(this)
       !-----------------------------------------------------------------
       !<  Get the column index of the entry of the matrix
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       import :: ip
       class(matrix_iterator_t), intent(in) :: this
       integer(ip) :: matrix_iterator_get_column
     end function matrix_iterator_get_column

     function matrix_iterator_get_entry(this)
       !-----------------------------------------------------------------
       !< Get the value of the entry of the matrix
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       import :: rp
       class(matrix_iterator_t), intent(in) :: this
       real(rp) :: matrix_iterator_get_entry
     end function matrix_iterator_get_entry

     subroutine matrix_iterator_set_entry(this,new_value)
       !-----------------------------------------------------------------
       !< Set the value of the entry of the matrix
       !-----------------------------------------------------------------
       import :: matrix_iterator_t
       import :: rp
       class(matrix_iterator_t), intent(inout) :: this
       real(rp)                , intent(in)    :: new_value
     end subroutine matrix_iterator_set_entry
  end interface
  
  ! Data types
  public :: matrix_t
  public :: matrix_iterator_t

contains
  function matrix_is_linear ( this )
    implicit none
    class(matrix_t), intent(in) :: this
    logical :: matrix_is_linear
    matrix_is_linear = .true. 
  end function matrix_is_linear
  
  subroutine matrix_reallocate_after_remesh(this)
    implicit none
    class(matrix_t), intent(inout) :: this
  end subroutine matrix_reallocate_after_remesh

  subroutine matrix_free_template_method ( this )
    implicit none
    class(matrix_t), intent(inout) :: this
    call this%free_in_stages(free_numerical_setup)
    call this%free_in_stages(free_symbolic_setup)
    call this%free_in_stages(free_clean)
  end subroutine matrix_free_template_method
end module matrix_names
