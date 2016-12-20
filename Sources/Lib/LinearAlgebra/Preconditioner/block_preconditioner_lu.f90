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
module block_preconditioner_lu_names
  use types_names
  use memor_names
  use operator_names
  use vector_names
  use vector_space_names
  use block_vector_names
  use block_preconditioner_l_names
  use block_preconditioner_u_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
# include "debug.i90"

  private
  


  ! Pointer to operator
  type p_abs_operator_t
     type(lvalue_operator_t), pointer :: p_op => null()
  end type p_abs_operator_t

  ! Lower block triangular preconditioner 
  type, extends(operator_t) :: block_preconditioner_lu_t
     private
     type(block_preconditioner_l_t) :: L
     type(block_preconditioner_u_t) :: U
  contains
     procedure  :: create             => block_preconditioner_lu_create
     procedure  :: set_block          => block_preconditioner_lu_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_lu_set_block_to_zero
     procedure  :: free               => block_preconditioner_lu_free
     procedure  :: apply              => block_preconditioner_lu_apply
     procedure  :: apply_add          => block_preconditioner_lu_apply_add
     procedure  :: is_linear          => block_preconditioner_lu_is_linear
  end type block_preconditioner_lu_t

  integer(ip), parameter :: lower = 0
  integer(ip), parameter :: upper = 1 

  ! Types
  public :: block_preconditioner_lu_t

  ! Functions
  ! public :: 

contains

  function block_preconditioner_lu_is_linear(this)
    implicit none
    class(block_preconditioner_lu_t), intent(in) :: this
    logical :: block_preconditioner_lu_is_linear
    block_preconditioner_lu_is_linear = .false.
  end function block_preconditioner_lu_is_linear


  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_lu_apply (this,x,y)
    implicit none
    class(block_preconditioner_lu_t)     , intent(in)   :: this
    class(vector_t)      , intent(in)    :: x
    class(vector_t)      , intent(inout) :: y
    class(vector_t), allocatable :: z
    class(vector_space_t), pointer :: range_vector_space_L
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    range_vector_space_L => this%L%get_range_vector_space()
    call range_vector_space_L%create_vector(z)
    call x%GuardTemp()
    call this%L%apply(x,z)
    call this%U%apply(z,y)
    call x%CleanTemp()
    call z%free()
    deallocate(z)
  end subroutine block_preconditioner_lu_apply
  
  ! op%apply_add(x,y) <=> y <- op*x + y
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_lu_apply_add (this,x,y)
    implicit none
    class(block_preconditioner_lu_t)     , intent(in)   :: this
    class(vector_t)      , intent(in)    :: x
    class(vector_t)      , intent(inout) :: y
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call x%GuardTemp()
    call this%L%apply_add(x,y)
    call this%U%apply_add(x,y)
    call x%CleanTemp()
  end subroutine block_preconditioner_lu_apply_add

  subroutine block_preconditioner_lu_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_lu_t)   , intent(inout) :: bop
    integer(ip)               , intent(in)    :: nblocks

    call bop%free()
    call bop%L%create(nblocks)
    call bop%U%create(nblocks)
  end subroutine block_preconditioner_lu_create

  subroutine block_preconditioner_lu_set_block (bop, factor, ib, jb, op)
    implicit none
    ! Parameters
    class(block_preconditioner_lu_t)    , intent(inout) :: bop
    integer(ip)                         , intent(in)    :: factor, ib, jb
    class(operator_t)              , intent(in)    :: op 

    assert(factor==lower .or. factor==upper)

    call op%GuardTemp()
    if ( factor == lower ) then
      call bop%L%set_block(ib,jb,op)
    else if ( factor == upper ) then
      call bop%U%set_block(ib,jb,op)
    end if
    call op%CleanTemp()
  end subroutine block_preconditioner_lu_set_block

  subroutine block_preconditioner_lu_set_block_to_zero (bop, factor, ib, jb)
    implicit none
    ! Parameters
    class(block_preconditioner_lu_t)   , intent(inout) :: bop
    integer(ip)               , intent(in)    :: factor,ib,jb

    assert(factor==lower .or. factor==upper)
    if ( factor == lower ) then
      call bop%L%set_block_to_zero(ib,jb)
    else if ( factor == upper ) then
      call bop%U%set_block_to_zero(ib,jb)
    end if
  end subroutine block_preconditioner_lu_set_block_to_zero


  subroutine block_preconditioner_lu_free (this)
    implicit none
    class(block_preconditioner_lu_t), intent(inout) :: this

    ! Locals
    integer(ip) :: iblk, jblk

    call this%L%free()
    call this%U%free()
  end subroutine block_preconditioner_lu_free

end module block_preconditioner_lu_names
