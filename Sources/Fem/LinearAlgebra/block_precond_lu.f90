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
# include "debug.i90"
module block_precond_lu_names
use types_names
use memor_names
  use base_operator_names
  use base_operand_names
  use block_operand_names
  use block_precond_l_names
  use block_precond_u_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
  private

  ! Pointer to operator
  type p_abs_operator
     type(abs_operator), pointer :: p_op => null()
  end type p_abs_operator

  ! Lower block triangular preconditioner 
  type, extends(base_operator) :: block_precond_lu
     private
     type(block_precond_l) :: L
     type(block_precond_u) :: U
  contains
     procedure  :: create             => block_precond_lu_create
     procedure  :: set_block          => block_precond_lu_set_block
     procedure  :: set_block_to_zero  => block_precond_lu_set_block_to_zero
     procedure  :: destroy            => block_precond_lu_destroy

     procedure  :: apply          => block_precond_lu_apply
     procedure  :: apply_fun      => block_precond_lu_apply_fun
     procedure  :: free           => block_precond_lu_free_tbp
  end type block_precond_lu

  integer(ip), parameter :: lower = 0
  integer(ip), parameter :: upper = 1 

  ! Types
  public :: block_precond_lu

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_precond_lu_apply (op,x,y)
    implicit none
    class(block_precond_lu)     , intent(in)   :: op
    class(base_operand)      , intent(in)    :: x
    class(base_operand)      , intent(inout) :: y

    class(base_operand), allocatable :: z
    allocate(z, mold=y)
    call z%clone(y)
    call x%GuardTemp()
    call op%L%apply(x,z)
    call op%U%apply(z,y)
    call x%CleanTemp()
    call z%free()
    deallocate(z)
  end subroutine block_precond_lu_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_precond_lu_apply_fun(op,x) result(y)
    implicit none
    class(block_precond_lu), intent(in)  :: op
    class(base_operand) , intent(in)   :: x
    class(base_operand) , allocatable  :: y

    class(base_operand), allocatable :: z

    allocate(z, mold=x)
    call x%GuardTemp()
    z = op%L * x
    call op%U%apply(z,y)
    call x%CleanTemp()
    call z%free()
    deallocate(z)
  end function block_precond_lu_apply_fun

  subroutine block_precond_lu_free_tbp(this)
    implicit none
    class(block_precond_lu), intent(inout) :: this
  end subroutine block_precond_lu_free_tbp

  subroutine block_precond_lu_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_precond_lu)   , intent(inout) :: bop
    integer(ip)               , intent(in)    :: nblocks

    call bop%destroy()
    call bop%L%create(nblocks)
    call bop%U%create(nblocks)
  end subroutine block_precond_lu_create

  subroutine block_precond_lu_set_block (bop, factor, ib, jb, op)
    implicit none
    ! Parameters
    class(block_precond_lu)             , intent(inout) :: bop
    integer(ip)                         , intent(in)    :: factor, ib, jb
    type(abs_operator)                  , intent(in)    :: op 

    assert(factor==lower .or. factor==upper)

    call op%GuardTemp()
    if ( factor == lower ) then
      call bop%L%set_block(ib,jb,op)
    else if ( factor == upper ) then
      call bop%U%set_block(ib,jb,op)
    end if
    call op%CleanTemp()
  end subroutine block_precond_lu_set_block

  subroutine block_precond_lu_set_block_to_zero (bop, factor, ib, jb)
    implicit none
    ! Parameters
    class(block_precond_lu)   , intent(inout) :: bop
    integer(ip)               , intent(in)    :: factor,ib,jb

    assert(factor==lower .or. factor==upper)
    if ( factor == lower ) then
      call bop%L%set_block_to_zero(ib,jb)
    else if ( factor == upper ) then
      call bop%U%set_block_to_zero(ib,jb)
    end if
  end subroutine block_precond_lu_set_block_to_zero


  subroutine block_precond_lu_destroy (bop)
    implicit none
    class(block_precond_lu), intent(inout) :: bop

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%L%destroy()
    call bop%U%destroy()
  end subroutine block_precond_lu_destroy

end module block_precond_lu_names
