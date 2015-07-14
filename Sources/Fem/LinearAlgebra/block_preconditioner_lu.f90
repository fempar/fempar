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
module block_preconditioner_lu_names
  use types_names
  use memor_names
  use base_operator_names
  use base_operand_names
  use block_operand_names
  use block_preconditioner_l_names
  use block_preconditioner_u_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
  private

  ! Pointer to operator
  type p_abs_operator_t
     type(abs_operator_t), pointer :: p_op => null()
  end type p_abs_operator_t

  ! Lower block triangular preconditioner 
  type, extends(base_operator_t) :: block_preconditioner_lu_t
     private
     type(block_preconditioner_l_t) :: L
     type(block_preconditioner_u_t) :: U
  contains
     procedure  :: create             => block_preconditioner_lu_create
     procedure  :: set_block          => block_preconditioner_lu_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_lu_set_block_to_zero
     procedure  :: destroy            => block_preconditioner_lu_destroy

     procedure  :: apply          => block_preconditioner_lu_apply
     procedure  :: apply_fun      => block_preconditioner_lu_apply_fun
     procedure  :: fill_values    => block_preconditioner_lu_fill_values
     procedure  :: free_values    => block_preconditioner_lu_free_values
     procedure  :: free           => block_preconditioner_lu_free_tbp
  end type block_preconditioner_lu_t

  integer(ip), parameter :: lower = 0
  integer(ip), parameter :: upper = 1 

  ! Types
  public :: block_preconditioner_lu_t

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_lu_apply (op,x,y)
    implicit none
    class(block_preconditioner_lu_t)     , intent(in)   :: op
    class(base_operand_t)      , intent(in)    :: x
    class(base_operand_t)      , intent(inout) :: y
    class(base_operand_t), allocatable :: z
    allocate(z, mold=y); call z%default_initialization()
    call z%clone(y)
    call x%GuardTemp()
    call op%L%apply(x,z)
    call op%U%apply(z,y)
    call x%CleanTemp()
    call z%free()
    deallocate(z)
  end subroutine block_preconditioner_lu_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_preconditioner_lu_apply_fun(op,x) result(y)
    implicit none
    class(block_preconditioner_lu_t), intent(in)  :: op
    class(base_operand_t) , intent(in)   :: x
    class(base_operand_t) , allocatable  :: y
    class(base_operand_t), allocatable :: z

    allocate(z, mold=x); call z%default_initialization()
    call x%GuardTemp()
    z = op%L * x
    call op%U%apply(z,y)
    call x%CleanTemp()
    call z%free()
    deallocate(z)
  end function block_preconditioner_lu_apply_fun

  ! op1%fill_values(op2)
  ! Fill preconditioner values
  subroutine block_preconditioner_lu_fill_values (op, stage)
    implicit none
    class(block_preconditioner_lu_t), intent(inout) :: op
    integer(ip), optional           , intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    if(op%fill_values_stage==stage_) then
       call op%L%fill_values(stage_)
       call op%U%fill_values(stage_)
    end if

  end subroutine block_preconditioner_lu_fill_values

  ! op1%free_values(op2)
  ! Free preconditioner values
  subroutine block_preconditioner_lu_free_values (op)
    implicit none
    class(block_preconditioner_lu_t), intent(inout) :: op

    if(op%do_free_values) then
       call op%L%free_values()
       call op%U%free_values()
    end if

  end subroutine block_preconditioner_lu_free_values


  subroutine block_preconditioner_lu_free_tbp(this)
    implicit none
    class(block_preconditioner_lu_t), intent(inout) :: this
  end subroutine block_preconditioner_lu_free_tbp

  subroutine block_preconditioner_lu_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_lu_t)   , intent(inout) :: bop
    integer(ip)               , intent(in)    :: nblocks

    call bop%destroy()
    call bop%L%create(nblocks)
    call bop%U%create(nblocks)
  end subroutine block_preconditioner_lu_create

  subroutine block_preconditioner_lu_set_block (bop, factor, ib, jb, op)
    implicit none
    ! Parameters
    class(block_preconditioner_lu_t)    , intent(inout) :: bop
    integer(ip)                         , intent(in)    :: factor, ib, jb
    class(base_operator_t)              , intent(in)    :: op 

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


  subroutine block_preconditioner_lu_destroy (bop)
    implicit none
    class(block_preconditioner_lu_t), intent(inout) :: bop

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%L%destroy()
    call bop%U%destroy()
  end subroutine block_preconditioner_lu_destroy

end module block_preconditioner_lu_names
