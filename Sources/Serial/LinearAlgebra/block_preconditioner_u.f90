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
module block_preconditioner_u_names
  use types_names
  use memor_names
  use base_operator_names
  use base_operand_names

  use block_operand_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
  private

  ! Pointer to operator
  type p_abs_operator_t
     type(abs_operator_t), pointer :: p_op => null()
  end type p_abs_operator_t


  ! Upper block triangular preconditioner 
  type, extends(base_operator_t) :: block_preconditioner_u_t
     private
     integer(ip)                       :: nblocks
     type(p_abs_operator_t), allocatable :: blocks(:,:)
   contains
     procedure  :: create             => block_preconditioner_u_create
     procedure  :: set_block          => block_preconditioner_u_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_u_set_block_to_zero
     procedure  :: destroy            => block_preconditioner_u_destroy
     procedure  :: get_block          => block_preconditioner_u_get_block

     procedure  :: apply          => block_preconditioner_u_apply
     procedure  :: apply_fun      => block_preconditioner_u_apply_fun
     procedure  :: fill_values    => block_preconditioner_u_fill_values
     procedure  :: free_values    => block_preconditioner_u_free_values
     procedure  :: free           => block_preconditioner_u_free_tbp
  end type block_preconditioner_u_t


  ! Types
  public :: block_preconditioner_u_t

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_u_apply (op,x,y)
    implicit none
    class(block_preconditioner_u_t)     , intent(in)   :: op
    class(base_operand_t)      , intent(in)    :: x
    class(base_operand_t)      , intent(inout) :: y

    ! Locals
    integer(ip) :: iblk, jblk
    class(base_operand_t), allocatable :: aux1, aux2

    call x%GuardTemp()
    select type(x)
    class is (block_operand_t)
       select type(y)
       class is(block_operand_t)
          allocate(aux1, mold=x%blocks(1)%p_op); call aux1%default_initialization()
          allocate(aux2, mold=x%blocks(1)%p_op); call aux2%default_initialization()
          do iblk=op%nblocks, 1, -1
             call aux1%clone(x%blocks(iblk)%p_op)
             call aux1%copy(x%blocks(iblk)%p_op)
             call aux2%clone(x%blocks(iblk)%p_op)
             do jblk=op%nblocks, iblk+1,-1
                if (associated(op%blocks(iblk,jblk)%p_op)) then
                   call op%blocks(iblk,jblk)%p_op%apply(y%blocks(jblk)%p_op,aux2)
                   call aux1%axpby(-1.0,aux2,1.0)
                end if
             end do
             call op%blocks(iblk,iblk)%p_op%apply(aux1,y%blocks(iblk)%p_op)
             call aux1%free()
             call aux2%free()
          end do
          deallocate(aux1, aux2)
       class default
          write(0,'(a)') 'block_preconditioner_u_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'block_preconditioner_u_t%apply: unsupported x class'
       check(1==0)
    end select
    call x%CleanTemp()
  end subroutine block_preconditioner_u_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_preconditioner_u_apply_fun(op,x) result(y) 
    implicit none
    class(block_preconditioner_u_t), intent(in)  :: op
    class(base_operand_t) , intent(in)   :: x
    class(base_operand_t) , allocatable  :: y

    type(block_operand_t), allocatable :: local_y
    class(base_operand_t), allocatable :: aux1, aux2
    integer(ip)                      :: iblk, jblk

    call x%GuardTemp()
    select type(x)
    class is (block_operand_t)
       allocate(local_y)
       call local_y%create(op%nblocks)
       allocate(aux1, mold=x%blocks(1)%p_op); call aux1%default_initialization()
       allocate(aux2, mold=x%blocks(1)%p_op); call aux2%default_initialization()
       do iblk=op%nblocks, 1, -1
          call aux1%clone(x%blocks(iblk)%p_op)
          call aux1%copy(x%blocks(iblk)%p_op)
          call aux2%clone(x%blocks(iblk)%p_op)
          do jblk=op%nblocks, iblk+1,-1
             if (associated(op%blocks(iblk,jblk)%p_op)) then
                call op%blocks(iblk,jblk)%p_op%apply(x%blocks(jblk)%p_op,aux2)
                call aux1%axpby(-1.0,aux2,1.0)
             end if
          end do
          allocate(local_y%blocks(iblk)%p_op, mold=aux1); call local_y%blocks(iblk)%p_op%default_initialization()
          local_y%blocks(iblk)%allocated = .true.
          local_y%blocks(iblk)%p_op = op%blocks(iblk,iblk)%p_op*aux1
          call aux1%free()
          call aux2%free()
       end do
       deallocate(aux1, aux2)
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'block_operand_t%apply_fun: unsupported x class'
       check(1==0)
    end select
    call x%CleanTemp()
  end function block_preconditioner_u_apply_fun

  ! op1%fill_values(op2)
  ! Fill preconditioner values
  subroutine block_preconditioner_u_fill_values (op,stage)
    implicit none
    class(block_preconditioner_u_t), intent(inout) :: op
    integer(ip), optional          , intent(in)    :: stage

    ! Locals
    integer(ip) :: iblk,stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    if(op%fill_values_stage==stage_) then
       do iblk=1, op%nblocks
          call op%blocks(iblk,iblk)%p_op%fill_values(stage_)
       end do
    end if

  end subroutine block_preconditioner_u_fill_values

  ! op1%free_values(op2)
  ! Free preconditioner values
  subroutine block_preconditioner_u_free_values (op)
    implicit none
    class(block_preconditioner_u_t), intent(inout) :: op

    ! Locals
    integer(ip) :: iblk

    if(op%do_free_values) then
       do iblk=1, op%nblocks
          call op%blocks(iblk,iblk)%p_op%free_values()
       end do
    end if

  end subroutine block_preconditioner_u_free_values

  subroutine block_preconditioner_u_free_tbp(this)
    implicit none
    class(block_preconditioner_u_t), intent(inout) :: this
  end subroutine block_preconditioner_u_free_tbp

  subroutine block_preconditioner_u_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)     :: nblocks

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%destroy()

    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks,nblocks) )
    do iblk=1, nblocks
       do jblk=iblk, bop%nblocks
          call bop%set_block_to_zero(iblk, jblk)
       end do
    end do
          
  end subroutine block_preconditioner_u_create


  subroutine block_preconditioner_u_set_block (bop, ib, jb, op)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t), intent(inout) :: bop
    integer(ip)           , intent(in)    :: ib, jb
    class(base_operator_t)    , intent(in)    :: op 

    assert ( ib <= jb )

    call op%GuardTemp()
    if ( .not. associated(bop%blocks(ib,jb)%p_op) ) then
       allocate(bop%blocks(ib,jb)%p_op)
    end if
    bop%blocks(ib,jb)%p_op = op
    call op%CleanTemp()

  end subroutine block_preconditioner_u_set_block

  subroutine block_preconditioner_u_set_block_to_zero (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t), intent(inout) :: bop
    integer(ip)           , intent(in)    :: ib,jb
    
    assert ( ib <= jb )

    if (associated(bop%blocks(ib,jb)%p_op)) then
       call bop%blocks(ib,jb)%p_op%free()
       deallocate(bop%blocks(ib,jb)%p_op)
    end if
    
    nullify ( bop%blocks(ib,jb)%p_op )
  end subroutine block_preconditioner_u_set_block_to_zero


  subroutine block_preconditioner_u_destroy (bop)
    implicit none
    class(block_preconditioner_u_t), intent(inout) :: bop

    ! Locals
    integer(ip) :: iblk, jblk

    do iblk=1, bop%nblocks
       do jblk=iblk, bop%nblocks
          if (associated(bop%blocks(iblk,jblk)%p_op)) then
             call bop%blocks(iblk,jblk)%p_op%free()
             deallocate(bop%blocks(iblk,jblk)%p_op)
          end if
       end do
    end do
    bop%nblocks = 0
    if(allocated(bop%blocks)) deallocate ( bop%blocks )
  end subroutine block_preconditioner_u_destroy

  function block_preconditioner_u_get_block (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t), target, intent(in) :: bop
    integer(ip)                            , intent(in) :: ib,jb
    class(base_operator_t)                 , pointer    :: block_preconditioner_u_get_block

    block_preconditioner_u_get_block =>  bop%blocks(ib,jb)%p_op
  end function block_preconditioner_u_get_block

end module block_preconditioner_u_names
