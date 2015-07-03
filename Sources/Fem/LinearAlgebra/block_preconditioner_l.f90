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
module block_preconditioner_l_names
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

  ! Lower block triangular preconditioner 
  type, extends(base_operator_t) :: block_preconditioner_l_t
     private
     integer(ip)                       :: nblocks
     type(p_abs_operator_t), allocatable :: blocks(:,:)
   contains
     procedure  :: create             => block_preconditioner_l_create
     procedure  :: set_block          => block_preconditioner_l_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_l_set_block_to_zero
     procedure  :: destroy            => block_preconditioner_l_destroy
     procedure  :: get_block          => block_preconditioner_l_get_block

     procedure  :: apply          => block_preconditioner_l_apply
     procedure  :: apply_fun      => block_preconditioner_l_apply_fun
     procedure  :: free           => block_preconditioner_l_free_tbp
  end type block_preconditioner_l_t


  ! Types
  public :: block_preconditioner_l_t

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_l_apply (op,x,y)
    implicit none
    class(block_preconditioner_l_t)     , intent(in)   :: op
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
          allocate(aux1, aux2, mold=x%blocks(1)%p_op)
          do iblk=1, op%nblocks
             call aux1%clone(x%blocks(iblk)%p_op)
             call aux1%copy(x%blocks(iblk)%p_op)
             call aux2%clone(x%blocks(iblk)%p_op)
             do jblk=1, iblk-1
                if (associated(op%blocks(iblk,jblk)%p_op)) then
                   call op%blocks(iblk,jblk)%p_op%apply(x%blocks(jblk)%p_op,aux2)
                   call aux1%axpby(-1.0,aux2,1.0)
                end if
             end do
             call op%blocks(iblk,iblk)%p_op%apply(aux1,y%blocks(iblk)%p_op)
             call aux1%free()
             call aux2%free()
          end do
          deallocate(aux1, aux2)
       class default
          write(0,'(a)') 'block_preconditioner_l_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'block_preconditioner_l_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine block_preconditioner_l_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_preconditioner_l_apply_fun(op,x) result(y)
    implicit none
    class(block_preconditioner_l_t), intent(in)  :: op
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
       allocate(aux1, aux2, mold=x%blocks(1)%p_op)
       do iblk=1, op%nblocks
          call aux1%clone(x%blocks(iblk)%p_op)
          call aux1%copy(x%blocks(iblk)%p_op)
          call aux2%clone(x%blocks(iblk)%p_op)
          do jblk=1, iblk-1
             if (associated(op%blocks(iblk,jblk)%p_op)) then
                call op%blocks(iblk,jblk)%p_op%apply(x%blocks(jblk)%p_op,aux2)
                call aux1%axpby(-1.0,aux2,1.0)
             end if
          end do
          allocate(local_y%blocks(iblk)%p_op, mold=aux1)
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
  end function block_preconditioner_l_apply_fun

  subroutine block_preconditioner_l_free_tbp(this)
    implicit none
    class(block_preconditioner_l_t), intent(inout) :: this
  end subroutine block_preconditioner_l_free_tbp

  subroutine block_preconditioner_l_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)     :: nblocks

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%destroy()

    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks,nblocks) )
    do iblk=1, nblocks
       do jblk=1, iblk
          call bop%set_block_to_zero(iblk, jblk)
       end do
    end do
          
  end subroutine block_preconditioner_l_create


  subroutine block_preconditioner_l_set_block (bop, ib, jb, op)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t)               , intent(inout) :: bop
    integer(ip)                         , intent(in)    :: ib, jb
    class(base_operator_t)                  , intent(in)    :: op 

    assert ( ib >= jb )

    call op%GuardTemp()
    if ( .not. associated(bop%blocks(ib,jb)%p_op) ) then
       allocate(bop%blocks(ib,jb)%p_op)
    end if
    bop%blocks(ib,jb)%p_op = op
    call op%CleanTemp()

  end subroutine block_preconditioner_l_set_block

  subroutine block_preconditioner_l_set_block_to_zero (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)    :: ib,jb
    
    assert ( ib >= jb )

    if (associated(bop%blocks(ib,jb)%p_op)) then
       call bop%blocks(ib,jb)%p_op%free()
       deallocate(bop%blocks(ib,jb)%p_op)
    end if
    
    nullify ( bop%blocks(ib,jb)%p_op )
  end subroutine block_preconditioner_l_set_block_to_zero


  subroutine block_preconditioner_l_destroy (bop)
    implicit none
    class(block_preconditioner_l_t), intent(inout) :: bop

    ! Locals
    integer(ip) :: iblk, jblk

    do iblk=1, bop%nblocks
       do jblk=1, iblk
          if (associated(bop%blocks(iblk,jblk)%p_op)) then
             call bop%blocks(iblk,jblk)%p_op%free()
             deallocate(bop%blocks(iblk,jblk)%p_op)
          end if
       end do
    end do
    bop%nblocks = 0
    if(allocated( bop%blocks)) deallocate ( bop%blocks )
  end subroutine block_preconditioner_l_destroy

  function block_preconditioner_l_get_block (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t), target, intent(in) :: bop
    integer(ip)                            , intent(in) :: ib,jb
    class(base_operator_t)                 , pointer    :: block_preconditioner_l_get_block

    block_preconditioner_l_get_block =>  bop%blocks(ib,jb)%p_op
  end function block_preconditioner_l_get_block

end module block_preconditioner_l_names
