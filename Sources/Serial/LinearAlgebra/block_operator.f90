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
module block_operator_names
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

  ! Added type(block_operator_t) as a new implementation of class(base_operator_t).
  ! This new kind of linear operator provides the functionality to express
  ! recursively block operators where in turn their blocks can be block operators.
  ! If one aims at using the block LU recursive preconditioning machinery,
  ! then it is a MUST that the linear coefficient matrix is
  ! provided as a type(block_operator_t) instance. type(block_operator_t) provides the
  ! necessary (concrete) interface to build type(block_operator_t) instances as views, 
  ! e.g., of the components of a type(block_matrix_t) or type(par_block_matrix_t).

  ! Pointer to operator
  type p_abs_operator_t
     type(abs_operator_t), pointer :: p_op => null()
  end type p_abs_operator_t


  ! Block operator
  type, extends(base_operator_t) :: block_operator_t
!     private ! IBM XLF 14.1 bug
     integer(ip)                       :: nblocks, mblocks
     type(p_abs_operator_t), allocatable :: blocks(:,:)
   contains
     procedure  :: create             => block_operator_create
     procedure  :: set_block          => block_operator_set_block
     procedure  :: set_block_to_zero  => block_operator_set_block_to_zero
     procedure  :: destroy            => block_operator_destroy
     procedure  :: get_block          => block_operator_get_block

     procedure  :: apply          => block_operator_apply
     procedure  :: apply_fun      => block_operator_apply_fun
     procedure  :: free           => block_operator_free_tbp
  end type block_operator_t


  ! Types
  public :: block_operator_t

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_operator_apply (op,x,y)
    implicit none
    class(block_operator_t)     , intent(in)   :: op
    class(base_operand_t)      , intent(in)    :: x
    class(base_operand_t)      , intent(inout) :: y

    ! Locals
    integer(ip) :: iblk, jblk
    class(base_operand_t), allocatable :: aux

    call x%GuardTemp()

    select type(x)
    class is (block_operand_t)
       select type(y)
       class is(block_operand_t)
          allocate(aux, mold=y%blocks(1)%p_op); call aux%default_initialization()
          do iblk=1, op%mblocks
             call y%blocks(iblk)%p_op%init(0.0_rp)
             call aux%clone(y%blocks(iblk)%p_op)
             do jblk=1, op%nblocks
                if (associated(op%blocks(iblk,jblk)%p_op)) then
                    call op%blocks(iblk,jblk)%p_op%apply(x%blocks(jblk)%p_op,aux)
                    call y%blocks(iblk)%p_op%axpby(1.0,aux,1.0)
                 end if
              end do
             call aux%free()
          end do
          deallocate(aux)
       class default
          write(0,'(a)') 'block_operator_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'block_operator_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine block_operator_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_operator_apply_fun(op,x) result(y)
    implicit none
    class(block_operator_t), intent(in)  :: op
    class(base_operand_t) , intent(in)   :: x
    class(base_operand_t) , allocatable  :: y

    type(block_operand_t), allocatable :: local_y
    class(base_operand_t), allocatable :: aux
    integer(ip)                      :: iblk, jblk, first_block_in_row

    call x%GuardTemp()

    select type(x)
    class is (block_operand_t)
       allocate(local_y)
       call local_y%create(op%mblocks)
       allocate(aux, mold=x%blocks(1)%p_op); call aux%default_initialization()
       do iblk=1, op%mblocks
          first_block_in_row = 1
          do jblk=1, op%nblocks
             if (associated(op%blocks(iblk,jblk)%p_op)) then
                if ( first_block_in_row == 1 ) then
                   aux = op%blocks(iblk,jblk)%p_op * x%blocks(jblk)%p_op
                   allocate(local_y%blocks(iblk)%p_op, mold=aux); call local_y%blocks(iblk)%p_op%default_initialization()
                   local_y%blocks(iblk)%allocated = .true.
                   call local_y%blocks(iblk)%p_op%clone(aux)
                   call local_y%blocks(iblk)%p_op%copy(aux)
                   first_block_in_row = 0
                else
                   call op%blocks(iblk,jblk)%p_op%apply(x%blocks(jblk)%p_op,aux)
                   call local_y%blocks(iblk)%p_op%axpby(1.0,aux,1.0)
                end if
             end if
          end do
          call aux%free()
       end do
       deallocate(aux)
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'block_operator_t%apply_fun: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end function block_operator_apply_fun

  subroutine block_operator_free_tbp(this)
    implicit none
    class(block_operator_t), intent(inout) :: this
  end subroutine block_operator_free_tbp

  subroutine block_operator_create (bop, mblocks, nblocks)
    implicit none
    ! Parameters
    class(block_operator_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)    :: mblocks, nblocks

    
    ! Locals
    integer(ip) :: iblk, jblk

    call bop%destroy()

    bop%nblocks = nblocks
    bop%mblocks = mblocks
    allocate ( bop%blocks(mblocks,nblocks) )
    do jblk=1, nblocks
       do iblk=1, mblocks
          call bop%set_block_to_zero(iblk, jblk)
       end do
    end do
          
  end subroutine block_operator_create

  subroutine block_operator_set_block (bop, ib, jb, op)
    implicit none
    ! Parameters
    class(block_operator_t)               , intent(inout) :: bop
    integer(ip)                         , intent(in)    :: ib, jb
    class(base_operator_t)              , intent(in)    :: op 

    call op%GuardTemp()
    if ( .not. associated(bop%blocks(ib,jb)%p_op) ) then
       allocate(bop%blocks(ib,jb)%p_op)
    end if
    bop%blocks(ib,jb)%p_op = op
    call op%CleanTemp()

  end subroutine block_operator_set_block

  subroutine block_operator_set_block_to_zero (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_operator_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)    :: ib,jb

    if (associated(bop%blocks(ib,jb)%p_op)) then
       call bop%blocks(ib,jb)%p_op%free()
       deallocate(bop%blocks(ib,jb)%p_op)
    end if
    
    nullify ( bop%blocks(ib,jb)%p_op )
  end subroutine block_operator_set_block_to_zero


  subroutine block_operator_destroy (bop)
    implicit none
    class(block_operator_t), intent(inout) :: bop

    ! Locals
    integer(ip) :: iblk, jblk
    
    do jblk=1, bop%nblocks
       do iblk=1, bop%mblocks
          if (associated(bop%blocks(iblk,jblk)%p_op)) then
             call bop%blocks(iblk,jblk)%p_op%free()
             deallocate(bop%blocks(iblk,jblk)%p_op)
          end if
       end do
    end do
    bop%nblocks = 0
    bop%mblocks = 0
    if(allocated(bop%blocks)) deallocate ( bop%blocks )
  end subroutine block_operator_destroy

  function block_operator_get_block (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_operator_t), target, intent(in) :: bop
    integer(ip)                    , intent(in) :: ib,jb
    class(base_operator_t)         , pointer    :: block_operator_get_block

    block_operator_get_block =>  bop%blocks(ib,jb)%p_op
  end function block_operator_get_block

end module block_operator_names
