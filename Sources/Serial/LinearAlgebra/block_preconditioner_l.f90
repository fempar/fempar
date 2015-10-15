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
module block_preconditioner_l_names
  use types_names
  use memor_names
  use operator_names
  use vector_names
  use block_vector_names

#ifdef memcheck
use iso_c_binding
#endif

  implicit none
  private
# include "debug.i90"

  ! Pointer to operator
  type p_abs_operator_t
     type(dynamic_state_operator_t), pointer :: p_op => null()
  end type p_abs_operator_t

  ! Lower block triangular preconditioner 
  type, extends(operator_t) :: block_preconditioner_l_t
!     private ! IBM XLF 14.1 bug
     integer(ip)                       :: nblocks
     type(p_abs_operator_t), allocatable :: blocks(:,:)
   contains
     procedure  :: create             => block_preconditioner_l_create
     procedure  :: set_block          => block_preconditioner_l_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_l_set_block_to_zero
     procedure  :: free               => block_preconditioner_l_free
     procedure  :: get_block          => block_preconditioner_l_get_block
     procedure  :: apply          => block_preconditioner_l_apply
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
    class(vector_t)      , intent(in)    :: x
    class(vector_t)      , intent(inout) :: y

    ! Locals
    integer(ip) :: iblk, jblk
    class(vector_t), allocatable :: aux1, aux2

    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call x%GuardTemp()
    select type(x)
    class is (block_vector_t)
       select type(y)
       class is(block_vector_t)
          allocate(aux1, mold=x%blocks(1)%vector); call aux1%default_initialization()
          allocate(aux2, mold=x%blocks(1)%vector); call aux2%default_initialization()
          do iblk=1, op%nblocks
             call aux1%clone(x%blocks(iblk)%vector)
             call aux1%copy(x%blocks(iblk)%vector)
             call aux2%clone(x%blocks(iblk)%vector)
             do jblk=1, iblk-1
                if (associated(op%blocks(iblk,jblk)%p_op)) then
                   call op%blocks(iblk,jblk)%p_op%apply(y%blocks(jblk)%vector,aux2)
                   call aux1%axpby(-1.0,aux2,1.0)
                end if
             end do
             call op%blocks(iblk,iblk)%p_op%apply(aux1,y%blocks(iblk)%vector)
             call aux1%free()
             call aux2%free()
          end do
          deallocate(aux1, aux2)
       end select
    end select
    call x%CleanTemp()
  end subroutine block_preconditioner_l_apply

  subroutine block_preconditioner_l_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)     :: nblocks

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%free()

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
    class(operator_t)                  , intent(in)    :: op 

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


  subroutine block_preconditioner_l_free (this)
    implicit none
    class(block_preconditioner_l_t), intent(inout) :: this

    ! Locals
    integer(ip) :: iblk, jblk

    do iblk=1, this%nblocks
       do jblk=1, iblk
          if (associated(this%blocks(iblk,jblk)%p_op)) then
             call this%blocks(iblk,jblk)%p_op%free()
             deallocate(this%blocks(iblk,jblk)%p_op)
          end if
       end do
    end do
    this%nblocks = 0
    if(allocated( this%blocks)) deallocate ( this%blocks )
  end subroutine block_preconditioner_l_free

  function block_preconditioner_l_get_block (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_l_t), target, intent(in) :: bop
    integer(ip)                            , intent(in) :: ib,jb
    class(operator_t)                 , pointer    :: block_preconditioner_l_get_block

    block_preconditioner_l_get_block =>  bop%blocks(ib,jb)%p_op
  end function block_preconditioner_l_get_block

end module block_preconditioner_l_names
