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
module block_preconditioner_u_names
  use types_names
  use memor_names
  use operator_names
  use vector_names

  use block_vector_names

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


  ! Upper block triangular preconditioner 
  type, extends(operator_t) :: block_preconditioner_u_t
     private
     integer(ip)                       :: nblocks
     type(p_abs_operator_t), allocatable :: blocks(:,:)
   contains
     procedure  :: create             => block_preconditioner_u_create
     procedure  :: set_block          => block_preconditioner_u_set_block
     procedure  :: set_block_to_zero  => block_preconditioner_u_set_block_to_zero
     procedure  :: free               => block_preconditioner_u_free
     procedure  :: get_block          => block_preconditioner_u_get_block
     procedure  :: apply              => block_preconditioner_u_apply
     procedure  :: apply_add          => block_preconditioner_u_apply_add
     procedure  :: is_linear          => block_preconditioner_u_is_linear
     procedure  :: update_matrix      => block_preconditioner_u_update_matrix
   end type block_preconditioner_u_t


  ! Types
  public :: block_preconditioner_u_t

  ! Functions
  ! public :: 

contains

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_u_apply (this,x,y)
    implicit none
    class(block_preconditioner_u_t)     , intent(inout)   :: this
    class(vector_t)      , intent(in)    :: x
    class(vector_t)      , intent(inout) :: y

    ! Locals
    integer(ip) :: iblk, jblk
    class(vector_t), allocatable :: aux

    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
    class is (block_vector_t)
       select type(y)
       class is(block_vector_t)          
          allocate(aux, mold=x%blocks(1)%vector); call aux%default_initialization()
          do iblk=this%nblocks, 1, -1
             call aux%clone(x%blocks(iblk)%vector)
             call aux%scal(-1.0_rp,x%blocks(iblk)%vector)
             do jblk=this%nblocks, iblk+1,-1
                if (associated(this%blocks(iblk,jblk)%p_op)) then
                   call this%blocks(iblk,jblk)%p_op%apply_add(y%blocks(jblk)%vector,aux)
                end if
             end do
             call aux%scal(-1.0_rp,aux)
             call this%blocks(iblk,iblk)%p_op%apply(aux,y%blocks(iblk)%vector)
          end do
          deallocate(aux)
       end select
    end select
    call x%CleanTemp()
  end subroutine block_preconditioner_u_apply
  
  ! op%apply_add(x,y) <=> y <- op*x + y
  ! Implicitly assumes that y is already allocated
  subroutine block_preconditioner_u_apply_add (this,x,y)
    implicit none
    class(block_preconditioner_u_t)     , intent(inout)   :: this
    class(vector_t)      , intent(in)    :: x
    class(vector_t)      , intent(inout) :: y

    ! Locals
    integer(ip) :: iblk, jblk
    class(vector_t), allocatable :: aux

    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call x%GuardTemp()
    select type(x)
    class is (block_vector_t)
       select type(y)
       class is(block_vector_t)
          allocate(aux, mold=x%blocks(1)%vector); call aux%default_initialization()
          do iblk=this%nblocks, 1, -1
             call aux%clone(x%blocks(iblk)%vector)
             call aux%scal(-1.0_rp,x%blocks(iblk)%vector)
             do jblk=this%nblocks, iblk+1,-1
                if (associated(this%blocks(iblk,jblk)%p_op)) then
                   call this%blocks(iblk,jblk)%p_op%apply_add(y%blocks(jblk)%vector,aux)
                end if
             end do
             call aux%scal(-1.0_rp,aux)
             call this%blocks(iblk,iblk)%p_op%apply_add(aux,y%blocks(iblk)%vector)
          end do
          deallocate(aux)
       end select
    end select
    call x%CleanTemp()
  end subroutine block_preconditioner_u_apply_add
  
  function block_preconditioner_u_is_linear(this)
    implicit none
    class(block_preconditioner_u_t), intent(in) :: this
    logical :: block_preconditioner_u_is_linear
    block_preconditioner_u_is_linear = .false.
  end function block_preconditioner_u_is_linear
  
  subroutine block_preconditioner_u_update_matrix(this, same_nonzero_pattern)
    implicit none
    class(block_preconditioner_u_t), intent(inout)    :: this
    logical                        , intent(in)       :: same_nonzero_pattern
    integer(ip) :: iblk, jblk
    do iblk=this%nblocks, 1, -1         
       do jblk=this%nblocks, iblk+1,-1
         if (associated(this%blocks(iblk,jblk)%p_op)) then
            call this%blocks(iblk,jblk)%p_op%update_matrix(same_nonzero_pattern)
         end if
       end do       
    end do
  end subroutine block_preconditioner_u_update_matrix

  subroutine block_preconditioner_u_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t)   , intent(inout) :: bop
    integer(ip)             , intent(in)     :: nblocks

    ! Locals
    integer(ip) :: iblk, jblk

    call bop%free()

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
    class(operator_t)    , intent(in)    :: op 

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


  subroutine block_preconditioner_u_free (this)
    implicit none
    class(block_preconditioner_u_t), intent(inout) :: this

    ! Locals
    integer(ip) :: iblk, jblk

    do iblk=1, this%nblocks
       do jblk=iblk, this%nblocks
          if (associated(this%blocks(iblk,jblk)%p_op)) then
             call this%blocks(iblk,jblk)%p_op%free()
             deallocate(this%blocks(iblk,jblk)%p_op)
          end if
       end do
    end do
    this%nblocks = 0
    if(allocated(this%blocks)) deallocate ( this%blocks )
  end subroutine block_preconditioner_u_free

  function block_preconditioner_u_get_block (bop,ib,jb)
    implicit none
    ! Parameters
    class(block_preconditioner_u_t), target, intent(in) :: bop
    integer(ip)                            , intent(in) :: ib,jb
    class(operator_t)                 , pointer    :: block_preconditioner_u_get_block

    block_preconditioner_u_get_block =>  bop%blocks(ib,jb)%p_op
  end function block_preconditioner_u_get_block

end module block_preconditioner_u_names
