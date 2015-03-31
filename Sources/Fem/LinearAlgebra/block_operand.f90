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
module block_operand_names
  use types
  use memor
  use base_operand_names

  implicit none
# include "debug.i90"

  private

  ! Pointer to operator
  type p_abs_operand
     logical                      :: allocated  = .false.
     class(base_operand), pointer :: p_op => null()
  end type p_abs_operand

  ! fem_vector
  type, extends(base_operand) :: block_operand
     integer(ip)                      :: nblocks
     type(p_abs_operand), allocatable :: blocks(:)
   contains
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => block_operand_dot_tbp
     procedure :: copy => block_operand_copy_tbp
     procedure :: init => block_operand_init_tbp
     procedure :: scal => block_operand_scal_tbp
     procedure :: axpby => block_operand_axpby_tbp
     procedure :: nrm2 => block_operand_nrm2_tbp
     procedure :: clone => block_operand_clone_tbp
     procedure :: comm  => block_operand_comm_tbp
     procedure :: free  => block_operand_free_tbp
  end type block_operand

  ! Types
  public :: block_operand

contains

  subroutine block_operand_alloc (nblocks, bop)
    implicit none
    ! Parameters
    integer(ip)            , intent(in)  :: nblocks
    type(block_operand)    , intent(out) :: bop
    
    ! Locals
    integer(ip) :: iblk

    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks) )
    do iblk=1, nblocks
       bop%blocks(iblk)%allocated = .false.
       nullify(bop%blocks(iblk)%p_op)
    end do
          
  end subroutine block_operand_alloc


  subroutine block_operand_set_block (ib, op, bop)
    implicit none
    ! Parameters
    integer(ip)                         , intent(in)    :: ib
    class(base_operand), target, intent(in)    :: op 
    type(block_operand)                , intent(inout) :: bop
    
    ! A base operand to be associated to a block cannot be temporary
    assert( .not. op%IsTemp() )
    
    if ( bop%blocks(ib)%allocated ) then
       deallocate(bop%blocks(ib)%p_op)
    end if

    bop%blocks(ib)%allocated = .false.    
    bop%blocks(ib)%p_op => op
  end subroutine block_operand_set_block


  subroutine block_operand_free_tbp (this)
    implicit none
    class(block_operand), intent(inout) :: this 

    ! Locals
    integer(ip) :: iblk

    do iblk=1, this%nblocks
       if ( this%blocks(iblk)%allocated ) then
          call this%blocks(iblk)%p_op%free()
          deallocate(this%blocks(iblk)%p_op) 
       end if 
       nullify(this%blocks(iblk)%p_op)
       this%blocks(iblk)%allocated = .false. 
    end do
    this%nblocks = 0 
  end subroutine block_operand_free_tbp

 ! alpha <- op1^T * op2
 function block_operand_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(block_operand), intent(in) :: op1
   class(base_operand) , intent(in) :: op2
   real(rp) :: alpha
   ! Locals
   integer(ip) :: iblk

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (block_operand)
      assert ( op1%nblocks == op2%nblocks )
      alpha = 0.0_rp
      do iblk=1, op1%nblocks
        alpha = alpha + op1%blocks(iblk)%p_op%dot(op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function block_operand_dot_tbp

 ! op1 <- op2 
 subroutine block_operand_copy_tbp(op1,op2)
   implicit none
   class(block_operand), intent(inout) :: op1
   class(base_operand), intent(in)  :: op2
   ! Locals
   integer(ip) :: iblk
   
   call op2%GuardTemp()
   select type(op2)
   class is (block_operand)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
        call op1%blocks(iblk)%p_op%copy(op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_copy_tbp

 ! op1 <- alpha * op2
 subroutine block_operand_scal_tbp(op1,alpha,op2)
   implicit none
   class(block_operand), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand), intent(in) :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
        call op1%blocks(iblk)%p_op%scal(alpha,op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_scal_tbp
 ! op <- alpha
 subroutine block_operand_init_tbp(op,alpha)
   implicit none
   class(block_operand), intent(inout) :: op
   real(rp), intent(in) :: alpha
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      call op%blocks(iblk)%p_op%init(alpha) 
   end do
 end subroutine block_operand_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine block_operand_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(block_operand), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand), intent(in) :: op2
   real(rp), intent(in) :: beta
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand)
     do iblk=1, op1%nblocks
        call op1%blocks(iblk)%p_op%axpby(alpha, op2%blocks(iblk)%p_op, beta) 
     end do
   class default
      write(0,'(a)') 'block_operand%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_axpby_tbp

 ! alpha <- nrm2(op)
 function block_operand_nrm2_tbp(op) result(alpha)
   implicit none
   class(block_operand), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()
   alpha = op%dot(op)
   alpha = sqrt(alpha)
   call op%CleanTemp()
 end function block_operand_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine block_operand_clone_tbp(op1,op2)
   implicit none
   class(block_operand), intent(inout) :: op1
   class(base_operand), intent(in)  :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand)
     call op1%free()
     call block_operand_alloc(op2%nblocks,op1) 
     do iblk=1, op2%nblocks
       allocate(op1%blocks(iblk)%p_op, mold=op2%blocks(iblk)%p_op)
       op1%blocks(iblk)%allocated = .true. 
       call op1%blocks(iblk)%p_op%clone(op2%blocks(iblk)%p_op) 
     end do
   class default
      write(0,'(a)') 'block_operand%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_clone_tbp

 ! op <- comm(op)
 subroutine block_operand_comm_tbp(op)
   implicit none
   class(block_operand), intent(inout) :: op
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      call op%blocks(iblk)%p_op%comm() 
   end do
 end subroutine block_operand_comm_tbp


end module block_operand_names
