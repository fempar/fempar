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
  use types_names
  use memor_names
  use vector_names

  implicit none
# include "debug.i90"

  private

  ! Pointer to abstract_vector
  type p_abs_operand_t
     logical                      :: allocated  = .false.
     class(vector_t), pointer :: p_op => null()
  end type p_abs_operand_t

  ! Added type(block_operand_t) as a new implementation of class(vector_t).
  ! type(block_operand_t) is the only type compatible with type(block_operator_t).
  ! Therefore, if one aims to solve a linear system by means of an, e.g., block
  ! LU recursive preconditioned GMRES, then both the right hand side, and sought-after
  ! solution have to be provided to abstract_gmres as type(block_operand_t) instances.
  ! type(block_operand_t) provides the necessary (concrete) interface to build
  ! type(block_operand_t) instances as views, e.g., of the components of a
  ! type(block_vector_t) or type(par_block_vector_t).

  ! block_operand
  type, extends(vector_t) :: block_operand_t
     integer(ip)                      :: nblocks
     type(p_abs_operand_t), allocatable :: blocks(:)
   contains
     
     procedure :: create    => block_operand_create
     procedure :: set_block => block_operand_set_block
     procedure :: destroy   => block_operand_free_tbp

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
  end type block_operand_t

  ! Types
  public :: block_operand_t

contains

  subroutine block_operand_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_operand_t)    , intent(inout) :: bop
    integer(ip)             , intent(in)    :: nblocks

    ! Locals
    integer(ip) :: iblk

    call bop%free()

    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks) )
    do iblk=1, nblocks
       bop%blocks(iblk)%allocated = .false.
       nullify(bop%blocks(iblk)%p_op)
    end do
          
  end subroutine block_operand_create

  subroutine block_operand_set_block (bop, ib, op)
    implicit none
    ! Parameters
    class(block_operand_t)        , intent(inout) :: bop
    integer(ip)                 , intent(in)    :: ib
    class(vector_t), target , intent(in)    :: op 
    
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
    class(block_operand_t), intent(inout) :: this 

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
    if(allocated( this%blocks )) deallocate( this%blocks )
  end subroutine block_operand_free_tbp

 ! alpha <- op1^T * op2
 function block_operand_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(block_operand_t), intent(in) :: op1
   class(vector_t) , intent(in) :: op2
   real(rp) :: alpha
   ! Locals
   integer(ip) :: iblk

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (block_operand_t)
      assert ( op1%nblocks == op2%nblocks )
      alpha = 0.0_rp
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%p_op))
         assert(associated(op2%blocks(iblk)%p_op))   
         alpha = alpha + op1%blocks(iblk)%p_op%dot(op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function block_operand_dot_tbp

 ! op1 <- op2 
 subroutine block_operand_copy_tbp(op1,op2)
   implicit none
   class(block_operand_t), intent(inout) :: op1
   class(vector_t), intent(in)  :: op2
   ! Locals
   integer(ip) :: iblk
   
   call op2%GuardTemp()
   select type(op2)
   class is (block_operand_t)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%p_op))
         assert(associated(op2%blocks(iblk)%p_op)) 
        call op1%blocks(iblk)%p_op%copy(op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_copy_tbp

 ! op1 <- alpha * op2
 subroutine block_operand_scal_tbp(op1,alpha,op2)
   implicit none
   class(block_operand_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand_t)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%p_op))
         assert(associated(op2%blocks(iblk)%p_op)) 
        call op1%blocks(iblk)%p_op%scal(alpha,op2%blocks(iblk)%p_op) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_scal_tbp
 ! op <- alpha
 subroutine block_operand_init_tbp(op,alpha)
   implicit none
   class(block_operand_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      assert(associated(op%blocks(iblk)%p_op))
      call op%blocks(iblk)%p_op%init(alpha) 
   end do
 end subroutine block_operand_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine block_operand_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(block_operand_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2
   real(rp), intent(in) :: beta
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand_t)
     do iblk=1, op1%nblocks
        assert(associated(op1%blocks(iblk)%p_op))
        assert(associated(op2%blocks(iblk)%p_op))
        call op1%blocks(iblk)%p_op%axpby(alpha, op2%blocks(iblk)%p_op, beta) 
     end do
   class default
      write(0,'(a)') 'block_operand_t%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_axpby_tbp

 ! alpha <- nrm2(op)
 function block_operand_nrm2_tbp(op) result(alpha)
   implicit none
   class(block_operand_t), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()
   alpha = op%dot(op)
   alpha = sqrt(alpha)
   call op%CleanTemp()
 end function block_operand_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine block_operand_clone_tbp(op1,op2)
   implicit none
   class(block_operand_t)         , intent(inout) :: op1
   class(vector_t) , target , intent(in)    :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_operand_t)
     call op1%free()
     call op1%create(op2%nblocks) 
     do iblk=1, op2%nblocks
        assert(associated(op2%blocks(iblk)%p_op))
       allocate(op1%blocks(iblk)%p_op, mold=op2%blocks(iblk)%p_op); call op1%blocks(iblk)%p_op%default_initialization()
       op1%blocks(iblk)%allocated = .true. 
       call op1%blocks(iblk)%p_op%clone(op2%blocks(iblk)%p_op) 
     end do
   class default
      write(0,'(a)') 'block_operand_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_operand_clone_tbp

 ! op <- comm(op)
 subroutine block_operand_comm_tbp(op)
   implicit none
   class(block_operand_t), intent(inout) :: op
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      assert(associated(op%blocks(iblk)%p_op))
      call op%blocks(iblk)%p_op%comm() 
   end do
 end subroutine block_operand_comm_tbp


end module block_operand_names
