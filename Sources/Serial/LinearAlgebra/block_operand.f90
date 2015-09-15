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
module block_vector_names
  use types_names
  use memor_names
  use vector_names

  implicit none
# include "debug.i90"

  private

  type p_vector_t
     logical :: allocated  = .false.
     class(vector_t), pointer :: vector => null()
  end type p_vector_t

  ! Added type(block_vector_t) as a new implementation of class(vector_t).
  ! type(block_vector_t) is the only type compatible with type(block_vector_t).
  ! Therefore, if one aims to solve a linear system by means of an, e.g., block
  ! LU recursive preconditioned GMRES, then both the right hand side, and sought-after
  ! solution have to be provided to abstract_gmres as type(block_vector_t) instances.
  ! type(block_vector_t) provides the necessary (concrete) interface to build
  ! type(block_vector_t) instances as views, e.g., of the components of a
  ! type(block_array_t) or type(par_block_array_t).

  ! block_operand
  type, extends(vector_t) :: block_vector_t
     integer(ip)                   :: nblocks
     type(p_vector_t), allocatable :: blocks(:)
   contains     
     procedure :: create    => block_vector_create
     procedure :: set_block => block_vector_set_block
     procedure :: dot  => block_vector_dot
     procedure :: copy => block_vector_copy
     procedure :: init => block_vector_init
     procedure :: scal => block_vector_scal
     procedure :: axpby => block_vector_axpby
     procedure :: nrm2 => block_vector_nrm2
     procedure :: clone => block_vector_clone
     procedure :: comm  => block_vector_comm
     procedure :: free  => block_vector_free
  end type block_vector_t

  ! Types
  public :: block_vector_t

contains

  subroutine block_vector_create (bop, nblocks)
    implicit none
    ! Parameters
    class(block_vector_t)    , intent(inout) :: bop
    integer(ip)             , intent(in)    :: nblocks

    ! Locals
    integer(ip) :: iblk

    call bop%free()

    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks) )
    do iblk=1, nblocks
       bop%blocks(iblk)%allocated = .false.
       nullify(bop%blocks(iblk)%vector)
    end do
          
  end subroutine block_vector_create

  subroutine block_vector_set_block (bop, ib, op)
    implicit none
    ! Parameters
    class(block_vector_t)        , intent(inout) :: bop
    integer(ip)                 , intent(in)    :: ib
    class(vector_t), target , intent(in)    :: op 
    
    ! A vector_t to be associated to a block cannot be temporary
    assert( .not. op%IsTemp() )
    
    if ( bop%blocks(ib)%allocated ) then
       deallocate(bop%blocks(ib)%vector)
    end if

    bop%blocks(ib)%allocated = .false.    
    bop%blocks(ib)%vector => op
  end subroutine block_vector_set_block


  subroutine block_vector_free (this)
    implicit none
    class(block_vector_t), intent(inout) :: this 

    ! Locals
    integer(ip) :: iblk

    do iblk=1, this%nblocks
       if ( this%blocks(iblk)%allocated ) then
          call this%blocks(iblk)%vector%free()
          deallocate(this%blocks(iblk)%vector) 
       end if 
       nullify(this%blocks(iblk)%vector)
       this%blocks(iblk)%allocated = .false. 
    end do
    this%nblocks = 0 
    if(allocated( this%blocks )) deallocate( this%blocks )
  end subroutine block_vector_free

 ! alpha <- op1^T * op2
 function block_vector_dot(op1,op2) result(alpha)
   implicit none
   class(block_vector_t), intent(in) :: op1
   class(vector_t) , intent(in) :: op2
   real(rp) :: alpha
   ! Locals
   integer(ip) :: iblk

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert ( op1%nblocks == op2%nblocks )
      alpha = 0.0_rp
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%vector))
         assert(associated(op2%blocks(iblk)%vector))   
         alpha = alpha + op1%blocks(iblk)%vector%dot(op2%blocks(iblk)%vector) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function block_vector_dot

 ! op1 <- op2 
 subroutine block_vector_copy(op1,op2)
   implicit none
   class(block_vector_t), intent(inout) :: op1
   class(vector_t), intent(in)  :: op2
   ! Locals
   integer(ip) :: iblk
   
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%vector))
         assert(associated(op2%blocks(iblk)%vector)) 
        call op1%blocks(iblk)%vector%copy(op2%blocks(iblk)%vector) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_vector_copy

 ! op1 <- alpha * op2
 subroutine block_vector_scal(op1,alpha,op2)
   implicit none
   class(block_vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert ( op1%nblocks == op2%nblocks )
      do iblk=1, op1%nblocks
         assert(associated(op1%blocks(iblk)%vector))
         assert(associated(op2%blocks(iblk)%vector)) 
        call op1%blocks(iblk)%vector%scal(alpha,op2%blocks(iblk)%vector) 
      end do
   class default
      write(0,'(a)') 'block_operand_t%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_vector_scal
 ! op <- alpha
 subroutine block_vector_init(op,alpha)
   implicit none
   class(block_vector_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      assert(associated(op%blocks(iblk)%vector))
      call op%blocks(iblk)%vector%init(alpha) 
   end do
 end subroutine block_vector_init

 ! op1 <- alpha*op2 + beta*op1
 subroutine block_vector_axpby(op1, alpha, op2, beta)
   implicit none
   class(block_vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2
   real(rp), intent(in) :: beta
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
     do iblk=1, op1%nblocks
        assert(associated(op1%blocks(iblk)%vector))
        assert(associated(op2%blocks(iblk)%vector))
        call op1%blocks(iblk)%vector%axpby(alpha, op2%blocks(iblk)%vector, beta) 
     end do
   class default
      write(0,'(a)') 'block_operand_t%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_vector_axpby

 ! alpha <- nrm2(op)
 function block_vector_nrm2(op) result(alpha)
   implicit none
   class(block_vector_t), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()
   alpha = op%dot(op)
   alpha = sqrt(alpha)
   call op%CleanTemp()
 end function block_vector_nrm2

 ! op1 <- clone(op2) 
 subroutine block_vector_clone(op1,op2)
   implicit none
   class(block_vector_t)         , intent(inout) :: op1
   class(vector_t) , target , intent(in)    :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
     call op1%free()
     call op1%create(op2%nblocks) 
     do iblk=1, op2%nblocks
        assert(associated(op2%blocks(iblk)%vector))
       allocate(op1%blocks(iblk)%vector, mold=op2%blocks(iblk)%vector); call op1%blocks(iblk)%vector%default_initialization()
       op1%blocks(iblk)%allocated = .true. 
       call op1%blocks(iblk)%vector%clone(op2%blocks(iblk)%vector) 
     end do
   class default
      write(0,'(a)') 'block_operand_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_vector_clone

 ! op <- comm(op)
 subroutine block_vector_comm(op)
   implicit none
   class(block_vector_t), intent(inout) :: op
   ! Locals
   integer(ip) :: iblk
   do iblk=1, op%nblocks
      assert(associated(op%blocks(iblk)%vector))
      call op%blocks(iblk)%vector%comm() 
   end do
 end subroutine block_vector_comm
 
end module block_vector_names
