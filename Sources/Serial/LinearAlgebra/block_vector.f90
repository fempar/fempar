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


  integer(ip), parameter :: not_created              = 0
  integer(ip), parameter :: blocks_container_created = 1
  integer(ip), parameter :: blocks_being_built       = 2 
  integer(ip), parameter :: assembled                = 3

  ! State transition diagram for type(serial_scalar_array_t)
  ! -------------------------------------------------
  ! Input State | Action               | Output State 
  ! -------------------------------------------------
  ! not_created | free_values              | not_created
  ! not_created | free_clean               | not_created
  ! not_created | free                     | not_created
  ! not_created | create                   | blocks_container_created 

  ! blocks_container_created | free_clean  | not_created
  ! blocks_container_created | free_values | blocks_container_created
  ! blocks_container_created | set_block   | blocks_being_built 

  ! blocks_being_built | free_values   | blocks_container_created 
  ! blocks_being_built | set_block     | blocks_being_built 
  ! blocks_being_built | assemble      | assembled 

  ! assembled | allocate    | assembled 
  ! assembled | free_values | blocks_container_created 


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
     integer(ip)                   :: state = not_created
     integer(ip)                   :: nblocks
     type(p_vector_t), allocatable :: blocks(:)
   contains     
     ! TBPs related to life cycle control
     procedure :: create         => block_vector_create
     procedure :: set_block      => block_vector_set_block
     procedure :: assemble       => block_vector_assemble
     procedure :: allocate       => block_vector_allocate
     procedure :: free           => block_vector_free
     procedure :: free_in_stages => block_vector_free_in_stages

     procedure :: dot  => block_vector_dot
     procedure :: copy => block_vector_copy
     procedure :: init => block_vector_init
     procedure :: scal => block_vector_scal
     procedure :: axpby => block_vector_axpby
     procedure :: nrm2 => block_vector_nrm2
     procedure :: clone => block_vector_clone
     procedure :: same_vector_space => block_vector_same_vector_space
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
    assert(bop%state == not_created)
    bop%nblocks = nblocks
    allocate ( bop%blocks(nblocks) )
    do iblk=1, nblocks
       bop%blocks(iblk)%allocated = .false.
       nullify(bop%blocks(iblk)%vector)
    end do
    bop%state = blocks_container_created
  end subroutine block_vector_create

  subroutine block_vector_set_block (bop, ib, op)
    implicit none
    ! Parameters
    class(block_vector_t)        , intent(inout) :: bop
    integer(ip)                 , intent(in)    :: ib
    class(vector_t), target , intent(in)    :: op 
    
    ! A vector_t to be associated to a block cannot be temporary
    assert(.not. op%IsTemp())
    assert(bop%state == blocks_container_created .or. bop%state == blocks_being_built)    
    if ( bop%blocks(ib)%allocated ) then
       deallocate(bop%blocks(ib)%vector)
    end if
    bop%blocks(ib)%allocated = .false.    
    bop%blocks(ib)%vector => op
    bop%state = blocks_being_built 
  end subroutine block_vector_set_block


  subroutine block_vector_free (this)
    implicit none
    class(block_vector_t), intent(inout) :: this 

    ! Locals
    integer(ip) :: iblk
    call this%free_in_stages(free_values)
    call this%free_in_stages(free_struct)
    call this%free_in_stages(free_clean)
  end subroutine block_vector_free

  subroutine block_vector_assemble (this)
    implicit none
    class(block_vector_t), intent(inout) :: this
    integer(ip) :: iblk
    assert ( this%state == blocks_being_built)
    do iblk=1, this%nblocks
      assert(associated(this%blocks(iblk)%vector))
    end do
    this%state = assembled
  end subroutine block_vector_assemble

  subroutine block_vector_allocate (this)
    implicit none
    class(block_vector_t), intent(inout) :: this
    integer(ip) :: iblk
    assert ( this%state == assembled)
    do iblk=1, this%nblocks
       if ( this%blocks(iblk)%allocated ) then
          call this%blocks(iblk)%vector%allocate()
       end if
    end do
  end subroutine block_vector_allocate

  subroutine block_vector_free_in_stages (this,action)
    implicit none
    class(block_vector_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: action

    ! Locals
    integer(ip) :: iblk

    assert ( action==free_values.or.action==free_struct.or.action==free_clean) 
    if  ( action == free_values ) then
      if ( this%state == assembled .or. this%state == blocks_being_built ) then
        do iblk=1, this%nblocks
           if ( this%blocks(iblk)%allocated ) then
             call this%blocks(iblk)%vector%free()
             deallocate(this%blocks(iblk)%vector)
           end if
           nullify(this%blocks(iblk)%vector)
           this%blocks(iblk)%allocated = .false.
        end do
        this%state = blocks_container_created
      end if
    else if ( action == free_clean ) then
       assert ( this%state == not_created .or. this%state == blocks_container_created)
       if ( this%state == blocks_container_created ) then
         this%nblocks = 0
         deallocate( this%blocks )
         this%state = not_created
       end if
    end if
  end subroutine block_vector_free_in_stages

 ! alpha <- op1^T * op2
 function block_vector_dot(op1,op2) result(alpha)
   implicit none
   class(block_vector_t), intent(in) :: op1
   class(vector_t) , intent(in) :: op2
   real(rp) :: alpha
   ! Locals
   integer(ip) :: iblk

   assert(op1%state==assembled)
   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert(op2%state==assembled)
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
   
   assert(op1%state==assembled)
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert(op2%state==assembled)
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
   assert(op1%state==assembled)
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
      assert(op2%state==assembled)
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
   assert(op%state==assembled)
   do iblk=1, op%nblocks
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
   assert(op1%state==assembled)
   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
     assert(op2%state==assembled)
     do iblk=1, op1%nblocks
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
   assert(op%state==assembled)
   call op%GuardTemp()
   alpha = op%dot(op)
   alpha = sqrt(alpha)
   call op%CleanTemp()
 end function block_vector_nrm2

 ! op1 <- clone(op2) 
 subroutine block_vector_clone(op1,op2)
   implicit none
   class(block_vector_t), intent(inout) :: op1
   class(vector_t), intent(in)    :: op2
   ! Locals
   integer(ip) :: iblk

   call op2%GuardTemp()
   select type(op2)
   class is (block_vector_t)
     assert(op2%state==assembled)
     call op1%free()
     call op1%create(op2%nblocks) 
     do iblk=1, op2%nblocks
       allocate(op1%blocks(iblk)%vector, mold=op2%blocks(iblk)%vector)
       call op1%blocks(iblk)%vector%default_initialization()
       op1%blocks(iblk)%allocated = .true. 
       call op1%blocks(iblk)%vector%clone(op2%blocks(iblk)%vector) 
     end do
   class default
      write(0,'(a)') 'block_operand_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine block_vector_clone
 
 function block_vector_same_vector_space(this,vector)
   implicit none
   class(block_vector_t), intent(in) :: this
   class(vector_t), intent(in) :: vector
   logical :: block_vector_same_vector_space
   integer(ip) :: iblk
   assert(this%state==assembled)
   block_vector_same_vector_space = .false.
   select type(vector)
   class is (block_vector_t)
     assert(vector%state==assembled)
     block_vector_same_vector_space = (this%nblocks == vector%nblocks)
     if ( block_vector_same_vector_space ) then
       do iblk=1, this%nblocks
          assert(associated(this%blocks(iblk)%vector))
          assert(associated(vector%blocks(iblk)%vector))
          block_vector_same_vector_space = this%blocks(iblk)%vector%same_vector_space(vector%blocks(iblk)%vector)
          if ( .not. block_vector_same_vector_space ) then
            exit
          end if
       end do
     end if
   end select
 end function block_vector_same_vector_space
 
end module block_vector_names
