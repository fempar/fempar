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
module serial_block_array_names
  use types_names
  use memor_names
  use serial_scalar_array_names
  use vector_names
  use graph_names
  implicit none
# include "debug.i90"

  private

  ! vector
  type, extends(vector_t) :: serial_block_array_t
     integer(ip)                              :: nblocks = 0
     type(serial_scalar_array_t), allocatable :: blocks(:)
   contains
     procedure, private :: serial_block_array_create_only_blocks_container
	 procedure, private :: serial_block_array_create_blocks_container_and_all_blocks
     generic :: create => serial_block_array_create_only_blocks_container, &
	                      serial_block_array_create_blocks_container_and_all_blocks
						  
     procedure :: create_view => serial_block_array_create_view
     procedure :: print => 	serial_block_array_print			

     procedure :: dot  => serial_block_array_dot
     procedure :: copy => serial_block_array_copy
     procedure :: init => serial_block_array_init
     procedure :: scal => serial_block_array_scal
     procedure :: axpby => serial_block_array_axpby
     procedure :: nrm2 => serial_block_array_nrm2
     procedure :: clone => serial_block_array_clone
     procedure :: comm  => serial_block_array_comm
     procedure :: free  => serial_block_array_free
  end type serial_block_array_t

  ! Types
  public :: serial_block_array_t

contains

!=============================================================================
  subroutine serial_block_array_create_only_blocks_container(this,nblocks)
    implicit none
    class(serial_block_array_t), intent(out) :: this
    integer(ip)                , intent(in)  :: nblocks
	
	integer(ip) :: istat

    this%nblocks = nblocks
    allocate ( this%blocks(nblocks), stat=istat )
	check ( istat == 0 )
  end subroutine serial_block_array_create_only_blocks_container
  
  !=============================================================================
  subroutine serial_block_array_create_blocks_container_and_all_blocks(this,nblocks,size_of_blocks)
    implicit none
    class(serial_block_array_t), intent(out) :: this
    integer(ip)                , intent(in)  :: nblocks
	integer(ip)                , intent(in)  :: size_of_blocks(nblocks)
	
	integer(ip) :: istat
    integer(ip) :: ib

    call this%create(nblocks)
	
    do ib=1, this%nblocks
      call this%blocks(ib)%create(size_of_blocks(ib))
    end do 
	
  end subroutine serial_block_array_create_blocks_container_and_all_blocks

  !=============================================================================
  subroutine serial_block_array_create_view (this, start, end, tvec)
    implicit none
    class(serial_block_array_t), intent(in) :: this
    integer(ip)                , intent(in) :: start
    integer(ip)                , intent(in) :: end
    type(serial_block_array_t), intent(out) :: tvec
    integer(ip) :: istat, ib

    call tvec%create(this%nblocks)

    do ib=1, this%nblocks
       call this%blocks(ib)%create_view (start, end, tvec%blocks(ib))
    end do
  end subroutine serial_block_array_create_view

  subroutine serial_block_array_print (this,luout)
    implicit none
    class(serial_block_array_t), intent(in) :: this
    integer(ip)                , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, this%nblocks
      call this%blocks(ib)%print(luout)
    end do 
  end subroutine serial_block_array_print
  
  ! alpha <- op1^T * op2
  function serial_block_array_dot(op1,op2) result(alpha)
    implicit none
    ! Parameters
    class(serial_block_array_t), intent(in)  :: op1
    class(vector_t)   , intent(in)  :: op2
    real(rp) :: alpha

    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
    class is (serial_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       alpha = 0.0_rp
       do ib=1,op1%nblocks
          aux = op1%blocks(ib)%dot(op2%blocks(ib))
          alpha = alpha + aux
       end do
    class default
       write(0,'(a)') 'serial_block_array_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function serial_block_array_dot

  ! op1 <- op2 
  subroutine serial_block_array_copy(op1,op2)
    implicit none
    ! Parameters
    class(serial_block_array_t), intent(inout) :: op1
    class(vector_t)   , intent(in)    :: op2

    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (serial_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%copy(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'serial_block_array_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_copy

  ! op <- alpha
  subroutine serial_block_array_init(op,alpha)
    implicit none
    class(serial_block_array_t), intent(inout) :: op 
    real(rp)             , intent(in)    :: alpha  
    ! Locals
    integer(ip) :: ib

    do ib=1, op%nblocks
       call op%blocks(ib)%init(alpha)
    end do    
  end subroutine serial_block_array_init
  
  ! op1 <- alpha * op2
  subroutine serial_block_array_scal(op1,alpha,op2)
    implicit none
    ! Parameters 
    class(serial_block_array_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(vector_t), intent(in)    :: op2
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (serial_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%scal(alpha,op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'serial_block_array_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_scal

  ! op1 <- alpha*op2 + beta*op1
  subroutine serial_block_array_axpby(op1,alpha,op2,beta)
    implicit none
    class(serial_block_array_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(vector_t), intent(in)    :: op2
    real(rp)             , intent(in)    :: beta
    ! Locals
    integer(ip)                           :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (serial_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%axpby(alpha,op2%blocks(ib),beta)
       end do
    class default
       write(0,'(a)') 'serial_block_array_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_axpby

  ! alpha <- nrm2(op)
  function serial_block_array_nrm2(op) result(alpha)
    implicit none
    class(serial_block_array_t), intent(in) :: op
    real(rp) :: alpha

    call op%GuardTemp()
    alpha = op%dot(op)
    alpha = sqrt(alpha)
   call op%CleanTemp()
  end function serial_block_array_nrm2

  ! op1 <- clone(op2) 
  subroutine serial_block_array_clone(op1,op2)
    implicit none
    ! Parameters
    class(serial_block_array_t)     , intent(inout) :: op1
    class(vector_t), target, intent(in)    :: op2
 
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (serial_block_array_t)
       op1%nblocks = op2%nblocks
       if(allocated(op1%blocks)) deallocate ( op1%blocks )
       allocate(op1%blocks(op1%nblocks))
       do ib=1,op1%nblocks
          call op1%blocks(ib)%clone(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'serial_block_array_t%clone: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_clone

  ! op <- comm(op)
  subroutine serial_block_array_comm(op)
    implicit none
    class(serial_block_array_t), intent(inout) :: op 
  end subroutine serial_block_array_comm

  subroutine serial_block_array_free(this)
    implicit none
    class(serial_block_array_t), intent(inout) :: this
    integer(ip)  :: ib
   
    do ib=1, this%nblocks
       call this%blocks(ib)%free()
    end do
    this%nblocks = 0
    deallocate( this%blocks )
  end subroutine serial_block_array_free

end module serial_block_array_names
