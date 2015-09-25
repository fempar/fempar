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
module serial_scalar_array_names
  use types_names
  use memor_names
#ifdef ENABLE_BLAS
  use blas77_interfaces_names
#endif
  use array_names
  use vector_names

  implicit none
# include "debug.i90"

  integer(ip), parameter :: not_created = 0 ! vector%b points to null
  integer(ip), parameter :: allocated   = 1 ! vector%b has been allocated
  integer(ip), parameter :: reference   = 2 ! vector%b references external memory

  private

  type, extends(array_t) :: serial_scalar_array_t
     integer(ip)                :: &
        neq = 0                       ! Number of equations
   
     integer(ip)                :: & 
        mode = not_created           

     real(rp), pointer          :: &
        b(:) => NULL()
   contains
   
	 procedure :: create_and_allocate => serial_scalar_array_create_and_allocate
	 procedure :: create              => serial_scalar_array_create
	 procedure :: allocate            => serial_scalar_array_allocate
	 procedure :: create_view         => serial_scalar_array_create_view
	 procedure :: print               => serial_scalar_array_print
	 procedure :: print_matrix_market => serial_scalar_array_print_matrix_market
   
     procedure :: dot  => serial_scalar_array_dot
     procedure :: copy => serial_scalar_array_copy
     procedure :: init => serial_scalar_array_init
     procedure :: scal => serial_scalar_array_scal
     procedure :: axpby => serial_scalar_array_axpby
     procedure :: nrm2 => serial_scalar_array_nrm2
     procedure :: clone => serial_scalar_array_clone
     procedure :: comm  => serial_scalar_array_comm
	 procedure :: free_in_stages  => serial_scalar_array_free_in_stages
     procedure :: default_initialization => serial_scalar_array_default_init
  end type serial_scalar_array_t

  ! Types
  public :: serial_scalar_array_t

  ! Constants 
  public :: reference

contains

  !=============================================================================
  subroutine serial_scalar_array_default_init (this)
    implicit none
    class(serial_scalar_array_t), intent(inout) :: this
    this%neq     = 0          ! Number of equations
    this%mode = not_created
    nullify(this%b)
    call this%NullifyTemporary()
  end subroutine serial_scalar_array_default_init
  
  !=============================================================================
  subroutine serial_scalar_array_create_and_allocate(this,size)
    implicit none
	class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: size
    call this%create(size)
	call this%allocate()
  end subroutine serial_scalar_array_create_and_allocate
  
  !=============================================================================
  subroutine serial_scalar_array_create(this,size)
    implicit none
	class(serial_scalar_array_t), intent(inout) :: this
    integer(ip)                 , intent(in)    :: size
    assert ( this%mode == not_created )
    this%neq  = size  
	this%mode = allocated
  end subroutine serial_scalar_array_create
  
  !=============================================================================
  subroutine serial_scalar_array_allocate(this)
    implicit none
	class(serial_scalar_array_t), intent(inout) :: this
    call memallocp(this%neq,this%b,__FILE__,__LINE__)
    this%b    = 0.0_rp
  end subroutine serial_scalar_array_allocate
  
  subroutine serial_scalar_array_create_view (this, start, end, tvec)
    implicit none
    class(serial_scalar_array_t), intent(in), target  :: this
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(serial_scalar_array_t), intent(inout)       :: tvec

    assert ( tvec%mode == not_created )

    tvec%neq =  end-start+1 ! Number of equations
    tvec%b => this%b(start:end)
    tvec%mode =  reference
  end subroutine serial_scalar_array_create_view

  subroutine serial_scalar_array_print (this, luout)
    implicit none
    class(serial_scalar_array_t), intent(in) :: this
    integer(ip)                , intent(in) :: luout

    write (luout, '(a)')     '*** begin serial_scalar_array data structure ***'
    write(luout,'(a,i10)') 'size', this%neq
    write (luout,'(e25.16)') this%b
	write (luout, '(a)')     '*** end serial_scalar_array data structure ***'
  end subroutine serial_scalar_array_print

  subroutine serial_scalar_array_print_matrix_market ( this, luout )
   implicit none
   class(serial_scalar_array_t), intent(in) :: this
   integer(ip)                , intent(in) :: luout
   integer (ip) :: i
   write (luout,'(a)') '%%MatrixMarket matrix array real general'
   write (luout,*) this%neq , 1
   do i=1,this%neq 
     write (luout,*) this%b( i )
   end do
 end subroutine serial_scalar_array_print_matrix_market

 ! alpha <- op1^T * op2
 function serial_scalar_array_dot(op1,op2) result(alpha)
   implicit none
   class(serial_scalar_array_t), intent(in)    :: op1
   class(vector_t), intent(in)  :: op2
   real(rp) :: alpha

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      assert ( op1%neq == op2%neq )
#ifdef ENABLE_BLAS
      alpha = ddot( op1%neq, op1%b, 1, op2%b, 1 )
#else
    check(1==0)
#endif
   class default
      write(0,'(a)') 'serial_scalar_array_t%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function serial_scalar_array_dot

 ! op1 <- op2 
 subroutine serial_scalar_array_copy(op1,op2)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   class(vector_t), intent(in)  :: op2
   
   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      assert ( op2%neq == op1%neq )
#ifdef ENABLE_BLAS
      call dcopy ( op2%neq, op2%b, 1, op1%b, 1 ) 
#else
      op1%b=op2%b
#endif
   class default
      write(0,'(a)') 'serial_scalar_array_t%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_copy

 ! op1 <- alpha * op2
 subroutine serial_scalar_array_scal(op1,alpha,op2)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      assert ( op2%neq == op1%neq )
#ifdef ENABLE_BLAS
      call dcopy ( op2%neq, op2%b, 1, op1%b, 1)
      call dscal ( op1%neq, alpha, op1%b, 1)
#else
      op1%b=alpha*op2%b
#endif
   class default
      write(0,'(a)') 'serial_scalar_array_t%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_scal
 ! op <- alpha
 subroutine serial_scalar_array_init(op,alpha)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   op%b=alpha
 end subroutine serial_scalar_array_init

 ! op1 <- alpha*op2 + beta*op1
 subroutine serial_scalar_array_axpby(op1, alpha, op2, beta)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(vector_t), intent(in) :: op2
   real(rp), intent(in) :: beta

   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      assert ( op2%neq == op1%neq )
      if ( beta == 0.0_rp ) then
         call op1%scal(alpha, op2)
      else if ( beta == 1.0_rp ) then
         ! AXPY
#ifdef ENABLE_BLAS
         call daxpy ( op2%neq, alpha, op2%b, 1, op1%b, 1 )
#else
         op1%b=op1%b+alpha*op2%b
#endif
      else
         ! SCAL + AXPY
         call op1%scal(beta, op1)
#ifdef ENABLE_BLAS
         call daxpy ( op2%neq, alpha, op2%b, 1, op1%b, 1 )    
#else
         op1%b=op1%b+alpha*op2%b
#endif  
      end if
   class default
      write(0,'(a)') 'serial_scalar_array_t%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_axpby

 ! alpha <- nrm2(op)
 function serial_scalar_array_nrm2(op) result(alpha)
   implicit none
   class(serial_scalar_array_t), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()

#ifdef ENABLE_BLAS
    alpha = dnrm2( op%neq, op%b, 1 )
#else
    alpha = op%dot(op)
    alpha = sqrt(alpha)
#endif

   call op%CleanTemp()
 end function serial_scalar_array_nrm2

 ! op1 <- clone(op2) 
 subroutine serial_scalar_array_clone(op1,op2)
   implicit none
   class(serial_scalar_array_t)           ,intent(inout) :: op1
   class(vector_t), target ,intent(in)    :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      if (op1%mode == allocated) call memfreep(op1%b,__FILE__,__LINE__)
      op1%neq     =  op2%neq
      call memallocp(op1%neq,op1%b,__FILE__,__LINE__)
      op1%mode = allocated
   class default
      write(0,'(a)') 'serial_scalar_array_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_clone

 ! op <- comm(op)
 subroutine serial_scalar_array_comm(op)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op
 end subroutine serial_scalar_array_comm
 
  subroutine serial_scalar_array_free_in_stages(this,action)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: this
   integer(ip)                 , intent(in)    :: action
   assert ( action == free_clean .or. action == free_struct .or. action == free_values )
   if ( action == free_clean ) then
     ! Undo create
     this%neq = 0
	 this%mode = not_created
   else if ( action == free_values ) then
     ! Undo allocate
     if (this%mode == allocated) call memfreep(this%b,__FILE__,__LINE__)
     nullify(this%b)
   end if
 end subroutine serial_scalar_array_free_in_stages

end module serial_scalar_array_names
