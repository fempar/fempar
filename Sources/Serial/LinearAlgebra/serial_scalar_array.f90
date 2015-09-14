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
  use abstract_vector_names

  implicit none
# include "debug.i90"

  integer(ip), parameter :: not_created = 0 ! vector%b points to null
  integer(ip), parameter :: allocated   = 1 ! vector%b has been allocated
  integer(ip), parameter :: reference   = 2 ! vector%b references external memory

  private

  ! vector
  type, extends(abstract_vector_t) :: serial_scalar_array_t
     integer(ip)                :: &
        neq = 0                       ! Number of equations
   
     integer(ip)                :: & 
        mode = not_created           

     real(rp), pointer          :: &
        b(:) => NULL()
   contains
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => serial_scalar_array_dot_tbp
     procedure :: copy => serial_scalar_array_copy_tbp
     procedure :: init => serial_scalar_array_init_tbp
     procedure :: scal => serial_scalar_array_scal_tbp
     procedure :: axpby => serial_scalar_array_axpby_tbp
     procedure :: nrm2 => serial_scalar_array_nrm2_tbp
     procedure :: clone => serial_scalar_array_clone_tbp
     procedure :: comm  => serial_scalar_array_comm_tbp
     procedure :: free  => serial_scalar_array_free_tbp
     procedure :: default_initialization => serial_scalar_array_default_init
  end type serial_scalar_array_t

  ! Types
  public :: serial_scalar_array_t

  ! Constants 
  public :: reference

  ! Functions
  public :: serial_scalar_array_free, &
            serial_scalar_array_alloc, serial_scalar_array_clone, serial_scalar_array_create_view ,           & 
            serial_scalar_array_dot, serial_scalar_array_copy,                      & 
            serial_scalar_array_zero, serial_scalar_array_init, serial_scalar_array_scale, serial_scalar_array_mxpy,   & 
            serial_scalar_array_axpy, serial_scalar_array_aypx, serial_scalar_array_pxpy, serial_scalar_array_pxmy,    & 
            serial_scalar_array_print, serial_scalar_array_print_matrix_market
contains

  !=============================================================================
  subroutine serial_scalar_array_free (vec)
    implicit none
    type(serial_scalar_array_t), intent(inout) :: vec
    assert (vec%mode == allocated .or. vec%mode == reference)
    vec%neq     = 0          ! Number of equations
    if (vec%mode == allocated) call memfreep(vec%b,__FILE__,__LINE__)
    vec%mode = not_created
  end subroutine serial_scalar_array_free

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
  subroutine serial_scalar_array_alloc(neq,vec)
    implicit none
    integer(ip)     , intent(in)    :: neq
    type(serial_scalar_array_t), intent(inout) :: vec
    assert ( vec%mode == not_created )
    vec%neq     = neq  ! Number of equations
    call memallocp(vec%neq,vec%b,__FILE__,__LINE__)
    vec%b    = 0.0_rp
    vec%mode = allocated
  end subroutine serial_scalar_array_alloc
  
  subroutine serial_scalar_array_create_view (svec, start, end, tvec)
    implicit none
    type(serial_scalar_array_t), intent(in), target  :: svec
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(serial_scalar_array_t), intent(inout)       :: tvec

    assert ( tvec%mode == not_created )

    tvec%neq     =  end-start+1           ! Number of equations
    tvec%b => svec%b(start:end)
    tvec%mode =  reference
  end subroutine serial_scalar_array_create_view

  subroutine serial_scalar_array_clone ( svec, tvec )
    implicit none
    type(serial_scalar_array_t), intent( in ) :: svec
    type(serial_scalar_array_t), intent(out) :: tvec
    tvec%neq     =  svec%neq       ! Number of equations
    call memallocp(tvec%neq,tvec%b,__FILE__,__LINE__)
    tvec%b = 0.0_rp
    tvec%mode = allocated 
  end subroutine serial_scalar_array_clone

  !=============================================================================
  subroutine serial_scalar_array_dot (x, y, t)
    implicit none
    type(serial_scalar_array_t), intent(in)  :: x
    type(serial_scalar_array_t), intent(in)  :: y
    real(rp)        , intent(out) :: t

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    t = ddot( x%neq, x%b, 1, y%b, 1 )
#else
!!$    AFM: A non BLAS-based implementation of the
!!$    dot product should go here
    check(1==0)
#endif

  end subroutine serial_scalar_array_dot

  !=============================================================================
  subroutine serial_scalar_array_copy(x,y)
    implicit none
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    call dcopy ( x%neq, x%b, 1, y%b, 1 ) 
#else
    y%b=x%b
#endif
  end subroutine serial_scalar_array_copy
  subroutine serial_scalar_array_zero(y)
    implicit none
    type(serial_scalar_array_t), intent(inout) :: y
    y%b=0.0_rp
  end subroutine serial_scalar_array_zero

  subroutine serial_scalar_array_init(alpha, y)
    implicit none
    type(serial_scalar_array_t), intent(inout) :: y 
    real(rp), intent(in)            :: alpha  
    y%b=alpha
  end subroutine serial_scalar_array_init
  
  subroutine serial_scalar_array_scale(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dcopy ( x%neq, x%b, 1, y%b, 1)
    call dscal ( y%neq, t, y%b, 1)
#else
    y%b=t*x%b
#endif
  end subroutine serial_scalar_array_scale

  subroutine serial_scalar_array_mxpy(x,y)
    implicit none
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, -1.0, x%b, 1, y%b, 1 )
#else
    y%b=y%b-x%b
#endif
  end subroutine serial_scalar_array_mxpy
  subroutine serial_scalar_array_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, t, x%b, 1, y%b, 1 )
#else
    y%b=y%b+t*x%b
#endif
  end subroutine serial_scalar_array_axpy

  subroutine serial_scalar_array_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dscal ( y%neq, t, y%b, 1)
    call daxpy ( x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=x%b+t*y%b
#endif
  end subroutine serial_scalar_array_aypx

  subroutine serial_scalar_array_pxpy(x,y)
    implicit none
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=y%b+x%b
#endif
  end subroutine serial_scalar_array_pxpy

  subroutine serial_scalar_array_pxmy(x,y)
    implicit none
    type(serial_scalar_array_t), intent(in)    :: x
    type(serial_scalar_array_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dscal ( y%neq, -1.0, y%b, 1)
    call daxpy ( x%neq, 1.0, x%b, 1, y%b, 1 ) 
#else
    y%b=x%b-y%b
#endif
  end subroutine serial_scalar_array_pxmy

  subroutine serial_scalar_array_print (luout, x)
    implicit none
    type(serial_scalar_array_t), intent(in) :: x
    integer(ip),      intent(in) :: luout
    
    ! Locals
    integer(ip) :: i

    write (luout, '(a)')     '*** begin vector data structure ***'
    write(luout,'(a,i10)') 'size', x%neq
    write (luout,'(e25.16)') x%b
	
  end subroutine serial_scalar_array_print

  subroutine serial_scalar_array_print_matrix_market ( luout, x )
   implicit none
   ! Parameters
   type(serial_scalar_array_t), intent(in) :: x
   integer(ip),      intent(in) :: luout

   ! Locals
   integer (ip) :: i, id

   write (luout,'(a)') '%%MatrixMarket matrix array real general'
   write (luout,*) x%neq , 1
   do i=1,x%neq 
     write (luout,*) x%b( i )
   end do

 end subroutine serial_scalar_array_print_matrix_market

 ! alpha <- op1^T * op2
 function serial_scalar_array_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(serial_scalar_array_t), intent(in)    :: op1
   class(abstract_vector_t), intent(in)  :: op2
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
      write(0,'(a)') 'vector_t%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function serial_scalar_array_dot_tbp

 ! op1 <- op2 
 subroutine serial_scalar_array_copy_tbp(op1,op2)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   class(abstract_vector_t), intent(in)  :: op2
   
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
      write(0,'(a)') 'vector_t%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_copy_tbp

 ! op1 <- alpha * op2
 subroutine serial_scalar_array_scal_tbp(op1,alpha,op2)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(abstract_vector_t), intent(in) :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      assert ( op2%neq == op1%neq )
#ifdef ENABLE_BLAS
      ! I guess that two calls to the level 1
      ! BLAS can not be competitive against 
      ! just one F90 vector operation. I have to
      ! measure the difference among these two
      ! options. 
      call dcopy ( op2%neq, op2%b, 1, op1%b, 1)
      call dscal ( op1%neq, alpha, op1%b, 1)
#else
      op1%b=alpha*op2%b
#endif
   class default
      write(0,'(a)') 'vector_t%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_scal_tbp
 ! op <- alpha
 subroutine serial_scalar_array_init_tbp(op,alpha)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   op%b=alpha
 end subroutine serial_scalar_array_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine serial_scalar_array_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(abstract_vector_t), intent(in) :: op2
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
      write(0,'(a)') 'vector_t%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_axpby_tbp

 ! alpha <- nrm2(op)
 function serial_scalar_array_nrm2_tbp(op) result(alpha)
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
 end function serial_scalar_array_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine serial_scalar_array_clone_tbp(op1,op2)
   implicit none
   class(serial_scalar_array_t)           ,intent(inout) :: op1
   class(abstract_vector_t), target ,intent(in)    :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (serial_scalar_array_t)
      if (op1%mode == allocated) call memfreep(op1%b,__FILE__,__LINE__)
      op1%neq     =  op2%neq       ! Number of equations
      call memallocp(op1%neq,op1%b,__FILE__,__LINE__)
      ! AFM: I think that clone should NOT init the memory just allocated.
      ! The code that surrounds clone (e.g., Krylov solvers) should not
      ! rely on vector_clone_tbp initializing the memory. I will comment
      ! it out, and expect that the codes continue working.
      ! op1%b = 0.0_rp 
      op1%mode = allocated
   class default
      write(0,'(a)') 'vector_t%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine serial_scalar_array_clone_tbp

 ! op <- comm(op)
 subroutine serial_scalar_array_comm_tbp(op)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: op
 end subroutine serial_scalar_array_comm_tbp

 subroutine serial_scalar_array_free_tbp(this)
   implicit none
   class(serial_scalar_array_t), intent(inout) :: this

   this%neq     = 0          ! Number of equations
   if (this%mode == allocated) call memfreep(this%b,__FILE__,__LINE__)
   nullify(this%b)
   this%mode = not_created
 end subroutine serial_scalar_array_free_tbp

end module serial_scalar_array_names
