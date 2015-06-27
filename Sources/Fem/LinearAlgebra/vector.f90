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
module vector_names
  use types_names
  use memor_names
#ifdef ENABLE_BLAS
  use blas77_interfaces_names
#endif
  use base_operand_names

!!$#ifdef memcheck
!!$  use iso_c_binding
!!$#endif
  implicit none
# include "debug.i90"

  !=============================================================
  ! TODO:
  ! 
  ! x Call to BLAS double precision or single precision 
  !   subroutines depending on the value of the rp parameter. 
  !   Currently we are always calling double precision variants 
  !   of the BLAS subroutines.
  ! 
  !=============================================================
  integer(ip), parameter :: not_created = 0 ! vector%b points to null
  integer(ip), parameter :: allocated   = 1 ! vector%b has been allocated
  integer(ip), parameter :: reference   = 2 ! vector%b references external memory

  private

  ! vector
  type, extends(base_operand_t) :: vector_t
     integer(ip)                :: &
        neq = 0                       ! Number of equations
   
     integer(ip)                :: & 
        mode = not_created           

     real(rp), pointer          :: &
        b(:) => NULL()
   contains
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => vector_dot_tbp
     procedure :: copy => vector_copy_tbp
     procedure :: init => vector_init_tbp
     procedure :: scal => vector_scal_tbp
     procedure :: axpby => vector_axpby_tbp
     procedure :: nrm2 => vector_nrm2_tbp
     procedure :: clone => vector_clone_tbp
     procedure :: comm  => vector_comm_tbp
     procedure :: free  => vector_free_tbp
  end type vector_t

  ! interface vector_assembly
  !    module procedure vector_assembly_w_dof_handler
  !    module procedure vector_assembly_nosq_w_dof_handler
  ! end interface vector_assembly

!  interface vector_create_view
!     module procedure vector_create_view_same_ndof, & 
!                      vector_create_view_new_ndof
!  end interface vector_create_view

!!$  ! Memalloc interfaces
!!$  interface memalloc
!!$     module procedure memalloc_vector1, memalloc_vector2, memalloc_vector3, &
!!$          &           memalloc_vector4
!!$  end interface memalloc
!!$  interface memrealloc
!!$     module procedure memrealloc_vector1, memrealloc_vector2, memrealloc_vector3, &
!!$          &           memrealloc_vector4
!!$  end interface memrealloc
!!$  interface memfree
!!$     module procedure memfree_vector1, memfree_vector2, memfree_vector3, &
!!$          &           memfree_vector4
!!$  end interface memfree
!!$  interface memmovealloc
!!$     module procedure memmovealloc_vector1, memmovealloc_vector2, memmovealloc_vector3, &
!!$          &           memmovealloc_vector4
!!$  end interface memmovealloc

  ! Types
  public :: vector_t

  ! Constants 
  public :: reference

  ! Functions
  public :: vector_free, &
            vector_alloc, vector_clone, vector_create_view ,           & 
            vector_comm,                                                       &
            vector_dot, vector_nrm2, vector_copy,                      & 
            vector_zero, vector_init, vector_scale, vector_mxpy,   & 
            vector_axpy, vector_aypx, vector_pxpy, vector_pxmy,    & 
            vector_print, vector_print_matrix_market!,                      &
            !memalloc, memrealloc, memfreep, memmovealloc
contains

!!$  !***********************************************************************
!!$  !***********************************************************************
!!$  ! Specialization to allocate type(vector_t)
!!$  !***********************************************************************
!!$  !***********************************************************************
!!$# define var_type type(vector_t)
!!$# define var_size 24
!!$# define var_attr allocatable,target
!!$# define point(a,b) call move_alloc(a,b)
!!$# define bound_kind ip
!!$
!!$# define generic_status_test     allocated
!!$# define generic_memalloc_1      memalloc_vector1    
!!$# define generic_memalloc_2      memalloc_vector2    
!!$# define generic_memalloc_3      memalloc_vector3    
!!$# define generic_memalloc_4      memalloc_vector4    
!!$# define generic_memrealloc_1    memrealloc_vector1  
!!$# define generic_memrealloc_2    memrealloc_vector2  
!!$# define generic_memrealloc_3    memrealloc_vector3  
!!$# define generic_memrealloc_4    memrealloc_vector4  
!!$# define generic_memfree_1       memfree_vector1     
!!$# define generic_memfree_2       memfree_vector2     
!!$# define generic_memfree_3       memfree_vector3     
!!$# define generic_memfree_4       memfree_vector4     
!!$# define generic_memmovealloc_1  memmovealloc_vector1
!!$# define generic_memmovealloc_2  memmovealloc_vector2
!!$# define generic_memmovealloc_3  memmovealloc_vector3
!!$# define generic_memmovealloc_4  memmovealloc_vector4
!!$# include "memor.i90"

  !=============================================================================
  subroutine vector_free (vec)
    implicit none
    type(vector_t), intent(inout) :: vec
    assert (vec%mode == allocated .or. vec%mode == reference)
    vec%neq     = 0          ! Number of equations
    if (vec%mode == allocated) call memfreep(vec%b,__FILE__,__LINE__)
    vec%mode = not_created
  end subroutine vector_free

  !=============================================================================
  subroutine vector_alloc(neq,vec)
    implicit none
    integer(ip)     , intent(in)    :: neq
    type(vector_t), intent(inout) :: vec
    assert ( vec%mode == not_created )
    vec%neq     = neq  ! Number of equations
    call memallocp(vec%neq,vec%b,__FILE__,__LINE__)
    vec%b    = 0.0_rp
    vec%mode = allocated
  end subroutine vector_alloc
  
  subroutine vector_create_view (svec, start, end, tvec)
    implicit none
    type(vector_t), intent(in), target  :: svec
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(vector_t), intent(inout)       :: tvec

    assert ( tvec%mode == not_created )

    tvec%neq     =  end-start+1           ! Number of equations
    tvec%b => svec%b(start:end)
    tvec%mode =  reference
  end subroutine vector_create_view

  ! ! This subroutine is required in order to create a view which has a 
  ! ! different ndof of an existing vector with a given ndof
  ! subroutine vector_create_view_new_ndof (svec, ndstart, ndend, start, end, tvec)
  !   implicit none
  !   ! Parameters
  !   type(vector_t), intent(in), target  :: svec
  !   integer(ip)     , intent(in)          :: ndstart
  !   integer(ip)     , intent(in)          :: ndend
  !   integer(ip)     , intent(in)          :: start
  !   integer(ip)     , intent(in)          :: end
  !   type(vector_t), intent(inout)       :: tvec

  !   ! Locals
  !   integer(ip) :: off, start_aux, end_aux

  !   assert ( tvec%mode == not_created )
      
  !   tvec%neq     =  end-start+1           ! Number of equations

  !   start_aux = mod(start-1, svec%neq)+1
  !   end_aux   = mod(end-1  , svec%neq)+1

  !     ! Compute initial offset
  !     off = (ndstart-1) * svec%neq
  !     tvec%b => svec%b(:, off+start_aux:off+end_aux)
    
  !   tvec%mode =  reference
  ! end subroutine vector_create_view_new_ndof



  subroutine vector_clone ( svec, tvec )
    implicit none
    type(vector_t), intent( in ) :: svec
    type(vector_t), intent(out) :: tvec
    tvec%neq     =  svec%neq       ! Number of equations
      call memallocp(tvec%neq,tvec%b,__FILE__,__LINE__)
    tvec%b = 0.0_rp
    tvec%mode = allocated 
  end subroutine vector_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine vector_comm ( vec )
    implicit none
    type(vector_t), intent( inout ) :: vec 
  end subroutine vector_comm

  !=============================================================================
  subroutine vector_dot (x, y, t)
    implicit none
    type(vector_t), intent(in)  :: x
    type(vector_t), intent(in)  :: y
    real(rp)        , intent(out) :: t

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    t = ddot( x%neq, x%b, 1, y%b, 1 )
#else
!!$    AFM: A non BLAS-based implementation of the
!!$    dot product should go here
    check(1==0)
#endif

  end subroutine vector_dot
  !=============================================================================
  subroutine vector_nrm2(x,t)
    implicit none
    type(vector_t), intent(in)  :: x
    real(rp)    , intent(out)     :: t

#ifdef ENABLE_BLAS
    t = dnrm2( x%neq, x%b, 1 )
#else
    call vector_dot (x, x, t)
    t = sqrt(t)
#endif
  end subroutine vector_nrm2
  !=============================================================================
  subroutine vector_copy(x,y)
    implicit none
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    call dcopy ( x%neq, x%b, 1, y%b, 1 ) 
#else
    y%b=x%b
#endif
  end subroutine vector_copy
  subroutine vector_zero(y)
    implicit none
    type(vector_t), intent(inout) :: y
    y%b=0.0_rp
  end subroutine vector_zero

  subroutine vector_init(alpha, y)
    implicit none
    type(vector_t), intent(inout) :: y 
    real(rp), intent(in)            :: alpha  
    y%b=alpha
  end subroutine vector_init
  
  subroutine vector_scale(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y

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
  end subroutine vector_scale

  subroutine vector_mxpy(x,y)
    implicit none
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, -1.0, x%b, 1, y%b, 1 )
#else
    y%b=y%b-x%b
#endif
  end subroutine vector_mxpy
  subroutine vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, t, x%b, 1, y%b, 1 )
#else
    y%b=y%b+t*x%b
#endif
  end subroutine vector_axpy

  subroutine vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y
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
  end subroutine vector_aypx

  subroutine vector_pxpy(x,y)
    implicit none
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=y%b+x%b
#endif
  end subroutine vector_pxpy

  subroutine vector_pxmy(x,y)
    implicit none
    type(vector_t), intent(in)    :: x
    type(vector_t), intent(inout) :: y
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
  end subroutine vector_pxmy

  subroutine vector_print (luout, x)
    implicit none
    type(vector_t), intent(in) :: x
    integer(ip),      intent(in) :: luout
    
    ! Locals
    integer(ip) :: i

    write (luout, '(a)')     '*** begin vector data structure ***'
    write(luout,'(a,i10)') 'size', x%neq
       write (luout,'(e25.16)') x%b
!!$    do i=1,x%neq
!!$!       write (luout,'(e14.7)') x%b(1,i)
!!$       write (luout,'(e25.16,1x)') x%b(1,i)
!!$    end do

  end subroutine vector_print

  subroutine vector_print_matrix_market ( luout, x )
   implicit none
   ! Parameters
   type(vector_t), intent(in) :: x
   integer(ip),      intent(in) :: luout

   ! Locals
   integer (ip) :: i, id

   write (luout,'(a)') '%%MatrixMarket matrix array real general'
   write (luout,*) x%neq , 1
   do i=1,x%neq 
            write (luout,*) x%b( i )
   end do

 end subroutine vector_print_matrix_market

 ! alpha <- op1^T * op2
 function vector_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(vector_t), intent(in)    :: op1
   class(base_operand_t), intent(in)  :: op2
   real(rp) :: alpha

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (vector_t)
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
 end function vector_dot_tbp

 ! op1 <- op2 
 subroutine vector_copy_tbp(op1,op2)
   implicit none
   class(vector_t), intent(inout) :: op1
   class(base_operand_t), intent(in)  :: op2
   
   call op2%GuardTemp()
   select type(op2)
   class is (vector_t)
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
 end subroutine vector_copy_tbp

 ! op1 <- alpha * op2
 subroutine vector_scal_tbp(op1,alpha,op2)
   implicit none
   class(vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand_t), intent(in) :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (vector_t)
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
 end subroutine vector_scal_tbp
 ! op <- alpha
 subroutine vector_init_tbp(op,alpha)
   implicit none
   class(vector_t), intent(inout) :: op
   real(rp), intent(in) :: alpha
   op%b=alpha
 end subroutine vector_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine vector_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(vector_t), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand_t), intent(in) :: op2
   real(rp), intent(in) :: beta

   call op2%GuardTemp()
   select type(op2)
   class is (vector_t)
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
 end subroutine vector_axpby_tbp

 ! alpha <- nrm2(op)
 function vector_nrm2_tbp(op) result(alpha)
   implicit none
   class(vector_t), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()

#ifdef ENABLE_BLAS
    alpha = dnrm2( op%neq, op%b, 1 )
#else
    alpha = op%dot(op)
    alpha = sqrt(alpha)
#endif

   call op%CleanTemp()
 end function vector_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine vector_clone_tbp(op1,op2)
   implicit none
   class(vector_t)           ,intent(inout) :: op1
   class(base_operand_t), target ,intent(in)    :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (vector_t)
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
 end subroutine vector_clone_tbp

 ! op <- comm(op)
 subroutine vector_comm_tbp(op)
   implicit none
   class(vector_t), intent(inout) :: op
 end subroutine vector_comm_tbp

 subroutine vector_free_tbp(this)
   implicit none
   class(vector_t), intent(inout) :: this

   this%neq     = 0          ! Number of equations
   if (this%mode == allocated) call memfreep(this%b,__FILE__,__LINE__)
   nullify(this%b)
   this%mode = not_created
 end subroutine vector_free_tbp

end module vector_names
