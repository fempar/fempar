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
module fem_vector_names
  use types
  use memor
#ifdef ENABLE_BLAS
  use blas77_interfaces
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
  integer(ip), parameter :: not_created = 0 ! fem_vector%b points to null
  integer(ip), parameter :: allocated   = 1 ! fem_vector%b has been allocated
  integer(ip), parameter :: reference   = 2 ! fem_vector%b references external memory

  private

  ! fem_vector
  type, extends(base_operand) :: fem_vector
     integer(ip)                :: &
        neq = 0                       ! Number of equations
   
     integer(ip)                :: & 
        mode = not_created           

     real(rp), pointer          :: &
        b(:) => NULL()
   contains
     ! Provide type bound procedures (tbp) implementors
     procedure :: dot  => fem_vector_dot_tbp
     procedure :: copy => fem_vector_copy_tbp
     procedure :: init => fem_vector_init_tbp
     procedure :: scal => fem_vector_scal_tbp
     procedure :: axpby => fem_vector_axpby_tbp
     procedure :: nrm2 => fem_vector_nrm2_tbp
     procedure :: clone => fem_vector_clone_tbp
     procedure :: comm  => fem_vector_comm_tbp
     procedure :: free  => fem_vector_free_tbp
  end type fem_vector

  ! interface fem_vector_assembly
  !    module procedure fem_vector_assembly_w_dof_handler
  !    module procedure fem_vector_assembly_nosq_w_dof_handler
  ! end interface fem_vector_assembly

!  interface fem_vector_create_view
!     module procedure fem_vector_create_view_same_ndof, & 
!                      fem_vector_create_view_new_ndof
!  end interface fem_vector_create_view

!!$  ! Memalloc interfaces
!!$  interface memalloc
!!$     module procedure memalloc_fem_vector1, memalloc_fem_vector2, memalloc_fem_vector3, &
!!$          &           memalloc_fem_vector4
!!$  end interface memalloc
!!$  interface memrealloc
!!$     module procedure memrealloc_fem_vector1, memrealloc_fem_vector2, memrealloc_fem_vector3, &
!!$          &           memrealloc_fem_vector4
!!$  end interface memrealloc
!!$  interface memfree
!!$     module procedure memfree_fem_vector1, memfree_fem_vector2, memfree_fem_vector3, &
!!$          &           memfree_fem_vector4
!!$  end interface memfree
!!$  interface memmovealloc
!!$     module procedure memmovealloc_fem_vector1, memmovealloc_fem_vector2, memmovealloc_fem_vector3, &
!!$          &           memmovealloc_fem_vector4
!!$  end interface memmovealloc

  ! Types
  public :: fem_vector

  ! Constants 
  public :: reference

  ! Functions
  public :: fem_vector_free, &
            fem_vector_alloc, fem_vector_clone, fem_vector_create_view ,           & 
            fem_vector_comm,                                                       &
            fem_vector_dot, fem_vector_nrm2, fem_vector_copy,                      & 
            fem_vector_zero, fem_vector_init, fem_vector_scale, fem_vector_mxpy,   & 
            fem_vector_axpy, fem_vector_aypx, fem_vector_pxpy, fem_vector_pxmy,    & 
            fem_vector_print, fem_vector_print_matrix_market!,                      &
            !memalloc, memrealloc, memfreep, memmovealloc
contains

!!$  !***********************************************************************
!!$  !***********************************************************************
!!$  ! Specialization to allocate type(fem_vector)
!!$  !***********************************************************************
!!$  !***********************************************************************
!!$# define var_type type(fem_vector)
!!$# define var_size 24
!!$# define var_attr allocatable,target
!!$# define point(a,b) call move_alloc(a,b)
!!$# define bound_kind ip
!!$
!!$# define generic_status_test     allocated
!!$# define generic_memalloc_1      memalloc_fem_vector1    
!!$# define generic_memalloc_2      memalloc_fem_vector2    
!!$# define generic_memalloc_3      memalloc_fem_vector3    
!!$# define generic_memalloc_4      memalloc_fem_vector4    
!!$# define generic_memrealloc_1    memrealloc_fem_vector1  
!!$# define generic_memrealloc_2    memrealloc_fem_vector2  
!!$# define generic_memrealloc_3    memrealloc_fem_vector3  
!!$# define generic_memrealloc_4    memrealloc_fem_vector4  
!!$# define generic_memfree_1       memfree_fem_vector1     
!!$# define generic_memfree_2       memfree_fem_vector2     
!!$# define generic_memfree_3       memfree_fem_vector3     
!!$# define generic_memfree_4       memfree_fem_vector4     
!!$# define generic_memmovealloc_1  memmovealloc_fem_vector1
!!$# define generic_memmovealloc_2  memmovealloc_fem_vector2
!!$# define generic_memmovealloc_3  memmovealloc_fem_vector3
!!$# define generic_memmovealloc_4  memmovealloc_fem_vector4
!!$# include "memor.i90"

  !=============================================================================
  subroutine fem_vector_free (vec)
    implicit none
    type(fem_vector), intent(inout) :: vec
    assert (vec%mode == allocated .or. vec%mode == reference)
    vec%neq     = 0          ! Number of equations
    if (vec%mode == allocated) call memfreep(vec%b,__FILE__,__LINE__)
    vec%mode = not_created
  end subroutine fem_vector_free

  !=============================================================================
  subroutine fem_vector_alloc(neq,vec)
    implicit none
    integer(ip)     , intent(in)    :: neq
    type(fem_vector), intent(inout) :: vec
    assert ( vec%mode == not_created )
    vec%neq     = neq  ! Number of equations
    call memallocp(vec%neq,vec%b,__FILE__,__LINE__)
    vec%b    = 0.0_rp
    vec%mode = allocated
  end subroutine fem_vector_alloc
  
  subroutine fem_vector_create_view (svec, start, end, tvec)
    implicit none
    type(fem_vector), intent(in), target  :: svec
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(fem_vector), intent(inout)       :: tvec

    assert ( tvec%mode == not_created )

    tvec%neq     =  end-start+1           ! Number of equations
    tvec%b => svec%b(start:end)
    tvec%mode =  reference
  end subroutine fem_vector_create_view

  ! ! This subroutine is required in order to create a view which has a 
  ! ! different ndof of an existing fem_vector with a given ndof
  ! subroutine fem_vector_create_view_new_ndof (svec, ndstart, ndend, start, end, tvec)
  !   implicit none
  !   ! Parameters
  !   type(fem_vector), intent(in), target  :: svec
  !   integer(ip)     , intent(in)          :: ndstart
  !   integer(ip)     , intent(in)          :: ndend
  !   integer(ip)     , intent(in)          :: start
  !   integer(ip)     , intent(in)          :: end
  !   type(fem_vector), intent(inout)       :: tvec

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
  ! end subroutine fem_vector_create_view_new_ndof



  subroutine fem_vector_clone ( svec, tvec )
    implicit none
    type(fem_vector), intent( in ) :: svec
    type(fem_vector), intent(out) :: tvec
    tvec%neq     =  svec%neq       ! Number of equations
      call memallocp(tvec%neq,tvec%b,__FILE__,__LINE__)
    tvec%b = 0.0_rp
    tvec%mode = allocated 
  end subroutine fem_vector_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine fem_vector_comm ( vec )
    implicit none
    type(fem_vector), intent( inout ) :: vec 
  end subroutine fem_vector_comm

  ! !=============================================================================
  ! subroutine fem_vector_assembly_w_dof_handler(nn, dofh, el, l2g, ev, vec)
  !   implicit none
  !   type(dof_handler), intent(in)    :: dofh
  !   type(fem_element), intent(in)    :: el
  !   type(elvec) ,      intent(in)    :: ev
  !   integer(ip) ,      intent(in)    :: nn(:)
  !   type(fem_vector),  intent(inout) :: vec
  !   integer(ip),       intent(inout) :: l2g(:)
  !   ! Locals
  !   integer(ip) :: nd,nl,ne,i,inode,j

  !   ! Asserts
  !   nd=size(ev%data,dim=1)
  !   assert(nd==vec%nd)
  !   assert(nd == 1)
  !   assert(size(nn,1)==el%nint) 

  !   l2g = 0
  !   do j = dofh%i_varsxprob(1)%a(el%prob),dofh%i_varsxprob(1)%a(el%prob+1)-1
  !      do inode = 1,nn(el%iv(dofh%j_varsxprob(1)%a(j)))
  !         i = el%p_nod(inode) + j - dofh%i_varsxprob(1)%a(el%prob)
  !         l2g(i) = el%elem2dof(inode,dofh%j_varsxprob(1)%a(j))
  !      end do
  !   end do

  !   call ass_vec_w_dof_handler ( ev%data, el%ldof, l2g, vec%neq, vec%b, dofh )

  ! end subroutine fem_vector_assembly_w_dof_handler

  !=============================================================================
 ! subroutine fem_vector_assembly_nosq_w_dof_handler(nn, dofh, el, l2g, bl, ev, vec)
 !    implicit none
 !    integer(ip) ,      intent(in)    :: nn(:)
 !    type(dof_handler), intent(in)    :: dofh
 !    type(fem_element), intent(in)    :: el
 !    type(fem_blocks),  intent(in)    :: bl
 !    type(elvec) ,      intent(in)    :: ev
 !    type(fem_vector),  intent(inout) :: vec
 !    integer(ip),       intent(inout) :: l2g(:)
 !    ! Locals
 !    integer(ip) :: nd,nl,ne,i,inode,j,maxnode


 !  end subroutine fem_vector_assembly_nosq_w_dof_handler

 ! !=============================================================================
 !  subroutine fem_face_vector_assembly_w_dof_handler(nn, l2g, dofh, el1, el2, ev, vec)
 !    implicit none
 !    integer(ip) ,      intent(in)        :: nn(2)
 !    type(dof_handler), intent(in)        :: dofh
 !    type(fem_element), intent(in)        :: el1,el2
 !    type(elvec) ,      intent(in)        :: ev
 !    integer(ip),       intent(inout)     :: l2g((nn(1)+nn(2))*dofh%nvarsxprob(el1%prob))
 !    type(fem_vector),  intent(inout)     :: vec
 !    integer(ip)                          :: i,inode,j


 !  end subroutine fem_face_vector_assembly_w_dof_handler

  ! !=============================================================================
  ! subroutine ass_vec_w_dof_handler (ev,nl,l2g,nv,b, dofhandler)
  !   implicit none
  !   integer(ip) , intent(in)    :: nv,nl
  !   integer(ip) , intent(in)    :: l2g(nl)
  !   real(rp)    , intent(in)    :: ev(1,nl)
  !   real(rp)    , intent(inout) :: b(1,nv)
  !   type(dof_handler), intent(in), optional :: dofhandler
    
  !   integer(ip)                 :: il,ig
  !   integer(ip)                 :: block_i, l_cing_i, cing_node, num_dofs, cing_glob_dof_i
  !   type(adapt_constraint), pointer :: constraint_ptr

  !   ! todo: find it out properly!
  !   ! todo: find it out properly!
  !   block_i = 1
        
  !   do il = 1,nl
  !      ig = l2g(il)
  !      if (ig /= 0) then
  !         if(ig > 0) then
  !            ! regular node
  !            !write(*,*) "adding regular rhs ", ig, il, ev(1,il)
  !            b(1, ig) = b(1, ig) + ev(1,il)
  !         else
  !            ! hanging node
  !            assert(present(dofhandler))
  !            constraint_ptr => dofhandler%constraint_list%list(-ig)
  !            do l_cing_i = 1, constraint_ptr%num_cing_nodes
  !               cing_node = constraint_ptr%cing_nodes(l_cing_i)
  !               num_dofs = dofhandler%ob2dof_p(block_i)%a(cing_node + 1) - dofhandler%ob2dof_p(block_i)%a(cing_node)
  !               assert(num_dofs <= 1)
  !               if(num_dofs == 0) then
  !                  cycle
  !               end if
  !               ! todo: second column in ob2dof_l
  !               cing_glob_dof_i = dofhandler%ob2dof_l(block_i)%a(dofhandler%ob2dof_p(block_i)%a(cing_node), 1)

  !               ! add the value
  !               b(1, cing_glob_dof_i) = b(1, cing_glob_dof_i) + ev(1, il) * constraint_ptr%cing_coefs(l_cing_i)
  !            end do
  !         end if
  !      end if
  !   end do
  ! end subroutine ass_vec_w_dof_handler

 ! !============================================================================
 !  subroutine ass_vec_nosq_w_dof_handler(nint,nn,nd,id,ld,ib,jb,nva,iv,pn,l2g,ev,nv,mn,jbn,b)

 !    implicit none
 !    integer(ip) , intent(in)      :: nint, nv, nd, nva, id, ld, mn
 !    integer(ip) , intent(in)      :: nn(nint),pn(nd),iv(nva)
 !    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
 !    integer(ip) , intent(in)      :: l2g(id)
 !    integer(ip) , intent(in)      :: jbn(mn,nva)
 !    real(rp)    , intent(in)      :: ev(1,ld)
 !    real(rp)    , intent(inout)   :: b(1,nv)

 !    ! local variables
 !    integer(ip)                   :: in, ip, il, ig, ie

 !    il = 0
 !    do ip = ib(1),ib(2)-1
 !       do in = 1,nn(iv(jb(ip)))
 !          il = il + 1
 !          ig = l2g(il)
 !          if (ig /= 0) then
 !             ie = pn(in)-1 + jbn(in,jb(ip))
 !             b(1,ig) = b(1,ig) + ev(1,ie)
 !          end if
 !       end do
 !    end do

 !  end subroutine ass_vec_nosq_w_dof_handler

  !=============================================================================
  subroutine fem_vector_dot (x, y, t)
    implicit none
    type(fem_vector), intent(in)  :: x
    type(fem_vector), intent(in)  :: y
    real(rp)        , intent(out) :: t

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    t = ddot( x%neq, x%b, 1, y%b, 1 )
#else
!!$    AFM: A non BLAS-based implementation of the
!!$    dot product should go here
    check(1==0)
#endif

  end subroutine fem_vector_dot
  !=============================================================================
  subroutine fem_vector_nrm2(x,t)
    implicit none
    type(fem_vector), intent(in)  :: x
    real(rp)    , intent(out)     :: t

#ifdef ENABLE_BLAS
    t = dnrm2( x%neq, x%b, 1 )
#else
    call fem_vector_dot (x, x, t)
    t = sqrt(t)
#endif
  end subroutine fem_vector_nrm2
  !=============================================================================
  subroutine fem_vector_copy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y

    assert ( x%neq == y%neq )

#ifdef ENABLE_BLAS
    call dcopy ( x%neq, x%b, 1, y%b, 1 ) 
#else
    y%b=x%b
#endif
  end subroutine fem_vector_copy
  subroutine fem_vector_zero(y)
    implicit none
    type(fem_vector), intent(inout) :: y
    y%b=0.0_rp
  end subroutine fem_vector_zero

  subroutine fem_vector_init(alpha, y)
    implicit none
    type(fem_vector), intent(inout) :: y 
    real(rp), intent(in)            :: alpha  
    y%b=alpha
  end subroutine fem_vector_init
  
  subroutine fem_vector_scale(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y

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
  end subroutine fem_vector_scale

  subroutine fem_vector_mxpy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, -1.0, x%b, 1, y%b, 1 )
#else
    y%b=y%b-x%b
#endif
  end subroutine fem_vector_mxpy
  subroutine fem_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, t, x%b, 1, y%b, 1 )
#else
    y%b=y%b+t*x%b
#endif
  end subroutine fem_vector_axpy

  subroutine fem_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
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
  end subroutine fem_vector_aypx

  subroutine fem_vector_pxpy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%neq == y%neq )
#ifdef ENABLE_BLAS
    call daxpy ( x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=y%b+x%b
#endif
  end subroutine fem_vector_pxpy

  subroutine fem_vector_pxmy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
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
  end subroutine fem_vector_pxmy

  subroutine fem_vector_print (luout, x)
    implicit none
    type(fem_vector), intent(in) :: x
    integer(ip),      intent(in) :: luout
    
    ! Locals
    integer(ip) :: i

    write (luout, '(a)')     '*** begin fem_vector data structure ***'
    write(luout,'(a,i10)') 'size', x%neq
       write (luout,'(e25.16)') x%b
!!$    do i=1,x%neq
!!$!       write (luout,'(e14.7)') x%b(1,i)
!!$       write (luout,'(e25.16,1x)') x%b(1,i)
!!$    end do

  end subroutine fem_vector_print

  subroutine fem_vector_print_matrix_market ( luout, x )
   implicit none
   ! Parameters
   type(fem_vector), intent(in) :: x
   integer(ip),      intent(in) :: luout

   ! Locals
   integer (ip) :: i, id

   write (luout,'(a)') '%%MatrixMarket matrix array real general'
   write (luout,*) x%neq , 1
   do i=1,x%neq 
            write (luout,*) x%b( i )
   end do

 end subroutine fem_vector_print_matrix_market

 ! alpha <- op1^T * op2
 function fem_vector_dot_tbp(op1,op2) result(alpha)
   implicit none
   class(fem_vector), intent(in)    :: op1
   class(base_operand), intent(in)  :: op2
   real(rp) :: alpha

   call op1%GuardTemp()
   call op2%GuardTemp()
   select type(op2)
   class is (fem_vector)
      assert ( op1%neq == op2%neq )
#ifdef ENABLE_BLAS
      alpha = ddot( op1%neq, op1%b, 1, op2%b, 1 )
#else
    check(1==0)
#endif
   class default
      write(0,'(a)') 'fem_vector%dot: unsupported op2 class'
      check(1==0)
   end select
   call op1%CleanTemp()
   call op2%CleanTemp()
 end function fem_vector_dot_tbp

 ! op1 <- op2 
 subroutine fem_vector_copy_tbp(op1,op2)
   implicit none
   class(fem_vector), intent(inout) :: op1
   class(base_operand), intent(in)  :: op2
   
   call op2%GuardTemp()
   select type(op2)
   class is (fem_vector)
      assert ( op2%neq == op1%neq )
#ifdef ENABLE_BLAS
      call dcopy ( op2%neq, op2%b, 1, op1%b, 1 ) 
#else
      op1%b=op2%b
#endif
   class default
      write(0,'(a)') 'fem_vector%copy: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine fem_vector_copy_tbp

 ! op1 <- alpha * op2
 subroutine fem_vector_scal_tbp(op1,alpha,op2)
   implicit none
   class(fem_vector), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand), intent(in) :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (fem_vector)
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
      write(0,'(a)') 'fem_vector%scal: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine fem_vector_scal_tbp
 ! op <- alpha
 subroutine fem_vector_init_tbp(op,alpha)
   implicit none
   class(fem_vector), intent(inout) :: op
   real(rp), intent(in) :: alpha
   op%b=alpha
 end subroutine fem_vector_init_tbp

 ! op1 <- alpha*op2 + beta*op1
 subroutine fem_vector_axpby_tbp(op1, alpha, op2, beta)
   implicit none
   class(fem_vector), intent(inout) :: op1
   real(rp), intent(in) :: alpha
   class(base_operand), intent(in) :: op2
   real(rp), intent(in) :: beta

   call op2%GuardTemp()
   select type(op2)
   class is (fem_vector)
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
      write(0,'(a)') 'fem_vector%axpby: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine fem_vector_axpby_tbp

 ! alpha <- nrm2(op)
 function fem_vector_nrm2_tbp(op) result(alpha)
   implicit none
   class(fem_vector), intent(in)  :: op
   real(rp) :: alpha
   call op%GuardTemp()

#ifdef ENABLE_BLAS
    alpha = dnrm2( op%neq, op%b, 1 )
#else
    alpha = op%dot(op)
    alpha = sqrt(alpha)
#endif

   call op%CleanTemp()
 end function fem_vector_nrm2_tbp

 ! op1 <- clone(op2) 
 subroutine fem_vector_clone_tbp(op1,op2)
   implicit none
   class(fem_vector), intent(inout) :: op1
   class(base_operand), intent(in)  :: op2

   call op2%GuardTemp()
   select type(op2)
   class is (fem_vector)
      if (op1%mode == allocated) call memfreep(op1%b,__FILE__,__LINE__)
      op1%neq     =  op2%neq       ! Number of equations
         call memallocp(op1%neq,op1%b,__FILE__,__LINE__)
      ! AFM: I think that clone should NOT init the memory just allocated.
      ! The code that surrounds clone (e.g., Krylov solvers) should not
      ! rely on fem_vector_clone_tbp initializing the memory. I will comment
      ! it out, and expect that the codes continue working.
      ! op1%b = 0.0_rp 
      op1%mode = allocated
   class default
      write(0,'(a)') 'fem_vector%clone: unsupported op2 class'
      check(1==0)
   end select
   call op2%CleanTemp()
 end subroutine fem_vector_clone_tbp

 ! op <- comm(op)
 subroutine fem_vector_comm_tbp(op)
   implicit none
   class(fem_vector), intent(inout) :: op
 end subroutine fem_vector_comm_tbp

 subroutine fem_vector_free_tbp(this)
   implicit none
   class(fem_vector), intent(inout) :: this

   this%neq     = 0          ! Number of equations
   if (this%mode == allocated) call memfreep(this%b,__FILE__,__LINE__)
   nullify(this%b)
   this%mode = not_created
 end subroutine fem_vector_free_tbp

end module fem_vector_names
