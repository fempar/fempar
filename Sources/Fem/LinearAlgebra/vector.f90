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
module fem_vector_class
  use types
  use memor
  use vector_dof
  use elvec_class
  use blas77_interfaces
  use dof_handler_class
  use fem_space_class
  use fem_blocks_class
  use adaptivity
!!$#ifdef memcheck
!!$  use iso_c_binding
!!$#endif
  implicit none
# include "debug.i90"

  ! *** IMPORTANT NOTE: This cpp macro should go to a 
  ! common include file or it should be a program 
  ! subroutine otherwise
# define blk2scal(iv,idof,ndof) (((iv)-1)*(ndof)+(idof))

  !=============================================================
  ! TODO:
  ! 
  ! x Call to BLAS double precision or single precision 
  !   subroutines depending on the value of the rp parameter. 
  !   Currently we are always calling double precision variants 
  !   of the BLAS subroutines.
  ! 
  !=============================================================

  integer(ip), parameter :: allocated  = 0 ! fem_vector%b has been allocated
  integer(ip), parameter :: reference  = 1 ! fem_vector%b references external memory

  private

  ! fem_vector
  type fem_vector
     integer(ip)                :: &
        nd  = 0,                   &  ! Number of degrees of freedom, ndof1
        neq = 0                       ! Number of equations
   
     integer(ip)                :: & 
        mode = reference              ! Creation mode (by default references to NULL)

     integer(ip)                :: &  ! Storage layout (blk: block; scal: scalar)
        storage = undef_sto

     real(rp), pointer          :: &
        b(:,:) => NULL()
  end type fem_vector

  interface fem_vector_assembly
     module procedure fem_vector_assembly_w_dof_handler
     module procedure fem_vector_assembly_nosq_w_dof_handler
  end interface fem_vector_assembly

  interface fem_vector_create_view
     module procedure fem_vector_create_view_same_ndof, & 
                      fem_vector_create_view_new_ndof
  end interface fem_vector_create_view

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
            fem_vector_alloc, fem_vector_create_view, fem_vector_clone,            & 
            fem_vector_assembly, fem_face_vector_assembly_w_dof_handler,           &
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
    vec%nd      = 0          ! Number of degrees of freedom
    vec%neq     = 0          ! Number of equations
    vec%storage = undef_sto
    if (vec%mode == allocated) call memfreep(vec%b,__FILE__,__LINE__)
  end subroutine fem_vector_free

  !=============================================================================
  subroutine fem_vector_alloc(storage,nd,neq,vec)
    implicit none
    integer(ip)     , intent(in)  :: storage, nd, neq
    type(fem_vector), intent(out) :: vec
    assert ( storage == blk .or. storage == scal )
    vec%nd      = nd   ! Number of degrees of freedom
    vec%neq     = neq  ! Number of equations
    vec%storage = storage
    if ( vec%storage == blk ) then
      call memallocp(vec%nd,vec%neq,vec%b,__FILE__,__LINE__)
    else if ( vec%storage == scal ) then
      call memallocp(1,vec%nd*vec%neq,vec%b,__FILE__,__LINE__)
    end if
    vec%b    = 0.0_rp
    vec%mode = allocated
  end subroutine fem_vector_alloc
  
  subroutine fem_vector_create_view_same_ndof (svec, start, end, tvec)
    implicit none
    type(fem_vector), intent(in), target  :: svec
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(fem_vector), intent(out)         :: tvec

    tvec%nd      =  svec%nd               ! Number of degrees of freedom
    tvec%neq     =  end-start+1           ! Number of equations
    tvec%storage =  svec%storage
         
    if ( tvec%storage == blk ) then
      tvec%b => svec%b(:, start:end)
    else if ( tvec%storage == scal ) then
      tvec%b => svec%b(:, blk2scal(start,1,tvec%nd):blk2scal(end,tvec%nd,tvec%nd))
    end if
    
    tvec%mode =  reference
  end subroutine fem_vector_create_view_same_ndof

  ! This subroutine is required in order to create a view which has a 
  ! different ndof of an existing fem_vector with a given ndof
  subroutine fem_vector_create_view_new_ndof (svec, ndstart, ndend, start, end, tvec)
    implicit none
    ! Parameters
    type(fem_vector), intent(in), target  :: svec
    integer(ip)     , intent(in)          :: ndstart
    integer(ip)     , intent(in)          :: ndend
    integer(ip)     , intent(in)          :: start
    integer(ip)     , intent(in)          :: end
    type(fem_vector), intent(out)         :: tvec

    ! Locals
    integer(ip) :: off, start_aux, end_aux

    tvec%nd      =  ndend-ndstart+1       ! Number of degrees of freedom
    tvec%neq     =  end-start+1           ! Number of equations
    tvec%storage =  svec%storage
         
    start_aux = mod(start-1, svec%neq)+1
    end_aux   = mod(end-1  , svec%neq)+1

    if ( tvec%storage == blk ) then
      tvec%b => svec%b(ndstart:ndend, start:end)
    else if ( tvec%storage == scal ) then

      ! Compute initial offset
      off = (ndstart-1) * svec%neq
      tvec%b => svec%b(:, off+blk2scal(start_aux,1,tvec%nd):off+blk2scal(end_aux,tvec%nd,tvec%nd))
    end if
    
    tvec%mode =  reference
  end subroutine fem_vector_create_view_new_ndof



  subroutine fem_vector_clone ( svec, tvec )
    implicit none
    type(fem_vector), intent( in ) :: svec
    type(fem_vector), intent(out) :: tvec
    tvec%nd      =  svec%nd        ! Number of degrees of freedom
    tvec%neq     =  svec%neq       ! Number of equations
    tvec%storage =  svec%storage
    if ( tvec%storage == blk ) then
      call memallocp(tvec%nd,tvec%neq,tvec%b,__FILE__,__LINE__)
    else if ( tvec%storage == scal ) then
      call memallocp(1,tvec%nd*tvec%neq,tvec%b,__FILE__,__LINE__)
    end if
    tvec%b = 0.0_rp
    tvec%mode = allocated 
  end subroutine fem_vector_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine fem_vector_comm ( vec )
    implicit none
    type(fem_vector), intent( inout ) :: vec 
  end subroutine fem_vector_comm

  !=============================================================================
  subroutine fem_vector_assembly_w_dof_handler(nn, dofh, el, l2g, ev, vec)
    implicit none
    type(dof_handler), intent(in)    :: dofh
    type(fem_element), intent(in)    :: el
    type(elvec) ,      intent(in)    :: ev
    integer(ip) ,      intent(in)    :: nn(:)
    type(fem_vector),  intent(inout) :: vec
    integer(ip),       intent(inout) :: l2g(:)
    ! Locals
    integer(ip) :: nd,nl,ne,i,inode,j

    ! Asserts
    nd=size(ev%data,dim=1)
    assert(nd==vec%nd)
    assert(nd == 1)
    assert(size(nn,1)==el%nint) 

    l2g = 0
    do j = dofh%i_varsxprob(1)%a(el%prob),dofh%i_varsxprob(1)%a(el%prob+1)-1
       do inode = 1,nn(el%iv(dofh%j_varsxprob(1)%a(j)))
          i = el%p_nod(inode) + j - dofh%i_varsxprob(1)%a(el%prob)
          l2g(i) = el%elem2dof(inode,dofh%j_varsxprob(1)%a(j))
       end do
    end do

    call ass_vec_w_dof_handler ( ev%data, el%ldof, l2g, vec%neq, vec%b, dofh )

  end subroutine fem_vector_assembly_w_dof_handler

  !=============================================================================
  subroutine fem_vector_assembly_nosq_w_dof_handler(nn, dofh, el, l2g, bl, ev, vec)
    implicit none
    integer(ip) ,      intent(in)    :: nn(:)
    type(dof_handler), intent(in)    :: dofh
    type(fem_element), intent(in)    :: el
    type(fem_blocks),  intent(in)    :: bl
    type(elvec) ,      intent(in)    :: ev
    type(fem_vector),  intent(inout) :: vec
    integer(ip),       intent(inout) :: l2g(:)
    ! Locals
    integer(ip) :: nd,nl,ne,i,inode,j,maxnode

    ! Asserts
    nd=size(ev%data,dim=1)
    assert(nd==vec%nd)
    assert(nd == 1)
    assert(size(nn,1)==el%nint) 

    maxnode=size(el%jvars,1)       

    i = 0
    do j =  bl%ib(1),bl%ib(2)-1
       do inode = 1,nn(el%iv(bl%jb(j)))
          i = i+1
          l2g(i) = el%elem2dof(inode,bl%jb(j))
       end do
    end do

    call ass_vec_nosq_w_dof_handler(el%nint,nn,size(el%p_nod,1),i,el%ldof,bl%ib,bl%jb,size(el%iv,1), &
         &                          el%iv,el%p_nod,l2g,ev%data,vec%neq,maxnode,el%jvars,vec%b)

  end subroutine fem_vector_assembly_nosq_w_dof_handler

 !=============================================================================
  subroutine fem_face_vector_assembly_w_dof_handler(nn, l2g, dofh, el1, el2, ev, vec)
    implicit none
    integer(ip) ,      intent(in)        :: nn(2)
    type(dof_handler), intent(in)        :: dofh
    type(fem_element), intent(in)        :: el1,el2
    type(elvec) ,      intent(in)        :: ev
    integer(ip),       intent(inout)     :: l2g((nn(1)+nn(2))*dofh%nvarsxprob(el1%prob))
    type(fem_vector),  intent(inout)     :: vec
    integer(ip)                          :: i,inode,j

    assert(size(ev%data,dim=1) == 1); assert(vec%nd == 1); assert(el1%prob == el2%prob)

    ! Construct l2g
    i = 0
    do inode = 1,nn(1)
       do j = 1,dofh%nvarsxprob(el1%prob)
          i = i+1
          l2g(i) = el1%elem2dof(inode,j)
       end do
    end do
    do inode = 1,nn(2)
       do j = 1,dofh%nvarsxprob(el2%prob)
          i = i+1
          l2g(i) = el2%elem2dof(inode,j)
       end do
    end do

    ! Assemble
    call ass_vec_w_dof_handler ( ev%data, (nn(1)+nn(2))*dofh%nvarsxprob(el1%prob), l2g, vec%neq, vec%b, dofh )

  end subroutine fem_face_vector_assembly_w_dof_handler

  !=============================================================================
  subroutine ass_vec_w_dof_handler (ev,nl,l2g,nv,b, dofhandler)
    implicit none
    integer(ip) , intent(in)    :: nv,nl
    integer(ip) , intent(in)    :: l2g(nl)
    real(rp)    , intent(in)    :: ev(1,nl)
    real(rp)    , intent(inout) :: b(1,nv)
    type(dof_handler), intent(in), optional :: dofhandler
    
    integer(ip)                 :: il,ig
    integer(ip)                 :: block_i, l_cing_i, cing_node, num_dofs, cing_glob_dof_i
    type(adapt_constraint), pointer :: constraint_ptr

    ! todo: find it out properly!
    ! todo: find it out properly!
    block_i = 1
        
    do il = 1,nl
       ig = l2g(il)
       if (ig /= 0) then
          if(ig > 0) then
             ! regular node
             !write(*,*) "adding regular rhs ", ig, il, ev(1,il)
             b(1, ig) = b(1, ig) + ev(1,il)
          else
             ! hanging node
             assert(present(dofhandler))
             constraint_ptr => dofhandler%constraint_list%list(-ig)
             do l_cing_i = 1, constraint_ptr%num_cing_nodes
                cing_node = constraint_ptr%cing_nodes(l_cing_i)
                num_dofs = dofhandler%ob2dof_p(block_i)%a(cing_node + 1) - dofhandler%ob2dof_p(block_i)%a(cing_node)
                assert(num_dofs <= 1)
                if(num_dofs == 0) then
                   cycle
                end if
                ! todo: second column in ob2dof_l
                cing_glob_dof_i = dofhandler%ob2dof_l(block_i)%a(dofhandler%ob2dof_p(block_i)%a(cing_node), 1)

                ! add the value
                b(1, cing_glob_dof_i) = b(1, cing_glob_dof_i) + ev(1, il) * constraint_ptr%cing_coefs(l_cing_i)
             end do
          end if
       end if
    end do
  end subroutine ass_vec_w_dof_handler

 !============================================================================
  subroutine ass_vec_nosq_w_dof_handler(nint,nn,nd,id,ld,ib,jb,nva,iv,pn,l2g,ev,nv,mn,jbn,b)

    implicit none
    integer(ip) , intent(in)      :: nint, nv, nd, nva, id, ld, mn
    integer(ip) , intent(in)      :: nn(nint),pn(nd),iv(nva)
    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
    integer(ip) , intent(in)      :: l2g(id)
    integer(ip) , intent(in)      :: jbn(mn,nva)
    real(rp)    , intent(in)      :: ev(1,ld)
    real(rp)    , intent(inout)   :: b(1,nv)

    ! local variables
    integer(ip)                   :: in, ip, il, ig, ie

    il = 0
    do ip = ib(1),ib(2)-1
       do in = 1,nn(iv(jb(ip)))
          il = il + 1
          ig = l2g(il)
          if (ig /= 0) then
             ie = pn(in)-1 + jbn(in,jb(ip))
             b(1,ig) = b(1,ig) + ev(1,ie)
          end if
       end do
    end do

  end subroutine ass_vec_nosq_w_dof_handler

  !=============================================================================
  subroutine fem_vector_dot (x, y, t)
    implicit none
    type(fem_vector), intent(in)  :: x
    type(fem_vector), intent(in)  :: y
    real(rp)        , intent(out) :: t

    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )

#ifdef ENABLE_BLAS
    t = ddot( x%nd * x%neq, x%b, 1, y%b, 1 )
#else
    if ( x%storage == blk ) then
       call dot_vec( x%nd, x%neq, x%b, y%b, t )
    else if ( x%storage == scal ) then
       call dot_vec_scal( x%nd, x%neq, x%b, y%b, t )
    end if 
#endif

  end subroutine fem_vector_dot
  !=============================================================================
  subroutine fem_vector_nrm2(x,t)
    implicit none
    type(fem_vector), intent(in)  :: x
    real(rp)    , intent(out)     :: t

#ifdef ENABLE_BLAS
    t = dnrm2( x%nd*x%neq, x%b, 1 )
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

    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )

#ifdef ENABLE_BLAS
    call dcopy ( x%nd*x%neq, x%b, 1, y%b, 1 ) 
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

    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )

#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dcopy ( x%nd*x%neq, x%b, 1, y%b, 1)
    call dscal ( y%nd*y%neq, t, y%b, 1)
#else
    y%b=t*x%b
#endif
  end subroutine fem_vector_scale

  subroutine fem_vector_mxpy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )
#ifdef ENABLE_BLAS
    call daxpy ( x%nd*x%neq, -1.0, x%b, 1, y%b, 1 )
#else
    y%b=y%b-x%b
#endif
  end subroutine fem_vector_mxpy
  subroutine fem_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )
#ifdef ENABLE_BLAS
    call daxpy ( x%nd*x%neq, t, x%b, 1, y%b, 1 )
#else
    y%b=y%b+t*x%b
#endif
  end subroutine fem_vector_axpy

  subroutine fem_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )
#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dscal ( y%nd*y%neq, t, y%b, 1)
    call daxpy ( x%nd*x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=x%b+t*y%b
#endif
  end subroutine fem_vector_aypx

  subroutine fem_vector_pxpy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )
#ifdef ENABLE_BLAS
    call daxpy ( x%nd*x%neq, 1.0, x%b, 1, y%b, 1 )    
#else
    y%b=y%b+x%b
#endif
  end subroutine fem_vector_pxpy

  subroutine fem_vector_pxmy(x,y)
    implicit none
    type(fem_vector), intent(in)    :: x
    type(fem_vector), intent(inout) :: y
    assert ( x%nd  == y%nd  )
    assert ( x%neq == y%neq )
    assert ( x%storage == y%storage )
#ifdef ENABLE_BLAS
    ! I guess that two calls to the level 1
    ! BLAS can not be competitive against 
    ! just one F90 vector operation. I have to
    ! measure the difference among these two
    ! options. 
    call dscal ( y%nd*y%neq, -1.0, y%b, 1)
    call daxpy ( x%nd*x%neq, 1.0, x%b, 1, y%b, 1 ) 
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
    write(luout,'(a,i10)') 'size', x%neq, 'ndof', x%nd
    if ( x%storage == blk ) then
       write (luout,'(e15.7)') x%b
    else if ( x%storage == scal ) then
       write (luout,'(e25.16)') x%b
    end if
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
   write (luout,*) x%neq * x%nd, 1
   do i=1,x%neq 
      do id=1,x%nd
         if ( x%storage == blk ) then
            write (luout,*) x%b(id,i)
         else if ( x%storage == scal ) then
            write (luout,*) x%b(1, blk2scal(i, id, x%nd) )
         end if
      end do
   end do

 end subroutine fem_vector_print_matrix_market

end module fem_vector_class
