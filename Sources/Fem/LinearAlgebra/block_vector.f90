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
module fem_block_vector_class
  use types
  use memor
  use fem_vector_class
  use block_elvec_class
  use fem_space_class
  use elvec_class
  use dof_handler_class
  use fem_blocks_class
  implicit none
# include "debug.i90"

  private

  ! fem_vector
  type fem_block_vector
     integer(ip)                   :: nblocks = 0
     type(fem_vector), allocatable :: blocks(:)
  end type fem_block_vector

  interface fem_block_vector_assembly
     module procedure fem_block_vector_assembly_w_dof_handler
     module procedure fem_block_vector_assembly_nosq_w_dof_handler
  end interface fem_block_vector_assembly

  ! Types
  public :: fem_block_vector

  ! Functions
  public :: fem_block_vector_free, fem_block_vector_alloc,        & 
            fem_block_vector_create_view, fem_block_vector_clone, & 
            fem_block_vector_comm,                                &
            fem_block_vector_assembly, fem_block_vector_dot,      & 
            fem_block_vector_nrm2, fem_block_vector_copy,         & 
            fem_block_vector_zero, fem_block_vector_init,         & 
            fem_block_vector_scale, fem_block_vector_mxpy,        & 
            fem_block_vector_axpy, fem_block_vector_aypx,         & 
            fem_block_vector_pxpy, fem_block_vector_pxmy,         & 
            fem_block_vector_print
contains

  !=============================================================================
  subroutine fem_block_vector_free (bvec)
    implicit none
    type(fem_block_vector), intent(inout) :: bvec
    integer(ip)  :: ib
   
    do ib=1, bvec%nblocks
       call fem_vector_free ( bvec%blocks(ib) )
    end do
    
    bvec%nblocks = 0
    deallocate( bvec%blocks )
  end subroutine fem_block_vector_free

  !=============================================================================
  subroutine fem_block_vector_alloc(nblocks, bvec)
    implicit none
    integer(ip)           , intent(in)  :: nblocks
    type(fem_block_vector), intent(out) :: bvec
    bvec%nblocks = nblocks
    allocate ( bvec%blocks(nblocks) )
  end subroutine fem_block_vector_alloc

  !=============================================================================
  subroutine fem_block_vector_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(fem_block_vector), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(fem_block_vector), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call fem_block_vector_alloc ( svec%nblocks, tvec )

    do ib=1, svec%nblocks
       call fem_vector_create_view (svec%blocks(ib), start, end, tvec%blocks(ib))
    end do
  end subroutine fem_block_vector_create_view

  !=============================================================================
  subroutine fem_block_vector_clone ( svec, tvec )
    implicit none
    ! Parameters
    type(fem_block_vector), intent( in ) :: svec
    type(fem_block_vector), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib
   
    call fem_block_vector_alloc ( svec%nblocks, tvec )
   
    do ib=1, svec%nblocks
       call fem_vector_clone (svec%blocks(ib), tvec%blocks(ib))
    end do

  end subroutine fem_block_vector_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine fem_block_vector_comm ( vec )
    implicit none
    type(fem_block_vector), intent( inout ) :: vec 
  end subroutine fem_block_vector_comm

  !=============================================================================
  subroutine fem_block_vector_assembly_w_dof_handler (nn, dofh, el, l2g, ev, vec)
    implicit none
    integer(ip) ,           intent(in)    :: nn(:)
    type(dof_handler),      intent(in)    :: dofh
    type(fem_element),      intent(in)    :: el
    type(elvec) ,           intent(in)    :: ev
    type(fem_block_vector), intent(inout) :: vec
    integer(ip),            intent(inout) :: l2g(:)
    ! Locals
    integer(ip)              :: nd,nl,ne,i,inode,j,ibloc,ib(3),maxnode
    integer(ip), allocatable :: jb(:)


    maxnode = size(el%jvars,1)
    
    do ibloc = 1,dofh%blocks%nb

       ! Assertsb
       nd=size(ev%data,dim=1)
       assert(nd==vec%blocks(ibloc)%nd)
       assert(nd == 1)
       assert(size(nn,1)==el%nint) ! Needed to modify to accomodate new fem info

       ! Construct block structures
       ib(1) = 1
       ib(2) = ib(1) + dofh%blknvarsxprob(ibloc)%a(el%prob) ! Rows
       ib(3) = ib(2)
       call memalloc(ib(3)-1,jb,__FILE__,__LINE__)
       jb(1:(ib(2)-1)) = dofh%j_varsxprob(ibloc)%a(dofh%i_varsxprob(ibloc)%a(el%prob) : &
            &                                     (dofh%i_varsxprob(ibloc)%a(el%prob+1)-1))
       jb(ib(2):(ib(3)-1)) = 0
       
       ! ROW DOFS
       l2g=0; i = 0
       do j = dofh%i_varsxprob(ibloc)%a(el%prob),dofh%i_varsxprob(ibloc)%a(el%prob+1)-1
          do inode = 1,nn(el%iv(j))
             i = i+1
             l2g(i) = el%elem2dof(inode,dofh%j_varsxprob(ibloc)%a(j))
          end do
       end do
       
       call ass_blkvec_w_dof_handler(el%nint,nn,size(el%p_nod,1),i,el%ldof,ib,jb,size(el%iv,1), &
            &                        el%iv,el%p_nod,l2g,ev%data,vec%blocks(ibloc)%neq,          &
            &                        maxnode,el%jvars,vec%blocks(ibloc)%b)

       call memfree(jb)

    end do

  end subroutine fem_block_vector_assembly_w_dof_handler
  
  !=============================================================================
  subroutine fem_block_vector_assembly_nosq_w_dof_handler (nn, dofh, el, l2g, bl, ev, vec)
    implicit none
    integer(ip) ,           intent(in)    :: nn(:)
    type(dof_handler),      intent(in)    :: dofh
    type(fem_element),      intent(in)    :: el
    type(fem_blocks),       intent(in)    :: bl
    type(elvec) ,           intent(in)    :: ev
    type(fem_block_vector), intent(inout) :: vec
    integer(ip),            intent(inout) :: l2g(:)

    write (*,*) 'fem_block_vector_assembly_nosq_w_dof_handler TO BE DONE'

  end subroutine fem_block_vector_assembly_nosq_w_dof_handler

  !=============================================================================
  subroutine fem_block_vector_dot (x, y, t)
    implicit none
    ! Parameters
    type(fem_block_vector), intent(in)  :: x
    type(fem_block_vector), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      call fem_vector_dot ( x%blocks(ib), y%blocks(ib), aux )
      t = t + aux
    end do 
  end subroutine fem_block_vector_dot
  !=============================================================================
  subroutine fem_block_vector_nrm2(x,t)
    implicit none
    type(fem_block_vector), intent(in)  :: x
    real(rp)    , intent(out)     :: t

    call fem_block_vector_dot (x, x, t)
    t = sqrt(t)
  end subroutine fem_block_vector_nrm2
  !=============================================================================
  subroutine fem_block_vector_copy(x,y)
    implicit none
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call fem_vector_copy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine fem_block_vector_copy

  subroutine fem_block_vector_zero(y)
    implicit none
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call fem_vector_zero ( y%blocks(ib) )
    end do 

  end subroutine fem_block_vector_zero

  subroutine fem_block_vector_init(alpha, y)
    implicit none
    type(fem_block_vector), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call fem_vector_init ( alpha, y%blocks(ib) )
    end do    
  end subroutine fem_block_vector_init
  
  subroutine fem_block_vector_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_scale ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine fem_block_vector_scale

  subroutine fem_block_vector_mxpy(x,y)
    implicit none
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_mxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine fem_block_vector_mxpy
  subroutine fem_block_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_axpy ( t, x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine fem_block_vector_axpy

  subroutine fem_block_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_aypx ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine fem_block_vector_aypx

  subroutine fem_block_vector_pxpy(x,y)
    implicit none
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_pxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine fem_block_vector_pxpy

  subroutine fem_block_vector_pxmy(x,y)
    implicit none
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call fem_vector_pxmy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine fem_block_vector_pxmy

  subroutine fem_block_vector_print (luout, x)
    implicit none
    type(fem_block_vector), intent(in) :: x
    integer(ip)           , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, x%nblocks
      call fem_vector_print ( luout, x%blocks(ib) )
    end do 
  end subroutine fem_block_vector_print

  subroutine ass_blkvec_w_dof_handler(nint,nn,nd,id,ld,ib,jb,nva,iv,pn,l2g,ev,nv,mn,jbn,b)
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

  end subroutine ass_blkvec_w_dof_handler

end module fem_block_vector_class
