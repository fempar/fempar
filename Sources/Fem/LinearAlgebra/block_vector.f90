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

  ! vector
  type block_vector_t
     integer(ip)                   :: nblocks = 0
     type(vector_t), allocatable :: blocks(:)
  end type block_vector_t

  ! Types
  public :: block_vector_t

  ! Functions
  public :: block_vector_free, block_vector_alloc,        & 
            block_vector_create_view, block_vector_clone, & 
            block_vector_comm,                                &
            block_vector_dot,                                 & 
            block_vector_nrm2, block_vector_copy,         & 
            block_vector_zero, block_vector_init,         & 
            block_vector_scale, block_vector_mxpy,        & 
            block_vector_axpy, block_vector_aypx,         & 
            block_vector_pxpy, block_vector_pxmy,         & 
            block_vector_print
contains

  !=============================================================================
  subroutine block_vector_free (bvec)
    implicit none
    type(block_vector_t), intent(inout) :: bvec
    integer(ip)  :: ib
   
    do ib=1, bvec%nblocks
       call vector_free ( bvec%blocks(ib) )
    end do
    
    bvec%nblocks = 0
    deallocate( bvec%blocks )
  end subroutine block_vector_free

  !=============================================================================
  subroutine block_vector_alloc(nblocks, bvec)
    implicit none
    integer(ip)           , intent(in)  :: nblocks
    type(block_vector_t), intent(out) :: bvec
    bvec%nblocks = nblocks
    allocate ( bvec%blocks(nblocks) )
  end subroutine block_vector_alloc

  !=============================================================================
  subroutine block_vector_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(block_vector_t), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(block_vector_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call block_vector_alloc ( svec%nblocks, tvec )

    do ib=1, svec%nblocks
       call vector_create_view (svec%blocks(ib), start, end, tvec%blocks(ib))
    end do
  end subroutine block_vector_create_view

  !=============================================================================
  subroutine block_vector_clone ( svec, tvec )
    implicit none
    ! Parameters
    type(block_vector_t), intent( in ) :: svec
    type(block_vector_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib
   
    call block_vector_alloc ( svec%nblocks, tvec )
   
    do ib=1, svec%nblocks
       call vector_clone (svec%blocks(ib), tvec%blocks(ib))
    end do

  end subroutine block_vector_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine block_vector_comm ( vec )
    implicit none
    type(block_vector_t), intent( inout ) :: vec 
  end subroutine block_vector_comm

  !=============================================================================
  subroutine block_vector_dot (x, y, t)
    implicit none
    ! Parameters
    type(block_vector_t), intent(in)  :: x
    type(block_vector_t), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      call vector_dot ( x%blocks(ib), y%blocks(ib), aux )
      t = t + aux
    end do 
  end subroutine block_vector_dot
  !=============================================================================
  subroutine block_vector_nrm2(x,t)
    implicit none
    type(block_vector_t), intent(in)  :: x
    real(rp)    , intent(out)     :: t

    call block_vector_dot (x, x, t)
    t = sqrt(t)
  end subroutine block_vector_nrm2
  !=============================================================================
  subroutine block_vector_copy(x,y)
    implicit none
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call vector_copy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine block_vector_copy

  subroutine block_vector_zero(y)
    implicit none
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call vector_zero ( y%blocks(ib) )
    end do 

  end subroutine block_vector_zero

  subroutine block_vector_init(alpha, y)
    implicit none
    type(block_vector_t), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call vector_init ( alpha, y%blocks(ib) )
    end do    
  end subroutine block_vector_init
  
  subroutine block_vector_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_scale ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine block_vector_scale

  subroutine block_vector_mxpy(x,y)
    implicit none
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_mxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine block_vector_mxpy
  subroutine block_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_axpy ( t, x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine block_vector_axpy

  subroutine block_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_aypx ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine block_vector_aypx

  subroutine block_vector_pxpy(x,y)
    implicit none
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_pxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine block_vector_pxpy

  subroutine block_vector_pxmy(x,y)
    implicit none
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call vector_pxmy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine block_vector_pxmy

  subroutine block_vector_print (luout, x)
    implicit none
    type(block_vector_t), intent(in) :: x
    integer(ip)           , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, x%nblocks
      call vector_print ( luout, x%blocks(ib) )
    end do 
  end subroutine block_vector_print

  subroutine ass_blkvec_w_dof_descriptor(nint,nn,nd,id,ld,ib,jb,nva,iv,pn,l2g,ev,nv,mn,jbn,b)
    implicit none
    integer(ip) , intent(in)      :: nint, nv, nd, nva, id, ld, mn
    integer(ip) , intent(in)      :: nn(nint),pn(nd),iv(nva)
    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
    integer(ip) , intent(in)      :: l2g(id)
    integer(ip) , intent(in)      :: jbn(mn,nva)
    real(rp)    , intent(in)      :: ev(1,ld)
    real(rp)    , intent(inout)   :: b(1,nv)

    ! local variables
    integer(ip)                   :: in, ipp, il, ig, ie

    il = 0
    do ipp = ib(1),ib(2)-1
       do in = 1,nn(iv(jb(ipp)))
          il = il + 1
          ig = l2g(il)
          if (ig /= 0) then
             ie = pn(in)-1 + jbn(in,jb(ipp))
             b(1,ig) = b(1,ig) + ev(1,ie)
          end if
       end do
    end do

  end subroutine ass_blkvec_w_dof_descriptor

end module block_vector_names
