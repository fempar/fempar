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
  use abstract_vector_names
  use block_graph_names
  use graph_names
  implicit none
# include "debug.i90"

  private

  ! vector
  type, extends(abstract_vector_t) :: serial_block_array_t
     integer(ip)                 :: nblocks = 0
     type(serial_scalar_array_t), allocatable :: blocks(:)
   contains
     procedure :: dot  => serial_block_array_dot_tbp
     procedure :: copy => serial_block_array_copy_tbp
     procedure :: init => serial_block_array_init_tbp
     procedure :: scal => serial_block_array_scal_tbp
     procedure :: axpby => serial_block_array_axpby_tbp
     procedure :: nrm2 => serial_block_array_nrm2_tbp
     procedure :: clone => serial_block_array_clone_tbp
     procedure :: comm  => serial_block_array_comm_tbp
     procedure :: free  => serial_block_array_free_tbp
     procedure :: serial_block_array_blocks
     procedure :: serial_block_array_alloc_all
     generic :: alloc => serial_block_array_blocks, serial_block_array_alloc_all
  end type serial_block_array_t

  ! Types
  public :: serial_block_array_t

  ! Functions
  public :: serial_block_array_free,        & 
            serial_block_array_create_view, serial_block_array_clone, & 
            serial_block_array_comm,                                &
            serial_block_array_dot,                                 & 
            serial_block_array_nrm2, serial_block_array_copy,         & 
            serial_block_array_zero, serial_block_array_init,         & 
            serial_block_array_scale, serial_block_array_mxpy,        & 
            serial_block_array_axpy, serial_block_array_aypx,         & 
            serial_block_array_pxpy, serial_block_array_pxmy,         & 
            serial_block_array_print
contains

  !=============================================================================
  subroutine serial_block_array_free (bvec)
    implicit none
    type(serial_block_array_t), intent(inout) :: bvec
    integer(ip)  :: ib
   
    do ib=1, bvec%nblocks
       call serial_scalar_array_free ( bvec%blocks(ib) )
    end do
    
    bvec%nblocks = 0
    deallocate( bvec%blocks )
  end subroutine serial_block_array_free

  !=============================================================================
  subroutine serial_block_array_alloc_all(bvec, bgraph)
    implicit none
    class(serial_block_array_t), intent(out) :: bvec
    type(block_graph_t) , intent(in)  :: bgraph
    integer(ip)  :: ib
    type(graph_t), pointer :: f_graph
    
    bvec%nblocks = bgraph%get_nblocks()
    allocate ( bvec%blocks(bvec%nblocks) )
    do ib=1, bvec%nblocks
       f_graph => bgraph%get_block(ib,ib)
       call serial_scalar_array_alloc ( f_graph%nv, bvec%blocks(ib) )
    end do

  end subroutine serial_block_array_alloc_all

  !=============================================================================
  subroutine serial_block_array_blocks(bvec,nblocks)
    implicit none
    class(serial_block_array_t), intent(out) :: bvec
    integer(ip)         , intent(in)  :: nblocks
    bvec%nblocks = nblocks
    allocate ( bvec%blocks(nblocks) )
  end subroutine serial_block_array_blocks

  !=============================================================================
  subroutine serial_block_array_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(serial_block_array_t), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(serial_block_array_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call tvec%serial_block_array_blocks(svec%nblocks)

    do ib=1, svec%nblocks
       call serial_scalar_array_create_view (svec%blocks(ib), start, end, tvec%blocks(ib))
    end do
  end subroutine serial_block_array_create_view

  !=============================================================================
  subroutine serial_block_array_clone ( svec, tvec )
    implicit none
    ! Parameters
    type(serial_block_array_t), intent( in ) :: svec
    type(serial_block_array_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib
   
    call tvec%serial_block_array_blocks(svec%nblocks)
   
    do ib=1, svec%nblocks
       call serial_scalar_array_clone (svec%blocks(ib), tvec%blocks(ib))
    end do

  end subroutine serial_block_array_clone

  !=============================================================================
  ! Dummy method required to specialize Krylov subspace methods
  subroutine serial_block_array_comm ( vec )
    implicit none
    type(serial_block_array_t), intent( inout ) :: vec 
  end subroutine serial_block_array_comm

  !=============================================================================
  subroutine serial_block_array_dot (x, y, t)
    implicit none
    ! Parameters
    type(serial_block_array_t), intent(in)  :: x
    type(serial_block_array_t), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      call serial_scalar_array_dot ( x%blocks(ib), y%blocks(ib), aux )
      t = t + aux
    end do 
  end subroutine serial_block_array_dot
  !=============================================================================
  subroutine serial_block_array_nrm2(x,t)
    implicit none
    type(serial_block_array_t), intent(in)  :: x
    real(rp)    , intent(out)     :: t

    call serial_block_array_dot (x, x, t)
    t = sqrt(t)
  end subroutine serial_block_array_nrm2
  !=============================================================================
  subroutine serial_block_array_copy(x,y)
    implicit none
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call serial_scalar_array_copy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine serial_block_array_copy

  subroutine serial_block_array_zero(y)
    implicit none
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call serial_scalar_array_zero ( y%blocks(ib) )
    end do 

  end subroutine serial_block_array_zero

  subroutine serial_block_array_init(alpha, y)
    implicit none
    type(serial_block_array_t), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call serial_scalar_array_init ( alpha, y%blocks(ib) )
    end do    
  end subroutine serial_block_array_init
  
  subroutine serial_block_array_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_scale ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine serial_block_array_scale

  subroutine serial_block_array_mxpy(x,y)
    implicit none
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_mxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine serial_block_array_mxpy
  subroutine serial_block_array_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_axpy ( t, x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine serial_block_array_axpy

  subroutine serial_block_array_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_aypx ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine serial_block_array_aypx

  subroutine serial_block_array_pxpy(x,y)
    implicit none
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_pxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine serial_block_array_pxpy

  subroutine serial_block_array_pxmy(x,y)
    implicit none
    type(serial_block_array_t), intent(in)    :: x
    type(serial_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call serial_scalar_array_pxmy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine serial_block_array_pxmy

  subroutine serial_block_array_print (luout, x)
    implicit none
    type(serial_block_array_t), intent(in) :: x
    integer(ip)           , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, x%nblocks
      call serial_scalar_array_print ( luout, x%blocks(ib) )
    end do 
  end subroutine serial_block_array_print

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
  
  ! alpha <- op1^T * op2
  function serial_block_array_dot_tbp(op1,op2) result(alpha)
    implicit none
    ! Parameters
    class(serial_block_array_t), intent(in)  :: op1
    class(abstract_vector_t), intent(in)  :: op2
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
       write(0,'(a)') 'block_vector_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function serial_block_array_dot_tbp

  ! op1 <- op2 
  subroutine serial_block_array_copy_tbp(op1,op2)
    implicit none
    ! Parameters
    class(serial_block_array_t), intent(inout) :: op1
    class(abstract_vector_t), intent(in)    :: op2

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
       write(0,'(a)') 'block_vector_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_copy_tbp

  ! op <- alpha
  subroutine serial_block_array_init_tbp(op,alpha)
    implicit none
    class(serial_block_array_t), intent(inout) :: op 
    real(rp)             , intent(in)    :: alpha  
    ! Locals
    integer(ip) :: ib

    do ib=1, op%nblocks
       call op%blocks(ib)%init(alpha)
    end do    
  end subroutine serial_block_array_init_tbp
  
  ! op1 <- alpha * op2
  subroutine serial_block_array_scal_tbp(op1,alpha,op2)
    implicit none
    ! Parameters 
    class(serial_block_array_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(abstract_vector_t), intent(in)    :: op2
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
       write(0,'(a)') 'block_vector_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_scal_tbp

  ! op1 <- alpha*op2 + beta*op1
  subroutine serial_block_array_axpby_tbp(op1,alpha,op2,beta)
    implicit none
    class(serial_block_array_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(abstract_vector_t), intent(in)    :: op2
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
       write(0,'(a)') 'block_vector_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_axpby_tbp

  ! alpha <- nrm2(op)
  function serial_block_array_nrm2_tbp(op) result(alpha)
    implicit none
    class(serial_block_array_t), intent(in) :: op
    real(rp) :: alpha

    call op%GuardTemp()
    alpha = op%dot(op)
    alpha = sqrt(alpha)
   call op%CleanTemp()
  end function serial_block_array_nrm2_tbp

  ! op1 <- clone(op2) 
  subroutine serial_block_array_clone_tbp(op1,op2)
    implicit none
    ! Parameters
    class(serial_block_array_t)        , intent(inout) :: op1
    class(abstract_vector_t), target, intent(in)    :: op2
 
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
       write(0,'(a)') 'block_vector_t%clone: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine serial_block_array_clone_tbp

  ! op <- comm(op)
  subroutine serial_block_array_comm_tbp(op)
    implicit none
    class(serial_block_array_t), intent(inout) :: op 
  end subroutine serial_block_array_comm_tbp

  subroutine serial_block_array_free_tbp(this)
    implicit none
    class(serial_block_array_t), intent(inout) :: this
    integer(ip)  :: ib
   
    do ib=1, this%nblocks
       call this%blocks(ib)%free()
    end do
    this%nblocks = 0
    deallocate( this%blocks )
  end subroutine serial_block_array_free_tbp

end module serial_block_array_names
