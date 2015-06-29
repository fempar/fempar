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
  use base_operand_names
  use block_graph_names
  use graph_names
  implicit none
# include "debug.i90"

  private

  ! vector
  type, extends(base_operand_t) :: block_vector_t
     integer(ip)                 :: nblocks = 0
     type(vector_t), allocatable :: blocks(:)
   contains
     procedure :: dot  => block_vector_dot_tbp
     procedure :: copy => block_vector_copy_tbp
     procedure :: init => block_vector_init_tbp
     procedure :: scal => block_vector_scal_tbp
     procedure :: axpby => block_vector_axpby_tbp
     procedure :: nrm2 => block_vector_nrm2_tbp
     procedure :: clone => block_vector_clone_tbp
     procedure :: comm  => block_vector_comm_tbp
     procedure :: free  => block_vector_free_tbp
     procedure :: alloc => block_vector_alloc_all
  end type block_vector_t

  interface block_vector_alloc
     module procedure block_vector_alloc_blocks, block_vector_alloc_all
  end interface block_vector_alloc

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
  subroutine block_vector_alloc_all(bvec, bgraph)
    implicit none
    class(block_vector_t), intent(out) :: bvec
    type(block_graph_t) , intent(in)  :: bgraph
    integer(ip)  :: ib
    type(graph_t), pointer :: f_graph
    
    bvec%nblocks = bgraph%get_nblocks()
    allocate ( bvec%blocks(bvec%nblocks) )
    do ib=1, bvec%nblocks
       f_graph => bgraph%get_block(ib,ib)
       call vector_alloc ( f_graph%nv, bvec%blocks(ib) )
    end do

  end subroutine block_vector_alloc_all

  !=============================================================================
  subroutine block_vector_alloc_blocks(nblocks, bvec)
    implicit none
    integer(ip)           , intent(in)  :: nblocks
    type(block_vector_t), intent(out) :: bvec
    bvec%nblocks = nblocks
    allocate ( bvec%blocks(nblocks) )
  end subroutine block_vector_alloc_blocks

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
  
  ! alpha <- op1^T * op2
  function block_vector_dot_tbp(op1,op2) result(alpha)
    implicit none
    ! Parameters
    class(block_vector_t), intent(in)  :: op1
    class(base_operand_t), intent(in)  :: op2
    real(rp) :: alpha

    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
    class is (block_vector_t)
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
  end function block_vector_dot_tbp

  ! op1 <- op2 
  subroutine block_vector_copy_tbp(op1,op2)
    implicit none
    ! Parameters
    class(block_vector_t), intent(inout) :: op1
    class(base_operand_t), intent(in)    :: op2

    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%copy(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'block_vector_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine block_vector_copy_tbp

  ! op <- alpha
  subroutine block_vector_init_tbp(op,alpha)
    implicit none
    class(block_vector_t), intent(inout) :: op 
    real(rp)             , intent(in)    :: alpha  
    ! Locals
    integer(ip) :: ib

    do ib=1, op%nblocks
       call op%blocks(ib)%init(alpha)
    end do    
  end subroutine block_vector_init_tbp
  
  ! op1 <- alpha * op2
  subroutine block_vector_scal_tbp(op1,alpha,op2)
    implicit none
    ! Parameters 
    class(block_vector_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(base_operand_t), intent(in)    :: op2
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%scal(alpha,op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'block_vector_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine block_vector_scal_tbp

  ! op1 <- alpha*op2 + beta*op1
  subroutine block_vector_axpby_tbp(op1,alpha,op2,beta)
    implicit none
    class(block_vector_t), intent(inout) :: op1
    real(rp)             , intent(in)    :: alpha
    class(base_operand_t), intent(in)    :: op2
    real(rp)             , intent(in)    :: beta
    ! Locals
    integer(ip)                           :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%axpby(alpha,op2%blocks(ib),beta)
       end do
    class default
       write(0,'(a)') 'block_vector_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine block_vector_axpby_tbp

  ! alpha <- nrm2(op)
  function block_vector_nrm2_tbp(op) result(alpha)
    implicit none
    class(block_vector_t), intent(in) :: op
    real(rp) :: alpha

    call op%GuardTemp()
    alpha = op%dot(op)
    alpha = sqrt(alpha)
   call op%CleanTemp()
  end function block_vector_nrm2_tbp

  ! op1 <- clone(op2) 
  subroutine block_vector_clone_tbp(op1,op2)
    implicit none
    ! Parameters
    class(block_vector_t)        , intent(inout) :: op1
    class(base_operand_t), target, intent(in)    :: op2
 
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (block_vector_t)
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
  end subroutine block_vector_clone_tbp

  ! op <- comm(op)
  subroutine block_vector_comm_tbp(op)
    implicit none
    class(block_vector_t), intent(inout) :: op 
  end subroutine block_vector_comm_tbp

  subroutine block_vector_free_tbp(this)
    implicit none
    class(block_vector_t), intent(inout) :: this
    integer(ip)  :: ib
   
    do ib=1, this%nblocks
       call this%blocks(ib)%free()
    end do
    this%nblocks = 0
    deallocate( this%blocks )
  end subroutine block_vector_free_tbp


end module block_vector_names
