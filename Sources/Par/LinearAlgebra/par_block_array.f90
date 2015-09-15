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
module par_block_array_names
  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_array_names
  use abstract_vector_names
    
  ! Parallel modules
  use par_scalar_array_names
  use par_block_graph_names
  use par_graph_names

  implicit none
# include "debug.i90"

  private

  ! par_vector
  type, extends(abstract_vector_t) :: par_block_array_t
     integer(ip)                     :: nblocks = 0
     type(par_scalar_array_t), allocatable :: blocks(:)
   contains
     procedure :: dot   => par_block_array_dot_tbp
     procedure :: copy  => par_block_array_copy_tbp
     procedure :: init  => par_block_array_init_tbp
     procedure :: scal  => par_block_array_scal_tbp
     procedure :: axpby => par_block_array_axpby_tbp
     procedure :: nrm2  => par_block_array_nrm2_tbp
     procedure :: clone => par_block_array_tbp
     procedure :: comm  => par_block_array_comm_tbp
     procedure :: free  => par_block_array_free_tbp
     procedure :: par_block_array_alloc_blocks
     procedure :: par_block_array_alloc_all
     generic   :: alloc => par_block_array_alloc_blocks, par_block_array_alloc_all
  end type par_block_array_t

  ! Types
  public :: par_block_array_t

  ! Functions
  public :: par_block_array_free,                                              &
            par_block_array_create_view,       &
            par_block_array_comm,              &
            par_block_array_weight,                                            &
            par_block_array_dot,           par_block_array_nrm2,              &
            par_block_array_copy,          par_block_array_zero,              &
            par_block_array_init,          par_block_array_scale,             & 
            par_block_array_mxpy,          par_block_array_axpy,              &
            par_block_array_aypx,          par_block_array_pxpy,              &
            par_block_array_pxmy,          par_block_array_print     
contains

  !=============================================================================
  subroutine par_block_array_free (bvec)
    implicit none
    type(par_block_array_t), intent(inout) :: bvec

    ! Locals
    integer(ip) :: ib

    do ib=1, bvec%nblocks
       call bvec%blocks(ib)%free()
    end do

    bvec%nblocks = 0
    deallocate( bvec%blocks )
  end subroutine par_block_array_free

  !=============================================================================
  subroutine par_block_array_alloc_all(bvec, p_bgraph)
    implicit none
    class(par_block_array_t), intent(out) :: bvec
    type(par_block_graph_t) , intent(in)  :: p_bgraph
    integer(ip)  :: ib
    type(par_graph_t), pointer :: p_graph
    
    bvec%nblocks = p_bgraph%get_nblocks()
    allocate ( bvec%blocks(bvec%nblocks) )
    do ib=1, bvec%nblocks
       p_graph => p_bgraph%get_block(ib,ib)
       call bvec%blocks(ib)%create ( p_graph%dof_dist, p_graph%p_env )
    end do

  end subroutine par_block_array_alloc_all

  !=============================================================================
  subroutine par_block_array_alloc_blocks ( bvec, nblocks)
    implicit none
    class(par_block_array_t), intent(out) :: bvec
    integer(ip)           , intent(in)  :: nblocks

    bvec%nblocks        = nblocks
    allocate ( bvec%blocks(nblocks) )
  end subroutine par_block_array_alloc_blocks

  !=============================================================================
  subroutine par_block_array_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(par_block_array_t), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(par_block_array_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call tvec%par_block_array_alloc_blocks(svec%nblocks)

    do ib=1, svec%nblocks
       call svec%blocks(ib)%create_view(start, end, tvec%blocks(ib))
    end do
  end subroutine par_block_array_create_view

  subroutine par_block_array_comm ( p_vec )
    implicit none

    ! Parameters
    type(par_block_array_t), intent(inout)         :: p_vec

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call p_vec%blocks(ib)%comm()
    end do

  end subroutine par_block_array_comm

  subroutine par_block_array_weight ( p_vec )
    implicit none

    ! Parameters
    type(par_block_array_t), intent(inout)         :: p_vec

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call p_vec%blocks(ib)%weight()
    end do

  end subroutine par_block_array_weight

  !=============================================================================
  subroutine par_block_array_dot (x, y, t)
    implicit none
    ! Parameters
    type(par_block_array_t), intent(in)  :: x
    type(par_block_array_t), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      aux = x%blocks(ib)%dot(y%blocks(ib))
      t = t + aux
    end do 
  end subroutine par_block_array_dot

  !=============================================================================
  subroutine par_block_array_nrm2(x,t)
    implicit none
    type(par_block_array_t), intent(in)  :: x
    real(rp)              , intent(out) :: t

    ! p_part%p_context is required within this subroutine
    assert ( associated(x%blocks(1)%p_env%p_context) )
    assert ( x%blocks(1)%p_env%p_context%created .eqv. .true.)

    if(x%blocks(1)%p_env%p_context%iam<0) return

    call par_block_array_dot (x, x, t)
    t = sqrt(t)
  end subroutine par_block_array_nrm2

  !=============================================================================
  subroutine par_block_array_copy(x,y)
    implicit none
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call y%blocks(ib)%copy(x%blocks(ib))
    end do 
  end subroutine par_block_array_copy

  subroutine par_block_array_zero(y)
    implicit none
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call y%blocks(ib)%init(0.0_rp)
    end do 

  end subroutine par_block_array_zero

  subroutine par_block_array_init(alpha, y)
    implicit none
    type(par_block_array_t), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call y%blocks(ib)%init(alpha)
    end do    
  end subroutine par_block_array_init
  
  subroutine par_block_array_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%scal(t, x%blocks(ib))
    end do 

  end subroutine par_block_array_scale

  subroutine par_block_array_mxpy(x,y)
    implicit none
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%axpby (-1.0_rp, x%blocks(ib), 1.0_rp)
    end do 
  end subroutine par_block_array_mxpy
  subroutine par_block_array_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%axpby(t, x%blocks(ib), 1.0_rp)
    end do 
  end subroutine par_block_array_axpy

  subroutine par_block_array_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%axpby ( 1.0_rp, x%blocks(ib),  t)
    end do 

  end subroutine par_block_array_aypx

  subroutine par_block_array_pxpy(x,y)
    implicit none
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%axpby ( 1.0_rp, x%blocks(ib), 1.0_rp )
    end do 
  end subroutine par_block_array_pxpy

  subroutine par_block_array_pxmy(x,y)
    implicit none
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib
       
    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call y%blocks(ib)%axpby ( 1.0_rp, x%blocks(ib), -1.0_rp )
    end do 
  end subroutine par_block_array_pxmy

  subroutine par_block_array_print (this,luout)
    implicit none
    class(par_block_array_t), intent(in) :: this
    integer(ip)             , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, this%nblocks
       write (*,*) 'Block-vector ', ib
       call this%blocks(ib)%print(luout)
    end do 
  end subroutine par_block_array_print

  ! alpha <- op1^T * op2
  function par_block_array_dot_tbp(op1,op2) result(alpha)
    implicit none
    ! Parameters
    class(par_block_array_t), intent(in)  :: op1
    class(abstract_vector_t)    , intent(in)  :: op2
    real(rp) :: alpha

    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
    class is (par_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       alpha = 0.0_rp
       do ib=1,op1%nblocks
          aux = op1%blocks(ib)%dot(op2%blocks(ib))
          alpha = alpha + aux
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%dot: unsupported op2 class'
       check(1==0)
    end select
    call op1%CleanTemp()
    call op2%CleanTemp()
  end function par_block_array_dot_tbp

  ! op1 <- op2 
  subroutine par_block_array_copy_tbp(op1,op2)
    implicit none
    ! Parameters
    class(par_block_array_t), intent(inout) :: op1
    class(abstract_vector_t)    , intent(in)    :: op2

    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%copy(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_array_copy_tbp

  ! op <- alpha
  subroutine par_block_array_init_tbp(op,alpha)
    implicit none
    class(par_block_array_t), intent(inout) :: op 
    real(rp)                 , intent(in)    :: alpha  
    ! Locals
    integer(ip) :: ib

    do ib=1, op%nblocks
       call op%blocks(ib)%init(alpha)
    end do
  end subroutine par_block_array_init_tbp

  ! op1 <- alpha * op2
  subroutine par_block_array_scal_tbp(op1,alpha,op2)
    implicit none
    ! Parameters 
    class(par_block_array_t), intent(inout) :: op1
    real(rp)                 , intent(in)    :: alpha
    class(abstract_vector_t)    , intent(in)    :: op2
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
       class is (par_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%scal(alpha,op2%blocks(ib))
       end do
       class default
       write(0,'(a)') 'par_block_vector_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_array_scal_tbp

  ! op1 <- alpha*op2 + beta*op1
  subroutine par_block_array_axpby_tbp(op1,alpha,op2,beta)
    implicit none
    class(par_block_array_t), intent(inout) :: op1
    real(rp)                 , intent(in)    :: alpha
    class(abstract_vector_t)    , intent(in)    :: op2
    real(rp)                 , intent(in)    :: beta
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_array_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%axpby(alpha,op2%blocks(ib),beta)
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_array_axpby_tbp

  ! alpha <- nrm2(op)
  function par_block_array_nrm2_tbp(op) result(alpha)
    implicit none
    class(par_block_array_t), intent(in) :: op
    real(rp) :: alpha

    call op%GuardTemp()
    alpha = op%dot(op)
    alpha = sqrt(alpha)
    call op%CleanTemp()
  end function par_block_array_nrm2_tbp

  ! op1 <- clone(op2) 
  subroutine par_block_array_tbp(op1,op2)
    implicit none
    ! Parameters
    class(par_block_array_t)    , intent(inout) :: op1
    class(abstract_vector_t), target, intent(in)    :: op2
 
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_array_t)
       op1%nblocks = op2%nblocks
       if(allocated(op1%blocks)) deallocate ( op1%blocks )
       allocate(op1%blocks(op1%nblocks))
       do ib=1,op1%nblocks
          call op1%blocks(ib)%clone(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%clone: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_array_tbp

  ! op <- comm(op)
  subroutine par_block_array_comm_tbp(op)
    implicit none
    class(par_block_array_t), intent(inout) :: op 
 
    ! Locals
    integer(ip) :: ib

    do ib=1,op%nblocks
       call op%blocks(ib)%comm()
    end do
    
  end subroutine par_block_array_comm_tbp

  subroutine par_block_array_free_tbp(this)
    implicit none
    class(par_block_array_t), intent(inout) :: this
    integer(ip)  :: ib
   
    do ib=1, this%nblocks
       call this%blocks(ib)%free()
    end do
    this%nblocks = 0
    deallocate( this%blocks )
  end subroutine par_block_array_free_tbp

end module par_block_array_names
