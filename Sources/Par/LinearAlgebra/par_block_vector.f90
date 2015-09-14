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
module par_block_vector_names
  ! Serial modules
  use types_names
  use memor_names
  use block_vector_names
  use serial_scalar_array_names
  use abstract_vector_names
    
  ! Parallel modules
  use par_vector_names
  use par_block_graph_names
  use par_graph_names

  implicit none
# include "debug.i90"

  private
 
  !=============================================================
  ! TO-CONSIDER:
  ! 
  ! x par_block_vector_dot and par_block_vector_nrm2 are not 
  !   optimal as a global communication (i.e., allreduce) is 
  !   required for each block. A more clever strategy could 
  !   replicate the body of par_vector_dot on each block of 
  !   par_block_vector_dot, and reduce the number of global 
  !   communications to one. However, with this strategy the code 
  !   would become significantly dirty. Does it really pay off ?
  ! 
  ! x Support for this parallel data structure in integrate.i90 
  !   would eliminate f_blk_vector member of par_block_vector and 
  !   par_block_vector_fill_complete method.
  !=============================================================

  ! par_vector
  type, extends(abstract_vector_t) :: par_block_vector_t
     integer(ip)                     :: nblocks = 0
     type(par_vector_t), allocatable :: blocks(:)

     ! **IMPORTANT NOTE**: This is an auxiliary data 
     ! structure provided in order to use SERIAL 
     ! block_vector assembly routines. The blocks of this 
     ! data structure are just VIEWS to the corresponding 
     ! counterparts in type(par_vector_t), allocatable :: blocks(:).
     ! This is required because currently integrate.i90 only
     ! accepts fem* data structures. If we provided support for 
     ! par* data structures in integrate.i90 we would not require 
     ! this aux. data structure
     type(block_vector_t)          :: f_blk_vector
     logical                       :: fill_completed
   contains
     procedure :: dot   => par_block_vector_dot_tbp
     procedure :: copy  => par_block_vector_copy_tbp
     procedure :: init  => par_block_vector_init_tbp
     procedure :: scal  => par_block_vector_scal_tbp
     procedure :: axpby => par_block_vector_axpby_tbp
     procedure :: nrm2  => par_block_vector_nrm2_tbp
     procedure :: clone => par_block_vector_clone_tbp
     procedure :: comm  => par_block_vector_comm_tbp
     procedure :: free  => par_block_vector_free_tbp
     procedure :: par_block_vector_alloc_blocks
     procedure :: par_block_vector_alloc_all
     generic   :: alloc => par_block_vector_alloc_blocks, par_block_vector_alloc_all
  end type par_block_vector_t

  ! Types
  public :: par_block_vector_t

  ! Functions
  public :: par_block_vector_free,                                              &
            par_block_vector_fill_complete, par_block_vector_create_view,       &
            par_block_vector_clone,         par_block_vector_comm,              &
            par_block_vector_weight,                                            &
            par_block_vector_dot,           par_block_vector_nrm2,              &
            par_block_vector_copy,          par_block_vector_zero,              &
            par_block_vector_init,          par_block_vector_scale,             & 
            par_block_vector_mxpy,          par_block_vector_axpy,              &
            par_block_vector_aypx,          par_block_vector_pxpy,              &
            par_block_vector_pxmy,          par_block_vector_print     
contains

  !=============================================================================
  subroutine par_block_vector_free (bvec)
    implicit none
    type(par_block_vector_t), intent(inout) :: bvec

    ! Locals
    integer(ip) :: ib

    do ib=1, bvec%nblocks
       call par_vector_free ( bvec%blocks(ib) )
    end do

    bvec%nblocks = 0
    deallocate( bvec%blocks )
    if ( bvec%fill_completed ) then
      call block_vector_free ( bvec%f_blk_vector )
      bvec%fill_completed = .false.
    end if
  end subroutine par_block_vector_free

  !=============================================================================
  subroutine par_block_vector_alloc_all(bvec, p_bgraph)
    implicit none
    class(par_block_vector_t), intent(out) :: bvec
    type(par_block_graph_t)  , intent(in)  :: p_bgraph
    integer(ip)  :: ib
    type(par_graph_t), pointer :: p_graph
    
    bvec%nblocks = p_bgraph%get_nblocks()
    allocate ( bvec%blocks(bvec%nblocks) )
    do ib=1, bvec%nblocks
       p_graph => p_bgraph%get_block(ib,ib)
       call par_vector_alloc ( p_graph%dof_dist, p_graph%p_env, bvec%blocks(ib) )
    end do

  end subroutine par_block_vector_alloc_all

  !=============================================================================
  subroutine par_block_vector_alloc_blocks ( bvec, nblocks)
    implicit none
    class(par_block_vector_t), intent(out) :: bvec
    integer(ip)           , intent(in)  :: nblocks

    bvec%nblocks        = nblocks
    bvec%fill_completed = .false.
    allocate ( bvec%blocks(nblocks) )
  end subroutine par_block_vector_alloc_blocks

  !=============================================================================
  subroutine par_block_vector_fill_complete (bvec)
    implicit none
    ! Parameters
    type(par_block_vector_t), intent(inout) :: bvec
    
    ! Locals
    integer(ip) :: ib  
  
    assert ( .not. bvec%fill_completed )
  
    call bvec%f_blk_vector%block_vector_alloc_blocks(bvec%nblocks)

    do ib=1, bvec%nblocks
       call serial_scalar_array_create_view ( bvec%blocks(ib)%f_vector,        &  
                                   & 1, bvec%blocks(ib)%f_vector%neq, &
                                   & bvec%f_blk_vector%blocks(ib))
    end do
    bvec%fill_completed = .true.
  end subroutine par_block_vector_fill_complete

  !=============================================================================
  subroutine par_block_vector_create_view (svec, start, end, tvec)
    implicit none
    ! Parameters
    type(par_block_vector_t), intent(in)  :: svec
    integer(ip)     , intent(in)        :: start
    integer(ip)     , intent(in)        :: end
    type(par_block_vector_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call tvec%par_block_vector_alloc_blocks(svec%nblocks)

    do ib=1, svec%nblocks
       call par_vector_create_view (svec%blocks(ib), start, end, tvec%blocks(ib))
    end do
  end subroutine par_block_vector_create_view

  !=============================================================================
  subroutine par_block_vector_clone ( svec, tvec )
    implicit none
    ! Parameters
    type(par_block_vector_t), intent( in ) :: svec
    type(par_block_vector_t), intent(out) :: tvec
 
    ! Locals
    integer(ip) :: ib

    call tvec%par_block_vector_alloc_blocks(svec%nblocks)

    do ib=1, svec%nblocks
       call par_vector_clone (svec%blocks(ib), tvec%blocks(ib))
    end do

  end subroutine par_block_vector_clone

  subroutine par_block_vector_comm ( p_vec )
    implicit none

    ! Parameters
    type(par_block_vector_t), intent(inout)         :: p_vec

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call par_vector_comm ( p_vec%blocks(ib) )
    end do

  end subroutine par_block_vector_comm

  subroutine par_block_vector_weight ( p_vec )
    implicit none

    ! Parameters
    type(par_block_vector_t), intent(inout)         :: p_vec

    ! Local variables
    integer(ip) :: ib
    do ib=1, p_vec%nblocks
       call par_vector_weight ( p_vec%blocks(ib) )
    end do

  end subroutine par_block_vector_weight

  !=============================================================================
  subroutine par_block_vector_dot (x, y, t)
    implicit none
    ! Parameters
    type(par_block_vector_t), intent(in)  :: x
    type(par_block_vector_t), intent(in)  :: y
    real(rp)              , intent(out) :: t
     
    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )

    t = 0.0_rp
    do ib=1,x%nblocks
      call par_vector_dot ( x%blocks(ib), y%blocks(ib), aux )
      t = t + aux
    end do 
  end subroutine par_block_vector_dot

  !=============================================================================
  subroutine par_block_vector_nrm2(x,t)
    implicit none
    type(par_block_vector_t), intent(in)  :: x
    real(rp)              , intent(out) :: t

    ! p_part%p_context is required within this subroutine
    assert ( associated(x%blocks(1)%p_env%p_context) )
    assert ( x%blocks(1)%p_env%p_context%created .eqv. .true.)

    if(x%blocks(1)%p_env%p_context%iam<0) return

    call par_block_vector_dot (x, x, t)
    t = sqrt(t)
  end subroutine par_block_vector_nrm2

  !=============================================================================
  subroutine par_block_vector_copy(x,y)
    implicit none
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y

    ! Locals
    integer(ip) :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, x%nblocks
      call par_vector_copy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_copy

  subroutine par_block_vector_zero(y)
    implicit none
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip) :: ib

    do ib=1, y%nblocks
      call par_vector_zero ( y%blocks(ib) )
    end do 

  end subroutine par_block_vector_zero

  subroutine par_block_vector_init(alpha, y)
    implicit none
    type(par_block_vector_t), intent(inout) :: y 
    real(rp), intent(in)                  :: alpha  
    ! Locals
    integer(ip)                           :: ib

    do ib=1, y%nblocks
      call par_vector_init ( alpha, y%blocks(ib) )
    end do    
  end subroutine par_block_vector_init
  
  subroutine par_block_vector_scale(t, x, y)
    implicit none
    ! Parameters 
    real(rp)              , intent(in)    :: t
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_scale ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine par_block_vector_scale

  subroutine par_block_vector_mxpy(x,y)
    implicit none
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_mxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_mxpy
  subroutine par_block_vector_axpy(t,x,y)
    implicit none
    real(rp)   , intent(in)         :: t
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_axpy ( t, x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_axpy

  subroutine par_block_vector_aypx(t,x,y)
    implicit none
    real(rp)        , intent(in)    :: t
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y

    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_aypx ( t, x%blocks(ib), y%blocks(ib) )
    end do 

  end subroutine par_block_vector_aypx

  subroutine par_block_vector_pxpy(x,y)
    implicit none
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib

    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_pxpy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_pxpy

  subroutine par_block_vector_pxmy(x,y)
    implicit none
    type(par_block_vector_t), intent(in)    :: x
    type(par_block_vector_t), intent(inout) :: y
    ! Locals
    integer(ip)                           :: ib
       
    assert ( x%nblocks == y%nblocks )
    do ib=1, y%nblocks
      call par_vector_pxmy ( x%blocks(ib), y%blocks(ib) )
    end do 
  end subroutine par_block_vector_pxmy

  subroutine par_block_vector_print (luout, x)
    implicit none
    type(par_block_vector_t), intent(in) :: x
    integer(ip)           , intent(in) :: luout
    
    ! Locals
    integer(ip) :: ib

    do ib=1, x%nblocks
       write (*,*) 'Block-vector ', ib
       call par_vector_print ( luout, x%blocks(ib) )
    end do 
  end subroutine par_block_vector_print

  ! alpha <- op1^T * op2
  function par_block_vector_dot_tbp(op1,op2) result(alpha)
    implicit none
    ! Parameters
    class(par_block_vector_t), intent(in)  :: op1
    class(abstract_vector_t)    , intent(in)  :: op2
    real(rp) :: alpha

    ! Locals
    real(rp)    :: aux
    integer(ip) :: ib

    call op1%GuardTemp()
    call op2%GuardTemp()
    select type(op2)
    class is (par_block_vector_t)
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
  end function par_block_vector_dot_tbp

  ! op1 <- op2 
  subroutine par_block_vector_copy_tbp(op1,op2)
    implicit none
    ! Parameters
    class(par_block_vector_t), intent(inout) :: op1
    class(abstract_vector_t)    , intent(in)    :: op2

    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%copy(op2%blocks(ib))
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%copy: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_vector_copy_tbp

  ! op <- alpha
  subroutine par_block_vector_init_tbp(op,alpha)
    implicit none
    class(par_block_vector_t), intent(inout) :: op 
    real(rp)                 , intent(in)    :: alpha  
    ! Locals
    integer(ip) :: ib

    do ib=1, op%nblocks
       call op%blocks(ib)%init(alpha)
    end do
  end subroutine par_block_vector_init_tbp

  ! op1 <- alpha * op2
  subroutine par_block_vector_scal_tbp(op1,alpha,op2)
    implicit none
    ! Parameters 
    class(par_block_vector_t), intent(inout) :: op1
    real(rp)                 , intent(in)    :: alpha
    class(abstract_vector_t)    , intent(in)    :: op2
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
       class is (par_block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%scal(alpha,op2%blocks(ib))
       end do
       class default
       write(0,'(a)') 'par_block_vector_t%scal: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_vector_scal_tbp

  ! op1 <- alpha*op2 + beta*op1
  subroutine par_block_vector_axpby_tbp(op1,alpha,op2,beta)
    implicit none
    class(par_block_vector_t), intent(inout) :: op1
    real(rp)                 , intent(in)    :: alpha
    class(abstract_vector_t)    , intent(in)    :: op2
    real(rp)                 , intent(in)    :: beta
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_vector_t)
       assert ( op1%nblocks == op2%nblocks )
       do ib=1,op1%nblocks
          call op1%blocks(ib)%axpby(alpha,op2%blocks(ib),beta)
       end do
    class default
       write(0,'(a)') 'par_block_vector_t%axpby: unsupported op2 class'
       check(1==0)
    end select
    call op2%CleanTemp()
  end subroutine par_block_vector_axpby_tbp

  ! alpha <- nrm2(op)
  function par_block_vector_nrm2_tbp(op) result(alpha)
    implicit none
    class(par_block_vector_t), intent(in) :: op
    real(rp) :: alpha

    call op%GuardTemp()
    alpha = op%dot(op)
    alpha = sqrt(alpha)
    call op%CleanTemp()
  end function par_block_vector_nrm2_tbp

  ! op1 <- clone(op2) 
  subroutine par_block_vector_clone_tbp(op1,op2)
    implicit none
    ! Parameters
    class(par_block_vector_t)    , intent(inout) :: op1
    class(abstract_vector_t), target, intent(in)    :: op2
 
    ! Locals
    integer(ip) :: ib

    call op2%GuardTemp()
    select type(op2)
    class is (par_block_vector_t)
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
  end subroutine par_block_vector_clone_tbp

  ! op <- comm(op)
  subroutine par_block_vector_comm_tbp(op)
    implicit none
    class(par_block_vector_t), intent(inout) :: op 
 
    ! Locals
    integer(ip) :: ib

    do ib=1,op%nblocks
       call op%blocks(ib)%comm()
    end do
    
  end subroutine par_block_vector_comm_tbp

  subroutine par_block_vector_free_tbp(this)
    implicit none
    class(par_block_vector_t), intent(inout) :: this
    integer(ip)  :: ib
   
    do ib=1, this%nblocks
       call this%blocks(ib)%free()
    end do
    this%nblocks = 0
    deallocate( this%blocks )
  end subroutine par_block_vector_free_tbp

end module par_block_vector_names
