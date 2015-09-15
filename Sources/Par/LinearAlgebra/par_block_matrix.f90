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
module par_block_matrix_names
  ! Serial modules
  use types_names
  use memor_names
  use graph_names
  use matrix_names
  use vector_names
  use abstract_operator_names

  ! Parallel modules
  use par_matrix_names
  use par_scalar_array_names
  use par_graph_names
  use par_block_graph_names
  use par_block_array_names

  implicit none
# include "debug.i90"

  private

  ! Pointer to matrix
  type p_par_matrix_t
    type(par_matrix_t), pointer :: p_p_matrix
  end type p_par_matrix_t


  ! Block Matrix
  type, extends(abstract_operator_t) :: par_block_matrix_t
!    private ! IBM XLF 14.1 bug
    integer(ip)                     :: nblocks
    type(p_par_matrix_t), allocatable :: blocks(:,:)
  contains
    procedure :: alloc                   => par_block_matrix_alloc
    generic   :: alloc_block             => alloc_diagonal_block, alloc_offdiagonal_block
	procedure :: alloc_diagonal_block    => par_block_matrix_alloc_diagonal_block
    procedure :: alloc_offdiagonal_block => par_block_matrix_alloc_offdiagonal_block
    procedure :: set_block_to_zero       => par_block_matrix_set_block_to_zero
    procedure :: free                    => par_block_matrix_free_tbp
    procedure :: get_block               => par_block_matrix_get_block
    procedure :: get_nblocks             => par_block_matrix_get_nblocks
    procedure :: apply                   => par_block_matrix_apply
    procedure :: apply_fun               => par_block_matrix_apply_fun
  end type par_block_matrix_t

  ! Types
  public :: par_block_matrix_t

contains

  !=============================================================================
  subroutine par_block_matrix_alloc (bmat,bgraph,diagonal_blocks_symmetric,sign_diagonal_blocks)
    implicit none
    ! Parameters
    class(par_block_matrix_t), intent(inout) :: bmat
    type(par_block_graph_t)  , intent(in)    :: bgraph
	logical                  , intent(in)    :: diagonal_blocks_symmetric(:)
    integer(ip)              , intent(in)    :: sign_diagonal_blocks(:)

    ! Locals
    integer(ip) :: ib,jb
    type(par_graph_t), pointer :: p_graph

    assert ( size(diagonal_blocks_symmetric) == bgraph%get_nblocks() )
    assert ( size(sign_diagonal_blocks) == bgraph%get_nblocks() )

    bmat%nblocks = bgraph%get_nblocks()
    allocate ( bmat%blocks(bmat%nblocks,bmat%nblocks) )

    do ib=1, bmat%nblocks
      do jb=1, bmat%nblocks
           p_graph => bgraph%blocks(ib,jb)%p_p_graph
           if (associated(p_graph)) then
             allocate ( bmat%blocks(ib,jb)%p_p_matrix )
             if ( ib == jb ) then
               call par_matrix_alloc ( diagonal_blocks_symmetric(ib), p_graph, bmat%blocks(ib,jb)%p_p_matrix, sign_diagonal_blocks(ib) )
             else
               call par_matrix_alloc ( .false., p_graph, bmat%blocks(ib,jb)%p_p_matrix )
             end if
          else
             nullify ( bmat%blocks(ib,jb)%p_p_matrix )
          end if
       end do
    end do
  end subroutine par_block_matrix_alloc

  subroutine par_block_matrix_alloc_diagonal_block (bmat,ib,p_graph,diagonal_block_symmetric,diagonal_block_sign)
    implicit none
    ! Parameters
    class(par_block_matrix_t),   target  , intent(inout) :: bmat
    integer(ip)                          , intent(in)    :: ib
    type(par_graph_t)                    , intent(in)    :: p_graph
	logical                              , intent(in)    :: diagonal_block_symmetric
    integer(ip)                          , intent(in)    :: diagonal_block_sign
    assert ( associated ( bmat%blocks(ib,ib)%p_p_matrix ) )
    if ( .not. associated( bmat%blocks(ib,ib)%p_p_matrix) ) then
       allocate ( bmat%blocks(ib,ib)%p_p_matrix )
       call par_matrix_alloc ( diagonal_block_symmetric, p_graph, bmat%blocks(ib,ib)%p_p_matrix, diagonal_block_sign )
    end if
  end subroutine par_block_matrix_alloc_diagonal_block
  
    subroutine par_block_matrix_alloc_offdiagonal_block (bmat,ib,jb,p_graph)
    implicit none
    ! Parameters
    class(par_block_matrix_t),   target  , intent(inout) :: bmat
    integer(ip)                          , intent(in)    :: ib,jb
    type(par_graph_t)                    , intent(in)    :: p_graph
    assert ( associated ( bmat%blocks(ib,ib)%p_p_matrix ) )
    if ( .not. associated( bmat%blocks(ib,ib)%p_p_matrix) ) then
       allocate ( bmat%blocks(ib,ib)%p_p_matrix )
       call par_matrix_alloc ( .false., p_graph, bmat%blocks(ib,ib)%p_p_matrix )
    end if
  end subroutine par_block_matrix_alloc_offdiagonal_block

  subroutine par_block_matrix_set_block_to_zero (bmat,ib,jb)
    implicit none
    ! Parameters
    class(par_block_matrix_t), intent(inout) :: bmat
    integer(ip)           , intent(in)  :: ib,jb

    if ( associated(bmat%blocks(ib,jb)%p_p_matrix) ) then
       call par_matrix_free( bmat%blocks(ib,jb)%p_p_matrix )
       deallocate (bmat%blocks(ib,jb)%p_p_matrix)
       nullify ( bmat%blocks(ib,jb)%p_p_matrix )
    end if

  end subroutine par_block_matrix_set_block_to_zero

  function par_block_matrix_get_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(par_block_matrix_t), target, intent(in) :: bmat
    integer(ip)                    , intent(in) :: ib,jb
    type(par_matrix_t)               , pointer    :: par_block_matrix_get_block

    par_block_matrix_get_block =>  bmat%blocks(ib,jb)%p_p_matrix
  end function par_block_matrix_get_block

  function par_block_matrix_get_nblocks (bmat)
    implicit none
    ! Parameters
    class(par_block_matrix_t), target, intent(in) :: bmat
    integer(ip)                                :: par_block_matrix_get_nblocks
    par_block_matrix_get_nblocks = bmat%nblocks
  end function par_block_matrix_get_nblocks


  subroutine par_block_matrix_print (lunou, p_b_matrix)
    implicit none
    type(par_block_matrix_t), intent(in)    :: p_b_matrix
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i

    check(.false.)
  end subroutine par_block_matrix_print

  !=============================================================================
  subroutine par_block_matrix_free (p_b_matrix)
    implicit none
    class(par_block_matrix_t), intent(inout) :: p_b_matrix
    integer(ip) :: ib,jb

    do ib=1, p_b_matrix%nblocks
       do jb=1, p_b_matrix%nblocks
          if ( associated(p_b_matrix%blocks(ib,jb)%p_p_matrix) ) then
             call par_matrix_free( p_b_matrix%blocks(ib,jb)%p_p_matrix )
             deallocate (p_b_matrix%blocks(ib,jb)%p_p_matrix) 
          end if
       end do
    end do

    p_b_matrix%nblocks = 0
    deallocate ( p_b_matrix%blocks )
  
  end subroutine par_block_matrix_free

  subroutine par_block_matvec (a, x, y)
    implicit none
    ! Parameters
    type(par_block_matrix_t), intent(in)    :: a
    type(par_block_array_t), intent(in)    :: x
    type(par_block_array_t), intent(inout) :: y

    ! Locals
    type(par_scalar_array_t)       :: aux
    integer(ip)            :: ib, jb

    assert ( a%nblocks == x%nblocks )
    assert ( a%nblocks == y%nblocks )     


    do ib=1, a%nblocks
       y%blocks(ib)%state = part_summed
       call y%blocks(ib)%init(0.0_rp)
       call aux%clone ( y%blocks(ib) ) 
       do jb=1, a%nblocks
          if ( associated(a%blocks(ib,jb)%p_p_matrix) ) then
             ! aux <- A(ib,jb) * x(jb)
             call par_matvec ( a%blocks(ib,jb)%p_p_matrix, x%blocks(jb), aux ) 

             ! y(ib) <- y(ib) + aux 
             call y%blocks(ib)%axpby ( 1.0_rp, aux, 1.0_rp )
          end if
       end do
       call aux%free()
    end do
  end subroutine par_block_matvec

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine par_block_matrix_apply(op,x,y)
    implicit none
    class(par_block_matrix_t), intent(in)    :: op
    class(vector_t)    , intent(in)    :: x
    class(vector_t)    , intent(inouT) :: y
    ! Locals
    integer(ip)        :: ib,jb
    type(par_scalar_array_t) :: aux

    call x%GuardTemp()

    call y%init(0.0_rp)
    select type(x)
    class is (par_block_array_t)
       select type(y)
       class is(par_block_array_t)
          do ib=1,op%nblocks
             call aux%clone(y%blocks(ib))
             do jb=1,op%nblocks
                if ( associated(op%blocks(ib,jb)%p_p_matrix) ) then
                   ! aux <- A(ib,jb) * x(jb)
                   call par_matvec(op%blocks(ib,jb)%p_p_matrix,x%blocks(jb),aux)
                   ! y(ib) <- y(ib) + aux
                   call y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
                end if
             end do
             call aux%free()
          end do
       class default
          write(0,'(a)') 'par_block_matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'par_block_matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine par_block_matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function par_block_matrix_apply_fun(op,x) result(y)
    implicit none
    class(par_block_matrix_t), intent(in)  :: op
    class(vector_t)    , intent(in)  :: x
    class(vector_t)    , allocatable :: y 
    ! Locals
    integer(ip) :: ib,jb
    type(par_block_array_t), allocatable :: local_y
    type(par_scalar_array_t) :: aux

    select type(x)
    class is (par_block_array_t)
       allocate(local_y)
       call local_y%create(op%nblocks)
       do ib=1,op%nblocks
          call aux%clone(local_y%blocks(ib))
          do jb=1,op%nblocks
             ! aux <- A(ib,jb) * x(jb)
             call par_matvec(op%blocks(ib,jb)%p_p_matrix,x%blocks(jb),aux)
             ! y(ib) <- y(ib) + aux
             call local_y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
          end do
          call aux%free()
       end do
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'par_block_matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function par_block_matrix_apply_fun

  subroutine par_block_matrix_free_tbp(this)
    implicit none
    class(par_block_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb

    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          if ( associated(this%blocks(ib,jb)%p_p_matrix) ) then
             call par_matrix_free( this%blocks(ib,jb)%p_p_matrix )
             deallocate (this%blocks(ib,jb)%p_p_matrix) 
          end if
       end do
    end do

    this%nblocks = -1 
    deallocate ( this%blocks ) 
  end subroutine par_block_matrix_free_tbp

end module par_block_matrix_names
