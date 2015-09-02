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
module block_matrix_names
  use types_names
  use memor_names
  use graph_names
  use block_graph_names
  use matrix_names
  use block_vector_names
  use vector_names
  
  ! Abstract types
  use base_operand_names
  use base_operator_names
  implicit none
# include "debug.i90"


  private

  ! Pointer to matrix
  type p_matrix_t
    type(matrix_t), pointer :: p_f_matrix
  end type p_matrix_t

  ! Block Matrix
  type, extends(base_operator_t):: block_matrix_t
!    private ! IBM XLF 14.1 bug
    integer(ip)                     :: nblocks = -1
    type(p_matrix_t), allocatable :: blocks(:,:)
  contains
    procedure :: alloc             => block_matrix_alloc
    procedure :: alloc_block       => block_matrix_alloc_block
    procedure :: set_block_to_zero => block_matrix_set_block_to_zero
    procedure :: free              => block_matrix_free_tbp
    procedure :: get_block         => block_matrix_get_block
    procedure :: get_nblocks       => block_matrix_get_nblocks
    procedure :: apply             => block_matrix_apply
    procedure :: apply_fun         => block_matrix_apply_fun
  end type block_matrix_t

  ! Types
  public :: block_matrix_t

  ! Functions
  public :: block_matrix_alloc, block_matrix_alloc_block,       & 
            block_matrix_set_block_to_zero, block_matrix_print, & 
            block_matrix_free, block_matrix_info

contains

  !=============================================================================
  subroutine block_matrix_alloc(bmat, bgraph, sign)
    implicit none
    ! Parameters
    class(block_matrix_t), intent(inout) :: bmat
    type(block_graph_t)  , intent(in)    :: bgraph
    integer(ip), optional  , intent(in)    :: sign(:)

    integer(ip)                             :: ib,jb
    type(graph_t), pointer                :: f_graph

    if ( present(sign) ) then
      assert ( size(sign) == bgraph%get_nblocks() )
    end if

    bmat%nblocks = bgraph%get_nblocks()
    allocate ( bmat%blocks(bmat%nblocks,bmat%nblocks) )
    do ib=1, bmat%nblocks 
      do jb=1, bmat%nblocks
           f_graph => bgraph%blocks(ib,jb)%p_f_graph
           if (associated(f_graph)) then
              allocate ( bmat%blocks(ib,jb)%p_f_matrix )
              if ( (ib == jb) .and. present(sign) ) then
                if ( f_graph%type == csr ) then
                   call matrix_alloc ( csr_mat, symm_false, f_graph, bmat%blocks(ib,jb)%p_f_matrix, sign(ib) )
                else if ( f_graph%type == csr_symm ) then
                   call matrix_alloc ( csr_mat, symm_true, f_graph, bmat%blocks(ib,jb)%p_f_matrix, sign(ib) )
                end if 
              else
                if ( ib == jb ) then
                  if ( f_graph%type == csr ) then
                     call matrix_alloc(csr_mat, symm_false, f_graph, bmat%blocks(ib,jb)%p_f_matrix)
                  else if ( f_graph%type == csr_symm ) then
                     call matrix_alloc(csr_mat, symm_true, f_graph, bmat%blocks(ib,jb)%p_f_matrix)
                  end if
                else
                  call matrix_alloc(csr_mat, symm_false, f_graph, bmat%blocks(ib,jb)%p_f_matrix)
                end if
              end if
           else
              nullify ( bmat%blocks(ib,jb)%p_f_matrix )
           end if
      end do
    end do
  end subroutine block_matrix_alloc

  subroutine block_matrix_alloc_block (bmat,ib,jb,f_graph,sign)
    implicit none
    ! Parameters
    class(block_matrix_t), intent(inout) :: bmat
    integer(ip)           , intent(in)     :: ib,jb
    type(graph_t)       , intent(in)     :: f_graph
    integer(ip), optional , intent(in)     :: sign

    assert ( associated ( bmat%blocks(ib,jb)%p_f_matrix ) )
    if ( .not. associated( bmat%blocks(ib,jb)%p_f_matrix) ) then
       allocate ( bmat%blocks(ib,jb)%p_f_matrix )
       if ( (ib == jb) ) then
          if ( f_graph%type == csr ) then
             call matrix_alloc ( csr_mat, symm_false, f_graph, bmat%blocks(ib,jb)%p_f_matrix, sign )
          else if ( f_graph%type == csr_symm ) then
             call matrix_alloc ( csr_mat, symm_true, f_graph, bmat%blocks(ib,jb)%p_f_matrix, sign )
          end if
       else
          call matrix_alloc ( csr_mat, symm_false, f_graph, bmat%blocks(ib,jb)%p_f_matrix )
       end if
    end if
  end subroutine block_matrix_alloc_block

  subroutine block_matrix_set_block_to_zero (bmat,ib,jb)
    implicit none
    ! Parameters
    class(block_matrix_t), intent(inout) :: bmat
    integer(ip)           , intent(in)    :: ib,jb

    if ( associated(bmat%blocks(ib,jb)%p_f_matrix) ) then
       call matrix_free( bmat%blocks(ib,jb)%p_f_matrix )
       deallocate (bmat%blocks(ib,jb)%p_f_matrix)
       nullify ( bmat%blocks(ib,jb)%p_f_matrix )
    end if
  end subroutine

  subroutine block_matrix_print (lunou, f_b_matrix)
    implicit none
    type(block_matrix_t), intent(in)    :: f_b_matrix
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i
    check(.false.)
  end subroutine block_matrix_print

  subroutine block_matrix_free(bmat)
    implicit none
    class(block_matrix_t), intent(inout) :: bmat
    integer(ip) :: ib,jb

    do ib=1, bmat%nblocks 
       do jb=1, bmat%nblocks
          if ( associated(bmat%blocks(ib,jb)%p_f_matrix) ) then
             call matrix_free( bmat%blocks(ib,jb)%p_f_matrix )
             deallocate (bmat%blocks(ib,jb)%p_f_matrix) 
          end if
       end do
    end do

    bmat%nblocks = -1 
    deallocate ( bmat%blocks ) 
  end subroutine block_matrix_free

  function block_matrix_get_block (bmat,ib,jb)
    implicit none
    ! Parameters
    class(block_matrix_t), target, intent(in) :: bmat
    integer(ip)                    , intent(in) :: ib,jb
    type(matrix_t)               , pointer    :: block_matrix_get_block

    block_matrix_get_block =>  bmat%blocks(ib,jb)%p_f_matrix
  end function block_matrix_get_block

  function block_matrix_get_nblocks (bmat)
    implicit none
    ! Parameters
    class(block_matrix_t), target, intent(in) :: bmat
    integer(ip)                                :: block_matrix_get_nblocks
    block_matrix_get_nblocks = bmat%nblocks
  end function block_matrix_get_nblocks

  subroutine block_matrix_info ( f_blk_mat, me, np )
    implicit none

    ! Parameters 
    type(block_matrix_t) , intent(in)    :: f_blk_mat
    integer                , intent(out)   :: me
    integer                , intent(out)   :: np
    
    me = 0
    np = 1 
  end subroutine block_matrix_info

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine block_matrix_apply(op,x,y)
    implicit none
    class(block_matrix_t), intent(in)    :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inouT) :: y
    ! Locals
    integer(ip) :: ib,jb
    type(vector_t) :: aux

    call x%GuardTemp()

    call y%init(0.0_rp)
    select type(x)
    class is (block_vector_t)
       select type(y)
       class is(block_vector_t)
          do ib=1,op%nblocks
             call aux%clone(y%blocks(ib))
             do jb=1,op%nblocks
                if ( associated(op%blocks(ib,jb)%p_f_matrix) ) then
                   ! aux <- A(ib,jb) * x(jb)
                   call matvec(op%blocks(ib,jb)%p_f_matrix,x%blocks(jb),aux)
                   ! y(ib) <- y(ib) + aux
                   call y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
                end if
             end do
             call aux%free()
          end do
       class default
          write(0,'(a)') 'block_matrix_t%apply: unsupported y class'
          check(1==0)
       end select
    class default
       write(0,'(a)') 'block_matrix_t%apply: unsupported x class'
       check(1==0)
    end select

    call x%CleanTemp()

  end subroutine block_matrix_apply

  ! op%apply(x)
  ! Allocates room for (temporary) y
  function block_matrix_apply_fun(op,x) result(y)
    implicit none
    class(block_matrix_t), intent(in)  :: op
    class(base_operand_t), intent(in)  :: x
    class(base_operand_t), allocatable :: y 
    ! Locals
    integer(ip) :: ib,jb
    type(block_vector_t), allocatable :: local_y
    type(vector_t) :: aux

    select type(x)
    class is (block_vector_t)
       allocate(local_y)
       call local_y%block_vector_alloc_blocks(op%nblocks)
       do ib=1,op%nblocks
          call aux%clone(local_y%blocks(ib))
          do jb=1,op%nblocks
             ! aux <- A(ib,jb) * x(jb)
             call matvec(op%blocks(ib,jb)%p_f_matrix,x%blocks(jb),aux)
             ! y(ib) <- y(ib) + aux
             call local_y%blocks(ib)%axpby(1.0_rp,aux,1.0_rp)
          end do
          call aux%free()
       end do
       call move_alloc(local_y, y)
       call y%SetTemp()
    class default
       write(0,'(a)') 'block_matrix_t%apply_fun: unsupported x class'
       check(1==0)
    end select
  end function block_matrix_apply_fun

    subroutine block_matrix_free_tbp(this)
    implicit none
    class(block_matrix_t), intent(inout) :: this
    integer(ip) :: ib,jb

    do ib=1, this%nblocks 
       do jb=1, this%nblocks
          if ( associated(this%blocks(ib,jb)%p_f_matrix) ) then
             call matrix_free( this%blocks(ib,jb)%p_f_matrix )
             deallocate (this%blocks(ib,jb)%p_f_matrix) 
          end if
       end do
    end do

    this%nblocks = -1 
    deallocate ( this%blocks ) 
  end subroutine block_matrix_free_tbp

end module block_matrix_names
