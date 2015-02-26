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
module fem_block_matrix_class
  use types
  use memor
  use fem_matrix_class
  implicit none
# include "debug.i90"

  private

  ! Pointer to matrix
  type p_fem_matrix
    type(fem_matrix), pointer :: p_f_matrix
  end type p_fem_matrix

  ! Block Matrix
  type fem_block_matrix
    integer(ip)                     :: nblocks
    type(p_fem_matrix), allocatable :: blocks(:,:)
  end type fem_block_matrix

  ! Types
  public :: fem_block_matrix

  ! Functions
  public :: fem_block_matrix_alloc, fem_block_matrix_alloc_block,       & 
            fem_block_matrix_set_block_to_zero, fem_block_matrix_print, & 
            fem_block_matrix_free,                                      & 
            fem_block_matrix_info, fem_block_matrix_zero

contains

  !=============================================================================
  subroutine fem_block_matrix_alloc(nblocks, bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: nblocks
    type(fem_block_matrix), intent(out) :: bmat
    integer(ip) :: ib,jb

    bmat%nblocks = nblocks
    allocate ( bmat%blocks(nblocks,nblocks) )
    do ib=1, nblocks 
      do jb=1, nblocks
           allocate ( bmat%blocks(ib,jb)%p_f_matrix )
      end do
    end do
  end subroutine fem_block_matrix_alloc

  subroutine fem_block_matrix_alloc_block (ib,jb,bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: ib,jb
    type(fem_block_matrix), intent(inout) :: bmat
    if ( .not. associated( bmat%blocks(ib,jb)%p_f_matrix)) then
       allocate ( bmat%blocks(ib,jb)%p_f_matrix )
    end if
  end subroutine

  subroutine fem_block_matrix_set_block_to_zero (ib,jb,bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: ib,jb
    type(fem_block_matrix), intent(inout) :: bmat

    if ( associated(bmat%blocks(ib,jb)%p_f_matrix) ) then
       deallocate (bmat%blocks(ib,jb)%p_f_matrix)
       nullify ( bmat%blocks(ib,jb)%p_f_matrix )
    end if
  end subroutine

  subroutine fem_block_matrix_zero(bmat)
    implicit none
    ! Parameters
    type(fem_block_matrix), intent(inout) :: bmat
    integer(ip)                           :: ib,jb

    do ib=1,bmat%nblocks
       do jb=1,bmat%nblocks
          call fem_matrix_zero (bmat%blocks(ib,jb)%p_f_matrix)
       end do
    end do

  end subroutine fem_block_matrix_zero


  subroutine fem_block_matrix_print (lunou, f_b_matrix)
    implicit none
    type(fem_block_matrix), intent(in)    :: f_b_matrix
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i

  end subroutine fem_block_matrix_print

  subroutine fem_block_matrix_info ( f_blk_mat, me, np )
    implicit none

    ! Parameters 
    type(fem_block_matrix) , intent(in)    :: f_blk_mat
    integer                , intent(out)   :: me
    integer                , intent(out)   :: np
    
    me = 0
    np = 1 
  end subroutine fem_block_matrix_info

  subroutine fem_block_matrix_free(f_b_matrix)
    implicit none
    type(fem_block_matrix), intent(inout) :: f_b_matrix
    integer(ip) :: ib,jb

    do ib=1, f_b_matrix%nblocks 
       do jb=1, f_b_matrix%nblocks
          if ( associated(f_b_matrix%blocks(ib,jb)%p_f_matrix) ) then
             deallocate (f_b_matrix%blocks(ib,jb)%p_f_matrix) 
          end if
       end do
    end do

    f_b_matrix%nblocks = 0
    deallocate ( f_b_matrix%blocks ) 
  end subroutine fem_block_matrix_free

end module fem_block_matrix_class
