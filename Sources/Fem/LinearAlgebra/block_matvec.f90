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
module fem_block_matrix_vector
  use types
  use fem_block_matrix_names
  use fem_block_vector_names
  use fem_matrix_names
  use fem_vector_names
  use fem_matrix_vector
  implicit none
# include "debug.i90"

  private

  public :: fem_block_matvec

contains

  subroutine fem_block_matvec (a, x, y)
    implicit none
    ! Parameters
    type(fem_block_matrix), intent(in)    :: a
    type(fem_block_vector), intent(in)    :: x
    type(fem_block_vector), intent(inout) :: y

    ! Locals
    type(fem_vector)       :: aux
    integer(ip)            :: ib, jb

    assert ( a%nblocks == x%nblocks )
    assert ( a%nblocks == y%nblocks )
    
    do ib=1, a%nblocks
       call fem_vector_zero  ( y%blocks(ib) )
       call fem_vector_clone ( y%blocks(ib), aux ) 
       do jb=1, a%nblocks
          if ( associated(a%blocks(ib,jb)%p_f_matrix) ) then
             ! aux <- A(ib,jb) * x(jb)
             call fem_matvec ( a%blocks(ib,jb)%p_f_matrix, x%blocks(jb), aux )
             
             ! y(ib) <- y(ib) + aux 
             call fem_vector_pxpy ( aux, y%blocks(ib) ) 
          end if
       end do
       call fem_vector_free ( aux )
    end do
  end subroutine fem_block_matvec

end module fem_block_matrix_vector
