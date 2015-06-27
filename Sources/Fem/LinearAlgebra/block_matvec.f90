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
module block_matrix_vector_names
use types_names
  use block_matrix_names
  use block_vector_names
  use matrix_names
  use vector_names
  implicit none
# include "debug.i90"

  private

  public :: block_matvec

contains

  subroutine block_matvec (a, x, y)
    implicit none
    ! Parameters
    type(block_matrix_t), intent(in)    :: a
    type(block_vector_t), intent(in)    :: x
    type(block_vector_t), intent(inout) :: y

    ! Locals
    type(matrix_t), pointer :: f_matrix
    type(vector_t)          :: aux
    integer(ip)               :: ib, jb

    assert ( a%get_nblocks() == x%nblocks )
    assert ( a%get_nblocks() == y%nblocks )
    
    do ib=1, a%get_nblocks()
       call vector_zero  ( y%blocks(ib) )
       call vector_clone ( y%blocks(ib), aux ) 
       do jb=1, a%get_nblocks()
          f_matrix => a%get_block(ib,jb)
          if ( associated(f_matrix) ) then
             ! aux <- A(ib,jb) * x(jb)
             call matvec ( f_matrix, x%blocks(jb), aux )
             
             ! y(ib) <- y(ib) + aux 
             call vector_pxpy ( aux, y%blocks(ib) ) 
          end if
       end do
       call vector_free ( aux )
    end do
  end subroutine block_matvec

end module block_matrix_vector_names
