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
module matrix_block_renum
  use types
  use memor
  use fem_matrix_class
  use fem_graph_class
  use fem_blocks_class
  use renum_class
  use graph_block_renum
  use psb_sort_mod

# include "debug.i90"
  implicit none
  private

  ! Functions
  public :: fem_matrix_split_into_blocks

contains

  subroutine fem_matrix_split_into_blocks ( matrix, blocks, ren, sign, graphs, output, matrices )
    implicit none
    ! Parameters
    type (fem_matrix),  intent(in) :: matrix
    type (fem_blocks), intent(in) :: blocks
    type (renum)     , intent(in) :: ren
    type (fem_graph), intent(in), target :: graphs(blocks%nb)
    integer(ip)     , intent(in)         :: sign(blocks%nb)
    logical         , intent(in)         :: output(blocks%nb)
    type (fem_matrix), intent(inout) :: matrices(blocks%nb)

    ! Locals
    integer (ip) :: ib, j, offset, k, row_orig, col_new, l

    assert ( matrix%gr%type == csr_symm )
    assert ( matrix%type == csr_mat )
    assert ( matrix%storage == scal )


    offset = 0
    do ib=1, blocks%nb

       if ( output(ib) ) then
          matrices(ib)%nd1     = matrix%nd1 ! With the old version of the codes (wo dof_handler), this does not work
          ! In particular, the value of nd for each block should be provided
          matrices(ib)%nd2     = matrix%nd2
          matrices(ib)%type    = csr_mat
          if ( graphs(ib)%type == csr_symm  ) then
             matrices(ib)%symm = symm_true
          else
             matrices(ib)%symm = symm_false
          end if
          matrices(ib)%storage = matrix%storage
          matrices(ib)%sign    = sign(ib)
          matrices(ib)%gr => graphs(ib)

          call memalloc ( 1, 1, graphs(ib)%ia(graphs(ib)%nv+1)-1, matrices(ib)%a, __FILE__,__LINE__)

          do j=blocks%ib(ib), blocks%ib(ib+1)-1
             row_orig=ren%iperm(j)
             do k=matrix%gr%ia(row_orig),matrix%gr%ia(row_orig+1)-1
                col_new=ren%lperm(matrix%gr%ja(k))
                if ( j-offset > 0 .and. col_new-offset > 0 & 
                     .and. j-offset <= blocks%ib(ib+1)-offset-1 .and. &
                     col_new-offset <= blocks%ib(ib+1)-offset-1 ) then ! Entry belongs to current diagonal block
                   ! New entry on (j-offset, col_new-offset)
                   if ( matrix%gr%type == csr_symm ) then
                      if (graphs(ib)%type == csr_symm ) then
                         if ( j-offset <= col_new-offset ) then ! U-entry
                            l=graphs(ib)%ia(j-offset)
                            do while (graphs(ib)%ja(l) < col_new-offset)
                               l=l+1
                            end do
                            matrices(ib)%a(1,1,l) = matrix%a(1,1,k)
                         else ! L-entry
                            l=graphs(ib)%ia(col_new-offset)
                            do while (graphs(ib)%ja(l) < j-offset)
                               l=l+1
                            end do
                            matrices(ib)%a(1,1,l) = matrix%a(1,1,k)
                         end if
                      else ! csr
                         l=graphs(ib)%ia(j-offset)
                         do while (graphs(ib)%ja(l) < col_new-offset)
                            l=l+1
                         end do
                         matrices(ib)%a(1,1,l) = matrix%a(1,1,k)

                         if ( j-offset /= col_new-offset) then
                            l=graphs(ib)%ia(col_new-offset)
                            do while (graphs(ib)%ja(l) < j-offset)
                               l=l+1
                            end do
                            matrices(ib)%a(1,1,l) = matrix%a(1,1,k)
                         end if
                      end if
                   else ! csr
                      stop ! Not implemented yet
                   end if
                end if
             end do
          end do

          ! do j=blocks%ib(ib+1)-blocks%ib(ib)+1,2,-1 
          !   graphs(ib)%ia(j) = graphs(ib)%ia(j-1)   
          ! end do
          ! graphs(ib)%ia(1)=1
       end if

       offset = offset + (blocks%ib(ib+1)-blocks%ib(ib))
    end do

  
  end subroutine fem_matrix_split_into_blocks

end module matrix_block_renum
