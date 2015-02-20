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
module graph_block_renum
  use types
  use memor
  use fem_graph_class
  use fem_blocks_class
  use renum_class
  use psb_sort_mod

# include "debug.i90"
  implicit none
  private

  ! Functions
  public :: fem_graph_split_into_blocks

contains

  subroutine fem_graph_split_into_blocks ( graph, blocks, ren, graphs_gtype, output, graphs )
    implicit none
    ! Parameters
    type (fem_graph),  intent(in) :: graph
    type (fem_blocks), intent(in) :: blocks
    type (renum)     , intent(in) :: ren
    integer(ip)      , intent(in) :: graphs_gtype(blocks%nb)
    logical          , intent(in) :: output(blocks%nb)
    type (fem_graph) , intent(inout) :: graphs(blocks%nb)

    ! Locals
    integer (ip) :: ib, j, offset, k, row_orig, col_new

    assert ( graph%type == csr_symm )
    
    offset = 0
    do ib=1, blocks%nb
       if ( output(ib) ) then 
          graphs(ib)%nv   = blocks%ib(ib+1)-blocks%ib(ib)
          graphs(ib)%nv2  = graphs(ib)%nv
          graphs(ib)%type = graphs_gtype(ib)
          call memalloc ( graphs(ib)%nv+1, graphs(ib)%ia,__FILE__,__LINE__)
          graphs(ib)%ia = 0
          do j=blocks%ib(ib), blocks%ib(ib+1)-1
             row_orig=ren%iperm(j)
             do k=graph%ia(row_orig),graph%ia(row_orig+1)-1
                col_new=ren%lperm(graph%ja(k))
                if ( j-offset > 0 .and. col_new-offset > 0 & 
                     .and. j-offset <= blocks%ib(ib+1)-offset-1 .and. &
                     col_new-offset <= blocks%ib(ib+1)-offset-1 ) then ! Entry belongs to current diagonal block
                   ! New entry on (j-offset, col_new-offset)
                   if ( graph%type == csr_symm ) then
                      if (graphs(ib)%type == csr_symm ) then
                         if ( j-offset <= col_new-offset ) then ! U-entry
                            graphs(ib)%ia(j-offset+1) = graphs(ib)%ia(j-offset+1)+1
                         else ! L-entry
                            graphs(ib)%ia(col_new-offset+1) = graphs(ib)%ia(col_new-offset+1)+1
                         end if
                      else ! csr
                         graphs(ib)%ia(j-offset+1) = graphs(ib)%ia(j-offset+1)+1
                         if ( j-offset /= col_new-offset) then
                            graphs(ib)%ia(col_new-offset+1) = graphs(ib)%ia(col_new-offset+1)+1
                         end if
                      end if
                   else ! csr
                      stop ! Not implemented yet
                   end if
                end if

             end do
          end do

          graphs(ib)%ia(1)=1 
          do j=1, blocks%ib(ib+1)-blocks%ib(ib)
             graphs(ib)%ia(j+1) = graphs(ib)%ia(j) + graphs(ib)%ia(j+1)  
          end do

          call memalloc ( graphs(ib)%ia(graphs(ib)%nv+1)-1, graphs(ib)%ja,__FILE__,__LINE__)

       end if
       offset = offset + (blocks%ib(ib+1)-blocks%ib(ib))

    end do


    offset = 0
    do ib=1, blocks%nb
       if ( output(ib) ) then 

          do j=blocks%ib(ib), blocks%ib(ib+1)-1
             row_orig=ren%iperm(j)
             do k=graph%ia(row_orig),graph%ia(row_orig+1)-1
                col_new=ren%lperm(graph%ja(k))
                if ( j-offset > 0 .and. col_new-offset > 0 & 
                     .and. j-offset <= blocks%ib(ib+1)-offset-1 .and. &
                     col_new-offset <= blocks%ib(ib+1)-offset-1 ) then ! Entry belongs to current diagonal block
                   ! New entry on (j-offset, col_new-offset)
                   if ( graph%type == csr_symm ) then
                      if (graphs(ib)%type == csr_symm ) then
                         if ( j-offset <= col_new-offset ) then ! U-entry
                            graphs(ib)%ja(graphs(ib)%ia(j-offset)) = col_new-offset
                            graphs(ib)%ia(j-offset) = graphs(ib)%ia(j-offset)+1
                         else ! L-entry
                            graphs(ib)%ja(graphs(ib)%ia(col_new-offset)) = j-offset
                            graphs(ib)%ia(col_new-offset) = graphs(ib)%ia(col_new-offset)+1
                         end if
                      else ! csr
                         graphs(ib)%ja(graphs(ib)%ia(j-offset)) = col_new-offset
                         graphs(ib)%ia(j-offset) = graphs(ib)%ia(j-offset)+1
                         if ( j-offset /= col_new-offset) then
                            graphs(ib)%ja(graphs(ib)%ia(col_new-offset)) = j-offset
                            graphs(ib)%ia(col_new-offset) = graphs(ib)%ia(col_new-offset)+1
                         end if
                      end if
                   else ! csr
                      stop ! Not implemented yet
                   end if
                end if

             end do
          end do

          do j=blocks%ib(ib+1)-blocks%ib(ib)+1,2,-1 
             graphs(ib)%ia(j) = graphs(ib)%ia(j-1)   
          end do
          graphs(ib)%ia(1)=1

          do j=1, graphs(ib)%nv
             call psb_hsort( graphs(ib)%ja(graphs(ib)%ia(j):graphs(ib)%ia(j+1)-1) )
          end do

       end if
       offset = offset + (blocks%ib(ib+1)-blocks%ib(ib))

    end do
  
  end subroutine fem_graph_split_into_blocks

end module graph_block_renum
