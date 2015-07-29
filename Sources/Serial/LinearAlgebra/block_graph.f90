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
module block_graph_names
  ! Serial modules
use types_names
  use graph_names

  implicit none
# include "debug.i90"

  private

  ! Pointer to graph
  type p_graph_t
    type(graph_t), pointer :: p_f_graph
  end type p_graph_t


  ! Block Graph 
  type block_graph_t
!     private ! IBM XLF 14.1 bug
     integer(ip)                    :: nblocks = -1
     type(p_graph_t), allocatable :: blocks(:,:)
  contains
    procedure :: alloc             => block_graph_alloc
    procedure :: alloc_block       => block_graph_alloc_block 
    procedure :: set_block_to_zero => block_graph_set_block_to_zero
    procedure :: free              => block_graph_free
    procedure :: get_block         => block_graph_get_block
    procedure :: get_nblocks       => block_graph_get_nblocks
  end type block_graph_t

  ! Types
  public :: block_graph_t

  ! Functions
  public :: block_graph_alloc, block_graph_alloc_block,       & 
            block_graph_set_block_to_zero, block_graph_print, & 
            block_graph_free, block_graph_get_block, block_graph_get_nblocks                                 

contains

  !=============================================================================
  subroutine block_graph_alloc ( p_f_graph, nblocks)
    implicit none
    ! Parameters
    class(block_graph_t)              , intent(inout) :: p_f_graph
    integer(ip)                         , intent(in)    :: nblocks

    ! Locals
    integer(ip) :: istat
    integer(ip) :: ib,jb

    p_f_graph%nblocks = nblocks 
    allocate ( p_f_graph%blocks(nblocks,nblocks), stat=istat )
    check(istat==0)
    do ib=1, nblocks
      do jb=1, nblocks
           allocate ( p_f_graph%blocks(ib,jb)%p_f_graph, stat=istat )
           check(istat==0)
      end do
    end do
  end subroutine block_graph_alloc

  subroutine block_graph_alloc_block (p_f_graph,ib,jb)
    implicit none
    ! Parameters
    class(block_graph_t), target, intent(inout) :: p_f_graph
    integer(ip)                  , intent(in)    :: ib,jb
    ! Locals
    integer(ip) :: istat

    assert ( .not. associated( p_f_graph%blocks(ib,jb)%p_f_graph) )
    if ( .not. associated( p_f_graph%blocks(ib,jb)%p_f_graph)) then
       allocate ( p_f_graph%blocks(ib,jb)%p_f_graph, stat=istat )
       check(istat==0)
    end if
  end subroutine block_graph_alloc_block

  subroutine block_graph_set_block_to_zero (p_f_graph,ib,jb)
    implicit none
    ! Parameters
    class(block_graph_t), intent(inout) :: p_f_graph
    integer(ip)           , intent(in)   :: ib,jb
    ! Locals
    integer(ip) :: istat

    if ( associated(p_f_graph%blocks(ib,jb)%p_f_graph) ) then
       deallocate (p_f_graph%blocks(ib,jb)%p_f_graph, stat=istat)
       check(istat==0)
       ! AFM: to address this scenario. The graph might be partially or fully created!!!
       ! call graph_free ( p_f_graph%blocks(ib,jb)%p_f_graph, free_clean)
       nullify    (p_f_graph%blocks(ib,jb)%p_f_graph)
    end if
  end subroutine block_graph_set_block_to_zero
  
  function block_graph_get_block (p_f_graph,ib,jb)
    implicit none
    ! Parameters
    class(block_graph_t), target, intent(in) :: p_f_graph
    integer(ip)                   , intent(in) :: ib,jb
    type(graph_t)               , pointer    :: block_graph_get_block

    block_graph_get_block =>  p_f_graph%blocks(ib,jb)%p_f_graph
  end function block_graph_get_block

  function block_graph_get_nblocks (p_f_graph)
    implicit none
    ! Parameters
    class(block_graph_t), target, intent(in) :: p_f_graph
    integer(ip)                                :: block_graph_get_nblocks
    block_graph_get_nblocks = p_f_graph%nblocks
  end function block_graph_get_nblocks

  subroutine block_graph_print (lunou, p_f_graph)
    implicit none
    class(block_graph_t), intent(in)    :: p_f_graph
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i

    check(.false.)
  end subroutine block_graph_print

  !=============================================================================
  subroutine block_graph_free (p_f_graph)
    implicit none
    class(block_graph_t), intent(inout) :: p_f_graph
    integer(ip) :: ib,jb
    ! Locals
    integer(ip) :: istat

    do ib=1, p_f_graph%nblocks
       do jb=1, p_f_graph%nblocks
          if ( associated(p_f_graph%blocks(ib,jb)%p_f_graph) ) then
             call graph_free ( p_f_graph%blocks(ib,jb)%p_f_graph )
             deallocate (p_f_graph%blocks(ib,jb)%p_f_graph, stat=istat) 
             ! AFM: At this point the graph MUST BE fully created
             check(istat==0)
          end if
       end do
    end do

    deallocate ( p_f_graph%blocks, stat=istat )
    check(istat==0)
    p_f_graph%nblocks=-1 
 end subroutine block_graph_free

end module block_graph_names
