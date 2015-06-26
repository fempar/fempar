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
module par_block_graph_names
  ! Serial modules
use types_names
use memor_names

  ! Parallel modules 
  use par_environment_names
  use dof_distribution_names
  use par_graph_names
  use block_dof_distribution_names

  implicit none
# include "debug.i90"

  private

  ! Pointer to graph
  type p_par_graph_t
    type(par_graph_t), pointer :: p_p_graph
  end type p_par_graph_t


  ! Block Graph 
  type par_block_graph_t
    private
    integer(ip)                    :: nblocks = -1
    type(p_par_graph_t), allocatable :: blocks(:,:)
  contains
    procedure :: alloc             => par_block_graph_alloc
    procedure :: alloc_block       => par_block_graph_alloc_block 
    procedure :: set_block_to_zero => par_block_graph_set_block_to_zero
    procedure :: free              => par_block_graph_free
    procedure :: get_block         => par_block_graph_get_block
    procedure :: get_nblocks       => par_block_graph_get_nblocks
  end type par_block_graph_t

  ! Types
  public :: par_block_graph_t

  ! Functions
  public :: par_block_graph_alloc, par_block_graph_alloc_block,       & 
            par_block_graph_set_block_to_zero, par_block_graph_print, & 
            par_block_graph_free, par_block_graph_get_block                                 

contains

  !=============================================================================
  subroutine par_block_graph_alloc (p_b_graph, blk_dof_dist)
    implicit none
    ! Parameters
    class(par_block_graph_t)              , intent(inout) :: p_b_graph
    type(block_dof_distribution_t),target, intent(in)     :: blk_dof_dist 

    ! Locals
    integer(ip) :: istat
    integer(ip) :: ib,jb

    p_b_graph%nblocks = blk_dof_dist%nblocks   
    allocate ( p_b_graph%blocks(p_b_graph%nblocks,p_b_graph%nblocks), stat=istat )
    check(istat==0)
    do ib=1, p_b_graph%nblocks
      do jb=1, p_b_graph%nblocks
           allocate ( p_b_graph%blocks(ib,jb)%p_p_graph, stat=istat )
           check(istat==0)
           call par_graph_create ( blk_dof_dist%get_block(ib), & 
                                   blk_dof_dist%get_block(jb), & 
                                   blk_dof_dist%p_env, & 
                                   p_b_graph%blocks(ib,jb)%p_p_graph )
      end do
    end do
  end subroutine par_block_graph_alloc

  subroutine par_block_graph_alloc_block (p_b_graph,dof_dist_rows,dof_dist_cols,p_env,ib,jb)
    implicit none
    ! Parameters
    class(par_block_graph_t), target, intent(inout) :: p_b_graph
    type(dof_distribution_t)        , intent(in)    :: dof_dist_rows
    type(dof_distribution_t)        , intent(in)    :: dof_dist_cols
    type(par_environment_t)         , intent(in)    :: p_env
    integer(ip)                   , intent(in)    :: ib,jb
    ! Locals
    integer(ip) :: istat

    if ( .not. associated( p_b_graph%blocks(ib,jb)%p_p_graph)) then
       allocate ( p_b_graph%blocks(ib,jb)%p_p_graph, stat=istat )
       check(istat==0)
       call par_graph_create ( dof_dist_rows, & 
                               dof_dist_cols, & 
                               p_env, & 
                               p_b_graph%blocks(ib,jb)%p_p_graph )
    end if
  end subroutine par_block_graph_alloc_block

  subroutine par_block_graph_set_block_to_zero (p_b_graph,ib,jb)
    implicit none
    ! Parameters
    class(par_block_graph_t), intent(inout) :: p_b_graph
    integer(ip)           , intent(in)   :: ib,jb
    ! Locals
    integer(ip) :: istat

    if ( associated(p_b_graph%blocks(ib,jb)%p_p_graph) ) then
       deallocate (p_b_graph%blocks(ib,jb)%p_p_graph, stat=istat)
       check(istat==0)
       ! AFM: to address this scenario. The graph might be partially or fully created!!!
       ! call par_graph_free ( p_b_graph%blocks(ib,jb)%p_p_graph, free_clean)
       nullify    (p_b_graph%blocks(ib,jb)%p_p_graph)
    end if
  end subroutine par_block_graph_set_block_to_zero
  
  function par_block_graph_get_block (p_b_graph,ib,jb)
    implicit none
    ! Parameters
    class(par_block_graph_t), target, intent(in) :: p_b_graph
    integer(ip)                   , intent(in) :: ib,jb
    type(par_graph_t)               , pointer    :: par_block_graph_get_block
    ! Locals
    integer(ip) :: istat

    par_block_graph_get_block =>  p_b_graph%blocks(ib,jb)%p_p_graph
  end function par_block_graph_get_block

  function par_block_graph_get_nblocks (p_p_graph)
    implicit none
    ! Parameters
    class(par_block_graph_t), target, intent(in) :: p_p_graph
    integer(ip)                                :: par_block_graph_get_nblocks
    par_block_graph_get_nblocks = p_p_graph%nblocks
  end function par_block_graph_get_nblocks

  subroutine par_block_graph_print (lunou, p_b_graph)
    implicit none
    class(par_block_graph_t), intent(in)    :: p_b_graph
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i

    check(.false.)
  end subroutine par_block_graph_print

  !=============================================================================
  subroutine par_block_graph_free (p_b_graph)
    implicit none
    class(par_block_graph_t), intent(inout) :: p_b_graph
    integer(ip) :: ib,jb
    ! Locals
    integer(ip) :: istat

    do ib=1, p_b_graph%nblocks
       do jb=1, p_b_graph%nblocks
          if ( associated(p_b_graph%blocks(ib,jb)%p_p_graph) ) then
             call par_graph_free ( p_b_graph%blocks(ib,jb)%p_p_graph)
             deallocate (p_b_graph%blocks(ib,jb)%p_p_graph, stat=istat) 
             ! AFM: At this point the graph MUST BE fully created
             check(istat==0)
          end if
       end do
    end do

    deallocate ( p_b_graph%blocks, stat=istat )
    check(istat==0)
    p_b_graph%nblocks = -1
  end subroutine par_block_graph_free

end module par_block_graph_names
