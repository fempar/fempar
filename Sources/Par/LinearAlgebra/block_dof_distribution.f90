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
module block_dof_distribution_names
  ! Serial modules
  use types
  
  ! Parallel modules 
  use dof_distribution_names
  use par_environment_names
 
  implicit none
# include "debug.i90"
  private
  
  type block_dof_distribution
     integer(ip)                          :: nblocks = -1 
     type (dof_distribution), allocatable :: blocks(:)
     type (par_environment)     , pointer :: p_env
  contains
     procedure :: alloc     => block_dof_distribution_alloc
     procedure :: free      => block_dof_distribution_free
     procedure :: get_block => block_dof_distribution_get_block
  end type block_dof_distribution

  ! Types
  public :: block_dof_distribution
  
contains

  subroutine block_dof_distribution_alloc(blk_dof_dist,p_env,nblocks)
    implicit none
    class(block_dof_distribution)     , intent(inout) :: blk_dof_dist 
    type(par_environment)     , target, intent(in)    :: p_env
    integer(ip)                       , intent(in)    :: nblocks

    ! Locals
    integer(ip) :: istat

    ! Parallel environment MUST BE already created
    assert ( p_env%created )
    
    blk_dof_dist%nblocks = nblocks
    allocate ( blk_dof_dist%blocks(nblocks), stat=istat )
    check ( istat == 0 )
    blk_dof_dist%p_env => p_env
  end subroutine block_dof_distribution_alloc

  subroutine block_dof_distribution_free(blk_dof_dist)
    implicit none
    class(block_dof_distribution), intent(inout) :: blk_dof_dist

    ! Locals
    integer(ip) :: istat, iblock

    ! Parallel environment MUST BE already created
    assert ( associated(blk_dof_dist%p_env) )
    assert ( blk_dof_dist%p_env%created )

    if( blk_dof_dist%p_env%p_context%iam >= 0 ) then
      do iblock = 1, blk_dof_dist%nblocks
        call dof_distribution_free ( blk_dof_dist%blocks(iblock) )
      end do 
    end if 
    
    blk_dof_dist%nblocks = -1 
    deallocate (blk_dof_dist%blocks,stat=istat)
    check(istat == 0)
    nullify(blk_dof_dist%p_env)
  end subroutine block_dof_distribution_free
  
  function block_dof_distribution_get_block(blk_dof_dist,iblock)
    implicit none
    class(block_dof_distribution), target, intent(in) :: blk_dof_dist
    integer(ip)                          , intent(in) :: iblock
    type(dof_distribution)      , pointer             :: block_dof_distribution_get_block

    ! Parallel environment MUST BE already created
    assert ( associated(blk_dof_dist%p_env) )
    assert ( blk_dof_dist%p_env%created )

    block_dof_distribution_get_block => blk_dof_dist%blocks(iblock) 
  end function block_dof_distribution_get_block
  
end module block_dof_distribution_names
