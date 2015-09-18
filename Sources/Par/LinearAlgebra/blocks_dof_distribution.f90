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
module blocks_dof_distribution_names
  ! Serial modules
  use types_names
  
  ! Parallel modules 
  use dof_distribution_names
  use par_environment_names
 
  implicit none
# include "debug.i90"
  private
  
  type blocks_dof_distribution_t
     integer(ip) :: nblocks = -1 
     type (dof_distribution_t), allocatable :: blocks(:)
     type (par_environment_t) , pointer :: p_env
  contains
     procedure :: create    => blocks_dof_distribution_create
     procedure :: free      => blocks_dof_distribution_free
     procedure :: get_block => blocks_dof_distribution_get_block
  end type blocks_dof_distribution_t

  ! Types
  public :: blocks_dof_distribution_t
  
contains

  subroutine blocks_dof_distribution_create(this,nblocks,p_env)
    implicit none
    class(blocks_dof_distribution_t)  , intent(inout) :: this 
    integer(ip)                       , intent(in)    :: nblocks
	type(par_environment_t), target   , intent(in)    :: p_env


    ! Locals
    integer(ip) :: istat

    ! Parallel environment MUST BE already created
    assert ( p_env%created )
    
    this%nblocks = nblocks
    allocate ( this%blocks(nblocks), stat=istat )
    check ( istat == 0 )
    this%p_env => p_env
  end subroutine blocks_dof_distribution_create

  subroutine blocks_dof_distribution_free(this)
    implicit none
    class(blocks_dof_distribution_t), intent(inout) :: this

    ! Locals
    integer(ip) :: istat, iblock

    ! Parallel environment MUST BE already created
    assert ( associated(this%p_env) )
    assert ( this%p_env%created )

    if( this%p_env%p_context%iam >= 0 ) then
      do iblock = 1, this%nblocks
        call dof_distribution_free ( this%blocks(iblock) )
      end do 
    end if 
    
    this%nblocks = -1 
    deallocate (this%blocks,stat=istat)
    check(istat == 0)
    nullify(this%p_env)
  end subroutine blocks_dof_distribution_free
  
  function blocks_dof_distribution_get_block(this,iblock)
    implicit none
    class(blocks_dof_distribution_t), target, intent(in) :: this
    integer(ip)                             , intent(in) :: iblock
    type(dof_distribution_t), pointer                    :: blocks_dof_distribution_get_block

    ! Parallel environment MUST BE already created
    assert ( associated(this%p_env) )
    assert ( this%p_env%created )

    blocks_dof_distribution_get_block => this%blocks(iblock) 
  end function blocks_dof_distribution_get_block
  
end module blocks_dof_distribution_names
