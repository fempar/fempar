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
module par_update_names
  ! Serial modules
  use types_names
  use memor_names
  use fem_update_names

  ! Parallel modules
  use par_fe_space_names
  use par_conditions_names
  use par_vector_names
  use par_block_vector_names
  
  implicit none
# include "debug.i90"
  private

  interface par_update_solution
     module procedure par_update_solution_mono, par_update_solution_block
  end interface par_update_solution

  ! Functions
  public :: par_update_strong_dirichlet_bcond, par_update_analytical_bcond, par_update_solution, &
       &    par_update_nonlinear
  
contains

  !==================================================================================================
  subroutine par_update_strong_dirichlet_bcond(p_fe_space, p_cond)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from par_conditions values.  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(par_fe_space_t) , intent(inout) :: p_fe_space
    type(par_conditions_t), intent(in)    :: p_cond

    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )

    ! If fine task call serial subroutine
    if( p_fe_space%p_trian%p_env%am_i_fine_task() ) then
       call fem_update_strong_dirichlet_bcond( p_fe_space%fe_space, p_cond%f_conditions )
    end if

  end subroutine par_update_strong_dirichlet_bcond

  !==================================================================================================
  subroutine par_update_analytical_bcond(vars_of_unk,case,ctime,p_fe_space,caset,t)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine updates Dirichlet boundary conditions in unkno from an analytical solution. !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)          , intent(in)    :: vars_of_unk(:)
    integer(ip)          , intent(in)    :: case
    real(rp)             , intent(in)    :: ctime
    type(par_fe_space_t)  , intent(inout) :: p_fe_space
    integer(ip), optional, intent(in)    :: caset,t

    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )

    ! If fine task call serial subroutine
    if( p_fe_space%p_trian%p_env%am_i_fine_task() ) then
       call fem_update_analytical_bcond( vars_of_unk,case,ctime,p_fe_space%fe_space,caset,t)
    end if

  end subroutine par_update_analytical_bcond
  
  !==================================================================================================
  subroutine par_update_solution_mono(p_vec,p_fe_space,iblock)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a fem_vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(par_vector_t)     , intent(in)    :: p_vec   
    type(par_fe_space_t)  , intent(inout) :: p_fe_space
    integer(ip), optional, intent(in)    :: iblock

    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )

    ! If fine task call serial subroutine
    if( p_fe_space%p_trian%p_env%am_i_fine_task() ) then
       call fem_update_solution(p_vec%f_vector,p_fe_space%fe_space,iblock)
    end if

  end subroutine par_update_solution_mono

  !==================================================================================================
  subroutine par_update_solution_block(blk_p_vec,p_fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the solution from a fem_vector into unkno.                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(par_block_vector_t), intent(in)    :: blk_p_vec   
    type(par_fe_space_t)   , intent(inout) :: p_fe_space
    ! Locals
    integer(ip) :: iblock

    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )

    ! If fine task call serial subroutine
    if( p_fe_space%p_trian%p_env%am_i_fine_task() ) then 

       ! Loop over blocks
       do iblock = 1,blk_p_vec%nblocks

          ! Call monolithic update
          call fem_update_solution(blk_p_vec%blocks(iblock)%f_vector,p_fe_space%fe_space,iblock)

       end do

    end if

  end subroutine par_update_solution_block 

  !==================================================================================================
  subroutine par_update_nonlinear(p_fe_space)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine stores the previous nonlinear solution.                                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(par_fe_space_t), intent(inout) :: p_fe_space

    ! Parallel environment MUST BE already created
    assert ( associated(p_fe_space%p_trian) )
    assert ( p_fe_space%p_trian%p_env%created )

    ! If fine task call serial subroutine
    if( p_fe_space%p_trian%p_env%am_i_fine_task() ) then
       call fem_update_nonlinear(p_fe_space%fe_space)
    end if

  end subroutine par_update_nonlinear

end module par_update_names
    
