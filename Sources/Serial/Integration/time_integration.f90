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
module time_integration_names
  use types_names
  implicit none
# include "debug.i90"
  private
  
  type :: time_integration_t
     real(rp)    :: ctime = 0.0_rp  ! Current time
     real(rp)    :: ftime = 0.0_rp  ! Final time
     real(rp)    :: itime = 0.0_rp  ! Initial time
     real(rp)    :: dtinv = 0.0_rp  ! Inverse of the time step
     integer(ip) :: istep = 1       ! Step id
   contains
     procedure :: update_solution => time_integration_update_solution
  end type time_integration_t
  
  ! Types
  public :: time_integration_t

contains

  !=================================================================================================
  subroutine  time_integration_update_solution(this,current_unkno,prev_step_unkno)
    implicit none
    class(time_integration_t), intent(in)    :: this
    real(rp)                 , intent(inout) :: current_unkno(:)
    real(rp)                 , intent(inout) :: prev_step_unkno(:)

    prev_step_unkno = current_unkno
  end subroutine time_integration_update_solution

end module time_integration_names
