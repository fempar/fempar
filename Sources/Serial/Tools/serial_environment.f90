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
module serial_environment_names
use types_names
  use abstract_environment_names
  implicit none
  
  private
  ! Abstract environment
  type, extends(abstract_environment_t) :: serial_environment_t
   contains
     procedure :: info                => serial_environment_info
     procedure :: am_i_fine_task      => serial_environment_am_i_fine_task
     procedure :: bcast               => serial_environment_bcast
     procedure :: first_level_barrier => serial_environment_first_level_barrier
  end type serial_environment_t
  
  public :: serial_environment_t
  
contains

  subroutine serial_environment_first_level_barrier(env) 
    implicit none
    class(serial_environment_t),intent(in)  :: env
  end subroutine serial_environment_first_level_barrier

  subroutine serial_environment_info(env,me,np) 
    implicit none
    class(serial_environment_t),intent(in)  :: env
    integer(ip)              ,intent(out) :: me
    integer(ip)              ,intent(out) :: np
    me = 0
    np = 1 
  end subroutine serial_environment_info
  
  function serial_environment_am_i_fine_task(env) 
    implicit none
    class(serial_environment_t) ,intent(in)  :: env
    logical                                :: serial_environment_am_i_fine_task 
    serial_environment_am_i_fine_task = .true.
  end function serial_environment_am_i_fine_task
  
  subroutine serial_environment_bcast (env, condition)
    implicit none
    class(serial_environment_t) ,intent(in)    :: env
    logical                   ,intent(inout) :: condition
  end subroutine serial_environment_bcast
  
end module serial_environment_names
