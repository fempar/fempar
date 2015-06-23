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
module abstract_environment_names
  use types
  implicit none

  private
  ! Abstract environment
  type, abstract :: abstract_environment
   contains
     procedure (info_interface)                , deferred  :: info
     procedure (am_i_fine_task_interface)      , deferred  :: am_i_fine_task
     procedure (bcast_interface)               , deferred  :: bcast
     procedure (first_level_barrier_interface) , deferred  :: first_level_barrier
  end type abstract_environment

  ! Abstract interfaces
  abstract interface
     subroutine info_interface(env,me,np) 
       import :: abstract_environment, ip
       implicit none
       class(abstract_environment),intent(in)  :: env
       integer(ip)                ,intent(out) :: me
       integer(ip)                ,intent(out) :: np
     end subroutine info_interface

     function am_i_fine_task_interface(env) 
       import :: abstract_environment, ip
       implicit none
       class(abstract_environment) ,intent(in)  :: env
       logical                                  :: am_i_fine_task_interface 
     end function am_i_fine_task_interface

     subroutine bcast_interface (env, condition)
       import :: abstract_environment
       implicit none
       class(abstract_environment) ,intent(in)    :: env
       logical                     ,intent(inout) :: condition
     end subroutine bcast_interface

     subroutine first_level_barrier_interface (env)
       import :: abstract_environment
       implicit none
       class(abstract_environment) ,intent(in)    :: env
     end subroutine first_level_barrier_interface
  end interface

  public :: abstract_environment

end module abstract_environment_names
