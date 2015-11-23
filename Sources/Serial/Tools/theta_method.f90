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
module theta_method_names
  use types_names
  use memor_names
  use time_integration_names
  implicit none
# include "debug.i90"
  private

  type, extends(time_integration_t) :: theta_method_t
     real(rp)    :: theta           ! Theta method to be used
     real(rp)    :: real_time       ! Time in wich u_{n+1} is computed (ctime is for u_{n+theta}
     real(rp)    :: time_step       ! The time_Step value
     integer(ip) :: number_of_steps ! Total amount of steps
     logical     :: finished        ! Boolean to decide whether to finish the time integration
   contains
     procedure   :: create     => theta_method_create
     procedure   :: initialize => theta_method_initialize
     procedure   :: update     => theta_method_update
     procedure   :: print      => theta_method_print
     procedure   :: update_solution => theta_method_update_solution
  end type theta_method_t

  ! Type
  public :: theta_method_t
  
contains

  !==================================================================================================
  subroutine theta_method_create(theta_method)
    implicit none
    class(theta_method_t), intent(out) :: theta_method

    theta_method%theta           = 0.5_rp
    theta_method%real_time       = 0.0_rp
    theta_method%time_step       = 0.0_rp
    theta_method%number_of_steps = 0
    theta_method%finished        = .false.
  end subroutine theta_method_create

  !==================================================================================================
  subroutine theta_method_initialize(theta_method)
    implicit none
    class(theta_method_t), intent(inout) :: theta_method

    ! First step
    theta_method%istep = 1

    if (theta_method%time_step == 0.0_rp) then
       ! Steady Problem
       theta_method%dtinv = 0.0_rp
       ! Calculate the total amount of steps
       theta_method%number_of_steps = 1 
    else
       ! Transient problem
       theta_method%dtinv = 1.0_rp/(theta_method%time_step*theta_method%theta)
       ! Step id
       theta_method%istep = 0
       ! Calculate the total amount of steps
       theta_method%number_of_steps = int((theta_method%ftime-theta_method%itime)/                  &
            &                              theta_method%time_step)
       ! Seek the time values for the first step
       theta_method%real_time = theta_method%itime
       theta_method%ctime = theta_method%itime + (theta_method%theta-1.0_rp)*theta_method%time_step
    end if

    ! At least we will perform 1 step
    theta_method%finished = .false.

   
  end subroutine theta_method_initialize

  !==================================================================================================
  subroutine theta_method_update(theta_method)
    implicit none
    class(theta_method_t), intent(inout) :: theta_method
    
    if (theta_method%time_step == 0.0_rp) then
       theta_method%finished = .true.
    else
       ! Transient problem
       theta_method%istep = theta_method%istep + 1
       ! Seek the time values for the first step
       theta_method%real_time = theta_method%real_time + theta_method%time_step
       theta_method%ctime = theta_method%ctime + theta_method%time_step
       ! Check if finished
       if (theta_method%istep > theta_method%number_of_steps) theta_method%finished = .true.
    end if
  end subroutine theta_method_update

  !==================================================================================================
  subroutine theta_method_print(theta_method,luout)
    implicit none
    class(theta_method_t), intent(in) :: theta_method
    integer(ip)          , intent(in) :: luout

    write(luout,*)         '========================================================================'
    write(luout,'(a10,i3,a1,i3)') 'Time step ', theta_method%istep ,'/',theta_method%number_of_steps

  end subroutine theta_method_print

  !=================================================================================================
  subroutine  theta_method_update_solution(this,current_unkno,prev_step_unkno)
    implicit none
    class(theta_method_t), intent(in)    :: this
    real(rp)             , intent(inout) :: current_unkno(:)
    real(rp)             , intent(inout) :: prev_step_unkno(:)

    real(rp) :: scale_prev, scale_current

    if (this%dtinv > 0.0_rp) then
       scale_prev    = (this%theta-1.0_rp)/this%theta
       scale_current = 1.0_rp/this%theta
       
       prev_step_unkno = scale_prev * prev_step_unkno + scale_current * current_unkno
       current_unkno   = prev_step_unkno
    end if
  end subroutine theta_method_update_solution
  
end module theta_method_names
