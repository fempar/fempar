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
# include "debug.i90"
module cdr_names
use types_names
use memor_names
  use problem_names
  implicit none
  private 

  type, extends(physical_problem) :: cdr_problem_t
     integer(ip) ::   & 
          kfl_conv,   & ! Flag for enabling advection
          kfl_tder,   & ! Flag for time derivative computation
          case_space, & ! Exact solution (space function)
          case_tempo    ! Exact solution (temporal function)
     real(rp) ::         &
          react,         & ! Reaction
          diffu            ! Diffusion
   contains
     procedure :: create => cdr_create
     procedure :: free => cdr_free
  end type cdr_problem_t

  public :: cdr_problem_t

contains

  !=================================================================================================
  subroutine cdr_create( prob, ndime)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(cdr_problem_t), intent(inout) :: prob
    integer(ip)       , intent(in)  :: ndime

    integer(ip) :: i

    ! Fill default problem data
    prob%ndime = ndime
    !prob%ntens = prob%ndime*(prob%ndime+1)/2 ! Number of tensor components 
    prob%nunks = 1     
    prob%nvars = 1     

    call memalloc(prob%nunks,prob%vars_of_unk,__FILE__,__LINE__)
    prob%vars_of_unk(1) = 1

    ! Flags
    prob%kfl_conv = 1 ! Enabling advection
    prob%kfl_tder = 0 ! Time derivative not computed 

    ! Problem variables
    prob%react  = 0.0_rp  ! Reaction
    prob%diffu  = 1.0_rp  ! Diffusion

    ! Analytical field variables
    prob%case_space   = 0      ! Exact solution (in space)
    prob%case_tempo   = 0      ! Exact solution (in time)

  end subroutine cdr_create



  !=================================================================================================
  subroutine cdr_free( prob )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(cdr_problem_t), intent(inout) :: prob

    call memfree ( prob%vars_of_unk,__FILE__,__LINE__)

  end subroutine cdr_free


end module cdr_names
