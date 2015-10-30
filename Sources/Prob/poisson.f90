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
module poisson_names
  use types_names
  use memor_names
  use problem_names
  use finite_element_names
  use element_fields_names
  use element_tools_names
  use analytical_function_names
  implicit none
# include "debug.i90"
  private 

  type, extends(physical_problem_t) :: poisson_problem_t
     integer(ip) ::   & 
          case_space, & ! Exact solution (space function)
          case_tempo    ! Exact solution (temporal function)
     real(rp) ::         &
          diffu            ! Diffusion
     
   contains
     procedure :: create => poisson_create
     procedure :: free => poisson_free
  end type poisson_problem_t

  public :: poisson_problem_t, poisson_analytical_force

contains

  !=================================================================================================
  subroutine poisson_create( prob, ndime)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(poisson_problem_t), intent(inout) :: prob
    integer(ip)       , intent(in)  :: ndime

    integer(ip) :: i

    ! Fill default problem data
    prob%ndime = ndime
    !prob%ntens = prob%ndime*(prob%ndime+1)/2 ! Number of tensor components 
    prob%nunks = 1     
    prob%nvars = 1     

    call memalloc(prob%nunks,prob%vars_of_unk,__FILE__,__LINE__)
    prob%vars_of_unk(1) = 1

    ! Problem variables
    prob%diffu  = 1.0_rp  ! Diffusion

    ! Analytical field variables
    prob%case_space   = 0      ! Exact solution (in space)
    prob%case_tempo   = 0      ! Exact solution (in time)

  end subroutine poisson_create

  !=================================================================================================
  subroutine poisson_free( prob )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(poisson_problem_t), intent(inout) :: prob

    call memfree ( prob%vars_of_unk,__FILE__,__LINE__)

  end subroutine poisson_free

  !==================================================================================================
  subroutine poisson_analytical_force(physics,finite_element,ctime,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine computes the elemental force needed to impose an analytical solution        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(poisson_problem_t)   , intent(inout) :: physics
    type(finite_element_t), intent(in)    :: finite_element
    real(rp)              , intent(in)    :: ctime
    type(scalar_t)        , intent(inout) :: force
    ! Locals
    integer(ip)    :: igaus,ivars,idime,conve
    real(rp)       :: gpvno
    real(rp)       :: params(11,physics%nvars)
    
    ! Loop over gauss points
    do igaus=1,finite_element%integ(1)%p%quad%ngaus
       
       ! Evaluate unknowns and derivatives
       params = 0.0_rp
       do ivars=1,physics%nvars
          call evaluate_analytical(finite_element%p_analytical_code%a(ivars,1),                  &
               &                   finite_element%p_analytical_code%a(ivars,2),                  &
               &                   physics%ndime,finite_element%integ(1)%p%femap%clocs(:,igaus), &
               &                   ctime,params(:,ivars))
       end do
          
       ! Evaluate force
       call poisson_force(physics%ndime,physics%diffu,params,force%a(igaus))

    end do
  end subroutine poisson_analytical_force

  !==================================================================================================
  subroutine  poisson_force(ndime,diffu,params,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the elemental force needed to impose an analytical solution       !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip), intent(in)    :: ndime
    real(rp)   , intent(in)    :: diffu,params(11,1)
    real(rp)   , intent(inout) :: force

    real(rp)    :: d2udx,d2udy,d2udz

    d2udx = params(5,1)
    d2udy = params(6,1)
    d2udz = params(7,1)

    force =  diffu*(d2udx+d2udy+d2udz)    

  end subroutine poisson_force

end module poisson_names

