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
module cdr_names
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

  type, extends(physical_problem_t) :: cdr_problem_t
     integer(ip) ::   & 
          kfl_conv,   & ! Flag for enabling advection
          kfl_tder,   & ! Flag for time derivative computation
          kfl_react,  & ! Flag for analytical reaction
          case_space, & ! Exact solution (space function)
          case_tempo    ! Exact solution (temporal function)
     real(rp) ::         &
          react,         & ! Reaction
          diffu            ! Diffusion
     real(rp), allocatable :: &
          convect(:)       ! Convection
     
   contains
     procedure :: create => cdr_create
     procedure :: free => cdr_free
  end type cdr_problem_t

  public :: cdr_problem_t, cdr_analytical_force, cdr_analytical_reaction,  cdr_analytical_convection

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
    prob%kfl_conv = 0 ! Enabling advection
    prob%kfl_tder = 0 ! Time derivative not computed 
    prob%kfl_react = 0 ! Non analytical reaction

    ! Problem variables
    prob%react  = 0.0_rp  ! Reaction
    prob%diffu  = 1.0_rp  ! Diffusion
    call memalloc(ndime,prob%convect,__FILE__,__LINE__)
    prob%convect = 0

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
    call memfree ( prob%convect,    __FILE__,__LINE__)

  end subroutine cdr_free

  !==================================================================================================
  subroutine cdr_analytical_force(physics,finite_element,ctime,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine computes the elemental force needed to impose an analytical solution        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(cdr_problem_t)   , intent(inout) :: physics
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

       ! Compute analytical reaction
       if(physics%kfl_react>0) then
          call cdr_analytical_reaction(physics%ndime,finite_element%integ(1)%p%femap%clocs(:,igaus), &
               &                       params(1,1),physics%kfl_react,physics%react)
       end if

       ! Compute analytical convection
       if (physics%kfl_conv>0) then
          call cdr_analytical_convection(physics%ndime,                                             &
               &                         finite_element%integ(1)%p%femap%clocs(:,igaus),            &
               &                         physics%kfl_conv,ctime,physics%convect)
       end if
          
       ! Evaluate force
       call cdr_force(physics%ndime,physics%diffu,physics%react,physics%convect,params,force%a(igaus))

    end do
  end subroutine cdr_analytical_force

  !==================================================================================================
  subroutine  cdr_force(ndime,diffu,react,convect,params,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the elemental force needed to impose an analytical solution       !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip), intent(in)    :: ndime
    real(rp)   , intent(in)    :: diffu,react,params(11,1)
    real(rp)   , intent(in)    :: convect(ndime)
    real(rp)   , intent(inout) :: force

    real(rp)    :: u,d2udx,d2udy,d2udz,d2udxy,d2udxz,d2udyz
    real(rp)    :: dudx,dudy,dudz,dudt

    !
    dudt = params(11,1)
    !
    u = params(1,1)
    !
    dudx = params(2,1)
    dudy = params(3,1)
    dudz = params(4,1)
    !
    d2udx = params(5,1)
    d2udy = params(6,1)
    d2udz = params(7,1)
    !
    d2udxy = params(8,1)
    d2udxz = params(9,1)
    d2udyz = params(10,1)

    force =  dudt - diffu*(d2udx+d2udy+d2udz) + react*u 
    force = force + convect(1)*dudx + convect(2)*dudy
    if (ndime == 3) force = force + convect(3)*dudz
   

  end subroutine cdr_force

  !==================================================================================================
  subroutine  cdr_analytical_reaction(ndime,clocs,unkno,kfl_react,react)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine computes the analytical expression of the reaction constant                 !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip), intent(in)    :: ndime,kfl_react
    real(rp)   , intent(in)    :: clocs(ndime),unkno
    real(rp)   , intent(inout) :: react

    react = 0.0_rp
    
    if(kfl_react==1) then
       react = 1.0_rp
    elseif(kfl_react==2) then
       react = clocs(1) + clocs(2)
    elseif(kfl_react==3) then
       react = unkno
    end if

  end subroutine cdr_analytical_reaction

 !===================================================================================================
  subroutine cdr_analytical_convection(ndime,coord,kfl_conv,ctime,convect, divergence_convect)
    !----------------------------------------------------------------------------------------
    ! This subroutine gives the value of the convection on coord depending on the kfl_conv
    !----------------------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)            :: kfl_conv,ndime
    real(rp)   , intent(in)            :: ctime
    real(rp)   , intent(in)            :: coord(ndime)
    real(rp)   , intent(out)           :: convect(ndime)
    real(rp)   , intent(out), optional :: divergence_convect

    real(rp)                 :: x,y,pi=3.14159265, divfi

    convect=0.0_rp
    divfi = 0.0_rp
    if(kfl_conv== 1) then
       ! b = (0,1)
       convect(1)= 0.0_rp 
       convect(2)= 1.0_rp
       divfi = 0.0_rp
    else if(kfl_conv== 3) then
       x = coord(1)-0.5
       y = coord(2)-0.5
       convect(1)=  2*pi*y
       convect(2)= -2*pi*x 
       divfi = 0.0_rp 
    else if(kfl_conv== 4) then
       x = coord(1)-0.5
       y = coord(2)-0.5
       convect(1)= -2*pi*y
       convect(2)= 2*pi*x 
       divfi = 0.0_rp 
    else if(kfl_conv== 5) then
       convect(1)= 0.0
       convect(2)= 0.0
       convect(3) = 1.0
       divfi = 0.0_rp 
    else if(kfl_conv== 6) then
       x = coord(1)
       y = coord(2)
       convect(1)= -2*pi*y
       convect(2)= 2*pi*x 
       divfi = 0.0_rp 
    else if(kfl_conv== 7) then
       convect(1)= -cos(55.0_rp/360.0_rp*2_rp*pi)
       convect(2)=   -sin(55.0_rp/360.0_rp*2_rp*pi)
    else if(kfl_conv== 8) then
       convect(1)= cos(-pi/3.0_rp)
       convect(2)= sin(-pi/3.0_rp)
       divfi = 0.0_rp 
    else if(kfl_conv== 9) then
       convect(1)= sin(pi/3.0_rp)
       convect(2)= cos(pi/3.0_rp)
       divfi = 0.0_rp 
    else if (kfl_conv == 10) then
       ! b = (1,0)
       convect(1)= 1.0_rp 
       convect(2)= 0.0_rp
       divfi = 0.0_rp
    else if(kfl_conv== 11) then
       ! b = (1,1)
       convect(1)= 1.0_rp 
       convect(2)= 1.0_rp
       divfi = 0.0_rp 
    else if(kfl_conv== 12) then
       ! b = 1/sqrt(5)*(1,-2)
       convect(1)= 1.0_rp/sqrt(5.0_rp) 
       convect(2)= -2.0_rp/sqrt(5.0_rp) 
       divfi = 0.0_rp 
    else if(kfl_conv== 13) then
       convect(1)= cos(pi/3.0_rp)
       convect(2)= sin(pi/3.0_rp)
       divfi = 0.0_rp 
    else if(kfl_conv== 14) then
       convect(1)= cos(7.0_rp*pi/20.0_rp)
       convect(2)= sin(7.0_rp*pi/20.0_rp)
       divfi = 0.0_rp 
    else if(kfl_conv== 15) then
       convect(1)= cos(pi/60.0_rp)
       convect(2)= sin(pi/60.0_rp)
       divfi = 0.0_rp 
    else if(kfl_conv== 16) then
       x = coord(1)
       y = coord(2)
       convect(1)= y
       convect(2)= -x 
       divfi = 0.0_rp 
    else if(kfl_conv== 17) then
       x = coord(1)
       y = coord(2)
       convect(1)= sin(pi*x)
       convect(2)= sin(pi*y) 
       divfi = pi*(cos(pi*x)+cos(pi*y)) 
    else if(kfl_conv== 18) then
       x = coord(1)
       y = coord(2)
       convect(1) = 2.0_rp*y*(1.0_rp-x**2)
       convect(2) = -2.0_rp*x*(1.0_rp-y**2)
       divfi = 0.0_rp
    else if(kfl_conv== 100) then
       ! b = (10,0)
       convect(1)= 10.0_rp 
       convect(2)= 0.0_rp
       divfi = 0.0_rp 
    end if

    if (present( divergence_convect))  divergence_convect = divfi
  end subroutine cdr_analytical_convection

end module cdr_names

