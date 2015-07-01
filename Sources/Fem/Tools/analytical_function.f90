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
module analytical_function_names
  use types_names
  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: evaluate_analytical, evaluate_analytical_spatial, evaluate_analytical_temporal

contains

  !==================================================================================================
  subroutine evaluate_analytical(case_spatial,case_temporal,ndime,coord,time,values,alpha,beta,gamma)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the scalar components for an analytical solution.                !
    !----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)       , intent(in)    :: case_spatial,case_temporal,ndime
    real(rp)          , intent(inout) :: values(11)
    real(rp)          , intent(in)    :: coord(ndime),time
    real(rp), optional, intent(in)    :: alpha,beta,gamma
    ! Locals
    real(rp) :: spatial_values(10)
    real(rp) :: temporal_values(3)
    
    call evaluate_analytical_spatial(case_spatial,ndime,coord,time,spatial_values,alpha,beta)
    call evaluate_analytical_temporal(case_temporal,time,temporal_values,gamma)
    values(1:10) = spatial_values(:)*temporal_values(1)
    values(11) = spatial_values(1)*temporal_values(2)  ! temporal derivative

  end subroutine evaluate_analytical

  !==================================================================================================
  subroutine evaluate_analytical_spatial(case,ndime,coord,time,values,alpha,beta)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the scalar components for an analytical spatial solution.        !
    !----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)      , intent(in)    :: case,ndime
    real(rp),optional, intent(in)    :: time
    real(rp)         , intent(inout) :: values(10)
    real(rp)         , intent(in)    :: coord(ndime)
    real(rp),optional, intent(in)    :: alpha,beta
    ! Local variables
    real(rp)    :: pi
    real(rp)    :: x,y,z,u,dudx,dudy,dudz
    real(rp)    :: d2udx,d2udy,d2udz,d2udxy,d2udxz,d2udyz

    u=0.0_rp; dudx = 0.0_rp; dudy = 0.0_rp; dudz = 0.0_rp
    d2udx = 0.0_rp; d2udy = 0.0_rp; d2udz = 0.0_rp
    d2udxy = 0.0_rp; d2udxz = 0.0_rp; d2udyz = 0.0_rp
    pi = 4.0_rp*atan(1.0_rp)
    
    ! Point coordinates
    x = coord(1)
    y = coord(2)
    if(ndime==3) then
       z = coord(3)
    else
       z = 0.0_rp
    end if

    select case(case)
    case(1)                ! u = x
       u = x
       dudx = 1.0_rp
    case(2)                ! u = -y
       u = -y       
       dudy = -1.0_rp
    case(3)                ! u = x+y
       u = 1.0_rp*(x+y)
       dudx = 1.0_rp
       dudy = 1.0_rp 
    case(4)                ! u = x^2*y
       u = x*x*y
       dudx = 2.0_rp*x*y
       dudy = x*x
       d2udx = 2.0_rp*y
       d2udxy = 2.0_rp*x
    case(5)                ! u = -y^2*x
       u = -x*y*y
       dudx = -y*y
       dudy = -2.0_rp*x*y
       d2udy = -2.0_rp*x
       d2udxy = -2.0_rp*y  
    end select
    
    ! Component assignment
    values(1) = u
    values(2) = dudx; values(3) = dudy; values(4) = dudz
    values(5) = d2udx; values(6) = d2udy; values(7) = d2udz
    values(8) = d2udxy; values(9) = d2udxz; values(10) = d2udyz

  end subroutine evaluate_analytical_spatial

  !==================================================================================================
  subroutine evaluate_analytical_temporal(case,t,values,gamma)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the scalar components for an analytical temporal solution.       !
    !----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)       , intent(in)    :: case
    real(rp)          , intent(inout) :: values(3)
    real(rp)          , intent(in)    :: t
    real(rp), optional, intent(in)    :: gamma
    ! Local variables
    real(rp)    :: u,dudt,d2udt,a,b,gamma_

    u=1.0_rp; dudt=0.0_rp; d2udt=0.0_rp
    
    gamma_ = 1.0_rp
    if(present(gamma)) gamma_ = gamma

    select case(case)
    case(1)
       u = t
       dudt = 1.0_rp
    case(2)
       u = t*t
       dudt = 2.0_rp*t
       d2udt = 2.0_rp
    case(3)
       a = pi/10.0_rp
       b = 1.0_rp/25.0_rp
       u = sin(a*t)*exp(b*t)
       dudt = a*cos(a*t)*exp(b*t) + b*sin(a*t)*exp(b*t)
       d2udt = (b*b-a*a)*sin(a*t)*exp(b*t) + 2.0_rp*a*b*cos(a*t)*exp(b*t)
    case(4)
       a = 2.0_rp*pi*pi*gamma_
       u = exp(-a*t)
       dudt = -a*exp(-a*t)
       d2udt = a*a*exp(-a*t)
    case(5)
       a = 4.0_rp*pi*pi*gamma_
       u = exp(-a*t)
       dudt = -a*exp(-a*t)
       d2udt = a*a*exp(-a*t)
    case(6)
       u = t*t*t
       dudt = 3.0_rp*t*t
       d2udt = 6.0_rp*t
    case(7)
       u = t*t*t*t
       dudt = 4.0_rp*t*t*t
       d2udt = 12.0_rp*t*t
    end select

    ! Component assignment
    values(1) = u
    values(2) = dudt
    values(3) = d2udt

  end subroutine evaluate_analytical_temporal
  
end module analytical_function_names
