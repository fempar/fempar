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
  subroutine evaluate_analytical(case_spatial,case_temporal,ndime,coord,time,values,alpha,beta,gamma, &
       &                         tvar)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the scalar components for an analytical solution.                !
    !----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)          , intent(in)    :: case_spatial,case_temporal,ndime
    real(rp)             , intent(inout) :: values(11)
    real(rp)             , intent(in)    :: coord(ndime),time
    real(rp), optional   , intent(in)    :: alpha,beta,gamma
    integer(ip), optional, intent(in)    :: tvar
    ! Locals
    real(rp) :: spatial_values(10)
    real(rp) :: temporal_values(3)
    
    call evaluate_analytical_spatial(case_spatial,ndime,coord,time,spatial_values,alpha,beta)
    call evaluate_analytical_temporal(case_temporal,time,temporal_values,gamma)
    if(present(tvar)) then       
       values(1:10) = spatial_values(:)*temporal_values(tvar)
    else
       values(1:10) = spatial_values(:)*temporal_values(1)
    end if
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
    real(rp)    :: pi,a,b,c,d,e
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
    case(7)                ! Taylor-Green Vortex, Vx
       a = 1.0_rp!2.0_rp*pi
       if(ndime==2) then
          u=sin(a*x)*cos(a*y)
          dudx = a*cos(a*x)*cos(a*y)
          dudy = -a*sin(a*x)*sin(a*y)
          d2udx = -a*a*u
          d2udy = -a*a*u 
          d2udxy = -a*a*cos(a*x)*sin(a*y)  
       elseif(ndime==3) then
          u=a*cos(x)*sin(y)*sin(z)
          dudx = -a*sin(x)*sin(y)*sin(z)
          dudy = a*cos(x)*cos(y)*sin(z)
          dudz = a*cos(x)*sin(y)*cos(z)
          d2udx = -u
          d2udy = -u
          d2udz = -u
          d2udxy = -cos(a*x)*sin(a*y)
          d2udxz = -a*sin(x)*sin(y)*cos(z)
          d2udyz = a*cos(x)*cos(y)*cos(z)    
       end if
    case(8)                ! Taylor-Green Vortex, Vy
       a = 1.0_rp!2.0_rp*pi
       if(ndime==2) then
          u=-cos(a*x)*sin(a*y)
          dudx = -a*cos(a*x)*cos(a*y)
          dudy = a*sin(a*x)*sin(a*y)
          d2udx = -a*a*u
          d2udy = -a*a*u
          d2udxy = a*a*sin(a*x)*cos(a*y)      
       elseif(ndime==3) then
          u=-a*sin(x)*cos(y)*sin(z)
          dudx = -a*cos(x)*cos(y)*sin(z) 
          dudy = a*sin(x)*sin(y)*sin(z)
          dudz = -a*sin(x)*cos(y)*cos(z)
          d2udx = -u; d2udy = -u; d2udz = -u
          d2udx = -u
          d2udy = -u
          d2udz = -u
          d2udxy = a*cos(x)*sin(y)*sin(z)
          d2udxz = -a*cos(x)*cos(y)*cos(z)
          d2udyz = a*sin(x)*sin(y)*cos(z)
       end if 
    case(9)                ! Channel flow (Re_tau=395), Vx
       check(ndime==3)
       a = 7.764_rp*(395.0_rp**(1/7))    ! C = 7.764*(Re_tau^(1/7))
       b = 0.1_rp*a                      ! eps = 0.1*C
       c = b*2.0_rp*pi/2.0_rp            ! eps*Lx/2
       d = 4.0_rp*pi/(2.0_rp*pi)         ! 4*pi/Lx
       e = 2.0_rp*pi/(2.0_rp*pi/3.0_rp)  ! 2*pi/Lz
       u = a*(1-y**8)+c*sin(pi*y)*cos(d*x)*sin(e*z)
       dudx = -d*c*sin(pi*y)*sin(d*x)*sin(e*z)
       dudy = -8.0_rp*a*y**7+c*pi*cos(pi*y)*cos(d*x)*sin(e*z)
       dudz = e*c*sin(pi*y)*cos(d*x)*cos(e*z)
       d2udx  = -d*d*c*sin(pi*y)*cos(d*x)*sin(e*z)
       d2udxy = -d*c*pi*cos(pi*y)*sin(d*x)*sin(e*z)
       d2udxz = -e*d*c*sin(pi*y)*sin(d*x)*cos(e*z)
       d2udy  = -56.0_rp*a*y**6-c*pi*pi*sin(pi*y)*cos(d*x)*sin(e*z)
       d2udyz = e*c*pi*cos(pi*y)*cos(d*x)*cos(e*z)
       d2udz  = -e*e*c*sin(pi*y)*cos(d*x)*sin(e*z)
    case(10)               ! Channel flow (Re_tau=395), Vy
       check(ndime==3)
       a = 7.764_rp*(395.0_rp**(1/7))    ! C = 7.764*(Re_tau^(1/7))
       b = 0.1_rp*a                      ! eps = 0.1*C
       d = 4.0_rp*pi/(2.0_rp*pi)         ! 4*pi/Lx
       e = 2.0_rp*pi/(2.0_rp*pi/3.0_rp)  ! 2*pi/Lz
       u = -b*(1+cos(pi*y))*sin(d*x)*sin(e*z)
       dudx = -d*b*(1+cos(pi*y))*cos(d*x)*sin(e*z)
       dudy = b*pi*sin(pi*y)*sin(d*x)*sin(e*z)
       dudz = -e*b*(1+cos(pi*y))*sin(d*x)*cos(e*z)
       d2udx  = d*d*b*(1+cos(pi*y))*sin(d*x)*sin(e*z)
       d2udxy = d*b*pi*sin(pi*y)*cos(d*x)*sin(e*z)
       d2udxz = -e*d*b*(1+cos(pi*y))*cos(d*x)*cos(e*z)
       d2udy  = b*pi*pi*cos(pi*y)*sin(d*x)*sin(e*z)
       d2udyz = e*b*pi*sin(pi*y)*sin(d*x)*cos(e*z)
       d2udz  = e*e*b*(1+cos(pi*y))*sin(d*x)*sin(e*z)
    case(11)                ! Channel flow (Re_tau=395), Vz
       check(ndime==3)
       a = 7.764_rp*(395.0_rp**(1/7))    ! C = 7.764*(Re_tau^(1/7))
       b = 0.1_rp*a                      ! eps = 0.1*C
       c = b*(2.0_rp*pi/3.0_rp)/2.0_rp   ! eps*Lz/2
       d = 4.0_rp*pi/(2.0_rp*pi)         ! 4*pi/Lx
       e = 2.0_rp*pi/(2.0_rp*pi/3.0_rp)  ! 2*pi/Lz
       u = -c*sin(pi*y)*sin(d*x)*cos(e*z)
       dudx = -d*c*sin(pi*y)*cos(d*x)*cos(e*z)
       dudy = -c*pi*cos(pi*y)*sin(d*x)*cos(e*z)
       dudz = e*c*sin(pi*y)*sin(d*x)*sin(e*z)
       d2udx  = d*d*c*sin(pi*y)*sin(d*x)*cos(e*z)
       d2udxy = -d*c*pi*cos(pi*y)*cos(d*x)*cos(e*z)
       d2udxz = e*d*c*sin(pi*y)*cos(d*x)*sin(e*z)
       d2udy  = c*pi*pi*sin(pi*y)*sin(d*x)*cos(e*z)
       d2udyz = e*c*pi*cos(pi*y)*sin(d*x)*sin(e*z)
       d2udz  = e*e*c*sin(pi*y)*sin(d*x)*cos(e*z)
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
