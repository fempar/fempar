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
module analytical_names
use types_names
  implicit none
  private

  public :: analytical_field

contains

  !=============================================================================
  ! Analytical (exact) solutions for vectorial & scalar fields
  ! fcomp contains the several components & derivatives
  !=============================================================================
  subroutine analytical_field(ncase,ndime,coord,time,fpara,alpha)
    implicit none
    integer(ip), intent(in)  :: ncase,ndime
    real(rp)   , intent(inout)  :: fpara(:)
    real(rp)   , intent(in)  :: coord(ndime),time
    real(rp), optional, intent(in) :: alpha

    real(rp) :: alpha_

    if(present(alpha)) then
       alpha_ = alpha
    else
       alpha_ = 0.0_rp
    end if
    
    if (size(fpara) == 10) then
       call analytical_scalar_field(ncase,ndime,coord,fpara,time) 
    elseif (size(fpara) == 30) then
       call analytical_vector_field(ncase,ndime,coord,fpara,alpha_)
    elseif (size(fpara) == 3) then
       call analytical_temporal_field(ncase,time,fpara,alpha_)   
    end if
      
  end subroutine analytical_field

  !=============================================================================
  subroutine analytical_vector_field(ncase,ndime,coord,fpara,alpha)
    implicit none
    integer(ip), intent(in)  :: ncase,ndime
    real(rp)   , intent(inout)  :: fpara(30)
    real(rp)   , intent(in)  :: coord(ndime),alpha
    ! Local variables
    real(rp)    :: x,y,z,u,v,w,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp)    :: d2udx,d2udy,d2udz,d2udxy,d2udxz,d2udyz
    real(rp)    :: d2vdx,d2vdy,d2vdz,d2vdxy,d2vdxz,d2vdyz
    real(rp)    :: d2wdx,d2wdy,d2wdz,d2wdxy,d2wdxz,d2wdyz
    real(rp)    :: a,b,c
    real(rp)    :: lambda,r,n,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,expd5,expd6,expp1,expp2,expp3
    real(rp)    :: Ha,db,visco,condu,dpdx,Vhunt(12),K
    real(rp)    :: cneg, uo, L, yy, zz

    u=0.0_rp; v=0.0_rp; w=0.0_rp
    dudx = 0.0_rp; dudy = 0.0_rp; dudz = 0.0_rp
    dvdx = 0.0_rp; dvdy = 0.0_rp; dvdz = 0.0_rp
    dwdx = 0.0_rp; dwdy = 0.0_rp; dwdz = 0.0_rp
    d2udx = 0.0_rp; d2udy = 0.0_rp; d2udz = 0.0_rp
    d2udxy = 0.0_rp; d2udxz = 0.0_rp; d2udyz = 0.0_rp
    d2vdx = 0.0_rp; d2vdy = 0.0_rp; d2vdz = 0.0_rp
    d2vdxy = 0.0_rp; d2vdxz = 0.0_rp; d2vdyz = 0.0_rp
    d2wdx = 0.0_rp; d2wdy = 0.0_rp; d2wdz = 0.0_rp
    d2wdxy = 0.0_rp; d2wdxz = 0.0_rp; d2wdyz = 0.0_rp

    ! Point coordinates
    x = coord(1)
    y = coord(2)
    if(ndime==3) then
       z = coord(3)
    else
       z = 0.0_rp
    end if

    if(ncase==1) then          ! u=(x,-y)
       u = x
       v = -y
       dudx = 1.0_rp
       dvdy = -1.0_rp
    elseif(ncase==2) then      
       if(ndime==2) then       ! u=(x^2*y, -y^2*x)
          u = x*x*y
          v = -x*y*y
          dudx = 2.0_rp*x*y
          dudy = x*x
          dvdx = -y*y
          dvdy = -2.0_rp*x*y
          !
          d2udx = 2.0_rp*y
          d2udxy = 2.0_rp*x
          d2vdy = -2.0_rp*x
          d2vdxy = -2.0_rp*y  
       elseif(ndime==3) then   ! u=(1/2*x^2*y*z, -y^2*x*z, 1/2*x*y*z^2)
          u = 0.5_rp*x*x*y*z
          v = -x*y*y*z
          w = 0.5_rp*x*y*z*z
          dudx = x*y*z
          dudy = 0.5_rp*x*x*z
          dudz = 0.5_rp*x*x*y

          dvdx = -y*y*z
          dvdy = -2.0_rp*x*y*z
          dvdz = -y*y*x

          dwdx = 0.5_rp*y*z*z
          dwdy = 0.5_rp*x*z*z
          dwdz = x*y*z
          !
          d2udx = y*z
          d2udxy = x*z
          d2udxz = x*y
          d2udyz = 0.5_rp*x*x
          !
          d2vdxy = -2.0_rp*y*z
          d2vdxz = -y*y
          d2vdy = -2.0_rp*x*z
          d2vdyz = -2.0_rp*x*y
          !
          d2wdxy = 0.5_rp*z*z
          d2wdxz = y*z
          d2wdyz = x*z
          d2wdz = x*y
       end if
    elseif(ncase==3) then      ! u=(x^2, -y*2*x)
       u = x*x
       v = -x*2.0_rp*y
       dudx = 2.0_rp*x
       dvdx = -2.0_rp*y
       dvdy = -2.0_rp*x
       !
       d2udx = 2.0_rp
       d2vdxy = -2.0_rp
    elseif(ncase==4) then      ! u=(sin(2*pi*x),-2*pi*cos(2*pi*x)*y)
       u = sin(2.0_rp*pi*x)
       v = -2.0_rp*pi*cos(2.0_rp*pi*x)*y
       dudx = 2.0_rp*pi*cos(2.0_rp*pi*x)
       dvdx = 4.0_rp*pi*pi*sin(2.0_rp*pi*x)*y
       dvdy = -2.0_rp*pi*cos(2.0_rp*pi*x)
       !
       d2udx = -4.0_rp*pi*pi*sin(2.0_rp*pi*x)
       d2vdx = 8.0_rp*pi*pi*pi*cos(2.0_rp*pi*x)*y
       d2vdxy = 4.0_rp*pi*pi*sin(2.0_rp*pi*x)    
    elseif(ncase==5) then      ! u=(1,1)
       u = 1.0_rp
       v = 1.0_rp
    elseif(ncase==6) then      ! u=(y,2x)
       u = y
       v = 2.0_rp*x
       dudy = 1.0_rp
       dvdx = 2.0_rp
    elseif(ncase==7) then      ! u=(x,-2y,z)
       u = x
       v = -2.0_rp*y
       w = z
       dudx = 1.0_rp
       dvdy = -2.0_rp
       dwdz = 1.0_rp
    elseif(ncase==8) then      
       if(ndime==2) then       ! u=(0.01*cos(16*pi*y),-0.01*cos(16*pi*x))
          a = 2.0_rp*pi
          b = 64.0_rp
          c = 63.0_rp
          u = sin(a*(x+y))!+0.005_rp*1.0_rp/b*cos(a*(b*x+c*y))   
          v = -sin(a*(x+y))!-0.005_rp*1.0_rp/c*cos(a*(b*x+c*y))
          dudx = a*(cos(a*(x+y)))!-0.005_rp*sin(a*(b*x+c*y)))
          dudy = a*(cos(a*(x+y)))!-0.005_rp*c/b*sin(a*(b*x+c*y)))
          dvdx = a*(-cos(a*(x+y)))!+0.005_rp*b/c*sin(a*(b*x+c*y)))
          dvdy = a*(-cos(a*(x+y)))!+0.005_rp*sin(a*(b*x+c*y)))
       else if(ndime==3) then  
          ! u=1/(20*sqrt(3))*(cos(2*pi*(x+y-z)),cos(2*pi*(x+y-z)),2*cos(2*pi*(x+y-z)))
          a = 1.0_rp/(20.0_rp*sqrt(3.0_rp))
          b = 2.0_rp*pi*(x+y-z)
          u = a*cos(b)
          v = a*cos(b)
          w = 2.0_rp*a*cos(b)
          !
          dudx = -2.0_rp*pi*a*sin(b); dvdx = dudx; dwdx = 2.0_rp*dudx
          dudy = -2.0_rp*pi*a*sin(b); dvdy = dudy; dwdy = 2.0_rp*dudy
          dudz = 2.0_rp*pi*a*sin(b); dvdz = dudz; dwdz = 2.0_rp*dudz
          !
          d2udx = -(2.0_rp*pi)**2*u; d2vdx = d2udx; d2wdx = 2.0_rp*d2udx
          d2udy = -(2.0_rp*pi)**2*u; d2vdy = d2udy; d2wdy = 2.0_rp*d2udy
          d2udz = -(2.0_rp*pi)**2*u; d2vdz = d2udz; d2wdz = 2.0_rp*d2udz
          !
          d2udxy = d2udx; d2udxz = -d2udx; d2udyz  = -d2udy
          d2vdxy = d2vdx; d2vdxz = -d2vdx; d2vdyz  = -d2vdy
          d2wdxy = d2wdx; d2wdxz = -d2wdx; d2wdyz  = -d2wdy
       end if
    elseif(ncase==9) then      ! u=(z,-2x,y)
       u = z
       v = -2.0_rp*x
       w = y
       dudz = 1.0_rp
       dvdx = -2.0_rp
       dwdy = 1.0_rp
    elseif(ncase==10) then      ! u=(x^2, -y*2*x)
       u = x*x
       v = -2.0_rp*x*y
       dudx = 2.0_rp*x
       dvdx = -2.0_rp*y
       dvdy = -2.0_rp*x
       !
       d2udx = 2.0_rp
       d2vdxy = -2.0_rp
    elseif(ncase==11) then      ! u=(x^2*y^2, -y^2*x^2)
       u = x*x*y*y
       v = -x*y*y*x
       dudx = 2.0_rp*y*y*x
       dudy = 2.0_rp*x*x*y
       dvdx = -2.0_rp*y*y*x
       dvdy = -2.0_rp*x*y*x
       !
       d2udx = 2.0_rp*y*y
       d2udxy = 4.0_rp*x*y
       d2udy = 2.0_rp*x*x
       d2vdx = -2.0_rp*y*y
       d2vdy = -2.0_rp*x*x
       d2vdxy = -4.0_rp*y*x 
    elseif(ncase==12) then      ! Taylor-Green vortex
       if(ndime==2) then
          a = 1.0_rp!2.0_rp*pi
          u=sin(a*x)*cos(a*y)
          v=-cos(a*x)*sin(a*y)
          !
          dudx = a*cos(a*x)*cos(a*y)
          dudy = -a*sin(a*x)*sin(a*y)
          dvdx = -dudy; dvdy = -dudx
          !
          d2udx = -a*a*u; d2udy = -a*a*u; 
          d2vdx = -a*a*v; d2vdy = -a*a*v;
          d2udxy = a*a*v; d2vdxy = a*a*u      
       elseif(ndime==3) then
          a = 1.0_rp!/2.0_rp*pi
          u=a*cos(x)*sin(y)*sin(z)
          v=-a*sin(x)*cos(y)*sin(z)
          w=0.0_rp
          !
          dudx = -a*sin(x)*sin(y)*sin(z)
          dudy = a*cos(x)*cos(y)*sin(z)
          dudz = a*cos(x)*sin(y)*cos(z)
          dvdx = -dudy; dvdy = -dudx
          dvdz = -a*sin(x)*cos(y)*cos(z)
          !
          d2udx = -u; d2udy = -u; d2udz = -u
          d2vdx = -v; d2vdy = -v; d2vdz = -v
          d2udxy = v; d2vdxy = u
          d2udxz = -a*sin(x)*sin(y)*cos(z)
          d2udyz = a*cos(x)*cos(y)*cos(z)
          d2vdxz = -d2udyz; d2vdyz = -d2udxz          
       end if
    elseif(ncase==13) then      
       if(ndime==2) then        ! u=(x^3,-3*x^2*y)
          u=x**3
          v=-3.0_rp*(x**2)*y
          !
          dudx = 3.0_rp*(x**2)
          dvdx = -6.0_rp*x*y
          dvdy = -dudx
          !
          d2udx = 6.0_rp*x
          d2vdxy = -d2udx
          d2vdx = -6.0_rp*y
       else if(ndime==3) then ! u=(x^3,y^2,-(3*x^2+3*y^2)*z)
          u=x**3
          v=y**3
          w=-(3.0_rp*(x**2)+3.0_rp*(y**2))*z
          !
          dudx = 3.0_rp*(x**2)
          dvdy = 3.0_rp*(y**2)
          dwdz = -3.0_rp*((x**2)+(y**2))
          dwdx = -6.0_rp*x*z
          dwdy = -6.0_rp*y*z
          !
          d2udx = 6.0_rp*x; d2wdx = -6.0_rp*z
          d2vdy = 6.0_rp*y; d2wdy = -6.0_rp*z
          d2wdxz = -6.0_rp*x
          d2wdyz = -6.0_rp*y    
       end if
!!$    elseif(ncase==14) then      ! Stokes' operator singular solution (u)
!!$       ! WARNING: 2nd derivatives not computed!!!!
!!$       lambda = 0.54448373678246_rp
!!$       n = 1.0_rp
!!$       call polar_values(x,y,z,lambda,n,r,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,&
!!$            expd5,expd6,expp1,expp2,expp3)
!!$       !
!!$       u = (r**lambda)*expr1
!!$       v = (r**lambda)*expr2
!!$       dudx = (r**(lambda-1.0_rp))*expd1
!!$       dudy = (r**(lambda-1.0_rp))*expd2
!!$       dvdx = (r**(lambda-1.0_rp))*expd3
!!$       dvdy = (r**(lambda-1.0_rp))*expd4
!!$       !
!!$       d2udx = (r**(lambda-2.0_rp))*expd5
!!$       d2udy = 0.0_rp ! The laplacian is fully calculated in d2udx
!!$       d2vdx = (r**(lambda-2.0_rp))*expd6
!!$       d2vdy = 0.0_rp ! The laplacian is fully calculated in d2vdx
!!$    elseif(ncase==15) then      ! Maxwell's operator singular solution (b)
!!$       ! WARNING: 2nd derivatives not computed!!!!
!!$       lambda = 0.54448373678246_rp
!!$       n = 1.0_rp
!!$       call polar_values(x,y,z,lambda,n,r,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,&
!!$            expd5,expd6,expp1,expp2,expp3)
!!$       !
!!$       u = (2.0_rp*n/3.0_rp)*r**((2.0_rp*n-3.0_rp)/3.0_rp)*ex1
!!$       v = (2.0_rp*n/3.0_rp)*r**((2.0_rp*n-3.0_rp)/3.0_rp)*ex2
!!$       !
!!$       dudx = ((4.0_rp*n*n-6.0_rp*n)/9.0_rp)*r**((2.0_rp*n-6.0_rp)/3.0_rp)*&
!!$            (COS(theta)*ex1 - SIN(theta)*ex2)
!!$       dudy = ((4.0_rp*n*n-6.0_rp*n)/9.0_rp)*r**((2.0_rp*n-6.0_rp)/3.0_rp)*&
!!$            (SIN(theta)*ex1 + COS(theta)*ex2)
!!$       dvdx = ((4.0_rp*n*n-6.0_rp*n)/9.0_rp)*r**((2.0_rp*n-6.0_rp)/3.0_rp)*&
!!$            (COS(theta)*ex2 + SIN(theta)*ex1)
!!$       dvdy = ((4.0_rp*n*n-6.0_rp*n)/9.0_rp)*r**((2.0_rp*n-6.0_rp)/3.0_rp)*&
!!$            (SIN(theta)*ex2 - COS(theta)*ex1)
!!$       !
!!$       if(x==0.0_rp .and. y==0.0_rp) then
!!$          u = 0.0_rp; v = 0.0_rp
!!$          dudx = 0.0_rp; dudy = 0.0_rp
!!$          dvdx = 0.0_rp; dvdy = 0.0_rp
!!$       end if
    elseif(ncase==16) then     ! Chanel flow perturbation
       u = alpha*sin(x)*(y**2-1.0_rp)
       dudx = alpha*cos(x)*(y**2-1.0_rp)
       dudy = 2.0_rp*alpha*sin(x)*y
       d2udx = -alpha*sin(x)*(y**2-1.0_rp)
       d2udy = 2.0_rp*alpha*sin(x)
       d2udxy = 2.0_rp*alpha*cos(x)*y
    elseif(ncase==17) then    ! Cavity
       if(ndime==2) then
          if(abs(y).lt.1.0e-8_rp.and.(abs(x).gt.1.0e-8_rp.and.abs(x-1.0_rp).gt.1.0e-8_rp)) then
             u = 1.0_rp
          end if
       end if      
    elseif(ncase==18) then     ! u = (x^2, y^2, -z(2x+2y))
       u = x*x
       v = y*y
       w = -z*(2.0_rp*x + 2.0_rp*y)
       dudx = 2.0_rp*x
       dvdy = 2.0_rp*y
       dwdx = -2.0_rp*z
       dwdy = -2.0_rp*z
       dwdz = -2.0_rp*x - 2.0_rp*y
       !
       d2udx = 2.0_rp
       d2vdy = 2.0_rp
       d2wdxz = -2.0_rp
       d2wdyz = -2.0_rp
    elseif(ncase==19) then   ! constant advection field
        u = 1.0_rp;   v = 1.0_rp
          
          dudx = 0.0_rp
          dudy = 0.0_rp
          dvdx = 0.0_rp
          dvdy = 0.0_rp
          !
          d2udx = 0.0_rp
          d2vdxy = 0.0_rp
          d2vdx = 0.0_rp
    elseif(ncase==20) then  ! u=uo*((1-y^2)*(1-z^2),0,0) parabolic inlet
       uo=0.5_rp
       !L=1.0_rp/2.0_rp
       L=2.0_rp/2.0_rp
       yy=(y-L)/L
       zz=(z-L)/L

       u = uo*((1-yy*yy)*(1-zz*zz))
       v = 0.0_rp
       w = 0.0_rp
       !
       dudx = 0.0_rp
       dudy = -uo*2.0_rp**yy*(1-zz*zz)
       dudz = -uo*2.0_rp*zz*(1-yy*yy)
       !
       d2udyz = uo*4.0_rp*yy*zz
    elseif(ncase==21.or.ncase==22) then   ! Hunt's velocity or current (insulated side walls), finite dB, centered velocity
       Ha=10.0_rp     !Hartmann number
       db=0.0_rp      !Hartmann wall conductivity ratio
       !db=1.0e30_rp  !Hartmann wall conductivity ratio
       visco=1.0_rp   !nu*rho
       condu=1.0_rp   !electrical conductivity
       a=1.0_rp       !dimension along B
       b=1.0_rp
       dpdx=15.6467_rp !pressure gradient
       !dpdx=63.56_rp  !pressure gradient
       if(alpha>0.0_rp) dpdx=alpha !used in the analytical solution for a given Uavg
       !
       if(ncase==21) then !velocity and velocity gradient
          call hunt_V(x,y,z,Ha,db,condu,visco,a,b,dpdx,Vhunt,1,1)
          K = a*a*dpdx/visco
          u = K*Vhunt(1)
          v = 0.0_rp
          w = 0.0_rp 
          dudy = K*Vhunt(2)
          dudz = K*Vhunt(3)
          d2udy = K*Vhunt(4)
          d2udyz = K*Vhunt(5)
          d2udz = K*Vhunt(6) 
       else !current, divergence free
          call hunt_V(x,y,z,Ha,db,condu,visco,a,b,dpdx,Vhunt,2,1)
          K = a*a*sqrt(condu)*dpdx/sqrt(visco)
          u = 0.0_rp
          v = K*Vhunt(1)
          w = K*Vhunt(2)
          dvdy = K*Vhunt(3)
          dvdz = K*Vhunt(4)
          dwdy = K*Vhunt(5)
          dwdz = K*Vhunt(6)
          d2vdy = K*Vhunt(7)
          d2vdyz = K*Vhunt(8)
          d2vdz = K*Vhunt(9)
          d2wdy = K*Vhunt(10)
          d2wdyz = K*Vhunt(11)
          d2wdz = K*Vhunt(12)
       end if  
    elseif(ncase==23) then ! Cylinder inlet
       if(x.le.1.0e-8_rp) then
          if(ndime==2) then
             u = 4.0_rp*1.5_rp*y*(0.41-y)/(0.41**2)
          elseif(ndime==3) then
             u = 16.0_rp*2.25_rp*y*z*(0.41-y)*(0.41-z)/(0.41**4)
          end if
       end if
       a = sqrt((x-0.5_rp)**2+(y-0.2_rp)**2)
       if(abs(a-0.05_rp).le.1.0e-4_rp) then ! Cylinder perturbation
          !if((x-0.5_rp).le.0.0_rp) then
          !   theta = atan((y-0.2_rp)/abs(x-0.5_rp))
          !   u = 1.0_rp*sin(theta)
          !   v = 1.0_rp*cos(theta)
          !else
          !   theta = atan((y-0.2_rp)/(x-0.5_rp))
          !   u = 1.0_rp*sin(theta)
          !   v = -1.0_rp*cos(theta)
          !end if
       end if
    elseif(ncase==24) then ! Cylinder parabolic
       if(((x-0.4_rp).le.1.0e-8_rp).or.((x-0.6_rp).ge.1.0e-8_rp)) then
          u = 16.0_rp*2.25_rp*y*z*(0.41-y)*(0.41-z)/(0.41**4)
       end if
    elseif(ncase==25) then ! Double-glazing problem (advection - 2D)
       if(ndime==2) then
          u = 2.0_rp*y*(1.0_rp-x**2)
          v = -2.0_rp*x*(1.0_rp-y**2)
          dudx = -4.0_rp*x*y
          dvdx = -2.0_rp*(1.0_rp-y**2)
          dudy = 2.0_rp*(1-x**2)
          dvdy = 4.0_rp*x*y
          d2udx = -4.0_rp*y
          d2vdy = 4.0_rp*x
          d2udxy = -4.0_rp*x
          d2vdxy = 4.0_rp*y
       end if
    elseif(ncase==26) then ! Colliding flow (2D)
       if(ndime==2) then
          u = 20.0_rp*x*y**3
          v = 5.0_rp*x**4-5.0_rp*y**4
          dudx = 20.0_rp*y**3
          dvdx = 20.0_rp*x**3
          dudy = 60.0_rp*x*y**2
          dvdy = -20.0_rp*y**3
          d2udy = 120.0_rp*x*y
          d2vdx = 60.0_rp*x**2
          d2vdy = -60.0_rp*y**2
          d2udxy = 60.0_rp*y**2
       end if 
    elseif(ncase==27) then ! Backward-facing step inlet
       if(x.le.1.0e-8_rp) then
          if((y-1.0_rp).ge.1.0e-8_rp) then
             !u=1.5_rp*4.0_rp*y*(1.0_rp-y)
             u=1.5_rp*4.0_rp*(2.0_rp-y)*(y-1.0_rp)
          else
             u=0.0_rp
          end if
       end if
    elseif(ncase==28) then ! NACA0012 inlet
       a = 10.0_rp*pi/180.0_rp
       if(abs(y).le.0.6_rp) then
          if(x.le.100.0_rp) then
             u = 1.0_rp*cos(a)
             v = 1.0_rp*sin(a)
          end if
       else
          u = 1.0_rp*cos(a)
          v = 1.0_rp*sin(a)
       end if
    end if


    ! Component assignment
    fpara(1) = u; fpara(2) = v; fpara(3) = w
    !
    fpara(4) = dudx; fpara(5) = dudy; fpara(6) = dudz
    fpara(7) = dvdx; fpara(8) = dvdy; fpara(9) = dvdz
    fpara(10) = dwdx; fpara(11) = dwdy; fpara(12) = dwdz
    !
    fpara(13) = d2udx; fpara(14) = d2udy; fpara(15) = d2udz
    fpara(16) = d2vdx; fpara(17) = d2vdy; fpara(18) = d2vdz
    fpara(19) = d2wdx; fpara(20) = d2wdy; fpara(21) = d2wdz
    fpara(22) = d2udxy; fpara(23) = d2udxz; fpara(24) = d2udyz
    fpara(25) = d2vdxy; fpara(26) = d2vdxz; fpara(27) = d2vdyz
    fpara(28) = d2wdxy; fpara(29) = d2wdxz; fpara(30) = d2wdyz

    return

  end subroutine analytical_vector_field

  !=============================================================================
  subroutine analytical_scalar_field(ncase,ndime,coord,fpara,time)
    implicit none
    integer(ip)      , intent(in)    :: ncase,ndime
    real(rp),optional, intent(in)    :: time
    real(rp)         , intent(inout) :: fpara(10)
    real(rp)         , intent(in)    :: coord(ndime)
    ! Local variables
    real(rp)    :: pi
    real(rp)    :: x,y,z,u,dudx,dudy,dudz
    real(rp)    :: d2udx,d2udy,d2udz,d2udxy,d2udxz,d2udyz
    real(rp)    :: drdx,drdy,dsdx,dsdy,d2rdxx,d2rdyy,d2rdxy,d2sdxx,d2sdyy,d2sdxy
    real(rp)    :: lambda,n,r,theta,tanh,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,expd5,expd6
    real(rp)    :: x0,y0,r0,expp1,expp2,expp3
   real(rp)    :: slope, t, t2, t_u_1, t_u_2

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

    if(ncase==1) then         ! u=x+y
       u = 1.0_rp*(x+y)
       dudx = 1.0_rp
       dudy = 1.0_rp 
    elseif(ncase==2) then     ! u=x^2+y^2
       u = x*x+y*y
       dudx = 2.0_rp*x
       dudy = 2.0_rp*y
       !
       d2udx = 2.0_rp
       d2udy = 2.0_rp
    elseif(ncase==3) then     ! u=sin(2*pi*x)
       u = sin(2.0_rp*pi*x)
       dudx = 2.0_rp*pi*cos(2.0_rp*pi*x)
       !
       d2udx = -4.0_rp*pi*pi*sin(2.0_rp*pi*x)
    elseif(ncase==4) then     ! u=1
       u = 1.0_rp
    elseif(ncase==5) then     ! u=x+y+z
       u = x+y+z
       dudx = 1.0_rp
       dudy = 1.0_rp
       dudz = 1.0_rp
!!$    elseif(ncase==6) then     ! Stokes's operator singular solution (p)
!!$       !WARNING: 2nd derivatives not computed
!!$       lambda = 0.54448373678246_rp
!!$       n = 1.0_rp
!!$       call polar_values(x,y,z,lambda,n,r,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,&
!!$            expd5,expd6,expp1,expp2,expp3)
!!$       !
!!$       u    = - (r**(lambda-1.0_rp))*expp1
!!$       dudx = (r**(lambda-2.0_rp))*expp2
!!$       dudy = (r**(lambda-2.0_rp))*expp3
    elseif(ncase==7) then      ! Taylor-Green vortex initial pressure
       if(ndime==2) then
          lambda = 1.0_rp/4.0_rp
          theta  = 4.0_rp*pi
          ex1 = cos(theta*x)+cos(theta*y)
          !
          u      = lambda*ex1
          dudx   = -theta*lambda*sin(theta*x)
          dudy   = -theta*lambda*sin(theta*y)
          d2udx  = -theta*theta*lambda*cos(theta*x)
          d2udy  = -theta*theta*lambda*cos(theta*y)
          d2udxy = 0.0_rp
       elseif(ndime==3) then
          lambda = 1.0_rp/16.0_rp
          ex1 = cos(2.0_rp*x)+cos(2.0_rp*y)
          ex2 = cos(2.0_rp*z)+2.0_rp
          !
          u     = lambda*ex1*ex2
          dudx  = -2.0_rp*lambda*sin(2.0_rp*x)*ex2
          dudy  = -2.0_rp*lambda*sin(2.0_rp*y)*ex2
          dudz  = -2.0_rp*lambda*sin(2.0_rp*z)*ex1
          !
          d2udx  = -4.0_rp*lambda*cos(2.0_rp*x)*ex2
          d2udy  = -4.0_rp*lambda*cos(2.0_rp*y)*ex2
          d2udz  = -4.0_rp*lambda*cos(2.0_rp*z)*ex1 
          d2udxz = 4.0_rp*lambda*sin(2.0_rp*x)*sin(2.0_rp*z)
          d2udxy = 4.0_rp*lambda*sin(2.0_rp*y)*sin(2.0_rp*z)
       end if
    elseif(ncase==8) then     ! u=x^3+y^3+z^3
       u = x*x*x+y*y*y+z*z*z
       dudx = 3.0_rp*x*x
       dudy = 3.0_rp*y*y
       dudz = 3.0_rp*z*z
       !
       d2udx = 6.0_rp*x
       d2udy = 6.0_rp*y
       d2udz = 6.0_rp*z
    elseif(ncase == 9) then   ! Unidireccional steady sinusoidal
       u = sin(2*pi*x)
       dudx = 2*pi*cos(2*pi*x)
       !
       d2udx =- 4*pi*pi*sin(2*pi*x)
    elseif(ncase == 10) then   ! Unidireccional shock
       x = modulo(x-time,1.0)
       if (x>=0.25 .and. x<=0.75) u = 1.0
!!$    elseif(ncase == 11) then   ! Casoni 2D figures 
!!$       !WARNING: derivatives not computed
!!$       x=2.0*x-1.0
!!$       y=2.0*y-1.0
!!$       r = (x+0.5)**2 + y**2
!!$       if (r.le. 0.25) then
!!$          u = (cos(2*pi*r))**2
!!$       elseif (x.gt.0.1 .and. x.lt.0.6 .and. y.gt. -0.25 .and. y.lt.0.25) then
!!$          u = 1.0
!!$       end if
    elseif(ncase == 12) then  ! Unidireccional sinusoidal with time
       x = modulo(x-time,1.0)
       u = sin(2*pi*x)
       !dudx = 2*pi*cos(2*pi*x)
       !d2udx =- 4*pi*pi*sin(2*pi*x)
!!$    elseif(ncase == 13) then
!!$       !WARNING: derivatives not computed
!!$       r = (x-0.5)**2 + (y-0.65)**2
!!$       if (r.le. 0.25**2) then
!!$          u = 1
!!$          if (x.gt.0.45 .and. x.lt.0.55 .and. y.gt. 0.3 .and. y.lt.0.6) u = 0.0
!!$       end if
    elseif(ncase==14) then     ! u=x^3+y^3
       u = x*x*x+y*y*y
       dudx = 3.0_rp*x*x
       dudy = 3.0_rp*y*y
       !
       d2udx = 6.0_rp*x
       d2udy = 6.0_rp*y
    elseif(ncase == 15) then  ! 3 body rotation
       !WARNING: derivatives not computed
       r0 = 0.15_rp
       !First figure
       x0 = 0.5_rp
       y0 = 0.75_rp
       r = sqrt((x-x0)**2 + (y-y0)**2)/r0
       if ( r .le. 1.0_rp) then
          if ((abs(x-x0) .ge. 0.025_rp) .or. (y .ge. 0.85_rp) ) then
             u = 1_rp
          else
             u = 0_rp
          end if
       end if
       !Second figure
       x0 = 0.5_rp
       y0 = 0.25_rp
       r = sqrt((x-x0)**2 + (y-y0)**2)/r0
       if ( r  .le. 1.0_rp) then
         u = 1_rp-r
       end if
       !Third figure
       x0 = 0.25_rp
       y0 = 0.5_rp
       r = sqrt((x-x0)**2 + (y-y0)**2)/r0
       if ( r  .le. 1.0_rp) then
          u = 0.25_rp*(1 + cos(pi * r))
       end if
    elseif(ncase==16) then     ! u=x^2
       u = x*x
       dudx = 2.0_rp*x
       !
       d2udx = 2.0_rp
    elseif(ncase==17) then     ! u=y^2
       u = y*y
       dudy = 2.0_rp*y
       !
       d2udy = 2.0_rp
    elseif(ncase==18) then     ! u=z^2
       u = z*z
       dudz = 2.0_rp*z
       !
       d2udz = 2.0_rp
    elseif(ncase==19) then     ! u=x^2+y^2+z^2
       u = x*x+y*y+z*z
       dudx = 2.0_rp*x
       dudy = 2.0_rp*y
       dudz = 2.0_rp*z
       !
       d2udx = 2.0_rp
       d2udy = 2.0_rp
       d2udz = 2.0_rp
!!$    elseif(ncase == 20) then   ! Smooth convergence test
!!$       ! WARNING: derivatives not computed!!!!
!!$       r = (x+0.5)**2 + y**2
!!$       if (r.le. 0.25) then
!!$          u = (cos(2*pi*r))**2
!!$       end if
    elseif(ncase == 21) then   ! Unidireccional shock (No periodic boundaries)
       x = x-time
       if (x>=0.25 .and. x<=0.75) u = 1.0
!!$    elseif(ncase==22) then     ! u=tanh((x^2+y^2)/9 - 1)
!!$       ! WARNING: derivatives not computed!!!!
!!$      theta = ((x-0.4_rp)**2 + y**2)/0.3_rp**2 - 1.0_rp
!!$      tanh = (exp(theta) - exp(-theta)) / (exp(theta) + exp(-theta)) 
!!$      u = 0.5_rp*(1.0_rp-tanh)
!!$      dudx = 0.0_rp
!!$      d2udx = 0.0_rp
   elseif(ncase==23) then ! u = x except on 3 boundary segments  
      if ((x<1.0_rp) .and. (y>0.0_rp) .and. (y<1.0_rp)) u = x
      dudx = 1.0_rp
   elseif(ncase==24) then   
      r = x*sin(-pi/3.0_rp) - (y-0.7_rp)*cos(-pi/3.0_rp)
      if (abs(y) < 1.0e-8_rp) then
         u =0.0_rp
      elseif (abs(x-1.0_rp) <  1.0e-8_rp) then
         u = 0.0_rp
      elseif (abs(r) < 10e-9_rp) then
         u = 0.5_rp
      elseif (r > 0) then
         u = 0.0_rp
      else
         u = 1.0_rp
      end if
   elseif(ncase==25) then  
      if ((x>=0.3_rp) .or. (y==1.0_rp)) then
         u = 0.0_rp
      else
         u = 1.0_rp
      end if
   elseif(ncase==26) then   
      if ((x>0.25_rp) .and. (x<0.75_rp)) then
         if ((y>0.25_rp) .and. (y<0.75_rp)) then
            u = (4.0_rp*x-1.0_rp)*(3.0_rp-4.0_rp*x)
            dudx = 4.0_rp*(3.0_rp-4.0_rp*x) - (4.0_rp*x-1.0_rp)*4.0_rp
            !
            d2udx = -32.0_rp
         end if
      end if
   elseif(ncase==27) then            
      if ((x<=0.9_rp) .and. (y>=0.5+0.4_rp*sin(pi/3.0_rp))) then
         u = 1.0_rp
      else
         u = 0.0_rp
      end if
   elseif(ncase==28) then            
      if ((x<=0.5_rp) .and. (y<=0.5)) then
         u = 1.0_rp
      else
         u = 0.0_rp
      end if
   elseif(ncase==29) then            
      if ((x<=-0.5_rp) .and. (y>=-0.5_rp)) then
         u = 1.0_rp
      else
         u = 0.0_rp
      end if
!!$   elseif(ncase==30) then 
!!$       ! WARNING: derivatives not computed!!!!
!!$      r = sqrt(x**2.00_rp + y**2.00_rp)
!!$      lambda = 5.00_rp*pi/3.00_rp
!!$      if ((r<=0.65_rp) .and. (r>=0.35_rp)) then
!!$         u = (cos(lambda*(2.00_rp*r-1.00_rp)))
!!$      else
!!$         u = 0.0_rp
!!$      end if
   elseif(ncase==31) then 
      r = sqrt(x**2.00_rp + y**2.00_rp)
      if ((r<=0.65_rp) .and. (r>=0.35_rp)) then
         u = 1.00_rp
      else
         u = 0.0_rp
      end if
   elseif(ncase==32) then 
      r = 100.0_rp
      x = r * x
      u = (exp(x)-exp(-x))/(exp(r)-exp(-r))
      dudx = r * (exp(x)+exp(-x))/(exp(r)-exp(-r))
      !
      d2udx =  r**2.0_rp * (exp(x)-exp(-x))/(exp(r)-exp(-r)) 
   elseif(ncase==33) then      ! sin(2*\pi(x+y)) + 0.005*cos(2*\pi(64x+63y))
      u = sin(2.0_rp*pi*(x+y)) + 0.005_rp*cos(2.0_rp*pi*(64*x+63*y))
      !
      dudx = 2.0_rp*pi*cos(2.0_rp*pi*(x+y)) - 0.005_rp*2.0_rp*pi*64*sin(2.0_rp*pi*(64*x+63*y))
      dudy = 2.0_rp*pi*cos(2.0_rp*pi*(x+y)) - 0.005_rp*2.0_rp*pi*63*sin(2.0_rp*pi*(64*x+63*y))
      !
      d2udx = - 4.0_rp*pi*pi*sin(2.0_rp*pi*(x+y)) - 0.005*4.0_rp*pi*pi*64*64*cos(2.0_rp*pi*(64*x+63*y))
      d2udy = - 4.0_rp*pi*pi*sin(2.0_rp*pi*(x+y)) - 0.005*4.0_rp*pi*pi*63*63*cos(2.0_rp*pi*(64*x+63*y))
      d2udxy = - 4.0_rp*pi*pi*sin(2.0_rp*pi*(x+y)) - 0.005*4.0_rp*pi*pi*64*63*cos(2.0_rp*pi*(64*x+63*y))
    elseif(ncase==34) then     ! Stokes's operator singular solution (p)
      theta = ((x-0.4_rp)**2 + y**2)/0.3_rp**2 - 1.0_rp
      tanh = (exp(theta) - exp(-theta)) / (exp(theta) + exp(-theta)) 
      u = 0.5_rp*(1.0_rp-tanh)
      dudx = 0.0_rp
      d2udx = 0.0_rp
   elseif(ncase==35) then     ! sin(2*pi*x)*sin(2*pi*y)
      theta = 2.0_rp*pi
      u = sin(theta*x)*sin(theta*y)
      dudx = theta*cos(theta*x)*sin(theta*y)
      dudy = theta*sin(theta*x)*cos(theta*y)
      d2udx = - theta**2*u
      d2udy = d2udx
      d2udxy = theta**2*cos(theta*x)*cos(theta*y)
   elseif(ncase==36) then     ! sin(2*pi*x)*sin(2*pi*y)
      x = x-0.5_rp
      y = y-0.5_rp
      r = sqrt(x**2.0_rp+y**2.0_rp)
      lambda = 2.0_rp !\sigma = 1/2lambda = 0.025
      u = lambda*exp(-lambda*r)
      dudx = -lambda*x/r*u
      dudy = -lambda*y/r*u
      d2udx = lambda*u/r*(-1+(x/r)**2.0_rp+lambda*x**2.0_rp/r)
      d2udy = lambda*u/r*(-1+(y/r)**2.0_rp+lambda*y**2.0_rp/r)
      d2udxy = lambda*x*y*u/r**2.0_rp*(-1/r+lambda)

    ! Internal layer problem, 2 variants 
   elseif((ncase == 37) .or. (ncase == 38)) then
      if(ncase == 37) then 
         x = x - 1.25_rp
         y = y + 0.25_rp
      end if
      if(ncase == 38) then 
         x = x + 0.25_rp
         y = y + 0.25_rp
      end if
      slope = 60
      u =  atan(slope * ((x**2 + y**2)**0.5 - pi/3));

      t2 = x**2 + y**2
      t = (t2)**0.5;
      t_u_1 = t * (slope**2 * (t - pi/3)**2 + 1);
      dudx = slope * x / t_u_1;
      dudy = slope * y / t_u_1;

      ! these values are probably not correct (at least d2udxy), but the laplacian d2udx + d2udy should be correct
      ! todo: derive properly!

      t_u_2 = ((pi - 3.0_rp*t)**2 * slope**2 + 9.0);

      d2udx = 27.0 * 2.0 * y**2 * (pi - 3.0*t) * slope**3.0 / (t_u_2**2 * t2) &
           - 9.0 * y**2 * slope / (t_u_2 * (t ** 3.0)) + 9.0 * slope / (t_u_2 * t)

      d2udy = 27.0 * 2.0 * x**2 * (pi - 3.0*t) * slope**3.0 / (t_u_2**2 * t2) &
           - 9.0 * x**2 * slope / (t_u_2 * (t ** 3.0)) + 9.0 * slope / (t_u_2 * t)

      d2udxy = 0

   ! Fichera corner domain
   elseif(ncase == 39) then
      lambda = 2.0_rp/3.0_rp
      theta = atan2(y,x)+pi/2.0_rp
      r = x*x + y*y
      u = r**(1.0_rp/3.0_rp)*sin(lambda*theta)
      !u = (x*x + y*y)**(1.0/3.0) * sin((2.0*atan2(y, x) + pi) /3.0)

      if (r>1.0e-14_rp) then
         drdx = lambda*x*r**(-lambda)
         drdy = lambda*y*r**(-lambda)
         dsdx = -lambda*y/r*cos(lambda*theta)
         dsdy =  lambda*x/r*cos(lambda*theta)

         dudx = drdx*sin(lambda*theta) + r**(1.0_rp/3.0_rp)*dsdx
         dudy = drdy*sin(lambda*theta) + r**(1.0_rp/3.0_rp)*dsdy

         d2rdxx = lambda*r**(-lambda) - 2.0_rp*lambda**2.0_rp*x**2.0_rp*r**(-5.0_rp/3.0_rp)
         d2rdxy = - 2.0_rp*lambda**2.0_rp*x*y*r**(-5.0_rp/3.0_rp)
         d2rdyy = lambda*r**(-lambda) - 2.0_rp*lambda**2.0_rp*y**2.0_rp*r**(-5.0_rp/3.0_rp)

         d2sdxx =  2.0_rp*lambda*x*y*r**(-2.0_rp)*cos(lambda*theta) - (lambda*y/r)**2.0_rp*sin(lambda*theta)
         d2sdyy = -2.0_rp*lambda*x*y*r**(-2.0_rp)*cos(lambda*theta) - (lambda*x/r)**2.0_rp*sin(lambda*theta)
         d2sdxy = lambda*(y**2.0_rp-x**2.0_rp)*r**(-2.0_rp)*cos(lambda*theta) -                     &
              &   (lambda/r)**2.0_rp*x*y*sin(lambda*theta)

         d2udx  = d2rdxx*sin(lambda*theta) + 2.0_rp*drdx*dsdx + r**(1.0_rp/3.0_rp)*d2sdxx
         d2udxy = d2rdxy*sin(lambda*theta) + drdy*dsdx + drdx*dsdy + r**(1.0_rp/3.0_rp)*d2sdxy
         d2udy  = d2rdyy*sin(lambda*theta) + 2.0_rp*drdy*dsdy + r**(1.0_rp/3.0_rp)*d2sdyy
      end if
      
      ! derivatives not calculated
      ! however, the laplacian of the function u should be 0
      ! todo: derive properly
      dudx = 0.0
      dudy = 0.0
      d2udx = 0.0
      d2udy = 0.0
      d2udxy = 0.0
   elseif (ncase == 111) then 
      u=1.0_rp
      dudx = 0.0_rp
      dudy = 0.0_rp

      
   elseif(ncase==40) then     ! u=x^7+y^7+z^7
      r = 7.0_rp
      u = x**r+y**r+z**r
      dudx = r*x**(r-1.0_rp)
      dudy = r*y**(r-1.0_rp)
      dudz = r*z**(r-1.0_rp)
 
       d2udx = r*(r-1.0_rp)*x**(r-2.0_rp)
       d2udy = r*(r-1.0_rp)*y**(r-2.0_rp)
       d2udz = r*(r-1.0_rp)*z**(r-2.0_rp)
    elseif(ncase==41) then ! Colliding flow (2D pressure field)
       if(ndime==2) then
          u = 60.0_rp*x**2*y-20.0_rp*y**3+40
          dudx = 120.0_rp*x*y
          dudy = 60.0_rp*x**2-60.0_rp*y**2
          d2udx = 120.0_rp*y
          d2udy = -120.0_rp*y
          d2udxy = 120.0_rp*x
       end if
   end if

    ! Component assignment
    fpara(1) = u
    fpara(2) = dudx; fpara(3) = dudy; fpara(4) = dudz
    !
    fpara(5) = d2udx; fpara(6) = d2udy; fpara(7) = d2udz
    fpara(8) = d2udxy; fpara(9) = d2udxz; fpara(10) = d2udyz


    return

  end subroutine analytical_scalar_field

  !=============================================================================
  subroutine analytical_temporal_field(ncase,t,fpara,alpha)
    implicit none
    integer(ip), intent(in)  :: ncase
    real(rp)   , intent(inout)  :: fpara(3)
    real(rp)   , intent(in)  :: t,alpha
    ! Local variables
    real(rp)    :: u,dudt,d2udt,a,b

    u=1.0_rp; dudt=0.0_rp; d2udt=0.0_rp

    if(ncase==1) then
       u = t
       dudt = 1.0_rp
    elseif(ncase==2) then
       u = t*t
       dudt = 2.0_rp*t
       d2udt = 2.0_rp
    elseif(ncase==3) then
       a = pi/10.0_rp
       b = 1.0_rp/25.0_rp
       u = sin(a*t)*exp(b*t)
       dudt = a*cos(a*t)*exp(b*t) + b*sin(a*t)*exp(b*t)
       d2udt = (b*b-a*a)*sin(a*t)*exp(b*t) + 2.0_rp*a*b*cos(a*t)*exp(b*t)
    elseif(ncase==4) then
       a = 2.0_rp*pi*pi*alpha
       u = exp(-a*t)
       dudt = -a*exp(-a*t)
       d2udt = a*a*exp(-a*t)
    elseif(ncase==5) then
       a = 4.0_rp*pi*pi*alpha
       u = exp(-a*t)
       dudt = -a*exp(-a*t)
       d2udt = a*a*exp(-a*t)
    elseif(ncase==6) then
       u = t*t*t
       dudt = 3.0_rp*t*t
       d2udt = 6.0_rp*t
    elseif(ncase==7) then
       u = t*t*t*t
       dudt = 4.0_rp*t*t*t
       d2udt = 12.0_rp*t*t
    end if

    ! Component assignment
    fpara(1) = u
    fpara(2) = dudt
    fpara(3) = d2udt

  end subroutine analytical_temporal_field

  !=============================================================================
  subroutine polar_values(x,y,z,lambda,n,r,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,&
       expd5,expd6,expp1,expp2,expp3)
    implicit none
    real(rp)   , intent(inout):: x,y,z,lambda,n
    real(rp)   , intent(out)  :: r,theta,ex1,ex2,expr1,expr2,expd1,expd2,expd3,expd4,expd5,expd6,expp1,expp2,expp3
    ! Local variables
    real(rp)    :: sinpt,sinmt,cospt,cosmt,coslo,fi,d1fi,d2fi,d3fi,d4fi,omega

    omega = 1.50_rp*pi
    !
    if (ABS(x).lt.1.0e-8) then
       x=0.0_rp
    end if
    if (ABS(y).lt.1.0e-8) then
       y=0.0_rp
    end if
    if (ABS(z).lt.1.0e-8) then
       z=0.0_rp
    end if
    !
    if (y==0.0_rp.and.ABS(x).gt.1.0e-8) then
       r=ABS(x)
       if (x.gt.0.0_rp) then
          theta=0.0_rp !!! 2*pi
       else
          theta=pi
       end if
    elseif (x==0.0_rp.and.ABS(y).gt.1.0e-8) then
       r=ABS(y)
       if (y.gt.0.0_rp) then
          theta=pi/2
       else
          theta=3*pi/2
       end if
    elseif (x==0.0_rp.and.y==0.0_rp) then
       r=0.0_rp
       theta=0.0_rp
    else
       r = SQRT(x*x + y*y)
       if (y.gt.0.0_rp) then
          theta = ATAN2(y,x)
       else
          theta = pi + ATAN2(-1.0_rp*y,-1.0_rp*x)
       end if
    end if
    !
    ex1 = SIN(2.0_rp*n*theta/3.0_rp)*COS(theta) - COS(2.0_rp*n*theta/3.0_rp)*SIN(theta)
    ex2 = SIN(2.0_rp*n*theta/3.0_rp)*SIN(theta) + COS(2.0_rp*n*theta/3.0_rp)*COS(theta) 
    !
    sinpt = SIN(theta*(1.0_rp+lambda))
    sinmt = SIN(theta*(1.0_rp-lambda))
    cospt = COS(theta*(1.0_rp+lambda))
    cosmt = COS(theta*(1.0_rp-lambda))
    coslo = COS(lambda*omega)
    !
    fi    = coslo*sinpt/(1.0_rp+lambda) - cospt - coslo*sinmt/(1.0_rp-lambda) + cosmt
    d1fi  = coslo*cospt + (1.0_rp+lambda)*sinpt - coslo*cosmt - (1.0_rp-lambda)*sinmt
    d2fi  = -(1.0_rp+lambda)*coslo*sinpt + ((1.0_rp+lambda)**2)*cospt + &
         (1.0_rp-lambda)*coslo*sinmt - ((1.0_rp-lambda)**2)*cosmt
    d3fi  = -((1.0_rp+lambda)**2)*coslo*cospt - ((1.0_rp+lambda)**3)*sinpt + &
         ((1.0_rp-lambda)**2)*coslo*cosmt + ((1.0_rp-lambda)**3)*sinmt
    d4fi  = ((1.0_rp+lambda)**3)*coslo*sinpt - ((1.0_rp+lambda)**4)*cospt - &
         ((1.0_rp-lambda)**3)*coslo*sinmt + ((1.0_rp-lambda)**4)*cosmt
    !
    expr1 = (1.0_rp+lambda)*SIN(theta)*fi + COS(theta)*d1fi
    expr2 = SIN(theta)*d1fi - (1.0_rp+lambda)*COS(theta)*fi
    !
    expd1 = lambda*COS(theta)*expr1 - SIN(theta)*((1.0_rp+lambda)*COS(theta)*fi + &
         lambda*SIN(theta)*d1fi + COS(theta)*d2fi)
    expd2 = lambda*SIN(theta)*expr1 + COS(theta)*((1.0_rp+lambda)*COS(theta)*fi + &
         lambda*SIN(theta)*d1fi + COS(theta)*d2fi)
    expd3 = lambda*COS(theta)*expr2 - SIN(theta)*((1.0_rp+lambda)*SIN(theta)*fi - &
         lambda*COS(theta)*d1fi + SIN(theta)*d2fi)
    expd4 = lambda*SIN(theta)*expr2 + COS(theta)*((1.0_rp+lambda)*SIN(theta)*fi - &
         lambda*COS(theta)*d1fi + SIN(theta)*d2fi)
    !
    expd5 = (lambda**2.0_rp)*expr1 - (1.0_rp+lambda)*SIN(theta)*fi + &
         (1.0_rp+2.0_rp*lambda)*COS(theta)*d1fi + (lambda-1.0_rp)*SIN(theta)*d2fi + &
         COS(theta)*d3fi
    expd6 = (lambda**2.0_rp)*expr2 + (1.0_rp+lambda)*COS(theta)*fi + &
         (1.0_rp+2.0_rp*lambda)*SIN(theta)*d1fi + (1.0_rp-lambda)*COS(theta)*d2fi + &
         SIN(theta)*d3fi
    !
    expp1 = (((1.0_rp+lambda)**2.0_rp)*d1fi + d3fi)/(1.0_rp-lambda)
    expp2 = (((1.0_rp+lambda)**2.0_rp)*d1fi*COS(theta) + &
         (((1.0_rp+lambda)**2.0_rp)*d2fi*SIN(theta)/(1.0_rp-lambda)) + d3fi*COS(theta) + &
         (d4fi*SIN(theta)/(1.0_rp-lambda)))
    expp3 = (((1.0_rp+lambda)**2.0_rp)*d1fi*SIN(theta) - &
         ((1.0_rp+lambda)**2.0_rp)*d2fi*COS(theta)/(1.0_rp-lambda) + d3fi*SIN(theta) - &
         d4fi*COS(theta)/(1.0_rp-lambda))


  end subroutine polar_values

  !=============================================================================
  subroutine hunt_V(x,y,z,Ha,db,condu,visco,a,b,dpdx,V,ncase,mesh)
    implicit none
    real(rp)   , intent(inout):: x,y,z,Ha,db,condu,visco,a,b
    real(rp)   , intent(out)  :: V(12)
    real(rp), optional  , intent(in) :: dpdx
    integer(ip), intent(in)  :: ncase,mesh
    ! Local variables
    integer(ip) :: k,kmax, Bdir, geo
    real(rp)    :: alphak,l,r1k,r2k,N,xi,eta,kr,C
    real(rp)    :: V2,V3,dV2deta,dV3deta,d2V2deta2,d2V3deta2
    real(rp)    :: H2,H3,dH2deta,dH3deta,d2H2deta2,d2H3deta2,d3H2deta3,d3H3deta3
    real(rp)    :: aux,num2,num3,den2,den3,coef2,coef3,coef2neg,coef3neg
    real(rp)    :: La,Lb
    integer(ip), parameter ::Bz = 1, By = 2

    kmax=1000
    l=b/a 
    Bdir=Bz !direction of B
    geo=mesh   !dimensions, 0: centered at zero (struct, nonuniform), 1: full duct section (GiD or struct uniform)
    if(geo==0) then
       La=0.0_rp
       Lb=0.0_rp
    else
       La=a
       Lb=b
    end if
    
    ! translate geometry in order to get the correct coordinate
    if(Bdir==Bz) then !Bz
       xi = (y-Lb)/a   !perpendicular to B    
       eta = (z-La)/a  !aligned with B
    else !By
       xi = (z-Lb)/a  !direction perpendicular to the magnetic field direction  
       eta = (y-La)/a !direction aligned with the magnetic field
    end if

    ! initial values By
    V(1) = 0.0_rp  !V         jy
    V(2) = 0.0_rp  !dVdy      jz
    V(3) = 0.0_rp  !dVdz      djydy 
    V(4) = 0.0_rp  !d2Vdydy   djydz
    V(5) = 0.0_rp  !d2Vdydz   djzdy
    V(6) = 0.0_rp  !d2Vdzdz   djzdz
    V(7) = 0.0_rp  !          d2jydyy
    V(8) = 0.0_rp  !          d2jydyz
    V(9) = 0.0_rp  !          d2jydzz
    V(10) = 0.0_rp !          d2jzdyy
    V(11) = 0.0_rp !          d2jzdyz
    V(12) = 0.0_rp !          d2jzdzz
  
    do k=0,kmax
       kr = REAL(k)
       alphak = (kr+0.5_rp)*pi/l
       N = sqrt(Ha*Ha + 4.0_rp*alphak*alphak)
       r1k = 0.5_rp*(Ha + N)
       r2k = 0.5_rp*(-Ha + N)
       
       aux = (1.0_rp-exp(-2.0_rp*(r1k+r2k)))

       den2 = db*N*0.5_rp*(1.0_rp+exp(-2.0_rp*r1k))+aux/(1.0_rp+exp(-2.0_rp*r2k))
       den3 = db*N*0.5_rp*(1.0_rp+exp(-2.0_rp*r2k))+aux/(1.0_rp+exp(-2.0_rp*r1k))

       num2 = db*r2k+(1.0_rp-exp(-2.0_rp*r2k))/(1.0_rp+exp(-2.0_rp*r2k))
       num3 = db*r1k+(1.0_rp-exp(-2.0_rp*r1k))/(1.0_rp+exp(-2.0_rp*r1k))

       coef2 = 0.5_rp*(exp(-r1k*(1.0_rp-eta)) + exp(-r1k*(1.0_rp+eta)))
       coef3 = 0.5_rp*(exp(-r2k*(1.0_rp-eta)) + exp(-r2k*(1.0_rp+eta)))
       coef2neg = 0.5_rp*(exp(-r1k*(1.0_rp-eta)) - exp(-r1k*(1.0_rp+eta)))
       coef3neg = 0.5_rp*(exp(-r2k*(1.0_rp-eta)) - exp(-r2k*(1.0_rp+eta)))

       V2 = num2*coef2/den2
       V3 = num3*coef3/den3

       H2 = num2*coef2neg/den2
       H3 = num3*coef3neg/den3

       C = 2.0_rp*((-1.0_rp)**kr)/l
       if(ncase==1) then !velocity
          V(1) = V(1) + C*(cos(alphak*xi)/(alphak**3.0_rp))*(1.0_rp-V2-V3)
          dV2deta = r1k*H2
          dV3deta = r2k*H3
          d2V2deta2 = (r1k**2.0_rp)*V2
          d2V3deta2 = (r2k**2.0_rp)*V3
          if(Bdir==Bz) then ! Bz
             V(2) = V(2) - (C/a)*(sin(alphak*xi)/(alphak**2.0_rp))*(1.0_rp-V2-V3)             !dVdy 
             V(3) = V(3) + (C/a)*(cos(alphak*xi)/(alphak**3.0_rp))*(-dV2deta-dV3deta)         !dVdz
             V(4) = V(4) - (C/(a*a))*(cos(alphak*xi)/alphak)*(1.0_rp-V2-V3)                   !d2Vdyy
             V(5) = V(5) - (C/(a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(-dV2deta-dV3deta)     !d2Vdyz  
             V(6) = V(6) + (C/(a*a))*(cos(alphak*xi)/(alphak**3.0_rp))*(-d2V2deta2-d2V3deta2) !d2Vdzz
          else ! By
             V(2) = V(2) + (C/a)*(cos(alphak*xi)/(alphak**3.0_rp))*(-dV2deta-dV3deta)         !dVdy
             V(3) = V(3) - (C/a)*(sin(alphak*xi)/(alphak**2.0_rp))*(1.0_rp-V2-V3)             !dVdz 
             V(4) = V(4) + (C/(a*a))*(cos(alphak*xi)/(alphak**3.0_rp))*(-d2V2deta2-d2V3deta2) !d2Vdyy
             V(5) = V(5) - (C/(a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(-dV2deta-dV3deta)     !d2Vdyz
             V(6) = V(6) - (C/(a*a))*(cos(alphak*xi)/alphak)*(1.0_rp-V2-V3)                   !d2Vdzz
          end if
       else !current
         ! V(1) = V(1) + C*(cos(alphak*xi)/(alphak**3.0_rp))*(H2-H3)
          dH2deta = r1k*V2
          dH3deta = r2k*V3
          d2H2deta2 = (r1k**2.0_rp)*H2
          d2H3deta2 = (r2k**2.0_rp)*H3
          d3H2deta3 = (r1k**3.0_rp)*V2
          d3H3deta3 = (r2k**3.0_rp)*V3
          if(Bdir==Bz) then !Bz
             V(1) = V(1) + (C/a)*(cos(alphak*xi)/(alphak**3.0_rp))*(dH2deta-dH3deta)           !jy
             V(2) = V(2) + (C/a)*(sin(alphak*xi)/(alphak**2.0_rp))*(H2-H3)                     !jz

             V(3) = V(3) - (C/(a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(dH2deta-dH3deta)         !djydy
             V(4) = V(4) + (C/(a*a))*(cos(alphak*xi)/(alphak**3.0_rp))*(d2H2deta2-d2H3deta2)     !djydz
             V(5) = V(5) + (C/(a*a))*(cos(alphak*xi)/alphak)*(H2-H3)                             !djzdy
             V(6) = V(6) + (C/(a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(dH2deta-dH3deta)         !djzdz

             V(7) = V(7) - (C/(a*a*a))*(cos(alphak*xi)/alphak)*(dH2deta-dH3deta)                 !d2jydyy
             V(8) = V(8) - (C/(a*a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(d2H2deta2-d2H3deta2)   !d2jydyz
             V(9) = V(9) + (C/(a*a*a))*(cos(alphak*xi)/(alphak**3.0_rp))*(d3H2deta3-d3H3deta3)   !d2jydzz
             V(10) = V(10) - (C/(a*a*a))*sin(alphak*xi)*(H2-H3)                                  !d2jzdyy
             V(11) = V(11) + (C/(a*a*a))*(cos(alphak*xi)/alphak)*(dH2deta-dH3deta)               !d2jzdyz
             V(12) = V(12) + (C/(a*a*a))*(sin(alphak*xi)/(alphak**2.0_rp))*(d2H2deta2-d2H3deta2) !d2jzdzz
          else !By
             V(1) = V(1) - (C/a)*(sin(alphak*xi)/(alphak**2.0_rp))*(H2-H3)                     !jy
             V(2) = V(2) - (C/a)*(cos(alphak*xi)/(alphak**3.0_rp))*(dH2deta-dH3deta)           !jz
             V(3) = V(3) - (C/a*a)*(sin(alphak*xi)/(alphak**2.0_rp))*(dH2deta-dH3deta)         !djydy
             V(4) = V(4) - (C/a*a)*(cos(alphak*xi)/alphak)*(H2-H3)                             !djydz
             V(5) = V(5) - (C/a*a)*(cos(alphak*xi)/(alphak**3.0_rp))*(d2H2deta2-d2H3deta2)     !djzdy
             V(6) = V(6) + (C/a*a)*(sin(alphak*xi)/(alphak**2.0_rp))*(dH2deta-dH3deta)         !djzdz
             V(7) = V(7) - (C/a*a*a)*(sin(alphak*xi)/(alphak**2.0_rp))*(d2H2deta2-d2H3deta2)   !d2jydyy
             V(8) = V(8) - (C/a*a*a)*(cos(alphak*xi)/alphak)*(dH2deta-dH3deta)                 !d2jydyz
             V(9) = V(9) + (C/a*a*a)*sin(alphak*xi)*(H2-H3)                                    !d2jydzz
             V(10) = V(10) - (C/a*a*a)*(cos(alphak*xi)/(alphak**3.0_rp))*(d3H2deta3-d3H3deta3) !d2jzdyy
             V(11) = V(11) + (C/a*a*a)*(sin(alphak*xi)/(alphak**2.0_rp))*(d2H2deta2-d2H3deta2) !d2jzdyz
             V(12) = V(12) + (C/a*a*a)*(cos(alphak*xi)/alphak)*(dH2deta-dH3deta)               !d2jzdzz
          end if
       end if
    end do

  end subroutine hunt_V
end module analytical_names
