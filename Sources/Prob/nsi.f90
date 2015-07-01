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
module nsi_names
  use types_names
  use memor_names
  use problem_names
  use finite_element_names
  use element_fields_names
  use analytical_names
  implicit none
# include "debug.i90"
  private 

  type, extends(physical_problem_t) :: nsi_problem_t
     integer(ip) ::   & 
          ksnsi,      & ! Symmetry flag (+-1) (1: no simetrica definida, -1: simetrica no definida)
          kfl_conv,   & ! Flag for enabling advection (Stokes=0; NavSto=1)
          kfl_symg,   & ! Flag for symmetric grad-grad term
          kfl_tder,   & ! Flag for time derivative computation
          kfl_skew,   & ! Flag for enabling skewsymmetric conv_tective terms (Off=0; type1=1, type2=2)
          kfl_vort,   & ! Flag for vorticity computation
          case_veloc, & ! Exact velocity
          case_press, & ! Exact pressure
          case_tempo, & ! Exact temporal solution
          case_t_pre    ! Exact temporal solution for pressure
     real(rp) ::         &
          react,         & ! Reaction
          diffu,         & ! Diffusion
          par_veloc(30), & ! Exact velocity field components (vector)
          par_press(10), & ! Exact pressure field components (scalar)
          par_tempo(3),  & ! Exact temporal field components (scalar)
          gravi(3)         ! Gravity field
   contains
     procedure :: create => nsi_create
     procedure :: free => nsi_free
  end type nsi_problem_t

  ! Types
  public :: nsi_problem_t

  ! Functions
  public :: nsi_elmchl, nsi_analytical_force

contains

  !=================================================================================================
  subroutine nsi_create(prob,ndime)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem                          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_problem_t), intent(inout) :: prob
    integer(ip)       , intent(in)    :: ndime
    integer(ip) :: i,istat

    ! Fill default problem data
    prob%ndime = ndime
    prob%ntens = prob%ndime*(prob%ndime+1)/2 ! Number of tensor components 
    prob%nunks = 2                           ! Velocity and pressure
    prob%nvars = prob%ndime+1                ! Number of degrees of freedom

    call memalloc(prob%nunks,prob%physical_problem_t%vars_of_unk,__FILE__,__LINE__)
    prob%vars_of_unk(1) = ndime
    prob%vars_of_unk(2) = 1

    !prob%vars_of_unk(2) = prob%vars_of_unk(1) + 1
    !prob%vars_of_unk(3) = prob%vars_of_unk(2) + 1

    ! allocate(prob%unkno_names(nunks),stat=istat)
    ! check(istat==0)
    ! prob%unkno_names(1)='velocity'
    ! prob%unkno_names(2)='pressure'

    ! Flags
    prob%ksnsi = 1    ! Symmetry flag
    prob%kfl_conv = 1 ! Enabling advection
    prob%kfl_vort = 0 ! Vorticity not computed
    prob%kfl_symg = 0 ! Symmetric grad-grad term (On=1, Off=0)
    prob%kfl_tder = 0 ! Time derivative not computed 
    prob%kfl_skew = 0 ! Enabling skewsymmetric convective terms: Off

    ! Problem variables
    prob%react  = 0.0_rp  ! Reaction
    prob%diffu  = 1.0_rp  ! Diffusion
    prob%gravi  = 0.0_rp  ! Gravity field

    ! Analytical field variables
    prob%case_veloc   = 0      ! Velocity field (see exact.f90)
    prob%case_press   = 0      ! Pressure field
    prob%case_tempo   = 0      ! Temporal field
    prob%case_t_pre   = 0      ! Temporal field for pressure
    prob%par_veloc(:) = 0.0_rp ! Exact velocity field
    prob%par_press(:) = 0.0_rp ! Exact pressure field
    prob%par_tempo(:) = 0.0_rp ! Exact temporal field

  end subroutine nsi_create

  !=================================================================================================
  subroutine nsi_free(prob)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates definitions of the Navier-Stokes problem                       !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_problem_t), intent(inout) :: prob

    call memfree ( prob%vars_of_unk,__FILE__,__LINE__)

  end subroutine nsi_free

  !==================================================================================================
  subroutine nsi_elmchl(jainv,hleng,elvel,ndime,pnode,iadve,chave,chale)
    !-----------------------------------------------------------------------------------------------!
    !    This routine computes the characteristic element lengths. The first one is the length in   ! 
    !    the flow direction if there are convective terms and the minimum length otherwise. The     ! 
    !    second one is the square root in 2D or cubic root in 3D of the volume.                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip), intent(in)    :: ndime,pnode,iadve
    real(rp),    intent(in)    :: jainv(ndime,ndime),hleng(ndime)
    real(rp),    intent(in)    :: elvel(pnode,ndime)
    real(rp),    intent(out)   :: chave(ndime,2)
    real(rp),    intent(out)   :: chale(2)
    integer(ip)                :: idime,ievab,inode,ivepo,kdime
    real(rp)                   :: elno1,elno2,hnatu

    ! Initialization
    chale(1)=hleng(ndime)                          ! Minimum element length
    !chale(2)=hleng(ndime)
    chale(2)=hleng(1)
    do idime=2,ndime
       chale(2)=chale(2)*hleng(idime)
    end do
    chale(2)=chale(2)**(1.0_rp/real(ndime))

    ! Natural element length
    if (ndime==2) then
       if(pnode==3.or.pnode==6.or.pnode==10) then  ! P1, P2 or P3
          hnatu=1.0_rp
       else if(pnode==4.or.pnode==9.or.pnode==16) then  ! Q1, Q2 or Q3
          hnatu=2.0_rp
       end if
    else if(ndime==3) then
       if(pnode==4.or.pnode==10.or.pnode==20) then  ! P1, P2 or P3
          hnatu=1.0_rp
       else if(pnode==8.or.pnode==27.or.pnode==64) then  ! Q1, Q2 or Q3
          hnatu=2.0_rp
       end if
    end if

    ! Length in the flow direction
    if(iadve==1) then
       chave=0.0_rp                                ! Characteristic element velocity
       do idime=1,ndime
          do inode=1,pnode
             chave(idime,1)=chave(idime,1)+elvel(inode,idime)
          end do
          chave(idime,1)=chave(idime,1)/real(pnode)
       end do
       do idime=1,ndime
          chave(idime,2)=0.0_rp
          do kdime=1,ndime
             chave(idime,2)=chave(idime,2)+jainv(idime,kdime)*chave(kdime,1)
          end do
       end do
       elno1=0.0_rp
       elno2=0.0_rp
       do idime=1,ndime
          elno1=elno1+chave(idime,1)*chave(idime,1)
       end do
       do idime=1,ndime
          elno2=elno2+chave(idime,2)*chave(idime,2)
       end do
       elno1=sqrt(elno1)
       elno2=sqrt(elno2)
       if(elno2>1.0e-16_rp.and.elno1>1.0e-16_rp) &
            chale(1)=elno1/elno2*hnatu              ! Characteristic element length
    end if

    !!! Only for CHANEL & NACA (h_min) 
    !chale(1)=hleng(ndime)
    !chale(2)=hleng(ndime)

    !!! (h_max) 
    !chale(1)=hleng(1)
    !chale(2)=hleng(1)

    ! Divide h by p^2
    if(ndime==2) then
       if(pnode>9) then                            ! Cubic elements
          chale = chale/(3.0_rp**2)
       else if(pnode>5) then                       ! Quadratic elements
          chale = (0.5_rp**2)*chale
       end if
    else if(ndime==3) then
       if(pnode==10.or.pnode==27) then
          chale = (0.5_rp**2)*chale                     ! Quadratic elements
       else if(pnode==20.or.pnode==64) then
          chale = chale/(3.0_rp**2)                     ! Cubic elements
       end if
    end if

  end subroutine nsi_elmchl

  !==================================================================================================
  subroutine nsi_analytical_force(physics,finite_element,ctime,gpvel,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine computes the elemental force needed to impose an analytical solution        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(nsi_problem_t)   , intent(in)    :: physics
    type(finite_element_t), intent(in)    :: finite_element
    real(rp)              , intent(in)    :: ctime
    type(vector_t)        , intent(in)    :: gpvel
    type(vector_t)        , intent(inout) :: force
    ! Locals
    integer(ip) :: igaus,idime,conve
    real(rp)    :: gpvno
    real(rp)    :: parv(30),parp(10),part(3),part_p(3)

    ! Loop over gauss points
    do igaus=1,finite_element%integ(1)%p%quad%ngaus
    
       ! Compute velocity norm
       gpvno=0.0_rp
       do idime = 1,physics%ndime
          gpvno = gpvno + gpvel%a(idime,igaus)*gpvel%a(idime,igaus)
       end do
       gpvno = sqrt(gpvno)

       ! Set convection flag
       if(gpvno.lt.1e-4) then
          conve = 0
       else
          conve = physics%kfl_conv
       end if
       
       ! Evaluate unknowns and derivatives
       parv   = 0.0_rp
       parp   = 0.0_rp
       part   = 0.0_rp
       part_p = 0.0_rp
       call analytical_field(physics%case_veloc,physics%ndime, &
            &                finite_element%integ(1)%p%femap%clocs(:,igaus),ctime,parv)
       call analytical_field(physics%case_press,physics%ndime, &
            &                finite_element%integ(1)%p%femap%clocs(:,igaus),ctime,parp)
       call analytical_field(physics%case_tempo,physics%ndime, &
            &                finite_element%integ(1)%p%femap%clocs(:,igaus),ctime,part)
       call analytical_field(physics%case_t_pre,physics%ndime, &
            &                finite_element%integ(1)%p%femap%clocs(:,igaus),ctime,part_p)

       ! Evaluate force
       call nsi_force(physics%ndime,physics%diffu,physics%react, &
            &         conve,physics%kfl_symg,physics%case_tempo,        &
            &         parv,parp,part,physics%gravi,force%a(:,igaus),part_p)

    end do

  end subroutine nsi_analytical_force

  !=================================================================================================
  subroutine nsi_force(ndime,diffu,react,kfl_conv,kfl_symg,caset,parv,parp,part,gravi,force,part_p)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the elemental force needed to impose an analytical solution       !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)       , intent(in)    :: ndime,kfl_conv,kfl_symg,caset
    real(rp)          , intent(in)    :: diffu,react,parv(:),parp(:),part(:),gravi(3)
    real(rp)          , intent(inout) :: force(ndime)
    real(rp), optional, intent(in)    :: part_p(:)

    real(rp)   :: u,v,w,dpdx,dpdy,dpdz,d2udx,d2udy,d2udz
    real(rp)   :: d2vdx,d2vdy,d2vdz,d2wdx,d2wdy,d2wdz
    real(rp)   :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp)   :: d2udxy,d2udxz,d2vdxy,d2vdyz,d2wdxz,d2wdyz
    real(rp)   :: dudt,dvdt,dwdt

    !
    if(caset>0) then
       dudt = part(2)*parv(1); dvdt = part(2)*parv(2); dwdt = part(2)*parv(3)
    else
       dudt=0.0_rp; dvdt=0.0_rp; dwdt=0.0_rp
    end if
    !
    u = parv(1)*part(1); v = parv(2)*part(1); w = parv(3)*part(1)
    dudx = parv(4)*part(1); dudy = parv(5)*part(1); dudz = parv(6)*part(1)
    dvdx = parv(7)*part(1); dvdy = parv(8)*part(1); dvdz = parv(9)*part(1)
    dwdx = parv(10)*part(1); dwdy = parv(11)*part(1); dwdz = parv(12)*part(1)
    !
    d2udx = parv(13)*part(1); d2udy = parv(14)*part(1); d2udz = parv(15)*part(1)
    d2vdx = parv(16)*part(1); d2vdy = parv(17)*part(1); d2vdz = parv(18)*part(1)
    d2wdx = parv(19)*part(1); d2wdy = parv(20)*part(1); d2wdz = parv(21)*part(1)
    !
    d2udxy = parv(22)*part(1); d2udxz = parv(23)*part(1)
    d2vdxy = parv(25)*part(1); d2vdyz = parv(27)*part(1)
    d2wdxz = parv(29)*part(1); d2wdyz = parv(30)*part(1)
    !
    if(present(part_p)) then
       dpdx = parp(2)*part_p(1); dpdy = parp(3)*part_p(1); dpdz = parp(4)*part_p(1)
    else
       dpdx = parp(2); dpdy = parp(3); dpdz = parp(4)
    end if
    !

    force(1) = dudt - diffu*(d2udx+d2udy+d2udz) + dpdx + react*u + &
         REAL(kfl_conv)*(u*dudx + v*dudy + w*dudz) - gravi(1) -    &
         REAL(kfl_symg)*diffu*(d2udx+d2vdxy+d2wdxz)
    force(2) = dvdt - diffu*(d2vdx+d2vdy+d2vdz) + dpdy + react*v + &
         REAL(kfl_conv)*(u*dvdx + v*dvdy + w*dvdz) - gravi(2) -    &
         REAL(kfl_symg)*diffu*(d2udxy+d2vdy+d2wdyz)
    if(ndime==3) force(3) = dwdt - diffu*(d2wdx+d2wdy+d2wdz) + dpdz + react*w + &
         REAL(kfl_conv)*(u*dwdx + v*dwdy + w*dwdz) - gravi(3) -    &
         REAL(kfl_symg)*diffu*(d2udxz+d2vdyz+d2wdz)

  end subroutine nsi_force


end module nsi_names
