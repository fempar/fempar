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
  use analytical_function_names
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
    real(rp)    :: params(11,physics%nvars)

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
       params = 0.0_rp
       do idime=1,physics%nvars
          call evaluate_analytical(finite_element%p_analytical_code%a(idime,1),                  &
               &                   finite_element%p_analytical_code%a(idime,2),                  &
               &                   physics%ndime,finite_element%integ(1)%p%femap%clocs(:,igaus), &
               &                   ctime,params(:,idime))
       end do

       ! Evaluate force
       call nsi_force(physics%ndime,physics%diffu,physics%react,conve,physics%kfl_symg,params, &
            &         physics%gravi,force%a(:,igaus))

    end do

  end subroutine nsi_analytical_force

  !=================================================================================================
  subroutine nsi_force(ndime,diffu,react,kfl_conv,kfl_symg,params,gravi,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the elemental force needed to impose an analytical solution       !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)       , intent(in)    :: ndime,kfl_conv,kfl_symg
    real(rp)          , intent(in)    :: diffu,react,params(11,ndime+1),gravi(3)
    real(rp)          , intent(inout) :: force(ndime)

    real(rp)   :: u,v,w,dpdx,dpdy,dpdz,d2udx,d2udy,d2udz
    real(rp)   :: d2vdx,d2vdy,d2vdz,d2wdx,d2wdy,d2wdz
    real(rp)   :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp)   :: d2udxy,d2udxz,d2udyz,d2vdxy,d2vdxz,d2vdyz,d2wdxy,d2wdxz,d2wdyz
    real(rp)   :: dudt,dvdt,dwdt

    !
    dudt = params(11,1); dvdt = params(11,2); dwdt = params(11,3)
    !
    u = params(1,1); v = params(1,2)
    !
    dudx = params(2,1); dvdx = params(2,2)
    dudy = params(3,1); dvdy = params(3,2)
    dudz = params(4,1); dvdz = params(4,2)
    !
    d2udx = params(5,1); d2vdx = params(5,2)
    d2udy = params(6,1); d2vdy = params(6,2)
    d2udz = params(7,1); d2vdz = params(7,2)
    !
    d2udxy = params(8,1);  d2vdxy = params(8,2)  
    d2udxz = params(9,1);  d2vdxz = params(9,2)  
    d2udyz = params(10,1); d2vdyz = params(10,2)
    !
    if(ndime==3) then
       w = params(1,3)
       dwdx = params(2,3); d2wdx = params(5,3); d2wdxy = params(8,3)    
       dwdy = params(3,3); d2wdy = params(6,3); d2wdxz = params(9,3) 
       dwdz = params(4,3); d2wdz = params(7,3); d2wdyz = params(10,3)
    else
       w = 0.0_rp
       dwdx = 0.0_rp; d2wdx = 0.0_rp; d2wdxy = 0.0_rp   
       dwdy = 0.0_rp; d2wdy = 0.0_rp; d2wdxz = 0.0_rp
       dwdz = 0.0_rp; d2wdz = 0.0_rp; d2wdyz = 0.0_rp       
    end if
    !
    dpdx = params(2,ndime+1); dpdy = params(3,ndime+1); dpdz = params(4,ndime+1)
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
