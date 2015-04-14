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
module nsi
 use types
 use memor
 use fem_space_names
 use problem_names
 use element_fields_names
 use element_tools_names
 implicit none
 private 

 type, extends(physical_problem) :: nsi_problem

    integer(ip) ::   & 
         ksnsi,      & ! Symmetry flag (+-1) (1: no simetrica definida, -1: simetrica no definida)
         kfl_conv,   & ! Flag for enabling advection (Stokes=0; NavSto=1)
         kfl_1vec,   & ! Flag for 1vec integration
         kfl_mtvc,   & ! Flag for matvec integration
         kfl_matr,   & ! Flag for matrix integration
         kfl_real,   & ! Flag for real integration
         kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
         kfl_vort,   & ! Flag for vorticity computation
         kfl_lump,   & ! Flag for lumped mass submatrix
         kfl_symg,   & ! Flag for symmetric grad-grad term
         kfl_tder,   & ! Flag for time derivative computation
         kfl_skew,   & ! Flag for enabling skewsymmetric convective terms (Off=0; type1=1, type2=2)
         kfl_stab,   & ! Flag for stabilization of the convective term (Off=0; OSS=2)
         kfl_proj,   & ! Flag for Projections weighted with tau's (On=1, Off=0)
         case_veloc, & ! Exact velocity
         case_press, & ! Exact pressure
         case_tempo, & ! Exact temporal solution
         case_t_pre, & ! Exact temporal solution for pressure
         tdimv,      & ! Number of temporal steps stored for velocity
         tdimp         ! Number of temporal steps stored for pressure   

    real(rp) ::         &
         react,         & ! Reaction
         diffu,         & ! Diffusion
         dtinv,         & ! Inverse of time step
         ctime,         & ! Current time
         par_veloc(30), & ! Exact velocity field components (vector)
         par_press(10), & ! Exact pressure field components (scalar)
         par_tempo(3),  & ! Exact temporal field components (scalar)
         gravi(3),      & ! Gravity field
         ktauc,         & ! Constant multiplying stabilization parameter tau_c 
         k1tau,         & ! C1 constant on stabilization parameter tau_m
         k2tau            ! C2 constant on stabilization parameter tau_m

    ! type(array_rp2), allocatable :: &
    !      vorti(:),      & ! Vorticity field
    !      dtunk(:)         ! Time derivative
    contains
      procedure :: create => nsi_create
      procedure :: matvec => nsi_matvec

 end type nsi_problem

 public :: nsi_problem

contains

  !=================================================================================================
  subroutine nsi_create(prob,ndime)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_problem), intent(out) :: prob
    integer(ip)      , intent(in)  :: ndime
    integer(ip) :: i

    ! Fill default problem data
    prob%ndime = ndime
    prob%ntens = prob%ndime*(prob%ndime+1)/2 ! Number of tensor components 
    prob%nunks = 2                           ! Velocity and pressure
    prob%nvars = prob%ndime+1                ! Number of degrees of freedom

    call memalloc(prob%nunks,prob%vars_of_unk,__FILE__,__LINE__)
    prob%vars_of_unk(1) = ndime
    prob%vars_of_unk(2) = 1

    !prob%vars_of_unk(2) = prob%vars_of_unk(1) + 1
    !prob%vars_of_unk(3) = prob%vars_of_unk(2) + 1

    ! Should be overwritten by the driver in a multiphyiscs context
    call memalloc(prob%nvars,prob%l2g_var,__FILE__,__LINE__)
    do i = 1,prob%nvars
       prob%l2g_var(i) = i
    end do


    ! Flags
    prob%ksnsi = 1    ! Symmetry flag
    prob%kfl_conv = 1 ! Enabling advection
    prob%kfl_1vec = 1 ! Integrate Full OSS
    prob%kfl_mtvc = 1 ! Integrate nsi_element_matvec
    prob%kfl_matr = 1 ! Integrate nsi_element_mat
    prob%kfl_real = 0 ! Integrate errornorm
    prob%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)
    prob%kfl_vort = 0 ! Vorticity not computed
    prob%kfl_symg = 0 ! Symmetric grad-grad term (On=1, Off=0)
    prob%kfl_tder = 0 ! Time derivative not computed 
    prob%kfl_skew = 0 ! Enabling skewsymmetric convective terms: Off
    prob%kfl_stab = 0 ! Stabilization of convective term (0: Off, 2: OSS)
    prob%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    prob%react  = 0.0_rp  ! Reaction
    prob%diffu  = 1.0_rp  ! Diffusion
    prob%gravi  = 0.0_rp  ! Gravity field
    prob%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter tau_m
    prob%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    prob%ktauc  = 0.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Analytical field variables
    prob%case_veloc   = 0      ! Velocity field (see exact.f90)
    prob%case_press   = 0      ! Pressure field
    prob%case_tempo   = 0      ! Temporal field
    prob%case_t_pre   = 0      ! Temporal field for pressure
    prob%par_veloc(:) = 0.0_rp ! Exact velocity field
    prob%par_press(:) = 0.0_rp ! Exact pressure field
    prob%par_tempo(:) = 0.0_rp ! Exact temporal field
    
    ! Time integration variables
    prob%dtinv  = 1.0_rp ! Inverse of time step
    prob%ctime  = 0.0_rp ! Current time
    prob%tdimv  =  2     ! Number of temporal steps stored for velocity
    prob%tdimp  =  2     ! Number of temporal steps stored for pressure
    
  end subroutine nsi_create

  subroutine nsi_matvec(prob,K)
    ! Thinking about e.g. iss + su
    implicit none
    class(nsi_problem), intent(inout) :: prob
    type(fem_element) , intent(inout) :: K
    !real(rp), intent(inout) :: mat(:,:)
    !real(rp), intent(inout) :: vec(:)

    type(basis_function) :: u ! Trial
    type(basis_function) :: p 
    type(basis_function) :: v ! Test
    type(basis_function) :: q 
    !type(given_function) :: a_h
    type(vector) :: a
    type(vector) :: u_n
    type(scalar) :: h,tau

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu
    ndime = prob%ndime

    u = basis_function(prob,1,K)
    p = basis_function(prob,2,K)
    v = basis_function(prob,1,K) 
    q = basis_function(prob,2,K)

    ! With a_h declared as given_function we could do:
    ! a_h   = given_function(prob,1,1,K)
    ! call interpolation(grad(a_h),grad_a) 
    ! with grad_a declared as tensor and, of course
    ! call interpolation(grad(given_function(prob,1,1,K)),grad_a)
    ! is also possible.
    call interpolation(given_function(prob,1,1,K),a)
    call interpolation(given_function(prob,1,3,K),u_n)

    ! Dirty, isnt'it?
    !h%a=K%integ(u%ivar)%p%femap%hleng(1,:)     ! max
    h%a=K%integ(u%ivar)%p%femap%hleng(ndime,:) ! min

    dtinv  = prob%dtinv
    mu = prob%diffu
    c1 = prob%k1tau
    c2 = prob%k1tau

    tau = c1*mu*inv(h*h) + c2*norm(a)*inv(h)
    tau = inv(tau)
    
    K%p_mat = integral(v,dtinv*u)
    mat = integral(v,dtinv*u) + integral(v, a*grad(u)) + integral(grad(v),mu*grad(u)) + integral(a*grad(v),tau*a*grad(u)) + integral(div(v),p) + integral(q,div(u))

    ! This will not work right now becaus + of basis_functions and gradients is not defined.
    !mat = integral(v,dtinv*u+a*grad(u)) + mu*integral(grad(v),grad(u)) + integral(a*grad(v),tau*a*grad(u) ) + integral(div(v),p) + integral(q,div(u))

  end subroutine nsi_matvec

end module nsi
