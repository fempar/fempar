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
module nsi_cg
 use types
 use memor
 use array_names
 use problem_names
 use nsi_names
 use element_fields_names
 use element_tools_names
 use volume_integration_tools_names
 implicit none
 private 

 type, extends(discrete_problem) :: nsi_approximation
    type(nsi_problem), pointer :: nsi
    integer(ip) ::   & 
         kfl_1vec,   & ! Flag for 1vec integration
         kfl_mtvc,   & ! Flag for matvec integration
         kfl_matr,   & ! Flag for matrix integration
         kfl_real,   & ! Flag for real integration
         kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
         kfl_lump,   & ! Flag for lumped mass submatrix
         kfl_stab,   & ! Flag for stabilization of the convective term (Off=0; OSS=2)
         kfl_proj,   & ! Flag for Projections weighted with tau's (On=1, Off=0)
         tdimv,      & ! Number of temporal steps stored for velocity
         tdimp         ! Number of temporal steps stored for pressure   
    real(rp) ::      &
         dtinv,      & ! Inverse of time step
         ctime,      & ! Current time
         ktauc,      & ! Constant multiplying stabilization parameter tau_c 
         k1tau,      & ! C1 constant on stabilization parameter tau_m
         k2tau         ! C2 constant on stabilization parameter tau_m
    contains
      procedure :: create => nsi_create
      procedure :: matvec => nsi_matvec
   end type nsi_approximation

 public :: nsi_approximation

contains

  !=================================================================================================
  subroutine nsi_create(aprox,prob)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_approximation)       , intent(out) :: aprox
    class(physical_problem), target, intent(in)  :: prob
    integer(ip) :: i

    select type (prob)
    type is(nsi_problem)
       aprox%nsi => prob
    class default
       check(.false.)
    end select

    ! Flags
    aprox%kfl_1vec = 1 ! Integrate Full OSS
    aprox%kfl_mtvc = 1 ! Integrate nsi_element_matvec
    aprox%kfl_matr = 1 ! Integrate nsi_element_mat
    aprox%kfl_real = 0 ! Integrate errornorm
    aprox%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)
    aprox%kfl_stab = 0 ! Stabilization of convective term (0: Off, 2: OSS)
    aprox%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    aprox%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter tau_m
    aprox%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    aprox%ktauc  = 0.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Time integration variables
    aprox%dtinv  = 1.0_rp ! Inverse of time step
    aprox%ctime  = 0.0_rp ! Current time
    aprox%tdimv  =  2     ! Number of temporal steps stored for velocity
    aprox%tdimp  =  2     ! Number of temporal steps stored for pressure
    
  end subroutine nsi_create

  subroutine nsi_matvec(aprox,integ,unkno,mat,vec)
    implicit none
    class(nsi_approximation)       , intent(in) :: aprox
    type(volume_integrator_pointer), intent(in) :: integ(:)
    real(rp)                       , intent(in) :: unkno(:,:,:)
    type(array_rp2), intent(inout) :: mat
    type(array_rp1), intent(inout) :: vec

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

    ndime = aprox%nsi%ndime

    u = basis_function(aprox%nsi,1,integ)
    p = basis_function(aprox%nsi,2,integ)
    v = basis_function(aprox%nsi,1,integ) 
    q = basis_function(aprox%nsi,2,integ)

    ! With a_h declared as given_function we could do:
    ! a_h   = given_function(aprox,1,1,integ)
    ! call interpolation(grad(a_h),grad_a) 
    ! with grad_a declared as tensor and, of course
    ! call interpolation(grad(given_function(aprox,1,1,integ)),grad_a)
    ! is also possible.
    !call interpolation(given_function(aprox,1,1,integ),a)
    !call interpolation(given_function(aprox,1,3,integ),u_n)

    ! The fields can be created once and reused on each element
    ! To do that we require an initial loop over elements and
    ! an initialization call.
    call create_vector (aprox%nsi, 1, integ, a)
    call create_vector (aprox%nsi, 1, integ, u_n)
    ! Then for each element fill values
    call interpolation (unkno, 1, 1, integ, a)
    call interpolation (unkno, 1, 3, integ, u_n)

    ! Dirty, isnt'it?
    !h%a=integ(1)%p%femap%hleng(1,:)     ! max
    h%a=integ(1)%p%femap%hleng(ndime,:) ! min

    mu = aprox%nsi%diffu

    dtinv  = aprox%dtinv
    c1 = aprox%k1tau
    c2 = aprox%k1tau

    tau = c1*mu*inv(h*h) + c2*norm(a)*inv(h)
    tau = inv(tau)
    
    mat = integral(v,dtinv*u)
    !mat = integral(v,dtinv*u) + integral(v, a*grad(u)) + integral(grad(v),mu*grad(u)) + integral(a*grad(v),tau*a*grad(u)) + integral(div(v),p) + integral(q,div(u))

    ! This will not work right now becaus + of basis_functions and gradients is not defined.
    !mat = integral(v,dtinv*u+a*grad(u)) + mu*integral(grad(v),grad(u)) + integral(a*grad(v),tau*a*grad(u) ) + integral(div(v),p) + integral(q,div(u))

  end subroutine nsi_matvec

end module nsi_cg
