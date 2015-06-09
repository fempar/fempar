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
module cdr_stabilized_continuous_Galerkin_names
 use types
 use memor
 use array_names
 use problem_names
 use cdr_names
 use element_fields_names
 use element_tools_names
 use volume_integration_tools_names
 implicit none
 private 

 type, extends(discrete_problem) :: cdr_approximation
    type(cdr_problem), pointer :: physics
    integer(ip) ::   & 
         kfl_1vec,   & ! Flag for 1vec integration
         kfl_mtvc,   & ! Flag for matvec integration
         kfl_matr,   & ! Flag for matrix integration
         kfl_real,   & ! Flag for real integration
         kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
         kfl_lump,   & ! Flag for lumped mass submatrix
         kfl_stab,   & ! Flag for stabilization of the convective term (Off=0; OSS=2)
         kfl_proj,   & ! Flag for Projections weighted with tau's (On=1, Off=0)
         tdimv         ! Number of temporal steps stored
    real(rp) ::      &
         dtinv,      & ! Inverse of time step
         ctime,      & ! Current time 
         k1tau,      & ! C1 constant on stabilization parameter 
         k2tau         ! C2 constant on stabilization parameter 
    contains
      procedure :: create => cdr_stabilized_continuous_Galerkin_create
      procedure :: matvec => cdr_stabilized_continuous_Galerkin_matvec
   end type cdr_approximation

 public :: cdr_approximation

contains

  !=================================================================================================
  subroutine cdr_stabilized_continuous_Galerkin_create(approx,prob)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(cdr_approximation)       , intent(out) :: approx
    class(physical_problem), target, intent(in)  :: prob
    !type(cdr_problem), target, intent(in)  :: prob
    integer(ip) :: i


    select type (prob)
    type is(cdr_problem)
       approx%physics => prob
       !approx%cdr => prob
    class default
       check(.false.)
    end select

    ! Flags
    approx%kfl_1vec = 1 ! Integrate Full OSS
    approx%kfl_mtvc = 1 ! Integrate cdr_element_matvec
    approx%kfl_matr = 1 ! Integrate cdr_element_mat
    approx%kfl_real = 0 ! Integrate errornorm
    approx%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)
    approx%kfl_stab = 0 ! Stabilization of convective term (0: Off, 2: OSS)
    approx%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    approx%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter
    approx%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter

    ! Time integration variables
    approx%dtinv  = 1.0_rp ! Inverse of time step
    approx%ctime  = 0.0_rp ! Current time
    approx%tdimv  =  2     ! Number of temporal steps stored
    
  end subroutine cdr_stabilized_continuous_Galerkin_create

  subroutine cdr_stabilized_continuous_Galerkin_matvec(approx,integ,unkno,start,mat,vec)
    implicit none
    class(cdr_approximation)       , intent(in) :: approx
    type(volume_integrator_pointer), intent(in) :: integ(:)
    real(rp)                       , intent(in) :: unkno(:,:,:)
    integer(ip)                    , intent(in) :: start(:)
    type(array_rp2), intent(inout) :: mat
    type(array_rp1), intent(inout) :: vec

    type(basis_function) :: u ! Trial
    type(basis_function) :: p 
    type(basis_function) :: v ! Test
    type(basis_function) :: q 

    type(vector) :: a
    type(vector) :: u_n
    type(scalar) :: h,tau

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus,nnode
    real(rp) :: factor


    mat%a = 0.0_rp
    vec%a = 0.0_rp
    ndime = approx%physics%ndime

    !u = basis_function(approx%physics,1,start,integ)
    !p = basis_function(approx%physics,2,start,integ)
    !v = basis_function(approx%physics,1,start,integ) 
    !q = basis_function(approx%physics,2,start,integ)

    ! The fields can be created once and reused on each element
    ! To do that we require an initial loop over elements and
    ! an initialization call.
    !call create_vector (approx%physics, 1, integ, a)
    !call create_vector (approx%physics, 1, integ, u_n)
    ! Then for each element fill values
    !call interpolation (unkno, 1, 1, integ, a)
    !call interpolation (unkno, 1, 3, integ, u_n)

    ! Dirty, isnt'it?
    !h%a=integ(1)%p%femap%hleng(1,:)     ! max
    !h%a=integ(1)%p%femap%hleng(ndime,:) ! min

    !mu = approx%physics%diffu

    !dtinv  = approx%dtinv
    !c1 = approx%k1tau
    !c2 = approx%k1tau

    ! tau = c1*mu*inv(h*h) + c2*norm(a)*inv(h)
    !tau = inv(tau)

    !mat = integral(grad(v),grad(u))

    nnode = integ(1)%p%uint_phy%nnode
    ngaus = integ(1)%p%uint_phy%nlocs

    do igaus = 1,ngaus
       factor = integ(1)%p%femap%detjm(igaus) * integ(1)%p%quad%weight(igaus)
       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                mat%a(inode,jnode) = mat%a(inode,jnode) + factor * &
                     & integ(1)%p%uint_phy%deriv(idime,inode,igaus) * &
                     & integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
          end do
       end do
    end do
 
    !mat = integral(v,dtinv*u) + integral(grad(v),grad(u))
    !mat = integral(v,dtinv*u) + integral(v, a*grad(u)) + integral(grad(v),mu*grad(u)) + integral(a*grad(v),tau*a*grad(u)) + integral(div(v),p) + integral(q,div(u))

    ! This will not work right now becaus + of basis_functions and gradients is not defined.
    !mat = integral(v,dtinv*u+a*grad(u)) + mu*integral(grad(v),grad(u)) + integral(a*grad(v),tau*a*grad(u) ) + integral(div(v),p) + integral(q,div(u))

  end subroutine cdr_stabilized_continuous_Galerkin_matvec

end module cdr_stabilized_continuous_Galerkin_names
