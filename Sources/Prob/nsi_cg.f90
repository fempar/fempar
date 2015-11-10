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
module nsi_cg_asgs_names
 use types_names
 use memor_names
 use problem_names
 use nsi_names
 use element_fields_names
 use element_tools_names
 use finite_element_names
 implicit none
 private 

 type, extends(discrete_problem_t) :: nsi_cg_asgs_discrete_t
    integer(ip) ::   & 
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
      procedure :: create => nsi_create_discrete
   end type nsi_cg_asgs_discrete_t

 type, extends(discrete_integration_t) :: nsi_cg_asgs_approximation_t
      type(nsi_cg_asgs_discrete_t), pointer :: discret
      type(nsi_problem_t)         , pointer :: physics
    contains
      procedure :: create  => nsi_matvec_create
      procedure :: compute => nsi_matvec
      procedure :: free    => nsi_matvec_free
   end type nsi_cg_asgs_approximation_t

 public :: nsi_cg_asgs_approximation_t, nsi_cg_asgs_discrete_t

contains

  !=================================================================================================
  subroutine nsi_create_discrete(discret,physics,l2g)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_asgs_discrete_t), intent(out) :: discret
    class(physical_problem_t)    , intent(in)  :: physics
    integer(ip), optional      , intent(in)  :: l2g(:)
    ! Locals
    integer(ip) :: i

    ! Flags
    discret%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)
    discret%kfl_stab = 0 ! Stabilization of convective term (0: Off, 2: OSS)
    discret%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    discret%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter tau_m
    discret%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    discret%ktauc  = 0.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Time integration variables
    discret%dtinv  = 1.0_rp ! Inverse of time step
    discret%ctime  = 0.0_rp ! Current time
    discret%tdimv  =  2     ! Number of temporal steps stored for velocity
    discret%tdimp  =  2     ! Number of temporal steps stored for pressure

    discret%nvars = physics%ndime+1
    call memalloc(discret%nvars,discret%l2g_var,__FILE__,__LINE__)
    if ( present(l2g) ) then
       assert ( size(l2g) == discret%nvars )
       discret%l2g_var = l2g
    else
       do i = 1,discret%nvars
          discret%l2g_var(i) = i
       end do
    end if
    
  end subroutine nsi_create_discrete

  !=================================================================================================
  subroutine nsi_matvec_create( this, physics, discret )
    implicit none
    class(nsi_cg_asgs_approximation_t), intent(inout) :: this
    class(physical_problem_t), target , intent(in)    :: physics
    class(discrete_problem_t), target , intent(in)    :: discret

    select type (physics)
    type is(nsi_problem_t)
       this%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_asgs_discrete_t)
       this%discret => discret
       class default
       check(.false.)
    end select

	this%domain_dimension = 3
  end subroutine nsi_matvec_create

  !=================================================================================================
  subroutine nsi_matvec_free(this)
    implicit none
    class(nsi_cg_asgs_approximation_t)  , intent(inout) :: this

    this%physics => null()
    this%discret => null()

  end subroutine nsi_matvec_free

  !=================================================================================================
  subroutine nsi_matvec(this,finite_element)
    implicit none
    class(nsi_cg_asgs_approximation_t), intent(inout) :: this
    type(finite_element_t)            , intent(inout) :: finite_element

    type(basis_function_t) :: u ! Trial
    type(basis_function_t) :: p 
    type(basis_function_t) :: v ! Test
    type(basis_function_t) :: q 

    type(vector_t) :: a
    type(vector_t) :: u_n
    type(scalar_t) :: h,tau

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu

    ndime = this%physics%ndime
       
    u = basis_function_t(this%physics,1,finite_element%start%a,finite_element%integ)
    p = basis_function_t(this%physics,2,finite_element%start%a,finite_element%integ)
    v = basis_function_t(this%physics,1,finite_element%start%a,finite_element%integ) 
    q = basis_function_t(this%physics,2,finite_element%start%a,finite_element%integ)

    ! With a_h declared as given_function we could do:
    ! a_h   = given_function(approx,1,1,finite_element%integ)
    ! call interpolation(grad(a_h),grad_a) 
    ! with grad_a declared as tensor and, of course
    ! call interpolation(grad(given_function(approx,1,1,finite_element%integ)),grad_a)
    ! is also possible.
    !call interpolation(given_function(approx,1,1,finite_element%integ),a)
    !call interpolation(given_function(approx,1,3,finite_element%integ),u_n)

    ! The fields can be created once and reused on each element
    ! To do that we require an initial loop over elements and
    ! an initialization call.
    call create_vector (this%physics, 1, finite_element%integ, a)
    call create_vector (this%physics, 1, finite_element%integ, u_n)
    ! Then for each element fill values
    call interpolation (finite_element%unkno, 1, 1, finite_element%integ, a)
    call interpolation (finite_element%unkno, 1, 3, finite_element%integ, u_n)

    ! Dirty, isnt'it?
    !h%a=finite_element%integ(1)%p%femap%hleng(1,:)     ! max
    h%a=finite_element%integ(1)%p%femap%hleng(ndime,:) ! min

    
    mu = this%physics%diffu

    dtinv  = this%discret%dtinv
    c1 = this%discret%k1tau
    c2 = this%discret%k1tau

    ! tau = c1*mu*inv(h*h) + c2*norm(a)*inv(h)
    tau = inv(tau)
    
    finite_element%p_mat = integral(v,dtinv*u) + integral(grad(v),grad(u))
    !finite_element%mat = integral(v,dtinv*u) + integral(v, a*grad(u)) + integral(grad(v),mu*grad(u)) + integral(a*grad(v),tau*a*grad(u)) + integral(div(v),p) + integral(q,div(u))

    ! This will not work right now becaus + of basis_functions and gradients is not defined.
    !finite_element%mat = integral(v,dtinv*u+a*grad(u)) + mu*integral(grad(v),grad(u)) + integral(a*grad(v),tau*a*grad(u) ) + integral(div(v),p) + integral(q,div(u))

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine nsi_matvec

end module nsi_cg_asgs_names
