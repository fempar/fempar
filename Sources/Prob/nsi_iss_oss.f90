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
module nsi_cg_iss_oss_names
  use types_names
  use memor_names
  use problem_names
  use nsi_names
  use finite_element_names
  use eltrm_gen_names
  use element_fields_names
  use element_tools_names
  use rungekutta_names
  use time_integration_names
  implicit none
# include "debug.i90"
  
  private
  
  ! INF-SUP STABLE (iss) with Orthogonal SubScales (oss) NAVIER-STOKES types 
  ! Problem data
  type, extends(discrete_problem_t) :: nsi_cg_iss_oss_discrete_t
     integer(ip) ::   & 
          kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
          kfl_lump,   & ! Flag for lumped mass submatrix
          kfl_proj,   & ! Flag for Projections weighted with tau's (On=1, Off=0)
          tdimv,      & ! Number of temporal steps stored for velocity
          tdimp         ! Number of temporal steps stored for pressure   
     real(rp) ::         &
          ktauc,         & ! Constant multiplying stabilization parameter tau_c 
          k1tau,         & ! C1 constant on stabilization parameter tau_m
          k2tau            ! C2 constant on stabilization parameter tau_m
   contains
     procedure :: create  => nsi_create_discrete
     procedure :: vars_block => nsi_vars_block
     procedure :: dof_coupling => nsi_dof_coupling
  end type nsi_cg_iss_oss_discrete_t

  ! Matvec
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_matvec_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(time_integration_t)       , pointer :: tinteg
   contains
     procedure :: create  => nsi_matvec_create
     procedure :: compute => nsi_matvec 
     procedure :: free    => nsi_matvec_free
  end type nsi_cg_iss_oss_matvec_t

  ! Mass_velocity
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_massu_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
   contains
     procedure :: create  => nsi_massu_create
     procedure :: compute => nsi_massu 
     procedure :: free    => nsi_massu_free
  end type nsi_cg_iss_oss_massu_t

  ! Mass_pressure
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_massp_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
   contains
     procedure :: create  => nsi_massp_create
     procedure :: compute => nsi_massp 
     procedure :: free    => nsi_massp_free
  end type nsi_cg_iss_oss_massp_t

  ! Mass_projection
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_massx_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
   contains
     procedure :: create  => nsi_massx_create
     procedure :: compute => nsi_massx 
     procedure :: free    => nsi_massx_free
  end type nsi_cg_iss_oss_massx_t

  ! Laplacian_pressure
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_lapla_p_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
   contains
     procedure :: create  => nsi_lapla_p_create
     procedure :: compute => nsi_lapla_p 
     procedure :: free    => nsi_lapla_p_free
  end type nsi_cg_iss_oss_lapla_p_t

  ! Runge-Kutta momentum equation (LHS)
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_rk_momentum_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(rungekutta_integrator_t)  , pointer :: rkinteg
   contains
     procedure :: create  => nsi_rk_momentum_create
     procedure :: compute => nsi_rk_momentum
     procedure :: free    => nsi_rk_momentum_free
  end type nsi_cg_iss_oss_rk_momentum_t

  ! Runge-Kutta momentum equation (RHS)
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_rk_momentum_rhs_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(rungekutta_integrator_t)  , pointer :: rkinteg
   contains
     procedure :: create  => nsi_rk_momentum_rhs_create
     procedure :: compute => nsi_rk_momentum_rhs
     procedure :: free    => nsi_rk_momentum_rhs_free
  end type nsi_cg_iss_oss_rk_momentum_rhs_t

  ! Runge-Kutta pressure equation
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_rk_pressure_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(time_integration_t)       , pointer :: tinteg
   contains
     procedure :: create  => nsi_rk_pressure_create
     procedure :: compute => nsi_rk_pressure
     procedure :: free    => nsi_rk_pressure_free
  end type nsi_cg_iss_oss_rk_pressure_t

  ! Runge-Kutta momentum equation update (U_n+1)
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_rk_momentum_update_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(rungekutta_integrator_t)  , pointer :: rkinteg
   contains
     procedure :: create  => nsi_rk_momentum_update_create
     procedure :: compute => nsi_rk_momentum_update
     procedure :: free    => nsi_rk_momentum_update_free
  end type nsi_cg_iss_oss_rk_momentum_update_t

  ! Runge-Kutta momentum equation update (U_n+1)
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_rk_projection_update_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
   contains
     procedure :: create  => nsi_rk_projection_update_create
     procedure :: compute => nsi_rk_projection_update
     procedure :: free    => nsi_rk_projection_update_free
  end type nsi_cg_iss_oss_rk_projection_update_t

  ! Dissipations
  type, extends(discrete_integration_t) :: nsi_cg_iss_oss_dissipation_t
     type(nsi_cg_iss_oss_discrete_t), pointer :: discret
     type(nsi_problem_t)            , pointer :: physics
     type(time_integration_t)       , pointer :: tinteg
     real(rp)                                 :: values(10)
   contains
     procedure :: create  => nsi_dissipation_create
     procedure :: compute => nsi_dissipation
     procedure :: free    => nsi_dissipation_free
  end type nsi_cg_iss_oss_dissipation_t

  ! Unkno components parameter definition
  integer(ip), parameter :: current   = 1
  integer(ip), parameter :: prev_iter = 2
  integer(ip), parameter :: prev_step = 3
  integer(ip), parameter :: explicit = 0, implicit=1

  ! Types
  public :: nsi_cg_iss_oss_matvec_t, nsi_cg_iss_oss_discrete_t, nsi_cg_iss_oss_massp_t, &
       &    nsi_cg_iss_oss_lapla_p_t, nsi_cg_iss_oss_rk_momentum_t, nsi_cg_iss_oss_rk_pressure_t, &
       &    nsi_cg_iss_oss_rk_momentum_update_t, nsi_cg_iss_oss_rk_projection_update_t, &
       &    nsi_cg_iss_oss_massu_t, nsi_cg_iss_oss_massx_t, nsi_cg_iss_oss_rk_momentum_rhs_t, &
       &    nsi_cg_iss_oss_dissipation_t
  
contains

  !=================================================================================================
  subroutine nsi_create_discrete(discret,physics,l2g)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets and orthogonal subscales.           !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_discrete_t), intent(out) :: discret
    class(physical_problem_t)       , intent(in)  :: physics
    integer(ip), optional           , intent(in)  :: l2g(:)
    ! Locals
    integer(ip) :: i

    ! Flags
    discret%kfl_lump = 0 ! Flag for lumped mass submatrix (Off=0, On=1)
    discret%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=1)
    discret%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    discret%k1tau = 12.0_rp ! C1 constant on stabilization parameter tau_m
    discret%k2tau = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    discret%ktauc = 4.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Time integration variables
    discret%tdimv = 2     ! Number of temporal steps stored for velocity
    discret%tdimp = 2     ! Number of temporal steps stored for pressure

    discret%nvars = 2*physics%ndime+1
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
  subroutine nsi_matvec_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_matvec_t)   , intent(inout) :: approx
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(discret%nvars,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,discret%nvars
       approx%working_vars(ivar) = ivar
    end do

  end subroutine nsi_matvec_create

  !=================================================================================================
  subroutine nsi_matvec_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_matvec_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_matvec_free

  !=================================================================================================
  subroutine nsi_massu_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massu_t)    , intent(inout) :: approx
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = ivar
    end do

  end subroutine nsi_massu_create

  !=================================================================================================
  subroutine nsi_massu_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massu_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_massu_free

  !=================================================================================================
  subroutine nsi_massp_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massp_t)    , intent(inout) :: approx
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(1,approx%working_vars,__FILE__,__LINE__)
    approx%working_vars(1) = physics%ndime+1

  end subroutine nsi_massp_create

  !=================================================================================================
  subroutine nsi_massp_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massp_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_massp_free

  !=================================================================================================
  subroutine nsi_massx_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massx_t)    , intent(inout) :: approx
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = physics%ndime + 1 + ivar
    end do

  end subroutine nsi_massx_create

  !=================================================================================================
  subroutine nsi_massx_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massx_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_massx_free

  !=================================================================================================
  subroutine nsi_lapla_p_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_lapla_p_t)  , intent(inout) :: approx
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(1,approx%working_vars,__FILE__,__LINE__)
    approx%working_vars(1) = physics%ndime+1

  end subroutine nsi_lapla_p_create

  !=================================================================================================
  subroutine nsi_lapla_p_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_lapla_p_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_lapla_p_free

  !=================================================================================================
  subroutine nsi_rk_momentum_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_t), intent(inout) :: approx
    class(physical_problem_t)  , target, intent(in)    :: physics
    class(discrete_problem_t)  , target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select 

    ! Domain dimension
    approx%domain_dimension = 3

   ! Allocate working variables
    call memalloc(2*physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = ivar
    end do
    do ivar=1,physics%ndime
       approx%working_vars( physics%ndime + ivar) = physics%ndime + 1 + ivar
    end do

  end subroutine nsi_rk_momentum_create

  !=================================================================================================
  subroutine nsi_rk_momentum_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_rk_momentum_free

  !=================================================================================================
  subroutine nsi_rk_momentum_rhs_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_rhs_t), intent(inout) :: approx
    class(physical_problem_t)      , target, intent(in)    :: physics
    class(discrete_problem_t)      , target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

   ! Allocate working variables
    call memalloc(2*physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = ivar
    end do
    do ivar=1,physics%ndime
       approx%working_vars( physics%ndime + ivar) = physics%ndime + 1 + ivar
    end do

  end subroutine nsi_rk_momentum_rhs_create

  !=================================================================================================
  subroutine nsi_rk_momentum_rhs_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_rhs_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_rk_momentum_rhs_free

  !=================================================================================================
  subroutine nsi_rk_pressure_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_pressure_t), intent(inout) :: approx
    class(physical_problem_t)  , target, intent(in)    :: physics
    class(discrete_problem_t)  , target, intent(in)    :: discret

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(1,approx%working_vars,__FILE__,__LINE__)
    approx%working_vars(1) = physics%ndime+1

  end subroutine nsi_rk_pressure_create

  !=================================================================================================
  subroutine nsi_rk_pressure_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_pressure_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_rk_pressure_free

  !=================================================================================================
  subroutine nsi_rk_momentum_update_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_update_t), intent(inout) :: approx
    class(physical_problem_t)  , target, intent(in)    :: physics
    class(discrete_problem_t)  , target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

   ! Allocate working variables
    call memalloc(physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = ivar
    end do

  end subroutine nsi_rk_momentum_update_create

  !=================================================================================================
  subroutine nsi_rk_momentum_update_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_update_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_rk_momentum_update_free

  !=================================================================================================
  subroutine nsi_rk_projection_update_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_projection_update_t), intent(inout) :: approx
    class(physical_problem_t)  , target, intent(in)    :: physics
    class(discrete_problem_t)  , target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(physics%ndime,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,physics%ndime
       approx%working_vars(ivar) = physics%ndime + 1 + ivar
    end do

  end subroutine nsi_rk_projection_update_create

  !=================================================================================================
  subroutine nsi_rk_projection_update_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_projection_update_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_rk_projection_update_free

  !=================================================================================================
  subroutine nsi_dissipation_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_dissipation_t), intent(inout) :: approx
    class(physical_problem_t)  , target, intent(in)    :: physics
    class(discrete_problem_t)  , target, intent(in)    :: discret
    ! Locals
    integer(ip) :: ivar

    select type (physics)
    type is(nsi_problem_t)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_oss_discrete_t)
       approx%discret => discret
       class default
       check(.false.)
    end select

    ! Domain dimension
    approx%domain_dimension = 3

    ! Allocate working variables
    call memalloc(discret%nvars,approx%working_vars,__FILE__,__LINE__)
    do ivar=1,discret%nvars
       approx%working_vars(ivar) = ivar
    end do

  end subroutine nsi_dissipation_create

  !=================================================================================================
  subroutine nsi_dissipation_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_dissipation_t), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

    ! Deallocate working variables
    call memfree(approx%working_vars,__FILE__,__LINE__)

  end subroutine nsi_dissipation_free

  !=================================================================================================
  subroutine nsi_matvec(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental matrix-vector integration selection.                !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_matvec_t), intent(inout) :: approx
    type(finite_element_t)        , intent(inout) :: finite_element
    ! Locals
    real(rp), allocatable :: elmat_vu(:,:,:,:)
    real(rp), allocatable :: elmat_vu_diag(:,:)
    real(rp), allocatable :: elmat_vp(:,:,:,:)
    real(rp), allocatable :: elmat_qu(:,:,:,:)
    real(rp), allocatable :: elmat_vx_diag(:,:)
    real(rp), allocatable :: elmat_wu_diag(:,:)
    real(rp), allocatable :: elmat_wx_diag(:,:)
    real(rp), allocatable :: elvec_u(:,:) 
    integer(ip)           :: igaus,idime,inode,jdime,jnode,idof,jdof
    integer(ip)           :: ngaus,ndime,nnodu,nnodp
    real(rp)              :: ctime,dtinv,dvolu,diffu,react
    real(rp)              :: work(4)
    real(rp)              :: agran(finite_element%integ(1)%p%uint_phy%nnode)
    real(rp)              :: testf(finite_element%integ(1)%p%uint_phy%nnode,finite_element%integ(1)%p%quad%ngaus)
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    type(vector_t)        :: gpvel, gpveln, force

    ! Checks
    !check(finite_element%reference_element_vars(1)%p%order > finite_element%reference_element_vars(approx%physics%ndime+1)%p%order)
    check(finite_element%integ(1)%p%quad%ngaus == finite_element%integ(approx%physics%ndime+1)%p%quad%ngaus)
    do idime=2,approx%physics%ndime
       check(finite_element%integ(1)%p%uint_phy%nnode == finite_element%integ(idime)%p%uint_phy%nnode)
    end do

    ! Initialize matrix and vector
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    react = approx%physics%react
    dtinv = approx%tinteg%dtinv

    ! Allocate auxiliar matrices and vectors
    call memalloc(ndime,ndime,nnodu,nnodu,elmat_vu,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)
    call memalloc(ndime,1,nnodu,nnodp,elmat_vp,__FILE__,__LINE__)
    call memalloc(1,ndime,nnodp,nnodu,elmat_qu,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_vx_diag,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_wu_diag,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_wx_diag,__FILE__,__LINE__)
    call memalloc(ndime,nnodu,elvec_u,__FILE__,__LINE__)

    ! Initialize to zero
    elmat_vu      = 0.0_rp
    elmat_vu_diag = 0.0_rp
    elmat_vp      = 0.0_rp
    elmat_qu      = 0.0_rp
    elmat_vx_diag = 0.0_rp
    elmat_wu_diag = 0.0_rp
    elmat_wx_diag = 0.0_rp
    elvec_u       = 0.0_rp

    ! Interpolation operations for velocity
    call create_vector (approx%physics, 1, finite_element%integ, gpvel)
    call create_vector (approx%physics, 1, finite_element%integ, gpveln)
    call interpolation (finite_element%unkno, 1, prev_iter, finite_element%integ, gpvel)
    if(dtinv == 0.0_rp) then
       call interpolation (finite_element%unkno, 1, prev_step, finite_element%integ, gpveln)
    else
       gpveln%a = 0.0_rp
    end if

    ! Set real time
    if(approx%discret%kfl_thet==0) then        ! BE
       ctime = approx%tinteg%ctime
    elseif(approx%discret%kfl_thet==1) then    ! CN
       ctime = approx%tinteg%ctime - 1.0_rp/dtinv
    end if

    ! Set force term
    call create_vector(approx%physics,1,finite_element%integ,force)
    force%a=0.0_rp
    ! Impose analytical solution
    if(finite_element%p_analytical_code%a(1,1)>0.and.finite_element%p_analytical_code%a(ndime+1,1)>0) then 
       call nsi_analytical_force(approx%physics,finite_element,ctime,gpvel,force)
    end if

    ! Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)
    
    ! Initializations
    work  = 0.0_rp
    agran = 0.0_rp 
    testf = 0.0_rp

    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + &
                     &         gpvel%a(idime,igaus)*finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
             testf(inode,igaus) = tau(1,igaus)*agran(inode)
          end do
       end if

       ! Add external force term
       do idime=1,approx%physics%ndime    
          force%a(idime,igaus) = force%a(idime,igaus) + approx%physics%gravi(idime)
       end do

       ! Computation of elemental terms
       ! ------------------------------
       ! Block U-V
       ! mu * ( grad u, grad v )
       call elmvis_gal(dvolu,diffu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elmat_vu_diag, &
            &          work)
       ! Add cross terms for symmetric grad
       if(approx%physics%kfl_symg==1) then
          call elmvis_gal_sym(dvolu,diffu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
               &              elmat_vu,work)
       end if
       if(approx%physics%kfl_skew==0) then
          ! (v, a·grad u) + s*(v,u) + (v, u/dt)
          call elmbuv_gal(dvolu,react,dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &          elmat_vu_diag,work)
       elseif(approx%physics%kfl_skew==1) then
          ! 1/2(v, a·grad u) - 1/2(u,a·grad v) + s*(v,u) + (v, u/dt)
          call elmbuv_gal_skew1(dvolu,react,dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &                elmat_vu_diag,work)
       end if
       ! tau*(a·grad u, a·grad v)
       call elmbuv_oss(dvolu,testf(:,igaus),agran,nnodu,elmat_vu_diag,work)
       ! tauc*(div v, div u)
       if(approx%discret%ktauc>0.0_rp) then
          call elmdiv_stab(tau(2,igaus),dvolu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
               &           elmat_vu,work)
       end if

       ! Block X-V
       ! -tau*(proj(a·grad u), a·grad v)
       if(approx%discret%kfl_proj==1) then
          work(2) = -tau(1,igaus)*dvolu
       else
          work(2) = -dvolu
       end if
       call elmbvu_gal(work(2),finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu,elmat_vx_diag, &
            &          work)

       ! Block P-V
       ! - ( div v, p )
       call elmbpv_gal_div_iss(dvolu,finite_element%integ(ndime+1)%p%uint_phy%shape(:,igaus),         &
            &                  finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_vp,work)

       ! Block U-W
       ! -tau*(a·grad u, v)
       call elmbuv_gal(-dvolu*tau(1,igaus),0.0_rp,0.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran, &
            &          nnodu,elmat_wu_diag,work)

       ! Block X-W
       ! tau*(proj(a·grad u),v)
       if(approx%discret%kfl_proj==1) then
          call elmmss_gal(dvolu,tau(1,igaus),finite_element%integ(1)%p%uint_phy%shape(:,igaus),nnodu, &
               &          elmat_wx_diag,work)
       else
          call elmmss_gal(dvolu,1.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus),nnodu,elmat_wx_diag, &
               &          work)
       end if
       
       ! Block U-Q
       ! ( div u, q )
       call elmbuq_gal_div_iss(dvolu,finite_element%integ(ndime+1)%p%uint_phy%shape(:,igaus),               &
            &                  finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_qu,work)

       ! RHS: Block U
       ! ( v, f ) + ( v, u_n/dt )
       call elmrhu_gal(dvolu,dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),gpveln%a(:,igaus), &
            &          force%a(:,igaus),nnodu,ndime,elvec_u,work)

    end do

    call memfree(gpvel%a,__FILE__,__LINE__)
    call memfree(gpveln%a,__FILE__,__LINE__)
    call memfree(force%a,__FILE__,__LINE__)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             do jdime=1,ndime
                ! Block V-U
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(jdime)+jnode-1
                finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vu(idime,jdime,inode,jnode)
             end do    
             ! Block V-U (diag)
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) +  elmat_vu_diag(inode,jnode)
             ! Block V-X
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(ndime+1+idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vx_diag(inode,jnode)
             ! Block W-U
             idof = finite_element%start%a(ndime+1+idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_wu_diag(inode,jnode)
             ! Block W-X
             idof = finite_element%start%a(ndime+1+idime)+inode-1
             jdof = finite_element%start%a(ndime+1+idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_wx_diag(inode,jnode)
          end do
          do jnode=1,nnodp
             ! Block V-P
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(ndime+1)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vp(idime,1,inode,jnode)
             ! Block Q-U
             idof = finite_element%start%a(ndime+1)+jnode-1
             jdof = finite_element%start%a(idime)+inode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) +  elmat_qu(1,idime,jnode,inode)
          end do
          ! Block U
          idof = finite_element%start%a(idime)+inode-1
          finite_element%p_vec%a(idof) = finite_element%p_vec%a(idof) + elvec_u(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu,__FILE__,__LINE__)
    call memfree(elmat_vu_diag,__FILE__,__LINE__)
    call memfree(elmat_vp,__FILE__,__LINE__)
    call memfree(elmat_qu,__FILE__,__LINE__)
    call memfree(elmat_vx_diag,__FILE__,__LINE__)
    call memfree(elmat_wu_diag,__FILE__,__LINE__)
    call memfree(elmat_wx_diag,__FILE__,__LINE__)
    call memfree(elvec_u,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 
    
  end subroutine nsi_matvec

  !=================================================================================================
  subroutine nsi_massu(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental mass matrix integration for pressure dofs.          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massu_t), intent(inout) :: approx
    type(finite_element_t)       , intent(inout) :: finite_element
    ! Locals
    integer(ip) :: igaus,inode,jnode,idof,jdof,ndime,nnodu,ngaus,idime
    real(rp)    :: dvolu
    real(rp), allocatable :: elmat_vu_diag(:,:)
    
    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp
    
    ! Allocate auxiliar matrices and vectors
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)

    ! Initialize to zero
    elmat_vu_diag = 0.0_rp

    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       do inode=1,nnodu
          do jnode=1,nnodu
             elmat_vu_diag(inode,jnode) = elmat_vu_diag(inode,jnode) + dvolu * &
                  &    finite_element%integ(1)%p%uint_phy%shape(inode,igaus) * &
                  &    finite_element%integ(1)%p%uint_phy%shape(jnode,igaus)
          end do
       end do

    end do

    ! Assembly to elemental p_mat
    if (approx%discret%kfl_lump==0) then
       do inode=1,nnodu
          do jnode=1,nnodu
             do idime=1,ndime
                ! Block V-U (diag)
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(idime)+jnode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + &
                     &                               elmat_vu_diag(inode,jnode)
             end do
          end do
       end do
    else
       do inode=1,nnodu
          do jnode=1,nnodu
             do idime=1,ndime
                ! Block V-U (diag)
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(idime)+inode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + &
                     &                               elmat_vu_diag(inode,jnode)
             end do
          end do
       end do
    end if

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu_diag,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 
    
  end subroutine nsi_massu

  !=================================================================================================
  subroutine nsi_massp(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental mass matrix integration for pressure dofs.          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massp_t), intent(inout) :: approx
    type(finite_element_t)       , intent(inout) :: finite_element
    ! Locals
    integer(ip) :: igaus,inode,jnode,idof,jdof,ndime,nnodp,ngaus
    real(rp)    :: dvolu
    
    ! Unpack variables
    ndime = approx%physics%ndime
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp
    
    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)
       if (approx%discret%kfl_lump==0) then
          do inode=1,nnodp
             do jnode=1,nnodp
                ! Block Q-P
                idof = finite_element%start%a(ndime+1)+inode-1
                jdof = finite_element%start%a(ndime+1)+jnode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + dvolu * &
                     &  finite_element%integ(ndime+1)%p%uint_phy%shape(inode,igaus) * &
                     &  finite_element%integ(ndime+1)%p%uint_phy%shape(jnode,igaus)
             end do
          end do
       else
          do inode=1,nnodp
             do jnode=1,nnodp
                ! Block Q-P
                idof = finite_element%start%a(ndime+1)+inode-1
                jdof = finite_element%start%a(ndime+1)+inode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + dvolu * &
                     &  finite_element%integ(ndime+1)%p%uint_phy%shape(inode,igaus) * &
                     &  finite_element%integ(ndime+1)%p%uint_phy%shape(jnode,igaus)
             end do
          end do
       end if
    end do
    
  end subroutine nsi_massp

  !=================================================================================================
  subroutine nsi_massx(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental mass matrix integration for pressure dofs.          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_massx_t), intent(inout) :: approx
    type(finite_element_t)       , intent(inout) :: finite_element
    ! Locals
    integer(ip)           :: igaus,inode,jnode,idof,jdof,ndime,nnodu,ngaus,idime
    real(rp)              :: dvolu,work
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    type(vector_t)        :: gpvel
    real(rp), allocatable :: elmat_wx_diag(:,:)
    
    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp
    
    ! Allocate auxiliar matrices and vectors
    call memalloc(nnodu,nnodu,elmat_wx_diag,__FILE__,__LINE__)

    ! Interpolation operations for velocity
    call create_vector (approx%physics, 1, finite_element%integ, gpvel)
    call interpolation (finite_element%unkno, 1, prev_iter, finite_element%integ, gpvel)

    ! Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)

    ! Initialize to zero
    elmat_wx_diag = 0.0_rp

    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       if(approx%discret%kfl_proj==1) then
          work = dvolu*tau(1,igaus)
       else
          work = dvolu
       end if

       do inode=1,nnodu
          do jnode=1,nnodu
             elmat_wx_diag(inode,jnode) = elmat_wx_diag(inode,jnode) + work  * &
                  &    finite_element%integ(1)%p%uint_phy%shape(inode,igaus) * &
                  &    finite_element%integ(1)%p%uint_phy%shape(jnode,igaus)
          end do
       end do

    end do

    ! Deallocate
    call memfree(gpvel%a,__FILE__,__LINE__)

    ! Assembly to elemental p_mat
    if (approx%discret%kfl_lump==0) then
       do inode=1,nnodu
          do jnode=1,nnodu
             do idime=ndime+2,2*ndime+1
                ! Block W-X (diag)
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(idime)+jnode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + &
                     &                               elmat_wx_diag(inode,jnode)
             end do
          end do
       end do
    else
       do inode=1,nnodu
          do jnode=1,nnodu
             do idime=ndime+2,2*ndime+1
                ! Block W-X (diag)
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(idime)+inode-1
                finite_element%p_mat%a(idof,jdof) =  finite_element%p_mat%a(idof,jdof) + &
                     &                               elmat_wx_diag(inode,jnode)
             end do
          end do
       end do
    end if

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_wx_diag,__FILE__,__LINE__)
    
  end subroutine nsi_massx

  !=================================================================================================
  subroutine nsi_lapla_p(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental laplacian matrix integration for pressure dofs.     !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_lapla_p_t), intent(inout) :: approx
    type(finite_element_t)         , intent(inout) :: finite_element
    ! Locals
    real(rp), allocatable :: elmat_qp(:,:)
    integer(ip)           :: igaus,inode,jnode,idof,jdof,ndime,nnodp,ngaus
    real(rp)              :: dvolu
    real(rp)              :: work(4)
    
    ! Unpack variables
    ndime = approx%physics%ndime
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Allocate auxiliar matrix
    call memalloc(nnodp,nnodp,elmat_qp,__FILE__,__LINE__)
    elmat_qp = 0.0_rp
    
    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)
       
       ! (grad p, grad q)
       call elmvis_gal(dvolu,1.0_rp,finite_element%integ(ndime+1)%p%uint_phy%deriv(:,:,igaus),ndime, &
            &          nnodp,elmat_qp,work)
    end do

    ! Assembly to elemental p_mat 
    do inode=1,nnodp
       do jnode=1,nnodp
          idof = finite_element%start%a(ndime+1)+inode-1
          jdof = finite_element%start%a(ndime+1)+jnode-1
          finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_qp(inode,jnode)
       end do
    end do

    call memfree(elmat_qp,__FILE__,__LINE__)

  end subroutine nsi_lapla_p

  !=================================================================================================
  subroutine nsi_rk_momentum(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental momentum matrix integration for velocity dofs.      !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_t), intent(inout) :: approx
    type(finite_element_t)             , intent(inout) :: finite_element
    ! Locals
    integer(ip)           :: ndime,nnodu,nnodp,ngaus
    integer(ip)           :: igaus,inode,jnode,idof,jdof,istage,jstge,idime,jdime
    real(rp)              :: dvolu,dtinv,diffu,ctime,prevtime
    real(rp)              :: work(4),beta
    real(rp)              :: agran(finite_element%integ(1)%p%uint_phy%nnode)
    real(rp)              :: testf(finite_element%integ(1)%p%uint_phy%nnode,finite_element%integ(1)%p%quad%ngaus)
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    real(rp), allocatable :: elmat_vu(:,:,:,:)
    real(rp), allocatable :: elmat_vu_diag(:,:)
    real(rp), allocatable :: elmat_vx_diag(:,:)
    real(rp), allocatable :: elmat_wu_diag(:,:)
    real(rp), allocatable :: elmat_wx_diag(:,:)
    type(vector_t)        :: gpvel
    type(rungekutta_integrator_t), pointer :: rkinteg

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    rkinteg => approx%rkinteg
    dtinv = rkinteg%dtinv
    istage = rkinteg%istage
    prevtime = rkinteg%ctime - 1.0_rp/dtinv

    ! Asserts
    !check(istage>1) ! First stage skiped (a_11=0)

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    call create_vector(approx%physics,1,finite_element%integ,gpvel)
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpvel) ! U^k,i

    ! Allocate auxiliar matrices and vectors
    call memalloc(ndime,ndime,nnodu,nnodu,elmat_vu,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_vx_diag,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_wu_diag,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_wx_diag,__FILE__,__LINE__)

    ! Allocate & compute Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)

    ! Initialize to zero
    elmat_vu      = 0.0_rp
    elmat_vu_diag = 0.0_rp
    elmat_vx_diag = 0.0_rp
    elmat_wu_diag = 0.0_rp
    elmat_wx_diag = 0.0_rp
    work          = 0.0_rp
    agran         = 0.0_rp
    testf         = 0.0_rp

    ! Loop on Gauss points
    gauss: do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + gpvel%a(idime,igaus) * &
                     &         finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
             testf(inode,igaus) = tau(1,igaus)*agran(inode)
          end do
       end if

       ! Construct LHS
       ! =============

       ! Diffusion
       if(rkinteg%rk_terms(1)%hsite == implicit.and.(approx%integration_stage==rkinteg%rk_terms(1)%ltype)) then
          call elmvis_gal(dvolu,diffu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime, &
               &          nnodu,elmat_vu_diag,work)
          ! Add cross terms for symmetric grad
          if(approx%physics%kfl_symg==1) then
             call elmvis_gal_sym(dvolu,diffu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus), &
                  &              ndime,nnodu,elmat_vu,work)
          end if
       end if

       ! Convection
       if(approx%physics%kfl_conv.ne.0) then
          if(rkinteg%rk_terms(2)%hsite == implicit.and.(approx%integration_stage==rkinteg%rk_terms(2)%ltype)) then
             if(approx%physics%kfl_skew==0) then
                ! (v, a·grad u)
                call elmbuv_gal(dvolu,0.0_rp,0.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                     &          agran,nnodu,elmat_vu_diag,work)
             elseif(approx%physics%kfl_skew==1) then
                ! 1/2(v, a·grad u) - 1/2(u,a·grad v)
                call elmbuv_gal_skew1(dvolu,0.0_rp,0.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                     &                agran,nnodu,elmat_vu_diag,work)
             end if
          end if
       end if

       ! OSS_vu
       if(rkinteg%rk_terms(4)%hsite == implicit.and.(approx%integration_stage==rkinteg%rk_terms(4)%ltype)) then
          ! tau*(a·grad u, a·grad v)
          call elmbuv_oss(dvolu,testf(:,igaus),agran,nnodu,elmat_vu_diag,work)
          ! tauc*(div v, div u)
          if(approx%discret%ktauc>0.0_rp) then
             call elmdiv_stab(tau(2,igaus),dvolu,finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus), &
                  &           ndime,nnodu,elmat_vu,work)
          end if
       end if

       ! OSS_vx
       if(rkinteg%rk_terms(5)%hsite == implicit.and.(approx%integration_stage==rkinteg%rk_terms(5)%ltype)) then
          ! -tau*(proj(a·grad u), a·grad v)
          if(approx%discret%kfl_proj==1) then
             beta = -tau(1,igaus)*dvolu
          else
             beta = -dvolu
          end if
          call elmbvu_gal(beta,finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &          elmat_vx_diag,work)
       end if

       if(approx%integration_stage==update_nonlinear) then
          ! OSS_wu
          ! -tau*(a·grad u, w)
          beta = -tau(1,igaus)*dvolu
          call elmbuv_gal(beta,0.0_rp,0.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus),agran, &
               &          nnodu,elmat_wu_diag,work)

          ! OSS_wx
          ! tau*(proj(a·grad u), w) / (proj(tau*a·grad u), w)
          if(approx%discret%kfl_proj==1) then
             call elmmss_gal(dvolu,tau(1,igaus),finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                  &          nnodu,elmat_wx_diag,work)
          else
             call elmmss_gal(dvolu,1.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                  &          nnodu,elmat_wx_diag,work)
          end if
       end if

    end do gauss

    ! Deallocate
    call memfree(gpvel%a,__FILE__,__LINE__)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             do jdime=1,ndime
                ! Block V-U
                idof = finite_element%start%a(idime)+inode-1
                jdof = finite_element%start%a(jdime)+jnode-1
                finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vu(idime,jdime,inode,jnode)
             end do    
             ! Block V-U (diag)
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vu_diag(inode,jnode)
             ! Block V-X
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(ndime+1+idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vx_diag(inode,jnode)
             ! Block W-U
             idof = finite_element%start%a(ndime+1+idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_wu_diag(inode,jnode)
             ! Block W-X
             idof = finite_element%start%a(ndime+1+idime)+inode-1
             jdof = finite_element%start%a(ndime+1+idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_wx_diag(inode,jnode)
          end do
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu,__FILE__,__LINE__)
    call memfree(elmat_vu_diag,__FILE__,__LINE__)
    call memfree(elmat_vx_diag,__FILE__,__LINE__)
    call memfree(elmat_wu_diag,__FILE__,__LINE__)
    call memfree(elmat_wx_diag,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine nsi_rk_momentum

  !=================================================================================================
  subroutine nsi_rk_momentum_rhs(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental momentum RHS integration for velocity dofs.         !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_rhs_t), intent(inout) :: approx
    type(finite_element_t)                 , intent(inout) :: finite_element
    ! Locals
    integer(ip)                 :: ndime,nnodu,nnodp,ngaus
    integer(ip)                 :: igaus,inode,jnode,idof,jdof,istage,jstge,idime,jdime
    real(rp)                    :: dvolu,dtinv,diffu,alpha,ctime,prevtime
    real(rp)                    :: work(4),beta
    real(rp)                    :: agran(finite_element%integ(1)%p%uint_phy%nnode)
    real(rp)      , allocatable :: testf(:,:,:)
    real(rp)      , allocatable :: tau(:,:,:)
    real(rp)      , allocatable :: elvec_u(:,:) 
    type(vector_t), allocatable :: gpvel(:),gposs(:),force(:),grpre(:)
    type(tensor_t), allocatable :: grvel(:)
    type(rungekutta_integrator_t), pointer :: rkinteg

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    rkinteg => approx%rkinteg
    dtinv = rkinteg%dtinv
    istage = rkinteg%istage
    prevtime = rkinteg%ctime - 1.0_rp/dtinv

    ! Asserts
    !check(istage>1) ! First stage skiped (a_11=0)

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    allocate(gpvel(istage+1))
    allocate(gposs(istage-1))
    allocate(grvel(istage-1))
    allocate(grpre(istage-1))
    do jstge=1,istage-1
       call create_vector(approx%physics,1,finite_element%integ,gpvel(jstge))
       call create_vector(approx%physics,1,finite_element%integ,gposs(jstge))
       call create_tensor(approx%physics,1,ndime,finite_element%integ,grvel(jstge))
       call create_vector(approx%physics,1,finite_element%integ,grpre(jstge))
    end do
    call create_vector (approx%physics,1,finite_element%integ,gpvel(istage))
    call create_vector (approx%physics,1,finite_element%integ,gpvel(istage+1))
    call interpolation(finite_element%unkno,1,prev_step,finite_element%integ,gpvel(1))                    ! U_n
    do jstge=2,istage
       call interpolation(finite_element%unkno,1,jstge+2,finite_element%integ,gpvel(jstge))               ! U_j
       call interpolation(finite_element%unkno,ndime+2,jstge+2,finite_element%integ,gposs(jstge-1))       ! X_j
       call interpolation(finite_element%unkno,1,jstge+2,finite_element%integ,grvel(jstge-1))             ! GradU_j
       call interpolation(finite_element%unkno,ndime+1,ndime,jstge+2,finite_element%integ,grpre(jstge-1)) ! GradP_j
    end do
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpvel(istage+1))             ! U^k_i

    ! Allocate auxiliar matrices and vectors
    call memalloc(ndime,nnodu,elvec_u,__FILE__,__LINE__)

    ! Allocate & compute Stabilization parameters
    call memalloc(nnodu,ngaus,istage,testf,__FILE__,__LINE__)
    call memalloc(2,ngaus,istage,tau,__FILE__,__LINE__)
    do jstge=2,istage+1
       call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel(jstge)%a,tau(:,:,jstge-1))
    end do

    ! Allocate & compute force term
    allocate(force(istage))
    do jstge=1,istage
       ctime = prevtime + 1.0_rp/dtinv*rkinteg%rk_table(1)%p%c(jstge)
       call create_vector(approx%physics,1,finite_element%integ,force(jstge))
       force(jstge)%a=0.0_rp
       ! Impose analytical solution
       if(finite_element%p_analytical_code%a(1,1)>0.and.finite_element%p_analytical_code%a(ndime+1,1)>0) then 
          call nsi_analytical_force(approx%physics,finite_element,ctime,gpvel(jstge+1),force(jstge))
       end if
       ! Add external force term
       do igaus=1,ngaus
          force(jstge)%a(:,igaus) = force(jstge)%a(:,igaus) + approx%physics%gravi(:)
       end do
    end do

    ! Initialize to zero
    elvec_u       = 0.0_rp
    work          = 0.0_rp
    agran         = 0.0_rp
    testf         = 0.0_rp

    ! Loop on Gauss points
    gauss: do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + gpvel(istage+1)%a(idime,igaus) * &
                     &         finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
             testf(inode,igaus,istage) = tau(1,igaus,istage)*agran(inode)
          end do
       end if

       ! Construct RHS
       ! =============
       ! Mass & force(istage)
       ! a_ii*( v, f ) + ( v, u_n/dt )
       force(istage)%a(:,igaus) = force(istage)%a(:,igaus)*rkinteg%rk_table(6)%p%A(istage,istage)
       call elmrhu_gal(dvolu,dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),gpvel(1)%a(:,igaus), &
            &           force(istage)%a(:,igaus),nnodu,ndime,elvec_u,work)

       ! Loop over previous stages
       stages: do jstge=1,istage-1

          ! Auxiliar variables
          if(approx%physics%kfl_conv.ne.0) then
             do inode = 1,nnodu
                agran(inode) = 0.0_rp
                do idime = 1,ndime
                   agran(inode) = agran(inode) + gpvel(jstge+1)%a(idime,igaus) * &
                        &         finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
                end do
                testf(inode,igaus,jstge) = tau(1,igaus,jstge)*agran(inode)
             end do
          end if

          ! Diffusion
          alpha = rkinteg%rk_table(1)%p%A(istage,jstge)
          beta = -dvolu*alpha
          ! nu( grad u, grad v )
          call elmrhu_gradgrad(beta,diffu,grvel(jstge)%a(:,:,igaus), &
               &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
          ! Add cross terms for symmetric grad
          if(approx%physics%kfl_symg==1) then
             call elmrhu_gradgrad_sym(beta,diffu,grvel(jstge)%a(:,:,igaus), &
                  &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
          end if

          ! Convection
          if(approx%physics%kfl_conv.ne.0) then
             alpha = rkinteg%rk_table(2)%p%A(istage,jstge)
             beta = -dvolu*alpha
             ! (u · grad u, v) 
             if(approx%physics%kfl_skew==0) then
                call oss_convu_arg_chk(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus), &
                     &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_u,work)
             elseif(approx%physics%kfl_skew==1) then
                call oss_convu_arg_chk_skew(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus),agran, &
                     &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                     &                 finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
             end if
          end if

          ! Pressure gradient
          alpha = rkinteg%rk_table(3)%p%A(istage,jstge)
          beta = -dvolu*alpha
          ! (grad p, v)
          call oss_gradp_arg_chk(beta,grpre(jstge)%a(:,igaus), &
               &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_u,work)

          ! OSS_vu
          alpha = rkinteg%rk_table(4)%p%A(istage,jstge)
          beta = -dvolu*alpha*tau(1,igaus,jstge)
          ! tau * ( u · grad u, u· grad v)
          call oss_convu_arg_chk(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus),agran, &
               &                 ndime,nnodu,elvec_u,work)
          ! tauc*(div v, div u)
          if(approx%discret%ktauc>0.0_rp) then
             beta = -dvolu*alpha*tau(2,igaus,jstge)
             call elmrhu_divudivv(beta,grvel(jstge)%a(:,:,igaus), &
                  &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
                  &               elvec_u,work)
          end if

          ! OSS_vx
          alpha = rkinteg%rk_table(5)%p%A(istage,jstge)
          beta = dvolu*alpha
          ! - tau * ( proj(u·grad u), u· grad v)
          if(approx%discret%kfl_proj==1) then
             call elmrhs_fce(beta,testf(:,igaus,jstge),gposs(jstge)%a(:,igaus),nnodu,ndime,elvec_u, &
                  &          work)
          else
             call elmrhs_fce(beta,agran,gposs(jstge)%a(:,igaus),nnodu,ndime,elvec_u, &
                  &          work)
          end if

          ! Force
          alpha = rkinteg%rk_table(6)%p%A(istage,jstge)
          beta = dvolu*alpha
          call elmrhs_fce(beta,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
               &          force(jstge)%a(:,igaus),nnodu,ndime,elvec_u,work)

       end do stages

    end do gauss

    ! Deallocate
    do jstge=1,istage-1
       call memfree(gpvel(jstge)%a,__FILE__,__LINE__)
       call memfree(gposs(jstge)%a,__FILE__,__LINE__)
       call memfree(grvel(jstge)%a,__FILE__,__LINE__)
       call memfree(grpre(jstge)%a,__FILE__,__LINE__)
       call memfree(force(jstge)%a,__FILE__,__LINE__)
    end do
    call memfree(gpvel(istage)%a,__FILE__,__LINE__)
    call memfree(gpvel(istage+1)%a,__FILE__,__LINE__)
    call memfree(force(istage)%a,__FILE__,__LINE__)
    deallocate(gpvel)
    deallocate(gposs)
    deallocate(grvel)
    deallocate(grpre)
    deallocate(force)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          ! Block U
          idof = finite_element%start%a(idime)+inode-1
          finite_element%p_vec%a(idof) = finite_element%p_vec%a(idof) + elvec_u(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elvec_u,__FILE__,__LINE__)
    call memfree(testf,__FILE__,__LINE__)
    call memfree(tau,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine nsi_rk_momentum_rhs

  !=================================================================================================
  subroutine nsi_rk_pressure(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental darcy matrix integration for pressure dofs.         !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_pressure_t), intent(inout) :: approx
    type(finite_element_t)             , intent(inout) :: finite_element
    ! Locals
    integer(ip)           :: ndime,nnodu,nnodp,ngaus
    integer(ip)           :: igaus,inode,jnode,idof,jdof,idime
    real(rp)              :: dvolu,diffu,dtinv,ctime 
    real(rp)              :: work(4),beta
    real(rp)              :: agran(finite_element%integ(1)%p%uint_phy%nnode)
    real(rp)              :: testf(finite_element%integ(1)%p%uint_phy%nnode,finite_element%integ(1)%p%quad%ngaus)
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    real(rp), allocatable :: elmat_vu_diag(:,:)
    real(rp), allocatable :: elmat_vp(:,:,:,:)
    real(rp), allocatable :: elmat_qu(:,:,:,:)
    real(rp), allocatable :: elvec_u(:,:) 
    type(vector_t)        :: gpvel,gposs,force
    type(tensor_t)        :: grvel
    
    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    dtinv = approx%tinteg%dtinv
    ctime = approx%tinteg%ctime

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    call create_vector(approx%physics,1,finite_element%integ,gpvel)
    call create_vector(approx%physics,1,finite_element%integ,gposs)
    call create_tensor(approx%physics,1,ndime,finite_element%integ,grvel)
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpvel)       ! U_j
    call interpolation(finite_element%unkno,ndime+2,prev_iter,finite_element%integ,gposs) ! X_j
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,grvel)       ! GradU_j

    ! Allocate auxiliar matrices and vectors
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)
    call memalloc(ndime,1,nnodu,nnodp,elmat_vp,__FILE__,__LINE__)
    call memalloc(1,ndime,nnodp,nnodu,elmat_qu,__FILE__,__LINE__)
    call memalloc(ndime,nnodu,elvec_u,__FILE__,__LINE__)

    ! Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)

    ! Set force term
    call create_vector(approx%physics,1,finite_element%integ,force)
    force%a=0.0_rp
    ! Impose analytical solution
    if(finite_element%p_analytical_code%a(1,1)>0.and.finite_element%p_analytical_code%a(ndime+1,1)>0) then 
       call nsi_analytical_force(approx%physics,finite_element,ctime,gpvel,force)
    end if
    ! Add external force term
    do igaus=1,ngaus
       force%a(:,igaus) = force%a(:,igaus) + approx%physics%gravi(:)
    end do  

    ! Initialize to zero
    elmat_vu_diag = 0.0_rp
    elmat_vp      = 0.0_rp
    elmat_qu      = 0.0_rp
    elvec_u       = 0.0_rp
    work          = 0.0_rp
    agran         = 0.0_rp
    testf         = 0.0_rp

    ! Loop on Gauss points
    do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + gpvel%a(idime,igaus) * &
                     &         finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
             testf(inode,igaus) = tau(1,igaus)*agran(inode)
          end do
       end if

       ! Construct LHS
       ! =============
       ! Mass matrix (M)
       call elmmss_gal(dvolu,1.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus),nnodu, &
            &          elmat_vu_diag,work)

       ! Pressure Gradient
       ! - ( div v, p )
       call elmbpv_gal_div_iss(dvolu,finite_element%integ(ndime+1)%p%uint_phy%shape(:,igaus),         &
            &                  finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_vp,work)

       ! Divergence
       ! ( div u, q )
       call elmbuq_gal_div_iss(dvolu,finite_element%integ(ndime+1)%p%uint_phy%shape(:,igaus),               &
            &                  finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_qu,work)
       
       ! Construct RHS
       ! =============
       ! Diffusion
       beta = -dvolu
       ! nu( grad u, grad v )
       call elmrhu_gradgrad(beta,diffu,grvel%a(:,:,igaus), &
            &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
       ! Add cross terms for symmetric grad
       if(approx%physics%kfl_symg==1) then
          call elmrhu_gradgrad_sym(beta,diffu,grvel%a(:,:,igaus), &
               &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
       end if

       ! Convection
       beta = -dvolu
       ! (u · grad u, v) 
       if(approx%physics%kfl_conv.ne.0) then
          if(approx%physics%kfl_skew==0) then
             call oss_convu_arg_chk(beta,gpvel%a(:,igaus),grvel%a(:,:,igaus), &
                  &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_u,work)
          elseif(approx%physics%kfl_skew==1) then
             call oss_convu_arg_chk_skew(beta,gpvel%a(:,igaus),grvel%a(:,:,igaus),agran, &
                  &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                  &                 finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
          end if
       end if
       
       ! OSS_vu
       beta = -dvolu*tau(1,igaus)
       ! tau * ( u · grad u, u· grad v)
       call oss_convu_arg_chk(beta,gpvel%a(:,igaus),grvel%a(:,:,igaus),agran,ndime,nnodu,elvec_u, &
            &                 work)
       ! tauc*(div v, div u)
       if(approx%discret%ktauc>0.0_rp) then
          beta = -dvolu*tau(2,igaus)
          call elmrhu_divudivv(beta,grvel%a(:,:,igaus), &
               &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
               &               elvec_u,work)
       end if

       ! OSS_vx
       beta = dvolu
       ! - tau * ( proj(u·grad u), u· grad v)
       if(approx%discret%kfl_proj==1) then
          call elmrhs_fce(beta,testf(:,igaus),gposs%a(:,igaus),nnodu,ndime,elvec_u,work)
       else
          call elmrhs_fce(beta,agran,gposs%a(:,igaus),nnodu,ndime,elvec_u,work)
       end if

       ! Force
       beta = dvolu
       call elmrhs_fce(beta,finite_element%integ(1)%p%uint_phy%shape(:,igaus),force%a(:,igaus), &
            &          nnodu,ndime,elvec_u,work)

    end do

    ! Deallocate
    call memfree(gpvel%a,__FILE__,__LINE__)
    call memfree(gposs%a,__FILE__,__LINE__)
    call memfree(grvel%a,__FILE__,__LINE__)
    call memfree(force%a,__FILE__,__LINE__)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             ! Block V-U (diag)
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) +  elmat_vu_diag(inode,jnode)
          end do
          do jnode=1,nnodp
             ! Block V-P
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(ndime+1)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_vp(idime,1,inode,jnode)
             ! Block Q-U
             idof = finite_element%start%a(ndime+1)+jnode-1
             jdof = finite_element%start%a(idime)+inode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) +  elmat_qu(1,idime,jnode,inode)
          end do
          ! Block U
          idof = finite_element%start%a(idime)+inode-1
          finite_element%p_vec%a(idof) = finite_element%p_vec%a(idof) + elvec_u(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu_diag,__FILE__,__LINE__)
    call memfree(elmat_vp,__FILE__,__LINE__)
    call memfree(elmat_qu,__FILE__,__LINE__)
    call memfree(elvec_u,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 
    
  end subroutine nsi_rk_pressure

  !=================================================================================================
  subroutine nsi_rk_momentum_update(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental momentum matrix integration for velocity dofs.      !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_momentum_update_t), intent(inout) :: approx
    type(finite_element_t)                    , intent(inout) :: finite_element
    ! Locals
    integer(ip)                 :: ndime,nnodu,nnodp,ngaus
    integer(ip)                 :: igaus,inode,jnode,idof,jdof,nstge,jstge,idime,jdime
    real(rp)                    :: dvolu,dtinv,diffu,alpha,ctime,prevtime,beta
    real(rp)                    :: work(4),aux
    real(rp)                    :: agran(finite_element%integ(1)%p%uint_phy%nnode)
    real(rp)      , allocatable :: testf(:,:,:)
    real(rp)      , allocatable :: tau(:,:,:)
    real(rp)      , allocatable :: elmat_vu_diag(:,:)
    real(rp)      , allocatable :: elvec_u(:,:) 
    type(vector_t), allocatable :: gpvel(:),gposs(:),force(:),grpre(:)
    type(tensor_t), allocatable :: grvel(:)
    type(rungekutta_integrator_t), pointer :: rkinteg

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    nnodp = finite_element%integ(ndime+1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    rkinteg => approx%rkinteg
    nstge = rkinteg%rk_table(1)%p%stage
    dtinv = rkinteg%dtinv
    ctime = rkinteg%ctime
    prevtime = rkinteg%ctime - 1.0_rp/dtinv

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    allocate(gpvel(nstge+1))
    allocate(gposs(nstge))
    allocate(grvel(nstge))
    allocate(grpre(nstge))
    do jstge=1,nstge
       call create_vector(approx%physics,1,finite_element%integ,gpvel(jstge))
       call create_vector(approx%physics,1,finite_element%integ,gposs(jstge))
       call create_tensor(approx%physics,1,ndime,finite_element%integ,grvel(jstge))
       call create_vector(approx%physics,1,finite_element%integ,grpre(jstge))
    end do
    call create_vector (approx%physics,1,finite_element%integ,gpvel(nstge+1))
    call interpolation(finite_element%unkno,1,prev_step,finite_element%integ,gpvel(1))                    ! U_n
    do jstge=2,nstge+1
       call interpolation(finite_element%unkno,1,jstge+2,finite_element%integ,gpvel(jstge))               ! U_j
       call interpolation(finite_element%unkno,ndime+2,jstge+2,finite_element%integ,gposs(jstge-1))       ! X_j
       call interpolation(finite_element%unkno,1,jstge+2,finite_element%integ,grvel(jstge-1))             ! GradU_j
       call interpolation(finite_element%unkno,ndime+1,ndime,jstge+2,finite_element%integ,grpre(jstge-1)) ! GradP_j
    end do

    ! Allocate auxiliar matrices and vectors
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)
    call memalloc(ndime,nnodu,elvec_u,__FILE__,__LINE__)

    ! Allocate & compute Stabilization parameters
    call memalloc(nnodu,ngaus,nstge,testf,__FILE__,__LINE__)
    call memalloc(2,ngaus,nstge,tau,__FILE__,__LINE__)
    do jstge=2,nstge+1
       call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel(jstge)%a,tau(:,:,jstge-1))
    end do

    ! Allocate & compute force term
    allocate(force(nstge))
    do jstge=1,nstge
       ctime = prevtime + 1.0_rp/dtinv*rkinteg%rk_table(1)%p%c(jstge)
       call create_vector(approx%physics,1,finite_element%integ,force(jstge))
       force(jstge)%a=0.0_rp
       ! Impose analytical solution
       if(finite_element%p_analytical_code%a(1,1)>0.and.finite_element%p_analytical_code%a(ndime+1,1)>0) then 
          call nsi_analytical_force(approx%physics,finite_element,ctime,gpvel(jstge+1),force(jstge))
       end if
       ! Add external force term
       do igaus=1,ngaus
          force(jstge)%a(:,igaus) = force(jstge)%a(:,igaus) + approx%physics%gravi(:)
       end do
    end do

    ! Initialize to zero
    elmat_vu_diag = 0.0_rp
    elvec_u       = 0.0_rp
    work          = 0.0_rp
    agran         = 0.0_rp
    testf         = 0.0_rp

    ! Loop on Gauss points
    gauss: do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Construct LHS
       ! =============
       ! Mass matrix (M/dt)
       call elmmss_gal(dvolu,dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),nnodu, &
            &          elmat_vu_diag,work)

       ! Construct RHS
       ! =============
       ! Mass
       ! ( v, u_n/dt )
       beta = dvolu*dtinv
       call elmrhs_fce(beta,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
            &          gpvel(1)%a(:,igaus),nnodu,ndime,elvec_u,work)

       ! Loop over previous stages
       stages: do jstge=1,nstge

          ! Auxiliar variables
          if(approx%physics%kfl_conv.ne.0) then
             do inode = 1,nnodu
                agran(inode) = 0.0_rp
                do idime = 1,ndime
                   agran(inode) = agran(inode) + gpvel(jstge+1)%a(idime,igaus) * &
                        &         finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
                end do
                testf(inode,igaus,jstge) = tau(1,igaus,jstge)*agran(inode)
             end do
          end if

          ! Diffusion
          alpha = rkinteg%rk_table(1)%p%b(jstge)
          beta = -dvolu*alpha
          ! nu( grad u, grad v )
          call elmrhu_gradgrad(beta,diffu,grvel(jstge)%a(:,:,igaus), &
               &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
          ! Add cross terms for symmetric grad
          if(approx%physics%kfl_symg==1) then
             call elmrhu_gradgrad_sym(beta,diffu,grvel(jstge)%a(:,:,igaus), &
                  &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
          end if

          ! Convection
          if(approx%physics%kfl_conv.ne.0) then
             alpha = rkinteg%rk_table(2)%p%b(jstge)
             beta = -dvolu*alpha
             ! (u · grad u, v) 
             if(approx%physics%kfl_skew==0) then
                call oss_convu_arg_chk(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus), &
                     &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_u,work)
             elseif(approx%physics%kfl_skew==1) then
                call oss_convu_arg_chk_skew(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus),agran, &
                     &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
                     &                 finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elvec_u,work)
             end if
          end if

          ! Pressure gradient
          alpha = rkinteg%rk_table(3)%p%b(jstge)
          beta = -dvolu*alpha
          ! (grad p, v)
          call oss_gradp_arg_chk(beta,grpre(jstge)%a(:,igaus), &
               &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_u,work)

          ! OSS_vu
          alpha = rkinteg%rk_table(4)%p%b(jstge)
          beta = -dvolu*alpha*tau(1,igaus,jstge)
          ! tau * ( u · grad u, u· grad v)
          call oss_convu_arg_chk(beta,gpvel(jstge+1)%a(:,igaus),grvel(jstge)%a(:,:,igaus),agran, &
               &                 ndime,nnodu,elvec_u,work)
          ! tauc*(div v, div u)
          if(approx%discret%ktauc>0.0_rp) then
             beta = -dvolu*alpha*tau(2,igaus,jstge)
             call elmrhu_divudivv(beta,grvel(jstge)%a(:,:,igaus), &
                  &               finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
                  &               elvec_u,work)
          end if

          ! OSS_vx
          alpha = rkinteg%rk_table(5)%p%b(jstge)
          beta = dvolu*alpha
          ! - tau * ( proj(u·grad u), u· grad v)
          if(approx%discret%kfl_proj==1) then
             call elmrhs_fce(beta,testf(:,igaus,jstge),gposs(jstge)%a(:,igaus),nnodu,ndime,elvec_u, &
                  &          work)
          else
             call elmrhs_fce(beta,agran,gposs(jstge)%a(:,igaus),nnodu,ndime,elvec_u,work)
          end if

          ! Force
          alpha = rkinteg%rk_table(6)%p%b(jstge)
          beta = dvolu*alpha
          call elmrhs_fce(beta,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
               &          force(jstge)%a(:,igaus),nnodu,ndime,elvec_u,work)

       end do stages

    end do gauss

    ! Deallocate
    do jstge=1,nstge
       call memfree(gpvel(jstge)%a,__FILE__,__LINE__)
       call memfree(gposs(jstge)%a,__FILE__,__LINE__)
       call memfree(grvel(jstge)%a,__FILE__,__LINE__)
       call memfree(grpre(jstge)%a,__FILE__,__LINE__)
       call memfree(force(jstge)%a,__FILE__,__LINE__)
    end do
    call memfree(gpvel(nstge+1)%a,__FILE__,__LINE__)
    deallocate(gpvel)
    deallocate(gposs)
    deallocate(grvel)
    deallocate(grpre)
    deallocate(force)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             ! Block V-U (diag)
             idof = finite_element%start%a(idime)+inode-1
             jdof = finite_element%start%a(idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) +  elmat_vu_diag(inode,jnode)
          end do
          ! Block U
          idof = finite_element%start%a(idime)+inode-1
          finite_element%p_vec%a(idof) = finite_element%p_vec%a(idof) + elvec_u(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu_diag,__FILE__,__LINE__)
    call memfree(elvec_u,__FILE__,__LINE__)
    call memfree(testf,__FILE__,__LINE__)
    call memfree(tau,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine nsi_rk_momentum_update

  !=================================================================================================
  subroutine nsi_rk_projection_update(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental projection matrix integration for projection dofs.  !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_rk_projection_update_t), intent(inout) :: approx
    type(finite_element_t)                      , intent(inout) :: finite_element
    ! Locals
    integer(ip)           :: ndime,nnodu,ngaus
    integer(ip)           :: igaus,inode,jnode,idof,jdof,idime,jdime
    real(rp)              :: dvolu,beta
    real(rp)              :: work(4)
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    real(rp), allocatable :: elmat_wx_diag(:,:)
    real(rp), allocatable :: elvec_x(:,:) 
    type(vector_t)        :: gpvel
    type(tensor_t)        :: grvel
    type(rungekutta_integrator_t), pointer :: rkinteg

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    call create_vector (approx%physics,1,finite_element%integ,gpvel)
    call create_tensor(approx%physics,1,ndime,finite_element%integ,grvel)
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpvel) 
    call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,grvel)

    ! Allocate auxiliar matrices and vectors
    call memalloc(nnodu,nnodu,elmat_wx_diag,__FILE__,__LINE__)
    call memalloc(ndime,nnodu,elvec_x,__FILE__,__LINE__)

    ! Allocate & compute Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)

    ! Initialize to zero
    elmat_wx_diag = 0.0_rp
    elvec_x       = 0.0_rp
    work          = 0.0_rp

    ! Loop on Gauss points
    gauss: do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Construct LHS
       ! =============
       ! OSS_wx
       ! tau*(proj(a·grad u), w)
       if(approx%discret%kfl_proj==1) then
          call elmmss_gal(dvolu,tau(1,igaus),finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
               &          nnodu,elmat_wx_diag,work)
       else
          call elmmss_gal(dvolu,1.0_rp,finite_element%integ(1)%p%uint_phy%shape(:,igaus), &
               &          nnodu,elmat_wx_diag,work)
       end if


       ! Construct RHS
       ! =============
       ! OSS_wu
       ! -tau*(a·grad u, w)
       beta = tau(1,igaus)*dvolu
       call oss_convu_arg_chk(beta,gpvel%a(:,igaus),grvel%a(:,:,igaus), &
            &                 finite_element%integ(1)%p%uint_phy%shape(:,igaus),ndime,nnodu,elvec_x,work)

    end do gauss

    ! Deallocate
    call memfree(gpvel%a,__FILE__,__LINE__)
    call memfree(grvel%a,__FILE__,__LINE__)

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             ! Block W-X
             idof = finite_element%start%a(ndime+1+idime)+inode-1
             jdof = finite_element%start%a(ndime+1+idime)+jnode-1
             finite_element%p_mat%a(idof,jdof) = finite_element%p_mat%a(idof,jdof) + elmat_wx_diag(inode,jnode)
          end do
          ! Block X
          idof = finite_element%start%a(ndime+1+idime)+inode-1
          finite_element%p_vec%a(idof) = finite_element%p_vec%a(idof) + elvec_x(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_wx_diag,__FILE__,__LINE__)
    call memfree(elvec_x,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine nsi_rk_projection_update

  !=================================================================================================
  subroutine nsi_dissipation(approx,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental dissipation parameters integration.                 !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_dissipation_t), intent(inout) :: approx
    type(finite_element_t)             , intent(inout) :: finite_element
    ! Locals
    integer(ip)           :: ndime,nnodu,ngaus
    integer(ip)           :: igaus,inode,jnode,idime,jdime
    real(rp)              :: dvolu,dtinv,diffu,gpvno
    real(rp)              :: tau(2,finite_element%integ(1)%p%quad%ngaus),auxva(3,2)
    real(rp)              :: divis,dicon,dinum,ditim,difex,energ,enstr,uxdx3,uxdx2,diver
    real(rp)              :: testu( approx%physics%ndime),ladjo( approx%physics%ndime)
    real(rp)              :: divel,strat,cnvec,extfo,dtime,rot_u
    type(vector_t)        :: gpvel,gpveln,grpre,gposs,force
    type(tensor_t)        :: grvel

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    dtinv = approx%tinteg%dtinv
    diffu = approx%physics%diffu

    ! Initialize to zero
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp

    ! Interpolation operations
    call create_vector(approx%physics,1,finite_element%integ,gpvel)
    call create_vector(approx%physics,1,finite_element%integ,gpveln)
    call create_vector(approx%physics,1,finite_element%integ,grpre)
    call create_vector(approx%physics,1,finite_element%integ,gposs)
    call create_tensor(approx%physics,1,ndime,finite_element%integ,grvel)
    call interpolation(finite_element%unkno,1,current,finite_element%integ,gpvel) 
    call interpolation(finite_element%unkno,1,current,finite_element%integ,grvel)
    call interpolation(finite_element%unkno,ndime+1,ndime,current,finite_element%integ,grpre)
    call interpolation(finite_element%unkno,ndime+2,current,finite_element%integ,gposs)
    if(dtinv == 0.0_rp) then
       gpveln%a = 0.0_rp
    else
       call interpolation (finite_element%unkno, 1, prev_step, finite_element%integ, gpveln)
    end if

    ! Allocate & compute Stabilization parameters
    call nsi_elmvsg(approx%physics,approx%discret,finite_element,gpvel%a,tau)

    ! Set force term
    call create_vector(approx%physics,1,finite_element%integ,force)
    force%a=0.0_rp
    ! Impose analytical solution
    if(finite_element%p_analytical_code%a(1,1)>0.and.finite_element%p_analytical_code%a(ndime+1,1)>0) then 
       call nsi_analytical_force(approx%physics,finite_element,approx%tinteg%ctime,gpvel,force)
    end if
    ! Add external force term
    do igaus=1,ngaus
       force%a(:,igaus) = force%a(:,igaus) + approx%physics%gravi(:)
    end do   

    ! Auxiliar vector for vorticity
    if(ndime==3) then
       auxva(1,:) = (/2,3/)
       auxva(2,:) = (/1,3/)
       auxva(3,:) = (/1,2/)
    end if

    ! Initialize to zero
    divis = 0.0_rp       ! Viscous term 
    dicon = 0.0_rp       ! Convective term
    dinum = 0.0_rp       ! Numerical dissipation
    ditim = 0.0_rp       ! Time derivative term (u_h) 
    difex = 0.0_rp       ! External force term
    energ = 0.0_rp       ! Global energy
    enstr = 0.0_rp       ! Enstrophy
    uxdx3 = 0.0_rp       
    uxdx2 = 0.0_rp
    diver = 0.0_rp       ! Divergence norm

    ! Loop on Gauss points
    gauss: do igaus = 1,ngaus
       dvolu = finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

       ! Velocity norm at gauss point
       gpvno=0.0_rp
       if(approx%physics%kfl_conv.ne.0) then
          do idime=1,ndime
             gpvno = gpvno + gpvel%a(idime,igaus)*gpvel%a(idime,igaus)
          end do
          gpvno = sqrt(gpvno)
       end if

       ! Initializtion
       testu = 0.0_rp
       divel = 0.0_rp
       ladjo = 0.0_rp
       strat = 0.0_rp
       cnvec = 0.0_rp
       extfo = 0.0_rp
       dtime = 0.0_rp
       rot_u = 0.0_rp

       ! Test operator and its adjoint acting on u
       ! =========================================
       ! Convective term (a·grad u)
       do idime = 1,ndime 
          do jdime = 1,ndime
             testu(idime) = testu(idime) + gpvel%a(jdime,igaus)*grvel%a(jdime,idime,igaus)
             ladjo(idime) = ladjo(idime) - gpvel%a(jdime,igaus)*grvel%a(jdime,idime,igaus) 
          end do
          divel = divel + grvel%a(idime,idime,igaus)
       end do

       ! Compute Auxiliars
       ! =================
       do idime = 1,ndime 
          do jdime = 1,ndime
             ! (grad u, grad u)
             strat = strat + grvel%a(idime,jdime,igaus)*(grvel%a(idime,jdime,igaus))
             if(approx%physics%kfl_symg==1) then
                strat = strat + grvel%a(idime,jdime,igaus)*grvel%a(jdime,idime,igaus)
             end if
             ! (a · grad u, u)
             if(approx%physics%kfl_skew==0) then
                cnvec = cnvec + gpvel%a(idime,igaus)*grvel%a(idime,jdime,igaus)*gpvel%a(jdime,igaus)
             end if
          end do
          ! - (f,u)
          extfo = extfo - gpvel%a(idime,igaus)*force%a(idime,igaus)
          ! (du/dt,u) ~ 1/(0.5*dt)*(u^2 - un^2)
          dtime = dtime + 0.5_rp*dtinv*((gpvel%a(idime,igaus))**2-(gpveln%a(idime,igaus))**2)
          ! Vorticity (|rot(u)|^2)
          if(ndime==3) then
             rot_u = rot_u + (grvel%a(auxva(idime,1),auxva(idime,2),igaus) - &
                  &           grvel%a(auxva(idime,2),auxva(idime,1),igaus))**2
          end if
          ! Numerical dissipation ((a·grad u)-Proj(a·grad u), a·grad u)
          if(approx%discret%kfl_proj==1) then
             dinum = dinum + tau(1,igaus)*(testu(idime)-gposs%a(idime,igaus))*testu(idime)*dvolu
          else
             dinum = dinum + (tau(1,igaus)*testu(idime)-gposs%a(idime,igaus))*testu(idime)*dvolu
          end if
       end do
       ! Vorticity (|rot(u)|^2)
       if(ndime==2) then
          rot_u = rot_u + (grvel%a(1,2,igaus)-grvel%a(2,1,igaus))**2
       end if
       ! Add tauc*(div u, div u)
       if(approx%discret%ktauc>0.0_rp) then
          dinum = dinum + tau(2,igaus)*divel*divel*dvolu
       end if

       ! Compute dissipations
       ! ====================
       divis = divis + diffu*strat*dvolu   
       dicon = dicon + cnvec*dvolu
       ditim = ditim + dtime*dvolu
       difex = difex + extfo*dvolu
       enstr = enstr + 0.5_rp*rot_u*dvolu
       do idime=1,ndime
          energ = energ + 0.5_rp*gpvel%a(idime,igaus)*gpvel%a(idime,igaus)*dvolu
       end do
       uxdx3 = uxdx3 + (grvel%a(1,1,igaus)**3)*dvolu
       uxdx2 = uxdx2 + (grvel%a(1,1,igaus)**2)*dvolu
       diver = diver + divel**2*dvolu

    end do gauss

    ! Deallocate
    call memfree(gpvel%a,__FILE__,__LINE__)
    call memfree(gpveln%a,__FILE__,__LINE__)
    call memfree(gposs%a,__FILE__,__LINE__)
    call memfree(grvel%a,__FILE__,__LINE__)
    call memfree(grpre%a,__FILE__,__LINE__)
    call memfree(force%a,__FILE__,__LINE__)

    ! Add contributions
    ! =================
    ! Energy balance quantities
    approx%values(1) = approx%values(1) + dinum      ! tau * (a · grad u) * ( (a · grad u) - Proj(a · grad u) )
    approx%values(2) = approx%values(2) + divis      ! nu * ||grad u||^2
    approx%values(3) = approx%values(3) + dicon      ! b (a , u, u)
    approx%values(4) = approx%values(4) + ditim      ! 1/2 * d/dt(||u||^2)
    approx%values(5) = approx%values(5) + difex      ! - (f,u)
    approx%values(6) = approx%values(6) + energ      ! 1/2 * (u,u)
    ! Turbulent quantities
    approx%values(7) = approx%values(7) + enstr      ! 1/2 * (w,w)
    approx%values(8) = approx%values(8) + uxdx3      ! < (du/dx)^3 >
    approx%values(9) = approx%values(9) + uxdx2      ! < (du/dx)^2 >
    ! Other quantities
    approx%values(10) = approx%values(10) + diver    ! (div u, div u)

  end subroutine nsi_dissipation

  !==================================================================================================
  subroutine nsi_elmvsg(physics,discret,finite_element,gpvel,tau)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine computes the stabilization parameters.                                     !
    !----------------------------------------------------------------------------------------------!
    !nnode,ndime,approx,hleng,ngaus,jainv,shape,elvel,gpvel,tau,difsma) 
    implicit none
    type(nsi_problem_t)            , intent(in)  :: physics
    type(nsi_cg_iss_oss_discrete_t), intent(in)  :: discret
    type(finite_element_t)         , intent(in)  :: finite_element
    real(rp)                       , intent(in)  :: gpvel(physics%ndime,finite_element%integ(1)%p%quad%ngaus)
    real(rp)                       , intent(out) :: tau(2,finite_element%integ(1)%p%quad%ngaus)
    ! Locals
    integer(ip) :: ngaus,ndime,nnodu
    integer(ip) :: igaus,idime
    real(rp)    :: alpha,gpvno,diffu
    real(rp)    :: chave(physics%ndime,2),chale(2)
    real(rp)    :: veloc(finite_element%integ(1)%p%uint_phy%nnode,physics%ndime)

    ! Unpack variables
    ndime = physics%ndime
    nnodu = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%quad%ngaus
    diffu = physics%diffu

    ! Initialize
    tau = 0.0_rp

    do igaus=1,ngaus

       ! Velocity norm at gauss point
       gpvno=0.0_rp
       do idime=1,ndime
          gpvno = gpvno + gpvel(idime,igaus)*gpvel(idime,igaus)
       end do
       gpvno = sqrt(gpvno)

       ! Compute the characteristic length chale
       veloc = finite_element%unkno(1:nnodu,1:ndime,1)
       call nsi_elmchl(physics%kfl_chale,finite_element%integ(1)%p%femap%jainv,          &
            &          finite_element%integ(1)%p%femap%hleng(:,igaus),veloc,ndime,nnodu, &
            &          physics%kfl_conv,chave,chale)
       
       ! Auxiliar computations
       alpha  = discret%k1tau*diffu/(chale(2)*chale(2)) + discret%k2tau*gpvno/chale(1) + &
            &   1.0_rp*physics%react  ! Old
       !alpha  = k1tau*diffu/(chale(2)*chale(2)) + k2tau*facto*gpvno/chale(1) + 1.0_rp*react   ! Codina-Guasch

       ! NS parameters
       if(alpha.gt.1e-8) tau(1,igaus) = 1.0_rp/(alpha)
       if(discret%k1tau.gt.1e-8) then
          tau(2,igaus) = discret%ktauc*(diffu+discret%k2tau/discret%k1tau*gpvno*chale(2))
       end if
    end do

  end subroutine nsi_elmvsg

  !==================================================================================================
  subroutine nsi_vars_block(discret,physics,vars_block)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generate the vars per block array needed for dof_handler creation.          !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_discrete_t), intent(in)  :: discret
    type(nsi_problem_t)             , intent(in)  :: physics
    integer(ip), allocatable        , intent(out) :: vars_block(:)
    ! Locals
    integer(ip) :: idime

    call memalloc(discret%nvars,vars_block,__FILE__,__LINE__)

    ! Block U
    do idime=1,physics%ndime
       vars_block(idime) = 1
    end do
    ! Block P
    vars_block(physics%ndime+1) = 2
    ! Block X
    do idime=1,physics%ndime
       vars_block(physics%ndime+1+idime) = 3
    end do

  end subroutine nsi_vars_block

  !==================================================================================================
  subroutine nsi_dof_coupling(discret,physics,dof_coupling)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generate the dof coupling array needed for dof_handler creation.            !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_oss_discrete_t), intent(in)  :: discret
    type(nsi_problem_t)             , intent(in)  :: physics
    integer(ip), allocatable        , intent(out) :: dof_coupling(:,:)
    ! Locals
    integer(ip) :: idime,jdime

    call memalloc(discret%nvars,discret%nvars,dof_coupling,__FILE__,__LINE__)
    dof_coupling = 0

    do idime=1,physics%ndime
       do jdime=1,physics%ndime
          ! Block V-U (all)
          dof_coupling(idime,jdime) = 1
          ! Block V-X (all)
          dof_coupling(idime,physics%ndime+1+jdime) = 1
          ! Block W-U (all)
          dof_coupling(physics%ndime+1+idime,jdime) = 1
       end do
       ! Block W-X (diag)
       dof_coupling(physics%ndime+1+idime,physics%ndime+1+idime) = 1
       ! Block V-P
       dof_coupling(idime,physics%ndime+1) = 1
       ! Block Q-U
       dof_coupling(physics%ndime+1,idime) = 1
    end do  
    ! Block Q-P
    dof_coupling(physics%ndime+1,physics%ndime+1) = 1  ! Needed to integrate Mass_p matrix (preconditioner)
    
  end subroutine nsi_dof_coupling    
  
end module nsi_cg_iss_oss_names
