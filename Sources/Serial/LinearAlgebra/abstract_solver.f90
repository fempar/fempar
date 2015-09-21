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
!=============================================================================
module abstract_solver_names
  use types_names
  use memor_names
  use operator_names
  use vector_names
  use abstract_environment_names

# include "debug.i90"

  implicit none


  ! List of convergence criteria available for iterative solvers 
  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration


  integer (ip), parameter :: mgs  = 1 ! mgs : Modified Gram-Schmidt 
                                      !       (appropriate for serial GMRES)
  integer (ip), parameter :: icgs = 2 ! icgs: Iterative Classical Gram-Schmidt 
                                      !       (appropriate for distributed GMRES)


  logical     , parameter :: default_track_conv_his  = .false.
  real    (rp), parameter :: default_relax           = 1.0_rp
  real    (rp), parameter :: default_atol            = 0.0_rp
  real    (rp), parameter :: default_rtol            = 1.0e-08_rp
  integer (ip), parameter :: default_luout           = 6
  integer (ip), parameter :: default_stopc           = res_nrmgiven_res_nrmgiven
  integer (ip), parameter :: default_trace           = 0
  integer (ip), parameter :: default_itmax           = 1000
  integer (ip), parameter :: default_dkrymax         = 30
  integer (ip), parameter :: default_orto            = icgs
  integer (ip), parameter :: default_nrhs            = 1

  ! List of Krylov subspace methods available
  integer (ip), parameter :: cg               = 1
  integer (ip), parameter :: lgmres           = 2
  integer (ip), parameter :: rgmres           = 3
  integer (ip), parameter :: fgmres           = 4
  integer (ip), parameter :: icg              = 7
  integer (ip), parameter :: lfom             = 8  ! Left preconditioned Full Orthogonalization Method
  integer (ip), parameter :: minres           = 9  ! Preconditioned MINimal RESidual method

  ! Not actually Krylov methods
  integer (ip), parameter :: richard          = 5  ! Richardson (fixed-point iteration)
  integer (ip), parameter :: direct           = 6  ! Apply preconditioner directly


  integer (ip), parameter :: default_kry_meth = direct

  character(len=*), parameter  :: methdname(9)=(/ &
       & 'PCG       ', &
       & 'LGMRES    ', &
       & 'RGMRES    ', &
       & 'FGMRES    ', &
       & 'RICHARDSON', &
       & 'DIRECT    ', &
       & 'IPCG      ', &
       & 'LFOM      ', &
       & 'MINRES    ' /)

  ! Control data (always initialized to default values)
  type solver_control_t
     ! Is the solver machinery going to track the convergence history?
     logical :: track_conv_his = default_track_conv_his

     ! Input parameters
     real(rp)      :: relax    = default_relax              ! Relaxation parameter
     real(rp)      :: rtol     = default_rtol               ! Relative tolerance
     real(rp)      :: atol     = default_atol               ! Absolute tolerance
     integer(ip)   :: dkrymax  = default_dkrymax            ! Maximum dimension of the Krylov subspace
     integer(ip)   :: orto     = default_orto               ! Ortonormalization strategy (mgs,icgs)
     integer(ip)   :: luout    = default_luout              ! Logical unit to output
     integer(ip)   :: stopc    = default_stopc              ! Stopping criteria
     integer(ip)   :: trace    = default_trace              ! Message every trace iterations
     integer(ip)   :: itmax    = default_itmax              ! Max. # of iterations
     integer(ip)   :: method   = default_kry_meth           ! Krylov subspace method
     integer(ip)   :: nrhs     = default_nrhs               ! Number of simultaneous RHS
     ! Output parameters
     logical               :: converged = .false.           ! Converged?
     integer(ip)           :: it        = 0                 ! # of iterations to converge
     real(rp)              :: bn2       = 0.0_rp            ! RHS 2-norm
     real(rp)              :: rn2       = 0.0_rp            ! Residual 2-norm
     real(rp)              :: dxn2      = 0.0_rp            ! dx 2-norm
     real(rp)              :: tol1      = default_atol      ! Tolerance 1
     real(rp)              :: tol2      = default_atol      ! Tolerance 2
     real(rp)              :: err1      = 0.0_rp            ! Current error estimate 1
     real(rp)              :: err2      = 0.0_rp            ! Current error estimate 2
     real(rp), allocatable :: err1h(:)                      ! Error estimates 1 history
     real(rp), allocatable :: err2h(:)                      ! Error estimates 2 history
  end type solver_control_t


  private
  
  ! Subroutines
  public :: abstract_solve

  ! Data types
  public :: solver_control_t
  public :: solver_control_log_conv_his, solver_control_free_conv_his

  ! Named constants
  public :: cg, lgmres, rgmres, fgmres, icg, lfom, minres, richard, direct
  public :: mgs, icgs
  public :: res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven, delta_rhs, delta_delta
  public :: res_res, res_rhs, delta_rhs_and_res_res, delta_rhs_and_res_rhs, delta_delta_and_res_res, delta_delta_and_res_rhs

contains

  !=============================================================================
  ! Abstract Solve (public driver)
  !=============================================================================
  subroutine abstract_solve(A,M,b,x,ctrl,env)
    implicit none
    class(operator_t), intent(in)        :: A        ! Matrix
    class(operator_t), intent(in)        :: M        ! Preconditioner
    class(vector_t) , intent(in)        :: b        ! RHS
    class(vector_t) , intent(inout)     :: x        ! Approximate solution
    type(solver_control_t), intent(inout)     :: ctrl     ! Solver parameters
    class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 
    
    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()
    
    call solver_control_allocate_conv_his ( ctrl )
    ctrl%converged = .false.

    ! Invoke Krylov Subspace Solver
    select case( ctrl%method )
    case ( cg )
        call abstract_pcg ( A, M, b, x, ctrl, env )
    case ( lgmres )
        call abstract_plgmres ( A, M, b, x, ctrl, env )
    case ( rgmres )
        call abstract_prgmres ( A, M, b, x, ctrl, env )
    case ( fgmres )
        call abstract_pfgmres ( A, M, b, x, ctrl, env )
    case ( richard )
        call abstract_prichard ( A, M, b, x, ctrl, env )
    case( direct )
        call M%apply(b, x)
    case ( icg )
        call abstract_ipcg ( A, M, b, x, ctrl, env )
    case ( lfom )
        call abstract_plfom ( A, M, b, x, ctrl, env )
    case ( minres )
        call abstract_pminres ( A, M, b, x, ctrl, env )
    case default
       ! Write an error message and stop ?      
    end select

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
  end subroutine abstract_solve

  !=============================================================================
  ! Abstract Preconditioned Conjugate Gradient
  !=============================================================================
  subroutine abstract_pcg( A, M, b, x, ctrl, env )
    !-----------------------------------------------------------------------------
    ! This routine performs pcg iterations on Ax=b with preconditioner M. 
    !-----------------------------------------------------------------------------
    implicit none

    ! Parameters
    class(operator_t) , intent(in)       :: A     ! Matrix
    class(operator_t) , intent(in)       :: M     ! Preconditioner
    class(vector_t)  , intent(in)       :: b     ! RHS
    class(vector_t)  , intent(inout)    :: x     ! Approximate solution
    type(solver_control_t) , intent(inout)    :: ctrl  ! Control data
    class(abstract_environment_t), intent(in) :: env   ! Serial/parallel environment 


    ! Locals
    real(rp)                         :: r_nrm_M      ! |r|_inv(M) 
    real(rp)                         :: b_nrm_M      ! |b|_inv(M)
    real(rp)                         :: r_z, Ap_p, alpha, beta
    integer                          :: me, np
    class(vector_t), allocatable :: r,p,Ap,z     ! Working vector_ts

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    allocate ( r, mold=x ); call r%default_initialization()
    allocate ( Ap, mold=x ); call Ap%default_initialization()
    allocate ( z, mold=x ); call z%default_initialization()
    allocate ( p, mold=x ); call p%default_initialization()

    call r%clone ( x  )
    call Ap%clone ( x )
    call z%clone ( x  )
    call p%clone ( x  )

    call env%info(me,np)

    ! Evaluate |b|_inv(M) if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
       ! r = inv(M) b
       call M%apply(b, r)
       b_nrm_M = b%dot(r)
       if ( env%am_i_fine_task() ) then ! Am I a fine task ?
          b_nrm_M = sqrt(b_nrm_M)
       end if
    else
       b_nrm_M = 0.0_rp
    endif

    ! 1) Compute initial residual
    ! 1.a) r=Ax
    call A%apply(x,r)

    ! 1.b) r=b-r
    call r%axpby(1.0, b ,-1.0)

    ! 2) z=inv(M)r
    call M%apply(r,z)


    ! 3) <r,z>
    r_z = r%dot(z) 

    if ( env%am_i_fine_task() ) then ! Am I a fine task ?
       r_nrm_M = sqrt( r_z )
    end if

    ! 4) Initializations:
    ! p=z
    call p%copy(z) 

    if ( env%am_i_fine_task() ) then
       ! Init and log convergence 
       call pcg_conv_init (  b, r, b_nrm_M, r_nrm_M, ctrl )
       if ((me == 0).and.(ctrl%trace/=0))  call solver_control_log_header(ctrl)
    end if

    ! 5) Iteration
    ctrl%it = 0
    loop_pcg: do
       ctrl%it = ctrl%it + 1

       ! Ap = A*p
       call A%apply(p,Ap)

       ! <Ap,p>
       Ap_p = Ap%dot(p)

       ! Is this correct/appropriate ?
       if (Ap_p /= 0.0_rp) then
          alpha = r_z / Ap_p
       else
          alpha = 0.0_rp
       end if

       ! x = x + alpha*p
       call x%axpby(alpha, p, 1.0)

       ! r = r - alpha*Ap
       call r%axpby (-alpha, Ap, 1.0 )

       if (env%am_i_fine_task()) then
          ! Check and log convergence
          call pcg_conv_check(r, r_nrm_M, alpha, p, ctrl )
          if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
       end if
       ! Send converged to coarse-grid tasks
       call env%bcast(ctrl%converged)
       if(ctrl%converged.or.(ctrl%it>=ctrl%itmax)) exit loop_pcg

       ! z = inv(M) r
       call M%apply(r,z)

       if ( env%am_i_fine_task()) then ! Am I a fine task ?
          beta = 1.0_rp/r_z
       else
          beta = 0.0_rp 
       end if

       r_z=r%dot(z)
       beta = beta*r_z

       if (env%am_i_fine_task()) then ! Am I a fine task ?
          r_nrm_M = sqrt( r_z )
       end if

       ! p = z + beta*p
       call p%axpby(1.0,z,beta)
    end do loop_pcg

    call r%free()
    call Ap%free()
    call z%free()
    call p%free()

    deallocate ( r )
    deallocate ( Ap )
    deallocate ( z )
    deallocate ( p )


    if ( env%am_i_fine_task() ) then
       if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
  end subroutine abstract_pcg

  subroutine pcg_conv_init ( b, r , nrm_b_given, nrm_r_given, ctrl )
    implicit none 
    class(vector_t) , intent(in)   :: b, r
    real(rp)            , intent(in)    :: nrm_b_given, nrm_r_given 
    type(solver_control_t), intent(inout) :: ctrl
    

    select case(ctrl%stopc)
    case ( delta_rhs, res_rhs, delta_rhs_and_res_rhs, delta_delta_and_res_rhs )
       ! call generic_nrm2 ( b, ctrl%bn2 )
       ctrl%bn2 = b%nrm2()
       ctrl%tol1 = ctrl%rtol * ctrl%bn2 + ctrl%atol
       ctrl%tol2 = ctrl%tol1
    case ( res_res, delta_delta_and_res_res )
       ! call generic_nrm2 ( r, ctrl%rn2 )
       ctrl%rn2 = r%nrm2()
       ctrl%tol1 = ctrl%rtol* ctrl%rn2 + ctrl%atol
       ctrl%tol2 = ctrl%tol1
    case ( delta_rhs_and_res_res )
       ! call generic_nrm2 ( b, ctrl%bn2 )
       ctrl%bn2 = b%nrm2()
       ! call generic_nrm2 ( r, ctrl%rn2 )
       ctrl%rn2 = r%nrm2() 
       ctrl%tol1 = ctrl%rtol* ctrl%bn2 + ctrl%atol
       ctrl%tol2 = ctrl%rtol* ctrl%rn2 + ctrl%atol
    case ( res_nrmgiven_rhs_nrmgiven )
       ctrl%tol1 = ctrl%rtol* nrm_b_given + ctrl%atol
    case ( res_nrmgiven_res_nrmgiven )
       ctrl%tol1 = ctrl%rtol* nrm_r_given + ctrl%atol
    case default
       ! Write an error message and stop ?      
    end select

    return

  end subroutine pcg_conv_init

  subroutine pcg_conv_check ( r, nrm_r_given, alpha, p, ctrl )
    implicit none 
    ! Parameters
    class(vector_t) , intent(in)    :: r, p
    real(rp)            , intent(in)    :: nrm_r_given
    real(rp)            , intent(in)    :: alpha
    type(solver_control_t), intent(inout) :: ctrl

    ! Compute 1st iteration error estimate and upper bound 
    ! for convergence criteria depending on ||dx(i)||
    if ( ctrl%it == 1 ) then
       select case( ctrl%stopc )
       case ( delta_rhs, delta_delta,  delta_rhs_and_res_res, delta_rhs_and_res_rhs, &
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          ! call generic_nrm2 ( p, ctrl%err1 )
          ctrl%err1 = p%nrm2()
          ctrl%err1 = alpha * ctrl%err1
          select case( ctrl%stopc )
          case ( delta_delta, delta_delta_and_res_res, delta_delta_and_res_rhs ) 
             ctrl%dxn2 = ctrl%err1
             ctrl%tol1 = ctrl%rtol*ctrl%dxn2 + ctrl%atol
          end select
       end select
    end if

    ctrl%converged = .false.

    ! Evaluate 1st convergence criterion
    select case( ctrl%stopc )

    case ( res_res, res_rhs )
       ! Compute || r(i) ||
       !call generic_nrm2 ( r, ctrl%err1 )
       ctrl%err1 = r%nrm2()
    case ( delta_rhs, delta_delta, delta_rhs_and_res_res, delta_rhs_and_res_rhs, &
         & delta_delta_and_res_res,delta_delta_and_res_rhs )
       if ( ctrl%it /= 1 ) then ! if false, no need to evaluate ||dx(i)|| again
          ! Compute || dx(i) ||
          ! call generic_nrm2 ( p, ctrl%err1 )
          ctrl%err1 = p%nrm2()
          ctrl%err1 = alpha * ctrl%err1
       end if

    case ( res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
       ctrl%err1 = nrm_r_given
    end select
    if ( ctrl%track_conv_his ) ctrl%err1h(ctrl%it) = ctrl%err1
    ctrl%converged = (ctrl%err1 <= ctrl%tol1 )

    ! Evaluate 2nd convergence criterion
    select case( ctrl%stopc )
    case ( delta_rhs_and_res_res,   delta_rhs_and_res_rhs, &
         & delta_delta_and_res_res, delta_delta_and_res_rhs )
       if ( ctrl%converged ) then
          ! Compute || r(i) ||
          ! call generic_nrm2 ( r, ctrl%err2 )
          ctrl%err2 = r%nrm2()
          ctrl%converged = (ctrl%err2 <= ctrl%tol2 )
          if ( ctrl%track_conv_his ) ctrl%err2h(ctrl%it) = ctrl%err2
       end if
    end select

    return
  end subroutine pcg_conv_check


!=============================================================================
! Abstract Inexact Preconditioned Conjugate Gradient
! (Taken from Golub et. al paper 
!  Inexact Preconditioned Conjugate Gradient Method with Inner-Outer Iteration)
!=============================================================================
subroutine abstract_ipcg( A, M, b, x, ctrl, env )
  !-----------------------------------------------------------------------------
  ! This routine performs pcg iterations on Ax=b with preconditioner M. 
  !-----------------------------------------------------------------------------
  implicit none

  ! Parameters
  class(operator_t), intent(in)    :: A        ! Matrix
  class(operator_t), intent(in)    :: M        ! Preconditioner
  class(vector_t) , intent(in)    :: b        ! RHS
  class(vector_t) , intent(inout) :: x        ! Approximate solution
  type(solver_control_t), intent(inout) :: ctrl     ! Control data
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 

  ! Locals
  real(rp)           :: r_nrm_M      ! |r|_inv(M) 
  real(rp)           :: b_nrm_M      ! |b|_inv(M)
  real(rp)           :: r_z, r2_z, Ap_p, alpha, beta
  integer            :: me, np
  class(vector_t), allocatable  :: r,r2,p,Ap,z     ! Working vector_ts

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    allocate(r, mold=x); call r%default_initialization()
    allocate(r2, mold=x); call r2%default_initialization()
    allocate(Ap, mold=x); call Ap%default_initialization()
    allocate(z, mold=x); call z%default_initialization()
    allocate(p, mold=x); call p%default_initialization()
    call r%clone(x)
    call r2%clone(x)
    call Ap%clone(x)
    call z%clone(x)
    call p%clone(x)

    call env%info(me, np)

    ! Evaluate |b|_inv(M) if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
        ! r = inv(M) b
        call M%apply(b, r)
        b_nrm_M = b%dot(r)

        if ( env%am_i_fine_task() ) then ! Am I a fine task ?
            b_nrm_M = sqrt(b_nrm_M)
        end if
    else
     b_nrm_M = 0.0_rp
  endif

    ! 1) Compute initial residual
    ! 1.a) r=Ax
    call A%apply(x, r)

    ! 1.b) r=b-r
    call r%axpby(1.0_rp,b,-1.0_rp)

    ! 2) z=inv(M)r
    call M%apply(r, z)

    ! 3) <r,z>
    r_z = r%dot(z)

    if ( env%am_i_fine_task() ) then ! Am I a fine task ?
        r_nrm_M = sqrt( r_z )
    end if

    ! 4) Initializations:
    ! p=z
    call p%copy(z)

    if ( env%am_i_fine_task() ) then
     ! Init and log convergence 
        call pcg_conv_init (  b, r, b_nrm_M, r_nrm_M, ctrl )
        if ((me == 0).and.(ctrl%trace/=0))  call solver_control_log_header(ctrl)
    end if

    ! 5) Iteration
    ctrl%it = 0
    loop_pcg: do
        ctrl%it = ctrl%it + 1

        ! Ap = A*p
        call A%apply(p, Ap)

        ! <Ap,p>
        Ap_p = Ap%dot(p)

        ! Is this correct/appropriate ?
        if (Ap_p /= 0.0_rp) then
            alpha = r_z / Ap_p
        else
            alpha = 0.0_rp
        end if

        ! x = x + alpha*p
        call x%axpby(alpha,p,1.0_rp)

        call r2%copy(r)

        ! r = r - alpha*Ap
        call r%axpby(-alpha,Ap,1.0_rp)

        if ( env%am_i_fine_task() ) then
            ! Check and log convergence
            call pcg_conv_check(r, r_nrm_M, alpha, p, ctrl )
            if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
        end if
        ! Send converged to coarse-grid tasks
        call env%bcast(ctrl%converged)

        if(ctrl%converged.or.(ctrl%it>=ctrl%itmax)) exit loop_pcg

        ! z = inv(M) r
        call M%apply(r, z)

        if ( env%am_i_fine_task() ) then ! Am I a fine task ?
            beta = 1.0_rp/r_z
            r_z = r%dot(z)
            r2_z = r2%dot(z)
            beta = beta*(r_z-r2_z)
        else
            beta = 0.0_rp
        end if
     
        if (env%am_i_fine_task() ) then ! Am I a fine task ?
            r_nrm_M = sqrt( r_z )
        end if

        ! p = z + beta*p
        call p%axpby(1.0_rp,z,beta)
    end do loop_pcg

    call r%free()
    call r2%free()
    call Ap%free()
    call z%free()
    call p%free()

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()

end subroutine abstract_ipcg


!=============================================================================
! Abstract Left Preconditioned GMRES
!=============================================================================
!
! TODO:  1. Review orthogonalization procedures. (Done: MGS or ICGS)
!        2. Implement the whole list of convergence criteria available
!           for PCG. Only two are currently implemented.
!        3. Set err1/err2 accordingly to the work performed in
!           point 2 above.
!        4. Use level-2-BLAS for backward/forward substitution + 
!           ICGS. This last option would imply defining a generic
!           data structure for storing Krylov subspace bases.  (DONE)
!        5. Breakdown detection ? (DONE: Taken from Y. Saad's SPARSKIT)
! 
! Help to point 1: from Trilinos AztecOO 3.6 user's guide:
!
! As mentioned above, GMRES is by far the most robust general-purpose Krylov
! method available. Part of its robustness comes from the ability to tune two parameters,
! namely options[AZ kspace] and options[AZ orthog]. options[AZ kspace] determines
! the number of Krylov vectors that will be used as part of the least-squares
! problem to generate the next approximate solution. Generally setting this value
! larger improves the robustness, decreases iteration count, but increases costs.
! This value is set to 30 by default, but for challenging problems one should set it
! (much) higher, especially if memory is available on the computer. For very difficult
! problems, we recommend setting options[AZ kspace] equal to the maximum
! number of iterations.

! The parameter options[AZ orthog] can be used to select the type of Gram-Schmidt
! algorithm to used. We provide two options:
!   1. Two steps of classical Gram-Schmidt (Double CGS).
!   2. One step of modified Gram-Schmidt. (Single MGS).
! For many years, (single) MGS was the preferred option for GMRES because of
! its superior numerical accuracy over single CGS. However, in the past several
! years it has been recognized that double CGS is more effective than single MGS,
! as effective as double MGS and has superior parallel performance. As a result,
! AztecOO uses double CGS by default. However, there may be instances where
! single MGS would be sufficient for robustness and it can have a lower cost in some
! situations.
  subroutine abstract_plgmres ( A, M, b, x, ctrl, env ) 
    !--------------------------------------------------------------------
    ! This routine performs plgmres iterations on Ax=b with preconditioner M. 
    !--------------------------------------------------------------------
#ifdef ENABLE_BLAS       
use blas77_interfaces_names
#endif
    implicit none
    ! Parameters
    class(operator_t)  , intent(in) :: A ! Matrix
    class(operator_t)  , intent(in) :: M ! Preconditioner
    class(vector_t)   , intent(inout) :: x ! Solution
    class(vector_t)   , intent(in)    :: b ! RHS
    type(solver_control_t), intent(inout) :: ctrl
    class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 


    integer(ip)                    :: ierrc
    integer(ip)                    :: kloc, kloc_aux, max_kloc, i, j, k_hh, id
    real(rp)                       :: res_norm, res_2_norm, rhs_norm
    real(rp)                       :: alpha, c, s 
    real(rp)   , allocatable       :: hh(:,:), g(:), g_aux(:), cs(:,:)
    integer                        :: me, np
    logical                        :: exit_loop

    class(vector_t), allocatable  :: r, z     ! Working vector_ts
    class(vector_t), allocatable  :: bkry(:)  ! Krylov basis

    assert(ctrl%stopc==res_res.or.ctrl%stopc==res_rhs.or.ctrl%stopc==res_nrmgiven_res_nrmgiven.or.ctrl%stopc==res_nrmgiven_rhs_nrmgiven)

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    ! Clone x in order to allocate working vectors
    allocate(r, mold=b); call r%default_initialization()
    allocate(z, mold=x); call z%default_initialization()
    call r%clone(b)
    call z%clone(x)

    allocate(bkry(ctrl%dkrymax+1), mold=x)

    ! Allocate working vectors
    call memalloc(ctrl%dkrymax+1,ctrl%dkrymax+1,hh,__FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1,g,__FILE__,__LINE__)
    if ( ctrl%stopc == res_res .or. ctrl%stopc == res_rhs ) then
       call memalloc(ctrl%dkrymax+1,g_aux,__FILE__,__LINE__)
    end if
    call memalloc(2,ctrl%dkrymax+1,cs,__FILE__,__LINE__)

    call env%info(me, np)

    ! Evaluate ||b||_2 if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven .or. ctrl%stopc == res_rhs  ) then
       rhs_norm = b%nrm2() 
    endif
    
    ! r = Ax
    call A%apply(x,r)

    ! r = b-r
    call r%axpby(1.0_rp,b,-1.0_rp)

    if ( ctrl%stopc == res_res .or. ctrl%stopc == res_rhs  ) then
       res_2_norm = r%nrm2()
    end if

    ! z=inv(M)r
    call M%apply(r,z)

    ! Evaluate ||z||_2
    res_norm = z%nrm2()

    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
       ctrl%tol1  = ctrl%rtol * rhs_norm + ctrl%atol
       ctrl%err1 = res_norm 
    else if ( ctrl%stopc == res_nrmgiven_res_nrmgiven ) then
       ctrl%tol1 = ctrl%rtol * res_norm + ctrl%atol
       ctrl%err1 = res_norm
    else if ( ctrl%stopc == res_res ) then
       ctrl%tol1 = ctrl%rtol * res_2_norm + ctrl%atol
       ctrl%err1 = res_2_norm
    else if ( ctrl%stopc == res_rhs ) then
       ctrl%tol1 = ctrl%rtol * rhs_norm + ctrl%atol
       ctrl%err1 = res_2_norm
    end if
    exit_loop = (ctrl%err1 <= ctrl%tol1)

    ! Send converged to coarse-grid tasks
    call env%bcast(exit_loop )


    if ( env%am_i_fine_task() ) then
       if ( (ctrl%trace > 0) .and. (me == 0) ) call solver_control_log_header(ctrl)
    end if

    max_kloc = 0
    ctrl%it = 0
    outer: do while ( (.not.exit_loop).and.(ctrl%it<ctrl%itmax) )

       ! Compute preconditioned residual from scratch (only if ctrl%it/=0)
       if ( ctrl%it /= 0 ) then 
          ! r = Ax
          call A%apply(x,r)

          ! r = b-r
          call r%axpby(1.0_rp,b,-1.0_rp)

          ! z=inv(M)r
          call M%apply(r,z)  

          ! Evaluate ||z||_2
          res_norm = z%nrm2()
       end if

       ! Normalize preconditioned residual direction (i.e., v_1 = z/||z||_2)
       call bkry(1)%default_initialization()
       call bkry(1)%clone(x)
       if (res_norm /= 0.0_rp) call bkry(1)%scal(1.0_rp/res_norm, z)

       ! residual in the krylov basis
       g(1)            = res_norm
       g(2:ctrl%dkrymax+1) = 0.0_rp

       ! start iterations
       kloc = 0
       kloc = 0
       inner: do while ( (.not.exit_loop) .and. &
            &            (ctrl%it < ctrl%itmax) .and. &
            &            (kloc < ctrl%dkrymax))
          kloc  = kloc  + 1
          ctrl%it = ctrl%it + 1

          ! Generate new basis vector
          call A%apply( bkry(kloc), r )
          call bkry(kloc+1)%default_initialization()
          call bkry(kloc+1)%clone(x)
          call M%apply(r, bkry(kloc+1) )


          if ( env%am_i_fine_task() ) then ! Am I a fine task ?
             ! Orthogonalize
             select case( ctrl%orto )
             case ( mgs )
                call mgsro  ( ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
             case ( icgs )
                call icgsro ( ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
             case default
                ! Write an error message and stop ?      
             end select
             if ( ierrc < 0 ) then
                ! The coarse-grid task should exit 
                ! the inner-do loop. Send signal.
                exit_loop = .true.
                call env%bcast(exit_loop) 
                exit inner ! Exit inner do-loop
             end if

             ! Apply previous given's rotations to kth column of hessenberg matrix
             k_hh = 1
             do j = 1,kloc-1    
                alpha = hh(k_hh,kloc)
                c = cs(1,j)
                s = cs(2,j)
                hh(k_hh,kloc) = c*alpha + s*hh(k_hh+1,kloc)
                hh(k_hh+1,kloc) = c*hh(k_hh+1,kloc) - s*alpha
                k_hh = k_hh +1
             enddo

             ! Compute (and apply) new given's rotation
             call givens(hh(k_hh,kloc), hh(k_hh+1,kloc), c, s)
             cs(1,kloc) = c
             cs(2,kloc) = s

             ! Update preconditioned residual vector 
             ! (expressed in terms of the preconditioned Krylov basis)
             alpha = -s*g(kloc)
             g(kloc) = c*g(kloc)
             g(kloc+1) = alpha

             ! Error norm
             res_norm = abs(alpha)

             if ( ctrl%stopc == res_res .or. ctrl%stopc == res_rhs ) then
                kloc_aux          = kloc 
                g_aux(1:kloc_aux) = g(1:kloc_aux) 
                if ( kloc_aux > 0 ) then
                   ! Compute the solution
                   ! If zero on the diagonal, 
                   ! solve a reduced linear system
                   do while ( kloc_aux > 0 ) ! .and. hh(kloc_aux,kloc_aux) == 0.0_rp  )
                      if(hh(kloc_aux,kloc_aux) /= 0.0_rp) exit
                      kloc_aux = kloc_aux - 1
                   end do

                   if ( kloc_aux <= 0 ) then
                      write (ctrl%luout,*) '** Warning: LGMRES: triangular system in GMRES has null rank'
                      ! The coarse-grid task should exit 
                      ! the outer-do loop. Send signal. 
                      exit_loop = .true.
                      call env%bcast(exit_loop)
                      call env%bcast(exit_loop)
                      exit outer ! Exit outer do loop     
                   end if

#ifdef ENABLE_BLAS       
                   !N    !A  !LDA        !X !INCX
                   call DTRSV ( 'U', 'N', 'N', kloc_aux, hh, ctrl%dkrymax+1, g_aux, 1)
#else
                   ! Solve the system hh*y = g_aux
                   ! Solution stored on g_aux itself
                   do j = kloc_aux,1,-1
                      g_aux(j) = g_aux(j)/hh(j,j)
                      do i = j-1,1,-1
                         g_aux(i) = g_aux(i) - hh(i,j) * g_aux(j)
                      end do
                   end do
#endif

                   ! Now g contains the solution in the krylov basis
                   ! Compute the solution in the global space
                   call z%copy(x)

                   ! z <-z +  g_aux_1 * v_1 + g_aux_2 * v_2 + ... + g_aux_kloc_aux * v_kloc_aux
                   ! call generic_krylov_basis_multiaxpy ( kloc_aux, 1.0_rp, bkry, g_aux, z )
                   do i=1, kloc_aux
                      call z%axpby(g_aux(i),bkry(i),1.0_rp)
                   end do

                   ! r = Az
                   call A%apply(z,r)

                   ! r = b-r
                   ! call generic_pxmy(b,r)
                   call r%axpby(1.0_rp,b,-1.0_rp)

                   res_2_norm = r%nrm2()
                end if
                ctrl%err1 = res_2_norm
             else
                ctrl%err1 = res_norm
             end if
          end if

          if ( ctrl%track_conv_his ) ctrl%err1h(ctrl%it) = ctrl%err1
          exit_loop = (ctrl%err1 <= ctrl%tol1)
          ! Send converged to coarse-grid tasks
          call env%bcast(exit_loop)

          if ( env%am_i_fine_task() ) then
             if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
          end if
       end do inner

       max_kloc = max(kloc+1,max_kloc)

       if ( env%am_i_fine_task() ) then ! Am I a fine task ?
          if ( ierrc == -2 ) then
             write (ctrl%luout,*) '** Warning: LGMRES: ortho failed due to abnormal numbers, no way to proceed'
             ! The coarse-grid task should exit 
             ! the outer-do loop. Send signal. 
             exit_loop = .true.
             call env%bcast(exit_loop)
             exit outer ! Exit outer do-loop
          end if

          if ( ctrl%stopc /= res_res .and. ctrl%stopc /= res_rhs ) then
             if ( kloc > 0 ) then
                ! Compute the solution
                ! If zero on the diagonal, 
                ! solve a reduced linear system
                do while ( kloc > 0 ) !.and. hh(kloc,kloc) == 0.0_rp  )
                   if(hh(kloc,kloc) /= 0.0_rp) exit
                   kloc = kloc - 1
                end do

                if ( kloc <= 0 ) then
                   write (ctrl%luout,*) '** Warning: LGMRES: triangular system in GMRES has null rank'
                   ! The coarse-grid task should exit 
                   ! the outer-do loop. Send signal. 
                   exit_loop = .true.
                   call env%bcast(exit_loop)
                   exit outer ! Exit outer do loop     
                end if
                
#ifdef ENABLE_BLAS       
                !N    !A  !LDA        !X !INCX
                call DTRSV ( 'U', 'N', 'N', kloc, hh, ctrl%dkrymax+1, g, 1)
#else
                ! Solve the system hh*y = g
                ! Solution stored on g itself
                do j = kloc,1,-1
                   g(j) = g(j)/hh(j,j)
                   do i = j-1,1,-1
                      g(i) = g(i) - hh(i,j) * g(j)
                   end do
                end do
#endif

                ! Now g contains the solution in the krylov basis
                ! Compute the solution in the global space
                call z%init(0.0_rp)

                ! z <- z +  g_1 * v_1 + g_2 * v_2 + ... + g_kloc * v_kloc
                do i=1, kloc
                   call z%axpby(g(i),bkry(i),1.0_rp)
                end do

                ! x <- x + z 
                call x%axpby(1.0_rp,z,1.0_rp)
             endif
          else
             ! x <- z
             call x%copy(z)
          end if
       end if

       exit_loop = (ctrl%err1 <= ctrl%tol1)
       ! Send converged to coarse-grid tasks
       call env%bcast(exit_loop)
    end do outer

    ctrl%converged = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast( ctrl%converged )

    call r%free()
    call z%free()
    deallocate(r)
    deallocate(z)

    ! Deallocate Krylov basis
    do i=1, max_kloc
       call bkry(i)%free()
    end do
    deallocate ( bkry )

    ! Deallocate working vectors
    call memfree(hh,__FILE__,__LINE__)
    call memfree(g,__FILE__,__LINE__)
    call memfree(cs,__FILE__,__LINE__)
    if ( ctrl%stopc == res_res .or. ctrl%stopc == res_rhs ) then
       call memfree(g_aux,__FILE__,__LINE__)
    end if

    if ( env%am_i_fine_task() ) then
       if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
  end subroutine abstract_plgmres




!=============================================================================
!
! Abstract Right Preconditioned GMRES
!
! IMPORTANT NOTE: it does not work for element-based data distributions.
! Indeed, if it were properly adapted for element-based data distributions,
! it would be terrible in term of communications (because the krylov
! basis must be composed of partially summed vectors). 
!
! IMPORTANT NOTE 2: The above note is not that severe as long as icgs is used. 
! I have already implemented PRGMRES for element-based data distributions. 
! Up to 6 nearest neighbour comms are required per Arnoldi iteration. This
! is indeed too much 6 can be reduced to 2 by having the krylov basis replicated
! in both states (i.e., part/full summed). 
! The expense comes from more memory consumption and more update vector operations.
! Does it pay off ? (TO CONSIDER)
!
! IMPORTANT NOTE 3: We have been able to develop a communication equivalent
! rgmres-like method by the introduction of the new generic vector routine 
! generic_comm.
!
!=============================================================================
!
! TODO:  1. Review orthogonalization procedures. (Done: MGS or ICGS)
!        2. Implement the whole list of convergence criteria available
!           for PCG. Only two are currently implemented.
!        3. Set err1/err2 accordingly to the work performed in
!           point 2 above.
!        4. Use level-2-BLAS for backward/forward substitution + 
!           ICGS. This last option would imply defining a generic
!           data structure for storing Krylov subspace bases.  (DONE)
!        5. Breakdown detection ? (DONE: Taken from Y. Saad's SPARSKIT)
subroutine abstract_prgmres ( A, M, b, x, ctrl, env)
  !--------------------------------------------------------------------
  ! This routine performs prgmres iterations on Ax=b with preconditioner M. 
  !--------------------------------------------------------------------
#ifdef ENABLE_BLAS       
use blas77_interfaces_names
#endif
  implicit none
  class(operator_t)   , intent(in)    :: A              ! Matrix
  class(operator_t)   , intent(in)    :: M              ! Preconditioner
  class(vector_t)    , intent(inout) :: x              ! Solution
  class(vector_t)    , intent(in)    :: b              ! RHS
  type(solver_control_t)   , intent(inout) :: ctrl
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 
      

  integer(ip)                :: ierrc
  integer(ip)                :: kloc, max_kloc, i, j, k_hh, id
  real(rp)                   :: res_norm, rhs_norm
  real(rp)                   :: alpha, c, s
  real(rp)   , allocatable   :: hh(:,:), g(:), cs(:,:)
  integer                    :: me, np
  logical                    :: exit_loop

  class(vector_t), allocatable  :: r, z       ! Working vector_ts
  class(vector_t), allocatable  :: bkry(:)    ! Krylov basis

    assert(ctrl%stopc==res_nrmgiven_rhs_nrmgiven.or.ctrl%stopc==res_nrmgiven_res_nrmgiven)


    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    ! Clone x in order to allocate working vectors
    allocate(r, mold=b); call r%default_initialization()
    allocate(z, mold=x); call z%default_initialization()
    call r%clone(b)
    call z%clone(x)

    allocate(bkry(ctrl%dkrymax+1), mold=x)

    ! Allocate working vectors
    call memalloc(ctrl%dkrymax+1,ctrl%dkrymax+1,hh,__FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1,g,__FILE__,__LINE__)
    call memalloc(2,ctrl%dkrymax+1,cs,__FILE__,__LINE__)

    call env%info(me, np)

    ! Evaluate ||b||_2 if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x,r)

    ! r = b-r
    call r%axpby(1.0_rp,b,-1.0_rp)

    ! part_summed to full_summed
    call r%comm() 

    ! res_norm = ||r||_2
    res_norm = r%nrm2()

    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
        ctrl%tol1  = ctrl%rtol * rhs_norm + ctrl%atol
        ctrl%err1 = res_norm 
    else if ( ctrl%stopc == res_nrmgiven_res_nrmgiven ) then
        ctrl%tol1  = ctrl%rtol * res_norm + ctrl%atol
        ctrl%err1 = res_norm
    end if
    exit_loop = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(exit_loop)

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_header(ctrl)
    end if

    max_kloc = 0

    ctrl%it = 0
  outer: do while ( (.not.exit_loop) .and. &
       &            (ctrl%it < ctrl%itmax))

        ! Compute residual from scratch (only if ctrl%it /= 0)
        if ( ctrl%it /= 0 ) then 
            ! r = Ax
            call A%apply(x,r)
    
            ! r = b-r
            call r%axpby(1.0_rp,b,-1.0_rp)
            
            ! part_summed to full_summed
            call r%comm()

            ! res_norm = ||r||_2
            res_norm = r%nrm2()
        end if

        ! Normalize residual direction (i.e., v_1 = r/||r||_2)
        call bkry(1)%default_initialization()
        call bkry(1)%clone(x)
        if ( env%am_i_fine_task() ) then 
          if (res_norm /= 0.0_rp) then
            call bkry(1)%scal(1.0_rp/res_norm, r)
          end if
        end if
         ! residual in the krylov basis
         g(1) = res_norm
         g(2:ctrl%dkrymax+1) = 0.0_rp

         ! start iterations
        kloc = 0
         inner: do while ( (.not.exit_loop) .and. &
              &            (ctrl%it < ctrl%itmax) .and. &
              &            (kloc < ctrl%dkrymax))
            kloc  = kloc  + 1
            ctrl%it = ctrl%it + 1

            ! Generate new basis vector
            call M%apply( bkry(kloc), z )
            call bkry(kloc+1)%default_initialization()
            call bkry(kloc+1)%clone(x)
            call A%apply(z, bkry(kloc+1))

            ! part_summed to full_summed
            call bkry(kloc+1)%comm()

            if ( env%am_i_fine_task() ) then
                ! Orthogonalize
                select case( ctrl%orto )
                    case ( mgs )
                        call mgsro  (ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
                    case ( icgs )
                        call icgsro (ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
                    case default
                        check(.false.)
                end select
           
                if ( ierrc < 0 ) then
                    ! The coarse-grid task should exit 
                    ! the inner-do loop. Send signal.
                    exit_loop = .true.
                    call env%bcast(exit_loop) 
                    exit inner ! Exit inner do-loop
                end if
           
                ! Apply previous given's rotations to kth column of hessenberg matrix
                k_hh = 1
                do j = 1,kloc-1    
                    alpha = hh(k_hh,kloc)
                    c = cs(1,j)
                    s = cs(2,j)
                    hh(k_hh,kloc) = c*alpha + s*hh(k_hh+1,kloc)
                    hh(k_hh+1,kloc) = c*hh(k_hh+1,kloc) - s*alpha
                    k_hh = k_hh +1
                enddo

                ! Compute (and apply) new given's rotation
                call givens(hh(k_hh,kloc), hh(k_hh+1,kloc), c, s)
                cs(1,kloc) = c
                cs(2,kloc) = s

                ! Update residual vector 
                ! (expressed in terms of the Krylov basis)
                alpha = -s*g(kloc)
                g(kloc) = c*g(kloc)
                g(kloc+1) = alpha

                ! Error norm
                res_norm = abs(alpha)
                ctrl%err1  = res_norm
                if ( ctrl%track_conv_his ) ctrl%err1h(ctrl%it) = ctrl%err1
            end if
            exit_loop = (ctrl%err1 <= ctrl%tol1) 
            ! Send converged to coarse-grid tasks
            call env%bcast(exit_loop)
        
            if ( env%am_i_fine_task() ) then
                if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
            end if
        end do inner

        max_kloc = max(kloc+1, max_kloc)

        if ( kloc > 0 ) then
            if ( env%am_i_fine_task() ) then
                if ( ierrc == -2 ) then
                    write (ctrl%luout,*) '** Warning: RGMRES: ortho failed due to abnormal numbers, no way to proceed'
                    ! The coarse-grid task should exit 
                    ! the outer-do loop. Send signal. 
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit outer ! Exit main do-loop
                end if

                ! Compute the solution
                ! If zero on the diagonal, 
                ! solve a reduced linear system
                do while ( kloc > 0 ) ! .and. hh(kloc,kloc) == 0.0_rp  )
                    if(hh(kloc,kloc) /= 0.0_rp) exit
                    kloc = kloc - 1
                end do

                if ( kloc <= 0 ) then
                    write (ctrl%luout,*) '** Warning: RGMRES: triangular system in GMRES has null rank'
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit outer ! Exit main do loop     
                end if

#ifdef ENABLE_BLAS       
                !N    !A  !LDA        !X !INCX
                call DTRSV ( 'U', 'N', 'N', kloc, hh, ctrl%dkrymax+1, g, 1)
#else
                ! Solve the system hh*y = g
                ! Solution stored on g itself
                do j = kloc,1,-1
                    g(j) = g(j)/hh(j,j)
                    do i = j-1,1,-1
                        g(i) = g(i) - hh(i,j) * g(j)
                    end do
                end do
#endif       
            end if

            ! Now g contains the solution in the krylov basis
            ! Compute the solution in the real space
            call r%init(0.0_rp)

            ! r <- g_1 * v_1 + g_2 * v_2 + ... + g_kloc * v_kloc
            do i=1, kloc
                call r%axpby(g(i),bkry(i),1.0_rp)
            end do

            ! Solve Mz = r
            call M%apply(r,z)

            ! x <- x + z
            call x%axpby(1.0_rp,z,1.0_rp)

            exit_loop = (ctrl%err1 <= ctrl%tol1)
            ! Send converged to coarse-grid tasks
            call env%bcast(exit_loop)
        end if
    end do outer

    ctrl%converged = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(ctrl%converged)

    ! Deallocate working vectors
    call memfree(hh,__FILE__,__LINE__)
    call memfree(g,__FILE__,__LINE__)
    call memfree(cs,__FILE__,__LINE__)
    
    call r%free()
    call z%free()
    deallocate(r)
    deallocate(z)

    ! Deallocate Krylov basis
    do i=1, max_kloc
       call bkry(i)%free()
    end do
    deallocate ( bkry )

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
end subroutine abstract_prgmres


!=============================================================================
!
! Abstract Flexible GMRES
!
subroutine abstract_pfgmres ( A, M, b, x, ctrl, env)
  !--------------------------------------------------------------------
  ! This routine performs pfgmres iterations on Ax=b with (right)-preconditioner M. 
  !--------------------------------------------------------------------
#ifdef ENABLE_BLAS       
use blas77_interfaces_names
#endif
  implicit none
  class(operator_t)   , intent(in)    :: A              ! Matrix
  class(operator_t)   , intent(in)    :: M              ! Preconditioner
  class(vector_t)    , intent(inout) :: x              ! Solution
  class(vector_t)    , intent(in)    :: b              ! RHS
  type(solver_control_t)  , intent(inout) :: ctrl
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 


  integer(ip)                :: ierrc
  integer(ip)                :: kloc, max_kloc, i, j, k_hh, id
  real(rp)                   :: res_norm, rhs_norm
  real(rp)                   :: alpha, c, s
  real(rp)   , allocatable   :: hh(:,:), g(:), cs(:,:)
  integer                    :: me, np
  logical                    :: exit_loop

  class(vector_t), allocatable :: r, z    ! Working vector_ts
  class(vector_t), allocatable :: bkry(:)    ! Krylov basis
  class(vector_t), allocatable :: bkryz(:)   ! (Right-)Preconditioned Krylov basis

    assert(ctrl%stopc==res_nrmgiven_rhs_nrmgiven.or.ctrl%stopc==res_nrmgiven_res_nrmgiven)

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    ! Clone x in order to allocate working vectors
    allocate(r, mold=b); call r%default_initialization()
    allocate(z, mold=x); call z%default_initialization()
    call r%clone(b)
    call z%clone(x)

    allocate(bkry(ctrl%dkrymax+1), mold=x)
    allocate(bkryz(ctrl%dkrymax), mold=x)

    ! Allocate working vectors
    call memalloc(ctrl%dkrymax+1,ctrl%dkrymax+1,hh,__FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1,g,__FILE__,__LINE__)
    call memalloc(2,ctrl%dkrymax+1,cs,__FILE__,__LINE__)

    call env%info(me, np)

    ! Evaluate ||b||_2 if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x, r)

    ! r = b-r
    call r%axpby(1.0_rp,b,-1.0_rp)

    ! part_summed to full_summed
    call r%comm()

    ! res_norm = ||r||_2
    res_norm = r%nrm2()

    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
        ctrl%tol1  = ctrl%rtol * rhs_norm + ctrl%atol
        ctrl%err1 = res_norm 
    else if ( ctrl%stopc == res_nrmgiven_res_nrmgiven ) then
        ctrl%tol1  = ctrl%rtol * res_norm + ctrl%atol
        ctrl%err1 = res_norm
    end if
    exit_loop = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(exit_loop)
    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_header(ctrl)
    end if

    max_kloc = 0
    ctrl%it = 0
    outer: do while ( (.not.exit_loop ) .and. &
        &            (ctrl%it < ctrl%itmax))
        ! Compute residual from scratch (only if ctrl%it /= 0)
        if ( ctrl%it /= 0 ) then 
            ! r = Ax
            call A%apply(x, r)

            ! r = b-r
            call r%axpby(1.0_rp,b,-1.0_rp)

            ! part_summed to full_summed
            call r%comm()

            ! res_norm = ||r||_2
            res_norm = r%nrm2()
        end if

        ! Normalize residual direction (i.e., v_1 = r/||r||_2)
        call bkry(1)%default_initialization()
        call bkry(1)%clone(x)
        if (res_norm /= 0.0_rp) call bkry(1)%scal(1.0_rp/res_norm, r)

        ! residual in the krylov basis
        g(1) = res_norm
        g(2:ctrl%dkrymax+1) = 0.0_rp

        ! start iterations
        kloc = 0
        inner: do while ( (.not.exit_loop) .and. &
            &            (ctrl%it < ctrl%itmax) .and. &
            &            (kloc < ctrl%dkrymax))
            kloc  = kloc  + 1
            ctrl%it = ctrl%it + 1

            ! Generate new basis vector
            call bkryz(kloc)%default_initialization()
            call bkryz(kloc)%clone(x)
            call M%apply(bkry(kloc),bkryz(kloc))
            call bkry(kloc+1)%default_initialization()
            call bkry(kloc+1)%clone(x)
            call A%apply(bkryz(kloc),bkry(kloc+1))

            ! part_summed to full_summed
            call bkry(kloc+1)%comm()
        
            if ( env%am_i_fine_task() ) then
                ! Orthogonalize
                select case( ctrl%orto )
                    case ( mgs )
                        call mgsro  (ctrl%luout,kloc+1, bkry, hh(1,kloc), ierrc )
                    case ( icgs )
                        call icgsro (ctrl%luout,kloc+1, bkry, hh(1,kloc), ierrc )
                    case default
                        ! Write an error message and stop ?      
                        check(.false.)
                end select

                if ( ierrc < 0 ) then
                    ! The coarse-grid task should exit 
                    ! the inner-do loop. Send signal.
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit inner ! Exit inner do-loop
                end if

                ! Apply previous given's rotations to kth column of hessenberg matrix
                k_hh = 1
                do j = 1,kloc-1    
                    alpha = hh(k_hh,kloc)
                    c = cs(1,j)
                    s = cs(2,j)
                    hh(k_hh,kloc) = c*alpha + s*hh(k_hh+1,kloc)
                    hh(k_hh+1,kloc) = c*hh(k_hh+1,kloc) - s*alpha
                    k_hh = k_hh +1
                enddo

                ! Compute (and apply) new given's rotation
                call givens(hh(k_hh,kloc), hh(k_hh+1,kloc), c, s)
                cs(1,kloc) = c
                cs(2,kloc) = s

                ! Update residual vector 
                ! (expressed in terms of the Krylov basis)
                alpha = -s*g(kloc)
                g(kloc) = c*g(kloc)
                g(kloc+1) = alpha

                ! Error norm
                res_norm = abs(alpha)
                ctrl%err1  = res_norm
                if (ctrl%track_conv_his) ctrl%err1h(ctrl%it) = ctrl%err1
            end if
            exit_loop = (ctrl%err1 <= ctrl%tol1)
            ! Send converged to coarse-grid tasks
            call env%bcast(exit_loop)

            if ( env%am_i_fine_task() ) then
                if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
            end if
        end do inner

        max_kloc = max(kloc+1, max_kloc)

        if ( kloc > 0 ) then
            if ( env%am_i_fine_task() ) then
                if ( ierrc == -2 ) then
                    write (ctrl%luout,*) '** Warning: FGMRES: ortho failed due to abnormal numbers, no way to proceed'
                    ! The coarse-grid task should exit 
                    ! the outer-do loop. Send signal. 
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit outer ! Exit main do-loop
                end if

                ! Compute the solution
                ! If zero on the diagonal, 
                ! solve a reduced linear system
                do while ( kloc > 0 ) ! .and. hh(kloc,kloc) == 0.0_rp  )
                    if(hh(kloc,kloc) /= 0.0_rp) exit
                    kloc = kloc - 1
                end do

                if ( kloc <= 0 ) then
                    write (ctrl%luout,*) '** Warning: FGMRES: triangular system in FGMRES has null rank'
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit outer    
                end if
#ifdef ENABLE_BLAS       
                !N    !A  !LDA        !X !INCX
                call DTRSV ( 'U', 'N', 'N', kloc, hh, ctrl%dkrymax+1, g, 1)
#else
                ! Solve the system hh*y = g
                ! Solution stored on g itself
                do j = kloc,1,-1
                    g(j) = g(j)/hh(j,j)
                    do i = j-1,1,-1
                        g(i) = g(i) - hh(i,j) * g(j)
                    end do
                end do
#endif       
            end if

            ! Now g contains the solution in the krylov basis
            ! Compute the solution in the real space
            call z%init(0.0_rp)

            ! z <- g_1 * M_1^-1v_1 + g_2 *  M_2^-1v_2  + ... + g_kloc *  M_kloc^-1v_kloc
            do j=1, kloc
                call z%axpby(g(j),bkryz(j),1.0_rp)
            enddo

            ! x <- x + z
            call x%axpby(1.0_rp,z,1.0_rp)
        
            exit_loop = (ctrl%err1 <= ctrl%tol1)

            ! Send converged to coarse-grid tasks
            call env%bcast(exit_loop)

        end if
    end do outer

    ctrl%converged = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(ctrl%converged)

    ! Deallocate working vectors
    call memfree(hh,__FILE__,__LINE__)
    call memfree(g,__FILE__,__LINE__)
    call memfree(cs,__FILE__,__LINE__)

    call r%free()
    call z%free()

    ! Deallocate Krylov basis
    do i=1, max_kloc
       call bkry(i)%free()
       call bkryz(i)%free()
    end do
    call bkry(ctrl%dkrymax+1)%free()
    deallocate ( bkry )
    deallocate ( bkryz )

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()

end subroutine abstract_pfgmres


!=============================================================================
!
! Abstract preconditioned RICHARDSON
!
subroutine abstract_prichard (A, M, b, x, ctrl, env )
  !-----------------------------------------------------------------------
  !
  ! This routine solves M^-1Ax = M^-1b where M is a given preconditioner
  ! using preconditioned Richardson fixed-point iterations, with relaxation
  ! parameter given by relax
  !
  !-----------------------------------------------------------------------
  implicit none
  class(operator_t)   , intent(in)    :: A              ! Matrix
  class(operator_t)   , intent(in)    :: M              ! Preconditioner
  class(vector_t)    , intent(inout) :: x              ! Solution
  class(vector_t)    , intent(in)    :: b              ! RHS
  type(solver_control_t)  , intent(inout) :: ctrl
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 


  integer                          :: me, np
  class(vector_t), allocatable :: r, z      ! Working vector_ts
  real(rp)                         :: res_norm, rhs_norm

  assert ( ctrl%stopc == res_res .or. ctrl%stopc == res_rhs ) 

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    allocate(r, mold=x); call r%default_initialization()
    allocate(z, mold=x); call z%default_initialization()

    call env%info(me, np)
    call r%clone(x)
    call z%clone(x)

    ! Evaluate ||b||_2 if required
    if ( ctrl%stopc == res_rhs ) then
        rhs_norm = b%nrm2()
    endif

    ctrl%converged = .false.
    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_header(ctrl)
    end if

    ctrl%it = 0
    loop_prichard: do while( (.not.ctrl%converged) .and. (ctrl%it < ctrl%itmax))
  
        ! r = Ax
        call A%apply(x, r)

        ! r = b-r
        call r%axpby(1.0_rp,b,-1.0_rp)

        ! Evaluate ||r||_L2
        res_norm = r%nrm2()

        ! Set upper bound (only in 1st iteration)
        if ( ctrl%it == 1 ) then
            if ( ctrl%stopc == res_rhs ) then
                ctrl%tol1  = ctrl%rtol * rhs_norm + ctrl%atol 
            else if ( ctrl%stopc == res_res ) then
                ctrl%tol1  = ctrl%rtol * res_norm + ctrl%atol
            end if
        end if
        ctrl%err1 = res_norm
        if (ctrl%it > 0 .and. ctrl%track_conv_his) ctrl%err1h(ctrl%it) = ctrl%err1
        ctrl%converged = (ctrl%err1 <= ctrl%tol1)

        ! Send converged to coarse-grid tasks
        call env%bcast(ctrl%converged)

        if ( env%am_i_fine_task() ) then
            if ((ctrl%it > 0).and.(me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
        end if

        ! z = inv(M) r
        call M%apply(r, z)

        ! x <- x + relax * z	
        call x%axpby(ctrl%relax,z,1.0_rp)

        ctrl%it = ctrl%it + 1
    end do loop_prichard

    call r%free()
    call z%free()

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if


  call A%CleanTemp()
  call M%CleanTemp()
  call b%CleanTemp()

end subroutine abstract_prichard


!====================================================
! Abstract Left Preconditioned FOM
!====================================================
subroutine abstract_plfom ( A, M, b, x, ctrl, env )
  !--------------------------------------------------------------------
  ! This routine performs plfom iterations on Ax=b with preconditioner M. 
  !--------------------------------------------------------------------
#ifdef ENABLE_LAPACK
use lapack77_interfaces_names
#endif

  implicit none
  class(operator_t), intent(in)     :: A      ! Matrix
  class(operator_t), intent(in)     :: M      ! Preconditioner
  class(vector_t),  intent(inout)  :: x      ! Solution
  class(vector_t),  intent(in)     :: b      ! RHS
  type(solver_control_t), intent(inout) :: ctrl
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 


  integer(ip)                    :: ierrc
  integer(ip)                    :: kloc, i, j, k_hh, id
  real(rp)                       :: res_norm, res_2_norm, rhs_norm, res_norm_initial
  real(rp)   , allocatable       :: hh(:,:), lu(:,:), g(:)
  integer(ip), allocatable       :: ipiv(:)
  real(rp)                       :: err1_it, ub1
  integer                        :: me, np, info
  logical                        :: exit_loop

  class(vector_t), allocatable :: r, z        ! Working vector_ts
  class(vector_t), allocatable :: bkry(:)     ! Krylov basis

  write(*,*) 'XXX', ctrl%stopc, res_res, res_rhs
 
    assert(ctrl%stopc==res_res.or.ctrl%stopc==res_rhs)

    

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    ! Clone x in order to allocate working vectors
    allocate(r, mold=b); call r%default_initialization()
    allocate(z, mold=x); call z%default_initialization()
    call r%clone(b)
    call z%clone(x)

    allocate(bkry(ctrl%dkrymax+1), mold=x)

    ! Allocate working vectors
    call memalloc(ctrl%dkrymax+1, ctrl%dkrymax+1,hh, __FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1, ctrl%dkrymax+1,lu, __FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1,                 g, __FILE__,__LINE__)
    call memalloc(ctrl%dkrymax+1,              ipiv, __FILE__,__LINE__)

    call env%info(me, np)

    ! Evaluate ||b||_2 if required
    if ( ctrl%stopc == res_rhs ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x, r)

    ! r = b-r
    call r%axpby(1.0_rp,b,-1.0_rp)

    res_2_norm = r%nrm2()

    ! z=inv(M)r
    call M%apply(r,z)

    ! Evaluate ||z||_2
    res_norm = z%nrm2()
    res_norm_initial = res_norm

    if ( ctrl%stopc == res_res ) then
        ctrl%tol1 = ctrl%rtol * res_2_norm + ctrl%atol
        ctrl%err1 = res_2_norm
    else if ( ctrl%stopc == res_rhs ) then
        ctrl%tol1 = ctrl%rtol * rhs_norm + ctrl%atol
        ctrl%err1 = res_2_norm
    end if

    exit_loop = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(exit_loop)

    if ( env%am_i_fine_task() ) then
        if ( (ctrl%trace > 0) .and. (me == 0) ) call solver_control_log_header(ctrl)
    end if

    ctrl%it = 0
    outer: do while ( (.not.exit_loop).and.(ctrl%it<ctrl%itmax) )
        hh = 0.0_rp
        ! Compute preconditioned residual from scratch (only if ctrl%it/=0)
        if ( ctrl%it /= 0 ) then 
            ! r = Ax
            call A%apply(x, r)

            ! r = b-r
            call r%axpby(1.0_rp,b,-1.0_rp)

            ! z=inv(M)r
            call M%apply(r, z)

            ! Evaluate ||z||_2
            res_norm = z%nrm2()
        end if

        ! Normalize preconditioned residual direction (i.e., v_1 = z/||z||_2)
        call bkry(1)%default_initialization()
        call bkry(1)%clone(x)
        if (res_norm/=0.0_rp) call bkry(1)%scal(1.0_rp/res_norm,z)

        ! start iterations
        kloc = 0
        inner: do while ( (.not.exit_loop) .and. &
            &            (ctrl%it < ctrl%itmax) .and. &
            &            (kloc < ctrl%dkrymax))
            kloc  = kloc  + 1
            ctrl%it = ctrl%it + 1

            ! Generate new basis vector
            call A%apply(bkry(kloc), r)
            call bkry(kloc+1)%default_initialization()
            call bkry(kloc+1)%clone(x)
            call M%apply(r,bkry(kloc+1))

            if ( env%am_i_fine_task() ) then ! Am I a fine task ?
                ! Orthogonalize
                select case( ctrl%orto )
                    case ( mgs )
                        call mgsro  ( ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
                    case ( icgs )
                        call icgsro ( ctrl%luout, kloc+1, bkry, hh(1,kloc), ierrc )
                    case default
                        check(.false.)
                        ! Write an error message and stop ?      
                end select

                if ( ierrc < 0 ) then
                    ! The coarse-grid task should exit 
                    ! the inner-do loop. Send signal.
                    exit_loop = .true.
                    call env%bcast(exit_loop)
                    exit inner ! Exit inner do-loop
                end if
    
                ! init right-hand-size to \beta*e_1, with \beta=||z_0||_2
                g(1)            = res_norm_initial
                g(2:kloc)       = 0.0_rp
    
                if ( kloc > 0 ) then
                    ! Compute the solution
    
                    lu(1:kloc,1:kloc) = hh(1:kloc,1:kloc)
                    ! write(*,*) 'XXX', lu(1:kloc,1:kloc)
#ifdef ENABLE_LAPACK
                    call dgetrf( kloc, kloc, lu, ctrl%dkrymax+1, ipiv, info )
                    if ( info /= 0 ) then
                        write (ctrl%luout,*) '** Warning: LFOM: dgetrf returned info /= 0'
                        exit_loop = .true.
                        call env%bcast(exit_loop)
                        call env%bcast(exit_loop)
                        exit outer ! Exit main do loop 
                    end if
    
                    call dgetrs( 'N' , kloc, 1, lu, ctrl%dkrymax+1, ipiv, g, ctrl%dkrymax+1, info )
                    if ( info /= 0 ) then
                        write (ctrl%luout,*) '** Warning: LFOM: dgetrs returned info /= 0'
                        exit_loop = .true.
                        call env%bcast(exit_loop)
                        call env%bcast(exit_loop)
                        exit outer ! Exit main do loop 
                    end if
#else
                    write (0,*) 'plfom ERROR :: dgetrf and dgetrs not available'
                    check(1==0)
#endif

                    ! Now g contains the solution in the krylov basis
                    ! Compute the solution in the global space
                    call z%copy(x)

                    ! z <-z +  g_1 * v_1 + g_2 * v_2 + ... + g_kloc * v_kloc
                    do j=1, kloc
                        call z%axpby(g(j),bkry(j),1.0_rp)
                    enddo
    
                    ! r = A(z+x)
                    call A%apply(z, r)
    
                    ! r = b-r
                    call r%axpby(1.0_rp,b,-1.0_rp)
    
                    res_2_norm = r%nrm2()
                end if
                ctrl%err1 = res_2_norm
            end if

            if (ctrl%track_conv_his) ctrl%err1h(ctrl%it) = ctrl%err1
            exit_loop = (ctrl%err1 <= ctrl%tol1)
            ! Send converged to coarse-grid tasks
            call env%bcast(exit_loop)
    
            if ( env%am_i_fine_task() ) then
                if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
            end if

        end do inner

        if ( env%am_i_fine_task() ) then ! Am I a fine task ?
            if ( ierrc == -2 ) then
                write (ctrl%luout,*) '** Warning: LFOM: ortho failed due to abnormal numbers, no way to proceed'
                ! The coarse-grid task should exit 
                ! the outer-do loop. Send signal. 
                exit_loop = .true.
                call env%bcast(exit_loop)
                exit outer ! Exit outer do-loop
            end if
        end if

        ! x <- z
        call x%copy(z)

        exit_loop = (ctrl%err1 <= ctrl%tol1)
        ! Send converged to coarse-grid tasks
        call env%bcast(exit_loop)

    end do outer

    ctrl%converged = (ctrl%err1 <= ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call env%bcast(ctrl%converged)

    call r%free()
    call z%free()

    do j=1,ctrl%dkrymax+1
        call bkry(j)%free()
    enddo
    deallocate(bkry)

    call memfree(hh,__FILE__,__LINE__)
    call memfree(lu,__FILE__,__LINE__)
    call memfree(g,__FILE__,__LINE__)
    call memfree(ipiv,__FILE__,__LINE__)

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()

end subroutine abstract_plfom



!=============================================================================
! Generic Preconditioned Minimal residual method
! Generic code extracted from F90 PMINRES implementation available at
! http://www.stanford.edu/group/SOL/software/minres/
!=============================================================================
!-------------------------------------------------------------------
!
! MINRES  is designed to solve the system of linear equations
!
!    Ax = b
!
! or the least-squares problem
!
!    min ||Ax - b||_2,
!
! where A is an n by n symmetric matrix and b is a given vector.
! The matrix A may be indefinite and/or singular.
!
! 1. If A is known to be positive definite, the Conjugate Gradient
! Method might be preferred, since it requires the same number
! of iterations as MINRES but less work per iteration.
!
! 2. If A is indefinite but Ax = b is known to have a solution
! (e.g. if A is nonsingular), SYMMLQ might be preferred,
! since it requires the same number of iterations as MINRES
! but slightly less work per iteration.
!
! The matrix A is intended to be large and sparse.  It is accessed
! by means of a subroutine call of the form
! SYMMLQ development:
!
!    call Aprod ( n, x, y )
!
! which must return the product y = Ax for any given vector x.
!
!
! More generally, MINRES is designed to solve the system
!
!    (A - shift*I) x = b
! or
!    min ||(A - shift*I) x - b||_2,
!
! where  shift  is a specified scalar value.  Again, the matrix
! (A - shift*I) may be indefinite and/or singular.
! The work per iteration is very slightly less if  shift = 0.
!
! Note: If  shift  is an approximate eigenvalue of  A
! and  b  is an approximate eigenvector,  x  might prove to be
! a better approximate eigenvector, as in the methods of
! inverse iteration and/or Rayleigh-quotient iteration.
! However, we're not yet sure on that -- it may be better to use SYMMLQ.
!
! A further option is that of preconditioning, which may reduce
! the number of iterations required.  If M = C C' is a positive
! definite matrix that is known to approximate  (A - shift*I)
! in some sense, and if systems of the form  My = x  can be
! solved efficiently, the parameters precon and Msolve may be
! used (see below).  When  precon = .true., MINRES will
! implicitly solve the system of equations
!
!    P (A - shift*I) P' xbar  =  P b,
!
! i.e.             Abar xbar  =  bbar
! where                    P  =  C**(-1),
!                       Abar  =  P (A - shift*I) P',
!                       bbar  =  P b,
!
! and return the solution       x  =  P' xbar.
! The associated residual is rbar  =  bbar - Abar xbar
!                                  =  P (b - (A - shift*I)x)
!                                  =  P r.
!
! In the discussion below, eps refers to the machine precision.
!
! Parameters
! ----------
!
! n       input      The dimension of the matrix A.
! b(n)    input      The rhs vector b.
! x(n)    output     Returns the computed solution x.
!
! Aprod   external   A subroutine defining the matrix A.
!                       call Aprod ( n, x, y )
!                    must return the product y = Ax
!                    without altering the vector x.
!
! Msolve  external   An optional subroutine defining a
!                    preconditioning matrix M, which should
!                    approximate (A - shift*I) in some sense.
!                    M must be positive definite.
!
!                       call Msolve( n, x, y )
!
!                    must solve the linear system My = x
!                    without altering the vector x.
!
!                    In general, M should be chosen so that Abar has
!                    clustered eigenvalues.  For example,
!                    if A is positive definite, Abar would ideally
!                    be close to a multiple of I.
!                    If A or A - shift*I is indefinite, Abar might
!                    be close to a multiple of diag( I  -I ).
!
! checkA  input      If checkA = .true., an extra call of Aprod will
!                    be used to check if A is symmetric.  Also,
!                    if precon = .true., an extra call of Msolve
!                    will be used to check if M is symmetric.
!
! precon  input      If precon = .true., preconditioning will
!                    be invoked.  Otherwise, subroutine Msolve
!                    will not be referenced; in this case the
!                    actual parameter corresponding to Msolve may
!                    be the same as that corresponding to Aprod.
!
! shift   input      Should be zero if the system Ax = b is to be
!                    solved.  Otherwise, it could be an
!                    approximation to an eigenvalue of A, such as
!                    the Rayleigh quotient b'Ab / (b'b)
!                    corresponding to the vector b.
!                    If b is sufficiently like an eigenvector
!                    corresponding to an eigenvalue near shift,
!                    then the computed x may have very large
!                    components.  When normalized, x may be
!                    closer to an eigenvector than b.
!
! nout    input      A file number.
!                    If nout > 0, a summary of the iterations
!                    will be printed on unit nout.
!
! itnlim  input      An upper limit on the number of iterations.
!
! rtol    input      A user-specified tolerance.  MINRES terminates
!                    if it appears that norm(rbar) is smaller than
!                       rtol * norm(Abar) * norm(xbar),
!                    where rbar is the transformed residual vector,
!                       rbar = bbar - Abar xbar.
!
!                    If shift = 0 and precon = .false., MINRES
!                    terminates if norm(b - A*x) is smaller than
!                       rtol * norm(A) * norm(x).
!
! istop   output     An integer giving the reason for termination...
!
!          -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                    second Lanczos vector is zero.  This means the
!                    rhs is very special.
!                    If there is no preconditioner, b is an
!                    eigenvector of A.
!                    Otherwise (if precon is true), let My = b.
!                    If shift is zero, y is a solution of the
!                    generalized eigenvalue problem Ay = lambda My,
!                    with lambda = alpha1 from the Lanczos vectors.
!
!                    In general, (A - shift*I)x = b
!                    has the solution         x = (1/alpha1) y
!                    where My = b.
!
!           0        b = 0, so the exact solution is x = 0.
!                    No iterations were performed.
!
!           1        Norm(rbar) appears to be less than
!                    the value  rtol * norm(Abar) * norm(xbar).
!                    The solution in  x  should be acceptable.
!
!           2        Norm(rbar) appears to be less than
!                    the value  eps * norm(Abar) * norm(xbar).
!                    This means that the residual is as small as
!                    seems reasonable on this machine.
!
!           3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                    which should indicate that x has essentially
!                    converged to an eigenvector of A
!                    corresponding to the eigenvalue shift.
!
!           4        Acond (see below) has exceeded 0.1/eps, so
!                    the matrix Abar must be very ill-conditioned.
!                    x may not contain an acceptable solution.
!
!           5        The iteration limit was reached before any of
!                    the previous criteria were satisfied.
!
!           6        The matrix defined by Aprod does not appear
!                    to be symmetric.
!                    For certain vectors y = Av and r = Ay, the
!                    products y'y and r'v differ significantly.
!
!           7        The matrix defined by Msolve does not appear
!                    to be symmetric.
!                    For vectors satisfying My = v and Mr = y, the
!                    products y'y and r'v differ significantly.
!
!           8        An inner product of the form  x' M**(-1) x
!                    was not positive, so the preconditioning matrix
!                    M does not appear to be positive definite.
!
!                    If istop >= 5, the final x may not be an
!                    acceptable solution.
!
! itn     output     The number of iterations performed.
!
! Anorm   output     An estimate of the norm of the matrix operator
!                    Abar = P (A - shift*I) P',   where P = C**(-1).
!
! Acond   output     An estimate of the condition of Abar above.
!                    This will usually be a substantial
!                    under-estimate of the true condition.
!
! rnorm   output     An estimate of the norm of the final
!                    transformed residual vector,
!                       P (b  -  (A - shift*I) x).
!
! ynorm   output     An estimate of the norm of xbar.
!                    This is sqrt( x'Mx ).  If precon is false,
!                    ynorm is an estimate of norm(x).
!-------------------------------------------------------------------
! MINRES is an implementation of the algorithm described in
! the following reference:
!
! C. C. Paige and M. A. Saunders (1975),
! Solution of sparse indefinite systems of linear equations,
! SIAM J. Numer. Anal. 12(4), pp. 617-629.
!-------------------------------------------------------------------
!
!
! MINRES development:
!    1972: First version, similar to original SYMMLQ.
!          Later lost @#%*!
!    Oct 1995: Tried to reconstruct MINRES from
!              1995 version of SYMMLQ.
! 30 May 1999: Need to make it more like LSQR.
!              In middle of major overhaul.
! 19 Jul 2003: Next attempt to reconstruct MINRES.
!              Seems to need two vectors more than SYMMLQ.  (w1, w2)
!              Lanczos is now at the top of the loop,
!              so the operator Aprod is called in just one place
!              (not counting the initial check for symmetry).
! 22 Jul 2003: Success at last.  Preconditioning also works.
!              minres.f added to http://www.stanford.edu/group/SOL/.
!
! 16 Oct 2007: Added a stopping rule for singular systems,
!              as derived in Sou-Cheng Choi's PhD thesis.
!              Note that ||Ar|| small => r is a null vector for A.
!              Subroutine minrestest2 in minresTestModule.f90
!              tests this option.  (NB: Not yet working.)
!-------------------------------------------------------------------
subroutine abstract_pminres(A, M, b, x, ctrl, env)
  !-----------------------------------------------------------------------------
  ! This routine performs pcg iterations on Ax=b with preconditioner M. 
  !-----------------------------------------------------------------------------
  implicit none

  ! Mandatory parameters
  class(operator_t), intent(in)    :: A        ! Matrix
  class(operator_t), intent(in)    :: M        ! Preconditioner
  class(vector_t),  intent(in)    :: b        ! RHS
  class(vector_t), intent(inout)  :: x        ! Approximate solution
  type(solver_control_t), intent(inout) :: ctrl
  class(abstract_environment_t), intent(in) :: env      ! Serial/parallel environment 


  !     Local arrays and variables
  class(vector_t), allocatable :: r1, r2, v1, v2, w, w1, w2, y
  real(rp)  :: alfa  , beta  , beta1 , cs    ,          &
       dbar  , delta , denom , diag  ,          &
       eps   , epsa  , epsln , epsr  , epsx  ,  &
       gamma , gbar  , gmax  , gmin  ,          &
       oldb  , oldeps, qrnorm, phi   , phibar,  &
       rhs1  , rhs2  , rnorml, rootl ,          &
       s     , sn    , t     , tnorm2, ynorm2, z

  integer(ip) :: i, istop, itn
  logical :: debug, prnt
  logical     :: beta1_lt_zero, beta1_eq_zero, beta_lt_zero, istop_neq_zero

  real(rp) ::   Anorm, Acond, rnorm, ynorm
  integer  :: me, np


  ! Local constants
  real(rp),         parameter :: zero =  0.0_rp,  one = 1.0_rp
  real(rp),         parameter :: ten  = 10.0_rp
  character(len=*), parameter :: msg(-1:8) =                  &
       (/ 'beta2 = 0.  If M = I, b and x are eigenvectors of A', & ! -1
       'beta1 = 0.  The exact solution is  x = 0           ', & !  0
       'Requested accuracy achieved, as determined by rtol ', & !  1
       'Reasonable accuracy achieved, given eps            ', & !  2
       'x has converged to an eigenvector                  ', & !  3
       'Acond has exceeded 0.1/eps                         ', & !  4
       'The iteration limit was reached                    ', & !  5
       'Aprod  does not define a symmetric matrix          ', & !  6
       'Msolve does not define a symmetric matrix          ', & !  7
       'Msolve does not define a pos-def preconditioner    ' /) !  8

  character(len=*), parameter    :: fmt11='(a,2x,es16.9,1x,a,1x,i4,1x,a)'
  character(len=*), parameter    :: fmt12='(a,3(2x,es16.9))'
  !-------------------------------------------------------------------

    assert ( ctrl%stopc == res_res )

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    call env%info(me, np)

    ! RHS space working vectors
    allocate(r1, mold=b); call r1%default_initialization()
    allocate(r2, mold=b); call r2%default_initialization()
    allocate(y, mold=b); call y%default_initialization()
    call r1%clone(b)
    call r2%clone(b)
    call y%clone(b)

    ! LHS space working vectors
    allocate(v1, mold=x); call v1%default_initialization()
    allocate(v2, mold=x); call v2%default_initialization()
    allocate(w , mold=x); call w%default_initialization()
    allocate(w1, mold=x); call w1%default_initialization()
    allocate(w2, mold=x); call w2%default_initialization()
    call v1%clone(x)
    call v2%clone(x)
    call w%clone( x)
    call w1%clone(x)
    call w2%clone(x)

    debug = .false.
    eps   = epsilon(eps)

    istop    = 0
    ctrl%it  = 0
    Anorm    = zero
    Acond    = zero
    rnorm    = zero
    ynorm    = zero

    !-------------------------------------------------------------------
    ! Set up y and v for the first Lanczos vector v1.
    ! y = beta1 P' v1, where P = C**(-1).
    ! v is really P' v1.
    !-------------------------------------------------------------------

    ! 1) Compute initial residual
    ! 1.a) r=Ax
    call A%apply(x, r1)

    ! 1.b) r=b-r
    call r1%axpby(1.0_rp,b,-1.0_rp)

    ! 2) y=inv(M)r1
    call M%apply(r1, v1)

    ! beta1 = r1 * M^{-1} * r1 = (||r1||_inv(M))^2
    beta1 = r1%dot(v1)

    beta1_lt_zero = ( beta1 < zero )
    call env%bcast(beta1_lt_zero)

    if (beta1_lt_zero) then     ! M must be indefinite.
        istop = 8
        go to 900
    end if

    beta1_eq_zero = ( beta1 == zero )
    call env%bcast(beta1_eq_zero)

    if (beta1_eq_zero) then    ! r1 = 0 exactly.  Stop with x
        istop = 0
        go to 900
    end if

    if ( env%am_i_fine_task() ) then ! Am I a fine task ?
        beta1  = sqrt( beta1 )         ! Normalize y to get v1 later.
    end if

    !-------------------------------------------------------------------
    ! Initialize other quantities.
    !-------------------------------------------------------------------
    oldb   = zero
    beta   = beta1
    dbar   = zero
    epsln  = zero
    qrnorm = beta1
    phibar = beta1
    rhs1   = beta1
    rhs2   = zero
    tnorm2 = zero
    ynorm2 = zero
    cs     = - one
    sn     = zero
    call w%init(0.0_rp)
    call w2%init(0.0_rp)
    call r2%copy(r1)

    if ( env%am_i_fine_task() ) then
        if ( (ctrl%trace > 0) .and. (me == 0) ) call solver_control_log_header(ctrl)
    end if

    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
        ctrl%it = ctrl%it + 1               ! k = itn = 1 first time through

        if ( env%am_i_fine_task() ) then
            !----------------------------------------------------------------
            ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
            ! The general iteration is similar to the case k = 1 with v0 = 0:
            !
            !   p1      = Operator * v1  -  beta1 * v0,
            !   alpha1  = v1'p1,
            !   q2      = p2  -  alpha1 * v1,
            !   beta2^2 = q2'q2,
            !   v2      = (1/beta2) q2.
            !
            ! Again, y = betak P vk,  where  P = C**(-1).
            ! .... more description needed.
            !----------------------------------------------------------------
            s      = one / beta            ! Normalize previous vector (in y).
            ! v      = s*y(1:n)            ! v = vk if P = I
            call v2%scal(s, v1)
            call A%apply(v2, y)

            if (ctrl%it >= 2) then
                ! y   = y - (beta/oldb)*r1    ! call daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
                call y%axpby(-(beta/oldb), r1, 1.0_rp)
            end if

            ! alfa   = dot_product(v,y)      ! alphak
            alfa = v2%dot(y)

            ! y      = y - (alfa/beta)*r2        ! call daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
            call y%axpby(-(alfa/beta), r2, 1.0_rp)

            ! r1     = r2
            call r1%copy(r2)

            ! r2     = y
            call r2%copy(y)

            ! v = inv(M) r2
            call M%apply(r2, v1)

            oldb   = beta                  ! oldb = betak

            ! beta   = dot_product(r2,v)     
            ! beta = betak+1^2
            beta = r2%dot(v1)

            beta_lt_zero = (beta < zero)

            call env%bcast(beta_lt_zero)
            if (beta_lt_zero) then
                istop = 6
                go to 900
            end if

            beta   = sqrt( beta )          ! beta = betak+1
            tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

            if (ctrl%it == 1) then                   ! Initialize a few things.
                if (beta/beta1 <= ten*eps) then   ! beta2 = 0 or ~ 0.
                    istop = -1                     ! Terminate later.
                end if
                !tnorm2 = alfa**2
                gmax   = abs( alfa )              ! alpha1
                gmin   = gmax                     ! alpha1
            end if

            ! Apply previous rotation Qk-1 to get
            !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
            !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

            oldeps = epsln
            delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
            gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
            epsln  =               sn * beta ! epsln2 = 0         epslnk+1
            dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

            ! Compute the next plane rotation Qk

            gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
            cs     = gbar / gamma                ! ck
            sn     = beta / gamma                ! sk
            phi    = cs * phibar                 ! phik
            phibar = sn * phibar                 ! phibark+1

            if (debug) then
                write(*,*) ' '
                write(*,*) 'alfa ', alfa
                write(*,*) 'beta ', beta
                write(*,*) 'gamma', gamma
                write(*,*) 'delta', delta
                write(*,*) 'gbar ', gbar
                write(*,*) 'epsln', epsln
                write(*,*) 'dbar ', dbar
                write(*,*) 'phi  ', phi
                write(*,*) 'phiba', phibar
                write(*,*) ' '
            end if

            ! Update  x.

            denom = one/gamma

!!$        do i = 1, n
!!$           w1(i) = w2(i)
!!$           w2(i) = w(i)
!!$           w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
!!$           x(i)  =   x(i) +   phi * w(i)
!!$        end do

            call w1%copy(w2)
            call w2%copy(w)
            call w%copy(v2)

            call w%axpby(-oldeps, w1, 1.0_rp)
            call w%axpby(-delta, w2, 1.0_rp)
            call w%scal(denom, w)
            call x%axpby(phi, w, 1.0_rp)

            ! Go round again.

            gmax   = max( gmax, gamma )
            gmin   = min( gmin, gamma )
            z      = rhs1 / gamma
            ynorm2 = z**2  +  ynorm2
            rhs1   = rhs2  -  delta * z
            rhs2   =       -  epsln * z

            ! Estimate various norms and test for convergence.

            Anorm  = sqrt( tnorm2 )
            ynorm  = sqrt( ynorm2 )
            epsa   = Anorm * eps
            epsx   = Anorm * ynorm * eps
            epsr   = Anorm * ynorm * ctrl%rtol + ctrl%atol
            ctrl%tol1 = epsr
            ctrl%err1 = phibar 
            diag   = gbar
            if (diag == zero) diag = epsa

            qrnorm = phibar
            rnorml = rnorm
            rnorm  = qrnorm
            rootl  = sqrt( gbar**2 +dbar**2  )  ! norm([gbar; dbar]);
            ! AFM Arnorml     = rnorml*rootl          ! ||A r_{k-1} ||
            ! AFM relArnorml  = rootl  /  Anorm;      ! ||Ar|| / (||A|| ||r||)     
            ! relArnorml = Arnorml / Anorm;           ! ||Ar|| / ||A|| 

            ! Estimate  cond(A).
            ! In this version we look at the diagonals of  R  in the
            ! factorization of the lower Hessenberg matrix,  Q * H = R,
            ! where H is the tridiagonal matrix from Lanczos with one
            ! extra row, beta(k+1) e_k^T.

            Acond  = gmax / gmin

            ! See if any of the stopping criteria are satisfied.
            ! In rare cases, istop is already -1 from above (Abar = const*I).

            if (istop == 0) then
                if (ctrl%it    >= ctrl%itmax    ) istop = 5
                if (Acond  >= 0.1d+0/eps) istop = 4
                if (epsx   >= beta1     ) istop = 3
                ! AFM if (qrnorm <= epsx  .or.  relArnorml <= epsx) istop = 2
                ! AFM if (qrnorm <= epsr  .or.  relArnorml <= epsr) istop = 1
                if (qrnorm <= epsx) istop = 2
                if (qrnorm <= epsr) istop = 1
            end if

            if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)

!!!$        ! See if it is time to print something.
!!!$        if ( me == 0 ) then
!!!$           prnt   = .false.
!!!$           if (ctrl%it    <= 10         ) prnt = .true.
!!!$           if (ctrl%it    >= ctrl%itmax - 10) prnt = .true.
!!!$           if (mod(ctrl%it,10)  ==     0) prnt = .true.
!!!$           if (qrnorm <=  ten * epsx) prnt = .true.
!!!$           if (qrnorm <=  ten * epsr) prnt = .true.
!!!$           if (Acond  >= 1.0d-2/eps ) prnt = .true.
!!!$           if (istop  /=  0         ) prnt = .true.
!!!$
!!!$           if ( prnt ) then
!!!$              if (    ctrl%it     == 1) write(luout_, 1200)
!!!$              write(luout_, 1300) ctrl%it,qrnorm, Anorm, Acond
!!!$              if (mod(ctrl%it,10) == 0) write(luout_, 1500)
!!!$           end if
!!!$        end if

            istop_neq_zero = (istop /= 0)
            call env%bcast(istop_neq_zero)

            if (istop_neq_zero) exit

        else
            ! y = inv(M) r2
            call M%apply(r2, y)
            call env%bcast(beta_lt_zero)

            if (beta_lt_zero) then
                istop = 6
                go to 900
            end if

            call env%bcast(istop_neq_zero)            

            if (istop_neq_zero) exit

        end if

    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

    900 call w2%free()
    call w1%free()
    call w%free()
    call v1%free()
    call v2%free()
    call r2%free()
    call r1%free()

    ! Check for convergence and output corresponding info. messages
    ctrl%converged = .true.

    if ( .not. (istop==0 .or. istop==-1 .or. istop==1 .or. istop==2) ) then
        ctrl%converged = .false.
    end if

    call env%bcast(ctrl%converged)

    if ( env%am_i_fine_task() ) then
        if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    if ( me == 0 ) then
        write(ctrl%luout, 2000) istop, ctrl%it,   &
            Anorm, Acond, &
            rnorm, ynorm
        write(ctrl%luout, 3000) msg(istop)
    end if

    return

    1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'    &
        / ' n      =', i7, 5x, 'checkA =', l4, 12x,         &
        'precon =', l4                                   &
        / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,       &
        'shift  =', e23.14)
    1200 format(// 5x, 'itn', 8x, 'x(1)', 10x,                       &
        'norm(r)', 3x, 'norm(A)', 3X, 'cond(A)')
    1300 format(1p, i8,3e10.2)
    1500 format(1x)
    2000 format(/ 1p, 5x, 'istop =', i3,   14x, 'itn   =', i8     &
        /       5x, 'Anorm =', e12.4, 5x, 'Acond =', e12.4  &
        /       5x, 'rnorm =', e12.4, 5x, 'ynorm =', e12.4, 5x)
    3000 format(      a )

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
end subroutine abstract_pminres




  !=============================================================================
  !
  ! Modified GS with re-orthogonalization
  ! (ideal for serial machines)
  !     ierrc   -- error code
  !                0 : successful return
  !               -1: zero input vector
  !               -2: input vector contains abnormal numbers
  !               -3: input vector is a linear combination of others
  subroutine mgsro(luout,m,bkry,hh,ierrc)
    implicit none
    ! Parameters
    integer(ip)               , intent(in)    :: luout
    integer(ip)               , intent(in)    :: m
    class(vector_t)       , intent(inout) :: bkry(m)
    real(rp)                  , intent(inout) :: hh(m)
    integer(ip)               , intent(out)   :: ierrc  

    ! Locals
    integer(ip)           :: i
    real(rp)              :: norm, thr, fct
    real(rp), parameter   :: reorth = 0.98, rzero = 0.0_rp, rone = 1.0_rp


    ! The last vector is orthogonalized against the others
    norm = bkry(m)%dot(bkry(m))

    if (norm <= rzero) then
       ierrc = -1
       write (luout,*) '** Warning: mgsro: zero input vector'
       return
    else if ( norm > rzero .and. rone/norm > rzero ) then
       ierrc =  0
    else
       ierrc = -2
       write (luout,*) '** Warning: mgsro: input vector contains abnormal numbers' 
       return
    endif

    thr = norm*reorth

    ! Orthogonalize against all the others 
    do i = 1,m-1
       fct =  bkry(m)%dot(bkry(i))
       hh(i) = fct
       call bkry(m)%axpby(-fct,bkry(i),rone)
       ! Reorthogonalization if it is 'too parallel' to the previous vector
       if (fct*fct > thr) then
          fct = bkry(m)%dot(bkry(i))
          hh(i) = hh(i) + fct
          call bkry(m)%axpby(-fct,bkry(i),rone)
       endif
       norm = norm - hh(i)*hh(i)
       if (norm < rzero) norm = rzero
       thr = norm*reorth
    enddo

    ! Normalize
    hh(m) = bkry(m)%nrm2() ! bkry(m)%dot(w)
    ! hh(m) = sqrt(hh(m))


    if ( hh(m) <= rzero ) then
       write (luout,*) '** Warning: mgsro: input vector is a linear combination of previous vectors'
       ierrc = -3
       return
    end if

    call bkry(m)%scal(rone/hh(m),bkry(m))
  end subroutine mgsro

  !=============================================================================
  !
  ! Iterative Classical Gram-Schmidt (ideal for distributed-memory machines)
  !
  ! Taken from Figure 5 of the following paper: 
  ! Efficient Gram-Schmidt orthonormalisation on parallel computers
  ! F. Commun. Numer. Meth. Engng. 2000; 16:57-66
  !
  ! *** IMPORTANT NOTE ***: It would be fine to have an implementation
  ! of icgs that allows to exploit level 2 BLAS for the following 2 operations
  ! (see below):
  !    x  p <- Q_k^T * Q(k)
  !    x  Q(k) <- Q(k) - Q_k * p
  ! This implementation would imply defining a generic
  ! data structure for storing Krylov subspace bases.
  !     ierrc   -- error code
  !                0 : successful return
  !               -1: zero input vector
  !               -2: input vector contains abnormal numbers
  !               -3: input vector is a linear combination of others
  subroutine icgsro (luout, k, Q, s, ierrc)
    implicit none
    ! Parameters  
    integer(ip), intent(in)                   :: luout
    integer(ip)               , intent(in)    :: k
    class(vector_t)       , intent(inout) :: Q(k)
    real(rp)                  , intent(inout) :: s(k)
    integer(ip)               , intent(out)   :: ierrc  


    ! Locals 
    real(rp), parameter              :: alpha = 0.5_rp, rone = 1.0_rp, rzero = 0.0_rp
    logical, parameter           :: debug = .false.
    integer(ip)                      :: i, j, m
    real(rp)                         :: p(k-1)
    real(rp)                         :: delta_i, delta_i_mone, beta_i


    s = 0.0_rp
    i = 1 

!    call generic_krylov_basis_extract_view ( k, Q, q_k )
!    Q(k)

    delta_i_mone = Q(k)%nrm2()

    if (delta_i_mone <= rzero) then
       ierrc = -1	
       write (luout,*) '** Warning: icgsro: zero input vector'
       return
    else if ( delta_i_mone > rzero .and. rone/delta_i_mone > rzero ) then
       ierrc =  0
    else
       ierrc = -2
       write (luout,*) '** Warning: icgsro: input vector contains abnormal numbers' 
       return
    endif

    loop_icgs: do


       do j=1, k-1
           ! p <- Q_k^T * Q(k), with Q_k = (Q(1), Q(2), .. Q(k-1))  
           p(j) = Q(k)%dot(Q(j))
           ! q_k <- q_k - Q_k * p 
           call Q(k)%axpby(-p(j),Q(j),1.0_rp)
       enddo


       s(1:k-1) = s(1:k-1) + p(1:k-1)

       ! OPTION 1: estimate ||q_k^{i+1}||
       ! ** IMPORTANT NOTE **: I have observed,
       ! for element-based data decompositions, when
       ! q_k^{i+1} is very close to zero, that this estimation
       ! becomes zero.
       ! beta_i = 1 - ||p||^2/(delta_i_mone^2)
       ! beta_i = 1 - (dot_product(p, p))/(delta_i_mone*delta_i_mone)
       ! delta_i = delta_i_mone * sqrt ( abs (beta_i) )

       ! OPTION 2: calculate ||q_k^{i+1}||
       ! write (*,*) i, k, 'xxx', aux,  delta_i      ! DBG:
       delta_i = Q(k)%nrm2()

       if ( delta_i > delta_i_mone * alpha .or. delta_i == rzero ) exit loop_icgs
       i = i + 1
       delta_i_mone = delta_i 
    end do loop_icgs

    if (debug) then
       write(luout,'(a,i1,a)') 'ICGS: ', i, ' ortho. iterations'
       write(luout,*) delta_i, delta_i_mone * alpha
    end if
    if (i>2)   write(luout,*) '** Warning: icgsro: required more than two iterations!'

    ! Normalize
    ! call generic_nrm2 ( q_k, s(k) )
    s(k) = delta_i 
    if ( s(k) <= rzero ) then
       write (luout,*) '** Warning: icgsro: input vector is a linear combination of previous vectors'
       ierrc = -3
       return
    end if

    call Q(k)%scal(rone/s(k),Q(k))

  end subroutine icgsro

  subroutine solver_control_allocate_conv_his( ctrl )
    implicit none
    type(solver_control_t), intent(inout) :: ctrl

    if (ctrl%track_conv_his) then
       call memalloc(ctrl%itmax,ctrl%err1h,__FILE__,__LINE__)
       call memalloc(ctrl%itmax,ctrl%err2h,__FILE__,__LINE__)
       ctrl%err1h=0.0_rp
       ctrl%err2h=0.0_rp
    end if
  end subroutine solver_control_allocate_conv_his

  subroutine solver_control_free_conv_his( ctrl )
    implicit none
    type(solver_control_t), intent(inout) :: ctrl

    if (ctrl%track_conv_his) then
       call memfree(ctrl%err1h,__FILE__,__LINE__)
       call memfree(ctrl%err2h,__FILE__,__LINE__)
    end if
  end subroutine solver_control_free_conv_his

  subroutine solver_control_log_header( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter    :: fmt1='(a18,1x,a4,3(2x,a15))'
    character(len=*), parameter    :: fmt2='(a18,1x,a4,3(2x,a15),3(2x,a15))'
    integer, parameter             :: outlen=18 
    character(len=len(methdname))  :: mname
    character(len=outlen)          :: outname

    mname = adjustl(trim(methdname(ctrl%method)))
    write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'

    select case(ctrl%stopc)

    case ( delta_rhs, delta_delta, res_res, res_rhs, res_nrmgiven_rhs_nrmgiven, &
         & res_nrmgiven_res_nrmgiven )

       write(ctrl%luout,fmt1) adjustl(outname),'Iteration','Error Estimate','Tolerance'

    case ( delta_rhs_and_res_res, delta_rhs_and_res_rhs,  & 
         delta_delta_and_res_res, delta_delta_and_res_rhs )

       write(ctrl%luout,fmt2) adjustl(outname), 'Iteration', 'Error Estimate', 'Tolerance', &
            & 'Error Estimate', 'Tolerance' 

    case default
       ! Write an error message and stop ?      
    end select
  end subroutine solver_control_log_header

  subroutine solver_control_log_conv ( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter   :: fmt1='(a18,1x,i4,3(2x,es16.9))'
    character(len=*), parameter   :: fmt2='(a18,1x,i4,3(2x,es16.9),3(2x,es16.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname

    if ( (mod(ctrl%it,ctrl%trace) == 0).or.ctrl%converged.or.(ctrl%it>=ctrl%itmax)) then 
       mname = adjustl(trim(methdname(ctrl%method)))
       write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
       select case(ctrl%stopc)
       case ( delta_rhs, delta_delta, res_res, res_rhs, & 
            & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
          write(ctrl%luout,fmt1) adjustl(outname), ctrl%it, ctrl%err1, ctrl%tol1
       case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          write(ctrl%luout,fmt2) adjustl(outname), ctrl%it, ctrl%err1, ctrl%tol1, &
               &  ctrl%err2, ctrl%tol2 
       case default
          ! Write an error message and stop ?
       end select
    endif

  end subroutine solver_control_log_conv

  subroutine solver_control_log_end ( ctrl )
    implicit none
    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    character(len=*), parameter  :: fmt11='(a,2x,es16.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt12='(a,3(2x,es16.9))'

    character(len=*), parameter  :: fmt21='(a,2x,es16.9,1x,es16.9,1x,a,1x,i4,1x,a)'
    character(len=*), parameter  :: fmt22='(a,3(2x,es16.9),3(2x,es16.9))'


    select case( ctrl%stopc )

    case ( delta_rhs,delta_delta,res_res,res_rhs,&
         & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven)

       if ( ctrl%converged ) then
          write(ctrl%luout,fmt11) trim(methdname(ctrl%method))//' converged to ', &
               & ctrl%tol1,' in ',ctrl%it,' iterations. '
          write(ctrl%luout,fmt12) 'Last iteration error estimate: ', ctrl%err1
       else
          write(ctrl%luout,fmt11) trim(methdname(ctrl%method))//' failed to converge to ', &
               & ctrl%tol1,' in ',ctrl%it,' iterations. '
          write(ctrl%luout,fmt12) 'Last iteration error estimate: ', ctrl%err1 
       end if

    case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
         & delta_delta_and_res_res, delta_delta_and_res_rhs )

       if ( ctrl%converged ) then
          write(ctrl%luout,fmt21) trim(methdname(ctrl%method))//' converged to ', &
               & ctrl%tol1, ctrl%tol2, ' in ', ctrl%it ,' iterations. '
          write(ctrl%luout,fmt22) 'Last iteration error estimates: ', ctrl%err1, ctrl%err2
       else
             write(ctrl%luout,fmt21) trim(methdname(ctrl%method))//' failed to converge to ', &
                  & ctrl%tol1, ctrl%tol2, ' in ', ctrl%it ,' iterations. '
             write(ctrl%luout,fmt22) 'Last iteration error estimates: ', ctrl%err1, ctrl%err2
       end if

    case default
       ! Write an error message and stop ?      

    end select

  end subroutine solver_control_log_end

  subroutine solver_control_log_conv_his ( ctrl )
    implicit none 

    ! Parameters
    type(solver_control_t), intent(in) :: ctrl

    ! Local variables
    character(len=*), parameter   :: fmt1='(a18,1x,i4,3(2x,es16.9))'
    character(len=*), parameter   :: fmt2='(a18,1x,i4,3(2x,es16.9),3(2x,es16.9))'
    integer, parameter            :: outlen=18 
    character(len=len(methdname)) :: mname
    character(len=outlen)         :: outname
    integer(ip)                   :: i

    if (ctrl%track_conv_his) then
       call solver_control_log_header(ctrl)

       mname = adjustl(trim(methdname(ctrl%method)))
       write(outname,'(a)') mname(1:min(len_trim(mname),outlen-1))//':'
       select case(ctrl%stopc)
       case ( delta_rhs, delta_delta, res_res, res_rhs, & 
            & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
          do i=1,ctrl%it
             write(ctrl%luout,fmt1) adjustl(outname), i, ctrl%err1h(i), ctrl%tol1
          end do
       case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, & 
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          do i=1,ctrl%it
             write(ctrl%luout,fmt2) adjustl(outname), i, ctrl%err1h(i), ctrl%tol1, &
                  &                                      ctrl%err2h(i), ctrl%tol2 
          end do
       case default
          ! Write an error message and stop ?
       end select

       call solver_control_log_end(ctrl)
    end if

  end subroutine solver_control_log_conv_his

  !Taken from SPARSKIT
  subroutine givens(x,y,c,s)
    real(rp) x,y,c,s

    !     Given x and y, this subroutine generates a Givens' rotation c, s.
    !     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
    !     (See P 202 of "matrix computations" by Golub and van Loan.)

    real(rp) t,one,rzero
    parameter (rzero=0.0_rp,one=1.0_rp)

    if (x.eq.rzero .and. y.eq.rzero) then
       c = one
       s = zero
    else if (abs(y).gt.abs(x)) then
       t = x / y
       x = sqrt(one+t*t)
       s = sign(one / x, y)
       c = t*s
    else if (abs(y).le.abs(x)) then
       t = y / x
       y = sqrt(one+t*t)
       c = sign(one / y, x)
       s = t*c
    else
       ! X or Y must be an invalid floating-point number, set both to zero
       x = rzero
       y = rzero
       c = one
       s = rzero
    endif
    x = abs(x*y)
    y = rzero
  end subroutine givens

  
end module abstract_solver_names
