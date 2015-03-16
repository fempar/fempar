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
module abstract_solver
  use solver_base 
  use base_operator_names
  use base_operand_names
  implicit none

  private
  
  public :: abstract_solve

contains

  !=============================================================================
  ! Abstract Solve (public driver)
  !=============================================================================
  subroutine abstract_solve(A,M,b,x,ctrl)
    implicit none
    class(base_operator), intent(inout) :: A        ! Matrix
    class(base_operator), intent(inout) :: M        ! Preconditioner
    class(base_operand) , intent(in)    :: b        ! RHS
    class(base_operand) , intent(inout) :: x        ! Approximate solution
    type(solver_control), intent(inout) :: ctrl
    
    call solver_control_allocate_conv_his ( ctrl )

    ! Invoke Krylov Subspace Solver
    select case( ctrl%method )
    case ( cg )
       call abstract_pcg ( A, M, b, x, ctrl )
    case ( lgmres )
       ! call abstract_plgmres ( A, M, b, x, ctrl )
    case ( rgmres )
       ! call abstract_prgmres ( A, M, b, x, ctrl )
    case ( fgmres )
       ! call abstract_pfgmres ( A, M, b, x, ctrl )
    case ( richard )
       ! call abstract_prichard ( A, M, b, x, ctrl )
    case( direct )
       call M%apply(b, x)
    case ( icg )
       ! call abstract_ipcg ( A, M, b, x, ctrl )
    case ( lfom )
       ! call abstract_plfom ( A, M, b, x, ctrl )
    case ( minres )
       ! call abstract_pminres ( A, M, b, x, ctrl )
    case default
       ! Write an error message and stop ?      
    end select
  end subroutine abstract_solve

  !=============================================================================
  ! Abstract Preconditioned Conjugate Gradient
  !=============================================================================
  subroutine abstract_pcg( A, M, b, x, ctrl )
    !-----------------------------------------------------------------------------
    ! This routine performs pcg iterations on Ax=b with preconditioner M. 
    !-----------------------------------------------------------------------------
    implicit none

    ! Parameters
    class(base_operator) , intent(inout) :: A        ! Matrix
    class(base_operator) , intent(inout) :: M        ! Preconditioner
    class(base_operand)  , intent(in)    :: b        ! RHS
    class(base_operand)  , intent(inout) :: x        ! Approximate solution
    type(solver_control) , intent(inout) :: ctrl     ! Control data

    ! Locals
    real(rp)                         :: r_nrm_M      ! |r|_inv(M) 
    real(rp)                         :: b_nrm_M      ! |b|_inv(M)
    real(rp)                         :: r_z, Ap_p, alpha, beta
    integer                          :: me, np
    class(base_operand), allocatable :: r,p,Ap,z     ! Working vectors

    call b%GuardTemp()

    allocate ( r, mold = x )
    allocate ( Ap, mold = x )
    allocate ( z, mold = x )
    allocate ( p, mold = x )

    call r%clone ( x  )
    call Ap%clone ( x )
    call z%clone ( x  )
    call p%clone ( x  )

    call A%info(me,np)

    ! Evaluate |b|_inv(M) if required
    if ( ctrl%stopc == res_nrmgiven_rhs_nrmgiven ) then
       ! r = inv(M) b
       call M%apply(b, r)
       b_nrm_M = b%dot(r)
       if ( M%am_i_fine_task() ) then ! Am I a fine task ?
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

    if ( M%am_i_fine_task() ) then ! Am I a fine task ?
       r_nrm_M = sqrt( r_z )
    end if

    ! 4) Initializations:
    ! p=z
    call p%copy(z) 

    if ( M%am_i_fine_task() ) then
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

       if (M%am_i_fine_task()) then
          ! Check and log convergence
          call pcg_conv_check(r, r_nrm_M, alpha, p, ctrl )
          if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
       end if
       ! Send converged to coarse-grid tasks
       call M%bcast(ctrl%converged)
       if(ctrl%converged.or.(ctrl%it>=ctrl%itmax)) exit loop_pcg

       ! z = inv(M) r
       call M%apply(r,z)

       if ( M%am_i_fine_task()) then ! Am I a fine task ?
          beta = 1.0_rp/r_z
       else
          beta = 0.0_rp 
       end if

       r_z=r%dot(z)
       beta = beta*r_z

       if (M%am_i_fine_task()) then ! Am I a fine task ?
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


    if ( M%am_i_fine_task() ) then
       if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call b%CleanTemp()
  end subroutine abstract_pcg

  subroutine pcg_conv_init ( b, r , nrm_b_given, nrm_r_given, ctrl )
    implicit none 
    class(base_operand) , intent(in)   :: b, r
    real(rp)            , intent(in)    :: nrm_b_given, nrm_r_given 
    type(solver_control), intent(inout) :: ctrl

    select case(ctrl%stopc)
    case ( delta_rhs, res_rhs, delta_rhs_and_res_rhs, delta_delta_and_res_rhs )
       ! call generic_nrm2 ( b, ctrl%bn2 )
       ctrl%bn2 = b%nrm2()
       ctrl%tol1 = ctrl%rtol* ctrl%bn2 + ctrl%atol
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
    class(base_operand) , intent(in)    :: r, p
    real(rp)            , intent(in)    :: nrm_r_given
    real(rp)            , intent(in)    :: alpha
    type(solver_control), intent(inout) :: ctrl

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
    ctrl%err1h(ctrl%it) = ctrl%err1
    ctrl%converged = (ctrl%err1 <= ctrl%tol1 )

    ! Evaluate 2nd convergence criterion
    select case( ctrl%stopc )
    case ( delta_rhs_and_res_res,   delta_rhs_and_res_rhs, &
         & delta_delta_and_res_res, delta_delta_and_res_rhs )
       if ( ctrl%converged ) then
          ! Compute || r(i) ||
          ! call generic_nrm2 ( r, ctrl%err2 )
          ctrl%err2 = r.nrm2()
          ctrl%converged = (ctrl%err2 <= ctrl%tol2 )
          ctrl%err2h(ctrl%it) = ctrl%err2
       end if
    end select

    return
  end subroutine pcg_conv_check

  
end module abstract_solver
