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

# include "debug.i90"

  implicit none

  private
  
  public :: abstract_solve

contains

  !=============================================================================
  ! Abstract Solve (public driver)
  !=============================================================================
  subroutine abstract_solve(A,M,b,x,ctrl)
    implicit none
    class(base_operator), intent(in) :: A        ! Matrix
    class(base_operator), intent(in) :: M        ! Preconditioner
    class(base_operand) , intent(in)    :: b        ! RHS
    class(base_operand) , intent(inout) :: x        ! Approximate solution
    type(solver_control), intent(inout) :: ctrl
    
    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()
    
    call solver_control_allocate_conv_his ( ctrl )

    ! Invoke Krylov Subspace Solver
    select case( ctrl%method )
    case ( cg )
       call abstract_pcg ( A, M, b, x, ctrl )
    case ( lgmres )
       call abstract_plgmres ( A, M, b, x, ctrl )
    case ( rgmres )
       check(1==0)
       ! call abstract_prgmres ( A, M, b, x, ctrl )
    case ( fgmres )
       check(1==0)
       ! call abstract_pfgmres ( A, M, b, x, ctrl )
    case ( richard )
       check(1==0)
       ! call abstract_prichard ( A, M, b, x, ctrl )
    case( direct )
       call M%apply(b, x)
    case ( icg )
       check(1==0)
       ! call abstract_ipcg ( A, M, b, x, ctrl )
    case ( lfom )
       check(1==0)
       ! call abstract_plfom ( A, M, b, x, ctrl )
    case ( minres )
       check(1==0)
       ! call abstract_pminres ( A, M, b, x, ctrl )
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
  subroutine abstract_pcg( A, M, b, x, ctrl )
    !-----------------------------------------------------------------------------
    ! This routine performs pcg iterations on Ax=b with preconditioner M. 
    !-----------------------------------------------------------------------------
    implicit none

    ! Parameters
    class(base_operator) , intent(in)    :: A        ! Matrix
    class(base_operator) , intent(in)    :: M        ! Preconditioner
    class(base_operand)  , intent(in)    :: b        ! RHS
    class(base_operand)  , intent(inout) :: x        ! Approximate solution
    type(solver_control) , intent(inout) :: ctrl     ! Control data

    ! Locals
    real(rp)                         :: r_nrm_M      ! |r|_inv(M) 
    real(rp)                         :: b_nrm_M      ! |b|_inv(M)
    real(rp)                         :: r_z, Ap_p, alpha, beta
    integer                          :: me, np
    class(base_operand), allocatable :: r,p,Ap,z     ! Working vectors

    call A%GuardTemp()
    call M%GuardTemp()
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

    call A%CleanTemp()
    call M%CleanTemp()
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
          ctrl%err2 = r%nrm2()
          ctrl%converged = (ctrl%err2 <= ctrl%tol2 )
          ctrl%err2h(ctrl%it) = ctrl%err2
       end if
    end select

    return
  end subroutine pcg_conv_check

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
  subroutine abstract_plgmres ( A, M, b, x, ctrl ) 
    !--------------------------------------------------------------------
    ! This routine performs plgmres iterations on Ax=b with preconditioner M. 
    !--------------------------------------------------------------------
#ifdef ENABLE_BLAS       
    use blas77_interfaces
#endif
    implicit none
    ! Parameters
    class(base_operator)  , intent(in) :: A ! Matrix
    class(base_operator)  , intent(in) :: M ! Preconditioner
    class(base_operand)   , intent(inout) :: x ! Solution
    class(base_operand)   , intent(in)    :: b ! RHS
    type(solver_control), intent(inout) :: ctrl

    integer(ip)                    :: ierrc
    integer(ip)                    :: kloc, kloc_aux, i, j, k_hh, id
    real(rp)                       :: res_norm, res_2_norm, rhs_norm
    real(rp)                       :: alpha, c, s 
    real(rp)   , allocatable       :: hh(:,:), g(:), g_aux(:), cs(:,:)
    integer                        :: me, np
    logical                        :: exit_loop

    class(base_operand), allocatable  :: r, z     ! Working vectors
    class(base_operand), allocatable  :: q_aux    
    class(base_operand), allocatable  :: bkry(:)  ! Krylov basis

    assert(ctrl%stopc==res_res.or.ctrl%stopc==res_rhs.or.ctrl%stopc==res_nrmgiven_res_nrmgiven.or.ctrl%stopc==res_nrmgiven_rhs_nrmgiven)

    call A%GuardTemp()
    call M%GuardTemp()
    call b%GuardTemp()

    ! Clone x in order to allocate working vectors
    allocate(r, mold=b)
    allocate(z, mold=x)
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

    call A%info(me, np)

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
    exit_loop = (ctrl%err1 < ctrl%tol1)

    ! Send converged to coarse-grid tasks
    call M%bcast(exit_loop )


    if ( M%am_i_fine_task() ) then
       if ( (ctrl%trace > 0) .and. (me == 0) ) call solver_control_log_header(ctrl)
    end if


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
       call bkry(1)%clone(x)
       call bkry(1)%scal(1.0_rp/res_norm, z)

       ! residual in the krylov basis
       g(1)            = res_norm
       g(2:ctrl%dkrymax+1) = 0.0_rp

       ! start iterations
       kloc = 0
       inner: do while ( (.not.exit_loop) .and. &
            &            (ctrl%it < ctrl%itmax) .and. &
            &            (kloc < ctrl%dkrymax))
          kloc  = kloc  + 1
          ctrl%it = ctrl%it + 1

          ! Generate new basis vector
          call A%apply( bkry(kloc), r )
          call bkry(kloc+1)%clone(x)
          call M%apply(r, bkry(kloc+1) )


          if ( M%am_i_fine_task() ) then ! Am I a fine task ?
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
                call M%bcast(exit_loop) 
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
                      call M%bcast(exit_loop)
                      call M%bcast(exit_loop)
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

          ctrl%err1h(ctrl%it) = ctrl%err1
          exit_loop = (ctrl%err1 < ctrl%tol1)
          ! Send converged to coarse-grid tasks
          call M%bcast(exit_loop)

          if ( M%am_i_fine_task() ) then
             if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_conv(ctrl)
          end if
       end do inner


       if ( M%am_i_fine_task() ) then ! Am I a fine task ?
          if ( ierrc == -2 ) then
             write (ctrl%luout,*) '** Warning: LGMRES: ortho failed due to abnormal numbers, no way to proceed'
             ! The coarse-grid task should exit 
             ! the outer-do loop. Send signal. 
             exit_loop = .true.
             call M%bcast(exit_loop)
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
                   call M%bcast(exit_loop)
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

       exit_loop = (ctrl%err1 < ctrl%tol1)
       ! Send converged to coarse-grid tasks
       call M%bcast(exit_loop)
    end do outer

    ctrl%converged = (ctrl%err1 < ctrl%tol1)
    ! Send converged to coarse-grid tasks
    call M%bcast( ctrl%converged )

    call r%free()
    call z%free()
    deallocate(r)
    deallocate(z)

    ! Deallocate Krylov basis
    do i=1, ctrl%dkrymax+1
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

    if ( M%am_i_fine_task() ) then
       if ((me == 0).and.(ctrl%trace/=0)) call solver_control_log_end(ctrl)
    end if

    call A%CleanTemp()
    call M%CleanTemp()
    call b%CleanTemp()
  end subroutine abstract_plgmres


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
    class(base_operand)       , intent(inout) :: bkry(m)
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
    class(base_operand)       , intent(inout) :: Q(k)
    real(rp)                  , intent(inout) :: s(k)
    integer(ip)               , intent(out)   :: ierrc  


    ! Locals 
    real(rp), parameter              :: alpha = 0.5_rp, rone = 1.0_rp, rzero = 0.0_rp
    logical(lg), parameter           :: debug = .false.
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

  
end module abstract_solver
