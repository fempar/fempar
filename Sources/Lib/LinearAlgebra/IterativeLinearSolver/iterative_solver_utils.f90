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

module iterative_linear_solver_utils_names

  use types_names
  use vector_names
  use multivector_names

contains

    !=============================================================================
  !
  ! Modified GS with re-orthogonalization
  ! (ideal for serial machines)
  !     ierrc   -- error code
  !                0 : successful return
  !               -1: zero input vector
  !               -2: input vector contains abnormal numbers
  !               -3: input vector is a linear combination of others
  subroutine modified_gs_reorthonorm(luout,m,bkry,hh,ierrc)
    implicit none
    ! Parameters
    integer(ip)               , intent(in)    :: luout
    integer(ip)               , intent(in)    :: m
    class(multivector_t)      , intent(inout) :: bkry
    real(rp)                  , intent(inout) :: hh(m)
    integer(ip)               , intent(out)   :: ierrc  

    ! Locals
    integer(ip)              :: i
    real(rp)                 :: norm, thr, fct
    real(rp), parameter      :: reorth = 0.98, rzero = 0.0_rp, rone = 1.0_rp
    class(vector_t), pointer :: bkrym, bkryi

    bkrym => bkry%get(m) 
    
    ! The last vector is orthogonalized against the others
    norm = bkrym%dot(bkrym)

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
       bkryi => bkry%get(i)
       fct =  bkrym%dot(bkryi)
       hh(i) = fct
       call bkrym%axpby(-fct,bkryi,rone)
       ! Reorthogonalization if it is 'too parallel' to the previous vector
       if (fct*fct > thr) then
          fct = bkrym%dot(bkryi)
          hh(i) = hh(i) + fct
          call bkrym%axpby(-fct,bkryi,rone)
       endif
       norm = norm - hh(i)*hh(i)
       if (norm < rzero) norm = rzero
       thr = norm*reorth
    enddo

    ! Normalize
    hh(m) = bkrym%nrm2()

    if ( hh(m) <= rzero ) then
       write (luout,*) '** Warning: mgsro: input vector is a linear combination of previous vectors'
       ierrc = -3
       return
    end if

    call bkrym%scal(rone/hh(m),bkrym)
  end subroutine modified_gs_reorthonorm

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
  subroutine iterative_gs_reorthonorm (luout, k, Q, s, ierrc)
    implicit none
    ! Parameters  
    integer(ip), intent(in)                   :: luout
    integer(ip)               , intent(in)    :: k
    class(multivector_t)      , intent(inout) :: Q
    real(rp)                  , intent(inout) :: s(k)
    integer(ip)               , intent(out)   :: ierrc  


    ! Locals 
    real(rp), parameter              :: alpha = 0.5_rp, rone = 1.0_rp, rzero = 0.0_rp
    logical, parameter               :: debug = .false.
    integer(ip)                      :: i, j, m
    real(rp)                         :: p(k-1)
    real(rp)                         :: delta_i, delta_i_mone, beta_i
    class(vector_t), pointer         :: q_k

    s = 0.0_rp
    i = 1 


    q_k => Q%get(k)
    delta_i_mone = q_k%nrm2()

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
       ! p <- Q_k^T * Q(k), with Q_k = (Q(1), Q(2), .. Q(k-1))  
       call Q%multidot(k-1,q_k,p)
       
       ! q_k <- q_k - Q_k * p 
       call Q%multiaxpy(k-1,q_k,-1.0_rp,p)

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
       delta_i = q_k%nrm2()

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
    s(k) = delta_i 
    if ( s(k) <= rzero ) then
       write (luout,*) '** Warning: icgsro: input vector is a linear combination of previous vectors'
       ierrc = -3
       return
    end if

    call q_k%scal(rone/s(k),q_k)
    
  end subroutine iterative_gs_reorthonorm
  
  ! Source: SPARSKIT
  ! Given x and y, this subroutine generates a Givens' rotation c, s.
  ! And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
  ! (See P 202 of "matrix computations" by Golub and van Loan)
  subroutine apply_givens_rotation(x,y,c,s)
    implicit none
    real(rp), intent(inout) :: x
    real(rp), intent(inout) :: y
    real(rp), intent(inout) :: c
    real(rp), intent(inout) :: s

    real(rp) :: t, one, rzero
    parameter (rzero=0.0_rp,one=1.0_rp)

    if (x.eq.rzero .and. y.eq.rzero) then
       c = one
       s = rzero
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
  end subroutine apply_givens_rotation


end module iterative_linear_solver_utils_names
