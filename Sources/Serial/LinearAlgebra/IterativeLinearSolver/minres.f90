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

module minres_names
  use types_names
  use stdio_names
  use memor_names
  
  ! Abstract modules
  use vector_names
  use vector_space_names
  use operator_names
  use environment_names
  use base_iterative_linear_solver_names
  use iterative_linear_solver_parameters_names
  use ParameterList

  implicit none
# include "debug.i90"
  private

  integer (ip), parameter :: default_minres_stopping_criteria = res_res

  type, extends(base_iterative_linear_solver_t) :: minres_t
    ! Working space vectors for type(minres_t)
    class(vector_t), allocatable :: r1, r2
    class(vector_t), allocatable :: v1, v2
    class(vector_t), allocatable :: w, w1, w2
    class(vector_t), allocatable :: y
  contains
    procedure          :: allocate_workspace            => minres_allocate_workspace
    procedure          :: free_workspace                => minres_free_workspace
    procedure          :: set_parameters_from_pl        => minres_set_parameters_from_pl
    procedure          :: solve_body                    => minres_solve_body
    procedure          :: supports_stopping_criteria    => minres_supports_stopping_criteria
    procedure          :: get_default_stopping_criteria => minres_get_default_stopping_criteria
  end type
  
  ! Data types
  public :: create_minres
  
contains
  subroutine minres_allocate_workspace(this)
    implicit none
    class(minres_t), intent(inout) :: this
    type(vector_space_t), pointer :: range
    type(dynamic_state_operator_t), pointer :: A, M
    A => this%get_A()
    range  => A%get_range_vector_space()
    call range%create_vector(this%r1)
    call range%create_vector(this%r2)
    call range%create_vector(this%y)

    M => this%get_M()
    range  => M%get_range_vector_space()
    call range%create_vector(this%v1)
    call range%create_vector(this%v2)
    call range%create_vector(this%w)
    call range%create_vector(this%w1)
    call range%create_vector(this%w2)
  end subroutine minres_allocate_workspace
  
  subroutine minres_free_workspace(this)
    implicit none
    class(minres_t), intent(inout) :: this
    call this%r1%free()
    call this%r2%free()
    call this%y%free()
    call this%v1%free() 
    call this%v2%free() 
    call this%w%free()
    call this%w1%free() 
    call this%w2%free() 
   
    deallocate(this%r1)
    deallocate(this%r2)
    deallocate(this%y)
    deallocate(this%v1)
    deallocate(this%v2)
    deallocate(this%w)
    deallocate(this%w1)
    deallocate(this%w2)
  end subroutine minres_free_workspace

  subroutine minres_set_parameters_from_pl(this, parameter_list) 
   implicit none
   class(minres_t),       intent(inout) :: this
   type(ParameterList_t), intent(in)    :: parameter_list
   call this%base_iterative_linear_solver_set_parameters_from_pl(parameter_list)
  end subroutine minres_set_parameters_from_pl
  
  subroutine minres_solve_body(this,b,x)
    implicit none
    class(minres_t), intent(inout) :: this
    class(vector_t)    , intent(in) :: b
    class(vector_t)    , intent(inout) :: x 
 
    ! Local variables to store a copy/reference of the corresponding member variables of base class
    class(environment_t), pointer :: environment
    class(operator_t)   , pointer :: A, M 
    class(vector_t)     , pointer :: initial_solution
    integer(ip)                   :: stopping_criteria, max_num_iterations, output_frequency, luout
    real(rp)                      :: atol, rtol
    logical                       :: track_convergence_history

    ! Pointers to freely modify/read private member variables of base class
    integer(ip), pointer :: num_iterations
    logical    , pointer :: did_converge
    real(rp)   , pointer :: rhs_convergence_test, error_estimate_convergence_test
    real(rp)   , pointer :: error_estimate_history_convergence_test(:)

    ! Local variables 
      real(rp)  :: alpha  , beta  , beta1 , cs    ,          &
       dbar  , delta , denom , diag  ,          &
       eps   , epsa  , epsln , epsr  , epsx  ,  &
       gamma , gbar  , gmax  , gmin  ,          &
       oldb  , oldeps, qrnorm, phi   , phibar,  &
       rhs1  , rhs2  , rnorml, rootl ,          &
       s     , sn    , t     , tnorm2, ynorm2, z

  integer(ip) :: i, istop, itn
  logical     :: debug, prnt
  logical     :: beta1_lt_zero, beta1_eq_zero, beta_lt_zero, istop_neq_zero

  real(rp) ::   Anorm, Acond, rnorm, ynorm
  integer  ::   me, np


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

! ------------------------------------------------------------------------------
   

    environment               => this%get_environment()
    A                         => this%get_A()
    M                         => this%get_M()
    initial_solution          => this%get_initial_solution()
    luout                     =  this%get_luout()
    stopping_criteria         =  this%get_stopping_criteria()
    max_num_iterations        =  this%get_max_num_iterations()
    atol                      =  this%get_atol()
    rtol                      =  this%get_rtol()
    output_frequency          =  this%get_output_frequency()
    track_convergence_history =  this%get_track_convergence_history()
     
    num_iterations                          => this%get_pointer_num_iterations()
    did_converge                            => this%get_pointer_did_converge()
    rhs_convergence_test                    => this%get_pointer_rhs_convergence_test()
    error_estimate_convergence_test         => this%get_pointer_error_estimate_convergence_test()
    error_estimate_history_convergence_test => this%get_pointer_error_estimate_history_convergence_test()

    debug = .false.
    eps   = epsilon(eps)

    istop           = 0
    num_iterations  = 0
    Anorm           = zero
    Acond           = zero
    rnorm           = zero
    ynorm           = zero

    call x%copy(initial_solution)
    !-------------------------------------------------------------------
    ! Set up y and v for the first Lanczos vector v1.
    ! y = beta1 P' v1, where P = C**(-1).
    ! v is really P' v1.
    !-------------------------------------------------------------------

    ! 1) Compute initial residual
    ! 1.a) r=Ax
    call A%apply(x, this%r1)

    ! 1.b) r=b-r
    call this%r1%axpby(1.0_rp,b,-1.0_rp)

    ! 2) y=inv(M)r1
    call M%apply(this%r1, this%v1)

    ! beta1 = r1 * M^{-1} * r1 = (||r1||_inv(M))^2
    beta1 = this%r1%dot(this%v1)

    beta1_lt_zero = ( beta1 < zero )
    call environment%bcast(beta1_lt_zero)

    if (beta1_lt_zero) then     ! M must be indefinite.
        istop = 8
        go to 900
    end if

    beta1_eq_zero = ( beta1 == zero )
    call environment%bcast(beta1_eq_zero)

    if (beta1_eq_zero) then    ! r1 = 0 exactly.  Stop with x
        istop = 0
        go to 900
    end if

    if ( environment%am_i_l1_task() ) then ! Am I a fine task ?
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
    call this%w%init(0.0_rp)
    call this%w2%init(0.0_rp)
    call this%r2%copy(this%r1)
 
    call this%print_convergence_history_header(luout)
 
    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
        num_iterations = num_iterations + 1               ! k = itn = 1 first time through

        if ( environment%am_i_l1_task() ) then
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
            call this%v2%scal(s, this%v1)
            call A%apply(this%v2, this%y)

            if (num_iterations >= 2) then
                ! y   = y - (beta/oldb)*r1    ! call daxpy ( n, (- beta/oldb), r1, 1, y, 1 )
                call this%y%axpby(-(beta/oldb), this%r1, 1.0_rp)
            end if

            ! alpha   = dot_product(v,y)      ! alphak
            alpha = this%v2%dot(this%y)

            ! y      = y - (alpha/beta)*r2        ! call daxpy ( n, (- alpha/beta), r2, 1, y, 1 )
            call this%y%axpby(-(alpha/beta), this%r2, 1.0_rp)

            ! r1     = r2
            call this%r1%copy(this%r2)

            ! r2     = y
            call this%r2%copy(this%y)

            ! v = inv(M) r2
            call M%apply(this%r2, this%v1)

            oldb   = beta                  ! oldb = betak

            ! beta   = dot_product(r2,v)     
            ! beta = betak+1^2
            beta = this%r2%dot(this%v1)

            beta_lt_zero = (beta < zero)

            call environment%bcast(beta_lt_zero)
            if (beta_lt_zero) then
                istop = 6
                go to 900
            end if

            beta   = sqrt( beta )          ! beta = betak+1
            tnorm2 = tnorm2 + alpha**2 + oldb**2 + beta**2

            if (num_iterations == 1) then                   ! Initialize a few things.
                if (beta/beta1 <= ten*eps) then   ! beta2 = 0 or ~ 0.
                    istop = -1                     ! Terminate later.
                end if
                !tnorm2 = alpha**2
                gmax   = abs( alpha )              ! alpha1
                gmin   = gmax                     ! alpha1
            end if

            ! Apply previous rotation Qk-1 to get
            !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
            !   [gbar k dbar k+1]   [sn -cs][alphak betak+1].

            oldeps = epsln
            delta  = cs * dbar  +  sn * alpha ! delta1 = 0         deltak
            gbar   = sn * dbar  -  cs * alpha ! gbar 1 = alpha1     gbar k
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
                write(*,*) 'alpha ', alpha
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

            call this%w1%copy(this%w2)
            call this%w2%copy(this%w)
            call this%w%copy(this%v2)

            call this%w%axpby(-oldeps, this%w1, 1.0_rp)
            call this%w%axpby(-delta, this%w2, 1.0_rp)
            call this%w%scal(denom, this%w)
            call x%axpby(phi, this%w, 1.0_rp)

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
            epsr   = Anorm * ynorm * rtol + atol
            rhs_convergence_test = epsr
            error_estimate_convergence_test = phibar 
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
                if (num_iterations   >= max_num_iterations    ) istop = 5
                if (Acond  >= 0.1d+0/eps) istop = 4
                if (epsx   >= beta1     ) istop = 3
                ! AFM if (qrnorm <= epsx  .or.  relArnorml <= epsx) istop = 2
                ! AFM if (qrnorm <= epsr  .or.  relArnorml <= epsr) istop = 1
                if (qrnorm <= epsx) istop = 2
                if (qrnorm <= epsr) istop = 1
            end if

             if (track_convergence_history) then 
                error_estimate_history_convergence_test(num_iterations) = error_estimate_convergence_test
             end if
            call this%print_convergence_history_new_line(luout)

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
            call environment%bcast(istop_neq_zero)

            if (istop_neq_zero) exit

        else
            ! y = inv(M) r2
            call M%apply(this%r2, this%y)
            call environment%bcast(beta_lt_zero)

            if (beta_lt_zero) then
                istop = 6
                go to 900
            end if

            call environment%bcast(istop_neq_zero)            

            if (istop_neq_zero) exit

        end if

    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

   ! Check for convergence and output corresponding info. messages
    900 did_converge = .true. 
    if ( .not. (istop==0 .or. istop==-1 .or. istop==1 .or. istop==2) ) then
        did_converge = .false.
    end if

    call environment%bcast(did_converge)
    call this%print_convergence_history_footer(luout)

    call environment%info(me,np)
    if ( me == 0 ) then
        write(luout, 2000) istop, num_iterations,   &
            Anorm, Acond, &
            rnorm, ynorm
        write(luout, 3000) msg(istop)
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

  end subroutine minres_solve_body

  function minres_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(minres_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: minres_supports_stopping_criteria
    minres_supports_stopping_criteria = ( stopping_criteria == res_res )
  end function minres_supports_stopping_criteria
  
  function minres_get_default_stopping_criteria(this)
    implicit none
    class(minres_t), intent(in) :: this
    integer(ip) :: minres_get_default_stopping_criteria
    minres_get_default_stopping_criteria = default_minres_stopping_criteria
  end function minres_get_default_stopping_criteria
  
  
  subroutine create_minres(environment, base_iterative_linear_solver)
    implicit none
    class(environment_t),                           intent(in)    :: environment
    class(base_iterative_linear_solver_t), pointer, intent(inout) :: base_iterative_linear_solver
    type(minres_t),                        pointer                :: minres
    assert(.not. associated(base_iterative_linear_solver))
    allocate(minres)
    call minres%set_environment(environment)
    call minres%set_name(minres_name)
    call minres%set_defaults()
    call minres%set_state(start)
    base_iterative_linear_solver => minres
  end subroutine create_minres
  
end module minres_names
