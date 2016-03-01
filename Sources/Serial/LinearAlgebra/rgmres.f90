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
module rgmres_names
  use types_names
  use stdio_names
  use memor_names
  
  ! Abstract modules
  use vector_names
  use vector_space_names
  use operator_names
  use environment_names
  use base_linear_solver_names  
  use multivector_names

  implicit none
# include "debug.i90"
  private
  
  character(len=*), parameter :: rgmres_name = 'RGMRES'
  character(len=*), parameter :: ls_dkrymax                   = 'linear_solver_dkrymax'
  character(len=*), parameter :: ls_orthonorm_strat           = 'linear_solver_orthonorm_strat'
  character(len=*), parameter :: orthonorm_strat_icgsro       = 'ICGSRO' 
  character(len=*), parameter :: orthonorm_strat_mgsro        = 'MGSRO'
  
  integer (ip), parameter :: mgsro  = 1 ! mgs : Modified Gram-Schmidt 
                                        !       (appropriate for serial GMRES)
  integer (ip), parameter :: icgsro = 2 ! icgs: Iterative Classical Gram-Schmidt 
                                        !       (appropriate for distributed GMRES)
  
  integer (ip)    , parameter :: default_rgmres_stopping_criteria = res_nrmgiven_res_nrmgiven
  integer (ip)    , parameter :: default_dkrymax           = 30
  integer (ip)    , parameter :: default_orthonorm_strat   = icgsro
  
  type, extends(base_linear_solver_t) :: rgmres_t
     ! Parameters
     integer(ip)                    :: dkrymax
     integer(ip)                    :: orthonorm_strat

     ! Working space data members
     type(multivector_t)            :: bkry
     class(vector_t), allocatable   :: r, z  
     real(rp)       , allocatable   :: hh(:,:), g(:), cs(:,:)
   contains
     procedure          :: allocate_workspace            => rgmres_allocate_workspace
     procedure          :: free_workspace                => rgmres_free_workspace
     procedure          :: set_parameters_from_pl        => rgmres_set_parameters_from_pl
     procedure          :: solve_body                    => rgmres_solve_body
     procedure          :: supports_stopping_criteria    => rgmres_supports_stopping_criteria
     procedure          :: get_default_stopping_criteria => rgmres_get_default_stopping_criteria
  end type rgmres_t
  
  ! Data types
  public :: rgmres_t, create_rgmres
  public :: ls_dkrymax, ls_orthonorm_strat, orthonorm_strat_icgsro, orthonorm_strat_mgsro
  public :: default_dkrymax, default_orthonorm_strat
  public :: mgsro, icgsro
  public :: modified_gs_reorthonorm, iterative_gs_reorthonorm, apply_givens_rotation
  
contains
  subroutine rgmres_allocate_workspace(this)
    implicit none
    class(rgmres_t), intent(inout) :: this
    type(vector_space_t), pointer :: range
    type(dynamic_state_operator_t), pointer :: A, M
    class(environment_t), pointer :: environment
    A => this%get_A()
    range  => A%get_range_vector_space()
    call range%create_vector(this%r)
    M => this%get_M()
    range  => M%get_range_vector_space()
    call range%create_vector(this%z)
    environment => this%get_environment()
    call this%bkry%create(environment,this%dkrymax+1,mold=this%z)
    call memalloc(this%dkrymax+1,this%dkrymax+1,this%hh,__FILE__,__LINE__)
    call memalloc(this%dkrymax+1,this%g,__FILE__,__LINE__)
    call memalloc(2,this%dkrymax+1,this%cs,__FILE__,__LINE__)
  end subroutine rgmres_allocate_workspace
  
  subroutine rgmres_free_workspace(this)
    implicit none
    class(rgmres_t), intent(inout) :: this
    call this%r%free()
    call this%z%free()
    call this%bkry%free()
    call memfree(this%hh,__FILE__,__LINE__)
    call memfree(this%g,__FILE__,__LINE__)
    call memfree(this%cs,__FILE__,__LINE__)
  end subroutine rgmres_free_workspace

  subroutine rgmres_set_parameters_from_pl(this) 
   implicit none
   class(rgmres_t), intent(inout) :: this
  end subroutine rgmres_set_parameters_from_pl
  
  subroutine rgmres_solve_body(this,x)
#ifdef ENABLE_BLAS
    use blas77_interfaces_names
#endif
    implicit none
    class(rgmres_t)    , intent(inout) :: this
    class(vector_t)    , intent(inout) :: x 

    ! Local variables to store a copy/reference of the corresponding member variables of base class
    class(environment_t), pointer :: environment
    class(operator_t)   , pointer :: A, M 
    class(vector_t)     , pointer :: initial_solution, b
    integer(ip)                   :: stopping_criteria, max_num_iterations, output_frequency, luout
    real(rp)                      :: atol, rtol
    logical                       :: track_convergence_history

    ! Pointers to freely modify/read private member variables of base class
    integer(ip), pointer :: num_iterations
    logical    , pointer :: did_converge
    real(rp)   , pointer :: rhs_convergence_test, error_estimate_convergence_test
    real(rp)   , pointer :: error_estimate_history_convergence_test(:)
    
    integer(ip)                :: ierrc
    integer(ip)                :: kloc, i, j, k_hh, id
    real(rp)                   :: res_norm, rhs_norm
    real(rp)                   :: alpha, c, s
    integer                    :: me, np
    logical                    :: exit_loop
    class(vector_t), pointer   :: bkryv


    environment               => this%get_environment()
    A                         => this%get_A()
    b                         => this%get_rhs()
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

    call x%copy(initial_solution)
    
    ! Evaluate ||b||_2 if required
    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x,this%r)

    ! r = b-r
    call this%r%axpby(1.0_rp,b,-1.0_rp)

    ! res_norm = ||r||_2
    res_norm = this%r%nrm2()
    
    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
        rhs_convergence_test = rtol * rhs_norm + atol
        error_estimate_convergence_test = res_norm 
    else if ( stopping_criteria == res_nrmgiven_res_nrmgiven ) then
        rhs_convergence_test = rtol * res_norm + atol
        error_estimate_convergence_test = res_norm
    end if
    
    exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
    ! Send converged to coarse-grid tasks
    call environment%bcast(exit_loop)
    
    call this%print_convergence_history_header(luout)
    
    num_iterations = 0
    outer: do while ( (.not.exit_loop) .and. &
         &            (num_iterations < max_num_iterations))

       ! Compute residual from scratch (only if num_iterations /= 0)
       if ( num_iterations /= 0 ) then 
          ! r = Ax
          call A%apply(x,this%r)

          ! r = b-r
          call this%r%axpby(1.0_rp,b,-1.0_rp)

          ! res_norm = ||r||_2
          res_norm = this%r%nrm2()
       end if

       ! Normalize residual direction (i.e., v_1 = r/||r||_2)
       bkryv => this%bkry%get(1)
       call bkryv%clone(x)
       if ( environment%am_i_fine_task() ) then 
          if (res_norm /= 0.0_rp) then
             call bkryv%scal(1.0_rp/res_norm, this%r)
          end if
       end if
       ! residual in the krylov basis
       this%g(1) = res_norm
       this%g(2:this%dkrymax+1) = 0.0_rp

       ! start iterations
       kloc = 0
       inner: do while ( (.not.exit_loop) .and. &
            &            (num_iterations < max_num_iterations) .and. &
            &            (kloc < this%dkrymax))
          kloc  = kloc  + 1
          num_iterations = num_iterations + 1

          ! Generate new basis vector
          bkryv => this%bkry%get(kloc)
          call M%apply( bkryv, this%z )
          bkryv => this%bkry%get(kloc+1)
          call bkryv%clone(x)
          call A%apply(this%z, bkryv)

          if ( environment%am_i_fine_task() ) then
             ! Orthogonalize
             select case( this%orthonorm_strat )
             case ( mgsro )
                call modified_gs_reorthonorm  (luout, kloc+1, this%bkry, this%hh(1,kloc), ierrc )
             case ( icgsro )
                call iterative_gs_reorthonorm (luout, kloc+1, this%bkry, this%hh(1,kloc), ierrc )
             case default
                check(.false.)
             end select

             if ( ierrc < 0 ) then
                ! The coarse-grid task should exit 
                ! the inner-do loop. Send signal.
                exit_loop = .true.
                call environment%bcast(exit_loop) 
                exit inner ! Exit inner do-loop
             end if

             ! Apply previous given's rotations to kth column of hessenberg matrix
             k_hh = 1
             do j = 1,kloc-1    
                alpha = this%hh(k_hh,kloc)
                c = this%cs(1,j)
                s = this%cs(2,j)
                this%hh(k_hh,kloc) = c*alpha + s*this%hh(k_hh+1,kloc)
                this%hh(k_hh+1,kloc) = c*this%hh(k_hh+1,kloc) - s*alpha
                k_hh = k_hh +1
             enddo

             ! Compute (and apply) new given's rotation
             call apply_givens_rotation(this%hh(k_hh,kloc), this%hh(k_hh+1,kloc), c, s)
             this%cs(1,kloc) = c
             this%cs(2,kloc) = s

             ! Update residual vector 
             ! (expressed in terms of the Krylov basis)
             alpha = -s*this%g(kloc)
             this%g(kloc) = c*this%g(kloc)
             this%g(kloc+1) = alpha

             ! Error norm
             res_norm = abs(alpha)
             error_estimate_convergence_test  = res_norm
             if (track_convergence_history) then 
                error_estimate_history_convergence_test(num_iterations) = error_estimate_convergence_test
             end if
          end if
          exit_loop = (error_estimate_convergence_test <= rhs_convergence_test) 
          ! Send converged to coarse-grid tasks
          call environment%bcast(exit_loop)

          call this%print_convergence_history_new_line(luout)
       end do inner

       if ( kloc > 0 ) then
          if ( environment%am_i_fine_task() ) then
             if ( ierrc == -2 ) then
                write (luout,*) '** Warning: RGMRES: ortho failed due to abnormal numbers, no way to proceed'
                ! The coarse-grid task should exit 
                ! the outer-do loop. Send signal. 
                exit_loop = .true.
                call environment%bcast(exit_loop)
                exit outer ! Exit main do-loop
             end if

             ! Compute the solution
             ! If zero on the diagonal, 
             ! solve a reduced linear system
             do while ( kloc > 0 ) ! .and. hh(kloc,kloc) == 0.0_rp  )
                if(this%hh(kloc,kloc) /= 0.0_rp) exit
                kloc = kloc - 1
             end do

             if ( kloc <= 0 ) then
                write (luout,*) '** Warning: RGMRES: triangular system in GMRES has null rank'
                exit_loop = .true.
                call environment%bcast(exit_loop)
                exit outer ! Exit main do loop     
             end if

#ifdef ENABLE_BLAS       
             !N    !A  !LDA        !X !INCX
             call DTRSV ( 'U', 'N', 'N', kloc, this%hh, this%dkrymax+1, this%g, 1)
#else
             ! Solve the system hh*y = g
             ! Solution stored on g itself
             do j = kloc,1,-1
                this%g(j) = this%g(j)/this%hh(j,j)
                do i = j-1,1,-1
                   this%g(i) = this%g(i) - this%hh(i,j) * this%g(j)
                end do
             end do
#endif       
          end if

          ! Now g contains the solution in the krylov basis
          ! Compute the solution in the real space
          call this%r%init(0.0_rp)

          ! r <- g_1 * v_1 + g_2 * v_2 + ... + g_kloc * v_kloc
          call this%bkry%multiaxpy(kloc, this%r, 1.0_rp, this%g)

          ! Solve Mz = r
          call M%apply(this%r,this%z)

          ! x <- x + z
          call x%axpby(1.0_rp,this%z,1.0_rp)

          exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
          ! Send converged to coarse-grid tasks
          call environment%bcast(exit_loop)
       end if
    end do outer
    did_converge = (error_estimate_convergence_test <= rhs_convergence_test)
    ! Send converged to coarse-grid tasks
    call environment%bcast(did_converge)
    call this%print_convergence_history_footer(luout)
  end subroutine rgmres_solve_body

  function rgmres_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(rgmres_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: rgmres_supports_stopping_criteria
    rgmres_supports_stopping_criteria = ( stopping_criteria == res_nrmgiven_rhs_nrmgiven .or. &
                                          stopping_criteria == res_nrmgiven_res_nrmgiven )
  end function rgmres_supports_stopping_criteria
  
  function rgmres_get_default_stopping_criteria(this)
    implicit none
    class(rgmres_t), intent(in) :: this
    integer(ip) :: rgmres_get_default_stopping_criteria
    rgmres_get_default_stopping_criteria = default_rgmres_stopping_criteria
  end function rgmres_get_default_stopping_criteria
  
  function create_rgmres(environment)
    implicit none
    class(environment_t), intent(in) :: environment
    class(base_linear_solver_t), pointer :: create_rgmres
    type(rgmres_t), pointer :: rgmres
    allocate(rgmres)
    call rgmres%set_environment(environment)
    call rgmres%set_name(rgmres_name)
    call rgmres%set_defaults()
    rgmres%dkrymax = default_dkrymax
    rgmres%orthonorm_strat = default_orthonorm_strat
    call rgmres%set_state(start)
    create_rgmres => rgmres
  end function create_rgmres
  
  
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
  
end module rgmres_names
