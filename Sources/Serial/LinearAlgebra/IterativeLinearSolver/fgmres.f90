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
module fgmres_names
  use types_names
  use stdio_names
  use memor_names
  
  ! Abstract modules
  use vector_names
  use vector_space_names
  use operator_names
  use environment_names
  use base_iterative_linear_solver_names  
  use iterative_linear_solver_utils_names
  use iterative_linear_solver_parameters_names
  use multivector_names

  implicit none
# include "debug.i90"
  private
  
  type, extends(base_iterative_linear_solver_t) :: fgmres_t
     ! Parameters
     integer(ip)                    :: dkrymax
     integer(ip)                    :: orthonorm_strat

     ! Working space data members
     type(multivector_t)            :: bkry  ! Krylov basis
     type(multivector_t)            :: bkryz ! (Right-)Preconditioned Krylov basis
     class(vector_t), allocatable   :: r, z  
     real(rp)       , allocatable   :: hh(:,:), g(:), cs(:,:)
   contains
     procedure          :: allocate_workspace            => fgmres_allocate_workspace
     procedure          :: free_workspace                => fgmres_free_workspace
     procedure          :: set_parameters_from_pl        => fgmres_set_parameters_from_pl
     procedure          :: solve_body                    => fgmres_solve_body
     procedure          :: supports_stopping_criteria    => fgmres_supports_stopping_criteria
     procedure          :: get_default_stopping_criteria => fgmres_get_default_stopping_criteria
  end type fgmres_t
  
  public :: create_fgmres
  
contains
  subroutine fgmres_allocate_workspace(this)
    implicit none
    class(fgmres_t), intent(inout) :: this
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
    call this%bkryz%create(environment,this%dkrymax+1,mold=this%z)
    call memalloc(this%dkrymax+1,this%dkrymax+1,this%hh,__FILE__,__LINE__)
    call memalloc(this%dkrymax+1,this%g,__FILE__,__LINE__)
    call memalloc(2,this%dkrymax+1,this%cs,__FILE__,__LINE__)
  end subroutine fgmres_allocate_workspace
  
  subroutine fgmres_free_workspace(this)
    implicit none
    class(fgmres_t), intent(inout) :: this
    call this%r%free()
    call this%z%free()
    call this%bkry%free()
    call this%bkryz%free()
    call memfree(this%hh,__FILE__,__LINE__)
    call memfree(this%g,__FILE__,__LINE__)
    call memfree(this%cs,__FILE__,__LINE__)
  end subroutine fgmres_free_workspace

  subroutine fgmres_set_parameters_from_pl(this) 
   implicit none
   class(fgmres_t), intent(inout) :: this
  end subroutine fgmres_set_parameters_from_pl
  
  subroutine fgmres_solve_body(this,b,x)
#ifdef ENABLE_BLAS
    use blas77_interfaces_names
#endif
    implicit none
    class(fgmres_t)    , intent(inout) :: this
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
    
    integer(ip)                :: ierrc
    integer(ip)                :: kloc, i, j, k_hh, id
    real(rp)                   :: res_norm, rhs_norm
    real(rp)                   :: alpha, c, s
    integer                    :: me, np
    logical                    :: exit_loop
    class(vector_t), pointer   :: bkryv, bkryzv

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

    call x%copy(initial_solution)
    
     ! Evaluate ||b||_2 if required
    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x, this%r)

    ! r = b-r
    call this%r%axpby(1.0_rp,b,-1.0_rp)

    ! res_norm = ||r||_2
    res_norm = this%r%nrm2()

    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
        rhs_convergence_test  = rtol * rhs_norm + atol
        error_estimate_convergence_test = res_norm 
    else if ( stopping_criteria == res_nrmgiven_res_nrmgiven ) then
        rhs_convergence_test  = rtol * res_norm + atol
        error_estimate_convergence_test = res_norm
    end if
    exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
    ! Send converged to coarse-grid tasks
    call environment%bcast(exit_loop)

    call this%print_convergence_history_header(luout)

    num_iterations = 0
    outer: do while ( (.not.exit_loop ) .and. &
        &            (num_iterations < max_num_iterations))
        ! Compute residual from scratch (only if num_iterations /= 0)
        if ( num_iterations /= 0 ) then 
            ! r = Ax
            call A%apply(x, this%r)

            ! r = b-r
            call this%r%axpby(1.0_rp,b,-1.0_rp)

            ! res_norm = ||r||_2
            res_norm = this%r%nrm2()
        end if

        ! Normalize residual direction (i.e., v_1 = r/||r||_2)
        bkryv => this%bkry%get(1)
        call bkryv%clone(x)
        if ( environment%am_i_fine_task() ) then 
        if (res_norm /= 0.0_rp) call bkryv%scal(1.0_rp/res_norm, this%r)
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
            bkryv  => this%bkry%get(kloc)
            bkryzv => this%bkryz%get(kloc)
            call bkryzv%clone(x)
            call M%apply(bkryv,bkryzv)
            bkryv  => this%bkry%get(kloc+1)
            call bkryv%clone(x)
            call A%apply(bkryzv,bkryv)
        
            if ( environment%am_i_fine_task() ) then
                ! Orthogonalize
                select case( this%orthonorm_strat )
                    case ( mgsro )
                        call modified_gs_reorthonorm  (luout,kloc+1, this%bkry, this%hh(1,kloc), ierrc )
                    case ( icgsro )
                        call iterative_gs_reorthonorm (luout,kloc+1, this%bkry, this%hh(1,kloc), ierrc )
                    case default
                        ! Write an error message and stop ?      
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
                    write (luout,*) '** Warning: FGMRES: ortho failed due to abnormal numbers, no way to proceed'
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
                    write (luout,*) '** Warning: FGMRES: triangular system in FGMRES has null rank'
                    exit_loop = .true.
                    call environment%bcast(exit_loop)
                    exit outer    
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
            call this%z%init(0.0_rp)

            ! z <- g_1 * M_1^-1v_1 + g_2 *  M_2^-1v_2  + ... + g_kloc *  M_kloc^-1v_kloc
            call this%bkryz%multiaxpy(kloc, this%z, 1.0_rp, this%g )
            
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
  end subroutine fgmres_solve_body

  function fgmres_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(fgmres_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: fgmres_supports_stopping_criteria
    fgmres_supports_stopping_criteria = ( stopping_criteria == res_nrmgiven_rhs_nrmgiven .or. &
                                          stopping_criteria == res_nrmgiven_res_nrmgiven )
  end function fgmres_supports_stopping_criteria
  
  function fgmres_get_default_stopping_criteria(this)
    implicit none
    class(fgmres_t), intent(in) :: this
    integer(ip) :: fgmres_get_default_stopping_criteria
    fgmres_get_default_stopping_criteria = default_fgmres_stopping_criteria
  end function fgmres_get_default_stopping_criteria
  
  function create_fgmres(environment)
    implicit none
    class(environment_t), intent(in) :: environment
    class(base_iterative_linear_solver_t), pointer :: create_fgmres
    type(fgmres_t), pointer :: fgmres
    allocate(fgmres)
    call fgmres%set_environment(environment)
    call fgmres%set_name(fgmres_name)
    call fgmres%set_defaults()
    fgmres%dkrymax = default_dkrymax
    fgmres%orthonorm_strat = default_orthonorm_strat
    call fgmres%set_state(start)
    create_fgmres => fgmres
  end function create_fgmres
  
end module fgmres_names
