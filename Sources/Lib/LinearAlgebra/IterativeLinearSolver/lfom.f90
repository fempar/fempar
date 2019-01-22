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
module lfom_names
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
  use rgmres_names
  use FPL

  implicit none
# include "debug.i90"
  private
  
  type, extends(base_iterative_linear_solver_t) :: lfom_t
     ! Parameters
     integer(ip)                    :: dkrymax
     integer(ip)                    :: orthonorm_strat

     ! Working space data members
     type(multivector_t)            :: bkry
     class(vector_t), allocatable   :: r, z  
     real(rp)       , allocatable   :: hh(:,:), lu(:,:), g(:)
     integer(ip)    , allocatable   :: ipiv(:)
   contains
     procedure          :: allocate_workspace            => lfom_allocate_workspace
     procedure          :: free_workspace                => lfom_free_workspace
     procedure          :: set_parameters_from_pl        => lfom_set_parameters_from_pl
     procedure          :: solve_body                    => lfom_solve_body
     procedure          :: supports_stopping_criteria    => lfom_supports_stopping_criteria
     procedure          :: get_default_stopping_criteria => lfom_get_default_stopping_criteria
  end type lfom_t
  
  public :: create_lfom
  
contains
  subroutine lfom_allocate_workspace(this)
    implicit none
    class(lfom_t), intent(inout) :: this
    type(vector_space_t), pointer :: range
    type(lvalue_operator_t), pointer :: A, M
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
    call memalloc(this%dkrymax+1,this%dkrymax+1,this%lu,__FILE__,__LINE__)
    call memalloc(this%dkrymax+1,this%g,    __FILE__,__LINE__)
    call memalloc(this%dkrymax+1,this%ipiv, __FILE__,__LINE__)
  end subroutine lfom_allocate_workspace
  
  subroutine lfom_free_workspace(this)
    implicit none
    class(lfom_t), intent(inout) :: this
    call this%r%free()
    call this%z%free()
    call this%bkry%free()
    call memfree(this%hh,__FILE__,__LINE__)
    call memfree(this%lu,__FILE__,__LINE__)
    call memfree(this%g,__FILE__,__LINE__)
    call memfree(this%ipiv,__FILE__,__LINE__)
  end subroutine lfom_free_workspace

  subroutine lfom_set_parameters_from_pl(this, parameter_list) 
   implicit none
   class(lfom_t),         intent(inout) :: this
   type(ParameterList_t), intent(in)    :: parameter_list
   integer(ip)                          :: FPLError
   call this%base_iterative_linear_solver_set_parameters_from_pl(parameter_list)
   ! Dkrymax
   if(parameter_list%isPresent(ils_max_dim_krylov_basis_key)) then
       assert(parameter_list%isAssignable(ils_max_dim_krylov_basis_key, this%dkrymax))
       FPLError   = parameter_list%Get(Key=ils_max_dim_krylov_basis_key, Value=this%dkrymax)
       assert(FPLError == 0)
   endif
   ! Orthonorm strat
   if(parameter_list%isPresent(ils_orthonorm_strategy_key)) then
       assert(parameter_list%isAssignable(ils_orthonorm_strategy_key, this%orthonorm_strat))
       FPLError   = parameter_list%Get(Key=ils_orthonorm_strategy_key, Value=this%orthonorm_strat)
       assert(FPLError == 0)
   endif
  end subroutine lfom_set_parameters_from_pl
  
  subroutine lfom_solve_body(this,b,x)
#ifdef ENABLE_LAPACK
    use lapack77_interfaces_names
#endif
    implicit none
    class(lfom_t)      , intent(inout) :: this
    class(vector_t)    , intent(in)    :: b
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
    real(rp)                   :: res_norm, res_2_norm, rhs_norm, res_norm_initial
    real(rp)                   :: err1_it, ub1
    integer                    :: me, np, info
    logical                    :: exit_loop
    class(vector_t), pointer   :: bkryv


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
    if ( stopping_criteria == res_rhs ) then
        rhs_norm = b%nrm2()
    endif

    ! r = Ax
    call A%apply(x, this%r)

    ! r = b-r
    call this%r%axpby(1.0_rp,b,-1.0_rp)

    res_2_norm = this%r%nrm2()

    ! z=inv(M)r
    call M%apply(this%r,this%z)

    ! Evaluate ||z||_2
    res_norm = this%z%nrm2()
    res_norm_initial = res_norm

    if ( stopping_criteria == res_res ) then
        rhs_convergence_test = rtol * res_2_norm + atol
        error_estimate_convergence_test = res_2_norm
    else if ( stopping_criteria == res_rhs ) then
        rhs_convergence_test = rtol * rhs_norm + atol
        error_estimate_convergence_test = res_2_norm
    end if

    exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
    ! Send converged to coarse-grid tasks
    call environment%l1_lgt1_bcast(exit_loop)

    call this%print_convergence_history_header(luout)

    num_iterations = 0
    did_converge = .false.
    outer: do while ( (.not.exit_loop).and.(num_iterations<max_num_iterations) )
        this%hh = 0.0_rp
        ! Compute preconditioned residual from scratch (only if num_iterations/=0)
        if ( num_iterations /= 0 ) then 
            ! r = Ax
            call A%apply(x, this%r)

            ! r = b-r
            call this%r%axpby(1.0_rp,b,-1.0_rp)

            ! z=inv(M)r
            call M%apply(this%r, this%z)

            ! Evaluate ||z||_2
            res_norm = this%z%nrm2()
        end if

        ! Normalize preconditioned residual direction (i.e., v_1 = z/||z||_2)
        bkryv => this%bkry%get(1)
        call bkryv%clone(x)
        if ( environment%am_i_l1_task() ) then 
          if (res_norm /= 0.0_rp) call bkryv%scal(1.0_rp/res_norm, this%z)
        end if
        ! start iterations
        kloc = 0
        inner: do while ( (.not.exit_loop) .and. &
            &            (num_iterations < max_num_iterations) .and. &
            &            (kloc < this%dkrymax))
            kloc  = kloc  + 1
            num_iterations = num_iterations + 1

            ! Generate new basis vector
            bkryv => this%bkry%get(kloc)
            call A%apply(bkryv,this%r)
            bkryv => this%bkry%get(kloc+1)
            call bkryv%clone(x)
            call M%apply(this%r, bkryv)
            
            if ( environment%am_i_l1_task() ) then ! Am I a fine task ?
                ! Orthogonalize
                select case( this%orthonorm_strat )
                    case ( mgsro )
                        call modified_gs_reorthonorm  ( luout, kloc+1, this%bkry, this%hh(1,kloc), ierrc )
                    case ( icgsro )
                        call iterative_gs_reorthonorm ( luout, kloc+1, this%bkry, this%hh(1,kloc), ierrc )
                    case default
                        check(.false.)
                        ! Write an error message and stop ?      
                end select

                if ( ierrc < 0 ) then
                    ! The coarse-grid task should exit 
                    ! the inner-do loop. Send signal.
                    exit_loop = .true.
                    call environment%l1_lgt1_bcast(exit_loop)
                    exit inner ! Exit inner do-loop
                end if
    
                ! init right-hand-size to \beta*e_1, with \beta=||z_0||_2
                this%g(1)            = res_norm_initial
                this%g(2:kloc)       = 0.0_rp
    
                if ( kloc > 0 ) then
                    ! Compute the solution
    
                    this%lu(1:kloc,1:kloc) = this%hh(1:kloc,1:kloc)
                    ! write(*,*) 'XXX', lu(1:kloc,1:kloc)
#ifdef ENABLE_LAPACK
                    call dgetrf( kloc, kloc, this%lu, this%dkrymax+1, this%ipiv, info )
                    if ( info /= 0 ) then
                        write (luout,*) '** Warning: LFOM: dgetrf returned info /= 0'
                        exit_loop = .true.
                        call environment%l1_lgt1_bcast(exit_loop)
                        call environment%l1_lgt1_bcast(exit_loop)
                        exit outer ! Exit main do loop 
                    end if
    
                    call dgetrs( 'N' , kloc, 1, this%lu, this%dkrymax+1, this%ipiv, this%g, this%dkrymax+1, info )
                    if ( info /= 0 ) then
                        write (luout,*) '** Warning: LFOM: dgetrs returned info /= 0'
                        exit_loop = .true.
                        call environment%l1_lgt1_bcast(exit_loop)
                        call environment%l1_lgt1_bcast(exit_loop)
                        exit outer ! Exit main do loop 
                    end if
#else
                    write (0,*) 'plfom ERROR :: dgetrf and dgetrs not available'
                    check(1==0)
#endif

                    ! Now g contains the solution in the krylov basis
                    ! Compute the solution in the global space
                    call this%z%copy(x)

                    ! z <-z +  g_1 * v_1 + g_2 * v_2 + ... + g_kloc * v_kloc
                    call this%bkry%multiaxpy(kloc, this%z, 1.0_rp, this%g )
    
                    ! r = A(z+x)
                    call A%apply(this%z, this%r)
    
                    ! r = b-r
                    call this%r%axpby(1.0_rp,b,-1.0_rp)
    
                    res_2_norm = this%r%nrm2()
                end if
                error_estimate_convergence_test = res_2_norm
            end if

            if (track_convergence_history) then  
            error_estimate_history_convergence_test(num_iterations) = error_estimate_convergence_test
            end if 
            exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
            ! Send converged to coarse-grid tasks
            call environment%l1_lgt1_bcast(exit_loop)
    
            call this%print_convergence_history_new_line(luout)
        end do inner

        if ( environment%am_i_l1_task() ) then ! Am I a fine task ?
            if ( ierrc == -2 ) then
                write (luout,*) '** Warning: LFOM: ortho failed due to abnormal numbers, no way to proceed'
                ! The coarse-grid task should exit 
                ! the outer-do loop. Send signal. 
                exit_loop = .true.
                call environment%l1_lgt1_bcast(exit_loop)
                exit outer ! Exit outer do-loop
            end if
        end if

        ! x <- z
        call x%copy(this%z)

        exit_loop = (error_estimate_convergence_test <= rhs_convergence_test)
        ! Send converged to coarse-grid tasks
        call environment%l1_lgt1_bcast(exit_loop)

    end do outer

    did_converge = (error_estimate_convergence_test <= rhs_convergence_test)
    ! Send converged to coarse-grid tasks
    call environment%l1_lgt1_bcast(did_converge)
    call this%print_convergence_history_footer(luout)
  end subroutine lfom_solve_body

  function lfom_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(lfom_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: lfom_supports_stopping_criteria
    lfom_supports_stopping_criteria = ( stopping_criteria == res_res .or. &
                                          stopping_criteria == res_rhs )
  end function lfom_supports_stopping_criteria
  
  function lfom_get_default_stopping_criteria(this)
    implicit none
    class(lfom_t), intent(in) :: this
    integer(ip) :: lfom_get_default_stopping_criteria
    lfom_get_default_stopping_criteria = default_lfom_stopping_criteria
  end function lfom_get_default_stopping_criteria
  
  subroutine create_lfom(environment, base_iterative_linear_solver)
    implicit none
    class(environment_t),                           intent(in)    :: environment
    class(base_iterative_linear_solver_t), pointer, intent(inout) :: base_iterative_linear_solver
    type(lfom_t),                          pointer                :: lfom
    assert(.not. associated(base_iterative_linear_solver))
    allocate(lfom)
    call lfom%set_environment(environment)
    call lfom%set_name(lfom_name)
    call lfom%set_defaults()
    lfom%dkrymax = default_dkrymax
    lfom%orthonorm_strat = default_orthonorm_strat
    call lfom%set_state(start)
    base_iterative_linear_solver => lfom
  end subroutine create_lfom
  
end module lfom_names
