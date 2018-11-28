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
module cg_names
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
  use FPL

  implicit none
# include "debug.i90"
  private

  type, extends(base_iterative_linear_solver_t) :: cg_t
    ! Working space vectors for type(cg_t)
    class(vector_t), allocatable :: r
    class(vector_t), allocatable :: z 
    class(vector_t), allocatable :: Ap 
    class(vector_t), allocatable :: p 
  contains
    procedure          :: allocate_workspace            => cg_allocate_workspace
    procedure          :: free_workspace                => cg_free_workspace
    procedure          :: set_parameters_from_pl        => cg_set_parameters_from_pl
    procedure          :: solve_body                    => cg_solve_body
    procedure          :: supports_stopping_criteria    => cg_supports_stopping_criteria
    procedure          :: get_default_stopping_criteria => cg_get_default_stopping_criteria
    procedure,private  :: init_convergence_data
    procedure,private  :: update_convergence_data_and_evaluate_stopping_criteria 
  end type
  
  public :: create_cg
  
contains
  subroutine cg_allocate_workspace(this)
    implicit none
    class(cg_t), intent(inout) :: this
    type(vector_space_t), pointer :: range
    type(lvalue_operator_t), pointer :: A, M
    A => this%get_A()
    range  => A%get_range_vector_space()
    call range%create_vector(this%r)
    call range%create_vector(this%Ap)
    M => this%get_M()
    range  => M%get_range_vector_space()
    call range%create_vector(this%z)
    call range%create_vector(this%p)
  end subroutine cg_allocate_workspace
  
  subroutine cg_free_workspace(this)
    implicit none
    class(cg_t), intent(inout) :: this
    call this%r%free()
    call this%Ap%free()
    call this%z%free()
    call this%p%free()
    deallocate(this%r)
    deallocate(this%Ap)
    deallocate(this%z)
    deallocate(this%p)
  end subroutine cg_free_workspace

  subroutine cg_set_parameters_from_pl(this, parameter_list) 
   implicit none
   class(cg_t),           intent(inout) :: this
   type(ParameterList_t), intent(in)    :: parameter_list
   call this%base_iterative_linear_solver_set_parameters_from_pl(parameter_list)
  end subroutine cg_set_parameters_from_pl
  
  subroutine cg_solve_body(this,b,x)
    implicit none
    class(cg_t), intent(inout) :: this
    class(vector_t)    , intent(in) :: b 
    class(vector_t)    , intent(inout) :: x 

    real(rp) :: r_nrm_M      ! |r|_inv(M) 
    real(rp) :: b_nrm_M      ! |b|_inv(M)
    real(rp) :: r_z, Ap_p, alpha, beta
 
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

    ! Evaluate |b|_inv(M) if required
    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
       ! r = inv(M) b
       call M%apply(b, this%r)
       b_nrm_M = b%dot(this%r)
       if ( environment%am_i_l1_task() ) then ! Am I a fine task ?
          b_nrm_M = sqrt(b_nrm_M)
       end if
    else
       b_nrm_M = 0.0_rp
    endif

    ! 1) Compute initial residual
    ! 1.a) r=Ax
    call A%apply(x,this%r)

    ! 1.b) r=b-r
    call this%r%axpby(1.0, b ,-1.0)
    
    ! 2) z=inv(M)r
    call M%apply(this%r,this%z)

    ! 3) <r,z>
    r_z = this%r%dot(this%z)

    if ( environment%am_i_l1_task() ) then ! Am I a fine task ?
       wassert(r_z>=0.0_rp, 'Square root of a negative number!')
       r_nrm_M = sqrt( r_z )
    end if

    ! 4) Initializations:
    ! p=z
    call this%p%copy(this%z)

    ! Init and log convergence 
    call this%init_convergence_data(b, this%r, b_nrm_M, r_nrm_M)
    call this%print_convergence_history_header(luout)

    ! 5) Iteration
    num_iterations = 0
    did_converge = .false.
    loop_cg: do
       num_iterations = num_iterations + 1

       ! Ap = A*p
       call A%apply(this%p,this%Ap)

       ! <Ap,p>
       Ap_p = this%Ap%dot(this%p)

       ! Is this correct/appropriate ?
       if (Ap_p /= 0.0_rp) then
          alpha = r_z / Ap_p
       else
          alpha = 0.0_rp
       end if

       ! x = x + alpha*p
       call x%axpby(alpha, this%p, 1.0_rp)

       ! r = r - alpha*Ap
       call this%r%axpby (-alpha, this%Ap, 1.0_rp)

       ! Check and log convergence
       call this%update_convergence_data_and_evaluate_stopping_criteria(this%r, r_nrm_M, alpha, this%p)
       call this%print_convergence_history_new_line(luout)
       
       ! Send did_converge to coarse-grid tasks
       call environment%l1_lgt1_bcast(did_converge)

       if(did_converge .or.(num_iterations>=max_num_iterations)) exit loop_cg

       ! z = inv(M) r
       call M%apply(this%r,this%z)

       if ( environment%am_i_l1_task()) then ! Am I a fine task ?
          beta = 1.0_rp/r_z
       else
          beta = 0.0_rp
       end if

       r_z=this%r%dot(this%z)
       beta=beta*r_z

       if (environment%am_i_l1_task()) then ! Am I a fine task ?
          wassert(r_z>=0.0, 'Square root of a negative number!')
          r_nrm_M = sqrt( r_z )
       end if

       ! p = z + beta*p
       call this%p%axpby(1.0,this%z,beta)
    end do loop_cg

    call this%print_convergence_history_footer(luout)
  end subroutine cg_solve_body

  function cg_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(cg_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: cg_supports_stopping_criteria
    cg_supports_stopping_criteria = (stopping_criteria == res_res .or. & 
                                     stopping_criteria == res_rhs .or. &
                                     stopping_criteria == res_nrmgiven_rhs_nrmgiven .or. &
                                     stopping_criteria == res_nrmgiven_res_nrmgiven .or. &
                                     stopping_criteria == delta_rhs .or. &
                                     stopping_criteria == delta_delta .or. &
                                     stopping_criteria == delta_rhs_and_res_res .or. &
                                     stopping_criteria == delta_rhs_and_res_rhs .or. &
                                     stopping_criteria == delta_delta_and_res_res .or. &
                                     stopping_criteria == delta_delta_and_res_rhs )
  end function cg_supports_stopping_criteria
  
  function cg_get_default_stopping_criteria(this)
    implicit none
    class(cg_t), intent(in) :: this
    integer(ip) :: cg_get_default_stopping_criteria
    cg_get_default_stopping_criteria = default_cg_stopping_criteria
  end function cg_get_default_stopping_criteria
  
  
  subroutine create_cg(environment, base_iterative_linear_solver)
    implicit none
    class(environment_t),                           intent(in)    :: environment
    class(base_iterative_linear_solver_t), pointer, intent(inout) :: base_iterative_linear_solver
    type(cg_t),                            pointer                :: cg
    assert(.not. associated(base_iterative_linear_solver))
    allocate(cg)
    call cg%set_environment(environment)
    call cg%set_name(cg_name)
    call cg%set_defaults()
    call cg%set_state(start)
    base_iterative_linear_solver => cg
  end subroutine create_cg

  subroutine init_convergence_data ( this, b, r, nrm_b_given, nrm_r_given )
    implicit none
    class(cg_t)     , intent(inout) :: this
    class(vector_t) , intent(in)    :: b, r 
    real(rp)        , intent(in)    :: nrm_b_given, nrm_r_given
    real(rp)                        :: atol, rtol
    integer(ip)                     :: stopping_criteria
    real(rp)            , pointer   :: rhs_convergence_test
    real(rp)            , pointer   :: rhs_extra_convergence_test
    class(environment_t), pointer   :: environment

    environment => this%get_environment()
    if ( environment%am_i_l1_task() ) then
       atol                        = this%get_atol()
       rtol                        = this%get_rtol()
       stopping_criteria           = this%get_stopping_criteria()
       rhs_convergence_test       => this%get_pointer_rhs_convergence_test()
       rhs_extra_convergence_test => this%get_pointer_rhs_extra_convergence_test()

       select case(stopping_criteria)
       case ( delta_rhs, res_rhs, delta_rhs_and_res_rhs, delta_delta_and_res_rhs )
          rhs_convergence_test = rtol * b%nrm2() + atol
          if (stopping_criteria == delta_rhs_and_res_rhs .or. &
             stopping_criteria == delta_delta_and_res_rhs ) then      
             rhs_extra_convergence_test = rhs_convergence_test 
          end if
       case ( res_res, delta_delta_and_res_res )
          rhs_convergence_test = rtol * r%nrm2() + atol
          if (stopping_criteria == delta_delta_and_res_res ) then
             rhs_extra_convergence_test = rhs_convergence_test 
          end if
       case ( delta_rhs_and_res_res )
          rhs_convergence_test       = rtol*b%nrm2() + atol
          rhs_extra_convergence_test = rtol*r%nrm2() + atol
       case ( res_nrmgiven_rhs_nrmgiven )
          rhs_convergence_test = rtol* nrm_b_given + atol
       case ( res_nrmgiven_res_nrmgiven )
          rhs_convergence_test = rtol* nrm_r_given + atol
       case default
          ! Write an error message and stop ?      
       end select
    end if
  end subroutine init_convergence_data

  subroutine update_convergence_data_and_evaluate_stopping_criteria ( this, r, nrm_r_given, alpha, p )
    implicit none
    class(cg_t)     , intent(inout) :: this
    class(vector_t) , intent(in)    :: r, p 
    real(rp)        , intent(in)    :: nrm_r_given
    real(rp)        , intent(in)    :: alpha
    real(rp)                        :: atol, rtol
    integer(ip)                     :: stopping_criteria
    logical                         :: track_convergence_history
    logical             , pointer   :: did_converge
    integer(ip)         , pointer   :: num_iterations
    real(rp)            , pointer   :: rhs_convergence_test
    real(rp)            , pointer   :: rhs_extra_convergence_test
    real(rp)            , pointer   :: error_estimate_convergence_test
    real(rp)            , pointer   :: error_estimate_history_convergence_test(:)
    real(rp)            , pointer   :: error_estimate_extra_convergence_test
    real(rp)            , pointer   :: error_estimate_history_extra_convergence_test(:)
    class(environment_t), pointer   :: environment

    environment => this%get_environment()
    if ( environment%am_i_l1_task() ) then
       atol                                     = this%get_atol()
       rtol                                     = this%get_rtol()
       track_convergence_history                = this%get_track_convergence_history()
       stopping_criteria                        = this%get_stopping_criteria()
       did_converge                             => this%get_pointer_did_converge()
       num_iterations                           => this%get_pointer_num_iterations()
       rhs_convergence_test                     => this%get_pointer_rhs_convergence_test()
       rhs_extra_convergence_test               => this%get_pointer_rhs_extra_convergence_test()
       error_estimate_convergence_test          => this%get_pointer_error_estimate_convergence_test()
       error_estimate_history_convergence_test  => this%get_pointer_error_estimate_history_convergence_test()
       error_estimate_extra_convergence_test    => this%get_pointer_error_estimate_extra_convergence_test()
       error_estimate_history_extra_convergence_test => this%get_pointer_error_estimate_history_extra_convergence_test()

       ! Compute 1st iteration error estimate and upper bound 
       ! for convergence criteria depending on ||dx(i)||
       if ( num_iterations == 1 ) then
          select case( stopping_criteria )
          case ( delta_rhs, delta_delta,  delta_rhs_and_res_res, delta_rhs_and_res_rhs, &
               & delta_delta_and_res_res, delta_delta_and_res_rhs )
             error_estimate_convergence_test = alpha*p%nrm2()
             select case( stopping_criteria )
             case ( delta_delta, delta_delta_and_res_res, delta_delta_and_res_rhs )
                rhs_convergence_test = rtol*error_estimate_convergence_test + atol
             end select
          end select
       end if

       did_converge = .false.

       ! Evaluate 1st convergence criterion
       select case( stopping_criteria )
       case ( res_res, res_rhs )
          ! Compute || r(i) ||
          error_estimate_convergence_test = r%nrm2()
       case ( delta_rhs, delta_delta, delta_rhs_and_res_res, delta_rhs_and_res_rhs, &
            & delta_delta_and_res_res,delta_delta_and_res_rhs )
          if ( num_iterations /= 1 ) then ! if false, no need to evaluate ||dx(i)|| again
             ! Compute || dx(i) ||
             error_estimate_convergence_test = alpha*p%nrm2()
          end if
       case ( res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven )
          error_estimate_convergence_test = nrm_r_given
       end select

       if ( track_convergence_history ) then
          error_estimate_history_convergence_test(num_iterations) = error_estimate_convergence_test
       end if
       did_converge = (error_estimate_convergence_test <= rhs_convergence_test)

       ! Evaluate 2nd convergence criterion
       select case( stopping_criteria )
       case ( delta_rhs_and_res_res,   delta_rhs_and_res_rhs, &
            & delta_delta_and_res_res, delta_delta_and_res_rhs )
          if ( did_converge ) then
             ! Compute || r(i) ||
             error_estimate_extra_convergence_test = r%nrm2()
             did_converge = (error_estimate_extra_convergence_test <= rhs_extra_convergence_test )
             if ( track_convergence_history ) then
                error_estimate_history_extra_convergence_test(num_iterations) = error_estimate_extra_convergence_test
             end if
          end if
       end select
    end if
  end subroutine update_convergence_data_and_evaluate_stopping_criteria
  
end module cg_names
