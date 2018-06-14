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
module nonlinear_solver_names
  use types_names
  use memor_names
  use vector_names
  use matrix_names 
  use environment_names
  
  !use base_sparse_matrix_names
  use linear_solver_names
  use fe_nonlinear_operator_names
  implicit none
# include "debug.i90"
  private

  character(len=*), parameter :: abs_res_norm                  = 'abs_res_norm'                  ! |r(x_i)| <= abs_tol
  character(len=*), parameter :: rel_inc_norm                  = 'rel_inc_norm'                  ! |dx_i|   <= rel_norm*|x_i|
  character(len=*), parameter :: rel_r0_res_norm               = 'rel_r0_res_norm'               ! |r(x_i)| <= rel_tol*|r(x_0)|
  character(len=*), parameter :: abs_res_norm_and_rel_inc_norm = 'abs_res_norm_and_rel_inc_norm' ! |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|

  character(len=*), parameter :: static              = 'static' 
  character(len=*), parameter :: cubic_backtracking  = 'cubic_backtracking' 
  character(len=*), parameter :: L2_iterative_secant = 'L2_iterative_secant'
    
  type, abstract :: line_search_t 
     private 
     real(rp)                     :: step_length 
     integer(ip)                  :: current_iterate 
     integer(ip)                  :: max_num_iterates 
   contains 
     procedure(create_interface)               , private, deferred :: create 
     procedure(free_interface)                 , private, deferred :: free 
     procedure(determine_step_length_interface), private, deferred :: determine_step_length
     procedure                                 , private           :: get_step_length => line_search_get_step_length 
  end type line_search_t

  abstract interface
     subroutine create_interface(this, nonlinear_operator) 
       import :: line_search_t, fe_nonlinear_operator_t 
       implicit none 
       class(line_search_t)          , intent(inout) :: this 
       class(fe_nonlinear_operator_t), intent(in)    :: nonlinear_operator 
     end subroutine create_interface

     subroutine free_interface(this) 
       import :: line_search_t
       implicit none 
       class(line_search_t)          , intent(inout) :: this 
     end subroutine free_interface

     subroutine determine_step_length_interface(this, increment_dof_values, current_dof_values, nonlinear_operator) 
       import :: line_search_t, vector_t, fe_nonlinear_operator_t 
       implicit none 
       class(line_search_t)          , intent(inout) :: this  
       class(vector_t)               , intent(in)    :: increment_dof_values 
       class(vector_t)               , intent(inout) :: current_dof_values
       class(fe_nonlinear_operator_t), intent(inout) :: nonlinear_operator 
     end subroutine determine_step_length_interface
  end interface

  type, extends(line_search_t) :: static_line_search_t 
     private 
   contains 
     procedure, private :: create                => static_line_search_create 
     procedure, private :: free                  => static_line_search_free 
     procedure, private :: determine_step_length => static_line_search_determine_step_length
  end type static_line_search_t

  type, extends(line_search_t) :: cubic_backtracking_line_search_t 
     private 
     real(ip)                     :: alpha 
     class(vector_t), allocatable :: initial_dof_values 
     class(vector_t), allocatable :: jacobian_incr_dof_values 
   contains 
     procedure, private :: create                => cubic_backtracking_line_search_create 
     procedure, private :: free                  => cubic_backtracking_line_search_free 
     procedure, private :: determine_step_length => cubic_backtracking_line_search_determine_step_length
  end type cubic_backtracking_line_search_t
  
    type, extends(line_search_t) :: L2_iterative_secant_line_search_t 
     private 
     class(vector_t), allocatable :: initial_dof_values 
   contains 
     procedure, private :: create                => L2_iterative_secant_line_search_create 
     procedure, private :: free                  => L2_iterative_secant_line_search_free 
     procedure, private :: determine_step_length => L2_iterative_secant_line_search_determine_step_length
  end type L2_iterative_secant_line_search_t

public :: static, cubic_backtracking, L2_iterative_secant 
  
  type :: nonlinear_solver_t
    private 
      integer(ip)                              :: current_iteration
      integer(ip)                              :: max_number_iterations
      real(rp)                                 :: absolute_tolerance
      real(rp)                                 :: relative_tolerance
      character(len=:),           allocatable  :: convergence_criteria
      class(line_search_t),       allocatable  :: line_search 
      class(vector_t),            pointer      :: current_dof_values => null()
      class(vector_t),            pointer      :: current_residual   => null()
      class(vector_t),            allocatable  :: minus_current_residual
      class(vector_t),            allocatable  :: initial_residual
      class(vector_t),            allocatable  :: increment_dof_values
      class(linear_solver_t),     pointer      :: linear_solver
      class(environment_t),       pointer      :: environment
    contains

      ! Public TBPs. They MUST be called in this order.
      procedure :: create                           => nonlinear_solver_create
      procedure :: solve                            => nonlinear_solver_solve
      procedure :: free                             => nonlinear_solver_free

      ! Getters and checkers
      procedure :: has_converged                     => nonlinear_solver_has_converged
      procedure :: has_finished                      => nonlinear_solver_has_finished
      procedure :: get_current_iteration             => nonlinear_solver_get_current_iteration

      ! Private TBPs
      !procedure, private :: initialize                       => nonlinear_solver_initialize
      procedure, private :: update_solution                  => nonlinear_solver_update_solution
      procedure, private :: print_iteration_output_header    => nonlinear_solver_print_iteration_output_header
      procedure, private :: print_current_iteration_output   => nonlinear_solver_print_current_iteration_output
      procedure, private :: print_final_output               => nonlinear_solver_print_final_output

  end type nonlinear_solver_t
  
  public :: nonlinear_solver_t
  public :: abs_res_norm
  public :: rel_inc_norm
  public :: rel_r0_res_norm
  public :: abs_res_norm_and_rel_inc_norm

contains

!==============================================================================
subroutine nonlinear_solver_create(this, &
                                   convergence_criteria, &
                                   abs_tol, &
                                   rel_tol, &
                                   max_iters, &
                                   linear_solver, &
                                   environment, &
                                   line_search_type )
  implicit none
  class(nonlinear_solver_t)        , intent(inout) :: this
  character(len=*)                 , intent(in)    :: convergence_criteria
  real(rp)                         , intent(in)    :: abs_tol
  real(rp)                         , intent(in)    :: rel_tol
  integer(ip)                      , intent(in)    :: max_iters
  class(linear_solver_t)   , target, intent(in)    :: linear_solver
  class(environment_t)     , target, intent(in)    :: environment
  character(len=*)       , optional, intent(in)    :: line_search_type 
  
  character(:), allocatable :: line_search_type_ 
    
  call this%free()
  
  ! Initialize options 
  this%current_iteration     = 1
  this%convergence_criteria  = convergence_criteria
  this%absolute_tolerance    = abs_tol
  this%relative_tolerance    = rel_tol
  this%max_number_iterations = max_iters
  this%linear_solver  => linear_solver
  this%environment    => environment
  
  if (present(line_search_type)) then  
     line_search_type_ = line_search_type 
  else 
     line_search_type_ = static 
  end if
  call make_line_search(line_search_type_, this%line_search) 
  
end subroutine nonlinear_solver_create

!==============================================================================
subroutine nonlinear_solver_solve(this,nonlinear_operator,unknown)

  implicit none
  class(nonlinear_solver_t)     , intent(inout) :: this
  class(fe_nonlinear_operator_t), intent(inout) :: nonlinear_operator
  class(vector_t), target       , intent(inout) :: unknown
  
  integer(ip) :: istat
  
  
  ! Initialize work data
  if(allocated(this%increment_dof_values)) then
    call this%increment_dof_values%free()
    deallocate(this%increment_dof_values,stat=istat); check(istat==0)
  end if
  call nonlinear_operator%create_domain_vector(this%increment_dof_values)
  
  ! Create line search 
  call this%line_search%create(nonlinear_operator) 
  
  ! Initialize nonlinear operator
  this%current_dof_values => unknown
  this%current_iteration = 0  
  call this%increment_dof_values%init(0.0) 
  
  ! Compute initial residual
  call nonlinear_operator%set_evaluation_point(unknown)
  call nonlinear_operator%compute_residual()
  this%current_residual => nonlinear_operator%get_translation()
  
  ! Store initial residual for the stopping criterium that needs it
  if (this%convergence_criteria == rel_r0_res_norm) then
    call nonlinear_operator%create_range_vector(this%initial_residual)
    this%initial_residual = this%current_residual
  end if
  call nonlinear_operator%create_range_vector(this%minus_current_residual)

  ! Print initial residual
  call this%print_iteration_output_header()
  call this%print_current_iteration_output()
  
  do while (.not. this%has_finished())
    this%current_iteration = this%current_iteration + 1
    call nonlinear_operator%compute_tangent()    

    ! Force numerical set-up
    call this%linear_solver%update_matrix(same_nonzero_pattern=.true.)
    
    ! sbadia : An extra copy is required to store "-residual" !!!!
    ! sbadia : To be optimized (Why don't we change the concept of residual?)
    ! ifort 18.0.0 (I do not know about other versions) had issues with the
    ! expression -this%current_residual passed to linear_solver%apply(...) below. 
    ! Thus, I had to explicitly change the sign of the residual on explicitly 
    ! allocated memory (handled by create/free)
    call this%minus_current_residual%scal(-1.0_rp, this%current_residual)
    
    call this%linear_solver%apply( this%minus_current_residual, this%increment_dof_values )
    call this%line_search%determine_step_length(this%increment_dof_values, this%current_dof_values, nonlinear_operator) 
    call this%update_solution(this%line_search%get_step_length()) ! x + lambda*dx
    call nonlinear_operator%set_evaluation_point(this%current_dof_values)
    call nonlinear_operator%compute_residual()
    call this%print_iteration_output_header()
    call this%print_current_iteration_output()
  end do
  call this%print_final_output()
end subroutine nonlinear_solver_solve

!==============================================================================
subroutine nonlinear_solver_free(this)
  implicit none
  class(nonlinear_solver_t), intent(inout)  :: this
  integer(ip) :: istat
  if(allocated(this%initial_residual)) then
    call this%initial_residual%free()
    deallocate(this%initial_residual,stat=istat); check(istat==0)
  end if
  if(allocated(this%increment_dof_values)) then
    call this%increment_dof_values%free()
    deallocate(this%increment_dof_values,stat=istat); check(istat==0)
  end if
  if(allocated(this%minus_current_residual)) then
    call this%minus_current_residual%free()
    deallocate(this%minus_current_residual,stat=istat); check(istat==0)
  end if
  this%current_residual   => null()
  this%current_dof_values => null()
  this%linear_solver      => null()
  this%environment        => null()
  if ( allocated(this%line_search)) call this%line_search%free() 
end subroutine nonlinear_solver_free

!==============================================================================
function nonlinear_solver_has_converged(this)

  implicit none
  class(nonlinear_solver_t), intent(in) :: this
  logical                                   :: nonlinear_solver_has_converged

  select case (this%convergence_criteria)
  case (abs_res_norm) !  |r(x_i)| <= abs_tol
    nonlinear_solver_has_converged = ( this%current_residual%nrm2() <= this%absolute_tolerance )
  case (rel_inc_norm) ! |dx_i| <= rel_norm*|x_i|
    nonlinear_solver_has_converged = ( this%increment_dof_values%nrm2() <= this%relative_tolerance*this%current_dof_values%nrm2() )
  case (rel_r0_res_norm) ! |r(x_i)| <= rel_tol*|r(x_0)|
    nonlinear_solver_has_converged = ( this%current_residual%nrm2() <= this%relative_tolerance*this%initial_residual%nrm2() )
  case (abs_res_norm_and_rel_inc_norm) !  |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|
    nonlinear_solver_has_converged = ( this%current_residual%nrm2() <= this%absolute_tolerance ) .and. &
                                        ( this%increment_dof_values%nrm2() <= this%relative_tolerance*this%current_dof_values%nrm2() )
  case default
    mcheck(.false., 'Unknown convergence criterium: '//this%convergence_criteria )
  end select
  call this%environment%l1_lgt1_bcast(nonlinear_solver_has_converged)
end function nonlinear_solver_has_converged

!==============================================================================
function nonlinear_solver_has_finished(this)
  implicit none
  class(nonlinear_solver_t), intent(in) :: this
  logical                                   :: nonlinear_solver_has_finished
  nonlinear_solver_has_finished = this%has_converged()
  nonlinear_solver_has_finished = nonlinear_solver_has_finished .or.  (this%current_iteration > this%max_number_iterations)
  nonlinear_solver_has_finished = nonlinear_solver_has_finished .or.  (this%current_residual%nrm2() > 1.0e10_rp)
  nonlinear_solver_has_finished = nonlinear_solver_has_finished .and. (this%current_iteration > 0)
  call this%environment%l1_lgt1_bcast(nonlinear_solver_has_finished)
end function nonlinear_solver_has_finished

!==============================================================================
function nonlinear_solver_get_current_iteration(this)
  implicit none
  class(nonlinear_solver_t), intent(in) :: this
  integer(ip)                               :: nonlinear_solver_get_current_iteration
  nonlinear_solver_get_current_iteration = this%current_iteration
end function nonlinear_solver_get_current_iteration

!==============================================================================
subroutine nonlinear_solver_update_solution(this, step_length)
  implicit none
  class(nonlinear_solver_t), intent(inout) :: this
  real(rp)                 , intent(in)    :: step_length
  
  call this%current_dof_values%axpby(step_length, this%increment_dof_values, 1.0_rp)
end subroutine nonlinear_solver_update_solution

!==============================================================================
subroutine nonlinear_solver_print_iteration_output_header(this)

  implicit none
  class(nonlinear_solver_t), intent(inout)  :: this

  if ( this%environment%am_i_l1_root() ) then
  select case (this%convergence_criteria)
  case (abs_res_norm) !  |r(x_i)| <= abs_tol
    write(*,'(a10,2a21)') 'NL iter', '|r(x_i)|', 'abs_tol'
  case (rel_inc_norm) ! |dx_i| <= rel_norm*|x_i|
    write(*,'(a10,3a21)') 'NL iter', '|dx_i|', 'rel_tol*|x_i|', 'rel_tol'
  case (rel_r0_res_norm) ! |r(x_i)| <= rel_tol*|r(x_0)|
    write(*,'(a10,3a21)') 'NL iter', '|r(x_i)|', 'rel_tol*|r(x_0)|', 'rel_tol'
  case (abs_res_norm_and_rel_inc_norm) !  |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|
    write(*,'(a8,5a15)') 'NL iter', '|r(x_i)|', 'abs_tol', '|dx_i|', 'rel_tol*|x_i|', 'rel_tol'
  case default
    mcheck(.false., 'Unknown convergence criterium: '//this%convergence_criteria )
  end select
  end if
  
end subroutine nonlinear_solver_print_iteration_output_header

!==============================================================================
subroutine nonlinear_solver_print_current_iteration_output(this)

  implicit none
  class(nonlinear_solver_t), intent(inout)  :: this
  real(rp) :: current_residual_nrm2, initial_residual_nrm2
  real(rp) :: current_dof_values_nrm2, increment_dof_values_nrm2
    
  select case (this%convergence_criteria)
  case (abs_res_norm) !  |r(x_i)| <= abs_tol
    current_residual_nrm2 = this%current_residual%nrm2()
    if ( this%environment%am_i_l1_root()) write(*,'(i10,2es21.10)') this%current_iteration, current_residual_nrm2, this%absolute_tolerance

  case (rel_inc_norm) ! |dx_i| <= rel_norm*|x_i|
    current_dof_values_nrm2 = this%current_dof_values%nrm2()
    increment_dof_values_nrm2 = this%increment_dof_values%nrm2()
    if ( this%environment%am_i_l1_root()) write(*,'(i10,3es21.10)') this%current_iteration, increment_dof_values_nrm2, &
                                          this%relative_tolerance*current_dof_values_nrm2, this%relative_tolerance

  case (rel_r0_res_norm) ! |r(x_i)| <= rel_tol*|r(x_0)|
    current_residual_nrm2 = this%current_residual%nrm2()
    initial_residual_nrm2 = this%initial_residual%nrm2()
    if ( this%environment%am_i_l1_root()) write(*,'(i10,3es21.10)') this%current_iteration, current_residual_nrm2, &
                                          this%relative_tolerance*initial_residual_nrm2, this%relative_tolerance

  case (abs_res_norm_and_rel_inc_norm) !  |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|
    current_residual_nrm2 = this%current_residual%nrm2()
    current_dof_values_nrm2 = this%current_dof_values%nrm2()
    increment_dof_values_nrm2 = this%increment_dof_values%nrm2()
    if (this%environment%am_i_l1_root() ) write(*,'(i8,5es15.4)') this%current_iteration, current_residual_nrm2, this%absolute_tolerance, &
                                     increment_dof_values_nrm2, this%relative_tolerance*current_dof_values_nrm2, this%relative_tolerance
  case default
    mcheck(.false., 'Unknown convergence criterium: '//this%convergence_criteria )
  end select
end subroutine nonlinear_solver_print_current_iteration_output

!==============================================================================
subroutine nonlinear_solver_print_final_output(this)
  implicit none
  class(nonlinear_solver_t), intent(inout)  :: this
  if ( this%has_converged() ) then
   if ( this%environment%am_i_l1_root() ) write(*,*) ' ---- Newton NL solver converged in', this%current_iteration, 'iterations ---- '
  else if (.not. this%has_converged() ) then
   assert (this%has_finished() .eqv. .true.)
   if ( this%environment%am_i_l1_root() ) write(*,*) ' ---- Newton NL failed to converge after', this%current_iteration, 'iterations --- '
  end if
end subroutine nonlinear_solver_print_final_output

! =============================================================================
subroutine make_line_search(line_search_type, line_search) 
  implicit none 
  character(len=*)                 , intent(in)      :: line_search_type 
  class(line_search_t), allocatable, intent(inout)   :: line_search 

  select case ( line_search_type ) 
  case ( static ) 
     allocate ( static_line_search_t :: line_search )
  case ( cubic_backtracking ) 
     allocate ( cubic_backtracking_line_search_t :: line_search )
  case ( L2_iterative_secant ) 
     allocate ( L2_iterative_secant_line_search_t :: line_search )
  case DEFAULT 
     massert(.false., 'Specified line search technique is not valid') 
  end select

end subroutine  make_line_search

! =============================================================================
function line_search_get_step_length (this) 
  implicit none 
  class(line_search_t)     , intent(inout) :: this 
  real(rp)  :: line_search_get_step_length 
  line_search_get_step_length = this%step_length 
end function  line_search_get_step_length 

! =============================================================================
subroutine static_line_search_free(this) 
  class(static_line_search_t)     , intent(inout) :: this 
  ! Do-nothing 
end subroutine static_line_search_free

! =============================================================================
subroutine static_line_search_create(this, nonlinear_operator) 
  class(static_line_search_t)     , intent(inout) :: this 
  class(fe_nonlinear_operator_t)  , intent(in)    :: nonlinear_operator 
  ! Do-nothing 
end subroutine static_line_search_create

! =============================================================================
subroutine cubic_backtracking_line_search_free(this) 
  class(cubic_backtracking_line_search_t), intent(inout) :: this 
  integer(ip) :: istat 

  if(allocated(this%initial_dof_values)) then
     call this%initial_dof_values%free()
     deallocate(this%initial_dof_values,stat=istat); check(istat==0)
  end if
  if(allocated(this%jacobian_incr_dof_values)) then
     call this%jacobian_incr_dof_values%free()
     deallocate(this%jacobian_incr_dof_values,stat=istat); check(istat==0)
  end if

end subroutine cubic_backtracking_line_search_free

! =============================================================================
subroutine cubic_backtracking_line_search_create(this, nonlinear_operator) 
  class(cubic_backtracking_line_search_t) , intent(inout) :: this 
  class(fe_nonlinear_operator_t)          , intent(in)    :: nonlinear_operator 
  call nonlinear_operator%create_domain_vector(this%initial_dof_values)
  call nonlinear_operator%create_domain_vector(this%jacobian_incr_dof_values) 
end subroutine cubic_backtracking_line_search_create

!==============================================================================
subroutine static_line_search_determine_step_length(this, increment_dof_values, current_dof_values, nonlinear_operator)  
  implicit none 
  class(static_line_search_t)      , intent(inout) :: this 
  class(vector_t)                  , intent(in)    :: increment_dof_values
  class(vector_t)                  , intent(inout) :: current_dof_values
  class(fe_nonlinear_operator_t)   , intent(inout) :: nonlinear_operator 

  ! Take full step 
  this%step_length = 1.0_rp 
end subroutine static_line_search_determine_step_length

!==============================================================================
! Algorithm from 'Brune, P. Knepley, M. Smith, B. and Tu, X. Composing scalable 
! nonlinear algebraic solvers. SIAM Review, 57(4) p. 540'
subroutine cubic_backtracking_line_search_determine_step_length(this, increment_dof_values, current_dof_values, nonlinear_operator) 
  implicit none 
  class(cubic_backtracking_line_search_t) , intent(inout) :: this 
  class(vector_t)                         , intent(in)    :: increment_dof_values
  class(vector_t)                         , intent(inout) :: current_dof_values
  class(fe_nonlinear_operator_t)          , intent(inout) :: nonlinear_operator 

  class(matrix_t), pointer :: jacobian   
  class(vector_t), pointer :: current_residual
  real(rp)                 :: initial_residual_norm 
  real(rp)                 :: previous_residual_norm 
  real(rp)                 :: previous_step_length  
  real(rp)                 :: trial_step_length 
  real(rp)                 :: t1, t2
  real(rp)                 :: a, b, d, s 

  ! Initialize parameters  
  this%step_length      = 1.0_rp 
  this%current_iterate  = 0
  this%max_num_iterates = 2
  this%alpha            = 1.0e-4_rp 

  current_residual      => nonlinear_operator%get_translation()
  initial_residual_norm =  current_residual%nrm2() 

  ! Copy of previous solution 
  call this%initial_dof_values%copy(current_dof_values) 
  ! Evaluate s = R(x)'*J(x)*dx  
  jacobian  => nonlinear_operator%get_matrix() 
  call jacobian%apply(increment_dof_values, this%jacobian_incr_dof_values)
  s = current_residual%dot(this%jacobian_incr_dof_values) 
  
  call current_dof_values%axpby(this%step_length, increment_dof_values, 1.0_rp) 
  call nonlinear_operator%set_evaluation_point(current_dof_values)
  call nonlinear_operator%compute_residual()
  if ( current_residual%nrm2()**2 <= initial_residual_norm + 2.0_rp*this%alpha*this%step_length*s ) then 
     ! Initial time steps fulfills the Armijo condition   
     call nonlinear_operator%set_evaluation_point(this%initial_dof_values) 
     return 
  end if

  previous_residual_norm = current_residual%nrm2()
  previous_step_length   = this%step_length 
  trial_step_length      = -s/(current_residual%nrm2()**2 - initial_residual_norm**2 - 2.0_rp*this%step_length*s)
  if ( trial_step_length < 0.1_rp*this%step_length ) then 
     this%step_length = 0.1_rp*this%step_length 
  elseif ( trial_step_length > 0.5_rp*this%step_length ) then 
     this%step_length = 0.5_rp*this%step_length 
  else 
     this%step_length = trial_step_length 
  end if

  call current_dof_values%axpby(this%step_length - previous_step_length, increment_dof_values, 1.0_rp) 
  call nonlinear_operator%set_evaluation_point(current_dof_values)
  call nonlinear_operator%compute_residual()
  if ( current_residual%nrm2()**2 <= initial_residual_norm + 2.0_rp*this%alpha*this%step_length*s ) then 
     ! Quadratic interpolation of F(step_length) is enough  
     call nonlinear_operator%set_evaluation_point(this%initial_dof_values) 
     return 
  end if

  do while ( this%current_iterate < this%max_num_iterates) 
     this%current_iterate=this%current_iterate+1

     t1 = 0.5_rp*( previous_residual_norm**2 - initial_residual_norm**2) - this%step_length*s 
     t2 = 0.5_rp*(current_residual%nrm2()**2 - initial_residual_norm**2) - this%step_length*s

     a  = ( t1/this%step_length**2 - t2/previous_step_length**2)/(this%step_length - previous_step_length)
     b  = this%step_length*( t2/previous_step_length**2 - t1/this%step_length**2)/(this%step_length - previous_step_length)
     d  = max ( b**2.0_rp - 3.0_rp*a*s, 0.0_rp )

     if (abs(a) < 1.0e-14_rp) then 
        trial_step_length = -s/(2.0_rp*b)
     else 
        trial_step_length = (sqrt(d)-b)/(3.0_rp*a) 
     end if
     
     previous_step_length   = this%step_length 
     previous_residual_norm = current_residual%nrm2() 
     ! Update step length 
     if ( trial_step_length < 0.1_rp*this%step_length ) then 
        this%step_length = 0.1_rp*this%step_length 
     elseif ( trial_step_length > 0.5_rp*this%step_length ) then 
        this%step_length = 0.5_rp*this%step_length 
     else 
        this%step_length = trial_step_length 
     end if

     call current_dof_values%axpby(this%step_length - previous_step_length, increment_dof_values, 1.0_rp) 
     call nonlinear_operator%set_evaluation_point(current_dof_values)
     call nonlinear_operator%compute_residual()
     if ( current_residual%nrm2()**2 <= initial_residual_norm + 2.0_rp*this%alpha*this%step_length*s) then 
        ! Cubic interpolation of F(step_length) is enough
        call nonlinear_operator%set_evaluation_point(this%initial_dof_values) 
        return 
     end if

  end do

  ! If it did not achieve to satisfy the condition, go back to the previous state 
  call nonlinear_operator%set_evaluation_point(this%initial_dof_values)
end subroutine cubic_backtracking_line_search_determine_step_length

! =============================================================================
subroutine L2_iterative_secant_line_search_free(this) 
  class(L2_iterative_secant_line_search_t), intent(inout) :: this 
  integer(ip) :: istat 

  if(allocated(this%initial_dof_values)) then
     call this%initial_dof_values%free()
     deallocate(this%initial_dof_values,stat=istat); check(istat==0)
  end if

end subroutine L2_iterative_secant_line_search_free

! =============================================================================
subroutine L2_iterative_secant_line_search_create(this, nonlinear_operator) 
  class(L2_iterative_secant_line_search_t) , intent(inout) :: this 
  class(fe_nonlinear_operator_t)          , intent(in)    :: nonlinear_operator 
  call nonlinear_operator%create_domain_vector(this%initial_dof_values)
end subroutine L2_iterative_secant_line_search_create

!==============================================================================
! Algorithm from 'Brune, P. Knepley, M. Smith, B. and Tu, X. Composing scalable 
! nonlinear algebraic solvers. SIAM Review, 57(4) p. 540'
subroutine L2_iterative_secant_line_search_determine_step_length(this, increment_dof_values, current_dof_values, nonlinear_operator) 
  implicit none 
  class(L2_iterative_secant_line_search_t) , intent(inout) :: this 
  class(vector_t)                          , intent(in)    :: increment_dof_values
  class(vector_t)                          , intent(inout) :: current_dof_values
  class(fe_nonlinear_operator_t)           , intent(inout) :: nonlinear_operator 

  class(vector_t), pointer :: current_residual
  real(rp)                 :: previous_step_length
  real(rp)                 :: trial_step_length 
  real(rp)                 :: previous_res_norm 
  real(rp)                 :: current_res_norm  
  real(rp)                 :: half_step_res_norm 
  real(rp)                 :: curr_grad_res_norm 
  real(rp)                 :: prev_grad_res_norm 
  real(rp)                 :: step_length_update

  ! Initialize parameters  
  this%step_length      = 1.0_rp 
  this%current_iterate  = 0
  this%max_num_iterates = 3

  previous_step_length = 0.0_rp 
  current_residual     => nonlinear_operator%get_translation()
  previous_res_norm    =  current_residual%nrm2() 
  call this%initial_dof_values%copy(current_dof_values)

  do while ( this%current_iterate < this%max_num_iterates )
     this%current_iterate = this%current_iterate + 1

     ! Evaluate F( (lambda_i + lambda_{i-1})/2 )
     call current_dof_values%axpby(0.5_rp*(previous_step_length+this%step_length), increment_dof_values, 1.0_rp)
     call nonlinear_operator%set_evaluation_point(current_dof_values)
     call nonlinear_operator%compute_residual()
     half_step_res_norm = current_residual%nrm2() 

     ! Evaluate F(lambda_i) 
     call current_dof_values%axpby(0.5_rp*(this%step_length-previous_step_length), increment_dof_values, 1.0_rp)
     call nonlinear_operator%set_evaluation_point(current_dof_values)
     call nonlinear_operator%compute_residual()
     current_res_norm = current_residual%nrm2() 

     curr_grad_res_norm = (3.0_rp*current_res_norm**2 - 4.0_rp*half_step_res_norm**2 +        previous_res_norm**2)/(this%step_length-previous_step_length)
     prev_grad_res_norm = (       current_res_norm**2 - 4.0_rp*half_step_res_norm**2 + 3.0_rp*previous_res_norm**2)/(this%step_length-previous_step_length)

     ! Update step length 
     step_length_update = -curr_grad_res_norm*(this%step_length-previous_step_length)/(curr_grad_res_norm-prev_grad_res_norm)
     if ( abs(step_length_update) > 1e-5_rp ) then 
        this%step_length  = this%step_length  + step_length_update 
     else 
        call nonlinear_operator%set_evaluation_point(this%initial_dof_values)
        return 
     end if

     ! Return nonlinear operator to the starting point  
     previous_step_length = this%step_length 
     previous_res_norm    = current_residual%nrm2()
  end do

end subroutine L2_iterative_secant_line_search_determine_step_length

end module nonlinear_solver_names
