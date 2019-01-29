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
  use solver_names
  use vector_space_names
  use vector_names
  use environment_names
  use linear_solver_names
  use fe_operator_names
  use FPL

  implicit none
# include "debug.i90"
  private

  !-------------------------------------------------------------------
  ! FPL keys of the parameters for nonlinear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: nls_rtol_key                   = 'nonlinear_solver_rtol'
  character(len=*), parameter :: nls_atol_key                   = 'nonlinear_solver_atol'
  character(len=*), parameter :: nls_stopping_criterium_key     = 'nonlinear_solver_stopping_criterium'
  character(len=*), parameter :: nls_max_num_iterations_key     = 'nonlinear_solver_max_num_iterations'
  character(len=*), parameter :: nls_print_iteration_output_key = 'nonlinear_solver_print_iteration_output'

  character(len=*), parameter :: abs_res_norm                   = 'abs_res_norm'                  ! |r(x_i)| <= abs_tol
  character(len=*), parameter :: rel_inc_norm                   = 'rel_inc_norm'                  ! |dx_i|   <= rel_norm*|x_i|
  character(len=*), parameter :: rel_r0_res_norm                = 'rel_r0_res_norm'               ! |r(x_i)| <= rel_tol*|r(x_0)|+abs_tol
  character(len=*), parameter :: abs_res_norm_and_rel_inc_norm  = 'abs_res_norm_and_rel_inc_norm' ! |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|

  !------------------------------------------------------------------
  ! Default values of the parameters for nonlinear solvers 
  !------------------------------------------------------------------
  real    (rp), parameter :: default_nls_rtol                     = 1.0e-09_rp
  real    (rp), parameter :: default_nls_atol                     = 1.0e-14_rp
  character(*), parameter :: default_nls_stopping_criterium       = rel_r0_res_norm
  integer (ip), parameter :: default_nls_max_iter                 = 10 
  logical     , parameter :: default_nls_print_iteration_output   = .true. 

  type, extends(solver_t) :: nonlinear_solver_t
  private
  integer(ip)                              :: current_iteration
  integer(ip)                              :: max_number_iterations
  real(rp)                                 :: absolute_tolerance
  real(rp)                                 :: relative_tolerance
  real(rp)                                 :: initial_residual_norm
  logical                                  :: print_iteration_output
  character(len=:),               allocatable  :: convergence_criteria
  class(vector_t),                pointer      :: current_dof_values => null()
  class(vector_t),                allocatable  :: initial_solution 
  class(vector_t),                allocatable  :: current_residual   
  class(vector_t),                allocatable  :: minus_current_residual
  class(vector_t),                allocatable  :: increment_dof_values
  class(linear_solver_t),         pointer      :: linear_solver
  class(environment_t),           pointer      :: environment
  class(fe_operator_t),           pointer      :: fe_operator

contains

  ! Public TBPs. 
  procedure :: create                           => nonlinear_solver_create
  procedure :: apply                            => nonlinear_solver_apply
  procedure :: apply_add                        => nonlinear_solver_apply_add
  procedure :: is_linear                        => nonlinear_solver_is_linear
  procedure :: update_matrix                    => nonlinear_solver_update_matrix  
  procedure :: free                             => nonlinear_solver_free

  ! Getters and checkers
  procedure :: has_converged                    => nonlinear_solver_has_converged
  procedure :: has_finished                     => nonlinear_solver_has_finished
  procedure :: get_current_iteration            => nonlinear_solver_get_current_iteration

  ! Setters
  procedure :: set_initial_solution             => nonlinear_solver_set_initial_solution

  ! Others

  ! Private TBPs
  procedure, private :: process_params                   => nonlinear_solver_process_params
  procedure, private :: update_solution                  => nonlinear_solver_update_solution
  procedure, private :: compute_residual                 => nonlinear_solver_compute_residual  
  procedure, private :: compute_tangent                  => nonlinear_solver_compute_tangent
  procedure, private :: set_current_dof_values           => nonlinear_solver_set_current_dof_values
  procedure, private :: print_iteration_output_header    => nonlinear_solver_print_iteration_output_header
  procedure, private :: print_current_iteration_output   => nonlinear_solver_print_current_iteration_output
  procedure, private :: print_final_output               => nonlinear_solver_print_final_output

end type nonlinear_solver_t

public :: nonlinear_solver_t

public :: nls_rtol_key
public :: nls_atol_key
public :: nls_stopping_criterium_key
public :: nls_max_num_iterations_key
public :: nls_print_iteration_output_key

public :: abs_res_norm
public :: rel_inc_norm
public :: rel_r0_res_norm
public :: abs_res_norm_and_rel_inc_norm

public :: default_nls_rtol
public :: default_nls_atol
public :: default_nls_stopping_criterium
public :: default_nls_max_iter 
public :: default_nls_print_iteration_output 

contains

!==============================================================================
subroutine nonlinear_solver_create(this, &
                                   parameters, &
                                   linear_solver, &
                                   fe_operator )
 implicit none
 class(nonlinear_solver_t)             , intent(inout) :: this
 type(parameterlist_t)                 , intent(in)    :: parameters
 class(linear_solver_t)        , target, intent(in)    :: linear_solver
 class(fe_operator_t)          , target, intent(in)    :: fe_operator

 call this%free()

 this%current_iteration     = 1
 call this%process_params(parameters)
 this%linear_solver  => linear_solver
 this%environment    => fe_operator%get_environment()
 this%fe_operator => fe_operator
 call this%fe_operator%create_range_vector(this%initial_solution)
 call this%initial_solution%init(0.0_rp)
end subroutine nonlinear_solver_create

!==============================================================================
subroutine nonlinear_solver_apply(this,x,y)
  !-----------------------------------------------------------------
  ! Find y such that A(y) = x
  !-----------------------------------------------------------------
 implicit none
 class(nonlinear_solver_t)     , intent(inout) :: this
 class(vector_t)               , intent(in)    :: x
 class(vector_t)               , intent(inout) :: y
 !x is in this context is the rhs of the nonlinear_operator_t not the evaluation_point, 
 !e.g. A(y) = x
 integer(ip) :: istat

 ! Initialize work data
 call this%fe_operator%create_domain_vector(this%increment_dof_values)

 ! Initialize nonlinear operator
 call this%set_current_dof_values(y)
 call this%current_dof_values%copy(this%initial_solution)
 this%current_iteration = 0  
 call this%increment_dof_values%init(0.0) 

 ! Compute initial residual
 call this%fe_operator%create_range_vector(this%current_residual)
 call this%compute_residual(x)

 ! Store initial residual norm for the stopping criterium that needs it
 if (this%convergence_criteria == rel_r0_res_norm) then
    this%initial_residual_norm = this%current_residual%nrm2()
 end if
 call this%fe_operator%create_range_vector(this%minus_current_residual)

 ! Print initial residual
 call this%print_iteration_output_header()
 call this%print_current_iteration_output()

 do while (.not. this%has_finished())
    this%current_iteration = this%current_iteration + 1
    call this%compute_tangent()    

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
    call this%update_solution() ! x + dx
    call this%compute_residual(x)
    call this%print_iteration_output_header()
    call this%print_current_iteration_output()
 end do
 call this%print_final_output()
end subroutine nonlinear_solver_apply

!==============================================================================
subroutine nonlinear_solver_apply_add(this, x, y)
  !-----------------------------------------------------------------
  ! Find y such that y = y + w with A(w) = x
  !-----------------------------------------------------------------
 class(nonlinear_solver_t), intent(inout)    :: this
 class(vector_t),        intent(in)    :: x
 class(vector_t),        intent(inout) :: y
 class(vector_t), allocatable          :: w
 type(vector_space_t), pointer         :: range_vector_space
 integer(ip)                           :: istat
 !-----------------------------------------------------------------
 range_vector_space => this%get_range_vector_space()
 call range_vector_space%create_vector(w)
 call this%apply(x,w)
 call y%axpby(1.0, w, 1.0)
 call w%free()
 deallocate(w, stat=istat); check(istat==0)
end subroutine nonlinear_solver_apply_add

!==============================================================================
function nonlinear_solver_is_linear(this) result(is_linear)
 class(nonlinear_solver_t), intent(in) :: this
 logical                            :: is_linear
 is_linear = .false.
end function nonlinear_solver_is_linear

!==============================================================================
subroutine nonlinear_solver_compute_residual(this,rhs)
 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 class(vector_t)          , intent(in)     :: rhs
 call this%fe_operator%apply(this%current_dof_values,this%current_residual)
 call this%current_residual%axpby(-1.0_rp, rhs,1.0_rp)
end subroutine nonlinear_solver_compute_residual

!==============================================================================
subroutine nonlinear_solver_compute_tangent(this)
 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 call this%fe_operator%compute_tangent() 
end subroutine nonlinear_solver_compute_tangent

!==============================================================================
subroutine nonlinear_solver_update_matrix(this, same_nonzero_pattern)
 implicit none 
 class(nonlinear_solver_t),        intent(inout) :: this
 logical,                       intent(in)    :: same_nonzero_pattern
 !-----------------------------------------------------------------
 call this%linear_solver%update_matrix(same_nonzero_pattern)  
end subroutine nonlinear_solver_update_matrix

!==============================================================================
subroutine nonlinear_solver_free(this)
 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 integer(ip) :: istat
 if(allocated(this%current_residual)) then
    call this%current_residual%free()
    deallocate(this%current_residual,stat=istat); check(istat==0)
 end if
 if(allocated(this%increment_dof_values)) then
    call this%increment_dof_values%free()
    deallocate(this%increment_dof_values,stat=istat); check(istat==0)
 end if
 if(allocated(this%minus_current_residual)) then
    call this%minus_current_residual%free()
    deallocate(this%minus_current_residual,stat=istat); check(istat==0)
 end if
 if(allocated(this%initial_solution)) then
    call this%initial_solution%free()
    deallocate(this%initial_solution,stat=istat); check(istat==0)
 end if
 this%current_dof_values => null()
 this%linear_solver      => null()
 this%environment        => null()
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
 case (rel_r0_res_norm) ! |r(x_i)| <= rel_tol*|r(x_0)| + abs_tol
    nonlinear_solver_has_converged = ( this%current_residual%nrm2() <= this%relative_tolerance*this%initial_residual_norm + this%absolute_tolerance )
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
 nonlinear_solver_has_finished = nonlinear_solver_has_finished .or.  (this%current_residual%nrm2() > 1.0e15_rp)
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
subroutine nonlinear_solver_set_initial_solution( this, initial_solution )
 implicit none
 class(nonlinear_solver_t), intent(inout) :: this
 class(vector_t)                 , intent(in)    :: initial_solution
 type(vector_space_t)            , pointer       :: nonlinear_operator_range 
 call initial_solution%GuardTemp()
 nonlinear_operator_range  => this%fe_operator%get_range_vector_space()
 if (.not. nonlinear_operator_range%belongs_to(initial_solution)) then
    write(0,'(a)') 'Warning: nonlinear_solver_t%set_initial_solution :: Ignoring initial solution; it does not belong to range(A)'
 else
    call this%initial_solution%copy(initial_solution)
 end if
 call initial_solution%CleanTemp()
end subroutine nonlinear_solver_set_initial_solution

!==============================================================================
subroutine nonlinear_solver_set_current_dof_values( this, current_dof_values )
 implicit none
 class(nonlinear_solver_t),        intent(inout) :: this
 class(vector_t),         target,  intent(in)    :: current_dof_values
 this%current_dof_values => current_dof_values
end subroutine nonlinear_solver_set_current_dof_values

!==============================================================================
subroutine nonlinear_solver_process_params(this,parameter_list)
 implicit none
 class(nonlinear_solver_t), intent(inout) :: this
 type(parameterlist_t)    , intent(in)    :: parameter_list
 integer(ip) :: FPLError
 ! Rtol
 if(parameter_list%isPresent(nls_rtol_key)) then
   assert(parameter_list%isAssignable(nls_rtol_key, this%relative_tolerance))
   FPLError   = parameter_list%Get(Key=nls_rtol_key, Value=this%relative_tolerance)
   assert(FPLError == 0)
 else
   this%relative_tolerance = default_nls_rtol
 endif
 ! Atol
 if(parameter_list%isPresent(nls_atol_key)) then
   assert(parameter_list%isAssignable(nls_atol_key, this%absolute_tolerance))
   FPLError   = parameter_list%Get(Key=nls_atol_key, Value=this%absolute_tolerance)
   assert(FPLError == 0)
 else
   this%absolute_tolerance = default_nls_atol
 endif
 ! Stopping criterium
 if(parameter_list%isPresent(nls_stopping_criterium_key)) then
   assert(parameter_list%isAssignable(nls_stopping_criterium_key, this%convergence_criteria))
   FPLError   = parameter_list%Get(Key=nls_stopping_criterium_key, Value=this%convergence_criteria)
   assert(FPLError == 0)
 else
   this%convergence_criteria = default_nls_stopping_criterium
 endif
 ! Max num iterations
 if(parameter_list%isPresent(nls_max_num_iterations_key)) then
   assert(parameter_list%isAssignable(nls_max_num_iterations_key, this%max_number_iterations))
   FPLError   = parameter_list%Get(Key=nls_max_num_iterations_key, Value=this%max_number_iterations)
   assert(FPLError == 0)
 else
   this%max_number_iterations = default_nls_max_iter
 endif
 ! Print iteration output 
 if(parameter_list%isPresent(nls_print_iteration_output_key)) then
   assert(parameter_list%isAssignable(nls_print_iteration_output_key, this%print_iteration_output))
   FPLError   = parameter_list%Get(Key=nls_print_iteration_output_key, Value=this%print_iteration_output)
   assert(FPLError == 0)
 else
   this%print_iteration_output = default_nls_print_iteration_output
 endif
end subroutine nonlinear_solver_process_params

!==============================================================================
subroutine nonlinear_solver_update_solution(this)
 implicit none
 class(nonlinear_solver_t), intent(inout) :: this
 call this%current_dof_values%axpby(1.0_rp, this%increment_dof_values, 1.0_rp)
end subroutine nonlinear_solver_update_solution

!==============================================================================
subroutine nonlinear_solver_print_iteration_output_header(this)

 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 
 if ( this%print_iteration_output ) then
   if ( this%environment%am_i_l1_root() ) then
      select case (this%convergence_criteria)
      case (abs_res_norm) !  |r(x_i)| <= abs_tol
         write(*,'(a10,2a21)') 'NL iter', '|r(x_i)|', 'abs_tol'
      case (rel_inc_norm) ! |dx_i| <= rel_norm*|x_i|
         write(*,'(a10,3a21)') 'NL iter', '|dx_i|', 'rel_tol*|x_i|', 'rel_tol'
      case (rel_r0_res_norm) ! |r(x_i)| <= rel_tol*|r(x_0)| + abs_tol
         write(*,'(a10,3a21)') 'NL iter', '|r(x_i)|', 'rel_tol*|r(x_0)|', 'rel_tol'
      case (abs_res_norm_and_rel_inc_norm) !  |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|
         write(*,'(a8,5a15)') 'NL iter', '|r(x_i)|', 'abs_tol', '|dx_i|', 'rel_tol*|x_i|', 'rel_tol'
      case default
         mcheck(.false., 'Unknown convergence criterium: '//this%convergence_criteria )
      end select
   end if
 end if
end subroutine nonlinear_solver_print_iteration_output_header

!==============================================================================
subroutine nonlinear_solver_print_current_iteration_output(this)

 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 real(rp) :: current_residual_nrm2
 real(rp) :: current_dof_values_nrm2, increment_dof_values_nrm2
 
 if ( this%print_iteration_output ) then
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
      if ( this%environment%am_i_l1_root()) write(*,'(i10,3es21.10)') this%current_iteration, current_residual_nrm2, &
           this%relative_tolerance*this%initial_residual_norm+this%absolute_tolerance, this%relative_tolerance

   case (abs_res_norm_and_rel_inc_norm) !  |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|
      current_residual_nrm2 = this%current_residual%nrm2()
      current_dof_values_nrm2 = this%current_dof_values%nrm2()
      increment_dof_values_nrm2 = this%increment_dof_values%nrm2()
      if (this%environment%am_i_l1_root() ) write(*,'(i8,5es15.4)') this%current_iteration, current_residual_nrm2, this%absolute_tolerance, &
           increment_dof_values_nrm2, this%relative_tolerance*current_dof_values_nrm2, this%relative_tolerance
   case default
      mcheck(.false., 'Unknown convergence criterium: '//this%convergence_criteria )
   end select
 end if
end subroutine nonlinear_solver_print_current_iteration_output

!==============================================================================
subroutine nonlinear_solver_print_final_output(this)
 implicit none
 class(nonlinear_solver_t), intent(inout)  :: this
 if ( this%has_converged() ) then
    if ( this%environment%am_i_l1_root()  .and.  this%print_iteration_output ) write(*,*) ' ---- Newton NL solver converged in', this%current_iteration, 'iterations ---- '
 else if (.not. this%has_converged() ) then
    assert (this%has_finished() .eqv. .true.)
    if ( this%environment%am_i_l1_root() ) write(*,*) ' ---- Newton NL failed to converge after', this%current_iteration, 'iterations --- '
 end if
end subroutine nonlinear_solver_print_final_output


end module nonlinear_solver_names
