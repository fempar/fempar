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

  type :: nonlinear_solver_t
    private
      integer(ip)                              :: current_iteration
      integer(ip)                              :: max_number_iterations
      real(rp)                                 :: absolute_tolerance
      real(rp)                                 :: relative_tolerance
      character(len=:),           allocatable  :: convergence_criteria
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
                                   environment)
  implicit none
  class(nonlinear_solver_t)        , intent(inout) :: this
  character(len=*)                 , intent(in)    :: convergence_criteria
  real(rp)                         , intent(in)    :: abs_tol
  real(rp)                         , intent(in)    :: rel_tol
  integer(ip)                      , intent(in)    :: max_iters
  class(linear_solver_t)   , target, intent(in)    :: linear_solver
  class(environment_t)     , target, intent(in)    :: environment
  
  call this%free()
  
  ! Initialize options
  this%current_iteration     = 1
  this%convergence_criteria  = convergence_criteria
  this%absolute_tolerance    = abs_tol
  this%relative_tolerance    = rel_tol
  this%max_number_iterations = max_iters
  this%linear_solver  => linear_solver
  this%environment    => environment
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
    call this%update_solution() ! x + dx
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
subroutine nonlinear_solver_update_solution(this)
  implicit none
  class(nonlinear_solver_t), intent(inout) :: this
  call this%current_dof_values%axpby(1.0_rp, this%increment_dof_values, 1.0_rp)
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

end module nonlinear_solver_names
