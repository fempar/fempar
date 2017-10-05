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
module hts_nonlinear_solver_names 
  use fempar_names
  use hts_discrete_integration_names
  use base_sparse_matrix_names 
  implicit none
# include "debug.i90"
  private
     
  type :: hts_line_search_t 
  real(rp)                        :: step_length 
  real(rp)                        :: alpha
  integer(ip)                     :: n, current_iterate 
  class(vector_t), allocatable    :: initial_dof_values 
  class(vector_t), allocatable    :: initial_residual 
  contains
  procedure :: initialize         => line_search_initialize 
  procedure :: cubic_backtracking => line_search_cubic_backtracking  
  procedure :: free               => line_search_free 
  
  end type hts_line_search_t 
  
  type :: hts_nonlinear_solver_t
     integer(ip)                                :: current_iteration 
     integer(ip)                                :: ideal_num_iterations 
     integer(ip)                                :: max_num_iterations
	    real(rp)                                   :: absolute_tolerance
     real(rp)                                   :: relative_tolerance
     type(hts_line_search_t)                    :: line_search  
     character(len=:) , allocatable             :: convergence_criteria 
     type(fe_affine_operator_t) , pointer       :: fe_affine_operator 
     class(vector_t), pointer                   :: current_dof_values 
     real(rp)                                   :: current_lagrange_multiplier 
     class(vector_t), pointer                   :: current_rhs 
     real(rp)                                   :: constraint_value 
     type(serial_scalar_array_t), pointer       :: constraint_vector 
     class(vector_t), allocatable               :: increment_dof_values 
     real(rp)                                   :: increment_lagrange_multiplier 
	    class(vector_t), allocatable               :: residual 
     class(vector_t), allocatable               :: initial_residual 
     type(direct_solver_t)                      :: direct_solver

    contains
     procedure :: create                           => hts_nonlinear_solver_create 
     procedure :: start_new_iteration              => hts_nonlinear_solver_start_new_iteration 
     procedure :: initialize                       => hts_nonlinear_solver_initialize 
     procedure :: compute_residual                 => hts_nonlinear_solver_compute_residual 
     procedure :: compute_Jacobian                 => hts_nonlinear_solver_compute_Jacobian 
     procedure :: update_solution                  => hts_nonlinear_solver_update_solution
     procedure :: converged                        => hts_nonlinear_solver_converged 
     procedure :: finished                         => hts_nonlinear_solver_finished 
     procedure :: get_current_iteration            => hts_nonlinear_solver_get_current_iteration 
     procedure :: get_ideal_num_iterations         => hts_nonlinear_solver_get_ideal_num_iterations 
     procedure :: print_current_iteration_output   => hts_nonlinear_solver_print_current_iteration_output
     procedure :: print_final_output               => hts_nonlinear_solver_print_final_output 
     procedure :: free                             => hts_nonlinear_solver_free 
					! Getters 
					procedure :: get_residual                     => hts_nonlinear_solver_get_residual 
					procedure :: get_increment_dof_values         => hts_nonlinear_solver_get_increment_dof_values 
  end type hts_nonlinear_solver_t 
  
  public :: hts_nonlinear_solver_t
  
contains

subroutine hts_nonlinear_solver_create(this, convergence_criteria, abs_tol, rel_tol, max_iters, ideal_iters, & 
                                       fe_affine_operator, current_dof_values )
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
character(len=*)                      , intent(in)     :: convergence_criteria 
real(rp)                              , intent(in)     :: abs_tol
real(rp)                              , intent(in)     :: rel_tol 
integer(ip)                           , intent(in)     :: max_iters 
integer(ip)                           , intent(in)     :: ideal_iters 
type(fe_affine_operator_t) , target   , intent(in)     :: fe_affine_operator 
class(vector_t)            , target   , intent(in)     :: current_dof_values

! Initialize values 
this%current_iteration     = 0
this%convergence_criteria  = convergence_criteria 
this%absolute_tolerance    = abs_tol
this%relative_tolerance    = rel_tol
this%max_num_iterations = max_iters 
this%ideal_num_iterations  = ideal_iters 

! Point current dof values 
this%current_dof_values => current_dof_values 
this%current_rhs        => fe_affine_operator%get_translation() 
this%fe_affine_operator => fe_affine_operator 
this%current_lagrange_multiplier = 0.0_rp

! Create auxiliar structures 
call fe_affine_operator%create_range_vector(this%residual) 
call fe_affine_operator%create_range_vector(this%initial_residual)
call fe_affine_operator%create_range_vector(this%increment_dof_values) 
   
end subroutine hts_nonlinear_solver_create

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_start_new_iteration(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

assert(.not. this%finished()  )
this%current_iteration = this%current_iteration + 1
call this%line_search%initialize() 

end subroutine hts_nonlinear_solver_start_new_iteration 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_initialize(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

this%current_iteration = 0 
call this%compute_residual()
this%initial_residual   = this%residual 

end subroutine hts_nonlinear_solver_initialize  

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_residual(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

! Compute current residual 
call this%residual%init(0.0_rp)
call this%fe_affine_operator%apply( this%current_dof_values, this%residual ) 

end subroutine hts_nonlinear_solver_compute_residual 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_Jacobian(this, discrete_integration )
implicit none 
class(hts_nonlinear_solver_t)           , intent(inout)   :: this 
type(hts_CG_discrete_integration_t)     , intent(inout)   :: discrete_integration

! This subroutine is in charge of adding tangent terms to the current operator 
! J(u) = A(u) + A'(u)·u, stored in fe_affine_operator 
call discrete_integration%set_integration_type('add_tangent_terms')
call this%fe_affine_operator%numerical_setup() 
call discrete_integration%set_integration_type('regular')

end subroutine hts_nonlinear_solver_compute_Jacobian 

  ! ------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_update_solution(this) 
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)     :: this

this%current_dof_values          = this%current_dof_values + this%line_search%step_length*this%increment_dof_values
end subroutine hts_nonlinear_solver_update_solution 

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_converged(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical         :: hts_nonlinear_solver_converged

select case (this%convergence_criteria) 
case ('abs_res_norm') !  |R| < abs_tol & |dH|/|H| < rel_tol 
hts_nonlinear_solver_converged = ( (this%residual%nrm2() .lt. this%absolute_tolerance ) .and. &
                                   (this%increment_dof_values%nrm2()/this%current_dof_values%nrm2() .lt. this%relative_tolerance) )
case ('rel_r0_res_norm') ! |R|/|Ro| < rel_tol 
hts_nonlinear_solver_converged = (this%residual%nrm2()/this%initial_residual%nrm2() .lt. this%relative_tolerance ) 
case ('rel_rhs_res_norm') ! |R|/|b| < rel_tol 
hts_nonlinear_solver_converged = (this%residual%nrm2()/this%current_rhs%nrm2() .lt. this%relative_tolerance ) 
case DEFAULT
assert(.false.) 
end select 

end function hts_nonlinear_solver_converged

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_finished(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical                                                :: hts_nonlinear_solver_finished

hts_nonlinear_solver_finished = ( ( this%converged() .or. (this%current_iteration .gt. this%max_num_iterations) .or. (this%residual%nrm2()>1e15_rp) )  &
                                   .and. this%current_iteration .gt. 0 )

end function hts_nonlinear_solver_finished

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_get_current_iteration(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(in)     :: this 
integer(ip)                                            :: hts_nonlinear_solver_get_current_iteration 

hts_nonlinear_solver_get_current_iteration = this%current_iteration 

end function hts_nonlinear_solver_get_current_iteration 

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_get_ideal_num_iterations(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(in)     :: this 
integer(ip)                                            :: hts_nonlinear_solver_get_ideal_num_iterations

hts_nonlinear_solver_get_ideal_num_iterations = this%ideal_num_iterations

end function hts_nonlinear_solver_get_ideal_num_iterations 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_print_current_iteration_output(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

select case (this%convergence_criteria) 
case ('abs_res_norm')         !  |R| < abs_tol & |dH|/|H| < rel_tol 
write(6,'(a14,i3,a10,es21.15, a14, es21.15)')  'NL iteration ', this%current_iteration, '  |R| ', this%residual%nrm2(), &
                                               '  |dH|/|H| ', this%increment_dof_values%nrm2()/this%current_dof_values%nrm2()
case ('rel_r0_res_norm') ! |R|/|Ro| < rel_tol 
write(6,'(a14,i3,a16,es21.15, a10, es21.15)')  'NL iteration ', this%current_iteration, '  |R|/|Ro| ', this%residual%nrm2()/this%initial_residual%nrm2(), &
                                               '  |R| ', this%residual%nrm2()
case ('rel_rhs_res_norm') ! |R|/|b| < rel_tol 
write(6,'(a14,i3,a16,es21.15, a10, es21.15, a5, f8.4)')  'NL iteration ', this%current_iteration, '  |R|/|b| ', this%residual%nrm2()/this%current_rhs%nrm2(), &
                                               '  |R| ', this%residual%nrm2(), 'BT ', this%line_search%step_length
case DEFAULT
assert(.false.) 
end select 

end subroutine hts_nonlinear_solver_print_current_iteration_output

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_print_final_output(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

if ( this%converged() ) then 
write(*,*) ' ---- Newton NL solver converged in', this%current_iteration, 'iterations ---- ' 
else if (.not. this%converged() ) then 
 assert (this%finished() .eqv. .true.) 
 write(*,*) ' ---- Newton NL failed to converge after', this%current_iteration, 'iterations --- ' 
end if 

end subroutine hts_nonlinear_solver_print_final_output

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_free(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

if (allocated(this%residual)) call this%residual%free()
call this%line_search%free()
if (allocated(this%initial_residual)) call this%initial_residual%free()
if (allocated(this%increment_dof_values)) call this%increment_dof_values%free()

end subroutine hts_nonlinear_solver_free

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_get_residual(this) 
implicit none 
class(hts_nonlinear_solver_t) , target   , intent(in)     :: this
class(vector_t), pointer  :: hts_nonlinear_solver_get_residual

hts_nonlinear_solver_get_residual => this%residual 
end function hts_nonlinear_solver_get_residual

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_get_increment_dof_values(this) 
implicit none 
class(hts_nonlinear_solver_t), target , intent(in)     :: this
class(vector_t), pointer  :: hts_nonlinear_solver_get_increment_dof_values

hts_nonlinear_solver_get_increment_dof_values => this%increment_dof_values 
end function hts_nonlinear_solver_get_increment_dof_values

! **********************    LINE SEARCH SUBROUTINES ***************************************

subroutine line_search_initialize(this)  
  implicit none 
  class(hts_line_search_t)  , intent(inout)   :: this
  
  this%step_length     = 1.0_rp
  this%alpha           = 1.0e-4_rp 
  this%n               = 4 
  this%current_iterate = 0
  
  end subroutine line_search_initialize 
  
  subroutine line_search_cubic_backtracking(this, newton_raphson_solver ) 
  implicit none 
  class(hts_line_search_t)       , intent(inout)  :: this
  type(hts_nonlinear_solver_t)   , intent(inout)  :: newton_raphson_solver 
  
  real(rp)    :: s, mu, mu_residual, lambda_residual 
  real(rp)    :: lambda, lambda_0, lambda_q, lambda_s, lambda_c 
  real(rp)    :: t1, t2, a, b, d 
  integer(ip) :: i

  call newton_raphson_solver%fe_affine_operator%create_range_vector(this%initial_residual)
  call newton_raphson_solver%fe_affine_operator%create_range_vector(this%initial_dof_values)

      ! Line search Cubic Backtracking algorithm 
  call newton_raphson_solver%fe_affine_operator%numerical_setup() 
  call newton_raphson_solver%compute_residual()
  
  call this%initial_dof_values%copy( newton_raphson_solver%current_dof_values )  
  call this%initial_residual%copy(newton_raphson_solver%residual)
   s = -newton_raphson_solver%residual%nrm2()**2.0_rp 
   lambda_0 = this%step_length 
   
  newton_raphson_solver%current_dof_values = this%initial_dof_values + lambda_0*newton_raphson_solver%increment_dof_values  
  call newton_raphson_solver%fe_affine_operator%numerical_setup()
  call newton_raphson_solver%compute_residual()
 
  ! Check if linear is sufficient 
  if ( newton_raphson_solver%residual%nrm2()**2 .le. (this%initial_residual%nrm2()**2.0_rp + 2.0_rp*this%alpha*lambda_0*s) ) then 
    this%step_length     = lambda_0
    this%current_iterate = 1
    newton_raphson_solver%current_dof_values     = this%initial_dof_values 
    return 
  end if 
  
  mu          = lambda_0
  mu_residual = newton_raphson_solver%residual%nrm2()
  lambda_q = -s/(newton_raphson_solver%residual%nrm2()**2.0_rp - this%initial_residual%nrm2()**2.0_rp - 2.0_rp*lambda_0*s )
  
  ! Update lambda
  if ( lambda_q .lt. 0.1_rp*lambda_0 ) then 
  lambda = 0.1_rp*lambda_0
  elseif ( lambda_q .gt. 0.5_rp*lambda_0 ) then 
  lambda = 0.5_rp*lambda_0
  else
  lambda = lambda_q 
  end if 
  
  !  Update solution with new length step 
    newton_raphson_solver%current_dof_values = this%initial_dof_values + lambda*newton_raphson_solver%increment_dof_values 
    call newton_raphson_solver%fe_affine_operator%numerical_setup() 
    call newton_raphson_solver%compute_residual() 

  ! Check if quadratic is sufficient 
    if ( newton_raphson_solver%residual%nrm2()**2 .le. (this%initial_residual%nrm2()**2 + 2.0_rp*this%alpha*lambda*s) ) then 
    this%step_length     = lambda
    this%current_iterate = 2
    newton_raphson_solver%current_dof_values     = this%initial_dof_values 
    return 
  end if 

      lambda_residual = newton_raphson_solver%residual%nrm2()
  do i = 1,this%n
  
  t1 = 0.5_rp* ( mu_residual**2.0_rp     - this%initial_residual%nrm2()**2.0_rp - lambda*s ) 
  t2 = 0.5_rp* ( lambda_residual**2.0_rp - this%initial_residual%nrm2()**2.0_rp - lambda*s ) 

  a = (t1/lambda**2.0_rp - t2/mu**2.0_rp)/(lambda - mu)
  b = ( lambda*t2/mu**2.0_rp - lambda*t1/lambda**2.0_rp )/(lambda - mu)
  d = max ( b**2.0_rp - 3.0_rp*a*s, 0.0_rp ) 
  
  if ( abs(a) .gt. 1e-10 ) then 
  lambda_c = (sqrt(d)-b)/(3.0_rp*a)
  else 
  lambda_c = -s/(2.0_rp*b)
  end if 
 
    mu = lambda
    mu_residual = lambda_residual
  
    ! Update lambda
  if ( lambda_c .lt. 0.1_rp*lambda ) then 
  lambda = 0.1_rp*lambda
  elseif ( lambda_c .gt. 0.5_rp*lambda ) then 
  lambda = 0.5_rp*lambda 
  else
  lambda = lambda_c 
  end if 
  
    newton_raphson_solver%current_dof_values = this%initial_dof_values + lambda*newton_raphson_solver%increment_dof_values 
    call newton_raphson_solver%fe_affine_operator%numerical_setup() 
    call newton_raphson_solver%compute_residual() 
    ! Check if cubic is sufficient 
    if ( newton_raphson_solver%residual%nrm2()**2 .le. ( this%initial_residual%nrm2()**2 + 2.0_rp*this%alpha*lambda*s)  ) then 
        this%step_length = lambda
        this%current_iterate = 2 + i
        newton_raphson_solver%current_dof_values     = this%initial_dof_values 
    return 
  end if 
  lambda_residual = newton_raphson_solver%residual%nrm2()

  end do 
  
  this%step_length                             = lambda 
  newton_raphson_solver%current_dof_values     = this%initial_dof_values 
  
  end subroutine line_search_cubic_backtracking 
  
    subroutine line_search_free(this)  
  implicit none 
  class(hts_line_search_t)  , intent(inout)   :: this
  
  if ( allocated(this%initial_dof_values) ) call this%initial_dof_values%free()
  if ( allocated(this%initial_residual)   ) call this%initial_residual%free() 
  
  end subroutine line_search_free
  
end module hts_nonlinear_solver_names 
