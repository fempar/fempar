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
  use hts_nedelec_discrete_integration_names
  implicit none
# include "debug.i90"
  private
  
  type :: hts_nonlinear_solver_t
     private
     integer(ip)                     :: current_iteration 
     integer(ip)                     :: ideal_num_iterations 
     integer(ip)                     :: max_number_iterations
	    real(rp)                        :: absolute_tolerance
     real(rp)                        :: relative_tolerance
     real(rp)                        :: step_length 
     class(vector_t), pointer        :: current_dof_values 
     class(vector_t), allocatable    :: increment_dof_values 
	    class(vector_t), allocatable    :: residual 
     type(direct_solver_t)           :: direct_solver
   contains
     procedure :: create                          => hts_nonlinear_solver_create 
     procedure :: start_new_iteration             => hts_nonlinear_solver_start_new_iteration 
     procedure :: initialize                      => hts_nonlinear_solver_initialize 
     procedure :: compute_residual                => hts_nonlinear_solver_compute_residual 
     procedure :: compute_Jacobian                => hts_nonlinear_solver_compute_Jacobian 
     procedure :: solve_tangent_system            => hts_nonlinear_solver_solve_tangent_system 
     procedure :: update_solution                 => hts_nonlinear_solver_update_solution
     procedure :: converged                       => hts_nonlinear_solver_converged 
     procedure :: finished                        => hts_nonlinear_solver_finished 
     procedure :: get_current_iteration           => hts_nonlinear_solver_get_current_iteration 
     procedure :: get_ideal_num_iterations        => hts_nonlinear_solver_get_ideal_num_iterations 
     procedure :: print_current_iteration_output  => hts_nonlinear_solver_print_current_iteration_output
     procedure :: print_final_output              => hts_nonlinear_solver_print_final_output 
     procedure :: free                            => hts_nonlinear_solver_free 
  end type hts_nonlinear_solver_t 
  
  public :: hts_nonlinear_solver_t
  
contains

subroutine hts_nonlinear_solver_create(this, abs_tol, rel_tol, max_iters, ideal_iters, fe_affine_operator, current_dof_values)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
real(rp)                              , intent(in)     :: abs_tol
real(rp)                              , intent(in)     :: rel_tol 
integer(ip)                           , intent(in)     :: max_iters 
integer(ip)                           , intent(in)     :: ideal_iters 
type(fe_affine_operator_t)            , intent(in)     :: fe_affine_operator 
class(vector_t)            , target   , intent(in)     :: current_dof_values 
integer                      :: FPLError
type(parameterlist_t)        :: parameter_list
integer                      :: iparm(64)
class(matrix_t), pointer     :: matrix

! Initialize values 
this%current_iteration     = 0
this%absolute_tolerance    = abs_tol
this%relative_tolerance    = rel_tol
this%max_number_iterations = max_iters 
this%ideal_num_iterations  = ideal_iters 
this%step_length           = 1.0_rp 

! Point current dof values 
this%current_dof_values => current_dof_values 

! Create auxiliar structures 
call fe_affine_operator%create_range_vector(this%residual) 
call fe_affine_operator%create_range_vector(this%increment_dof_values) 

! Create direct solver to update iterates 
    call parameter_list%init()
   
    FPLError =            parameter_list%set(key = direct_solver_type     ,   value = pardiso_mkl)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_uns)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_message_level, value = 0)
    iparm = 0
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_iparm,         value = iparm)
    assert(FPLError == 0)
    
    call this%direct_solver%set_type_from_pl(parameter_list)
    call this%direct_solver%set_parameters_from_pl(parameter_list)
    
     matrix => fe_affine_operator%get_matrix()
    select type(matrix)
    class is (sparse_matrix_t)  
       call this%direct_solver%set_matrix(matrix)
    class DEFAULT
       assert(.false.) 
    end select
  
    call parameter_list%free()
    
end subroutine hts_nonlinear_solver_create

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_start_new_iteration(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

assert(.not. this%finished()  )
this%current_iteration = this%current_iteration + 1
end subroutine hts_nonlinear_solver_start_new_iteration 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_initialize(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

this%current_iteration = 0
end subroutine hts_nonlinear_solver_initialize  

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_residual(this, fe_affine_operator, current_dof_values )
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
type(fe_affine_operator_t)            , intent(in)     :: fe_affine_operator
class(vector_t)                       , intent(in)     :: current_dof_values 

call this%residual%init(0.0_rp)
call fe_affine_operator%apply( current_dof_values, this%residual ) 

end subroutine hts_nonlinear_solver_compute_residual 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_Jacobian(this, discrete_integration, fe_affine_operator )
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)   :: this 
class(hts_nedelec_discrete_integration_t)   , intent(inout)   :: discrete_integration
type(fe_affine_operator_t)                  , intent(inout)   :: fe_affine_operator

! This subroutine is in charge of adding tangent terms to the current operator 
! J(u) = A(u) + A'(u)·u, stored in fe_affine_operator 
call discrete_integration%set_integration_type('add_tangent_terms')
call fe_affine_operator%numerical_setup() 
call discrete_integration%set_integration_type('regular')

end subroutine hts_nonlinear_solver_compute_Jacobian 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_solve_tangent_system(this, fe_affine_operator )
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)     :: this  
type(fe_affine_operator_t)                  , intent(in)        :: fe_affine_operator
class(vector_t), pointer       :: current_dof_values 
class(matrix_t), pointer       :: Jacobian
class(vector_t), pointer       :: residual 

Jacobian => fe_affine_operator%get_matrix() 
 
select type (Jacobian) 
class is (sparse_matrix_t) 
  call this%direct_solver%update_matrix( Jacobian, same_nonzero_pattern=.false.) 
  call this%direct_solver%solve( -this%residual , this%increment_dof_values )
class DEFAULT 
assert(.false.) 
end select 

end subroutine hts_nonlinear_solver_solve_tangent_system   

subroutine hts_nonlinear_solver_update_solution(this, current_solution) 
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)     :: this
class(fe_function_t)                        , intent(inout)     :: current_solution
class(vector_t), pointer   :: current_dof_values

current_dof_values => current_solution%get_dof_values()
current_dof_values = current_dof_values + this%step_length*this%increment_dof_values

end subroutine hts_nonlinear_solver_update_solution 

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_converged(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical         :: hts_nonlinear_solver_converged

hts_nonlinear_solver_converged = ( (this%residual%nrm2() .lt. this%absolute_tolerance ) .and. &
                                   (this%increment_dof_values%nrm2()/this%current_dof_values%nrm2() .lt. this%relative_tolerance) )

end function hts_nonlinear_solver_converged

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_finished(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical                                                :: hts_nonlinear_solver_finished

hts_nonlinear_solver_finished = ( ( this%converged() .or. (this%current_iteration .gt. this%max_number_iterations) .or. (this%residual%nrm2()>1e12_rp) )  &
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

  ! Screen output 
write(6,'(a14,i3,a10,es21.15, a14, es21.15)')  'NL iteration ', this%current_iteration, '  |R| ', this%residual%nrm2(), &
                                               '  |dH|/|H| ', this%increment_dof_values%nrm2()/this%current_dof_values%nrm2()

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
if (allocated(this%increment_dof_values)) call this%increment_dof_values%free()
! call this%direct_solver%free() 

end subroutine hts_nonlinear_solver_free



end module hts_nonlinear_solver_names 
