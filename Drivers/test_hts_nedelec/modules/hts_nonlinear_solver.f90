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
  use base_sparse_matrix_names 
  implicit none
# include "debug.i90"
  private
  
  type :: hts_nonlinear_solver_t
     private
     integer(ip)                                :: current_iteration 
     integer(ip)                                :: ideal_num_iterations 
     integer(ip)                                :: max_number_iterations
	    real(rp)                                   :: absolute_tolerance
     real(rp)                                   :: relative_tolerance
     real(rp)                                   :: step_length 
     character(len=:) , allocatable             :: convergence_criteria 
     type(fe_affine_operator_t) , pointer       :: fe_affine_operator 
     class(vector_t), pointer                   :: current_dof_values 
     real(rp)                                   :: current_lagrange_multiplier 
     class(vector_t), pointer                   :: current_rhs 
     class(vector_t), allocatable               :: increment_dof_values 
     real(rp)                                   :: increment_lagrange_multiplier 
	    class(vector_t), allocatable               :: residual 
     type(serial_scalar_array_t)                :: constrained_residual 
     class(vector_t), allocatable               :: initial_residual 
     type(direct_solver_t)                      :: direct_solver
     logical                                    :: apply_current_constraint 
    contains
     procedure :: create                           => hts_nonlinear_solver_create 
     procedure :: start_new_iteration              => hts_nonlinear_solver_start_new_iteration 
     procedure :: initialize                       => hts_nonlinear_solver_initialize 
     procedure :: compute_residual                 => hts_nonlinear_solver_compute_residual 
     procedure :: compute_constrained_residual     => hts_nonlinear_solver_compute_constrained_residual 
     procedure :: compute_Jacobian                 => hts_nonlinear_solver_compute_Jacobian 
     procedure :: solve_tangent_system             => hts_nonlinear_solver_solve_tangent_system 
     procedure :: solve_constrained_tangent_system => hts_nonlinear_solver_solve_constrained_tangent_system 
     procedure :: update_solution                  => hts_nonlinear_solver_update_solution
     procedure :: converged                        => hts_nonlinear_solver_converged 
     procedure :: finished                         => hts_nonlinear_solver_finished 
     procedure :: get_current_iteration            => hts_nonlinear_solver_get_current_iteration 
     procedure :: get_ideal_num_iterations         => hts_nonlinear_solver_get_ideal_num_iterations 
     procedure :: get_apply_current_constraint     => hts_nonlinear_solver_get_apply_current_constraint 
     procedure :: print_current_iteration_output   => hts_nonlinear_solver_print_current_iteration_output
     procedure :: print_final_output               => hts_nonlinear_solver_print_final_output 
     procedure :: free                             => hts_nonlinear_solver_free 
  end type hts_nonlinear_solver_t 
  
  public :: hts_nonlinear_solver_t
  
contains

subroutine hts_nonlinear_solver_create(this, convergence_criteria, abs_tol, rel_tol, max_iters, ideal_iters, fe_affine_operator, current_dof_values, apply_constraint)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
character(len=*)                      , intent(in)     :: convergence_criteria 
real(rp)                              , intent(in)     :: abs_tol
real(rp)                              , intent(in)     :: rel_tol 
integer(ip)                           , intent(in)     :: max_iters 
integer(ip)                           , intent(in)     :: ideal_iters 
type(fe_affine_operator_t) , target   , intent(in)     :: fe_affine_operator 
class(vector_t)            , target   , intent(in)     :: current_dof_values
logical                               , intent(in)     :: apply_constraint 
integer                      :: FPLError
type(parameterlist_t)        :: parameter_list
integer                      :: iparm(64)
class(matrix_t), pointer     :: matrix

! Initialize values 
this%current_iteration     = 0
this%convergence_criteria  = convergence_criteria 
this%absolute_tolerance    = abs_tol
this%relative_tolerance    = rel_tol
this%max_number_iterations = max_iters 
this%ideal_num_iterations  = ideal_iters 
this%step_length           = 1.0_rp  
this%apply_current_constraint = apply_constraint 

! Point current dof values 
this%current_dof_values => current_dof_values 
this%current_rhs        => fe_affine_operator%get_translation() 
this%fe_affine_operator => fe_affine_operator 
this%current_lagrange_multiplier = 0.0_rp

! Create auxiliar structures 
call fe_affine_operator%create_range_vector(this%residual) 
call fe_affine_operator%create_range_vector(this%initial_residual)
call fe_affine_operator%create_range_vector(this%increment_dof_values) 

    select type(residual=>this%residual)
    class is (serial_scalar_array_t)  
     call this%constrained_residual%create_and_allocate(residual%get_size() + 1)
    class DEFAULT
       assert(.false.) 
    end select

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
call this%compute_residual()
this%initial_residual  = this%residual 

end subroutine hts_nonlinear_solver_initialize  

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_residual(this)
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 

call this%residual%init(0.0_rp)
call this%fe_affine_operator%apply( this%current_dof_values, this%residual ) 

end subroutine hts_nonlinear_solver_compute_residual 

! -------------------------------------------------------------------------------------------------
  subroutine hts_nonlinear_solver_compute_constrained_residual(this, constraint_vector, constraint_value ) 
  implicit none 
  class(hts_nonlinear_solver_t), intent(inout)    :: this
  type(serial_scalar_array_t) ,  intent(in)       :: constraint_vector 
  real(rp)                    ,  intent(in)       :: constraint_value 
  ! Locals 
    integer(ip)                                      :: num_dofs 
    type(serial_scalar_array_t)      , pointer       :: original_residual 
    real(rp)                         , allocatable   :: constrained_residual_entries(:) 
    class(vector_t)                  , pointer       :: dof_values_current
    type(serial_scalar_array_t)      , pointer       :: original_sol
    integer(ip)                      , allocatable   :: original_indices(:)
    integer(ip) :: idof 
	 
       select type (residual=>this%residual)
    type is (serial_scalar_array_t)
       original_residual => residual
       num_dofs = residual%get_size() 
    class DEFAULT
       assert(.false.) 
    end select
	
	   select type (current_dof_values=>this%current_dof_values)
    type is (serial_scalar_array_t)
       original_sol => current_dof_values 
    class DEFAULT
       assert(.false.) 
    end select
	
    call this%constrained_residual%init(0.0_rp) 
    call memalloc ( num_dofs, original_indices, __FILE__,__LINE__) 
    do idof = 1,num_dofs
    original_indices(idof) = idof 
    end do 
    
  ! Build full constrained residual: 
	   !       [ R - lag_mult · C^t ]
    !  R* = [ J_o - C · H          ]
    		    
    call this%constrained_residual%init(0.0_rp) 
    call this%constrained_residual%insert_subvector(iblock=1,                          &
                                                    size_indices = num_dofs,           &
                                                    indices = original_indices,        &
                                                    values = original_residual%get_entries() + this%current_lagrange_multiplier*constraint_vector%get_entries() )
    call this%constrained_residual%insert(num_dofs+1, constraint_vector%dot(original_sol)  - constraint_value )
    call memfree( original_indices, __FILE__, __LINE__ )

  end subroutine hts_nonlinear_solver_compute_constrained_residual 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_compute_Jacobian(this, discrete_integration )
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)   :: this 
class(hts_nedelec_discrete_integration_t)   , intent(inout)   :: discrete_integration

! This subroutine is in charge of adding tangent terms to the current operator 
! J(u) = A(u) + A'(u)·u, stored in fe_affine_operator 
call discrete_integration%set_integration_type('add_tangent_terms')
call this%fe_affine_operator%numerical_setup() 
call discrete_integration%set_integration_type('regular')

end subroutine hts_nonlinear_solver_compute_Jacobian 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_solve_tangent_system(this)
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)     :: this  
class(matrix_t), pointer       :: Jacobian

Jacobian => this%fe_affine_operator%get_matrix() 
 
select type (Jacobian) 
class is (sparse_matrix_t) 
  call this%direct_solver%update_matrix( Jacobian, same_nonzero_pattern=.false.) 
  call this%direct_solver%solve( -this%residual , this%increment_dof_values )
class DEFAULT 
assert(.false.) 
end select 

end subroutine hts_nonlinear_solver_solve_tangent_system   

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_solve_constrained_tangent_system(this, constraint_matrix, constraint_vector, constraint_value )
  implicit none
  class(hts_nonlinear_solver_t) , intent(inout)  :: this
  type(coo_sparse_matrix_t)     , intent(in)     :: constraint_matrix 
  type(serial_scalar_array_t)   , intent(in)     :: constraint_vector 
  real(rp)                      , intent(in)     :: constraint_value 
  
  ! Locals 
  class(matrix_t)                  , pointer       :: matrix

  type(sparse_matrix_t)      , pointer :: original_coefficient_matrix
  type(sparse_matrix_t)                :: constrained_coefficient_matrix

  type(serial_scalar_array_t), pointer :: original_residual
  type(serial_scalar_array_t), pointer :: original_sol
  type(serial_scalar_array_t)          :: minus_constrained_residual
  type(serial_scalar_array_t)          :: constrained_sol
  type(serial_scalar_array_t)          :: constraints_vector  
  real(rp), allocatable                :: minus_constrained_residual_entries(:)
  real(rp), allocatable                :: constrained_sol_entries(:)
  integer(ip)                          :: num_rows
  integer(ip), allocatable             :: original_indices(:) 
  integer(ip)                          :: FPLError, idof 
  type(parameterlist_t)                :: parameter_list
  integer(ip)                          :: iparm(64)
  type(direct_solver_t)                :: direct_solver

  matrix  => this%fe_affine_operator%get_matrix()
  this%increment_lagrange_multiplier = 0.0_rp 

  select type (matrix)
  type is (sparse_matrix_t)
     original_coefficient_matrix => matrix
     class DEFAULT
     assert(.false.) 
  end select

  select type (residual=>this%residual)
  type is (serial_scalar_array_t)
     original_residual => residual
     class DEFAULT
     assert(.false.) 
  end select

  select type (current_dof_values => this%current_dof_values)
  type is (serial_scalar_array_t)
     original_sol => current_dof_values 
     class DEFAULT
     assert(.false.) 
  end select

  ! Expand original_coefficient_matrix with this%constraint_matrix, and
  ! store the result in constrained_coefficient_matrix 
  call original_coefficient_matrix%expand_matrix_numeric(C_T = constraint_matrix, &
                                                    to  = constrained_coefficient_matrix)
  num_rows = original_coefficient_matrix%get_num_rows()
  call memalloc( num_rows, original_indices, __FILE__, __LINE__ )
  do idof=1,num_rows
     original_indices(idof) = idof
  end do

  ! Set-up constrained Residual  
  !       [ R - lag_mult · C^t ]
  !  R⁺ = [ J - C · H          ]

  call minus_constrained_residual%create(num_rows+1)
  call memalloc ( num_rows+1, minus_constrained_residual_entries, __FILE__,__LINE__)   
  minus_constrained_residual_entries(1:num_rows)  = - original_residual%get_entries() - this%current_lagrange_multiplier*constraint_vector%get_entries()   
  minus_constrained_residual_entries(num_rows+1)  = - constraint_vector%dot(original_sol) + constraint_value ! To Add I(t) 
  call minus_constrained_residual%set_view_entries( minus_constrained_residual_entries(1:num_rows+1) )

  ! Set-up constrained sol
  call constrained_sol%create(num_rows+1)
  call memalloc ( num_rows+1, constrained_sol_entries, __FILE__,__LINE__)  
  call constrained_sol%set_view_entries(constrained_sol_entries)

  ! The following lines involving the set_type + set_parameters MUST BE removed
  ! in the final implementation
  call parameter_list%init()
  FPLError = 0
  FPLError = FPLError + parameter_list%set(key = direct_solver_type,        value = pardiso_mkl)
  FPLError = FPLError + parameter_list%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_uss)
  FPLError = FPLError + parameter_list%set(key = pardiso_mkl_message_level, value = 0)
  iparm = 0
  FPLError = FPLError + parameter_list%set(key = pardiso_mkl_iparm,         value = iparm)
  check(FPLError == 0)

  call direct_solver%set_type_from_pl(parameter_list)
  call direct_solver%set_parameters_from_pl(parameter_list)

  call parameter_list%free()

  call direct_solver%update_matrix(constrained_coefficient_matrix, same_nonzero_pattern=.false.)
  call direct_solver%solve( minus_constrained_residual, constrained_sol )
  call direct_solver%free()
  
  ! Insert constrained solution DOFs into original vector, get increment lagrange multiplier 
  constrained_sol_entries = constrained_sol%get_entries() 
  call this%increment_dof_values%insert_subvector( 1, num_rows, original_indices, constrained_sol_entries(1:num_rows) ) 
  this%increment_lagrange_multiplier = constrained_sol_entries(num_rows+1) 

  call memfree ( minus_constrained_residual_entries, __FILE__,__LINE__)  
  call memfree ( constrained_sol_entries,            __FILE__,__LINE__) 
  call memfree ( original_indices       ,            __FILE__,__LINE__)

  call constrained_coefficient_matrix%free()
end subroutine hts_nonlinear_solver_solve_constrained_tangent_system 
  
  ! ------------------------------------------------------------------------------------

subroutine hts_nonlinear_solver_update_solution(this, current_solution) 
implicit none 
class(hts_nonlinear_solver_t)               , intent(inout)     :: this
class(fe_function_t)                        , intent(inout)     :: current_solution
class(vector_t), pointer   :: current_dof_values

this%current_dof_values          = this%current_dof_values + this%step_length*this%increment_dof_values
this%current_lagrange_multiplier = this%current_lagrange_multiplier + this%step_length*this%increment_lagrange_multiplier 

end subroutine hts_nonlinear_solver_update_solution 

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_converged(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical         :: hts_nonlinear_solver_converged

if ( .not. this%get_apply_current_constraint() ) then 

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

else 

select case (this%convergence_criteria) 
case ('abs_res_norm') !  |R*| < abs_tol & |dH|/|H| < rel_tol 
hts_nonlinear_solver_converged = ( (this%constrained_residual%nrm2() .lt. this%absolute_tolerance ) )
case ('rel_r0_res_norm') ! |R*|/|Ro*| < rel_tol 
assert(.false.) 
case ('rel_rhs_res_norm') ! |R*|/|b*| < rel_tol 
hts_nonlinear_solver_converged = (this%constrained_residual%nrm2()/this%current_rhs%nrm2() .lt. this%relative_tolerance ) 
case DEFAULT
assert(.false.) 
end select 

end if 

end function hts_nonlinear_solver_converged

! ---------------------------------------------------------------------------------------
function hts_nonlinear_solver_finished(this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
logical                                                :: hts_nonlinear_solver_finished

hts_nonlinear_solver_finished = ( ( this%converged() .or. (this%current_iteration .gt. this%max_number_iterations) .or. (this%residual%nrm2()>1e15_rp) )  &
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
function hts_nonlinear_solver_get_apply_current_constraint (this) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(in)     :: this 
logical                                                :: hts_nonlinear_solver_get_apply_current_constraint

hts_nonlinear_solver_get_apply_current_constraint = this%apply_current_constraint

end function hts_nonlinear_solver_get_apply_current_constraint 

! ---------------------------------------------------------------------------------------
subroutine hts_nonlinear_solver_print_current_iteration_output(this, constraint_value) 
implicit none 
class(hts_nonlinear_solver_t)         , intent(inout)  :: this 
real(rp) , optional                   , intent(in)     :: constraint_value 

real(rp) :: constrained_rhs_norm 

  ! Screen output
if ( .not. this%get_apply_current_constraint() ) then 

select case (this%convergence_criteria) 
case ('abs_res_norm')         !  |R| < abs_tol & |dH|/|H| < rel_tol 
write(6,'(a14,i3,a10,es21.15, a14, es21.15)')  'NL iteration ', this%current_iteration, '  |R| ', this%residual%nrm2(), &
                                               '  |dH|/|H| ', this%increment_dof_values%nrm2()/this%current_dof_values%nrm2()
case ('rel_r0_res_norm') ! |R|/|Ro| < rel_tol 
write(6,'(a14,i3,a16,es21.15, a10, es21.15)')  'NL iteration ', this%current_iteration, '  |R|/|Ro| ', this%residual%nrm2()/this%initial_residual%nrm2(), &
                                               '  |R| ', this%residual%nrm2()
case ('rel_rhs_res_norm') ! |R|/|b| < rel_tol 
write(6,'(a14,i3,a16,es21.15, a10, es21.15)')  'NL iteration ', this%current_iteration, '  |R|/|b| ', this%residual%nrm2()/this%current_rhs%nrm2(), &
                                               '  |R| ', this%residual%nrm2()
case DEFAULT
assert(.false.) 
end select 

else 

select case (this%convergence_criteria) 
case ('abs_res_norm')         !  |R*| < abs_tol 
write(6,'(a14,i3,a10,es21.15, a10, es21.15)')  'NL iteration ', this%current_iteration, '  |R*| ', this%constrained_residual%nrm2(), ' |R| ', this%residual%nrm2() 
case ('rel_r0_res_norm') ! |R*|/|Ro*| < rel_tol 
assert(.false.) ! Not implemented 
case ('rel_rhs_res_norm') ! |R*|/|b*| < rel_tol 
constrained_rhs_norm = sqrt(this%current_rhs%nrm2()**2 + constraint_value**2)
write(6,'(a14,i3,a16,es21.15, a10, es21.15)')  'NL iteration ', this%current_iteration, '  |R*|/|b*| ', this%constrained_residual%nrm2()/constrained_rhs_norm, &
                                               '  |R*| ', this%constrained_residual%nrm2()
case DEFAULT
assert(.false.) 
end select 

end if 

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
call this%constrained_residual%free()
if (allocated(this%initial_residual)) call this%initial_residual%free()
if (allocated(this%increment_dof_values)) call this%increment_dof_values%free()
 
! call this%direct_solver%free() 

end subroutine hts_nonlinear_solver_free



end module hts_nonlinear_solver_names 
