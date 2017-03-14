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

module poisson_unfitted_analytical_functions_names
  use fempar_names
  use poisson_unfitted_exact_solutions_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: base_scalar_function_t
    class(scalar_exact_solution_t), pointer :: exact_solution => null()
    contains
      procedure, non_overridable :: create => base_scalar_function_create
      procedure, non_overridable :: free   => base_scalar_function_free
  end type base_scalar_function_t

  type, extends(base_scalar_function_t) :: source_term_t
    private
    contains
      procedure :: get_value_space    => source_term_get_value_space
  end type source_term_t

  type, extends(base_scalar_function_t) :: boundary_function_t
    private
    contains
      procedure :: get_value_space => boundary_function_get_value_space
  end type boundary_function_t

  type, extends(base_scalar_function_t) :: solution_function_t
    private
    contains
      procedure :: get_value_space    => solution_function_get_value_space
      procedure :: get_gradient_space => solution_function_get_gradient_space
  end type solution_function_t

  type poisson_unfitted_analytical_functions_t
    private
    class(scalar_exact_solution_t), pointer :: exact_solution => null()
    type(source_term_t)       :: source_term
    type(boundary_function_t) :: boundary_function
    type(solution_function_t) :: solution_function
    contains
      procedure :: create                  => poisson_unfitted_analytical_functions_create
      procedure :: free                    => poisson_unfitted_analytical_functions_free
      procedure :: get_source_term         => poisson_unfitted_analytical_functions_get_source_term
      procedure :: get_boundary_function   => poisson_unfitted_analytical_functions_get_boundary_function
      procedure :: get_solution_function   => poisson_unfitted_analytical_functions_get_solution_function
  end type poisson_unfitted_analytical_functions_t

  public :: poisson_unfitted_analytical_functions_t, base_scalar_function_t

contains  

  !===============================================================================================
  subroutine base_scalar_function_create ( this, exact_solution )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    class(scalar_exact_solution_t), pointer, intent(in) :: exact_solution
    this%exact_solution => exact_solution
  end subroutine base_scalar_function_create

  !===============================================================================================
  subroutine base_scalar_function_free ( this )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    this%exact_solution => null()
  end subroutine base_scalar_function_free

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    call this%exact_solution%lapl_u(point,result)
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine boundary_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    call this%exact_solution%u(point,result)
  end subroutine boundary_function_get_value_space

  !===============================================================================================
  subroutine solution_function_get_value_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    call this%exact_solution%u(point,result)
  end subroutine solution_function_get_value_space

  !===============================================================================================
  subroutine solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    call this%exact_solution%grad_u(point,result)
  end subroutine solution_function_get_gradient_space

  !===============================================================================================
  subroutine poisson_unfitted_analytical_functions_create ( this, solution_name)
    implicit none
    class(poisson_unfitted_analytical_functions_t), intent(inout)    :: this
    character(len=*), intent(in) :: solution_name
    this%exact_solution => scalar_exact_solution_t(solution_name)
    call this%source_term%create(this%exact_solution)
    call this%boundary_function%create(this%exact_solution)
    call this%solution_function%create(this%exact_solution)
  end subroutine poisson_unfitted_analytical_functions_create

  !===============================================================================================
  subroutine poisson_unfitted_analytical_functions_free ( this, solution_name)
    implicit none
    class(poisson_unfitted_analytical_functions_t), intent(inout)    :: this
    character(len=*), intent(in) :: solution_name
    this%exact_solution => null()
    call this%source_term%free()
    call this%boundary_function%free()
    call this%solution_function%free()
  end subroutine poisson_unfitted_analytical_functions_free

  !===============================================================================================
  function poisson_unfitted_analytical_functions_get_source_term ( this )
    implicit none
    class(poisson_unfitted_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_unfitted_analytical_functions_get_source_term
    poisson_unfitted_analytical_functions_get_source_term => this%source_term
  end function poisson_unfitted_analytical_functions_get_source_term

  !===============================================================================================
  function poisson_unfitted_analytical_functions_get_boundary_function ( this )
    implicit none
    class(poisson_unfitted_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_unfitted_analytical_functions_get_boundary_function
    poisson_unfitted_analytical_functions_get_boundary_function => this%boundary_function
  end function poisson_unfitted_analytical_functions_get_boundary_function

  !===============================================================================================
  function poisson_unfitted_analytical_functions_get_solution_function ( this )
    implicit none
    class(poisson_unfitted_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_unfitted_analytical_functions_get_solution_function
    poisson_unfitted_analytical_functions_get_solution_function => this%solution_function
  end function poisson_unfitted_analytical_functions_get_solution_function

end module poisson_unfitted_analytical_functions_names
!***************************************************************************************************
