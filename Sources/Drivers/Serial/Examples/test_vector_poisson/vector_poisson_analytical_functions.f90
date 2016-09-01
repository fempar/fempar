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

module vector_poisson_analytical_functions_names
  use serial_names
  implicit none
# include "debug.i90"
  private

  type, extends(vector_function_t) :: source_term_t
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t

  type, extends(scalar_function_t) :: boundary_values_t
   contains
     procedure :: get_value_space => boundary_values_get_value_space
  end type boundary_values_t

  type, extends(vector_function_t) :: solution_values_t
   contains
     procedure :: get_value_space => solution_values_get_value_space
  end type solution_values_t
  
  type, extends(tensor_function_t) :: solution_gradient_t
   contains
     procedure :: get_value_space => solution_gradient_get_value_space
  end type solution_gradient_t

  type vector_poisson_analytical_functions_t
     private
     type(source_term_t)       :: source_term
     type(boundary_values_t)   :: boundary_values
     type(solution_values_t)   :: solution_values
     type(solution_gradient_t) :: solution_gradient
   contains
     procedure :: get_source_term       => poisson_analytical_functions_get_source_term
     procedure :: get_boundary_values   => poisson_analytical_functions_get_boundary_values
     procedure :: get_solution_values   => poisson_analytical_functions_get_solution_values
     procedure :: get_solution_gradient => poisson_analytical_functions_get_solution_gradient
  end type vector_poisson_analytical_functions_t

  public :: vector_poisson_analytical_functions_t

contains  

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    call result%set(1,-2.0_rp)
    call result%set(2,-2.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine boundary_values_get_value_space ( this, point, result )
    implicit none
    class(boundary_values_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    result = point%get(1)**2 !sin ( pi * point%get(1) ) * sin ( pi * point%get(2) ) + point%get(1)
  end subroutine boundary_values_get_value_space

  !===============================================================================================
  subroutine solution_values_get_value_space ( this, point, result )
    implicit none
    class(solution_values_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    call result%set(1, point%get(1)**2 ) !result%set(1, sin ( pi * point%get(1) ) * sin ( pi * point%get(2) ) + point%get(1) )
    call result%set(2, point%get(1)**2 ) !result%set(2, sin ( pi * point%get(1) ) * sin ( pi * point%get(2) ) + point%get(1) )
  end subroutine solution_values_get_value_space
  
  !===============================================================================================
  subroutine solution_gradient_get_value_space ( this, point, result )
    implicit none
    class(solution_gradient_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    call result%set( 1, 1, 0.0_rp ) 
    call result%set( 2, 1, 0.0_rp )
    call result%set( 1, 2, 0.0_rp ) 
    call result%set( 2, 2, 0.0_rp )
    !call result%set( 1, 1, pi * cos ( pi * point%get(1) ) * sin ( pi * point%get(2) ) + 1.0_rp ) 
    !call result%set( 2, 1, pi * sin ( pi * point%get(1) ) * cos ( pi * point%get(2) ) )
    !call result%set( 1, 2, pi * cos ( pi * point%get(1) ) * sin ( pi * point%get(2) ) + 1.0_rp ) 
    !call result%set( 2, 2, pi * sin ( pi * point%get(1) ) * cos ( pi * point%get(2) ) )
  end subroutine solution_gradient_get_value_space
  
  !===============================================================================================
  function poisson_analytical_functions_get_source_term ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: poisson_analytical_functions_get_source_term
    poisson_analytical_functions_get_source_term => this%source_term
  end function poisson_analytical_functions_get_source_term
  
  !===============================================================================================
  function poisson_analytical_functions_get_boundary_values ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_analytical_functions_get_boundary_values
    poisson_analytical_functions_get_boundary_values => this%boundary_values
  end function poisson_analytical_functions_get_boundary_values
  
  !===============================================================================================
  function poisson_analytical_functions_get_solution_values ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: poisson_analytical_functions_get_solution_values
    poisson_analytical_functions_get_solution_values => this%solution_values
  end function poisson_analytical_functions_get_solution_values
  
  !===============================================================================================
  function poisson_analytical_functions_get_solution_gradient ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(tensor_function_t), pointer :: poisson_analytical_functions_get_solution_gradient
    poisson_analytical_functions_get_solution_gradient => this%solution_gradient
  end function poisson_analytical_functions_get_solution_gradient

end module vector_poisson_analytical_functions_names
!***************************************************************************************************


