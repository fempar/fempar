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

module poisson_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims = -1  
  contains
    procedure :: set_num_dims    => base_scalar_function_set_num_dims 
  end type base_scalar_function_t
  
    type, extends(base_scalar_function_t) :: source_term_t
    private 
   contains
     procedure :: get_value_space       => source_term_get_value_space
     procedure :: get_value_space_time  => source_term_get_value_space_time
  end type source_term_t

  type, extends(base_scalar_function_t) :: boundary_function_t
    private
   contains
     procedure :: get_value_space      => boundary_function_get_value_space
     procedure :: get_value_space_time => boundary_function_get_value_space_time
  end type boundary_function_t
  
type, extends(base_scalar_function_t) :: boundary_derivative_function_t
    private
  contains
     procedure :: get_value_space      => boundary_derivative_function_get_value_space
     procedure :: get_value_space_time => boundary_derivative_function_get_value_space_time
end type boundary_derivative_function_t
   
type, extends(base_scalar_function_t) :: solution_function_t
    private 
   contains
     procedure :: get_value_space      => solution_function_get_value_space
     procedure :: get_gradient_space   => solution_function_get_gradient_space
     procedure :: get_value_space_time => solution_function_get_value_space_time
     procedure :: get_gradient_space_time   => solution_function_get_gradient_space_time
  end type solution_function_t

  type poisson_analytical_functions_t
     private
     type(source_term_t)       :: source_term
     type(boundary_function_t) :: boundary_function
     type(boundary_derivative_function_t) :: boundary_derivative_function
     type(solution_function_t) :: solution_function
   contains
     procedure :: set_num_dims            => poisson_analytical_functions_set_num_dims
     procedure :: get_source_term         => poisson_analytical_functions_get_source_term
     procedure :: get_boundary_function   => poisson_analytical_functions_get_boundary_function
     procedure :: get_boundary_derivative_function => poisson_analytical_functions_get_boundary_derivative_function
     procedure :: get_solution_function   => poisson_analytical_functions_get_solution_function
  end type poisson_analytical_functions_t

  public :: poisson_analytical_functions_t, base_scalar_function_t

contains  

  !===============================================================================================
  subroutine base_scalar_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_scalar_function_set_num_dims

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    mcheck(.false., "Implementation pending")
  end subroutine source_term_get_value_space
  
  !===============================================================================================
  subroutine source_term_get_value_space_time( this, point, time, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(in)    :: time
    real(rp)            , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      !result = 0 ! f = 0
      !result = (point%get(1)+point%get(2)) !f = x + y 
      result = (point%get(1)+point%get(2)) * 2.0_rp * time !f = ( x + y ) * t
      ! (x**2-x)*(y**2-y) - 2*t*((x**2-x)+(y**2-y))
      !result = (point%get(1)**2-point%get(1))*(point%get(2)**2-point%get(2)) - & 
      !         2.0_rp * time * ((point%get(1)**2-point%get(1)) +  (point%get(2)**2-point%get(2)))
    else if ( this%num_dims == 3 ) then
      mcheck(.false., "Implementation pending")
    end if
  end subroutine source_term_get_value_space_time

  !===============================================================================================
  subroutine boundary_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    mcheck(.false., "Implementation pending")
  end subroutine boundary_function_get_value_space 
  
  !===============================================================================================
  subroutine boundary_derivative_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_derivative_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    mcheck(.false., "Implementation pending")
  end subroutine boundary_derivative_function_get_value_space 
  
  !===============================================================================================
  subroutine boundary_function_get_value_space_time( this, point, time, result )
    implicit none
    class(boundary_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(in)    :: time
    real(rp)                  , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      !result = (point%get(1)+point%get(2)) ! uD = x + y
      !result = ( point%get(1) + point%get(2) ) * time ! uD = ( x + y ) * t
      result = ( point%get(1) + point%get(2) ) * time ** 2 ! uD = ( x + y ) * t ^ 2
      ! (x**2-x)*(y**2-y)*t
      !result = (point%get(1)**2-point%get(1))*(point%get(2)**2-point%get(2))*time
    else if ( this%num_dims == 3 ) then
      mcheck(.false., "Implementation pending")
    end if
  end subroutine boundary_function_get_value_space_time
  
    !===============================================================================================
  subroutine boundary_derivative_function_get_value_space_time( this, point, time, result )
    implicit none
    class(boundary_derivative_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(in)    :: time
    real(rp)                  , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      !result = 0 ! du/dt = 0
      !result = point%get(1) + point%get(2) !du/dt = x + y 
      result = ( point%get(1) + point%get(2) ) * 2.0_rp * time !du/dt = ( x + y ) * 2 t
    else if ( this%num_dims == 3 ) then
      mcheck(.false., "Implementation pending")
    end if
  end subroutine boundary_derivative_function_get_value_space_time

  !===============================================================================================
  subroutine solution_function_get_value_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    mcheck(.false., "Implementation pending")  
  end subroutine solution_function_get_value_space

  
  !===============================================================================================
  subroutine solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    mcheck(.false., "Implementation pending")
  end subroutine solution_function_get_gradient_space
  
  !===============================================================================================
  subroutine solution_function_get_value_space_time( this, point, time, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(in)    :: time
    real(rp)                  , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      !result = (point%get(1)+point%get(2)) ! u = x + y
      !result = ( point%get(1) + point%get(2) ) * time ! u = ( x + y ) * t
      result = ( point%get(1) + point%get(2) ) * time ** 2 ! u = ( x + y ) * t ^ 2
      !result = (point%get(1)**2-point%get(1))*(point%get(2)**2-point%get(2))*time   ! (x**2-x)*(y**2-y)*t
    else if ( this%num_dims == 3 ) then
      mcheck(.false., "Implementation pending")
    end if
  end subroutine solution_function_get_value_space_time
  
    subroutine solution_function_get_gradient_space_time ( this, point, time, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(in)    :: time
    type(vector_field_t)      , intent(inout) :: result
    mcheck(.false., "Implementation pending")
  end subroutine solution_function_get_gradient_space_time
  
  !===============================================================================================
  subroutine poisson_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(poisson_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%boundary_function%set_num_dims(num_dims)
    call this%boundary_derivative_function%set_num_dims(num_dims)
    call this%solution_function%set_num_dims(num_dims)
  end subroutine poisson_analytical_functions_set_num_dims 
  
  !===============================================================================================
  function poisson_analytical_functions_get_source_term ( this )
    implicit none
    class(poisson_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_analytical_functions_get_source_term
    poisson_analytical_functions_get_source_term => this%source_term
  end function poisson_analytical_functions_get_source_term
  
  !===============================================================================================
  function poisson_analytical_functions_get_boundary_function ( this )
    implicit none
    class(poisson_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_analytical_functions_get_boundary_function
    poisson_analytical_functions_get_boundary_function => this%boundary_function
  end function poisson_analytical_functions_get_boundary_function
 
  !===============================================================================================
  function poisson_analytical_functions_get_boundary_derivative_function ( this )
    implicit none
    class(poisson_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_analytical_functions_get_boundary_derivative_function
    poisson_analytical_functions_get_boundary_derivative_function => this%boundary_derivative_function
  end function poisson_analytical_functions_get_boundary_derivative_function
  !===============================================================================================
  function poisson_analytical_functions_get_solution_function ( this )
    implicit none
    class(poisson_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: poisson_analytical_functions_get_solution_function
    poisson_analytical_functions_get_solution_function => this%solution_function
  end function poisson_analytical_functions_get_solution_function

end module poisson_analytical_functions_names
!***************************************************************************************************
