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

module fempar_sm_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: base_scalar_function_t
    private
    integer(ip) :: num_dimensions = -1  
  contains
    procedure :: set_num_dimensions    => base_scalar_function_set_num_dimensions
  end type base_scalar_function_t

  type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dimensions = -1  
  contains
    procedure :: set_num_dimensions    => base_vector_function_set_num_dimensions
  end type base_vector_function_t

  !===============================================================================================

  type, extends(base_scalar_function_t) :: boundary_function_t
    private
   contains
     procedure :: get_value_space => boundary_function_get_value_space
  end type boundary_function_t
  
  !===============================================================================================
  type, extends(base_scalar_function_t) :: scalar_source_term_t
    private 
   contains
     procedure :: get_value_space  => scalar_source_term_get_value_space
  end type scalar_source_term_t
  
  type, extends(base_vector_function_t) :: vector_source_term_t
   contains
     procedure :: get_value_space => vector_source_term_get_value_space
  end type vector_source_term_t

  !===============================================================================================
  type, extends(base_scalar_function_t) :: scalar_solution_function_t
    private 
   contains
     procedure :: get_value_space    => scalar_solution_function_get_value_space
     procedure :: get_gradient_space => scalar_solution_function_get_gradient_space
  end type scalar_solution_function_t

  type, extends(base_vector_function_t) :: vector_solution_function_t
   contains
     procedure :: get_value_space    => vector_solution_function_get_value_space
     procedure :: get_gradient_space => vector_solution_function_get_gradient_space
  end type vector_solution_function_t
  !===============================================================================================


  type fempar_sm_analytical_functions_t
     private
     type(scalar_source_term_t)       :: scalar_source_term
     type(vector_source_term_t)       :: vector_source_term
     type(boundary_function_t)        :: boundary_function
     type(scalar_solution_function_t) :: scalar_solution_function
     type(vector_solution_function_t) :: vector_solution_function
   contains
     procedure :: set_num_dimensions      => fempar_sm_analytical_functions_set_num_dimensions
     procedure :: get_boundary_function   => fempar_sm_analytical_functions_get_boundary_function
     procedure :: get_scalar_source_term  => fempar_sm_analytical_functions_get_scalar_source_term
     procedure :: get_vector_source_term  => fempar_sm_analytical_functions_get_vector_source_term
     procedure :: get_scalar_solution_function   => fempar_sm_analytical_functions_get_scalar_solution_function
     procedure :: get_vector_solution_function   => fempar_sm_analytical_functions_get_vector_solution_function
  end type fempar_sm_analytical_functions_t

  public :: fempar_sm_analytical_functions_t

contains  

  subroutine base_scalar_function_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine base_scalar_function_set_num_dimensions

  subroutine base_vector_function_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine base_vector_function_set_num_dimensions

  !===============================================================================================
  subroutine boundary_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    if ( this%num_dimensions == 2 ) then
      result = point%get(1) + point%get(2) ! x+y
    else if ( this%num_dimensions == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3) ! x+y+z
    end if  
  end subroutine boundary_function_get_value_space 

  !===============================================================================================
  subroutine scalar_source_term_get_value_space ( this, point, result )
    implicit none
    class(scalar_source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    result = 0.0_rp 
  end subroutine scalar_source_term_get_value_space

  subroutine vector_source_term_get_value_space ( this, point, result )
    implicit none
    class(vector_source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set(1,0.0_rp)
      call result%set(2,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
    else
      call result%set(1,0.0_rp)
      call result%set(2,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
      call result%set(3,0.0_rp) !2 * ( pi**2 ) * sin ( pi * point%get(1) ) * sin ( pi * point%get(2) )
    end if  
  end subroutine vector_source_term_get_value_space

  !===============================================================================================
  subroutine scalar_solution_function_get_value_space ( this, point, result )
    implicit none
    class(scalar_solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    if ( this%num_dimensions == 2 ) then
      result = point%get(1) + point%get(2) ! x+y 
    else if ( this%num_dimensions == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3) ! x+y+z
    end if  
      
  end subroutine scalar_solution_function_get_value_space

    subroutine vector_solution_function_get_value_space ( this, point, result )
    implicit none
    class(vector_solution_function_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set(1, point%get(1)+point%get(2) ) 
      call result%set(2, point%get(1)+point%get(2) ) 
    else
      call result%set(1, point%get(1)+point%get(2)+point%get(3) ) 
      call result%set(2, point%get(1)+point%get(2)+point%get(3) ) 
      call result%set(3, point%get(1)+point%get(2)+point%get(3) )
    end if
  end subroutine vector_solution_function_get_value_space

  !===============================================================================================
  subroutine scalar_solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(scalar_solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
    if ( this%num_dimensions == 2 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, 1.0_rp )
    else if ( this%num_dimensions == 3 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, 1.0_rp )
      call result%set( 3, 1.0_rp ) 
    end if
  end subroutine scalar_solution_function_get_gradient_space

  subroutine vector_solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(vector_solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    if ( this%num_dimensions == 2 ) then
      call result%set( 1, 1, 1.0_rp ) 
      call result%set( 2, 1, 1.0_rp )
      call result%set( 1, 2, 1.0_rp ) 
      call result%set( 2, 2, 1.0_rp )
    else
      call result%init(1.0_rp)
    end if
    !if ( this%num_dimensions == 2 ) then
    !  call result%set( 1, 1, 1.0_rp ) 
    !  call result%set( 2, 1, 0.0_rp )
    !  call result%set( 1, 2, 1.0_rp ) 
    !  call result%set( 2, 2, 0.0_rp )
    !else
    !  call result%set( 1, 1, 1.0_rp ) 
    !  call result%set( 2, 1, 0.0_rp )
    !  call result%set( 3, 1, 0.0_rp )
    !  call result%set( 1, 2, 1.0_rp ) 
    !  call result%set( 2, 2, 0.0_rp )
    !  call result%set( 3, 2, 0.0_rp )
    !  call result%set( 1, 3, 1.0_rp ) 
    !  call result%set( 2, 3, 0.0_rp )
    !  call result%set( 3, 3, 0.0_rp )
    !end if
  end subroutine vector_solution_function_get_gradient_space
  
  !===============================================================================================
  subroutine fempar_sm_analytical_functions_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(fempar_sm_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    call this%scalar_source_term%set_num_dimensions(num_dimensions)
    call this%vector_source_term%set_num_dimensions(num_dimensions)
    call this%boundary_function%set_num_dimensions(num_dimensions)
    call this%scalar_solution_function%set_num_dimensions(num_dimensions)
    call this%vector_solution_function%set_num_dimensions(num_dimensions)
  end subroutine fempar_sm_analytical_functions_set_num_dimensions 
  
  !===============================================================================================
  function fempar_sm_analytical_functions_get_scalar_source_term ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_scalar_source_term
    fempar_sm_analytical_functions_get_scalar_source_term => this%scalar_source_term
  end function fempar_sm_analytical_functions_get_scalar_source_term

  function fempar_sm_analytical_functions_get_vector_source_term ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: fempar_sm_analytical_functions_get_vector_source_term
    fempar_sm_analytical_functions_get_vector_source_term => this%vector_source_term
  end function fempar_sm_analytical_functions_get_vector_source_term
  
  !===============================================================================================
  function fempar_sm_analytical_functions_get_boundary_function ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_boundary_function
    fempar_sm_analytical_functions_get_boundary_function => this%boundary_function
  end function fempar_sm_analytical_functions_get_boundary_function
 
  !===============================================================================================
  function fempar_sm_analytical_functions_get_scalar_solution_function ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: fempar_sm_analytical_functions_get_scalar_solution_function
    fempar_sm_analytical_functions_get_scalar_solution_function => this%scalar_solution_function
  end function fempar_sm_analytical_functions_get_scalar_solution_function

  function fempar_sm_analytical_functions_get_vector_solution_function ( this )
    implicit none
    class(fempar_sm_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: fempar_sm_analytical_functions_get_vector_solution_function
    fempar_sm_analytical_functions_get_vector_solution_function => this%vector_solution_function
  end function fempar_sm_analytical_functions_get_vector_solution_function

end module fempar_sm_analytical_functions_names
!***************************************************************************************************
