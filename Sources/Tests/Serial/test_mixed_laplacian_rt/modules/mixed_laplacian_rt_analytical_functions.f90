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

module mixed_laplacian_rt_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

 type, extends(scalar_function_t) :: base_scalar_function_t
    private
    integer(ip) :: num_dims = -1  
  contains
    procedure :: set_num_dims    => base_scalar_function_set_num_dims
  end type base_scalar_function_t
  
 type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dims = -1  
  contains
    procedure :: set_num_dims    => base_vector_function_set_num_dims
  end type base_vector_function_t
  
  type, extends(base_scalar_function_t) :: pressure_source_term_t
   contains
     procedure :: get_value_space => pressure_source_term_get_value_space
  end type pressure_source_term_t
 
  type, extends(base_scalar_function_t) :: pressure_boundary_function_t
   contains
     procedure :: get_value_space => pressure_boundary_function_get_value_space
  end type pressure_boundary_function_t
  
  type, extends(base_scalar_function_t) :: pressure_solution_t
   contains
     procedure :: get_value_space    => pressure_solution_get_value_space
     procedure :: get_gradient_space => pressure_solution_get_gradient_space
  end type pressure_solution_t
  
  type, extends(base_vector_function_t) :: velocity_solution_t
   contains
     procedure :: get_value_space    => velocity_solution_get_value_space
     procedure :: get_gradient_space => velocity_solution_get_gradient_space
  end type velocity_solution_t

  type mixed_laplacian_rt_analytical_functions_t
     private
     type(pressure_source_term_t)         :: pressure_source_term
     type(pressure_boundary_function_t)   :: pressure_boundary_function
     type(pressure_solution_t)            :: pressure_solution
     type(velocity_solution_t)            :: velocity_solution
   contains
     procedure :: set_num_dims               => mlrtaf_set_num_dims
     procedure :: get_pressure_source_term         => mlrtaf_get_pressure_source_term
     procedure :: get_pressure_boundary_function   => mlrtaf_get_pressure_boundary_function
     procedure :: get_pressure_solution            => mlrtaf_get_pressure_solution
     procedure :: get_velocity_solution            => mlrtaf_get_velocity_solution
  end type mixed_laplacian_rt_analytical_functions_t

  public :: mixed_laplacian_rt_analytical_functions_t

contains  

  !===============================================================================================
  subroutine base_scalar_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_scalar_function_set_num_dims
 
  
  !===============================================================================================
  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims

  !===============================================================================================
  subroutine pressure_source_term_get_value_space ( this, point, result )
    implicit none
    class(pressure_source_term_t), intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    real(rp)                     , intent(inout) :: result
    result = 0.0_rp
  end subroutine pressure_source_term_get_value_space

  !===============================================================================================
  subroutine pressure_boundary_function_get_value_space ( this, point, result )
    implicit none
    class(pressure_boundary_function_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      result = point%get(1) + point%get(2) ! x+y
    else if ( this%num_dims == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3) ! x+y+z
    end if  
  end subroutine pressure_boundary_function_get_value_space
  
  !===============================================================================================
  subroutine pressure_solution_get_value_space ( this, point, result )
    implicit none
    class(pressure_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      result = point%get(1) + point%get(2) ! x+y
    else if ( this%num_dims == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3) ! x+y+z
    end if  
  end subroutine pressure_solution_get_value_space
  
  !===============================================================================================
  subroutine pressure_solution_get_gradient_space ( this, point, result )
    implicit none
    class(pressure_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      call result%set(1,1.0_rp)
      call result%set(2,1.0_rp)
    else if ( this%num_dims == 3 ) then
      call result%set(1,1.0_rp)
      call result%set(2,1.0_rp)
      call result%set(3,1.0_rp)
    end if  
    
  end subroutine pressure_solution_get_gradient_space
  
 !===============================================================================================
  subroutine velocity_solution_get_value_space ( this, point, result )
    implicit none
    class(velocity_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    if ( this%num_dims == 2 ) then
      call result%set(1,-1.0_rp)
      call result%set(2,-1.0_rp)
    else if ( this%num_dims == 3 ) then
      call result%set(1,-1.0_rp)
      call result%set(2,-1.0_rp)
      call result%set(3,-1.0_rp)
    end if  
  end subroutine velocity_solution_get_value_space
  
  !===============================================================================================
  subroutine velocity_solution_get_gradient_space ( this, point, result )
    implicit none
    class(velocity_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
    call result%init(0.0_rp)
  end subroutine velocity_solution_get_gradient_space 
  
  !===============================================================================================
  subroutine mlrtaf_set_num_dims ( this, num_dims )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%pressure_source_term%set_num_dims(num_dims)
    call this%pressure_boundary_function%set_num_dims(num_dims)
    call this%pressure_solution%set_num_dims(num_dims)
    call this%velocity_solution%set_num_dims(num_dims)
  end subroutine mlrtaf_set_num_dims 
  
  !===============================================================================================
  function mlrtaf_get_pressure_source_term ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mlrtaf_get_pressure_source_term
    mlrtaf_get_pressure_source_term => this%pressure_source_term
  end function mlrtaf_get_pressure_source_term
  
  !===============================================================================================
  function mlrtaf_get_pressure_boundary_function ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mlrtaf_get_pressure_boundary_function
    mlrtaf_get_pressure_boundary_function => this%pressure_boundary_function
  end function mlrtaf_get_pressure_boundary_function
  
  !===============================================================================================
  function mlrtaf_get_pressure_solution ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mlrtaf_get_pressure_solution
    mlrtaf_get_pressure_solution => this%pressure_solution
  end function mlrtaf_get_pressure_solution
  
    !===============================================================================================
  function mlrtaf_get_velocity_solution ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mlrtaf_get_velocity_solution
    mlrtaf_get_velocity_solution => this%velocity_solution
  end function mlrtaf_get_velocity_solution
  
end module mixed_laplacian_rt_analytical_functions_names
!***************************************************************************************************


