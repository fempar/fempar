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

module projections_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
 
 type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dims = -1  
  contains
    procedure :: set_num_dims    => base_vector_function_set_num_dims
  end type base_vector_function_t
		
		type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims = -1  
  contains
		procedure :: set_num_dims => base_scalar_function_set_num_dims 
  end type base_scalar_function_t
  
  type, extends(base_vector_function_t) :: source_term_t
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t
  
  type, extends(base_vector_function_t) :: magnetic_field_solution_t
   contains
     procedure :: get_value_space    => magnetic_field_solution_get_value_space
     procedure :: get_gradient_space => magnetic_field_solution_get_gradient_space
  end type magnetic_field_solution_t
		
		  type, extends(base_scalar_function_t) :: pressure_solution_t
   contains
     procedure :: get_value_space    => pressure_solution_get_value_space
     procedure :: get_gradient_space => pressure_solution_get_gradient_space
  end type pressure_solution_t
  
  type, extends(base_scalar_function_t) :: boundary_function_Hx_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hx_get_value_space
  end type boundary_function_Hx_t 
  
    type, extends(base_scalar_function_t) :: boundary_function_Hy_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hy_get_value_space
  end type boundary_function_Hy_t 
  
   type, extends(base_scalar_function_t) :: boundary_function_Hz_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_Hz_get_value_space
  end type boundary_function_Hz_t 
		
		   type, extends(base_scalar_function_t) :: boundary_function_pressure_t
    private 
   contains
     procedure :: get_value_space    => boundary_function_pressure_get_value_space 
  end type boundary_function_pressure_t 
  
  
  type projections_analytical_functions_t
     private
  type(source_term_t)                     :: source_term
  type(magnetic_field_solution_t)         :: magnetic_field_solution
		type(pressure_solution_t)               :: pressure_solution
	 type(boundary_function_Hx_t)            :: boundary_function_Hx
	 type(boundary_function_Hy_t)            :: boundary_function_Hy
	 type(boundary_function_Hz_t)            :: boundary_function_Hz
		type(boundary_function_pressure_t)      :: boundary_function_pressure 
   contains
  procedure :: set_num_dims                     => mn_set_num_dims
  procedure :: get_source_term                  => mn_get_source_term
  procedure :: get_magnetic_field_solution      => mn_get_magnetic_field_solution
		procedure :: get_pressure_solution            => mn_get_pressure_solution
	 procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
	 procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
	 procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
		procedure :: get_boundary_function_pressure   => mn_get_boundary_function_pressure
  end type projections_analytical_functions_t

  public :: projections_analytical_functions_t

contains  

  !===============================================================================================
  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims
		
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
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )

    call result%set(1, -point%get(2) )
    call result%set(2,  point%get(1))
    if ( this%num_dims == 3 ) then
       call result%set(3, 0.0_rp) 
    end if
											
  end subroutine source_term_get_value_space
		
		  !===============================================================================================
  subroutine source_term_get_value_space_time ( this, point, time, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
				real(rp)                , intent(in)    :: time
    type(vector_field_t)    , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )

    call result%set(1, -point%get(2) * time )
    call result%set(2,  point%get(1) * time )
    if ( this%num_dims == 3 ) then
       call result%set(3, 0.0_rp) 
    end if
											
  end subroutine source_term_get_value_space_time
		
  !===============================================================================================
  subroutine magnetic_field_solution_get_value_space ( this, point, result )
    implicit none
    class(magnetic_field_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )

    call result%set(1, -point%get(2))
    call result%set(2,  point%get(1))  
    if ( this%num_dims == 3 ) then
       call result%set(3, 0.0_rp)
    end if

  end subroutine magnetic_field_solution_get_value_space
		
		  !===============================================================================================
  subroutine magnetic_field_solution_get_value_space_time ( this, point, time, result )
    implicit none
    class(magnetic_field_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
				real(rp)                , intent(in)    :: time 
    type(vector_field_t)    , intent(inout) :: result
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )

    call result%set(1, -point%get(2) * time)
    call result%set(2,  point%get(1) * time)  
    if ( this%num_dims == 3 ) then
       call result%set(3, 0.0_rp)
    end if

  end subroutine magnetic_field_solution_get_value_space_time

  !===============================================================================================
  subroutine magnetic_field_solution_get_gradient_space ( this, point, result )
    implicit none
    class(magnetic_field_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
				call result%init(0.0_rp) 
    call result%set(2,1, -1.0_rp)
	   call result%set(1,2,  1.0_rp)
  end subroutine magnetic_field_solution_get_gradient_space
		
		  !===============================================================================================
  subroutine magnetic_field_solution_get_gradient_space_time ( this, point, time, result )
    implicit none
    class(magnetic_field_solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
				real(rp)                , intent(in)    :: time 
    type(tensor_field_t), intent(inout) :: result
				call result%init(0.0_rp) 
    call result%set(2,1, -1.0_rp * time)
	   call result%set(1,2,  1.0_rp * time)
  end subroutine magnetic_field_solution_get_gradient_space_time
		
		  !===============================================================================================
  subroutine pressure_solution_get_value_space ( this, point, result )
    implicit none
    class(pressure_solution_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    result = point%get(1)+point%get(2)
  end subroutine pressure_solution_get_value_space

  !===============================================================================================
  subroutine pressure_solution_get_gradient_space ( this, point, result )
    implicit none
    class(pressure_solution_t)   , intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    type(vector_field_t)         , intent(inout) :: result
				call result%set(1, 1.0_rp) 
				call result%set(2, 1.0_rp) 
  end subroutine pressure_solution_get_gradient_space
  
  !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    result = -point%get(2)
  end subroutine boundary_function_Hx_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
     result = point%get(1)
  end subroutine boundary_function_Hy_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hz_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    result = 0.0_rp          
  end subroutine boundary_function_Hz_get_value_space
		
		  !===============================================================================================
  subroutine boundary_function_pressure_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_pressure_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    result = point%get(1) + point%get(2) 
  end subroutine boundary_function_pressure_get_value_space

  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(projections_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%magnetic_field_solution%set_num_dims(num_dims)
				call this%pressure_solution%set_num_dims(num_dims)
  end subroutine mn_set_num_dims

  !===============================================================================================
  function mn_get_magnetic_field_solution ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_magnetic_field_solution
    mn_get_magnetic_field_solution => this%magnetic_field_solution
  end function mn_get_magnetic_field_solution
		
		  !===============================================================================================
  function mn_get_pressure_solution ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_pressure_solution
    mn_get_pressure_solution => this%pressure_solution
  end function mn_get_pressure_solution

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term

  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz 
  end function mn_get_boundary_function_Hz
		
		  !===============================================================================================
  function mn_get_boundary_function_pressure ( this )
    implicit none
    class(projections_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_pressure
    mn_get_boundary_function_pressure => this%boundary_function_pressure 
  end function mn_get_boundary_function_pressure

end module projections_analytical_functions_names
!***************************************************************************************************


