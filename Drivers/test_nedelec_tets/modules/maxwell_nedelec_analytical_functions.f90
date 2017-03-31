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

module maxwell_nedelec_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
 
 type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dimensions = -1  
  contains
    procedure :: set_num_dimensions    => base_vector_function_set_num_dimensions
  end type base_vector_function_t
  
  type, extends(base_vector_function_t) :: source_term_t
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t
  
  type, extends(base_vector_function_t) :: solution_t
   contains
     procedure :: get_value_space    => solution_get_value_space
     procedure :: get_gradient_space => solution_get_gradient_space
  end type solution_t
  
    type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dimensions = -1  
  contains
  end type base_scalar_function_t
  
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
  
  
  type maxwell_nedelec_analytical_functions_t
     private
     type(source_term_t)                     :: source_term
     type(solution_t)                        :: solution
	 type(boundary_function_Hx_t)            :: boundary_function_Hx
	 type(boundary_function_Hy_t)            :: boundary_function_Hy
	 type(boundary_function_Hz_t)            :: boundary_function_Hz
   contains
     procedure :: set_num_dimensions               => mn_set_num_dimensions
     procedure :: get_source_term                  => mn_get_source_term
     procedure :: get_solution                     => mn_get_solution
	 procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
	 procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
	 procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
  end type maxwell_nedelec_analytical_functions_t

  public :: maxwell_nedelec_analytical_functions_t

contains  

  !===============================================================================================
  subroutine base_vector_function_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine base_vector_function_set_num_dimensions

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
 
	real(rp) :: x,y,z 
	
	assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
	x = point%get(1); y=point%get(2); z=point%get(3)  
      
	 call result%init(0.0_rp) 
	 !call result%set(1, -y) 
	 !call result%set(2, 0.0_rp) 
	 !if (this%num_dimensions == 3) then 
	 !call result%set(3, 0.0_rp) 
	 !end if 
     call result%set(1, pi*pi*(sin(pi*x)*sin(pi*y))  + sin(pi*x)*sin(pi*y) )
     call result%set(2, pi*pi*(cos(pi*x)*cos(pi*y)))
     if (this%num_dimensions == 3) then 
	 call result%set(3, 0.0_rp) 
	 end if 
     !call result%set(1,  sin(pi*x)*sin(pi*y) )
     !call result%set(2, 0.0_rp )
	
	! call result%set(1, 2.0_rp*x*(1.0_rp-x) + x*(1-x)*y*(1-y) ) 
	! call result%set(2, (1.0_rp-2.0_rp*x)*(1.0_rp-2.0_rp*y) ) 
	! call result%set(3, 0.0_rp) 
	
	!call result%set(1, -y ) 
	!call result%set(2,  x ) 
	!if (this%num_dimensions == 3) then 
	!call result%set(3,  0.0_rp) 
	!end if 
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine source_term_get_gradient_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)        :: this
    type(point_t)       , intent(in)        :: point
    type(tensor_field_t), intent(inout)     :: result
    check(.false.)
  end subroutine source_term_get_gradient_space

  !===============================================================================================
  subroutine solution_get_value_space ( this, point, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
	
    real(rp) :: x,y,z 
	assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 )
	x = point%get(1); y=point%get(2); z=point%get(3) 
	
	 call result%init(0.0_rp) 
	 !call result%set(1, -y) 
	 !call result%set(2,  0.0_rp) 
	 !if (this%num_dimensions == 3) then 
	 !call result%set(3, 0.0_rp) 
	 !end if 
     call result%set(1, sin(pi*x)*sin(pi*y) )
	 call result%set(2, 0.0_rp) 
     if (this%num_dimensions == 3) then 
	 call result%set(3, 0.0_rp) 
	 end if 
	 ! call result%set(1, x*(1-x)*y*(1-y) ) 
	 !call result%set(2, 0.0_rp) 
	 !call result%set(3, 0.0_rp) 
	!call result%set(1, -y ) 
	!call result%set(2,  x ) 
	!if (this%num_dimensions == 3) then 
	!call result%set(3,  0.0_rp) 
	!end if 
  end subroutine solution_get_value_space

  !===============================================================================================
  subroutine solution_get_gradient_space ( this, point, result )
    implicit none
    class(solution_t)   , intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
	
	real(rp) :: x,y,z 
	x = point%get(1); y=point%get(2); z=point%get(3)

	call result%init(0.0_rp) 
 	call result%set(1, 1, pi*cos(pi*x)*sin(pi*y) )
  	call result%set(2, 1, pi*sin(pi*x)*cos(pi*y) ) 
	!call result%set(1, 1, y*(1.0_rp-y)*(1.0_rp-2.0_rp*x) )
	!call result%set(2, 1, x*(1.0_rp-x)*(1.0_rp-2.0_rp*y) )
	
	!call result%set(2,1, -1.0_rp)
	!call result%set(1,2,  1.0_rp)
  end subroutine solution_get_gradient_space
  
  !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
		real(rp) :: x,y,z 
	x = point%get(1); y=point%get(2); z=point%get(3)
	 result = sin(pi*x)*sin(pi*y) 
  end subroutine boundary_function_Hx_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
		real(rp) :: x,y,z 
	x = point%get(1); y=point%get(2); z=point%get(3)
    result = 0.0_rp        
  end subroutine boundary_function_Hy_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hz_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
		real(rp) :: x,y,z 
	x = point%get(1); y=point%get(2); z=point%get(3)
    result = 0.0_rp         
  end subroutine boundary_function_Hz_get_value_space

  !===============================================================================================
  subroutine mn_set_num_dimensions ( this, num_dimensions )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dimensions
    call this%source_term%set_num_dimensions(num_dimensions)
    call this%solution%set_num_dimensions(num_dimensions)
  end subroutine mn_set_num_dimensions

  !===============================================================================================
  function mn_get_solution ( this )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_solution
    mn_get_solution => this%solution
  end function mn_get_solution

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term

  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz 
  end function mn_get_boundary_function_Hz

end module maxwell_nedelec_analytical_functions_names
!***************************************************************************************************


