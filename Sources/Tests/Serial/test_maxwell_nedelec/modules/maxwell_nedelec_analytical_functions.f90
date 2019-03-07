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

  character(len=*), parameter :: in_fe_space       ='in_fe_space'
  character(len=*), parameter :: sinusoidal_2D     ='sinusoidal_2D'
  character(len=*), parameter :: sinusoidal_3D     ='sinusoidal_3D' 

  type, extends(vector_function_t) :: base_vector_function_t
     character(len=:), allocatable :: function_case
   contains
     procedure :: set_case        => base_vector_function_set_case
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
     character(len=:), allocatable :: function_case
   contains
     procedure :: set_case        => base_scalar_function_set_case
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
     procedure :: set_case                         => mn_set_case 
     procedure :: set_num_dims                     => mn_set_num_dims
     procedure :: get_source_term                  => mn_get_source_term
     procedure :: get_solution                     => mn_get_solution
     procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
     procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
     procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
  end type maxwell_nedelec_analytical_functions_t

  public :: in_fe_space, sinusoidal_2D, sinusoidal_3D 
  public :: maxwell_nedelec_analytical_functions_t

contains  

  !===============================================================================================
  subroutine mn_set_case ( this, function_case )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), intent(inout) :: this
    character(len=*)                     , intent(in)    :: function_case 
    call this%source_term%set_case(function_case)
    call this%boundary_function_Hx%set_case(function_case)
    call this%boundary_function_Hy%set_case(function_case)
    call this%boundary_function_Hz%set_case(function_case)
    call this%solution%set_case(function_case)
  end subroutine mn_set_case

  subroutine base_scalar_function_set_case ( this, function_case )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    character(len=*)             , intent(in)       :: function_case
    this%function_case = function_case
  end subroutine base_scalar_function_set_case

  subroutine base_vector_function_set_case ( this, function_case )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    character(len=*)             , intent(in)       :: function_case
    this%function_case = function_case
  end subroutine base_vector_function_set_case

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    real(rp) :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x = point%get(1); y = point%get(2); z = point%get(3)
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       call result%set(1, -point%get(2) )
       call result%set(2,  point%get(1))
       if ( this%get_num_dims() == 3 ) then
          call result%set(3, 0.0_rp) 
       end if
    case ( sinusoidal_2D ) 
       call result%set(1, (2.0_rp*pi**2+1.0_rp)*cos(pi*x)*cos(pi*y) ) 
       call result%set(2, (2.0_rp*pi**2+1.0_rp)*sin(pi*x)*sin(pi*y) ) 
       call result%set(3, 0.0_rp) 
    case ( sinusoidal_3D ) 
       call result%set(1, (pi**2+1.0_rp)*cos(pi*x)*cos(pi*y) + pi**2*sin(pi*x)*sin(pi*z) )
       call result%set(2, (pi**2+1.0_rp)*sin(pi*y)*sin(pi*z) + pi**2*sin(pi*x)*sin(pi*y) ) 
       call result%set(3, (pi**2+1.0_rp)*cos(pi*x)*cos(pi*z) + pi**2*cos(pi*y)*cos(pi*z) ) 
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select

  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine source_term_get_gradient_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
    check(.false.)
  end subroutine source_term_get_gradient_space

  !===============================================================================================
  subroutine solution_get_value_space ( this, point, result )
    implicit none
    class(solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    real(rp) :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x = point%get(1); y = point%get(2); z = point%get(3)
    call result%init(0.0_rp) 
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       call result%set(1, -point%get(2))
       call result%set(2,  point%get(1))  
    case ( sinusoidal_2D ) 
       call result%set(1, cos(pi*x)*cos(pi*y) ) 
       call result%set(2, sin(pi*x)*sin(pi*y) ) 
    case ( sinusoidal_3D ) 
       call result%set(1, cos(pi*x)*cos(pi*y) ) 
       call result%set(2, sin(pi*y)*sin(pi*z) ) 
       call result%set(3, cos(pi*x)*cos(pi*z) )
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select
  end subroutine solution_get_value_space

  !===============================================================================================
  subroutine solution_get_gradient_space ( this, point, result )
    implicit none
    class(solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t), intent(inout) :: result
    real(rp) :: x,y,z
    x = point%get(1); y = point%get(2); z = point%get(3)
    call result%init(0.0_rp) 
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       call result%set(2,1, -1.0_rp)
       call result%set(1,2,  1.0_rp)
    case ( sinusoidal_2D ) 
       call result%set(1, 1, -pi*sin(pi*x)*cos(pi*y) ) 
       call result%set(2, 1, -pi*cos(pi*x)*sin(pi*y) )
       call result%set(1, 2, pi*cos(pi*x)*sin(pi*y)  ) 
       call result%set(2, 2, pi*sin(pi*x)*cos(pi*y)  ) 
    case ( sinusoidal_3D ) 
       call result%set(1, 1, -pi*sin(pi*x)*cos(pi*y) ) 
       call result%set(2, 1, -pi*cos(pi*x)*sin(pi*y) )

       call result%set(2, 2, pi*cos(pi*y)*sin(pi*z) ) 
       call result%set(3, 2, pi*sin(pi*y)*cos(pi*z) ) 

       call result%set(1, 3, -pi*sin(pi*x)*cos(pi*z) ) 
       call result%set(3, 3, -pi*cos(pi*x)*sin(pi*z) ) 
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select

  end subroutine solution_get_gradient_space

  !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z
    x = point%get(1); y = point%get(2); z = point%get(3)
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       result = -point%get(2)
    case ( sinusoidal_2D, sinusoidal_3D )
       result = cos(pi*x)*cos(pi*y)
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select
  end subroutine boundary_function_Hx_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z
    x = point%get(1); y = point%get(2); z = point%get(3)
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       result = point%get(1)  
    case ( sinusoidal_2D)
       result = sin(pi*x)*sin(pi*y) 
    case ( sinusoidal_3D) 
       result = sin(pi*y)*sin(pi*z)
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select

  end subroutine boundary_function_Hy_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hz_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    real(rp) :: x,y,z
    x = point%get(1); y = point%get(2); z = point%get(3)
    select case ( this%function_case ) 
    case ( in_fe_space, sinusoidal_2D ) 
       result = 0.0_rp   
    case ( sinusoidal_3D )      
       result = cos(pi*x)*cos(pi*z)
    case DEFAULT 
       massert( .false., 'Unknown analytical function case' ) 
    end select
  end subroutine boundary_function_Hz_get_value_space

  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(maxwell_nedelec_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%solution%set_num_dims(num_dims)
  end subroutine mn_set_num_dims

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


