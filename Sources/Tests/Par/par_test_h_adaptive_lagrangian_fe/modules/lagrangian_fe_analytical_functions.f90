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

module lagrangian_fe_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: base_scalar_function_t
    private
    integer(ip) :: order 
  contains
    procedure :: set_polynomial_order  => base_scalar_function_set_polynomial_order
  end type base_scalar_function_t
  
  type, extends(base_scalar_function_t) :: source_term_t
    private 
   contains
     procedure :: get_value_space      => source_term_get_value_space
     procedure :: get_values_set_space => source_term_get_values_set_space
  end type source_term_t

  type, extends(base_scalar_function_t) :: boundary_function_t
    private
   contains
     procedure :: get_value_space    => boundary_function_get_value_space
     procedure :: get_gradient_space => boundary_function_get_gradient_space
  end type boundary_function_t
  
  type, extends(tensor_function_t) :: solution_stress_function_t
   contains
     procedure :: get_value_space      => solution_stress_function_get_value_space
     !procedure :: get_gradient_space   => solution_stress_function_get_gradient_space
  end type solution_stress_function_t  
  
  type, extends(vector_function_t) :: solution_displacement_function_t
   contains
     procedure :: get_value_space      => solution_displacement_function_get_value_space
     procedure :: get_gradient_space   => solution_displacement_function_get_gradient_space
  end type solution_displacement_function_t  

  type, extends(base_scalar_function_t) :: solution_temperature_function_t
    private 
   contains
     procedure :: get_value_space    => solution_temperature_function_get_value_space
     procedure :: get_gradient_space => solution_temperature_function_get_gradient_space
  end type solution_temperature_function_t
  
  type, extends(base_scalar_function_t) :: solution_pressure_function_t
    private 
   contains
     procedure :: get_value_space    => solution_pressure_function_get_value_space
     procedure :: get_gradient_space => solution_pressure_function_get_gradient_space
  end type solution_pressure_function_t  

  type lagrangian_fe_analytical_functions_t
     private
     type(source_term_t)       :: source_term
     type(boundary_function_t) :: boundary_function
     type(solution_temperature_function_t)  :: solution_temperature_function
     type(solution_pressure_function_t)     :: solution_pressure_function
     type(solution_displacement_function_t) :: solution_displacement_function
     type(solution_stress_function_t)       :: solution_stress_function
   contains
     procedure :: set_num_dims                     => lagrangian_fe_analytical_functions_set_num_dims
     procedure :: set_solution_polynomial_order    => lfe_analytical_functions_set_solution_polynomial_order
     procedure :: get_source_term                  => lagrangian_fe_analytical_functions_get_source_term
     procedure :: get_boundary_function            => lagrangian_fe_analytical_functions_get_boundary_function
     procedure :: get_solution_temperature_function  => lagrangian_af_get_solution_temperature_function
     procedure :: get_solution_pressure_function     => lagrangian_af_get_solution_pressure_function
     procedure :: get_solution_displacement_function => lagrangian_af_get_solution_displacement_function
     procedure :: get_solution_stress_function       => lagrangian_af_get_solution_stress_function
  end type lagrangian_fe_analytical_functions_t

  public :: lagrangian_fe_analytical_functions_t

contains  
  
    subroutine base_scalar_function_set_polynomial_order ( this, order )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  order
    this%order = order
  end subroutine base_scalar_function_set_polynomial_order
  
  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    real(rp) :: n
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result =  point%get(1) + point%get(2)
    else if ( this%get_num_dims() == 3 ) then
      result =  point%get(1) + point%get(2) + point%get(3) 
    end if
    !n = real(this%order, rp)
    !if ( this%get_num_dims() == 2 ) then
    !  result =  -n*(n-1.0_rp)*(point%get(1)**(n-2.0_rp) + point%get(2)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2})
    !else if ( this%get_num_dims() == 3 ) then
    !  result =  -n*(n-1.0_rp)*(point%get(1)**(n-2.0_rp) + & 
    !            point%get(2)**(n-2.0_rp) + point%get(3)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2}+z^{n-2})
    !end if
  end subroutine source_term_get_value_space
  
  !===============================================================================================
  subroutine source_term_get_values_set_space( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point(:)
    real(rp)                , intent(inout) :: result(:)
    integer(ip) :: n, i, num_points
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    num_points = size(point)
    do i = 1, num_points
      if ( this%get_num_dims() == 2 ) then
       result(i) = 0.0_rp ! x+y
      else if ( this%get_num_dims() == 3 ) then
       result(i) = 0.0_rp ! x+y+z
      end if  
    end do
    !num_points = size(point)
    !n = real(this%order, rp)
    !do i = 1, num_points
    !  if ( this%get_num_dims() == 2 ) then
    !    result =  -n*(n-1.0_rp)*(point(i)%get(1)**(n-2.0_rp) + point(i)%get(2)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2})
    !  else if ( this%get_num_dims() == 3 ) then
    !    result =  -n*(n-1.0_rp)*(point(i)%get(1)**(n-2.0_rp) + & 
    !              point(i)%get(2)**(n-2.0_rp) + point(i)%get(3)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2}+z^{n-2})
    !  end if 
    !end do
  end subroutine source_term_get_values_set_space
  
  !===============================================================================================
  subroutine boundary_function_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1) + point%get(2)
      !result = point%get(1)**this%order + point%get(2)**this%order ! x^n+y^n
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3) 
      !result = point%get(1)**this%order + point%get(2)**this%order + point%get(3)**this%order ! x^n+y^n+z^n
    end if  
  end subroutine boundary_function_get_value_space 

  !===============================================================================================
  subroutine boundary_function_get_gradient_space ( this, point, result )
    implicit none
    class(boundary_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    real(rp) :: n 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    n = real(this%order, rp) 
    if ( this%get_num_dims() == 2 ) then
      call result%set( 1, 1.0_rp ) ! nx^{n-1}
      call result%set( 2, 1.0_rp ) ! ny^{n-1}
      ! call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      ! call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
    else if ( this%get_num_dims() == 3 ) then
      call result%set( 1, 1.0_rp ) ! nx^{n-1}
      call result%set( 2, 1.0_rp ) ! ny^{n-1}
      call result%set( 3, 1.0_rp ) ! nz^{n-1}
      ! call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      ! call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
      ! call result%set( 3, n*point%get(3)**(n-1.0_rp) ) ! nz^{n-1}
    end if
  end subroutine boundary_function_get_gradient_space
  
  !===============================================================================================
  subroutine solution_temperature_function_get_value_space ( this, point, result )
    implicit none
    class(solution_temperature_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1) + point%get(2) 
      ! result = point%get(1)**this%order + point%get(2)**this%order ! x^n+y^n 
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1) + point%get(2) + point%get(3)
      ! result = point%get(1)**this%order + point%get(2)**this%order + point%get(3)**this%order ! x^n+y^n+z^n 
    end if  
      
  end subroutine solution_temperature_function_get_value_space
  
  !===============================================================================================
  subroutine solution_temperature_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_temperature_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    real(rp) :: n 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    n = real(this%order, rp) 
    if ( this%get_num_dims() == 2 ) then
      call result%set( 1, 1.0_rp ) ! nx^{n-1}
      call result%set( 2, 1.0_rp ) ! ny^{n-1}
      ! call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      ! call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
    else if ( this%get_num_dims() == 3 ) then
      call result%set( 1, 1.0_rp ) ! nx^{n-1}
      call result%set( 2, 1.0_rp ) ! ny^{n-1}
      call result%set( 3, 1.0_rp ) ! nz^{n-1}
      ! call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      ! call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
      ! call result%set( 3, n*point%get(3)**(n-1.0_rp) ) ! nz^{n-1}
    end if
  end subroutine solution_temperature_function_get_gradient_space
  
  !===============================================================================================
  subroutine solution_pressure_function_get_value_space ( this, point, result )
    implicit none
    class(solution_pressure_function_t), intent(in) :: this
    type(point_t)             , intent(in)          :: point
    real(rp)                  , intent(inout)       :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1) - point%get(2) 
      ! result = point%get(1)**this%order + point%get(2)**this%order ! x^n+y^n 
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1) - point%get(2) - point%get(3)
      ! result = point%get(1)**this%order + point%get(2)**this%order + point%get(3)**this%order ! x^n+y^n+z^n 
    end if  
      
  end subroutine solution_pressure_function_get_value_space   
  
  !===============================================================================================
  subroutine solution_pressure_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_pressure_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    real(rp) :: n 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, -1.0_rp ) 
    else if ( this%get_num_dims() == 3 ) then
      call result%set( 1, 1.0_rp ) 
      call result%set( 2, -1.0_rp ) 
      call result%set( 3, -1.0_rp ) 
    end if
  end subroutine solution_pressure_function_get_gradient_space  
  
  !===============================================================================================  
  subroutine solution_displacement_function_get_value_space ( this, point, result )
    implicit none
    class(solution_displacement_function_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    if ( this%get_num_dims() == 2 ) then
      call result%set(1, point%get(1)-point%get(2) ) 
      call result%set(2, -point%get(1)+point%get(2) ) 
      call result%set(3,0.0_rp) 
    else
      call result%set(1, point%get(1)-point%get(2)-point%get(3) ) 
      call result%set(2, -point%get(1)+point%get(2)-point%get(3) ) 
      call result%set(3, -point%get(1)-point%get(2)+point%get(3) )
    end if
  end subroutine solution_displacement_function_get_value_space  
  
  !===============================================================================================  
  subroutine solution_displacement_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_displacement_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    if ( this%get_num_dims() == 2 ) then
      call result%init(0.0_rp)    
      call result%set( 1, 1, 1.0_rp ) 
      call result%set( 2, 1, -1.0_rp )
      call result%set( 1, 2, -1.0_rp ) 
      call result%set( 2, 2, 1.0_rp )
    else
      call result%init(-1.0_rp)
      call result%set( 1, 1, 1.0_rp ) 
      call result%set( 2, 2, 1.0_rp )
      call result%set( 3, 3, 1.0_rp ) 
    end if
  end subroutine solution_displacement_function_get_gradient_space  
  
  !===============================================================================================  
  subroutine solution_stress_function_get_value_space ( this, point, result )
    implicit none
    class(solution_stress_function_t), intent(in)  :: this
    type(point_t)           , intent(in)           :: point
    type(tensor_field_t)    , intent(inout)        :: result
    call result%init(0.0_rp)
    if ( this%get_num_dims() > 1) then
      call result%set(1, 1,  point%get(1)) 
      call result%set(2, 1,  point%get(1)-point%get(2)) 
      call result%set(1, 2, -point%get(1)+point%get(2))
      call result%set(2, 2,  point%get(2))    
    end if
    if ( this%get_num_dims() > 2) then
      call result%set(3, 1,  point%get(1)-point%get(3)) 
      call result%set(3, 2,  point%get(2)-point%get(3)) 
      call result%set(1, 3, -point%get(1)+point%get(3))
      call result%set(2, 3, -point%get(2)+point%get(3))
      call result%set(3, 3,  point%get(3))
    end if
  end subroutine solution_stress_function_get_value_space  
  
  !===============================================================================================
  subroutine lagrangian_fe_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(lagrangian_fe_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%boundary_function%set_num_dims(num_dims)
    call this%solution_temperature_function%set_num_dims(num_dims)
    call this%solution_pressure_function%set_num_dims(num_dims)
    call this%solution_displacement_function%set_num_dims(num_dims)
    call this%solution_stress_function%set_num_dims(num_dims)
  end subroutine lagrangian_fe_analytical_functions_set_num_dims 
  
  !===============================================================================================
  subroutine lfe_analytical_functions_set_solution_polynomial_order ( this, order )
    implicit none
    class(lagrangian_fe_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  order
    call this%source_term%set_polynomial_order(order)
    call this%boundary_function%set_polynomial_order(order)
    call this%solution_temperature_function%set_polynomial_order(order)
  end subroutine lfe_analytical_functions_set_solution_polynomial_order 
  
  !===============================================================================================
  function lagrangian_fe_analytical_functions_get_source_term ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: lagrangian_fe_analytical_functions_get_source_term
    lagrangian_fe_analytical_functions_get_source_term => this%source_term
  end function lagrangian_fe_analytical_functions_get_source_term
  
  !===============================================================================================
  function lagrangian_fe_analytical_functions_get_boundary_function ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: lagrangian_fe_analytical_functions_get_boundary_function
    lagrangian_fe_analytical_functions_get_boundary_function => this%boundary_function
  end function lagrangian_fe_analytical_functions_get_boundary_function
  
  !===============================================================================================
  function lagrangian_af_get_solution_temperature_function ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: lagrangian_af_get_solution_temperature_function
    lagrangian_af_get_solution_temperature_function => this%solution_temperature_function
  end function lagrangian_af_get_solution_temperature_function
  
  !===============================================================================================
  function lagrangian_af_get_solution_displacement_function ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: lagrangian_af_get_solution_displacement_function
    lagrangian_af_get_solution_displacement_function => this%solution_displacement_function
  end function lagrangian_af_get_solution_displacement_function  
  
  !===============================================================================================
  function lagrangian_af_get_solution_pressure_function ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: lagrangian_af_get_solution_pressure_function 
    lagrangian_af_get_solution_pressure_function  => this%solution_pressure_function 
  end function lagrangian_af_get_solution_pressure_function     
  
  !===============================================================================================
  function lagrangian_af_get_solution_stress_function ( this )
    implicit none
    class(lagrangian_fe_analytical_functions_t), target, intent(in)    :: this
    class(tensor_function_t), pointer :: lagrangian_af_get_solution_stress_function
    lagrangian_af_get_solution_stress_function => this%solution_stress_function
  end function lagrangian_af_get_solution_stress_function    

end module lagrangian_fe_analytical_functions_names
!***************************************************************************************************
