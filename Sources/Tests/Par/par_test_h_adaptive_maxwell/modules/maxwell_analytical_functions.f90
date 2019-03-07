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

module maxwell_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  character(len=*), parameter :: in_fe_space    ='in_fe_space'
  character(len=*), parameter :: fichera_2D     ='fichera_2D'
  character(len=*), parameter :: fichera_3D     ='fichera_3D'
  character(len=*), parameter :: sinusoidal     ='sinusoidal'

  public :: in_fe_space, fichera_2D, fichera_3D, sinusoidal 

  real(rp), parameter :: n     = 1.0_rp 
  real(rp), parameter :: alpha = 2.0_rp*n/3.0_rp 

  type, extends(vector_function_t) :: base_vector_function_t
     private
     character(len=:), allocatable :: function_case 
     integer(ip)                   :: degree
   contains
     procedure :: set_case        => base_vector_function_set_case 
     procedure :: set_degree      => base_vector_function_set_degree
  end type base_vector_function_t

  type, extends(scalar_function_t) :: base_scalar_function_t
     private
     character(len=:), allocatable :: function_case 
     integer(ip)                   :: degree
   contains
     procedure :: set_case        => base_scalar_function_set_case 
     procedure :: set_degree      => base_scalar_function_set_degree
  end type base_scalar_function_t

  type, extends(base_vector_function_t) :: source_term_t
     private 
   contains
     procedure :: get_value_space    => source_term_get_value_space
  end type source_term_t

  type, extends(base_scalar_function_t) :: boundary_function_Hx_t
     private
   contains
     procedure :: get_value_space => boundary_function_Hx_get_value_space
  end type boundary_function_Hx_t

  type, extends(base_scalar_function_t) :: boundary_function_Hy_t
     private
   contains
     procedure :: get_value_space => boundary_function_Hy_get_value_space
  end type boundary_function_Hy_t

  type, extends(base_scalar_function_t) :: boundary_function_Hz_t
     private
   contains
     procedure :: get_value_space => boundary_function_Hz_get_value_space
  end type boundary_function_Hz_t

  type, extends(base_vector_function_t) :: solution_function_t
     private 
   contains
     procedure :: get_value_space    => solution_function_get_value_space
     procedure :: get_gradient_space => solution_function_get_gradient_space
  end type solution_function_t

  type maxwell_analytical_functions_t
     private
     type(source_term_t)          :: source_term
     type(boundary_function_Hx_t) :: boundary_function_Hx
     type(boundary_function_Hy_t) :: boundary_function_Hy
     type(boundary_function_Hz_t) :: boundary_function_Hz
     type(solution_function_t)    :: solution_function
   contains
     procedure :: set_case                   => maxwell_analytical_functions_set_case 
     procedure :: set_degree                 => maxwell_analytical_functions_set_degree
     procedure :: set_num_dims               => maxwell_analytical_functions_set_num_dims
     procedure :: get_source_term            => maxwell_analytical_functions_get_source_term
     procedure :: get_boundary_function_Hx   => maxwell_analytical_functions_get_boundary_function_Hx
     procedure :: get_boundary_function_Hy   => maxwell_analytical_functions_get_boundary_function_Hy
     procedure :: get_boundary_function_Hz   => maxwell_analytical_functions_get_boundary_function_Hz
     procedure :: get_solution_function      => maxwell_analytical_functions_get_solution_function
  end type maxwell_analytical_functions_t

  public :: maxwell_analytical_functions_t

contains  

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
  
  subroutine base_scalar_function_set_degree ( this, degree )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip)                  , intent(in)       :: degree
    this%degree = degree
  end subroutine base_scalar_function_set_degree

  subroutine base_vector_function_set_degree ( this, degree )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip)                  , intent(in)       :: degree
    this%degree = degree
  end subroutine base_vector_function_set_degree

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    real(rp) :: x,y,z,k
    real(rp) :: r,theta 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )

    x = point%get(1); y = point%get(2); z = point%get(3)  
    call result%init(0.0_rp)  
           
    select case ( this%function_case ) 
    case ( in_fe_space )
       k = real(this%degree,kind=rp) 
       call result%set(1,   k*y**(k-2.0_rp)*(k-1.0_rp) - y**k )
       call result%set(2, - k*x**(k-2.0_rp)*(k-1.0_rp) + x**k ) 
       if ( this%get_num_dims() == 3 ) then
       call result%add(1, k*z**(k-2.0_rp)*(k-1.0_rp) - z**k)
       call result%add(2, k*z**(k-2.0_rp)*(k-1.0_rp) - z**k) 
       call result%set(3, -k*x**(k-2.0_rp)*(k-1.0_rp) - k*y**(k-2.0_rp)*(k-1.0_rp) + x**k + y**k) 
       end if
       
    case ( fichera_2D ) 
       r = sqrt(x*x+y*y) 
       theta = compute_polar_angle(x,y) 
       call result%set(1, alpha*r**(alpha-2.0_rp)*(x*sin(alpha*theta)-y*cos(alpha*theta)) ) 
       call result%set(2, alpha*r**(alpha-2.0_rp)*(y*sin(alpha*theta)+x*cos(alpha*theta)) ) 

    case ( fichera_3D )       
       call result%set(1, alpha*x*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x**2.0_rp*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
       call result%set(2, alpha*y*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y**2.0_rp*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
       call result%set(3, alpha*z*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*y)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
    
    case ( sinusoidal ) 
       call result%set(1, (pi**2+1.0_rp)*cos(pi*x)*cos(pi*y) + pi**2*sin(pi*x)*sin(pi*z) )
       call result%set(2, (pi**2+1.0_rp)*sin(pi*y)*sin(pi*z) + pi**2*sin(pi*x)*sin(pi*y) ) 
       call result%set(3, (pi**2+1.0_rp)*cos(pi*x)*cos(pi*z) + pi**2*cos(pi*y)*cos(pi*z) ) 
       
    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine boundary_function_Hx_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_Hx_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    real(rp) :: x,y,z,k
    real(rp) :: r,theta 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x = point%get(1); y = point%get(2); z = point%get(3) 

    select case ( this%function_case ) 
    case (in_fe_space) 
       k = real(this%degree,kind=rp) 
       result = -y**k 
       if ( this%get_num_dims() == 3 ) result = result - z**k  
    case ( fichera_2D ) 
       r = sqrt(x*x+y*y) 
       theta = compute_polar_angle(x,y) 
       result = alpha*r**(alpha-2.0_rp )*(x*sin(alpha*theta)-y*cos(alpha*theta)) 
    case ( fichera_3D ) 
       result = alpha*x*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x**2.0_rp*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) 
    case ( sinusoidal ) 
       result = cos(pi*x)*cos(pi*y) 

    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select

  end subroutine boundary_function_Hx_get_value_space
  
  subroutine boundary_function_Hy_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_Hy_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    real(rp) :: x,y,z, k
    real(rp) :: r,theta 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x = point%get(1); y = point%get(2); z = point%get(3) 

    select case ( this%function_case ) 
    case (in_fe_space) 
       k = real(this%degree,kind=rp) 
       result = x**k
       if ( this%get_num_dims()==3 ) result = result - z**k
    case ( fichera_2D ) 
       r = sqrt(x*x+y*y) 
       theta = compute_polar_angle(x,y) 
       result = alpha*r**(alpha-2.0_rp)*(y*sin(alpha*theta)+x*cos(alpha*theta))
    case ( fichera_3D ) 
       result = alpha*y*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y**2.0_rp*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp)  
    case ( sinusoidal ) 
       result = sin(pi*y)*sin(pi*z) 
    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select
  end subroutine boundary_function_Hy_get_value_space
  
    !===============================================================================================
  subroutine boundary_function_Hz_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_Hz_t), intent(in)  :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    real(rp) :: x,y,z, k
    real(rp) :: r,theta 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x = point%get(1); y = point%get(2); z = point%get(3) 
    select case ( this%function_case ) 
    case (in_fe_space) 
       k = real(this%degree,kind=rp)
       if ( this%get_num_dims()==3) result = x**k + y**k
    case ( fichera_2D ) 
       r = sqrt(x*x+y*y) 
       theta = compute_polar_angle(x,y) 
       result = 0.0_rp 
    case ( fichera_3D ) 
       result = alpha*z*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*y)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp)
    case ( sinusoidal ) 
       result = cos(pi*x)*cos(pi*z) 
    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select
  end subroutine boundary_function_Hz_get_value_space

  !===============================================================================================
  subroutine solution_function_get_value_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    real(rp) :: x,y,z, k
    real(rp) :: r,theta 
    real(rp) :: t,A,C
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )

    x = point%get(1); y = point%get(2); z = point%get(3)  
    call result%init(0.0_rp) 
    
    select case ( this%function_case ) 
    case ( in_fe_space ) ! u = [ -y^k - z^k, x^k - z^k, x^k + y^k ]
       k = real(this%degree, kind=rp) 
       call result%set(1, -y**k )
       call result%set(2,  x**k ) 
       if ( this%get_num_dims() == 3 ) then
       call result%add(1, - z**k )
       call result%add(2, - z**k )
       call result%set(3, x**k + y**k) 
       end if

    case ( fichera_2D ) ! u = grad(w); w = r^alpha*sin(alpha*theta)
       r = sqrt(x*x+y*y) 
       theta = compute_polar_angle(x,y) 
       call result%init(0.0_rp) 
       call result%set(1, alpha*r**(alpha-2.0_rp)*(x*sin(alpha*theta)-y*cos(alpha*theta)) ) 
       call result%set(2, alpha*r**(alpha-2.0_rp)*(y*sin(alpha*theta)+x*cos(alpha*theta)) ) 
       if ( this%get_num_dims() == 3 ) then 
          call result%set(3, 0.0_rp) 
       end if
       
    case ( fichera_3D ) ! u = grad(w); w = r^alpha*sin(alpha*t), t = arccos(x*y*z/r)          
       call result%set(1, alpha*x*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x**2.0_rp*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
       call result%set(2, alpha*y*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y**2.0_rp*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
       call result%set(3, alpha*z*sin(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp - 1.0_rp) - (alpha*cos(alpha*acos((x*y*z)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp)))*((x*y)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(1.0_rp/2.0_rp) - (x*y*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(3.0_rp/2.0_rp))*(x**2.0_rp + y**2.0_rp + z**2.0_rp)**(alpha/2.0_rp))/(1.0_rp - (x**2.0_rp*y**2.0_rp*z**2.0_rp)/(x**2.0_rp + y**2.0_rp + z**2.0_rp))**(1.0_rp/2.0_rp) )
   
    case ( sinusoidal ) 
       call result%set(1, cos(pi*x)*cos(pi*y) ) 
       call result%set(2, sin(pi*y)*sin(pi*z) ) 
       call result%set(3, cos(pi*x)*cos(pi*z) )

    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select
  end subroutine solution_function_get_value_space

  !===============================================================================================
  subroutine solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    real(rp) :: x,y,z,k,r,t  
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )

    x = point%get(1); y = point%get(2); z = point%get(3)  
    call result%init(0.0_rp) 
    
    select case ( this%function_case ) 
    case ( in_fe_space ) 
       k = real(this%degree, kind=rp) 
       call result%set( 1, 2,  k*x**(k-1.0_rp) )
       call result%set( 2, 1, -k*y**(k-1.0_rp) )
       if ( this%get_num_dims()==3) then 
       call result%set( 1, 3,  k*x**(k-1.0_rp) )
       call result%set( 2, 3,  k*y**(k-1.0_rp) )
       call result%set( 3, 1, -k*z**(k-1.0_rp) ) 
       call result%set( 3, 2, -k*z**(k-1.0_rp) )
       end if 

    case ( fichera_2D, fichera_3D ) 
       ! Gradients are not null, but are needed to evaluate curl norm 
       ! and the curl of any gradient is curl(u) = 0; given u = grad(w);
    case ( sinusoidal ) 
      call result%init(0.0_rp) 
      call result%set(1, 1, -pi*sin(pi*x)*cos(pi*y) ) 
      call result%set(2, 1, -pi*cos(pi*x)*sin(pi*y) )
      
      call result%set(2, 2, pi*cos(pi*y)*sin(pi*z)  ) 
      call result%set(3, 2, pi*sin(pi*y)*cos(pi*z) ) 

      call result%set(1, 3, -pi*sin(pi*x)*cos(pi*z) ) 
      call result%set(3, 3, -pi*cos(pi*x)*sin(pi*z) ) 
     
    case DEFAULT 
       massert(.false., 'maxwell_analytical_funcions: case is not defined')
    end select
  end subroutine solution_function_get_gradient_space

  !===============================================================================================
  subroutine maxwell_analytical_functions_set_case ( this, function_case )
    implicit none
    class(maxwell_analytical_functions_t), intent(inout) :: this
    character(len=*)                     , intent(in)    :: function_case 
    call this%source_term%set_case(function_case)
    call this%boundary_function_Hx%set_case(function_case)
    call this%boundary_function_Hy%set_case(function_case)
    call this%boundary_function_Hz%set_case(function_case)
    call this%solution_function%set_case(function_case)
  end subroutine maxwell_analytical_functions_set_case
  
    !===============================================================================================
  subroutine maxwell_analytical_functions_set_degree ( this, degree )
    implicit none
    class(maxwell_analytical_functions_t), intent(inout) :: this
    integer(ip)                          , intent(in)    :: degree
    call this%source_term%set_degree(degree)
    call this%boundary_function_Hx%set_degree(degree)
    call this%boundary_function_Hy%set_degree(degree)
    call this%boundary_function_Hz%set_degree(degree)
    call this%solution_function%set_degree(degree)
  end subroutine maxwell_analytical_functions_set_degree

  !===============================================================================================
  subroutine maxwell_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(maxwell_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%boundary_function_Hx%set_num_dims(num_dims)
    call this%boundary_function_Hy%set_num_dims(num_dims)
    call this%boundary_function_Hz%set_num_dims(num_dims)
    call this%solution_function%set_num_dims(num_dims)
  end subroutine maxwell_analytical_functions_set_num_dims

  !===============================================================================================
  function maxwell_analytical_functions_get_source_term ( this )
    implicit none
    class(maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: maxwell_analytical_functions_get_source_term
    maxwell_analytical_functions_get_source_term => this%source_term
  end function maxwell_analytical_functions_get_source_term

  !===============================================================================================
  function maxwell_analytical_functions_get_boundary_function_Hx ( this )
    implicit none
    class(maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: maxwell_analytical_functions_get_boundary_function_Hx
    maxwell_analytical_functions_get_boundary_function_Hx => this%boundary_function_Hx
  end function maxwell_analytical_functions_get_boundary_function_Hx

  !===============================================================================================
  function maxwell_analytical_functions_get_boundary_function_Hy ( this )
    implicit none
    class(maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: maxwell_analytical_functions_get_boundary_function_Hy
    maxwell_analytical_functions_get_boundary_function_Hy => this%boundary_function_Hy
  end function maxwell_analytical_functions_get_boundary_function_Hy

  !===============================================================================================
  function maxwell_analytical_functions_get_boundary_function_Hz ( this )
    implicit none
    class(maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: maxwell_analytical_functions_get_boundary_function_Hz
    maxwell_analytical_functions_get_boundary_function_Hz => this%boundary_function_Hz
  end function maxwell_analytical_functions_get_boundary_function_Hz

  !===============================================================================================
  function maxwell_analytical_functions_get_solution_function ( this )
    implicit none
    class(maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: maxwell_analytical_functions_get_solution_function
    maxwell_analytical_functions_get_solution_function => this%solution_function
  end function maxwell_analytical_functions_get_solution_function

  ! Auxiliar functions 
  ! ==============================================================================================
  function compute_polar_angle(x,y) 
    implicit none 
    real(rp)  , intent(in) :: x 
    real(rp)  , intent(in) :: y 
    real(rp) :: compute_polar_angle 
    real(rp) :: beta

    beta = atan(y/x) 
    if (x > 0.0_rp .and. y > 0.0_rp ) then 
       compute_polar_angle = beta
    elseif ( x < 0.0_rp .and. y > 0.0_rp ) then
       compute_polar_angle= beta + pi
    elseif ( x < 0.0_rp .and. y == 0.0_rp ) then 
       compute_polar_angle = pi 
    elseif ( x < 0.0_rp .and. y < 0.0_rp ) then 
       compute_polar_angle = beta + pi  
    elseif ( x> 0.0_rp .and. y< 0.0_rp ) then 
       compute_polar_angle = 2.0_rp*pi + beta
    elseif ( x> 0.0_rp .and. y==0.0_rp ) then 
       compute_polar_angle = 0.0_rp 
    elseif ( x==0.0_rp ) then 
       if ( y > 0.0_rp ) then 
          compute_polar_angle = pi/2.0_rp 
       elseif (y < 0.0_rp ) then 
          compute_polar_angle = 3.0_rp*pi/2.0_rp 
       elseif ( y==0.0_rp) then 
          compute_polar_angle = 0.0_rp 
          wassert(.false., 'Trying to evaluate theta at x=0, y=0: Set it to 0.0') 
       end if
    end if

  end function compute_polar_angle

end module maxwell_analytical_functions_names
!***************************************************************************************************
