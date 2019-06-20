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
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(vector_function_t) :: base_vector_function_t
    logical     :: in_fe_space = .false.
    integer(ip) :: degree      = 1
  contains
    procedure :: set_is_in_fe_space    => base_vector_function_set_is_in_fe_space
    procedure :: set_degree            => base_vector_function_set_degree
  end type base_vector_function_t
  
  type, extends(base_vector_function_t) :: source_term_t
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t

  type, extends(base_vector_function_t) :: solution_function_t
   contains
     procedure :: get_value_space    => solution_function_get_value_space
     procedure :: get_gradient_space => solution_function_get_gradient_space
  end type solution_function_t

  type vector_poisson_analytical_functions_t
     private
     type(source_term_t)         :: source_term
     type(solution_function_t)   :: solution_function
   contains
     procedure :: set_num_dims            => poisson_analytical_functions_set_num_dims
     procedure :: set_is_in_fe_space      => poisson_analytical_functions_set_is_in_fe_space
     procedure :: set_degree              => poisson_analytical_functions_set_degree
     procedure :: get_source_term         => poisson_analytical_functions_get_source_term
     procedure :: get_solution_function   => poisson_analytical_functions_get_solution_function
  end type vector_poisson_analytical_functions_t

  public :: vector_poisson_analytical_functions_t

contains  

  subroutine sol_ex001_2d_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1, x2 + x1*x2 )
    call val%set(2, x1 + x1*x2 )
  end subroutine sol_ex001_2d_u

  subroutine sol_ex001_2d_grad_u(point,val,q)
    implicit none
    type(point_t),        intent(in)    :: point
    type(tensor_field_t), intent(inout) :: val
    integer(ip),          intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1,1, x2)
    call val%set(1,2, x2 + 1.0)
    call val%set(2,1, x1 + 1.0)
    call val%set(2,2, x1)
  end subroutine sol_ex001_2d_grad_u

  subroutine sol_ex001_2d_lapl_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1, 0.0 )
    call val%set(2, 0.0 )
  end subroutine sol_ex001_2d_lapl_u

  subroutine base_vector_function_set_is_in_fe_space(this,val)
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    logical,                       intent(in)       :: val
    this%in_fe_space = val
  end subroutine base_vector_function_set_is_in_fe_space

  subroutine base_vector_function_set_degree(this,degree)
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip),                   intent(in)       :: degree
    this%degree = degree
  end subroutine base_vector_function_set_degree

  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    type(vector_field_t), intent(inout) :: result
    if ( this%get_num_dims() == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_lapl_u(point,result,this%degree)
        result = (-1.0)*result
      else
        check(.false.)
      end if
    else if ( this%get_num_dims() == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine solution_function_get_value_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    if ( this%get_num_dims() == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_u(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%get_num_dims() == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_get_value_space
  
  !===============================================================================================
  subroutine solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    if ( this%get_num_dims() == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_grad_u(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%get_num_dims() == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_get_gradient_space
  
  !===============================================================================================
  subroutine poisson_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(vector_poisson_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%solution_function%set_num_dims(num_dims)
  end subroutine poisson_analytical_functions_set_num_dims 

  !===============================================================================================
  subroutine poisson_analytical_functions_set_is_in_fe_space(this,val)
    implicit none
    class(vector_poisson_analytical_functions_t), intent(inout)    :: this
    logical, intent(in) :: val
    call this%source_term      %set_is_in_fe_space(val)
    call this%solution_function%set_is_in_fe_space(val)
  end subroutine poisson_analytical_functions_set_is_in_fe_space

  !===============================================================================================
  subroutine poisson_analytical_functions_set_degree(this,degree)
    implicit none
    class(vector_poisson_analytical_functions_t), intent(inout) :: this
    integer(ip),                                    intent(in)    :: degree
    call this%source_term      %set_degree(degree)
    call this%solution_function%set_degree(degree)
  end subroutine poisson_analytical_functions_set_degree
  
  !===============================================================================================
  function poisson_analytical_functions_get_source_term ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: poisson_analytical_functions_get_source_term
    poisson_analytical_functions_get_source_term => this%source_term
  end function poisson_analytical_functions_get_source_term
  
  !===============================================================================================
  function poisson_analytical_functions_get_solution_function ( this )
    implicit none
    class(vector_poisson_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: poisson_analytical_functions_get_solution_function
    poisson_analytical_functions_get_solution_function => this%solution_function
  end function poisson_analytical_functions_get_solution_function

end module vector_poisson_analytical_functions_names
!***************************************************************************************************


