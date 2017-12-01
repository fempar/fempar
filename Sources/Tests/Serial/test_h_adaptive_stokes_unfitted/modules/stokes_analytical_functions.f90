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

module stokes_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  !===================================================

  type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims    = -1  
    logical     :: in_fe_space = .true.
    integer(ip) :: degree      = 1
  contains
    procedure :: set_num_dims          => base_scalar_function_set_num_dims 
    procedure :: set_is_in_fe_space    => base_scalar_function_set_is_in_fe_space
    procedure :: set_degree            => base_scalar_function_set_degree
  end type base_scalar_function_t

  type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dims    = -1  
    logical     :: in_fe_space = .true.
    integer(ip) :: degree      = 1
  contains
    procedure :: set_num_dims          => base_vector_function_set_num_dims 
    procedure :: set_is_in_fe_space    => base_vector_function_set_is_in_fe_space
    procedure :: set_degree            => base_vector_function_set_degree
  end type base_vector_function_t

  !===================================================

  type, extends(base_vector_function_t) :: source_term_u_t
    private 
   contains
     procedure :: get_value_space    => source_term_u_get_value_space
  end type source_term_u_t

  type, extends(base_vector_function_t) :: solution_function_u_t
    private 
   contains
     procedure :: get_value_space    => solution_function_u_get_value_space
     procedure :: get_gradient_space => solution_function_u_get_gradient_space
  end type solution_function_u_t
  
  type, extends(base_scalar_function_t) :: source_term_p_t
    private 
   contains
     procedure :: get_value_space    => source_term_p_get_value_space
  end type source_term_p_t

  type, extends(base_scalar_function_t) :: solution_function_p_t
    private 
   contains
     procedure :: get_value_space    => solution_function_p_get_value_space
     procedure :: get_gradient_space => solution_function_p_get_gradient_space
  end type solution_function_p_t

  !===================================================

  type stokes_analytical_functions_t
     private
     type(source_term_u_t)         :: source_term_u
     type(source_term_p_t)         :: source_term_p
     type(solution_function_u_t)   :: solution_function_u
     type(solution_function_p_t)   :: solution_function_p
   contains
     procedure :: set_num_dims            => stokes_analytical_functions_set_num_dims
     procedure :: set_is_in_fe_space      => stokes_analytical_functions_set_is_in_fe_space
     procedure :: set_degree              => stokes_analytical_functions_set_degree
     procedure :: get_degree              => stokes_analytical_functions_get_degree
     procedure :: get_source_term_u       => stokes_analytical_functions_get_source_term_u
     procedure :: get_source_term_p       => stokes_analytical_functions_get_source_term_p
     procedure :: get_solution_function_u => stokes_analytical_functions_get_solution_function_u
     procedure :: get_solution_function_p => stokes_analytical_functions_get_solution_function_p
  end type stokes_analytical_functions_t

  public :: stokes_analytical_functions_t, base_scalar_function_t

contains  

  !===================================================

  subroutine sol_ex001_2d_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1, x1**q + x1**q*x2**q )
    call val%set(2, x2**q + x1**q*x2**q )
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
    call val%set(1,1, q*x1**(q - 1) + q*x1**(q - 1)*x2**q  )
    call val%set(1,2, q*x1**(q - 1)*x2**q                 )
    call val%set(2,1, q*x1**q*x2**(q - 1)                 )
    call val%set(2,2, q*x2**(q - 1) + q*x1**q*x2**(q - 1)  )
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
    call val%set(1, q*x1**(q - 2)*(q - 1) + q**2*x1**(q - 1)*x2**(q - 1) + q*x1**(q - 2)*x2**q*(q - 1) )
    call val%set(2, q*x2**(q - 2)*(q - 1) + q**2*x1**(q - 1)*x2**(q - 1) + q*x1**q*x2**(q - 2)*(q - 1) )
  end subroutine sol_ex001_2d_lapl_u

  subroutine sol_ex001_2d_div_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    val = q*x1**(q - 1) + q*x2**(q - 1) + q*x1**q*x2**(q - 1) + q*x1**(q - 1)*x2**q
  end subroutine sol_ex001_2d_div_u

  !===================================================

  subroutine sol_ex001_2d_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    val = x1**(q - 1)*x2**(q - 1)
  end subroutine sol_ex001_2d_p

  subroutine sol_ex001_2d_grad_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1, x1**(q - 2)*x2**(q - 1)*(q - 1) )
    call val%set(2, x1**(q - 1)*x2**(q - 2)*(q - 1) )
  end subroutine sol_ex001_2d_grad_p

  !===================================================

  subroutine base_scalar_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_scalar_function_set_num_dims

  subroutine base_scalar_function_set_is_in_fe_space(this,val)
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    logical,                       intent(in)       :: val
    this%in_fe_space = val
  end subroutine base_scalar_function_set_is_in_fe_space

  subroutine base_scalar_function_set_degree(this,degree)
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip),                   intent(in)       :: degree
    this%degree = degree
  end subroutine base_scalar_function_set_degree

  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims

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

  !===================================================

  subroutine source_term_u_get_value_space ( this, point, result )
    implicit none
    class(source_term_u_t), intent(in)    :: this
    type(point_t)         , intent(in)    :: point
    type(vector_field_t)  , intent(inout) :: result
    type(vector_field_t) :: val
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_lapl_u(point,val,this%degree)
        result = (-1.0)*val
        call sol_ex001_2d_grad_p(point,val,this%degree)
        result = result + val
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine source_term_u_get_value_space

  subroutine solution_function_u_get_value_space ( this, point, result )
    implicit none
    class(solution_function_u_t), intent(in)    :: this
    type(point_t)         , intent(in)    :: point
    type(vector_field_t)  , intent(inout) :: result
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_u(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_u_get_value_space
  
  subroutine solution_function_u_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_u_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(tensor_field_t)      , intent(inout) :: result
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_grad_u(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_u_get_gradient_space

  subroutine source_term_p_get_value_space ( this, point, result )
    implicit none
    class(source_term_p_t), intent(in)    :: this
    type(point_t)         , intent(in)    :: point
    real(rp)              , intent(inout) :: result
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_div_u(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine source_term_p_get_value_space

  subroutine solution_function_p_get_value_space ( this, point, result )
    implicit none
    class(solution_function_p_t), intent(in)    :: this
    type(point_t)         , intent(in)    :: point
    real(rp)              , intent(inout) :: result
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_p(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_p_get_value_space
  
  subroutine solution_function_p_get_gradient_space ( this, point, result )
    implicit none
    class(solution_function_p_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    if ( this%num_dims == 2 ) then
      if (this%in_fe_space) then
        call sol_ex001_2d_grad_p(point,result,this%degree)
      else
        check(.false.)
      end if
    else if ( this%num_dims == 3 ) then
      if (this%in_fe_space) then
        check(.false.)
      else
        check(.false.)
      end if
    else
        check(.false.)
    end if  
  end subroutine solution_function_p_get_gradient_space

  !===================================================
  
  subroutine stokes_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(stokes_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term_u%set_num_dims(num_dims)
    call this%source_term_p%set_num_dims(num_dims)
    call this%solution_function_u%set_num_dims(num_dims)
    call this%solution_function_p%set_num_dims(num_dims)
  end subroutine stokes_analytical_functions_set_num_dims 

  subroutine stokes_analytical_functions_set_is_in_fe_space(this,val)
    implicit none
    class(stokes_analytical_functions_t), intent(inout)    :: this
    logical, intent(in) :: val
    call this%source_term_u%set_is_in_fe_space(val)
    call this%source_term_p%set_is_in_fe_space(val)
    call this%solution_function_u%set_is_in_fe_space(val)
    call this%solution_function_p%set_is_in_fe_space(val)
  end subroutine stokes_analytical_functions_set_is_in_fe_space

  subroutine stokes_analytical_functions_set_degree(this,degree)
    implicit none
    class(stokes_analytical_functions_t), intent(inout) :: this
    integer(ip),                                    intent(in)    :: degree
    call this%source_term_u%set_degree(degree)
    call this%source_term_p%set_degree(degree)
    call this%solution_function_u%set_degree(degree)
    call this%solution_function_p%set_degree(degree)
  end subroutine stokes_analytical_functions_set_degree

  function stokes_analytical_functions_get_degree(this)
    implicit none
    class(stokes_analytical_functions_t), intent(in)    :: this
    integer(ip) :: stokes_analytical_functions_get_degree
    stokes_analytical_functions_get_degree = this%source_term_u%degree
  end function stokes_analytical_functions_get_degree

  function stokes_analytical_functions_get_source_term_u ( this )
    implicit none
    class(stokes_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: stokes_analytical_functions_get_source_term_u
    stokes_analytical_functions_get_source_term_u => this%source_term_u
  end function stokes_analytical_functions_get_source_term_u
  
  function stokes_analytical_functions_get_source_term_p ( this )
    implicit none
    class(stokes_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: stokes_analytical_functions_get_source_term_p
    stokes_analytical_functions_get_source_term_p => this%source_term_p
  end function stokes_analytical_functions_get_source_term_p

  function stokes_analytical_functions_get_solution_function_u ( this )
    implicit none
    class(stokes_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: stokes_analytical_functions_get_solution_function_u
    stokes_analytical_functions_get_solution_function_u => this%solution_function_u
  end function stokes_analytical_functions_get_solution_function_u
  
  function stokes_analytical_functions_get_solution_function_p ( this )
    implicit none
    class(stokes_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: stokes_analytical_functions_get_solution_function_p
    stokes_analytical_functions_get_solution_function_p => this%solution_function_p
  end function stokes_analytical_functions_get_solution_function_p
  

end module stokes_analytical_functions_names
!***************************************************************************************************
