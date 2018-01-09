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

  real(rp), parameter :: offset_x = -0.3
  real(rp), parameter :: offset_y = 0.5
  real(rp), parameter :: offset_z = 0.0

  !===================================================

  type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims    = -1  
    integer(ip) :: case_id = 1
    integer(ip) :: degree      = 1
  contains
    procedure :: set_num_dims          => base_scalar_function_set_num_dims 
    procedure :: set_case_id    => base_scalar_function_set_case_id
    procedure :: set_degree            => base_scalar_function_set_degree
  end type base_scalar_function_t

  type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dims    = -1  
    integer(ip) :: case_id = 1
    integer(ip) :: degree      = 1
  contains
    procedure :: set_num_dims          => base_vector_function_set_num_dims 
    procedure :: set_case_id    => base_vector_function_set_case_id
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
     procedure :: set_case_id             => stokes_analytical_functions_set_case_id
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

  subroutine sol_ex001_2d_div_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    val = x1 + x2
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
    val = x1*x2
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
    call val%set(1, x2 )
    call val%set(2, x1 )
  end subroutine sol_ex001_2d_grad_p

  !===================================================

  subroutine sol_ex002_2d_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    call val%init(0.0)
    call val%set(1, -x2/(x1**2 + x2**2)**(1.0/2.0) )
    call val%set(2,  x1/(x1**2 + x2**2)**(1.0/2.0) )
  end subroutine sol_ex002_2d_u

  subroutine sol_ex002_2d_grad_u(point,val,q)
    implicit none
    type(point_t),        intent(in)    :: point
    type(tensor_field_t), intent(inout) :: val
    integer(ip),          intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    call val%init(0.0)
    call val%set(1,1,  (x1*x2)/(x1**2 + x2**2)**(3.0/2.0)                                 )
    call val%set(1,2, 1.0/(x1**2 + x2**2)**(1.0/2.0) - x1**2/(x1**2 + x2**2)**(3.0/2.0)   )
    call val%set(2,1,   x2**2/(x1**2 + x2**2)**(3.0/2.0) - 1.0/(x1**2 + x2**2)**(1.0/2.0) )
    call val%set(2,2, -(x1*x2)/(x1**2 + x2**2)**(3.0/2.0)                                 )
  end subroutine sol_ex002_2d_grad_u

  subroutine sol_ex002_2d_lapl_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    call val%init(0.0)
    call val%set(1, (4.0*x2)/(x1**2 + x2**2)**(3.0/2.0) - (3.0*x2**3.0)/(x1**2 + x2**2)**(5.0/2.0) - (3.0*x1**2*x2)/(x1**2 + x2**2)**(5.0/2.0) )
    call val%set(2, (3.0*x1**3.0)/(x1**2 + x2**2)**(5.0/2.0) - (4.0*x1)/(x1**2 + x2**2)**(3.0/2.0) + (3.0*x1*x2**2)/(x1**2 + x2**2)**(5.0/2.0) )
  end subroutine sol_ex002_2d_lapl_u

  subroutine sol_ex002_2d_div_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    val = 0.0
  end subroutine sol_ex002_2d_div_u

  !===================================================

  subroutine sol_ex002_2d_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    val = x1**3*x2**3
  end subroutine sol_ex002_2d_p

  subroutine sol_ex002_2d_grad_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2
    x1 = point%get(1)
    x2 = point%get(2)
    call val%init(0.0)
    call val%set(1, 3.0*x1**2*x2**3 )
    call val%set(2, 3.0*x1**3*x2**2 )
  end subroutine sol_ex002_2d_grad_p

  !===================================================

  subroutine sol_ex001_3d_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    call val%init(0.0)
    call val%set(1, x2*x3 )
    call val%set(2, x1*x3 )
    call val%set(3, x1*x2 )
  end subroutine sol_ex001_3d_u

  subroutine sol_ex001_3d_grad_u(point,val,q)
    implicit none
    type(point_t),        intent(in)    :: point
    type(tensor_field_t), intent(inout) :: val
    integer(ip),          intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    call val%init(0.0)
    call val%set(1,1, 0.0)
    call val%set(1,2, x3 )
    call val%set(1,3, x2 )
    call val%set(2,1, x3 )
    call val%set(2,2, 0.0)
    call val%set(2,3, x1 )
    call val%set(3,1, x2 )
    call val%set(3,2, x1 )
    call val%set(3,3, 0.0)
  end subroutine sol_ex001_3d_grad_u

  subroutine sol_ex001_3d_lapl_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    call val%init(0.0)
    call val%set(1, 0.0 )
    call val%set(2, 0.0 )
    call val%set(3, 0.0 )
  end subroutine sol_ex001_3d_lapl_u

  subroutine sol_ex001_3d_div_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    val = 0.0
  end subroutine sol_ex001_3d_div_u

  !===================================================

  subroutine sol_ex001_3d_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    val = x1*x2
  end subroutine sol_ex001_3d_p

  subroutine sol_ex001_3d_grad_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    call val%init(0.0)
    call val%set(1, x2 )
    call val%set(2, x1 )
    call val%set(3, 0.0)
  end subroutine sol_ex001_3d_grad_p

  !===================================================

  subroutine sol_ex002_3d_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    x3 = point%get(3) - offset_z
    call val%init(0.0)
    call val%set(1,  x2/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0)         )
    call val%set(2,  -(x1 + x3)/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0) )
    call val%set(3,  x2/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0)         )
  end subroutine sol_ex002_3d_u

  subroutine sol_ex002_3d_grad_u(point,val,q)
    implicit none
    type(point_t),        intent(in)    :: point
    type(tensor_field_t), intent(inout) :: val
    integer(ip),          intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    x3 = point%get(3) - offset_z
    call val%init(0.0)
    call val%set(1,1,  -(x2*(2.0*x1 + 2.0*x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0))                                                   )
    call val%set(1,2,  ((2.0*x1 + 2.0*x3)*(x1 + x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0)) - 1.0/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0) )
    call val%set(1,3,  -(x2*(2.0*x1 + 2.0*x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0))                                                   )
    call val%set(2,1,  1.0/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0) - (2.0*x2**2)/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0)                         )
    call val%set(2,2,  (2.0*x2*(x1 + x3))/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0)                                                              )
    call val%set(2,3,  1.0/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0) - (2.0*x2**2)/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0)                         )
    call val%set(3,1,  -(x2*(2.0*x1 + 2.0*x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0))                                                   )
    call val%set(3,2,  ((2.0*x1 + 2.0*x3)*(x1 + x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0)) - 1.0/((x1 + x3)**2 + 2.0*x2**2)**(1.0/2.0) )
    call val%set(3,3,  -(x2*(2.0*x1 + 2.0*x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0))                                                   )
  end subroutine sol_ex002_3d_grad_u

  subroutine sol_ex002_3d_lapl_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    x3 = point%get(3) - offset_z
    call val%init(0.0)
    call val%set(1, (12*x2**3)/((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0) - (8*x2)/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0) + (3.0*x2*(2.0*x1 + 2.0*x3)**2)/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0))                                                                                           )
    call val%set(2, (2.0*(2.0*x1 + 2.0*x3))/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0) + (4.0*(x1 + x3))/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0) - (12*x2**2*(x1 + x3))/((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0) - (3.0*(2.0*x1 + 2.0*x3)**2*(x1 + x3))/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0)) )
    call val%set(3, (12*x2**3)/((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0) - (8*x2)/((x1 + x3)**2 + 2.0*x2**2)**(3.0/2.0) + (3.0*x2*(2.0*x1 + 2.0*x3)**2)/(2.0*((x1 + x3)**2 + 2.0*x2**2)**(5.0/2.0))                                                                                           )
  end subroutine sol_ex002_3d_lapl_u

  subroutine sol_ex002_3d_div_u(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1) - offset_x
    x2 = point%get(2) - offset_y
    x3 = point%get(3) - offset_z
    val = 0.0
  end subroutine sol_ex002_3d_div_u

  !===================================================

  subroutine sol_ex002_3d_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    val = x1**3*x2**3
  end subroutine sol_ex002_3d_p

  subroutine sol_ex002_3d_grad_p(point,val,q)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    integer(ip),   intent(in)    :: q
    real(rp) :: x1, x2, x3
    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)
    call val%init(0.0)
    call val%set(1, 3.0*x1**2*x2**3 )
    call val%set(2, 3.0*x1**3*x2**2 )
    call val%set(3, 0.0 )
  end subroutine sol_ex002_3d_grad_p

  !===================================================

  subroutine sol_ex003_3d_u(point,val)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    real(rp) :: x, y, z
    real(rp), parameter :: tol = 1.0e-10
    x = point%get(1)
    y = point%get(2)
    z = point%get(3)
    call val%init(0.0)
    if (abs(x)<=tol .or. abs(x-1)<=tol) then
      call val%set(1, 10*y*(y-1)*z*(z-1) )
    end if
  end subroutine sol_ex003_3d_u

  subroutine sol_ex003_3d_grad_u(point,val)
    implicit none
    type(point_t),        intent(in)    :: point
    type(tensor_field_t), intent(inout) :: val
    call val%init(0.0)
  end subroutine sol_ex003_3d_grad_u

  subroutine sol_ex003_3d_lapl_u(point,val)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    call val%init(0.0)
  end subroutine sol_ex003_3d_lapl_u

  subroutine sol_ex003_3d_div_u(point,val)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    val = 0.0
  end subroutine sol_ex003_3d_div_u

  !===================================================

  subroutine sol_ex003_3d_p(point,val)
    implicit none
    type(point_t), intent(in)    :: point
    real(rp),      intent(inout) :: val
    val = 0.0
  end subroutine sol_ex003_3d_p

  subroutine sol_ex003_3d_grad_p(point,val)
    implicit none
    type(point_t), intent(in)    :: point
    type(vector_field_t), intent(inout) :: val
    call val%init(0.0)
  end subroutine sol_ex003_3d_grad_p

  !===================================================

  subroutine base_scalar_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_scalar_function_set_num_dims

  subroutine base_scalar_function_set_case_id(this,val)
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip),                   intent(in)       :: val
    this%case_id = val
  end subroutine base_scalar_function_set_case_id

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

  subroutine base_vector_function_set_case_id(this,val)
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip),                   intent(in)       :: val
    this%case_id = val
  end subroutine base_vector_function_set_case_id

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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_lapl_u(point,val,this%degree)
          result = (-1.0)*val
          call sol_ex001_2d_grad_p(point,val,this%degree)
          result = result + (-1.0)*val
        case (2)
          call sol_ex002_2d_lapl_u(point,val,this%degree)
          result = (-1.0)*val
          call sol_ex002_2d_grad_p(point,val,this%degree)
          result = result + (-1.0)*val
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_lapl_u(point,val,this%degree)
          result = (-1.0)*val
          call sol_ex001_3d_grad_p(point,val,this%degree)
          result = result + (-1.0)*val
        case (2)
          call sol_ex002_3d_lapl_u(point,val,this%degree)
          result = (-1.0)*val
          call sol_ex002_3d_grad_p(point,val,this%degree)
          result = result + (-1.0)*val
        case (3)
          call sol_ex003_3d_lapl_u(point,val)
          result = (-1.0)*val
          call sol_ex003_3d_grad_p(point,val)
          result = result + (-1.0)*val
        case default
          check(.false.)
      end select
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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_u(point,result,this%degree)
        case (2)
          call sol_ex002_2d_u(point,result,this%degree)
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_u(point,result,this%degree)
        case (2)
          call sol_ex002_3d_u(point,result,this%degree)
        case (3)
          call sol_ex003_3d_u(point,result)
        case default
          check(.false.)
      end select
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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_grad_u(point,result,this%degree)
        case (2)
          call sol_ex002_2d_grad_u(point,result,this%degree)
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_grad_u(point,result,this%degree)
        case (2)
          call sol_ex002_3d_grad_u(point,result,this%degree)
        case (3)
          call sol_ex003_3d_grad_u(point,result)
        case default
          check(.false.)
      end select
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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_div_u(point,result,this%degree)
        case (2)
          call sol_ex002_2d_div_u(point,result,this%degree)
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_div_u(point,result,this%degree)
        case (2)
          call sol_ex002_3d_div_u(point,result,this%degree)
        case (3)
          call sol_ex003_3d_div_u(point,result)
        case default
          check(.false.)
      end select
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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_p(point,result,this%degree)
        case (2)
          call sol_ex002_2d_p(point,result,this%degree)
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_p(point,result,this%degree)
        case (2)
          call sol_ex002_3d_p(point,result,this%degree)
        case (3)
          call sol_ex003_3d_p(point,result)
        case default
          check(.false.)
      end select
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
      select case (this%case_id)
        case (1)
          call sol_ex001_2d_grad_p(point,result,this%degree)
        case (2)
          call sol_ex002_2d_grad_p(point,result,this%degree)
        case default
          check(.false.)
      end select
    else if ( this%num_dims == 3 ) then
      select case (this%case_id)
        case (1)
          call sol_ex001_3d_grad_p(point,result,this%degree)
        case (2)
          call sol_ex002_3d_grad_p(point,result,this%degree)
        case (3)
          call sol_ex003_3d_grad_p(point,result)
        case default
          check(.false.)
      end select
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

  subroutine stokes_analytical_functions_set_case_id(this,val)
    implicit none
    class(stokes_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) :: val
    call this%source_term_u%set_case_id(val)
    call this%source_term_p%set_case_id(val)
    call this%solution_function_u%set_case_id(val)
    call this%solution_function_p%set_case_id(val)
  end subroutine stokes_analytical_functions_set_case_id

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
