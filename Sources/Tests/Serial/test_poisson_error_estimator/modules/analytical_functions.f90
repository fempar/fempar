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

module analytical_functions_names
  use fempar_names
  use test_poisson_error_estimator_params_names
  implicit none
# include "debug.i90"
  private

  real(rp), parameter :: alpha = 60.0_rp
  real(rp), parameter :: x0(3) = (/ 1.25_rp, -0.25_rp, 0.00_rp /)

  type, extends(scalar_function_t) :: base_scalar_function_t
    private
    integer(ip) :: order    =  5
  contains
    procedure :: set_num_dims          => base_scalar_function_set_num_dims
    procedure :: set_polynomial_order  => base_scalar_function_set_polynomial_order
  end type base_scalar_function_t
  
  type, extends(base_scalar_function_t) :: polynomial_source_term_t
    private 
   contains
     procedure :: get_value_space    => polynomial_source_term_get_value_space
  end type polynomial_source_term_t

  type, extends(base_scalar_function_t) :: polynomial_boundary_function_t
    private
   contains
     procedure :: get_value_space => polynomial_boundary_function_get_value_space
  end type polynomial_boundary_function_t

  type, extends(base_scalar_function_t) :: polynomial_solution_function_t
    private 
   contains
     procedure :: get_value_space    => polynomial_solution_function_get_value_space
     procedure :: get_gradient_space => polynomial_solution_function_get_gradient_space
  end type polynomial_solution_function_t

  type, extends(base_scalar_function_t) :: shock_source_term_t
    private 
   contains
     procedure :: get_value_space    => shock_source_term_get_value_space
  end type shock_source_term_t

  type, extends(base_scalar_function_t) :: shock_boundary_function_t
    private
   contains
     procedure :: get_value_space    => shock_boundary_function_get_value_space
  end type shock_boundary_function_t

  type, extends(base_scalar_function_t) :: shock_solution_function_t
    private 
   contains
     procedure :: get_value_space    => shock_solution_function_get_value_space
     procedure :: get_gradient_space => shock_solution_function_get_gradient_space
  end type shock_solution_function_t

  type, abstract :: base_analytical_functions_t
     private
   contains
     procedure ( set_num_dims_interface )         , deferred :: set_num_dims
     procedure ( get_source_term_interface )      , deferred :: get_source_term
     procedure ( get_boundary_function_interface ), deferred :: get_boundary_function
     procedure ( get_solution_function_interface ), deferred :: get_solution_function
  end type base_analytical_functions_t

  abstract interface
  
    subroutine set_num_dims_interface(this,num_dims)
      import :: base_analytical_functions_t, ip
      class(base_analytical_functions_t), intent(inout) :: this
      integer(ip)                       , intent(in)    :: num_dims
    end subroutine set_num_dims_interface
  
    function get_source_term_interface(this)
      import :: base_analytical_functions_t, scalar_function_t
      class(base_analytical_functions_t), target, intent(in) :: this
      class(scalar_function_t), pointer :: get_source_term_interface
    end function get_source_term_interface
  
    function get_boundary_function_interface(this)
      import :: base_analytical_functions_t, scalar_function_t
      class(base_analytical_functions_t), target, intent(in) :: this
      class(scalar_function_t), pointer :: get_boundary_function_interface
    end function get_boundary_function_interface
  
    function get_solution_function_interface(this)
      import :: base_analytical_functions_t, scalar_function_t
      class(base_analytical_functions_t), target, intent(in) :: this
      class(scalar_function_t), pointer :: get_solution_function_interface
    end function get_solution_function_interface
  
  end interface

  type, extends(base_analytical_functions_t) :: polynomial_analytical_functions_t
     private
     type(polynomial_source_term_t)       :: source_term
     type(polynomial_boundary_function_t) :: boundary_function
     type(polynomial_solution_function_t) :: solution_function
   contains
     procedure :: set_num_dims          => polynomial_analytical_functions_set_num_dims
     procedure :: get_source_term       => polynomial_analytical_functions_get_source_term
     procedure :: get_boundary_function => polynomial_analytical_functions_get_boundary_function
     procedure :: get_solution_function => polynomial_analytical_functions_get_solution_function
  end type polynomial_analytical_functions_t

  type, extends(base_analytical_functions_t) :: shock_analytical_functions_t
     private
     type(shock_source_term_t)       :: source_term
     type(shock_boundary_function_t) :: boundary_function
     type(shock_solution_function_t) :: solution_function
   contains
     procedure :: set_num_dims          => shock_analytical_functions_set_num_dims
     procedure :: get_source_term       => shock_analytical_functions_get_source_term
     procedure :: get_boundary_function => shock_analytical_functions_get_boundary_function
     procedure :: get_solution_function => shock_analytical_functions_get_solution_function
  end type shock_analytical_functions_t
  
  type analytical_functions_t
     private
     class(base_analytical_functions_t), allocatable :: analytical_functions
   contains
     procedure :: create                   => analytical_functions_create
     procedure :: free                     => analytical_functions_free
     procedure :: set_num_dims             => analytical_functions_set_num_dims
     procedure :: get_source_term          => analytical_functions_get_source_term
     procedure :: get_boundary_function    => analytical_functions_get_boundary_function
     procedure :: get_solution_function    => analytical_functions_get_solution_function
     procedure :: get_analytical_functions => analytical_functions_get_analytical_functions
  end type analytical_functions_t
  
  public :: base_analytical_functions_t, analytical_functions_t

contains  

  !===============================================================================================
  subroutine base_scalar_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%scalar_function_t%set_num_dims(num_dims)
  end subroutine base_scalar_function_set_num_dims
  
  !===============================================================================================
  subroutine base_scalar_function_set_polynomial_order ( this, order )
    implicit none
    class(base_scalar_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  order
    this%order = order
  end subroutine base_scalar_function_set_polynomial_order
  
  !===============================================================================================
  subroutine polynomial_source_term_get_value_space ( this, point, result )
    implicit none
    class(polynomial_source_term_t), intent(in)    :: this
    type(point_t)                  , intent(in)    :: point
    real(rp)                       , intent(inout) :: result
    real(rp) :: n
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    n = real(this%order, rp)
    if ( this%get_num_dims() == 2 ) then
      result =  -n*(n-1.0_rp)*(point%get(1)**(n-2.0_rp) + point%get(2)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2})
    else if ( this%get_num_dims() == 3 ) then
      result =  -n*(n-1.0_rp)*(point%get(1)**(n-2.0_rp) + & 
                point%get(2)**(n-2.0_rp) + point%get(3)**(n-2.0_rp)) ! -n(n-1)(x^{n-2}+y^{n-2}+z^{n-2})
    end if 
  end subroutine polynomial_source_term_get_value_space

  !===============================================================================================
  subroutine polynomial_boundary_function_get_value_space ( this, point, result )
    implicit none
    class(polynomial_boundary_function_t), intent(in)  :: this
    type(point_t)                        , intent(in)    :: point
    real(rp)                             , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1)**this%order + point%get(2)**this%order ! x^n+y^n
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1)**this%order + point%get(2)**this%order + point%get(3)**this%order ! x^n+y^n+z^n
    end if  
  end subroutine polynomial_boundary_function_get_value_space 

  !===============================================================================================
  subroutine polynomial_solution_function_get_value_space ( this, point, result )
    implicit none
    class(polynomial_solution_function_t), intent(in)    :: this
    type(point_t)                        , intent(in)    :: point
    real(rp)                             , intent(inout) :: result
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    if ( this%get_num_dims() == 2 ) then
      result = point%get(1)**this%order + point%get(2)**this%order ! x^n+y^n 
    else if ( this%get_num_dims() == 3 ) then
      result = point%get(1)**this%order + point%get(2)**this%order + point%get(3)**this%order ! x^n+y^n+z^n 
    end if  
      
  end subroutine polynomial_solution_function_get_value_space
  
  !===============================================================================================
  subroutine polynomial_solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(polynomial_solution_function_t), intent(in)    :: this
    type(point_t)                        , intent(in)    :: point
    type(vector_field_t)                 , intent(inout) :: result
    real(rp) :: n 
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    n = real(this%order, rp) 
    if ( this%get_num_dims() == 2 ) then
      call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
    else if ( this%get_num_dims() == 3 ) then
      call result%set( 1, n*point%get(1)**(n-1.0_rp) ) ! nx^{n-1}
      call result%set( 2, n*point%get(2)**(n-1.0_rp) ) ! ny^{n-1}
      call result%set( 3, n*point%get(3)**(n-1.0_rp) ) ! nz^{n-1}
    end if
  end subroutine polynomial_solution_function_get_gradient_space
  
  !===============================================================================================
  subroutine polynomial_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(polynomial_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%boundary_function%set_num_dims(num_dims)
    call this%solution_function%set_num_dims(num_dims)
  end subroutine polynomial_analytical_functions_set_num_dims 
  
  !===============================================================================================
  function polynomial_analytical_functions_get_source_term ( this )
    implicit none
    class(polynomial_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: polynomial_analytical_functions_get_source_term
    polynomial_analytical_functions_get_source_term => this%source_term
  end function polynomial_analytical_functions_get_source_term
  
  !===============================================================================================
  function polynomial_analytical_functions_get_boundary_function ( this )
    implicit none
    class(polynomial_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: polynomial_analytical_functions_get_boundary_function
    polynomial_analytical_functions_get_boundary_function => this%boundary_function
  end function polynomial_analytical_functions_get_boundary_function
  
  !===============================================================================================
  function polynomial_analytical_functions_get_solution_function ( this )
    implicit none
    class(polynomial_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: polynomial_analytical_functions_get_solution_function
    polynomial_analytical_functions_get_solution_function => this%solution_function
  end function polynomial_analytical_functions_get_solution_function
  
  !===============================================================================================
  subroutine shock_source_term_get_value_space ( this, point, result )
    implicit none
    class(shock_source_term_t), intent(in)    :: this
    type(point_t)       , intent(in)    :: point
    real(rp)            , intent(inout) :: result
    type(vector_field_t) :: r_vec
    real(rp)             :: r_nrm
    r_vec = point - vector_field_t(x0)
    r_nrm = r_vec%nrm2()
    result = 2.0_rp*(alpha**3.0_rp)*(r_nrm-1)/((1+(alpha**2.0_rp)*(r_nrm-1)**2.0_rp)**2.0_rp) - &
             (this%get_num_dims()-1)*alpha/(r_nrm*(1+(alpha**2.0_rp)*(r_nrm-1)**2.0_rp))
  end subroutine shock_source_term_get_value_space

  !===============================================================================================
  subroutine shock_boundary_function_get_value_space ( this, point, result )
    implicit none
    class(shock_boundary_function_t), intent(in)  :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    type(vector_field_t) :: r_vec
    real(rp)             :: r_nrm
    r_vec = point - vector_field_t(x0)
    r_nrm = r_vec%nrm2()
    result = atan(alpha*(r_nrm-1.0_rp))
  end subroutine shock_boundary_function_get_value_space 

  !===============================================================================================
  subroutine shock_solution_function_get_value_space ( this, point, result )
    implicit none
    class(shock_solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    real(rp)                  , intent(inout) :: result
    type(vector_field_t) :: r_vec
    real(rp)             :: r_nrm
    r_vec = point - vector_field_t(x0)
    r_nrm = r_vec%nrm2()
    result = atan(alpha*(r_nrm-1.0_rp))
  end subroutine shock_solution_function_get_value_space
  
  !===============================================================================================
  subroutine shock_solution_function_get_gradient_space ( this, point, result )
    implicit none
    class(shock_solution_function_t), intent(in)    :: this
    type(point_t)             , intent(in)    :: point
    type(vector_field_t)      , intent(inout) :: result
    type(vector_field_t) :: r_vec
    real(rp)             :: r_nrm
    integer(ip)          :: idime
    r_vec = point - vector_field_t(x0)
    r_nrm = r_vec%nrm2()
    do idime = 1,this%get_num_dims()
      call result%set(idime,(alpha*r_vec%get(idime))/(r_nrm*(1+(alpha**2.0_rp)*(r_nrm-1)**2.0_rp)))
    end do
  end subroutine shock_solution_function_get_gradient_space
  
  !===============================================================================================
  subroutine shock_analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(shock_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%boundary_function%set_num_dims(num_dims)
    call this%solution_function%set_num_dims(num_dims)
  end subroutine shock_analytical_functions_set_num_dims 
  
  !===============================================================================================
  function shock_analytical_functions_get_source_term ( this )
    implicit none
    class(shock_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: shock_analytical_functions_get_source_term
    shock_analytical_functions_get_source_term => this%source_term
  end function shock_analytical_functions_get_source_term
  
  !===============================================================================================
  function shock_analytical_functions_get_boundary_function ( this )
    implicit none
    class(shock_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: shock_analytical_functions_get_boundary_function
    shock_analytical_functions_get_boundary_function => this%boundary_function
  end function shock_analytical_functions_get_boundary_function
  
  !===============================================================================================
  function shock_analytical_functions_get_solution_function ( this )
    implicit none
    class(shock_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: shock_analytical_functions_get_solution_function
    shock_analytical_functions_get_solution_function => this%solution_function
  end function shock_analytical_functions_get_solution_function
  
  !===============================================================================================
  subroutine analytical_functions_set_num_dims ( this, num_dims )
    implicit none
    class(analytical_functions_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: num_dims
    call this%analytical_functions%set_num_dims(num_dims)
  end subroutine analytical_functions_set_num_dims 
  
  !===============================================================================================
  function analytical_functions_get_source_term ( this )
    implicit none
    class(analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: analytical_functions_get_source_term
    analytical_functions_get_source_term => this%analytical_functions%get_source_term()
  end function analytical_functions_get_source_term
  
  !===============================================================================================
  function analytical_functions_get_boundary_function ( this )
    implicit none
    class(analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: analytical_functions_get_boundary_function
    analytical_functions_get_boundary_function => this%analytical_functions%get_boundary_function()
  end function analytical_functions_get_boundary_function
  
  !===============================================================================================
  function analytical_functions_get_solution_function ( this )
    implicit none
    class(analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: analytical_functions_get_solution_function
    analytical_functions_get_solution_function => this%analytical_functions%get_solution_function()
  end function analytical_functions_get_solution_function
  
  !===============================================================================================
  function analytical_functions_get_analytical_functions ( this )
    implicit none
    class(analytical_functions_t), target, intent(inout) :: this
    class(base_analytical_functions_t), pointer :: analytical_functions_get_analytical_functions
    analytical_functions_get_analytical_functions => this%analytical_functions
  end function analytical_functions_get_analytical_functions
  
  !===============================================================================================
  subroutine analytical_functions_create ( this, parameters )
    implicit none
    class(analytical_functions_t)              , intent(inout) :: this
    type(test_poisson_error_estimator_params_t), intent(in)    :: parameters
    integer(ip) :: istat
    call this%free()
    if ( parameters%get_analytical_functions_type() == 'polynomial' ) then
      allocate ( polynomial_analytical_functions_t :: this%analytical_functions, stat = istat )
      check(istat == 0)
    else if ( parameters%get_analytical_functions_type() == 'shock' ) then
      allocate ( shock_analytical_functions_t :: this%analytical_functions, stat = istat )
      check(istat == 0)
    else
      mcheck(.false.,'this family of analytical functions is not implemented')
    end if
  end subroutine analytical_functions_create
  
  !===============================================================================================
  subroutine analytical_functions_free ( this )
    implicit none
    class(analytical_functions_t), intent(inout) :: this
    integer(ip) :: istat
    if ( allocated(this%analytical_functions) ) then
      deallocate( this%analytical_functions, stat = istat ); check( istat == 0 )
    end if
  end subroutine analytical_functions_free

end module analytical_functions_names
!***************************************************************************************************
