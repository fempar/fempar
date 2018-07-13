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

module par_pb_bddc_maxwell_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  ! Scalar functions 
  type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims = -1  
  contains
  end type base_scalar_function_t
    
  type, extends(base_scalar_function_t) :: resistivity_function_t
   private 
    real(rp) :: default_value 
   contains
     procedure :: set_value          => resistivity_function_set_default_value
  end type resistivity_function_t 
  
  type :: resistivity_holder_t 
     class(resistivity_function_t), pointer :: p => NULL()  
  end type resistivity_holder_t 
  
  type, extends(resistivity_function_t) :: resistivity_function_white_t
    private 
   contains
     procedure :: get_value_space    => resistivity_function_white_get_value_space
  end type resistivity_function_white_t 
  
  type, extends(resistivity_function_t) :: resistivity_function_black_t
    private 
   contains
     procedure :: get_value_space    => resistivity_function_black_get_value_space
  end type resistivity_function_black_t 
  
  type, extends(base_scalar_function_t) :: permeability_function_t
   private 
    real(rp) :: default_value 
   contains
     procedure :: set_value          => permeability_function_set_default_value
  end type permeability_function_t 
  
  type :: permeability_holder_t 
     class(permeability_function_t), pointer :: p => NULL()  
  end type permeability_holder_t 
  
  type, extends(permeability_function_t) :: permeability_function_white_t
    private 
   contains
     procedure :: get_value_space    => permeability_function_white_get_value_space
  end type permeability_function_white_t 
  
  type, extends(permeability_function_t) :: permeability_function_black_t
    private 
   contains
     procedure :: get_value_space    => permeability_function_black_get_value_space
  end type permeability_function_black_t 
  
  ! Boundary scalar functions definition 
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
  
  ! Vector functions 
   type, extends(vector_function_t) :: base_vector_function_t
    integer(ip) :: num_dims = -1  
  contains
    procedure :: set_num_dims    => base_vector_function_set_num_dims
  end type base_vector_function_t
  
    type, extends(base_vector_function_t) :: source_term_t
     real(rp) :: permeability 
     real(rp) :: resistivity 
   contains
     procedure :: get_value_space => source_term_get_value_space
  end type source_term_t

   type, extends(base_vector_function_t) :: solution_t
   contains
     procedure :: get_value_space    => solution_get_value_space
     procedure :: get_gradient_space => solution_get_gradient_space
  end type solution_t
  
  type par_pb_bddc_maxwell_analytical_functions_t
     private
   type(resistivity_holder_t), pointer     :: resistivity(:)
   type(permeability_holder_t), pointer    :: permeability(:) 
	  type(boundary_function_Hx_t)            :: boundary_function_Hx
	  type(boundary_function_Hy_t)            :: boundary_function_Hy
	  type(boundary_function_Hz_t)            :: boundary_function_Hz
   type(source_term_t)                     :: source_term
   type(solution_t)                        :: solution
   contains
   procedure :: set_num_dims                     => mn_set_num_dims
   procedure :: set_resistivity                  => mn_set_resistivity_holder
   procedure :: set_permeability                 => mn_set_permeability_holder 
   ! Getters 
   procedure :: get_resistivity                  => mn_get_resistivity
   procedure :: get_permeability                 => mn_get_permeability 
	  procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
	  procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
	  procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
   procedure :: get_source_term                  => mn_get_source_term
   procedure :: get_solution_function            => mn_get_solution_function
  end type par_pb_bddc_maxwell_analytical_functions_t

  public :: par_pb_bddc_maxwell_analytical_functions_t
  public :: resistivity_holder_t, resistivity_function_white_t, resistivity_function_black_t 
  public :: permeability_holder_t, permeability_function_white_t, permeability_function_black_t 

contains  
  !===============================================================================================
  subroutine resistivity_function_set_default_value ( this, default_value )
    implicit none
    class(resistivity_function_t), intent(inout)    :: this
    real(rp), intent(in) ::  default_value
    this%default_value = default_value
  end subroutine resistivity_function_set_default_value
  
    !===============================================================================================
  subroutine permeability_function_set_default_value ( this, default_value )
    implicit none
    class(permeability_function_t), intent(inout)    :: this
    real(rp), intent(in) ::  default_value
    this%default_value = default_value
  end subroutine permeability_function_set_default_value
  
  !===============================================================================================
  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims
  
        !===============================================================================================
  subroutine resistivity_function_white_get_value_space( this, point, result )
    implicit none 
    class(resistivity_function_white_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
	   result = this%default_value 
    ! result = exp( 4.0_rp*sin(5*pi*(point%get(1)+point%get(2)+point%get(3))))
  end subroutine resistivity_function_white_get_value_space
  
    !===============================================================================================
  subroutine resistivity_function_black_get_value_space( this, point, result )
    implicit none 
    class(resistivity_function_black_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
    result = this%default_value 
    ! result = 10**(point%get(1) + point%get(2) + point%get(3) ) 
  end subroutine resistivity_function_black_get_value_space
  
          !===============================================================================================
  subroutine permeability_function_white_get_value_space( this, point, result )
    implicit none 
    class(permeability_function_white_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
	   result = this%default_value  
    !result = 10**(point%get(1) + point%get(2) + point%get(3) )
  end subroutine permeability_function_white_get_value_space
  
    !===============================================================================================
  subroutine permeability_function_black_get_value_space( this, point, result )
    implicit none 
    class(permeability_function_black_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
     result = this%default_value  
    ! result = 10**(-(point%get(1)+point%get(2)+point%get(3)))
  end subroutine permeability_function_black_get_value_space
  
      !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
		real(rp) :: x,y,z 
	x = point%get(1); y=point%get(2); z=point%get(3)
	 result = -y
  end subroutine boundary_function_Hx_get_value_space
  
  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
		  real(rp) :: x,y,z 
	   x = point%get(1); y=point%get(2); z=point%get(3)
     result = x
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
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
 
	real(rp) :: x,y,z 

	x = point%get(1); y=point%get(2); z=point%get(3)     
	 ! call result%init(1.0_rp) 
	 call result%set(1, -y ) 
	 call result%set(2,  x ) 
	 call result%set(3,  0.0_rp )
  
  !call result%init(1.0_rp)
  end subroutine source_term_get_value_space

  !===============================================================================================
  subroutine solution_get_value_space ( this, point, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
	
    real(rp) :: x,y,z 
	   x = point%get(1); y=point%get(2); z=point%get(3) 
	 ! call result%init(0.0_rp) 
	 call result%set(1, -y ) 
	 call result%set(2,  x ) 
	 call result%set(3,  0.0_rp )
	 
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
	call result%set(2,1, -1.0_rp)
	call result%set(1,2,  1.0_rp)

  end subroutine solution_get_gradient_space
  
  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%solution%set_num_dims(num_dims)
  end subroutine mn_set_num_dims
  
    !===============================================================================================
  subroutine mn_set_resistivity_holder ( this, resistivity_holder )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    type(resistivity_holder_t), target   , intent(in)       :: resistivity_holder(:)
    this%resistivity => resistivity_holder 
  end subroutine mn_set_resistivity_holder
  
      !===============================================================================================
  subroutine mn_set_permeability_holder ( this, permeability_holder )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), intent(inout)    :: this
    type(permeability_holder_t), target   , intent(in)      :: permeability_holder(:)
    this%permeability => permeability_holder 
  end subroutine mn_set_permeability_holder
  
      !===============================================================================================
  function mn_get_resistivity ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    type(resistivity_holder_t), pointer     :: mn_get_resistivity(:)
    mn_get_resistivity => this%resistivity 
  end function mn_get_resistivity
  
        !===============================================================================================
  function mn_get_permeability ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    type(permeability_holder_t), pointer     :: mn_get_permeability(:)
    mn_get_permeability => this%permeability 
  end function mn_get_permeability
  
  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz 
  end function mn_get_boundary_function_Hz
  
  !===============================================================================================
  function mn_get_solution_function ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_solution_function
    mn_get_solution_function => this%solution
  end function mn_get_solution_function

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(par_pb_bddc_maxwell_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term

end module par_pb_bddc_maxwell_analytical_functions_names
!***************************************************************************************************
