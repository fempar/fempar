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

module hts_nedelec_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(vector_function_t) :: base_vector_function_t
     integer(ip) :: num_dims = -1  
   contains
     procedure :: set_num_dims    => base_vector_function_set_num_dims
  end type base_vector_function_t

  type, extends(base_vector_function_t) :: source_term_t
   contains
     procedure :: get_value_space      => source_term_get_value_space
     procedure :: get_value_space_time => source_term_get_value_space_time 
  end type source_term_t

  type, extends(base_vector_function_t) :: solution_t
   contains
     procedure :: get_value_space         => solution_get_value_space 
     procedure :: get_value_space_time    => solution_get_value_space_time 
     procedure :: get_gradient_space      => solution_get_gradient_space
     procedure :: get_gradient_space_time => solution_get_gradient_space_time 
     procedure :: get_curl_space          => solution_get_curl_space
     procedure :: get_curl_space_time     => solution_get_curl_space_time 
  end type solution_t

  type, extends(scalar_function_t) :: base_scalar_function_t
     integer(ip) :: num_dims = -1  
   contains
  end type base_scalar_function_t

  type, extends(base_scalar_function_t) :: boundary_function_Hx_t
   private 
     real(rp) :: amplitude
     real(rp) :: frequency 
   contains
     procedure :: get_value_space        => boundary_function_Hx_get_value_space
     procedure :: get_value_space_time   => boundary_function_Hx_get_value_space_time
  end type boundary_function_Hx_t

  type, extends(base_scalar_function_t) :: boundary_function_Hy_t
    private
     real(rp) :: amplitude
     real(rp) :: frequency  
   contains
     procedure :: get_value_space        => boundary_function_Hy_get_value_space
     procedure :: get_value_space_time   => boundary_function_Hy_get_value_space_time
  end type boundary_function_Hy_t

  type, extends(base_scalar_function_t) :: boundary_function_Hz_t
    private 
     real(rp) :: amplitude 
     real(rp) :: frequency 
   contains
     procedure :: get_value_space        => boundary_function_Hz_get_value_space
     procedure :: get_value_space_time   => boundary_function_Hz_get_value_space_time
  end type boundary_function_Hz_t
  
    type, extends(base_scalar_function_t) :: boundary_function_p_t
     private 
   contains
     procedure :: get_value_space        => boundary_function_p_get_value_space
     procedure :: get_value_space_time   => boundary_function_p_get_value_space_time
  end type boundary_function_p_t

   type, extends(base_scalar_function_t) :: constraint_value_t 
    private 
     real(rp)  :: amplitude 
     real(rp)  :: frequency 
   contains
     procedure :: get_constraint_value  => constraint_value_get_constraint_value
  end type constraint_value_t

  type hts_nedelec_analytical_functions_t
     private
     type(source_term_t)                     :: source_term
     type(solution_t)                        :: solution
     type(boundary_function_Hx_t)            :: boundary_function_Hx
     type(boundary_function_Hy_t)            :: boundary_function_Hy
     type(boundary_function_Hz_t)            :: boundary_function_Hz 
     type(boundary_function_p_t)             :: boundary_function_p
     type(constraint_value_t)                :: constraint_value 
   contains
     procedure :: set_num_dims               => mn_set_num_dims
     procedure :: initialize                       => mn_initialize 
     procedure :: get_source_term                  => mn_get_source_term
     procedure :: get_solution                     => mn_get_solution
     procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
     procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
     procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
     procedure :: get_boundary_function_p          => mn_get_boundary_function_p
     procedure :: get_constraint_value             => mn_get_constraint_value 
  end type hts_nedelec_analytical_functions_t

  public :: hts_nedelec_analytical_functions_t, constraint_value_t

contains  

  !===============================================================================================
  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims
  
  !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z, n  
    x = point%get(1)
    y = point%get(2) 
    n = 3.0_rp 
    
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    !call result%set(1, 0.0_rp) 
    !call result%set(2, -(n+1.0_rp)*((3.0_rp*x*x)**n)*6.0_rp*x + x*x*x )
    !call result%set(3, 0.0_rp) 
  end subroutine source_term_get_value_space
  
    !===============================================================================================
  subroutine source_term_get_value_space_time ( this, point, time, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(in)    :: time 
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z, n  
    x = point%get(1)
    y = point%get(2) 
    n = 3.0_rp 
    
    call result%init(0.0_rp) 
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    !call result%set(1, 0.0_rp ) 
    !call result%set(2, -(n+1.0_rp)*((time*3.0_rp*x*x)**n)*(time*6.0_rp*x) + x*x*x )
    !call result%set(3, 0.0_rp) 

  end subroutine source_term_get_value_space_time 

  !===============================================================================================
  subroutine source_term_get_gradient_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t)    , intent(inout) :: result
    check(.false.)
  end subroutine source_term_get_gradient_space

  !===============================================================================================
  subroutine solution_get_value_space ( this, point, result )
    implicit none
    class(solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z
    x = point%get(1)
    y = point%get(2) 
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    call result%set(1, 0.0_rp ) 
    call result%set(2,  x*x*x  ) 
    call result%set(3, 0.0_rp) 
 
  end subroutine solution_get_value_space
  
    !===============================================================================================
  subroutine solution_get_value_space_time ( this, point, time, result )
    implicit none
    class(solution_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(in)    :: time 
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z
    x = point%get(1)
    y = point%get(2) 
    assert ( this%num_dims == 2 .or. this%num_dims == 3 )
    call result%set(1,  0.0_rp ) 
    call result%set(2,  time*x*x*x ) 
    call result%set(3, 0.0_rp) 
 
  end subroutine solution_get_value_space_time
    
  !===============================================================================================
  subroutine solution_get_gradient_space ( this, point, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(tensor_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z  
    x = point%get(1)
    y = point%get(2) 
    z = point%get(3) 
    
    call result%set(1, 2, 3.0_rp*x*x )

  end subroutine solution_get_gradient_space
  
   !===============================================================================================
  subroutine solution_get_gradient_space_time ( this, point, time, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(in)    :: time
    type(tensor_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z  
    x = point%get(1)
    y = point%get(2) 
    z = point%get(3) 

    call result%set(1,2, time*3.0_rp*x*x) 

  end subroutine solution_get_gradient_space_time

  !===============================================================================================
  subroutine solution_get_curl_space ( this, point, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z  
    x = point%get(1)
    y = point%get(2) 

    !call result%set(1, 0.0_rp) 
    !call result%set(2, 0.0_rp) 
    !call result%set(3, 3.0_rp*x*x ) 
    
  end subroutine solution_get_curl_space
  
   subroutine solution_get_curl_space_time ( this, point, time, result )
    implicit none
    class(solution_t)       , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(in)    :: time 
    type(vector_field_t)    , intent(inout) :: result
    ! Locals 
    real(rp)  :: x, y, z  
    x = point%get(1)
    y = point%get(2) 

    !call result%set(1, 0.0_rp) 
    !call result%set(2, 0.0_rp) 
    !call result%set(3, time*3.0_rp*x*x ) 
    
  end subroutine solution_get_curl_space_time
    
  !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 

    result = 0.0_rp
    
  end subroutine boundary_function_Hx_get_value_space
  
    !===============================================================================================
  subroutine boundary_function_Hx_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time
    real(rp)                       , intent(inout) :: result 

    result = this%amplitude*sin(2.0_rp*this%frequency*pi*time)

  end subroutine boundary_function_Hx_get_value_space_time 

  !===============================================================================================
  subroutine boundary_function_Hy_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result    

     result = 0.0_rp  

  end subroutine boundary_function_Hy_get_value_space
  
    !===============================================================================================
  subroutine boundary_function_Hy_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time
    real(rp)                       , intent(inout) :: result  

    result = this%amplitude*sin(2.0_rp*this%frequency*pi*time) 
    
  end subroutine boundary_function_Hy_get_value_space_time

    ! ============================================================================================
  subroutine boundary_function_Hz_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result    

    result = 0.0_rp
    
  end subroutine boundary_function_Hz_get_value_space
  
  ! ============================================================================================
  subroutine boundary_function_Hz_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result    

    result = this%amplitude*sin(2.0_rp*this%frequency*pi*time) 
 
	! Benchmark ramp 
	if ( time .le. 5 ) then 
	result = this%amplitude*(time/5.0_rp) 
	elseif ( time .le. 10 ) then 
	result = this%amplitude 
	elseif ( time .le. 15 ) then 
	result = this%amplitude*(15.0_rp-time)/(5.0_rp)
	else 
	result = 0.0_rp
	end if 
    
  end subroutine boundary_function_Hz_get_value_space_time
  
      ! ============================================================================================
  subroutine boundary_function_p_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_p_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result    

    result = 0.0_rp    
  end subroutine boundary_function_p_get_value_space
  
    ! ============================================================================================
  subroutine boundary_function_p_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_p_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result    

    result = 0.0_rp    
  end subroutine boundary_function_p_get_value_space_time
  
      ! ============================================================================================
  subroutine constraint_value_get_constraint_value( this, time, result )
    implicit none 
    class(constraint_value_t)      , intent(in)    :: this 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result    

    result = this%amplitude*sin(2.0_rp*pi*this%frequency*time) 

  end subroutine constraint_value_get_constraint_value
  
  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(hts_nedelec_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%solution%set_num_dims(num_dims)
  end subroutine mn_set_num_dims
  
    !===============================================================================================
  subroutine mn_initialize ( this, H, wH, J, wJ )
    implicit none
    class(hts_nedelec_analytical_functions_t), intent(inout)    :: this
    real(rp)      , intent(in)   :: H(0:SPACE_DIM-1), J(0:SPACE_DIM-1) 
    real(rp)      , intent(in)   :: wH, wJ 
    
    this%boundary_function_Hx%amplitude = H(0) 
    this%boundary_function_Hx%frequency = wH
    this%boundary_function_Hy%amplitude = H(1) 
    this%boundary_function_Hy%frequency = wH
    this%boundary_function_Hz%amplitude = H(2) 
    this%boundary_function_Hz%frequency = wH
   
    this%constraint_value%amplitude = J(2) 
    this%constraint_value%frequency = wJ 

  end subroutine mn_initialize 

  !===============================================================================================
  function mn_get_solution ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_solution
    mn_get_solution => this%solution
  end function mn_get_solution

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term

  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz
  end function mn_get_boundary_function_Hz
  
    !===============================================================================================
  function mn_get_boundary_function_p ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_p
    mn_get_boundary_function_p => this%boundary_function_p
  end function mn_get_boundary_function_p
  
     !===============================================================================================
  function mn_get_constraint_value ( this )
    implicit none
    class(hts_nedelec_analytical_functions_t), target, intent(in)    :: this
    class(constraint_value_t), pointer :: mn_get_constraint_value
    mn_get_constraint_value=> this%constraint_value
  end function mn_get_constraint_value

end module hts_nedelec_analytical_functions_names
!***************************************************************************************************


