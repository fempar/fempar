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

module hts_analytical_functions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  ! Scalar functions 
  type, extends(scalar_function_t) :: base_scalar_function_t
    integer(ip) :: num_dims = -1  
  contains
  end type base_scalar_function_t
  
    type, extends(base_scalar_function_t) :: boundary_function_Hx_t
    private 
    real(rp) :: amplitude
    real(rp) :: frequency  
   contains
     procedure :: get_value_space         => boundary_function_Hx_get_value_space
     procedure :: get_value_space_time    => boundary_function_Hx_get_value_space_time
  end type boundary_function_Hx_t 
  
    type, extends(base_scalar_function_t) :: boundary_function_Hy_t
    private 
    real(rp) :: amplitude
    real(rp) :: frequency 
   contains
     procedure :: get_value_space          => boundary_function_Hy_get_value_space
      procedure :: get_value_space_time    => boundary_function_Hy_get_value_space_time 
  end type boundary_function_Hy_t 
  
   type, extends(base_scalar_function_t) :: boundary_function_Hz_t
    private 
    real(rp) :: amplitude
    real(rp) :: frequency 
   contains
     procedure :: get_value_space         => boundary_function_Hz_get_value_space
     procedure :: get_value_space_time    => boundary_function_Hz_get_value_space_time
  end type boundary_function_Hz_t 
  
  ! Vector functions 
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
  
      type, extends(base_vector_function_t) :: initial_conditions_t
   contains
     procedure :: get_value_space_time => initial_conditions_get_value_space_time
  end type initial_conditions_t
 
  type hts_analytical_functions_t
     private
   type(boundary_function_Hx_t)            :: boundary_function_Hx
   type(boundary_function_Hy_t)            :: boundary_function_Hy
   type(boundary_function_Hz_t)            :: boundary_function_Hz
   type(source_term_t)                     :: source_term
   type(initial_conditions_t)              :: initial_conditions 
   contains
   procedure :: get_parameter_values             => mn_get_parameter_values 
   procedure :: set_num_dims                     => mn_set_num_dims
   procedure :: get_boundary_function_Hx         => mn_get_boundary_function_Hx
   procedure :: get_boundary_function_Hy         => mn_get_boundary_function_Hy
   procedure :: get_boundary_function_Hz         => mn_get_boundary_function_Hz
   procedure :: get_source_term                  => mn_get_source_term
   procedure :: get_initial_conditions           => mn_get_initial_conditions 
  end type hts_analytical_functions_t

  public :: hts_analytical_functions_t

contains  
  !===============================================================================================
  subroutine base_vector_function_set_num_dims ( this, num_dims )
    implicit none
    class(base_vector_function_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    this%num_dims = num_dims
  end subroutine base_vector_function_set_num_dims
  
    !===============================================================================================
  subroutine boundary_function_Hx_get_value_space( this, point, result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(inout) :: result 
  real(rp) :: x,y,z 
 x = point%get(1); y=point%get(2); z=point%get(3)
 result = 0.0_rp 

  end subroutine boundary_function_Hx_get_value_space
  
      !===============================================================================================
  subroutine boundary_function_Hx_get_value_space_time( this, point, time , result )
    implicit none 
    class(boundary_function_Hx_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result 
    
  real(rp) :: x,y,z 
  x = point%get(1); y=point%get(2); z=point%get(3)

 result = this%amplitude*sin(2.0_rp*pi*this%frequency*time) 

  end subroutine boundary_function_Hx_get_value_space_time

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
  subroutine boundary_function_Hy_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_Hy_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result 
    
  real(rp) :: x,y,z 
 x = point%get(1); y=point%get(2); z=point%get(3)
 result = this%amplitude*sin(2.0_rp*pi*this%frequency*time)

  end subroutine boundary_function_Hy_get_value_space_time

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
  subroutine boundary_function_Hz_get_value_space_time( this, point, time, result )
    implicit none 
    class(boundary_function_Hz_t)  , intent(in)    :: this 
    type(point_t)                  , intent(in)    :: point 
    real(rp)                       , intent(in)    :: time 
    real(rp)                       , intent(inout) :: result 
    
   real(rp) :: x,y,z 
   x = point%get(1); y=point%get(2); z=point%get(3)
   result = this%amplitude*sin(2.0_rp*pi*this%frequency*time)

  end subroutine boundary_function_Hz_get_value_space_time

 !===============================================================================================
  subroutine source_term_get_value_space ( this, point, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    type(vector_field_t)    , intent(inout) :: result
 
 real(rp) :: x,y,z 
 
 assert ( this%num_dims == 2 .or. this%num_dims == 3 )
 x = point%get(1); y=point%get(2); z=point%get(3)     
  call result%init(0.0_rp) 
  call result%set(1,  0.0_rp ) 
  call result%set(2,  0.0_rp ) 
  call result%set(3,  0.0_rp )

  end subroutine source_term_get_value_space
  
   !===============================================================================================
  subroutine source_term_get_value_space_time ( this, point, time, result )
    implicit none
    class(source_term_t)    , intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(in)    :: time 
    type(vector_field_t)    , intent(inout) :: result
 
 real(rp) :: x,y,z 
 
 assert ( this%num_dims == 2 .or. this%num_dims == 3 )
 x = point%get(1); y=point%get(2); z=point%get(3)     
  call result%init(0.0_rp) 
  call result%set(1,  0.0_rp ) 
  call result%set(2,  0.0_rp ) 
  call result%set(3,  0.0_rp )

  end subroutine source_term_get_value_space_time 
  
     !===============================================================================================
  subroutine initial_conditions_get_value_space_time ( this, point, time, result )
    implicit none
    class(initial_conditions_t)    , intent(in)    :: this
    type(point_t)                  , intent(in)    :: point
    real(rp)                       , intent(in)    :: time 
    type(vector_field_t)           , intent(inout) :: result
 
 real(rp) :: x,y,z 
 
 assert ( this%num_dims == 2 .or. this%num_dims == 3 )
 x = point%get(1); y=point%get(2); z=point%get(3)     
  call result%init(0.0_rp) 
  call result%set(1,  0.0_rp ) 
  call result%set(2,  0.0_rp ) 
  call result%set(3,  0.0_rp )

  end subroutine initial_conditions_get_value_space_time 
  
      !===============================================================================================
  subroutine mn_get_parameter_values ( this, H, wH  )
    implicit none
    class(hts_analytical_functions_t), intent(inout) :: this
    real(rp)                         , intent(in)    :: H(0:SPACE_DIM-1) 
    real(rp)                         , intent(in)    :: wH
    
    this%boundary_function_Hx%amplitude = H(0) 
    this%boundary_function_Hx%frequency = wH
    this%boundary_function_Hy%amplitude = H(1) 
    this%boundary_function_Hy%frequency = wH
    this%boundary_function_Hz%amplitude = H(2) 
    this%boundary_function_Hz%frequency = wH

  end subroutine mn_get_parameter_values 

  !===============================================================================================
  subroutine mn_set_num_dims ( this, num_dims )
    implicit none
    class(hts_analytical_functions_t), intent(inout)    :: this
    integer(ip), intent(in) ::  num_dims
    call this%source_term%set_num_dims(num_dims)
    call this%initial_conditions%set_num_dims(num_dims)
  end subroutine mn_set_num_dims
  
  !===============================================================================================
  function mn_get_boundary_function_Hx ( this )
    implicit none
    class(hts_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hx
    mn_get_boundary_function_Hx => this%boundary_function_Hx 
  end function mn_get_boundary_function_Hx

  !===============================================================================================
  function mn_get_boundary_function_Hy ( this )
    implicit none
    class(hts_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hy
    mn_get_boundary_function_Hy => this%boundary_function_Hy 
  end function mn_get_boundary_function_Hy

  !===============================================================================================
  function mn_get_boundary_function_Hz ( this )
    implicit none
    class(hts_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mn_get_boundary_function_Hz
    mn_get_boundary_function_Hz => this%boundary_function_Hz 
  end function mn_get_boundary_function_Hz

  !===============================================================================================
  function mn_get_source_term ( this )
    implicit none
    class(hts_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_source_term
    mn_get_source_term => this%source_term
  end function mn_get_source_term
  
    !===============================================================================================
  function mn_get_initial_conditions ( this )
    implicit none
    class(hts_analytical_functions_t), target, intent(in)    :: this
    class(vector_function_t), pointer :: mn_get_initial_conditions
    mn_get_initial_conditions => this%initial_conditions
  end function mn_get_initial_conditions

end module hts_analytical_functions_names
!***************************************************************************************************
