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
module function_library_names
  use types_names
  use field_names
  use function_names

  implicit none
# include "debug.i90"

  private

  type, extends(scalar_function_t) :: constant_scalar_function_t
     private
     real(rp) :: function_value
   contains
     procedure :: constant_scalar_function_create
     generic   :: create                    => constant_scalar_function_create
     procedure :: get_value_space           => constant_scalar_function_get_value_space
     procedure :: get_value_space_time      => constant_scalar_function_get_value_space_time
     procedure :: get_values_set_space      => constant_scalar_function_get_values_set_space
     procedure :: get_values_set_space_time => constant_scalar_function_get_values_set_space_time
  end type constant_scalar_function_t
  
  interface constant_scalar_function_t
    module procedure constant_scalar_function_constructor
  end interface constant_scalar_function_t
  
  type, extends(vector_function_t) :: constant_vector_function_t
     private
     type(vector_field_t) :: function_value
   contains
     procedure :: constant_vector_function_create
     generic   :: create                    => constant_vector_function_create
     procedure :: get_value_space           => constant_vector_function_get_value_space
     procedure :: get_value_space_time      => constant_vector_function_get_value_space_time
     procedure :: get_values_set_space      => constant_vector_function_get_values_set_space
     procedure :: get_values_set_space_time => constant_vector_function_get_values_set_space_time
  end type constant_vector_function_t
  
  interface constant_vector_function_t
    module procedure constant_vector_function_constructor
  end interface constant_vector_function_t
  
  public :: constant_scalar_function_t, constant_vector_function_t

contains
  subroutine constant_scalar_function_create( this, number_dimensions, function_value )
    class(constant_scalar_function_t), intent(inout)     :: this
    integer(ip)                      , intent(in)        :: number_dimensions
    real(rp)                         , intent(in)        :: function_value
    call this%create(number_dimensions)
    this%function_value = function_value
  end subroutine constant_scalar_function_create

  subroutine constant_scalar_function_get_value_space( this, point, result )
    class(constant_scalar_function_t), intent(in)     :: this
    type(point_t)                    , intent(in)     :: point
    real(rp)                         , intent(inout)  :: result
    result = this%function_value
  end subroutine constant_scalar_function_get_value_space

  subroutine constant_scalar_function_get_value_space_time( this, point, time, result )
    class(constant_scalar_function_t), intent(in)    :: this
    type(point_t)                    , intent(in)    :: point
    real(rp)                         , intent(in)    :: time
    real(rp)                         , intent(inout) :: result
    result = this%function_value
  end subroutine constant_scalar_function_get_value_space_time

  subroutine constant_scalar_function_get_values_set_space( this, point, result )
    class(constant_scalar_function_t), intent(in)     :: this
    type(point_t)                    , intent(in)     :: point(:)
    real(rp)                         , intent(inout)  :: result(:)
    result = this%function_value
  end subroutine constant_scalar_function_get_values_set_space

  subroutine constant_scalar_function_get_values_set_space_time( this, point, time, result )
    class(constant_scalar_function_t), intent(in)    :: this
    type(point_t)                    , intent(in)    :: point(:)
    real(rp)                         , intent(in)    :: time(:)
    real(rp)                         , intent(inout) :: result(:,:)
    result = this%function_value
  end subroutine constant_scalar_function_get_values_set_space_time
  
  function constant_scalar_function_constructor ( number_dimensions, function_value ) result(new_scalar_constant_function)
    integer(ip), intent(in) :: number_dimensions
    real(rp)   , intent(in) :: function_value
    type(constant_scalar_function_t) :: new_scalar_constant_function
    call new_scalar_constant_function%create(number_dimensions, function_value)
  end function constant_scalar_function_constructor
  
  subroutine constant_vector_function_create( this, number_dimensions, function_value )
    class(constant_vector_function_t), intent(inout)     :: this
    integer(ip)                      , intent(in)        :: number_dimensions
    type(vector_field_t)             , intent(in)        :: function_value
    integer(ip)                                          :: idime
    
    call this%create(number_dimensions)
    call this%function_value%init(0.0_rp)
    do idime = 1, this%get_number_dimensions()
      call this%function_value%set(idime,function_value%get(idime))
    end do
  end subroutine constant_vector_function_create

  subroutine constant_vector_function_get_value_space( this, point, result )
    class(constant_vector_function_t), intent(in)     :: this
    type(point_t)                    , intent(in)     :: point
    type(vector_field_t)             , intent(inout)  :: result
    result = this%function_value
  end subroutine constant_vector_function_get_value_space

  subroutine constant_vector_function_get_value_space_time( this, point, time, result )
    class(constant_vector_function_t), intent(in)    :: this
    type(point_t)                    , intent(in)    :: point
    real(rp)                         , intent(in)    :: time
    type(vector_field_t)             , intent(inout) :: result
    result = this%function_value
  end subroutine constant_vector_function_get_value_space_time

  subroutine constant_vector_function_get_values_set_space( this, point, result )
    class(constant_vector_function_t), intent(in)     :: this
    type(point_t)                    , intent(in)     :: point(:)
    type(vector_field_t)             , intent(inout)  :: result(:)
    result = this%function_value
  end subroutine constant_vector_function_get_values_set_space

  subroutine constant_vector_function_get_values_set_space_time( this, point, time, result )
    class(constant_vector_function_t), intent(in)    :: this
    type(point_t)                    , intent(in)    :: point(:)
    real(rp)                         , intent(in)    :: time(:)
    type(vector_field_t)             , intent(inout) :: result(:,:)
    result = this%function_value
  end subroutine constant_vector_function_get_values_set_space_time
  
  function constant_vector_function_constructor ( number_dimensions, function_value ) result(new_vector_constant_function)
    integer(ip)         , intent(in) :: number_dimensions
    type(vector_field_t), intent(in) :: function_value
    type(constant_vector_function_t) :: new_vector_constant_function
    call new_vector_constant_function%create(number_dimensions, function_value)
  end function constant_vector_function_constructor
  
end module function_library_names
