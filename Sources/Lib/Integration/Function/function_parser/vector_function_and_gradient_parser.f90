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
module vector_function_and_gradient_parser_names
    use types_names
    use field_names
    use function_names
    use vector_function_parser_names
    use tensor_function_parser_names

    implicit none
# include "debug.i90"

    private

    type, extends(vector_function_t) :: vector_function_and_gradient_parser_t
    private
        type(vector_function_parser_t), pointer :: function => null()
        type(tensor_function_parser_t), pointer :: gradient => null()
    contains
        procedure :: create                       => vector_function_and_gradient_parser_create
        procedure :: get_value_space              => vector_function_and_gradient_parser_get_value_space
        procedure :: get_value_space_time         => vector_function_and_gradient_parser_get_value_space_time
        procedure :: get_gradient_space           => vector_function_and_gradient_parser_get_gradient_space
        procedure :: get_gradient_space_time      => vector_function_and_gradient_parser_get_gradient_space_time
    end type vector_function_and_gradient_parser_t

    public :: vector_function_and_gradient_parser_t

contains

    subroutine vector_function_and_gradient_parser_create(this, function, gradient)
    !-----------------------------------------------------------------
    !< Initialize the time independant vector analytical function
    !-----------------------------------------------------------------
        class(vector_function_and_gradient_parser_t), intent(inout) :: this
        type(vector_function_parser_t), target,       intent(in)    :: function
        type(tensor_function_parser_t), target,       intent(in)    :: gradient
    !----------------------------------------------------------------- 
        assert(function%get_num_dims() == gradient%get_num_dims())
        this%function => function
        this%gradient => gradient
    end subroutine vector_function_and_gradient_parser_create


    subroutine vector_function_and_gradient_parser_get_value_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant vector analytical function in a given point
    !-----------------------------------------------------------------
        class(vector_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        type(vector_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%function%get_value_space(point, result)
    end subroutine vector_function_and_gradient_parser_get_value_space


    subroutine vector_function_and_gradient_parser_get_value_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant vector function in a given point and time step
    !-----------------------------------------------------------------
        class(vector_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(in)    :: time
        type(vector_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%function%get_value_space_time(point, time, result)
    end subroutine vector_function_and_gradient_parser_get_value_space_time


    subroutine vector_function_and_gradient_parser_get_gradient_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant gradient in a given point
    !-----------------------------------------------------------------
        class(vector_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        type(tensor_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%gradient%get_value_space(point, result)
    end subroutine vector_function_and_gradient_parser_get_gradient_space


    subroutine vector_function_and_gradient_parser_get_gradient_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant gradient in a given point and time step
    !-----------------------------------------------------------------
        class(vector_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(in)    :: time
        type(tensor_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%gradient%get_value_space_time(point, time, result)
    end subroutine vector_function_and_gradient_parser_get_gradient_space_time

end module vector_function_and_gradient_parser_names
