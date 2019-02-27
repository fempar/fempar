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
module scalar_function_and_gradient_parser_names
    use types_names
    use field_names
    use function_names
    use scalar_function_parser_names
    use vector_function_parser_names

    implicit none
# include "debug.i90"

    private

    type, extends(scalar_function_t) :: scalar_function_and_gradient_parser_t
    private
        type(scalar_function_parser_t), pointer :: function => null()
        type(vector_function_parser_t), pointer :: gradient => null()
        type(p_scalar_function_parser_t)        :: time_derivatives(0:1)
    contains
        procedure :: create                        => scalar_function_and_gradient_parser_create
        procedure :: get_value_space               => scalar_function_and_gradient_parser_get_value_space
        procedure :: get_value_space_time          => scalar_function_and_gradient_parser_get_value_space_time
        procedure :: get_gradient_space            => scalar_function_and_gradient_parser_get_gradient_space
        procedure :: get_gradient_space_time       => scalar_function_and_gradient_parser_get_gradient_space_time
        procedure :: get_value_temporal_derivative => scalar_function_and_gradient_parser_get_temporal_derivative
        procedure :: free                          => scalar_function_and_gradient_parser_free
    end type scalar_function_and_gradient_parser_t

    public :: scalar_function_and_gradient_parser_t

contains

    subroutine scalar_function_and_gradient_parser_create(this, function, gradient, first_order_time_derivative)
    !-----------------------------------------------------------------
    !< Initialize the time independant scalar analytical function
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t),     intent(inout) :: this
        type(scalar_function_parser_t), target,           intent(in)    :: function
        type(vector_function_parser_t), target,           intent(in)    :: gradient
        type(scalar_function_parser_t), target, optional, intent(in)    :: first_order_time_derivative
    !-----------------------------------------------------------------
        assert(function%get_num_dims() == gradient%get_num_dims())
        call this%free() 
        this%function => function
        this%gradient => gradient
        this%time_derivatives(0)%function => function
        if(present(first_order_time_derivative)) then
            assert(function%get_num_dims() == first_order_time_derivative%get_num_dims()) 
            this%time_derivatives(1)%function => first_order_time_derivative
        endif 
    end subroutine scalar_function_and_gradient_parser_create


    subroutine scalar_function_and_gradient_parser_get_value_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant scalar analytical function in a given point
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(inout) :: result
    !-----------------------------------------------------------------
        call this%function%get_value_space(point, result)
    end subroutine scalar_function_and_gradient_parser_get_value_space


    subroutine scalar_function_and_gradient_parser_get_value_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant scalar function in a given point and time step
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(in)    :: time
        real(rp),                                     intent(inout) :: result
    !-----------------------------------------------------------------
        call this%function%get_value_space_time(point, time, result)
    end subroutine scalar_function_and_gradient_parser_get_value_space_time


    subroutine scalar_function_and_gradient_parser_get_gradient_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant gradient in a given point
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        type(vector_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%gradient%get_value_space(point, result)
    end subroutine scalar_function_and_gradient_parser_get_gradient_space


    subroutine scalar_function_and_gradient_parser_get_gradient_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant gradient in a given point and time step
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(in)    :: time
        type(vector_field_t),                         intent(inout) :: result
    !-----------------------------------------------------------------
        call this%gradient%get_value_space_time(point, time, result)
    end subroutine scalar_function_and_gradient_parser_get_gradient_space_time


    subroutine scalar_function_and_gradient_parser_get_temporal_derivative( this, point, time, time_derivative_order, result )
    !-----------------------------------------------------------------
    !< Evaluate the time derivative in a given point and time step
    !< Only zero and first time derivative order supported
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t), intent(in)    :: this
        type(point_t),                                intent(in)    :: point
        real(rp),                                     intent(in)    :: time
        integer(ip),                                  intent(in)    :: time_derivative_order
        real(rp),                                     intent(inout) :: result
    !-----------------------------------------------------------------
        massert( time_derivative_order >= 0 .or. time_derivative_order <= 1, "Only order 0 and 1 time derivatives supported"  )
        assert( associated(this%time_derivatives(time_derivative_order)%function)  )
        call this%time_derivatives(time_derivative_order)%function%get_value_space_time(point, time, result)
    end subroutine scalar_function_and_gradient_parser_get_temporal_derivative


    subroutine scalar_function_and_gradient_parser_free(this)
    !-----------------------------------------------------------------
    !< Free scalar function, gradient and time derivatives
    !-----------------------------------------------------------------
        class(scalar_function_and_gradient_parser_t),     intent(inout) :: this
    !-----------------------------------------------------------------
        nullify(this%function)
        nullify(this%gradient)
        nullify(this%time_derivatives(0)%function)
        nullify(this%time_derivatives(1)%function) 
    end subroutine scalar_function_and_gradient_parser_free


end module scalar_function_and_gradient_parser_names
