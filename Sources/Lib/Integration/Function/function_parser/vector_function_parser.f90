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
module vector_function_parser_names
    use types_names
    use field_names
    use function_names
    use scalar_function_parser_names

    implicit none
# include "debug.i90"

    private

    type, extends(vector_function_t) :: vector_function_parser_t
    private
        type(p_scalar_function_parser_t)  :: components(SPACE_DIM)
        integer                           :: ncomponents = -1
    contains
        procedure :: create_2D            => vector_function_parser_create_2D
        procedure :: create_3D            => vector_function_parser_create_3D
        procedure :: get_value_space      => vector_function_parser_get_value_space
        procedure :: get_value_space_time => vector_function_parser_get_value_space_time
        generic   :: create               => create_2D, create_3D        
    end type vector_function_parser_t

    public :: vector_function_parser_t
        
contains

    subroutine vector_function_parser_create_2D( this, component_1, component_2)
    !-----------------------------------------------------------------
    !< Initialize a time independant vector analytical function
    !-----------------------------------------------------------------
        class(vector_function_parser_t),        intent(inout) :: this
        type(scalar_function_parser_t), target, intent(in)    :: component_1
        type(scalar_function_parser_t), target, intent(in)    :: component_2
    !----------------------------------------------------------------- 
        this%ncomponents = 2
        call this%set_num_dims(2)
        this%components(1)%function => component_1
        this%components(2)%function => component_2
    end subroutine vector_function_parser_create_2D


    subroutine vector_function_parser_create_3D( this, component_1, component_2, component_3)
    !-----------------------------------------------------------------
    !< Initialize a time independant vector analytical function
    !-----------------------------------------------------------------
        class(vector_function_parser_t),        intent(inout) :: this
        type(scalar_function_parser_t), target, intent(in)    :: component_1
        type(scalar_function_parser_t), target, intent(in)    :: component_2
        type(scalar_function_parser_t), target, intent(in)    :: component_3
    !----------------------------------------------------------------- 
        this%ncomponents = 3
        call this%set_num_dims(3)
        this%components(1)%function => component_1
        this%components(2)%function => component_2
        this%components(3)%function => component_3
    end subroutine vector_function_parser_create_3D


    subroutine vector_function_parser_get_value_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant vector analytical function in a given point
    !-----------------------------------------------------------------
        class(vector_function_parser_t),     intent(in)    :: this
        type(point_t),                       intent(in)    :: point
        type(vector_field_t),                intent(inout) :: result
        real(rp)                                           :: tmp_result
        integer                                            :: i, dims
    !-----------------------------------------------------------------
        dims = this%get_num_dims()
        do i=1, dims
            call this%components(i)%function%get_value_space(point, tmp_result)
            call result%set(i, tmp_result)
        enddo
    end subroutine vector_function_parser_get_value_space


    subroutine vector_function_parser_get_value_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant vector analytical function in a given point and time step
    !-----------------------------------------------------------------
        class(vector_function_parser_t), intent(in)    :: this
        type(point_t),                   intent(in)    :: point
        real(rp),                        intent(in)    :: time
        type(vector_field_t),            intent(inout) :: result
        real(rp)                                       :: tmp_result
        integer                                        :: i, dims
    !-----------------------------------------------------------------
        dims = this%get_num_dims()
        do i=1, dims
            call this%components(i)%function%get_value_space_time(point, time,  tmp_result)
            call result%set(i, tmp_result)
        enddo
    end subroutine vector_function_parser_get_value_space_time

end module vector_function_parser_names
