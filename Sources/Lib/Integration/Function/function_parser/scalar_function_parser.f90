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
module scalar_function_parser_names
    use types_names
    use field_names
    use function_names
    use FortranParser, only: EquationParser

    implicit none
# include "debug.i90"

    private

    character(len=1), dimension(4), parameter :: variables = ["x", "y", "z", "t"]

    type, extends(scalar_function_t) :: scalar_function_parser_t
    private
        type(EquationParser) :: equation
    contains
        procedure :: create               => scalar_function_parser_create
        procedure :: get_value_space      => scalar_function_parser_get_value_space
        procedure :: get_value_space_time => scalar_function_parser_get_value_space_time
        procedure :: free                 => scalar_function_parser_free
    end type scalar_function_parser_t

    type :: p_scalar_function_parser_t
        type(scalar_function_parser_t), pointer    :: function => null()
    end type

    public :: scalar_function_parser_t, p_scalar_function_parser_t

contains

    subroutine scalar_function_parser_create(this, expression, num_dims)
    !-----------------------------------------------------------------
    !< Initialize the time independant scalar analytical function   
    !-----------------------------------------------------------------
        class(scalar_function_parser_t),     intent(inout) :: this
        character(len=*),                    intent(in)    :: expression
        integer,                             intent(in)    :: num_dims
    !-----------------------------------------------------------------
        assert(num_dims==2 .or. num_dims==3) 
        call this%free()
        call this%set_num_dims(num_dims)
        this%equation = EquationParser(expression, variables)
        assert(this%equation%Error == 0)
    end subroutine scalar_function_parser_create


    subroutine scalar_function_parser_get_value_space( this, point, result )
    !-----------------------------------------------------------------
    !< Evaluate the time independant scalar analytical function in a given point
    !-----------------------------------------------------------------
        class(scalar_function_parser_t),     intent(in)    :: this
        type(point_t),                       intent(in)    :: point
        real(rp),                            intent(inout) :: result
        real(rp), dimension(4)                             :: values
    !-----------------------------------------------------------------
        values(1:SPACE_DIM) = point%get_value()
        values(SPACE_DIM+1:4) = 0._rp
        result = this%equation%evaluate(values)
        assert(this%equation%Error == 0)
    end subroutine scalar_function_parser_get_value_space


    subroutine scalar_function_parser_get_value_space_time( this, point, time, result )
    !-----------------------------------------------------------------
    !< Evaluate the time dependant scalar function in a given point and time step
    !-----------------------------------------------------------------
        class(scalar_function_parser_t),    intent(in)    :: this
        type(point_t),                      intent(in)    :: point
        real(rp),                           intent(in)    :: time
        real(rp),                           intent(inout) :: result
        real(rp), dimension(4)                            :: values
    !-----------------------------------------------------------------
        values(1:SPACE_DIM) = point%get_value()
        values(SPACE_DIM+1:3) = 0._rp
        values(4) = time
        result = this%equation%evaluate(values)
        assert(this%equation%Error == 0)
    end subroutine scalar_function_parser_get_value_space_time


    subroutine scalar_function_parser_free(this)
    !-----------------------------------------------------------------
    !< Free scalar analytical function
    !-----------------------------------------------------------------
        class(scalar_function_parser_t),     intent(inout) :: this
    !----------------------------------------------------------------- 
        !call this%equation%finalize()
        call this%set_num_dims(-1)
    end subroutine scalar_function_parser_free

end module scalar_function_parser_names
