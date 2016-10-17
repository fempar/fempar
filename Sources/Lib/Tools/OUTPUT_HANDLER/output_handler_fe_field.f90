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

module output_handler_fe_field_names

USE types_names
USE fe_function_names,           only: fe_function_t

implicit none
#include "debug.i90"
private

    type :: output_handler_fe_field_t
    private
        character(len=:),  allocatable :: name
        integer(ip)                    :: field_id = 0
        type(fe_function_t), pointer   :: fe_function => NULL()
    contains
        procedure         ::                    output_handler_fe_field_assign
        procedure, public :: set             => output_handler_fe_field_set
        procedure, public :: get_name        => output_handler_fe_field_get_name
        procedure, public :: get_field_id    => output_handler_fe_field_get_field_id
        procedure, public :: get_fe_function => output_handler_fe_field_get_fe_function
        procedure, public :: free            => output_handler_fe_field_free
        generic,   public :: assignment(=)   => output_handler_fe_field_assign
    end type

public :: output_handler_fe_field_t

contains

!---------------------------------------------------------------------
!< output_handler_fe_field_t PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_fe_field_free(this)
    !-----------------------------------------------------------------
    !< Free output_handler_fe_field_t derived type
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%name)) deallocate(this%name)
        nullify(this%fe_function)
        this%field_id = 0
    end subroutine output_handler_fe_field_free


    subroutine output_handler_fe_field_assign(this, output_handler_fe_field)
    !----------------------------------------------------------------- 
    !< Assign operator overloading for output_handler_fe_field_t
    !----------------------------------------------------------------- 
        class(output_handler_fe_field_t), intent(inout) :: this
        type(output_handler_fe_field_t),  intent(in)    :: output_handler_fe_field
        type(fe_function_t), pointer                    :: fe_function
        integer(ip),                                    :: field_id
        character(len=:), allocatable                   :: name
    !----------------------------------------------------------------- 
        call this%free()
        fe_function => output_handler_fe_field%get_fe_function() 
        field_id    =  output_handler_fe_field%get_field_id()
        name        =  output_handler_fe_field%get_name()
        call this%set(fe_function, field_id, name)
    end subroutine output_handler_fe_field_assign


    subroutine output_handler_fe_field_set(this, fe_function, field_id, name)
    !-----------------------------------------------------------------
    !< Associate a fe_function with a field name
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t),  intent(inout) :: this
        type(fe_function_t), target,       intent(in)    :: fe_function
        integer(ip),                       intent(in)    :: field_id
        character(len=*),                  intent(in)    :: name
    !-----------------------------------------------------------------
        call this%free()
        this%fe_function => fe_function
        this%field_id    = field_id
        this%name        = name
    end subroutine output_handler_fe_field_set


    function output_handler_fe_field_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t), intent(in) :: this
        character(len=:), allocatable                :: name
    !-----------------------------------------------------------------
        assert(allocated(this%name))
        name = this%name
    end function output_handler_fe_field_get_name


    function output_handler_fe_field_get_field_id(this) result(field_id)
    !-----------------------------------------------------------------
    !< Return the field_id of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t), intent(in) :: this
        integer(ip)                                  :: field_id
    !-----------------------------------------------------------------
        assert(this%field_id /= 0)
        field_id = this%field_id
    end function output_handler_fe_field_get_field_id


    function output_handler_fe_field_get_fe_function(this) result(fe_function)
    !-----------------------------------------------------------------
    !< Return a fe_function pointer
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t), target, intent(in) :: this
        type(fe_function_t), pointer                         :: fe_function
    !-----------------------------------------------------------------
        assert(associated(this%fe_function))
        fe_function => this%fe_function
    end function output_handler_fe_field_get_fe_function

end module output_handler_fe_field_names
