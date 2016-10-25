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
USE memor_names
USE fe_function_names,           only: fe_function_t
USE output_handler_parameters_names

implicit none
#include "debug.i90"
private

    type :: output_handler_fe_field_t
    private
        character(len=:),  allocatable :: name
        integer(ip)                    :: field_id = 0
        type(fe_function_t), pointer   :: fe_function => NULL()
        character(len=:),  allocatable :: diff_operator
    contains
        procedure, non_overridable         ::                          output_handler_fe_field_assign
        procedure, non_overridable, public :: set                   => output_handler_fe_field_set
        procedure, non_overridable, public :: get_name              => output_handler_fe_field_get_name
        procedure, non_overridable, public :: get_diff_operator     => output_handler_fe_field_get_diff_operator
        procedure, non_overridable, public :: get_field_id          => output_handler_fe_field_get_field_id
        procedure, non_overridable, public :: get_fe_function       => output_handler_fe_field_get_fe_function
        procedure, non_overridable, public :: free                  => output_handler_fe_field_free
        generic,                    public :: assignment(=)         => output_handler_fe_field_assign
    end type


    type :: output_handler_cell_vector_t
    private
        character(len=:),  allocatable :: name
        real(rp), pointer              :: cell_vector(:) => NULL()
    contains
        procedure, non_overridable         ::                          output_handler_cell_vector_assign
        procedure, non_overridable, public :: set                   => output_handler_cell_vector_set
        procedure, non_overridable, public :: get_name              => output_handler_cell_vector_get_name
        procedure, non_overridable, public :: get_cell_vector       => output_handler_cell_vector_get_cell_vector
        procedure, non_overridable, public :: free                  => output_handler_cell_vector_free
        generic,                    public :: assignment(=)         => output_handler_cell_vector_assign
    end type

    type :: output_handler_fe_field_1D_value_t
    private
        integer(ip)                    :: number_components = 0
        real(rp), allocatable          :: value(:)
    contains
        procedure, non_overridable, public :: value_is_allocated    => output_handler_fe_field_1D_value_value_is_allocated
        procedure, non_overridable, public :: allocate_value        => output_handler_fe_field_1D_value_allocate_value
        procedure, non_overridable, public :: get_value             => output_handler_fe_field_1D_value_get_value
        procedure, non_overridable, public :: get_number_components => output_handler_fe_field_1D_value_get_number_components
        procedure, non_overridable, public :: free                  => output_handler_fe_field_1D_value_free
    end type

    type :: output_handler_fe_field_2D_value_t
    private
        integer(ip)                    :: number_components = 0
        real(rp), allocatable          :: value(:,:)
    contains
        procedure, non_overridable, public :: value_is_allocated    => output_handler_fe_field_2D_value_value_is_allocated
        procedure, non_overridable, public :: allocate_value        => output_handler_fe_field_2D_value_allocate_value
        procedure, non_overridable, public :: get_value             => output_handler_fe_field_2D_value_get_value
        procedure, non_overridable, public :: get_number_components => output_handler_fe_field_2D_value_get_number_components
        procedure, non_overridable, public :: free                  => output_handler_fe_field_2D_value_free
    end type

public :: output_handler_fe_field_t
public :: output_handler_cell_vector_t
public :: output_handler_fe_field_1D_value_t
public :: output_handler_fe_field_2D_value_t

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
        this%field_id          = 0
    end subroutine output_handler_fe_field_free


    subroutine output_handler_fe_field_assign(this, output_handler_fe_field)
    !----------------------------------------------------------------- 
    !< Assign operator overloading for output_handler_fe_field_t
    !----------------------------------------------------------------- 
        class(output_handler_fe_field_t), intent(inout) :: this
        type(output_handler_fe_field_t),  intent(in)    :: output_handler_fe_field
        type(fe_function_t), pointer                    :: fe_function
        real(rp),            pointer                    :: value(:,:)
        integer(ip)                                     :: field_id
        character(len=:), allocatable                   :: name
    !----------------------------------------------------------------- 
        call this%free()
        fe_function => output_handler_fe_field%get_fe_function() 
        field_id    =  output_handler_fe_field%get_field_id()
        name        =  output_handler_fe_field%get_name()
        call this%set(fe_function, field_id, name)
    end subroutine output_handler_fe_field_assign


    subroutine output_handler_fe_field_set(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Associate a fe_function with a field name
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t),  intent(inout) :: this
        type(fe_function_t), target,       intent(in)    :: fe_function
        integer(ip),                       intent(in)    :: field_id
        character(len=*),                  intent(in)    :: name
        character(len=*), optional,        intent(in)    :: diff_operator
    !-----------------------------------------------------------------
        call this%free()
        this%fe_function   => fe_function
        this%field_id      = field_id
        this%name          = name
        this%diff_operator = no_diff_operator
        if(present(diff_operator)) this%diff_operator = diff_operator
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


    function output_handler_fe_field_get_diff_operator(this) result(diff_operator)
    !-----------------------------------------------------------------
    !< Return the diff_operator of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_fe_field_t), intent(in) :: this
        character(len=:), allocatable                :: diff_operator
    !-----------------------------------------------------------------
        assert(allocated(this%diff_operator))
        diff_operator = this%diff_operator
    end function output_handler_fe_field_get_diff_operator


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


!---------------------------------------------------------------------
!< output_handler_CELL_VECTOR_t PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_cell_vector_free(this)
    !-----------------------------------------------------------------
    !< Free output_handler_cell_vector_t derived type
    !-----------------------------------------------------------------
        class(output_handler_cell_vector_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%name)) deallocate(this%name)
        nullify(this%cell_vector)
    end subroutine output_handler_cell_vector_free


    subroutine output_handler_cell_vector_assign(this, output_handler_cell_vector)
    !----------------------------------------------------------------- 
    !< Assign operator overloading for output_handler_cell_vector_t
    !----------------------------------------------------------------- 
        class(output_handler_cell_vector_t), intent(inout) :: this
        type(output_handler_cell_vector_t),  intent(in)    :: output_handler_cell_vector
        real(rp), pointer                                  :: cell_vector(:)
        character(len=:), allocatable                      :: name
    !----------------------------------------------------------------- 
        call this%free()
        cell_vector => output_handler_cell_vector%get_cell_vector() 
        name        =  output_handler_cell_vector%get_name()
        call this%set(cell_vector, name)
    end subroutine output_handler_cell_vector_assign


    subroutine output_handler_cell_vector_set(this, cell_Vector, name)
    !-----------------------------------------------------------------
    !< Associate a fe_function with a field name
    !-----------------------------------------------------------------
        class(output_handler_cell_vector_t),  intent(inout) :: this
        real(rp),               target,       intent(in)    :: cell_vector(:)
        character(len=*),                     intent(in)    :: name
    !-----------------------------------------------------------------
        call this%free()
        this%cell_vector   => cell_vector
        this%name          = name
    end subroutine output_handler_cell_vector_set


    function output_handler_cell_vector_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(output_handler_cell_vector_t), intent(in) :: this
        character(len=:), allocatable                :: name
    !-----------------------------------------------------------------
        assert(allocated(this%name))
        name = this%name
    end function output_handler_cell_vector_get_name


    function output_handler_cell_vector_get_cell_vector(this) result(cell_vector)
    !-----------------------------------------------------------------
    !< Return a fe_function pointer
    !-----------------------------------------------------------------
        class(output_handler_cell_vector_t), target, intent(in) :: this
        real(rp), pointer                                       :: cell_vector(:)
    !-----------------------------------------------------------------
        assert(associated(this%cell_vector))
        cell_vector => this%cell_vector
    end function output_handler_cell_vector_get_cell_vector


!---------------------------------------------------------------------
!< output_handler_fe_field_1D_value_t PROCEDURES
!---------------------------------------------------------------------


    subroutine output_handler_fe_field_1D_value_free(this)
    !-----------------------------------------------------------------
    !< Check if value is allocated
    !-----------------------------------------------------------------
        class(output_handler_fe_field_1D_value_t), intent(inout) :: this
        logical                                                  :: value_is_allocated
    !-----------------------------------------------------------------
        if(allocated(this%value)) call memfree(this%value, __FILE__, __LINE__)
        this%number_components = 0
    end subroutine output_handler_fe_field_1D_value_free


    function output_handler_fe_field_1D_value_value_is_allocated(this) result(value_is_allocated)
    !-----------------------------------------------------------------
    !< Check if value is allocated
    !-----------------------------------------------------------------
        class(output_handler_fe_field_1D_value_t), intent(in) :: this
        logical                                               :: value_is_allocated
    !-----------------------------------------------------------------
        value_is_allocated =  allocated(this%value)
    end function output_handler_fe_field_1D_value_value_is_allocated


    subroutine output_handler_fe_field_1D_value_allocate_value(this, number_components, number_nodes) 
    !-----------------------------------------------------------------
    !< Allocate value
    !-----------------------------------------------------------------
        class(output_handler_fe_field_1D_value_t), intent(inout) :: this
        integer(ip)                                              :: number_components
        integer(ip)                                              :: number_nodes
    !-----------------------------------------------------------------
        assert(.not. allocated(this%value))
        call memalloc(number_components*number_nodes, this%value, __FILE__, __LINE__)
        this%number_components = number_components
    end subroutine output_handler_fe_field_1D_value_allocate_value


    function output_handler_fe_field_1D_value_get_value(this) result(value)
    !-----------------------------------------------------------------
    !< Return a pointer to the raw values
    !-----------------------------------------------------------------
        class(output_handler_fe_field_1D_value_t), target, intent(in) :: this
        real(rp),                                  pointer            :: value(:)
    !-----------------------------------------------------------------
        value => this%value
    end function output_handler_fe_field_1D_value_get_value


    function output_handler_fe_field_1D_value_get_number_components(this) result(number_components)
    !-----------------------------------------------------------------
    !< Return the number of components
    !-----------------------------------------------------------------
        class(output_handler_fe_field_1D_value_t), target, intent(in) :: this
        integer(ip)                                                   :: number_components
    !-----------------------------------------------------------------
        number_components = this%number_components
    end function output_handler_fe_field_1D_value_get_number_components


!---------------------------------------------------------------------
!< output_handler_fe_field_2D_value_t PROCEDURES
!---------------------------------------------------------------------


    subroutine output_handler_fe_field_2D_value_free(this)
    !-----------------------------------------------------------------
    !< Check if value is allocated
    !-----------------------------------------------------------------
        class(output_handler_fe_field_2D_value_t), intent(inout) :: this
        logical                                                  :: value_is_allocated
    !-----------------------------------------------------------------
        if(allocated(this%value)) call memfree(this%value, __FILE__, __LINE__)
        this%number_components = 0
    end subroutine output_handler_fe_field_2D_value_free


    function output_handler_fe_field_2D_value_value_is_allocated(this) result(value_is_allocated)
    !-----------------------------------------------------------------
    !< Check if value is allocated
    !-----------------------------------------------------------------
        class(output_handler_fe_field_2D_value_t), intent(in) :: this
        logical                                               :: value_is_allocated
    !-----------------------------------------------------------------
        value_is_allocated =  allocated(this%value)
    end function output_handler_fe_field_2D_value_value_is_allocated


    subroutine output_handler_fe_field_2D_value_allocate_value(this, number_components, number_nodes) 
    !-----------------------------------------------------------------
    !< Allocate value
    !-----------------------------------------------------------------
        class(output_handler_fe_field_2D_value_t), intent(inout) :: this
        integer(ip)                                              :: number_components
        integer(ip)                                              :: number_nodes
    !-----------------------------------------------------------------
        assert(.not. allocated(this%value))
        call memalloc(number_components, number_nodes, this%value, __FILE__, __LINE__)
        this%number_components = number_components
    end subroutine output_handler_fe_field_2D_value_allocate_value


    function output_handler_fe_field_2D_value_get_value(this) result(value)
    !-----------------------------------------------------------------
    !< Return a pointer to the raw values
    !-----------------------------------------------------------------
        class(output_handler_fe_field_2D_value_t), target, intent(in) :: this
        real(rp),                                  pointer            :: value(:,:)
    !-----------------------------------------------------------------
        value => this%value
    end function output_handler_fe_field_2D_value_get_value


    function output_handler_fe_field_2D_value_get_number_components(this) result(number_components)
    !-----------------------------------------------------------------
    !< Return the number of components
    !-----------------------------------------------------------------
        class(output_handler_fe_field_2D_value_t), target, intent(in) :: this
        integer(ip)                                                   :: number_components
    !-----------------------------------------------------------------
        number_components = this%number_components
    end function output_handler_fe_field_2D_value_get_number_components

end module output_handler_fe_field_names
