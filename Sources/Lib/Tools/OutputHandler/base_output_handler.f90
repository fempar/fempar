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

module base_output_handler_names

USE FPL
USE types_names
USE memor_names
USE fe_space_names,              only: serial_fe_space_t, fe_iterator_t, fe_accessor_t
USE fe_function_names,           only: fe_function_t
USE output_handler_fe_field_names
USE output_handler_patch_names
USE output_handler_cell_fe_function_names
USE output_handler_fe_iterator_names

implicit none
#include "debug.i90"
private

    integer(ip), parameter :: cell_node_field_array_size = 10

    type, abstract :: base_output_handler_t
    private
        class(serial_fe_space_t),            pointer    :: fe_space => NULL()
        class(output_handler_fe_iterator_t), pointer    :: iterator => NULL()
        type(output_handler_fe_iterator_t)              :: default_iterator
        type(output_handler_fe_field_t),    allocatable :: fe_fields(:)
        type(output_handler_cell_vector_t), allocatable :: cell_vectors(:)
        integer(ip)                                     :: number_cells        = 0
        integer(ip)                                     :: number_nodes        = 0
        integer(ip)                                     :: number_fields       = 0
        integer(ip)                                     :: number_dimensions   = 0
        integer(ip)                                     :: number_cell_vectors = 0
    contains
    private
        procedure, non_overridable, public :: free                         => output_handler_base_free
        procedure, non_overridable, public :: get_number_nodes             => output_handler_base_get_number_nodes
        procedure, non_overridable, public :: get_number_cells             => output_handler_base_get_number_cells
        procedure, non_overridable, public :: get_number_dimensions        => output_handler_base_get_number_dimensions
        procedure, non_overridable, public :: get_number_fields            => output_handler_base_get_number_fields
        procedure, non_overridable, public :: get_number_cell_vectors      => output_handler_base_get_number_cell_vectors
        procedure, non_overridable, public :: get_fe_field                 => output_handler_base_get_fe_field
        procedure, non_overridable, public :: get_cell_vector              => output_handler_base_get_cell_vector
        procedure, non_overridable, public :: attach_fe_space              => output_handler_base_attach_fe_space
        procedure, non_overridable, public :: get_fe_space                 => output_handler_base_get_fe_space
        procedure, non_overridable, public :: set_iterator                 => output_handler_base_set_iterator
        procedure, non_overridable         :: resize_fe_fields_if_needed   => output_handler_base_resize_fe_fields_if_needed
        procedure, non_overridable         :: resize_cell_vectors_if_needed=> output_handler_base_resize_cell_vectors_if_needed
        procedure, non_overridable, public :: add_fe_function              => output_handler_base_add_fe_function
        procedure, non_overridable, public :: add_cell_vector              => output_handler_base_add_cell_vector
        procedure, non_overridable, public :: fill_data                    => output_handler_base_fill_data
        procedure(output_handler_base_open),                           public, deferred :: open
        procedure(output_handler_base_append_time_step),               public, deferred :: append_time_step
        procedure(output_handler_base_allocate_cell_and_nodal_arrays),         deferred :: allocate_cell_and_nodal_arrays
        procedure(output_handler_base_append_cell),                            deferred :: append_cell
        procedure(output_handler_base_write),                          public, deferred :: write
        procedure(output_handler_base_close),                          public, deferred :: close
        procedure(output_handler_base_free_body),                              deferred :: free_body
    end type

    abstract interface
        subroutine output_handler_base_open(this, dir_path, prefix, parameter_list)
            import base_output_handler_t
            import ParameterList_t
            class(base_output_handler_t),    intent(inout) :: this
            character(len=*),                intent(in)    :: dir_path
            character(len=*),                intent(in)    :: prefix
            type(ParameterList_t), optional, intent(in)    :: parameter_list
        end subroutine

        subroutine output_handler_base_append_time_step(this, value)
            import base_output_handler_t
            import rp
            class(base_output_handler_t), intent(inout) :: this
            real(rp),                     intent(in)    :: value
        end subroutine

        subroutine output_handler_base_write(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine output_handler_base_allocate_cell_and_nodal_arrays(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine output_handler_base_append_cell(this, subcell_accessor)
            import base_output_handler_t
            import patch_subcell_accessor_t
            class(base_output_handler_t),   intent(inout) :: this
            type(patch_subcell_accessor_t), intent(in)    :: subcell_accessor
        end subroutine

        subroutine output_handler_base_close(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine output_handler_base_free_body(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine
    end interface

public :: base_output_handler_t

contains

!---------------------------------------------------------------------
!< output_handler_BASE_T PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_base_free(this)
    !-----------------------------------------------------------------
    !< Free output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(inout) :: this
        integer(ip)                                 :: i
    !-----------------------------------------------------------------
        call this%free_body()
        if(allocated(this%fe_fields)) then
            do i=1, size(this%fe_fields)
                call this%fe_fields(i)%free()
            enddo
            deallocate(this%fe_fields)
        endif
        nullify(this%iterator)
        call this%default_iterator%free()
        nullify(this%fe_space)
        this%number_cell_vectors = 0
        this%number_fields       = 0
        this%number_nodes        = 0
        this%number_cells        = 0
        this%number_dimensions   = 0
    end subroutine output_handler_base_free

    subroutine output_handler_base_set_iterator(this, iterator)
    !-----------------------------------------------------------------
    !< Set output handler fe_iterator
    !-----------------------------------------------------------------
        class(base_output_handler_t),                intent(inout) :: this
        class(output_handler_fe_iterator_t), target, intent(in)  :: iterator
    !-----------------------------------------------------------------
        this%iterator => iterator
    end subroutine output_handler_base_set_iterator


    function output_handler_base_get_number_nodes(this) result(number_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_nodes
    !-----------------------------------------------------------------
        number_nodes = this%number_nodes
    end function output_handler_base_get_number_nodes


    function output_handler_base_get_number_cells(this) result(number_cells)
    !-----------------------------------------------------------------
    !< Return the number of cells
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_cells
    !-----------------------------------------------------------------
        number_cells = this%number_cells
    end function output_handler_base_get_number_cells


    function output_handler_base_get_number_dimensions(this) result(number_dimensions)
    !-----------------------------------------------------------------
    !< Return the number of dimensions
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_dimensions
    !-----------------------------------------------------------------
        number_dimensions = this%number_dimensions
    end function output_handler_base_get_number_dimensions


    function output_handler_base_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_fields
    !-----------------------------------------------------------------
        number_fields = this%number_fields
    end function output_handler_base_get_number_fields


    function output_handler_base_get_number_cell_vectors(this) result(number_cell_vectors)
    !-----------------------------------------------------------------
    !< Return the number of cell_vectors
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_cell_vectors
    !-----------------------------------------------------------------
        number_cell_vectors = this%number_cell_vectors
    end function output_handler_base_get_number_cell_vectors


    function output_handler_base_get_fe_field(this, field_id) result(field)
    !-----------------------------------------------------------------
    !< Return a fe field given its id
    !-----------------------------------------------------------------
        class(base_output_handler_t),    target, intent(in) :: this
        integer(ip),                             intent(in) :: field_id
        type(output_handler_fe_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(field_id <= this%number_fields)
        field => this%fe_fields(field_id)
    end function output_handler_base_get_fe_field


    function output_handler_base_get_cell_vector(this, cell_vector_id) result(cell_vector)
    !-----------------------------------------------------------------
    !< Return a cell vector given its id
    !-----------------------------------------------------------------
        class(base_output_handler_t),       target, intent(in) :: this
        integer(ip),                                intent(in) :: cell_vector_id
        type(output_handler_cell_vector_t), pointer            :: cell_vector
    !-----------------------------------------------------------------
        assert(cell_vector_id <= this%number_cell_vectors)
        cell_vector => this%cell_vectors(cell_vector_id)
    end function output_handler_base_get_cell_vector


    subroutine output_handler_base_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t),          intent(inout) :: this
        class(serial_fe_space_t), target,      intent(in)    :: fe_space
    !-----------------------------------------------------------------
        this%fe_space => fe_space
    end subroutine output_handler_base_attach_fe_space


    function output_handler_base_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Return a fe_space pointer
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        class(serial_fe_space_t), pointer        :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function output_handler_base_get_fe_space


    subroutine output_handler_base_resize_fe_fields_if_needed(this, number_fields)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
        integer(ip)                                       :: current_size
        type(output_handler_fe_field_t), allocatable      :: temp_fe_functions(:)
    !-----------------------------------------------------------------
        if(.not. allocated(this%fe_fields)) then
            allocate(this%fe_fields(cell_node_field_array_size))
        elseif(number_fields > size(this%fe_fields)) then
            current_size = size(this%fe_fields)
            call move_alloc(from=this%fe_fields, to=temp_fe_functions)
            allocate(this%fe_fields(int(1.5*current_size)))
            this%fe_fields(1:current_size) = temp_fe_functions(1:current_size)
            deallocate(temp_fe_functions)
        endif
    end subroutine output_handler_base_resize_fe_fields_if_needed


    subroutine output_handler_base_resize_cell_vectors_if_needed(this, number_fields)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
        integer(ip)                                       :: current_size
        type(output_handler_cell_vector_t), allocatable   :: temp_cell_vectors(:)
    !-----------------------------------------------------------------
        if(.not. allocated(this%cell_vectors)) then
            allocate(this%cell_vectors(cell_node_field_array_size))
        elseif(number_fields > size(this%cell_vectors)) then
            current_size = size(this%cell_vectors)
            call move_alloc(from=this%cell_vectors, to=temp_cell_vectors)
            allocate(this%cell_vectors(int(1.5*current_size)))
            this%cell_vectors(1:current_size) = temp_cell_vectors(1:current_size)
            deallocate(temp_cell_vectors)
        endif
    end subroutine output_handler_base_resize_cell_vectors_if_needed


    subroutine output_handler_base_add_fe_function(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Add a fe_function to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        type(fe_function_t),                intent(in)    :: fe_function
        integer(ip),                        intent(in)    :: field_id
        character(len=*),                   intent(in)    :: name
        character(len=*), optional,         intent(in)    :: diff_operator
    !-----------------------------------------------------------------
        call this%resize_fe_fields_if_needed(this%number_fields+1)
        this%number_fields = this%number_fields + 1
        call this%fe_fields(this%number_fields)%set(fe_function, field_id, name, diff_operator)
    end subroutine output_handler_base_add_fe_function


    subroutine output_handler_base_add_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Add a fe_function to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        real(rp), allocatable,              intent(in)    :: cell_vector(:)
        character(len=*),                   intent(in)    :: name
    !-----------------------------------------------------------------
        call this%resize_cell_vectors_if_needed(this%number_cell_vectors+1)
        this%number_cell_vectors = this%number_cell_vectors + 1
        call this%cell_vectors(this%number_cell_vectors)%set(cell_vector, name)
    end subroutine output_handler_base_add_cell_Vector


    subroutine output_handler_base_fill_data(this)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(base_output_handler_t), target, intent(inout) :: this
        type(fe_accessor_t)                                 :: fe
        type(output_handler_cell_fe_function_t)             :: output_handler_cell_function
        type(output_handler_patch_t)                        :: patch
        type(patch_subcell_iterator_t)                      :: subcell_iterator
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))

        ! Ouput_handler_FE_iterator 
        if(.not. associated (this%iterator)) then
            this%iterator => this%default_iterator
            call this%iterator%create(this%fe_space)
        endif
        
        ! Create Output Cell Handler and allocate patch fields
        call output_handler_cell_function%create(this%fe_space, this%iterator, this%number_fields, this%fe_fields(1:this%number_fields))

        ! Allocate geometry and connectivity arrays
        this%number_nodes      = output_handler_cell_function%get_number_nodes()
        this%number_cells      = output_handler_cell_function%get_number_cells()
        this%number_dimensions = output_handler_cell_function%get_number_dimensions()
        call this%allocate_cell_and_nodal_arrays()

        call patch%create(this%number_fields, this%number_cell_vectors)
        ! Translate coordinates and connectivities to VTK format for every subcell
        call this%iterator%init()
        do while ( .not. this%iterator%has_finished())
            ! Get Finite element
            call this%iterator%current(fe)
            if ( fe%is_local() ) then
                call output_handler_cell_function%fill_patch(fe, &
                                                             this%number_fields, &
                                                             this%fe_fields(1:this%number_fields), &
                                                             this%number_cell_vectors, &
                                                             this%cell_vectors(1:this%number_cell_vectors), &
                                                             patch)
                subcell_iterator = patch%get_subcells_iterator()
!               ! Fill data
                do while(.not. subcell_iterator%has_finished())
                    call this%append_cell(subcell_iterator%get_accessor())
                    call subcell_iterator%next()
                enddo
            endif
            call this%iterator%next()
        end do

        call patch%free()
        call output_handler_cell_function%free()
        call fe%free()
    end subroutine output_handler_base_fill_data

end module base_output_handler_names


