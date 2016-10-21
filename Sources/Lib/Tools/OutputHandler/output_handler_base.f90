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

module output_handler_base_names

USE FPL
USE types_names
USE fe_space_names,              only: serial_fe_space_t, fe_iterator_t, fe_accessor_t
USE fe_function_names,           only: fe_function_t
USE output_handler_fe_field_names
USE output_handler_patch_names
USE output_handler_cell_fe_function_names
USE output_handler_fe_iterator_names

implicit none
#include "debug.i90"
private

    integer(ip), parameter :: fe_functions_initial_size = 10

    type, abstract :: output_handler_base_t
    private
        class(serial_fe_space_t),            pointer :: fe_space => NULL()
        class(output_handler_fe_iterator_t), pointer :: iterator => NULL()
        type(output_handler_fe_field_t), allocatable :: fields(:)
        integer(ip)                                  :: number_cells  = 0
        integer(ip)                                  :: number_nodes  = 0
        integer(ip)                                  :: number_fields = 0
    contains
        procedure,                  public :: free                       => output_handler_base_free
        procedure, non_overridable, public :: get_number_nodes           => output_handler_base_get_number_nodes
        procedure, non_overridable, public :: get_number_cells           => output_handler_base_get_number_cells
        procedure, non_overridable, public :: get_number_fields          => output_handler_base_get_number_fields
        procedure, non_overridable, public :: get_field                  => output_handler_base_get_field
        procedure, non_overridable, public :: attach_fe_space            => output_handler_base_attach_fe_space
        procedure, non_overridable, public :: get_fe_space               => output_handler_base_get_fe_space
        procedure, non_overridable, public :: set_iterator               => output_handler_base_set_iterator
        procedure, non_overridable, public :: set_default_iterator       => output_handler_base_set_default_iterator
        procedure, non_overridable         :: get_default_iterator       => output_handler_base_get_default_iterator
        procedure, non_overridable         :: resize_fields_if_needed    => output_handler_base_resize_fields_if_needed
        procedure, non_overridable, public :: add_fe_function            => output_handler_base_add_fe_function
        procedure, non_overridable, public :: fill_data                  => output_handler_base_fill_data
        procedure(output_handler_base_open),                           public, deferred :: open
        procedure(output_handler_base_append_time_step),               public, deferred :: append_time_step
        procedure(output_handler_base_allocate_cell_and_nodal_arrays), public, deferred :: allocate_cell_and_nodal_arrays
        procedure(output_handler_base_append_cell),                    public, deferred :: append_cell
        procedure(output_handler_base_write),                          public, deferred :: write
        procedure(output_handler_base_close),                          public, deferred :: close
    end type

    abstract interface
        subroutine output_handler_base_open(this, dir_path, prefix, parameter_list)
            import output_handler_base_t
            import ParameterList_t
            class(output_handler_base_t),    intent(inout) :: this
            character(len=*),                intent(in)    :: dir_path
            character(len=*),                intent(in)    :: prefix
            type(ParameterList_t), optional, intent(in)    :: parameter_list
        end subroutine

        subroutine output_handler_base_append_time_step(this, value)
            import output_handler_base_t
            import rp
            class(output_handler_base_t), intent(inout) :: this
            real(rp),                     intent(in)    :: value
        end subroutine

        subroutine output_handler_base_write(this)
            import output_handler_base_t
            class(output_handler_base_t), intent(inout) :: this
        end subroutine

        subroutine output_handler_base_allocate_cell_and_nodal_arrays(this)
            import output_handler_base_t
            class(output_handler_base_t), intent(inout) :: this
        end subroutine

        subroutine output_handler_base_append_cell(this, subcell_accessor)
            import output_handler_base_t
            import patch_subcell_accessor_t
            class(output_handler_base_t),   intent(inout) :: this
            type(patch_subcell_accessor_t), intent(in)    :: subcell_accessor
        end subroutine

        subroutine output_handler_base_close(this)
            import output_handler_base_t
            class(output_handler_base_t), intent(inout) :: this
        end subroutine
    end interface

    class(output_handler_fe_iterator_t), allocatable, target, save :: default_output_handler_fe_iterator

public :: output_handler_base_t

contains

!---------------------------------------------------------------------
!< output_handler_BASE_T PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_base_free(this)
    !-----------------------------------------------------------------
    !< Free output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(output_handler_base_t), intent(inout) :: this
        integer(ip)                                 :: i
    !-----------------------------------------------------------------
        if(allocated(this%fields)) then
            do i=1, size(this%fields)
                call this%fields(i)%free()
            enddo
            deallocate(this%fields)
        endif
        if(associated(this%iterator)) then
            call this%iterator%free()
        endif
        nullify(this%iterator)
        nullify(this%fe_space)
        this%number_fields = 0
        this%number_nodes  = 0
        this%number_cells  = 0
    end subroutine output_handler_base_free


    subroutine output_handler_base_set_default_iterator(this, iterator)
    !-----------------------------------------------------------------
    !< Set default output handler
    !-----------------------------------------------------------------
        class(output_handler_base_t),        intent(in) :: this
        class(output_handler_fe_iterator_t), intent(in) :: iterator
        integer                                         :: error
    !-----------------------------------------------------------------
        if (allocated(default_output_handler_fe_iterator)) deallocate(default_output_handler_fe_iterator) 
        allocate(default_output_handler_fe_iterator, mold=iterator, stat=error)
        check(error==0)
    end subroutine output_handler_base_set_default_iterator


    function output_handler_base_get_default_iterator(this) result(iterator)
    !-----------------------------------------------------------------
    !< Return default output_handler_fe_iterator
    !-----------------------------------------------------------------
        class(output_handler_base_t),        intent(in) :: this
        class(output_handler_fe_iterator_t), pointer    :: iterator
    !-----------------------------------------------------------------
        if (.not. allocated(default_output_handler_fe_iterator)) then 
            allocate(output_handler_fe_iterator_t :: default_output_handler_fe_iterator)
        end if
        iterator => default_output_handler_fe_iterator
    end function output_handler_base_get_default_iterator


    subroutine output_handler_base_set_iterator(this, iterator)
    !-----------------------------------------------------------------
    !< Set output handler fe_iterator
    !-----------------------------------------------------------------
        class(output_handler_base_t),        intent(inout) :: this
        class(output_handler_fe_iterator_t), intent(in)    :: iterator
        integer                                            :: error
    !-----------------------------------------------------------------
        if (associated(this%iterator)) deallocate(this%iterator) 
        allocate(this%iterator, mold=iterator, stat=error)
        check(error==0)
    end subroutine output_handler_base_set_iterator


    function output_handler_base_get_number_nodes(this) result(number_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !-----------------------------------------------------------------
        class(output_handler_base_t), intent(in) :: this
        integer(ip)                              :: number_nodes
    !-----------------------------------------------------------------
        number_nodes = this%number_nodes
    end function output_handler_base_get_number_nodes


    function output_handler_base_get_number_cells(this) result(number_cells)
    !-----------------------------------------------------------------
    !< Return the number of cells
    !-----------------------------------------------------------------
        class(output_handler_base_t), intent(in) :: this
        integer(ip)                              :: number_cells
    !-----------------------------------------------------------------
        number_cells = this%number_cells
    end function output_handler_base_get_number_cells


    function output_handler_base_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields
    !-----------------------------------------------------------------
        class(output_handler_base_t), intent(in) :: this
        integer(ip)                              :: number_fields
    !-----------------------------------------------------------------
        number_fields = this%number_fields
    end function output_handler_base_get_number_fields


    function output_handler_base_get_field(this, field_id) result(field)
    !-----------------------------------------------------------------
    !< Return the number of fields
    !-----------------------------------------------------------------
        class(output_handler_base_t),    target, intent(in) :: this
        integer(ip),                             intent(in) :: field_id
        type(output_handler_fe_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(field_id <= this%number_fields)
        field => this%fields(field_id)
    end function output_handler_base_get_field


    subroutine output_handler_base_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(output_handler_base_t),          intent(inout) :: this
        class(serial_fe_space_t), target,      intent(in)    :: fe_space
    !-----------------------------------------------------------------
        this%fe_space => fe_space
    end subroutine output_handler_base_attach_fe_space


    function output_handler_base_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Return a fe_space pointer
    !-----------------------------------------------------------------
        class(output_handler_base_t), intent(in) :: this
        class(serial_fe_space_t), pointer        :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function output_handler_base_get_fe_space


    subroutine output_handler_base_resize_fields_if_needed(this, number_fields)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(output_handler_base_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
        integer(ip)                                       :: current_size
        type(output_handler_fe_field_t), allocatable      :: temp_fe_functions(:)
    !-----------------------------------------------------------------
        if(.not. allocated(this%fields)) then
            allocate(this%fields(fe_functions_initial_size))
        elseif(number_fields > size(this%fields)) then
            current_size = size(this%fields)
            allocate(temp_fe_functions(current_size))
            temp_fe_functions(1:current_size) = this%fields(1:current_size)
            deallocate(this%fields)
            allocate(temp_fe_functions(int(1.5*current_size)))
            this%fields(1:current_size) = temp_fe_functions(1:current_size)
            deallocate(temp_fe_functions)
        endif
    end subroutine output_handler_base_resize_fields_if_needed


    subroutine output_handler_base_add_fe_function(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Add a fe_function to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(output_handler_base_t),       intent(inout) :: this
        type(fe_function_t),                intent(in)    :: fe_function
        integer(ip),                        intent(in)    :: field_id
        character(len=*),                   intent(in)    :: name
        character(len=*), optional,         intent(in)    :: diff_operator
    !-----------------------------------------------------------------
        call this%resize_fields_if_needed(this%number_fields+1)
        this%number_fields = this%number_fields + 1
        call this%fields(this%number_fields)%set(fe_function, field_id, name, diff_operator)
    end subroutine output_handler_base_add_fe_function


    subroutine output_handler_base_fill_data(this)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the output_handler_base_t derived type
    !-----------------------------------------------------------------
        class(output_handler_base_t),     intent(inout) :: this
        type(fe_accessor_t)                             :: fe
        type(output_handler_cell_fe_function_t)         :: output_handler_cell_function
        type(output_handler_patch_t)                    :: patch
        type(patch_subcell_iterator_t)                  :: subcell_iterator
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))

        ! Ouput_handler_FE_iterator 
        if(.not. associated (this%iterator)) then
            this%iterator => this%get_default_iterator()
        endif
        call this%iterator%set_fe_iterator(this%fe_space%create_fe_iterator())

        ! Create Output Cell Handler and allocate patch fields
        call output_handler_cell_function%create(this%fe_space, this%iterator, this%number_fields, this%fields(1:this%number_fields))

        ! Allocate geometry and connectivity arrays
        this%number_nodes = output_handler_cell_function%get_number_nodes()
        this%number_cells = output_handler_cell_function%get_number_cells()
        call this%allocate_cell_and_nodal_arrays()

        call patch%create(this%number_fields)
        ! Translate coordinates and connectivities to VTK format for every subcell
        call this%iterator%init()
        do while ( .not. this%iterator%has_finished())
            ! Get Finite element
            call this%iterator%current(fe)
            if ( fe%is_local() ) then
                call output_handler_cell_function%fill_patch(fe, this%number_fields, this%fields(1:this%number_fields), patch)
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

end module output_handler_base_names

