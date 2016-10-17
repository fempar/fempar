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

module data_out_base_names

USE types_names
USE fe_space_names,              only: serial_fe_space_t, fe_iterator_t, fe_accessor_t
USE fe_function_names,           only: fe_function_t
USE data_out_patch_names
USE data_out_cell_fe_function_names

implicit none
#include "debug.i90"
private

    integer(ip), parameter :: fe_functions_initial_size = 10

    type :: data_out_fe_function_t
    private
        character(len=:),  allocatable :: name
        integer(ip)                    :: field_id = 0
        type(fe_function_t), pointer   :: fe_function => NULL()
    contains
        procedure         ::                    data_out_fe_function_assign
        procedure, public :: set             => data_out_fe_function_set
        procedure, public :: get_name        => data_out_fe_function_get_name
        procedure, public :: get_field_id    => data_out_fe_function_get_field_id
        procedure, public :: get_fe_function => data_out_fe_function_get_fe_function
        procedure, public :: free            => data_out_fe_function_free
        generic,   public :: assignment(=)   => data_out_fe_function_assign
    end type

    type :: data_out_base_t
    private
        class(serial_fe_space_t),     pointer     :: fe_space => NULL()
        type(data_out_fe_function_t), allocatable :: fe_functions(:)
        integer(ip)                               :: number_fe_functions = 0
    contains
        procedure, public :: free                          => data_out_base_free
        procedure, public :: attach_fe_space               => data_out_base_attach_fe_space
        procedure         :: resize_fe_functions_if_needed => data_out_base_resize_fe_functions_if_needed
        procedure, public :: add_fe_function               => data_out_base_add_fe_function
        procedure, public :: fill_data                     => data_out_base_fill_data
!        procedure(data_out_base_write),       public, deferred :: write
!        procedure(data_out_base_append_cell), public, deferred :: append_cell
    end type

    abstract interface
        subroutine data_out_base_write(this)
            import data_out_base_t
            class(data_out_base_t), intent(in) :: this
        end subroutine

        subroutine data_out_base_append_cell(this)
            import data_out_base_t
            class(data_out_base_t), intent(in) :: this
        end subroutine
    end interface

public :: data_out_base_t

contains

!---------------------------------------------------------------------
!< DATA_OUT_FE_FUNCTION_T PROCEDURES
!---------------------------------------------------------------------

    subroutine data_out_fe_function_free(this)
    !-----------------------------------------------------------------
    !< Free data_out_fe_function_t derived type
    !-----------------------------------------------------------------
        class(data_out_fe_function_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%name)) deallocate(this%name)
        nullify(this%fe_function)
        this%field_id = 0
    end subroutine data_out_fe_function_free


    subroutine data_out_fe_function_assign(this, data_out_fe_function)
    !----------------------------------------------------------------- 
    !< Assign operator overloading for data_out_fe_function_t
    !----------------------------------------------------------------- 
        class(data_out_fe_function_t), intent(inout) :: this
        type(data_out_fe_function_t),  intent(in)    :: data_out_fe_function
    !----------------------------------------------------------------- 
        call this%free()
        call this%set(data_out_fe_function%get_fe_function(), data_out_fe_function%get_field_id(), data_out_fe_function%get_name())
    end subroutine data_out_fe_function_assign


    subroutine data_out_fe_function_set(this, fe_function, field_id, name)
    !-----------------------------------------------------------------
    !< Associate a fe_function with a field name
    !-----------------------------------------------------------------
        class(data_out_fe_function_t), intent(inout) :: this
        type(fe_function_t), target,   intent(in)    :: fe_function
        integer(ip),                   intent(in)    :: field_id
        character(len=*),              intent(in)    :: name
    !-----------------------------------------------------------------
        call this%free()
        this%fe_function => fe_function
        this%field_id    = field_id
        this%name        = name
    end subroutine data_out_fe_function_set


    function data_out_fe_function_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(data_out_fe_function_t), intent(in) :: this
        character(len=:), allocatable             :: name
    !-----------------------------------------------------------------
        assert(allocated(this%name))
        name = this%name
    end function data_out_fe_function_get_name


    function data_out_fe_function_get_field_id(this) result(field_id)
    !-----------------------------------------------------------------
    !< Return the field_id of a field associated with a fe_function
    !-----------------------------------------------------------------
        class(data_out_fe_function_t), intent(in) :: this
        integer(ip)                               :: field_id
    !-----------------------------------------------------------------
        assert(this%field_id /= 0)
        field_id = this%field_id
    end function data_out_fe_function_get_field_id


    function data_out_fe_function_get_fe_function(this) result(fe_function)
    !-----------------------------------------------------------------
    !< Return a fe_function pointer
    !-----------------------------------------------------------------
        class(data_out_fe_function_t), target, intent(in) :: this
        type(fe_function_t), pointer                      :: fe_function
    !-----------------------------------------------------------------
        assert(associated(this%fe_function))
        fe_function => this%fe_function
    end function data_out_fe_function_get_fe_function

!---------------------------------------------------------------------
!< DATA_OUT_BASE_T PROCEDURES
!---------------------------------------------------------------------

    subroutine data_out_base_free(this)
    !-----------------------------------------------------------------
    !< Free data_out_base_t derived type
    !-----------------------------------------------------------------
        class(data_out_base_t), intent(inout) :: this
        integer(ip)                           :: i
    !-----------------------------------------------------------------
        if(allocated(this%fe_functions)) then
            do i=1, size(this%fe_functions)
                call this%fe_functions(i)%free()
            enddo
            deallocate(this%fe_functions)
        endif
        this%number_fe_functions = 0
        nullify(this%fe_space)
    end subroutine data_out_base_free


    subroutine data_out_base_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the data_out_base_t derived type
    !-----------------------------------------------------------------
        class(data_out_base_t),          intent(inout) :: this
        class(serial_fe_space_t), target, intent(in)   :: fe_space
    !-----------------------------------------------------------------
        this%fe_space => fe_space
    end subroutine data_out_base_attach_fe_space


    subroutine data_out_base_resize_fe_functions_if_needed(this, number_fe_functions)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the data_out_base_t derived type
    !-----------------------------------------------------------------
        class(data_out_base_t),       intent(inout) :: this
        integer(ip),                  intent(in)    :: number_fe_functions
        integer(ip)                                 :: current_size
        type(data_out_fe_function_t), allocatable   :: temp_fe_functions(:)
    !-----------------------------------------------------------------
        if(.not. allocated(this%fe_functions)) then
            allocate(this%fe_functions(fe_functions_initial_size))
        elseif(number_fe_functions > size(this%fe_functions)) then
            current_size = size(this%fe_functions)
            allocate(temp_fe_functions(current_size))
            temp_fe_functions(1:current_size) = this%fe_functions(1:current_size)
            deallocate(this%fe_functions)
            allocate(temp_fe_functions(int(1.5*current_size)))
            this%fe_functions(1:current_size) = temp_fe_functions(1:current_size)
            deallocate(temp_fe_functions)
        endif
    end subroutine data_out_base_resize_fe_functions_if_needed


    subroutine data_out_base_add_fe_function(this, fe_function, field_id, name)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the data_out_base_t derived type
    !-----------------------------------------------------------------
        class(data_out_base_t),       intent(inout) :: this
        type(fe_function_t),          intent(in)    :: fe_function
        integer(ip),                  intent(in)    :: field_id
        character(len=*),             intent(in)    :: name
    !-----------------------------------------------------------------
        call this%resize_fe_functions_if_needed(this%number_fe_functions+1)
        this%number_fe_functions = this%number_fe_functions + 1
        call this%fe_functions(this%number_fe_functions)%set(fe_function, field_id, name)
    end subroutine data_out_base_add_fe_function


    subroutine data_out_base_fill_data(this)
    !-----------------------------------------------------------------
    !< Attach a fe_space to the data_out_base_t derived type
    !-----------------------------------------------------------------
        class(data_out_base_t),     intent(inout) :: this
        type(fe_iterator_t)                       :: fe_iterator
        type(fe_accessor_t)                       :: fe
        type(data_out_cell_fe_function_t)         :: data_out_cell_function
        type(data_out_patch_t)                    :: data_out_patch
        type(data_out_patch_subcell_iterator_t)   :: subcell_iterator
        real(rp), allocatable                     :: coordinates(:,:)
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))

        ! Create FE iterator 
        fe_iterator             = this%fe_space%create_fe_iterator()

        ! Allocate VTK geometry and connectivity data
!        call this%allocate_elemental_arrays()
!        call this%allocate_nodal_arrays()
!        call this%initialize_coordinates()

        ! Create Output Cell Handler
        call data_out_cell_function%create(this%fe_space)
        allocate(data_out_patch%fields(1))
        data_out_patch%fields(1)%fe_function => this%fe_functions(1)%get_fe_function()
        data_out_patch%fields(1)%field_id = 1
        data_out_patch%fields(1)%name = 'solution'
        ! Translate coordinates and connectivities to VTK format for every subcell
        do while ( .not. fe_iterator%has_finished())
            ! Get Finite element
            call fe_iterator%current(fe)
            if ( fe%is_local() ) then
                call data_out_cell_function%build_patch(fe, data_out_patch)
                subcell_iterator = data_out_patch%get_subcells_iterator()
!               ! Fill data
                do while(.not. subcell_iterator%has_finished())
                    call subcell_iterator%get_coordinates(coordinates)
print*, coordinates
                    call subcell_iterator%next()
!                   this%append_cell(SUBCELL)
                enddo
            endif
            call fe_iterator%next()
        end do

        call data_out_patch%free()
        call data_out_cell_function%free()
        call fe_iterator%free()
        call fe%free()
    end subroutine data_out_base_fill_data

end module data_out_base_names

