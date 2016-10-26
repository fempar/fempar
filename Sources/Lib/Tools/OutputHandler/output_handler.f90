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

module output_handler_names

USE FPL
USE types_names
USE fe_space_names,              only: serial_fe_space_t, fe_iterator_t, fe_accessor_t
USE fe_function_names,           only: fe_function_t
USE base_output_handler_names
USE output_handler_fe_iterator_names
USE vtk_output_handler_names
USE xh5_output_handler_names

implicit none
#include "debug.i90"
private

    character(len=*), parameter :: VTK = 'VTK'
    character(len=*), parameter :: XH5 = 'XH5'

    type :: output_handler_t
    private
        class(base_output_handler_t), pointer :: state => NULL()
    contains
        procedure, non_overridable         ::                                   output_handler_create
        procedure, non_overridable         ::                                   output_handler_create_string
        procedure, non_overridable         :: get_default_output_handler     => output_handler_get_default_output_handler
        procedure, non_overridable, public :: set_default_output_handler     => output_handler_set_default_output_handler
        procedure, non_overridable, public :: attach_fe_space                => output_handler_attach_fe_space
        procedure, non_overridable, public :: set_iterator                   => output_handler_set_iterator
        procedure, non_overridable, public :: add_fe_function                => output_handler_add_fe_function
        procedure, non_overridable, public :: add_cell_vector                => output_handler_add_cell_vector
        procedure, non_overridable, public :: open                           => output_handler_open
        procedure, non_overridable, public :: append_time_step               => output_handler_append_time_step
        procedure, non_overridable, public :: write                          => output_handler_write
        procedure, non_overridable, public :: close                          => output_handler_close
        procedure, non_overridable, public :: free                           => output_handler_free
        generic,                    public :: create                         => output_handler_create, &
                                                                                output_handler_create_string
    end type

    class(base_output_handler_t), allocatable, target, save :: default_output_handler

public :: output_handler_t
public :: VTK
public :: XH5

contains

    subroutine output_handler_create(this)
    !-----------------------------------------------------------------
    !< Create the state output handler by means of the default output handler
    !-----------------------------------------------------------------
        class(output_handler_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state => this%get_default_output_handler()
    end subroutine output_handler_create


    subroutine output_handler_create_string(this, descriptor)
    !-----------------------------------------------------------------
    !< Create the state output handler given its descriptor string
    !-----------------------------------------------------------------
        class(output_handler_t), intent(inout) :: this
        character(len=*),        intent(in)    :: descriptor
    !-----------------------------------------------------------------
        if(associated(this%state)) then
            call this%state%free()
            deallocate(this%state)
        endif
        select case (descriptor)
            case (VTK)
                allocate(vtk_output_handler_t :: this%state)
            case (XH5)
                allocate(xh5_output_handler_t :: this%state)
            case DEFAULT
                check(.false.)
        end select
    end subroutine output_handler_create_string


    subroutine output_handler_set_default_output_handler(this, output_handler)
    !-----------------------------------------------------------------
    !< Set default output handler
    !-----------------------------------------------------------------
        class(output_handler_t),      intent(in) :: this
        class(base_output_handler_t), intent(in) :: output_handler
        integer                                  :: error
    !-----------------------------------------------------------------
        if (allocated(default_output_handler)) deallocate(default_output_handler) 
        allocate(default_output_handler, mold=output_handler, stat=error)
        check(error==0)
    end subroutine output_handler_set_default_output_handler


    function output_handler_get_default_output_handler(this) result(output_handler)
    !-----------------------------------------------------------------
    !< Return default output handler
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
        class(base_output_handler_t), pointer           :: output_handler
    !-----------------------------------------------------------------
        if (.not. allocated(default_output_handler)) then 
#ifdef ENABLE_HDF5
            allocate(xh5_output_handler_t :: default_output_handler)
#else
            allocate(vtk_output_handler_t :: default_output_handler)
#endif
        end if
        output_handler => default_output_handler
    end function output_handler_get_default_output_handler


    subroutine output_handler_set_iterator(this, iterator)
    !-----------------------------------------------------------------
    !< Set output handler fe_iterator
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
        class(output_handler_fe_iterator_t), intent(in) :: iterator
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%set_iterator(iterator)
    end subroutine output_handler_set_iterator

    subroutine output_handler_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a fe_space 
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
        class(serial_fe_space_t), target, intent(in)    :: fe_space
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%attach_fe_space(fe_space)
    end subroutine output_handler_attach_fe_space


    subroutine output_handler_add_fe_function(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Add fe function
    !-----------------------------------------------------------------
        class(output_handler_t),    intent(inout) :: this
        type(fe_function_t),        intent(in)    :: fe_function
        integer(ip),                intent(in)    :: field_id
        character(len=*),           intent(in)    :: name
        character(len=*), optional, intent(in)    :: diff_operator
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%add_fe_function(fe_function, field_id, name, diff_operator)
    end subroutine output_handler_add_fe_function


    subroutine output_handler_add_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Add fe function
    !-----------------------------------------------------------------
        class(output_handler_t),    intent(inout) :: this
        real(rp), allocatable,      intent(in)    :: cell_vector(:)
        character(len=*),           intent(in)    :: name
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%add_cell_vector(cell_vector, name)
    end subroutine output_handler_add_cell_vector


    subroutine output_handler_open(this, dir_path, prefix, parameter_list)
    !-----------------------------------------------------------------
    !< Open procedure
    !-----------------------------------------------------------------
        class(output_handler_t),         intent(inout) :: this
        character(len=*),                intent(in)    :: dir_path
        character(len=*),                intent(in)    :: prefix
        type(ParameterList_t), optional, intent(in)    :: parameter_list
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%open(dir_path, prefix, parameter_list)
    end subroutine output_handler_open


    subroutine output_handler_append_time_step(this, value)
    !-----------------------------------------------------------------
    !< Open procedure
    !-----------------------------------------------------------------
        class(output_handler_t), intent(inout) :: this
        real(rp),                intent(in)    :: value
    !-----------------------------------------------------------------
        assert(associated(this%state))    
        call this%state%append_time_step(value)
    end  subroutine output_handler_append_time_step


    subroutine output_handler_write(this)
    !-----------------------------------------------------------------
    !< Write procedure
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%write()
    end subroutine output_handler_write


    subroutine output_handler_close(this)
    !-----------------------------------------------------------------
    !< Close procedure
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%state))
        call this%state%close()
    end subroutine output_handler_close


    subroutine output_handler_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        if(associated(this%state)) then
            call this%State%Free()
            deallocate(this%state)
        endif
    end subroutine output_handler_free

end module output_handler_names
