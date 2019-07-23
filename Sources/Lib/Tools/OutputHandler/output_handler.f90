!   Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
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
!---------------------------------------------------------------------
!*Author: Víctor Sande
! Date: 2016-11-28
! Version: 0.0.1
! Category: IO
!
!---------------------------------------------------------------------
!### Software subsystem implementing the OO IO layer.
!
! Contains the following public entities:
! [[output_handler_names(module)]]
! 
! @note Look at [[base_output_handler_names(module)]] to see the
! *state transition diagram*
!---------------------------------------------------------------------
module output_handler_names
!---------------------------------------------------------------------
!*Author: Víctor Sande
! Date: 2016-11-28
! Version: 0.0.1
! Category: IO
!
!---------------------------------------------------------------------
!### Software subsystem implementing the OO IO layer.
!
! Contains the following public entities:
! [[output_handler_t(type)]], 
! [[output_handler_prototype_reset(subroutine)]], 
! [[output_handler_prototype_free(subroutine)]], 
! [[output_handler_names(module):VTK(variable)]] and  
! [[output_handler_names(module):XH5(variable)]]
! 
! @note Look at [[base_output_handler_t(type)]] to see the
! *state transition diagram*
!---------------------------------------------------------------------

USE FPL
USE types_names
USE fe_space_names,              only: serial_fe_space_t, fe_cell_iterator_t, fe_function_t
USE std_vector_real_rp_names
USE output_handler_field_generator_names
USE output_handler_parameters_names
USE base_output_handler_names
USE vtk_output_handler_names
USE xh5_output_handler_names

implicit none
#include "debug.i90"
private

    type :: output_handler_t
    !-----------------------------------------------------------------
    !*Author: Víctor Sande
    ! Date: 2016-11-28
    ! Version: 0.0.1
    ! Category: IO
    ! 
    !-----------------------------------------------------------------
    !### Reading/writting mesh files and results
    !
    ! The main task of this container is to host the chosen IO *strategy*
    ! as the current **state**.
    !
    ! Container of a polymorphic entity of type [[base_output_handler_t(type)]].
    ! 
    !-----------------------------------------------------------------
    !### Usage
    !```fortran
    !...
    !    type(output_handler_t) :: oh
    !...
    !    call oh%create()
    !    call oh%attach_fe_space(fe_space)
    !    call oh%add_fe_function(fe_function, field_id, 'field_name')
    !    call oh%add_fe_function(fe_function, field_id, 'grad_field_name', grad_diff_operator)
    !    call oh%open(Path, Prefix)
    !    call oh%write()
    !    call oh%close()
    !    call oh%free()
    !...
    !```
    !-----------------------------------------------------------------
    private
        character(len=:),             allocatable :: format
        character(len=:),             allocatable :: prefix
        character(len=:),             allocatable :: dir_path
        class(base_output_handler_t), allocatable :: state
    contains
    private
        procedure, non_overridable, public :: create                         => output_handler_create
        procedure, non_overridable, public :: attach_fe_space                => output_handler_attach_fe_space
        procedure, non_overridable, public :: set_create_fe_cell_iterator    => output_handler_set_create_fe_cell_iterator
        procedure, non_overridable, public :: set_free_fe_cell_iterator      => output_handler_set_free_fe_cell_iterator
        procedure, non_overridable, public :: add_fe_function                => output_handler_add_fe_function
        procedure, non_overridable, public :: add_field_generator            => output_handler_add_field_generator
        procedure, non_overridable, public :: add_cell_vector                => output_handler_add_cell_vector
        procedure, non_overridable, public :: update_cell_vector             => output_handler_update_cell_vector
        procedure, non_overridable         :: open_noargs                    => output_handler_open_noargs
        procedure, non_overridable         :: open_dir_path_prefix           => output_handler_open_dir_path_prefix
        procedure, non_overridable, public :: append_time_step               => output_handler_append_time_step
        procedure, non_overridable, public :: write                          => output_handler_write
        procedure, non_overridable, public :: close                          => output_handler_close
        procedure, non_overridable, public :: free                           => output_handler_free
        generic,                    public :: open                           => open_dir_path_prefix, open_noargs
    end type

    public :: output_handler_t

contains

    subroutine output_handler_create(this, parameterlist)
    !-----------------------------------------------------------------
    !< Create the output handler **state** given a parameter list
    !-----------------------------------------------------------------
        class(output_handler_t), intent(inout) :: this
        type(ParameterList_t),   intent(in)    :: parameterlist
        integer(ip)                            :: error
    !-----------------------------------------------------------------
        this%dir_path = output_handler_default_dir_path
        this%prefix   = output_handler_default_prefix
        this%format   = output_handler_default_format

        ! Mandatory parameters
        ! Get fomat value from parameterlist
        assert(parameterlist%isPresent(output_handler_format_key))
        assert(parameterlist%isAssignable(output_handler_format_key, 'string'))
        error = parameterlist%getasstring(key = output_handler_format_key, string = this%format)
        check(error==0)

        ! Get dir_path value from parameterlist
        assert(parameterlist%isPresent(output_handler_dir_path_key))
        assert(parameterlist%isAssignable(output_handler_dir_path_key, 'string'))
        error = parameterlist%GetAsString(Key=output_handler_dir_path_key, String=this%dir_path)
        assert(error == 0)

        ! Get prefix value from parameterlist
        assert(parameterlist%isPresent(output_handler_prefix_key))
        assert(parameterlist%isAssignable(output_handler_prefix_key, 'string'))
        error = parameterlist%GetAsString(Key=output_handler_prefix_key, String=this%prefix)
        assert(error == 0)

        ! Allocate state according to the given format
        call this%free()
        select case (this%format)
            case (VTK)
                allocate(vtk_output_handler_t :: this%state)
            case (XH5)
                allocate(xh5_output_handler_t :: this%state)
            case DEFAULT
                check(.false.)
        end select

        ! Call concrete implementation of create procedure
        call this%state%create(parameterlist)
    end subroutine output_handler_create


    subroutine output_handler_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a **fe_space** to work with. 
    !< The user is responsible for ensuring compatibility 
    !< between the [[serial_fe_space_t(type)]] attached and the 
    !< [[fe_function_t(type)]] and **cell_vector_t** added.
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
        class(serial_fe_space_t), target, intent(in)    :: fe_space
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%attach_fe_space(fe_space)
    end subroutine output_handler_attach_fe_space
    
    subroutine output_handler_set_create_fe_cell_iterator(this, create_fe_cell_iterator)
      class(output_handler_t), intent(inout) :: this
       procedure(create_fe_cell_iterator_interface) :: create_fe_cell_iterator
       call this%state%set_create_fe_cell_iterator(create_fe_cell_iterator)
    end subroutine output_handler_set_create_fe_cell_iterator

    subroutine output_handler_set_free_fe_cell_iterator(this, free_fe_cell_iterator)
      class(output_handler_t), intent(inout) :: this
      procedure(free_fe_cell_iterator_interface) :: free_fe_cell_iterator
      call this%state%set_free_fe_cell_iterator(free_fe_cell_iterator)
    end subroutine output_handler_set_free_fe_cell_iterator

    subroutine output_handler_add_fe_function(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Add [[fe_function_t(type)]].
    !<
    !< The user is responsible for ensuring compatibility 
    !< between the [[serial_fe_space_t(type)]] attached and the 
    !< [[fe_function_t(type)]] added.
    !-----------------------------------------------------------------
        class(output_handler_t),    intent(inout) :: this
        type(fe_function_t),        intent(in)    :: fe_function
        integer(ip),                intent(in)    :: field_id
        character(len=*),           intent(in)    :: name
        character(len=*), optional, intent(in)    :: diff_operator
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%add_fe_function(fe_function, field_id, name, diff_operator)
    end subroutine output_handler_add_fe_function
    
    subroutine output_handler_add_field_generator(this, name, output_handler_field_generator)
        class(output_handler_t)                ,    intent(inout) :: this
        character(len=*)                       ,    intent(in)    :: name
        class(output_handler_field_generator_t),    intent(in)    :: output_handler_field_generator
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%add_field_generator(name,output_handler_field_generator)
    end subroutine output_handler_add_field_generator

    subroutine output_handler_add_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Adds a new **cell_vector** entry to output_handler_t identified 
    !< by **name**. The user is responsible for ensuring compatibility 
    !< between the [[serial_fe_space_t(type)]] attached and the newly
    !< added **cell_vector**.
    !-----------------------------------------------------------------
        class(output_handler_t),    intent(inout) :: this
        type(std_vector_real_rp_t), intent(in)    :: cell_vector
        character(len=*),           intent(in)    :: name
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%add_cell_vector(cell_vector, name)
    end subroutine output_handler_add_cell_vector
    
    subroutine output_handler_update_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Update the pointer to which an existing cell_vector entry already
    !< registered in output_handler_t is pointing to s.t. it points to
    !< **cell_vector**. The provided **name** is used to locate the cell_vector
    !< within output_handler_t. If **name** is not found, a warning message is issued. 
    !< The user is responsible for ensuring compatibility between the [[serial_fe_space_t(type)]] 
    !< attached and the updated **cell_vector**.
    !-----------------------------------------------------------------
        class(output_handler_t),    intent(inout) :: this
        type(std_vector_real_rp_t), intent(in)    :: cell_vector
        character(len=*),           intent(in)    :: name
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%update_cell_vector(cell_vector, name)
    end subroutine output_handler_update_cell_vector


    subroutine output_handler_open_noargs(this)
    !-----------------------------------------------------------------
    !< Open procedure. 
    !< First procedure to be called after attaching 
    !< the [[serial_fe_space_t(type)]] and adding [[fe_function_t(type)]] and 
    !< **cell_vector** and before starting to write data to disk.
    !< This procedure must be called only once and out-of-the-loop
    !< in transient simulations
    !-----------------------------------------------------------------
        class(output_handler_t),         intent(inout) :: this
    !-----------------------------------------------------------------
        call this%open(this%dir_path, this%prefix)
    end subroutine output_handler_open_noargs


    subroutine output_handler_open_dir_path_prefix(this, dir_path, prefix)
    !-----------------------------------------------------------------
    !< Open procedure. 
    !< First procedure to be called after attaching 
    !< the [[serial_fe_space_t(type)]] and adding [[fe_function_t(type)]] and 
    !< **cell_vector** and before starting to write data to disk.
    !< This procedure must be called only once and out-of-the-loop
    !< in transient simulations
    !-----------------------------------------------------------------
        class(output_handler_t),         intent(inout) :: this
        character(len=*),                intent(in)    :: dir_path
        character(len=*),                intent(in)    :: prefix
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%open(dir_path, prefix)
    end subroutine output_handler_open_dir_path_prefix


    subroutine output_handler_append_time_step(this, value)
    !-----------------------------------------------------------------
    !< Set the current time step value and append it to a vector.
    !-----------------------------------------------------------------
        class(output_handler_t), intent(inout) :: this
        real(rp),                intent(in)    :: value
    !-----------------------------------------------------------------
        assert(allocated(this%state))    
        call this%state%append_time_step(value)
    end  subroutine output_handler_append_time_step


    subroutine output_handler_write(this)
    !-----------------------------------------------------------------
    !< Converts and writes to a mesh file the attached/added data.
    !< [[serial_fe_space_t(type)]] to mesh file, 
    !< the added [[fe_function_T(type)]] to nodal fields and 
    !< the added **cell_vector** as cell fields.
    !< This procedure must be called one time each time step
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%write()
    end subroutine output_handler_write


    subroutine output_handler_close(this)
    !-----------------------------------------------------------------
    !< Close procedure. Last procedure to call when finished the data
    !< writing to write.
    !< This procedure must be called only once and out-of-the-loop
    !< in transient simulations
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        assert(allocated(this%state))
        call this%state%close()
    end subroutine output_handler_close


    subroutine output_handler_free(this)
    !-----------------------------------------------------------------
    !< Free and deallocate the **State**
    !-----------------------------------------------------------------
        class(output_handler_t),          intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%state)) then
            call this%State%Free()
            deallocate(this%state)
        endif
    end subroutine output_handler_free

end module output_handler_names
