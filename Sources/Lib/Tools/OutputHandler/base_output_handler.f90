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
!---------------------------------------------------------------------
!*Author: Víctor Sande
! Date: 2016-11-28
! Version: 0.0.1
! Category: IO
!
!---------------------------------------------------------------------
!### Abstract base class implementing IO operations.
!
! Contains the following public entities: 
! [[base_output_handler_names(module)]]
!---------------------------------------------------------------------

module base_output_handler_names
!---------------------------------------------------------------------
!*Author: Víctor Sande
! Date: 2016-11-28
! Version: 0.0.1
! Category: IO
!
!---------------------------------------------------------------------
!### Abstract base class implementing IO operations.
!
! Contains the following public entities: 
! [[base_output_handler_t(type)]]
!---------------------------------------------------------------------

USE FPL
USE types_names
USE memor_names
USE reference_fe_names,          only: field_type_scalar, field_type_vector, field_type_tensor, field_type_symmetric_tensor
USE fe_space_names,              only: serial_fe_space_t, fe_iterator_t
USE fe_function_names,           only: fe_function_t
USE output_handler_fe_field_names
USE output_handler_patch_names
USE output_handler_cell_fe_function_names
USE output_handler_parameters_names

implicit none
#include "debug.i90"
private

    integer(ip), parameter :: cell_node_field_array_size = 10

    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_INIT   = 0
    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_OPEN   = 1
    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_FILL   = 2

    type, abstract :: base_output_handler_t
    !-----------------------------------------------------------------
    !*Author: Víctor Sande
    ! Date: 2016-11-28
    ! Version: 0.0.1
    ! Category: IO
    ! 
    !-----------------------------------------------------------------
    !### Base class to perform mesh and fields IO.
    !
    ! Abstract entity [[base_output_handler_t(type)]] defines the common
    ! variables and procedures for all classes extending it.
    ! 
    ! ----------------------------------------------------------------
    ! **State transition diagram for [[base_output_handler_t(type)]] **
    ! --------------------------------------------------------------------------------------------
    ! Input State         | Action                                           | Output State 
    ! :------------------:|:------------------------------------------------:|:------------------:
    ! Init                | [[base_output_handler_t(type):Free(bound)]]      | Init
    ! Init                | [[base_output_handler_t(type):Open(bound)]]      | Open
    ! Open                | [[base_output_handler_t(type):Free(bound)]]      | Init
    ! Open                | [[base_output_handler_t(type):Fill_data(bound)]] | Fill
    ! Open                | [[base_output_handler_t(type):Close(bound)]]     | Init
    ! Fill                | [[base_output_handler_t(type):Free(bound)]]      | Init
    ! Fill                | [[base_output_handler_t(type):Close(bound)]]     | Init
    ! 
    ! @note it is desirable that the state management occurs only
    !       inside this class to get a cleaner implementation
    !       of the son classes
    ! 
    ! @note
    ! All setters must be called before OPEN.
    ! All getters (except get_fe_space) must be called after OPEN.
    ! FILLED occurs when metadata is filled from [[output_handler_cell_fe_function_t(type)]].
    !-----------------------------------------------------------------
    private
        class(serial_fe_space_t),            pointer             :: fe_space              => NULL()
        type(output_handler_cell_fe_function_t)                  :: ohcff
        type(output_handler_fe_field_t),    allocatable          :: fe_fields(:)
        type(output_handler_cell_vector_t), allocatable          :: cell_vectors(:)
        procedure(create_fe_iterator_interface), nopass, pointer :: create_fe_iterator    => NULL()
        procedure(free_fe_iterator_interface)  , nopass, pointer :: free_fe_iterator      => NULL()
        integer(ip)                                              :: state                 = BASE_OUTPUT_HANDLER_STATE_INIT
        integer(ip)                                              :: number_fields         = 0
        integer(ip)                                              :: number_cell_vectors   = 0
    contains
    private
        procedure, non_overridable, public :: free                          => base_output_handler_free
        procedure, non_overridable, public :: get_number_nodes              => base_output_handler_get_number_nodes
        procedure, non_overridable, public :: get_number_cells              => base_output_handler_get_number_cells
        procedure, non_overridable, public :: has_mixed_cell_topologies     => base_output_handler_has_mixed_cell_topologies
        procedure, non_overridable, public :: get_number_dimensions         => base_output_handler_get_number_dimensions
        procedure, non_overridable, public :: get_number_fields             => base_output_handler_get_number_fields
        procedure, non_overridable, public :: get_number_cell_vectors       => base_output_handler_get_number_cell_vectors
        procedure, non_overridable, public :: get_fe_field                  => base_output_handler_get_fe_field
        procedure, non_overridable, public :: get_cell_vector               => base_output_handler_get_cell_vector
        procedure, non_overridable, public :: get_fe_space                  => base_output_handler_get_fe_space
        procedure, non_overridable         :: resize_fe_fields_if_needed    => base_output_handler_resize_fe_fields_if_needed
        procedure, non_overridable         :: resize_cell_vectors_if_needed => base_output_handler_resize_cell_vectors_if_needed
        procedure, non_overridable, public :: attach_fe_space               => base_output_handler_attach_fe_space
        procedure, non_overridable, public :: set_create_fe_iterator        => base_output_handler_set_create_fe_iterator
        procedure, non_overridable, public :: set_free_fe_iterator          => base_output_handler_set_free_fe_iterator
        procedure, non_overridable, private:: create_fe_iterator_wrapper    => base_output_handler_create_fe_iterator_wrapper
        procedure, non_overridable, private:: free_fe_iterator_wrapper      => base_output_handler_free_fe_iterator_wrapper
        procedure, non_overridable, public :: add_fe_function               => base_output_handler_add_fe_function
        procedure, non_overridable, public :: add_cell_vector               => base_output_handler_add_cell_vector
        procedure, non_overridable, public :: fill_data                     => base_output_handler_fill_data
        procedure, non_overridable, public :: open                          => base_output_handler_open
        procedure, non_overridable, public :: close                         => base_output_handler_close
        procedure(base_output_handler_open_body),                      public, deferred :: open_body
        procedure(base_output_handler_append_time_step),               public, deferred :: append_time_step
        procedure(base_output_handler_allocate_cell_and_nodal_arrays),         deferred :: allocate_cell_and_nodal_arrays
        procedure(base_output_handler_append_cell),                            deferred :: append_cell
        procedure(base_output_handler_write),                          public, deferred :: write
        procedure(base_output_handler_close_body),                     public, deferred :: close_body
        procedure(base_output_handler_free_body),                              deferred :: free_body
    end type

    interface
      subroutine create_fe_iterator_interface(fe)
        import :: fe_iterator_t
        class(fe_iterator_t), allocatable, intent(inout) :: fe
      end subroutine create_fe_iterator_interface
      
      subroutine free_fe_iterator_interface(fe)
        import :: fe_iterator_t
        class(fe_iterator_t), allocatable, intent(inout) :: fe
      end subroutine free_fe_iterator_interface
    end interface
    
    
    abstract interface
        subroutine base_output_handler_open_body(this, dir_path, prefix, parameter_list)
            import base_output_handler_t
            import ParameterList_t
            class(base_output_handler_t),    intent(inout) :: this
            character(len=*),                intent(in)    :: dir_path
            character(len=*),                intent(in)    :: prefix
            type(ParameterList_t), optional, intent(in)    :: parameter_list
        end subroutine

        subroutine base_output_handler_append_time_step(this, value)
            import base_output_handler_t
            import rp
            class(base_output_handler_t), intent(inout) :: this
            real(rp),                     intent(in)    :: value
        end subroutine

        subroutine base_output_handler_write(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine base_output_handler_allocate_cell_and_nodal_arrays(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine base_output_handler_append_cell(this, subcell_accessor)
            import base_output_handler_t
            import patch_subcell_accessor_t
            class(base_output_handler_t),   intent(inout) :: this
            type(patch_subcell_accessor_t), intent(in)    :: subcell_accessor
        end subroutine

        subroutine base_output_handler_close_body(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine

        subroutine base_output_handler_free_body(this)
            import base_output_handler_t
            class(base_output_handler_t), intent(inout) :: this
        end subroutine
    end interface

public :: base_output_handler_t, create_fe_iterator_interface, free_fe_iterator_interface 

contains

!---------------------------------------------------------------------
! base_output_handler_t PROCEDURES
!---------------------------------------------------------------------

    subroutine base_output_handler_free(this)
    !-----------------------------------------------------------------
    !< Free [[base_output_handler_t(type)]] derived type
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
        call this%ohcff%free()
        nullify(this%fe_space)
        nullify(this%create_fe_iterator)
        nullify(this%free_fe_iterator)
        this%number_cell_vectors   = 0
        this%number_fields         = 0
        this%state                 = BASE_OUTPUT_HANDLER_STATE_INIT
    end subroutine base_output_handler_free

    function base_output_handler_get_number_nodes(this) result(number_nodes)
    !-----------------------------------------------------------------
    !< Return the number of visualization nodes
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_nodes
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        number_nodes = this%ohcff%get_number_nodes()
    end function base_output_handler_get_number_nodes


    function base_output_handler_get_number_cells(this) result(number_cells)
    !-----------------------------------------------------------------
    !< Return the number of visualization cells
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_cells
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        number_cells = this%ohcff%get_number_cells()
    end function base_output_handler_get_number_cells


    function base_output_handler_get_number_dimensions(this) result(number_dimensions)
    !-----------------------------------------------------------------
    !< Return the number of dimensions
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_dimensions
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        number_dimensions = this%ohcff%get_number_dimensions()
    end function base_output_handler_get_number_dimensions

    function base_output_handler_has_mixed_cell_topologies(this) result(mixed_cell_topologies)
    !-----------------------------------------------------------------
    !< Return if the mesh has mixed cell topologies
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        logical                                  :: mixed_cell_topologies
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        mixed_cell_topologies = this%ohcff%has_mixed_cell_topologies()
    end function base_output_handler_has_mixed_cell_topologies


    function base_output_handler_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_fields
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        number_fields = this%number_fields
    end function base_output_handler_get_number_fields


    function base_output_handler_get_number_cell_vectors(this) result(number_cell_vectors)
    !-----------------------------------------------------------------
    !< Return the number of cell_vectors
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: number_cell_vectors
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        number_cell_vectors = this%number_cell_vectors
    end function base_output_handler_get_number_cell_vectors


    function base_output_handler_get_fe_field(this, field_id) result(field)
    !-----------------------------------------------------------------
    !< Return a [[output_handler_fe_field_t(type)]] given its **ID**
    !-----------------------------------------------------------------
        class(base_output_handler_t),    target, intent(in) :: this
        integer(ip),                             intent(in) :: field_id
        type(output_handler_fe_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        assert(field_id <= this%number_fields)
        field => this%fe_fields(field_id)
    end function base_output_handler_get_fe_field  


    function base_output_handler_get_cell_vector(this, cell_vector_id) result(cell_vector)
    !-----------------------------------------------------------------
    !< Return a cell [[output_handler_cell_vector_t(type)]] given its **ID*
    !-----------------------------------------------------------------
        class(base_output_handler_t),       target, intent(in) :: this
        integer(ip),                                intent(in) :: cell_vector_id
        type(output_handler_cell_vector_t), pointer            :: cell_vector
    !-----------------------------------------------------------------
        assert(cell_vector_id <= this%number_cell_vectors)
        cell_vector => this%cell_vectors(cell_vector_id)
    end function base_output_handler_get_cell_vector


    subroutine base_output_handler_attach_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Attach a [[serial_fe_space_t(type)]]. 
    !< Only a single [[serial_fe_space_t(type)]] is managed by the handler.
    !< The attached [[serial_fe_space_t(type)]] must agree with
    !< all [[fe_function_t(type)]] and **cell_vector** added.
    !-----------------------------------------------------------------
        class(base_output_handler_t),          intent(inout) :: this
        class(serial_fe_space_t), target,      intent(in)    :: fe_space
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        this%fe_space => fe_space
    end subroutine base_output_handler_attach_fe_space
    
    subroutine base_output_handler_set_create_fe_iterator(this, create_fe_iterator)
      class(base_output_handler_t), intent(inout) :: this
       procedure(create_fe_iterator_interface) :: create_fe_iterator
       this%create_fe_iterator => create_fe_iterator
    end subroutine base_output_handler_set_create_fe_iterator

    subroutine base_output_handler_set_free_fe_iterator(this, free_fe_iterator)
      class(base_output_handler_t), intent(inout) :: this
      procedure(free_fe_iterator_interface) :: free_fe_iterator
      this%free_fe_iterator => free_fe_iterator
    end subroutine base_output_handler_set_free_fe_iterator

    subroutine base_output_handler_create_fe_iterator_wrapper(this, fe)
      class(base_output_handler_t)     , intent(inout) :: this
      class(fe_iterator_t), allocatable, intent(inout) :: fe
      if (associated(this%create_fe_iterator)) then
        call this%create_fe_iterator(fe)
      else 
        call this%fe_space%create_fe_iterator(fe)
      end if
    end subroutine base_output_handler_create_fe_iterator_wrapper
    
    subroutine base_output_handler_free_fe_iterator_wrapper(this, fe)
      class(base_output_handler_t)     , intent(inout) :: this
      class(fe_iterator_t), allocatable, intent(inout) :: fe
      if (associated(this%free_fe_iterator)) then
        call this%free_fe_iterator(fe)
      else 
        call this%fe_space%free_fe_iterator(fe)
      end if
    end subroutine base_output_handler_free_fe_iterator_wrapper
    
    function base_output_handler_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Return a [[serial_fe_space_t(type)]] pointer.
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        class(serial_fe_space_t), pointer        :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function base_output_handler_get_fe_space


    subroutine base_output_handler_resize_fe_fields_if_needed(this, number_fields)
    !-----------------------------------------------------------------
    !< Resize the array of added [[output_handler_fe_field_t(type)]]
    !< only if it's needed
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
        integer(ip)                                       :: current_size
        type(output_handler_fe_field_t), allocatable      :: temp_fe_functions(:)
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        if(.not. allocated(this%fe_fields)) then
            allocate(this%fe_fields(cell_node_field_array_size))
        elseif(number_fields > size(this%fe_fields)) then
            current_size = size(this%fe_fields)
            call move_alloc(from=this%fe_fields, to=temp_fe_functions)
            allocate(this%fe_fields(int(1.5*current_size)))
            this%fe_fields(1:current_size) = temp_fe_functions(1:current_size)
            deallocate(temp_fe_functions)
        endif
    end subroutine base_output_handler_resize_fe_fields_if_needed


    subroutine base_output_handler_resize_cell_vectors_if_needed(this, number_fields)
    !-----------------------------------------------------------------
    !< Resize the array of added **cell_vector** only if it's needed
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
        integer(ip)                                       :: current_size
        type(output_handler_cell_vector_t), allocatable   :: temp_cell_vectors(:)
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        if(.not. allocated(this%cell_vectors)) then
            allocate(this%cell_vectors(cell_node_field_array_size))
        elseif(number_fields > size(this%cell_vectors)) then
            current_size = size(this%cell_vectors)
            call move_alloc(from=this%cell_vectors, to=temp_cell_vectors)
            allocate(this%cell_vectors(int(1.5*current_size)))
            this%cell_vectors(1:current_size) = temp_cell_vectors(1:current_size)
            deallocate(temp_cell_vectors)
        endif
    end subroutine base_output_handler_resize_cell_vectors_if_needed


    subroutine base_output_handler_add_fe_function(this, fe_function, field_id, name, diff_operator)
    !-----------------------------------------------------------------
    !< Add the field stored in the [[fe_function_t(type)]] with the
    !< given **field_id**
    !< The attached [[serial_fe_space_t(type)]] must agree with
    !< all [[fe_function_t(type)]] added.
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        type(fe_function_t),                intent(in)    :: fe_function
        integer(ip),                        intent(in)    :: field_id
        character(len=*),                   intent(in)    :: name
        character(len=*), optional,         intent(in)    :: diff_operator
        class(serial_fe_space_t), pointer                 :: fe_space
        character(:), allocatable                         :: field_type
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        assert(associated(this%fe_space))
        assert(field_id >=1 .and. field_id <= this%fe_space%get_number_fields())
        field_type = this%fe_space%get_field_type(field_id)
        call this%resize_fe_fields_if_needed(this%number_fields+1)
        this%number_fields = this%number_fields + 1
        call this%fe_fields(this%number_fields)%set(fe_function, field_id, name, field_type, diff_operator)
    end subroutine base_output_handler_add_fe_function


    subroutine base_output_handler_add_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Add a **cell_vector**. 
    !< The attached [[serial_fe_space_t(type)]] must agree with
    !< all **cell_vector** added.
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        real(rp), allocatable,              intent(in)    :: cell_vector(:)
        character(len=*),                   intent(in)    :: name
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        call this%resize_cell_vectors_if_needed(this%number_cell_vectors+1)
        this%number_cell_vectors = this%number_cell_vectors + 1
        call this%cell_vectors(this%number_cell_vectors)%set(cell_vector, name)
    end subroutine base_output_handler_add_cell_Vector


    subroutine base_output_handler_open(this, dir_path, prefix, parameter_list)
    !-----------------------------------------------------------------
    !< Open procedure. 
    !< Only manages the *State transition diagram* and
    !< call [[base_output_handler_t(type):Open_body(bound)]] of the
    !< concrete object
    !-----------------------------------------------------------------
        class(base_output_handler_t),    intent(inout) :: this
        character(len=*),                intent(in)    :: dir_path
        character(len=*),                intent(in)    :: prefix
        type(ParameterList_t), optional, intent(in)    :: parameter_list
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        call this%open_body(dir_path, prefix, parameter_list)
        this%state = BASE_OUTPUT_HANDLER_STATE_OPEN
    end subroutine base_output_handler_open


    subroutine base_output_handler_close(this)
    !-----------------------------------------------------------------
    !< Close procedure. 
    !< Only manages the *State transition diagram* and
    !< call [[base_output_handler_t(type):Close_body(bound)]] of the
    !< concrete object
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        call this%close_body()
        this%state = BASE_OUTPUT_HANDLER_STATE_INIT
    end subroutine base_output_handler_close


    subroutine base_output_handler_fill_data(this, update_mesh)
    !-----------------------------------------------------------------
    !< Translation of the data contained in [[serial_fe_space_t(type)]],
    !< [[fe_function_t(type)]] and **cell_vector**.
    !< This is the kernel of the implemented strategy for writting to disk. 
    !< It uses [[output_handler_cell_fe_function_t(type)]] in order to fill
    !< the [[output_handler_patch_t(type)]] with local view of the data delimited 
    !< in each cell. 
    !< This procedure is supported by the deferred **append_cell**
    !< procedure implemented in all extended classes.
    !-----------------------------------------------------------------
        class(base_output_handler_t), target, intent(inout) :: this
        logical,                              intent(in)    :: update_mesh
        class(fe_iterator_t), allocatable                   :: fe
        type(output_handler_patch_t)                        :: patch
        type(patch_subcell_iterator_t)                      :: subcell_iterator
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        assert(associated(this%fe_space))
        
        call this%create_fe_iterator_wrapper(fe)
        if(update_mesh) then
            ! Create Output Cell Handler and allocate patch fields
            call this%ohcff%create(fe, this%number_fields, this%fe_fields(1:this%number_fields))
            this%state = BASE_OUTPUT_HANDLER_STATE_FILL

            ! Allocate geometry and connectivity arrays
            call this%allocate_cell_and_nodal_arrays()
        endif

        ! Otherwise ifc issues an error in the fill_patch call below
        if(.not.allocated(this%cell_vectors)) allocate(this%cell_vectors(1))
        
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        call patch%create(this%number_fields, this%number_cell_vectors)
        call fe%first()
        ! Translate coordinates and connectivities to VTK format for every subcell
        do while ( .not. fe%has_finished())
            ! Get Finite element
            if ( fe%is_local() ) then
                call this%ohcff%fill_patch(fe, &
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
            call fe%next()
        end do
        call patch%free()
        call this%free_fe_iterator_wrapper(fe)
    end subroutine base_output_handler_fill_data
end module base_output_handler_names
