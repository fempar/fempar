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
USE field_names
USE environment_names
USE reference_fe_names
USE fe_space_names,              only: serial_fe_space_t, fe_cell_iterator_t, fe_function_t
USE allocatable_array_names
USE std_vector_real_rp_names
USE output_handler_field_generator_names
USE output_handler_fe_field_names
USE output_handler_patch_names
USE output_handler_fe_cell_function_names
USE output_handler_parameters_names

implicit none
#include "debug.i90"
private

    integer(ip), parameter :: cell_node_field_array_size = 10

    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_INIT   = 0
    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_OPEN   = 1
    integer(ip), parameter :: BASE_OUTPUT_HANDLER_STATE_FILL   = 2

    type :: fill_patch_field_procedure_t
    !-----------------------------------------------------------------
    !*Author: Víctor Sande
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    ! 
    !-----------------------------------------------------------------
    !### Derived type containing a pointer to a procedure interface
    !-----------------------------------------------------------------
        procedure(fill_patch_field_interface), nopass, pointer :: p => NULL()
    end type
    
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
    ! FILLED occurs when metadata is filled from [[output_handler_fe_cell_function_t(type)]].
    !-----------------------------------------------------------------
    private
        class(serial_fe_space_t),            pointer             :: fe_space              => NULL()
        class(fe_cell_iterator_t), allocatable                   :: fe
        type(output_handler_fe_cell_function_t)                  :: ohcff
        type(output_handler_fe_field_t),    allocatable          :: fe_fields(:)
        type(output_handler_field_generator_info_t), allocatable :: field_generators(:)
        type(output_handler_cell_vector_t), allocatable          :: cell_vectors(:)
        procedure(create_fe_cell_iterator_interface)   , pointer :: create_fe_cell_iterator    => NULL()
        procedure(free_fe_cell_iterator_interface)     , pointer :: free_fe_cell_iterator      => NULL()
        integer(ip)                                              :: state                 = BASE_OUTPUT_HANDLER_STATE_INIT
        integer(ip)                                              :: num_fields           = 0
        integer(ip)                                              :: num_field_generators = 0
        integer(ip)                                              :: num_cell_vectors     = 0
        type(fill_patch_field_procedure_t), allocatable:: fill_patch_field(:)
    contains
    private
        procedure, non_overridable, public :: free                          => base_output_handler_free
        procedure, non_overridable, public :: get_num_nodes              => base_output_handler_get_num_nodes
        procedure, non_overridable, public :: get_num_cells              => base_output_handler_get_num_cells
        procedure, non_overridable, public :: has_mixed_cell_topologies     => base_output_handler_has_mixed_cell_topologies
        procedure, non_overridable, public :: get_num_dims         => base_output_handler_get_num_dims
        procedure, non_overridable, public :: get_num_fields             => base_output_handler_get_num_fields
        procedure, non_overridable, public :: get_num_field_generators   => base_output_handler_get_num_field_generators
        procedure, non_overridable, public :: get_num_cell_vectors       => base_output_handler_get_num_cell_vectors
        procedure, non_overridable, public :: get_fe_field                  => base_output_handler_get_fe_field
        procedure, non_overridable, public :: get_field_generator          => base_output_handler_get_field_generator
        procedure, non_overridable, public :: get_cell_vector               => base_output_handler_get_cell_vector
        procedure, non_overridable, public :: get_fe_space                  => base_output_handler_get_fe_space
        procedure, non_overridable         :: resize_fe_fields_if_needed    => base_output_handler_resize_fe_fields_if_needed
        procedure, non_overridable         :: resize_field_generators_if_needed    => base_output_handler_resize_field_generators_if_needed
        procedure, non_overridable         :: resize_cell_vectors_if_needed => base_output_handler_resize_cell_vectors_if_needed
        procedure, non_overridable, public :: attach_fe_space               => base_output_handler_attach_fe_space
        procedure, non_overridable, public :: set_create_fe_cell_iterator        => base_output_handler_set_create_fe_cell_iterator
        procedure, non_overridable, public :: set_free_fe_cell_iterator          => base_output_handler_set_free_fe_cell_iterator
        procedure, non_overridable, private:: create_fe_cell_iterator_wrapper    => base_output_handler_create_fe_cell_iterator_wrapper
        procedure, non_overridable, private:: free_fe_cell_iterator_wrapper      => base_output_handler_free_fe_cell_iterator_wrapper
        procedure, non_overridable, public :: add_fe_function               => base_output_handler_add_fe_function
        procedure, non_overridable, public :: add_field_generator           => base_output_handler_add_field_generator
        procedure, non_overridable, public :: add_cell_vector               => base_output_handler_add_cell_vector
        procedure, non_overridable, public :: update_cell_vector            => base_output_handler_update_cell_vector
        procedure, non_overridable, public :: fill_data                     => base_output_handler_fill_data
        procedure, non_overridable, public :: fill_patch                    => base_output_handler_fill_patch
        procedure, non_overridable, public :: open                          => base_output_handler_open
        procedure, non_overridable, public :: close                         => base_output_handler_close
        
        ! Strategy procedures to fill patch field data
        procedure, non_overridable :: apply_fill_patch_field_strategy => &
                                                            base_output_handler_apply_fill_patch_field_strategy
        
        procedure(base_output_handler_open_body),                      public, deferred :: open_body
        procedure(base_output_handler_append_time_step),               public, deferred :: append_time_step
        procedure(base_output_handler_allocate_cell_and_nodal_arrays),         deferred :: allocate_cell_and_nodal_arrays
        procedure(base_output_handler_append_cell),                            deferred :: append_cell
        procedure(base_output_handler_write),                          public, deferred :: write
        procedure(base_output_handler_close_body),                     public, deferred :: close_body
        procedure(base_output_handler_free_body),                              deferred :: free_body
    end type

    interface
      subroutine create_fe_cell_iterator_interface(this,fe)
        import base_output_handler_t
        import :: fe_cell_iterator_t
        class(base_output_handler_t)             , intent(in)    :: this
        class(fe_cell_iterator_t)   , allocatable, intent(inout) :: fe
      end subroutine create_fe_cell_iterator_interface
      
      subroutine free_fe_cell_iterator_interface(this,fe)
        import base_output_handler_t
        import :: fe_cell_iterator_t
        class(base_output_handler_t)             , intent(in)    :: this
        class(fe_cell_iterator_t)   , allocatable, intent(inout) :: fe
      end subroutine free_fe_cell_iterator_interface
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
    
    interface 
        subroutine fill_patch_field_interface(this,  fe_function, field_id, patch_field)
            import ip
            import fe_function_t
            import output_handler_patch_field_t
            import output_handler_fe_cell_function_t
            type(output_handler_fe_cell_function_t), intent(in)     :: this
            type(fe_function_t),                      intent(in)    :: fe_function
            integer(ip),                              intent(in)    :: field_id
            type(output_handler_patch_field_t),       intent(inout) :: patch_field
        end subroutine fill_patch_field_interface
    end interface

public :: base_output_handler_t, create_fe_cell_iterator_interface, free_fe_cell_iterator_interface 

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
        call this%free_fe_cell_iterator_wrapper(this%fe)
        nullify(this%fe_space)
        nullify(this%create_fe_cell_iterator)
        nullify(this%free_fe_cell_iterator)
        this%num_cell_vectors      = 0
        this%num_fields            = 0
        this%num_field_generators  = 0
        this%state                 = BASE_OUTPUT_HANDLER_STATE_INIT
        if(allocated(this%fill_patch_field)) deallocate(this%fill_patch_field)
    end subroutine base_output_handler_free

    function base_output_handler_get_num_nodes(this) result(num_nodes)
    !-----------------------------------------------------------------
    !< Return the number of visualization nodes
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_nodes
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_nodes = this%ohcff%get_num_nodes()
    end function base_output_handler_get_num_nodes


    function base_output_handler_get_num_cells(this) result(num_cells)
    !-----------------------------------------------------------------
    !< Return the number of visualization cells
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_cells
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_cells = this%ohcff%get_num_cells()
    end function base_output_handler_get_num_cells


    function base_output_handler_get_num_dims(this) result(num_dims)
    !-----------------------------------------------------------------
    !< Return the number of dimensions
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_dims
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_dims = this%ohcff%get_num_dims()
    end function base_output_handler_get_num_dims

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


    function base_output_handler_get_num_fields(this) result(num_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_fields
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_fields = this%num_fields
    end function base_output_handler_get_num_fields
    
    function base_output_handler_get_num_field_generators(this) result(num_field_generators)
    !-----------------------------------------------------------------
    !< Return the number of field generators
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_field_generators
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_field_generators = this%num_field_generators
    end function base_output_handler_get_num_field_generators

    function base_output_handler_get_num_cell_vectors(this) result(num_cell_vectors)
    !-----------------------------------------------------------------
    !< Return the number of cell_vectors
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        integer(ip)                              :: num_cell_vectors
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        num_cell_vectors = this%num_cell_vectors
    end function base_output_handler_get_num_cell_vectors


    function base_output_handler_get_fe_field(this, field_id) result(field)
    !-----------------------------------------------------------------
    !< Return a [[output_handler_fe_field_t(type)]] given its **ID**
    !-----------------------------------------------------------------
        class(base_output_handler_t),    target, intent(in) :: this
        integer(ip),                             intent(in) :: field_id
        type(output_handler_fe_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        assert(field_id <= this%num_fields)
        field => this%fe_fields(field_id)
    end function base_output_handler_get_fe_field  
    
    function base_output_handler_get_field_generator(this, field_generator_id) result(field_generator)
    !-----------------------------------------------------------------
    !< Return a [[output_handler_fe_field_t(type)]] given its **ID**
    !-----------------------------------------------------------------
        class(base_output_handler_t),    target, intent(in) :: this
        integer(ip),                             intent(in) :: field_generator_id
        type(output_handler_field_generator_info_t), pointer    :: field_generator
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        assert(field_generator_id <= this%get_num_field_generators())
        field_generator => this%field_generators(field_generator_id)
    end function base_output_handler_get_field_generator


    function base_output_handler_get_cell_vector(this, cell_vector_id) result(cell_vector)
    !-----------------------------------------------------------------
    !< Return a cell [[output_handler_cell_vector_t(type)]] given its **ID*
    !-----------------------------------------------------------------
        class(base_output_handler_t),       target, intent(in) :: this
        integer(ip),                                intent(in) :: cell_vector_id
        type(output_handler_cell_vector_t), pointer            :: cell_vector
    !-----------------------------------------------------------------
        assert(cell_vector_id <= this%num_cell_vectors)
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
    
    subroutine base_output_handler_set_create_fe_cell_iterator(this, create_fe_cell_iterator)
      class(base_output_handler_t), intent(inout) :: this
       procedure(create_fe_cell_iterator_interface) :: create_fe_cell_iterator
       this%create_fe_cell_iterator => create_fe_cell_iterator
    end subroutine base_output_handler_set_create_fe_cell_iterator

    subroutine base_output_handler_set_free_fe_cell_iterator(this, free_fe_cell_iterator)
      class(base_output_handler_t), intent(inout) :: this
      procedure(free_fe_cell_iterator_interface) :: free_fe_cell_iterator
      this%free_fe_cell_iterator => free_fe_cell_iterator
    end subroutine base_output_handler_set_free_fe_cell_iterator

    subroutine base_output_handler_create_fe_cell_iterator_wrapper(this, fe)
      class(base_output_handler_t)     , intent(inout) :: this
      class(fe_cell_iterator_t), allocatable, intent(inout) :: fe
      if (associated(this%create_fe_cell_iterator)) then
        call this%create_fe_cell_iterator(fe)
      else 
        call this%fe_space%create_fe_cell_iterator(fe)
      end if
    end subroutine base_output_handler_create_fe_cell_iterator_wrapper
    
    subroutine base_output_handler_free_fe_cell_iterator_wrapper(this, fe)
      class(base_output_handler_t)     , intent(inout) :: this
      class(fe_cell_iterator_t), allocatable, intent(inout) :: fe
      if (associated(this%free_fe_cell_iterator)) then
        call this%free_fe_cell_iterator(fe)
      else 
        call this%fe_space%free_fe_cell_iterator(fe)
      end if
    end subroutine base_output_handler_free_fe_cell_iterator_wrapper
    
    function base_output_handler_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Return a [[serial_fe_space_t(type)]] pointer.
    !-----------------------------------------------------------------
        class(base_output_handler_t), intent(in) :: this
        class(serial_fe_space_t), pointer        :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function base_output_handler_get_fe_space


    subroutine base_output_handler_resize_fe_fields_if_needed(this, num_fields)
    !-----------------------------------------------------------------
    !< Resize the array of added [[output_handler_fe_field_t(type)]]
    !< only if it's needed
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: num_fields
        integer(ip)                                       :: current_size
        type(output_handler_fe_field_t), allocatable      :: temp_fe_functions(:)
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        if(.not. allocated(this%fe_fields)) then
            allocate(this%fe_fields(cell_node_field_array_size))
        elseif(num_fields > size(this%fe_fields)) then
            current_size = size(this%fe_fields)
            call move_alloc(from=this%fe_fields, to=temp_fe_functions)
            allocate(this%fe_fields(int(1.5*current_size)))
            this%fe_fields(1:current_size) = temp_fe_functions(1:current_size)
            deallocate(temp_fe_functions)
        endif
    end subroutine base_output_handler_resize_fe_fields_if_needed
    
   subroutine base_output_handler_resize_field_generators_if_needed(this, num_field_generators)
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: num_field_generators
        integer(ip)                                       :: current_size
        type(output_handler_field_generator_info_t), allocatable :: temp_field_generators(:)
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        if(.not. allocated(this%field_generators)) then
            allocate(this%field_generators(cell_node_field_array_size))
        elseif(num_field_generators > size(this%field_generators)) then
            current_size = size(this%field_generators)
            call move_alloc(from=this%field_generators, to=temp_field_generators)
            allocate(this%field_generators(int(1.5*current_size)))
            this%field_generators(1:current_size) = temp_field_generators(1:current_size)
            deallocate(temp_field_generators)
        endif
    end subroutine base_output_handler_resize_field_generators_if_needed
    
    


    subroutine base_output_handler_resize_cell_vectors_if_needed(this, num_fields)
    !-----------------------------------------------------------------
    !< Resize the array of added **cell_vector** only if it's needed
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        integer(ip),                        intent(in)    :: num_fields
        integer(ip)                                       :: current_size
        type(output_handler_cell_vector_t), allocatable   :: temp_cell_vectors(:)
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        if(.not. allocated(this%cell_vectors)) then
            allocate(this%cell_vectors(cell_node_field_array_size))
        elseif(num_fields > size(this%cell_vectors)) then
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
        assert(field_id >=1 .and. field_id <= this%fe_space%get_num_fields())
        field_type = this%fe_space%get_field_type(field_id)
        call this%resize_fe_fields_if_needed(this%num_fields+1)
        this%num_fields = this%num_fields + 1
        call this%fe_fields(this%num_fields)%set(fe_function, field_id, name, field_type, diff_operator)
    end subroutine base_output_handler_add_fe_function
    
    subroutine base_output_handler_add_field_generator(this, name, output_handler_field_generator)
        class(base_output_handler_t)           ,    intent(inout) :: this
        character(len=*)                       ,    intent(in)    :: name
        class(output_handler_field_generator_t),    intent(in)    :: output_handler_field_generator
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        assert(associated(this%fe_space))
        call this%resize_field_generators_if_needed(this%num_field_generators+1)
        this%num_field_generators = this%num_field_generators + 1
        call this%field_generators(this%num_field_generators)%create(name, output_handler_field_generator)
    end subroutine base_output_handler_add_field_generator

    subroutine base_output_handler_add_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Adds a new **cell_vector** entry to base_output_handler_t identified 
    !< by **name**. The user is responsible for ensuring compatibility 
    !< between the [[serial_fe_space_t(type)]] attached and the newly
    !< added **cell_vector**.
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        real(rp)                    ,       intent(in)    :: cell_vector(:)
        character(len=*),                   intent(in)    :: name
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        call this%resize_cell_vectors_if_needed(this%num_cell_vectors+1)
        this%num_cell_vectors = this%num_cell_vectors + 1
        call this%cell_vectors(this%num_cell_vectors)%set(cell_vector, name)
    end subroutine base_output_handler_add_cell_Vector
    
    subroutine base_output_handler_update_cell_vector(this, cell_vector, name)
    !-----------------------------------------------------------------
    !< Update the pointer to which an existing cell_vector entry already
    !< registered in output_handler_t is pointing to s.t. it points to
    !< **cell_vector**. The provided **name** is used to locate the cell_vector
    !< within output_handler_t. If **name** is not found, a warning message is issued. 
    !< The user is responsible for ensuring compatibility between the [[serial_fe_space_t(type)]] 
    !< attached and the updated **cell_vector**.
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        real(rp)                    ,       intent(in)    :: cell_vector(:)
        character(len=*),                   intent(in)    :: name
    !-----------------------------------------------------------------
        integer(ip) :: i
        logical :: found
        found = .false. 
        do i=1, this%num_cell_vectors
           if ( this%cell_vectors(i)%get_name() == name ) then
             found = .true. 
             exit
           end if
        end do
        wassert( found, "base_output_handler_update_cell_vector :: name not found, cell_vector cannot be updated")
        if ( found ) then
          call this%cell_vectors(i)%set_cell_vector(cell_vector)
        end if
    end subroutine base_output_handler_update_cell_vector
    
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
        character(len=:), allocatable :: field_type, diff_operator
        integer(ip) :: num_field
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_INIT)
        
        ! Configure fill_patch_field strategy for each field
        if(allocated(this%fill_patch_field)) deallocate(this%fill_patch_field)
          allocate(this%fill_patch_field(this%num_fields))
          do num_field = 1, this%num_fields
              field_type    = this%fe_fields(num_field)%get_field_type()
              diff_operator = this%fe_fields(num_field)%get_diff_operator()
              call this%apply_fill_patch_field_strategy(field_type, diff_operator, this%fill_patch_field(num_field)%p)
          end do
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
    
    
    subroutine base_output_handler_apply_fill_patch_field_strategy(this, field_type, diff_operator, proc)
    !-----------------------------------------------------------------
    !< Choose strategy to fill a patch field.
    !< Patch field calculation is distributed in several procedures
    !< in order to apply some diff operators.
    !-----------------------------------------------------------------
        class(base_output_handler_t),       intent(inout) :: this
        character(len=:), allocatable,                  intent(in)    :: field_type
        character(len=:), allocatable,                  intent(in)    :: diff_operator
        procedure(fill_patch_field_interface), pointer, intent(inout) :: proc
    !-----------------------------------------------------------------
        select case(field_type)
            ! Select procedures to fill patch field from a scalar field  (Value or Grad)
            case ( field_type_scalar )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => fill_patch_scalar_field_val
                    case (grad_diff_operator)
                        proc => fill_patch_scalar_field_grad
                    case DEFAULT
                        check(.false.)
                end select
            ! Select procedures to fill patch field from a vector field (Value or Grad or Div or Curl)
            case ( field_type_vector )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => fill_patch_vector_field_val
                    case (grad_diff_operator)
                        proc => fill_patch_vector_field_grad
                    case (div_diff_operator)
                        proc => fill_patch_vector_field_div
                    case (curl_diff_operator)
                        proc => fill_patch_vector_field_curl
                    case DEFAULT
                        check(.false.)
                end select
            ! Select procedures to fill patch field from a tensor field (only Value)
            case ( field_type_tensor )
                select case (diff_operator)
                    case (no_diff_operator)
                        proc => fill_patch_tensor_field_val
                    case DEFAULT
                        check(.false.)
                end select
        end select
    end subroutine base_output_handler_apply_fill_patch_field_strategy


    subroutine fill_patch_scalar_field_val(ohcff, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a scalar fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        
        ! Get reference_Fe
        reference_fe => current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_scalar)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! Calculate scalar field values
        call patch_field%set_field_type(field_type_scalar)
        patch_field_scalar_function_values => patch_field%get_scalar_function_values()
        call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
        call cell_integrator%evaluate_fe_function(nodal_values, scalar_function_values)
        call patch_field_scalar_function_values%move_alloc_in(scalar_function_values) 

    end subroutine fill_patch_scalar_field_val


    subroutine fill_patch_scalar_field_grad(ohcff, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a scalar fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)     :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(cell_integrator_t),                  pointer       :: cell_integrator
        type(std_vector_real_rp_t),               pointer       :: patch_field_nodal_values
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe => current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_scalar)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! Calculate scalar field gradients
        call patch_field%set_field_type(field_type_vector)
        patch_field_vector_function_values => patch_field%get_vector_function_values()
        call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
        call cell_integrator%evaluate_gradient_fe_function(nodal_values, vector_function_values)
        call patch_field_vector_function_values%move_alloc_in(vector_function_values) 

    end subroutine fill_patch_scalar_field_grad


    subroutine fill_patch_vector_field_val(ohcff, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a vector fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe => current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! Calculate vector field values
        call patch_field%set_field_type(field_type_vector)
        patch_field_vector_function_values => patch_field%get_vector_function_values()
        call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
        call cell_integrator%evaluate_fe_function(nodal_values, vector_function_values)
        call patch_field_vector_function_values%move_alloc_in(vector_function_values) 

    end subroutine fill_patch_vector_field_val


    subroutine fill_patch_vector_field_grad(ohcff, fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a vector fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
     !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe => current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! Calculate vector field gradients
        call patch_field%set_field_type(field_type_tensor)
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 
        call cell_integrator%evaluate_gradient_fe_function(nodal_values, tensor_function_values)
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values) 

    end subroutine fill_patch_vector_field_grad


    subroutine fill_patch_vector_field_div(ohcff,  fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field divergence given a vector fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(quadrature_t),                       pointer       :: quadrature
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        integer(ip)                                             :: qpoint
        integer(ip)                                             :: dim
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe    => current_fe%get_reference_fe(field_id)
        quadrature      => ohcff%get_quadrature()
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        call patch_field%set_field_type(field_type_scalar)

        ! get scalar and tensor function values
        patch_field_scalar_function_values => patch_field%get_scalar_function_values()
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 

        ! Calculate gradients
        call cell_integrator%evaluate_gradient_fe_function(nodal_values, tensor_function_values)

        ! Allocate scalar function values
        if ( allocated(scalar_function_values) ) then
           call memrealloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
        else
           call memalloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
        end if
        
        ! Calculate divergence
        scalar_function_values = 0._rp
        do qpoint = 1, quadrature%get_num_quadrature_points()
            do dim = 1, ohcff%get_num_dims()
                scalar_function_values(qpoint) = scalar_function_values(qpoint) + tensor_function_values(qpoint)%get(dim, dim)
            enddo
        enddo 

        ! return scalar and tensor function values
        call patch_field_scalar_function_values%move_alloc_in(scalar_function_values)
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine fill_patch_vector_field_div


    subroutine fill_patch_vector_field_curl(ohcff,  fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field gradients given a vector fe_field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(quadrature_t),                       pointer       :: quadrature
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        real(rp),              allocatable                      :: scalar_function_values(:)
        type(vector_field_t),  allocatable                      :: vector_function_values(:)
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_rp1_t),            pointer       :: patch_field_scalar_function_values
        type(allocatable_array_vector_field_t),   pointer       :: patch_field_vector_function_values
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        integer(ip)                                             :: qpoint
        integer(ip)                                             :: istat
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe    => current_fe%get_reference_fe(field_id)
        quadrature      => ohcff%get_quadrature()
        assert(reference_fe%get_field_type() == field_type_vector)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! get tensor function values
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 

        ! Calculate gradients
        call cell_integrator%evaluate_gradient_fe_function(nodal_values, tensor_function_values)

        if(ohcff%get_num_dims() == 2) then
            call patch_field%set_field_type(field_type_scalar)
            patch_field_scalar_function_values => patch_field%get_scalar_function_values()
            call patch_field_scalar_function_values%move_alloc_out(scalar_function_values) 
            ! Allocate scalar function values
            if ( allocated(scalar_function_values) ) then
               call memrealloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
            else
               call memalloc(quadrature%get_num_quadrature_points(), scalar_function_values, __FILE__, __LINE__)
            end if
            ! Calculate curl
            do qpoint = 1, quadrature%get_num_quadrature_points()
                scalar_function_values(qpoint) = tensor_function_values(qpoint)%get(1,2)-tensor_function_values(qpoint)%get(2,1)
            enddo
            call patch_field_scalar_function_values%move_alloc_in(scalar_function_values)     

        elseif(ohcff%get_num_dims() == 3) then
            call patch_field%set_field_type(field_type_vector)
            patch_field_vector_function_values => patch_field%get_vector_function_values()
            call patch_field_vector_function_values%move_alloc_out(vector_function_values) 
            ! Allocate vector function values
            if ( allocated(vector_function_values) ) then
               if ( size(vector_function_values) < quadrature%get_num_quadrature_points() ) then
                  deallocate(vector_function_values, stat=istat); check(istat==0)
                  allocate(vector_function_values(quadrature%get_num_quadrature_points()), stat=istat); check(istat==0)
               endif
            else
               allocate(vector_function_values(quadrature%get_num_quadrature_points()), stat=istat); check(istat==0)
            end if
            ! Calculate curl
            do qpoint = 1, quadrature%get_num_quadrature_points()
                call vector_function_values(qpoint)%set(1, &
                        tensor_function_values(qpoint)%get(2,3)-tensor_function_values(qpoint)%get(3,2))
                call vector_function_values(qpoint)%set(2, &
                        tensor_function_values(qpoint)%get(3,1)-tensor_function_values(qpoint)%get(1,3))
                call vector_function_values(qpoint)%set(3, &
                        tensor_function_values(qpoint)%get(1,2)-tensor_function_values(qpoint)%get(2,1))
            enddo
            call patch_field_vector_function_values%move_alloc_in(vector_function_values) 
        endif

        ! return tensor function values
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine fill_patch_vector_field_curl


    subroutine fill_patch_tensor_field_val(ohcff,  fe_function, field_id, patch_field)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a tensor field
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)    :: ohcff
        type(fe_function_t),                      intent(in)    :: fe_function
        integer(ip),                              intent(in)    :: field_id
        type(output_handler_patch_field_t),       intent(inout) :: patch_field
        class(reference_fe_t),                    pointer       :: reference_fe
        type(cell_integrator_t),                pointer       :: cell_integrator
        type(std_vector_real_rp_t),            pointer       :: patch_field_nodal_values
        type(tensor_field_t),  allocatable                      :: tensor_function_values(:)
        type(allocatable_array_tensor_field_t),   pointer       :: patch_field_tensor_function_values
        class(fe_cell_iterator_t),                pointer       :: current_fe
        real(rp), pointer :: nodal_values(:)
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Get reference_Fe
        reference_fe => current_fe%get_reference_fe(field_id)
        assert(reference_fe%get_field_type() == field_type_tensor)

        ! Get cell integrator
        cell_integrator => ohcff%get_cell_integrator(field_id) 

        ! Gather DoFs of current cell + field_id on nodal_values 
        patch_field_nodal_values => patch_field%get_nodal_values()
        call patch_field_nodal_values%resize(reference_fe%get_num_shape_functions())
        nodal_values => patch_field_nodal_values%get_pointer()
        call fe_function%gather_nodal_values(current_fe, field_id, nodal_values)

        ! Calculate tensor field values
        call patch_field%set_field_type(field_type_tensor)
        patch_field_tensor_function_values => patch_field%get_tensor_function_values()
        call patch_field_tensor_function_values%move_alloc_out(tensor_function_values) 
        call cell_integrator%evaluate_fe_function(nodal_values, tensor_function_values )
        call patch_field_tensor_function_values%move_alloc_in(tensor_function_values)

    end subroutine fill_patch_tensor_field_val

    subroutine fill_patch_cell_vector( ohcff,cell_vector, num_subcells, patch_cell_vector)
    !-----------------------------------------------------------------
    !< Fill the [[output_handler_patch_field_t(type)]] with field values given a **cell_vector**
    !-----------------------------------------------------------------
        type(output_handler_fe_cell_function_t), intent(in)     :: ohcff
        type(output_handler_cell_vector_t),       intent(in)    :: cell_vector
        integer(ip),                              intent(in)    :: num_subcells
        type(std_vector_real_rp_t),               intent(inout) :: patch_cell_vector
        real(rp), pointer                                       :: cell_vector_values(:)
        real(rp), pointer                                       :: patch_cell_vector_values(:)
        class(fe_cell_iterator_t),                pointer       :: current_fe
    !-----------------------------------------------------------------
        current_fe => ohcff%get_fe()
        ! Gather DoFs of current cell + field_id on nodal_values 
        cell_vector_values => cell_vector%get_cell_vector()
        call patch_cell_vector%resize(num_subcells)
        patch_cell_vector_values => patch_cell_vector%get_pointer()
        patch_cell_vector_values = cell_vector_values(current_fe%get_gid())
    end subroutine fill_patch_cell_vector

    subroutine base_output_handler_fill_data(this, update_mesh)
    !-----------------------------------------------------------------
    !< Translation of the data contained in [[serial_fe_space_t(type)]],
    !< [[fe_function_t(type)]] and **cell_vector**.
    !< This is the kernel of the implemented strategy for writting to disk. 
    !< It uses [[output_handler_fe_cell_function_t(type)]] in order to fill
    !< the [[output_handler_patch_t(type)]] with local view of the data delimited 
    !< in each cell. 
    !< This procedure is supported by the deferred **append_cell**
    !< procedure implemented in all extended classes.
    !-----------------------------------------------------------------
        class(base_output_handler_t), target, intent(inout) :: this
        logical,                              intent(in)    :: update_mesh
        type(output_handler_patch_t)                        :: patch
        type(patch_subcell_iterator_t)                      :: subcell_iterator
    !-----------------------------------------------------------------
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_OPEN .or. this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        assert(associated(this%fe_space))
        
        if(update_mesh) then
            ! Create Output Cell Handler and allocate patch fields
            call this%create_fe_cell_iterator_wrapper(this%fe)
            call this%ohcff%create(this%fe)
            this%state = BASE_OUTPUT_HANDLER_STATE_FILL

            ! Allocate geometry and connectivity arrays
            call this%allocate_cell_and_nodal_arrays()
        endif

        ! Otherwise ifc issues an error in the fill_patch call below
        if(.not.allocated(this%cell_vectors)) allocate(this%cell_vectors(1))
        
        assert(this%state == BASE_OUTPUT_HANDLER_STATE_FILL)
        call patch%create(this%num_fields+this%num_field_generators, this%num_cell_vectors)
        call this%fe%first()
        ! Translate coordinates and connectivities to VTK format for every subcell
        do while ( .not. this%fe%has_finished())
            ! Get Finite element
            if ( this%fe%is_local() ) then
                call this%fill_patch(patch)
                subcell_iterator = patch%get_subcells_iterator()
!               ! Fill data
                do while(.not. subcell_iterator%has_finished())
                    call this%append_cell(subcell_iterator%get_accessor())
                    call subcell_iterator%next()
                enddo
            endif
            call this%fe%next()
        end do
        call patch%free()
    end subroutine base_output_handler_fill_data
    
    subroutine base_output_handler_fill_patch(this,  patch)
    !-----------------------------------------------------------------
    !< Fill a [[output_handler_patch_t(type)]] from a given [[fe_cell_iterator_t(type)]].
    !< The **pach** contains a local view of the coordinates, connectivities 
    !< and field data per cell.
    !-----------------------------------------------------------------
        class(base_output_handler_t),              intent(inout) :: this
        type(output_handler_patch_t),              intent(inout) :: patch
        integer(ip)                                              :: reference_fe_id
        integer(ip)                                              :: idx
        integer(ip)                                              :: field_id
        integer(ip)                                              :: max_order_within_fe
        class(serial_fe_space_t),          pointer               :: fe_space
        type(fe_function_t),               pointer               :: fe_function
        class(reference_fe_t),             pointer               :: reference_fe_geo
        class(environment_t),              pointer               :: environment
        type(point_t),                     pointer               :: coordinates(:)
        type(cell_map_t),                  pointer               :: cell_map
        type(cell_integrator_t),           pointer               :: cell_integrator
        type(quadrature_t),                pointer               :: quadrature
        type(output_handler_patch_field_t),pointer               :: patch_field
        type(std_vector_real_rp_t)        ,pointer               :: patch_cell_vector
        type(allocatable_array_ip2_t),     pointer               :: patch_subcells_connectivity
        character(len=:), allocatable                            :: field_type
        character(len=:), allocatable                            :: diff_operator
        integer(ip), allocatable                                 :: patch_subcells_connectivity_a(:,:)
        class(output_handler_field_generator_t), pointer         :: field_generator
        class(fe_cell_iterator_t)   ,  pointer                   :: fe_cell_iterator
        integer(ip)                                              :: fe_get_lev
    !-----------------------------------------------------------------
        fe_cell_iterator => this%ohcff%get_fe()
        fe_space => fe_cell_iterator%get_fe_space()
        environment => fe_space%get_environment()
        if (environment%am_i_l1_task()) then
            max_order_within_fe = max(fe_cell_iterator%get_max_order_all_fields(),1)
            reference_fe_geo    => fe_cell_iterator%get_reference_fe_geo()
            cell_map              => this%ohcff%get_cell_map()
            coordinates         => cell_map%get_coordinates()
            call fe_cell_iterator%get_nodes_coordinates(coordinates)
            fe_get_lev = fe_cell_iterator%get_level()
            quadrature => this%ohcff%get_quadrature()
            call cell_map%update(fe_get_lev,quadrature,no_ressemblance)
            
            do field_id=1, fe_cell_iterator%get_num_fields()
              cell_integrator => this%ohcff%get_cell_integrator(field_id)
              call cell_integrator%update(fe_get_lev,no_ressemblance,cell_map)
            end do 

            ! Set subcell information into patch
            call patch%set_cell_type(reference_fe_geo%get_topology())
            call patch%set_num_dims(reference_fe_geo%get_num_dims())
            call patch%set_num_vertices_x_subcell(quadrature%get_num_quadrature_points())
            call patch%set_num_subcells(reference_fe_geo%get_num_subcells(num_refinements=max_order_within_fe-1))
            call patch%set_num_vertices_x_subcell(reference_fe_geo%get_num_vertices())

            ! Set patch coordinates from cell_map
            call patch%set_coordinates(cell_map%get_quadrature_points_coordinates())

            ! Set patch connectivities from reference_fe_geo given num_refinements
            patch_subcells_connectivity => patch%get_subcells_connectivity()
            
            call patch_subcells_connectivity%resize(reference_fe_geo%get_num_vertices(), &
                                                    reference_fe_geo%get_num_subcells(num_refinements=max_order_within_fe-1))
           
            
            call reference_fe_geo%get_subcells_connectivity(num_refinements=max_order_within_fe-1, &
                                                            connectivity=patch_subcells_connectivity%a)

            ! Fill patch fe field data
            do idx = 1, this%num_fields
                fe_function  => this%fe_fields(idx)%get_fe_function()
                field_id     =  this%fe_fields(idx)%get_field_id()
                patch_field  => patch%get_field(idx)

                assert(associated(this%fill_patch_field(idx)%p))
                call this%fill_patch_field(idx)%p(this%ohcff, fe_function, field_id, patch_field)
            end do
            
            ! Fill patch fe field data using field generators
            do idx = 1, this%num_field_generators
                patch_field     => patch%get_field(this%num_fields+idx)
                field_generator => this%field_generators(idx)%get_field_generator()
                call field_generator%generate_patch_field(this%ohcff, cell_map%get_quadrature_points_coordinates(), patch_field)
            end do
            

            ! Fill patch cell vectors data
            do idx = 1, this%num_cell_vectors
                patch_cell_vector  => patch%get_cell_vector(idx)
                call fill_patch_cell_vector(this%ohcff, this%cell_vectors(idx), patch%get_num_subcells(), patch_cell_vector)
            end do
        end if
    end subroutine base_output_handler_fill_patch
    
    
end module base_output_handler_names
