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
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Patch-related derived types
!
! Contains the following public entities:
! [[output_handler_patch_names(module)]]
!---------------------------------------------------------------------
module output_handler_patch_names
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Patch-related derived types
!
! A **Patch** is and entity containing raw values (connectivities, 
! coordinates, cell and nodal fields restricted to a single
! finite element
!
! Contains the following public entities:
! [[output_handler_patch_t(type)]], 
! [[output_handler_patch_field_t(type)]], 
! [[patch_subcell_iterator_t(type)]], 
! [[patch_subcell_accessor_t(type)]]
!---------------------------------------------------------------------

USE types_names
USE field_names
USE allocatable_array_names
USE reference_fe_names, only: field_type_scalar, field_type_vector, field_type_tensor, field_type_symmetric_tensor

implicit none
#include "debug.i90"
private

    type :: output_handler_patch_field_t
    !-----------------------------------------------------------------
    !* Author: Víctor Sande Veiga
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    !
    !-----------------------------------------------------------------
    !### Patch-field derived type
    !
    ! A [[output_handler_patch_field_t(type)]] is and entity containing 
    ! raw values (cell and nodal fields) restricted to a single finite element
    !-----------------------------------------------------------------
    private
        character(len=:), allocatable          :: field_type
        type(allocatable_array_rp1_t)          :: nodal_values
        type(allocatable_array_rp1_t)          :: scalar_function_values
        type(allocatable_array_vector_field_t) :: vector_function_values
        type(allocatable_array_tensor_field_t) :: tensor_function_values
    contains
    private
        procedure, non_overridable, public :: free                          => output_handler_patch_field_free
        procedure, non_overridable, public :: set_field_type                => output_handler_patch_field_set_field_type
        procedure, non_overridable, public :: get_field_type                => output_handler_patch_field_get_field_type
        procedure, non_overridable, public :: get_nodal_values              => output_handler_patch_field_get_nodal_values
        procedure, non_overridable, public :: get_scalar_function_values    => output_handler_patch_field_get_scalar_function_values
        procedure, non_overridable, public :: get_vector_function_values    => output_handler_patch_field_get_vector_function_values
        procedure, non_overridable, public :: get_tensor_function_values    => output_handler_patch_field_get_tensor_function_values
        procedure, non_overridable, public :: get_num_components         => output_handler_patch_field_get_num_components
    end type

    type :: output_handler_patch_t
    !-----------------------------------------------------------------
    !* Author: Víctor Sande Veiga
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    !
    !-----------------------------------------------------------------
    !### Patch derived type
    !
    ! A [[output_handler_patch_t(type)]] is and entity containing raw values
    ! (connectivities, coordinates, cell and nodal fields restricted 
    ! to a single finite element
    !-----------------------------------------------------------------
    private
        character(len=:), allocatable              :: cell_type
        integer(ip)                                :: num_dims           = 0
        integer(ip)                                :: num_fields               = 0
        integer(ip)                                :: num_cell_vectors         = 0
        integer(ip)                                :: num_subcells             = 0
        integer(ip)                                :: num_vertices_x_subcell = 0
        type(point_t),                 pointer     :: coordinates(:)
        type(allocatable_array_ip2_t)              :: subcells_connectivity
        type(output_handler_patch_field_t), allocatable :: fields(:)
        type(allocatable_array_rp1_t),      allocatable :: cell_vectors(:)
    contains
    private
        procedure, non_overridable, public :: set_cell_type                   => output_handler_patch_set_cell_type
        procedure, non_overridable, public :: set_num_dims           => output_handler_patch_set_num_dims
        procedure, non_overridable, public :: set_num_points               => output_handler_patch_set_num_points
        procedure, non_overridable, public :: set_num_subcells             => output_handler_patch_set_num_subcells
        procedure, non_overridable, public :: set_num_vertices_x_subcell => output_handler_patch_set_num_vertices_x_subcell
        procedure, non_overridable, public :: set_coordinates                 => output_handler_patch_set_coordinates
        procedure, non_overridable, public :: create                          => output_handler_patch_create
        procedure, non_overridable, public :: free                            => output_handler_patch_free
        procedure, non_overridable, public :: get_cell_type                   => output_handler_patch_get_cell_type
        procedure, non_overridable, public :: get_num_dims           => output_handler_patch_get_num_dims
        procedure, non_overridable, public :: get_subcells_connectivity       => output_handler_patch_get_subcells_connectivity
        procedure, non_overridable, public :: get_num_fields               => output_handler_patch_get_num_fields
        procedure, non_overridable, public :: get_num_cell_vectors         => output_handler_patch_get_num_cell_vectors
        procedure, non_overridable, public :: get_num_subcells             => output_handler_patch_get_num_subcells
        procedure, non_overridable, public :: get_num_vertices_x_subcell => output_handler_patch_get_num_vertices_x_subcell
        procedure, non_overridable, public :: get_coordinates                 => output_handler_patch_get_coordinates
        procedure, non_overridable, public :: get_field                       => output_handler_patch_get_field
        procedure, non_overridable, public :: get_cell_vector                 => output_handler_patch_get_cell_vector
        procedure, non_overridable, public :: get_subcells_iterator           => output_handler_patch_get_subcells_iterator
    end type


    type :: patch_subcell_accessor_t
    !-----------------------------------------------------------------
    !* Author: Víctor Sande Veiga
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    !
    !-----------------------------------------------------------------
    !### Isolate a single linear subcell from a [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
    private
        type(output_handler_patch_t), pointer :: patch => NULL()
        integer(ip)                           :: current_subcell = 0
    contains
    private
        procedure, non_overridable, public :: create                      => patch_subcell_accessor_create
        procedure, non_overridable, public :: free                        => patch_subcell_accessor_free
        procedure, non_overridable, public :: get_cell_type               => patch_subcell_accessor_get_cell_type
        procedure, non_overridable, public :: get_num_dims       => patch_subcell_accessor_get_num_dims
        procedure, non_overridable, public :: get_num_vertices         => patch_subcell_accessor_get_num_vertices
        procedure, non_overridable, public :: get_connectivity            => patch_subcell_accessor_get_connectivity
        procedure, non_overridable, public :: get_num_fields           => patch_subcell_accessor_get_num_fields
        procedure, non_overridable, public :: get_cell_vector             => patch_subcell_accessor_get_cell_vector
        procedure, non_overridable         ::                                patch_subcell_accessor_get_coordinates_X_Y_Z
        procedure, non_overridable         ::                                patch_subcell_accessor_get_coordinates_XYZ
        procedure, non_overridable         ::                                patch_subcell_accessor_get_field_1D
        procedure, non_overridable         ::                                patch_subcell_accessor_get_field_2D
        generic,                    public :: get_coordinates             => patch_subcell_accessor_get_coordinates_X_Y_Z, &
                                                                             patch_subcell_accessor_get_coordinates_XYZ
        generic,                    public :: get_field                   => patch_subcell_accessor_get_field_1D, &
                                                                             patch_subcell_accessor_get_field_2D
    end type

    type :: patch_subcell_iterator_t
    !-----------------------------------------------------------------
    !* Author: Víctor Sande Veiga
    ! Date: 2016-11-29
    ! Version: 0.0.1
    ! Category: IO
    !
    !-----------------------------------------------------------------
    !### Iterate over the linear subcells contained in a [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
    private
        type(output_handler_patch_t), pointer :: patch => NULL()
        integer(ip)                           :: current_subcell = 0
        integer(ip)                           :: num_subcells = 0
    contains
    private
        procedure, non_overridable         :: create                      => patch_subcell_iterator_create
        procedure, non_overridable, public :: begin                       => patch_subcell_iterator_begin
        procedure, non_overridable, public :: next                        => patch_subcell_iterator_next
        procedure, non_overridable, public :: has_finished                => patch_subcell_iterator_has_finished
        procedure, non_overridable, public :: get_accessor                => patch_subcell_iterator_get_accessor
        procedure, non_overridable, public :: free                        => patch_subcell_iterator_free
    end type

public :: output_handler_patch_t
public :: output_handler_patch_field_t
public :: patch_subcell_iterator_t
public :: patch_subcell_accessor_t

contains

!---------------------------------------------------------------------
! output_handler_PATCH_FIELD_T PROCEDURES
!---------------------------------------------------------------------


    subroutine output_handler_patch_field_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(this%field_type)) deallocate(this%field_type)
        call this%nodal_values%free()
        call this%scalar_function_values%free()
        call this%vector_function_values%free()
        call this%tensor_function_values%free()
    end subroutine output_handler_patch_field_free


    subroutine output_handler_patch_field_set_field_type(this, field_type)
    !-----------------------------------------------------------------
    !< Set the field type
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(inout) :: this
        character(len=*),                    intent(in)    :: field_type
    !-----------------------------------------------------------------
        this%field_type = field_type
    end subroutine output_handler_patch_field_set_field_type


    function output_handler_patch_field_get_field_type(this) result(field_type)
    !-----------------------------------------------------------------
    !< Return the field type
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(inout) :: this
        character(len=:), allocatable                      :: field_type
    !-----------------------------------------------------------------
        field_type = this%field_type
    end function output_handler_patch_field_get_field_type


    function output_handler_patch_field_get_nodal_values(this) result(nodal_values)
    !-----------------------------------------------------------------
    !< Return a pointer to nodal_values
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_rp1_t),       pointer               :: nodal_values
    !-----------------------------------------------------------------
        nodal_values => this%nodal_values
    end function output_handler_patch_field_get_nodal_values


    function output_handler_patch_field_get_scalar_function_values(this) result(scalar_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to scalar_function_values
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_rp1_t),       pointer               :: scalar_function_values
    !-----------------------------------------------------------------
        scalar_function_values => this%scalar_function_values
    end function output_handler_patch_field_get_scalar_function_values


    function output_handler_patch_field_get_vector_function_values(this) result(vector_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_values
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_vector_field_t), pointer               :: vector_function_values
    !-----------------------------------------------------------------
        vector_function_values => this%vector_function_values
    end function output_handler_patch_field_get_vector_function_values


    function output_handler_patch_field_get_tensor_function_values(this) result(tensor_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_values
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_tensor_field_t), pointer               :: tensor_function_values
    !-----------------------------------------------------------------
        tensor_function_values => this%tensor_function_values
    end function output_handler_patch_field_get_tensor_function_values


    function output_handler_patch_field_get_num_components(this) result(num_components)
    !-----------------------------------------------------------------
    !< Return the number of components of the field
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(in) :: this
        integer(ip)                                     :: num_components       
    !-----------------------------------------------------------------
        select case(this%field_type)
            case (field_type_scalar)
                num_components = 1
            case (field_type_vector)
                num_components = SPACE_DIM
            case (field_type_tensor)
                num_components = SPACE_DIM*SPACE_DIM
            case (field_type_symmetric_tensor)
                num_components = SPACE_DIM*SPACE_DIM
            case DEFAULT
                check(.false.)
        end select
    end function output_handler_patch_field_get_num_components


!---------------------------------------------------------------------
! output_handler_PATCH_T PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_patch_create(this, num_fields, num_cell_vectors)
    !-----------------------------------------------------------------
    !< Create procedure. Allocate fields
    !-----------------------------------------------------------------
        class(output_handler_patch_t),      intent(inout) :: this
        integer(ip),                        intent(in)    :: num_fields
        integer(ip),                        intent(in)    :: num_cell_vectors
    !-----------------------------------------------------------------
        call this%free()
        this%num_fields       = num_fields
        this%num_cell_vectors = num_cell_vectors
        allocate(this%fields(this%num_fields))
        allocate(this%cell_vectors(this%num_cell_vectors))
    end subroutine output_handler_patch_create


    subroutine output_handler_patch_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip)                                  :: i
    !-----------------------------------------------------------------
        if(allocated(this%cell_type))  deallocate(this%cell_type)
        nullify(this%coordinates)
        call this%subcells_connectivity%free()
        if(allocated(this%fields)) then
            do i=1, this%num_fields
                call this%fields(i)%free()
            enddo
            deallocate(this%fields)
        endif
        if(allocated(this%cell_vectors)) then
            do i=1, this%num_cell_vectors
                call this%cell_vectors(i)%free()
            enddo
            deallocate(this%cell_vectors)
        endif
        this%num_points               = 0
        this%num_fields               = 0
        this%num_cell_vectors         = 0
        this%num_subcells             = 0
        this%num_vertices_x_subcell = 0
    end subroutine output_handler_patch_free


    subroutine output_handler_patch_set_cell_type(this, cell_type)
    !-----------------------------------------------------------------
    !< Set the cell type of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        character(len=*),              intent(in)    :: cell_type
    !-----------------------------------------------------------------
        this%cell_type = cell_type
    end subroutine output_handler_patch_set_cell_type


    subroutine output_handler_patch_set_num_dims(this, num_dims)
    !-----------------------------------------------------------------
    !< Set the number of dimensions of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: num_dims
    !-----------------------------------------------------------------
        this%num_dims = num_dims
    end subroutine output_handler_patch_set_num_dims


    subroutine output_handler_patch_set_num_points(this, num_points)
    !-----------------------------------------------------------------
    !< Set the number of points of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: num_points
    !-----------------------------------------------------------------
        this%num_points = num_points
    end subroutine output_handler_patch_set_num_points


    subroutine output_handler_patch_set_num_subcells(this, num_subcells)
    !-----------------------------------------------------------------
    !< Set the number of subcells of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: num_subcells
    !-----------------------------------------------------------------
        this%num_subcells = num_subcells
    end subroutine output_handler_patch_set_num_subcells


    subroutine output_handler_patch_set_num_vertices_x_subcell(this, num_vertices_x_subcell)
    !-----------------------------------------------------------------
    !< Set the number of vertices per subcell of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: num_vertices_x_subcell
    !-----------------------------------------------------------------
        this%num_vertices_x_subcell = num_vertices_x_subcell
    end subroutine output_handler_patch_set_num_vertices_x_subcell


    subroutine output_handler_patch_set_coordinates(this, coordinates)
    !-----------------------------------------------------------------
    !< Set the coordinates of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        type(point_t), pointer,        intent(in)    :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(coordinates))
        this%coordinates => coordinates
    end subroutine output_handler_patch_set_coordinates


    function output_handler_patch_get_cell_type(this) result(cell_type)
    !-----------------------------------------------------------------
    !< Return the cell type of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        character(len=:), allocatable             :: cell_type
    !-----------------------------------------------------------------
        cell_type = this%cell_type
    end function output_handler_patch_get_cell_type


    pure function output_handler_patch_get_num_dims(this) result(num_dims)
    !-----------------------------------------------------------------
    !< Return the number of dimensions of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: num_dims
    !-----------------------------------------------------------------
        num_dims = this%num_dims
    end function output_handler_patch_get_num_dims


    function output_handler_patch_get_num_subcells(this) result(num_subcells)
    !-----------------------------------------------------------------
    !< Return the number of subcells contained in the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: num_subcells
    !-----------------------------------------------------------------
        num_subcells = this%num_subcells
    end function output_handler_patch_get_num_subcells


    pure function output_handler_patch_get_num_vertices_x_subcell(this) result(num_vertices_x_subcell)
    !-----------------------------------------------------------------
    !< Return the number of vertices per subcell in the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: num_vertices_x_subcell
    !-----------------------------------------------------------------
        num_vertices_x_subcell = this%num_vertices_x_subcell
    end function output_handler_patch_get_num_vertices_x_subcell


    function output_handler_patch_get_num_fields(this) result(num_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields handled by the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: num_fields
    !-----------------------------------------------------------------
        num_fields = this%num_fields
    end function output_handler_patch_get_num_fields


    function output_handler_patch_get_num_cell_vectors(this) result(num_cell_vectors)
    !-----------------------------------------------------------------
    !< Return the number of cell vectors handled by the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: num_cell_vectors
    !-----------------------------------------------------------------
        num_cell_vectors = this%num_cell_vectors
    end function output_handler_patch_get_num_cell_vectors


    function output_handler_patch_get_field(this, num_field) result(field)
    !-----------------------------------------------------------------
    !< Return a fields handled by the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t),      target, intent(in) :: this
        integer(ip),                                intent(in) :: num_field
        type(output_handler_patch_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(num_field <= this%num_fields)
        field => this%fields(num_field)
    end function output_handler_patch_get_field


    function output_handler_patch_get_cell_vector(this, num_cell_vector) result(cell_vector)
    !-----------------------------------------------------------------
    !< Return a cell vector handled by the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t),      target, intent(in) :: this
        integer(ip),                                intent(in) :: num_cell_vector
        type(allocatable_array_rp1_t), pointer                 :: cell_vector
    !-----------------------------------------------------------------
        assert(num_cell_vector <= this%num_cell_vectors)
        cell_vector => this%cell_vectors(num_cell_vector)
    end function output_handler_patch_get_cell_vector


    function output_handler_patch_get_subcells_connectivity(this) result(subcells_connectivity)
    !-----------------------------------------------------------------
    !< Return a pointer to subcells connectivity 
    !-----------------------------------------------------------------
        class(output_handler_patch_t),       target, intent(in) :: this
        type(allocatable_array_ip2_t), pointer                  :: subcells_connectivity
    !-----------------------------------------------------------------
        subcells_connectivity => this%subcells_connectivity
    end function output_handler_patch_get_subcells_connectivity


    function output_handler_patch_get_coordinates(this) result(coordinates)
    !-----------------------------------------------------------------
    !< Set the number of points of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        type(point_t), pointer                       :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(this%coordinates))
        coordinates => this%coordinates
    end function output_handler_patch_get_coordinates


    function output_handler_patch_get_subcells_iterator(this) result(iterator)
    !-----------------------------------------------------------------
    !< Set the number of points of the [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in)     :: this
        type(patch_subcell_iterator_t) :: iterator
    !-----------------------------------------------------------------
        call iterator%create(this)
    end function output_handler_patch_get_subcells_iterator

!---------------------------------------------------------------------
! patch_subcell_iterator_T PROCEDURES
!---------------------------------------------------------------------

    subroutine patch_subcell_iterator_free(this)
    !-----------------------------------------------------------------
    !< Free the [[patch_subcell_iterator_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        nullify(this%patch)
        this%current_subcell = 0
        this%num_subcells = 0
    end subroutine patch_subcell_iterator_free


    subroutine patch_subcell_iterator_create(this, patch)
    !-----------------------------------------------------------------
    !< Create a [[patch_subcell_iterator_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t),      intent(inout) :: this
        type(output_handler_patch_t), target, intent(in)    :: patch
    !-----------------------------------------------------------------
        this%patch => patch
        this%num_subcells = this%patch%get_num_subcells()
        call this%begin()
    end subroutine patch_subcell_iterator_create


    subroutine patch_subcell_iterator_begin(this)
    !-----------------------------------------------------------------
    !< Rewind [[patch_subcell_iterator_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = 1
    end subroutine patch_subcell_iterator_begin


    subroutine patch_subcell_iterator_next(this)
    !-----------------------------------------------------------------
    !< Jump to the next subcell
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = this%current_subcell+1
    end subroutine patch_subcell_iterator_next


    function patch_subcell_iterator_has_finished(this) result(has_finished)
    !-----------------------------------------------------------------
    !< Ask if the [[patch_subcell_iterator_t(type)]] has reached the last position
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t), intent(inout) :: this
        logical                                        :: has_finished
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        has_finished = this%current_subcell > this%num_subcells
    end function patch_subcell_iterator_has_finished


    function patch_subcell_iterator_get_accessor(this) result(accessor)
    !-----------------------------------------------------------------
    !< Return the [[patch_subcell_accessor_t(type)]] for the current position
    !-----------------------------------------------------------------
        class(patch_subcell_iterator_t), intent(inout) :: this
        type(patch_subcell_accessor_t)                 :: accessor
    !-----------------------------------------------------------------
        assert(associated(this%patch) .and. this%current_subcell>0 .and. this%current_subcell<=this%num_subcells)
        call accessor%create(this%patch, this%current_subcell)
    end function patch_subcell_iterator_get_accessor


!---------------------------------------------------------------------
! patch_subcell_accessor_t PROCEDURES
!---------------------------------------------------------------------


    subroutine patch_subcell_accessor_free(this, patch, current_subcell)
    !-----------------------------------------------------------------
    !< Create a [[patch_subcell_accessor_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),      intent(inout) :: this
        type(output_handler_patch_t), target, intent(in)    :: patch
        integer(ip),                          intent(in)    :: current_subcell
    !-----------------------------------------------------------------
        nullify(this%patch)
        this%current_subcell = 0
    end subroutine patch_subcell_accessor_free


    subroutine patch_subcell_accessor_create(this, patch, current_subcell)
    !-----------------------------------------------------------------
    !< Create a [[patch_subcell_accessor_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),      intent(inout) :: this
        type(output_handler_patch_t), target, intent(in)    :: patch
        integer(ip),                          intent(in)    :: current_subcell
    !-----------------------------------------------------------------
        this%patch => patch
        this%current_subcell = current_subcell
    end subroutine patch_subcell_accessor_create


    function patch_subcell_accessor_get_cell_type(this) result(cell_type)
    !-----------------------------------------------------------------
    !< Return the topology type of the current [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t), intent(in) :: this
        character(len=:), allocatable               :: cell_type
    !-----------------------------------------------------------------
        cell_type = this%patch%get_cell_type()
    end function patch_subcell_accessor_get_cell_type


    pure function patch_subcell_accessor_get_num_dims(this) result(num_dims)
    !-----------------------------------------------------------------
    !< Return the number of dimensions of the current [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t), intent(in) :: this
        integer(ip)                                 :: num_dims
    !-----------------------------------------------------------------
        num_dims = this%patch%get_num_dims()
    end function patch_subcell_accessor_get_num_dims


    pure function patch_subcell_accessor_get_num_vertices(this) result(num_vertices)
    !-----------------------------------------------------------------
    !< Return number of vertices per subcell of the current [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t), intent(in)    :: this
        integer(ip)                                    :: num_vertices
    !-----------------------------------------------------------------
        num_vertices = this%patch%get_num_vertices_x_subcell()
    end function patch_subcell_accessor_get_num_vertices


    subroutine patch_subcell_accessor_get_coordinates_X_Y_Z(this, X, Y, Z)
    !-----------------------------------------------------------------
    !< Return [[patch_subcell_accessor_t(type)]] coordinates
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),        intent(in)    :: this
        real(rp),                               intent(inout) :: X(this%patch%get_num_vertices_x_subcell())
        real(rp),                               intent(inout) :: Y(this%patch%get_num_vertices_x_subcell())
        real(rp),                               intent(inout) :: Z(this%patch%get_num_vertices_x_subcell())
        type(point_t),                 pointer                :: patch_coordinates(:)
        type(allocatable_array_ip2_t), pointer                :: subcells_connectivity
        integer(ip)                                           :: num_dims
        integer(ip)                                           :: num_vertices
        integer(ip)                                           :: vertex
    !-----------------------------------------------------------------
        num_vertices       =  this%get_num_vertices()
        num_dims     =  this%get_num_dims()
        patch_coordinates     => this%patch%get_coordinates()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        do vertex = 1, num_vertices
            if(num_dims>=1) X(vertex) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(1)
            if(num_dims>=2) Y(vertex) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(2)
            if(num_dims>=3) Z(vertex) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(3)
        end do
    end subroutine patch_subcell_accessor_get_coordinates_X_Y_Z


    subroutine patch_subcell_accessor_get_coordinates_XYZ(this, XYZ)
    !-----------------------------------------------------------------
    !< Return [[patch_subcell_accessor_t(type)]] coordinates (x1,y1,z1,x2,y2,z2,...)
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),        intent(in)    :: this
        real(rp),                               intent(inout) :: XYZ(this%patch%get_num_vertices_x_subcell()*this%patch%get_num_dims())
        type(point_t),                 pointer                :: patch_coordinates(:)
        type(allocatable_array_ip2_t), pointer                :: subcells_connectivity
        integer(ip)                                           :: num_vertices
        integer(ip)                                           :: vertex
        integer(ip)                                           :: dim
        integer(ip)                                           :: counter
    !-----------------------------------------------------------------
        num_vertices       =  this%get_num_vertices()
        patch_coordinates     => this%patch%get_coordinates()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        counter = 1
        do vertex = 1, num_vertices
            do dim = 1, this%get_num_dims()
                XYZ(counter) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(dim)
                counter = counter + 1
            enddo
        end do
    end subroutine patch_subcell_accessor_get_coordinates_XYZ


    subroutine patch_subcell_accessor_get_connectivity(this, connectivity)
    !-----------------------------------------------------------------
    !< Return [[patch_subcell_accessor_t(type)]] connectivity
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),        intent(in)    :: this
        integer(ip),                            intent(inout) :: connectivity(this%patch%get_num_vertices_x_subcell())
        type(allocatable_array_ip2_t), pointer                :: subcells_connectivity
        integer(ip)                                           :: num_vertices
    !-----------------------------------------------------------------
        num_vertices       =  this%patch%get_num_vertices_x_subcell()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        connectivity(1:num_vertices) = subcells_connectivity%a(1:num_vertices, this%current_subcell)
    end subroutine patch_subcell_accessor_get_connectivity


    function patch_subcell_accessor_get_num_fields(this) result(num_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields of [[output_handler_patch_t(type)]]
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t), intent(in) :: this
        integer(ip)                                 :: num_fields
    !-----------------------------------------------------------------
        num_fields = this%patch%get_num_fields()
    end function patch_subcell_accessor_get_num_fields


    subroutine patch_subcell_accessor_get_cell_vector(this, cell_vector_id, cell_vector)
    !-----------------------------------------------------------------
    !< Return a cell field of the current [[patch_subcell_accessor_t(type)]] 
    !< given its **cell_vector_id**
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),             intent(in)    :: this
        integer(ip),                                 intent(in)    :: cell_vector_id
        real(rp),                                    intent(inout) :: cell_vector(1)
        type(allocatable_array_rp1_t), pointer                     :: cell_vector_Values
        type(output_handler_patch_field_t),     pointer            :: patch_field
        integer(ip)                                                :: num_components
        integer(ip)                                                :: i_comp, j_comp, counter
    !-----------------------------------------------------------------
        cell_vector_values                => this%patch%get_cell_vector(cell_vector_id)
        assert(size(cell_vector_values%a) >= this%current_subcell)

        cell_vector(1) = cell_vector_values%a(this%current_subcell)
    end subroutine patch_subcell_accessor_get_cell_vector


    subroutine patch_subcell_accessor_get_field_1D(this, field_id, field)
    !-----------------------------------------------------------------
    !< Return a *nodal field* of the current [[patch_subcell_accessor_t(type)]] 
    !< given its **field_id**
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),             intent(in)    :: this
        integer(ip),                                 intent(in)    :: field_id
        real(rp),                                    intent(inout) :: field(:)
        type(output_handler_patch_field_t),     pointer            :: patch_field
        type(allocatable_array_ip2_t),          pointer            :: subcells_connectivity
        type(allocatable_array_rp1_t),          pointer            :: scalar_function_values
        type(allocatable_array_vector_field_t), pointer            :: vector_function_values
        type(allocatable_array_tensor_field_t), pointer            :: tensor_function_values
        type(vector_field_t),                   pointer            :: vector_field(:)
        type(tensor_field_t),                   pointer            :: tensor_field(:)
        integer(ip)                                                :: num_components
        integer(ip)                                                :: num_vertices
        integer(ip)                                                :: vertex
        integer(ip)                                                :: i_comp, j_comp, counter
    !-----------------------------------------------------------------
        patch_field       => this%patch%get_field(field_id)
        num_vertices   =  this%get_num_vertices()
        num_components =  patch_field%get_num_components()
        subcells_connectivity => this%patch%get_subcells_connectivity()

        assert(size(field) == num_vertices*num_components)

        counter = 1
        select case(patch_field%get_field_type())
            case (field_type_scalar)
                scalar_function_values => patch_field%get_scalar_function_values()
                do vertex=1, num_vertices
                    field(vertex) = scalar_function_values%a(subcells_connectivity%a(vertex, this%current_subcell))
                enddo
            case (field_type_vector)
                vector_function_values => patch_field%get_vector_function_values()
                vector_field => patch_field%vector_function_values%get_array()
                do vertex=1, num_vertices
                    do i_comp=1, SPACE_DIM
                        field(counter) = vector_field(subcells_connectivity%a(vertex, this%current_subcell))%get(i_comp)
                        counter = counter+1
                    enddo
                enddo
            case (field_type_tensor, field_type_symmetric_tensor)
                tensor_function_values => patch_field%get_tensor_function_values()
                tensor_field => patch_field%tensor_function_values%get_array()
                do vertex=1, num_vertices
                    do i_comp=1, SPACE_DIM
                        do j_comp=1, SPACE_DIM
                            field(counter) = tensor_field(subcells_connectivity%a(vertex, this%current_subcell))%get(i_comp,j_comp)
                            counter = counter+1
                        enddo
                    enddo
                enddo
            case DEFAULT
                check(.false.)
        end select
    end subroutine patch_subcell_accessor_get_field_1D


    subroutine patch_subcell_accessor_get_field_2D(this, field_id, field)
    !-----------------------------------------------------------------
    !< Return a *field* corresponding with the given **field_id** as a 2D matrix
    !-----------------------------------------------------------------
        class(patch_subcell_accessor_t),             intent(in)    :: this
        integer(ip),                                 intent(in)    :: field_id
        real(rp),                                    intent(inout) :: field(:,:)
        type(output_handler_patch_field_t),     pointer            :: patch_field
        type(allocatable_array_ip2_t),          pointer            :: subcells_connectivity
        type(allocatable_array_rp1_t),          pointer            :: scalar_function_values
        type(allocatable_array_vector_field_t), pointer            :: vector_function_values
        type(allocatable_array_tensor_field_t), pointer            :: tensor_function_values
        type(vector_field_t),                   pointer            :: vector_field(:)
        type(tensor_field_t),                   pointer            :: tensor_field(:)
        integer(ip)                                                :: num_components
        integer(ip)                                                :: num_vertices
        integer(ip)                                                :: vertex
        integer(ip)                                                :: i_comp, j_comp
    !-----------------------------------------------------------------
        patch_field       => this%patch%get_field(field_id)
        num_vertices   =  this%get_num_vertices()
        num_components =  patch_field%get_num_components()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        
       ! assert(size(field,1)==num_components)
       ! assert(size(field,2)==num_vertices)
        
        select case(patch_field%get_field_type())
            case (field_type_scalar)
                scalar_function_values => patch_field%get_scalar_function_values()
                do vertex=1, num_vertices
                    field(1,vertex) = scalar_function_values%a(subcells_connectivity%a(vertex, this%current_subcell))
                enddo
            case (field_type_vector)
                vector_function_values => patch_field%get_vector_function_values()
                vector_field => patch_field%vector_function_values%get_array()
                do vertex=1, num_vertices
                    do i_comp=1, SPACE_DIM
                        field(i_comp,vertex) = vector_field(subcells_connectivity%a(vertex, this%current_subcell))%get(i_comp)
                    enddo
                enddo
            case (field_type_tensor, field_type_symmetric_tensor)
                tensor_function_values => patch_field%get_tensor_function_values()
                tensor_field => patch_field%tensor_function_values%get_array()
                do vertex=1, num_vertices
                    do i_comp=1, SPACE_DIM
                        do j_comp=1, SPACE_DIM
                            field(((i_comp-1)*SPACE_DIM)+j_comp,vertex) = tensor_field(subcells_connectivity%a(vertex, this%current_subcell))%get(i_comp,j_comp)
                        enddo
                    enddo
                enddo
            case DEFAULT
                check(.false.)
        end select
    end subroutine patch_subcell_accessor_get_field_2D

end module output_handler_patch_names
