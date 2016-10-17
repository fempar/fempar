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

module data_out_patch_names

USE types_names
USE field_names
USE allocatable_array_names
USE fe_function_names,           only: fe_function_t

implicit none
#include "debug.i90"
private

    type :: data_out_patch_field_t
!    private
        type(fe_function_t),  pointer          :: fe_function => NULL()
        character(len=:), allocatable          :: name
        integer(ip)                            :: field_id = 0
        type(allocatable_array_rp1_t)          :: nodal_values
        ! Values + gradients for scalar fields
        type(allocatable_array_rp1_t)          :: scalar_function_values
        type(allocatable_array_vector_field_t) :: scalar_function_gradients
        
        ! Values + gradients for vector fields
        type(allocatable_array_vector_field_t) :: vector_function_values
        type(allocatable_array_tensor_field_t) :: vector_function_gradients
        
        ! Values for tensor fields (gradients not supported yet)
        type(allocatable_array_tensor_field_t) :: tensor_function_values
    contains
        procedure, public :: free                          => data_out_patch_field_free
        procedure, public :: get_name                      => data_out_patch_field_get_name
        procedure, public :: get_field_id                  => data_out_patch_field_get_field_id
        procedure, public :: get_fe_function               => data_out_patch_field_get_fe_function
        procedure, public :: get_nodal_values              => data_out_patch_field_get_nodal_values
        procedure, public :: get_scalar_function_values    => data_out_patch_field_get_scalar_function_values
        procedure, public :: get_scalar_function_gradients => data_out_patch_field_get_scalar_function_gradients
        procedure, public :: get_vector_function_values    => data_out_patch_field_get_vector_function_values
        procedure, public :: get_vector_function_gradients => data_out_patch_field_get_vector_function_gradients
        procedure, public :: get_tensor_function_values    => data_out_patch_field_get_tensor_function_values
    end type

    type :: data_out_patch_t
    private
        integer(ip)                                :: number_dimensions = 0
        integer(ip)                                :: number_points = 0
        integer(ip)                                :: number_fields = 0
        integer(ip)                                :: number_subcells = 0
        integer(ip)                                :: number_vertices_per_subcell = 0
        type(point_t),                 pointer     :: coordinates(:)
        type(allocatable_array_ip2_t)              :: subcells_connectivity
        type(data_out_patch_field_t),  public,allocatable :: fields(:)
    contains
        procedure, public :: set_number_dimensions           => data_out_patch_set_number_dimensions
        procedure, public :: set_number_points               => data_out_patch_set_number_points
        procedure, public :: set_number_subcells             => data_out_patch_set_number_subcells
        procedure, public :: set_number_vertices_per_subcell => data_out_patch_set_number_vertices_per_subcell
        procedure, public :: set_coordinates                 => data_out_patch_set_coordinates
        procedure, public :: free                            => data_out_patch_free
        procedure, public :: get_number_dimensions           => data_out_patch_get_number_dimensions
        procedure, public :: get_subcells_connectivity       => data_out_patch_get_subcells_connectivity
        procedure, public :: get_number_fields               => data_out_patch_get_number_fields
        procedure, public :: get_number_subcells             => data_out_patch_get_number_subcells
        procedure, public :: get_number_vertices_per_subcell => data_out_patch_get_number_vertices_per_subcell
        procedure, public :: get_coordinates                 => data_out_patch_get_coordinates
        procedure, public :: get_field                       => data_out_patch_get_field
        procedure, public :: get_subcells_iterator           => data_out_patch_get_subcells_iterator
    end type

    type :: data_out_patch_subcell_iterator_t
    private
        type(data_out_patch_t), pointer :: patch => NULL()
        integer(ip)                     :: current_subcell = 0
        integer(ip)                     :: number_subcells = 0
    contains
        procedure, public :: create          => data_out_patch_subcell_iterator_create
        procedure, public :: begin           => data_out_patch_subcell_iterator_begin
        procedure, public :: next            => data_out_patch_subcell_iterator_next
        procedure, public :: has_finished    => data_out_patch_subcell_iterator_has_finished
        procedure, public :: free            => data_out_patch_subcell_iterator_free
        procedure, public :: get_coordinates => data_out_patch_subcell_iterator_get_coordinates
    end type

public :: data_out_patch_t
public :: data_out_patch_field_t
public :: data_out_patch_subcell_iterator_t

contains

!---------------------------------------------------------------------
!< DATA_OUT_PATCH_FIELD_T PROCEDURES
!---------------------------------------------------------------------


    subroutine data_out_patch_field_free(this)
        class(data_out_patch_field_t), intent(inout) :: this
        nullify(this%fe_function)
        if(allocated(this%name)) deallocate(this%name)
        this%field_id = 0
        call this%nodal_values%free()
        call this%scalar_function_values%free()
        call this%scalar_function_gradients%free()
        call this%vector_function_values%free()
        call this%vector_function_gradients%free()
        call this%tensor_function_values%free()
    end subroutine data_out_patch_field_free


    function data_out_patch_field_get_name(this) result(name)
    !-----------------------------------------------------------------
    !< Return the name of the field 
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), intent(inout) :: this
        character(len=:), allocatable                :: name
    !-----------------------------------------------------------------
        assert(allocated(name))
        name = this%name
    end function data_out_patch_field_get_name


    function data_out_patch_field_get_field_id(this) result(field_id)
    !-----------------------------------------------------------------
    !< Return the associated field_id
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), intent(inout) :: this
        integer(ip)                                  :: field_id
    !-----------------------------------------------------------------
        assert(this%field_id > 0)
        field_id = this%field_id
    end function data_out_patch_field_get_field_id


    function data_out_patch_field_get_fe_function(this) result(fe_function)
    !-----------------------------------------------------------------
    !< Return a pointer to the associated fe_function
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), intent(inout) :: this
        type(fe_function_t), pointer                 :: fe_function
    !-----------------------------------------------------------------
        assert(associated(this%fe_function))
        fe_function => this%fe_function
    end function data_out_patch_field_get_fe_function


    function data_out_patch_field_get_nodal_values(this) result(nodal_values)
    !-----------------------------------------------------------------
    !< Return a pointer to nodal_values
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_rp1_t), pointer               :: nodal_values
    !-----------------------------------------------------------------
        nodal_values => this%nodal_values
    end function data_out_patch_field_get_nodal_values


    function data_out_patch_field_get_scalar_function_values(this) result(scalar_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to scalar_function_values
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_rp1_t), pointer               :: scalar_function_values
    !-----------------------------------------------------------------
        scalar_function_values => this%scalar_function_values
    end function data_out_patch_field_get_scalar_function_values


    function data_out_patch_field_get_scalar_function_gradients(this) result(scalar_function_gradients)
    !-----------------------------------------------------------------
    !< Return a pointer to scalar_function_gradients
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_vector_field_t), pointer      :: scalar_function_gradients
    !-----------------------------------------------------------------
        scalar_function_gradients => this%scalar_function_gradients
    end function data_out_patch_field_get_scalar_function_gradients


    function data_out_patch_field_get_vector_function_values(this) result(vector_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_values
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_vector_field_t), pointer      :: vector_function_values
    !-----------------------------------------------------------------
        vector_function_values => this%vector_function_values
    end function data_out_patch_field_get_vector_function_values


    function data_out_patch_field_get_vector_function_gradients(this) result(vector_function_gradients)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_gradients
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_tensor_field_t), pointer      :: vector_function_gradients
    !-----------------------------------------------------------------
        vector_function_gradients => this%vector_function_gradients
    end function data_out_patch_field_get_vector_function_gradients


    function data_out_patch_field_get_tensor_function_values(this) result(tensor_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_gradients
    !-----------------------------------------------------------------
        class(data_out_patch_field_t), target, intent(inout) :: this
        type(allocatable_array_tensor_field_t), pointer      :: tensor_function_values
    !-----------------------------------------------------------------
        tensor_function_values => this%tensor_function_values
    end function data_out_patch_field_get_tensor_function_values

!---------------------------------------------------------------------
!< DATA_OUT_PATCH_T PROCEDURES
!---------------------------------------------------------------------

    subroutine data_out_patch_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        integer(ip)                            :: i
    !-----------------------------------------------------------------
        nullify(this%coordinates)
        call this%subcells_connectivity%free()
        if(allocated(this%fields)) then
            do i=1, this%number_fields
                call this%fields(i)%free()
            enddo
            deallocate(this%fields)
        endif
        this%number_points = 0
        this%number_fields = 0
        this%number_subcells = 0
        this%number_vertices_per_subcell = 0
    end subroutine data_out_patch_free


    subroutine data_out_patch_set_number_dimensions(this, number_dimensions)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        integer(ip),             intent(in)    :: number_dimensions
    !-----------------------------------------------------------------
        this%number_dimensions = number_dimensions
    end subroutine data_out_patch_set_number_dimensions


    subroutine data_out_patch_set_number_points(this, number_points)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        integer(ip),             intent(in)    :: number_points
    !-----------------------------------------------------------------
        this%number_points = number_points
    end subroutine data_out_patch_set_number_points


    subroutine data_out_patch_set_number_subcells(this, number_subcells)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        integer(ip),             intent(in)    :: number_subcells
    !-----------------------------------------------------------------
        this%number_subcells = number_subcells
    end subroutine data_out_patch_set_number_subcells


    subroutine data_out_patch_set_number_vertices_per_subcell(this, number_vertices_per_subcell)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        integer(ip),             intent(in)    :: number_vertices_per_subcell
    !-----------------------------------------------------------------
        this%number_vertices_per_subcell = number_vertices_per_subcell
    end subroutine data_out_patch_set_number_vertices_per_subcell


    subroutine data_out_patch_set_coordinates(this, coordinates)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        type(point_t), pointer,  intent(in)    :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(coordinates))
        this%coordinates => coordinates
    end subroutine data_out_patch_set_coordinates


    function data_out_patch_get_number_dimensions(this) result(number_dimensions)
    !-----------------------------------------------------------------
    !< Return the number of dimensions of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(in) :: this
        integer(ip)                         :: number_dimensions
    !-----------------------------------------------------------------
        number_dimensions = this%number_dimensions
    end function data_out_patch_get_number_dimensions


    function data_out_patch_get_number_subcells(this) result(number_subcells)
    !-----------------------------------------------------------------
    !< Return the number of subcells in the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(in) :: this
        integer(ip)                         :: number_subcells
    !-----------------------------------------------------------------
        number_subcells = this%number_subcells
    end function data_out_patch_get_number_subcells


    function data_out_patch_get_number_vertices_per_subcell(this) result(number_vertices_per_subcell)
    !-----------------------------------------------------------------
    !< Return the number of subcells in the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(in) :: this
        integer(ip)                         :: number_vertices_per_subcell
    !-----------------------------------------------------------------
        number_vertices_per_subcell = this%number_vertices_per_subcell
    end function data_out_patch_get_number_vertices_per_subcell


    function data_out_patch_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields handled by the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(in) :: this
        integer(ip)                         :: number_fields
    !-----------------------------------------------------------------
        number_fields = this%number_fields
    end function data_out_patch_get_number_fields


    function data_out_patch_get_field(this, number_field) result(field)
    !-----------------------------------------------------------------
    !< Return the number of fields handled by the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t),      target, intent(in) :: this
        integer(ip),                          intent(in) :: number_field
        type(data_out_patch_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(number_field <= this%number_fields)
        field => this%fields(number_field)
    end function data_out_patch_get_field


    function data_out_patch_get_subcells_connectivity(this) result(subcells_connectivity)
    !-----------------------------------------------------------------
    !< Return a pointer to subcells connectivity 
    !-----------------------------------------------------------------
        class(data_out_patch_t),       target, intent(in) :: this
        type(allocatable_array_ip2_t), pointer            :: subcells_connectivity
    !-----------------------------------------------------------------
        subcells_connectivity => this%subcells_connectivity
    end function data_out_patch_get_subcells_connectivity


    function data_out_patch_get_coordinates(this) result(coordinates)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(inout) :: this
        type(point_t), pointer                 :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(this%coordinates))
        coordinates => this%coordinates
    end function data_out_patch_get_coordinates


    function data_out_patch_get_subcells_iterator(this) result(iterator)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(data_out_patch_t), intent(in)     :: this
        type(data_out_patch_subcell_iterator_t) :: iterator
    !-----------------------------------------------------------------
        call iterator%create(this)
    end function data_out_patch_get_subcells_iterator

!---------------------------------------------------------------------
!< DATA_OUT_PATCH_SUBCELL_ITERATOR_T PROCEDURES
!---------------------------------------------------------------------

    subroutine data_out_patch_subcell_iterator_free(this)
    !-----------------------------------------------------------------
    !< Free a patch subcell iterator
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        nullify(this%patch)
        this%current_subcell = 0
        this%number_subcells = 0
    end subroutine data_out_patch_subcell_iterator_free


    subroutine data_out_patch_subcell_iterator_create(this, patch)
    !-----------------------------------------------------------------
    !< Create a patch subcell iterator
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(inout) :: this
        type(data_out_patch_t), target,           intent(in)    :: patch
    !-----------------------------------------------------------------
        this%patch => patch
        this%number_subcells = this%patch%get_number_subcells()
        call this%begin()
    end subroutine data_out_patch_subcell_iterator_create


    subroutine data_out_patch_subcell_iterator_begin(this)
    !-----------------------------------------------------------------
    !< Rewind patch subcell iterator
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = 1
    end subroutine data_out_patch_subcell_iterator_begin


    subroutine data_out_patch_subcell_iterator_next(this)
    !-----------------------------------------------------------------
    !< Next subcell
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = this%current_subcell+1
    end subroutine data_out_patch_subcell_iterator_next


    function data_out_patch_subcell_iterator_has_finished(this) result(has_finished)
    !-----------------------------------------------------------------
    !< Rewind patch subcell iterator
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(inout) :: this
        logical                                                 :: has_finished
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        has_finished = this%current_subcell > this%number_subcells
    end function data_out_patch_subcell_iterator_has_finished


    subroutine data_out_patch_subcell_iterator_get_coordinates(this, subcell_coordinates)
    !-----------------------------------------------------------------
    !< Return subcell coordinates
    !-----------------------------------------------------------------
        class(data_out_patch_subcell_iterator_t), intent(in)    :: this
        real(rp), allocatable,                    intent(inout) :: subcell_coordinates(:,:)
        type(point_t),                 pointer                  :: patch_coordinates(:)
        type(allocatable_array_ip2_t), pointer                  :: subcells_connectivity
        integer(ip)                                             :: dim
        integer(ip)                                             :: vertex
    !-----------------------------------------------------------------
        if(.not. allocated(subcell_coordinates)) then
            call memalloc(this%patch%get_number_dimensions(),this%patch%get_number_vertices_per_subcell(), subcell_coordinates, __FILE__, __LINE__)
        elseif(size(subcell_coordinates,1) /= this%patch%get_number_dimensions() .or. &
               size(subcell_coordinates,2) /= this%patch%get_number_vertices_per_subcell()) then
            call memfree(subcell_coordinates, __FILE__, __LINE__)
            call memalloc(this%patch%get_number_dimensions(),this%patch%get_number_vertices_per_subcell(), subcell_coordinates, __FILE__, __LINE__)
        endif
        patch_coordinates     => this%patch%get_coordinates()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        do vertex = 1, this%patch%get_number_vertices_per_subcell()
            do dim = 1, this%patch%get_number_dimensions()
                subcell_coordinates(dim, vertex) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(dim)
            enddo
        end do
    end subroutine data_out_patch_subcell_iterator_get_coordinates

end module data_out_patch_names
