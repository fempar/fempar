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

module output_handler_patch_names

USE types_names
USE field_names
USE allocatable_array_names
USE reference_fe_names, only: field_type_scalar, field_type_vector, field_type_tensor, field_type_symmetric_tensor

implicit none
#include "debug.i90"
private

    type :: output_handler_patch_field_t
    private
        character(len=:), allocatable          :: field_type
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
        procedure, public :: free                          => output_handler_patch_field_free
        procedure, public :: set_field_type                => output_handler_patch_field_set_field_type
        procedure, public :: get_field_type                => output_handler_patch_field_get_field_type
        procedure, public :: get_nodal_values              => output_handler_patch_field_get_nodal_values
        procedure, public :: get_scalar_function_values    => output_handler_patch_field_get_scalar_function_values
        procedure, public :: get_scalar_function_gradients => output_handler_patch_field_get_scalar_function_gradients
        procedure, public :: get_vector_function_values    => output_handler_patch_field_get_vector_function_values
        procedure, public :: get_vector_function_gradients => output_handler_patch_field_get_vector_function_gradients
        procedure, public :: get_tensor_function_values    => output_handler_patch_field_get_tensor_function_values
        procedure, public :: get_number_components         => output_handler_patch_field_get_number_components
        procedure, public :: get_raw_values                => output_handler_patch_field_get_raw_values
    end type

    type :: output_handler_patch_t
    private
        integer(ip)                                :: number_dimensions = 0
        integer(ip)                                :: number_points = 0
        integer(ip)                                :: number_fields = 0
        integer(ip)                                :: number_subcells = 0
        integer(ip)                                :: number_vertices_per_subcell = 0
        type(point_t),                 pointer     :: coordinates(:)
        type(allocatable_array_ip2_t)              :: subcells_connectivity
        type(output_handler_patch_field_t),  public,allocatable :: fields(:)
    contains
        procedure, public :: set_number_dimensions           => output_handler_patch_set_number_dimensions
        procedure, public :: set_number_points               => output_handler_patch_set_number_points
        procedure, public :: set_number_subcells             => output_handler_patch_set_number_subcells
        procedure, public :: set_number_vertices_per_subcell => output_handler_patch_set_number_vertices_per_subcell
        procedure, public :: set_coordinates                 => output_handler_patch_set_coordinates
        procedure, public :: create                          => output_handler_patch_create
        procedure, public :: free                            => output_handler_patch_free
        procedure, public :: get_number_dimensions           => output_handler_patch_get_number_dimensions
        procedure, public :: get_subcells_connectivity       => output_handler_patch_get_subcells_connectivity
        procedure, public :: get_number_fields               => output_handler_patch_get_number_fields
        procedure, public :: get_number_subcells             => output_handler_patch_get_number_subcells
        procedure, public :: get_number_vertices_per_subcell => output_handler_patch_get_number_vertices_per_subcell
        procedure, public :: get_coordinates                 => output_handler_patch_get_coordinates
        procedure, public :: get_field                       => output_handler_patch_get_field
        procedure, public :: get_subcells_iterator           => output_handler_patch_get_subcells_iterator
    end type

    type :: output_handler_patch_subcell_iterator_t
    private
        type(output_handler_patch_t), pointer :: patch => NULL()
        integer(ip)                     :: current_subcell = 0
        integer(ip)                     :: number_subcells = 0
    contains
        procedure, public :: create            => output_handler_patch_subcell_iterator_create
        procedure, public :: begin             => output_handler_patch_subcell_iterator_begin
        procedure, public :: next              => output_handler_patch_subcell_iterator_next
        procedure, public :: has_finished      => output_handler_patch_subcell_iterator_has_finished
        procedure, public :: free              => output_handler_patch_subcell_iterator_free
        procedure, public :: get_coordinates   => output_handler_patch_subcell_iterator_get_coordinates
        procedure, public :: get_connectivity  => output_handler_patch_subcell_iterator_get_connectivity
        procedure, public :: get_number_fields => output_handler_patch_subcell_iterator_get_number_fields
        procedure, public :: get_field         => output_handler_patch_subcell_iterator_get_field
    end type

public :: output_handler_patch_t
public :: output_handler_patch_field_t
public :: output_handler_patch_subcell_iterator_t

contains

!---------------------------------------------------------------------
!< output_handler_PATCH_FIELD_T PROCEDURES
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
        call this%scalar_function_gradients%free()
        call this%vector_function_values%free()
        call this%vector_function_gradients%free()
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


    function output_handler_patch_field_get_scalar_function_gradients(this) result(scalar_function_gradients)
    !-----------------------------------------------------------------
    !< Return a pointer to scalar_function_gradients
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_vector_field_t), pointer               :: scalar_function_gradients
    !-----------------------------------------------------------------
        scalar_function_gradients => this%scalar_function_gradients
    end function output_handler_patch_field_get_scalar_function_gradients


    function output_handler_patch_field_get_vector_function_values(this) result(vector_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_values
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_vector_field_t), pointer               :: vector_function_values
    !-----------------------------------------------------------------
        vector_function_values => this%vector_function_values
    end function output_handler_patch_field_get_vector_function_values


    function output_handler_patch_field_get_vector_function_gradients(this) result(vector_function_gradients)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_gradients
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_tensor_field_t), pointer               :: vector_function_gradients
    !-----------------------------------------------------------------
        vector_function_gradients => this%vector_function_gradients
    end function output_handler_patch_field_get_vector_function_gradients


    function output_handler_patch_field_get_tensor_function_values(this) result(tensor_function_values)
    !-----------------------------------------------------------------
    !< Return a pointer to vector_function_gradients
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t),    target, intent(inout) :: this
        type(allocatable_array_tensor_field_t), pointer               :: tensor_function_values
    !-----------------------------------------------------------------
        tensor_function_values => this%tensor_function_values
    end function output_handler_patch_field_get_tensor_function_values


    function output_handler_patch_field_get_number_components(this) result(number_components)
    !-----------------------------------------------------------------
    !< Return the number of components of the field
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(in) :: this
        integer(ip)                                     :: number_components       
    !-----------------------------------------------------------------
        select case(this%field_type)
            case (field_type_scalar)
                number_components = 1
            case (field_type_vector)
                number_components = SPACE_DIM
            case (field_type_tensor)
                number_components = SPACE_DIM*SPACE_DIM
            case (field_type_symmetric_tensor)
                number_components = SPACE_DIM*SPACE_DIM
            case DEFAULT
                check(.false.)
        end select
    end function output_handler_patch_field_get_number_components


    subroutine output_handler_patch_field_get_raw_values(this, number_vertices, raw_values)
    !-----------------------------------------------------------------
    !< Return the number of components of the field
    !-----------------------------------------------------------------
        class(output_handler_patch_field_t), intent(in)    :: this
        integer(ip),                         intent(in)    :: number_vertices
        real(rp),                            intent(inout) :: raw_values(:,:)
        integer(ip)                                        :: i, j, k
        type(vector_field_t), pointer                      :: vector_field(:)
        type(tensor_field_t), pointer                      :: tensor_field(:)
    !-----------------------------------------------------------------
        select case(this%field_type)
            case (field_type_scalar)
                raw_values(1,1:number_vertices) = this%scalar_function_values%a(1:number_vertices)
            case (field_type_vector)
                vector_field => this%vector_function_values%get_array()
                do i=1, number_vertices
                    do j=1, SPACE_DIM
                        raw_values(j,i) = vector_field(i)%get(j)
                    enddo
                enddo
            case (field_type_tensor, field_type_symmetric_tensor)
                tensor_field => this%tensor_function_values%get_array()
                do i=1, number_vertices
                    do j=1, SPACE_DIM
                        do k=1, SPACE_DIM
                            raw_values((j-1*SPACE_DIM)+k,i) = tensor_field(i)%get(j,k)
                        enddo
                    enddo
                enddo
            case DEFAULT
                check(.false.)
        end select
    end subroutine output_handler_patch_field_get_raw_values


!---------------------------------------------------------------------
!< output_handler_PATCH_T PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_patch_create(this, number_fields)
    !-----------------------------------------------------------------
    !< Create procedure. Allocate fields
    !-----------------------------------------------------------------
        class(output_handler_patch_t),      intent(inout) :: this
        integer(ip),                        intent(in)    :: number_fields
    !-----------------------------------------------------------------
        call this%free()
        this%number_fields = number_fields
        allocate(this%fields(this%number_fields))
    end subroutine output_handler_patch_create


    subroutine output_handler_patch_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip)                                  :: i
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
    end subroutine output_handler_patch_free


    subroutine output_handler_patch_set_number_dimensions(this, number_dimensions)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: number_dimensions
    !-----------------------------------------------------------------
        this%number_dimensions = number_dimensions
    end subroutine output_handler_patch_set_number_dimensions


    subroutine output_handler_patch_set_number_points(this, number_points)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: number_points
    !-----------------------------------------------------------------
        this%number_points = number_points
    end subroutine output_handler_patch_set_number_points


    subroutine output_handler_patch_set_number_subcells(this, number_subcells)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: number_subcells
    !-----------------------------------------------------------------
        this%number_subcells = number_subcells
    end subroutine output_handler_patch_set_number_subcells


    subroutine output_handler_patch_set_number_vertices_per_subcell(this, number_vertices_per_subcell)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        integer(ip),                   intent(in)    :: number_vertices_per_subcell
    !-----------------------------------------------------------------
        this%number_vertices_per_subcell = number_vertices_per_subcell
    end subroutine output_handler_patch_set_number_vertices_per_subcell


    subroutine output_handler_patch_set_coordinates(this, coordinates)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        type(point_t), pointer,        intent(in)    :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(coordinates))
        this%coordinates => coordinates
    end subroutine output_handler_patch_set_coordinates


    function output_handler_patch_get_number_dimensions(this) result(number_dimensions)
    !-----------------------------------------------------------------
    !< Return the number of dimensions of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: number_dimensions
    !-----------------------------------------------------------------
        number_dimensions = this%number_dimensions
    end function output_handler_patch_get_number_dimensions


    function output_handler_patch_get_number_subcells(this) result(number_subcells)
    !-----------------------------------------------------------------
    !< Return the number of subcells in the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: number_subcells
    !-----------------------------------------------------------------
        number_subcells = this%number_subcells
    end function output_handler_patch_get_number_subcells


    function output_handler_patch_get_number_vertices_per_subcell(this) result(number_vertices_per_subcell)
    !-----------------------------------------------------------------
    !< Return the number of subcells in the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: number_vertices_per_subcell
    !-----------------------------------------------------------------
        number_vertices_per_subcell = this%number_vertices_per_subcell
    end function output_handler_patch_get_number_vertices_per_subcell


    function output_handler_patch_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields handled by the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in) :: this
        integer(ip)                               :: number_fields
    !-----------------------------------------------------------------
        number_fields = this%number_fields
    end function output_handler_patch_get_number_fields


    function output_handler_patch_get_field(this, number_field) result(field)
    !-----------------------------------------------------------------
    !< Return the number of fields handled by the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t),      target, intent(in) :: this
        integer(ip),                                intent(in) :: number_field
        type(output_handler_patch_field_t), pointer            :: field
    !-----------------------------------------------------------------
        assert(number_field <= this%number_fields)
        field => this%fields(number_field)
    end function output_handler_patch_get_field


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
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(inout) :: this
        type(point_t), pointer                       :: coordinates(:)
    !-----------------------------------------------------------------
        assert(associated(this%coordinates))
        coordinates => this%coordinates
    end function output_handler_patch_get_coordinates


    function output_handler_patch_get_subcells_iterator(this) result(iterator)
    !-----------------------------------------------------------------
    !< Set the number of points of the patch
    !-----------------------------------------------------------------
        class(output_handler_patch_t), intent(in)     :: this
        type(output_handler_patch_subcell_iterator_t) :: iterator
    !-----------------------------------------------------------------
        call iterator%create(this)
    end function output_handler_patch_get_subcells_iterator

!---------------------------------------------------------------------
!< output_handler_PATCH_SUBCELL_ITERATOR_T PROCEDURES
!---------------------------------------------------------------------

    subroutine output_handler_patch_subcell_iterator_free(this)
    !-----------------------------------------------------------------
    !< Free a patch subcell iterator
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        nullify(this%patch)
        this%current_subcell = 0
        this%number_subcells = 0
    end subroutine output_handler_patch_subcell_iterator_free


    subroutine output_handler_patch_subcell_iterator_create(this, patch)
    !-----------------------------------------------------------------
    !< Create a patch subcell iterator
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(inout) :: this
        type(output_handler_patch_t), target,           intent(in)    :: patch
    !-----------------------------------------------------------------
        this%patch => patch
        this%number_subcells = this%patch%get_number_subcells()
        call this%begin()
    end subroutine output_handler_patch_subcell_iterator_create


    subroutine output_handler_patch_subcell_iterator_begin(this)
    !-----------------------------------------------------------------
    !< Rewind patch subcell iterator
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = 1
    end subroutine output_handler_patch_subcell_iterator_begin


    subroutine output_handler_patch_subcell_iterator_next(this)
    !-----------------------------------------------------------------
    !< Next subcell
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        this%current_subcell = this%current_subcell+1
    end subroutine output_handler_patch_subcell_iterator_next


    function output_handler_patch_subcell_iterator_has_finished(this) result(has_finished)
    !-----------------------------------------------------------------
    !< Rewind patch subcell iterator
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(inout) :: this
        logical                                                       :: has_finished
    !-----------------------------------------------------------------
        assert(associated(this%patch))
        has_finished = this%current_subcell > this%number_subcells
    end function output_handler_patch_subcell_iterator_has_finished


    subroutine output_handler_patch_subcell_iterator_get_coordinates(this, subcell_coordinates)
    !-----------------------------------------------------------------
    !< Return subcell coordinates
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(in)    :: this
        real(rp), allocatable,                          intent(inout) :: subcell_coordinates(:,:)
        type(point_t),                 pointer                        :: patch_coordinates(:)
        type(allocatable_array_ip2_t), pointer                        :: subcells_connectivity
        integer(ip)                                                   :: number_dimensions
        integer(ip)                                                   :: number_vertices
        integer(ip)                                                   :: dim
        integer(ip)                                                   :: vertex
    !-----------------------------------------------------------------
        number_vertices   = this%patch%get_number_vertices_per_subcell()
        number_dimensions = this%patch%get_number_dimensions()
        if(.not. allocated(subcell_coordinates)) then
            call memalloc(number_dimensions, number_vertices, subcell_coordinates, __FILE__, __LINE__)
        elseif(size(subcell_coordinates,1) /= number_dimensions .or. size(subcell_coordinates,2) /= number_vertices) then
            call memfree(subcell_coordinates, __FILE__, __LINE__)
            call memalloc(number_dimensions, number_vertices, subcell_coordinates, __FILE__, __LINE__)
        endif
        patch_coordinates     => this%patch%get_coordinates()
        subcells_connectivity => this%patch%get_subcells_connectivity()
        do vertex = 1, number_vertices
            do dim = 1, number_dimensions
                subcell_coordinates(dim, vertex) = patch_coordinates(subcells_connectivity%a(vertex, this%current_subcell))%get(dim)
            enddo
        end do
    end subroutine output_handler_patch_subcell_iterator_get_coordinates


    subroutine output_handler_patch_subcell_iterator_get_connectivity(this, subcell_connectivity)
    !-----------------------------------------------------------------
    !< Return subcell connectivity
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(in)    :: this
        integer(ip), allocatable,                       intent(inout) :: subcell_connectivity(:)
        type(allocatable_array_ip2_t), pointer                        :: subcells_connectivity
        integer(ip)                                                   :: number_vertices
    !-----------------------------------------------------------------
        number_vertices = this%patch%get_number_vertices_per_subcell()
        if(.not. allocated(subcell_connectivity)) then
            call memalloc(number_vertices, subcell_connectivity, __FILE__, __LINE__)
        elseif(size(subcell_connectivity) /= number_vertices) then
            call memfree(subcell_connectivity, __FILE__, __LINE__)
            call memalloc(number_vertices, subcell_connectivity, __FILE__, __LINE__)
        endif
        subcells_connectivity => this%patch%get_subcells_connectivity()
        subcell_connectivity(1:number_vertices) = subcells_connectivity%a(1:number_vertices, this%current_subcell)
    end subroutine output_handler_patch_subcell_iterator_get_connectivity


    function output_handler_patch_subcell_iterator_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return subcell number of fields
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(in)    :: this
        integer(ip)                                                   :: number_fields
    !-----------------------------------------------------------------
        number_fields = this%patch%get_number_fields()
    end function output_handler_patch_subcell_iterator_get_number_fields


    subroutine output_handler_patch_subcell_iterator_get_field(this, field_id, field)
    !-----------------------------------------------------------------
    !< Return subcell field corresponding to the field_id
    !-----------------------------------------------------------------
        class(output_handler_patch_subcell_iterator_t), intent(in)    :: this
        integer(ip),                                    intent(in)    :: field_id
        real(rp), allocatable,                          intent(inout) :: field(:,:)
        type(output_handler_patch_field_t), pointer                   :: patch_field
        type(allocatable_array_ip2_t),      pointer                   :: subcells_connectivity
        integer(ip)                                                   :: number_components
        integer(ip)                                                   :: number_vertices
        integer(ip)                                                   :: vertex
        integer(ip)                                                   :: comp
    !-----------------------------------------------------------------
        patch_field       => this%patch%get_field(field_id)
        number_vertices   =  this%patch%get_number_vertices_per_subcell()
        number_components =  patch_field%get_number_components()
        if(.not. allocated(field)) then
            call memalloc(number_components, number_vertices, field, __FILE__, __LINE__)
        elseif(size(field,1) /= number_components .or. size(field,2) /= number_vertices) then
            call memfree(field, __FILE__, __LINE__)
            call memalloc(number_components, number_vertices, field, __FILE__, __LINE__)
        endif
        subcells_connectivity => this%patch%get_subcells_connectivity()
        call patch_field%get_raw_values(number_vertices, field)
    end subroutine output_handler_patch_subcell_iterator_get_field

end module output_handler_patch_names
