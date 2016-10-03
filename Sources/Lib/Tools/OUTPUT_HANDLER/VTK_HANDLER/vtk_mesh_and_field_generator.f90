
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

module vtk_mesh_and_field_generator

USE types_names
USE memor_names
USE IR_Precision,                    only: I1P
USE list_types_names
USE iso_fortran_env,                 only: error_unit
USE field_names,                     only: point_t, vector_field_t, tensor_field_t
USE vector_names,                    only: vector_t
USE base_static_triangulation_names, only: base_static_triangulation_t, cell_iterator_t, cell_accessor_t
USE serial_scalar_array_names,       only: serial_scalar_array_t
USE fe_space_names,                  only: serial_fe_space_t, fe_iterator_t, fe_accessor_t
USE fe_function_names,               only: fe_function_t
USE output_handler_cell_fe_function_names
USE reference_fe_names,              only: reference_fe_t, lagrangian_reference_fe_t, fe_map_t,   &
                                           quadrature_t, interpolation_t, topology_hex, topology_tet, &
                                           fe_type_lagrangian, field_type_scalar, field_type_vector, field_type_tensor
USE vtk_parameters_names

implicit none

#include "debug.i90"

private

    ! Type for storing mesh data
    type vtk_mesh_and_field_generator_t
    private
        class(serial_fe_space_t), pointer :: fe_space      => NULL()                ! Points to fe_space_t
        real(rp),          allocatable    :: X(:)                                   ! Mesh X coordintates
        real(rp),          allocatable    :: Y(:)                                   ! Mesh Y coordintates
        real(rp),          allocatable    :: Z(:)                                   ! Mesh Z coordintates
        integer(ip),       allocatable    :: connectivities(:)                      ! Connectivity matrix
        integer(ip),       allocatable    :: offset(:)                              ! VTK element offset
        integer(I1P),      allocatable    :: cell_types(:)                          ! VTK element type
        integer(ip)                       :: number_of_nodes    = 0                 ! Number of nodes
        integer(ip)                       :: number_of_elements = 0                 ! Number of elements
        integer(ip)                       :: dimensions         = 0                 ! Dimensions of the mesh
        logical                           :: filled             = .false.           ! Mesh data was already filled
    contains
    private
        procedure, non_overridable, public :: set_fe_space                      => vtk_mesh_and_field_generator_set_fe_space
        procedure, non_overridable, public :: get_number_nodes                  => vtk_mesh_and_field_generator_get_number_nodes
        procedure, non_overridable, public :: get_number_elements               => vtk_mesh_and_field_generator_get_number_elements
        procedure, non_overridable, public :: get_number_fields                 => vtk_mesh_and_field_generator_get_number_fields
        procedure, non_overridable, public :: get_X_coordinates                 => vtk_mesh_and_field_generator_get_X_coordinates
        procedure, non_overridable, public :: get_Y_coordinates                 => vtk_mesh_and_field_generator_get_Y_coordinates
        procedure, non_overridable, public :: get_Z_coordinates                 => vtk_mesh_and_field_generator_get_Z_coordinates
        procedure, non_overridable, public :: get_connectivities                => vtk_mesh_and_field_generator_get_connectivities
        procedure, non_overridable, public :: get_offset                        => vtk_mesh_and_field_generator_get_offset
        procedure, non_overridable, public :: get_cell_types                    => vtk_mesh_and_field_generator_get_cell_types
        procedure, non_overridable, public :: generate_mesh                     => vtk_mesh_and_field_generator_generate_mesh
        procedure, non_overridable, public :: generate_field                    => vtk_mesh_and_field_generator_generate_field
        procedure, non_overridable         :: initialize_coordinates            => vtk_mesh_and_field_generator_initialize_coordinates
        procedure, non_overridable         :: allocate_nodal_arrays             => vtk_mesh_and_field_generator_allocate_nodal_arrays
        procedure, non_overridable         :: allocate_elemental_arrays         => vtk_mesh_and_field_generator_allocate_elemental_arrays
        procedure, non_overridable         :: topology_to_cell_type             => vtk_mesh_and_field_generator_topology_to_cell_type
        procedure, non_overridable         :: generate_max_order_mesh           => vtk_mesh_and_field_generator_generate_max_order_mesh
        procedure, non_overridable, public :: free                              => vtk_mesh_and_field_generator_free
    end type vtk_mesh_and_field_generator_t

public :: vtk_mesh_and_field_generator_t

contains

    subroutine vtk_mesh_and_field_generator_set_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Set the fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),                intent(INOUT) :: this
        class(serial_fe_space_t), target, intent(IN)    :: fe_space
    !-----------------------------------------------------------------
        this%fe_space => fe_space
    end subroutine vtk_mesh_and_field_generator_set_fe_space


    function vtk_mesh_and_field_generator_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Get the fe_space pointer
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),                intent(IN) :: this
        class(serial_fe_space_t), pointer            :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function vtk_mesh_and_field_generator_get_fe_space


    function vtk_mesh_and_field_generator_is_filled(this) result(filled)
    !-----------------------------------------------------------------
    !< Ask if the mesh data is filled
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(IN) :: this
        logical                       :: filled
    !-----------------------------------------------------------------
        filled = this%filled
    end function vtk_mesh_and_field_generator_is_filled


    function vtk_mesh_and_field_generator_get_number_nodes(this) result(number_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        integer(ip)                      :: number_nodes
    !-----------------------------------------------------------------
        assert(this%filled)
        number_nodes = this%number_of_nodes
    end function vtk_mesh_and_field_generator_get_number_nodes


    function vtk_mesh_and_field_generator_get_number_elements(this) result(number_elements)
    !-----------------------------------------------------------------
    !< Return the number of elements
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        integer(ip)                      :: number_elements
    !-----------------------------------------------------------------
        assert(this%filled)
        number_elements = this%number_of_elements
    end function vtk_mesh_and_field_generator_get_number_elements


    function vtk_mesh_and_field_generator_get_dimensions(this) result(dimensions)
    !-----------------------------------------------------------------
    !< Return the space dimensions of the mesh
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        integer(ip)                      :: dimensions
    !-----------------------------------------------------------------
        assert(this%filled)
        dimensions = this%dimensions
    end function vtk_mesh_and_field_generator_get_dimensions


    function vtk_mesh_and_field_generator_get_X_coordinates(this) result(X)
    !-----------------------------------------------------------------
    !< Return a pointer to the X coordinates array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        real(rp),          pointer               :: X(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        X => this%X
    end function vtk_mesh_and_field_generator_get_X_coordinates


    function vtk_mesh_and_field_generator_get_Y_coordinates(this) result(Y)
    !-----------------------------------------------------------------
    !< Return a pointer to the Y coordinates array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        real(rp),          pointer               :: Y(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        Y => this%Y
    end function vtk_mesh_and_field_generator_get_Y_coordinates


    function vtk_mesh_and_field_generator_get_Z_coordinates(this) result(Z)
    !-----------------------------------------------------------------
    !< Return a pointer to the Z coordinates array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        real(rp),          pointer               :: Z(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        Z => this%Z
    end function vtk_mesh_and_field_generator_get_Z_coordinates


    function vtk_mesh_and_field_generator_get_connectivities(this) result(connectivities)
    !-----------------------------------------------------------------
    !< Return a pointer to the connectivities array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        integer(ip),       pointer               :: connectivities(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        connectivities => this%connectivities
    end function vtk_mesh_and_field_generator_get_connectivities


    function vtk_mesh_and_field_generator_get_offset(this) result(offset)
    !-----------------------------------------------------------------
    !< Return a pointer to the offset array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        integer(ip),       pointer               :: offset(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        offset => this%offset
    end function vtk_mesh_and_field_generator_get_offset


    function vtk_mesh_and_field_generator_get_cell_types(this) result(cell_types)
    !-----------------------------------------------------------------
    !< Return a pointer to the cell_types array
    !< Generate the mesh if it is not filled yet
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), target, intent(INOUT) :: this
        integer(I1P),      pointer               :: cell_types(:)
    !-----------------------------------------------------------------
        assert(this%filled)
        cell_types => this%cell_types
    end function vtk_mesh_and_field_generator_get_cell_types


    subroutine vtk_mesh_and_field_generator_initialize_coordinates(this)
    !-----------------------------------------------------------------
    !< Set the z coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(allocated(this%X))
        assert(allocated(this%Y))
        assert(allocated(this%Z))
        this%X = 0.0_rp
        this%Y = 0.0_rp
        this%Z = 0.0_rp
    end subroutine vtk_mesh_and_field_generator_initialize_coordinates


    function vtk_mesh_and_field_generator_get_number_fields(this) result(number_fields)
    !-----------------------------------------------------------------
    !< Return the number of fields 
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(IN) :: this
        integer(ip)                   :: number_fields
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))
        number_fields =  this%fe_space%get_number_fields()
    end function vtk_mesh_and_field_generator_get_number_fields


    subroutine vtk_mesh_and_field_generator_allocate_elemental_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate all arrays of size number of elements
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(.not. allocated(this%offset))
        assert(.not. allocated(this%cell_types))
        call memalloc(this%number_of_elements, this%offset, __FILE__, __LINE__)
        call memalloc(this%number_of_elements, this%cell_types, __FILE__, __LINE__)
    end subroutine vtk_mesh_and_field_generator_allocate_elemental_arrays


    subroutine vtk_mesh_and_field_generator_allocate_nodal_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate all arrays with size number of nodes
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(.not. allocated(this%connectivities))
        assert(.not. allocated(this%X))
        assert(.not. allocated(this%Y))
        assert(.not. allocated(this%Z))
        call memalloc (this%number_of_nodes, this%connectivities, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%X, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%Y, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%Z, __FILE__,__LINE__)
    end subroutine vtk_mesh_and_field_generator_allocate_nodal_arrays


    function vtk_mesh_and_field_generator_topology_to_cell_type(this, topology, dimension) result(cell_type)
    !-----------------------------------------------------------------
    !< Translate the topology type of the reference_fe_geo into VTK cell type
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        character(len=*),                      intent(IN)    :: topology
        integer(ip),                           intent(IN)    :: dimension
        integer(I1P)                                         :: cell_type
    !-----------------------------------------------------------------
        if(topology == topology_hex) then 
            if(dimension == 2) then
                cell_type = vtk_pixel
            elseif(dimension == 3) then
                cell_type = vtk_voxel
            endif
        elseif(topology == topology_tet) then
            if(dimension == 2) then
                cell_type = vtk_triangle
            elseif(dimension == 3) then
                cell_type = vtk_tetra
            endif
        else
            write(error_unit,*) 'fill_mesh_from_triangulation: Topology not supported'
            check(.false.)    
        endif
    end function vtk_mesh_and_field_generator_topology_to_cell_type


    subroutine vtk_mesh_and_field_generator_generate_mesh(this, fe_function)
    !-----------------------------------------------------------------
    !< Generate the mesh data from fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        type(fe_function_t),                   intent(IN)    :: fe_function
    !-----------------------------------------------------------------
        if(.not. this%filled) then
            call this%generate_max_order_mesh(fe_function)
        endif
    end subroutine vtk_mesh_and_field_generator_generate_mesh


    subroutine vtk_mesh_and_field_generator_generate_max_order_mesh(this, fe_function)
    !-----------------------------------------------------------------
    !< Store a superlinear_order mesh in a from a fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        type(fe_function_t),                   intent(IN)    :: fe_function
        type(point_t),                    pointer            :: subcells_coordinates(:)
        integer(ip),                      pointer            :: subcells_connectivity(:,:)
        type(point_t),                    pointer            :: quadrature_coordinates(:)
        class(reference_fe_t),            pointer            :: reference_fe
        class(lagrangian_reference_fe_t), pointer            :: reference_fe_geo
        type(fe_iterator_t)                                  :: fe_iterator
        type(fe_accessor_t)                                  :: fe
        type(output_handler_cell_fe_function_t)              :: output_cell_handler
        integer(ip)                                          :: num_refinements
        integer(ip)                                          :: num_elements
        integer(ip)                                          :: num_nodes_per_element
        integer(ip)                                          :: num_vertices_per_element
        integer(ip)                                          :: num_subelements_per_element
        integer(ip)                                          :: elements_counter
        integer(ip)                                          :: vertex
        integer(ip)                                          :: reference_fe_id
        integer(ip)                                          :: dimensions
        integer(ip)                                          :: subelement_vertex
        integer(ip)                                          :: subelement_index
        integer(ip)                                          :: nodes_counter
        integer(ip)                                          :: istat
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))

        nullify(quadrature_coordinates)
        nullify(reference_fe)

        ! Create FE iterator 
        fe_iterator             = this%fe_space%create_fe_iterator()
        this%number_of_nodes    = 0
        this%number_of_elements = 0
        ! Count number of vtk elements and nodes
        do while ( .not. fe_iterator%has_finished())
            ! Get Finite element
            call fe_iterator%current(fe)
            if ( fe%is_local() ) then
                ! Create FE_MAP for current cell
                reference_fe            => fe%get_max_order_reference_fe()
                reference_fe_geo        => fe%get_reference_fe_geo()
                num_refinements         = reference_fe%get_max_order()-1
                this%number_of_elements = this%number_of_elements+reference_fe_geo%get_number_subcells(num_refinements)
                this%number_of_nodes    = this%number_of_nodes+(reference_fe_geo%get_number_subcells(num_refinements)*reference_fe_geo%get_number_vertices())
            endif
            call fe_iterator%next()
        enddo

        ! Allocate VTK geometry and connectivity data
        call this%allocate_elemental_arrays()
        call this%allocate_nodal_arrays()
        call this%initialize_coordinates()

        nodes_counter    = 0
        elements_counter = 0
        ! Create FE iterator
        fe_iterator      = this%fe_space%create_fe_iterator()
        ! Create Output Cell Handler
        call output_cell_handler%create(this%fe_space)
        ! Translate coordinates and connectivities to VTK format for every subcell
        do while ( .not. fe_iterator%has_finished())
            ! Get Finite element
            call fe_iterator%current(fe)

            if ( fe%is_local() ) then

                call output_cell_handler%update(fe, fe_function)

                reference_fe            => fe%get_max_order_reference_fe()
                reference_fe_geo        => fe%get_reference_fe_geo()
                num_refinements         =  reference_fe%get_max_order()-1

                subcells_coordinates  => output_cell_handler%get_subcells_vertex_coordinates()
                subcells_connectivity => output_cell_handler%get_subcells_connectivity()

                ! Fill VTK mesh
                dimensions = reference_fe_geo%get_number_dimensions()
                do subelement_index = 1, reference_fe_geo%get_number_subcells(num_refinements)
                    elements_counter = elements_counter + 1
                    do vertex = 1, reference_fe%get_number_vertices()
                        subelement_vertex = subcells_connectivity(vertex, subelement_index)
                        nodes_counter = nodes_counter + 1
                        if(dimensions>=1) this%X(nodes_counter) = subcells_coordinates(subelement_vertex)%get(1)
                        if(dimensions>=2) this%Y(nodes_counter) = subcells_coordinates(subelement_vertex)%get(2)
                        if(dimensions>=3) this%Z(nodes_counter) = subcells_coordinates(subelement_vertex)%get(3)
                        this%connectivities(nodes_counter) = nodes_counter-1
                    end do

                    ! Store the type of element
                    this%cell_types(elements_counter) = this%topology_to_cell_type(reference_fe%get_topology(), dimensions)

                    ! Fill offset
                    this%offset(elements_counter) = nodes_counter
                end do
            endif
            call fe_iterator%next()
        end do
        this%filled  = .true.

        call output_cell_handler%free()
        call fe_iterator%free()
        call fe%free()
    end subroutine vtk_mesh_and_field_generator_generate_max_order_mesh


    function vtk_mesh_and_field_generator_generate_field(this, fe_function, field_id, field_name, field, number_components) result(E_IO)
    !-----------------------------------------------------------------
    !< Generate the field data from fe_space and fe_function
    !-----------------------------------------------------------------
        implicit none
        class(vtk_mesh_and_field_generator_t),  intent(INOUT) :: this         !< this raw mesh
        type(fe_function_t),        intent(IN)     :: fe_function              !< Postprocess field structure to be written
        integer(ip),                intent(IN)     :: field_id                 !< Fe space index
        character(len=*),           intent(IN)     :: field_name               !< name of the field
        real(rp),    allocatable,   intent(INOUT)  :: field(:,:)               !< FIELD(ncomp,nnod)
        integer(ip),                intent(OUT)    :: number_components        !< number of components
        type(output_handler_cell_fe_function_t)    :: output_cell_handler
        type(fe_iterator_t)                        :: fe_iterator              !< finite element iterator
        type(fe_accessor_t)                        :: fe                       !< finite element accessor
        type(base_static_triangulation_t), pointer :: triangulation            !< triangulation
        class(vector_t),                   pointer :: fe_function_dof_values   !< dof values of the fe_function
        type(serial_scalar_array_t),       pointer :: strong_dirichlet_values  !< Strong dirichlet values
        class(lagrangian_reference_fe_t),  pointer :: reference_fe_geo         !< reference finite element
        class(reference_fe_t),             pointer :: reference_fe             !< reference finite element
        integer(ip),                       pointer :: subcells_connectivity(:,:) !< scalar nodal values
        type(vector_field_t),              pointer :: vector_values(:) 
        type(tensor_field_t),              pointer :: tensor_values(:) 
        real(rp),                          pointer :: values(:)                !< scalar nodal values
        integer(ip)                                :: reference_fe_id          !< reference_fe_id
        integer(ip)                                :: element_index            !< element index
        integer(ip)                                :: component_index          !< component index
        integer(ip)                                :: node_index               !< node index
        integer(ip)                                :: subelement_index         !< subelement index
        integer(ip)                                :: subnode_index            !< subelement node index
        integer(ip)                                :: pos                      !< array position
        integer(ip)                                :: E_IO                     !< IO Error
    !-----------------------------------------------------------------
        E_IO = 0

        nullify(triangulation)
        assert(associated(this%fe_space))
        triangulation => this%fe_space%get_triangulation()
        assert(associated(triangulation))
        assert(this%filled)

        nullify(reference_fe_geo)
        nullify(reference_fe)
        
        ! Create FE iterator and get number of components
        fe_iterator = this%fe_space%create_fe_iterator()
        call fe_iterator%current(fe)
        reference_fe         => fe%get_reference_fe(field_id)
        number_components    = reference_fe%get_number_field_components()

        ! Get number of field components (constant in all fes for field) and allocate VTK field array
        if(allocated(field)) call memfree(field, __FILE__, __LINE__)
!        call memalloc(number_components, this%number_of_nodes , field, __FILE__, __LINE__)
		if (number_components > 1) then 
            call memalloc(3, this%number_of_nodes , field, __FILE__, __LINE__)
            number_components = 3
    		field = 0.0_rp
		else 
    		call memalloc(number_components, this%number_of_nodes , field, __FILE__, __LINE__)
		end if 

        ! Create Output Cell Handler
        call output_cell_handler%create(this%fe_space)
        ! Translate fe_function to VTK field format
        ! Loop on elements
        subnode_index=1
        do while ( .not. fe_iterator%has_finished())
            ! Get Finite element
            call fe_iterator%current(fe)

            if ( fe%is_local() ) then

                call output_cell_handler%update(fe, fe_function)

                ! Get reference fes
                reference_fe_geo    => fe%get_reference_fe_geo()
                reference_fe        => fe%get_max_order_reference_fe()

                subcells_connectivity => output_cell_handler%get_subcells_connectivity()
                select case(fe%get_field_type(field_id))
                    case ( field_type_scalar )
                        call output_cell_handler%get_values_scalar(field_id, values)
                    case ( field_type_vector )
                        call output_cell_handler%get_values_vector(field_id, vector_values)
                    case ( field_type_tensor )
                        call output_cell_handler%get_values_tensor(field_id, tensor_values)
                    case default
                        assert(.false.)
                end select

                ! Loop on subelements: Build field in VTK-like format
                do subelement_index = 1, reference_fe_geo%get_number_subcells(num_refinements=reference_fe%get_max_order()-1)
                    ! Loop on geometrical nodes per subelement
                    do node_index=1, reference_fe_geo%get_number_vertices()
                        ! Loop on components
                        do component_index=1, number_components
                            select case(fe%get_field_type(field_id))
                                case ( field_type_scalar )
                                    field(component_index , subnode_index) = values(subcells_connectivity(node_index,subelement_index))
                                case ( field_type_vector )
                                    field(component_index , subnode_index) = vector_values(subcells_connectivity(node_index,subelement_index))%get(component_index)
                                case ( field_type_tensor )
!                                    call output_cell_handler%get_values_tensor(field_id, tensor_values)
                                case default
                                    assert(.false.)
                            end select
                        enddo
                        subnode_index=subnode_index+1
                    end do
                end do
            endif

            call fe_iterator%next()
        enddo

        call output_cell_handler%free()
    end function vtk_mesh_and_field_generator_generate_field


    subroutine vtk_mesh_and_field_generator_free(this) 
    !-----------------------------------------------------------------
    !< Free the vtk_mesh_and_field_generator_t derived type
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(inout) :: this
        integer(ip)                                          :: i
    !-----------------------------------------------------------------
        if(allocated(this%X))                        call memfree(this%X, __FILE__, __LINE__)
        if(allocated(this%Y))                        call memfree(this%Y, __FILE__, __LINE__)
        if(allocated(this%Z))                        call memfree(this%Z, __FILE__, __LINE__)
        if(allocated(this%connectivities))           call memfree(this%connectivities, __FILE__, __LINE__)
        if(allocated(this%offset))                   call memfree(this%offset, __FILE__, __LINE__)
        if(allocated(this%cell_types))               call memfree(this%cell_types, __FILE__, __LINE__)
        nullify(this%fe_space)
        this%number_of_nodes    = 0
        this%number_of_elements = 0
        this%filled             = .false.
    end subroutine

end module vtk_mesh_and_field_generator
