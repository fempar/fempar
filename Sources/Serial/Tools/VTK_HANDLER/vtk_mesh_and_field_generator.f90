
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
USE IR_Precision,                only: I1P
USE list_types_names
USE iso_fortran_env,             only: error_unit
USE field_names,                 only: point_t
USE vector_names,                only: vector_t
USE triangulation_names,         only: triangulation_t
USE serial_scalar_array_names,   only: serial_scalar_array_t
USE serial_fe_space_names,       only: serial_fe_space_t, finite_element_t, fe_function_t
USE reference_fe_names,          only: reference_fe_t, quad_lagrangian_reference_fe_t, fe_map_t,   &
                                       quadrature_t, interpolation_t, topology_quad, topology_tet, &
                                       fe_type_lagrangian

implicit none

#include "debug.i90"

private


    ! VTK cell type parameters
    integer(ip), parameter :: vtk_vertex               = 1_I1P
    integer(ip), parameter :: vtk_poly_vertex          = 2_I1P
    integer(ip), parameter :: vtk_line                 = 3_I1P
    integer(ip), parameter :: vtk_poly_line            = 4_I1P
    integer(ip), parameter :: vtk_triangle             = 5_I1P
    integer(ip), parameter :: vtk_triangle_strip       = 6_I1P
    integer(ip), parameter :: vtk_polygon              = 7_I1P
    integer(ip), parameter :: vtk_pixel                = 8_I1P
    integer(ip), parameter :: vtk_quad                 = 9_I1P
    integer(ip), parameter :: vtk_tetra                = 10_I1P
    integer(ip), parameter :: vtk_voxel                = 11_I1P
    integer(ip), parameter :: vtk_hexahedron           = 12_I1P
    integer(ip), parameter :: vtk_wedge                = 13_I1P
    integer(ip), parameter :: vtk_pyramid              = 14_I1P
    integer(ip), parameter :: vtk_quadratic_edge       = 21_I1P
    integer(ip), parameter :: vtk_quadratic_triangle   = 22_I1P
    integer(ip), parameter :: vtk_quadratic_quad       = 23_I1P
    integer(ip), parameter :: vtk_quadratic_tetra      = 24_I1P
    integer(ip), parameter :: vtk_quadratic_hexahedron = 25_I1P


    ! Type for storing mesh data
    type vtk_mesh_and_field_generator_t
    private
        class(serial_fe_space_t), pointer :: fe_space      => NULL()                ! Poins to fe_space_t
        real(rp),          allocatable    :: X(:)                                   ! Mesh X coordintates
        real(rp),          allocatable    :: Y(:)                                   ! Mesh Y coordintates
        real(rp),          allocatable    :: Z(:)                                   ! Mesh Z coordintates
        integer(ip),       allocatable    :: connectivities(:)                      ! Connectivity matrix
        integer(ip),       allocatable    :: offset(:)                              ! VTK element offset
        integer(I1P),      allocatable    :: cell_types(:)                          ! VTK element type
        integer(ip),       allocatable    :: subelements_connectivity(:,:)          ! Connectivities of subelements (If not linear_order)
        integer(ip)                       :: number_of_nodes    = 0                 ! Number of nodes
        integer(ip)                       :: number_of_elements = 0                 ! Number of elements
        integer(ip)                       :: dimensions         = 0                 ! Dimensions of the mesh
        logical                           :: linear_order       = .false.           ! Order 1 (.true.) or higher
        logical                           :: filled             = .false.           ! Mesh data was already filled
    contains
    private
        procedure, non_overridable, public :: set_fe_space                      => vtk_mesh_and_field_generator_set_fe_space
        procedure, non_overridable, public :: set_linear_order                  => vtk_mesh_and_field_generator_set_linear_order
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
        procedure, non_overridable         :: allocate_subelements_connectivity => vtk_mesh_and_field_generator_allocate_subelements_connectivity
        procedure, non_overridable         :: generate_linear_mesh              => vtk_mesh_and_field_generator_generate_linear_mesh
        procedure, non_overridable         :: generate_superlinear_mesh         => vtk_mesh_and_field_generator_generate_superlinear_mesh
        procedure, non_overridable         :: generate_linear_field             => vtk_mesh_and_field_generator_generate_linear_field
        procedure, non_overridable         :: generate_superlinear_field        => vtk_mesh_and_field_generator_generate_superlinear_field
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


    subroutine vtk_mesh_and_field_generator_set_linear_order(this, linear_order)
    !-----------------------------------------------------------------
    !< Set linear order
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
        logical,           intent(IN)    :: linear_order
    !-----------------------------------------------------------------
        assert(.not. this%filled)
        this%linear_order = linear_order
    end subroutine vtk_mesh_and_field_generator_set_linear_order


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
        number_fields =  this%fe_space%get_number_fe_spaces()
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


    subroutine vtk_mesh_and_field_generator_allocate_subelements_connectivity(this, number_vertices, number_subelements)
    !-----------------------------------------------------------------
    !< Allocate subelements connectivity array
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: number_vertices
        integer(ip),           intent(IN)    :: number_subelements
    !-----------------------------------------------------------------
        assert(.not. allocated(this%subelements_connectivity))
        call memalloc(number_vertices, number_subelements, this%subelements_connectivity, __FILE__, __LINE__)
    end subroutine vtk_mesh_and_field_generator_allocate_subelements_connectivity


    subroutine vtk_mesh_and_field_generator_generate_mesh(this)
    !-----------------------------------------------------------------
    !< Generate the mesh data from fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(INOUT) :: this
    !-----------------------------------------------------------------
        if(.not. this%filled) then
            if(this%linear_order) then
                call this%generate_linear_mesh()
            else
                call this%generate_superlinear_mesh()
            endif
        endif
    end subroutine vtk_mesh_and_field_generator_generate_mesh


    subroutine vtk_mesh_and_field_generator_generate_linear_mesh(this)
    !-----------------------------------------------------------------
    !< Store a linear_order mesh from a triangulation
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),                intent(INOUT) :: this
        type(triangulation_t), pointer                  :: triangulation
        integer(ip)                                     :: i
        integer(ip)                                     :: j
        integer(ip)                                     :: number_nodes
        integer(ip)                                     :: counter
    !-----------------------------------------------------------------
        triangulation => this%fe_space%get_triangulation()
        assert(associated(triangulation))

        this%dimensions = triangulation%num_dims
        this%number_of_elements = triangulation%num_elems
        call this%allocate_elemental_arrays()

        ! Fill VTK cell type and offset arrays and and count nodes
        number_nodes = 0
        do i=1, this%number_of_elements
            if(triangulation%elems(i)%reference_fe_geo%get_topology() == topology_quad) then
                this%cell_types(i) = vtk_pixel
            else
                write(error_unit,*) 'fill_mesh_from_triangulation: Topology not supported'
                check(.false.)
            endif
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
            this%offset(i) = number_nodes
        enddo
        this%number_of_nodes = number_nodes
        call this%allocate_nodal_arrays()
        call this%initialize_coordinates()

        ! Fill VTK coordinate arrays
        counter = 1
        do i=1, this%number_of_elements
            do j=1, triangulation%elems(i)%reference_fe_geo%get_number_vertices()
                this%connectivities(counter) = counter-1
                if (triangulation%num_dims >=1) this%X(counter) = triangulation%elems(i)%coordinates(1,j)
                if (triangulation%num_dims >=2) this%Y(counter) = triangulation%elems(i)%coordinates(2,j)
                if (triangulation%num_dims >=3) this%Z(counter) = triangulation%elems(i)%coordinates(3,j)
                counter = counter + 1
            enddo
            number_nodes = number_nodes + triangulation%elems(i)%reference_fe_geo%get_number_vertices()
        enddo
        this%filled  = .true.
    end subroutine vtk_mesh_and_field_generator_generate_linear_mesh


    subroutine vtk_mesh_and_field_generator_generate_superlinear_mesh(this)
    !-----------------------------------------------------------------
    !< Store a superlinear_order mesh in a from a fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),                intent(INOUT) :: this
        type(point_t),            pointer               :: vertex_coordinates(:)
        type(point_t),            pointer               :: nodal_coordinates(:)
        type(quadrature_t),       pointer               :: nodal_quadrature
        class(reference_fe_t),    pointer               :: reference_fe_geo
        type(finite_element_t),   pointer               :: fe
        class(reference_fe_t),    pointer               :: max_order_reference_fe_phy
        type(fe_map_t)                                  :: fe_map
        integer(ip)                                     :: num_elements
        integer(ip)                                     :: num_nodes_per_element
        integer(ip)                                     :: num_vertices_per_element
        integer(ip)                                     :: num_subelements_per_element
        integer(ip)                                     :: elements_counter
        integer(ip)                                     :: vertex
        integer(ip)                                     :: max_order
        integer(ip)                                     :: dimensions
        integer(ip)                                     :: subelement_vertex
        integer(ip)                                     :: fe_space_index
        integer(ip)                                     :: element_index
        integer(ip)                                     :: subelement_index
        integer(ip)                                     :: max_order_fe_space_index
        integer(ip)                                     :: nodes_counter
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))

        nullify(vertex_coordinates)
        nullify(nodal_coordinates)
        nullify(nodal_quadrature)
        nullify(reference_fe_geo)
        nullify(fe)
        nullify(max_order_reference_fe_phy)

        ! Get max order and its fe_space_index
        max_order                = this%fe_space%get_max_order()
        max_order_fe_space_index = this%fe_space%get_max_order_fe_space_component()
        
        ! Getters
        max_order_reference_fe_phy => this%fe_space%get_reference_fe_phy(max_order_fe_space_index)
        reference_fe_geo           => this%fe_space%get_reference_fe_geo()
        nodal_quadrature           => max_order_reference_fe_phy%get_nodal_quadrature()

        ! Get dimensions from the geometric reference finite element
        dimensions = reference_fe_geo%get_number_dimensions()
        this%dimensions = dimensions
     
        ! Calculate the number of subelems and points for the postprocess
        num_elements                = this%fe_space%get_number_elements()
        num_vertices_per_element    = reference_fe_geo%get_number_vertices()
        num_subelements_per_element = max_order_reference_fe_phy%get_number_subelements()
        this%number_of_elements = num_elements * num_subelements_per_element
        this%number_of_nodes = num_vertices_per_element * this%number_of_elements

        ! Allocate VTK geometry and connectivity data
        call this%allocate_elemental_arrays()
        call this%allocate_nodal_arrays()
        call this%initialize_coordinates()
        
        ! Get the connectivity of the subelements
        call this%allocate_subelements_connectivity(num_vertices_per_element, num_subelements_per_element)
        call max_order_reference_fe_phy%get_subelements_connectivity(this%subelements_connectivity)

        ! Create FE map
        call fe_map%create(nodal_quadrature, reference_fe_geo)

        ! Translate coordinates and connectivities to VTK format
        nodes_counter = 0
        elements_counter = 0
        do element_index = 1, num_elements
            ! Get Finite element
            fe => this%fe_space%get_finite_element(element_index)  

            ! Interpolate coordinates
            vertex_coordinates => fe_map%get_coordinates()
            call fe%get_cell_coordinates(vertex_coordinates)
            call fe_map%compute_quadrature_coordinates()
            nodal_coordinates => fe_map%get_quadrature_coordinates()

            ! Fill VTK mesh
            do subelement_index = 1, num_subelements_per_element
                elements_counter = elements_counter + 1
                do vertex = 1, num_vertices_per_element
                    subelement_vertex = this%subelements_connectivity(vertex, subelement_index)
                    nodes_counter = nodes_counter + 1
                    if(dimensions>=1) this%X(nodes_counter) = nodal_coordinates(subelement_vertex)%get(1)
                    if(dimensions>=2) this%Y(nodes_counter) = nodal_coordinates(subelement_vertex)%get(2)
                    if(dimensions>=3) this%Z(nodes_counter) = nodal_coordinates(subelement_vertex)%get(3)
                    this%connectivities(nodes_counter) = nodes_counter-1
                end do

                ! Store the type of element
                this%cell_types(elements_counter) = vtk_pixel

                ! Fill offset
                this%offset(elements_counter) = nodes_counter
            end do
        end do
        this%filled  = .true.
        call fe_map%free()
    end subroutine vtk_mesh_and_field_generator_generate_superlinear_mesh


    function vtk_mesh_and_field_generator_generate_field(this, fe_function, fe_space_index, field_name, field, number_components) result(E_IO)
    !-----------------------------------------------------------------
    !< Generate the field data from fe_space and fe_function
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t),          intent(INOUT) :: this             !< vtk_mesh_and_field_generator_t derived type
        type(fe_function_t),        intent(IN)    :: fe_function      !< Postprocess field structure to be written
        integer(ip),                intent(IN)    :: fe_space_index   !< Fe space index
        character(len=*),           intent(IN)    :: field_name       !< name of the field
        real(rp), allocatable,      intent(INOUT) :: field(:,:)       !< FIELD(ncomp,nnod)
        integer(ip),                intent(OUT)   :: number_components!< number of components
        integer(ip)                               :: E_IO             !< IO Error
    !-----------------------------------------------------------------
        assert(this%filled)
        if(this%linear_order) then
            E_IO = this%generate_linear_field(fe_function, fe_space_index, field_name, field, number_components)
        else
           E_IO = this%generate_superlinear_field(fe_function, fe_space_index, field_name, field, number_components)
        endif
    end function vtk_mesh_and_field_generator_generate_field


    function vtk_mesh_and_field_generator_generate_linear_field(this, fe_function, fe_space_index, field_name, field, number_components) result(E_IO)
    !-----------------------------------------------------------------
    !< Generate the linear field data from fe_space and fe_function
    !-----------------------------------------------------------------
        implicit none
        class(vtk_mesh_and_field_generator_t),     intent(INOUT) :: this                               !< vtk_mesh_and_field_generator_t derived type
        type(fe_function_t),        intent(IN)    :: fe_function                        !< Postprocess field structure to be written
        integer(ip),                intent(IN)    :: fe_space_index                     !< Fe space index
        character(len=*),           intent(IN)    :: field_name                         !< name of the field
        real(rp), allocatable,      intent(INOUT) :: field(:,:)                         !< FIELD(ncomp,nnod)
        integer(ip),                intent(OUT)   :: number_components                  !< number of components
        type(serial_scalar_array_t), pointer      :: strong_dirichlet_values            !< Strong dirichlet values
        type(finite_element_t),      pointer      :: fe                                 !< finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_origin            !< reference finite element
        class(vector_t),             pointer      :: fe_function_dof_values             !< dof values of the fe_function
        type(i1p_t),                 pointer      :: elem2dof(:)                        !< element 2 dof translator
        integer(ip),                 pointer      :: field_blocks(:)                    !< field blocks
        real(rp),                    pointer      :: strong_dirichlet_values_entries(:) !< strong dirichlet values
        type(list_t),                pointer      :: nodes_vef                          !< list of reference_fe_phy nodes
        real(rp), allocatable                     :: nodal_values(:)                    !< nodal values
        type(list_iterator_t)                     :: nodes_vertex_iterator              !< iterator on vertex nodes of the reference_fe_phy
        integer(ip)                               :: number_elements                    !< number of elements
        integer(ip)                               :: number_vertices                    !< number of geo vertex
        integer(ip)                               :: number_nodes                       !< number of nodes per fe space
        integer(ip)                               :: element_index                      !< element index
        integer(ip)                               :: component_index                    !< component index
        integer(ip)                               :: vertex_index                       !< vertex index
        integer(ip)                               :: node_index                         !< node index
        integer(ip)                               :: E_IO                               !< IO Error
    !-----------------------------------------------------------------
        assert(associated(this%fe_space))
        assert(this%filled)
        assert(this%linear_order)

        nullify(strong_dirichlet_values)
        nullify(fe)
        nullify(reference_fe_phy_origin)
        nullify(fe_function_dof_values)
        nullify(elem2dof)
        nullify(field_blocks)
        nullify(strong_dirichlet_values_entries)

        ! Point to some fe_space content
        reference_fe_phy_origin => this%fe_space%get_reference_fe_phy(fe_space_index)
        field_blocks            => this%fe_space%get_field_blocks()
        nodes_vef               => reference_fe_phy_origin%get_nodes_vef()


        ! Extract nodal values associated to dirichlet bcs and dof values
        strong_dirichlet_values         => fe_function%get_strong_dirichlet_values()
        strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
        fe_function_dof_values          => fe_function%get_dof_values()

        ! Get number components, elements, nodes and vertices
        number_elements   = this%fe_space%get_number_elements()
        number_components = reference_fe_phy_origin%get_number_field_components()
        number_nodes      = reference_fe_phy_origin%get_number_nodes()
        number_vertices   = reference_fe_phy_origin%get_number_vertices()

        ! Allocate nodal values per finite element and VTK field
        if(allocated(field)) call memfree(field, __FILE__, __LINE__)
        call memalloc(number_nodes, nodal_values, __FILE__, __LINE__)
        call memalloc(number_components, number_elements*number_vertices , field, __FILE__, __LINE__)

        do element_index=1, number_elements
            fe => this%fe_space%get_finite_element(element_index)  
            elem2dof => fe%get_elem2dof()

            ! Extract nodal values associated to dofs
            call fe_function_dof_values%extract_subvector ( field_blocks(fe_space_index), number_nodes, elem2dof(fe_space_index)%p, nodal_values )

            ! Build field in VTK-like format
            ! Loop on geometrical nodes in subelement
            do vertex_index=1, number_vertices
                nodes_vertex_iterator = nodes_vef%create_iterator(vertex_index)
                assert(nodes_vertex_iterator%get_size() == number_components)
                ! Loop on field components
                do while(.not. nodes_vertex_iterator%is_upper_bound())
                    node_index = nodes_vertex_iterator%get_current()
                    if ( elem2dof(fe_space_index)%p(node_index) < 0 ) then
                        ! Fill field with strong dirichlet values
                        field(reference_fe_phy_origin%get_component_node(node_index), (element_index-1)*number_vertices+vertex_index) = &
                            strong_dirichlet_values_entries(-elem2dof(fe_space_index)%p(node_index))
                    else
                        ! Fill field with nodal values
                        field(reference_fe_phy_origin%get_component_node(node_index), (element_index-1)*number_vertices+vertex_index) = nodal_values(node_index)
                    endif
                    call nodes_vertex_iterator%next()
                end do
            end do
        enddo
        call memfree(nodal_values, __FILE__, __LINE__)
    end function vtk_mesh_and_field_generator_generate_linear_field


    function vtk_mesh_and_field_generator_generate_superlinear_field(this, fe_function, fe_space_index, field_name, field, number_components) result(E_IO)
    !-----------------------------------------------------------------
    !< Write superlinear field to file
    !-----------------------------------------------------------------
        implicit none
        class(vtk_mesh_and_field_generator_t),          intent(INOUT) :: this                               !< this raw mesh
        type(fe_function_t),        intent(IN)    :: fe_function                        !< Postprocess field structure to be written
        integer(ip),                intent(IN)    :: fe_space_index                     !< Fe space index
        character(len=*),           intent(IN)    :: field_name                         !< name of the field
        real(rp), allocatable,      intent(INOUT) :: field(:,:)                         !< FIELD(ncomp,nnod)
        integer(ip),                intent(OUT)   :: number_components                  !< number of components
        type(quadrature_t),          pointer      :: nodal_quadrature_target            !< Nodal quadrature
        type(serial_scalar_array_t), pointer      :: strong_dirichlet_values            !< Strong dirichlet values
        type(finite_element_t),      pointer      :: fe                                 !< finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_origin            !< reference finite element
        class(reference_fe_t),       pointer      :: reference_fe_phy_target            !< reference finite element
        class(vector_t),             pointer      :: fe_function_dof_values             !< dof values of the fe_function
        type(i1p_t),                 pointer      :: elem2dof(:)                        !< element 2 dof translator
        integer(ip),                 pointer      :: field_blocks(:)                    !< field blocks
        real(rp),                    pointer      :: strong_dirichlet_values_entries(:) !< strong dirichlet values
        real(rp), allocatable                     :: nodal_values_origin(:)             !< nodal values of the origin fe_space
        real(rp), allocatable                     :: nodal_values_target(:)             !< nodal values for the interpolation
        type(interpolation_t)                     :: interpolation                      !< interpolator
        integer(ip)                               :: number_elements                    !< number of elements
        integer(ip)                               :: number_subelements                 !< number of subelements per element
        integer(ip)                               :: number_vertices                    !< number of geo vertex
        integer(ip)                               :: number_nodes_scalar                !< number of scalar nodes
        integer(ip)                               :: number_nodes                       !< number of nodes per fe space
        integer(ip)                               :: element_index                      !< element index
        integer(ip)                               :: component_index                    !< component index
        integer(ip)                               :: node_index                         !< node index
        integer(ip)                               :: subelement_index                   !< subelement index
        integer(ip)                               :: subnode_index                      !< subelement node index
        integer(ip)                               :: order                              !< current fe_space order
        integer(ip)                               :: max_order                          !< max fe_space order
        integer(ip)                               :: max_order_fe_space_index           !< index of the fe_space with max order
        integer(ip)                               :: idx                                !< Indices
        integer(ip)                               :: E_IO                               !< IO Error
        logical                                   :: ft                                 !< Fine Task
    !-----------------------------------------------------------------
        E_IO = 0

        assert(associated(this%fe_space))
        assert(this%filled)
        assert(.not. this%linear_order)

        nullify(nodal_quadrature_target)
        nullify(strong_dirichlet_values)
        nullify(fe)
        nullify(reference_fe_phy_origin)
        nullify(reference_fe_phy_target)
        nullify(fe_function_dof_values)
        nullify(elem2dof)
        nullify(field_blocks)
        nullify(strong_dirichlet_values_entries)

        ! Get max order and its fe_space_index
        max_order = this%fe_space%get_max_order()
        max_order_fe_space_index = this%fe_space%get_max_order_fe_space_component()

        ! Point physical referece fe, field blocks, nodal quadrature and subelements connectivity
        reference_fe_phy_origin  => this%fe_space%get_reference_fe_phy(fe_space_index)
        reference_fe_phy_target  => this%fe_space%get_reference_fe_phy(max_order_fe_space_index)
        field_blocks             => this%fe_space%get_field_blocks()
        nodal_quadrature_target  => reference_fe_phy_target%get_nodal_quadrature()

        ! Extract nodal values associated to dirichlet bcs and dof values
        strong_dirichlet_values         => fe_function%get_strong_dirichlet_values()
        strong_dirichlet_values_entries => strong_dirichlet_values%get_entries()
        fe_function_dof_values          => fe_function%get_dof_values()
        
        ! Calculate number nodes of target field given a certain order of interpolation, the fe_type and the number of components
        number_elements     = this%fe_space%get_number_elements()
        number_components   = reference_fe_phy_origin%get_number_field_components()
        number_subelements  = reference_fe_phy_target%get_number_subelements()
        number_nodes_scalar = reference_fe_phy_origin%get_number_nodes_scalar()
        number_nodes        = reference_fe_phy_origin%get_number_nodes()
        number_vertices     = reference_fe_phy_target%get_number_vertices()
        order               = reference_fe_phy_origin%get_order()

        ! Allocate VTK field, origin and target nodal values
        if(allocated(field)) call memfree(field, __FILE__, __LINE__)
        call memalloc(number_nodes, nodal_values_origin, __FILE__, __LINE__)
        call memalloc(number_components*nodal_quadrature_target%get_number_quadrature_points(), nodal_values_target, __FILE__, __LINE__)
        call memalloc(number_components, number_subelements*number_elements*number_vertices , field, __FILE__, __LINE__)

        ! Create interpolation from nodal quadrature
        call reference_fe_phy_origin%create_interpolation(nodal_quadrature_target, interpolation)
        
        ! Loop on elements
        subnode_index=1
        do element_index=1, number_elements
            fe => this%fe_space%get_finite_element(element_index)  
            elem2dof => fe%get_elem2dof()

            ! Extract nodal values associated to dofs
            call fe_function_dof_values%extract_subvector ( field_blocks(fe_space_index), number_nodes, elem2dof(fe_space_index)%p, nodal_values_origin )

            ! Fill nodal values with strong dirichlet values
            do node_index = 1, number_nodes
                if ( elem2dof(fe_space_index)%p(node_index) < 0 ) then
                    nodal_values_origin(node_index) = strong_dirichlet_values_entries(-elem2dof(fe_space_index)%p(node_index))
                end if
            end do

            ! interpolate nodal values if needed
            if(order /= max_order) then
                call reference_fe_phy_origin%interpolate_nodal_values( interpolation, nodal_values_origin, nodal_values_target ) 
            else
                nodal_values_target=nodal_values_origin
            endif

            ! Build field in VTK-like format
            ! Loop on subelements
            do subelement_index = 1, number_subelements
                ! Loop on geometrical nodes per subelement
                do node_index=1, number_vertices
                    ! Loop on components
                    do component_index=1, number_components
                        idx = this%subelements_connectivity(node_index,subelement_index) + (component_index-1)*number_nodes_scalar
                        field(reference_fe_phy_target%get_component_node(idx) , subnode_index) = nodal_values_target(idx)
                    enddo
                    subnode_index=subnode_index+1
                end do
            end do
        enddo
        call interpolation%free()
        call memfree(nodal_values_origin, __FILE__, __LINE__)
        call memfree(nodal_values_target, __FILE__, __LINE__)
      end function vtk_mesh_and_field_generator_generate_superlinear_field


    subroutine vtk_mesh_and_field_generator_free(this) 
    !-----------------------------------------------------------------
    !< Free the vtk_mesh_and_field_generator_t derived type
    !-----------------------------------------------------------------
        class(vtk_mesh_and_field_generator_t), intent(inout) :: this
        integer(ip)                      :: i
    !-----------------------------------------------------------------
        if(allocated(this%X))                        call memfree(this%X, __FILE__, __LINE__)
        if(allocated(this%Y))                        call memfree(this%Y, __FILE__, __LINE__)
        if(allocated(this%Z))                        call memfree(this%Z, __FILE__, __LINE__)
        if(allocated(this%connectivities))           call memfree(this%connectivities, __FILE__, __LINE__)
        if(allocated(this%offset))                   call memfree(this%offset, __FILE__, __LINE__)
        if(allocated(this%cell_types))               call memfree(this%cell_types, __FILE__, __LINE__)
        if(allocated(this%subelements_connectivity)) call memfree(this%subelements_connectivity, __FILE__, __LINE__)
        nullify(this%fe_space)
        this%number_of_nodes    = 0
        this%number_of_elements = 0
        this%linear_order       = .false.
        this%filled             = .false.
    end subroutine

end module vtk_mesh_and_field_generator
