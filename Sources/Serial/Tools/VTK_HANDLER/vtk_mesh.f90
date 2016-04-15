
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

module vtk_mesh

USE triangulation_names,         only: triangulation_t
USE serial_fe_space_names,       only: serial_fe_space_t
USE types_names
USE memor_names
USE Lib_VTK_IO

implicit none

#include "debug.i90"

private


    ! STATE PARAMETERS
    integer(ip), parameter :: VTK_STATE_UNKNOWN          = 0
    integer(ip), parameter :: VTK_STATE_WRITE_STARTED    = 1
    integer(ip), parameter :: VTK_STATE_POINTDATA_OPENED = 2
    integer(ip), parameter :: VTK_STATE_POINTDATA_CLOSED = 3
    integer(ip), parameter :: VTK_STATE_ENDED            = 4
    integer(ip), parameter :: VTK_STATE_READ_STARTED     = 5

    ! Type for storing mesh data
    type vtk_mesh_t
    private
        class(serial_fe_space_t), pointer :: fe_space      => NULL()                ! Poins to fe_space_t
        character(len=:),  allocatable    :: directory_path                         ! Directory where the results are going to be stored
        character(len=:),  allocatable    :: name_prefix                            ! Name prefix of the VTK files
        integer(ip)                       :: file_id = 0                            ! File descriptor (unit)
        real(rp),          allocatable    :: X(:)                                   ! Mesh X coordintates
        real(rp),          allocatable    :: Y(:)                                   ! Mesh Y coordintates
        real(rp),          allocatable    :: Z(:)                                   ! Mesh Z coordintates
        integer(ip),       allocatable    :: connectivities(:)                      ! Connectivity matrix
        integer(ip),       allocatable    :: offset(:)                              ! VTK element offset
        integer(ip),       allocatable    :: cell_types(:)                          ! VTK element type
        integer(ip),       allocatable    :: subelements_connectivity(:,:)          ! Connectivities of subelements (If not linear_order)
        integer(ip)                       :: number_of_nodes    = 0                 ! Number of nodes
        integer(ip)                       :: number_of_elements = 0                 ! Number of elements
        integer(ip)                       :: dimensions         = 0                 ! Dimensions of the mesh
        logical                           :: linear_order       = .false.           ! Order 1 (.true.) or higher
        logical                           :: filled             = .false.           ! Mesh data was already filled
        integer(ip)                       :: status             = VTK_STATE_UNKNOWN ! Status of the write process
    contains
    private
        procedure, non_overridable, public :: move_to                           => vtk_mesh_move_to
        procedure, non_overridable, public :: set_fe_space                      => vtk_mesh_set_fe_space
        procedure, non_overridable, public :: get_fe_space                      => vtk_mesh_get_fe_space
        procedure, non_overridable, public :: get_triangulation                 => vtk_mesh_get_triangulation
        procedure, non_overridable, public :: set_path                          => vtk_mesh_set_path
        procedure, non_overridable, public :: get_path                          => vtk_mesh_get_path
        procedure, non_overridable, public :: set_prefix                        => vtk_mesh_set_prefix
        procedure, non_overridable, public :: get_prefix                        => vtk_mesh_get_prefix
        procedure, non_overridable, public :: set_linear_order                  => vtk_mesh_set_linear_order
        procedure, non_overridable, public :: is_linear_order                   => vtk_mesh_is_linear_order
        procedure, non_overridable, public :: set_filled                        => vtk_mesh_set_filled
        procedure, non_overridable, public :: is_filled                         => vtk_mesh_is_filled
        procedure, non_overridable, public :: set_number_nodes                  => vtk_mesh_set_number_nodes
        procedure, non_overridable, public :: get_number_nodes                  => vtk_mesh_get_number_nodes
        procedure, non_overridable, public :: set_number_elements               => vtk_mesh_set_number_elements
        procedure, non_overridable, public :: get_number_elements               => vtk_mesh_get_number_elements
        procedure, non_overridable, public :: set_dimensions                    => vtk_mesh_set_dimensions
        procedure, non_overridable, public :: get_dimensions                    => vtk_mesh_get_dimensions
        procedure, non_overridable, public :: set_cell_type                     => vtk_mesh_set_cell_type
        procedure, non_overridable, public :: get_cell_type                     => vtk_mesh_get_cell_type
        procedure, non_overridable, public :: set_offset                        => vtk_mesh_set_offset
        procedure, non_overridable, public :: get_offset                        => vtk_mesh_get_offset
        procedure, non_overridable, public :: set_connectivity                  => vtk_mesh_set_connectivity
        procedure, non_overridable, public :: get_connectivity                  => vtk_mesh_get_connectivity
        procedure, non_overridable, public :: set_x_coordinate                  => vtk_mesh_set_x_coordinate
        procedure, non_overridable, public :: get_x_coordinate                  => vtk_mesh_get_x_coordinate
        procedure, non_overridable, public :: set_y_coordinate                  => vtk_mesh_set_y_coordinate
        procedure, non_overridable, public :: get_y_coordinate                  => vtk_mesh_get_y_coordinate
        procedure, non_overridable, public :: set_z_coordinate                  => vtk_mesh_set_z_coordinate
        procedure, non_overridable, public :: get_z_coordinate                  => vtk_mesh_get_z_coordinate
        procedure, non_overridable, public :: initialize_coordinates            => vtk_mesh_initialize_coordinates
        procedure, non_overridable, public :: get_subelements_connectivity      => vtk_mesh_get_subelements_connectivity
        procedure, non_overridable, public :: allocate_nodal_arrays             => vtk_mesh_allocate_nodal_arrays
        procedure, non_overridable, public :: allocate_elemental_arrays         => vtk_mesh_allocate_elemental_arrays
        procedure, non_overridable, public :: allocate_subelements_connectivity => vtk_mesh_allocate_subelements_connectivity
        procedure, non_overridable, public :: begin_write                       => vtk_mesh_begin_write
        procedure, non_overridable, public :: write_node_field                  => vtk_mesh_write_node_field
        procedure, non_overridable, public :: end_write                         => vtk_mesh_end_write
        procedure, non_overridable, public :: free                              => vtk_mesh_free
    end type vtk_mesh_t

public :: vtk_mesh_t

contains

    subroutine vtk_mesh_move_to(this, mesh)
    !-----------------------------------------------------------------
    !< Move this to other mesh host
    !-----------------------------------------------------------------
        class(vtk_mesh_t), intent(INOUT) :: this                      ! Input mesh
        type(vtk_mesh_t),  intent(INOUT) :: mesh                      ! Output mesh
    !-----------------------------------------------------------------
        call mesh%free()
        if(associated(this%fe_space)) mesh%fe_space => this%fe_space
        if(allocated(this%directory_path)) mesh%directory_path = this%directory_path
        if(allocated(this%name_prefix)) mesh%name_prefix = this%name_prefix
        if(allocated(this%X)) call memmovealloc(this%X, mesh%X, __FILE__, __LINE__)
        if(allocated(this%Y)) call memmovealloc(this%Y, mesh%Y, __FILE__, __LINE__)
        if(allocated(this%Z)) call memmovealloc(this%Z, mesh%Z, __FILE__, __LINE__)
        if(allocated(this%offset)) call memmovealloc(this%offset, mesh%offset, __FILE__, __LINE__)
        if(allocated(this%cell_types)) call memmovealloc(this%cell_types, mesh%cell_types, __FILE__, __LINE__)
        if(allocated(this%subelements_connectivity)) call memmovealloc(this%subelements_connectivity, mesh%subelements_connectivity, __FILE__, __LINE__)
        mesh%file_id = this%file_id
        mesh%number_of_nodes = this%number_of_nodes
        mesh%number_of_elements = this%number_of_elements
        mesh%dimensions = this%dimensions
        mesh%linear_order = this%linear_order
        mesh%filled = this%filled
        mesh%status = this%status
        call this%free()
    end subroutine vtk_mesh_move_to


    subroutine vtk_mesh_set_fe_space(this, fe_space)
    !-----------------------------------------------------------------
    !< Set the fe_space
    !-----------------------------------------------------------------
        class(vtk_mesh_t),                intent(INOUT) :: this
        class(serial_fe_space_t), target, intent(IN)    :: fe_space
    !-----------------------------------------------------------------
        this%fe_space => fe_space
    end subroutine vtk_mesh_set_fe_space


    function vtk_mesh_get_fe_space(this) result(fe_space)
    !-----------------------------------------------------------------
    !< Get the fe_space pointer
    !-----------------------------------------------------------------
        class(vtk_mesh_t),                intent(IN) :: this
        class(serial_fe_space_t), pointer            :: fe_space
    !-----------------------------------------------------------------
        fe_space => this%fe_space
    end function vtk_mesh_get_fe_space


    function vtk_mesh_get_triangulation(this) result(triangulation)
    !-----------------------------------------------------------------
    !< Get the fe_space pointer
    !-----------------------------------------------------------------
        class(vtk_mesh_t),             intent(IN) :: this
        type(triangulation_t), pointer            :: triangulation
    !-----------------------------------------------------------------
        triangulation => this%fe_space%get_triangulation()
    end function vtk_mesh_get_triangulation


    subroutine vtk_mesh_set_path(this, path)
    !-----------------------------------------------------------------
    !< Set the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        character(len=*),      intent(IN)    :: path
    !-----------------------------------------------------------------
        this%directory_path = path
    end subroutine vtk_mesh_set_path


    subroutine vtk_mesh_get_path(this, path)
    !-----------------------------------------------------------------
    !< Get the name of the output directory
    !-----------------------------------------------------------------
        class(vtk_mesh_t),             intent(IN)    :: this
        character(len=:), allocatable, intent(INOUT) :: path
    !-----------------------------------------------------------------
        assert(allocated(this%directory_path))
        path = this%directory_path
    end subroutine vtk_mesh_get_path


    subroutine vtk_mesh_set_prefix(this, prefix)
    !-----------------------------------------------------------------
    !< Set the file name prefix
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        character(len=*),      intent(IN)    :: prefix
    !-----------------------------------------------------------------
        this%name_prefix = prefix
    end subroutine vtk_mesh_set_prefix


    subroutine vtk_mesh_get_prefix(this, prefix)
    !-----------------------------------------------------------------
    !< Get the filename prefix
    !-----------------------------------------------------------------
        class(vtk_mesh_t),             intent(IN)    :: this
        character(len=:), allocatable, intent(INOUT) :: prefix
    !-----------------------------------------------------------------
        assert(allocated(this%name_prefix))
        prefix = this%name_prefix
    end subroutine vtk_mesh_get_prefix


    subroutine vtk_mesh_set_filled(this, filled)
    !-----------------------------------------------------------------
    !< True if the mesh data is already filled
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        logical,               intent(IN)    :: filled
    !-----------------------------------------------------------------
        this%filled = filled
    end subroutine vtk_mesh_set_filled


    function vtk_mesh_is_filled(this) result(filled)
    !-----------------------------------------------------------------
    !< Ask if the mesh data is filled
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        logical                              :: filled
    !-----------------------------------------------------------------
        filled = this%filled
    end function vtk_mesh_is_filled


    subroutine vtk_mesh_set_linear_order(this, linear_order)
    !-----------------------------------------------------------------
    !< Set linear order
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        logical,               intent(IN)    :: linear_order
    !-----------------------------------------------------------------
        this%linear_order = linear_order
    end subroutine vtk_mesh_set_linear_order


    function vtk_mesh_is_linear_order(this) result(linear_order)
    !-----------------------------------------------------------------
    !< Ask if the stored mesh data is from linear order mesh
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        logical                              :: linear_order
    !-----------------------------------------------------------------
        linear_order = this%linear_order
    end function vtk_mesh_is_linear_order


    subroutine vtk_mesh_set_number_nodes(this, number_nodes)
    !-----------------------------------------------------------------
    !< Set the number of nodes
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: number_nodes
    !-----------------------------------------------------------------
        this%number_of_nodes = number_nodes
    end subroutine vtk_mesh_set_number_nodes


    function vtk_mesh_get_number_nodes(this) result(number_nodes)
    !-----------------------------------------------------------------
    !< Return the number of nodes
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: number_nodes
    !-----------------------------------------------------------------
        number_nodes = this%number_of_nodes
    end function vtk_mesh_get_number_nodes


    subroutine vtk_mesh_set_number_elements(this, number_elements)
    !-----------------------------------------------------------------
    !< Set the number of elements
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: number_elements
    !-----------------------------------------------------------------
        this%number_of_elements = number_elements
    end subroutine vtk_mesh_set_number_elements


    function vtk_mesh_get_number_elements(this) result(number_elements)
    !-----------------------------------------------------------------
    !< Return the number of elements
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: number_elements
    !-----------------------------------------------------------------
        number_elements = this%number_of_elements
    end function vtk_mesh_get_number_elements


    subroutine vtk_mesh_set_dimensions(this, dimensions)
    !-----------------------------------------------------------------
    !< Set the space dimensions of the mesh
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: dimensions
    !-----------------------------------------------------------------
        this%dimensions = dimensions
    end subroutine vtk_mesh_set_dimensions


    function vtk_mesh_get_dimensions(this) result(dimensions)
    !-----------------------------------------------------------------
    !< Return the space dimensions of the mesh
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: dimensions
    !-----------------------------------------------------------------
        dimensions = this%dimensions
    end function vtk_mesh_get_dimensions


    subroutine vtk_mesh_set_cell_type(this, index, type)
    !-----------------------------------------------------------------
    !< Set the cell types given the index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        integer(ip),           intent(IN)    :: type
    !-----------------------------------------------------------------
        assert(allocated(this%cell_types))
        assert(index>0 .and. index<=this%number_of_elements)
        this%cell_types(index) = type
    end subroutine vtk_mesh_set_cell_type


    function vtk_mesh_get_cell_type(this, index) result(type)
    !-----------------------------------------------------------------
    !< Return cell type given the cell index of 
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        integer(ip)                          :: type
    !-----------------------------------------------------------------
        assert(allocated(this%cell_types))
        assert(index>0 .and. index<=this%number_of_elements)
        type = this%cell_types(index)
    end function vtk_mesh_get_cell_type


    subroutine vtk_mesh_set_offset(this, index, offset)
    !-----------------------------------------------------------------
    !< Set the offset given the cellindex
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        integer(ip),           intent(IN)    :: offset
    !-----------------------------------------------------------------
        assert(allocated(this%offset))
        assert(index>0 .and. index<=this%number_of_elements)
        this%offset(index) = offset
    end subroutine vtk_mesh_set_offset


    function vtk_mesh_get_offset(this, index) result(offset)
    !-----------------------------------------------------------------
    !< Return the offset given the cellindex
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        integer(ip)                          :: offset
    !-----------------------------------------------------------------
        assert(allocated(this%offset))
        assert(index>0 .and. index<=this%number_of_elements)
        offset = this%offset(index)
    end function vtk_mesh_get_offset


    subroutine vtk_mesh_set_connectivity(this, index, node)
    !-----------------------------------------------------------------
    !< Set element connectivities
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        integer(ip),           intent(IN)    :: node
    !-----------------------------------------------------------------
        assert(allocated(this%connectivities))
        assert(index>0 .and. index<=this%number_of_nodes)
        this%connectivities(index) = node
    end subroutine vtk_mesh_set_connectivity


    function vtk_mesh_get_connectivity(this, index) result(node)
    !-----------------------------------------------------------------
    !< Return element connectivities given index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        integer(ip)                          :: node
    !-----------------------------------------------------------------
        assert(allocated(this%connectivities))
        assert(index>0 .and. index<=this%number_of_nodes)
        node = this%connectivities(index)
    end function vtk_mesh_get_connectivity


    subroutine vtk_mesh_set_x_coordinate(this, index, coordinate)
    !-----------------------------------------------------------------
    !< Set the x coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        real(rp),              intent(IN)    :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%X))
        assert(index>0 .and. index<=this%number_of_nodes)
        this%X(index) = coordinate
    end subroutine vtk_mesh_set_x_coordinate


    function vtk_mesh_get_x_coordinate(this, index) result(coordinate)
    !-----------------------------------------------------------------
    !< Return x coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        real(rp)                             :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%X))
        assert(index>0 .and. index<=this%number_of_nodes)
        coordinate = this%X(index)
    end function vtk_mesh_get_x_coordinate


    subroutine vtk_mesh_set_y_coordinate(this, index, coordinate)
    !-----------------------------------------------------------------
    !< Set the y coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        real(rp),              intent(IN)    :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%Y))
        assert(index>0 .and. index<=this%number_of_nodes)
        this%Y(index) = coordinate
    end subroutine vtk_mesh_set_y_coordinate


    function vtk_mesh_get_y_coordinate(this, index) result(coordinate)
    !-----------------------------------------------------------------
    !< Return the y coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        real(rp)                             :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%Y))
        assert(index>0 .and. index<=this%number_of_nodes)
        coordinate = this%Y(index)
    end function vtk_mesh_get_y_coordinate


    subroutine vtk_mesh_set_z_coordinate(this, index, coordinate)
    !-----------------------------------------------------------------
    !< Set the z coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: index
        real(rp),              intent(IN)    :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%Z))
        assert(index>0 .and. index<=this%number_of_nodes)
        this%Z(index) = coordinate
    end subroutine vtk_mesh_set_z_coordinate


    function vtk_mesh_get_z_coordinate(this, index) result(coordinate)
    !-----------------------------------------------------------------
    !< Return the Z coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip)                          :: index
        real(rp)                             :: coordinate
    !-----------------------------------------------------------------
        assert(allocated(this%Z))
        assert(index>0 .and. index<=this%number_of_nodes)
        coordinate = this%Z(index)
    end function vtk_mesh_get_z_coordinate


    subroutine vtk_mesh_initialize_coordinates(this)
    !-----------------------------------------------------------------
    !< Set the z coordinate given the node index
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(allocated(this%X))
        assert(allocated(this%Y))
        assert(allocated(this%Z))
        this%X = 0.0_rp
        this%Y = 0.0_rp
        this%Z = 0.0_rp
    end subroutine vtk_mesh_initialize_coordinates


    function vtk_mesh_get_subelements_connectivity(this) result(subelements_connectivity)
    !-----------------------------------------------------------------
    !< Return a pointer to the subelements connectivity array
    !-----------------------------------------------------------------
        class(vtk_mesh_t), target, intent(INOUT) :: this
        integer(ip),       pointer               :: subelements_connectivity(:,:) 
    !-----------------------------------------------------------------
        assert(allocated(this%subelements_connectivity))
        nullify(subelements_connectivity)
        subelements_connectivity => this%subelements_connectivity
    end function vtk_mesh_get_subelements_connectivity


    subroutine vtk_mesh_allocate_elemental_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate all arrays of size number of elements
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(.not. allocated(this%offset))
        assert(.not. allocated(this%cell_types))
        call memalloc(this%number_of_elements, this%offset, __FILE__, __LINE__)
        call memalloc(this%number_of_elements, this%cell_types, __FILE__, __LINE__)
    end subroutine vtk_mesh_allocate_elemental_arrays


    subroutine vtk_mesh_allocate_nodal_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate all arrays with size number of nodes
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
    !-----------------------------------------------------------------
        assert(.not. allocated(this%connectivities))
        assert(.not. allocated(this%X))
        assert(.not. allocated(this%Y))
        assert(.not. allocated(this%Z))
        call memalloc (this%number_of_nodes, this%connectivities, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%X, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%Y, __FILE__,__LINE__)
        call memalloc (this%number_of_nodes, this%Z, __FILE__,__LINE__)
    end subroutine vtk_mesh_allocate_nodal_arrays


    subroutine vtk_mesh_allocate_subelements_connectivity(this, number_vertices, number_subelements)
    !-----------------------------------------------------------------
    !< Allocate subelements connectivity array
    !-----------------------------------------------------------------
        class(vtk_mesh_t),     intent(INOUT) :: this
        integer(ip),           intent(IN)    :: number_vertices
        integer(ip),           intent(IN)    :: number_subelements
    !-----------------------------------------------------------------
        assert(.not. allocated(this%subelements_connectivity))
        call memalloc(number_vertices, number_subelements, this%subelements_connectivity, __FILE__, __LINE__)
    end subroutine vtk_mesh_allocate_subelements_connectivity


    function vtk_mesh_begin_write(this, file_name, part_number, time_step, format) result(E_IO)
    !-----------------------------------------------------------------
    !< Start the writing of a single VTK file to disk (if I am fine MPI task)
    !< Writes connectivities and coordinates ( VTK_INI_XML, 
    !< VTK_GEO_XML, VTK_CON_XML )
    !-----------------------------------------------------------------
        class(vtk_mesh_t),          intent(INOUT) :: this        !< VTK_t derived type
        character(len=*),           intent(IN)    :: file_name   !< VTK File NAME
        integer(ip),                intent(IN)    :: part_number !< Number of the PART
        real(rp),                   intent(IN)    :: time_step   !< Time STEP value
        character(len=*), optional, intent(IN)    :: format      !< Ouput ForMaT
        character(len=:), allocatable             :: of          !< Real Output Format
        integer(ip)                               :: nnods       !< Number of NODeS
        integer(ip)                               :: nels        !< Number of ELementS
        integer(ip)                               :: E_IO        !< Error IO
      ! ----------------------------------------------------------------------------------
        if(this%status == VTK_STATE_UNKNOWN .or. this%status == VTK_STATE_ENDED) then
            of = 'raw'
            if(present(format)) of = trim(adjustl(format))

            E_IO = VTK_INI_XML(output_format = trim(adjustl(of)),   &
                               filename = trim(adjustl(file_name)), &
                               mesh_topology = 'UnstructuredGrid',  &
                               cf=this%file_id)
            E_IO = VTK_GEO_XML(NN = this%number_of_nodes,    &
                               NC = this%number_of_elements, &
                               X  = this%X,                  &
                               Y  = this%Y,                  &
                               Z  = this%Z,                  &
                               cf = this%file_id)
            E_IO = VTK_CON_XML(NC = this%number_of_elements,    &
                               connect   = this%connectivities, &
                               offset    = this%offset,         &
                               cell_type = int(this%cell_types,1),     &
                               cf        = this%file_id)

            this%status = VTK_STATE_WRITE_STARTED
        endif
    end function vtk_mesh_begin_write


    function vtk_mesh_write_node_field(this, field, field_name) result(E_IO)
    !-----------------------------------------------------------------
    !< Writes a VTK field ( VTK_DAT_XML )
    !-----------------------------------------------------------------
        class(vtk_mesh_t),          intent(INOUT) :: this        !< VTK_t derived type
        character(len=*),           intent(IN)    :: field_name  !< VTK field NAME
        real(rp),                   intent(IN)    :: field(:,:)  !< Field to write
        integer(ip)                               :: E_IO        !< Error IO
    !-----------------------------------------------------------------
        if((this%status == VTK_STATE_WRITE_STARTED) .or. (this%status == VTK_STATE_POINTDATA_OPENED)) then
            if(this%status == VTK_STATE_WRITE_STARTED) then
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=this%file_id)
                this%status = VTK_STATE_POINTDATA_OPENED
            endif
            E_IO = VTK_VAR_XML(NC_NN=this%number_of_nodes,N_COL=size(field,1), varname=field_name, var=field, cf=this%file_id)
        endif
    end function vtk_mesh_write_node_field


    function vtk_mesh_end_write(this) result(E_IO)
    !-----------------------------------------------------------------
    !< Ends the writing of a single VTK file to disk (if I am fine MPI task)
    !< Closes geometry ( VTK_END_XML, VTK_GEO_XML )
    !-----------------------------------------------------------------
        implicit none
        class(vtk_mesh_t),          intent(INOUT) :: this        !< VTK_t derived type
        integer(ip)                               :: nm          !< Real Number of Mesh
        integer(ip)                               :: E_IO        !< IO Error
        logical                                   :: ft          !< Fine Task
    ! -----------------------------------------------------------------
        if ((this%status >= VTK_STATE_WRITE_STARTED) .and. (this%status /= VTK_STATE_ENDED)) then
            if(this%status == VTK_STATE_POINTDATA_OPENED) then
                E_IO = VTK_DAT_XML(var_location='node',var_block_action='close', cf=this%file_id)
                this%status = VTK_STATE_POINTDATA_CLOSED
            endif
              
            E_IO = VTK_GEO_XML(cf=this%file_id)
            E_IO = VTK_END_XML(cf=this%file_id)
              
            this%status = VTK_STATE_ENDED
       endif
    end function vtk_mesh_end_write


    subroutine vtk_mesh_free(this) 
    !-----------------------------------------------------------------
    !< Free the vtk_mesh_t derived type
    !-----------------------------------------------------------------
        class(vtk_mesh_t), intent(inout) :: this
        integer(ip)                      :: i
        integer(ip)                      :: error
    !-----------------------------------------------------------------
        error = 0
        if(allocated(this%directory_path))           deallocate(this%directory_path, stat=error)
        assert(error==0)
        if(allocated(this%name_prefix))              deallocate(this%name_prefix, stat=error)
        assert(error==0)
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
        this%status             = VTK_STATE_UNKNOWN
    end subroutine

end module vtk_mesh
