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

module vtk_output_handler_names

USE iso_c_binding
USE types_names
USE memor_names
USE lib_vtk_io
USE environment_names
USE vtk_parameters_names
USE output_handler_base_names
USE output_handler_fe_field_names
USE fe_space_names,             only: serial_fe_space_t
USE output_handler_patch_names, only: patch_subcell_iterator_t
USE reference_fe_names,         only: topology_hex, topology_tet
USE iso_fortran_env,            only: error_unit
USE IR_Precision,               only: I1P


implicit none
#include "debug.i90"
private

    interface
        function mkdir_recursive(path) bind(c,name="mkdir_recursive")
            use iso_c_binding
            integer(kind=c_int) :: mkdir_recursive
            character(kind=c_char,len=1), intent(IN) :: path(*)
        end function mkdir_recursive
    end interface

    type, extends(output_handler_base_t) :: vtk_output_handler_t
    private 
        integer(ip)               :: file_id
        real(rp),     allocatable :: X(:)
        real(rp),     allocatable :: Y(:)
        real(rp),     allocatable :: Z(:)
        integer(ip),  allocatable :: Offset(:)
        integer(I1P), allocatable :: CellTypes(:)
        integer(ip),  allocatable :: Connectivities(:)
        integer(ip)               :: node_offset = 0
        integer(ip)               :: cell_offset = 0
    contains
        procedure, public :: write                          => vtk_output_handler_write
        procedure, public :: allocate_cell_and_nodal_arrays => vtk_output_handler_allocate_cell_and_nodal_arrays
        procedure, public :: topology_to_celltype           => vtk_output_handler_topology_to_celltype
        procedure, public :: append_cell                    => vtk_output_handler_append_cell
        procedure, public :: free                           => vtk_output_handler_free
    end type

public :: vtk_output_handler_t

contains


    subroutine vtk_output_handler_free(this)
        class(vtk_output_handler_t), intent(inout) :: this
        if(allocated(this%X))              call memfree(this%X,              __FILE__, __LINE__)
        if(allocated(this%Y))              call memfree(this%Y,              __FILE__, __LINE__)
        if(allocated(this%Z))              call memfree(this%Z,              __FILE__, __LINE__)
        if(allocated(this%Offset))         call memfree(this%Offset,         __FILE__, __LINE__)
        if(allocated(this%CellTypes))      call memfree(this%CellTypes,      __FILE__, __LINE__)
        if(allocated(this%Connectivities)) call memfree(this%Connectivities, __FILE__, __LINE__)
    end subroutine vtk_output_handler_free


    function vtk_output_handler_topology_to_celltype(this, topology, dimension) result(cell_type)
    !-----------------------------------------------------------------
    !< Translate the topology type of the reference_fe_geo into VTK cell type
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(inout) :: this
        character(len=*),            intent(in)    :: topology
        integer(ip),                 intent(in)    :: dimension
        integer(I1P)                               :: cell_type
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
            write(error_unit,*) 'Topology_to_CellType: Topology not supported ('//trim(adjustl(topology))//')'
            check(.false.)    
        endif
    end function vtk_output_handler_topology_to_celltype


    subroutine vtk_output_handler_allocate_cell_and_nodal_arrays(this)
        class(vtk_output_handler_t), intent(inout) :: this
        integer(ip)                                :: number_nodes
        integer(ip)                                :: number_cells
        assert(.not. allocated(this%X))
        assert(.not. allocated(this%Y))
        assert(.not. allocated(this%Z))
        assert(.not. allocated(this%Offset))
        assert(.not. allocated(this%CellTypes))
        assert(.not. allocated(this%Connectivities))
        number_nodes = this%get_number_nodes()
        number_cells = this%get_number_cells()
        call memalloc(number_nodes, this%X,              __FILE__, __LINE__)
        call memalloc(number_nodes, this%Y,              __FILE__, __LINE__)
        call memalloc(number_nodes, this%Z,              __FILE__, __LINE__)
        call memalloc(number_cells, this%Offset,         __FILE__, __LINE__)
        call memalloc(number_cells, this%CellTypes,      __FILE__, __LINE__)
        call memalloc(number_nodes, this%Connectivities, __FILE__, __LINE__)
        this%Z = 0_rp
    end subroutine vtk_output_handler_allocate_cell_and_nodal_arrays


    subroutine vtk_output_handler_append_cell(this, subcell_iterator)
        class(vtk_output_handler_t), intent(inout) :: this
        type(patch_subcell_iterator_t), intent(in) :: subcell_iterator
        type(output_handler_fe_field_t), pointer   :: field
        integer(ip)                                :: number_vertices
        integer(ip)                                :: number_dimensions
        integer(ip)                                :: number_fields
        integer(ip)                                :: number_components
        integer(ip)                                :: i
        real(rp),    allocatable :: Coordinates(:,:)
        integer(ip), allocatable :: Connectivity(:)
        real(rp),    pointer     :: Value(:,:)


        number_vertices   = subcell_iterator%get_number_vertices()
        number_dimensions = subcell_iterator%get_number_dimensions()
        number_fields     = this%get_number_fields()

        call subcell_iterator%get_coordinates(this%X(this%node_offset+1:this%node_offset+number_vertices), &
                                              this%Y(this%node_offset+1:this%node_offset+number_vertices), &
                                              this%Z(this%node_offset+1:this%node_offset+number_vertices))

        this%Connectivities(this%node_offset+1:this%node_offset+number_vertices) = &
                            (/(i, i=this%node_offset, this%node_offset+number_vertices-1)/)

        do i=1, number_fields
            number_components = subcell_iterator%get_number_field_components(i)
            field => this%get_field(i)
            if(.not. field%value_is_allocated()) call field%allocate_value(number_components, this%get_number_nodes())
            Value => field%get_value()
            call subcell_iterator%get_field(i, number_components, Value(1:number_components,this%node_offset+1:this%node_offset+number_vertices))
        enddo
        this%node_offset = this%node_offset + number_vertices
        this%cell_offset = this%cell_offset + 1

        this%Offset(this%cell_offset) = this%node_offset
        this%CellTypes(this%cell_offset) = this%Topology_to_cellType(subcell_iterator%get_cell_type(), number_dimensions)

    end subroutine vtk_output_handler_append_cell


    subroutine vtk_output_handler_write(this)
        class(vtk_output_handler_t), intent(inout) :: this
        type(output_handler_fe_field_t), pointer   :: field
        integer(ip)                                :: number_fields
        integer(ip)                                :: number_components
        real(rp), pointer                          :: Value(:,:)
        integer(ip)                                :: E_IO, i

        call this%fill_data()
        E_IO = VTK_INI_XML(output_format = 'ascii',     &
                   filename = 'output.vtu',             &
                   mesh_topology = 'UnstructuredGrid',  &
                   cf=this%file_id)

        E_IO = VTK_GEO_XML(NN = this%get_number_nodes(), &
                           NC = this%get_number_cells(), &
                           X  = this%X,                  &
                           Y  = this%Y,                  &
                           Z  = this%Z,                  &
                           cf = this%file_id)

        E_IO = VTK_CON_XML(NC        = this%get_number_cells(), &
                           connect   = this%Connectivities,     &
                           offset    = this%Offset,             &
                           cell_type = this%CellTypes,          &
                           cf        = this%file_id)

        E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=this%file_id)

        number_fields     = this%get_number_fields()
        do i=1, number_fields
            field => this%get_field(i)
            Value => field%get_value()
            number_components = size(Value,1)
            E_IO = VTK_VAR_XML(NC_NN=this%get_number_nodes(), N_COL=number_components, varname=field%get_name(), var=Value, cf=this%file_id)
        enddo

        E_IO = VTK_DAT_XML(var_location='node', var_block_action='close', cf=this%file_id)
        E_IO = VTK_GEO_XML(cf=this%file_id)
        E_IO = VTK_END_XML(cf=this%file_id)

    end subroutine vtk_output_handler_write


    function vtk_output_handler_create_directory(this, path, task, issue_final_barrier) result(error)
    !-----------------------------------------------------------------
    !< The root process create a hierarchy of directories
    !-----------------------------------------------------------------
        class(vtk_handler_t), intent(inout) :: this
        character(len=*),     intent(in)    :: path
        class(serial_fe_space_t), pointer   :: fe_space
        logical                             :: ifb
        integer(kind=c_int)                 :: error
        integer(ip)                         :: default_root_task
    !-----------------------------------------------------------------
        default_root_task =  0
        fe_space          => this%get_fe_space()
        mpi_environment   => fe_space%get_fe_space()
        ifb = .False.; if(present(issue_final_barrier)) ifb = issue_final_barrier

        error = 0
        if(me == default_root_task) then
            error = mkdir_recursive(path//C_NULL_CHAR)
            check ( error == 0 ) 
        end if

        if(ifb) call mpi_environment%l1_barrier()
    end function vtk_output_handler_create_directory

end module vtk_output_handler_names
