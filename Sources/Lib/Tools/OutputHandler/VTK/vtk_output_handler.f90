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

USE types_names
USE memor_names
USE lib_vtk_io
USE environment_names
USE vtk_utils_names
USE output_handler_base_names
USE output_handler_fe_field_names
USE vtk_parameters_names,       only: default_root_task
USE fe_space_names,             only: serial_fe_space_t
USE output_handler_patch_names, only: patch_subcell_iterator_t
USE IR_Precision,               only: I1P


implicit none
#include "debug.i90"
private

    type, extends(output_handler_base_t) :: vtk_output_handler_t
    private 
        real(rp),                                 allocatable :: X(:)
        real(rp),                                 allocatable :: Y(:)
        real(rp),                                 allocatable :: Z(:)
        integer(ip),                              allocatable :: Offset(:)
        integer(I1P),                             allocatable :: CellTypes(:)
        integer(ip),                              allocatable :: Connectivities(:)
        type(output_handler_fe_field_2D_value_t), allocatable :: FieldValues(:)
        integer(ip)                                           :: node_offset = 0
        integer(ip)                                           :: cell_offset = 0
    contains
        procedure, non_overridable        :: write_vtu                      => vtk_output_handler_write_vtu
        procedure, non_overridable        :: write_pvtu                     => vtk_output_handler_write_pvtu
        procedure,                 public :: write                          => vtk_output_handler_write
        procedure,                 public :: allocate_cell_and_nodal_arrays => vtk_output_handler_allocate_cell_and_nodal_arrays
        procedure,                 public :: append_cell                    => vtk_output_handler_append_cell
        procedure,                 public :: free                           => vtk_output_handler_free
    end type

public :: vtk_output_handler_t

contains


    subroutine vtk_output_handler_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(inout) :: this
        integer(ip)                                :: i
    !-----------------------------------------------------------------
        if(allocated(this%X))              call memfree(this%X,              __FILE__, __LINE__)
        if(allocated(this%Y))              call memfree(this%Y,              __FILE__, __LINE__)
        if(allocated(this%Z))              call memfree(this%Z,              __FILE__, __LINE__)
        if(allocated(this%Offset))         call memfree(this%Offset,         __FILE__, __LINE__)
        if(allocated(this%CellTypes))      call memfree(this%CellTypes,      __FILE__, __LINE__)
        if(allocated(this%Connectivities)) call memfree(this%Connectivities, __FILE__, __LINE__)
        if(allocated(this%FieldValues)) then
            do i=1, size(this%Fieldvalues)
                call this%FieldValues(i)%Free()
            enddo
            deallocate(this%FieldValues)
        endif
        this%node_offset = 0
        this%cell_offset = 0
    end subroutine vtk_output_handler_free


    subroutine vtk_output_handler_allocate_cell_and_nodal_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate cell and nodal arrays
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(inout) :: this
        integer(ip)                                :: number_nodes
        integer(ip)                                :: number_cells
    !-----------------------------------------------------------------
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
    !-----------------------------------------------------------------
    !< Append cell data to global arrays
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(inout) :: this
        type(patch_subcell_iterator_t), intent(in) :: subcell_iterator
        real(rp),    allocatable                   :: Coordinates(:,:)
        integer(ip), allocatable                   :: Connectivity(:)
        real(rp),                        pointer   :: Value(:,:)
        integer(ip)                                :: number_vertices
        integer(ip)                                :: number_dimensions
        integer(ip)                                :: number_fields
        integer(ip)                                :: number_components
        integer(ip)                                :: i
    !-----------------------------------------------------------------
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
            if(.not. this%FieldValues(i)%value_is_allocated()) call this%FieldValues(i)%allocate_value(number_components, this%get_number_nodes())
            Value => this%FieldValues(i)%get_value()
            call subcell_iterator%get_field(i, number_components, Value(1:number_components,this%node_offset+1:this%node_offset+number_vertices))
        enddo
        this%node_offset = this%node_offset + number_vertices
        this%cell_offset = this%cell_offset + 1

        this%Offset(this%cell_offset) = this%node_offset
        this%CellTypes(this%cell_offset) = topology_to_vtk_celltype(subcell_iterator%get_cell_type(), number_dimensions)

    end subroutine vtk_output_handler_append_cell


    subroutine vtk_output_handler_write(this)
    !-----------------------------------------------------------------
    !< Fill global arrays and Write VTU and PVTU files
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(inout) :: this
        class(serial_fe_space_t),        pointer   :: fe_space
        class(environment_t),            pointer   :: mpi_environment
        character(len=:), allocatable              :: path
        character(len=:), allocatable              :: prefix
        integer(ip)                                :: E_IO, i
        integer(ip)                                :: me, np
    !-----------------------------------------------------------------
        allocate(this%FieldValues(this%get_number_fields()))

        call this%fill_data()

        fe_space          => this%get_fe_space()
        assert(associated(fe_space))
        mpi_environment   => fe_space%get_environment()
        assert(associated(mpi_environment))
        call mpi_environment%info(me, np)

        prefix   = 'output'
        path     = get_vtk_output_directory('.', 0._rp)
        E_IO     = if_iam_root_create_directory(path, me)
        check(E_IO == 0)
        call mpi_environment%l1_barrier()

        ! Write VTU
        call this%write_vtu(path, prefix, me)

        if( mpi_environment%am_i_l1_task() .and. me == default_root_task) then
            call this%write_pvtu(path, prefix, np)
        endif

    end subroutine vtk_output_handler_write


    subroutine vtk_output_handler_write_vtu(this, dir_path, prefix, task_id)
    !-----------------------------------------------------------------
    !< Write the vtu file 
    !-----------------------------------------------------------------
        class(vtk_output_handler_t),  intent(in) :: this
        character(len=*),             intent(in) :: dir_path
        character(len=*),             intent(in) :: prefix
        integer(ip),                  intent(in) :: task_id
        character(len=:), allocatable            :: filename
        type(output_handler_fe_field_t), pointer :: field
        integer(ip)                              :: number_fields
        integer(ip)                              :: number_components
        integer(ip)                              :: file_id
        real(rp), pointer                        :: Value(:,:)
        integer(ip)                              :: E_IO, i
    !-----------------------------------------------------------------
        E_IO = 0
        filename = get_vtk_filename(prefix, task_id)
        ! Write VTU
        E_IO = VTK_INI_XML(output_format = 'ascii',    &
                   filename = dir_path//'/'//filename, &
                   mesh_topology = 'UnstructuredGrid', &
                   cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_GEO_XML(NN = this%get_number_nodes(), &
                           NC = this%get_number_cells(), &
                           X  = this%X,                  &
                           Y  = this%Y,                  &
                           Z  = this%Z,                  &
                           cf = file_id)
        assert(E_IO == 0)
        E_IO = VTK_CON_XML(NC        = this%get_number_cells(), &
                           connect   = this%Connectivities,     &
                           offset    = this%Offset,             &
                           cell_type = this%CellTypes,          &
                           cf        = file_id)
        assert(E_IO == 0)
        E_IO = VTK_DAT_XML(var_location='node',var_block_action='open', cf=file_id)

        number_fields     = this%get_number_fields()
        do i=1, number_fields
            field => this%get_field(i)
            Value => this%FieldValues(i)%get_value()
            number_components = size(Value,1)
            E_IO = VTK_VAR_XML(NC_NN=this%get_number_nodes(), N_COL=number_components, varname=field%get_name(), var=Value, cf=file_id)
            assert(E_IO == 0)
        enddo

        E_IO = VTK_DAT_XML(var_location='node', var_block_action='close', cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_GEO_XML(cf=file_id)
        assert(E_IO == 0)
        E_IO = VTK_END_XML(cf=file_id)
        assert(E_IO == 0)
    end subroutine vtk_output_handler_write_vtu


    subroutine vtk_output_handler_write_pvtu(this, dir_path, prefix, number_parts)
    !-----------------------------------------------------------------
    !< Write the pvtu file containing the number of parts
    !< (only root processor)
    !-----------------------------------------------------------------
        class(vtk_output_handler_t), intent(in)  :: this
        character(len=*),            intent(in)  :: dir_path
        character(len=*),             intent(in) :: prefix
        integer(ip),                  intent(in) :: number_parts
        character(len=:), allocatable            :: filename
        class(serial_fe_space_t),        pointer :: fe_space
        class(environment_t),            pointer :: mpi_environment
        type(output_handler_fe_field_t), pointer :: field
        integer(ip)                              :: file_id
        integer(ip)                              :: E_IO
        integer(ip)                              :: me, np, i
    !-----------------------------------------------------------------
        E_IO = 0
        filename = trim(adjustl(get_pvtu_filename(prefix)))

        E_IO = PVTK_INI_XML(filename      = dir_path//'/'//filename,  &
                            mesh_topology = 'PUnstructuredGrid',      &
                            tp            = 'Float64',                &
                            cf            = file_id)
        assert(E_IO == 0)
        do i=0, number_parts-1
            E_IO = PVTK_GEO_XML(source=trim(adjustl(get_vtk_filename(prefix, i))), cf=file_id)
            assert(E_IO == 0)
        enddo
        E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'OPEN', cf=file_id)
        ! Write point data fields
        do i=1, this%get_number_fields()
            field => this%get_field(i)
            E_IO = PVTK_VAR_XML(varname = trim(adjustl(field%get_name())),             &
                                tp      = 'Float64',                                   &
                                Nc      = this%FieldValues(i)%get_number_components(), &
                                cf      = file_id )
            assert(E_IO == 0)
        enddo
        E_IO = PVTK_DAT_XML(var_location = 'Node', var_block_action = 'CLOSE', cf=file_id)
        assert(E_IO == 0)
        E_IO = PVTK_END_XML(cf = file_id)
        assert(E_IO == 0)

    end subroutine vtk_output_handler_write_pvtu


end module vtk_output_handler_names
