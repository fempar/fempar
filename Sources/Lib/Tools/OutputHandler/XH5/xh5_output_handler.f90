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

module xh5_output_handler_names

USE types_names
USE memor_names
USE xh5for
USE FPL
USE xh5_utils_names
USE xh5_parameters_names
USE environment_names
USE execution_context_names
USE mpi_context_names
USE base_output_handler_names
USE output_handler_parameters_names
USE output_handler_fe_field_names
USE fe_space_names,             only: serial_fe_space_t
USE output_handler_patch_names, only: patch_subcell_accessor_t
USE reference_fe_names,         only: topology_hex, topology_tet


implicit none
#include "debug.i90"
private

    type, extends(base_output_handler_t) :: xh5_output_handler_t
    private 
        type(xh5for_t)                                        :: xh5
        character(:), allocatable                             :: FilePrefix
        character(:), allocatable                             :: Path
        character(:), allocatable                             :: CellType
        logical                                               :: StaticGrid = .false.
        integer(ip)                                           :: GridType
        integer(ip)                                           :: Strategy
        integer(ip)                                           :: Action
        real(rp),                                 allocatable :: XYZ(:)
        integer(ip),                              allocatable :: Connectivities(:)
        type(output_handler_fe_field_1D_value_t), allocatable :: FieldValues(:)
        type(output_handler_fe_field_1D_value_t), allocatable :: CellValues(:)
        integer(ip)                                           :: number_steps = 0
        integer(ip)                                           :: node_offset = 0
        integer(ip)                                           :: cell_offset = 0
    contains
    private
        procedure,                 public :: open_body                      => xh5_output_handler_open_body
        procedure,                 public :: append_time_step               => xh5_output_handler_append_time_step
        procedure                         :: allocate_cell_and_nodal_arrays => xh5_output_handler_allocate_cell_and_nodal_arrays
        procedure                         :: append_cell                    => xh5_output_handler_append_cell
        procedure,                 public :: write                          => xh5_output_handler_write
        procedure,                 public :: close_body                     => xh5_output_handler_close_body
        procedure                         :: free_body                      => xh5_output_handler_free_body
    end type

public :: xh5_output_handler_t

contains


    subroutine xh5_output_handler_free_body(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        integer(ip)                                :: i
    !-----------------------------------------------------------------
        call this%xh5%free()
        if(allocated(this%Path))           deallocate(this%Path)
        if(allocated(this%FilePrefix))     deallocate(this%FilePrefix)
        if(allocated(this%XYZ))            call memfree(this%XYZ,            __FILE__, __LINE__)
        if(allocated(this%Connectivities)) call memfree(this%Connectivities, __FILE__, __LINE__)
        if(allocated(this%FieldValues)) then
            do i=1, size(this%Fieldvalues)
                call this%FieldValues(i)%Free()
            enddo
            deallocate(this%FieldValues)
        endif
        if(allocated(this%CellValues)) then
            do i=1, size(this%CellValues)
                call this%CellValues(i)%Free()
            enddo
            deallocate(this%CellValues)
        endif
        this%number_steps = 0
        this%node_offset  = 0
        this%cell_offset  = 0
    end subroutine xh5_output_handler_free_body


    subroutine xh5_output_handler_open_body(this, dir_path, prefix, parameter_list)
    !-----------------------------------------------------------------
    !< Open xh5for_t derive dtype. Set parameters from parameter list
    !-----------------------------------------------------------------
        class(xh5_output_handler_t),     intent(inout) :: this
        character(len=*),                intent(in)    :: dir_path
        character(len=*),                intent(in)    :: prefix
        type(ParameterList_t), optional, intent(in)    :: parameter_list
        class(serial_fe_space_t),        pointer       :: fe_space
        type(environment_t),             pointer       :: environment
        class(execution_context_t),      pointer       :: cntxt
        integer(ip)                                    :: mpi_info
        integer(ip)                                    :: mpi_comm
        integer(ip)                                    :: FPLError
    !-----------------------------------------------------------------
        fe_space    => this%get_fe_space()
        assert(associated(fe_space))
        environment => fe_space%get_environment()
        assert(associated(environment))

        if( environment%am_i_l1_task()) then

            this%Path       = dir_path
            this%FilePrefix = prefix

            ! Set defaults
            this%StaticGrid = xh5_default_StaticGrid
            this%Strategy   = xh5_default_Strategy
            this%GridType   = xh5_default_GridType
            this%Action     = xh5_default_Action
            mpi_info        = xh5_default_Info
            mpi_comm        = xh5_default_Comm

            cntxt => environment%get_l1_context()

            select type (cntxt)
                type is (mpi_context_t)
                    mpi_comm  =  cntxt%get_icontxt()
            end select

            if(present(parameter_list)) then
                ! Get StaticGrid value from parameter_list
                if(parameter_list%isPresent(oh_staticgrid)) then
                    assert(parameter_list%isAssignable(oh_StaticGrid, this%StaticGrid))
                    FPLError   = parameter_list%Get(Key=oh_StaticGrid, Value=this%StaticGrid)
                    assert(FPLError == 0)
                endif

                ! Get Strategy value from parameter_list
                if(parameter_list%isPresent(xh5_Strategy)) then
                    assert(parameter_list%isAssignable(xh5_Strategy, this%Strategy))
                    FPLError   = parameter_list%Get(Key=xh5_Strategy, Value=this%Strategy)
                    assert(FPLError == 0)
                endif

                ! Get Info value from parameter_list
                if(parameter_list%isPresent(xh5_Info)) then
                    assert(parameter_list%isAssignable(xh5_info, mpi_info))
                    FPLError   = parameter_list%Get(Key=xh5_Info, Value=mpi_info)
                    assert(FPLError == 0)
                endif
            endif

            call this%xh5%Open(FilePrefix = this%FilePrefix, &
                               Path       = this%Path,       &
                               GridType   = this%GridType,   &
                               StaticGrid = this%StaticGrid, &
                               Strategy   = this%Strategy,   &
                               Action     = this%Action,     &
                               Comm       = mpi_comm,        &
                               Info       = mpi_info)
        endif

    end subroutine xh5_output_handler_open_body


    subroutine xh5_output_handler_append_time_step(this, value)
    !-----------------------------------------------------------------
    !< Append a new time step
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        real(rp),                    intent(in)    :: value
        class(serial_fe_space_t),        pointer   :: fe_space
        type(environment_t),             pointer   :: environment
    !-----------------------------------------------------------------
        fe_space          => this%get_fe_space()
        assert(associated(fe_space))
        environment   => fe_space%get_environment()
        assert(associated(environment))

        if( environment%am_i_l1_task()) call this%xh5%AppendStep(Value)
        this%number_steps = this%number_steps + 1
        this%node_offset  = 0
        this%cell_offset  = 0
    end subroutine xh5_output_handler_append_time_step


    subroutine xh5_output_handler_allocate_cell_and_nodal_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate cell and nodal arrays
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        integer(ip)                                :: number_nodes
        integer(ip)                                :: number_cells
        integer(ip)                                :: number_dimensions
    !-----------------------------------------------------------------
        number_nodes      = this%get_number_nodes()
        number_cells      = this%get_number_cells()
        number_dimensions = this%get_number_dimensions()
        if(allocated(this%XYZ)) then
            call memrealloc(number_nodes*number_dimensions, this%XYZ,            __FILE__, __LINE__)
        else
            call memalloc(  number_nodes*number_dimensions, this%XYZ,            __FILE__, __LINE__)
        endif
        if(this%has_mixed_cell_topologies()) then
            if(allocated(this%Connectivities)) then
                call memrealloc(number_nodes+number_cells,  this%Connectivities, __FILE__, __LINE__)
            else
                call memalloc(  number_nodes+number_cells,  this%Connectivities, __FILE__, __LINE__)
            endif
        else
            if(allocated(this%Connectivities)) then
                call memrealloc(number_nodes,               this%Connectivities, __FILE__, __LINE__)
            else
                call memalloc(  number_nodes,               this%Connectivities, __FILE__, __LINE__)
            endif
        endif
    end subroutine xh5_output_handler_allocate_cell_and_nodal_arrays


    subroutine xh5_output_handler_append_cell(this, subcell_accessor)
    !-----------------------------------------------------------------
    !< Append cell data to global arrays
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        type(patch_subcell_accessor_t), intent(in) :: subcell_accessor
        real(rp),    allocatable                   :: Coordinates(:,:)
        integer(ip), allocatable                   :: Connectivity(:)
        real(rp),                        pointer   :: Value(:)
        integer(ip)                                :: number_vertices
        integer(ip)                                :: number_dimensions
        integer(ip)                                :: number_fields
        integer(ip)                                :: number_cell_vectors
        integer(ip)                                :: number_components
        integer(ip)                                :: connectivities_offset
        integer(ip)                                :: i
    !-----------------------------------------------------------------
        number_vertices   = subcell_accessor%get_number_vertices()
        number_dimensions = subcell_accessor%get_number_dimensions()
        number_fields     = this%get_number_fields()
        number_cell_vectors = this%get_number_cell_vectors()


        if(.not. this%StaticGrid .or. this%number_steps <= 1) then
            call subcell_accessor%get_coordinates(this%XYZ( &
                                this%node_offset*number_dimensions+1:(this%node_offset+number_vertices)*number_dimensions) )

            if(this%has_mixed_cell_topologies()) then
                ! Add element topology ID before its connectivities
                connectivities_offset = this%node_offset+this%cell_offset
                this%Connectivities(connectivities_offset+1) = topology_to_xh5_celltype(subcell_accessor%get_cell_type(), number_dimensions)
                connectivities_offset = connectivities_offset+1
            else
                connectivities_offset = this%node_offset
                if(.not. allocated(this%CellType)) this%CellType = subcell_accessor%get_cell_type()
            endif

            select case (subcell_accessor%get_cell_type())
                case (topology_hex) 
                    select case (number_dimensions)
                        case (2)
                            this%Connectivities(connectivities_offset+1:connectivities_offset+number_vertices) = &
                                    [0,1,3,2]+this%node_offset
                        case (3)
                            this%Connectivities(connectivities_offset+1:connectivities_offset+number_vertices) = &
                                    [0,1,3,2,4,5,7,6]+this%node_offset
                        case DEFAULT
                            check(.false.)
                    end select
                case (topology_tet) 
                    this%Connectivities(connectivities_offset+1:connectivities_offset+number_vertices) = &
                                    [(i, i=this%node_offset, this%node_offset+number_vertices-1)]
                case DEFAULT
                    check(.false.)    
            end select
        endif

        do i=1, number_fields
            number_components = subcell_accessor%get_number_field_components(i)
            if(.not. this%FieldValues(i)%value_is_allocated()) call this%FieldValues(i)%allocate_value(number_components, this%get_number_nodes())
            Value => this%FieldValues(i)%get_value()
            call subcell_accessor%get_field(i, Value((number_components*this%node_offset)+1:number_components*(this%node_offset+number_vertices)))
        enddo

        this%node_offset = this%node_offset + number_vertices
        this%cell_offset = this%cell_offset + 1

        do i=1, number_cell_vectors
            if(.not. this%CellValues(i)%value_is_allocated()) call this%CellValues(i)%allocate_value(1, this%get_number_cells())
            Value => this%CellValues(i)%get_value()
            call subcell_accessor%get_cell_vector(i, Value(this%cell_offset:this%cell_offset))
        enddo

    end subroutine xh5_output_handler_append_cell


    subroutine xh5_output_handler_write(this)
    !-----------------------------------------------------------------
    !< Fill global arrays and Write VTU and PVTU files
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        class(serial_fe_space_t),           pointer   :: fe_space
        type(environment_t),                pointer   :: environment
        type(output_handler_fe_field_t),    pointer   :: field
        type(output_handler_cell_vector_t), pointer   :: cell_vector
        real(rp), pointer                             :: Value(:)
        integer(ip)                                   :: attribute_type
        integer(ip)                                   :: geometry_type
        integer(ip)                                   :: topology_type
        integer(ip)                                   :: E_IO, i
        integer(ip)                                   :: me, np
    !-----------------------------------------------------------------
        fe_space          => this%get_fe_space()
        assert(associated(fe_space))
        environment   => fe_space%get_environment()
        assert(associated(environment))

        if( environment%am_i_l1_task()) then
            if(.not. allocated(this%FieldValues)) allocate(this%FieldValues(this%get_number_fields()))
            if(.not. allocated(this%CellValues)) allocate(this%CellValues(this%get_number_cell_vectors()))

            call this%fill_data(update_mesh = (.not. this%StaticGrid .or. this%number_steps <= 1))

            call environment%info(me, np)

            if(.not. this%StaticGrid .or. this%number_steps <= 1) then                 
                geometry_type = dimensions_to_xh5_unstructured_GeometryType(this%get_number_dimensions())
                if(allocated(this%CellType) .and.  .not. this%has_mixed_cell_topologies()) then
                    topology_type = topology_to_xh5_topologytype(this%CellType, this%get_number_dimensions())
                else
                    topology_type = XDMF_TOPOLOGY_TYPE_MIXED
                endif

                call this%xh5%SetGrid(NumberOfNodes        = this%get_number_nodes(),  &
                                      NumberOfElements     = this%get_number_cells(),  &
                                      TopologyType         = topology_type,            &
                                      GeometryType         = geometry_type)

                call this%xh5%WriteGeometry(XYZ = this%XYZ)

                call this%xh5%WriteTopology(Connectivities = this%Connectivities)
            endif

            do i=1, this%get_number_fields()
                field => this%get_fe_field(i)
                Value => this%FieldValues(i)%get_value()
                attribute_type = number_components_to_xh5_AttributeType(this%FieldValues(i)%get_number_components())
                call this%xh5%WriteAttribute(Name   = field%get_name(),             &
                                             Type   = attribute_type,               &
                                             Center = XDMF_ATTRIBUTE_CENTER_NODE ,  &
                                             Values = Value)
            enddo

            do i=1, this%get_number_cell_vectors()
                cell_vector => this%get_cell_vector(i)
                Value => this%CellValues(i)%get_value()
                call this%xh5%WriteAttribute(Name   = cell_vector%get_name(),       &
                                             Type   = XDMF_ATTRIBUTE_TYPE_SCALAR,   &
                                             Center = XDMF_ATTRIBUTE_CENTER_CELL ,  &
                                             Values = Value)
            enddo

            call this%xh5%Serialize()
        endif
    end subroutine xh5_output_handler_write


    subroutine xh5_output_handler_close_body(this)
    !-----------------------------------------------------------------
    !< Close xh5for_t derived type
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        class(serial_fe_space_t),        pointer   :: fe_space
        type(environment_t),             pointer   :: environment
    !-----------------------------------------------------------------
        fe_space          => this%get_fe_space()
        assert(associated(fe_space))
        environment   => fe_space%get_environment()
        assert(associated(environment))

        if( environment%am_i_l1_task()) call this%xh5%close()
    end subroutine xh5_output_handler_close_body

end module xh5_output_handler_names
