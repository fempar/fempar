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
USE xh5_utils_names
USE environment_names
USE output_handler_base_names
USE output_handler_fe_field_names
USE fe_space_names,             only: serial_fe_space_t
USE output_handler_patch_names, only: patch_subcell_iterator_t
USE reference_fe_names,         only: topology_hex, topology_tet


implicit none
#include "debug.i90"
private

    type, extends(output_handler_base_t) :: xh5_output_handler_t
    private 
        type(xh5for_t)                                        :: xh5
        real(rp),                                 allocatable :: X(:)
        real(rp),                                 allocatable :: Y(:)
        real(rp),                                 allocatable :: Z(:)
        integer(ip),                              allocatable :: Connectivities(:)
        type(output_handler_fe_field_1D_value_t), allocatable :: FieldValues(:)
        integer(ip)                                           :: node_offset = 0
        integer(ip)                                           :: cell_offset = 0
    contains
        procedure,                 public :: write                          => xh5_output_handler_write
        procedure,                 public :: allocate_cell_and_nodal_arrays => xh5_output_handler_allocate_cell_and_nodal_arrays
        procedure,                 public :: append_cell                    => xh5_output_handler_append_cell
        procedure,                 public :: free                           => xh5_output_handler_free
    end type

public :: xh5_output_handler_t

contains


    subroutine xh5_output_handler_free(this)
    !-----------------------------------------------------------------
    !< Free procedure
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        integer(ip)                                :: i
    !-----------------------------------------------------------------
        call this%xh5%free()
        if(allocated(this%X))              call memfree(this%X,              __FILE__, __LINE__)
        if(allocated(this%Y))              call memfree(this%Y,              __FILE__, __LINE__)
        if(allocated(this%Z))              call memfree(this%Z,              __FILE__, __LINE__)
        if(allocated(this%Connectivities)) call memfree(this%Connectivities, __FILE__, __LINE__)
        if(allocated(this%FieldValues)) then
            do i=1, size(this%Fieldvalues)
                call this%FieldValues(i)%Free()
            enddo
            deallocate(this%FieldValues)
        endif
        this%node_offset = 0
        this%cell_offset = 0
    end subroutine xh5_output_handler_free


    subroutine xh5_output_handler_allocate_cell_and_nodal_arrays(this)
    !-----------------------------------------------------------------
    !< Allocate cell and nodal arrays
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        integer(ip)                                :: number_nodes
        integer(ip)                                :: number_cells
    !-----------------------------------------------------------------
        assert(.not. allocated(this%X))
        assert(.not. allocated(this%Y))
        assert(.not. allocated(this%Z))
        assert(.not. allocated(this%Connectivities))
        number_nodes = this%get_number_nodes()
        number_cells = this%get_number_cells()
        call memalloc(number_nodes,              this%X,              __FILE__, __LINE__)
        call memalloc(number_nodes,              this%Y,              __FILE__, __LINE__)
        call memalloc(number_nodes,              this%Z,              __FILE__, __LINE__)
        call memalloc(number_nodes+number_cells, this%Connectivities, __FILE__, __LINE__)
        this%Z = 0_rp
    end subroutine xh5_output_handler_allocate_cell_and_nodal_arrays


    subroutine xh5_output_handler_append_cell(this, subcell_iterator)
    !-----------------------------------------------------------------
    !< Append cell data to global arrays
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        type(patch_subcell_iterator_t), intent(in) :: subcell_iterator
        real(rp),    allocatable                   :: Coordinates(:,:)
        integer(ip), allocatable                   :: Connectivity(:)
        real(rp),                        pointer   :: Value(:)
        integer(ip)                                :: number_vertices
        integer(ip)                                :: number_dimensions
        integer(ip)                                :: number_fields
        integer(ip)                                :: number_components
        integer(ip)                                :: node_and_cell_offset
        integer(ip)                                :: i
    !-----------------------------------------------------------------
        number_vertices   = subcell_iterator%get_number_vertices()
        number_dimensions = subcell_iterator%get_number_dimensions()
        number_fields     = this%get_number_fields()

        call subcell_iterator%get_coordinates(this%X(this%node_offset+1:this%node_offset+number_vertices), &
                                              this%Y(this%node_offset+1:this%node_offset+number_vertices), &
                                              this%Z(this%node_offset+1:this%node_offset+number_vertices))

        node_and_cell_offset = this%node_offset+this%cell_offset
        this%Connectivities(node_and_cell_offset+1) = topology_to_xh5_celltype(subcell_iterator%get_cell_type(), number_dimensions)
        node_and_cell_offset = node_and_cell_offset+1

        select case (subcell_iterator%get_cell_type())
            case (topology_hex) 
                select case (number_dimensions)
                    case (2)
                        this%Connectivities(node_and_cell_offset+1:node_and_cell_offset+number_vertices) = &
                                (/0,1,3,2/)+this%node_offset
                    case (3)
                        this%Connectivities(node_and_cell_offset+1:node_and_cell_offset+number_vertices) = &
                                (/0,1,3,2,4,5,7,6/)+this%node_offset
                    case DEFAULT
                        check(.false.)
                end select
            case (topology_tet) 
                this%Connectivities(node_and_cell_offset+1:this%node_offset+number_vertices) = (/(i, i=this%node_offset, this%node_offset+number_vertices-1)/)
            case DEFAULT
                check(.false.)    
        end select

        do i=1, number_fields
            number_components = subcell_iterator%get_number_field_components(i)
            if(.not. this%FieldValues(i)%value_is_allocated()) call this%FieldValues(i)%allocate_value(number_components, this%get_number_nodes())
            Value => this%FieldValues(i)%get_value()
            call subcell_iterator%get_field(i, Value((number_components*this%node_offset)+1:number_components*(this%node_offset+number_vertices)))
        enddo
        this%node_offset = this%node_offset + number_vertices
        this%cell_offset = this%cell_offset + 1

    end subroutine xh5_output_handler_append_cell


    subroutine xh5_output_handler_write(this)
    !-----------------------------------------------------------------
    !< Fill global arrays and Write VTU and PVTU files
    !-----------------------------------------------------------------
        class(xh5_output_handler_t), intent(inout) :: this
        class(serial_fe_space_t),        pointer   :: fe_space
        class(environment_t),            pointer   :: mpi_environment
        type(output_handler_fe_field_t), pointer   :: field
        character(len=:), allocatable              :: path
        character(len=:), allocatable              :: prefix
        real(rp), pointer                          :: Value(:)
        integer(ip)                                :: number_fields
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

        prefix = 'output'

        call this%xh5%Open(FilePrefix = prefix, &
                      Strategy   = XDMF_STRATEGY_CONTIGUOUS_HYPERSLAB, &
                      Action     = XDMF_ACTION_WRITE)
        call this%xh5%SetGrid(NumberOfNodes=this%get_number_nodes(), &
                              NumberOfElements=this%get_number_cells(), &
                              TopologyType=XDMF_TOPOLOGY_TYPE_MIXED, &
                              GeometryType=XDMF_GEOMETRY_TYPE_X_Y_Z)
        call this%xh5%WriteTopology(Connectivities=this%Connectivities)
        call this%xh5%WriteGeometry(X=this%X, Y=this%Y, Z=this%Z)

        do i=1, this%get_number_fields()
            field => this%get_field(i)
            Value => this%FieldValues(i)%get_value()
            call this%xh5%WriteAttribute(Name=field%get_name(), &
                                         Type=XDMF_ATTRIBUTE_TYPE_SCALAR, &
                                         Center=XDMF_ATTRIBUTE_CENTER_NODE , &
                                         Values=Value)
        enddo

        call this%xh5%Close()
        call this%xh5%Free()

    end subroutine xh5_output_handler_write

end module xh5_output_handler_names
