
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
!### Public parameters of the [[xh5_output_handler_names(module)]] module
!
! Contains the following public entities:
! [[xh5_parameters_names(module)]]
!---------------------------------------------------------------------
module xh5_parameters_names
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Public parameters of the [[xh5_output_handler_names(module)]] module
!
! Contains the following public parameters:
! [[xh5_parameters_names(module):xh5_default_StaticGrid(variable)]], 
! [[xh5_parameters_names(module):xh5_default_Strategy(variable)]], 
! [[xh5_parameters_names(module):xh5_default_GridType(variable)]], 
! [[xh5_parameters_names(module):xh5_default_Action(variable)]], 
! [[xh5_parameters_names(module):xh5_default_Info(variable)]], 
! [[xh5_parameters_names(module):xh5_default_Comm(variable)]]
! and also several XH5For parameters
!--------------------------------------------------------------------- 
USE types_names
USE xh5for, only: XDMF_GRID_TYPE_CURVILINEAR, XDMF_GRID_TYPE_RECTILINEAR,               &
                  XDMF_GRID_TYPE_REGULAR, XDMF_GRID_TYPE_UNSTRUCTURED,                  &
                  XDMF_GEOMETRY_TYPE_XY, XDMF_GEOMETRY_TYPE_XYZ,                        &
                  XDMF_TOPOLOGY_TYPE_TRIANGLE, XDMF_TOPOLOGY_TYPE_QUADRILATERAL,        &
                  XDMF_TOPOLOGY_TYPE_TETRAHEDRON, XDMF_TOPOLOGY_TYPE_HEXAHEDRON,        &
                  XDMF_TOPOLOGY_TYPE_MIXED,                                             &
                  XDMF_ATTRIBUTE_TYPE_SCALAR, XDMF_ATTRIBUTE_TYPE_VECTOR,               &
                  XDMF_ATTRIBUTE_TYPE_TENSOR, XDMF_ATTRIBUTE_TYPE_TENSOR6,              &
                  XDMF_STRATEGY_CONTIGUOUS_HYPERSLAB,XDMF_STRATEGY_DATASET_PER_PROCESS, &
                  XDMF_ACTION_READ, XDMF_ACTION_WRITE


implicit none
private

    ! VTK cell type parameters
    integer(ip), parameter, public :: XDMF_NOTOPOLOGY               = 0   ! p(new XdmfTopologyType(0, 0, faces, 0, "NoTopology", NoCellType, 0x0));
    integer(ip), parameter, public :: XDMF_POLYVERTEX               = 1   ! p(new XdmfTopologyType(1, 0, faces, 0, "Polyvertex", Linear, 0x1));
    integer(ip), parameter, public :: XDMF_POLYLINE                 = 2   ! p(new XdmfTopologyType(nodesPerElement, 0, faces, nodesPerElement - 1,"Polyline", Linear, 0x2));
    integer(ip), parameter, public :: XDMF_POLYGON                  = 3   ! p(new XdmfTopologyType(nodesPerElement, 1, faces, nodesPerElement, "Polygon", Linear, 0x3));
    integer(ip), parameter, public :: XDMF_TRIANGLE                 = 4   ! p(new XdmfTopologyType(3, 1, faces, 3, "Triangle", Linear, 0x4));
    integer(ip), parameter, public :: XDMF_QUADRILATERAL            = 5   ! p(new XdmfTopologyType(4, 1, faces, 4, "Quadrilateral", Linear, 0x5));
    integer(ip), parameter, public :: XDMF_TETRAHEDRON              = 6   ! p(new XdmfTopologyType(4, 4, faces, 6, "Tetrahedron", Linear, 0x6));
    integer(ip), parameter, public :: XDMF_PYRAMID                  = 7   ! p(new XdmfTopologyType(5, 5, faces, 8, "Pyramid", Linear, 0x7));
    integer(ip), parameter, public :: XDMF_WEDGE                    = 8   ! p(new XdmfTopologyType(6, 5, faces, 9, "Wedge", Linear, 0x8));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON               = 9   ! p(new XdmfTopologyType(8, 6, faces, 12, "Hexahedron", Linear, 0x9));
    integer(ip), parameter, public :: XDMF_POLYHEDRON               = 16  ! p(new XdmfTopologyType(0, 0, faces, 0, "Polyhedron", Linear, 0x10));
    integer(ip), parameter, public :: XDMF_EDGE_3                   = 34  ! p(new XdmfTopologyType(3, 0, faces, 1, "Edge_3", Quadratic, 0x22));
    integer(ip), parameter, public :: XDMF_QUADRILATERAL_9          = 35  ! p(new XdmfTopologyType(9, 1, faces, 4, "Quadrilateral_9", Quadratic, 0x23));
    integer(ip), parameter, public :: XDMF_TRIANGLE_6               = 36  ! p(new XdmfTopologyType(6, 1, faces, 3, "Triangle_6", Quadratic, 0x24));
    integer(ip), parameter, public :: XDMF_QUADRILATERAL_8          = 37  ! p(new XdmfTopologyType(8, 1, faces, 4, "Quadrilateral_8", Quadratic, 0x25));
    integer(ip), parameter, public :: XDMF_TETRAHEDRON_10           = 38  ! p(new XdmfTopologyType(10, 4, faces, 6, "Tetrahedron_10", Quadratic, 0x26));
    integer(ip), parameter, public :: XDMF_PYRAMID_13               = 39  ! p(new XdmfTopologyType(13, 5, faces, 8, "Pyramid_13", Quadratic, 0x27));
    integer(ip), parameter, public :: XDMF_WEDGE_15                 = 40  ! p(new XdmfTopologyType(15, 5, faces, 9, "Wedge_15", Quadratic, 0x28));
    integer(ip), parameter, public :: XDMF_WEDGE_18                 = 41  ! p(new XdmfTopologyType(18, 5, faces, 9, "Wedge_18", Quadratic, 0x29));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_20            = 48  ! p(new XdmfTopologyType(20, 6, faces, 12, "Hexahedron_20", Quadratic, 0x30));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_24            = 49  ! p(new XdmfTopologyType(24, 6, faces, 12, "Hexahedron_24", Quadratic, 0x31));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_27            = 50  ! p(new XdmfTopologyType(27, 6, faces, 12, "Hexahedron_27", Quadratic, 0x32));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_64            = 51  ! p(new XdmfTopologyType(64, 6, faces, 12, "Hexahedron_64", Cubic, 0x33));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_125           = 52  ! p(new XdmfTopologyType(125, 6, faces, 12, "Hexahedron_125", Quartic, 0x34));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_216           = 53  ! p(new XdmfTopologyType(216, 6, faces, 12, "Hexahedron_216", Quintic, 0x35));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_343           = 54  ! p(new XdmfTopologyType(343, 6, faces, 12, "Hexahedron_343", Sextic, 0x36));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_512           = 55  ! p(new XdmfTopologyType(512, 6, faces, 12, "Hexahedron_512", Septic, 0x37));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_729           = 56  ! p(new XdmfTopologyType(729, 6, faces, 12, "Hexahedron_729", Octic, 0x38));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_1000          = 57  ! p(new XdmfTopologyType(1000, 6, faces, 12, "Hexahedron_1000", Nonic, 0x39));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_1331          = 64  ! p(new XdmfTopologyType(1331, 6, faces, 12, "Hexahedron_1331", Decic, 0x40));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_64   = 65  ! p(new XdmfTopologyType(64, 6, faces, 12, "Hexahedron_Spectral_64", Cubic, 0x41));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_125  = 66  ! p(new XdmfTopologyType(125, 6, faces, 12, "Hexahedron_Spectral_125", Quartic, 0x42));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_216  = 67  ! p(new XdmfTopologyType(216, 6, faces, 12, "Hexahedron_Spectral_216", Quintic, 0x43));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_343  = 68  ! p(new XdmfTopologyType(343, 6, faces, 12, "Hexahedron_Spectral_343", Sextic, 0x44));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_512  = 69  ! p(new XdmfTopologyType(512, 6, faces, 12, "Hexahedron_Spectral_512", Septic, 0x45));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_729  = 70  ! p(new XdmfTopologyType(729, 6, faces, 12, "Hexahedron_Spectral_729", Octic, 0x46));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_1000 = 71  ! p(new XdmfTopologyType(1000, 6, faces, 12, "Hexahedron_Spectral_1000", Nonic, 0x47));
    integer(ip), parameter, public :: XDMF_HEXAHEDRON_SPECTRAL_1331 = 72  ! p(new XdmfTopologyType(1331, 6, faces, 12, "Hexahedron_Spectral_1331", Decic, 0x48));
    integer(ip), parameter, public :: XDMF_MIXED                    = 112 ! p(new XdmfTopologyType(0, 0, faces, 0, "Mixed", Arbitrary, 0x70));

    ! PARAMETERS IDENTIFIERS
    character(*), parameter, public :: xh5_Strategy   = 'xh5_Strategy'
    character(*), parameter, public :: xh5_GridType   = 'xh5_GridType'
    character(*), parameter, public :: xh5_Action     = 'xh5_Action'
    character(*), parameter, public :: xh5_Info       = 'xh5_Info'

    ! GRID TYPE PARAMETERS (from xh5for)
    public :: XDMF_GRID_TYPE_CURVILINEAR
    public :: XDMF_GRID_TYPE_RECTILINEAR
    public :: XDMF_GRID_TYPE_REGULAR
    public :: XDMF_GRID_TYPE_UNSTRUCTURED

    ! GEOMETRY TYPE PARAMETERS (from xh5for)
    public :: XDMF_GEOMETRY_TYPE_XY
    public :: XDMF_GEOMETRY_TYPE_XYZ

    ! TOPOLOGY TYPE PARAMETERS (from xh5for)
    public :: XDMF_TOPOLOGY_TYPE_TRIANGLE
    public :: XDMF_TOPOLOGY_TYPE_QUADRILATERAL
    public :: XDMF_TOPOLOGY_TYPE_TETRAHEDRON
    public :: XDMF_TOPOLOGY_TYPE_HEXAHEDRON
    public :: XDMF_TOPOLOGY_TYPE_MIXED

    ! ATTRIBUTE TYPE PARAMETERS (from xh5for)
    public :: XDMF_ATTRIBUTE_TYPE_SCALAR
    public :: XDMF_ATTRIBUTE_TYPE_VECTOR
    public :: XDMF_ATTRIBUTE_TYPE_TENSOR
    public :: XDMF_ATTRIBUTE_TYPE_TENSOR6

    ! STRATEGY PARAMETERS (from xh5for)
    public :: XDMF_STRATEGY_CONTIGUOUS_HYPERSLAB
    public :: XDMF_STRATEGY_DATASET_PER_PROCESS

    ! ACTION PARAMETERS (from xh5for)
    public :: XDMF_ACTION_READ
    public :: XDMF_ACTION_WRITE

    ! DEFAULT PARAMETERS
    logical,     parameter, public :: xh5_default_StaticGrid = .true.
    integer(ip), parameter, public :: xh5_default_Strategy   = XDMF_STRATEGY_CONTIGUOUS_HYPERSLAB
    integer(ip), parameter, public :: xh5_default_GridType   = XDMF_GRID_TYPE_UNSTRUCTURED
    integer(ip), parameter, public :: xh5_default_Action     = XDMF_ACTION_WRITE
    integer(ip), parameter, public :: xh5_default_Info       = 0
    integer(ip), parameter, public :: xh5_default_Comm       = 0


end module xh5_parameters_names
