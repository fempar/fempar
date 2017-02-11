
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
!### Public parameters of the [[vtk_output_handler_names(module)]] module
!
! Contains the following public entities:
! [[vtk_parameters_names(module)]]
!---------------------------------------------------------------------
module vtk_parameters_names
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
! Contains the following public parameters:
! [[vtk_parameters_names(module):vtk_format(variable]], 
! [[vtk_parameters_names(module):vtk_format_ascii(variable]], 
! [[vtk_parameters_names(module):vtk_format_raw(variable]], 
! [[vtk_parameters_names(module):vtk_format_binary_appended(variable]], 
! [[vtk_parameters_names(module):vtk_format_binary(variable]], 
! [[vtk_parameters_names(module):vtk_default_root_task(variable]], 
! [[vtk_parameters_names(module):vtk_default_number_of_tasks(variable]], 
! [[vtk_parameters_names(module):vtk_default_guess_number_of_steps(variable]], 
! [[vtk_parameters_names(module):vtk_default_step_value(variable]], 
! [[vtk_parameters_names(module):vtk_default_format(variable]], 
! [[vtk_parameters_names(module):vtk_default_StaticGrid(variable]]
! and also several VTK parameters
!--------------------------------------------------------------------- 
!### Public parameters of the [[vtk_output_handler_names(module)]] module
USE types_names
USE IR_Precision, only: I1P

implicit none
private

    ! VTK cell type parameters
    integer(ip), parameter, public :: vtk_vertex               = 1_I1P
    integer(ip), parameter, public :: vtk_poly_vertex          = 2_I1P
    integer(ip), parameter, public :: vtk_line                 = 3_I1P
    integer(ip), parameter, public :: vtk_poly_line            = 4_I1P
    integer(ip), parameter, public :: vtk_triangle             = 5_I1P
    integer(ip), parameter, public :: vtk_triangle_strip       = 6_I1P
    integer(ip), parameter, public :: vtk_polygon              = 7_I1P
    integer(ip), parameter, public :: vtk_pixel                = 8_I1P
    integer(ip), parameter, public :: vtk_quad                 = 9_I1P
    integer(ip), parameter, public :: vtk_tetra                = 10_I1P
    integer(ip), parameter, public :: vtk_voxel                = 11_I1P
    integer(ip), parameter, public :: vtk_hexahedron           = 12_I1P
    integer(ip), parameter, public :: vtk_wedge                = 13_I1P
    integer(ip), parameter, public :: vtk_pyramid              = 14_I1P
    integer(ip), parameter, public :: vtk_quadratic_edge       = 21_I1P
    integer(ip), parameter, public :: vtk_quadratic_triangle   = 22_I1P
    integer(ip), parameter, public :: vtk_quadratic_quad       = 23_I1P
    integer(ip), parameter, public :: vtk_quadratic_tetra      = 24_I1P
    integer(ip), parameter, public :: vtk_quadratic_hexahedron = 25_I1P

    ! PARAMETERS IDENTIFIERS
    character(len=*), parameter, public :: vtk_format                  = 'vtk_format'

    ! VTK_FORMAT PARAMETERS
    character(len=*), parameter, public :: vtk_format_ascii            = 'ASCII'
    character(len=*), parameter, public :: vtk_format_raw              = 'RAW'
    character(len=*), parameter, public :: vtk_format_binary_appended  = 'BINARY-APPENDED'
    character(len=*), parameter, public :: vtk_format_binary           = 'BINARY'

    ! DEFAULT PARAMETERS
    integer(ip),      parameter, public :: vtk_default_root_task             = 0
    integer(ip),      parameter, public :: vtk_default_number_of_tasks       = 1
    integer(ip),      parameter, public :: vtk_default_guess_number_of_steps = 100
    real(rp),         parameter, public :: vtk_default_step_value            = 0.0_rp
    character(len=*), parameter, public :: vtk_default_format                = vtk_format_raw
    logical,          parameter, public :: vtk_default_StaticGrid            = .true.

end module vtk_parameters_names
