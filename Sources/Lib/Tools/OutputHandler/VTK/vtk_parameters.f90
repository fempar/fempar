
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

module vtk_parameters_names

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
    character(len=*), parameter, public :: vtk_format            = 'vtk_format'

    ! DEFAULT PARAMETERS
    integer(ip),      parameter, public :: vtk_default_root_task             = 0
    integer(ip),      parameter, public :: vtk_default_number_of_tasks       = 1
    integer(ip),      parameter, public :: vtk_default_guess_number_of_steps = 100
    real(rp),         parameter, public :: vtk_default_step_value            = 0.0_rp
    character(len=3), parameter, public :: vtk_default_format                = 'raw'
    logical,          parameter, public :: vtk_default_StaticGrid            = .false.

end module vtk_parameters_names
