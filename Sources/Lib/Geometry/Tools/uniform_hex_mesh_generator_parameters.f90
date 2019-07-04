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
module uniform_hex_mesh_generator_parameters_names
  use types_names
  implicit none  
  character(len=*), parameter :: struct_hex_mesh_generator_num_dims_key                = 'STRUCT_HEX_MESH_GENERATOR_NUM_DIMS'
  character(len=*), parameter :: struct_hex_mesh_generator_num_levels_key              = 'STRUCT_HEX_MESH_GENERATOR_NUM_LEVELS'
  character(len=*), parameter :: struct_hex_mesh_generator_num_cells_x_dim_key         = 'STRUCT_HEX_MESH_GENERATOR_NUM_CELLS_X_DIM'
  character(len=*), parameter :: struct_hex_mesh_generator_num_parts_x_dim_x_level_key = 'STRUCT_HEX_MESH_GENERATOR_NUM_PARTS_X_DIM_AND_X_LEVEL'
  character(len=*), parameter :: struct_hex_mesh_generator_is_dir_periodic_key         = 'STRUCT_HEX_MESH_GENERATOR_IS_DIR_PERIODIC'
  character(len=*), parameter :: struct_hex_mesh_generator_domain_limits_key           = 'STRUCT_HEX_MESH_GENERATOR_DOMAIN_LIMITS'
  
  character(len=*), parameter :: struct_hex_mesh_generator_num_dims_cla_name                = '--'//struct_hex_mesh_generator_num_dims_key
  character(len=*), parameter :: struct_hex_mesh_generator_num_levels_cla_name              = '--'//struct_hex_mesh_generator_num_levels_key
  character(len=*), parameter :: struct_hex_mesh_generator_num_cells_x_dim_cla_name         = '--'//struct_hex_mesh_generator_num_cells_x_dim_key
  character(len=*), parameter :: struct_hex_mesh_generator_num_parts_x_dim_x_level_cla_name = '--'//struct_hex_mesh_generator_num_parts_x_dim_x_level_key
  character(len=*), parameter :: struct_hex_mesh_generator_is_dir_periodic_cla_name         = '--'//struct_hex_mesh_generator_is_dir_periodic_key
  character(len=*), parameter :: struct_hex_mesh_generator_domain_limits_cla_name           = '--'//struct_hex_mesh_generator_domain_limits_key
  
  character(len=*), parameter :: struct_hex_mesh_generator_num_dims_cla_choices = '2,3'  

  integer(ip),               parameter :: default_struct_hex_mesh_generator_num_levels              = 1
  integer(ip),               parameter :: default_struct_hex_mesh_generator_num_dims                = 2
  real(rp),    dimension(*), parameter :: default_struct_hex_mesh_generator_domain_limits           = [0.0,1.0,0.0,1.0,0.0,1.0]
  integer(ip), dimension(*), parameter :: default_struct_hex_mesh_generator_is_dir_periodic         = [0,0,0]
  integer(ip), dimension(*), parameter :: default_struct_hex_mesh_generator_num_cells_x_dim         = [10,10,10]
  integer(ip), dimension(*), parameter :: default_struct_hex_mesh_generator_num_parts_x_dim_x_level = [1,1,1]

  character(len=*), parameter :: struct_hex_mesh_generator_num_levels_cla_help              = 'Number of levels in the triangulation hierarchy required for the MLBDDC preconditioner'
  character(len=*), parameter :: struct_hex_mesh_generator_num_dims_cla_help                = 'Number of space dimensions'
  character(len=*), parameter :: struct_hex_mesh_generator_domain_limits_cla_help           = 'Domain interval per dimension'
  character(len=*), parameter :: struct_hex_mesh_generator_is_dir_periodic_cla_help         = 'Is the mesh periodic per dimension'
  character(len=*), parameter :: struct_hex_mesh_generator_num_cells_x_dim_cla_help         = 'Number of cells per dimension per level'
  character(len=*), parameter :: struct_hex_mesh_generator_num_parts_x_dim_x_level_cla_help = 'Number of parts per dimension per level'
  
end module uniform_hex_mesh_generator_parameters_names
