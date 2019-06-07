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
module triangulation_parameters_names
  use types_names
  implicit none
  
  ! Parameters to define vef_type. Observe that
  ! mod(vef_type,10)     = dimension
  ! mod(vef_type/10,10)  = interior (0) or boundary (1)
  ! vef_type/100         = local (0), interface(1) or ghost(2)
  integer(ip), parameter :: local_interior_dim0 =  0
  integer(ip), parameter :: local_interior_dim1 =  1
  integer(ip), parameter :: local_interior_dim2 =  2
  integer(ip), parameter :: local_boundary_dim0 = 10
  integer(ip), parameter :: local_boundary_dim1 = 11
  integer(ip), parameter :: local_boundary_dim2 = 12

  integer(ip), parameter :: interface_interior_dim0 = 100
  integer(ip), parameter :: interface_interior_dim1 = 101
  integer(ip), parameter :: interface_interior_dim2 = 102
  integer(ip), parameter :: interface_boundary_dim0 = 110
  integer(ip), parameter :: interface_boundary_dim1 = 111
  integer(ip), parameter :: interface_boundary_dim2 = 112

  integer(ip), parameter :: ghost_dim0 = 200
  integer(ip), parameter :: ghost_dim1 = 201
  integer(ip), parameter :: ghost_dim2 = 202
  
  character(len=*), parameter :: static_triang_generate_from_mesh_data_files           = "mesh_data_files"
  character(len=*), parameter :: static_triang_generate_from_struct_hex_mesh_generator = "struct_hex_mesh_generator"
  
  character(len=*), parameter :: static_triang_generate_from_key                      = 'STATIC_TRIANG_GENERATE_FROM'
  character(len=*), parameter :: static_triang_geometric_interpolation_order_key      = 'STATIC_TRIANG_GEOMETRIC_INTERPOLATION_ORDER'
  character(len=*), parameter :: static_triang_geometric_interpolation_order_cla_name = '--'//static_triang_geometric_interpolation_order_key
  character(len=*), parameter :: static_triang_generate_cla_name                      = '--'//static_triang_generate_from_key
  character(len=*), parameter :: static_triang_generate_cla_choices                   = static_triang_generate_from_mesh_data_files // ' ' // &
                                                                                        static_triang_generate_from_struct_hex_mesh_generator
                                                                                        
  character(len=*), parameter :: static_triang_mesh_data_files_dir_path_key           = 'STATIC_TRIANG_MESH_DATA_FILES_DIR_PATH'
  character(len=*), parameter :: static_triang_mesh_data_files_dir_path_cla_name      = '--'//static_triang_mesh_data_files_dir_path_key
  character(len=*), parameter :: static_triang_mesh_data_files_prefix_key             = 'STATIC_TRIANG_MESH_DATA_FILES_PREFIX'
  character(len=*), parameter :: static_triang_mesh_data_files_prefix_cla_name        = '--'//static_triang_mesh_data_files_prefix_key
                                                                                        
    
  character(len=*), parameter :: subparts_coupling_criteria_key = 'subparts_coupling_criteria'
  character(len=*), parameter :: all_coupled                    = 'all_coupled'
  character(len=*), parameter :: loose_coupling                 = 'loose_coupling' 
  character(len=*), parameter :: strong_coupling                = 'strong_coupling' 
  
  integer(ip), parameter :: refinement = 1 
  integer(ip), parameter :: coarsening = -1 
  integer(ip), parameter :: do_nothing = 0   
  
end module triangulation_parameters_names
