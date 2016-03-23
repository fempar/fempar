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
module serial_names
  ! Tools
  use types_names
  use memor_names
  use allocatable_array_names
  use stdio_names
  use hash_table_names
  use postpro_names
  use environment_names
  use serial_environment_names
  use Data_Type_Command_Line_Interface
  use FPL

  ! Geometry
  use map_names
  use map_apply_names
  use triangulation_names
  use mesh_to_triangulation_names

  use triangulation_names
  use mesh_to_triangulation_names
  use mesh_distribution_names
  use create_mesh_distribution_names
  use generate_uniform_triangulation_names

  use element_import_names
  use element_import_create_names
  use conditions_names
  use conditions_io_names
  use materials_names
  use materials_io_names
  use elem_to_subset_id_names
  use graph_names
  use generate_vefs_mesh_conditions_names
  use renumbering_names
  use mesh_names
  use partitioning_params_names
  use create_mesh_distribution_names

  use mesh_io_names
  use migratory_element_names

  ! Linear Algebra
  use iterative_linear_solver_names
  use sparse_matrix_names
  use serial_scalar_array_names
  use serial_block_array_names
  use operator_names
  use vector_space_names
  use vector_names
  use matrix_names
  use array_names
  use block_preconditioner_l_names
  use block_preconditioner_u_names
  use block_preconditioner_lu_names
  use block_operator_names
  use block_vector_names
  use direct_solver_names
  use direct_solver_parameters_names
  use direct_solver_creational_methods_dictionary_names
  
  ! Integration 
  use reference_fe_names
  use reference_fe_factory_names
  use field_names
  use serial_fe_space_names
  use discrete_integration_names
  use assembler_names
  use fe_affine_operator_names
  use function_names
  use function_library_names
end module serial_names
