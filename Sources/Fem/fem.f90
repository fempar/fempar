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
module fem_names

  ! Tools
  use types_names
  use memor_names
  use array_names
  use stdio_names
  use hash_table_names
  use postpro_names
  use serial_environment_names
  use fem_update_names

  ! Geometry
  use maps_names
  use map_apply_names
  use fem_triangulation_names
  use mesh_triangulation_names

  !use fem_import_names
  !use partition_import
  use fem_triangulation_names
  use mesh_triangulation_names
  use fem_mesh_distribution_names
  use fem_mesh_gen_distribution_names

  use fem_element_import_names
  use fem_element_import_create_names
  use fem_conditions_names
  use fem_conditions_io_names
  use fem_materials_names
  use fem_materials_io_names
  use fem_graph_names
  use mesh_graph_names
  use geom2topo_names
  use renum_names
  use fem_mesh_names
  use fem_mesh_partition_base_names
  use fem_mesh_partition_distribution_names

  use fem_mesh_io_names
  use fem_mesh_gen_names
  use fem_mesh_refine_names
  !use fem_mesh_faces
  !use fem_mesh_lelpo
  use migratory_element_names
  use template_element_names
  use template_mesh_names

  ! Linear Algebra
  use fem_matrix_names
  use fem_block_graph_names
  use fem_precond_names
  use fem_vector_names
  use fem_vector_krylov_basis_names
  use fem_block_matrix_names
  use fem_block_precond_names
  use fem_block_vector_names
  use fem_block_vector_krylov_basis_names
  use fem_block_matrix_vector_names
  use solver_base_names
  use solver_names
  use abstract_solver_names
  use base_operand_names
  use base_operator_names
  
  ! Integration 
  use fem_space_names
  use fem_space_types_names
  use fem_space_faces_names
  use dof_handler_names
  use integration_names
  use quadrature_names
  use interpolation_names
  use femap_names
  use femap_interp_names
  use create_global_dof_info_names
  !use element_gather_tools
  use problem_names
  use integration_names

end module fem_names
