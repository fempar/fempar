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
module fem

  ! Tools
  use types
  use memor
  use array_names
  use stdio
  use hash_table_names
  use postpro_names

  ! Geometry
  use maps_names
  use map_apply
  use fem_partition_names

  !use fem_import_names
  !use partition_import
  use fem_mesh_distribution_names

  use fem_element_import_names
  use partition_element_import
  use fem_conditions_names
  use fem_conditions_io
  use fem_materials_names
  use fem_materials_io
  use fem_graph_names
  use fem_graph_partition
  use mesh_graph
  use mesh_graph_partition
  use geom2topo
  use renum_names
  use fem_mesh_names
  use fem_mesh_partition_base
  use fem_mesh_partition
  use fem_mesh_partition_distribution

  use fem_mesh_io
  use fem_mesh_gen
  use fem_mesh_gen_partition
  use fem_mesh_refine
  !use fem_mesh_faces
  !use fem_mesh_lelpo
  use fem_mesh_dual
  use migratory_element_names
  use template_element_names
  use template_mesh_names

  ! Linear Algebra
  use fem_matrix_names
  use fem_precond_names
  use fem_vector_names
  use fem_vector_krylov_basis_names
  use fem_block_matrix_names
  use fem_block_precond_names
  use fem_block_vector_names
  use fem_block_vector_krylov_basis_names
  use fem_block_matrix_vector
  use solver_base
  use solver
  use abstract_solver
  
  ! Integration 
  use fem_space_names
  use fem_space_types
  use fem_space_faces
  use dof_handler_names
  use integration_names
  use quadrature_names
  use interpolation_names
  use femap_names
  use femap_interp
  use create_global_dof_info_names
  !use element_gather_tools

end module fem
