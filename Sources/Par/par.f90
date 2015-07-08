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
module par_names
  ! Tools
  use par_context_names
  use psb_penv_mod_names
  use par_sparse_global_collectives_names
  use par_element_exchange_names
  use par_timer_names
  use par_io_names
  use par_environment_names
  use par_update_names

  ! Geometry
  use par_mesh_names
  use par_triangulation_names
  use par_mesh_to_triangulation_names
  use par_conditions_names
  use par_generate_uniform_triangulation_names
  use par_uniform_refinement_names

  ! Linear algebra
  use par_vector_names
  use par_matrix_names
  use par_graph_names
  use par_block_matrix_names
  use par_block_vector_names
  use par_block_graph_names
  use block_dof_distribution_names
  use par_dd_base_names
  use par_preconditioner_dd_diagonal_names
  use par_preconditioner_dd_mlevel_bddc_names
  use par_preconditioner_dd_identity_names

  ! Integration
  use par_fe_space_names
  use par_create_global_dof_info_names
    
end module par_names
