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
module mesh_partitioner_parameters_names
  use types_names
  use metis_names
  implicit none

  character(len=*), parameter :: mesh_partitioner_num_levels_key        = 'MESH_PARTITIONER_NUM_LEVELS'
  character(len=*), parameter :: mesh_partitioner_num_parts_x_level_key = 'MESH_PARTITIONER_NUM_PARTS_X_LEVEL'
  character(len=*), parameter :: mesh_partitioner_strategy_key          = 'MESH_PARTITIONER_STRATEGY'
  character(len=*), parameter :: mesh_partitioner_dir_path_key          = 'MESH_PARTITIONER_DIR_PATH'
  character(len=*), parameter :: mesh_partitioner_prefix_key            = 'MESH_PARTITIONER_PREFIX'
  
  character(len=*), parameter :: mesh_partitioner_metis_option_debug_key   = 'MESH_PARTITIONER_METIS_DEBUG'
  character(len=*), parameter :: mesh_partitioner_metis_option_ufactor_key = 'MESH_PARTITIONER_METIS_OPTION_UFACTOR'
  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_key = 'MESH_PARTITIONER_METIS_OPTION_MINCONN'
  character(len=*), parameter :: mesh_partitioner_metis_option_contig_key  = 'MESH_PARTITIONER_METIS_OPTION_CONFIG'
  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_key   = 'MESH_PARTITIONER_METIS_OPTION_CTYPE'
  character(len=*), parameter :: metis_option_iptype_key                   = 'MESH_PARTITIONER_METIS_OPTION_IPTYPE'

  character(len=*), parameter :: mesh_partitioner_num_levels_cla_name         = '--'//mesh_partitioner_num_levels_key
  character(len=*), parameter :: mesh_partitioner_num_parts_x_level_cla_name  = '--'//mesh_partitioner_num_parts_x_level_key
  character(len=*), parameter :: mesh_partitioner_strategy_cla_name           = '--'//mesh_partitioner_strategy_key
  character(len=*), parameter :: mesh_partitioner_dir_path_cla_name           = '--'//mesh_partitioner_dir_path_key
  character(len=*), parameter :: mesh_partitioner_prefix_cla_name             = '--'//mesh_partitioner_prefix_key

  character(len=*), parameter :: mesh_partitioner_metis_option_debug_cla_name   = '--'//mesh_partitioner_metis_option_debug_key
  character(len=*), parameter :: mesh_partitioner_metis_option_ufactor_cla_name = '--'//mesh_partitioner_metis_option_ufactor_key
  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_cla_name = '--'//mesh_partitioner_metis_option_minconn_key
  character(len=*), parameter :: mesh_partitioner_metis_option_contig_cla_name  = '--'//mesh_partitioner_metis_option_contig_key
  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_cla_name   = '--'//mesh_partitioner_metis_option_ctype_key
  character(len=*), parameter :: mesh_partitioner_metis_option_iptype_cla_name  = '--'//metis_option_iptype_key

  character(len=*), parameter :: metis_part_kway      = 'metis_part_kway'
  character(len=*), parameter :: metis_part_recursive = 'metis_part_recursive'
  character(len=*), parameter :: part_strip           = 'part_strip'
  character(len=*), parameter :: part_rcm_strip       = 'part_rcm_strip'
  
  character(len=*), parameter :: mesh_partitioner_default_dir_path_output = '.'  
  character(len=*), parameter :: mesh_partitioner_default_prefix_output   = 'mesh_partition'
  
  character(len=*), parameter :: mesh_partitioner_default_strat = metis_part_kway
  
  ! Use METIS defaults (i.e., == -1) 30 for part_kway, and 1 for part_recursive
  integer(ip), parameter :: mesh_partitioner_default_metis_option_ufactor = -1
  
  ! (Try to) Minimize maximum degree of subdomain graph
  integer(ip), parameter :: mesh_partitioner_default_metis_option_minconn = 1
  
  ! (Try to) Produce partitions that are contiguous
  integer(ip), parameter :: mesh_partitioner_default_metis_option_contig = 1
  
  ! Random matching
  integer(ip), parameter :: mesh_partitioner_default_metis_option_ctype = METIS_CTYPE_RM
  
  ! Grow bisection greedy
  integer(ip), parameter :: mesh_partitioner_default_metis_option_iptype = METIS_IPTYPE_GROW
  
  integer(ip), parameter :: mesh_partitioner_default_metis_option_debug = 0

end module mesh_partitioner_parameters_names
