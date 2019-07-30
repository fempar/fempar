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
  use vtk_parameters_names, only: output_handler_vtk_format_cla_choices, &
                                  output_handler_vtk_format_default,     &
                                  output_handler_vtk_format_cla_help
  implicit none

  ! ---------------------------------------------------------------------------
  ! CLI keys
  ! ---------------------------------------------------------------------------

  character(len=*), parameter :: mesh_partitioner_num_levels_key        = 'MESH_PARTITIONER_NUM_LEVELS'
  character(len=*), parameter :: mesh_partitioner_num_parts_x_level_key = 'MESH_PARTITIONER_NUM_PARTS_X_LEVEL'
  character(len=*), parameter :: mesh_partitioner_strategy_key          = 'MESH_PARTITIONER_STRATEGY'
  character(len=*), parameter :: mesh_partitioner_dir_path_key          = 'MESH_PARTITIONER_DIR_PATH'
  character(len=*), parameter :: mesh_partitioner_prefix_key            = 'MESH_PARTITIONER_PREFIX'
  character(len=*), parameter :: mesh_partitioner_vtk_format_key        = 'MESH_PARTITIONER_VTK_FORMAT'
  
  character(len=*), parameter :: mesh_partitioner_metis_option_debug_key   = 'MESH_PARTITIONER_METIS_DEBUG'
  character(len=*), parameter :: mesh_partitioner_metis_option_ufactor_key = 'MESH_PARTITIONER_METIS_OPTION_UFACTOR'
  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_key = 'MESH_PARTITIONER_METIS_OPTION_MINCONN'
  character(len=*), parameter :: mesh_partitioner_metis_option_contig_key  = 'MESH_PARTITIONER_METIS_OPTION_CONFIG'
  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_key   = 'MESH_PARTITIONER_METIS_OPTION_CTYPE'
  character(len=*), parameter :: metis_option_iptype_key                   = 'MESH_PARTITIONER_METIS_OPTION_IPTYPE'

  ! ---------------------------------------------------------------------------
  ! CLI arguments
  ! ---------------------------------------------------------------------------

  character(len=*), parameter :: mesh_partitioner_num_levels_cla_name         = '--'//mesh_partitioner_num_levels_key
  character(len=*), parameter :: mesh_partitioner_num_parts_x_level_cla_name  = '--'//mesh_partitioner_num_parts_x_level_key
  character(len=*), parameter :: mesh_partitioner_strategy_cla_name           = '--'//mesh_partitioner_strategy_key
  character(len=*), parameter :: mesh_partitioner_dir_path_cla_name           = '--'//mesh_partitioner_dir_path_key
  character(len=*), parameter :: mesh_partitioner_prefix_cla_name             = '--'//mesh_partitioner_prefix_key
  character(len=*), parameter :: mesh_partitioner_vtk_format_cla_name         = '--'//mesh_partitioner_vtk_format_key

  character(len=*), parameter :: mesh_partitioner_metis_option_debug_cla_name   = '--'//mesh_partitioner_metis_option_debug_key
  character(len=*), parameter :: mesh_partitioner_metis_option_ufactor_cla_name = '--'//mesh_partitioner_metis_option_ufactor_key
  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_cla_name = '--'//mesh_partitioner_metis_option_minconn_key
  character(len=*), parameter :: mesh_partitioner_metis_option_contig_cla_name  = '--'//mesh_partitioner_metis_option_contig_key
  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_cla_name   = '--'//mesh_partitioner_metis_option_ctype_key
  character(len=*), parameter :: mesh_partitioner_metis_option_iptype_cla_name  = '--'//metis_option_iptype_key

  character(len=*), parameter :: metis_part_kway      = 'METIS_PART_KWAY'
  character(len=*), parameter :: metis_part_recursive = 'METIS_PART_RECURSIVE'
  character(len=*), parameter :: part_strip           = 'PART_STRIP'
  character(len=*), parameter :: part_rcm_strip       = 'PART_RCM_STRIP'

  ! ---------------------------------------------------------------------------
  ! CLI choices
  ! ---------------------------------------------------------------------------

  character(len=*), parameter :: mesh_partitioner_strategy_cla_choices = metis_part_kway      // ',' // &
                                                                         metis_part_recursive // ',' // &
                                                                         part_strip           // ',' // &
                                                                         part_rcm_strip 

  character(len=*), parameter :: mesh_partitioner_vtk_format_cla_choices = output_handler_vtk_format_cla_choices

  character(len=*), parameter :: mesh_partitioner_metis_option_iptype_cla_choices = '0,1,2,3,4'
  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_cla_choices   = '0,1'
  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_cla_choices = '0,1'
  character(len=*), parameter :: mesh_partitioner_metis_option_contig_cla_choices  = '0,1'

  ! ---------------------------------------------------------------------------
  ! Default values
  ! ---------------------------------------------------------------------------
  
  character(len=*),          parameter :: mesh_partitioner_default_strat = metis_part_kway
  character(len=*),          parameter :: mesh_partitioner_default_dir_path_output = '.'  
  character(len=*),          parameter :: mesh_partitioner_default_prefix_output   = 'mesh_partition'
  character(len=*),          parameter :: mesh_partitioner_default_vtk_format = output_handler_vtk_format_default
  character(len=*),          parameter :: default_dir_path                              ='.'
  character(len=*),          parameter :: default_prefix                                ='A'
  integer(ip),               parameter :: default_mesh_partitioner_num_levels           = 1
  integer(ip), dimension(*), parameter :: default_mesh_partitioner_num_parts_x_level    = [1]
  
  ! Use METIS defaults (i.e., == -1) 30 for part_kway, and 1 for part_recursive
  integer(ip), parameter :: mesh_partitioner_default_metis_option_ufactor = -1
  
  ! (Try to) Minimize maximum degree of subdomain graph
  integer(ip), parameter :: mesh_partitioner_default_metis_option_minconn = 0
  
  ! (Try to) Produce partitions that are contiguous
  integer(ip), parameter :: mesh_partitioner_default_metis_option_contig = 1

  ! Debug option
  integer(ip), parameter :: mesh_partitioner_default_metis_option_debug = 0

  ! Random matching
  integer(ip), parameter :: mesh_partitioner_default_metis_option_ctype = METIS_CTYPE_RM
  
  ! Grow bisection greedy
  integer(ip), parameter :: mesh_partitioner_default_metis_option_iptype = METIS_IPTYPE_GROW

  ! ---------------------------------------------------------------------------
  ! CLI help messages
  ! ---------------------------------------------------------------------------

  character(len=*), parameter :: dir_path_cla_help                              = "Directory of the source files to read from GiD mesh" 
  character(len=*), parameter :: prefix_cla_help                                = "Name of the GiD files"
  character(len=*), parameter :: mesh_partitioner_dir_path_cla_help             = "The relative or full file system path to the folder where the data files of the" // BRK_LINE // &
                                                                                  "mesh partitioner are generated"
  character(len=*), parameter :: mesh_partitioner_prefix_cla_help               = "Token string which is used as a prefix to compose the names of the data files"   // BRK_LINE // &
                                                                                  "generated by the mesh partitioner as prefix.*"
  character(len=*), parameter :: mesh_partitioner_num_levels_cla_help           = "Number of levels in the triangulation hierarchy required for the "               // BRK_LINE // &
                                                                                  "MLBDDC preconditioner"
  character(len=*), parameter :: mesh_partitioner_num_parts_x_level_cla_help    = "Number of parts per level"
  character(len=*), parameter :: mesh_partitioner_vtk_format_cla_help           = output_handler_vtk_format_cla_help
  character(len=*), parameter :: mesh_partitioner_metis_option_debug_cla_help   = "METIS_OPTION_DBGLVL (see METIS users' manual for additional details)"
  character(len=*), parameter :: mesh_partitioner_metis_option_ufactor_cla_help = "METIS_OPTION_UFACTOR, (see METIS users' manual for additional details)."         // BRK_LINE // &
                                                                                  "The special value -1 sets the METIS defaults for part_kway and part_recursive"

  character(len=*), parameter :: mesh_partitioner_strategy_cla_help = &
                       'Mesh partitioning strategy/algorithm'                                   // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // metis_part_kway // &
                            ": Partitions the mesh (dual graph) into $k$ parts using"           // BRK_LINE // &
                            "  multilevel k-way partitioning as implemented in METIS"           // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // metis_part_recursive // &
                            ": Partitions the mesh (dual graph) into $k$ parts"                 // BRK_LINE // &
                            "  using multilevel recursive bisection as implemented in METIS"    // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // part_strip // &
                            ": Uses a trivial algorithm to partition the mesh (dual graph)"     // BRK_LINE // &
                            "  cells into equally-sized blocks, where each block is composed"   // BRK_LINE // &
                            "  by cells with consecutive cell identifiers (the mesh partition"  // BRK_LINE // &
                            "  thus depends on the numbering of mesh cells)"                    // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // part_rcm_strip // &
                            ": As "//part_strip//", but the mesh cell identifiers (the dual"    // BRK_LINE // &
                            "  graph vertex identifiers) are re-ordered using the Reverse"      // BRK_LINE // &
                            "  Cuthill-McKee algorithm in advance"

  character(len=*), parameter :: mesh_partitioner_metis_option_iptype_cla_help = &
                       "METIS_OPTION_IPTYPE (see METIS users' manual for additional details)"                                         // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "0: METIS_IPTYPE_GROW (Grows a bisection using a greedy strategy)"                 // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "1: METIS_IPTYPE_RANDOM (Computes a bisection at random followed by a refinement)" // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // "2: METIS_IPTYPE_EDGE (Derives a separator from an edge cut)"                      // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // "3: METIS_IPTYPE_NODE (Grow a bisection using a greedy node-based strategy)"       // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // "4: METIS_IPTYPE_METISRB"

  character(len=*), parameter :: mesh_partitioner_metis_option_ctype_cla_help = &
                       "METIS_OPTION_CTYPE (see METIS users' manual for additional details)" // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "0: METIS_CTYPE_RM (Random matching)"                       // BRK_LINE // &
                       BULLET_FLAP_HELP_MESSAGE // "1: METIS_CTYPE_SHEM (Sorted heavy-edge matching)" 

  character(len=*), parameter :: mesh_partitioner_metis_option_minconn_cla_help = &
                       "METIS_OPTION_MINCONN (see METIS users' manual for additional details)"                  // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "0: Does not explicitly minimize the maximum connectivity"   // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "1: Explicitly minimize the maximum connectivity"


  character(len=*), parameter :: mesh_partitioner_metis_option_contig_cla_help = &
                       "METIS_OPTION_CONTIG (see METIS users' manual for additional details)"  // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "0: Does not force contiguous partitions"   // BRK_LINE // & 
                       BULLET_FLAP_HELP_MESSAGE // "1: Forces contiguous partitions"

end module mesh_partitioner_parameters_names
