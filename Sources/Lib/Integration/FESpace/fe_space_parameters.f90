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
module fe_space_parameters_names
  use types_names
  use reference_fe_names
  implicit none
  
  ! These three parameter constants are thought be used as FPL keys. The corresponding pairs 
  ! <XXXkey,.true.|.false.> in FPL let the user to control whether or not coarse vertex, edge, or 
  ! face DoFs are to be included into coarse_fe_space_t. These parameters might be used during 
  ! coarse_fe_space_t set-up, as well as by the deferred TBP methods corresponding to 
  ! class(coarse_fe_handler_t).
  character(len=*), parameter :: coarse_space_use_vertices_key = 'FES_COARSE_SPACE_USE_VERTICES'
  character(len=*), parameter :: coarse_space_use_edges_key    = 'FES_COARSE_SPACE_USE_EDGES'
  character(len=*), parameter :: coarse_space_use_faces_key    = 'FES_COARSE_SPACE_USE_FACES'
  
  ! Keys being used by the FE space constructor that relies on the parameter handler
  character(len=*), parameter, public :: fes_num_fields_key = 'FES_NUM_FIELDS'
  character(len=*), parameter, public :: fes_num_ref_fes_key = 'FES_NUM_REF_FES'
  character(len=*), parameter, public :: fes_field_types_key = 'FES_FIELD_TYPES'
  character(len=*), parameter, public :: fe_space_field_ids_key = 'FE_SPACE_FIELD_IDS'
  character(len=*), parameter, public :: fes_field_blocks_key = 'FES_FIELD_BLOCKS'
  character(len=*), parameter, public :: fes_same_ref_fes_all_cells_key = 'FES_SAME_REFS_ALL_CELLS'
      
  character(len=*), parameter, public :: fes_ref_fe_conformities_key = 'FES_REF_FE_CONFORMITIES'
  character(len=*), parameter, public :: fes_ref_fe_continuities_key = 'FES_REF_FE_CONTINUITIES'
  character(len=*), parameter, public :: fes_ref_fe_orders_key = 'FES_REF_FE_ORDERS'
  character(len=*), parameter, public :: fes_ref_fe_types_key = 'FES_REF_FE_TYPES'
  character(len=*), parameter, public :: fes_set_ids_ref_fes_key = 'FES_SET_IDS_REF_FES'

  ! These three parameter constants are thought be used as CLA names. 

  character(len=*), parameter, public :: coarse_space_use_vertices_cla_name          = '--'//coarse_space_use_vertices_key
  character(len=*), parameter, public :: coarse_space_use_edges_cla_name             = '--'//coarse_space_use_edges_key
  character(len=*), parameter, public :: coarse_space_use_faces_cla_name             = '--'//coarse_space_use_faces_key
  character(len=*), parameter, public :: fes_num_fields_cla_name             = '--'//fes_num_fields_key
  character(len=*), parameter, public :: fes_num_ref_fes_cla_name            = '--'//fes_num_ref_fes_key
  character(len=*), parameter, public :: fes_field_types_cla_name            = '--'//fes_field_types_key
  character(len=*), parameter, public :: fe_space_field_ids_cla_name         = '--'//fe_space_field_ids_key
  character(len=*), parameter, public :: fes_field_blocks_cla_name           = '--'//fes_field_blocks_key
  character(len=*), parameter, public :: fes_same_ref_fes_all_cells_cla_name = '--'//fes_same_ref_fes_all_cells_key
  character(len=*), parameter, public :: fes_ref_fe_conformities_cla_name    = '--'//fes_ref_fe_conformities_key
  character(len=*), parameter, public :: fes_ref_fe_continuities_cla_name    = '--'//fes_ref_fe_continuities_key
  character(len=*), parameter, public :: fes_ref_fe_orders_cla_name          = '--'//fes_ref_fe_orders_key
  character(len=*), parameter, public :: fes_ref_fe_types_cla_name           = '--'//fes_ref_fe_types_key
  character(len=*), parameter, public :: fes_set_ids_ref_fes_cla_name        = '--'//fes_set_ids_ref_fes_key

  ! CLA defaults

  logical,                   parameter, public :: default_coarse_space_use_vertices  = .true.
  logical,                   parameter, public :: default_coarse_space_use_edges     = .true.
  logical,                   parameter, public :: default_coarse_space_use_faces     = .true.
  integer(ip),               parameter, public :: default_fes_num_fields             = 1
  integer(ip),               parameter, public :: default_fes_num_ref_fes            = 1
  integer(ip), dimension(*), parameter, public :: default_fes_set_ids_ref_fes        = [1]
  integer(ip), dimension(*), parameter, public :: default_fes_ref_fe_orders          = [1]
  logical,     dimension(*), parameter, public :: default_fes_ref_fe_conformities    = [ .true. ]
  logical,     dimension(*), parameter, public :: default_fes_ref_fe_continuities    = [ .true. ]
  integer(ip), dimension(*), parameter, public :: default_fes_field_blocks           = [1]
  logical,                   parameter, public :: default_fes_same_ref_fes_all_cells = .true.

  ! CLA choices
  
  character(len=*), parameter, public :: fes_field_types_cla_choices = field_type_scalar // "," // &
                                                                       field_type_vector // "," // &
                                                                       field_type_tensor

  ! CLA help  

  character(len=*), parameter, public :: coarse_space_use_vertices_cla_help  = 'Coarse-space shape functions on vertices'
  character(len=*), parameter, public :: coarse_space_use_edges_cla_help     = 'Coarse-space shape functions on edges'
  character(len=*), parameter, public :: coarse_space_use_faces_cla_help     = 'Coarse-space shape functions on faces'
  character(len=*), parameter, public :: fes_num_fields_cla_help             = 'Finite element space number of fields'
  character(len=*), parameter, public :: fes_num_ref_fes_cla_help            = 'Finite element space number of fields'
  character(len=*), parameter, public :: fes_set_ids_ref_fes_cla_help        = 'Set IDs to reference FEs for every field'
  character(len=*), parameter, public :: fes_ref_fe_orders_cla_help          = 'Reference finite element orders'
  character(len=*), parameter, public :: fes_ref_fe_conformities_cla_help    = 'Finite element space conformities'
  character(len=*), parameter, public :: fes_ref_fe_continuities_cla_help    = 'Finite element space continuities'
  character(len=*), parameter, public :: fes_field_blocks_cla_help           = 'Finite element space map field_id to block_id'
  character(len=*), parameter, public :: fes_same_ref_fes_all_cells_cla_help = 'Finite element space fixed reference fe logical'
                                                                       
  character(len=*), parameter, public :: fes_field_types_cla_help    = "Finite element space field types" // BRK_LINE // & 
                                                                       BULLET_FLAP_HELP_MESSAGE // field_type_scalar // ": Scalar-Valued Reference FE" // BRK_LINE // & 
                                                                       BULLET_FLAP_HELP_MESSAGE // field_type_vector // ": Vector-Valued Reference FE" // BRK_LINE // &
                                                                       BULLET_FLAP_HELP_MESSAGE // field_type_tensor // ": Tensor-Valued Reference FE"
                                                                       
  character(len=*), parameter, public :: fes_ref_fe_types_cla_choices = fe_type_lagrangian // "," // &
                                                                        fe_type_lagrangian_gp // "," // &
                                                                        fe_type_raviart_thomas // "," // &
                                                                        fe_type_nedelec // "," // &
                                                                        fe_type_void
                                                                       
  character(len=*), parameter, public :: fes_ref_fe_types_cla_help    = "Reference finite element types" // BRK_LINE // & 
                                                                       BULLET_FLAP_HELP_MESSAGE // fe_type_lagrangian // &
                                                                        ": Grad-conforming Lagrangian Reference FE of arbitrary degree $k$ on top of n-cubes (i.e., $Q_k$)" // BRK_LINE // &
                                                                        "  and n-simplices (i.e., $P_k$), for the discretization of either scalar-valued, vector-valued or tensor-valued fields." // BRK_LINE // & 
                                                                       BULLET_FLAP_HELP_MESSAGE // fe_type_lagrangian_gp  // &
                                                                       ": " // fe_type_lagrangian_gp // BRK_LINE // &
                                                                       BULLET_FLAP_HELP_MESSAGE // fe_type_raviart_thomas // &
                                                                       ": The vector-valued div-conforming Raviart-Thomas Reference FE of arbitrary degree $k$ on top of " // BRK_LINE // &
                                                                       "  n-cubes (n-simplices not supported yet) suitable for the mixed Laplacian problem and some fluid flow problems" // BRK_LINE // &
                                                                       BULLET_FLAP_HELP_MESSAGE // fe_type_nedelec // & 
                                                                       ": The vector-valued curl-conforming Nèdèlec Reference FE (of first kind) of arbitrary degree $k$ on " // BRK_LINE // &
                                                                       "  top of n-cubes and n-simplices suitable for electromagnetic problems" // BRK_LINE // &
                                                                       BULLET_FLAP_HELP_MESSAGE // fe_type_void // & 
                                                                       ": A software artifact that represents a Reference FE with no DOFs at all, neither at the cell interiors," // BRK_LINE // &
                                                                       "  nor at their boundary n-faces."
  
  ! These parameter constants are used in order to generate a unique (non-consecutive) 
  ! but consistent across MPI tasks global ID (integer(igp)) of a given DoF.
  ! See type(par_fe_space_t)%generate_non_consecutive_dof_ggid()
  integer(ip), parameter :: cell_ggid_shift              = 44
  integer(ip), parameter :: dofs_x_reference_fe_shift = 14
  integer(ip), parameter :: num_fields_shift         = igp*8-(cell_ggid_shift+dofs_x_reference_fe_shift)-1


  ! Node type 
  integer(ip), parameter :: num_node_types     = 2
  integer(ip), parameter :: edge_boundary_node = 1
  integer(ip), parameter :: interior_node      = 2

  ! Fine edge direction related to coarse edge direction 
  integer(ip), parameter :: opposite_to_coarse_edge = 0
  integer(ip), parameter :: same_as_coarse_edge     = 1  

  ! Coarse Edges enforced Continuity algorithm 
  character(len=*), parameter :: bddc_edge_continuity_algorithm_key        = 'BDDC_EDGE_CONTINUYRY_ALGORITHM'
  character(len=*), parameter :: bddc_edge_continuity_algorithm_cla_name   = '--'//bddc_edge_continuity_algorithm_key
  character(len=*), parameter :: tangential_average                        = 'tangential_average'
  character(len=*), parameter :: tangential_average_and_first_order_moment = 'tangential_average_and_first_order_moment'
  character(len=*), parameter :: all_dofs_in_coarse_edges                  = 'all_dofs_in_coarse_edges'

  ! BDDC Scaling function 
  character(len=*), parameter :: bddc_scaling_function_case_key      = 'BDDC_SCALING_FUNCTION_CASE'
  character(len=*), parameter :: bddc_scaling_function_case_cla_name = '--'//bddc_scaling_function_case_key
  character(len=*), parameter :: bddc_scaling_function_case_cla_help = 'Scaling type for the BDDC weighting operator'

  character(len=*), parameter :: average_mass_coeff_key         = 'average_mass_coeff' 
  character(len=*), parameter :: average_curl_curl_coeff_key    = 'average_curl_curl_coeff'
  character(len=*), parameter :: cardinality           = 'cardinality'
  character(len=*), parameter :: curl_curl_coeff       = 'curl_curl_coeff' 
  character(len=*), parameter :: mass_coeff            = 'mass_coeff' 
  character(len=*), parameter :: stiffness             = 'stiffness'
  character(len=*), parameter :: weighted_coefficients = 'weighted_coefficients'


end module fe_space_parameters_names
