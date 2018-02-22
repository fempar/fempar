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
module fe_space_names
  ! Serial modules
  use types_names
  use stdio_names
  use list_types_names
  use memor_names
  use sort_names
  use allocatable_array_names
  use hash_table_names
  use std_vector_names
  use FPL

  use environment_names
  use triangulation_names
  use conditions_names
  
  use reference_fe_names
  use field_names
  use function_names
  
  use block_layout_names
  use vector_names
  use array_names
  use assembler_names
  use base_sparse_matrix_names
  use sparse_matrix_names
  use block_sparse_matrix_names
  use serial_scalar_array_names
  use serial_block_array_names
  use sparse_assembler_names
  use block_sparse_assembler_names
  
  use direct_solver_names
  use direct_solver_parameters_names
  use iterative_linear_solver_names
  use iterative_linear_solver_parameters_names
  
  use piecewise_cell_map_names
  
  ! Adaptivity
  use std_vector_real_rp_names
  use std_vector_integer_ip_names
  use std_vector_integer_igp_names
  use p4est_triangulation_names

  ! Parallel modules
  use environment_names
  use dof_import_names
  use cell_import_names
  use par_sparse_matrix_names
  use par_scalar_array_names
  use par_sparse_assembler_names
  
  implicit none
# include "debug.i90"
  private
  
  ! Towards having the final hierarchy of FE spaces, I moved  to a root superclass
  ! those member variables which are in common by type(serial/par_fe_space_t) and 
  ! type(coarse_fe_space_t) and in turn which are required by either type(mlbddc_t) (L1) 
  ! or type(mlbddc_coarse_t) (L2-LN). This is not (by far) the best that can be done for 
  ! minimizing the amount of code replication within the hierarchy of FE spaces. In particular, 
  ! more member variables of sub-classes should we moved here, as well as the corresponding 
  ! TBPs in charge of handling those.
  type :: base_fe_space_t
    private
    integer(ip)                                 :: num_fields
    integer(ip) , allocatable                   :: fe_space_type_x_field(:)
    integer(ip) , allocatable                   :: num_dofs_x_field(:)
    type(dof_import_t)            , allocatable :: blocks_dof_import(:)
    
    ! Pointer to data structure which is in charge of coarse DoF handling.
    ! It will be a nullified pointer on L1 tasks, and associated via target 
    ! allocation in the case of L2-Ln tasks.
    type(coarse_fe_space_t)       , pointer :: coarse_fe_space => NULL()
  contains 
    procedure, non_overridable, private        :: free_blocks_dof_import                       => base_fe_space_free_blocks_dof_import
    
    procedure, non_overridable                 :: get_num_fields                               => base_fe_space_get_num_fields
    procedure, non_overridable                 :: set_num_fields                               => base_fe_space_set_num_fields
    procedure, non_overridable                 :: get_fe_space_type                               => base_fe_space_get_fe_space_type
    procedure                                  :: get_num_blocks                               => base_fe_space_get_num_blocks
    procedure                                  :: get_field_blocks                                => base_fe_space_get_field_blocks
    procedure                                  :: get_field_coupling                              => base_fe_space_get_field_coupling
    procedure, non_overridable                 :: get_total_num_dofs                           => base_fe_space_get_total_num_dofs
    procedure, non_overridable                 :: get_field_num_dofs                           => base_fe_space_get_field_num_dofs
    procedure, non_overridable                 :: set_field_num_dofs                           => base_fe_space_set_field_num_dofs
    procedure                                  :: get_block_num_dofs                           => base_fe_space_get_block_num_dofs
    procedure, non_overridable                 :: get_total_num_interior_dofs                  => base_fe_space_get_total_num_interior_dofs
    procedure, non_overridable                 :: get_total_num_interface_dofs                 => base_fe_space_get_total_num_interface_dofs
    procedure, non_overridable                 :: get_block_num_interior_dofs                  => base_fe_space_get_block_num_interior_dofs
    procedure, non_overridable                 :: get_block_num_interface_dofs                 => base_fe_space_get_block_num_interface_dofs
    procedure, non_overridable                 :: get_block_dof_import                            => base_fe_space_get_block_dof_import
    procedure, non_overridable                 :: get_coarse_fe_space                             => base_fe_space_get_coarse_fe_space
    
    ! These three TBPs are (currently) though to be overriden by subclasses.
    ! In the future, once we harmonize type(coarse_fe_space_t) and type(serial/par_fe_space_t),
    ! the code of these TBPs should be developed at this level, being valid for all subclasses.
    procedure                                  :: get_environment                                 => base_fe_space_get_environment
    procedure                                  :: get_total_num_coarse_dofs                    => base_fe_space_get_total_num_coarse_dofs
    procedure                                  :: get_block_num_coarse_dofs                    => base_fe_space_get_block_num_coarse_dofs
  end type base_fe_space_t
  
  public :: base_fe_space_t

  type :: base_fe_cell_iterator_t
    private
    class(cell_iterator_t), allocatable :: cell
  contains
    ! Methods exploiting ("inherited from") cell_t common to all descendants
    procedure                            :: next                    => base_fe_cell_iterator_next
    procedure                            :: first                   => base_fe_cell_iterator_first
    procedure                            :: set_gid                 => base_fe_cell_iterator_set_gid
    procedure, non_overridable           :: has_finished            => base_fe_cell_iterator_has_finished
    procedure, non_overridable           :: get_reference_fe_geo    => base_fe_cell_iterator_get_reference_fe_geo
    procedure, non_overridable           :: get_reference_fe_geo_id => base_fe_cell_iterator_get_reference_fe_geo_id
    procedure, non_overridable           :: get_num_nodes           => base_fe_cell_iterator_get_num_nodes
    procedure, non_overridable           :: get_nodes_coordinates   => base_fe_cell_iterator_get_nodes_coordinates
    procedure, non_overridable           :: get_gid                 => base_fe_cell_iterator_get_gid
    procedure, non_overridable           :: get_ggid                => base_fe_cell_iterator_get_ggid
    procedure, non_overridable           :: get_my_part             => base_fe_cell_iterator_get_mypart
    procedure, non_overridable           :: get_my_subpart          => base_fe_cell_iterator_get_mysubpart
    procedure, non_overridable           :: get_my_subpart_lid      => base_fe_cell_iterator_get_mysubpart_lid
    procedure, non_overridable           :: get_set_id              => base_fe_cell_iterator_get_set_id
    procedure, non_overridable           :: get_num_vefs            => base_fe_cell_iterator_get_num_vefs
    procedure, non_overridable           :: get_vef_gid             => base_fe_cell_iterator_get_vef_gid
    procedure, non_overridable           :: get_vefs_gid            => base_fe_cell_iterator_get_vefs_gid
    procedure, non_overridable           :: get_vef_ggid            => base_fe_cell_iterator_get_vef_ggid
    procedure, non_overridable           :: get_vef_lid_from_gid    => base_fe_cell_iterator_get_vef_lid_from_gid
    procedure, non_overridable           :: get_vef_lid_from_ggid   => base_fe_cell_iterator_get_vef_lid_from_ggid
    procedure, non_overridable           :: is_local                => base_fe_cell_iterator_is_local
    procedure, non_overridable           :: is_ghost                => base_fe_cell_iterator_is_ghost
    procedure, non_overridable, private  :: base_fe_cell_iterator_get_vef
    generic                              :: get_vef                 => base_fe_cell_iterator_get_vef
    procedure, non_overridable           :: get_triangulation       => base_fe_cell_iterator_get_triangulation
    procedure, non_overridable           :: get_permutation_index   => base_fe_cell_iterator_get_permutation_index
    
    procedure                            :: update_sub_triangulation    => base_fe_cell_iterator_update_sub_triangulation
    procedure                            :: get_mc_case                 => base_fe_cell_iterator_get_mc_case
    procedure                            :: get_num_subcells            => base_fe_cell_iterator_get_num_subcells
    procedure                            :: get_num_subcell_nodes       => base_fe_cell_iterator_get_num_subcell_nodes
    procedure                            :: get_phys_coords_of_subcell  => base_fe_cell_iterator_get_phys_coords_of_subcell
    procedure                            :: get_ref_coords_of_subcell   => base_fe_cell_iterator_get_ref_coords_of_subcell
    procedure                            :: get_num_subfacets           => base_fe_cell_iterator_get_num_subfacets
    procedure                            :: get_num_subfacet_nodes      => base_fe_cell_iterator_get_num_subfacet_nodes
    procedure                            :: get_phys_coords_of_subfacet => base_fe_cell_iterator_get_phys_coords_of_subfacet
    procedure                            :: get_ref_coords_of_subfacet  => base_fe_cell_iterator_get_ref_coords_of_subfacet
    procedure                            :: is_cut                      => base_fe_cell_iterator_is_cut
    procedure                            :: is_interior                 => base_fe_cell_iterator_is_interior
    procedure                            :: is_exterior                 => base_fe_cell_iterator_is_exterior
    procedure                            :: is_interior_subcell         => base_fe_cell_iterator_is_interior_subcell
    procedure                            :: is_exterior_subcell         => base_fe_cell_iterator_is_exterior_subcell
  end type base_fe_cell_iterator_t  
  
  type, extends(base_fe_cell_iterator_t) :: fe_cell_iterator_t
    private
    class(serial_fe_space_t) , pointer     :: fe_space => NULL()
    ! Scratch data to support FE integration
    integer(ip)                        , allocatable :: num_cell_dofs_x_field(:)
    type(i1p_t)                        , allocatable :: fe_dofs(:)
    type(cell_map_t)                   , pointer     :: cell_map => NULL()
    type(p_cell_integrator_t)          , allocatable :: cell_integrators(:)
    type(std_vector_integer_ip_t)      , allocatable :: extended_fe_dofs(:)
    type(allocatable_array_ip1_t)      , allocatable :: gid_to_lid_map(:)
    type(allocatable_array_rp2_t)                    :: extended_elmat
    type(allocatable_array_rp1_t)                    :: extended_elvec
  contains
  
    procedure                           :: create                                     => fe_cell_iterator_create
    procedure                           :: create_cell                                => fe_cell_iterator_create_cell
    procedure                           :: free_cell                                  => fe_cell_iterator_free_cell
    procedure                           :: free                                       => fe_cell_iterator_free
    final                               :: fe_cell_iterator_free_final
    procedure                           :: set_fe_space                               => fe_cell_iterator_set_fe_space
    procedure                           :: nullify_fe_space                           => fe_cell_iterator_nullify_fe_space
    procedure, non_overridable          :: allocate_scratch_data                      => fe_cell_iterator_allocate_scratch_data
    procedure, non_overridable          :: allocate_block_based_scratch_data          => fe_cell_iterator_allocate_block_based_scratch_data
    procedure, non_overridable          :: free_scratch_data                          => fe_cell_iterator_free_scratch_data
    procedure, non_overridable, private :: count_own_dofs_cell                        => fe_cell_iterator_count_own_dofs_cell
    procedure, non_overridable, private :: count_own_dofs_vef                         => fe_cell_iterator_count_own_dofs_vef
    procedure, non_overridable          :: generate_own_dofs_cell                     => fe_cell_iterator_generate_own_dofs_cell
    procedure, non_overridable          :: generate_own_dofs_vef                      => fe_cell_iterator_generate_own_dofs_vef
    procedure, non_overridable          :: generate_own_dofs_vef_component_wise       => fe_cell_iterator_generate_own_dofs_vef_component_wise
    procedure, non_overridable          :: fetch_own_dofs_vef_from_source_fe          => fe_cell_iterator_fetch_own_dofs_vef_from_source_fe
    procedure, non_overridable, private :: generate_dofs_facet_integration_coupling   => fe_cell_iterator_generate_dofs_facet_integration_coupling
    procedure, non_overridable, private :: renum_dofs_block                           => fe_cell_iterator_renum_dofs_block
    procedure, non_overridable, private :: renum_dofs_field                           => fe_cell_iterator_renum_dofs_field
    procedure, non_overridable, private :: update_scratch_data                        => fe_cell_iterator_update_scratch_data
    procedure                           :: update_integration                         => fe_cell_iterator_update_integration
    procedure                           :: update_cell_map                            => fe_cell_iterator_update_cell_map
    procedure                           :: update_cell_integrators                    => fe_cell_iterator_update_cell_integrators
    procedure                           :: set_cell_map                               => fe_cell_iterator_set_cell_map
    procedure                           :: set_cell_integrator                        => fe_cell_iterator_set_cell_integrator
    
    procedure, non_overridable :: get_quadrature_points_coordinates => fe_cell_iterator_get_quadrature_points_coordinates
    procedure, non_overridable :: get_det_jacobian                  => fe_cell_iterator_get_det_jacobian

    procedure, non_overridable          :: get_fe_space                               => fe_cell_iterator_get_fe_space
    procedure, non_overridable          :: get_num_fields                             => fe_cell_iterator_get_num_fields
    procedure, non_overridable, private :: get_fe_space_type                          => fe_cell_iterator_get_fe_space_type
    procedure, non_overridable          :: get_field_type                             => fe_cell_iterator_get_field_type

    procedure, non_overridable          :: get_field_blocks                           => fe_cell_iterator_get_field_blocks
    procedure, non_overridable          :: get_field_coupling                         => fe_cell_iterator_get_field_coupling
    procedure, non_overridable          :: get_num_dofs                               => fe_cell_iterator_get_num_dofs
    procedure, non_overridable          :: get_num_dofs_field                         => fe_cell_iterator_get_num_dofs_field
    procedure, non_overridable          :: get_field_fe_dofs                          => fe_cell_iterator_get_field_fe_dofs
    procedure, non_overridable          :: get_fe_dofs                                => fe_cell_iterator_get_fe_dofs
    procedure, non_overridable          :: get_order                                  => fe_cell_iterator_get_order
    
    procedure, non_overridable          :: get_max_order_single_field                 => fe_cell_iterator_get_max_order_single_field
    procedure, non_overridable          :: get_max_order_all_fields                   => fe_cell_iterator_get_max_order_all_fields
    generic                             :: get_max_order                              => get_max_order_single_field, &
                                                                                         get_max_order_all_fields
    procedure, non_overridable          :: at_strong_dirichlet_boundary               => fe_cell_iterator_at_strong_dirichlet_boundary
    procedure, non_overridable          :: has_fixed_dofs                             => fe_cell_iterator_has_fixed_dofs
    procedure, non_overridable          :: has_hanging_dofs                           => fe_cell_iterator_has_hanging_dofs
    procedure, non_overridable          :: fe_cell_iterator_set_at_strong_dirichlet_boundary_single_field
    procedure, non_overridable          :: fe_cell_iterator_set_at_strong_dirichlet_boundary_all_fields
    generic                             :: determine_at_strong_dirichlet_boundary     => fe_cell_iterator_set_at_strong_dirichlet_boundary_single_field, &
                                                                                         fe_cell_iterator_set_at_strong_dirichlet_boundary_all_fields
    procedure, non_overridable          :: determine_has_fixed_dofs                   => fe_cell_iterator_set_has_fixed_dofs
    procedure                           :: determine_has_hanging_dofs                 => fe_cell_iterator_set_has_hanging_dofs
    procedure, non_overridable          :: is_free_dof                                => fe_cell_iterator_is_free_dof
    procedure                           :: is_strong_dirichlet_dof                    => fe_cell_iterator_is_strong_dirichlet_dof   
    procedure, non_overridable          :: is_fixed_dof                               => fe_cell_iterator_is_fixed_dof
    procedure                           :: is_hanging_dof                             => fe_cell_iterator_is_hanging_dof
    procedure, non_overridable          :: compute_volume                             => fe_cell_iterator_compute_volume
    
    procedure, non_overridable          :: get_default_quadrature_degree              => fe_cell_iterator_get_default_quadrature_degree
    procedure, non_overridable          :: get_quadrature_degree                      => fe_cell_iterator_get_quadrature_degree
    procedure, non_overridable          :: set_quadrature_degree                      => fe_cell_iterator_set_quadrature_degree
    procedure                           :: get_quadrature                             => fe_cell_iterator_get_quadrature
    procedure                           :: get_cell_map                               => fe_cell_iterator_get_cell_map
    procedure                           :: get_cell_integrator                        => fe_cell_iterator_get_cell_integrator
    
    procedure, non_overridable, private :: fe_cell_iterator_get_fe_vef
    generic                             :: get_vef                                    => fe_cell_iterator_get_fe_vef
    procedure, non_overridable          :: get_reference_fe                           => fe_cell_iterator_get_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe                 => fe_cell_iterator_get_max_order_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe_id              => fe_cell_iterator_get_max_order_reference_fe_id
    procedure, non_overridable          :: get_reference_fe_id                        => fe_cell_iterator_get_reference_fe_id
    procedure, non_overridable          :: set_reference_fe_id                        => fe_cell_iterator_set_reference_fe_id
    procedure, non_overridable          :: is_void                                    => fe_cell_iterator_is_void
    procedure, non_overridable          :: create_own_dofs_on_vef_iterator            => fe_cell_iterator_create_own_dofs_on_vef_iterator
    procedure                 , private :: fe_cell_iterator_impose_strong_dirichlet_bcs_conforming
    procedure                 , private :: fe_cell_iterator_impose_strong_dirichlet_bcs_non_conforming
    generic                             :: impose_strong_dirichlet_bcs                => fe_cell_iterator_impose_strong_dirichlet_bcs_conforming, &
                                                                                         fe_cell_iterator_impose_strong_dirichlet_bcs_non_conforming
    procedure, private :: assembly_array => fe_cell_iterator_assembly_array
    procedure, private :: assembly_matrix => fe_cell_iterator_assembly_matrix
    procedure, private :: assembly_matrix_array => fe_cell_iterator_assembly_matrix_array
    procedure, private :: assembly_matrix_array_with_strong_bcs => fe_cell_iterator_assembly_matrix_array_with_strong_bcs
    generic            :: assembly                                   => assembly_array,        &
                                                                        assembly_matrix,       &
                                                                        assembly_matrix_array, &
                                                                        assembly_matrix_array_with_strong_bcs
    procedure, non_overridable          :: first_local_non_void                       => fe_cell_iterator_first_local_non_void

    ! Added by unfitted_fe_cell_iterator
    procedure                           :: get_boundary_quadrature                    => fe_cell_iterator_get_boundary_quadrature
    procedure                           :: get_boundary_piecewise_cell_map              => fe_cell_iterator_get_boundary_piecewise_cell_map
    procedure                           :: get_boundary_cell_map                        => fe_cell_iterator_get_boundary_cell_map
    procedure                           :: get_boundary_cell_integrator               => fe_cell_iterator_get_boundary_cell_integrator
    procedure                           :: update_boundary_integration                => fe_cell_iterator_update_boundary_integration
    
    procedure, non_overridable, private :: get_values_scalar           => fe_cell_iterator_get_values_scalar
    procedure, non_overridable, private :: get_values_vector           => fe_cell_iterator_get_values_vector
    generic                             :: get_values                  => get_values_scalar, get_values_vector 
    procedure, non_overridable, private :: get_gradients_scalar => fe_cell_iterator_get_gradients_scalar
    procedure, non_overridable, private :: get_gradients_vector => fe_cell_iterator_get_gradients_vector
    generic                             :: get_gradients        => get_gradients_scalar, get_gradients_vector 
    procedure, non_overridable, private :: get_divergences_vector => fe_cell_iterator_get_divergences_vector
    generic                             :: get_divergences        => get_divergences_vector
    procedure, non_overridable, private :: get_curls_vector => fe_cell_iterator_get_curls_vector
    generic                             :: get_curls        => get_curls_vector
    
        
    procedure, non_overridable, private :: fe_cell_iterator_evaluate_fe_function_scalar
    procedure, non_overridable, private :: fe_cell_iterator_evaluate_fe_function_vector
    procedure, non_overridable, private :: fe_cell_iterator_evaluate_fe_function_tensor
    generic :: evaluate_fe_function => fe_cell_iterator_evaluate_fe_function_scalar, &
    & fe_cell_iterator_evaluate_fe_function_vector, &
    & fe_cell_iterator_evaluate_fe_function_tensor

    procedure, non_overridable, private :: fe_cell_iterator_evaluate_gradient_fe_function_scalar
    procedure, non_overridable, private :: fe_cell_iterator_evaluate_gradient_fe_function_vector
    generic :: evaluate_gradient_fe_function => fe_cell_iterator_evaluate_gradient_fe_function_scalar, &
    & fe_cell_iterator_evaluate_gradient_fe_function_vector

    procedure, non_overridable, private :: apply_constraints                          => fe_cell_iterator_apply_constraints

  end type fe_cell_iterator_t
   
  type p_fe_cell_iterator_t
    class(fe_cell_iterator_t), pointer :: p 
  end type p_fe_cell_iterator_t
  
  type :: base_fe_vef_iterator_t
    private
    class(vef_iterator_t), allocatable :: vef
  contains
     ! Methods exploiting ("inherited from") vef_t
     procedure                           :: first                     => base_fe_vef_iterator_first
     procedure                           :: next                      => base_fe_vef_iterator_next
     procedure                           :: has_finished              => base_fe_vef_iterator_has_finished
     procedure                           :: set_gid                   => base_fe_vef_iterator_set_gid
     procedure                           :: get_gid                   => base_fe_vef_iterator_get_gid
     procedure, non_overridable          :: get_nodes_coordinates     => base_fe_vef_iterator_get_nodes_coordinates
     procedure, non_overridable          :: get_set_id                => base_fe_vef_iterator_get_set_id
     
     procedure, non_overridable          :: get_dim                   => base_fe_vef_iterator_get_dim
     procedure, non_overridable          :: is_at_boundary            => base_fe_vef_iterator_is_at_boundary
     procedure, non_overridable          :: is_local                  => base_fe_vef_iterator_is_local
     procedure, non_overridable          :: is_ghost                  => base_fe_vef_iterator_is_ghost
     procedure, non_overridable          :: is_at_interface           => base_fe_vef_iterator_is_at_interface
     procedure, non_overridable          :: is_facet                   => base_fe_vef_iterator_is_facet
     
     procedure                           :: get_num_cells_around      => base_fe_vef_iterator_get_num_cells_around
     procedure, non_overridable          :: base_fe_vef_iterator_get_cell_around
     generic                             :: get_cell_around           => base_fe_vef_iterator_get_cell_around
  end type base_fe_vef_iterator_t

  type , extends(base_fe_vef_iterator_t) :: fe_vef_iterator_t
    private
    class(serial_fe_space_t), pointer :: fe_space => NULL()
  contains
     procedure                           :: create                            => fe_vef_iterator_create
     procedure                           :: free                              => fe_vef_iterator_free
     final                               :: fe_vef_iterator_final
     
     procedure, non_overridable          :: is_proper                         => fe_vef_iterator_is_proper
     procedure, non_overridable          :: is_fe_facet                       => fe_vef_iterator_is_fe_facet
     procedure, non_overridable, private :: has_free_dofs                     => fe_vef_iterator_has_free_dofs
     
     procedure                 , private :: fe_vef_iterator_get_fe_around
     generic                             :: get_cell_around                   => fe_vef_iterator_get_fe_around
     
     procedure, non_overridable          :: get_num_improper_cells_around     => fe_vef_iterator_get_num_improper_cells_around
     procedure, non_overridable, private :: fe_vef_iterator_get_improper_cell_around
     procedure, non_overridable, private :: fe_vef_iterator_get_improper_fe_around
     generic                             :: get_improper_cell_around          => fe_vef_iterator_get_improper_cell_around, &
                                                                                 fe_vef_iterator_get_improper_fe_around         
     procedure                           :: get_improper_cell_around_ivef     => fe_vef_iterator_get_improper_cell_around_ivef
     procedure                           :: get_improper_cell_around_subvef   => fe_vef_iterator_get_improper_cell_around_subvef  
  end type fe_vef_iterator_t
    
  type, extends(fe_vef_iterator_t) :: fe_facet_iterator_t
    private
    integer(ip)                              :: facet_gid
    ! Scratch data to support FE face integration
    class(fe_cell_iterator_t) , allocatable  :: fe1
    class(fe_cell_iterator_t) , allocatable  :: fe2
    type(p_fe_cell_iterator_t)               :: fes_around(2)
    type(facet_maps_t)  , pointer           :: facet_maps => NULL()
    type(p_facet_integrator_t), allocatable :: facet_integrators(:)
   contains
    procedure                           :: create                         => fe_facet_iterator_create
    procedure                           :: free                           => fe_facet_iterator_free
    procedure, non_overridable          :: allocate_scratch_data          => fe_facet_iterator_allocate_scratch_data
    procedure, non_overridable          :: free_scratch_data              => fe_facet_iterator_free_scratch_data
    procedure                           :: first                          => fe_facet_iterator_first
    procedure                           :: next                           => fe_facet_iterator_next
    procedure                           :: has_finished                   => fe_facet_iterator_has_finished
    procedure, non_overridable          :: set_gid                        => fe_facet_iterator_set_gid
    procedure, non_overridable          :: get_gid                        => fe_facet_iterator_get_gid
    procedure, non_overridable          :: is_at_field_boundary          => fe_facet_iterator_is_at_field_boundary
    procedure, non_overridable          :: is_at_field_interior          => fe_facet_iterator_is_at_field_interior
    procedure, non_overridable          :: get_num_cells_around          => fe_facet_iterator_get_num_cells_around
    procedure, non_overridable, private :: fe_vef_iterator_get_fe_around => fe_facet_iterator_get_fe_around
    procedure, non_overridable, private :: update_fes_around             => fe_facet_iterator_update_fes_around
    procedure, non_overridable          :: update_integration            => fe_facet_iterator_update_integration
    procedure, non_overridable          :: get_num_dofs_field            => fe_facet_iterator_get_num_dofs_field
    procedure, non_overridable, private :: fe_facet_iterator_assembly_array
    procedure, non_overridable, private :: fe_facet_iterator_assembly_matrix
    procedure, non_overridable, private :: fe_facet_iterator_assembly_matrix_array
    generic                             :: assembly                      => fe_facet_iterator_assembly_array,  &
                                                                            fe_facet_iterator_assembly_matrix, &
                                                                            fe_facet_iterator_assembly_matrix_array
    procedure, non_overridable          :: get_fe_space                  => fe_facet_iterator_get_fe_space
    procedure, non_overridable          :: get_fe_dofs                   => fe_facet_iterator_get_fe_dofs
    procedure, non_overridable          :: get_default_quadrature_degree  => fe_facet_iterator_get_default_quadrature_degree
    procedure, non_overridable          :: get_quadrature_degree          => fe_facet_iterator_get_quadrature_degree
    procedure, non_overridable          :: set_quadrature_degree          => fe_facet_iterator_set_quadrature_degree
    procedure, non_overridable          :: get_quadrature                 => fe_facet_iterator_get_quadrature

    procedure, non_overridable, private :: get_facet_maps                => fe_facet_iterator_get_facet_map
    procedure, non_overridable          :: update_facet_maps             => fe_facet_iterator_update_facet_maps
    procedure, non_overridable          :: update_facet_integrators      => fe_facet_iterator_update_facet_integrators
    procedure, non_overridable, private :: get_facet_integrator          => fe_facet_iterator_get_facet_integrator
    procedure, non_overridable          :: compute_surface                => fe_facet_iterator_compute_surface
    procedure                 , private :: compute_fe_facet_permutation_index => fe_facet_iterator_compute_fe_facet_permutation_index
    procedure, non_overridable          :: get_lpos_within_cell_around   => fe_facet_iterator_get_lpos_within_cell_around
    procedure, non_overridable          :: get_facet_permutation_index   => fe_facet_iterator_get_fe_facet_permutation_index 
    procedure, non_overridable, private :: get_subfacet_lid_cell_around  => fe_facet_iterator_get_subfacet_lid_cell_around
    
    procedure, non_overridable :: get_quadrature_points_coordinates => fe_facet_iterator_get_quadrature_points_coordinates
    procedure, non_overridable :: get_normals                       => fe_facet_iterator_get_normals
    procedure, non_overridable :: get_det_jacobian                  => fe_facet_iterator_get_det_jacobian
    procedure, non_overridable :: compute_characteristic_length     => fe_facet_iterator_compute_characteristic_length
    
    
    procedure, non_overridable :: get_values_scalar     => fe_facet_iterator_get_values_scalar
    procedure, non_overridable :: get_values_vector     => fe_facet_iterator_get_values_vector
    generic                    :: get_values            => get_values_scalar, get_values_vector
    procedure, non_overridable :: get_gradients_scalar  => fe_facet_iterator_get_gradients_scalar
    generic                    :: get_gradients         => get_gradients_scalar
    procedure, non_overridable :: get_curls             => fe_facet_iterator_get_curls_vector 
    
    procedure, non_overridable :: get_active_cell_id    => fe_facet_iterator_get_active_cell_id
    
    procedure, non_overridable, private :: fe_facet_iterator_evaluate_fe_function_scalar
    procedure, non_overridable, private :: fe_facet_iterator_evaluate_fe_function_vector
    procedure, non_overridable, private :: fe_facet_iterator_evaluate_fe_function_tensor
    generic :: evaluate_fe_function => fe_facet_iterator_evaluate_fe_function_scalar, &
    & fe_facet_iterator_evaluate_fe_function_vector, &
    & fe_facet_iterator_evaluate_fe_function_tensor

    procedure, non_overridable, private :: fe_facet_iterator_evaluate_gradient_fe_function_scalar
    procedure, non_overridable, private :: fe_facet_iterator_evaluate_gradient_fe_function_vector
    generic :: evaluate_gradient_fe_function => fe_facet_iterator_evaluate_gradient_fe_function_scalar, &
    & fe_facet_iterator_evaluate_gradient_fe_function_vector
    
    procedure, non_overridable :: get_current_qpoints_perm => fe_facet_iterator_get_current_qpoints_perm
    
  end type fe_facet_iterator_t
      
  integer(ip), parameter :: fe_space_type_cg                        = 0 ! H^1 conforming FE space
  integer(ip), parameter :: fe_space_type_dg                        = 1 ! L^2 conforming FE space + .not. H^1 conforming (weakly imposed via face integration)
  integer(ip), parameter :: fe_space_type_dg_conforming             = 2 ! DG approximation of L^2 spaces (does not involve coupling by face)
  
  integer(ip), parameter :: fe_space_default_quadrature_degree_flag = -1000

  type, extends(base_fe_space_t) :: serial_fe_space_t 
     !private      ! UNDER QUARANTINE
     type(cell_integrator_duties_t) :: cell_integrator_duties

     ! Reference FE container
     integer(ip)                                 :: reference_fes_size
     type(p_reference_fe_t)        , allocatable :: reference_fes(:)

     ! Finite Element-related integration containers
     type(std_vector_quadrature_t)               :: cell_quadratures
     type(std_vector_cell_map_t)                 :: cell_maps
     type(std_vector_cell_integrator_t)          :: cell_integrators
     type(std_vector_integer_ip_t)               :: cell_quadratures_degree
     
     ! Mapping of FEs to reference FE and FEs-related integration containers
     integer(ip)                   , allocatable :: set_ids_to_reference_fes(:,:)         ! The (num_fields,num_cell_set_ids) table provided by the user to %create()
     type(std_vector_integer_ip_t) , allocatable :: field_cell_to_ref_fes(:)
     type(std_vector_integer_ip_t)               :: max_order_reference_fe_id_x_cell      ! Stores Key=max_order_reference_fe_id for all FEs
     type(hash_table_ip_ip_t)                    :: cell_quadratures_and_maps_position    ! Key = [geo_reference_fe_id,quadrature_degree]
     type(hash_table_ip_ip_t)                    :: cell_integrators_position             ! Key = [geo_reference_fe_id,quadrature_degree,reference_fe_id]
     
     ! Finite Face-related integration containers
     type(std_vector_quadrature_t)               :: facet_quadratures
     type(std_vector_facet_maps_t)               :: facet_maps
     type(std_vector_facet_integrator_t)         :: facet_integrators
     type(std_vector_integer_ip_t)               :: facet_quadratures_degree
     
     ! Mapping of Finite Faces and integration containers
     type(std_vector_integer_ip_t)               :: max_order_field_cell_to_ref_fes_face ! Stores max_order_field_cell_to_ref_fes_face for all faces
     type(hash_table_ip_ip_t)                    :: facet_quadratures_position  ! Key = [quadrature_degree, 
                                                                                  !       left_geo_reference_fe_id,
                                                                                  !       right_geo_reference_fe_id (with 0 for boundary faces)]
     type(hash_table_ip_ip_t)                    :: facet_integrators_position  ! Key = [quadrature_degree,
                                                                                  !       left_reference_fe_id,
                                                                                  !       left_reference_fe_id (with 0 for boundary faces)]
     
     ! Member variables to provide support to fe_facet_iterator_t
     type(std_vector_integer_ip_t)               :: facet_gids
     type(std_vector_integer_ip_t)               :: facet_permutation_indices
     
     ! DoF identifiers associated to each FE and field within FE
     type(std_vector_integer_ip_t) , allocatable :: ptr_dofs_x_fe(:)
     type(std_vector_integer_ip_t)               :: lst_dofs_gids
     
     ! Strong Dirichlet BCs-related member variables
     class(conditions_t)           , pointer     :: conditions    => NULL()
     integer(ip)                                 :: num_fixed_dofs
     integer(ip)                                 :: num_hanging_dofs
     integer(ip)                                 :: num_dirichlet_dofs
     type(std_vector_logical_t)    , allocatable :: at_strong_dirichlet_boundary_x_fe(:)
     type(std_vector_logical_t)    , allocatable :: has_fixed_dofs_x_fe(:)
     type(std_vector_logical_t)    , allocatable :: has_hanging_dofs_x_fe(:)
     
     ! Descriptor of the block layout selected for the PDE system at hand
     type(block_layout_t)                        :: block_layout
     
     ! ( Polymorphic ) pointer to a triangulation it was created from
     class(triangulation_t)            , pointer :: triangulation => NULL()
     
     ! Constraining DOFs arrays to give support to h-adaptivity
     type(std_vector_integer_ip_t)               :: ptr_constraining_free_dofs
     type(std_vector_integer_ip_t)               :: ptr_constraining_dirichlet_dofs
     type(std_vector_integer_ip_t)               :: constraining_free_dofs
     type(std_vector_real_rp_t)                  :: constraining_free_dofs_coefficients
     type(std_vector_integer_ip_t)               :: constraining_dirichlet_dofs
     type(std_vector_real_rp_t)                  :: constraining_dirichlet_dofs_coefficients
     type(std_vector_real_rp_t)                  :: constraints_independent_term
     
   contains
     procedure                           :: serial_fe_space_create_same_reference_fes_on_all_cells
     procedure                           :: serial_fe_space_create_different_ref_fes_between_cells
     generic                             :: create                                       => serial_fe_space_create_same_reference_fes_on_all_cells,&
                                                                                            serial_fe_space_create_different_ref_fes_between_cells
     procedure                           :: free                                         => serial_fe_space_free
     procedure                           :: print                                        => serial_fe_space_print
     procedure, non_overridable          :: allocate_and_fill_reference_fes              => serial_fe_space_allocate_and_fill_reference_fes
     procedure, non_overridable, private :: free_reference_fes                           => serial_fe_space_free_reference_fes
     
     procedure                           :: allocate_and_fill_set_ids_to_reference_fes   => sfs_allocate_and_fill_set_ids_to_reference_fes
     procedure                           :: free_set_ids_to_reference_fes                => sfs_free_set_ids_to_reference_fes
     procedure, non_overridable          :: allocate_field_cell_to_ref_fes               => serial_fe_space_allocate_field_cell_to_ref_fes
     procedure, non_overridable, private :: free_field_cell_to_ref_fes                   => serial_fe_space_free_field_cell_to_ref_fes
     procedure, non_overridable          :: fill_field_cell_to_ref_fes_same_on_all_cells => serial_fe_space_fill_field_cell_to_ref_fes_same_on_all_cells
     procedure, non_overridable          :: fill_field_cell_to_ref_fes_different_ref_fes_between_cells => sfes_fill_field_cell_to_ref_fes_different_ref_fes_between_cells
     
     procedure, non_overridable          :: check_cell_vs_fe_topology_consistency        => serial_fe_space_check_cell_vs_fe_topology_consistency
     
     procedure, non_overridable          :: allocate_and_fill_fe_space_type_x_field      => serial_fe_space_allocate_and_fill_fe_space_type_x_field
     procedure, non_overridable, private :: free_fe_space_type_x_field                   => serial_fe_space_free_fe_space_type_x_field
     
     procedure, non_overridable          :: allocate_and_init_ptr_lst_dofs_gids          => serial_fe_space_allocate_and_init_ptr_lst_dofs_gids
     procedure, non_overridable          :: copy_ptr_lst_dofs                            => serial_fe_space_copy_ptr_lst_dofs
     procedure, non_overridable, private :: free_ptr_lst_dofs                            => serial_fe_space_free_ptr_lst_dofs
     
     procedure, non_overridable          :: allocate_and_init_at_strong_dirichlet_bound  => serial_fe_space_allocate_and_init_at_strong_dirichlet_bound  
     procedure, non_overridable, private :: free_at_strong_dirichlet_boundary               => serial_fe_space_free_at_strong_dirichlet_boundary
     procedure, non_overridable          :: allocate_and_init_has_fixed_dofs             => serial_fe_space_allocate_and_init_has_fixed_dofs
     procedure, non_overridable, private :: free_has_fixed_dofs                          => serial_fe_space_free_has_fixed_dofs
     procedure                           :: set_up_strong_dirichlet_bcs                  => serial_fe_space_set_up_strong_dirichlet_bcs
     procedure, non_overridable, private :: set_up_strong_dirichlet_bcs_on_vef_and_field => serial_fe_space_set_up_strong_dirichlet_bcs_on_vef_and_field
     procedure                           :: interpolate_scalar_function                  => serial_fe_space_interpolate_scalar_function 
     procedure                           :: interpolate_vector_function                  => serial_fe_space_interpolate_vector_function
     generic                             :: interpolate                                  => interpolate_scalar_function, interpolate_vector_function     
     procedure                           :: interpolate_dirichlet_values                 => serial_fe_space_interpolate_dirichlet_values
     procedure                           :: project_dirichlet_values_curl_conforming     => serial_fe_space_project_dirichlet_values_curl_conforming
     procedure, non_overridable, private :: allocate_and_fill_fields_to_project_         => serial_fe_space_allocate_and_fill_fields_to_project_
     procedure, non_overridable, private :: allocate_and_fill_offset_component           => serial_fe_space_allocate_and_fill_offset_component
     procedure, non_overridable, private :: allocate_and_fill_global2subset_and_inverse  => serial_fe_space_allocate_and_fill_global2subset_and_inverse 
     procedure, non_overridable, private :: get_function_scalar_components               => serial_fe_space_get_function_scalar_components
     procedure, non_overridable, private :: evaluate_vector_function_scalar_components   => serial_fe_space_evaluate_vector_function_scalar_components
     procedure, non_overridable, private :: project_curl_conforming_compute_elmat_elvec  => serial_fe_space_project_curl_conforming_compute_elmat_elvec     
     procedure, non_overridable, private :: allocate_and_init_cell_quadratures_degree      => serial_fe_space_allocate_and_init_cell_quadratures_degree
     procedure, non_overridable, private :: free_cell_quadratures_degree                   => serial_fe_space_free_cell_quadratures_degree
     
     procedure, non_overridable, private :: free_max_order_reference_fe_id_x_cell        => serial_fe_space_free_max_order_reference_fe_id_x_cell
     procedure, non_overridable, private :: compute_max_order_reference_fe_id_x_cell     => serial_fe_space_compute_max_order_reference_fe_id_x_cell         
     
     procedure, non_overridable          :: set_up_cell_integration                    => serial_fe_space_set_up_cell_integration
     procedure, non_overridable, private :: free_fe_integration                          => serial_fe_space_free_fe_integration
     procedure, non_overridable, private :: generate_cell_quadratures_position_key         => serial_fe_space_generate_cell_quadratures_position_key
     procedure, non_overridable, private :: generate_cell_integrators_position_key  => serial_fe_space_generate_cell_integrators_position_key

     procedure, non_overridable, private :: allocate_and_init_facet_quadratures_degree => serial_fe_space_allocate_and_init_facet_quadratures_degree
     procedure, non_overridable, private :: free_facet_quadratures_degree              => serial_fe_space_free_facet_quadratures_degree
     
     procedure, non_overridable, private :: free_max_order_field_cell_to_ref_fes_face     => serial_fe_space_free_max_order_field_cell_to_ref_fes_face
     procedure, non_overridable, private :: compute_max_order_field_cell_to_ref_fes_face  => serial_fe_space_compute_max_order_field_cell_to_ref_fes_face   
     
     procedure                 , private :: fill_facet_gids                    => serial_fe_space_fill_facet_gids
     procedure, non_overridable, private :: free_facet_gids                    => serial_fe_space_free_facet_gids
     
     procedure, non_overridable, private :: compute_facet_permutation_indices             => serial_fe_space_compute_facet_permutation_indices
     procedure, non_overridable, private :: free_facet_permutation_indices                => serial_fe_space_free_facet_permutation_indices
     
     procedure, non_overridable          :: set_up_facet_integration                 => serial_fe_space_set_up_facet_integration
     procedure, non_overridable, private :: free_facet_integration                   => serial_fe_space_free_facet_integration
     procedure, non_overridable, private :: generate_facet_quadratures_position_key  => serial_fe_space_facet_quadratures_position_key
     procedure, non_overridable, private :: generate_facet_integrators_position_key  => serial_fe_space_facet_integrators_position_key

     procedure                           :: create_dof_values                        => serial_fe_space_create_dof_values
     procedure                           :: generate_global_dof_numbering            => serial_fe_space_generate_global_dof_numbering
     procedure                           :: allocate_num_dofs_x_field                =>  serial_fe_space_allocate_num_dofs_x_field
     procedure                 , private :: count_dofs                               => serial_fe_space_count_dofs
     procedure                 , private :: list_dofs                                => serial_fe_space_list_dofs
     procedure                 , private :: renum_dofs_block                         => serial_fe_space_renum_dofs_block
 
     ! Getters
     procedure                           :: get_num_dims                              => serial_fe_space_get_num_dims
     procedure, non_overridable          :: get_num_reference_fes                     => serial_fe_space_get_num_reference_fes
     procedure, non_overridable          :: get_reference_fe                          => serial_fe_space_get_reference_fe
     procedure, non_overridable          :: get_field_type                            => serial_fe_space_get_field_type 
     procedure, non_overridable, private :: determine_fe_space_type                   => serial_fe_space_determine_fe_space_type
     procedure, non_overridable          :: get_num_components                        => serial_fe_space_get_num_components
     procedure, non_overridable          :: get_max_num_shape_functions               => serial_fe_space_get_max_num_shape_functions
     procedure, non_overridable          :: get_max_num_dofs_on_a_cell                => serial_fe_space_get_max_num_dofs_on_a_cell
     procedure, non_overridable          :: get_max_num_quadrature_points             => serial_fe_space_get_max_num_quadrature_points
     procedure, non_overridable          :: get_max_num_nodal_quadrature_points       => serial_fe_space_get_max_num_nodal_quadrature_points
     procedure, non_overridable          :: get_max_num_facet_quadrature_points        => serial_fe_space_get_max_num_facet_quadrature_points     
     procedure, non_overridable          :: get_max_order                                => serial_fe_space_get_max_order
     procedure, non_overridable          :: get_triangulation                            => serial_fe_space_get_triangulation
     procedure, non_overridable          :: set_triangulation                            => serial_fe_space_set_triangulation
     procedure                           :: get_environment                              => serial_fe_space_get_environment
     procedure, non_overridable          :: get_conditions                               => serial_fe_space_get_conditions
     procedure, non_overridable          :: set_conditions                               => serial_fe_space_set_conditions
     procedure                           :: get_num_fixed_dofs                           => serial_fe_space_get_num_fixed_dofs
     procedure                           :: get_num_dirichlet_dofs                       => serial_fe_space_get_num_dirichlet_dofs
     procedure                           :: get_num_hanging_dofs                         => serial_fe_space_get_num_hanging_dofs
     procedure                           :: get_num_blocks                            => serial_fe_space_get_num_blocks
     procedure                           :: get_field_blocks                             => serial_fe_space_get_field_blocks
     procedure                           :: get_field_coupling                           => serial_fe_space_get_field_coupling
     procedure                           :: get_block_num_dofs                        => serial_fe_space_get_block_num_dofs
     procedure                           :: set_block_num_dofs                        => serial_fe_space_set_block_num_dofs
     procedure, non_overridable          :: get_block_layout                             => serial_fe_space_get_block_layout
     
     ! fes, fe_vefs and fe_faces traversals-related TBPs
     procedure                           :: create_fe_cell_iterator                           => serial_fe_space_create_fe_cell_iterator
     procedure                           :: free_fe_cell_iterator                             => serial_fe_space_free_fe_cell_iterator
     procedure, non_overridable          :: create_fe_vef_iterator                       => serial_fe_space_create_fe_vef_iterator     
     procedure, non_overridable          :: create_itfc_fe_vef_iterator                  => serial_fe_space_create_itfc_fe_vef_iterator     
     procedure, non_overridable          :: free_fe_vef_iterator                         => serial_fe_space_free_fe_vef_iterator
     procedure, non_overridable          :: create_fe_facet_iterator                      => serial_fe_space_create_fe_facet_iterator
     procedure, non_overridable          :: free_fe_facet_iterator                        => serial_fe_space_free_fe_facet_iterator

     procedure, non_overridable          :: allocate_and_init_has_hanging_dofs_x_fe       => serial_fe_space_allocate_and_init_has_hanging_dofs_x_fe
     procedure, non_overridable, private :: free_has_hanging_dofs_x_fe                    => serial_fe_space_free_has_hanging_dofs_x_fe
     
     procedure, non_overridable, private :: free_constraining_dofs_arrays                 => serial_fe_space_free_constraining_dofs_arrays
     procedure, non_overridable, private :: free_ptr_constraining_free_dofs               => serial_fe_space_free_ptr_constraining_free_dofs
     procedure, non_overridable, private :: free_constraining_free_dofs                   => serial_fe_space_free_constraining_free_dofs
     procedure, non_overridable, private :: free_constraining_free_dofs_coefficients      => serial_fe_space_free_constraining_free_dofs_coefficients
     procedure, non_overridable, private :: free_ptr_constraining_dirichlet_dofs          => serial_fe_space_free_ptr_constraining_dirichlet_dofs
     procedure, non_overridable, private :: free_constraining_dirichlet_dofs              => serial_fe_space_free_constraining_dirichlet_dofs
     procedure, non_overridable, private :: free_constraining_dirichlet_dofs_coefficients => serial_fe_space_free_constraining_dirichlet_dofs_coefficients
     procedure, non_overridable, private :: free_constraints_independent_term             => serial_fe_space_free_constraints_independent_term

     procedure, non_overridable          :: setup_hanging_node_constraints               => serial_fe_space_setup_hanging_node_constraints

     procedure, non_overridable, private :: project_field_cell_to_ref_fes                 => serial_fe_space_project_field_cell_to_ref_fes
     procedure, non_overridable, private :: project_fe_integration_arrays                 => serial_fe_space_project_fe_integration_arrays
     procedure, non_overridable, private :: project_facet_integration_arrays              => serial_fe_space_project_facet_integration_arrays

     procedure,                  private :: serial_fe_space_refine_and_coarsen_single_fe_function
     procedure,                  private :: serial_fe_space_refine_and_coarsen_fe_function_array
     generic                             :: refine_and_coarsen                            => serial_fe_space_refine_and_coarsen_single_fe_function, &
                                                                                             serial_fe_space_refine_and_coarsen_fe_function_array

     procedure, non_overridable          :: update_hanging_dof_values                     => serial_fe_space_update_hanging_dof_values
     
#ifndef ENABLE_P4EST
    procedure, non_overridable           :: not_enabled_error                             => serial_fe_space_not_enabled_error
#endif     

 end type serial_fe_space_t  

 public :: serial_fe_space_t, serial_fe_space_set_up_strong_dirichlet_bcs
 public :: fe_space_type_cg, fe_space_type_dg
 public :: fe_cell_iterator_t
 public :: fe_vef_iterator_t
 public :: fe_facet_iterator_t

 type base_fe_object_iterator_t
   private
   type(object_iterator_t) :: object
 contains
   ! Methods "inherited" from object_iterator_t
   procedure                            :: first                                 => base_fe_object_iterator_first
   procedure                            :: next                                  => base_fe_object_iterator_next
   procedure                            :: set_gid                               => base_fe_object_iterator_set_gid
   procedure                            :: has_finished                          => base_fe_object_iterator_has_finished
   procedure, non_overridable           :: get_gid                               => base_fe_object_iterator_get_gid
   procedure, non_overridable           :: get_ggid                              => base_fe_object_iterator_get_ggid
   procedure, non_overridable           :: get_dim                               => base_fe_object_iterator_get_dim
   procedure, non_overridable           :: get_num_parts_around                  => base_fe_object_iterator_get_num_parts_around
   procedure, non_overridable           :: get_num_subparts_around               => base_fe_object_iterator_get_num_subparts_around
   procedure, non_overridable           :: create_parts_around_iterator          => base_fe_object_iterator_create_parts_around_iterator
   procedure, non_overridable           :: create_subparts_around_iterator       => base_fe_object_iterator_create_subparts_around_iterator
   procedure, non_overridable           :: get_num_vefs                          => base_fe_object_iterator_get_num_vefs
   procedure, non_overridable, private  :: base_fe_object_iterator_get_vef
   generic                              :: get_vef                               => base_fe_object_iterator_get_vef
 end type base_fe_object_iterator_t
 
 
 type, extends(base_fe_object_iterator_t) :: fe_object_iterator_t
    private
    type(par_fe_space_t), pointer :: fe_space => NULL()
  contains
    procedure, non_overridable, private  :: create                                => fe_object_iterator_create
    procedure, non_overridable, private  :: free                                  => fe_object_iterator_free
    final                                :: fe_object_iterator_free_final
    ! Own methods of fe_object_iterator_t
    procedure, non_overridable, private  :: fe_object_iterator_get_fe_vef
    generic                              :: get_vef                               => fe_object_iterator_get_fe_vef 
    procedure, non_overridable, private  :: get_face                              => fe_object_iterator_get_fe_face
    procedure, non_overridable           :: get_num_facets                         => fe_object_iterator_get_num_facets
    procedure, non_overridable           :: get_num_coarse_dofs                => fe_object_iterator_get_num_coarse_dofs
    procedure, non_overridable           :: create_own_coarse_dofs_iterator       => fe_object_iterator_create_own_coarse_dofs_iterator
    procedure, non_overridable           :: create_faces_object_iterator          => fe_object_iterator_create_faces_object_iterator
  end type fe_object_iterator_t
    
  
  ! These parameter constants are used in order to generate a unique (non-consecutive) 
  ! but consistent across MPI tasks global ID (integer(igp)) of a given DoF.
  ! See type(par_fe_space_t)%generate_non_consecutive_dof_ggid()
  integer(ip), parameter :: cell_ggid_shift              = 44
  integer(ip), parameter :: dofs_x_reference_fe_shift = 14
  integer(ip), parameter :: num_fields_shift         = igp*8-(cell_ggid_shift+dofs_x_reference_fe_shift)-1
  
  ! These three parameter constants are thought be used as FPL keys. The corresponding pairs 
  ! <XXXkey,.true.|.false.> in FPL let the user to control whether or not coarse vertex, edge, or 
  ! face DoFs are to be included into coarse_fe_space_t. These parameters might be used during 
  ! coarse_fe_space_t set-up, as well as by the deferred TBP methods corresponding to 
  ! class(coarse_fe_handler_t).
  character(len=*), parameter :: coarse_space_use_vertices_key = 'coarse_space_use_vertices_key'
  character(len=*), parameter :: coarse_space_use_edges_key    = 'coarse_space_use_edges_key'
  character(len=*), parameter :: coarse_space_use_faces_key    = 'coarse_space_use_faces_key'

  type :: p_l1_coarse_fe_handler_t
    class(l1_coarse_fe_handler_t), pointer :: p
  end type p_l1_coarse_fe_handler_t
  
 type, extends(serial_fe_space_t) :: par_fe_space_t
   private   
   ! Multilevel fe space   
   ! It is the equivalent to the "element_to_dof" at the finer level
   ! Pointers to the start/end of coarse DoFs GIDs of each field (lst_coarse_dofs)
   integer(ip)                   , allocatable :: ptr_coarse_dofs_x_field(:)
   ! List of coarse DoFs GIDs
   integer(ip)                   , allocatable :: lst_coarse_dofs(:)
   
   ! Coarse DoFs GIDs on top of coarse n_faces per field
   ! It provides (for every field) what we get from the reference_fe_t at the
   ! finer level
   type(list_t), allocatable                   :: own_coarse_dofs_x_field(:)
   
   ! GIDs of the coarse DoFs. For their generation, we though to be a good idea 
   ! to take as basis the GIDs of the VEFs objects they are built from (instead of
   ! generating them from scratch, which in turn would imply further communication).
   type(allocatable_array_igp1_t), allocatable :: coarse_dof_gids_x_field(:)
   
   ! Polymorphic data type in charge of filling some of the member variables above
   ! (so far, lst_coarse_dofs + own_coarse_dofs_x_field)
   type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)

   ! Member variables to provide support to fe_facet_iterator_t 
   ! to iterate the faces of an object
   type(list_t)                                :: faces_object
   
   ! Scratch data required for non-conforming triangulations
   type(std_vector_integer_ip_t)               :: max_part_id_vefs
   type(std_vector_integer_ip_t)               :: my_part_id_vefs
   type(std_vector_integer_ip_t)               :: rcv_my_part_id_vefs 
   type(std_vector_integer_ip_t)               :: ptr_dofs_field
   type(std_vector_integer_igp_t)              :: lst_dofs_ggids
   type(std_vector_integer_ip_t)               :: ptr_ghosts_per_local_cell
   type(std_vector_integer_ip_t)               :: lst_ghosts_per_local_cell
   type(std_vector_integer_ip_t)               :: num_cells_to_send_x_local_cell
   type(std_vector_integer_ip_t)               :: snd_ptrs_complete_itfc_couplings
   type(std_vector_integer_ip_t)               :: snd_leids_complete_itfc_couplings
   type(std_vector_integer_ip_t)               :: rcv_ptrs_complete_itfc_couplings
   type(std_vector_integer_ip_t)               :: ptr_ghosts_per_ghost_cell
   type(std_vector_integer_ip_t)               :: lst_ghosts_per_ghost_cell
   type(std_vector_integer_ip_t)               :: rcv_my_part_id_vefs_complete_itfc_couplings 

   type(std_vector_real_rp_t)                  :: old_fe_function_nodal_values
   type(std_vector_real_rp_t)                  :: new_fe_function_nodal_values
 contains
   procedure          :: serial_fe_space_create_same_reference_fes_on_all_cells    => par_fe_space_serial_create_same_reference_fes_on_all_cells 
   procedure          :: serial_fe_space_create_different_ref_fes_between_cells    => par_fe_space_serial_create_different_ref_fes_between_cells 
   procedure          :: par_fe_space_create_same_reference_fes_on_all_cells 
   procedure          :: par_fe_space_create_different_ref_fes_between_cells
   generic                                     :: create                                          => par_fe_space_create_same_reference_fes_on_all_cells, &
                                                                                                     par_fe_space_create_different_ref_fes_between_cells
   procedure                         , private :: allocate_and_fill_coarse_fe_handlers            => par_fe_space_allocate_and_fill_coarse_fe_handlers
   procedure                         , private :: free_coarse_fe_handlers                         => par_fe_space_free_coarse_fe_handlers
   procedure                                   :: generate_global_dof_numbering                   => par_fe_space_generate_global_dof_numbering
   procedure                                   :: renum_dofs_first_interior_then_interface        => par_fe_space_renum_dofs_first_interior_then_interface
   procedure        , non_overridable          :: compute_num_global_dofs_and_their_ggids         => par_fe_space_compute_num_global_dofs_and_their_ggids
   
   procedure        , non_overridable, private :: compute_blocks_dof_import                                    => par_fe_space_compute_blocks_dof_import
   procedure        , non_overridable, private :: compute_dof_import                                           => par_fe_space_compute_dof_import
   procedure        , non_overridable, private :: compute_dof_import_non_conforming_mesh                       => par_fe_space_compute_dof_import_non_conforming_mesh
   procedure        , non_overridable, private :: compute_raw_interface_data_by_continuity                     => par_fe_space_compute_raw_interface_data_by_continuity
   procedure        , non_overridable, private :: compute_raw_interface_data_by_continuity_non_conforming_mesh => pfs_compute_raw_itfc_data_by_continuity_non_conforming_mesh
   procedure        , non_overridable, private :: raw_interface_data_by_continuity_decide_owner                => par_fe_space_raw_interface_data_by_continuity_decide_owner
   procedure        , non_overridable, private :: compute_max_part_id_my_part_id_and_dofs_ggids_field          => pfs_compute_max_part_id_my_part_id_and_dofs_ggids_field 
   procedure        , non_overridable, private :: compute_exchange_control_data_to_complete_itfc_couplings     => pfs_compute_exchange_control_data_to_complete_itfc_couplings 

   procedure        , non_overridable, private :: compute_raw_interface_data_by_facet_integ                => par_fe_space_compute_raw_interface_data_by_facet_integ
   procedure        , non_overridable, private :: compute_ubound_num_itfc_couplings_by_continuity          => pfs_compute_ubound_num_itfc_couplings_by_continuity
   procedure        , non_overridable, private :: compute_ubound_num_itfc_couplings_by_continuity_nc_mesh  => pfs_compute_ubound_num_itfc_couplings_by_continuity_nc_mesh
   procedure        , non_overridable, private :: compute_ubound_num_itfc_couplings_by_facet_integ         => pfs_compute_ubound_num_itfc_couplings_by_facet_integ
   procedure, nopass, non_overridable, private :: generate_non_consecutive_dof_ggid                        => par_fe_space_generate_non_consecutive_dof_ggid
   
   ! These set of three subroutines are in charge of generating a dof_import for the distributed-memory solution of boundary mass matrices
   procedure        , non_overridable, private :: compute_boundary_dof_import                              => par_fe_space_compute_boundary_dof_import
   procedure        , non_overridable, private :: compute_ubound_boundary_num_itfc_couplings_by_continuity => pfs_compute_ubound_boundary_num_itfc_couplings_by_continuity
   procedure        , non_overridable, private :: compute_boundary_raw_interface_data_by_continuity        => par_fe_space_compute_boundary_raw_interface_data_by_continuity
   
   procedure        , non_overridable, private :: compute_faces_object                            => par_fe_space_compute_faces_object
   procedure        , non_overridable, private :: free_faces_object                               => par_fe_space_free_faces_object
   
   procedure        , non_overridable          :: get_num_fe_objects                              => par_fe_space_get_num_fe_objects
   
   procedure                                   :: print                                           => par_fe_space_print
   procedure                                   :: free                                            => par_fe_space_free
   procedure                                   :: create_dof_values                               => par_fe_space_create_dof_values
   procedure                                   :: interpolate_scalar_function                     => par_fe_space_interpolate_scalar_function
   procedure                                   :: interpolate_vector_function                     => par_fe_space_interpolate_vector_function  
   procedure                                   :: interpolate_dirichlet_values                    => par_fe_space_interpolate_dirichlet_values
   procedure                                   :: project_dirichlet_values_curl_conforming        => par_fe_space_project_dirichlet_values_curl_conforming
   
   procedure        , non_overridable, private :: setup_coarse_dofs                               => par_fe_space_setup_coarse_dofs
   procedure, nopass, non_overridable, private :: generate_coarse_dof_ggid                        => par_fe_space_generate_coarse_dof_ggid
   procedure        , non_overridable, private :: free_coarse_dofs                                => par_fe_space_free_coarse_dofs
   procedure        , non_overridable          :: setup_coarse_fe_space                           => par_fe_space_setup_coarse_fe_space
   procedure        , non_overridable, private :: transfer_num_fields                             => par_fe_space_transfer_num_fields
   procedure        , non_overridable, private :: transfer_fe_space_type                          => par_fe_space_transfer_fe_space_type
   procedure        , non_overridable, private :: gather_ptr_dofs_x_fe                            => par_fe_space_gather_ptr_dofs_x_fe
   procedure        , non_overridable, private :: gather_coarse_dofs_ggids_rcv_counts_and_displs  => par_fe_space_gather_coarse_dofs_ggids_rcv_counts_and_displs
   procedure        , non_overridable, private :: gather_coarse_dofs_ggids                        => par_fe_space_gather_coarse_dofs_ggids
   procedure        , non_overridable, private :: gather_vefs_ggids_dofs_objects                  => par_fe_space_gather_vefs_ggids_dofs_objects
   procedure                                   :: get_total_num_coarse_dofs                       => par_fe_space_get_total_num_coarse_dofs
   procedure                                   :: get_block_num_coarse_dofs                       => par_fe_space_get_block_num_coarse_dofs
   procedure       , non_overridable           :: get_coarse_fe_handler                           => par_fe_space_get_coarse_fe_handler
   
   ! Transfer and redistribution of FE functions
   procedure,                          private :: serial_fe_space_refine_and_coarsen_single_fe_function   => par_fe_space_refine_and_coarsen_single_fe_function
   procedure,                          private :: serial_fe_space_refine_and_coarsen_fe_function_array    => par_fe_space_refine_and_coarsen_fe_function_array
   procedure                                   :: redistribute                                            => par_fe_space_redistribute
   procedure                                   :: update_after_redistribute                               => par_fe_space_update_after_redistribute
   procedure,                          private :: update_after_refine_coarsen                             => par_fe_space_update_after_refine_coarsen
   procedure,                          private :: migrate_field_cell_to_ref_fes                           => par_fe_space_migrate_field_cell_to_ref_fes
   procedure,                          private :: migrate_fe_integration_arrays                           => par_fe_space_migrate_fe_integration_arrays
   procedure,                          private :: migrate_facet_integration_arrays                        => par_fe_space_migrate_facet_integration_arrays
   
   ! Objects-related traversals
   procedure, non_overridable                  :: create_fe_object_iterator                       => par_fe_space_create_fe_object_iterator
   procedure, non_overridable                  :: free_fe_object_iterator                         => par_fe_space_free_fe_object_iterator
   end type par_fe_space_t
 
 public :: par_fe_space_t
 public :: fe_object_iterator_t
 
  type, abstract :: l1_coarse_fe_handler_t
  contains
    ! Deferred methods
    procedure (l1_get_num_coarse_dofs_interface), deferred :: get_num_coarse_dofs
    procedure (l1_setup_constraint_matrix)      , deferred :: setup_constraint_matrix
    procedure (l1_setup_weighting_operator)     , deferred :: setup_weighting_operator
  end type l1_coarse_fe_handler_t
 
  abstract interface
    ! Returns the number of coarse DoFs that the object customizing
    ! l1_coarse_fe_handler_t requires to introduce on the subdomain 
    ! interface
    subroutine l1_get_num_coarse_dofs_interface(this, field_id, par_fe_space, parameter_list, num_coarse_dofs) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, ip
      class(l1_coarse_fe_handler_t), intent(in)    :: this
      integer(ip)                  , intent(in)    :: field_id
      type(par_fe_space_t)         , intent(in)    :: par_fe_space 
      type(parameterlist_t)        , intent(in)    :: parameter_list
      integer(ip)                  , intent(inout) :: num_coarse_dofs(:)
    end subroutine l1_get_num_coarse_dofs_interface
   
    subroutine l1_setup_constraint_matrix(this, field_id, par_fe_space, parameter_list, constraint_matrix) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, coo_sparse_matrix_t, ip
      class(l1_coarse_fe_handler_t), intent(in)    :: this
      integer(ip)                  , intent(in)    :: field_id
      type(par_fe_space_t)         , intent(in)    :: par_fe_space
      type(parameterlist_t)        , intent(in)    :: parameter_list
      type(coo_sparse_matrix_t)    , intent(inout) :: constraint_matrix
    end subroutine l1_setup_constraint_matrix
  
    subroutine l1_setup_weighting_operator(this, field_id, par_fe_space, parameter_list, weighting_operator) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, rp, ip
      class(l1_coarse_fe_handler_t) , intent(in)    :: this
      integer(ip)                  , intent(in)    :: field_id
      type(par_fe_space_t)          , intent(in)    :: par_fe_space
      type(parameterlist_t)         , intent(in)    :: parameter_list
      real(rp)         , allocatable, intent(inout) :: weighting_operator(:)
    end subroutine l1_setup_weighting_operator
  end interface

  type, extends(l1_coarse_fe_handler_t) :: standard_l1_coarse_fe_handler_t
    private
  contains
    procedure             :: get_num_coarse_dofs                       => standard_l1_get_num_coarse_dofs
    procedure             :: setup_constraint_matrix                   => standard_l1_setup_constraint_matrix
    procedure             :: setup_weighting_operator                  => standard_l1_setup_weighting_operator
    procedure, nopass     :: get_coarse_space_use_vertices_edges_faces => standard_get_coarse_space_use_vertices_edges_faces
  end type standard_l1_coarse_fe_handler_t
  
  type, extends(standard_l1_coarse_fe_handler_t) :: H1_l1_coarse_fe_handler_t
    private
    real(rp), public :: diffusion_inclusion = 1.0_rp
  contains
    procedure :: setup_constraint_matrix  => H1_l1_setup_constraint_matrix_multiple
    procedure :: setup_weighting_operator => H1_l1_setup_weighting_operator
  end type H1_l1_coarse_fe_handler_t
  
  type, extends(standard_l1_coarse_fe_handler_t) :: vector_laplacian_pb_bddc_l1_coarse_fe_handler_t
    private
    real(rp), public :: diffusion_inclusion = 1.0_rp
  contains
    procedure :: setup_constraint_matrix  => vector_laplacian_l1_setup_constraint_matrix_multiple
    procedure :: setup_weighting_operator => vector_laplacian_l1_setup_weighting_operator
  end type vector_laplacian_pb_bddc_l1_coarse_fe_handler_t

  type, extends(standard_l1_coarse_fe_handler_t) :: elasticity_pb_bddc_l1_coarse_fe_handler_t
    private
    real(rp), public :: elastic_modulus = 1.0_rp
  contains
    procedure :: get_num_coarse_dofs      => elasticity_l1_get_num_coarse_dofs  
    procedure :: setup_constraint_matrix  => elasticity_l1_setup_constraint_matrix_multiple
    procedure :: setup_weighting_operator => elasticity_l1_setup_weighting_operator
  end type elasticity_pb_bddc_l1_coarse_fe_handler_t
  
  public :: p_l1_coarse_fe_handler_t, l1_coarse_fe_handler_t, standard_l1_coarse_fe_handler_t, H1_l1_coarse_fe_handler_t, vector_laplacian_pb_bddc_l1_coarse_fe_handler_t, elasticity_pb_bddc_l1_coarse_fe_handler_t 
  
  type, extends(vector_function_t) :: rigid_body_mode_t
    private
    integer(ip)  :: num_dims = -1
    integer(ip)  :: imode = -1  
  contains
    procedure :: set_num_dims         => rigid_body_mode_set_num_dims
    procedure :: set_mode             => rigid_body_mode_set_mode
    procedure :: get_values_set_space => rigid_body_mode_get_values_set_space 
  end type rigid_body_mode_t
    
  type , extends(base_fe_cell_iterator_t) :: coarse_fe_cell_iterator_t
    private
    type(coarse_fe_space_t), pointer :: coarse_fe_space => NULL()
  contains
    procedure, non_overridable, private :: create                                      => coarse_fe_cell_iterator_create
    procedure, non_overridable, private :: free                                        => coarse_fe_cell_iterator_free
    final                               :: coarse_fe_cell_iterator_free_final
    procedure, non_overridable          :: create_own_dofs_on_vef_iterator             => coarse_fe_cell_iterator_create_own_dofs_on_vef_iterator
    procedure, non_overridable, private :: generate_own_dofs_vef                       => coarse_fe_cell_iterator_generate_own_dofs_vef
    procedure, non_overridable, private :: generate_own_dofs_vef_from_source_coarse_fe => cfeci_generate_own_dofs_vef_from_source_coarse_fe
    procedure, non_overridable, private :: renum_dofs_block                            => coarse_fe_cell_iterator_renum_dofs_block
    procedure, non_overridable, private :: renum_dofs_field                            => coarse_fe_cell_iterator_renum_dofs_field
    
    procedure, non_overridable, private :: get_scan_sum_num_dofs                       => coarse_fe_cell_iterator_get_scan_sum_num_dofs
    procedure, non_overridable          :: get_num_fe_spaces                           => coarse_fe_cell_iterator_get_num_fe_spaces
    procedure, non_overridable          :: get_num_dofs                                => coarse_fe_cell_iterator_get_num_dofs
    procedure, non_overridable          :: get_fe_dofs                                 => coarse_fe_cell_iterator_get_fe_dofs
    procedure, non_overridable          :: get_field_fe_dofs                           => coarse_fe_cell_iterator_get_field_fe_dofs
    procedure, non_overridable          :: get_coarse_fe_vef                           => coarse_fe_cell_iterator_get_coarse_fe_vef
    
    procedure, non_overridable          :: last                                        => coarse_fe_cell_iterator_last
    procedure, non_overridable          :: scan_sum_num_vefs                           => coarse_fe_cell_iterator_scan_sum_num_vefs
  end type coarse_fe_cell_iterator_t
  
  type, extends(base_fe_vef_iterator_t) :: coarse_fe_vef_iterator_t
    private
    type(coarse_fe_space_t), pointer    :: coarse_fe_space => NULL()
  contains
     procedure, non_overridable, private :: create                    => coarse_fe_vef_iterator_create
     procedure, non_overridable, private :: free                      => coarse_fe_vef_iterator_free
     final                               :: coarse_fe_vef_iterator_final
     procedure, non_overridable          :: get_coarse_fe_around      => coarse_fe_vef_iterator_get_coarse_fe_around
     generic                             :: get_cell_around           => get_coarse_fe_around
  end type coarse_fe_vef_iterator_t
  
  type, extends(base_fe_object_iterator_t) :: coarse_fe_object_iterator_t
    private
    type(coarse_fe_space_t), pointer :: coarse_fe_space => NULL()
  contains
    procedure, non_overridable, private  :: create                                   => coarse_fe_object_iterator_create
    procedure, non_overridable, private  :: free                                     => coarse_fe_object_iterator_free
    
    ! Own methods of coarse_fe_object_iterator_t
    procedure, non_overridable, private  :: coarse_fe_object_iterator_get_fe_vef
    generic                              :: get_vef                                  => coarse_fe_object_iterator_get_fe_vef 
    procedure, non_overridable           :: get_num_coarse_dofs                   => coarse_fe_object_iterator_get_num_coarse_dofs
    procedure, non_overridable           :: create_own_coarse_dofs_iterator          => coarse_fe_object_iterator_create_own_coarse_dofs_iterator
  end type coarse_fe_object_iterator_t
    
  
  type, extends(base_fe_space_t) :: coarse_fe_space_t
    private
    
    ! Data related to block structure of the FE system + size of each block
    integer(ip)                                 :: num_blocks
    integer(ip)                   , allocatable :: field_blocks(:)
    logical                       , allocatable :: field_coupling(:,:)
    logical                       , allocatable :: blocks_coupling(:,:)
    integer(ip)                   , allocatable :: num_dofs_x_block(:)
    
    integer(ip) , allocatable                   :: ptr_dofs_x_fe(:)
    integer(ip) , allocatable                   :: lst_dofs_gids(:)
    type(list_t), allocatable                   :: own_dofs_vef_x_fe(:)
     
    ! Pointer to coarse triangulation this coarse_fe_space has been built from
    type(coarse_triangulation_t), pointer       :: coarse_triangulation => NULL()
 
    ! It is the equivalent to the "element_to_dof" at the finer level
    ! Pointers to the start/end of coarse DoFs gids of each field (lst_coarse_dofs)
    integer(ip)                   , allocatable :: ptr_coarse_dofs_x_field(:) 
    ! List of coarse DoFs gids
    integer(ip)                   , allocatable :: lst_coarse_dofs(:) 
    
    ! Coarse DoFs gids on top of coarse n_faces per field
    ! It provides (for every field) what we get from the reference_fe_t at the
    ! finer level
    type(list_t), allocatable                   :: own_coarse_dofs_x_field(:)
   
    
    ! GIDs of the coarse DoFs. For their generation, we though to be a good idea 
    ! to take as basis the GIDs of the VEFs objects they are built from (instead of
    ! generating them from scratch, which in turn would imply further communication).
    type(allocatable_array_igp1_t), allocatable :: coarse_dof_ggids_x_field(:)
  contains
    procedure                                   :: create                                          => coarse_fe_space_create
    procedure                                   :: free                                            => coarse_fe_space_free
    procedure, non_overridable                  :: print                                           => coarse_fe_space_print
    procedure, non_overridable, private         :: allocate_and_fill_field_blocks_and_coupling     => coarse_fe_space_allocate_and_fill_field_blocks_and_coupling
    procedure, non_overridable, private         :: free_field_blocks_and_coupling                  => coarse_fe_space_free_field_blocks_and_coupling
    procedure, non_overridable, private         :: allocate_and_fill_fe_space_type_x_field         => coarse_fe_space_allocate_and_fill_fe_space_type
    procedure, non_overridable, private         :: free_fe_space_type_x_field                      => coarse_fe_space_free_fe_space_type
    procedure, non_overridable, private         :: allocate_and_fill_ptr_dofs_x_fe                 => coarse_fe_space_allocate_and_fill_ptr_dofs_x_fe
    procedure, non_overridable, private         :: free_ptr_dofs_x_fe                              => coarse_fe_space_free_ptr_dofs_x_fe
    procedure, non_overridable, private         :: fetch_ghost_fes_data                            => coarse_fe_space_fetch_ghost_fes_data
    procedure, non_overridable, private         :: allocate_and_generate_own_dofs_vef_x_fe         => coarse_fe_space_allocate_and_generate_own_dofs_vef_x_fe
    procedure, non_overridable, private         :: generate_own_dofs_cell_x_fe_field               => coarse_fe_space_generate_own_dofs_cell_x_fe_field
    procedure, non_overridable, private         :: free_own_dofs_vef_x_fe                          => coarse_fe_space_free_free_own_dofs_vef_x_fe
    procedure, non_overridable, private         :: allocate_lst_dofs_gids                          => coarse_fe_space_allocate_lst_dofs_gids
    procedure, non_overridable, private         :: count_dofs_and_fill_lst_dof_gids                => coarse_fe_space_count_dofs_and_fill_lst_dof_gids
    procedure, non_overridable, private         :: count_dofs_and_fill_lst_dof_gids_field          => coarse_fe_space_count_dofs_and_fill_lst_dof_gids_field
    procedure, non_overridable, private         :: free_lst_dofs_gids                              => coarse_fe_space_free_lst_dofs_gids
    procedure, non_overridable, private         :: free_num_dofs_x_field_and_block                 => coarse_fe_space_free_num_dofs_x_field_and_block
    procedure, non_overridable, nopass, private :: coarse_fe_size                                  => coarse_fe_space_coarse_fe_size
    procedure, non_overridable, nopass, private :: coarse_fe_pack                                  => coarse_fe_space_coarse_fe_pack
    procedure, non_overridable, nopass, private :: coarse_fe_unpack                                => coarse_fe_space_coarse_fe_unpack
    procedure, non_overridable, private         :: compute_blocks_dof_import                       => coarse_fe_space_compute_blocks_dof_import
    procedure, non_overridable, private         :: compute_dof_import                              => coarse_fe_space_compute_dof_import
    procedure, non_overridable, private         :: compute_ubound_num_itfc_couplings_by_continuity => coarse_fe_space_compute_ubound_num_itfc_couplings_by_continuity
    procedure, non_overridable, private         :: compute_raw_interface_data_by_continuity        => coarse_fe_space_compute_raw_interface_data_by_continuity
    procedure, non_overridable, private         :: raw_interface_data_by_continuity_decide_owner   => coarse_fe_space_raw_interface_data_by_continuity_decide_owner
    
    procedure, non_overridable                  :: setup_coarse_dofs                               => coarse_fe_space_setup_coarse_dofs
    procedure, non_overridable, private         :: setup_num_coarse_dofs                           => coarse_fe_space_setup_num_coarse_dofs
    procedure, non_overridable                  :: setup_constraint_matrix                         => coarse_fe_space_setup_constraint_matrix
    procedure, non_overridable, private         :: setup_coarse_fe_space                           => coarse_fe_space_setup_coarse_fe_space
    procedure, non_overridable, private         :: transfer_num_fields                          => coarse_fe_space_transfer_num_fields
    procedure, non_overridable, private         :: transfer_fe_space_type                          => coarse_fe_space_transfer_fe_space_type
    procedure, non_overridable, private         :: gather_ptr_dofs_x_fe                => coarse_fe_space_gather_ptr_dofs_x_fe
    procedure, non_overridable, private         :: gather_coarse_dofs_gids_rcv_counts_and_displs   => coarse_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
    procedure, non_overridable, private         :: gather_coarse_dofs_gids                         => coarse_fe_space_gather_coarse_dofs_gids
    procedure, non_overridable, private         :: gather_vefs_gids_dofs_objects                   => coarse_fe_space_gather_vefs_gids_dofs_objects

    procedure                                   :: get_num_fe_objects                           => coarse_fe_space_get_num_fe_objects
    procedure                                   :: get_total_num_coarse_dofs                    => coarse_fe_space_get_total_num_coarse_dofs
    procedure                                   :: get_block_num_coarse_dofs                    => coarse_fe_space_get_block_num_coarse_dofs
    procedure, non_overridable, private         :: free_coarse_dofs                                => coarse_fe_space_free_coarse_dofs
    
    procedure                                   :: renum_dofs_first_interior_then_interface     => coarse_fe_space_renum_dofs_first_interior_then_interface
    procedure                         , private :: renum_dofs_block                             => coarse_fe_space_renum_dofs_block
    
    ! Coarse FE traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_cell_iterator                       => coarse_fe_space_create_coarse_fe_cell_iterator
    procedure, non_overridable                  :: free_coarse_fe_cell_iterator                         => coarse_fe_space_free_coarse_fe_cell_iterator

    
    ! Coarse FE VEFS traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_vef_iterator                   => coarse_fe_space_create_coarse_fe_vef_iterator
    procedure, non_overridable                  :: create_itfc_coarse_fe_vef_iterator              => coarse_fe_space_create_itfc_coarse_fe_vef_iterator
    procedure, non_overridable                  :: free_coarse_fe_vef_iterator                     => coarse_fe_space_free_coarse_fe_vef_iterator
    
    
     ! Objects-related traversals
     procedure, non_overridable                 :: create_coarse_fe_object_iterator                => coarse_fe_space_create_coarse_fe_object_iterator
     procedure, non_overridable                 :: free_coarse_fe_object_iterator                  => coarse_fe_space_free_coarse_fe_object_iterator
     
     ! Getters
     procedure, non_overridable                 :: get_num_local_coarse_fes                     => coarse_fe_space_get_num_local_coarse_fes
     procedure, non_overridable                 :: get_num_ghost_coarse_fes                     => coarse_fe_space_get_num_ghost_coarse_fes
     procedure, non_overridable                 :: get_num_coarse_fe_objects                    => coarse_fe_space_get_num_coarse_fe_objects
     procedure, non_overridable                 :: get_triangulation                            => coarse_fe_space_get_triangulation
     procedure                                  :: get_environment                              => coarse_fe_space_get_environment
     
     procedure                                  :: get_num_blocks                               => coarse_fe_space_get_num_blocks
     procedure                                  :: get_field_blocks                             => coarse_fe_space_get_field_blocks
     procedure                                  :: get_field_coupling                           => coarse_fe_space_get_field_coupling
     procedure                                  :: get_block_num_dofs                           => coarse_fe_space_get_block_num_dofs
 end type coarse_fe_space_t
 
 public :: coarse_fe_space_t, coarse_fe_cell_iterator_t, coarse_fe_vef_iterator_t
 public :: coarse_fe_object_iterator_t
 public :: coarse_space_use_vertices_key, coarse_space_use_edges_key, coarse_space_use_faces_key
 
 type, abstract :: lgt1_coarse_fe_handler_t
  contains
    ! Deferred methods
    procedure (lgt1_setup_coarse_dofs_interface), deferred :: setup_coarse_dofs
  end type lgt1_coarse_fe_handler_t

  abstract interface
    subroutine lgt1_setup_coarse_dofs_interface(this, coarse_fe_space) 
      import :: lgt1_coarse_fe_handler_t, coarse_fe_space_t
      implicit none
      class(lgt1_coarse_fe_handler_t), intent(in)    :: this
      type(coarse_fe_space_t)        , intent(inout) :: coarse_fe_space 
    end subroutine lgt1_setup_coarse_dofs_interface
  end interface
  
  type, extends(lgt1_coarse_fe_handler_t) :: standard_lgt1_coarse_fe_handler_t
    private
  contains
    procedure :: setup_coarse_dofs       => standard_lgt1_setup_coarse_dofs
 end type standard_lgt1_coarse_fe_handler_t
 
  type fe_function_t
   private
   class(vector_t), allocatable  :: free_dof_values
   type(serial_scalar_array_t)   :: fixed_dof_values
   type(serial_scalar_array_t)   :: constraining_x_fixed_dof_values ! C_D u_D
  contains
     procedure, non_overridable          :: create                         => fe_function_create
     procedure, non_overridable          :: gather_nodal_values_through_iterator => fe_function_gather_nodal_values_through_iterator
     procedure, non_overridable          :: gather_nodal_values_from_raw_data    => fe_function_gather_nodal_values_from_raw_data
     generic                             :: gather_nodal_values                  => gather_nodal_values_through_iterator, &
                                                                                    gather_nodal_values_from_raw_data
     procedure, non_overridable          :: insert_nodal_values            => fe_function_insert_nodal_values
     procedure, non_overridable          :: axpby                          => fe_function_axpby
     procedure, non_overridable          :: copy                           => fe_function_copy
     procedure, non_overridable          :: get_free_dof_values            => fe_function_get_free_dof_values
     procedure, non_overridable          :: get_fixed_dof_values           => fe_function_get_fixed_dof_values
     procedure, non_overridable          :: set_free_dof_values            => fe_function_set_free_dof_values
     procedure, non_overridable          :: set_fixed_dof_values           => fe_function_set_fixed_dof_values
     procedure, non_overridable          :: free                           => fe_function_free
     generic                             :: assignment(=)                  => copy
  end type fe_function_t 
  
  type :: p_fe_function_t
    class(fe_function_t), pointer :: p
  end type p_fe_function_t
  
  public :: fe_function_t, p_fe_function_t

 character(*), parameter :: interpolator_type_nodal = "interpolator_type_nodal"
 character(*), parameter :: interpolator_type_Hcurl = "interpolator_type_Hcurl"
 
 ! Abstract interpolator
 type, abstract :: interpolator_t
    private
    integer(ip) :: field_id
  contains    
    procedure, private  :: get_function_values_from_scalar_components => interpolator_t_get_function_values_from_scalar_components
    procedure, private  :: get_vector_function_values                 => interpolator_t_get_vector_function_values 
    generic             :: get_function_values                        => get_vector_function_values 
    procedure, private  :: reallocate_vector_function_values          => interpolator_t_reallocate_vector_function_values 
    procedure, private  :: reallocate_scalar_function_values          => interpolator_t_reallocate_scalar_function_values
    procedure, private  :: reallocate_array                           => interpolator_t_reallocate_array 
    generic             :: reallocate_if_needed => reallocate_vector_function_values, reallocate_scalar_function_values, reallocate_array 
    ! Deferred TBPs 
    procedure(interpolator_create_interface)                               , deferred :: create
    procedure(interpolator_evaluate_scalar_function_moments_interface)     , deferred :: evaluate_scalar_function_moments
    procedure(interpolator_evaluate_vector_function_moments_interface)     , deferred :: evaluate_vector_function_moments 
    procedure(interpolator_evaluate_function_components_moments_interface) , deferred :: evaluate_function_scalar_components_moments
    procedure(interpolator_free_interface)                                 , deferred :: free    
 end type interpolator_t

 type p_interpolator_t
    class(interpolator_t), allocatable :: p 
  contains
 end type p_interpolator_t

 abstract interface 

    subroutine interpolator_create_interface( this, fe_space, field_id )
      import :: interpolator_t, reference_fe_t, ip, serial_fe_space_t
      implicit none
      class(interpolator_t)    , intent(inout) :: this
      class(serial_fe_space_t) , intent(in)    :: fe_space
      integer(ip)              , intent(in)    :: field_id
    end subroutine interpolator_create_interface

    subroutine interpolator_evaluate_scalar_function_moments_interface( this, fe, scalar_function, dof_values, n_face_mask, time )
      import :: interpolator_t, scalar_function_t, fe_cell_iterator_t, rp   
      implicit none
      class(interpolator_t)           , intent(inout) :: this
      class(fe_cell_iterator_t)       , intent(in)    :: fe
      class(scalar_function_t)        , intent(in)    :: scalar_function
      real(rp) , allocatable          , intent(inout) :: dof_values(:) 
      logical  , optional             , intent(in)    :: n_face_mask(:)
      real(rp) , optional             , intent(in)    :: time 
    end subroutine interpolator_evaluate_scalar_function_moments_interface

    subroutine interpolator_evaluate_vector_function_moments_interface( this, fe, vector_function, dof_values, n_face_mask, time )
      import :: interpolator_t, vector_function_t, fe_cell_iterator_t, rp   
      implicit none
      class(interpolator_t)           , intent(inout) :: this
      class(fe_cell_iterator_t)       , intent(in)    :: fe
      class(vector_function_t)        , intent(in)    :: vector_function
      real(rp) , allocatable          , intent(inout) :: dof_values(:) 
      logical  , optional             , intent(in)    :: n_face_mask(:)
      real(rp) , optional             , intent(in)    :: time 
    end subroutine interpolator_evaluate_vector_function_moments_interface

    subroutine interpolator_evaluate_function_components_moments_interface( this, n_face_mask, fe, vector_function_scalar_components, dof_values, time )
      import :: interpolator_t, fe_cell_iterator_t, p_scalar_function_t, rp, ip 
      implicit none 
      class(interpolator_t)           , intent(inout) :: this
      logical                         , intent(in)    :: n_face_mask(:) 
      class(fe_cell_iterator_t)       , intent(in)    :: fe
      class(p_scalar_function_t)      , intent(in)    :: vector_function_scalar_components(:,:)
      real(rp) , allocatable          , intent(inout) :: dof_values(:) 
      real(rp) , optional             , intent(in)    :: time 
    end subroutine interpolator_evaluate_function_components_moments_interface

    subroutine interpolator_free_interface( this )
      import :: interpolator_t, cell_map_t  
      implicit none
      class(interpolator_t)           , intent(inout) :: this
    end subroutine interpolator_free_interface
 end interface

 type, extends(interpolator_t) :: nodal_interpolator_t 
 private 
 type(cell_map_t)     , allocatable  :: cell_maps(:) 
 type(p_quadrature_t) , allocatable  :: nodal_quadratures(:)   
 ! Boundary values array 
 real(rp), allocatable             :: scalar_function_values(:,:)
 type(vector_field_t), allocatable :: function_values(:,:)
contains 
 procedure :: create                                             => nodal_interpolator_create
 procedure :: evaluate_scalar_function_moments                   => nodal_interpolator_evaluate_scalar_function_moments 
 procedure :: evaluate_vector_function_moments                   => nodal_interpolator_evaluate_vector_function_moments  
 procedure :: evaluate_function_scalar_components_moments        => nodal_interpolator_evaluate_function_scalar_components_moments
 procedure :: free                                               => nodal_interpolator_free
end type nodal_interpolator_t

type, extends(interpolator_t), abstract :: Hcurl_interpolator_t 
private 
! Maps 
type(edge_map_t)  , allocatable  :: edge_maps(:) 
type(facet_map_t) , allocatable  :: facet_maps(:) 
type(cell_map_t)  , allocatable  :: cell_maps(:) 
! Quadratures 
type(quadrature_t)  , allocatable  :: edge_quadratures(:) 
type(quadrature_t)  , allocatable  :: facet_quadratures(:) 
type(quadrature_t)  , allocatable  :: cell_quadratures(:)  
! Interpolation 
type(interpolation_t) , allocatable :: edge_interpolations(:)
type(interpolation_t) , allocatable :: facet_interpolations(:) 
type(interpolation_t) , allocatable :: cell_interpolations(:)  
! Function values arrays 
type(vector_field_t), allocatable :: edge_function_values(:,:) 
type(vector_field_t), allocatable :: facet_function_values(:,:) 
type(vector_field_t), allocatable :: cell_function_values(:,:)
! Boundary values array 
real(rp), allocatable             :: scalar_function_values_on_edge(:,:)
real(rp), allocatable             :: scalar_function_values_on_facet(:,:)
contains   
procedure :: evaluate_scalar_function_moments    => Hcurl_interpolator_evaluate_scalar_function_moments 
procedure, private :: update_edge_map_coordinates    => Hcurl_interpolator_update_edge_map_coordinates
procedure, private :: update_facet_map_coordinates   => Hcurl_interpolator_update_facet_map_coordinates
generic            :: update_map_coordinates         => update_edge_map_coordinates, update_facet_map_coordinates
end type Hcurl_interpolator_t 

type, extends(Hcurl_interpolator_t) :: hex_Hcurl_interpolator_t 
private       
type(hex_lagrangian_reference_fe_t)       , allocatable   :: fes_1D(:) 
type(hex_nedelec_reference_fe_t)          , allocatable   :: fes_2D(:) 
type(hex_raviart_thomas_reference_fe_t )  , allocatable   :: fes_rt(:) 
type(interpolation_t)                     , allocatable   :: real_cell_interpolations(:) 
contains 
procedure :: create                                             => hex_Hcurl_interpolator_create
procedure :: evaluate_vector_function_moments                   => hex_Hcurl_interpolator_evaluate_vector_function_moments  
procedure :: evaluate_function_scalar_components_moments        => hex_Hcurl_interpolator_evaluate_function_components_moments
procedure :: free                                               => hex_Hcurl_interpolator_free
end type hex_Hcurl_interpolator_t

type, extends(Hcurl_interpolator_t) :: tet_Hcurl_interpolator_t 
private       
type(tet_lagrangian_reference_fe_t)    , allocatable      :: fes_1D(:)
type(tet_lagrangian_reference_fe_t)    , allocatable      :: fes_2D(:) 
type(tet_lagrangian_reference_fe_t )   , allocatable      :: fes_lagrangian(:) 
contains 
procedure :: create                                             => tet_Hcurl_interpolator_create
procedure :: evaluate_vector_function_moments                   => tet_Hcurl_interpolator_evaluate_vector_function_moments  
procedure :: evaluate_function_scalar_components_moments        => tet_Hcurl_interpolator_evaluate_function_components_moments
procedure :: free                                               => tet_Hcurl_interpolator_free
end type tet_Hcurl_interpolator_t

contains
!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_base_fe_space.i90"
#include "sbm_base_fe_cell_iterator.i90"
#include "sbm_serial_fe_space.i90"
#include "sbm_fe_cell_iterator.i90"
#include "sbm_fe_facet_iterator.i90"
#include "sbm_base_fe_vef_iterator.i90"
#include "sbm_fe_vef_iterator.i90"
#include "sbm_par_fe_space.i90"
#include "sbm_base_fe_object_iterator.i90"
#include "sbm_fe_object_iterator.i90"
#include "sbm_standard_coarse_fe_handler.i90"
#include "sbm_H1_coarse_fe_handler.i90"
#include "sbm_vector_laplacian_coarse_fe_handler.i90"
#include "sbm_elasticity_coarse_fe_handler.i90"

#include "sbm_coarse_fe_space.i90"
#include "sbm_coarse_fe_object_iterator.i90"
#include "sbm_coarse_fe_cell_iterator.i90"
#include "sbm_coarse_fe_vef_iterator.i90"

#include "sbm_fe_function.i90"

#include "sbm_interpolator.i90"
#include "sbm_nodal_interpolator.i90"
#include "sbm_Hcurl_interpolator.i90"
#include "sbm_hex_Hcurl_interpolator.i90" 
#include "sbm_tet_Hcurl_interpolator.i90"

end module fe_space_names
