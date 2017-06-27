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
  use list_types_names
  use memor_names
  use sort_names
  use allocatable_array_names
  use hash_table_names
  use FPL

  use environment_names
  use base_static_triangulation_names
  use conditions_names
  
  use reference_fe_names
  use field_names
  use function_names
  
  use operator_names
  use matrix_names
  use vector_names
  use array_names
  use matrix_array_assembler_names
  use base_sparse_matrix_names
  use sparse_matrix_names
  use block_sparse_matrix_names
  use serial_scalar_array_names
  use serial_block_array_names
  use sparse_matrix_array_assembler_names
  use block_sparse_matrix_array_assembler_names
  use direct_solver_names
  use direct_solver_parameters_names
  use iterative_linear_solver_names
  use iterative_linear_solver_parameters_names
  
  
 ! Parallel modules
  use environment_names
  !use par_context_names
  use dof_import_names
  use cell_import_names
  use par_sparse_matrix_names
  use par_scalar_array_names
  use par_sparse_matrix_array_assembler_names
  
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
    integer(ip)                                 :: number_fields
    integer(ip) , allocatable                   :: fe_space_type_per_field(:)
    
    ! Data related to block structure of the FE system + size of each block
    integer(ip)                                 :: number_blocks
    integer(ip)                   , allocatable :: field_blocks(:)
    logical                       , allocatable :: field_coupling(:,:)
    logical                       , allocatable :: blocks_coupling(:,:)
    integer(ip)                   , allocatable :: number_dofs_per_block(:) 
    integer(ip) , allocatable                   :: number_dofs_per_field(:)
    
    type(dof_import_t)            , allocatable :: blocks_dof_import(:)
    
    ! Pointer to data structure which is in charge of coarse DoF handling.
    ! It will be a nullified pointer on L1 tasks, and associated via target 
    ! allocation in the case of L2-Ln tasks.
    type(coarse_fe_space_t)       , pointer :: coarse_fe_space => NULL()
  contains 
    procedure, non_overridable                 :: get_number_fields                               => base_fe_space_get_number_fields
    procedure, non_overridable                 :: get_fe_space_type                               => base_fe_space_get_fe_space_type
    procedure, non_overridable                 :: get_number_blocks                               => base_fe_space_get_number_blocks
    procedure, non_overridable                 :: get_field_blocks                                => base_fe_space_get_field_blocks
    procedure, non_overridable                 :: get_field_coupling                              => base_fe_space_get_field_coupling
    procedure, non_overridable                 :: get_total_number_dofs                           => base_fe_space_get_total_number_dofs
    procedure, non_overridable                 :: get_field_number_dofs                           => base_fe_space_get_field_number_dofs
    procedure, non_overridable                 :: get_block_number_dofs                           => base_fe_space_get_block_number_dofs
    procedure, non_overridable                 :: get_total_number_interior_dofs                  => base_fe_space_get_total_number_interior_dofs
    procedure, non_overridable                 :: get_total_number_interface_dofs                 => base_fe_space_get_total_number_interface_dofs
    procedure, non_overridable                 :: get_block_number_interior_dofs                  => base_fe_space_get_block_number_interior_dofs
    procedure, non_overridable                 :: get_block_number_interface_dofs                 => base_fe_space_get_block_number_interface_dofs
    procedure, non_overridable                 :: get_block_dof_import                            => base_fe_space_get_block_dof_import
    procedure, non_overridable                 :: get_coarse_fe_space                             => base_fe_space_get_coarse_fe_space
    
    ! These three TBPs are (currently) though to be overriden by subclasses.
    ! In the future, once we harmonize type(coarse_fe_space_t) and type(serial/par_fe_space_t),
    ! the code of these TBPs should be developed at this level, being valid for all subclasses.
    procedure                                  :: get_environment                                 => base_fe_space_get_environment
    procedure                                  :: get_total_number_coarse_dofs                    => base_fe_space_get_total_number_coarse_dofs
    procedure                                  :: get_block_number_coarse_dofs                    => base_fe_space_get_block_number_coarse_dofs
  end type base_fe_space_t
  
  public :: base_fe_space_t

  type :: base_fe_iterator_t
    private
    class(cell_iterator_t), allocatable :: cell
  contains
    ! Methods exploiting ("inherited from") cell_t common to all descendants
    procedure, non_overridable           :: next                    => base_fe_iterator_next
    procedure, non_overridable           :: first                   => base_fe_iterator_first
    procedure, non_overridable           :: last                    => base_fe_iterator_last
    procedure, non_overridable           :: set_lid                 => base_fe_iterator_set_lid
    procedure, non_overridable           :: has_finished            => base_fe_iterator_has_finished
    procedure, non_overridable           :: get_reference_fe_geo    => base_fe_iterator_get_reference_fe_geo
    procedure, non_overridable           :: get_reference_fe_geo_id => base_fe_iterator_get_reference_fe_geo_id
    procedure, non_overridable           :: get_coordinates         => base_fe_iterator_get_coordinates
    procedure, non_overridable           :: set_coordinates         => base_fe_iterator_set_coordinates
    procedure, non_overridable           :: get_lid                 => base_fe_iterator_get_lid
    procedure, non_overridable           :: get_gid                 => base_fe_iterator_get_gid
    procedure, non_overridable           :: get_my_part             => base_fe_iterator_get_mypart
    procedure, non_overridable           :: get_my_subpart          => base_fe_iterator_get_mysubpart
    procedure, non_overridable           :: get_my_subpart_lid      => base_fe_iterator_get_mysubpart_lid
    procedure, non_overridable           :: get_set_id              => base_fe_iterator_get_set_id
    procedure, non_overridable           :: get_num_vefs            => base_fe_iterator_get_num_vefs
    procedure, non_overridable           :: get_vef_lid             => base_fe_iterator_get_vef_lid
    procedure, non_overridable           :: get_vef_lids            => base_fe_iterator_get_vef_lids
    procedure, non_overridable           :: get_vef_gid             => base_fe_iterator_get_vef_gid
    procedure, non_overridable           :: find_lpos_vef_lid       => base_fe_iterator_find_lpos_vef_lid
    procedure, non_overridable           :: find_lpos_vef_gid       => base_fe_iterator_find_lpos_vef_gid
    procedure, non_overridable           :: is_local                => base_fe_iterator_is_local
    procedure, non_overridable           :: is_ghost                => base_fe_iterator_is_ghost
    procedure, non_overridable           :: scan_sum_number_vefs    => base_fe_iterator_get_scan_sum_number_vefs
    procedure, non_overridable, private  :: base_fe_iterator_get_vef
    generic                              :: get_vef                 => base_fe_iterator_get_vef
  end type base_fe_iterator_t
  
  
  type, extends(base_fe_iterator_t) :: fe_iterator_t
    private
    class(serial_fe_space_t), pointer    :: fe_space => NULL()
  contains
    procedure                 , private :: create                                     => fe_iterator_create
    procedure                 , private :: free                                       => fe_iterator_free
    final                               :: fe_iterator_free_final
    procedure, non_overridable, private :: fill_own_dofs                              => fe_iterator_fill_own_dofs
    procedure, non_overridable, private :: fill_own_dofs_on_vef                       => fe_iterator_fill_own_dofs_on_vef
    procedure, non_overridable, private :: fill_own_dofs_on_vef_component_wise        => fe_iterator_fill_own_dofs_on_vef_component_wise
    procedure, non_overridable, private :: fill_own_dofs_on_vef_from_source_fe        => fe_iterator_fill_own_dofs_on_vef_from_source_fe
    procedure, non_overridable, private :: fill_dofs_face_integration_coupling        => fe_iterator_fill_dofs_face_integration_coupling
    procedure, non_overridable, private :: renumber_dofs_block                        => fe_iterator_renumber_dofs_block
    procedure, non_overridable, private :: renumber_dofs_field                        => fe_iterator_renumber_dofs_field
    procedure, non_overridable          :: update_integration                         => fe_iterator_update_integration

    procedure, non_overridable          :: get_fe_space                               => fe_iterator_get_fe_space
    procedure, non_overridable          :: get_number_fields                          => fe_iterator_get_number_fields
    procedure, non_overridable, private :: get_fe_space_type                          => fe_iterator_get_fe_space_type
    procedure, non_overridable          :: get_field_type                             => fe_iterator_get_field_type

    procedure, non_overridable          :: get_field_blocks                           => fe_iterator_get_field_blocks
    procedure, non_overridable          :: get_number_dofs                            => fe_iterator_get_number_dofs
    procedure, non_overridable          :: get_number_dofs_per_field                  => fe_iterator_get_number_dofs_per_field
    procedure, non_overridable          :: get_field_elem2dof                         => fe_iterator_get_field_elem2dof
    procedure, non_overridable          :: get_elem2dof                               => fe_iterator_get_elem2dof
    procedure, non_overridable          :: get_order                                  => fe_iterator_get_order
    
    procedure, non_overridable          :: get_max_order_single_field                 => fe_iterator_get_max_order_single_field
    procedure, non_overridable          :: get_max_order_all_fields                   => fe_iterator_get_max_order_all_fields
    generic                             :: get_max_order                              => get_max_order_single_field, &
                                                                                         get_max_order_all_fields

    procedure, non_overridable          :: at_strong_dirichlet_boundary               => fe_iterator_at_strong_dirichlet_boundary
    procedure, non_overridable          :: fe_iterator_determine_at_strong_dirichlet_boundary_single_field
    procedure, non_overridable          :: fe_iterator_determine_at_strong_dirichlet_boundary_all_fields
    generic                             :: determine_at_strong_dirichlet_boundary     => fe_iterator_determine_at_strong_dirichlet_boundary_single_field, &
                                                                                         fe_iterator_determine_at_strong_dirichlet_boundary_all_fields 
    procedure, non_overridable          :: compute_volume                             => fe_iterator_compute_volume
    
    procedure, non_overridable          :: get_quadrature                             => fe_iterator_get_quadrature
    procedure, non_overridable          :: get_fe_map                                 => fe_iterator_get_fe_map
    procedure, non_overridable          :: get_volume_integrator                      => fe_iterator_get_volume_integrator    
    
    procedure, non_overridable, private :: fe_iterator_get_fe_vef
    generic                             :: get_vef                                    => fe_iterator_get_fe_vef
    procedure, non_overridable          :: get_reference_fe                           => fe_iterator_get_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe                 => fe_iterator_get_max_order_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe_id              => fe_iterator_get_max_order_reference_fe_id
    procedure, non_overridable          :: get_reference_fe_id                        => fe_iterator_get_reference_fe_id
    procedure, non_overridable          :: create_own_dofs_on_vef_iterator            => fe_iterator_create_own_dofs_on_vef_iterator
    procedure, non_overridable          :: impose_strong_dirichlet_bcs                => fe_iterator_impose_strong_dirichlet_bcs
  end type fe_iterator_t
   
  type :: base_fe_vef_iterator_t
    private
    class(vef_iterator_t), allocatable :: vef
  contains
     ! Methods exploiting ("inherited from") vef_t
     procedure                           :: first                     => base_fe_vef_iterator_first
     procedure                           :: next                      => base_fe_vef_iterator_next
     procedure, non_overridable          :: set_lid                   => base_fe_vef_iterator_set_lid
     procedure, non_overridable          :: has_finished              => base_fe_vef_iterator_has_finished
     procedure, non_overridable          :: get_lid                   => base_fe_vef_iterator_get_lid
     procedure, non_overridable          :: get_set_id                => base_fe_vef_iterator_get_set_id

     procedure, non_overridable          :: get_dimension             => base_fe_vef_iterator_get_dimension
     procedure, non_overridable          :: is_at_boundary            => base_fe_vef_iterator_is_at_boundary
     procedure, non_overridable          :: is_local                  => base_fe_vef_iterator_is_local
     procedure, non_overridable          :: is_ghost                  => base_fe_vef_iterator_is_ghost
     procedure, non_overridable          :: is_at_interface           => base_fe_vef_iterator_is_at_interface
     procedure, non_overridable          :: is_face                   => base_fe_vef_iterator_is_face
     
     procedure, non_overridable          :: get_num_cells_around      => base_fe_vef_iterator_get_num_cells_around
     procedure, non_overridable          :: base_fe_vef_iterator_get_cell_around
     generic                             :: get_cell_around           => base_fe_vef_iterator_get_cell_around
  end type base_fe_vef_iterator_t

  type , extends(base_fe_vef_iterator_t) :: fe_vef_iterator_t
    private
    class(serial_fe_space_t), pointer :: fe_space => NULL()
  contains
     procedure                 , private :: create                    => fe_vef_iterator_create
     procedure                 , private :: free                      => fe_vef_iterator_free
     final                               :: fe_vef_iterator_final
     procedure, non_overridable, private :: fe_vef_iterator_get_fe_around
     generic                             :: get_cell_around           => fe_vef_iterator_get_fe_around
  end type fe_vef_iterator_t
    
  type, extends(fe_vef_iterator_t) :: fe_face_iterator_t
    private
    class(fe_iterator_t)  , allocatable :: fe
    class(face_iterator_t), pointer     :: face
   contains
    procedure         , private         :: create                                     => fe_face_iterator_create
    procedure         , private         :: free                                       => fe_face_iterator_free
    procedure         , non_overridable :: update_integration                         => fe_face_iterator_update_integration 
    procedure         , non_overridable :: get_fe_space                               => fe_face_iterator_get_fe_space
    procedure         , non_overridable :: get_elem2dof                               => fe_face_iterator_get_elem2dof
    procedure         , non_overridable :: get_quadrature                             => fe_face_iterator_get_quadrature
    procedure         , non_overridable :: get_face_map                               => fe_face_iterator_get_face_map
    procedure         , non_overridable :: get_face_integrator                        => fe_face_iterator_get_face_integrator
    procedure         , non_overridable :: impose_strong_dirichlet_bcs                => fe_face_iterator_impose_strong_dirichlet_bcs
    procedure         , non_overridable :: compute_surface                            => fe_face_iterator_compute_surface
     ! Methods exploiting ("inherited from") face_t
    procedure, non_overridable          :: get_coordinates                  => fe_face_iterator_get_coordinates
    procedure, non_overridable          :: get_face_lid                     => fe_face_iterator_get_face_lid
    procedure, non_overridable          :: get_face_lpos_within_cell_around => fe_face_iterator_get_face_lpos_within_cell_around
    procedure, non_overridable          :: get_face_orientation             => fe_face_iterator_get_face_orientation
    procedure, non_overridable          :: get_face_rotation                => fe_face_iterator_get_face_rotation
  end type fe_face_iterator_t
      
  integer(ip), parameter :: fe_space_type_cg            = 0 ! H^1 conforming FE space
  integer(ip), parameter :: fe_space_type_dg            = 1 ! L^2 conforming FE space + .not. H^1 conforming (weakly imposed via face integration)
  integer(ip), parameter :: fe_space_type_dg_conforming = 2 ! DG approximation of L^2 spaces (does not involve coupling by face)
  
  type, extends(base_fe_space_t) :: serial_fe_space_t 
     private     
     ! Reference FE container
     integer(ip)                                 :: reference_fes_size
     type(p_reference_fe_t)        , allocatable :: reference_fes(:)
     
     ! Finite Element-related integration containers
     type(quadrature_t)            , allocatable :: fe_quadratures(:)
     type(fe_map_t)                , allocatable :: fe_maps(:)
     type(volume_integrator_t)     , allocatable :: fe_volume_integrators(:)
     
     ! Mapping of FEs to reference FE and FEs-related integration containers
     integer(ip)                   , allocatable :: reference_fe_id_per_fe(:,:)         ! (number_fields, number_fes)
     integer(ip)                   , allocatable :: max_order_reference_fe_id_per_fe(:) ! Stores Key=max_order_reference_fe_id for all FEs
     type(hash_table_ip_ip_t)                    :: fe_quadratures_and_maps_position    ! Key = max_order_reference_fe_id
     type(hash_table_ip_ip_t)                    :: fe_volume_integrators_position      ! Key = [max_order_reference_fe_id,
                                                                                        !       reference_fe_id]
     
     ! Finite Face-related integration containers
     type(quadrature_t)            , allocatable :: fe_face_quadratures(:)
     type(face_map_t)              , allocatable :: fe_face_maps(:)
     type(face_integrator_t)       , allocatable :: fe_face_integrators(:)
     
     
     ! Mapping of Finite Faces and integration containers
     integer(ip)                   , allocatable :: max_order_per_fe_face(:)      ! Stores Key=max_order_fes_around_fe_face for all faces
     type(hash_table_ip_ip_t)                    :: fe_face_quadratures_position  ! Key=max_order_fes_around_fe_face
     type(hash_table_ip_ip_t)                    :: fe_face_maps_position         ! Key=[max_order_fes_around_fe_face, 
                                                                                  !     left_geo_reference_fe_id,
                                                                                  !     right_geo_reference_fe_id (with 0 for boundary faces)]
     type(hash_table_ip_ip_t)                    :: fe_face_integrators_position  ! Key = [max_order_fes_around_fe_face,
                                                                                  !       left_reference_fe_id,
                                                                                  !       left_reference_fe_id (with 0 for boundary faces)]
    
     ! DoF identifiers associated to each FE and field within FE
     integer(ip)                   , allocatable :: ptr_dofs_per_fe(:,:) ! (number_fields, number_fes+1)
     integer(ip)                   , allocatable :: lst_dofs_lids(:)
    
     ! Strong Dirichlet BCs-related member variables
     integer(ip)                                 :: number_strong_dirichlet_dofs
     type(serial_scalar_array_t)                 :: strong_dirichlet_values
     logical                       , allocatable :: at_strong_dirichlet_boundary_per_fe(:,:)
     
     class(base_static_triangulation_t), pointer :: triangulation =>  NULL()
   contains
     procedure,                  private :: serial_fe_space_create_same_reference_fes_on_all_cells
     procedure,                  private :: serial_fe_space_create_different_between_cells
     generic                             :: create                                       => serial_fe_space_create_same_reference_fes_on_all_cells,&
                                                                                            serial_fe_space_create_different_between_cells
     procedure                           :: free                                         => serial_fe_space_free
     procedure                           :: print                                        => serial_fe_space_print
     procedure, non_overridable, private :: allocate_and_fill_reference_fes              => serial_fe_space_allocate_and_fill_reference_fes
     procedure, non_overridable, private :: free_reference_fes                           => serial_fe_space_free_reference_fes
     procedure, non_overridable, private :: allocate_and_fill_field_blocks_and_coupling  => serial_fe_space_allocate_and_fill_field_blocks_and_coupling 
     procedure, non_overridable, private :: free_field_blocks_and_coupling               => serial_fe_space_free_field_blocks_and_coupling
     procedure, non_overridable, private :: allocate_ref_fe_id_per_fe                    => serial_fe_space_allocate_ref_fe_id_per_fe
     procedure, non_overridable, private :: free_ref_fe_id_per_fe                        => serial_fe_space_free_ref_fe_id_per_fe
     procedure, non_overridable, private :: fill_ref_fe_id_per_fe_same_on_all_cells      => serial_fe_space_fill_ref_fe_id_per_fe_same_on_all_cells
     procedure, non_overridable, private :: fill_ref_fe_id_per_fe_different_between_cells=> serial_fe_space_fill_ref_fe_id_per_fe_different_between_cells
     procedure, non_overridable, private :: check_cell_vs_fe_topology_consistency        => serial_fe_space_check_cell_vs_fe_topology_consistency
     procedure, non_overridable, private :: allocate_and_fill_fe_space_type_per_field    => serial_fe_space_allocate_and_fill_fe_space_type_per_field
     procedure, non_overridable, private :: free_fe_space_type_per_field                 => serial_fe_space_free_fe_space_type_per_field
     procedure, non_overridable, private :: allocate_and_init_ptr_lst_dofs               => serial_fe_space_allocate_and_init_ptr_lst_dofs
     procedure, non_overridable, private :: free_ptr_lst_dofs                            => serial_fe_space_free_ptr_lst_dofs
     procedure, non_overridable, private :: allocate_and_init_at_strong_dirichlet_bound  => serial_fe_space_allocate_and_init_at_strong_dirichlet_bound  
     procedure, non_overridable, private :: free_at_strong_dirichlet_bound               => serial_fe_space_free_at_strong_dirichlet_bound
     procedure, non_overridable, private :: set_up_strong_dirichlet_bcs                  => serial_fe_space_set_up_strong_dirichlet_bcs
     procedure, non_overridable, private :: set_up_strong_dirichlet_bcs_on_vef_and_field => serial_fe_space_set_up_strong_dirichlet_bcs_on_vef_and_field
     procedure                           :: interpolate_dirichlet_values                 => serial_fe_space_interpolate_dirichlet_values
     
     procedure                           :: project_dirichlet_values_curl_conforming     => serial_fe_space_project_dirichlet_values_curl_conforming
     procedure, non_overridable, private :: allocate_and_fill_fields_to_project_         => serial_fe_space_allocate_and_fill_fields_to_project_
     procedure, non_overridable, private :: allocate_and_fill_offset_component           => serial_fe_space_allocate_and_fill_offset_component
     procedure, non_overridable, private :: allocate_and_fill_global2subset_and_inverse  => serial_fe_space_allocate_and_fill_global2subset_and_inverse 
     procedure, non_overridable, private :: get_function_scalar_components               => serial_fe_space_get_function_scalar_components
     procedure, non_overridable, private :: evaluate_vector_function_scalar_components   => serial_fe_space_evaluate_vector_function_scalar_components
     procedure, non_overridable, private :: project_curl_conforming_compute_elmat_elvec  => serial_fe_space_project_curl_conforming_compute_elmat_elvec
     
     procedure, non_overridable          :: initialize_fe_integration                    => serial_fe_space_initialize_fe_integration
     procedure, non_overridable, private :: free_fe_integration                          => serial_fe_space_free_fe_integration
     procedure, non_overridable, private :: generate_fe_volume_integrators_position_key  => serial_fe_space_generate_fe_volume_integrators_position_key
     
     procedure, non_overridable          :: initialize_fe_face_integration               => serial_fe_space_initialize_fe_face_integration
     procedure, non_overridable, private :: free_fe_face_integration                     => serial_fe_space_free_fe_face_integration
     procedure, non_overridable, private :: generate_fe_face_maps_position_key           => serial_fe_space_fe_face_maps_position_key
     procedure, non_overridable, private :: generate_fe_face_integrators_position_key    => serial_fe_space_fe_face_integrators_position_key

     procedure                           :: create_assembler                             => serial_fe_space_create_assembler
     procedure                           :: symbolic_setup_assembler                     => serial_fe_space_symbolic_setup_assembler
     procedure                           :: create_dof_values                            => serial_fe_space_create_dof_values
     procedure                           :: fill_dof_info                                => serial_fe_space_fill_dof_info
     procedure                 , private :: fill_elem2dof_and_count_dofs                 => serial_fe_space_fill_elem2dof_and_count_dofs
     procedure                 , private :: renumber_dofs_block                          => serial_fe_space_renumber_dofs_block
 
     ! Getters
     procedure                           :: get_num_dimensions                           => serial_fe_space_get_num_dimensions
     procedure, non_overridable          :: get_number_reference_fes                     => serial_fe_space_get_number_reference_fes
     procedure, non_overridable          :: get_reference_fe                             => serial_fe_space_get_reference_fe
     procedure, non_overridable          :: get_field_type                               => serial_fe_space_get_field_type 
     procedure, non_overridable          :: get_number_components                        => serial_fe_space_get_number_components
     procedure, non_overridable          :: get_max_number_shape_functions               => serial_fe_space_get_max_number_shape_functions
     procedure, non_overridable          :: get_max_number_dofs_on_a_cell                => serial_fe_space_get_max_number_dofs_on_a_cell
     procedure, non_overridable          :: get_max_number_quadrature_points             => serial_fe_space_get_max_number_quadrature_points
     procedure, non_overridable          :: get_max_number_nodal_quadrature_points       => serial_fe_space_get_max_number_nodal_quadrature_points
     procedure, non_overridable          :: get_max_number_face_quadrature_points        => serial_fe_space_get_max_number_face_quadrature_points     
     procedure, non_overridable          :: get_triangulation                            => serial_fe_space_get_triangulation
     procedure                           :: get_environment                              => serial_fe_space_get_environment
     procedure                           :: get_strong_dirichlet_values                  => serial_fe_space_get_strong_dirichlet_values
     
     ! fes, fe_vefs and fe_faces traversals-related TBPs
     procedure                           :: create_fe_iterator                           => serial_fe_space_create_fe_iterator
     procedure                           :: free_fe_iterator                             => serial_fe_space_free_fe_iterator
     procedure, non_overridable          :: create_fe_vef_iterator                       => serial_fe_space_create_fe_vef_iterator     
     procedure, non_overridable          :: create_itfc_fe_vef_iterator                  => serial_fe_space_create_itfc_fe_vef_iterator     
     procedure, non_overridable          :: create_fe_face_iterator                      => serial_fe_space_create_fe_face_iterator
     procedure, non_overridable          :: free_fe_vef_iterator                         => serial_fe_space_free_fe_vef_iterator
 end type serial_fe_space_t  
 
 public :: serial_fe_space_t   
 public :: fe_iterator_t
 public :: fe_vef_iterator_t
 public :: fe_face_iterator_t
 
 type base_fe_object_iterator_t
   private
   type(object_iterator_t) :: object
 contains
   ! Methods "inherited" from object_iterator_t
   procedure                            :: first                                 => base_fe_object_iterator_first
   procedure                            :: next                                  => base_fe_object_iterator_next
   procedure                            :: set_lid                               => base_fe_object_iterator_set_lid
   procedure                            :: has_finished                          => base_fe_object_iterator_has_finished
   procedure, non_overridable           :: get_lid                               => base_fe_object_iterator_get_lid
   procedure, non_overridable           :: get_gid                               => base_fe_object_iterator_get_gid
   procedure, non_overridable           :: get_dimension                         => base_fe_object_iterator_get_dimension
   procedure, non_overridable           :: get_number_parts_around               => base_fe_object_iterator_get_number_parts_around
   procedure, non_overridable           :: get_number_subparts_around            => base_fe_object_iterator_get_number_subparts_around
   procedure, non_overridable           :: create_parts_around_iterator          => base_fe_object_iterator_create_parts_around_iterator
   procedure, non_overridable           :: create_subparts_around_iterator       => base_fe_object_iterator_create_subparts_around_iterator
   procedure, non_overridable           :: get_num_vefs                          => base_fe_object_iterator_get_num_vefs
   procedure, non_overridable           :: get_num_faces                         => base_fe_object_iterator_get_num_faces
   procedure, non_overridable, private  :: base_fe_object_iterator_get_vef
   procedure, non_overridable, private  :: base_fe_object_iterator_get_face
   generic                              :: get_vef                               => base_fe_object_iterator_get_vef
   generic                              :: get_face                              => base_fe_object_iterator_get_face
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
    procedure, non_overridable, private  :: fe_object_iterator_get_fe_face
    generic                              :: get_vef                               => fe_object_iterator_get_fe_vef 
    generic                              :: get_face                              => fe_object_iterator_get_fe_face
    procedure, non_overridable           :: get_number_coarse_dofs                => fe_object_iterator_get_number_coarse_dofs
    procedure, non_overridable           :: create_own_coarse_dofs_iterator       => fe_object_iterator_create_own_coarse_dofs_iterator
  end type fe_object_iterator_t
        
  ! These parameter constants are used in order to generate a unique (non-consecutive) 
  ! but consistent across MPI tasks global ID (integer(igp)) of a given DoF.
  ! See type(par_fe_space_t)%generate_non_consecutive_dof_gid()
  integer(ip), parameter :: cell_gid_shift              = 44
  integer(ip), parameter :: dofs_per_reference_fe_shift = 14
  integer(ip), parameter :: number_fields_shift         = igp*8-(cell_gid_shift+dofs_per_reference_fe_shift)-1
  
  ! These three parameter constants are thought be used as FPL keys. The corresponding pairs 
  ! <XXXkey,.true.|.false.> in FPL let the user to control whether or not coarse vertex, edge, or 
  ! face DoFs are to be included into coarse_fe_space_t. These parameters might be used during 
  ! coarse_fe_space_t set-up, as well as by the deferred TBP methods corresponding to 
  ! class(coarse_fe_handler_t).
  character(len=*), parameter :: coarse_space_use_vertices_key = 'coarse_space_use_vertices_key'
  character(len=*), parameter :: coarse_space_use_edges_key    = 'coarse_space_use_edges_key'
  character(len=*), parameter :: coarse_space_use_faces_key    = 'coarse_space_use_faces_key'

  
 type, extends(serial_fe_space_t) :: par_fe_space_t
   private   
   ! Multilevel fe space   
   ! It is the equivalent to the "element_to_dof" at the finer level
   ! Pointers to the start/end of coarse DoFs LIDs of each field (lst_coarse_dofs)
   integer(ip)                   , allocatable :: ptr_coarse_dofs_per_field(:)	
   ! List of coarse DoFs LIDs
   integer(ip)                   , allocatable :: lst_coarse_dofs(:)	
   
   ! Coarse DoFs LIDs on top of coarse n_faces per field
   ! It provides (for every field) what we get from the reference_fe_t at the
   ! finer level
   type(list_t), allocatable                   :: own_coarse_dofs_per_field(:)
   
   ! GIDs of the coarse DoFs. For their generation, we though to be a good idea 
   ! to take as basis the GIDs of the VEFs objects they are built from (instead of
   ! generating them from scratch, which in turn would imply further communication).
   type(allocatable_array_igp1_t), allocatable :: coarse_dof_gids_per_field(:)
   
   ! Polymorphic data type in charge of filling some of the member variables above
   ! (so far, lst_coarse_dofs + own_coarse_dofs_per_field)
   class(l1_coarse_fe_handler_t), pointer      :: coarse_fe_handler => NULL()
 contains
   procedure, private :: serial_fe_space_create_same_reference_fes_on_all_cells                   => par_fe_space_serial_create_same_reference_fes_on_all_cells 
   procedure, private :: serial_fe_space_create_different_between_cells                           => par_fe_space_serial_create_different_between_cells 
   procedure, private :: par_fe_space_create_same_reference_fes_on_all_cells 
   procedure, private :: par_fe_space_create_different_between_cells
   generic                                     :: create                                          => par_fe_space_create_same_reference_fes_on_all_cells, &
                                                                                                     par_fe_space_create_different_between_cells
   procedure                                   :: fill_dof_info                                   => par_fe_space_fill_dof_info
   procedure                         , private :: fill_elem2dof_and_count_dofs                    => par_fe_space_fill_elem2dof_and_count_dofs
   procedure                                   :: renumber_dofs_first_interior_then_interface     => par_fe_space_renumber_dofs_first_interior_then_interface
   procedure        , non_overridable, private :: renumber_interface_dofs_first_E_then_Ec         => par_fe_space_renumber_interface_dofs_first_E_then_Ec
   procedure        , non_overridable, private :: set_up_strong_dirichlet_bcs_ghost_fes           => par_fe_space_set_up_strong_dirichlet_bcs_ghost_fes

   procedure        , non_overridable, private :: compute_blocks_dof_import                       => par_fe_space_compute_blocks_dof_import
   procedure        , non_overridable, private :: compute_dof_import                              => par_fe_space_compute_dof_import
   procedure        , non_overridable, private :: compute_raw_interface_data_by_continuity        => par_fe_space_compute_raw_interface_data_by_continuity
   procedure        , non_overridable, private :: raw_interface_data_by_continuity_decide_owner   => par_fe_space_raw_interface_data_by_continuity_decide_owner
   procedure        , non_overridable, private :: compute_raw_interface_data_by_face_integ        => par_fe_space_compute_raw_interface_data_by_face_integ
   procedure        , non_overridable, private :: compute_ubound_num_itfc_couplings_by_continuity => pfs_compute_ubound_num_itfc_couplings_by_continuity
   procedure        , non_overridable, private :: compute_ubound_num_itfc_couplings_by_face_integ => pfs_compute_ubound_num_itfc_couplings_by_face_integ
   procedure, nopass, non_overridable, private :: generate_non_consecutive_dof_gid                => par_fe_space_generate_non_consecutive_dof_gid
   
   ! These set of three subroutines are in charge of generating a dof_import for the distributed-memory solution of boundary mass matrices
   procedure        , non_overridable, private :: compute_boundary_dof_import                              => par_fe_space_compute_boundary_dof_import
   procedure        , non_overridable, private :: compute_ubound_boundary_num_itfc_couplings_by_continuity => pfs_compute_ubound_boundary_num_itfc_couplings_by_continuity
   procedure        , non_overridable, private :: compute_boundary_raw_interface_data_by_continuity        => par_fe_space_compute_boundary_raw_interface_data_by_continuity


   procedure        , non_overridable          :: get_number_fe_objects                           => par_fe_space_get_number_fe_objects
   procedure                                   :: get_par_environment                             => par_fe_space_get_par_environment
   procedure                                   :: get_environment                                 => par_fe_space_get_environment
   
   procedure                                   :: print                                           => par_fe_space_print
   procedure                                   :: free                                            => par_fe_space_free
   procedure                                   :: create_assembler                                => par_fe_space_create_assembler
   procedure                                   :: symbolic_setup_assembler                        => par_fe_space_symbolic_setup_assembler
   procedure                                   :: create_dof_values                               => par_fe_space_create_dof_values
   procedure                                   :: interpolate_dirichlet_values                    => par_fe_space_interpolate_dirichlet_values
   procedure                                   :: project_dirichlet_values_curl_conforming        => par_fe_space_project_dirichlet_values_curl_conforming
   
   procedure        , non_overridable, private  :: setup_coarse_dofs                               => par_fe_space_setup_coarse_dofs
   procedure, nopass, non_overridable, private  :: generate_coarse_dof_gid                         => par_fe_space_generate_coarse_dof_gid
   procedure        , non_overridable, private  :: free_coarse_dofs                                => par_fe_space_free_coarse_dofs
   procedure        , non_overridable           :: setup_coarse_fe_space                           => par_fe_space_setup_coarse_fe_space
   procedure        , non_overridable, private  :: transfer_number_fields                          => par_fe_space_transfer_number_fields
   procedure        , non_overridable, private  :: transfer_fe_space_type                          => par_fe_space_transfer_fe_space_type
   procedure        , non_overridable, private  :: gather_ptr_dofs_per_fe_and_field                => par_fe_space_gather_ptr_dofs_per_fe_and_field
   procedure        , non_overridable, private  :: gather_coarse_dofs_gids_rcv_counts_and_displs   => par_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
   procedure        , non_overridable, private  :: gather_coarse_dofs_gids                         => par_fe_space_gather_coarse_dofs_gids
   procedure        , non_overridable, private  :: gather_vefs_gids_dofs_objects                   => par_fe_space_gather_vefs_gids_dofs_objects
   procedure                                    :: get_total_number_coarse_dofs                    => par_fe_space_get_total_number_coarse_dofs
   procedure                                    :: get_block_number_coarse_dofs                    => par_fe_space_get_block_number_coarse_dofs
   procedure       , non_overridable            :: get_coarse_fe_handler                           => par_fe_space_get_coarse_fe_handler

   ! Objects-related traversals
   procedure, non_overridable                   :: create_fe_object_iterator                       => par_fe_space_create_fe_object_iterator
   procedure, non_overridable                   :: free_fe_object_iterator                         => par_fe_space_free_fe_object_iterator
   end type par_fe_space_t
 
 public :: par_fe_space_t
 public :: fe_object_iterator_t
 
  type, abstract :: l1_coarse_fe_handler_t
  contains
    ! Deferred methods
       procedure (l1_compute_change_basis_matrix)           , deferred :: compute_change_basis_matrix 
       procedure (l1_get_num_coarse_dofs_interface)         , deferred :: get_num_coarse_dofs
	   procedure (l1_renumber_interface_dofs_first_E_then_Ec),deferred :: renumber_interface_dofs_first_E_then_Ec
	   procedure (l1_setup_constraint_matrix)               , deferred :: setup_constraint_matrix
	   procedure (l1_setup_weighting_operator)              , deferred :: setup_weighting_operator
	   procedure (l1_apply_weighting_operator_and_comm)     , deferred :: apply_weighting_operator_and_comm
	   procedure (l1_apply_transpose_weighting_operator)    , deferred :: apply_transpose_weighting_operator 
  end type l1_coarse_fe_handler_t
 
  abstract interface	
    ! In H(curl) conforming spaces a specific change of basis is needed in the 3D case. 
    ! In all other cases this subroutine will have no impact 
    subroutine l1_compute_change_basis_matrix( this, par_fe_space ) 
	import :: l1_coarse_fe_handler_t, par_fe_space_t
	class(l1_coarse_fe_handler_t), intent(inout) :: this
    type(par_fe_space_t)         , intent(inout) :: par_fe_space 
	end subroutine l1_compute_change_basis_matrix 
	
    ! Returns the number of coarse DoFs that the object customizing
    ! l1_coarse_fe_handler_t requires to introduce on the subdomain 
    ! interface
    subroutine l1_get_num_coarse_dofs_interface(this, par_fe_space, parameter_list, num_coarse_dofs) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, ip
      class(l1_coarse_fe_handler_t), intent(in)    :: this
      type(par_fe_space_t)         , intent(in)    :: par_fe_space 
      type(parameterlist_t)        , intent(in)    :: parameter_list
      integer(ip)                  , intent(inout) :: num_coarse_dofs(:)
    end subroutine l1_get_num_coarse_dofs_interface
	
	subroutine l1_renumber_interface_dofs_first_E_then_Ec(this, par_fe_space, iblock, perm_old2new ) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, ip
      class(l1_coarse_fe_handler_t), intent(inout) :: this
      type(par_fe_space_t)         , intent(in)    :: par_fe_space 
	  integer(ip)                  , intent(in)    :: iblock 
	  integer(ip)  , allocatable   , intent(inout) :: perm_old2new(:)
    end subroutine l1_renumber_interface_dofs_first_E_then_Ec
   
    subroutine l1_setup_constraint_matrix(this, par_fe_space, parameter_list, constraint_matrix) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, coo_sparse_matrix_t
      class(l1_coarse_fe_handler_t), intent(in)    :: this
      type(par_fe_space_t)         , intent(in)    :: par_fe_space
      type(parameterlist_t)        , intent(in)    :: parameter_list
	     type(coo_sparse_matrix_t)    , intent(inout) :: constraint_matrix
    end subroutine l1_setup_constraint_matrix
  
    subroutine l1_setup_weighting_operator(this, par_fe_space, parameter_list, weighting_operator) 
	     import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, operator_t, rp
      class(l1_coarse_fe_handler_t) , intent(in)    :: this
      type(par_fe_space_t)          , intent(in)    :: par_fe_space
      type(parameterlist_t)         , intent(in)    :: parameter_list
	     real(rp)         , allocatable, intent(inout) :: weighting_operator(:)
    end subroutine l1_setup_weighting_operator
	
	 subroutine l1_apply_weighting_operator_and_comm(this, W, x, y) 
	     import :: l1_coarse_fe_handler_t, par_scalar_array_t, rp
      class(l1_coarse_fe_handler_t) , intent(in)    :: this
	  real(rp)         , allocatable, intent(in)    :: W(:)
	  type(par_scalar_array_t)      , intent(inout) :: x
      type(par_scalar_array_t)      , intent(inout) :: y
    end subroutine l1_apply_weighting_operator_and_comm
	
	 subroutine l1_apply_transpose_weighting_operator(this, W, x, y) 
	     import :: l1_coarse_fe_handler_t, par_scalar_array_t, rp
      class(l1_coarse_fe_handler_t) , intent(in)    :: this
	  real(rp)         , allocatable, intent(in)    :: W(:)
	  type(par_scalar_array_t)      , intent(inout) :: x
      type(par_scalar_array_t)      , intent(inout) :: y
    end subroutine l1_apply_transpose_weighting_operator

  end interface
  
  type, extends(l1_coarse_fe_handler_t) :: standard_l1_coarse_fe_handler_t
    private
  contains
       procedure             :: compute_change_basis_matrix               => standard_l1_compute_change_basis_matrix    
       procedure             :: get_num_coarse_dofs                       => standard_l1_get_num_coarse_dofs
	   procedure             :: renumber_interface_dofs_first_E_then_Ec   => standard_l1_renumber_interface_dofs_first_E_then_Ec 
	   procedure             :: setup_constraint_matrix                   => standard_l1_setup_constraint_matrix
	   procedure             :: setup_weighting_operator                  => standard_l1_setup_weighting_operator
	   procedure             :: apply_weighting_operator_and_comm         => standard_l1_apply_weighting_operator_and_comm 
	   procedure             :: apply_transpose_weighting_operator        => standard_l1_apply_transpose_weighting_operator 
       procedure, nopass     :: get_coarse_space_use_vertices_edges_faces => standard_get_coarse_space_use_vertices_edges_faces
	   procedure             :: free                                      => standard_l1_free 
  end type standard_l1_coarse_fe_handler_t
  
  type, extends(standard_l1_coarse_fe_handler_t) :: H1_l1_coarse_fe_handler_t
    private
    real(rp), public :: diffusion_inclusion
  contains
	   procedure :: setup_constraint_matrix  => H1_l1_setup_constraint_matrix
	   procedure :: setup_weighting_operator => H1_l1_setup_weighting_operator
  end type H1_l1_coarse_fe_handler_t
  
     type :: edge_change_basis_matrix_t 
   private 
   	   ! Change Basis Matrix is stored as 
	   ! [ G 0      G = diag(G_Ei) block diagonal for every coarse edge  
	   !   B I ];   B = coupling interface dofs with coarse edge dofs 
	   type(sparse_matrix_t), allocatable  :: G(:) 
	   type(sparse_matrix_t)               :: B  
   end type edge_change_basis_matrix_t
   
    type, extends(standard_l1_coarse_fe_handler_t) :: Hcurl_l1_coarse_fe_handler_t
    private 
	   integer(ip)                             :: order
	   integer(ip)                             :: number_interior_dofs 
	   integer(ip)                             :: number_edge_wire_dofs
	   integer(ip)                             :: number_total_wire_dofs 	   
	   integer(ip)                             :: number_coarse_edges
	   integer(ip) , allocatable               :: number_fine_edges_per_coarse_edge(:) 
	   integer(ip) , allocatable               :: perm_sorted_edges(:,:) 
	   logical     , allocatable               :: follows_coarse_edge_orientation(:,:)
	   type(hash_table_ip_ip_t), allocatable   :: coupled_vefs_added(:)
	   type(hash_table_ip_ip_t), allocatable   :: coarse_edge_nodes_order(:)
       type(edge_change_basis_matrix_t)        :: change_basis_matrix 
	   
  contains
       ! Overriding procedures 
       procedure                           :: free                                          => Hcurl_l1_free 
	   procedure                           :: get_num_coarse_dofs                           => Hcurl_l1_get_num_coarse_dofs 
	   procedure                           :: setup_constraint_matrix                       => Hcurl_l1_setup_constraint_matrix
	   procedure                           :: apply_weighting_operator_and_comm             => Hcurl_l1_apply_weighting_operator_and_comm   
	   procedure                           :: apply_transpose_weighting_operator            => Hcurl_l1_apply_transpose_weighting_operator
	   procedure                           :: renumber_interface_dofs_first_E_then_Ec       => Hcurl_l1_renumber_interface_dofs_first_E_then_Ec	   
	   procedure                           :: compute_change_basis_matrix                   => Hcurl_l1_compute_change_basis_matrix 
	   ! Auxiliar procedures for renumbering 
	   procedure, non_overridable, private :: fill_dofs_in_coarse_edges_renumbering         => Hcurl_l1_fill_dofs_in_coarse_edges_renumbering 
	   procedure, non_overridable, private :: fill_dofs_coupled_to_coarse_edges_renumbering => Hcurl_l1_fill_dofs_coupled_to_coarse_edges_renumbering
	   procedure, non_overridable, private, nopass :: set_coarse_edge_orientation           => Hcurl_l1_set_coarse_edge_orientation
	   procedure, non_overridable, private :: sort_fine_edges_within_coarse_edge            => Hcurl_l1_sort_fine_edges_within_coarse_edge
	   ! Change of basis construction     
	   procedure, non_overridable, private :: compute_edge_discrete_gradient_elmat          => Hcurl_l1_compute_edge_discrete_gradient_elmat
	   procedure, non_overridable, private :: compute_face_discrete_gradient_elmat          => Hcurl_l1_compute_face_discrete_gradient_elmat
	   ! Change of basis application 
	   procedure, non_overridable, private :: apply_global_change_basis                     => Hcurl_l1_apply_global_change_basis
	   procedure, non_overridable, private :: apply_global_change_basis_transpose           => Hcurl_l1_apply_global_change_basis_transpose
	   procedure, non_overridable, private :: apply_inverse_local_change_basis              => Hcurl_l1_apply_inverse_local_change_basis 
	   procedure, non_overridable, private :: apply_inverse_local_change_basis_transpose    => Hcurl_l1_apply_inverse_local_change_basis_transpose
	   ! Private TBPs 	 
	   procedure, non_overridable, private :: compute_first_order_moment_in_edges           => Hcurl_l1_compute_first_order_moment_in_edges
	   procedure, non_overridable, private :: compute_edge_elvec                            => Hcurl_l1_compute_edge_elvec
	   procedure, non_overridable, private :: fill_edge_local_change_of_basis               => Hcurl_l1_fill_edge_local_change_of_basis  
	   procedure, non_overridable, private :: fill_edge_coupled_to_edges_local_change_of_basis => Hcurl_l1_fill_edge_coupled_to_edges_local_change_of_basis
	   procedure, non_overridable, private :: fill_face_coupled_to_edges_local_change_of_basis => Hcurl_l1_fill_face_coupled_to_edges_local_change_of_basis
	   procedure, non_overridable, private :: find_interface_fe_face_around_edge            => Hcurl_l1_find_interface_fe_face_around_edge 
	   procedure, non_overridable, private :: assemble_face_coupled_to_edges_B_elmat        => Hcurl_l1_assemble_face_coupled_to_edges_B_elmat
	   ! DoF numbering getters 
	   procedure, non_overridable, private :: get_wire_basis_dofs_from_vef     => Hcurl_l1_get_wire_basis_dofs_from_vef 
	   procedure, non_overridable, private :: get_new_basis_dof_from_node_id   => Hcurl_l1_get_new_basis_dof_from_node_id 
	   procedure, non_overridable, private :: get_new_basis_dofs_from_edge     => Hcurl_l1_get_new_basis_dofs_from_edge
	   ! Logical getters 
	   procedure, non_overridable, private :: is_first_edge                    => Hcurl_l1_is_first_edge
	   procedure, non_overridable, private :: is_last_edge                     => Hcurl_l1_is_last_edge 
	   ! Auxiliar basic procedures 
	   procedure, non_overridable, private :: compute_edge_length              => Hcurl_l1_compute_edge_length 
  end type  Hcurl_l1_coarse_fe_handler_t
     
  public :: l1_coarse_fe_handler_t, standard_l1_coarse_fe_handler_t, H1_l1_coarse_fe_handler_t, Hcurl_l1_coarse_fe_handler_t 
    
  type , extends(base_fe_iterator_t) :: coarse_fe_iterator_t
    private
    type(coarse_fe_space_t), pointer :: coarse_fe_space => NULL()
  contains
    procedure, non_overridable, private :: create                                     => coarse_fe_iterator_create
    procedure, non_overridable, private :: free                                       => coarse_fe_iterator_free
    final                               :: coarse_fe_iterator_free_final
    procedure, non_overridable          :: create_own_dofs_on_vef_iterator            => coarse_fe_iterator_create_own_dofs_on_vef_iterator
    procedure, non_overridable, private :: fill_own_dofs_on_vef                       => coarse_fe_iterator_fill_own_dofs_on_vef
    procedure, non_overridable, private :: fill_own_dofs_on_vef_from_source_coarse_fe => coarse_fe_iterator_fill_own_dofs_on_vef_from_source_coarse_fe
    procedure, non_overridable, private :: renumber_dofs_block                        => coarse_fe_iterator_renumber_dofs_block
    procedure, non_overridable, private :: renumber_dofs_field                        => coarse_fe_iterator_renumber_dofs_field
    
    procedure, non_overridable, private :: get_scan_sum_number_dofs                   => coarse_fe_iterator_get_scan_sum_number_dofs
    procedure, non_overridable          :: get_number_fe_spaces                       => coarse_fe_iterator_get_number_fe_spaces
    procedure, non_overridable          :: get_number_dofs                            => coarse_fe_iterator_get_number_dofs
    procedure, non_overridable          :: get_elem2dof                               => coarse_fe_iterator_get_elem2dof
    procedure, non_overridable          :: get_field_elem2dof                         => coarse_fe_iterator_get_field_elem2dof
    procedure, non_overridable          :: get_coarse_fe_vef                          => coarse_fe_iterator_get_coarse_fe_vef
  end type coarse_fe_iterator_t
  
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
    procedure, non_overridable           :: get_number_coarse_dofs                   => coarse_fe_object_iterator_get_number_coarse_dofs
    procedure, non_overridable           :: create_own_coarse_dofs_iterator          => coarse_fe_object_iterator_create_own_coarse_dofs_iterator
  end type coarse_fe_object_iterator_t
    
  
  type, extends(base_fe_space_t) :: coarse_fe_space_t
    private
    integer(ip) , allocatable                   :: ptr_dofs_per_fe_and_field(:)
    integer(ip) , allocatable                   :: lst_dofs_lids(:)
    type(list_t), allocatable                   :: own_dofs_vef_per_fe(:)
    	
	   ! Pointer to coarse triangulation this coarse_fe_space has been built from
    type(coarse_triangulation_t), pointer       :: coarse_triangulation => NULL()
	
    ! It is the equivalent to the "element_to_dof" at the finer level
    ! Pointers to the start/end of coarse DoFs LIDs of each field (lst_coarse_dofs)
    integer(ip)                   , allocatable :: ptr_coarse_dofs_per_field(:)	
    ! List of coarse DoFs LIDs
    integer(ip)                   , allocatable :: lst_coarse_dofs(:)	
    
    ! Coarse DoFs LIDs on top of coarse n_faces per field
    ! It provides (for every field) what we get from the reference_fe_t at the
    ! finer level
    type(list_t), allocatable                   :: own_coarse_dofs_per_field(:)
   
    
    ! GIDs of the coarse DoFs. For their generation, we though to be a good idea 
    ! to take as basis the GIDs of the VEFs objects they are built from (instead of
    ! generating them from scratch, which in turn would imply further communication).
    type(allocatable_array_igp1_t), allocatable :: coarse_dof_gids_per_field(:)
  contains
    procedure                                   :: create                                          => coarse_fe_space_create
    procedure                                   :: free                                            => coarse_fe_space_free
    procedure, non_overridable                  :: print                                           => coarse_fe_space_print
    procedure, non_overridable, private         :: allocate_and_fill_field_blocks_and_coupling     => coarse_fe_space_allocate_and_fill_field_blocks_and_coupling
    procedure, non_overridable, private         :: free_field_blocks_and_coupling                  => coarse_fe_space_free_field_blocks_and_coupling
    procedure, non_overridable, private         :: allocate_and_fill_fe_space_type_per_field       => coarse_fe_space_allocate_and_fill_fe_space_type
    procedure, non_overridable, private         :: free_fe_space_type_per_field                    => coarse_fe_space_free_fe_space_type
    procedure, non_overridable, private         :: allocate_and_fill_ptr_dofs_per_fe_and_field     => coarse_fe_space_allocate_and_fill_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: free_ptr_dofs_per_fe_and_field                  => coarse_fe_space_free_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: fetch_ghost_fes_data                            => coarse_fe_space_fetch_ghost_fes_data
    procedure, non_overridable, private         :: allocate_and_fill_own_dofs_vef_per_fe           => coarse_fe_space_allocate_and_fill_own_dofs_vef_per_fe
    procedure, non_overridable, private         :: fill_own_dofs_per_fe_field                      => coarse_fe_space_fill_own_dofs_per_fe_field
    procedure, non_overridable, private         :: free_own_dofs_vef_per_fe                        => coarse_fe_space_free_free_own_dofs_vef_per_fe
    procedure, non_overridable, private         :: allocate_lst_dofs_lids                          => coarse_fe_space_allocate_lst_dofs_lids
    procedure, non_overridable, private         :: count_dofs_and_fill_lst_dof_lids                => coarse_fe_space_count_dofs_and_fill_lst_dof_lids
    procedure, non_overridable, private         :: count_dofs_and_fill_lst_dof_lids_field          => coarse_fe_space_count_dofs_and_fill_lst_dof_lids_field
    procedure, non_overridable, private         :: free_lst_dofs_lids                              => coarse_fe_space_free_lst_dofs_lids
    procedure, non_overridable, private         :: free_number_dofs_per_field_and_block            => coarse_fe_space_free_number_dofs_per_field_and_block
    procedure, non_overridable, private         :: free_blocks_dof_import                          => coarse_fe_space_free_blocks_dof_import
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
    procedure, non_overridable, private         :: transfer_number_fields                          => coarse_fe_space_transfer_number_fields
    procedure, non_overridable, private         :: transfer_fe_space_type                          => coarse_fe_space_transfer_fe_space_type
    procedure, non_overridable, private         :: gather_ptr_dofs_per_fe_and_field                => coarse_fe_space_gather_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: gather_coarse_dofs_gids_rcv_counts_and_displs   => coarse_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
    procedure, non_overridable, private         :: gather_coarse_dofs_gids                         => coarse_fe_space_gather_coarse_dofs_gids
    procedure, non_overridable, private         :: gather_vefs_gids_dofs_objects                   => coarse_fe_space_gather_vefs_gids_dofs_objects

    procedure                                   :: get_number_fe_objects                           => coarse_fe_space_get_number_fe_objects
    procedure                                   :: get_total_number_coarse_dofs                    => coarse_fe_space_get_total_number_coarse_dofs
    procedure                                   :: get_block_number_coarse_dofs                    => coarse_fe_space_get_block_number_coarse_dofs
    procedure, non_overridable, private         :: free_coarse_dofs                                => coarse_fe_space_free_coarse_dofs
    
    procedure                                   :: renumber_dofs_first_interior_then_interface     => coarse_fe_space_renumber_dofs_first_interior_then_interface
    procedure                         , private :: renumber_dofs_block                             => coarse_fe_space_renumber_dofs_block
    
    ! Coarse FE traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_iterator                       => coarse_fe_space_create_coarse_fe_iterator
    procedure, non_overridable                  :: free_coarse_fe_iterator                         => coarse_fe_space_free_coarse_fe_iterator

    
    ! Coarse FE VEFS traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_vef_iterator                   => coarse_fe_space_create_coarse_fe_vef_iterator
    procedure, non_overridable                  :: create_itfc_coarse_fe_vef_iterator              => coarse_fe_space_create_itfc_coarse_fe_vef_iterator
    procedure, non_overridable                  :: free_coarse_fe_vef_iterator                     => coarse_fe_space_free_coarse_fe_vef_iterator
    
    
     ! Objects-related traversals
     procedure, non_overridable                 :: create_coarse_fe_object_iterator                => coarse_fe_space_create_coarse_fe_object_iterator
     procedure, non_overridable                 :: free_coarse_fe_object_iterator                  => coarse_fe_space_free_coarse_fe_object_iterator
     
     ! Getters
     procedure, non_overridable                 :: get_number_local_coarse_fes                     => coarse_fe_space_get_number_local_coarse_fes
     procedure, non_overridable                 :: get_number_ghost_coarse_fes                     => coarse_fe_space_get_number_ghost_coarse_fes
     procedure, non_overridable                 :: get_number_coarse_fe_objects                    => coarse_fe_space_get_number_coarse_fe_objects
     procedure, non_overridable                 :: get_triangulation                               => coarse_fe_space_get_triangulation
     procedure, non_overridable                 :: get_par_environment                             => coarse_fe_space_get_par_environment
     procedure                                  :: get_environment                                 => coarse_fe_space_get_environment
 end type coarse_fe_space_t
 
 public :: coarse_fe_space_t, coarse_fe_iterator_t, coarse_fe_vef_iterator_t
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
 
contains
!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_base_fe_space.i90"
#include "sbm_base_fe_iterator.i90"
#include "sbm_serial_fe_space.i90"
#include "sbm_fe_iterator.i90"
#include "sbm_fe_face_iterator.i90"
#include "sbm_base_fe_vef_iterator.i90"
#include "sbm_fe_vef_iterator.i90"
#include "sbm_par_fe_space.i90"
#include "sbm_base_fe_object_iterator.i90"
#include "sbm_fe_object_iterator.i90"
#include "sbm_standard_coarse_fe_handler.i90"
#include "sbm_H1_coarse_fe_handler.i90"
#include "sbm_Hcurl_coarse_fe_handler.i90"

#include "sbm_coarse_fe_space.i90"
#include "sbm_coarse_fe_object_iterator.i90"
#include "sbm_coarse_fe_iterator.i90"
#include "sbm_coarse_fe_vef_iterator.i90"

end module fe_space_names
