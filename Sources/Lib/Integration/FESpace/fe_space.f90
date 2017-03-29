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
  
  type, extends(cell_accessor_t) :: fe_accessor_t
    private
    class(serial_fe_space_t), pointer :: fe_space
  contains
    procedure, private, non_overridable :: fe_accessor_create
    generic                             :: create                                     => fe_accessor_create
    procedure                           :: cell_accessor_create                       => fe_accessor_cell_accessor_create
    procedure                           :: free                                       => fe_accessor_free
    procedure, non_overridable, private :: fill_own_dofs                              => fe_accessor_fill_own_dofs
    procedure, non_overridable, private :: fill_own_dofs_on_vef                       => fe_accessor_fill_own_dofs_on_vef
    procedure, non_overridable, private :: fill_own_dofs_on_vef_from_source_fe        => fe_accessor_fill_own_dofs_on_vef_from_source_fe
    procedure, non_overridable, private :: fill_own_strong_dirichlet_dofs_on_vef      => fe_accessor_fill_own_strong_dirichlet_dofs_on_vef
    procedure, non_overridable, private :: fill_dofs_face_integration_coupling        => fe_accessor_fill_dofs_face_integration_coupling
    procedure, non_overridable, private :: renumber_dofs_block                        => fe_accessor_renumber_dofs_block
    procedure, non_overridable, private :: renumber_dofs_field                        => fe_accessor_renumber_dofs_field
    procedure, non_overridable          :: update_integration                         => fe_accessor_update_integration

    procedure, non_overridable          :: get_fe_space                               => fe_accessor_get_fe_space
    procedure, non_overridable          :: get_number_fields                          => fe_accessor_get_number_fields
    procedure, non_overridable, private :: get_fe_space_type                          => fe_accessor_get_fe_space_type
    procedure, non_overridable          :: get_field_type                             => fe_accessor_get_field_type

    procedure, non_overridable          :: get_field_blocks                           => fe_accessor_get_field_blocks
    procedure, non_overridable          :: get_number_dofs                            => fe_accessor_get_number_dofs
    procedure, non_overridable          :: get_number_dofs_per_field                  => fe_accessor_get_number_dofs_per_field
    procedure, non_overridable          :: get_field_elem2dof                         => fe_accessor_get_field_elem2dof
    procedure, non_overridable          :: get_elem2dof                               => fe_accessor_get_elem2dof
    procedure, non_overridable          :: get_order                                  => fe_accessor_get_order
    
    procedure, non_overridable, private :: get_max_order_single_field                 => fe_accessor_get_max_order_single_field
    procedure, non_overridable, private :: get_max_order_all_fields                   => fe_accessor_get_max_order_all_fields
    generic                             :: get_max_order                              => get_max_order_single_field, &
                                                                                         get_max_order_all_fields

    procedure, non_overridable          :: at_strong_dirichlet_boundary               => fe_accessor_at_strong_dirichlet_boundary
    procedure, non_overridable          :: set_at_strong_dirichlet_boundary           => fe_accessor_set_at_strong_dirichlet_boundary
    procedure, non_overridable          :: unset_at_strong_dirichlet_boundary         => fe_accessor_unset_at_strong_dirichlet_boundary
    procedure, non_overridable          :: compute_volume                             => fe_accessor_compute_volume
    
    procedure, non_overridable          :: get_quadrature                             => fe_accessor_get_quadrature
    procedure, non_overridable          :: get_fe_map                                 => fe_accessor_get_fe_map
    procedure, non_overridable          :: get_volume_integrator                      => fe_accessor_get_volume_integrator    
    
    procedure, non_overridable, private :: fe_accessor_get_fe_vef
    generic                             :: get_vef                                    => fe_accessor_get_fe_vef
    procedure, non_overridable          :: get_reference_fe                           => fe_accessor_get_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe                 => fe_accessor_get_max_order_reference_fe
    procedure, non_overridable          :: get_max_order_reference_fe_id              => fe_accessor_get_max_order_reference_fe_id
    procedure, non_overridable          :: get_reference_fe_id                        => fe_accessor_get_reference_fe_id
    procedure, non_overridable          :: create_own_dofs_on_vef_iterator            => fe_accessor_create_own_dofs_on_vef_iterator
    procedure, non_overridable          :: impose_strong_dirichlet_bcs                => fe_accessor_impose_strong_dirichlet_bcs
  end type fe_accessor_t
  
  type fe_iterator_t
    private
    type(fe_accessor_t) :: fe_accessor
  contains
    procedure, non_overridable, private :: create       => fe_iterator_create
    procedure, non_overridable          :: free         => fe_iterator_free
    procedure, non_overridable          :: init         => fe_iterator_init
    procedure, non_overridable          :: next         => fe_iterator_next
    procedure, non_overridable          :: has_finished => fe_iterator_has_finished
    procedure, non_overridable, private :: fe_iterator_current
    generic                             :: current      => fe_iterator_current
  end type fe_iterator_t
  
  type, extends(vef_accessor_t) :: fe_vef_accessor_t
    private
    class(serial_fe_space_t), pointer :: fe_space
  contains
     procedure                 , private :: fe_vef_accessor_create
     procedure                           :: vef_accessor_create               => fe_vef_accessor_vef_accessor_create
     generic                             :: create                            => fe_vef_accessor_create
     procedure                           :: free                              => fe_vef_accessor_free
     procedure, non_overridable, private :: fe_vef_accessor_get_cell_around
     generic                             :: get_cell_around                   => fe_vef_accessor_get_cell_around
  end type fe_vef_accessor_t
  
  type fe_vef_iterator_t
    private
    type(fe_vef_accessor_t) :: current_fe_vef_accessor
  contains
     procedure, non_overridable, private :: create       => fe_vef_iterator_create
     procedure, non_overridable          :: free         => fe_vef_iterator_free
     procedure, non_overridable          :: init         => fe_vef_iterator_init
     procedure, non_overridable          :: next         => fe_vef_iterator_next
     procedure, non_overridable          :: has_finished => fe_vef_iterator_has_finished
     procedure, non_overridable, private :: fe_vef_iterator_current
     generic                             :: current      => fe_vef_iterator_current
  end type fe_vef_iterator_t  
  
  type :: itfc_fe_vef_iterator_t
    private
    type(itfc_vef_iterator_t) :: itfc_vef_iterator
    type(fe_vef_accessor_t)   :: current_vef_accessor
  contains
    procedure, non_overridable, private :: create       => itfc_fe_vef_iterator_create
    procedure, non_overridable          :: free         => itfc_fe_vef_iterator_free
    procedure, non_overridable          :: init         => itfc_fe_vef_iterator_init
    procedure, non_overridable          :: next         => itfc_fe_vef_iterator_next
    procedure, non_overridable          :: has_finished => itfc_fe_vef_iterator_has_finished
    procedure, non_overridable, private :: itfc_fe_vef_iterator_current
    generic                             :: current      => itfc_fe_vef_iterator_current
  end type itfc_fe_vef_iterator_t
  
  
  type, extends(face_accessor_t) :: fe_face_accessor_t
    private
    class(serial_fe_space_t), pointer :: fe_space
  contains
    procedure, private, non_overridable :: fe_face_accessor_create
    generic                             :: create                                     => fe_face_accessor_create
    procedure                           :: vef_accessor_create                        => fe_face_accessor_vef_accessor_create
    procedure                           :: free                                       => fe_face_accessor_free
    generic                             :: get_cell_around                            => fe_face_accessor_get_fe_around
    procedure, private, non_overridable :: fe_face_accessor_get_fe_around
    procedure         , non_overridable :: update_integration                         => fe_face_accessor_update_integration 
    procedure         , non_overridable :: get_fe_space                               => fe_face_accessor_get_fe_space
    procedure         , non_overridable :: get_elem2dof                               => fe_face_accessor_get_elem2dof
    procedure         , non_overridable :: get_quadrature                             => fe_face_accessor_get_quadrature
    procedure         , non_overridable :: get_face_map                               => fe_face_accessor_get_face_map
    procedure         , non_overridable :: get_face_integrator                        => fe_face_accessor_get_face_integrator
    procedure         , non_overridable :: impose_strong_dirichlet_bcs                => fe_face_accessor_impose_strong_dirichlet_bcs
    procedure         , non_overridable :: compute_surface                            => fe_face_accessor_compute_surface
  end type fe_face_accessor_t
  
  
  type :: fe_face_iterator_t
    private
    type(face_iterator_t)     :: face_iterator
    type(fe_face_accessor_t)  :: current_fe_face_accessor
  contains
    procedure, non_overridable, private :: create       => fe_face_iterator_create
    procedure, non_overridable          :: free         => fe_face_iterator_free
    procedure, non_overridable          :: init         => fe_face_iterator_init
    procedure, non_overridable          :: next         => fe_face_iterator_next
    procedure, non_overridable          :: has_finished => fe_face_iterator_has_finished
    procedure, non_overridable, private :: fe_face_iterator_current
    generic                             :: current      => fe_face_iterator_current
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
     
     ! fes and fe_faces traversals-related TBPs
     procedure, non_overridable          :: create_fe_iterator                           => serial_fe_space_create_fe_iterator
     procedure, non_overridable          :: create_fe_vef_iterator                       => serial_fe_space_create_fe_vef_iterator
     procedure, non_overridable          :: create_itfc_fe_vef_iterator                  => serial_fe_space_create_itfc_fe_vef_iterator
     procedure, non_overridable          :: create_fe_face_iterator                      => serial_fe_space_create_fe_face_iterator
 end type serial_fe_space_t  
 
 public :: serial_fe_space_t   
 public :: fe_iterator_t, fe_accessor_t
 public :: fe_vef_iterator_t, itfc_fe_vef_iterator_t, fe_vef_accessor_t
 public :: fe_face_iterator_t, fe_face_accessor_t
 
 type, extends(object_accessor_t) :: fe_object_accessor_t
    private
    type(par_fe_space_t), pointer :: fe_space
  contains
    procedure                            :: fe_object_accessor_create
    generic                              :: create                                => fe_object_accessor_create
    procedure                            :: free                                  => fe_object_accessor_free
    procedure, non_overridable           :: get_number_fe_vefs_on_object          => fe_object_accessor_get_number_fe_vefs_on_object
    procedure, non_overridable           :: create_fe_vefs_on_object_iterator     => fe_object_accessor_create_fe_vefs_on_object_iterator
    procedure, non_overridable           :: create_fe_faces_on_object_iterator    => fe_object_accessor_create_fe_faces_on_object_iterator
    procedure, non_overridable           :: get_number_coarse_dofs                => fe_object_accessor_get_number_coarse_dofs
    procedure, non_overridable           :: create_own_coarse_dofs_iterator       => fe_object_accessor_create_own_coarse_dofs_iterator
  end type fe_object_accessor_t
  
  type fe_object_iterator_t
    private
    type(fe_object_accessor_t) :: current_fe_object_accessor
  contains
    procedure, non_overridable          :: create       => fe_object_iterator_create
    procedure, non_overridable          :: free         => fe_object_iterator_free
    procedure, non_overridable          :: init         => fe_object_iterator_init
    procedure, non_overridable          :: next         => fe_object_iterator_next
    procedure, non_overridable          :: has_finished => fe_object_iterator_has_finished
    procedure, non_overridable          :: current      => fe_object_iterator_current
  end type fe_object_iterator_t
  
  type :: fe_vefs_on_object_iterator_t
    private
    type(vefs_on_object_iterator_t)     :: vefs_on_object_iterator
    type(fe_vef_accessor_t)             :: current_fe_vef_accessor
  contains
    procedure, non_overridable          :: create       => fe_vefs_on_object_iterator_create
    procedure, non_overridable          :: free         => fe_vefs_on_object_iterator_free
    procedure, non_overridable          :: init         => fe_vefs_on_object_iterator_init
    procedure, non_overridable          :: next         => fe_vefs_on_object_iterator_next
    procedure, non_overridable          :: has_finished => fe_vefs_on_object_iterator_has_finished
    procedure, non_overridable          :: current      => fe_vefs_on_object_iterator_current
  end type fe_vefs_on_object_iterator_t
  
  type :: fe_faces_on_object_iterator_t
    private
    type(vefs_on_object_iterator_t) :: vefs_on_object_iterator
    type(fe_face_accessor_t)        :: current_fe_face_accessor
  contains
    procedure, non_overridable          :: create       => fe_faces_on_object_iterator_create
    procedure, non_overridable          :: free         => fe_faces_on_object_iterator_free
    procedure, non_overridable          :: init         => fe_faces_on_object_iterator_init
    procedure, non_overridable          :: next         => fe_faces_on_object_iterator_next
    procedure, non_overridable          :: has_finished => fe_faces_on_object_iterator_has_finished
    procedure, non_overridable          :: current      => fe_faces_on_object_iterator_current
  end type fe_faces_on_object_iterator_t

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
   procedure, private :: par_fe_space_create_same_reference_fes_on_all_cells 
   generic                                     :: create                                          => par_fe_space_create_same_reference_fes_on_all_cells   
   procedure                                   :: fill_dof_info                                   => par_fe_space_fill_dof_info
   procedure                         , private :: fill_elem2dof_and_count_dofs                    => par_fe_space_fill_elem2dof_and_count_dofs
   procedure                                   :: renumber_dofs_first_interior_then_interface     => par_fe_space_renumber_dofs_first_interior_then_interface

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
   
   procedure       , non_overridable, private  :: setup_coarse_dofs                               => par_fe_space_setup_coarse_dofs
   procedure       , non_overridable, private  :: free_coarse_dofs                                => par_fe_space_free_coarse_dofs
   procedure       , non_overridable           :: setup_coarse_fe_space                           => par_fe_space_setup_coarse_fe_space
   procedure       , non_overridable, private  :: transfer_number_fields                          => par_fe_space_transfer_number_fields
   procedure       , non_overridable, private  :: transfer_fe_space_type                          => par_fe_space_transfer_fe_space_type
   procedure       , non_overridable, private  :: gather_ptr_dofs_per_fe_and_field                => par_fe_space_gather_ptr_dofs_per_fe_and_field
   procedure       , non_overridable, private  :: gather_coarse_dofs_gids_rcv_counts_and_displs   => par_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
   procedure       , non_overridable, private  :: gather_coarse_dofs_gids                         => par_fe_space_gather_coarse_dofs_gids
   procedure       , non_overridable, private  :: gather_vefs_gids_dofs_objects                   => par_fe_space_gather_vefs_gids_dofs_objects
   procedure                                   :: get_total_number_coarse_dofs                    => par_fe_space_get_total_number_coarse_dofs
   procedure                                   :: get_block_number_coarse_dofs                    => par_fe_space_get_block_number_coarse_dofs
   procedure       , non_overridable           :: get_coarse_fe_handler                           => par_fe_space_get_coarse_fe_handler

   ! Objects-related traversals
   procedure, non_overridable                  :: create_fe_object_iterator                       => par_fe_space_create_fe_object_iterator
   procedure, non_overridable                  :: create_fe_vefs_on_object_iterator               => par_fe_space_create_fe_vefs_on_object_iterator
 end type par_fe_space_t
 
 public :: par_fe_space_t
 public :: fe_object_accessor_t, fe_object_iterator_t, fe_vefs_on_object_iterator_t, fe_faces_on_object_iterator_t
 
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
    subroutine l1_get_num_coarse_dofs_interface(this, par_fe_space, parameter_list, num_coarse_dofs) 
      import :: l1_coarse_fe_handler_t, par_fe_space_t, parameterlist_t, ip
      class(l1_coarse_fe_handler_t), intent(in)    :: this
      type(par_fe_space_t)         , intent(in)    :: par_fe_space 
      type(parameterlist_t)        , intent(in)    :: parameter_list
      integer(ip)                  , intent(inout) :: num_coarse_dofs(:)
    end subroutine l1_get_num_coarse_dofs_interface
   
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
    real(rp), public :: diffusion_inclusion
  contains
	   procedure :: setup_constraint_matrix  => H1_l1_setup_constraint_matrix
	   procedure :: setup_weighting_operator => H1_l1_setup_weighting_operator
  end type H1_l1_coarse_fe_handler_t
  
  public :: l1_coarse_fe_handler_t, standard_l1_coarse_fe_handler_t, H1_l1_coarse_fe_handler_t
    
  type, extends(cell_accessor_t) :: coarse_fe_accessor_t
    private
    type(coarse_fe_space_t), pointer :: coarse_fe_space
  contains
    procedure, private, non_overridable :: coarse_fe_accessor_create
    generic                             :: create                                     => coarse_fe_accessor_create
    procedure                           :: cell_accessor_create                       => coarse_fe_accessor_cell_accessor_create
    procedure                           :: free                                       => coarse_fe_accessor_free
    procedure, non_overridable, private :: create_own_dofs_on_vef_iterator            => coarse_fe_accessor_create_own_dofs_on_vef_iterator
    procedure, non_overridable, private :: fill_own_dofs_on_vef                       => coarse_fe_accessor_fill_own_dofs_on_vef
    procedure, non_overridable, private :: fill_own_dofs_on_vef_from_source_coarse_fe => coarse_fe_accessor_fill_own_dofs_on_vef_from_source_coarse_fe
    procedure, non_overridable, private :: renumber_dofs_block                        => coarse_fe_accessor_renumber_dofs_block
    procedure, non_overridable, private :: renumber_dofs_field                        => coarse_fe_accessor_renumber_dofs_field
    
    procedure, non_overridable, private :: get_scan_sum_number_dofs                   => coarse_fe_accessor_get_scan_sum_number_dofs
    procedure, non_overridable          :: get_number_fe_spaces                       => coarse_fe_accessor_get_number_fe_spaces
    procedure, non_overridable          :: get_number_dofs                            => coarse_fe_accessor_get_number_dofs
    procedure, non_overridable          :: get_elem2dof                               => coarse_fe_accessor_get_elem2dof
    procedure, non_overridable          :: get_coarse_fe_vef                          => coarse_fe_accessor_get_coarse_fe_vef
  end type coarse_fe_accessor_t
  
  type coarse_fe_iterator_t
    private
    type(coarse_fe_accessor_t) :: current_coarse_fe_accessor
  contains
    procedure, non_overridable, private :: create       => coarse_fe_iterator_create
    procedure, non_overridable          :: free         => coarse_fe_iterator_free
    procedure, non_overridable          :: init         => coarse_fe_iterator_init
    procedure, non_overridable          :: next         => coarse_fe_iterator_next
    procedure, non_overridable          :: has_finished => coarse_fe_iterator_has_finished
    procedure, non_overridable          :: current      => coarse_fe_iterator_current
  end type coarse_fe_iterator_t
  
  type, extends(vef_accessor_t) :: coarse_fe_vef_accessor_t
    private
    type(coarse_fe_space_t), pointer            :: coarse_fe_space
  contains
     procedure :: coarse_fe_vef_accessor_create
     procedure :: vef_accessor_create                                 => coarse_fe_vef_accessor_vef_accessor_create
     generic   :: create                                              => coarse_fe_vef_accessor_create
     procedure :: free                                                => coarse_fe_vef_accessor_free
     procedure, non_overridable          :: get_num_coarse_fes_around => coarse_fe_vef_accessor_get_num_coarse_fes_around
     procedure, non_overridable          :: get_coarse_fe_around      => coarse_fe_vef_accessor_get_coarse_fe_around
  end type coarse_fe_vef_accessor_t
  
  type coarse_fe_vef_iterator_t
    private
    type(coarse_fe_vef_accessor_t) :: current_coarse_fe_vef_accessor
  contains
     procedure, non_overridable, private :: create       => coarse_fe_vef_iterator_create
     procedure, non_overridable          :: free         => coarse_fe_vef_iterator_free
     procedure, non_overridable          :: init         => coarse_fe_vef_iterator_init
     procedure, non_overridable          :: next         => coarse_fe_vef_iterator_next
     procedure, non_overridable          :: has_finished => coarse_fe_vef_iterator_has_finished
     procedure, non_overridable          :: current      => coarse_fe_vef_iterator_current
  end type coarse_fe_vef_iterator_t  
  
  type :: itfc_coarse_fe_vef_iterator_t
    private
    type(itfc_vef_iterator_t)           :: itfc_vef_iterator
    type(coarse_fe_vef_accessor_t)      :: current_coarse_vef_accessor
  contains
    procedure, non_overridable, private :: create       => itfc_coarse_fe_vef_iterator_create
    procedure, non_overridable          :: free         => itfc_coarse_fe_vef_iterator_free
    procedure, non_overridable          :: init         => itfc_coarse_fe_vef_iterator_init
    procedure, non_overridable          :: next         => itfc_coarse_fe_vef_iterator_next
    procedure, non_overridable          :: has_finished => itfc_coarse_fe_vef_iterator_has_finished
    procedure, non_overridable          :: current      => itfc_coarse_fe_vef_iterator_current
  end type itfc_coarse_fe_vef_iterator_t
  
  
  type, extends(object_accessor_t) :: coarse_fe_object_accessor_t
    private
    type(coarse_fe_space_t), pointer :: coarse_fe_space
  contains
    procedure                            :: coarse_fe_object_accessor_create
    generic                              :: create                                   => coarse_fe_object_accessor_create
    procedure                            :: free                                     => coarse_fe_object_accessor_free
    procedure, non_overridable           :: get_number_coarse_fe_vefs_on_object      => coarse_fe_object_accessor_get_number_coarse_fe_vefs_on_object
    procedure, non_overridable           :: create_coarse_fe_vefs_on_object_iterator => coarse_fe_object_accessor_get_coarse_fe_vefs_on_object_iterator
  end type coarse_fe_object_accessor_t
  
  type coarse_fe_object_iterator_t
    private
    type(coarse_fe_object_accessor_t) :: current_coarse_fe_object_accessor
  contains
     procedure, non_overridable          :: create       => coarse_fe_object_iterator_create
     procedure, non_overridable          :: free         => coarse_fe_object_iterator_free
     procedure, non_overridable          :: init         => coarse_fe_object_iterator_init
     procedure, non_overridable          :: next         => coarse_fe_object_iterator_next
     procedure, non_overridable          :: has_finished => coarse_fe_object_iterator_has_finished
     procedure, non_overridable          :: current      => coarse_fe_object_iterator_current
  end type coarse_fe_object_iterator_t
  
  type :: coarse_fe_vefs_on_object_iterator_t
    private
    type(vefs_on_object_iterator_t)     :: vefs_on_object_iterator
    type(coarse_fe_vef_accessor_t)      :: current_coarse_fe_vef_accessor
  contains
    procedure, non_overridable          :: create       => coarse_fe_vefs_on_object_iterator_create
    procedure, non_overridable          :: free         => coarse_fe_vefs_on_object_iterator_free
    procedure, non_overridable          :: init         => coarse_fe_vefs_on_object_iterator_init
    procedure, non_overridable          :: next         => coarse_fe_vefs_on_object_iterator_next
    procedure, non_overridable          :: has_finished => coarse_fe_vefs_on_object_iterator_has_finished
    procedure, non_overridable          :: current      => coarse_fe_vefs_on_object_iterator_current
  end type coarse_fe_vefs_on_object_iterator_t
  
 type, extends(coarse_fe_object_accessor_t) :: coarse_dof_object_accessor_t
   private
   integer(ip)                      :: field_id
   integer(ip)                      :: coarse_dof_object_lid
 contains
    procedure  :: coarse_dof_object_accessor_create
    generic    :: create                      => coarse_dof_object_accessor_create
    procedure  :: free                        => coarse_dof_object_accessor_free
    procedure  :: next                        => coarse_dof_object_accessor_next
    procedure  :: set_lid                     => coarse_dof_object_accessor_set_lid
    procedure  :: get_number_dofs_on_object   => coarse_dof_object_accessor_get_number_dofs_on_object
    procedure  :: get_dofs_on_object_iterator => coarse_dof_object_accessor_get_dofs_on_object_iterator
 end type coarse_dof_object_accessor_t
 
 type coarse_dof_object_iterator_t
    private
    type(coarse_dof_object_accessor_t) :: current_coarse_dof_object_accessor
  contains
     procedure, non_overridable          :: create       => coarse_dof_object_iterator_create
     procedure, non_overridable          :: free         => coarse_dof_object_iterator_free
     procedure, non_overridable          :: init         => coarse_dof_object_iterator_init
     procedure, non_overridable          :: next         => coarse_dof_object_iterator_next
     procedure, non_overridable          :: has_finished => coarse_dof_object_iterator_has_finished
     procedure, non_overridable          :: current      => coarse_dof_object_iterator_current
  end type coarse_dof_object_iterator_t
  
  
  type, extends(base_fe_space_t) :: coarse_fe_space_t
    private
    integer(ip) , allocatable                   :: ptr_dofs_per_fe_and_field(:)
    integer(ip) , allocatable                   :: lst_dofs_lids(:)
    type(list_t), allocatable                   :: own_dofs_vef_per_fe(:)
    	
	   ! Pointer to coarse triangulation this coarse_fe_space has been built from
    type(coarse_triangulation_t), pointer       :: coarse_triangulation => NULL()
	
	   ! Multilevel fe space   
    integer(ip)                   , allocatable :: num_coarse_dofs_per_field(:)
   	
    ! Aggregation of fine DoFs into coarse DoFs
    type(list_t)                  , allocatable :: fine_dofs_coarse_dofs_per_field(:) 
    ! GIDs of the coarse DoFs. For their generation, we though to be a good idea 
    ! to take as basis the GIDs of the VEFs objects they are built from (instead of
    ! generating them from scratch, which in turn would imply further communication).
    type(allocatable_array_igp1_t), allocatable :: coarse_dof_gids_per_field(:)
    ! LIDs of the VEFs objects the coarse DoFs are put on top of
    type(allocatable_array_ip1_t) , allocatable :: coarse_n_face_lids_coarse_dofs_per_field(:)
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
    procedure, non_overridable, private         :: setup_dofs_objects_by_continuity                => coarse_fe_space_setup_dofs_objects_by_continuity
    procedure, non_overridable                  :: setup_constraint_matrix                         => coarse_fe_space_setup_constraint_matrix
    procedure, non_overridable, private         :: setup_coarse_fe_space                           => coarse_fe_space_setup_coarse_fe_space
    procedure, non_overridable, private         :: transfer_number_fields                          => coarse_fe_space_transfer_number_fields
    procedure, non_overridable, private         :: transfer_fe_space_type                          => coarse_fe_space_transfer_fe_space_type
    procedure, non_overridable, private         :: gather_ptr_dofs_per_fe_and_field                => coarse_fe_space_gather_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: gather_coarse_dofs_gids_rcv_counts_and_displs   => coarse_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
    procedure, non_overridable, private         :: gather_coarse_dofs_gids                         => coarse_fe_space_gather_coarse_dofs_gids
    procedure, non_overridable, private         :: gather_vefs_gids_dofs_objects                   => coarse_fe_space_gather_vefs_gids_dofs_objects
    procedure                                   :: get_total_number_coarse_dofs                    => coarse_fe_space_get_total_number_coarse_dofs
    procedure                                   :: get_block_number_coarse_dofs                    => coarse_fe_space_get_block_number_coarse_dofs
	   procedure, non_overridable, private         :: free_coarse_dofs                                => coarse_fe_space_free_coarse_dofs
    
    procedure                                   :: renumber_dofs_first_interior_then_interface     => coarse_fe_space_renumber_dofs_first_interior_then_interface
    procedure                         , private :: renumber_dofs_block                             => coarse_fe_space_renumber_dofs_block
    
    ! Coarse FE traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_iterator                       => coarse_fe_space_create_coarse_fe_iterator
    
    ! Coarse FE VEFS traversals-related TBPs
    procedure, non_overridable                  :: create_coarse_fe_vef_iterator                   => coarse_fe_space_create_coarse_fe_vef_iterator
    procedure, non_overridable                  :: create_itfc_coarse_fe_vef_iterator              => coarse_fe_space_create_itfc_coarse_fe_vef_iterator
    
     ! Objects-related traversals
     procedure, non_overridable                 :: create_coarse_fe_object_iterator                => coarse_fe_space_create_coarse_fe_object_iterator
     procedure, non_overridable                 :: create_coarse_fe_vefs_on_object_iterator        => coarse_fe_space_create_coarse_fe_vefs_on_object_iterator
     procedure, non_overridable                 :: create_field_dofs_object_iterator               => coarse_fe_space_create_field_dofs_object_iterator
     
     ! Getters
     procedure, non_overridable                 :: get_number_local_coarse_fes                     => coarse_fe_space_get_number_local_coarse_fes
     procedure, non_overridable                 :: get_number_ghost_coarse_fes                     => coarse_fe_space_get_number_ghost_coarse_fes
     procedure, non_overridable                 :: get_number_coarse_fe_objects                    => coarse_fe_space_get_number_coarse_fe_objects
     procedure, non_overridable                 :: get_triangulation                               => coarse_fe_space_get_triangulation
     procedure, non_overridable                 :: get_par_environment                             => coarse_fe_space_get_par_environment
     procedure                                  :: get_environment                                 => coarse_fe_space_get_environment
 end type coarse_fe_space_t
 
 public :: coarse_fe_space_t, coarse_fe_iterator_t, coarse_fe_accessor_t
 public :: coarse_fe_object_accessor_t, coarse_dof_object_accessor_t, coarse_dof_object_iterator_t
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
#include "sbm_serial_fe_space.i90"
#include "sbm_fe_accessor.i90"
#include "sbm_fe_iterator.i90"
#include "sbm_fe_face_accessor.i90"
#include "sbm_fe_face_iterator.i90"
#include "sbm_fe_vef_accessor.i90"
#include "sbm_fe_vef_iterator.i90"
#include "sbm_par_fe_space.i90"
#include "sbm_fe_object_accessor.i90"
#include "sbm_fe_object_iterator.i90"
#include "sbm_fe_vefs_on_object_iterator.i90"
#include "sbm_fe_faces_on_object_iterator.i90"
#include "sbm_standard_coarse_fe_handler.i90"
#include "sbm_H1_coarse_fe_handler.i90"

#include "sbm_coarse_fe_space.i90"
#include "sbm_coarse_fe_object_accessor.i90"
#include "sbm_coarse_fe_object_iterator.i90"
#include "sbm_coarse_fe_vefs_on_object_iterator.i90"
#include "sbm_coarse_fe_accessor.i90"
#include "sbm_coarse_fe_iterator.i90"
#include "sbm_coarse_fe_vef_accessor.i90"
#include "sbm_coarse_fe_vef_iterator.i90"
#include "sbm_coarse_dof_object_accessor.i90"
#include "sbm_coarse_dof_object_iterator.i90"

end module fe_space_names
