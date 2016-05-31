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
module serial_fe_space_names
  
 ! Serial modules
  use types_names
  use list_types_names
  use memor_names
  use sort_names
  use allocatable_array_names
  use hash_table_names
  
  use triangulation_names
  use conditions_names
  
  use reference_fe_names
  use field_names
  use function_names
  
  use matrix_names
  use vector_names
  use array_names
  use matrix_array_assembler_names
  
  use sparse_matrix_names
  use block_sparse_matrix_names
  use serial_scalar_array_names
  use serial_block_array_names
  use sparse_matrix_array_assembler_names
  use block_sparse_matrix_array_assembler_names
  use base_static_triangulation_names 
  
  
 ! Parallel modules
  use par_environment_names
  use par_context_names
  use par_triangulation_names
  use par_conditions_names
  use dof_import_names
  use element_import_names
  use par_sparse_matrix_names
  use par_scalar_array_names
  use par_sparse_matrix_array_assembler_names
  
  implicit none
# include "debug.i90"
  private

  ! This module includes the physical FE and FE space structure. It includes the
  ! following types:
  !
  ! * finite_element_t: It includes pointers to the unknown/geometrical reference_fe_t
  !   and the corresponding volume integrator, the local-to-global mapping for DOFs in
  !   elem2dof and the boundary conditions restricted to the element
  ! * fe_space_t: The abstract fe_space_t, which only includes some deferred methods 
  !   that must be provided by any fe_space_t concretization
  ! * serial_fe_space_t: The serial version of the fe_space_t, which includes an array
  !   of finite_element_t, the array of all possible reference_fe_t and volume_integrator_t
  !   in the FE space, and the total number of DOFs.
  ! * fe_function_scalar/vector/tensor_t: 3 different field-dependent objects that contain
  !   two work arrays at the element level: The first stores the nodal values of the FE  
  !   approximation from a global dof-vector. The second gives the values of the FE
  !   approximation at the quadrature points.

  type :: finite_element_t
     private 
     integer(ip)                                 :: number_nodes
     integer(ip)                                 :: number_fe_spaces
     
     integer(ip)                                 :: number_blocks
     integer(ip), pointer                        :: field_blocks(:)         => NULL()
     
     type(elem_topology_t)         , pointer     :: cell                    => NULL()
     type(fe_map_t)                , pointer     :: fe_map                  => NULL()
     
     type(p_reference_fe_t)        , pointer     :: reference_fe_phy(:)     => NULL()
     type(quadrature_t)            , pointer     :: quadrature              => NULL()
     type(p_volume_integrator_t)   , pointer     :: volume_integrator(:)    => NULL()
     
     type(i1p_t)                   , allocatable :: elem2dof(:)
     
     logical                       , allocatable :: at_strong_dirichlet_boundary(:)
     ! The next member variable points to the counterpart one in type(serial_fe_space_t)
     ! We would have like to add a pointer to the whole type(serial_fe_space_t) instead
     ! of to a particular member variable, but we could not as this is not Fortran standard
     ! conforming. This member variable is required to be able to impose strongly Dirichlet
     ! boundary conditions in the local scope of type(finite_element_t).
     type(serial_scalar_array_t)   , pointer     :: strong_dirichlet_values => NULL()
   contains
     
     procedure, non_overridable, private :: create =>  finite_element_create
     procedure, non_overridable          :: update_integration => finite_element_update_integration
     procedure, non_overridable, private :: free => finite_element_free
     ! procedure :: print Pending
     procedure, non_overridable, private :: fill_own_dofs => finite_element_fill_own_dofs
     procedure, non_overridable, private :: fill_own_dofs_on_vef => finite_element_fill_own_dofs_on_vef
     procedure, non_overridable, private :: fill_own_dofs_on_vef_from_source_element => finite_element_fill_own_dofs_on_vef_from_source_element
     procedure, non_overridable, private :: fill_dofs_face_integration_coupling => finite_element_fill_dofs_face_integration_coupling
     
     procedure, non_overridable :: get_number_nodes => finite_element_get_number_nodes
     procedure, non_overridable :: get_fe_map => finite_element_get_fe_map
     procedure, non_overridable :: get_quadrature => finite_element_get_quadrature
     procedure, non_overridable :: get_volume_integrator => finite_element_get_volume_integrator
     procedure, non_overridable :: get_elem2dof  => finite_element_get_elem2dof
     procedure, non_overridable :: get_number_nodes_per_field => finite_element_get_number_nodes_per_field
     procedure, non_overridable :: get_subset_id => finite_element_get_subset_id
     procedure, non_overridable :: get_order     => finite_element_get_order     
     procedure, non_overridable, private :: is_at_strong_dirichlet_boundary => finite_element_is_at_strong_dirichlet_boundary

     procedure, non_overridable :: compute_volume     => finite_element_compute_volume
     
     procedure, non_overridable, private :: update_scalar_values => finite_element_update_scalar_values
     procedure, non_overridable, private :: update_vector_values => finite_element_update_vector_values
     procedure, non_overridable, private :: update_tensor_values => finite_element_update_tensor_values
     generic :: update_values => update_scalar_values, &
                               & update_vector_values, &
                               & update_tensor_values
                               
     procedure, non_overridable :: impose_strong_dirichlet_bcs => finite_element_impose_strong_dirichlet_bcs
     procedure, non_overridable :: get_cell_coordinates => finite_element_get_cell_coordinates
  end type finite_element_t

  type :: p_finite_element_t
     type(finite_element_t), pointer :: p => NULL()
  end type p_finite_element_t

  public :: finite_element_t, p_finite_element_t

  type :: finite_face_t
     private
     integer(ip)                            :: number_fe_spaces
     type(face_topology_t)    , pointer     :: face_topology    =>   NULL()
     type(p_finite_element_t)               :: neighbour_fe(2)
     type(face_map_t)         , pointer     :: map              =>   NULL()
     type(quadrature_t)       , pointer     :: quadrature       =>   NULL() 
     type(p_face_integrator_t), allocatable :: face_integrator(:)
   contains
     procedure, non_overridable :: create              => finite_face_create
     procedure, non_overridable :: update_integration  => finite_face_update_integration
     procedure, non_overridable :: free                => finite_face_free
     procedure, non_overridable :: is_boundary         => finite_face_is_boundary
     procedure, non_overridable :: number_neighbours   => finite_face_number_neighbours
     procedure, non_overridable :: get_elem2dof        => finite_face_get_elem2dof
     procedure, non_overridable :: get_map             => finite_face_get_map
     procedure, non_overridable :: get_quadrature      => finite_face_get_quadrature
     procedure, non_overridable :: get_face_integrator => finite_face_get_face_integrator
     procedure, non_overridable :: get_relative_orientation => finite_face_get_relative_orientation
     procedure, non_overridable :: get_relative_rotation => finite_face_get_relative_rotation
  end type finite_face_t

  public :: finite_face_t

  
  integer(ip), parameter :: fe_space_type_cg            = 0 ! H^1 conforming FE space
  integer(ip), parameter :: fe_space_type_dg            = 1 ! L^2 conforming FE space + .not. H^1 conforming (weakly imposed via face integration)
  integer(ip), parameter :: fe_space_type_dg_conforming = 2 ! DG approximation of L^2 spaces (does not involve coupling by face)
  
  type :: serial_fe_space_t
     private
     integer(ip)                                 :: number_fe_spaces   
     type(p_fe_map_t)              , allocatable :: fe_map(:)
     type(face_map_t)              , allocatable :: face_map(:)
     type(p_reference_fe_t)        , allocatable :: reference_fe_phy_list(:)
     type(p_quadrature_t)          , allocatable :: quadrature(:)
     type(quadrature_t)            , allocatable :: face_quadrature(:)
     type(p_volume_integrator_t)   , allocatable :: volume_integrator(:)
     type(p_face_integrator_t)     , allocatable :: face_integrator(:)
     integer(ip)                   , allocatable :: fe_space_type(:)
     
     ! Strong Dirichlet data
     integer(ip)                 :: number_strong_dirichlet_dofs
     integer(ip), allocatable    :: strong_dirichlet_codes(:)
     type(serial_scalar_array_t) :: strong_dirichlet_values
     
     type(triangulation_t)         , pointer     :: triangulation =>  NULL()
     type(finite_element_t)        , allocatable :: fe_array(:)
     type(finite_face_t)           , allocatable :: face_array(:)
     
     ! Data related to block structure of the FE system + size of each block
     integer(ip)                                 :: number_blocks
     integer(ip)                   , allocatable :: field_blocks(:)
     logical                       , allocatable :: field_coupling(:,:)
     logical                       , allocatable :: blocks_coupling(:,:)
     integer(ip)                   , allocatable :: number_dofs_per_block(:)
     integer(ip)                   , allocatable :: number_dofs_per_field(:)
   contains
     procedure, private         :: serial_fe_space_create
     generic                    :: create => serial_fe_space_create
     procedure, non_overridable, private :: set_up_strong_dirichlet_bcs => serial_fe_space_set_up_strong_dirichlet_bcs
     procedure                  :: fill_dof_info => serial_fe_space_fill_dof_info
     procedure                 , private :: fill_elem2dof_and_count_dofs => serial_fe_space_fill_elem2dof_and_count_dofs
     procedure                  :: free => serial_fe_space_free
     procedure                  :: print => serial_fe_space_print
     procedure, non_overridable :: initialize_integration => serial_fe_space_initialize_integration
     procedure, non_overridable, private :: initialize_quadrature => serial_fe_space_initialize_quadrature
     procedure, non_overridable, private :: initialize_volume_integrator => serial_fe_space_initialize_volume_integrator
     procedure, non_overridable, private :: initialize_fe_map => serial_fe_space_initialize_fe_map
     procedure                  :: create_assembler => serial_fe_space_create_assembler
     procedure                  :: symbolic_setup_assembler => serial_fe_space_symbolic_setup_assembler
     procedure, non_overridable :: get_number_elements => serial_fe_space_get_number_elements
     procedure, non_overridable :: get_number_interior_faces => serial_fe_space_get_number_interior_faces
     procedure, non_overridable :: get_number_boundary_faces => serial_fe_space_get_number_boundary_faces
     procedure, non_overridable :: get_number_fe_spaces => serial_fe_space_get_number_fe_spaces
     procedure                  :: get_finite_element => serial_fe_space_get_finite_element
     procedure, non_overridable :: get_finite_face => serial_fe_space_get_finite_face
     procedure, non_overridable :: get_number_blocks => serial_fe_space_get_number_blocks
     procedure, non_overridable :: get_field_blocks => serial_fe_space_get_field_blocks
     procedure, non_overridable :: get_field_coupling => serial_fe_space_get_field_coupling
     procedure, private         :: get_max_number_nodes_field => serial_fe_space_get_max_number_nodes_field
     procedure, private         :: get_max_number_nodes_fe_space => serial_fe_space_get_max_number_nodes_fe_space
     generic :: get_max_number_nodes => get_max_number_nodes_field, & 
                                      & get_max_number_nodes_fe_space
     procedure, non_overridable :: get_max_number_quadrature_points => serial_fe_space_get_max_number_quadrature_points
     procedure, non_overridable :: create_face_array => serial_fe_space_create_face_array

     procedure :: update_bc_value_scalar => serial_fe_space_update_bc_value_scalar
     procedure :: update_bc_value_vector => serial_fe_space_update_bc_value_vector
     procedure :: update_bc_value_tensor => serial_fe_space_update_bc_value_tensor
     generic :: update_bc_value => update_bc_value_scalar, &
                                 & update_bc_value_vector, &
                                 & update_bc_value_tensor 
     
     procedure                  :: create_global_fe_function => serial_fe_space_create_global_fe_function
     procedure                  :: update_global_fe_function_bcs => serial_fe_space_update_global_fe_function_bcs
     procedure, private         :: create_fe_function_scalar => serial_fe_space_create_fe_function_scalar
     procedure, private         :: create_fe_function_vector => serial_fe_space_create_fe_function_vector
     procedure, private         :: create_fe_function_tensor => serial_fe_space_create_fe_function_tensor
     generic :: create_fe_function => create_fe_function_scalar, &
                                    & create_fe_function_vector, &
                                    & create_fe_function_tensor

     procedure, non_overridable :: interpolate_fe_function_scalar => serial_fe_space_interpolate_fe_function_scalar
     procedure, non_overridable :: interpolate_fe_function_vector => serial_fe_space_interpolate_fe_function_vector
     procedure, non_overridable :: interpolate_fe_function_tensor => serial_fe_space_interpolate_fe_function_tensor
     generic :: interpolate_fe_function => interpolate_fe_function_scalar, &
                                         & interpolate_fe_function_vector, &
                                         & interpolate_fe_function_tensor    

     procedure, non_overridable :: get_max_order => serial_fe_space_get_max_order
     procedure, non_overridable :: get_max_order_fe_space_component => serial_fe_space_get_max_order_fe_space_component
     procedure, non_overridable :: get_triangulation => serial_fe_space_get_triangulation
     procedure, non_overridable :: get_reference_fe_geo => serial_fe_space_get_reference_fe_geo
     procedure, non_overridable :: get_reference_fe_phy => serial_fe_space_get_reference_fe_phy
  end type serial_fe_space_t

  public :: serial_fe_space_t
    
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
    procedure, non_overridable, private :: get_scan_sum_number_dofs                   => coarse_fe_accessor_get_scan_sum_number_dofs
    procedure, non_overridable          :: get_number_fe_spaces                       => coarse_fe_accessor_get_number_fe_spaces
    procedure, non_overridable          :: get_number_dofs                            => coarse_fe_accessor_get_number_dofs
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
  
 
  type :: coarse_fe_space_t
    private
    integer(ip)                                 :: number_fields
    integer(ip) , allocatable                   :: field_type(:)
    integer(ip) , allocatable                   :: number_dofs_per_field(:)
    
    integer(ip)                                 :: number_blocks
    integer(ip)                   , allocatable :: field_blocks(:)
    logical                       , allocatable :: field_coupling(:,:)
    logical                       , allocatable :: blocks_coupling(:,:)
    integer(ip)                   , allocatable :: number_dofs_per_block(:)
    
    integer(ip) , allocatable                   :: ptr_dofs_per_fe_and_field(:)
    integer(ip) , allocatable                   :: lst_dofs_lids(:)
    type(list_t), allocatable                   :: own_dofs_vef_per_fe(:)
    
    type(dof_import_t)            , allocatable :: blocks_dof_import(:)
    
    type(coarse_triangulation_t), pointer       :: coarse_triangulation => NULL()
  contains
    procedure                                   :: create                                          => coarse_fe_space_create
    procedure                                   :: free                                            => coarse_fe_space_free
    procedure, non_overridable, private         :: allocate_and_fill_field_blocks_and_coupling     => coarse_fe_space_allocate_and_fill_field_blocks_and_coupling
    procedure, non_overridable, private         :: free_field_blocks_and_coupling                  => coarse_fe_space_free_field_blocks_and_coupling
    procedure, non_overridable, private         :: allocate_and_fill_field_type                    => coarse_fe_space_allocate_and_fill_field_type
    procedure, non_overridable, private         :: free_field_type                                 => coarse_fe_space_free_field_type
    procedure, non_overridable, private         :: allocate_and_fill_ptr_dofs_per_fe_and_field     => coarse_fe_space_allocate_and_fill_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: free_ptr_dofs_per_fe_and_field                  => coarse_fe_space_free_ptr_dofs_per_fe_and_field
    procedure, non_overridable, private         :: fetch_ghost_fes_data                            => coarse_fe_space_fetch_ghost_fes_data
    procedure, non_overridable, private         :: dof_gid2vef_gid                                 => coarse_fe_space_dof_gid2vef_gid
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
    ! Coarse FE traversals-related TBPs
    procedure, non_overridable                 :: create_coarse_fe_iterator                        => coarse_fe_space_create_coarse_fe_iterator
    
    ! Coarse FE VEFS traversals-related TBPs
    procedure, non_overridable                 :: create_coarse_fe_vef_iterator                    => coarse_fe_space_create_coarse_fe_vef_iterator
    procedure, non_overridable                 :: create_itfc_coarse_fe_vef_iterator               => coarse_fe_space_create_itfc_coarse_fe_vef_iterator
  end type coarse_fe_space_t
  
  ! These parameter constants are used in order to generate a unique (non-consecutive) 
  ! but consistent across MPI tasks global ID (integer(igp)) of a given DoF.
  ! See type(par_fe_space_t)%generate_non_consecutive_dof_gid()
  integer(ip), parameter :: cell_gid_shift              = 44
  integer(ip), parameter :: dofs_per_reference_fe_shift = 14
  integer(ip), parameter :: number_fields_shift         = igp*8-(cell_gid_shift+dofs_per_reference_fe_shift)-1
  
  type, extends(serial_fe_space_t) :: par_fe_space_t
     private
     type(par_triangulation_t), pointer     :: par_triangulation => NULL()
     type(finite_element_t)   , allocatable :: ghost_fe_array(:)
     type(dof_import_t)       , allocatable :: blocks_dof_import(:)
     
     ! GIDs of the DoFs objects. For their generation, we though to be a good idea 
     ! to take as basis the GIDs of the VEFs objects they are built from (instead of
     ! generating them from scratch, which in turn would imply further communication).
     integer(igp), allocatable               :: dofs_objects_gids(:)
     
     ! GIDs of the VEFs objects the DoFs objects are put on top of.
     integer(igp), allocatable               :: vefs_gids_dofs_objects(:)
     
     integer(ip) , allocatable               :: num_dofs_objects_per_field(:)
     type(list_t), allocatable               :: dofs_objects_per_field(:)
     
     type(coarse_fe_space_t), pointer        :: coarse_fe_space
  contains
     
     procedure, non_overridable, private :: par_fe_space_create
     procedure, non_overridable, private :: serial_fe_space_create                          => par_fe_space_serial_fe_space_create
     generic                             :: create                                          => par_fe_space_create
     procedure                           :: fill_dof_info                                   => par_fe_space_fill_dof_info
     procedure                 , private :: fill_elem2dof_and_count_dofs                    => par_fe_space_fill_elem2dof_and_count_dofs
     procedure, non_overridable, private :: compute_blocks_dof_import                       => par_fe_space_compute_blocks_dof_import
     procedure, non_overridable, private :: compute_dof_import                              => par_fe_space_compute_dof_import
     procedure, non_overridable, private :: compute_raw_interface_data_by_continuity        => par_fe_space_compute_raw_interface_data_by_continuity
     procedure, non_overridable, private :: raw_interface_data_by_continuity_decide_owner   => par_fe_space_raw_interface_data_by_continuity_decide_owner
     procedure, non_overridable, private :: compute_raw_interface_data_by_face_integ        => par_fe_space_compute_raw_interface_data_by_face_integ
     procedure, non_overridable, private :: compute_ubound_num_itfc_couplings_by_continuity => par_fe_space_compute_ubound_num_itfc_couplings_by_continuity
     procedure, non_overridable, private :: compute_ubound_num_itfc_couplings_by_face_integ => par_fe_space_compute_ubound_num_itfc_couplings_by_face_integ
     procedure, nopass, non_overridable, private :: generate_non_consecutive_dof_gid        => par_fe_space_generate_non_consecutive_dof_gid
     
     procedure                           :: get_finite_element            => par_fe_space_get_finite_element
     procedure                           :: print                         => par_fe_space_print
     procedure                           :: free                          => par_fe_space_free
     procedure                           :: create_assembler              => par_fe_space_create_assembler
     procedure                           :: symbolic_setup_assembler      => par_fe_space_symbolic_setup_assembler
     procedure                           :: update_bc_value_scalar        => par_fe_space_update_bc_value_scalar
     procedure                           :: update_bc_value_vector        => par_fe_space_update_bc_value_vector
     procedure                           :: update_bc_value_tensor        => par_fe_space_update_bc_value_tensor
     procedure                           :: create_global_fe_function     => par_fe_space_create_global_fe_function
     procedure                           :: update_global_fe_function_bcs => par_fe_space_update_global_fe_function_bcs
     
     ! Private methods for coarser FE space set-up
     procedure, non_overridable, private :: compute_dofs_objects                          => par_fe_space_compute_dofs_objects
     procedure, non_overridable, private :: compute_dofs_objects_by_continuity            => par_fe_space_compute_dofs_objects_by_continuity
     procedure, non_overridable, private :: setup_coarse_fe_space                         => par_fe_space_setup_coarse_fe_space
     procedure, non_overridable, private :: transfer_field_type                           => par_fe_space_transfer_field_type
     procedure, non_overridable, private :: gather_ptr_dofs_per_fe_and_field              => par_fe_space_gather_ptr_dofs_per_fe_and_field
     procedure, non_overridable, private :: gather_coarse_dofs_gids_rcv_counts_and_displs => par_fe_space_gather_coarse_dofs_gids_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_dofs_gids                       => par_fe_space_gather_coarse_dofs_gids
     procedure, non_overridable, private :: gather_vefs_gids_dofs_objects                 => par_fe_space_gather_vefs_gids_dofs_objects
  end type

  public :: par_fe_space_t
  
   type fe_function_t
   private
   class(vector_t), allocatable  :: dof_values
   type(serial_scalar_array_t)   :: strong_dirichlet_values
  contains
     procedure, non_overridable, private :: create                      => fe_function_create
     procedure, non_overridable, private :: copy_bc_values              => fe_function_copy_bc_values
     procedure, non_overridable          :: copy                        => fe_function_copy
     procedure, non_overridable          :: get_dof_values              => fe_function_get_dof_values
     procedure, non_overridable          :: get_strong_dirichlet_values => fe_function_get_strong_dirichlet_values
     procedure, non_overridable          :: free                        => fe_function_free
     generic                             :: assignment(=)               => copy
  end type fe_function_t 
  
  type fe_function_scalar_t
   private
   integer(ip) :: fe_space_id
   
   integer(ip) :: current_number_nodes             
   integer(ip) :: current_number_quadrature_points
   
   integer(ip) :: max_number_nodes  
   integer(ip) :: max_number_quadrature_points          
   
   real(rp), allocatable :: nodal_values(:)  
   real(rp), allocatable :: quadrature_points_values(:)
   
  contains
     procedure, non_overridable, private :: create                      => fe_function_scalar_create
     procedure, non_overridable :: get_fe_space_id                      => fe_function_scalar_get_fe_space_id
     procedure, non_overridable :: get_nodal_values                     => fe_function_scalar_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values         => fe_function_scalar_get_quadrature_points_values
     procedure, non_overridable :: get_value                            => fe_function_scalar_get_value
     procedure, non_overridable :: set_current_number_nodes             => fe_function_scalar_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points => fe_function_scalar_set_current_number_quadrature_points
     procedure, non_overridable :: free                                 => fe_function_scalar_free
  end type fe_function_scalar_t
  
  type fe_function_vector_t
   private
   integer(ip) :: fe_space_id
   
   integer(ip) :: current_number_nodes             
   integer(ip) :: current_number_quadrature_points           
   
   integer(ip) :: max_number_quadrature_points
   integer(ip) :: max_number_nodes            
   
   real(rp)            , allocatable :: nodal_values(:)  
   type(vector_field_t), allocatable :: quadrature_points_values(:)
   
  contains
     procedure, non_overridable, private :: create                      => fe_function_vector_create
     procedure, non_overridable :: get_fe_space_id                      => fe_function_vector_get_fe_space_id
     procedure, non_overridable :: get_nodal_values                     => fe_function_vector_get_nodal_values      
     procedure, non_overridable :: get_quadrature_points_values         => fe_function_vector_get_quadrature_points_values
     procedure, non_overridable :: get_value                            => fe_function_vector_get_value
     procedure, non_overridable :: set_current_number_nodes             => fe_function_vector_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points => fe_function_vector_set_current_number_quadrature_points
     procedure, non_overridable :: free                                 => fe_function_vector_free
  end type fe_function_vector_t
  
  type fe_function_tensor_t
   private
   integer(ip) :: fe_space_id
   
   integer(ip) :: current_number_nodes            
   integer(ip) :: current_number_quadrature_points     
   
   integer(ip) :: max_number_nodes    
   integer(ip) :: max_number_quadrature_points        
   
   real(rp)            , allocatable :: nodal_values(:)
   type(tensor_field_t), allocatable :: quadrature_points_values(:)
   
  contains
     procedure, non_overridable, private :: create                      => fe_function_tensor_create
     procedure, non_overridable :: get_fe_space_id                      => fe_function_tensor_get_fe_space_id
     procedure, non_overridable :: get_nodal_values                     => fe_function_tensor_get_nodal_values  
     procedure, non_overridable :: get_quadrature_points_values         => fe_function_tensor_get_quadrature_points_values          
     procedure, non_overridable :: get_value                            => fe_function_tensor_get_value
     procedure, non_overridable :: set_current_number_nodes             => fe_function_tensor_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points => fe_function_tensor_set_current_number_quadrature_points
     procedure, non_overridable :: free                                 => fe_function_tensor_free
  end type fe_function_tensor_t
  
  public :: fe_function_t, fe_function_scalar_t, fe_function_vector_t, fe_function_tensor_t
  
contains

  ! Includes with all the TBP and supporting subroutines for the types above.
  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_finite_element.i90"
#include "sbm_finite_face.i90"
#include "sbm_serial_fe_space.i90"
#include "../../Par/Integration/sbm_par_fe_space.i90"
#include "sbm_coarse_fe_space.i90"
#include "sbm_coarse_fe_accessor.i90"
#include "sbm_coarse_fe_iterator.i90"
#include "sbm_coarse_fe_vef_accessor.i90"
#include "sbm_coarse_fe_vef_iterator.i90"
#include "sbm_serial_fe_space_faces.i90"
#include "sbm_fe_function.i90"

end module serial_fe_space_names
