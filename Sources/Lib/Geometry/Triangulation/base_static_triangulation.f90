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
module base_static_triangulation_names
  ! Serial modules
  use types_names
  use memor_names
  use sort_names
  use reference_fe_names
  use field_names
  use cell_import_names
  use hash_table_names
  use list_types_names
  use mesh_names
  use mesh_distribution_names
  use par_io_names
  use stdio_names
  use field_names
  use FPL
  use uniform_hex_mesh_generator_names
  
  ! Geometry modules
  use sisl_names
  use geometry_names
  use gid_geometry_reader_names

  ! Par modules
  use environment_names

  implicit none
# include "debug.i90"
  private

  ! Parameters to define vef_type. Observe that
  ! mod(vef_type,10)     = dimension
  ! mod(vef_type/10,10)  = interior (0) or boundary (1)
  ! vef_type/100         = local (0), interface(1) or ghost(2)
  integer(ip), parameter :: local_interior_dim0 =  0
  integer(ip), parameter :: local_interior_dim1 =  1
  integer(ip), parameter :: local_interior_dim2 =  2
  integer(ip), parameter :: local_boundary_dim0 = 10
  integer(ip), parameter :: local_boundary_dim1 = 11
  integer(ip), parameter :: local_boundary_dim2 = 12

  integer(ip), parameter :: interface_interior_dim0 = 100
  integer(ip), parameter :: interface_interior_dim1 = 101
  integer(ip), parameter :: interface_interior_dim2 = 102
  integer(ip), parameter :: interface_boundary_dim0 = 110
  integer(ip), parameter :: interface_boundary_dim1 = 111
  integer(ip), parameter :: interface_boundary_dim2 = 112

  integer(ip), parameter :: ghost_dim0 = 200
  integer(ip), parameter :: ghost_dim1 = 201
  integer(ip), parameter :: ghost_dim2 = 202
  
  integer(ip), parameter :: triangulation_generate_from_mesh  = 0
  integer(ip), parameter :: triangulation_generate_structured = 1
  public :: triangulation_generate_from_mesh
  public :: triangulation_generate_structured
  
  character(len=*), parameter :: geometry_interpolation_order_key    = 'geometry_interpolation_order'
  character(len=*), parameter :: triangulation_generate_key          = 'triangulation_generate'
  public :: geometry_interpolation_order_key
  public :: triangulation_generate_key
  
  type cell_iterator_t
    private
    integer(ip)                                 :: gid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
    ! create/free/next/set_gid CANNOT longer be private as type(coarse_fe_cell_iterator_t)
    ! extends type(cell_iterator_t), requires to call these TBPs and cannot be placed
    ! either within this module or within a submodule of base_static_triangulation_names
    ! due to module dependencies cycle
    procedure                            :: create                  => cell_iterator_create
    procedure                            :: free                    => cell_iterator_free
    final                                ::                            cell_iterator_free_final
    procedure, non_overridable           :: next                    => cell_iterator_next
    procedure, non_overridable           :: first                   => cell_iterator_first
    procedure, non_overridable           :: last                    => cell_iterator_last
    procedure, non_overridable           :: set_gid                 => cell_iterator_set_gid
    procedure, non_overridable, private  :: set_ggid                 => cell_iterator_set_ggid
    procedure, non_overridable, private  :: set_mypart              => cell_iterator_set_mypart
    procedure, non_overridable, private  :: get_triangulation       => cell_iterator_get_triangulation
    procedure, non_overridable, private  :: cell_iterator_get_vef
    procedure, non_overridable           :: has_finished            => cell_iterator_has_finished
    procedure, non_overridable           :: get_reference_fe_geo    => cell_iterator_get_reference_fe_geo
    procedure, non_overridable           :: get_reference_fe_geo_id => cell_iterator_get_reference_fe_geo_id
    procedure, non_overridable           :: get_coordinates         => cell_iterator_get_coordinates
    procedure, non_overridable           :: set_coordinates         => cell_iterator_set_coordinates
    procedure, non_overridable           :: get_coordinates_ref_space => cell_iterator_get_coordinates_ref_space
    procedure, non_overridable           :: get_gid                 => cell_iterator_get_gid
    procedure, non_overridable           :: get_ggid                 => cell_iterator_get_ggid
    procedure, non_overridable           :: get_my_part             => cell_iterator_get_mypart
    procedure, non_overridable           :: get_my_subpart          => cell_iterator_get_mysubpart
    procedure, non_overridable           :: get_my_subpart_lid      => cell_iterator_get_mysubpart_lid
    procedure, non_overridable           :: get_set_id              => cell_iterator_get_set_id
    procedure, non_overridable           :: get_num_vefs            => cell_iterator_get_num_vefs
    procedure, non_overridable           :: get_num_nodes           => cell_iterator_get_num_nodes
    procedure, non_overridable           :: get_node_gid            => cell_iterator_get_node_gid
    procedure, non_overridable           :: get_vef_gid             => cell_iterator_get_vef_gid
    procedure, non_overridable           :: get_vef_gids  => cell_iterator_get_vef_gids
    procedure, non_overridable           :: get_vef_ggid            => cell_iterator_get_vef_ggid
    procedure, non_overridable           :: get_vef_lid_from_gid    => cell_iterator_get_vef_lid_from_gid
    procedure, non_overridable           :: get_vef_lid_from_ggid   => cell_iterator_get_vef_lid_from_ggid
    procedure, non_overridable           :: get_vef                 => cell_iterator_get_vef
    procedure, non_overridable           :: is_local                => cell_iterator_is_local
    procedure, non_overridable           :: is_ghost                => cell_iterator_is_ghost
    procedure, non_overridable           :: scan_sum_num_vefs    => cell_iterator_get_scan_sum_num_vefs

    ! Declare dummy procedures to be implemented in the corresponding derived classes 
    procedure :: update_sub_triangulation    => cell_iterator_update_sub_triangulation
    procedure :: get_mc_case                 => cell_iterator_get_mc_case
    procedure :: get_num_subcells      => cell_iterator_get_num_subcells
    procedure :: get_num_subcell_nodes => cell_iterator_get_num_subcell_nodes
    procedure :: get_phys_coords_of_subcell  => cell_iterator_get_phys_coords_of_subcell
    procedure :: get_ref_coords_of_subcell   => cell_iterator_get_ref_coords_of_subcell
    procedure :: get_num_subfacets      => cell_iterator_get_num_subfacets
    procedure :: get_num_subfacet_nodes => cell_iterator_get_num_subfacet_nodes
    procedure :: get_phys_coords_of_subfacet  => cell_iterator_get_phys_coords_of_subfacet
    procedure :: get_ref_coords_of_subfacet   => cell_iterator_get_ref_coords_of_subfacet
    procedure :: is_cut                      => cell_iterator_is_cut
    procedure :: is_interior                 => cell_iterator_is_interior
    procedure :: is_exterior                 => cell_iterator_is_exterior
    procedure :: is_interior_subcell         => cell_iterator_is_interior_subcell
    procedure :: is_exterior_subcell         => cell_iterator_is_exterior_subcell

    procedure, non_overridable, private  :: fill_nodes_on_vertices        => cell_iterator_fill_nodes_on_vertices
    procedure, non_overridable, private  :: fill_nodes_on_vef_new         => cell_iterator_fill_nodes_on_vef_new
    procedure, non_overridable, private  :: fill_nodes_on_vef_from_source => cell_iterator_fill_nodes_on_vef_from_source
    procedure, non_overridable, private  :: fill_internal_nodes_new       => cell_iterator_fill_internal_nodes_new
  end type cell_iterator_t
  
  type vef_iterator_t
    private
    integer(ip)                                 :: gid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
     procedure                           :: create                    => vef_iterator_create
     procedure                           :: free                      => vef_iterator_free
     final                               ::                              vef_iterator_free_final
     procedure                           :: first                     => vef_iterator_first
     procedure                           :: next                      => vef_iterator_next
     procedure, non_overridable          :: set_gid                   => vef_iterator_set_gid
     procedure                           :: has_finished              => vef_iterator_has_finished
     procedure, non_overridable          :: get_coordinates           => vef_iterator_get_coordinates
     procedure, non_overridable          :: get_triangulation         => vef_iterator_get_triangulation
     procedure, non_overridable          :: get_gid                   => vef_iterator_get_gid
     procedure, non_overridable          :: get_ggid                  => vef_iterator_get_ggid

     procedure, non_overridable          :: set_geom_id               => vef_iterator_set_geom_id
     procedure, non_overridable          :: set_set_id                => vef_iterator_set_set_id
     procedure, non_overridable          :: get_geom_id               => vef_iterator_get_geom_id
     procedure, non_overridable          :: get_set_id                => vef_iterator_get_set_id

     procedure, non_overridable          :: set_dim                   => vef_iterator_set_dim
     procedure, non_overridable          :: set_it_at_boundary        => vef_iterator_set_it_at_boundary
     procedure, non_overridable          :: set_it_as_local           => vef_iterator_set_it_as_local
     procedure, non_overridable          :: set_it_as_ghost           => vef_iterator_set_it_as_ghost
     procedure, non_overridable          :: set_it_at_interface       => vef_iterator_set_it_at_interface

     procedure, non_overridable          :: get_dim             => vef_iterator_get_dim
     procedure, non_overridable          :: is_at_boundary            => vef_iterator_is_at_boundary
     procedure, non_overridable          :: is_local                  => vef_iterator_is_local
     procedure, non_overridable          :: is_ghost                  => vef_iterator_is_ghost
     procedure, non_overridable          :: is_at_interface           => vef_iterator_is_at_interface
     procedure, non_overridable          :: is_face                   => vef_iterator_is_face
     
     procedure, non_overridable          :: get_num_cells_around      => vef_iterator_get_num_cells_around
     procedure, non_overridable          :: vef_iterator_get_cell_around
     generic                             :: get_cell_around           => vef_iterator_get_cell_around
  end type vef_iterator_t

  type, extends(vef_iterator_t) :: itfc_vef_iterator_t
    private
    integer(ip)  :: itfc_gid = -1
    contains
     procedure          :: first          => itfc_vef_iterator_first
     procedure          :: next           => itfc_vef_iterator_next
     procedure          :: has_finished   => itfc_vef_iterator_has_finished
  end type itfc_vef_iterator_t


  type object_iterator_t
    private
    integer(ip)                                 :: gid = -1
    type(list_iterator_t)                       :: vefs_object_iterator
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
    procedure                           :: create                          => object_iterator_create
    procedure                           :: free                            => object_iterator_free
    final                               ::                                    object_iterator_free_final
    procedure, non_overridable, private :: update_vefs_object_iterator     => object_iterator_update_vefs_object_iterator
    procedure                           :: first                           => object_iterator_first
    procedure                           :: next                            => object_iterator_next
    procedure                           :: set_gid                         => object_iterator_set_gid
    procedure                           :: has_finished                    => object_iterator_has_finished
    procedure, non_overridable          :: get_gid                         => object_iterator_get_gid
    procedure, non_overridable          :: get_ggid                         => object_iterator_get_ggid
    procedure, non_overridable          :: get_dim                   => object_iterator_get_dim
    procedure, non_overridable          :: get_num_parts_around         => object_iterator_get_num_parts_around
    procedure, non_overridable          :: get_num_subparts_around      => object_iterator_get_num_subparts_around
    procedure, non_overridable          :: create_parts_around_iterator    => object_iterator_create_parts_around_iterator
    procedure, non_overridable          :: create_subparts_around_iterator => object_iterator_create_subparts_around_iterator
    procedure, non_overridable          :: get_num_vefs                    => object_iterator_get_num_vefs
    procedure, non_overridable          :: get_vef                         => object_iterator_get_vef
  end type object_iterator_t
       
   ! JP-TODO: implement states: discuss with Alberto and Victor.
   !
   ! State transition diagram for type(base_static_triangulation_t). The
   ! creation of reference elements must be performed together
   ! with reading from mesh when going from created to elements_filled.
   !
   ! ------------------------------------------------...--------
   ! Input State      | Action                    | Output State 
   ! -----------------------------------------------------------
   ! not_created      | create                    | created 
   ! not_created      | free                      | not_created
   ! created          | free                      | not_created
   ! created          | external (read from mesh) | elements_filled
   ! elements_filled  | create_dual               | vefs_filled
   ! elements_filled  | free                      | not_created
   ! vefs_filled      | free                      | not_created

   !integer(ip), parameter :: base_static_triangulation_not_created     = 0 
   !integer(ip), parameter :: base_static_triangulation_created         = 1 
   !integer(ip), parameter :: base_static_triangulation_elements_filled = 2 
   !integer(ip), parameter :: base_static_triangulation_vefs_filled     = 3 

  integer(ip), parameter      :: max_num_reference_fes_geo = 3

  type base_static_triangulation_t ! Base class for coarse_triangulation_t and fine_triangulation_t
     private

     ! Parallel environment describing MPI tasks among which the triangulation is distributed
     ! (NULL for serial_triangulation_t)
     type(environment_t), pointer      :: p_env => NULL()
     type(environment_t)               :: par_environment
    
     ! Sizes
     integer(ip)                           :: num_dims  = 0

     ! Data structures to store cell related information
     integer(ip)                           :: num_local_cells = 0
     integer(ip)                           :: num_ghost_cells = 0
     integer(ip)                           :: max_vefs_x_cell = 0
     integer(igp), allocatable             :: cells_gid(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_mypart(:)            ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_set(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: ptr_vefs_x_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable             :: lst_vefs_gids(:)

     ! Data type describing the layout in distributed-memory of the dual graph
     ! (It is required, e.g., for nearest neighbour comms on this graph)
     type(cell_import_t)                   :: cell_import   
     
     ! Data structures to store vef related information
     integer(ip)                           :: num_vefs = 0         ! = num_vertices + num_edges + num_faces
     integer(ip)                           :: num_vertices = 0
     integer(ip)                           :: num_edges = 0
     integer(ip)                           :: num_faces = 0
     integer(igp), allocatable             :: vefs_gid(:)          ! num_vefs
     integer(ip) , allocatable             :: vefs_set(:)          ! num_vefs
     integer(ip) , allocatable             :: vefs_geometry(:)     ! num_vefs
     integer(ip) , allocatable             :: vefs_type(:)         ! num_vefs, will replace vefs_dim
                                                                   ! above and vef_itfc_gid below (which is currently only accessed
                                                                   ! to check whether a vef is interface or not).
     integer(ip)                           :: num_itfc_vefs  = 0
     integer(ip), allocatable              :: lst_itfc_vefs(:)
     integer(ip), allocatable              :: ptrs_cells_around(:) ! num_itfc_vefs+1
     integer(ip), allocatable              :: lst_cells_around(:)  ! ptrs_cells_around(num_itfc_vefs+1)-1

     ! Data structures to create objects (coarse cell info)
     integer(ip)                             :: num_global_objects = 0
     integer(ip)                             :: num_objects = 0
     integer(igp), allocatable               :: objects_gids(:)
     integer(ip) , allocatable               :: objects_dim(:)
     type(list_t)                            :: vefs_object
     type(list_t)                            :: faces_object
     type(list_t)                            :: parts_object
     integer(ip)                             :: num_subparts        ! Number of subparts around part (including those subparts which are local)
     type(list_t)                            :: subparts_object        ! Number and list of subparts GIDs around each coarse n_face
     type(hash_table_ip_ip_t)                :: g2l_subparts           ! Translator among the GIDs of subparts and gids
     type(coarse_triangulation_t), pointer   :: coarse_triangulation => NULL()

     ! Data structures that should be defined in fine_triangulation_t (which requires extensive refactoring)     
     type(geometry_t)                        :: geometry
     type(p_lagrangian_reference_fe_t)       :: reference_fe_geo_list(max_num_reference_fes_geo)
     type(hash_table_ip_ip_t)                :: reference_fe_geo_index
     
     ! Geometry interpolation
     integer(ip)                             :: num_nodes
     integer(ip) , allocatable               :: ptr_nodes_x_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable               :: lst_nodes(:)
     type(point_t), allocatable              :: coordinates(:)
     
 contains  

     ! Getters
     procedure, non_overridable          :: get_par_environment                 => bst_get_par_environment
     procedure, non_overridable          :: get_cell_import                     => bst_get_cell_import
     procedure, non_overridable          :: get_coarse_triangulation            => bst_get_coarse_triangulation
     procedure, non_overridable          :: get_num_dims                  => bst_get_num_dims
     procedure, non_overridable          :: get_num_vefs                        => bst_get_num_vefs 
     procedure, non_overridable          :: get_num_facets                       => bst_get_num_facets
     procedure, non_overridable          :: get_num_cells                       => bst_get_num_cells
     procedure, non_overridable          :: get_num_local_cells                 => bst_get_num_local_cells
     procedure, non_overridable          :: get_num_ghost_cells                 => bst_get_num_ghost_cells
     procedure, non_overridable          :: get_num_objects                  => bst_get_num_objects

     ! Cell traversals-related TBPs
     procedure                           :: create_cell_iterator                => bst_create_cell_iterator
     procedure                           :: free_cell_iterator                  => bst_free_cell_iterator
     
     ! Vef traversals-related TBPs
     procedure                           :: create_vef_iterator              => bst_create_vef_iterator
     procedure                           :: create_itfc_vef_iterator         => bst_create_itfc_vef_iterator
     procedure                           :: free_vef_iterator                => bst_free_vef_iterator
     
     ! Objects-related traversals
     procedure, non_overridable          :: create_object_iterator              => bst_create_object_iterator
     procedure, non_overridable          :: free_object_iterator                => bst_free_object_iterator
       
     ! Other
     procedure                           :: print                               => bst_print    
     procedure                           :: free                                => bst_free
     
     ! Getters
     procedure                           :: get_num_reference_fes_geo        => bst_get_num_reference_fes_geo
     procedure                           :: get_max_num_shape_functions      => bst_get_max_num_shape_functions
     procedure                           :: is_tet_mesh                         => bst_is_tet_mesh
     procedure                           :: is_hex_mesh                         => bst_is_hex_mesh
     procedure                           :: is_mix_mesh                         => bst_is_mix_mesh
   

     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_ptr_vefs_x_cell => bst_allocate_and_fill_ptr_vefs_x_cell
     procedure, non_overridable, private :: allocate_cells_gid                  => bst_allocate_cells_gid
     procedure, non_overridable, private :: fill_local_cells_gid                => bst_fill_local_cells_gid
     procedure, non_overridable, private :: allocate_cells_mypart               => bst_allocate_cells_mypart
     procedure, non_overridable, private :: fill_local_cells_mypart             => bst_fill_local_cells_mypart
     procedure, non_overridable, private :: allocate_cells_set                  => bst_allocate_cells_set
     procedure, non_overridable          :: fill_cells_set                      => bst_fill_cells_set
     procedure, non_overridable, private :: free_ptr_vefs_x_cell              => bst_free_ptr_vefs_x_cell
     procedure, non_overridable, private :: free_lst_vefs_gids                  => bst_free_lst_vefs_gids 
     procedure, non_overridable, private :: free_cells_gid                      => bst_free_cells_gid
     procedure, non_overridable, private :: free_cells_mypart                   => bst_free_cells_mypart
     procedure, non_overridable, private :: free_cells_set                      => bst_free_cells_set
     procedure, non_overridable, private :: orient_tet_mesh                     => bst_orient_tet_mesh

     ! Private methods to perform nearest neighbor exchange
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_gids
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_gids_and_dim
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_gids_and_coordinates
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_gids
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_gids_and_dim
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_gids_and_coordinates
     procedure, non_overridable, private :: fetch_ghost_cells_data         => bst_fetch_ghost_cells_data
     procedure, non_overridable, nopass, private :: cell_size              => bst_cell_size
     generic,                            private :: cell_pack              => bst_cell_pack_vef_gids, &
                                                                              bst_cell_pack_vef_gids_and_dim, &
                                                                              bst_cell_pack_vef_gids_and_coordinates
     generic,                            private :: cell_unpack            => bst_cell_unpack_vef_gids, &
                                                                              bst_cell_unpack_vef_gids_and_dim, &
                                                                              bst_cell_unpack_vef_gids_and_coordinates

     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: compute_num_vefs                    => bst_compute_num_vefs
     procedure, non_overridable, private :: allocate_and_fill_vefs_gid          => bst_allocate_and_fill_vefs_gid
     procedure, non_overridable, private :: free_vefs_gid                       => bst_free_vefs_gid
     procedure, non_overridable, private :: free_vefs_type                      => bst_triangulation_free_vefs_type
     procedure, non_overridable, private :: allocate_and_fill_cells_around      => bst_allocate_and_fill_cells_around
     procedure, non_overridable, private :: free_cells_around                   => bst_free_cells_around
     procedure, non_overridable, private :: find_and_list_vefs_at_interfaces    => bst_find_and_list_vefs_at_interfaces
     procedure, non_overridable, private :: free_lst_itfc_vefs                  => bst_free_lst_itfc_vefs

     ! Private methods to compute objects
     procedure, non_overridable          :: get_num_subparts                            => bst_get_num_subparts
     procedure, non_overridable          :: get_subpart_lid                                => bst_get_subpart_lid
     procedure, non_overridable, private :: compute_vefs_and_parts_object                  => bst_compute_vefs_and_parts_object
     procedure, non_overridable, private :: compute_vefs_and_parts_object_body             => bst_compute_vefs_and_parts_object_body
     procedure, non_overridable, private :: compute_parts_itfc_vefs                        => bst_compute_parts_itfc_vefs
     procedure, non_overridable, private :: compute_subparts_itfc_vefs                     => bst_compute_subparts_itfc_vefs
     procedure, non_overridable, private :: compute_parts_object_from_subparts_object      => bst_compute_parts_object_from_subparts_object
     procedure, non_overridable, private :: compute_part_id_from_subpart_gid               => bst_compute_part_id_from_subpart_gid
     procedure, non_overridable, private :: compute_objects_dim                      => bst_compute_objects_dim
     procedure, non_overridable, private :: compute_objects_neighbours_exchange_data       => bst_compute_objects_neighbours_exchange_data
     procedure, non_overridable, private :: compute_num_global_objects_and_their_gids   => bst_compute_num_global_objs_and_their_gids
     procedure, non_overridable, private :: free_objects_gids_and_dim                      => bst_free_objects_gids_and_dim

     ! Private methods for coarser triangulation set-up
     procedure, non_overridable          :: setup_coarse_triangulation                     => bst_setup_coarse_triangulation
     procedure, non_overridable, private :: gather_coarse_cell_gids                        => bst_gather_coarse_cell_gids
     procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs       => bst_gather_coarse_vefs_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_vefs_gids                        => bst_gather_coarse_vefs_gids
     procedure, non_overridable, private :: gather_coarse_vefs_dim                   => bst_gather_coarse_vefs_dim
     procedure, non_overridable, private :: fetch_l2_part_id_neighbours                    => bst_fetch_l2_part_id_neighbours
     procedure, non_overridable, private :: gather_coarse_dgraph_rcv_counts_and_displs     => bst_gather_coarse_dgraph_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_dgraph_lextn_and_lextp           => bst_gather_coarse_dgraph_lextn_and_lextp
     procedure, non_overridable, private :: adapt_coarse_raw_arrays                        => bst_adapt_coarse_raw_arrays
     
  end type base_static_triangulation_t
  
  type, extends(base_static_triangulation_t) :: fine_triangulation_t
     private
  contains
     ! Private methods to create vefs (these functions make use of the reference fe and therefore are not bounded
     ! to the mother class)
     procedure, non_overridable, private :: fill_reference_fe_geo_list          => fine_triangulation_fill_reference_fe_geo_list
     procedure, non_overridable, private :: generate_vefs                       => fine_triangulation_generate_vefs
     procedure, non_overridable, private :: allocate_and_fill_geometry_and_set  => fine_triangulation_allocate_and_fill_geometry_and_set
     procedure, non_overridable, private :: free_geometry_and_set               => fine_triangulation_free_geometry_and_set
     procedure, non_overridable, private :: compute_vefs_dim              => fine_triangulation_compute_vefs_dim
     procedure, non_overridable, private :: find_vefs_at_boundary               => fine_triangulation_find_vefs_at_boundary
     ! Geometry interpolation 
     procedure, non_overridable, private :: allocate_and_fill_nodes             => fine_triangulation_allocate_and_fill_nodes
     procedure, non_overridable, private :: free_nodes                          => fine_triangulation_free_nodes
     procedure, non_overridable          :: allocate_and_fill_coordinates       => fine_triangulation_allocate_and_fill_coordinates
     procedure, non_overridable          :: free_coordinates                    => fine_triangulation_free_coordinates
     procedure                           :: free                                => fine_triangulation_free     
  end type fine_triangulation_t


  type, extends(fine_triangulation_t) :: serial_triangulation_t
  contains
     procedure                 , private :: serial_triangulation_create
     generic                             :: create                              => serial_triangulation_create
  end type serial_triangulation_t
  
  type, extends(fine_triangulation_t) :: par_triangulation_t
  contains
     generic                             :: create                              => par_triangulation_create
     procedure, private                  :: par_triangulation_create
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_gids     => par_triangulation_allocate_and_fill_lst_vefs_gids 
  end type par_triangulation_t

  
  type, extends(base_static_triangulation_t) :: coarse_triangulation_t
  contains
     procedure, non_overridable          :: create                              => coarse_triangulation_create
     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_gids     => coarse_triangulation_allocate_and_fill_lst_vefs_gids 
     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: allocate_and_fill_vefs_dim    => coarse_triangulation_allocate_and_fill_vefs_dim
     procedure                           :: print                               => coarse_triangulation_print
  end type coarse_triangulation_t

  public :: base_static_triangulation_t
  public :: serial_triangulation_t
  public :: coarse_triangulation_t 
  public :: par_triangulation_t
  public :: cell_iterator_t
  public :: vef_iterator_t
  public :: itfc_vef_iterator_t, object_iterator_t
  
contains

#include "sbm_cell_iterator.i90"
#include "sbm_vef_iterator.i90"
#include "sbm_object_iterator.i90"
#include "sbm_base_static_triangulation.i90"
#include "sbm_fine_triangulation.i90"
#include "sbm_serial_triangulation.i90"
#include "sbm_par_triangulation.i90"
#include "sbm_coarse_triangulation.i90"

end module base_static_triangulation_names
