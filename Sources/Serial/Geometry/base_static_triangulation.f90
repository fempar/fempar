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
  use reference_fe_factory_names
  use field_names
  use cell_import_names
  use hash_table_names
  use list_types_names
  use mesh_names
  use mesh_distribution_names
  use par_io_names
  use stdio_names
  use field_names

  ! Geometry modules
  use sisl_names
  use geometry_names

  ! Par modules
  use par_environment_names

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

  type cell_accessor_t
    private
    integer(ip)                                 :: lid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation
  contains
    ! create/free/next/set_lid CANNOT longer be private as type(coarse_fe_accessor_t)
    ! extends type(cell_accessor_t), requires to call these TBPs and cannot be placed
    ! either within this module or within a submodule of base_static_triangulation_names
    ! due to module dependencies cycle
    procedure                            :: cell_accessor_create
    generic                              :: create                  => cell_accessor_create
    procedure                            :: free                    => cell_accessor_free
    procedure, non_overridable           :: next                    => cell_accessor_next
    procedure, non_overridable           :: set_lid                 => cell_accessor_set_lid
    procedure, non_overridable, private  :: set_gid                 => cell_accessor_set_gid
    procedure, non_overridable, private  :: set_mypart              => cell_accessor_set_mypart
    procedure, non_overridable, private  :: get_triangulation       => cell_accessor_get_triangulation
    procedure, non_overridable, private  :: cell_accessor_get_vef
    procedure, non_overridable           :: past_the_end            => cell_accessor_past_the_end
    procedure, non_overridable           :: get_reference_fe_geo    => cell_accessor_get_reference_fe_geo
    procedure, non_overridable           :: get_reference_fe_geo_id => cell_accessor_get_reference_fe_geo_id
    procedure, non_overridable           :: get_coordinates         => cell_accessor_get_coordinates
    procedure, non_overridable           :: set_coordinates         => cell_accessor_set_coordinates
    procedure, non_overridable           :: get_lid                 => cell_accessor_get_lid
    procedure, non_overridable           :: get_gid                 => cell_accessor_get_gid
    procedure, non_overridable           :: get_mypart              => cell_accessor_get_mypart
    procedure, non_overridable           :: get_num_vefs            => cell_accessor_get_num_vefs
    procedure, non_overridable           :: get_num_nodes           => cell_accessor_get_num_nodes
    procedure, non_overridable           :: get_node_lid            => cell_accessor_get_node_lid
    procedure, non_overridable           :: get_vef_lid             => cell_accessor_get_vef_lid
    procedure, non_overridable           :: get_vef_lids            => cell_accessor_get_vef_lids
    procedure, non_overridable           :: get_vef_gid             => cell_accessor_get_vef_gid
    procedure, non_overridable           :: find_lpos_vef_lid       => cell_accessor_find_lpos_vef_lid
    procedure, non_overridable           :: find_lpos_vef_gid       => cell_accessor_find_lpos_vef_gid
    generic                              :: get_vef                 => cell_accessor_get_vef
    procedure, non_overridable           :: is_local                => cell_accessor_is_local
    procedure, non_overridable           :: is_ghost                => cell_accessor_is_ghost
    procedure, non_overridable           :: scan_sum_number_vefs    => cell_accessor_get_scan_sum_number_vefs
    procedure, non_overridable           :: fill_nodes_on_vertices  => cell_accessor_fill_nodes_on_vertices
  end type cell_accessor_t
  
  type cell_iterator_t
    private
    type(cell_accessor_t) :: current_cell_accessor
  contains
     procedure, non_overridable, private ::                 cell_iterator_current
     procedure, non_overridable, private :: create       => cell_iterator_create
     procedure, non_overridable          :: free         => cell_iterator_free
     procedure, non_overridable          :: init         => cell_iterator_init
     procedure, non_overridable          :: next         => cell_iterator_next
     procedure, non_overridable          :: has_finished => cell_iterator_has_finished
     generic                             :: current      => cell_iterator_current
  end type cell_iterator_t  
  
  type vef_accessor_t
    private
    integer(ip)                                 :: lid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
     ! create/free/next/set_lid/past_the_end/get_triangulation CANNOT longer be private as 
     ! type(coarse_fe_vef_accessor_t) extends type(vef_accessor_t), requires to call these TBPs
     procedure                           :: vef_accessor_create
     procedure, non_overridable          :: vef_accessor_get_cell_around
     generic                             :: create => vef_accessor_create
     procedure                           :: free                      => vef_accessor_free
     procedure, non_overridable          :: next                      => vef_accessor_next
     procedure, non_overridable          :: set_lid                   => vef_accessor_set_lid
     procedure, non_overridable          :: past_the_end              => vef_accessor_past_the_end
     procedure, non_overridable          :: get_triangulation         => vef_accessor_get_triangulation
     procedure, non_overridable          :: get_lid                   => vef_accessor_get_lid
     procedure, non_overridable          :: get_gid                   => vef_accessor_get_gid

     procedure, non_overridable          :: set_geom_id               => vef_accessor_set_geom_id
     procedure, non_overridable          :: set_set_id                => vef_accessor_set_set_id
     procedure, non_overridable          :: get_geom_id               => vef_accessor_get_geom_id
     procedure, non_overridable          :: get_set_id                => vef_accessor_get_set_id

     procedure, non_overridable          :: set_dimension             => vef_accessor_set_dimension
     procedure, non_overridable          :: set_it_at_boundary        => vef_accessor_set_it_at_boundary
     procedure, non_overridable          :: set_it_as_local           => vef_accessor_set_it_as_local
     procedure, non_overridable          :: set_it_as_ghost           => vef_accessor_set_it_as_ghost
     procedure, non_overridable          :: set_it_at_interface       => vef_accessor_set_it_at_interface

     procedure, non_overridable          :: get_dimension             => vef_accessor_get_dimension
     procedure, non_overridable          :: is_at_boundary            => vef_accessor_is_at_boundary
     procedure, non_overridable          :: is_local                  => vef_accessor_is_local
     procedure, non_overridable          :: is_ghost                  => vef_accessor_is_ghost
     procedure, non_overridable          :: is_at_interface           => vef_accessor_is_at_interface
     procedure, non_overridable          :: is_face                   => vef_accessor_is_face
     
     procedure, non_overridable          :: get_num_cells_around      => vef_accessor_get_num_cells_around
     generic                             :: get_cell_around           => vef_accessor_get_cell_around
     !procedure, non_overridable          :: get_vertices              => vef_accessor_get_vertices
  end type vef_accessor_t

  ! So far this is just a rename
  !type, extends(vef_accessor_t) :: vertex_accessor_t
  !end type vertex_accessor_t

  ! So far this is just a rename
  !type, extends(vef_accessor_t) :: edge_accessor_t
  !end type edge_accessor_t

  ! This will include functions to ask for orientation, rotation, etc.
  type, extends(vef_accessor_t) :: face_accessor_t
    private
  contains
    procedure, non_overridable          :: get_coordinates                  => face_accessor_get_coordinates
    procedure, non_overridable          :: get_face_lid                     => face_accessor_get_face_lid
    procedure, non_overridable          :: get_face_lpos_within_cell_around => face_accessor_get_face_lpos_within_cell_around
    procedure, non_overridable          :: get_face_orientation             => face_accessor_get_face_orientation
    procedure, non_overridable          :: get_face_rotation                => face_accessor_get_face_rotation
    procedure, non_overridable, private :: set_face_orientation             => face_accessor_set_face_orientation
    procedure, non_overridable, private :: set_face_rotation                => face_accessor_set_face_rotation
  end type face_accessor_t

  ! In order to define iterators over vertices, edges and faces as extensions of vef_iterator
  ! we need to overwrite init and next. The alternative is to repeat all TBPs.
  type vef_iterator_t
    private
    type(vef_accessor_t) :: current_vef_accessor
  contains
     procedure, non_overridable, private ::                 vef_iterator_current
     procedure, non_overridable          :: create       => vef_iterator_create
     procedure, non_overridable          :: free         => vef_iterator_free
     procedure                           :: init         => vef_iterator_init
     procedure                           :: next         => vef_iterator_next
     procedure, non_overridable          :: has_finished => vef_iterator_has_finished
     generic                             :: current      => vef_iterator_current
  end type vef_iterator_t
  
  type :: itfc_vef_iterator_t
    private
    integer(ip)          :: itfc_lid = -1
    type(vef_accessor_t) :: current_vef_accessor
  contains
    procedure, non_overridable, private ::                 itfc_vef_iterator_current
    procedure, non_overridable          :: create       => itfc_vef_iterator_create
    procedure, non_overridable          :: free         => itfc_vef_iterator_free
    procedure, non_overridable          :: init         => itfc_vef_iterator_init
    procedure, non_overridable          :: next         => itfc_vef_iterator_next
    procedure, non_overridable          :: has_finished => itfc_vef_iterator_has_finished
    generic                             :: current      => itfc_vef_iterator_current
  end type itfc_vef_iterator_t

  type, extends(vef_iterator_t) :: vertex_iterator_t
    private
  contains
    procedure, non_overridable          :: next         => vertex_iterator_next
  end type vertex_iterator_t

  type, extends(vef_iterator_t) :: edge_iterator_t
    private
  contains
    procedure, non_overridable          :: init         => edge_iterator_init
    procedure, non_overridable          :: next         => edge_iterator_next
  end type edge_iterator_t

  type, extends(vef_iterator_t) :: face_iterator_t
    private
  contains
    procedure, non_overridable          :: init         => face_iterator_init
    procedure, non_overridable, private :: face_iterator_current
    generic                             :: current      => face_iterator_current
  end type face_iterator_t


  type object_accessor_t
    private
    integer(ip)                                 :: lid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation
  contains
    procedure                   :: object_accessor_create
    generic                     :: create                         => object_accessor_create
    procedure                   :: free                           => object_accessor_free
    procedure                   :: next                           => object_accessor_next
    procedure                   :: set_lid                        => object_accessor_set_lid
    procedure                   :: past_the_end                   => object_accessor_past_the_end
    procedure, non_overridable  :: get_lid                        => object_accessor_get_lid
    procedure, non_overridable  :: get_gid                        => object_accessor_get_gid
    procedure, non_overridable  :: get_number_parts_around        => object_accessor_get_number_parts_around
    procedure, non_overridable  :: create_parts_around_iterator   => object_accessor_create_parts_around_iterator
    procedure, non_overridable  :: get_number_vefs_on_object      => object_accessor_get_number_vefs_on_object
    procedure, non_overridable  :: create_vefs_on_object_iterator => object_accessor_get_vefs_on_object_iterator
  end type object_accessor_t
  
  type object_iterator_t
    private
    type(object_accessor_t) :: current_object_accessor
  contains
     procedure, non_overridable, private ::                 object_iterator_current
     procedure, non_overridable          :: create       => object_iterator_create
     procedure, non_overridable          :: free         => object_iterator_free
     procedure, non_overridable          :: init         => object_iterator_init
     procedure, non_overridable          :: next         => object_iterator_next
     procedure, non_overridable          :: has_finished => object_iterator_has_finished
     generic                             :: current      => object_iterator_current
  end type object_iterator_t
  
  type :: vefs_on_object_iterator_t
    private
    type(list_iterator_t) :: vefs_lids_on_object_iterator
    type(vef_accessor_t)  :: current_vef_accessor
  contains
    procedure, non_overridable, private ::                 vefs_on_object_iterator_current
    procedure, non_overridable          :: create       => vefs_on_object_iterator_create
    procedure, non_overridable          :: free         => vefs_on_object_iterator_free
    procedure, non_overridable          :: init         => vefs_on_object_iterator_init
    procedure, non_overridable          :: next         => vefs_on_object_iterator_next
    procedure, non_overridable          :: has_finished => vefs_on_object_iterator_has_finished
    generic                             :: current      => vefs_on_object_iterator_current
  end type vefs_on_object_iterator_t
   
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
     type(par_environment_t),   pointer      :: p_env => NULL()
     
     ! Sizes
     integer(ip)                           :: num_dimensions  = -1

     ! Data structures to store cell related information
     integer(ip)                           :: num_local_cells = -1
     integer(ip)                           :: num_ghost_cells = -1
     integer(ip)                           :: max_vefs_per_cell = -1
     integer(igp), allocatable             :: cells_gid(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_mypart(:)            ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: ptr_vefs_per_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable             :: lst_vefs_lids(:)

     ! Data type describing the layout in distributed-memory of the dual graph
     ! (It is required, e.g., for nearest neighbour comms on this graph)
     type(cell_import_t)                   :: cell_import   
     
     ! Data structures to store vef related information
     integer(ip)                           :: num_vefs = -1        ! = num_vertices + num_edges + num_faces
     integer(ip)                           :: num_vertices = 0
     integer(ip)                           :: num_edges = 0
     integer(ip)                           :: num_faces = 0
     integer(igp), allocatable             :: vefs_gid(:)          ! num_vefs
     integer(ip) , allocatable             :: vefs_set(:)          ! num_vefs
     integer(ip) , allocatable             :: vefs_geometry(:)     ! num_vefs
     integer(ip) , allocatable             :: vefs_dimension(:)    ! num_vefs
     integer(ip) , allocatable             :: vefs_type(:)         ! num_vefs, will replace vefs_dimension
                                                                   ! above and vef_itfc_lid below (which is currently only accessed
                                                                   ! to check whether a vef is interface or not).
     integer(ip)                           :: num_itfc_vefs  = -1
     integer(ip), allocatable              :: lst_itfc_vefs(:)
     integer(ip), allocatable              :: ptrs_cells_around(:) ! num_itfc_vefs+1
     integer(ip), allocatable              :: lst_cells_around(:)  ! ptrs_cells_around(num_itfc_vefs+1)-1
     
     ! Data structures related to face integration
     integer(ip), allocatable              :: face_orientation(:)
     integer(ip), allocatable              :: face_rotation(:)

     ! Data structures to create objects
     integer(ip)                             :: number_global_objects = -1
     integer(ip)                             :: number_objects = -1
     integer(igp), allocatable               :: objects_gids(:)
     integer(ip) , allocatable               :: objects_dimension(:)
     type(list_t)                            :: vefs_object
     type(list_t)                            :: parts_object
     type(coarse_triangulation_t), pointer   :: coarse_triangulation

     ! Data structures that should be defined in fine_triangulation_t (which requires extensive refactoring)     
     type(geometry_t)                        :: geometry
     type(p_lagrangian_reference_fe_t)       :: reference_fe_geo_list(max_num_reference_fes_geo)
     type(hash_table_ip_ip_t)                :: reference_fe_geo_index
     
     ! Geometry interpolation
     integer(ip)                             :: num_nodes
     integer(ip) , allocatable               :: ptr_nodes_per_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable               :: lst_nodes(:)
     type(point_t), allocatable              :: coordinates(:)
     
 contains  

     ! Getters
     procedure, non_overridable          :: get_par_environment                 => bst_get_par_environment
     procedure, non_overridable          :: get_cell_import                     => bst_get_cell_import
     procedure, non_overridable          :: get_coarse_triangulation            => bst_get_coarse_triangulation
     procedure, non_overridable          :: get_num_dimensions                  => bst_get_num_dimensions
     procedure, non_overridable          :: get_num_vefs                        => bst_get_num_vefs 
     procedure, non_overridable          :: get_num_faces                       => bst_get_num_faces
     procedure, non_overridable          :: get_num_cells                       => bst_get_num_cells
     procedure, non_overridable          :: get_num_local_cells                 => bst_get_num_local_cells
     procedure, non_overridable          :: get_num_ghost_cells                 => bst_get_num_ghost_cells
     procedure, non_overridable          :: get_number_objects                  => bst_get_number_objects

     ! Cell traversals-related TBPs
     procedure, non_overridable          :: create_cell_iterator                => bst_create_cell_iterator
  
     ! Vef traversals-related TBPs
     procedure, non_overridable          :: create_vef_iterator                 => bst_create_vef_iterator
     procedure, non_overridable          :: create_vertex_iterator              => bst_create_vertex_iterator
     procedure, non_overridable          :: create_edge_iterator                => bst_create_edge_iterator
     procedure, non_overridable          :: create_face_iterator                => bst_create_face_iterator
     procedure, non_overridable          :: create_itfc_vef_iterator            => bst_create_itfc_vef_iterator

     ! Objects-related traversals
     procedure, non_overridable          :: create_object_iterator              => bst_create_object_iterator
     procedure, non_overridable          :: create_vefs_on_object_iterator      => bst_create_vefs_on_object_iterator
  
     ! Other
     procedure                           :: print                               => bst_print    
     procedure                           :: free                                => bst_free
     
     ! Getters
     procedure                           :: get_number_reference_fes_geo        => bst_get_number_reference_fes_geo
   

     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_ptr_vefs_per_cell => bst_allocate_and_fill_ptr_vefs_per_cell
     procedure, non_overridable, private :: allocate_cells_gid                  => bst_allocate_cells_gid
     procedure, non_overridable, private :: fill_local_cells_gid                => bst_fill_local_cells_gid
     procedure, non_overridable, private :: allocate_cells_mypart               => bst_allocate_cells_mypart
     procedure, non_overridable, private :: fill_local_cells_mypart             => bst_fill_local_cells_mypart
     procedure, non_overridable, private :: free_ptr_vefs_per_cell              => bst_free_ptr_vefs_per_cell
     procedure, non_overridable, private :: free_lst_vefs_lids                  => bst_free_lst_vefs_lids 
     procedure, non_overridable, private :: free_cells_gid                      => bst_free_cells_gid
     procedure, non_overridable, private :: free_cells_mypart                   => bst_free_cells_mypart

     ! Private methods to perform nearest neighbor exchange
     procedure, non_overridable, private :: fetch_ghost_cells_data              => bst_fetch_ghost_cells_data
     procedure, non_overridable, nopass, private :: cell_size                   => bst_cell_size
     procedure, non_overridable, nopass, private :: cell_pack                   => bst_cell_pack
     procedure, non_overridable, nopass, private :: cell_unpack                 => bst_cell_unpack

     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: compute_num_vefs                    => bst_compute_num_vefs
     procedure, non_overridable, private :: allocate_and_fill_vefs_gid          => bst_allocate_and_fill_vefs_gid
     procedure, non_overridable, private :: free_vefs_gid                       => bst_free_vefs_gid
     procedure, non_overridable, private :: free_vefs_type                      => bst_triangulation_free_vefs_type
     procedure, non_overridable, private :: allocate_and_fill_cells_around      => bst_allocate_and_fill_cells_around
     procedure, non_overridable, private :: free_cells_around                   => bst_free_cells_around
     procedure, non_overridable, private :: find_and_list_vefs_at_interfaces    => bst_find_and_list_vefs_at_interfaces
     procedure, non_overridable, private :: free_lst_itfc_vefs                  => bst_free_lst_itfc_vefs
     
     ! Private methods for creating face-related data (for face integration)
     procedure, non_overridable, private :: allocate_and_fill_face_orientation_rotation => bst_allocate_and_fill_face_orientation_rotation
     procedure, non_overridable, private :: free_face_orientation_rotation              => bst_free_face_orientation_rotation

     ! Private methods to compute objects
     procedure, non_overridable, private :: compute_vefs_and_parts_object                  => bst_compute_vefs_and_parts_object
     procedure, non_overridable, private :: compute_objects_dimension                      => bst_compute_objects_dimension
     procedure, non_overridable, private :: compute_objects_neighbours_exchange_data       => bst_compute_objects_neighbours_exchange_data
     procedure, non_overridable, private :: compute_number_global_objects_and_their_gids   => bst_compute_num_global_objs_and_their_gids
     procedure, non_overridable, private :: free_objects_gids_and_dim                      => bst_free_objects_gids_and_dim

     ! Private methods for coarser triangulation set-up
     procedure, non_overridable, private :: setup_coarse_triangulation                     => bst_setup_coarse_triangulation
     procedure, non_overridable, private :: gather_coarse_cell_gids                        => bst_gather_coarse_cell_gids
     procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs       => bst_gather_coarse_vefs_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_vefs_gids                        => bst_gather_coarse_vefs_gids
     procedure, non_overridable, private :: gather_coarse_vefs_dimension                   => bst_gather_coarse_vefs_dimension
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
     procedure, non_overridable, private :: generate_vefs                       => fine_triangulation_generate_vefs
     procedure, non_overridable, private :: allocate_and_fill_geometry_and_set  => fine_triangulation_allocate_and_fill_geometry_and_set
     procedure, non_overridable, private :: free_geometry_and_set               => fine_triangulation_free_geometry_and_set
     procedure, non_overridable, private :: compute_vefs_dimension              => fine_triangulation_compute_vefs_dimension
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
     procedure, non_overridable          :: create                              => serial_triangulation_create
  end type serial_triangulation_t
  
  type, extends(fine_triangulation_t) :: par_triangulation_t
  contains
     procedure, non_overridable          :: create                              => par_triangulation_create
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_lids     => par_triangulation_allocate_and_fill_lst_vefs_lids 
  end type par_triangulation_t

  
  type, extends(base_static_triangulation_t) :: coarse_triangulation_t
  contains
     procedure, non_overridable          :: create                              => coarse_triangulation_create
     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_lids     => coarse_triangulation_allocate_and_fill_lst_vefs_lids 
     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: allocate_and_fill_vefs_dimension    => coarse_triangulation_allocate_and_fill_vefs_dimension
     procedure                           :: print                               => coarse_triangulation_print
     !procedure, non_overridable, private :: compute_num_itfc_vefs               => coarse_triangulation_compute_num_itfc_vefs
     !procedure, non_overridable, private :: allocate_and_fill_lst_itfc_vefs     => coarse_triangulation_allocate_and_fill_lst_itfc_vefs
  end type coarse_triangulation_t

  public :: base_static_triangulation_t
  public :: serial_triangulation_t
  public :: coarse_triangulation_t 
  public :: par_triangulation_t
  public :: cell_iterator_t, vef_iterator_t, face_iterator_t, itfc_vef_iterator_t, object_iterator_t, vefs_on_object_iterator_t
  public :: cell_accessor_t, vef_accessor_t, face_accessor_t, object_accessor_t
  
contains

#include "sbm_cell_accessor.i90"
#include "sbm_cell_iterator.i90"
#include "sbm_vef_accessor.i90"
#include "sbm_vef_iterator.i90"
#include "sbm_object_accessor.i90"
#include "sbm_object_iterator.i90"
#include "sbm_vefs_on_object_iterator.i90"
#include "sbm_base_static_triangulation.i90"
#include "sbm_fine_triangulation.i90"
#include "sbm_serial_triangulation.i90"
#include "sbm_par_triangulation.i90"
#include "sbm_coarse_triangulation.i90"

end module base_static_triangulation_names
