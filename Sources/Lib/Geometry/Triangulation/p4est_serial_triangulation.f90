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
module p4est_serial_triangulation_names
  use, intrinsic :: iso_c_binding
  
  use types_names
  use list_types_names
  use stdio_names
  use memor_names
  use environment_names
  use uniform_hex_mesh_generator_names
  use p4est_bindings_names
  use reference_fe_names
  use base_static_triangulation_names
  use std_vector_integer_ip_names
  use field_names
  
  use FPL

  implicit none
# include "debug.i90"
  private
  
  ! For 2D
  integer(ip), parameter :: NUM_SUBCELLS_IN_TOUCH_FACE_2D = 2
  integer(ip), parameter :: NUM_CORNERS_2D                = 4
  integer(ip), parameter :: NUM_FACES_2D                  = 4
  integer(ip), parameter :: NUM_SUBFACES_FACE_2D          = 2
  integer(ip), parameter :: NUM_FACE_CORNERS_2D           = 2
  integer(ip), parameter :: NUM_FACES_AT_CORNER_2D        = 2
  integer(ip), parameter :: NUM_VEFS_2D              = NUM_CORNERS_2D+NUM_FACES_2D
  integer(ip), target :: P4EST_FACE_CORNERS_2D(NUM_FACE_CORNERS_2D,NUM_FACES_2D) = & 
                                                  reshape([1, 3,&
                                                           2, 4,&  
                                                           1, 2,&
                                                           3, 4], [NUM_FACE_CORNERS_2D,NUM_FACES_2D])
                                                  
  integer(ip), target :: P4EST_FACES_AT_CORNER_2D(NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D) = & 
                                                    reshape([1, 3,&
                                                             2, 3,&  
                                                             1, 4,&
                                                             2, 4], [NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D])
                                                    
  integer(ip), target :: FEMPAR_SUBCELLS_IN_TOUCH_FACE_2D(NUM_SUBCELLS_IN_TOUCH_FACE_2D,NUM_FACES_2D) = &
                                                    reshape([1, 2,&
                                                             3, 4,&  
                                                             1, 3,&
                                                             2, 4], [NUM_SUBCELLS_IN_TOUCH_FACE_2D, NUM_FACES_2D])

  integer(ip), target :: P4EST_CORNER_IN_FACE_2D(NUM_FACES_2D,NUM_CORNERS_2D) = & 
                                                  reshape([ 1,-1, 1,-1,&
                                                           -1, 1, 2,-1,&
                                                            2,-1,-1, 1,&
                                                           -1, 2,-1, 2],[NUM_FACES_2D,NUM_CORNERS_2D])
  
  integer(ip), target :: P4EST_OPPOSITE_CORNER_2D(NUM_CORNERS_2D) = [ 4, 3, 2, 1 ]
  integer(ip), target :: P4EST_2_FEMPAR_CORNER_2D(NUM_CORNERS_2D) = [ 1, 2, 3, 4 ]
  integer(ip), target :: P4EST_2_FEMPAR_FACE_2D  (NUM_FACES_2D)   = [ 3, 4, 1, 2 ]

  ! For 3D
  integer(ip), parameter :: NUM_SUBCELLS_IN_TOUCH_FACE_3D = 4
  integer(ip), parameter :: NUM_CORNERS_3D                = 8
  integer(ip), parameter :: NUM_FACES_3D                  = 6
  integer(ip), parameter :: NUM_EDGES_3D                  = 12
  integer(ip), parameter :: NUM_SUBFACES_FACE_3D          = 4
  integer(ip), parameter :: NUM_FACE_CORNERS_3D           = 4
  integer(ip), parameter :: NUM_EDGE_CORNERS_3D           = 2
  integer(ip), parameter :: NUM_FACES_AT_CORNER_3D        = 3
  integer(ip), parameter :: NUM_EDGES_AT_CORNER_3D        = 3
  integer(ip), parameter :: NUM_VEFS_3D              = NUM_CORNERS_3D+NUM_FACES_3D+NUM_EDGES_3D

  integer(ip), target :: P4EST_FACE_CORNERS_3D(NUM_FACE_CORNERS_3D,NUM_FACES_3D) = & 
                                                  reshape([1, 3, 5, 7,&
                                                           2, 4, 6, 8,&  
                                                           1, 2, 5, 6,&
                                                           3, 4, 7, 8,&
                                                           1, 2, 3, 4,&
                                                           5, 6, 7, 8], [NUM_FACE_CORNERS_3D,NUM_FACES_3D])

  integer(ip), target :: P4EST_EDGE_CORNERS_3D(NUM_EDGE_CORNERS_3D,NUM_EDGES_3D) = & 
                                                  reshape([ 1, 2,&
                                                            3, 4,&
                                                            5, 6,&
                                                            7, 8,&
                                                            1, 3,&
                                                            2, 4,&
                                                            5, 7,&
                                                            6, 8,&
                                                            1, 5,&
                                                            2, 6,&
                                                            3, 7,&
                                                            4, 8], [NUM_EDGE_CORNERS_3D,NUM_EDGES_3D])

  integer(ip), target :: P4EST_FACES_AT_CORNER_3D(NUM_FACES_AT_CORNER_3D, NUM_CORNERS_3D) = & 
                                                    reshape([ 1, 3, 5,&
                                                              2, 3, 5,&
                                                              1, 4, 5,&
                                                              2, 4, 5,&
                                                              1, 3, 6,&
                                                              2, 3, 6,&
                                                              1, 4, 6,&
                                                              2, 4, 6], [NUM_FACES_AT_CORNER_3D, NUM_CORNERS_3D])

  integer(ip), target :: P4EST_EDGES_AT_CORNER_3D(NUM_EDGES_AT_CORNER_3D, NUM_CORNERS_3D) = & 
                                                    reshape([ 1,  5,  9,&
                                                              1,  6, 10,&
                                                              2,  5, 11,&
                                                              2,  6, 12,&
                                                              3,  7,  9,&
                                                              3,  8, 10,&
                                                              4,  7, 11,&
                                                              4,  8, 12], [NUM_EDGES_AT_CORNER_3D, NUM_CORNERS_3D])

  integer(ip), target :: P4EST_CORNER_IN_FACE_3D(NUM_FACES_3D,NUM_CORNERS_3D) = & 
                                                  reshape([ 1,-1, 1,-1, 1,-1,&
                                                           -1, 1, 2,-1, 2,-1,&
                                                            2,-1,-1, 1, 3,-1,&
                                                           -1, 2,-1, 2, 4,-1,&
                                                            3,-1, 3,-1,-1, 1,&
                                                           -1, 3, 4,-1,-1, 2,&
                                                            4,-1,-1, 3,-1, 3,&
                                                           -1, 4,-1, 4,-1, 4 ],[NUM_FACES_3D,NUM_CORNERS_3D])

  integer(ip), target :: P4EST_OPPOSITE_CORNER_3D(NUM_CORNERS_3D) = [ 8, 7, 6, 5, 4, 3, 2, 1 ]
  integer(ip), target :: P4EST_2_FEMPAR_CORNER_3D(NUM_CORNERS_3D) = [ 1, 2, 3, 4, 5, 6, 7, 8 ]
  integer(ip), target :: P4EST_2_FEMPAR_FACE_3D  (NUM_FACES_3D)   = [ 5, 6, 3, 4, 1, 2 ]
  
  integer(ip), parameter :: refinement = 1 
  integer(ip), parameter :: coarsening = -1 
  integer(ip), parameter :: do_nothing = 0 
  
  
  type, extends(cell_iterator_t) :: p4est_cell_iterator_t
    private
    type(p4est_serial_triangulation_t), pointer :: p4est_triangulation => NULL()
  contains
    procedure                            :: create                  => p4est_cell_iterator_create
    procedure                            :: free                    => p4est_cell_iterator_free
    !final                                ::                            p4est_cell_iterator_free_final
    !procedure, non_overridable           :: next                    => p4est_cell_iterator_next
    !procedure, non_overridable           :: first                   => p4est_cell_iterator_first
    procedure                            :: last                    => p4est_cell_iterator_last
    !procedure, non_overridable           :: set_lid                 => p4est_cell_iterator_set_lid
    !procedure, non_overridable, private  :: set_gid                 => p4est_cell_iterator_set_gid
    !procedure, non_overridable, private  :: set_mypart              => p4est_cell_iterator_set_mypart
    !procedure, non_overridable, private  :: get_triangulation       => p4est_cell_iterator_get_triangulation
    procedure                            :: has_finished            => p4est_cell_iterator_has_finished
    procedure                            :: get_reference_fe_geo    => p4est_cell_iterator_get_reference_fe_geo
    procedure                            :: get_reference_fe_geo_id => p4est_cell_iterator_get_reference_fe_geo_id
    procedure                            :: get_coordinates         => p4est_cell_iterator_get_coordinates
    !procedure, non_overridable           :: set_coordinates         => p4est_cell_iterator_set_coordinates
    !procedure, non_overridable           :: get_lid                 => p4est_cell_iterator_get_lid
    !procedure, non_overridable           :: get_gid                 => p4est_cell_iterator_get_gid
    !procedure, non_overridable           :: get_my_part             => p4est_cell_iterator_get_mypart
    !procedure, non_overridable           :: get_my_subpart          => p4est_cell_iterator_get_mysubpart
    !procedure, non_overridable           :: get_my_subpart_lid      => p4est_cell_iterator_get_mysubpart_lid
    procedure                            :: get_set_id              => p4est_cell_iterator_get_set_id
    procedure                            :: get_level               => p4est_cell_iterator_get_level
    procedure                            :: get_num_vefs            => p4est_cell_iterator_get_num_vefs
    !procedure, non_overridable           :: get_num_nodes           => p4est_cell_iterator_get_num_nodes
    !procedure, non_overridable           :: get_node_lid            => p4est_cell_iterator_get_node_lid
    procedure                            :: get_vef_lid             => p4est_cell_iterator_get_vef_lid
    procedure                            :: get_vef_lids            => p4est_cell_iterator_get_vef_lids
    !procedure, non_overridable           :: get_vef_gid             => p4est_cell_iterator_get_vef_gid
    !procedure                            :: find_lpos_vef_lid       => p4est_cell_iterator_find_lpos_vef_lid
    !procedure, non_overridable           :: find_lpos_vef_gid       => p4est_cell_iterator_find_lpos_vef_gid
    !procedure, non_overridable           :: get_vef                 => p4est_cell_iterator_get_vef
    procedure                            :: is_local                => p4est_cell_iterator_is_local
    procedure                            :: is_ghost                => p4est_cell_iterator_is_ghost
    
    procedure                            :: is_ancestor             => p4est_cell_iterator_is_ancestor
    
    procedure                            :: set_for_coarsening      => p4est_cell_iterator_set_for_coarsening
    procedure                            :: set_for_refinement      => p4est_cell_iterator_set_for_refinement
    procedure                            :: set_for_do_nothing      => p4est_cell_iterator_set_for_do_nothing
    procedure                            :: get_transformation_flag => p4est_cell_iterator_get_transformation_flag
  end type p4est_cell_iterator_t
  
  type, extends(vef_iterator_t) :: p4est_vef_iterator_t
    private
    type(p4est_serial_triangulation_t), pointer :: p4est_triangulation => NULL()
  contains
     procedure                           :: create                    => p4est_vef_iterator_create
     procedure                           :: free                      => p4est_vef_iterator_free
     final                               ::                              p4est_vef_iterator_free_final
     !procedure                           :: first                     => p4est_vef_iterator_first
     procedure                           :: next                      => p4est_vef_iterator_next
     !procedure, non_overridable          :: set_lid                   => p4est_vef_iterator_set_lid
     procedure                           :: has_finished              => p4est_vef_iterator_has_finished
     !procedure, non_overridable          :: get_triangulation         => p4est_vef_iterator_get_triangulation
     !procedure, non_overridable          :: get_lid                   => p4est_vef_iterator_get_lid
     !procedure, non_overridable          :: get_gid                   => p4est_vef_iterator_get_gid

     !procedure, non_overridable          :: set_geom_id               => p4est_vef_iterator_set_geom_id
     procedure                           :: set_set_id                => p4est_vef_iterator_set_set_id
     !procedure, non_overridable          :: get_geom_id               => p4est_vef_iterator_get_geom_id
     procedure                           :: get_set_id                => p4est_vef_iterator_get_set_id

     !procedure, non_overridable          :: set_dimension             => p4est_vef_iterator_set_dimension
     !procedure, non_overridable          :: set_it_at_boundary        => p4est_vef_iterator_set_it_at_boundary
     !procedure, non_overridable          :: set_it_as_local           => p4est_vef_iterator_set_it_as_local
     !procedure, non_overridable          :: set_it_as_ghost           => p4est_vef_iterator_set_it_as_ghost
     !procedure, non_overridable          :: set_it_at_interface       => p4est_vef_iterator_set_it_at_interface

     procedure                           :: get_dimension             => p4est_vef_iterator_get_dimension
     procedure                           :: is_at_boundary            => p4est_vef_iterator_is_at_boundary
     procedure                           :: is_local                  => p4est_vef_iterator_is_local
     procedure                           :: is_ghost                  => p4est_vef_iterator_is_ghost
     procedure                           :: is_at_interface           => p4est_vef_iterator_is_at_interface
     !procedure, non_overridable          :: is_face                   => p4est_vef_iterator_is_face
     
     procedure                           :: is_proper                       => p4est_vef_iterator_is_proper
     procedure                           :: is_within_valid_range           => p4est_vef_iterator_is_within_valid_range
     
     procedure                           :: get_num_cells_around            => p4est_vef_iterator_get_num_cells_around
     procedure                           :: get_cell_around                 => p4est_vef_iterator_get_cell_around
     
     ! Improper VEFs-only TBPs
     procedure                           :: get_num_improper_cells_around   => p4est_vef_iterator_get_num_improper_cells_around
     procedure                           :: get_improper_cell_around        => p4est_vef_iterator_get_improper_cell_around
     procedure                           :: get_improper_cell_around_ivef   => p4est_vef_iterator_get_improper_cell_around_ivef
     procedure                           :: get_improper_cell_around_subvef => p4est_vef_iterator_get_improper_cell_around_subvef
  end type p4est_vef_iterator_t
  
  
  
  ! TODO: this data type should extend an abstract triangulation,
  !       and implement its corresponding accessors
  type, extends(serial_triangulation_t) ::  p4est_serial_triangulation_t
    private
    integer(ip) :: p4est_num_cells          = -1
    integer(ip) :: p4est_num_dimensions     = -1
    integer(ip) :: p4est_num_vefs           = -1
    integer(ip) :: num_proper_vefs          = -1 
    integer(ip) :: num_improper_vefs        = -1 
    
    type(hex_lagrangian_reference_fe_t) :: reference_fe_geo
    type(point_t), allocatable          :: per_cell_vertex_coordinates(:)
    
    ! p4est-related data
    type(c_ptr) :: p4est_connectivity = c_null_ptr
    type(c_ptr) :: p4est              = c_null_ptr
    type(c_ptr) :: p4est_mesh         = c_null_ptr
    
    ! TODO: I am pretty sure that a type(c_ptr) :: p4est_ghost
    !       member variable will be needed (at least in the parallel realization)
    
    ! p4est quadrant connectivity (1:NUM_FACES_2D/3D,1:nQuads) => neighbor quadrant
    integer(P4EST_F90_LOCIDX),pointer     :: quad_to_quad(:,:)   => NULL()
    ! p4est face connectivity (1:NUM_FACES_2D/3D,1:nQuads) => neighbor faceId + orientation + non-conform info
    integer(P4EST_F90_QLEVEL),pointer     :: quad_to_face(:,:)   => NULL()   
    ! p4est face connectivity for mortars NUM_SUBFACES_FACE_2D/3D,1:nHalfFaces), (~small sides)
    integer(P4EST_F90_LOCIDX),pointer     :: quad_to_half(:,:)   => NULL()
    integer(P4EST_F90_LOCIDX),pointer     :: quad_to_corner(:,:) => NULL()

    ! TODO: The following 2x member variables should be replaced by our F200X implementation of "std::vector<T>" 
    ! p4est Integer coordinates of first quadrant node (xy/xyz,nQuads)
    integer(P4EST_F90_LOCIDX), allocatable :: quad_coords(:,:)
    ! p4est Integer Level of quadrant
    integer(P4EST_F90_QLEVEL), allocatable :: quad_level(:)
    
    type(std_vector_integer_ip_t)          :: p4est_lst_vefs_lids
    
    type(std_vector_integer_ip_t)          :: p4est_ptr_cells_around_proper_vefs
    type(std_vector_integer_ip_t)          :: p4est_lst_cells_around_proper_vefs
    type(std_vector_integer_ip_t)          :: p4est_ptr_cells_around_improper_vefs
    type(std_vector_integer_ip_t)          :: p4est_lst_cells_around_improper_vefs
    type(std_vector_integer_ip_t)          :: p4est_ptr_improper_cells_around
    type(std_vector_integer_ip_t)          :: p4est_lst_improper_cells_around
    type(std_vector_integer_ip_t)          :: p4est_improper_vefs_improper_cell_around_ivef
    type(std_vector_integer_ip_t)          :: p4est_improper_vefs_improper_cell_around_subvef
    type(std_vector_integer_ip_t)          :: p4est_proper_vefs_dimension
    type(std_vector_integer_ip_t)          :: p4est_improper_vefs_dimension
    type(std_vector_integer_ip_t)          :: p4est_proper_vefs_at_boundary
    type(std_vector_integer_ip_t)          :: p4est_refinement_and_coarsening_flags
    type(std_vector_integer_ip_t)          :: p4est_cell_set_ids
    type(std_vector_integer_ip_t)          :: p4est_proper_vefs_set_ids
    type(std_vector_integer_ip_t)          :: p4est_improper_vefs_set_ids
  contains
  
    ! Getters
    procedure                                   :: get_num_dimensions                            => p4est_serial_triangulation_get_num_dimensions
    procedure                                   :: get_num_cells                                 => p4est_serial_triangulation_get_num_cells
    procedure                                   :: get_num_vefs                                  => p4est_serial_triangulation_get_num_vefs
    procedure                                   :: get_num_proper_vefs                           => p4est_serial_triangulation_get_num_proper_vefs
    procedure                                   :: get_num_improper_vefs                         => p4est_serial_triangulation_get_num_improper_vefs
    procedure                                   :: get_p4est_refinement_and_coarsening_flags     => p4est_st_get_p4est_refinement_and_coarsening_flags
    
    ! Set up related methods
    procedure                                   :: serial_triangulation_create                   => p4est_serial_triangulation_create
    procedure                                   :: free                                          => p4est_serial_triangulation_free
    procedure                 , non_overridable :: refine_and_coarsen                            => p4est_serial_triangulation_refine_and_coarsen
    procedure, private        , non_overridable :: update_p4est_mesh                             => p4est_serial_triangulation_update_p4est_mesh
    procedure, private        , non_overridable :: update_topology_from_p4est_mesh               => p4est_serial_triangulation_update_topology_from_p4est_mesh
    procedure, private        , non_overridable :: get_ptr_vefs_per_cell                         => p4est_serial_triangulation_get_ptr_vefs_per_cell
    procedure, private        , non_overridable :: update_lst_vefs_lids_and_cells_around         => p4est_st_update_lst_vefs_lids_and_cells_around
    procedure, private        , non_overridable :: update_cell_set_ids                           => p4est_st_update_cell_set_ids
    procedure, private        , non_overridable :: update_vef_set_ids                            => p4est_st_update_vef_set_ids
    procedure                 , non_overridable :: std_vector_transform_length_to_header         => p4est_st_std_vector_transform_length_to_header
    procedure, private        , non_overridable :: allocate_and_fill_per_cell_vertex_coordinates => p4est_st_allocate_and_fill_per_cell_vertex_coordinates
    procedure, private        , non_overridable :: free_per_cell_vertex_coordinates              => p4est_st_free_per_cell_vertex_coordinates
    procedure                 , non_overridable :: clear_refinement_and_coarsening_flags         => p4est_st_clear_refinement_and_coarsening_flags
    procedure                 , non_overridable :: clear_cell_set_ids                            => p4est_st_clear_cell_set_ids
    procedure                                   :: fill_cells_set                                => p4est_st_fill_cells_set
    procedure                 , non_overridable :: clear_vef_set_ids                             => p4est_st_clear_vef_set_ids

    ! Cell traversals-related TBPs
    procedure                                   :: create_cell_iterator                  => p4est_create_cell_iterator
    
    ! VEF traversals-related TBPs
    procedure                                   :: create_vef_iterator                   => p4est_create_vef_iterator

#ifndef ENABLE_P4EST
    procedure, non_overridable :: not_enabled_error => p4est_serial_triangulation_not_enabled_error
#endif
  end type p4est_serial_triangulation_t
  
  public :: p4est_serial_triangulation_t, p4est_cell_iterator_t, p4est_vef_iterator_t
  public :: refinement, coarsening, do_nothing
  
contains

#include "sbm_p4est_serial_triangulation.i90"
#include "sbm_p4est_cell_iterator.i90"
#include "sbm_p4est_vef_iterator.i90"
end module p4est_serial_triangulation_names
