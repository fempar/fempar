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
  use stdio_names
  use memor_names
  use p4est_bindings_names
  use base_static_triangulation_names
  use std_vector_integer_ip_names
  
  use FPL

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: NUM_CORNERS_2D           = 4
  integer(ip), parameter :: NUM_FACES_2D             = 4
  integer(ip), parameter :: NUM_SUBFACES_FACE_2D     = 2
  integer(ip), parameter :: NUM_FACE_CORNERS_2D      = 2
  integer(ip), parameter :: NUM_FACES_AT_CORNER_2D   = 2
  integer(ip), parameter :: NUM_VEFS_2D              = NUM_CORNERS_2D+NUM_FACES_2D
  integer(ip), parameter :: P4EST_FACE_CORNERS_2D(NUM_FACE_CORNERS_2D,NUM_FACES_2D) = & 
                                                  reshape([1, 3,&
                                                           2, 4,&  
                                                           1, 2,&
                                                           3, 4], [NUM_FACE_CORNERS_2D,NUM_FACES_2D])
                                                  
  integer(ip), parameter :: P4EST_FACES_AT_CORNER_2D(NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D) = & 
                                                    reshape([1, 3,&
                                                             2, 3,&  
                                                             1, 4,&
                                                             2, 4], [NUM_FACES_AT_CORNER_2D, NUM_CORNERS_2D])
                                                  
  integer(ip), parameter :: P4EST_OPPOSITE_CORNER(NUM_CORNERS_2D) = [ 4, 3, 2, 1 ]
  integer(ip), parameter :: P4EST_2_FEMPAR_CORNER(NUM_CORNERS_2D) = [ 1, 2, 3, 4 ]
  integer(ip), parameter :: P4EST_2_FEMPAR_FACE  (NUM_FACES_2D)   = [ 3, 4, 1, 2 ]
  
  
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
    !procedure, non_overridable           :: get_reference_fe_geo    => p4est_cell_iterator_get_reference_fe_geo
    !procedure, non_overridable           :: get_reference_fe_geo_id => p4est_cell_iterator_get_reference_fe_geo_id
    !procedure, non_overridable           :: get_coordinates         => p4est_cell_iterator_get_coordinates
    !procedure, non_overridable           :: set_coordinates         => p4est_cell_iterator_set_coordinates
    !procedure, non_overridable           :: get_lid                 => p4est_cell_iterator_get_lid
    !procedure, non_overridable           :: get_gid                 => p4est_cell_iterator_get_gid
    !procedure, non_overridable           :: get_my_part             => p4est_cell_iterator_get_mypart
    !procedure, non_overridable           :: get_my_subpart          => p4est_cell_iterator_get_mysubpart
    !procedure, non_overridable           :: get_my_subpart_lid      => p4est_cell_iterator_get_mysubpart_lid
    !procedure, non_overridable           :: get_set_id              => p4est_cell_iterator_get_set_id
    procedure                            :: get_num_vefs            => p4est_cell_iterator_get_num_vefs
    !procedure, non_overridable           :: get_num_nodes           => p4est_cell_iterator_get_num_nodes
    !procedure, non_overridable           :: get_node_lid            => p4est_cell_iterator_get_node_lid
    procedure                            :: get_vef_lid             => p4est_cell_iterator_get_vef_lid
    !procedure, non_overridable           :: get_vef_lids            => p4est_cell_iterator_get_vef_lids
    !procedure, non_overridable           :: get_vef_gid             => p4est_cell_iterator_get_vef_gid
    !procedure                            :: find_lpos_vef_lid       => p4est_cell_iterator_find_lpos_vef_lid
    !procedure, non_overridable           :: find_lpos_vef_gid       => p4est_cell_iterator_find_lpos_vef_gid
    !procedure, non_overridable           :: get_vef                 => p4est_cell_iterator_get_vef
    procedure                            :: is_local                => p4est_cell_iterator_is_local
    procedure                            :: is_ghost                => p4est_cell_iterator_is_ghost
  end type p4est_cell_iterator_t
  
  ! TODO: this data type should extend an abstract triangulation,
  !       and implement its corresponding accessors
  type, extends(serial_triangulation_t) ::  p4est_serial_triangulation_t
    private
    integer(ip) :: p4est_num_cells          = -1
    integer(ip) :: p4est_num_dimensions     = -1
    integer(ip) :: p4est_num_vefs           = -1
    integer(ip) :: num_proper_vefs    = -1 
    integer(ip) :: num_improper_vefs  = -1 
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
  contains
    procedure                                   :: create                                => p4est_serial_triangulation_create
    procedure                                   :: free                                  => p4est_serial_triangulation_free
    procedure                 , non_overridable :: refine_and_coarsen                    => p4est_serial_triangulation_refine_and_coarsen
    procedure, private        , non_overridable :: update_p4est_mesh                     => p4est_serial_triangulation_update_p4est_mesh
    procedure, private        , non_overridable :: update_topology_from_p4est_mesh       => p4est_serial_triangulation_update_topology_from_p4est_mesh
    procedure, private        , non_overridable :: get_ptr_vefs_per_cell                 => p4est_serial_triangulation_get_ptr_vefs_per_cell
    procedure, private        , non_overridable :: update_lst_vefs_lids_and_cells_around => p4est_st_update_lst_vefs_lids_and_cells_around
    procedure, private        , non_overridable :: std_vector_transform_length_to_header => p4est_st_std_vector_transform_length_to_header

    ! Cell traversals-related TBPs
    procedure                                   :: create_cell_iterator                  => p4est_create_cell_iterator
    

#ifndef ENABLE_P4EST
    procedure, non_overridable :: not_enabled_error => p4est_serial_triangulation_not_enabled_error
#endif
  end type p4est_serial_triangulation_t
  
  public :: p4est_serial_triangulation_t
  
contains

#include "sbm_p4est_serial_triangulation.i90"
#include "sbm_p4est_cell_iterator.i90"

end module p4est_serial_triangulation_names
