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
module uniform_hex_mesh_generator
  ! Serial modules
  use types_names
  use memor_names
  ! WE SHOULD TRY TO MINIMIZE THE MODULES/DATA TYPES ON WHICH THIS MODULE DEPENDS.
  ! IN PARTICULAR, IT MUST NOT DEPEND ON THE FEMPAR'S TRIANGULATION MODULES, BUT 
  ! PROVIDE THE DATA AS PLAIN ARRAYS THAT THE TRIANGULATION REQUIRES IN ORDER TO BE
  ! GENERATED

  ! Formely uniform_mesh_descriptor_t. Provided here for convenience (so that legacy 
  ! code can be re-used here AS IS). However, the user of FEMPAR should not be aware of this 
  ! data type. Instead, it would be better that he uses FPL to parametrize the creation of 
  ! the triangulation. I would like to avoid having multiple data types for algorithms parameters 
  ! (e.g., one per algorithm/module), using a single one (i.e., FPL) for all of them ...
  type uniform_hex_mesh_descriptor_t 
    private 
    
  end type uniform_hex_mesh_descriptor_t
  
  ! Convention used within this uniform hex mesh generator module 
  ! for numbering objects within an hexaedron (lexicographic, ZYX)
  ! 2D:
  ! Points code: 1=bottom left; 2=top left;   3=bottom right; 4=top right
  ! Lines code:  1=left line;   2=right line; 3=bottom line;  4=top line 

  ! 3D:
  ! Points code:  1=(0,0,0), 2=(0,0,1), 3=(0,1,0), 4=(0,1,1), 5=(1,0,0),
  !               6=(1,0,1), 7=(1,1,0), 8=(1,1,1)
  ! Lines code:   1=(0,0,z), 2=(0,1,z), 3=(1,0,z), 4=(1,1,z), 5=(0,y,0),
  !               6=(0,y,1), 7=(1,y,0), 8=(1,y,1), 9=(x,0,0), 10=(x,0,1),
  !               11=(x,1,0), 12=(x,1,1)
  ! Surface code: 1=(x,y,0), 2=(x,y,1), 3=(x,0,z), 4=(x,1,z), 5=(0,y,z),
  !               6=(1,y,z)
  ! Formely uniform_conditions_descriptor_t
  type hex_boundary_set_ids_descriptor_t 
    private
    integer(ip)              :: number_dimensions
    integer(ip), allocatable :: boundary_set_ids_vertex(:)
    integer(ip), allocatable :: boundary_set_ids_edges(:)
    integer(ip), allocatable :: boundary_set_ids_faces(:)
  !contains
  !  procedure                :: create                       => hex_boundary_set_ids_descriptor_create
  !  procedure                :: free                         => hex_boundary_set_ids_descriptor_free
  !  procedure                :: set_vertex_boundary_set_id   => hex_boundary_set_ids_descriptor_set_vertex_boundary_set_id
  !  procedure                :: set_edge_boundary_set_id     => hex_boundary_set_ids_descriptor_set_edge_boundary_set_id
  !  procedure                :: set_face_boundary_set_id     => hex_boundary_set_ids_descriptor_set_face_boundary_set_id
  !  procedure                :: get_vertex_boundary_set_id   => hex_boundary_set_ids_descriptor_get_vertex_boundary_set_id
  !  procedure                :: get_edge_boundary_set_id     => hex_boundary_set_ids_descriptor_get_edge_boundary_set_id
  !  procedure                :: get_face_boundary_set_id     => hex_boundary_set_ids_descriptor_get_face_boundary_set_id
  end type hex_boundary_set_ids_descriptor_t
  
  public :: uniform_hex_mesh_descriptor_t
  public :: hex_boundary_set_ids_descriptor_t
  public :: uniform_hex_mesh_generator_generate_fempar_triangulation_arrays
  
contains

  ! Main driver subroutine of this module
  subroutine uniform_hex_mesh_generator_generate_fempar_triangulation_arrays(part_id,                         &
                                                                             uniform_hex_mesh_descriptor,     &
                                                                             hex_boundary_set_ids_descriptor, &
                                                                             num_local_cells,                 &
                                                                             num_local_vefs,                  &
                                                                             ptr_vefs_per_cell,               &
                                                                             lst_vefs_lids,                   &
                                                                             cell_gids,                       &
                                                                             vefs_gids,                       &
                                                                             vefs_boundary_set_ids,           &
                                                                             vefs_geometry_ids,               &
                                                                             vertex_coordinates,              &  
                                                                             num_itfc_cells, &
                                                                             lst_itfc_cells, &
                                                                             ptr_ext_neighs_per_itfc_cell, &
                                                                             lst_ext_neighs_gids, &
                                                                             lst_ext_neighs_part_ids)
                                                                             
    implicit none
    integer(ip)                            , intent(in)    :: part_id
    type(uniform_hex_mesh_descriptor_t)    , intent(in)    :: uniform_hex_mesh_descriptor
    type(hex_boundary_set_ids_descriptor_t), intent(in)    :: hex_boundary_set_ids_descriptor
    integer(ip)                            , intent(out)   :: num_local_cells
    integer(ip)                            , intent(out)   :: num_local_vefs
    integer(ip)           , allocatable    , intent(inout) :: ptr_vefs_per_cell(:)            ! Size = num_local_cells + 1
    integer(ip)           , allocatable    , intent(inout) :: lst_vefs_lids(:)                ! Size = ptr_vefs_per_cell(num_local_cells+1)-1
    integer(igp)          , allocatable    , intent(inout) :: cell_gids(:)                    ! Size = num_local_cells 
    integer(igp)          , allocatable    , intent(inout) :: vefs_gids(:)                    ! Size = num_local_vefs
    integer(ip)           , allocatable    , intent(inout) :: vefs_boundary_set_ids(:)        ! Size = num_local_vefs
    integer(ip)           , allocatable    , intent(inout) :: vefs_geometry_ids(:)            ! Size = num_local_vefs
    real(rp)              , allocatable    , intent(inout) :: vertex_coordinates(:,:)         ! Size = (number_dimensions, num_local_vertices)
    integer(ip) , optional, allocatable    , intent(out)   :: num_itfc_cells                  ! NONE or ALL OPTIONAL ARGUMENTS MUST BE PRESENT LOGIC
    integer(ip) , optional, allocatable    , intent(inout) :: lst_itfc_cells(:)               ! Size = num_itfc_cells 
    integer(ip) , optional, allocatable    , intent(inout) :: ptr_ext_neighs_per_itfc_cell(:) ! Size = num_itfc_cells + 1
    integer(ip) , optional, allocatable    , intent(inout) :: lst_ext_neighs_gids(:)          ! Size = ptr_ext_neighs_per_itfc_cell(num_itfc_cells + 1)-1
    integer(ip) , optional, allocatable    , intent(inout) :: lst_ext_neighs_part_ids(:)      ! Size = ptr_ext_neighs_per_itfc_cell(num_itfc_cells + 1)-1
  end subroutine uniform_hex_mesh_generator_generate_fempar_triangulation_arrays


end module uniform_hex_mesh_generator
