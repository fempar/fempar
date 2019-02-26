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


! Module file taken "AS-IS" and adapted to FEMPAR needs from p4est_binding.f90 
! https://github.com/cburstedde/hopest 4feed803f0c61564203a7bc3f2ca1a6adb63d3cd

!> date:   10-May-2017                                                                           
!> author: Alberto F. Martin
!> summary: F90 <-> C wrapper subroutines for the p4est (2D) subroutines
!> This module provides F90 <-> C wrapper subroutines for the p4est (2D) subroutines
module p4est_bindings_names
  use types_names
  use, intrinsic :: iso_c_binding
  implicit none
  
    ! ONLY change these parameter values if p4est changes its types
  integer(ip), parameter :: P4EST_F90_TOPIDX = C_INT32_T
  integer(ip), parameter :: P4EST_F90_QCOORD = C_INT32_T
  integer(ip), parameter :: P4EST_F90_LOCIDX = C_INT32_T
  integer(ip), parameter :: P4EST_F90_GLOIDX = C_INT64_T
  integer(ip), parameter :: P4EST_F90_QLEVEL = C_INT8_T

#ifdef ENABLE_P4EST   
  interface 
     !=================================================================================================================================
     !> summary: Initializes (sc_)MPI, sc, and p4est
     !>
     !> @note (sc_)MPI is only initialized if p4est/sc was compiled without MPI support
     !=================================================================================================================================
     subroutine F90_p4est_init(Fcomm,sc_log_level) bind(c,name="F90_p4est_init")
       use, intrinsic :: iso_c_binding
       implicit none
       integer, value, intent(in) :: Fcomm
       integer, value, intent(in) :: sc_log_level
     end subroutine F90_p4est_init

     subroutine F90_p4est_finalize() bind(c,name="F90_p4est_finalize")
       implicit none
     end subroutine F90_p4est_finalize
     
     !=================================================================================================================================
     !> summary: Creates p4est connectivity corresponding to unit square domain
     !=================================================================================================================================
     subroutine F90_p4est_connectivity_new_unitsquare(p4est_connectivity) bind(c,name="F90_p4est_connectivity_new_unitsquare")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p4est_connectivity
     end subroutine F90_p4est_connectivity_new_unitsquare

     !=================================================================================================================================
     !> summary: Creates p8est connectivity corresponding to unit cube domain
     !=================================================================================================================================
     subroutine F90_p8est_connectivity_new_unitcube(p8est_connectivity) bind(c,name="F90_p8est_connectivity_new_unitcube")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p8est_connectivity
     end subroutine F90_p8est_connectivity_new_unitcube
     
     !=================================================================================================================================
     !> summary: Scales p4est connectivity with bounding box limits
     !=================================================================================================================================
     subroutine F90_p4est_connectivity_set_bounding_box_limits(p4est_connectivity,bounding_box_limits) bind(c,name="F90_p4est_connectivity_set_bounding_box_limits")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr)       , intent(inout)  :: p4est_connectivity
       type(c_ptr), value, intent(in)     :: bounding_box_limits
     end subroutine F90_p4est_connectivity_set_bounding_box_limits
     !=================================================================================================================================
     !> summary: Creates unrefined p4est (it will contain a single root octant)
     !=================================================================================================================================
     subroutine F90_p4est_new(mpi_comm, p4est_connectivity, p4est) bind(c,name="F90_p4est_new")
       use, intrinsic :: iso_c_binding
       implicit none
       integer    , value, intent(in)     :: mpi_comm
       type(c_ptr), value, intent(in)     :: p4est_connectivity
       type(c_ptr)       , intent(inout)  :: p4est
     end subroutine F90_p4est_new

     !=================================================================================================================================
     !> summary: Creates unrefined p8est (it will contain a single root octant)
     !=================================================================================================================================
     subroutine F90_p8est_new(mpi_comm, p8est_connectivity, p8est) bind(c,name="F90_p8est_new")
       use, intrinsic :: iso_c_binding
       implicit none
       integer    , value, intent(in)     :: mpi_comm
       type(c_ptr), value, intent(in)     :: p8est_connectivity
       type(c_ptr)       , intent(inout)  :: p8est
     end subroutine F90_p8est_new
     
     !=================================================================================================================================
     !> summary: set user_pointer member variable of p4est_t struct to input integer user_data array
     !=================================================================================================================================
     subroutine F90_p4est_set_user_pointer(user_data, p4est) bind(c,name="F90_p4est_set_user_pointer")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) , value, intent(in)  :: user_data
       type(c_ptr) , value, intent(in)  :: p4est
     end subroutine F90_p4est_set_user_pointer

     !=================================================================================================================================
     !> summary: set user_pointer member variable of p8est_t struct to input integer user_data array
     !=================================================================================================================================
     subroutine F90_p8est_set_user_pointer(user_data, p8est) bind(c,name="F90_p8est_set_user_pointer")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) , value, intent(in)  :: user_data
       type(c_ptr) , value, intent(in)  :: p8est
     end subroutine F90_p8est_set_user_pointer
     
     !=================================================================================================================================
     !> summary: Creates p4est_mesh from p4est
     !=================================================================================================================================
     subroutine F90_p4est_mesh_new(p4est, p4est_ghost, p4est_mesh) bind(c,name="F90_p4est_mesh_new")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in)     :: p4est
       type(c_ptr)       , intent(inout)  :: p4est_ghost
       type(c_ptr)       , intent(inout)  :: p4est_mesh
     end subroutine F90_p4est_mesh_new
     
     !=================================================================================================================================
     !> summary: Creates p8est_mesh from p8est
     !=================================================================================================================================
     subroutine F90_p8est_mesh_new(p8est, p8est_ghost, p8est_mesh) bind(c,name="F90_p8est_mesh_new")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value, intent(in)     :: p8est
       type(c_ptr)       , intent(inout)  :: p8est_ghost
       type(c_ptr)       , intent(inout)  :: p8est_mesh
     end subroutine F90_p8est_mesh_new
     
     !=================================================================================================================================
     !> summary: Gets mesh info from (p4est, p4est_mesh)
     !=================================================================================================================================
     subroutine F90_p4est_get_mesh_info (p4est, &
                                         p4est_mesh, &
                                         local_num_quadrants, &
                                         ghost_num_quadrants, &
                                         global_num_quadrants, &
                                         global_first_quadrant, &
                                         num_half_faces) bind(c,name="F90_p4est_get_mesh_info")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX, P4EST_F90_GLOIDX
       type(c_ptr), value       , intent(in)     :: p4est
       type(c_ptr), value       , intent(in)     :: p4est_mesh
       integer(P4EST_F90_LOCIDX), intent(out)    :: local_num_quadrants
       integer(P4EST_F90_LOCIDX), intent(out)    :: ghost_num_quadrants
       integer(P4EST_F90_GLOIDX), intent(out)    :: global_num_quadrants
       integer(P4EST_F90_GLOIDX), intent(out)    :: global_first_quadrant(*)
       integer(P4EST_F90_LOCIDX), intent(out)    :: num_half_faces
     end subroutine F90_p4est_get_mesh_info 
     
     !=================================================================================================================================
     !> summary: Gets mesh info from (p8est, p8est_mesh)
     !=================================================================================================================================
     subroutine F90_p8est_get_mesh_info (p8est, &
                                         p8est_mesh, &
                                         local_num_quadrants, &
                                         ghost_num_quadrants, &
                                         global_num_quadrants, &
                                         global_first_quadrant, &
                                         num_half_faces) bind(c,name="F90_p8est_get_mesh_info")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX, P4EST_F90_GLOIDX
       type(c_ptr), value       , intent(in)     :: p8est
       type(c_ptr), value       , intent(in)     :: p8est_mesh
       integer(P4EST_F90_LOCIDX), intent(out)    :: local_num_quadrants
       integer(P4EST_F90_LOCIDX), intent(out)    :: ghost_num_quadrants
       integer(P4EST_F90_GLOIDX), intent(out)    :: global_num_quadrants
       integer(P4EST_F90_GLOIDX), intent(out)    :: global_first_quadrant(*)
       integer(P4EST_F90_LOCIDX), intent(out)    :: num_half_faces
     end subroutine F90_p8est_get_mesh_info 
     
     !=================================================================================================================================
     !> summary: Gets mesh topology arrays from (p4est, p4est_mesh)
     !=================================================================================================================================
     subroutine F90_p4est_get_mesh_topology_arrays(p4est, &
                                                   p4est_mesh, &
                                                   p4est_ghost, &
                                                   quad_to_quad, &
                                                   quad_to_face, &
                                                   quad_to_half, &
                                                   quad_to_corner, &
                                                   quadcoords, &
                                                   quadlevel) bind(c,name="F90_p4est_get_mesh_topology_arrays")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est
       type(c_ptr), value       , intent(in)     :: p4est_mesh
       type(c_ptr), value       , intent(in)     :: p4est_ghost
       type(c_ptr)              , intent(out)    :: quad_to_quad
       type(c_ptr)              , intent(out)    :: quad_to_face
       type(c_ptr)              , intent(out)    :: quad_to_half
       type(c_ptr)              , intent(out)    :: quad_to_corner
       integer(P4EST_F90_QCOORD), intent(out)    :: quadcoords(2,*)
       integer(P4EST_F90_QLEVEL), intent(out)    :: quadlevel(*)
     end subroutine F90_p4est_get_mesh_topology_arrays

     !=================================================================================================================================
     !> summary: Gets mesh topology arrays from (p8est, p8est_mesh)
     !=================================================================================================================================
     subroutine F90_p8est_get_mesh_topology_arrays(p8est, &
                                                   p8est_mesh, &
                                                   p8est_ghost, &
                                                   quad_to_quad, &
                                                   quad_to_face, &
                                                   quad_to_half, &
                                                   quad_to_quad_by_edge,&
                                                   quad_to_edge,&
                                                   num_half_edges, &
                                                   quad_to_half_by_edge, & 
                                                   quad_to_corner, &
                                                   quadcoords, &
                                                   quadlevel) bind(c,name="F90_p8est_get_mesh_topology_arrays")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL, P4EST_F90_LOCIDX
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est
       type(c_ptr), value       , intent(in)     :: p8est_mesh
       type(c_ptr), value       , intent(in)     :: p8est_ghost
       type(c_ptr)              , intent(out)    :: quad_to_quad
       type(c_ptr)              , intent(out)    :: quad_to_face
       type(c_ptr)              , intent(out)    :: quad_to_half
       integer(P4EST_F90_LOCIDX), intent(out)    :: quad_to_quad_by_edge(12,*)
       integer(P4EST_F90_QLEVEL), intent(out)    :: quad_to_edge(12,*)
       integer(P4EST_F90_LOCIDX), intent(out)    :: num_half_edges
       type(c_ptr)              , intent(inout)  :: quad_to_half_by_edge
       type(c_ptr)              , intent(out)    :: quad_to_corner
       integer(P4EST_F90_QCOORD), intent(out)    :: quadcoords(3,*)
       integer(P4EST_F90_QLEVEL), intent(out)    :: quadlevel(*)
     end subroutine F90_p8est_get_mesh_topology_arrays
     
     
     !=================================================================================================================================
     !> summary: Refines in place the p4est data structure
     !=================================================================================================================================
     subroutine F90_p4est_refine(p4est) bind(c,name="F90_p4est_refine")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est
     end subroutine F90_p4est_refine

     !=================================================================================================================================
     !> summary: Refines in place the p8est data structure
     !=================================================================================================================================
     subroutine F90_p8est_refine(p8est) bind(c,name="F90_p8est_refine")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est
     end subroutine F90_p8est_refine
     
     !=================================================================================================================================
     !> summary: Coarsens in place the p4est data structure
     !=================================================================================================================================
     subroutine F90_p4est_coarsen(p4est) bind(c,name="F90_p4est_coarsen")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est
     end subroutine F90_p4est_coarsen

     !=================================================================================================================================
     !> summary: Coarsens in place the p8est data structure
     !=================================================================================================================================
     subroutine F90_p8est_coarsen(p8est) bind(c,name="F90_p8est_coarsen")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est
     end subroutine F90_p8est_coarsen
     
     !=================================================================================================================================
     !> summary: Makes a deep copy of p4est_input into p4est_output
     !=================================================================================================================================
     subroutine F90_p4est_copy(p4est_input, p4est_output) bind(c,name="F90_p4est_copy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est_input
       type(c_ptr)              , intent(inout)  :: p4est_output
     end subroutine F90_p4est_copy

     !=================================================================================================================================
     !> summary: Makes a deep copy of p8est_input into p8est_output
     !=================================================================================================================================
     subroutine F90_p8est_copy(p8est_input, p8est_output) bind(c,name="F90_p8est_copy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est_input
       type(c_ptr)              , intent(inout)  :: p8est_output
     end subroutine F90_p8est_copy
     
     !=================================================================================================================================
     !> summary: 2:1 (FULL) balance the size differences of neighboring elements in a forest
     !=================================================================================================================================
     subroutine F90_p4est_balance(p4est) bind(c,name="F90_p4est_balance")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est
     end subroutine F90_p4est_balance

     !=================================================================================================================================
     !> summary: 2:1 (FULL) balance the size differences of neighboring elements in a forest
     !=================================================================================================================================
     subroutine F90_p8est_balance(p8est) bind(c,name="F90_p8est_balance")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est
     end subroutine F90_p8est_balance
     
     !=================================================================================================================================
     !> summary: Equally partition the forest. The forest will be partitioned between processors such that they have an approximately 
     !>          equal number of quadrants (or sum of weights).
     !=================================================================================================================================
     subroutine F90_p4est_partition(p4est) bind(c,name="F90_p4est_partition")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est
     end subroutine F90_p4est_partition

     !=================================================================================================================================
     !> summary: Equally partition the forest. The forest will be partitioned between processors such that they have an approximately 
     !>          equal number of quadrants (or sum of weights).
     !=================================================================================================================================
     subroutine F90_p8est_partition(p8est) bind(c,name="F90_p8est_partition")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est
     end subroutine F90_p8est_partition
     
     
     !=================================================================================================================================
     !> summary: compares p4est_old with p4est_new and updates refinement+coarsening flags of the former accordingly to how  
     !>          the former has been transformed into the latter via refine+coarsen+balance
     !=================================================================================================================================
     subroutine F90_p4est_update_refinement_and_coarsening_flags(p4est_old, p4est_new) bind(c,name="F90_p4est_update_refinement_and_coarsening_flags")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est_old
       type(c_ptr), value       , intent(in)     :: p4est_new
     end subroutine F90_p4est_update_refinement_and_coarsening_flags

     !=================================================================================================================================
     !> summary: compares p8est_old with p8est_new and updates refinement+coarsening flags of the former accordingly to how  
     !>          the former has been transformed into the latter via refine+coarsen+balance
     !=================================================================================================================================
     subroutine F90_p8est_update_refinement_and_coarsening_flags(p8est_old, p8est_new) bind(c,name="F90_p8est_update_refinement_and_coarsening_flags")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est_old
       type(c_ptr), value       , intent(in)     :: p8est_new
     end subroutine F90_p8est_update_refinement_and_coarsening_flags
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p4est_connectivity_t
     !=================================================================================================================================
     subroutine F90_p4est_connectivity_destroy(p4est_connectivity) bind(c, name="F90_p4est_connectivity_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p4est_connectivity
     end subroutine F90_p4est_connectivity_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p8est_connectivity_t
     !=================================================================================================================================
     subroutine F90_p8est_connectivity_destroy(p8est_connectivity) bind(c, name="F90_p8est_connectivity_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p8est_connectivity
     end subroutine F90_p8est_connectivity_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p4est_t
     !=================================================================================================================================
     subroutine F90_p4est_destroy(p4est) bind(c, name="F90_p4est_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p4est
     end subroutine F90_p4est_destroy

     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p8est_t
     !=================================================================================================================================
     subroutine F90_p8est_destroy(p8est) bind(c, name="F90_p8est_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p8est
     end subroutine F90_p8est_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p4mesh_t
     !=================================================================================================================================
     subroutine F90_p4est_mesh_destroy(p4est_mesh) bind(c, name="F90_p4est_mesh_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p4est_mesh
     end subroutine F90_p4est_mesh_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p8mesh_t
     !=================================================================================================================================
     subroutine F90_p8est_mesh_destroy(p8est_mesh) bind(c, name="F90_p8est_mesh_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p8est_mesh
     end subroutine F90_p8est_mesh_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p4ghost_t
     !=================================================================================================================================
     subroutine F90_p4est_ghost_destroy(p4est_ghost) bind(c, name="F90_p4est_ghost_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p4est_ghost
     end subroutine F90_p4est_ghost_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory involved in p8ghost_t
     !=================================================================================================================================
     subroutine F90_p8est_ghost_destroy(p8est_ghost) bind(c, name="F90_p8est_ghost_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: p8est_ghost
     end subroutine F90_p8est_ghost_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory in a p4est_locidx_t work array
     !=================================================================================================================================
     subroutine F90_p4est_locidx_buffer_destroy(buffer) bind(c, name="F90_p4est_locidx_buffer_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: buffer
     end subroutine F90_p4est_locidx_buffer_destroy
     
     !=================================================================================================================================
     !> summary: Frees all dynamic memory in a int8_t work array
     !=================================================================================================================================
     subroutine F90_p4est_int8_buffer_destroy(buffer) bind(c, name="F90_p4est_int8_buffer_destroy")
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), intent(inout)  :: buffer
     end subroutine F90_p4est_int8_buffer_destroy
     
     !===========================================================================================================================
     !> summary: Provides in vxyz the coordinates in real space of a vertex given a quadrant(x,y,l) and corner ID within quadrant
     !===========================================================================================================================
     subroutine F90_p4est_get_quadrant_vertex_coordinates(connectivity, &
                                                          treeid, &
                                                          x, &
                                                          y, &
                                                          level, &
                                                          corner, &
                                                          vxyz) bind(c, name="F90_p4est_get_quadrant_vertex_coordinates")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_TOPIDX, P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)               , value, intent(in)     :: connectivity
       integer(P4EST_F90_TOPIDX) , value, intent(in)     :: treeid
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: y
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: level
       integer(c_int)            , value, intent(in)     :: corner
       real(c_double)                   , intent(inout)  :: vxyz(3)
     end subroutine F90_p4est_get_quadrant_vertex_coordinates

     !===========================================================================================================================
     !> summary: Provides in vxyz the coordinates in real space of a vertex given a quadrant(x,y,l) and corner ID within quadrant
     !===========================================================================================================================
     subroutine F90_p8est_get_quadrant_vertex_coordinates(connectivity, &
                                                          treeid, &
                                                          x, &
                                                          y, &
                                                          z, &
                                                          level, &
                                                          corner, &
                                                          vxyz) bind(c, name="F90_p8est_get_quadrant_vertex_coordinates")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_TOPIDX, P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)               , value, intent(in)     :: connectivity
       integer(P4EST_F90_TOPIDX) , value, intent(in)     :: treeid
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: y
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: z
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: level
       integer(c_int)            , value, intent(in)     :: corner
       real(c_double)                   , intent(inout)  :: vxyz(3)
     end subroutine F90_p8est_get_quadrant_vertex_coordinates
     
     !=================================================================================================================================
     !> summary: Return true if quadrant q1 is unequal to and an ancestor of another quadrant q2.
     !=================================================================================================================================
     function F90_p4est_is_ancestor(q1_x,q1_y,q1_level,q2_x,q2_y,q2_level) bind(c, name="F90_p4est_is_ancestor")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q1_x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q1_y
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: q1_level
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q2_x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q2_y
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: q2_level
       integer(c_int)                                    :: F90_p4est_is_ancestor
     end function F90_p4est_is_ancestor
     
     !=================================================================================================================================
     !> summary: Return true if quadrant q1 describes the same quadrant as q2.
     !=================================================================================================================================
     function F90_p4est_is_equal(q1_x,q1_y,q1_level,q2_x,q2_y,q2_level) bind(c, name="F90_p4est_is_equal")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q1_x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q1_y
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: q1_level
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q2_x
       integer(P4EST_F90_QCOORD) , value, intent(in)     :: q2_y
       integer(P4EST_F90_QLEVEL) , value, intent(in)     :: q2_level
       integer(c_int)                                    :: F90_p4est_is_equal
     end function F90_p4est_is_equal
     
     !=================================================================================================================================
     !> summary: Set quadrant q Morton indices based on linear position id in uniform grid of level l.
     !=================================================================================================================================
     subroutine F90_p4est_quadrant_set_morton(l,id,q_x,q_y) bind(c, name="F90_p4est_quadrant_set_morton")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD
       implicit none
       integer(c_int)           , value, intent(in)     :: l
       integer(c_int64_t)       , value, intent(in)     :: id
       integer(P4EST_F90_QCOORD)       , intent(out)    :: q_x
       integer(P4EST_F90_QCOORD)       , intent(out)    :: q_y
     end subroutine F90_p4est_quadrant_set_morton
     
     !=================================================================================================================================
     !> summary: Compute the position of this child within its sibling. Returns its child id in 0..3
     !=================================================================================================================================
     function F90_p4est_quadrant_child_id(q_x,q_y,level) bind(c, name="F90_p4est_quadrant_child_id")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD), value, intent(in)    :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)    :: q_y
       integer(P4EST_F90_QLEVEL), value, intent(in)    :: level
       integer(c_int) :: F90_p4est_quadrant_child_id
     end function F90_p4est_quadrant_child_id
     
     !=================================================================================================================================
     !> summary: Compute the position of this child within its sibling. Returns its child id in 0..7
     !=================================================================================================================================
     function F90_p8est_quadrant_child_id(q_x,q_y,q_z,level) bind(c, name="F90_p8est_quadrant_child_id")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD), value, intent(in)    :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)    :: q_y
       integer(P4EST_F90_QCOORD), value, intent(in)    :: q_z
       integer(P4EST_F90_QLEVEL), value, intent(in)    :: level
        integer(c_int) :: F90_p8est_quadrant_child_id
     end function F90_p8est_quadrant_child_id
     
     !=================================================================================================================================
     !> summary: Fills the remote processor identifiers of the ghost cells of this processor in 1..P
     !=================================================================================================================================
     subroutine F90_p4est_fill_ghost_procs(p4est_ghost,ghost_procs) bind(c, name="F90_p4est_fill_ghost_procs")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr)       , value, intent(in) :: p4est_ghost
       integer(P4EST_F90_LOCIDX), intent(in) :: ghost_procs(*)
     end subroutine F90_p4est_fill_ghost_procs
     
     !=================================================================================================================================
     !> summary: Fills the remote processor identifiers of the ghost cells of this processor in 1..P
     !=================================================================================================================================
     subroutine F90_p8est_fill_ghost_procs(p8est_ghost,ghost_procs) bind(c, name="F90_p8est_fill_ghost_procs")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr)       , value, intent(in) :: p8est_ghost
       integer(P4EST_F90_LOCIDX), intent(in) :: ghost_procs(*)
     end subroutine F90_p8est_fill_ghost_procs

     !=================================================================================================================================
     !> summary: Fills the global identifiers of the ghost cells of this processor
     !=================================================================================================================================
     subroutine F90_p4est_fill_ghost_ggids(p4est_ghost,first_global_quadrant,ghost_ggids) bind(c, name="F90_p4est_fill_ghost_ggids")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX, P4EST_F90_GLOIDX
       implicit none
       type(c_ptr)       , value, intent(in)  :: p4est_ghost
       integer(P4EST_F90_GLOIDX), intent(in)  :: first_global_quadrant(*)
       integer(P4EST_F90_GLOIDX), intent(out) :: ghost_ggids(*)
     end subroutine F90_p4est_fill_ghost_ggids
     
     !=================================================================================================================================
     !> summary: Fills the global identifiers of the ghost cells of this processor
     !=================================================================================================================================
     subroutine F90_p8est_fill_ghost_ggids(p8est_ghost,first_global_quadrant,ghost_ggids) bind(c, name="F90_p8est_fill_ghost_ggids")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX, P4EST_F90_GLOIDX
       implicit none
       type(c_ptr)       , value, intent(in)  :: p8est_ghost
       integer(P4EST_F90_GLOIDX), intent(in)  :: first_global_quadrant(*)
       integer(P4EST_F90_GLOIDX), intent(out) :: ghost_ggids(*)
     end subroutine F90_p8est_fill_ghost_ggids
     
     subroutine F90_p4est_allocate_and_fill_cell_import_raw_arrays(p4est, &
                                                                   p4est_ghost, &
                                                                   num_neighbours, &
                                                                   neighbour_ids, &
                                                                   rcv_ptrs, &
                                                                   rcv_leids, & 
                                                                   snd_ptrs, &
                                                                   snd_leids) bind(c, name="F90_p4est_allocate_and_fill_cell_import_raw_arrays")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr)              , value, intent(in)   :: p4est
       type(c_ptr)              , value, intent(in)   :: p4est_ghost
       integer(P4EST_F90_LOCIDX)       , intent(out)  :: num_neighbours
       type(c_ptr)                     , intent(out)  :: neighbour_ids
       type(c_ptr)                     , intent(out)  :: rcv_ptrs
       type(c_ptr)                     , intent(out)  :: rcv_leids
       type(c_ptr)                     , intent(out)  :: snd_ptrs
       type(c_ptr)                     , intent(out)  :: snd_leids
     end subroutine F90_p4est_allocate_and_fill_cell_import_raw_arrays
     
     subroutine F90_p8est_allocate_and_fill_cell_import_raw_arrays(p4est, &
                                                                   p4est_ghost, &
                                                                   num_neighbours, &
                                                                   neighbour_ids, &
                                                                   rcv_ptrs, &
                                                                   rcv_leids, & 
                                                                   snd_ptrs, &
                                                                   snd_leids) bind(c, name="F90_p8est_allocate_and_fill_cell_import_raw_arrays")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr)              , value, intent(in)   :: p4est
       type(c_ptr)              , value, intent(in)   :: p4est_ghost
       integer(P4EST_F90_LOCIDX)       , intent(out)  :: num_neighbours
       type(c_ptr)                     , intent(out)  :: neighbour_ids
       type(c_ptr)                     , intent(out)  :: rcv_ptrs
       type(c_ptr)                     , intent(out)  :: rcv_leids
       type(c_ptr)                     , intent(out)  :: snd_ptrs
       type(c_ptr)                     , intent(out)  :: snd_leids
     end subroutine F90_p8est_allocate_and_fill_cell_import_raw_arrays
     
     
     subroutine F90_p4est_compute_migration_control_data (p4est_old, p4est_new, num_ranks, lst_ranks, ptr_ranks, local_ids, old2new) &
          bind(c, name="F90_p4est_compute_migration_control_data")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est_old
       type(c_ptr), value       , intent(in)     :: p4est_new
       integer(P4EST_F90_LOCIDX), intent(out)    :: num_ranks
       type(c_ptr)              , intent(inout)  :: lst_ranks
       type(c_ptr)              , intent(inout)  :: ptr_ranks
       type(c_ptr)              , intent(inout)  :: local_ids
       type(c_ptr)              , intent(inout)  :: old2new
    end subroutine F90_p4est_compute_migration_control_data
    
    subroutine F90_p8est_compute_migration_control_data (p8est_old, p8est_new, num_ranks, lst_ranks, ptr_ranks, local_ids, old2new) &
          bind(c, name="F90_p8est_compute_migration_control_data")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est_old
       type(c_ptr), value       , intent(in)     :: p8est_new
       integer(P4EST_F90_LOCIDX), intent(out)    :: num_ranks
       type(c_ptr)              , intent(inout)  :: lst_ranks
       type(c_ptr)              , intent(inout)  :: ptr_ranks
       type(c_ptr)              , intent(inout)  :: local_ids
       type(c_ptr)              , intent(inout)  :: old2new
    end subroutine F90_p8est_compute_migration_control_data

    subroutine F90_p4est_fill_proc_offsets_and_ghost_gids_remote_neighbours ( p4est_ghost, proc_offsets, ghost_gids_remote_neighbours ) &
          bind(c, name="F90_p4est_fill_proc_offsets_and_ghost_gids_remote_neighbours")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr), value       , intent(in)     :: p4est_ghost
       integer(P4EST_F90_LOCIDX), intent(out)    :: proc_offsets(*) 
       integer(P4EST_F90_LOCIDX), intent(out)    :: ghost_gids_remote_neighbours(*) 
    end subroutine F90_p4est_fill_proc_offsets_and_ghost_gids_remote_neighbours 
    
    subroutine F90_p8est_fill_proc_offsets_and_ghost_gids_remote_neighbours ( p8est_ghost, proc_offsets, ghost_gids_remote_neighbours ) &
          bind(c, name="F90_p8est_fill_proc_offsets_and_ghost_gids_remote_neighbours")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_LOCIDX
       implicit none
       type(c_ptr), value       , intent(in)     :: p8est_ghost
       integer(P4EST_F90_LOCIDX), intent(out)    :: proc_offsets(*) 
       integer(P4EST_F90_LOCIDX), intent(out)    :: ghost_gids_remote_neighbours(*) 
    end subroutine F90_p8est_fill_proc_offsets_and_ghost_gids_remote_neighbours 
    
    !=================================================================================================================================
    !> summary: Computes the face neighbor of a quadrant
    !=================================================================================================================================
    subroutine F90_p4est_quadrant_face_neighbor(q_x, q_y, q_level, face, n_x, n_y, n_level) &
       bind(c, name="F90_p4est_quadrant_face_neighbor")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int)           , value, intent(in)  :: face
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_x
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_y
       integer(P4EST_F90_QLEVEL)       , intent(out) :: n_level
    end subroutine F90_p4est_quadrant_face_neighbor
    
    !=================================================================================================================================
    !> summary: Computes the face neighbor of a quadrant
    !=================================================================================================================================
    subroutine F90_p8est_quadrant_face_neighbor(q_x, q_y, q_z, q_level, face, n_x, n_y, n_z, n_level) &
       bind(c, name="F90_p8est_quadrant_face_neighbor")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_z
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int)           , value, intent(in)  :: face
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_x
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_y
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_z
       integer(P4EST_F90_QLEVEL)       , intent(out) :: n_level
    end subroutine F90_p8est_quadrant_face_neighbor
    
    !=================================================================================================================================
    !> summary: Computes the edge neighbor of a quadrant
    !=================================================================================================================================
    subroutine F90_p8est_quadrant_edge_neighbor(q_x, q_y, q_z, q_level, edge, n_x, n_y, n_z, n_level) &
       bind(c, name="F90_p8est_quadrant_edge_neighbor")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_z
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int)           , value, intent(in)  :: edge
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_x
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_y
       integer(P4EST_F90_QCOORD)       , intent(out) :: n_z
       integer(P4EST_F90_QLEVEL)       , intent(out) :: n_level
    end subroutine F90_p8est_quadrant_edge_neighbor
    
    !=================================================================================================================================
    !> summary: Conduct binary search for exact match on a range of the p4est layer
    !=================================================================================================================================
    function F90_p4est_bsearch (p4est, q_x, q_y, q_level) &
       bind(c, name="F90_p4est_bsearch")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)              , value, intent(in)  :: p4est
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int) :: F90_p4est_bsearch
    end function F90_p4est_bsearch
    
    !=================================================================================================================================
    !> summary: Conduct binary search for exact match on a range of the p8est layer
    !=================================================================================================================================
    function F90_p8est_bsearch (p8est, q_x, q_y, q_z, q_level) &
       bind(c, name="F90_p8est_bsearch")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)              , value, intent(in)  :: p8est
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_z
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int) :: F90_p8est_bsearch
    end function F90_p8est_bsearch
    
    !=================================================================================================================================
    !> summary: Conduct binary search for exact match on a range of the ghost layer
    !=================================================================================================================================
    function F90_p4est_ghost_bsearch (ghost, q_x, q_y, q_level) & 
      bind(c, name="F90_p4est_ghost_bsearch")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)              , value, intent(in)  :: ghost
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int) :: F90_p4est_ghost_bsearch
    end function F90_p4est_ghost_bsearch
    
    !=================================================================================================================================
    !> summary: Conduct binary search for exact match on a range of the ghost layer
    !=================================================================================================================================
    function F90_p8est_ghost_bsearch (ghost, q_x, q_y, q_z, q_level) & 
      bind(c, name="F90_p8est_ghost_bsearch")
       use, intrinsic :: iso_c_binding
       import :: P4EST_F90_QCOORD, P4EST_F90_QLEVEL
       implicit none
       type(c_ptr)              , value, intent(in)  :: ghost
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_x
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_y
       integer(P4EST_F90_QCOORD), value, intent(in)  :: q_z
       integer(P4EST_F90_QLEVEL), value, intent(in)  :: q_level
       integer(c_int) :: F90_p8est_ghost_bsearch
    end function F90_p8est_ghost_bsearch
    
     !subroutine p4_savemesh(filename, p4est) bind(c)
     !  !=================================================================================================================================
     !  ! save the p4est data  to a p4est state file 
     !  !=================================================================================================================================
     !  ! modules
     !  use, intrinsic :: iso_c_binding  
     !  ! implicit variable handling
     !  implicit none
     !  !---------------------------------------------------------------------------------------------------------------------------------
     !  ! input variables
     !  character(kind=c_char),dimension(*) :: filename
     !  type(c_ptr),value                   :: p4est
     !  !---------------------------------------------------------------------------------------------------------------------------------
     !  ! output variables
     !  !=================================================================================================================================
     !end subroutine p4_savemesh
  end interface
#endif
end module p4est_bindings_names 
