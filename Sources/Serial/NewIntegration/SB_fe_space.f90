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
module SB_fe_space_names
  use memor_names
  use allocatable_array_names
  use matrix_names
  use array_names
  use sparse_matrix_names
  use serial_scalar_array_names
  use block_sparse_matrix_names
  use serial_block_array_names
  use SB_matrix_array_assembler_names
  use SB_sparse_matrix_array_assembler_names
  use SB_block_sparse_matrix_array_assembler_names
  use types_names
  use list_types_names
  use reference_fe_names
  use triangulation_names
  use reference_fe_factory_names
  use integration_tools_names
  use migratory_element_names
  use conditions_names
  use hash_table_names
  use graph_names
  use sort_names
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
  !   in the FE space, the total number of DOFs and an array that provides the DOFs
  !   in a VEF.

  type :: SB_finite_element_t
     private 
     integer(ip)                                 :: number_nodes
     integer(ip)                                 :: number_fe_spaces
     
     type(elem_topology_t)         , pointer     :: cell 
     class(reference_fe_t)         , pointer     :: reference_fe_geo 
     type(fe_map_t)                , pointer     :: fe_map
     
     type(p_reference_fe_t)        , pointer     :: reference_fe_phy(:) 
     type(SB_quadrature_t)         , pointer     :: quadrature
     type(SB_p_volume_integrator_t), pointer     :: volume_integrator(:)
     
     type(i1p_t)                   , allocatable :: elem2dof(:)
     type(i1p_t)                   , allocatable :: bc_code(:)
     type(r1p_t)                   , allocatable :: bc_value(:)
   contains
     
     ! procedure :: create Pending
     procedure, non_overridable :: update_integration
     ! procedure :: free Pending
     ! procedure :: print Pending
     procedure, non_overridable :: get_number_nodes 
     procedure, non_overridable :: get_fe_map
     procedure, non_overridable :: get_quadrature	
     procedure, non_overridable :: get_volume_integrator
     procedure, non_overridable :: get_elem2dof 
     procedure, non_overridable :: get_bc_code 
     procedure, non_overridable :: get_bc_value 
     procedure, non_overridable :: get_number_nodes_per_field 
					procedure, non_overridable :: get_subset_id 
  end type SB_finite_element_t

  type :: p_SB_finite_element_t
     type(SB_finite_element_t), pointer :: p
  end type p_SB_finite_element_t

  public :: SB_finite_element_t, p_SB_finite_element_t

  type :: finite_face_t
     private
     integer(ip)                            :: number_fe_spaces
     type(face_topology_t)        , pointer :: face_topology
     type(p_SB_finite_element_t)            :: neighbour_fe(2)
     type(p_face_integrator_t), allocatable :: face_integrator(:)
   contains
     procedure, non_overridable :: create                  => finite_face_create
     procedure, non_overridable :: update_face_integration => finite_face_update_face_integration
     procedure                  :: free                    => finite_face_free
  end type finite_face_t

  public :: finite_face_t

  type :: SB_serial_fe_space_t
     private
     integer(ip)                                 :: number_fe_spaces   
     type(p_reference_fe_t)        , allocatable :: reference_fe_geo_list(:)
     type(p_fe_map_t)              , allocatable :: fe_map(:)
     type(p_reference_fe_t)        , allocatable :: reference_fe_phy_list(:)
     type(SB_p_quadrature_t)       , allocatable :: quadrature(:)
     type(SB_p_volume_integrator_t), allocatable :: volume_integrator(:)
     type(p_face_integrator_t)     , allocatable :: face_integrator(:)
     
     type(triangulation_t)         , pointer     :: triangulation
     type(SB_finite_element_t)     , allocatable :: fe_array(:)
     type(finite_face_t)           , allocatable :: face_array(:)
     
     ! Data related to block structure of the FE system + size of each block
     integer(ip)                                 :: number_blocks
     integer(ip)                   , allocatable :: field_blocks(:)
     logical                       , allocatable :: field_coupling(:,:)
     integer(ip)                   , allocatable :: number_dofs(:)
     
     ! Acceleration arrays
     type(list_2d_t)               , allocatable :: vef2dof(:)
   contains
     procedure, non_overridable :: create
     procedure, non_overridable :: fill_dof_info
     procedure, non_overridable :: free
     procedure, non_overridable :: print
     procedure, non_overridable :: initialize_integration
     procedure, non_overridable, private :: initialize_quadrature
     procedure, non_overridable, private :: initialize_volume_integrator
     procedure, non_overridable, private :: initialize_fe_map
     procedure, non_overridable :: create_assembler
     procedure, non_overridable :: symbolic_setup_assembler
     procedure, non_overridable :: get_number_elements
     procedure, non_overridable :: get_number_interior_faces
     procedure, non_overridable :: get_number_fe_spaces
     procedure, non_overridable :: get_finite_element
     procedure, non_overridable :: get_finite_face
     procedure, non_overridable :: get_number_blocks
     procedure, non_overridable :: get_field_blocks
     procedure, non_overridable :: get_field_coupling
     procedure, non_overridable :: get_max_number_nodes
     procedure, non_overridable :: create_face_array
  end type SB_serial_fe_space_t

  public :: SB_serial_fe_space_t

contains

  ! Includes with all the TBP and supporting subroutines for the types above.
  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_finite_element.i90"

#include "sbm_finite_face.i90"

#include "sbm_serial_fe_space.i90"

#include "sbm_serial_fe_space_faces.i90"

end module SB_fe_space_names
