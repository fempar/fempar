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
  use serial_scalar_matrix_names
  use serial_scalar_array_names
  use serial_block_matrix_names
  use serial_block_array_names
  use SB_matrix_array_assembler_names
  use SB_serial_scalar_matrix_array_assembler_names
  use SB_serial_block_matrix_array_assembler_names
  use types_names
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
     type(elem_topology_t), pointer :: cell 
     class(reference_fe_t), pointer :: geometry_reference_fe 
     type(p_reference_fe_t), pointer :: reference_fe(:) 
     type(SB_quadrature_t), pointer :: quadrature
     type(fe_map_t), pointer :: fe_map
     type(SB_p_volume_integrator_t), pointer :: volume_integrator(:)
     integer(ip) :: number_nodes
     type(i1p_t), pointer :: elem2dof(:)
     type(i1p_t), pointer :: bc_code(:)
     type(r1p_t), pointer :: bc_value(:)
   contains
     procedure :: get_volume_integrator
     procedure :: get_quadrature	
     procedure :: get_fe_map
     procedure :: update_integration

     procedure :: get_elem2dof 
     procedure :: get_bc_code 
     procedure :: get_bc_value 
     procedure :: get_number_nodes 
     procedure :: get_number_nodes_field 
  end type SB_finite_element_t

  public :: SB_finite_element_t

  type :: SB_serial_fe_space_t
     private
     type(triangulation_t), pointer :: triangulation
     type(p_reference_fe_t), allocatable :: reference_fe_geo_list(:)
     type(p_reference_fe_t), allocatable :: reference_fe_phy_list(:)
     type(SB_p_quadrature_t), allocatable :: quadrature(:)
     type(p_fe_map_t), allocatable :: fe_map(:)
     type(SB_p_volume_integrator_t), allocatable :: volume_integrator(:)
     type(SB_finite_element_t), allocatable :: fe_array(:)
     integer(ip) :: number_fe_spaces   
     integer(ip), allocatable :: field_blocks(:)
     integer(ip) :: number_blocks
     logical, allocatable :: fields_coupling(:,:)
     integer(ip), allocatable :: number_dofs(:)
     ! Acceleration arrays
     type(list_2d_t), allocatable       :: vef2dof(:)
   contains
     procedure :: get_number_elements
     procedure :: get_fe
     procedure :: initialize_integration
     procedure :: initialize_volume_integrator
     procedure :: initialize_quadrature
     procedure :: initialize_fe_map
     procedure :: create_assembler
     procedure :: symbolic_setup_assembler
     procedure :: get_blocks
     procedure :: get_max_number_nodes
     procedure :: get_fields_coupling
     procedure :: fill_dof_info
     procedure :: get_number_blocks
     procedure :: get_number_fe_spaces
     procedure :: create
     procedure :: print
  end type SB_serial_fe_space_t

  public :: SB_serial_fe_space_t

contains

  ! Includes with all the TBP and supporting subroutines for the types above.
  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_finite_element.i90"

#include "sbm_serial_fe_space.i90"

end module SB_fe_space_names
