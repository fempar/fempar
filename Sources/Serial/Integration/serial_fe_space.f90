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
  use allocatable_array_names
  use hash_table_names
  
  use triangulation_names
  use conditions_names
  
  use reference_fe_names
  use reference_fe_factory_names
  use field_names
  
  use matrix_names
  use vector_names
  use array_names
  use sparse_matrix_names
  use block_sparse_matrix_names
  
  use serial_scalar_array_names
  use serial_block_array_names
  use matrix_array_assembler_names
  use sparse_matrix_array_assembler_names
  use block_sparse_matrix_array_assembler_names
  
 ! Parallel modules
  use par_triangulation_names
  use par_conditions_names
  
  
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
     integer(ip), pointer                        :: field_blocks(:)
     
     type(elem_topology_t)         , pointer     :: cell 
     type(fe_map_t)                , pointer     :: fe_map
     
     type(p_reference_fe_t)        , pointer     :: reference_fe_phy(:) 
     type(quadrature_t)         , pointer     :: quadrature
     type(p_volume_integrator_t), pointer     :: volume_integrator(:)
     
     type(i1p_t)                   , allocatable :: elem2dof(:)
     type(i1p_t)                   , allocatable :: bc_code(:)
     type(r1p_t)                   , allocatable :: bc_value(:)
   contains
     
     procedure, non_overridable, private :: create =>  finite_element_create
     procedure, non_overridable          :: update_integration => finite_element_update_integration
     procedure, non_overridable, private :: free => finite_element_free
     ! procedure :: print Pending
     procedure, non_overridable, private :: fill_interior_dofs => finite_element_fill_interior_dofs
     procedure, non_overridable, private :: fill_interior_dofs_on_vef => finite_element_fill_interior_dofs_on_vef
     procedure, non_overridable, private :: fill_interior_dofs_on_vef_from_source_element => finite_element_fill_interior_dofs_on_vef_from_source_element
     procedure, non_overridable, private :: fill_dofs_on_vef => finite_element_fill_dofs_on_vef
     
     procedure, non_overridable :: get_number_nodes => finite_element_get_number_nodes
     procedure, non_overridable :: get_fe_map => finite_element_get_fe_map
     procedure, non_overridable :: get_quadrature => finite_element_get_quadrature
     procedure, non_overridable :: get_volume_integrator => finite_element_get_volume_integrator
     procedure, non_overridable :: get_elem2dof  => finite_element_get_elem2dof
     procedure, non_overridable :: get_bc_code  => finite_element_get_bc_code
     procedure, non_overridable :: get_bc_value => finite_element_get_bc_value
     procedure, non_overridable :: get_number_nodes_per_field => finite_element_get_number_nodes_per_field
     procedure, non_overridable :: get_subset_id => finite_element_get_subset_id
     
     procedure, non_overridable, private :: update_scalar_values => finite_element_update_scalar_values
     procedure, non_overridable, private :: update_vector_values => finite_element_update_vector_values
     procedure, non_overridable, private :: update_tensor_values => finite_element_update_tensor_values
     generic :: update_values => update_scalar_values, &
                               & update_vector_values, &
                               & update_tensor_values
     
  end type finite_element_t

  type :: p_finite_element_t
     type(finite_element_t), pointer :: p
  end type p_finite_element_t

  public :: finite_element_t, p_finite_element_t

  type :: finite_face_t
     private
     integer(ip)                            :: number_fe_spaces
     type(face_topology_t)    , pointer     :: face_topology
     type(p_finite_element_t)            :: neighbour_fe(2)
     type(face_map_t)         , pointer     :: map
     type(quadrature_t)    , pointer     :: quadrature
     type(p_face_integrator_t), allocatable :: face_integrator(:)
   contains
     procedure, non_overridable :: create              => finite_face_create
     procedure, non_overridable :: update_integration  => finite_face_update_integration
     procedure, non_overridable :: free                => finite_face_free
     procedure, non_overridable :: is_boundary         => finite_face_is_boundary
     procedure, non_overridable :: number_neighbours   => finite_face_number_neighbours
     procedure, non_overridable :: get_bc_code         => finite_face_get_bc_code
     procedure, non_overridable :: get_bc_value        => finite_face_get_bc_value
     procedure, non_overridable :: get_elem2dof        => finite_face_get_elem2dof
     procedure, non_overridable :: get_map             => finite_face_get_map
     procedure, non_overridable :: get_quadrature      => finite_face_get_quadrature
     procedure, non_overridable :: get_face_integrator => finite_face_get_face_integrator
  end type finite_face_t

  public :: finite_face_t

  type :: serial_fe_space_t
     private
     integer(ip)                                 :: number_fe_spaces   
     type(p_fe_map_t)              , allocatable :: fe_map(:)
     type(face_map_t)              , allocatable :: face_map(:)
     type(p_reference_fe_t)        , allocatable :: reference_fe_phy_list(:)
     type(p_quadrature_t)       , allocatable :: quadrature(:)
     type(quadrature_t)         , allocatable :: face_quadrature(:)
     type(p_volume_integrator_t), allocatable :: volume_integrator(:)
     type(p_face_integrator_t)     , allocatable :: face_integrator(:)
     
     type(triangulation_t)         , pointer     :: triangulation
     type(finite_element_t)        , allocatable :: fe_array(:)
     type(finite_face_t)           , allocatable :: face_array(:)
     
     ! Data related to block structure of the FE system + size of each block
     integer(ip)                                 :: number_blocks
     integer(ip)                   , allocatable :: field_blocks(:)
     logical                       , allocatable :: field_coupling(:,:)
     integer(ip)                   , allocatable :: number_dofs_per_block(:)
     integer(ip)                   , allocatable :: number_dofs_per_field(:)
   contains
     procedure, non_overridable :: create => serial_fe_space_create
     procedure, non_overridable :: fill_dof_info => serial_fe_space_fill_dof_info
     procedure, non_overridable, private :: fill_elem2dof_and_count_dofs => serial_fe_space_fill_elem2dof_and_count_dofs
     procedure, non_overridable :: free => serial_fe_space_free
     procedure, non_overridable :: print => serial_fe_space_print
     procedure, non_overridable :: initialize_integration => serial_fe_space_initialize_integration
     procedure, non_overridable, private :: initialize_quadrature => serial_fe_space_initialize_quadrature
     procedure, non_overridable, private :: initialize_volume_integrator => serial_fe_space_initialize_volume_integrator
     procedure, non_overridable, private :: initialize_fe_map => serial_fe_space_initialize_fe_map
     procedure, non_overridable :: create_assembler => serial_fe_space_create_assembler
     procedure, non_overridable :: symbolic_setup_assembler => serial_fe_space_symbolic_setup_assembler
     procedure, non_overridable :: get_number_elements => serial_fe_space_get_number_elements
     procedure, non_overridable :: get_number_interior_faces => serial_fe_space_get_number_interior_faces
     procedure, non_overridable :: get_number_boundary_faces => serial_fe_space_get_number_boundary_faces
     procedure, non_overridable :: get_number_fe_spaces => serial_fe_space_get_number_fe_spaces
     procedure, non_overridable :: get_finite_element => serial_fe_space_get_finite_element
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
     procedure, private         :: create_fe_function_scalar => serial_fe_space_create_fe_function_scalar
     procedure, private         :: create_fe_function_vector => serial_fe_space_create_fe_function_vector
     procedure, private         :: create_fe_function_tensor => serial_fe_space_create_fe_function_tensor
     generic :: create_fe_function => create_fe_function_scalar, &
                                    & create_fe_function_vector, &
                                    & create_fe_function_tensor
     
  end type serial_fe_space_t

  public :: serial_fe_space_t
  
  type, extends(serial_fe_space_t) :: par_fe_space_t
     private
     type(par_triangulation_t), pointer     :: par_triangulation
     type(finite_element_t)   , allocatable :: ghost_fe_array(:)
  contains
     procedure :: par_fe_space_create
     procedure :: par_fe_space_fill_dof_info
     procedure :: par_fe_space_fill_elem2dof_and_count_dofs
     procedure :: par_fe_space_print
     procedure :: par_fe_space_free
  end type

  public :: par_fe_space_t
  
  type fe_function_scalar_t
   private
   integer(ip) :: fe_space_id
   
   integer(ip) :: current_number_nodes             ! Not being used
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
   
   integer(ip) :: current_number_nodes             ! Not being used
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
   
   integer(ip) :: current_number_nodes             ! Not being used
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
  
  public :: fe_function_scalar_t, fe_function_vector_t, fe_function_tensor_t
  
contains

  ! Includes with all the TBP and supporting subroutines for the types above.
  ! In a future, we would like to use the submodule features of FORTRAN 2008.

#include "sbm_finite_element.i90"

#include "sbm_finite_face.i90"

#include "sbm_serial_fe_space.i90"

#include "../../Par/Integration/sbm_par_fe_space.i90"

#include "sbm_serial_fe_space_faces.i90"

#include "sbm_serial_fe_space_fe_function.i90"

#include "sbm_fe_function.i90"

end module serial_fe_space_names
