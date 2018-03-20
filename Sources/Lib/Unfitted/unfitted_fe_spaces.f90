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

module unfitted_fe_spaces_names
  use types_names
  use memor_names
  use field_names
  use list_types_names
  use block_layout_names
  use std_vector_names
  use hash_table_names
  
  use triangulation_names
  use p4est_triangulation_names
  use reference_fe_names
  use environment_names
  use fe_space_names
  use conditions_names
  use assembler_names
  
  use unfitted_triangulations_names
  use piecewise_cell_map_names
  use level_set_functions_gallery_names

  use ParameterList

  implicit none
# include "debug.i90"
  private
  
  
  integer(ip), parameter :: pos_map_in_domain = 1
  integer(ip), parameter :: pos_map_on_boundary = 2
  integer(ip), parameter :: pos_map_max_id = 2
  
 ! Types from unfitted branch that are *** UNDER QUARANTINE ***
 type, extends(fe_cell_iterator_t) :: unfitted_fe_cell_iterator_t
    private
    class(unfitted_integration_manager_t), pointer :: unfitted_integration_manager => NULL()
  contains

    ! Creation / deletion methods
    procedure :: create => unfitted_fe_cell_iterator_create
    procedure :: free   => unfitted_fe_cell_iterator_free

    ! Getters that override
    procedure          :: get_quadrature        => unfitted_fe_cell_iterator_get_quadrature
    
    ! Getters that extend
    procedure          :: get_boundary_quadrature          => unfitted_fe_cell_iterator_get_boundary_quadrature
    procedure          :: get_boundary_piecewise_cell_map    => unfitted_fe_cell_iterator_get_boundary_piecewise_cell_map
    procedure          :: get_boundary_cell_map              => unfitted_fe_cell_iterator_get_boundary_cell_map
    procedure          :: get_boundary_cell_integrator     => unfitted_fe_cell_iterator_get_boundary_cell_integrator

    ! Updater that overrides
    procedure :: update_integration      => unfitted_fe_cell_iterator_update_integration
    procedure :: update_cell_map         => unfitted_fe_cell_iterator_update_cell_map
    procedure :: update_cell_integrators => unfitted_fe_cell_iterator_update_cell_integrators
    
    ! Updater that extends
    procedure :: update_boundary_integration  => unfitted_fe_cell_iterator_update_boundary_integration

    ! Private TBPs
    procedure, non_overridable, private :: update_cut_quadratures => unfitted_fe_cell_iterator_update_cut_quadratures
    procedure, non_overridable, private :: update_cut_cell_maps     => unfitted_fe_cell_iterator_update_cut_cell_maps
    procedure, non_overridable, private :: update_cut_cell_integrators  => unfitted_fe_cell_iterator_update_cut_cell_integrators
    procedure, non_overridable, private :: update_cut_boundary_quadratures => unfitted_fe_cell_iterator_update_cut_boundary_quadratures
    procedure, non_overridable, private :: update_cut_boundary_cell_maps     => unfitted_fe_cell_iterator_update_cut_boundary_cell_maps
    procedure, non_overridable, private :: update_cut_boundary_cell_integrators  => &
    unfitted_fe_cell_iterator_update_cut_boundary_cell_integrators

  end type unfitted_fe_cell_iterator_t

 type, extends(fe_facet_iterator_t) :: unfitted_fe_facet_iterator_t
    private
    class(unfitted_integration_manager_t), pointer :: unfitted_integration_manager => NULL()
    type(facet_maps_t), pointer :: unfitted_facet_maps => NULL()
    type(p_facet_integrator_t), allocatable  :: unfitted_facet_integrators(:)
  contains

    procedure :: create                               => unfitted_fe_facet_iterator_create
    procedure :: free                                 => unfitted_fe_facet_iterator_free
    procedure :: update_integration                   => unfitted_fe_facet_iterator_update_integration
    procedure, private :: update_quadrature           => unfitted_fe_facet_iterator_update_quadrature
    procedure :: get_quadrature                       => unfitted_fe_facet_iterator_get_quadrature
    procedure, non_overridable :: update_facet_maps_interpolation  => unfitted_fe_facet_iterator_update_facet_maps_interpolation
    procedure, non_overridable :: update_facet_integrators_interpolation => uffi_update_facet_integrators_interpolation
    procedure :: get_quadrature_points_coordinates    => unfitted_fe_facet_iterator_get_quadrature_points_coordinates
    procedure :: get_normals                          => unfitted_fe_facet_iterator_get_normals
    procedure :: get_det_jacobian                     => unfitted_fe_facet_iterator_get_det_jacobian
    procedure :: compute_characteristic_length        => unfitted_fe_facet_iterator_compute_characteristic_length

    procedure :: get_values_scalar     => unfitted_fe_facet_iterator_get_values_scalar
    procedure :: get_values_vector     => unfitted_fe_facet_iterator_get_values_vector
    procedure :: get_gradients_scalar  => unfitted_fe_facet_iterator_get_gradients_scalar
    procedure :: get_curls             => unfitted_fe_facet_iterator_get_curls_vector 

  end type unfitted_fe_facet_iterator_t

  type :: unfitted_integration_manager_t

    private

    class(serial_fe_space_t), pointer :: fe_space       => NULL()

    ! All the machinery for integrating in subcells
    type(quadrature_t)                     :: quadrature_subelem
    type(tet_lagrangian_reference_fe_t)    :: geo_reference_subelem
    type(cell_map_t)                         :: cell_map_subelem
    type(quadrature_t),        allocatable :: cut_quadratures(:)
    type(cell_map_t),            allocatable :: cut_cell_maps(:)
    type(cell_integrator_t),   allocatable :: cut_cell_integrators(:,:)
    type(hash_table_ip_ip_t) :: num_sub_cells_to_pos

    ! All the machinery for integrating in subfacets
    type(quadrature_t)                     :: quadrature_subfacet
    type(quadrature_t),        allocatable :: cut_boundary_quadratures_cell_dim(:)
    type(piecewise_cell_map_t),  allocatable :: cut_boundary_piecewise_cell_maps(:)
    type(cell_map_t),            allocatable :: cut_boundary_cell_maps(:)
    type(cell_integrator_t),   allocatable :: cut_boundary_cell_integrators(:,:)    
    type(hash_table_ip_ip_t) :: num_unfitted_sub_facets_to_pos

    ! All the machinery to integrate in fitted subfacets
    type(tet_lagrangian_reference_fe_t)   :: geo_reference_subfacet
    type(quadrature_t)                     :: quadrature_fitted_subfacet
    type(cell_map_t)                     :: cell_map_subfacet
    type(quadrature_t),       allocatable :: cut_fitted_facet_quadratures(:)
    type(facet_maps_t),       allocatable :: cut_fitted_facet_maps(:,:)
    type(facet_integrator_t), allocatable :: cut_fitted_facet_integrators(:,:,:)
    type(hash_table_ip_ip_t) :: num_fitted_sub_facets_to_pos

    ! Auxiliary dummy empty quadratures
    type(quadrature_t)             :: empty_quadrature
    type(cell_map_t)                 :: empty_cell_map
    type(piecewise_cell_map_t)       :: empty_piecewise_cell_map
    type(cell_integrator_t), allocatable  :: empty_cell_integrator(:)
    type(quadrature_t)             :: empty_facet_quadrature
    type(facet_maps_t)             :: empty_facet_maps(pos_map_max_id)
    type(facet_integrator_t), allocatable :: empty_facet_integrators(:,:)
    
    contains

      ! Creation / Deletion methods
      procedure, non_overridable :: create => uim_create
      procedure, non_overridable :: free   => uim_free
      
      !Private TBPs
      procedure, non_overridable, private :: check_assumptions      => uim_check_assumptions
      procedure, non_overridable, private :: init_reference_subelem => uim_init_reference_subelem
      procedure, non_overridable, private :: free_reference_subelem => uim_free_reference_subelem
      procedure, non_overridable, private :: init_cut_integration   => uim_init_cut_integration
      procedure, non_overridable, private :: free_cut_integration   => uim_free_cut_integration
      procedure, non_overridable, private :: init_cut_boundary_integration   => uim_init_cut_boundary_integration
      procedure, non_overridable, private :: free_cut_boundary_integration   => uim_free_cut_boundary_integration
      procedure, non_overridable, private :: init_cut_fitted_facets_integration   => uim_init_cut_fitted_facets_integration
      procedure, non_overridable, private :: free_cut_fitted_facets_integration   => uim_free_cut_fitted_facets_integration

  end type unfitted_integration_manager_t

  type, extends(serial_fe_space_t) :: serial_unfitted_fe_space_t
    private
      type(unfitted_integration_manager_t) :: unfitted_integration
      integer(ip), allocatable :: aggregate_ids(:)
      logical    , allocatable :: is_in_aggregate_x_cell(:)
      real(rp)   , allocatable :: aggregate_size(:)
      logical :: use_constraints = .false.
      type(ParameterList_t), public :: debug_info
      integer(ip) :: num_hanging_dofs_full_cells
      integer(ip) :: num_hanging_dofs_other
      integer(ip) :: num_aggregated_dofs
    contains
      ! Creation / deletion methods
      procedure           :: serial_fe_space_create_same_reference_fes_on_all_cells => sufs_create_same_reference_fes_on_all_cells
      procedure           :: serial_fe_space_create_different_ref_fes_between_cells => sufs_space_create_different_ref_fes_between_cells
      procedure           :: free                                                   => sufs_free
      procedure           :: set_use_constraints                                    => sufs_set_use_constraints
      
      ! Creation of the iterator (overrides)
      procedure :: create_fe_cell_iterator       => sufs_create_fe_cell_iterator
      procedure :: create_fe_facet_iterator      => sufs_create_fe_facet_iterator

      ! Generation of Dofs (overrides)
      procedure :: count_dofs => sufs_count_dofs
      procedure :: list_dofs  => sufs_list_dofs
      procedure :: setup_hanging_node_constraints => sufs_setup_hanging_node_constraints
      
      ! Setup integration structures (overrides)
      procedure    :: set_up_cell_integration   => sufs_set_up_cell_integration
      procedure    :: set_up_facet_integration  => sufs_set_up_facet_integration

      ! Mesh refinement
      procedure :: refine_mesh_for_small_aggregates => sufs_refine_mesh_for_small_aggregates

      ! Getters
      procedure, non_overridable :: get_aggregate_ids          => sufs_get_aggregate_ids
      procedure, non_overridable :: get_aggregate_size         => sufs_get_aggregate_size
      procedure, non_overridable :: get_is_in_aggregate_x_cell => sufs_get_is_in_aggregate_x_cell

      ! Printers
      procedure, non_overridable :: print_debug_info => sufs_print_debug_info

      ! Private TBPs
      procedure, private, non_overridable :: allocate_and_fill_aggregate_ids         => sufs_allocate_and_fill_aggregate_ids
      procedure, private, non_overridable :: check_for_full_neighbors                => sufs_check_for_full_neighbors
      procedure, private, non_overridable :: fill_proper_vef_constrains_full_cell    => sufs_fill_proper_vef_constrains_full_cell
      procedure, private, non_overridable :: setup_only_hanging_node_constraints     => sufs_setup_only_hanging_node_constraints
      procedure, private, non_overridable :: setup_only_cell_aggregation_constraints => sufs_setup_only_cell_aggregation_constraints
      procedure, private, non_overridable :: compute_aggregate_sizes                 => sufs_compute_aggregate_sizes

  end type serial_unfitted_fe_space_t

  type, extends(par_fe_space_t) :: par_unfitted_fe_space_t
    private

      class(par_unfitted_triangulation_t), pointer :: unfitted_triangulation =>  NULL()
      type(unfitted_integration_manager_t) :: unfitted_integration

    contains

      ! Creation / deletion methods
      procedure           :: par_fe_space_create_same_reference_fes_on_all_cells => pufs_create_same_reference_fes_on_all_cells
      procedure           :: par_fe_space_create_different_ref_fes_between_cells         => pufs_create_different_ref_fes_between_cells
      procedure           :: free  => pufs_free

      ! Creation of the iterator
      procedure :: create_fe_cell_iterator           => pufs_create_fe_cell_iterator
      !procedure :: create_fe_facet_iterator          => pufs_create_fe_facet_iterator

      ! Setup integration structures (overrides)
      procedure    :: set_up_cell_integration   => pufs_set_up_cell_integration
      procedure    :: set_up_facet_integration  => pufs_set_up_facet_integration

  end type par_unfitted_fe_space_t


  public :: unfitted_fe_cell_iterator_t
  public :: serial_unfitted_fe_space_t
  public :: par_unfitted_fe_space_t

contains

#include "../Unfitted/sbm_unfitted_fe_cell_iterator.i90"
#include "../Unfitted/sbm_unfitted_fe_facet_iterator.i90"
#include "../Unfitted/sbm_unfitted_integration_manager.i90"
#include "../Unfitted/sbm_serial_unfitted_fe_space.i90"
#include "../Unfitted/sbm_par_unfitted_fe_space.i90"

end module unfitted_fe_spaces_names
