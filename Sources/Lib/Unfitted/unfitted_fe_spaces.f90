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
  use block_layout_names
  
  use triangulation_names
  use reference_fe_names
  use environment_names
  use fe_space_names
  use conditions_names
  use assembler_names
  
  use unfitted_triangulations_names
  use piecewise_cell_map_names

  use ParameterList

  implicit none
# include "debug.i90"
  private

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

  type, extends(unfitted_fe_cell_iterator_t) :: unfitted_hp_adaptive_fe_cell_iterator_t
      type (hp_adaptive_fe_cell_iterator_t) :: adaptive_fe
    contains

      procedure :: create               => uhpafeci_create
      procedure :: free                 => uhpafeci_free
      procedure :: next                 => uhpafeci_next
      procedure :: first                => uhpafeci_first
      procedure :: set_gid              => uhpafeci_set_gid

      procedure                  :: is_strong_dirichlet_dof    => uhpafeci_is_strong_dirichlet_dof
      procedure                  :: is_hanging_dof             => uhpafeci_is_hanging_dof
      !procedure, non_overridable :: has_hanging_dofs           => uhpafeci_has_hanging_dofs
      procedure                  :: determine_has_hanging_dofs => uhpafeci_set_has_hanging_dofs
      !procedure, non_overridable, private :: apply_constraints => uhpafeci_apply_constraints
      procedure, private :: assembly_array =>  uhpafeci_assembly_array
      procedure, private :: assembly_matrix => uhpafeci_assembly_matrix
      procedure, private :: assembly_matrix_array => uhpafeci_assembly_matrix_array
      procedure, private :: assembly_matrix_array_with_strong_bcs => uhpafeci_assembly_matrix_array_with_strong_bcs

  end type unfitted_hp_adaptive_fe_cell_iterator_t


  type :: unfitted_integration_manager_t

    private

    class(serial_fe_space_t), pointer :: fe_space       => NULL()
    class(marching_cubes_t),  pointer :: marching_cubes => NULL()

    ! All the machinery for integrating in subcells
    type(quadrature_t)                     :: quadrature_subelem
    type(tet_lagrangian_reference_fe_t)    :: geo_reference_subelem
    type(cell_map_t)                         :: cell_map_subelem
    type(quadrature_t),        allocatable :: cut_quadratures(:)
    type(cell_map_t),            allocatable :: cut_cell_maps(:)
    type(cell_integrator_t),   allocatable :: cut_cell_integrators(:,:)

    ! All the machinery for integrating in subfacets
    type(quadrature_t)                     :: quadrature_subfacet
    type(quadrature_t),        allocatable :: cut_boundary_quadratures_cell_dim(:)
    type(piecewise_cell_map_t),  allocatable :: cut_boundary_piecewise_cell_maps(:)
    type(cell_map_t),            allocatable :: cut_boundary_cell_maps(:)
    type(cell_integrator_t),   allocatable :: cut_boundary_cell_integrators(:,:)    

    ! Auxiliary dummy empty quadratures
    type(quadrature_t)             :: empty_quadrature
    type(cell_map_t)                 :: empty_cell_map
    type(piecewise_cell_map_t)       :: empty_piecewise_cell_map
    type(cell_integrator_t), allocatable  :: empty_cell_integrator(:)
    
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

  end type unfitted_integration_manager_t

  type, extends(serial_fe_space_t) :: serial_unfitted_fe_space_t
    private

      class(serial_unfitted_triangulation_t), pointer :: unfitted_triangulation =>  NULL()
      type(unfitted_integration_manager_t) :: unfitted_integration

    contains

      ! Creation / deletion methods
      procedure,  private :: serial_fe_space_create_same_reference_fes_on_all_cells => sufs_create_same_reference_fes_on_all_cells
      procedure,  private :: serial_fe_space_create_different_ref_fes_between_cells         => sufs_space_create_different_ref_fes_between_cells
      procedure           :: free  => sufs_free

      ! Creation of the iterator
      procedure :: create_fe_cell_iterator           => sufs_create_fe_cell_iterator

  end type serial_unfitted_fe_space_t

  type, extends(serial_hp_adaptive_fe_space_t) :: serial_unfitted_hp_adaptive_fe_space_t
    private
      class(unfitted_p4est_serial_triangulation_t), pointer :: unfitted_triangulation =>  NULL()
      type(unfitted_integration_manager_t) :: unfitted_integration
      integer(ip), allocatable :: aggregate_ids(:)
      real(rp) :: max_separation_from_root = -1.0_rp
      logical :: use_constraints = .true.
      type(ParameterList_t), public :: debug_info
    contains
      ! Creation / deletion methods
      procedure           :: serial_fe_space_create_same_reference_fes_on_all_cells => suhpafs_create_same_reference_fes_on_all_cells
      procedure           :: serial_fe_space_create_different_ref_fes_between_cells         => suhpafs_space_create_different_ref_fes_between_cells
      procedure           :: free                                                   => suhpafs_free
      procedure           :: set_use_constraints                                    => suhpafs_set_use_constraints
      
      ! Creation of the iterator
      procedure :: create_fe_cell_iterator                                               => suhpafs_create_fe_cell_iterator
      
      ! Creation of constrained degrees of freedom
      procedure          :: generate_global_dof_numbering                                           => suhpafs_generate_global_dof_numbering 
      procedure          :: fill_fe_dofs_and_count_dofs                            => suhpafs_procedure_fill_fe_dofs_and_count_dofs

      ! Getters
      procedure, non_overridable :: get_aggregate_ids                               => suhpafs_get_aggregate_ids
      procedure, non_overridable :: get_max_separation_from_root                    => suhpafs_get_max_separation_from_root

      ! Private TBPs
      procedure, private, non_overridable :: allocate_and_fill_aggregate_ids        => suhpafs_allocate_and_fill_aggregate_ids
      procedure, private, non_overridable :: compute_aggregate_size                 => suhpafs_compute_aggregate_size
      procedure, private, non_overridable :: setup_cut_cells_constraints            => suhpafs_setup_cut_cells_constraints
      

  end type serial_unfitted_hp_adaptive_fe_space_t

  type, extends(par_fe_space_t) :: par_unfitted_fe_space_t
    private

      class(par_unfitted_triangulation_t), pointer :: unfitted_triangulation =>  NULL()
      type(unfitted_integration_manager_t) :: unfitted_integration

    contains

      ! Creation / deletion methods
      procedure,  private :: par_fe_space_create_same_reference_fes_on_all_cells => pufs_create_same_reference_fes_on_all_cells
      procedure,  private :: par_fe_space_create_different_ref_fes_between_cells         => pufs_create_different_ref_fes_between_cells
      procedure           :: free  => pufs_free

      ! Creation of the iterator
      procedure :: create_fe_cell_iterator           => pufs_create_fe_cell_iterator

  end type par_unfitted_fe_space_t


  public :: unfitted_fe_cell_iterator_t
  public :: unfitted_hp_adaptive_fe_cell_iterator_t
  public :: serial_unfitted_fe_space_t
  public :: serial_unfitted_hp_adaptive_fe_space_t
  public :: par_unfitted_fe_space_t

contains

#include "../Unfitted/sbm_unfitted_fe_cell_iterator.i90"
#include "../Unfitted/sbm_unfitted_hp_adaptive_fe_cell_iterator.i90"
#include "../Unfitted/sbm_unfitted_integration_manager.i90"
#include "../Unfitted/sbm_serial_unfitted_fe_space.i90"
#include "../Unfitted/sbm_serial_unfitted_hp_adaptive_fe_space.i90"
#include "../Unfitted/sbm_par_unfitted_fe_space.i90"

end module unfitted_fe_spaces_names
