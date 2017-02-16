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

module serial_unfitted_triangulation_names
  use fempar_names
  use level_set_functions_gallery_names
  
  implicit none
# include "debug.i90"
  private

  ! Include the look-up tables 
# include "mc_tables_qua4.i90"
!# include "mc_tables_hex8.i90"

  type, extends(cell_accessor_t) :: unfitted_cell_accessor_t
    private
    class(serial_unfitted_triangulation_t), pointer :: serial_unfitted_triangulation
  contains
    procedure :: cell_accessor_create => unfitted_cell_accessor_cell_accessor_create
    procedure :: cell_accessor_free   => unfitted_cell_accessor_cell_accessor_free
    !procedure, non_overridable, private :: get_marching_cube_case => unfitted_cell_accessor_get_marching_cube_case
  end type unfitted_cell_accessor_t

  type :: unfitted_cell_iterator_t
    private
    type(unfitted_cell_accessor_t) :: current_cell_accessor
  contains
    procedure, non_overridable, private :: create       => unfitted_cell_iterator_create
    procedure, non_overridable          :: free         => unfitted_cell_iterator_free
    procedure, non_overridable          :: init         => unfitted_cell_iterator_init
    procedure, non_overridable          :: next         => unfitted_cell_iterator_next
    procedure, non_overridable          :: has_finished => unfitted_cell_iterator_has_finished
    procedure, non_overridable          :: current      => unfitted_cell_iterator_current
  end type unfitted_cell_iterator_t

  type, extends(serial_triangulation_t) :: serial_unfitted_triangulation_t
    private

    ! Look up-tables (precomputed off-line, for each cell type)
    integer(ip)                :: mc_table_num_cases
    integer(ip)                :: mc_table_max_num_subcells
    integer(ip)                :: mc_table_num_nodes_subcell
    integer(ip),   allocatable :: mc_table_num_subcells_per_case(:)
    integer(ip),   allocatable :: mc_table_num_cut_edges_per_case(:)
    integer(ip),   allocatable :: mc_table_inout_subcells_per_case(:,:)
    integer(ip),   allocatable :: mc_table_subcell_node_ids_per_case(:,:,:)

    ! Info related to cut cells on this triangulation (this is computed at runtime)
    integer(ip),   allocatable :: mc_case_per_cell(:)        ! Num local cells
    integer(ip),   allocatable :: mc_ptr_to_intersections(:) ! Num local cells + 1
    type(point_t), allocatable :: mc_intersection_points(:)  ! mc_ptr_to_intersections(num_local_cells+num_ghost_cells+1) - mc_ptr_to_intersections(1)

  contains
    procedure, private :: serial_triangulation_create                        => serial_unfitted_triangulation_serial_triangulation_create
    procedure, private :: serial_unfitted_triangulation_create
    generic            :: create                                          => serial_unfitted_triangulation_create
    procedure          :: free                                            => serial_unfitted_triangulation_free
    procedure, non_overridable          :: create_unfitted_cell_iterator  => serial_unfitted_triangulation_create_unfitted_cell_iterator
    procedure, non_overridable, private :: mc_tables_create               => serial_unfitted_triangulation_mc_tables_create
    procedure, non_overridable, private :: mc_tables_free                 => serial_unfitted_triangulation_mc_tables_free
    procedure, non_overridable, private :: mc_runtime_info_create         => serial_unfitted_triangulation_mc_runtime_info_create
    procedure, non_overridable, private :: mc_runtime_info_free           => serial_unfitted_triangulation_mc_runtime_info_free
    procedure, non_overridable, private :: determine_mc_case              => serial_unfitted_triangulation_determine_mc_case
    procedure, non_overridable, private :: compute_mc_intersection_points => serial_unfitted_triangulation_compute_mc_intersection_points

  end type serial_unfitted_triangulation_t

  ! Derived types
  !public :: unfitted_cell_accessor_t
  !public :: unfitted_cell_iterator_t
  public :: serial_unfitted_triangulation_t

contains

#include "sbm_unfitted_cell_accessor.i90"
#include "sbm_unfitted_cell_iterator.i90"
#include "sbm_serial_unfitted_triangulation.i90"

end module serial_unfitted_triangulation_names
  
