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
  use IR_Precision ! VTK_IO
  use Lib_VTK_IO ! VTK_IO
  
  implicit none
# include "debug.i90"
  private

  ! Include the look-up tables 
# include "mc_tables_qua4.i90"
# include "mc_tables_hex8.i90"
  
  type, extends(cell_accessor_t) :: unfitted_cell_accessor_t
  
    private
    class(marching_cubes_t), pointer :: marching_cubes => NULL()
    
  contains

    ! Creation /deletion methods
    generic :: create => unfitted_cell_accessor_create
    procedure, private :: unfitted_cell_accessor_create
    procedure :: cell_accessor_create  => unfitted_cell_accessor_cell_accessor_create
    procedure :: free   => unfitted_cell_accessor_cell_accessor_free
    
    ! To be called each time the lid changes
    procedure, non_overridable :: update_sub_triangulation    => unfitted_cell_accessor_update_sub_triangulation

    ! TODO this one should be private in the future
    ! Now it is public because the fe accessor uses it
    ! The goal is that the fe accessor does not assume that the mc algorithm is used to subdivide the elements
    procedure, non_overridable :: get_mc_case   => unfitted_cell_accessor_get_mc_case
    
    ! TBPs related with the subcells
    procedure, non_overridable :: get_number_of_subcells      => unfitted_cell_accessor_get_number_of_subcells
    procedure, non_overridable :: get_number_of_subcell_nodes => unfitted_cell_accessor_get_number_of_subcell_nodes
    procedure, non_overridable :: get_phys_coords_of_subcell  => unfitted_cell_accessor_get_phys_coords_of_subcell
    procedure, non_overridable :: get_ref_coords_of_subcell   => unfitted_cell_accessor_get_ref_coords_of_subcell
    
    ! TBPs related with the subfaces
    procedure, non_overridable :: get_number_of_subfaces      => unfitted_cell_accessor_get_number_of_subfaces
    procedure, non_overridable :: get_number_of_subface_nodes => unfitted_cell_accessor_get_number_of_subface_nodes
    procedure, non_overridable :: get_phys_coords_of_subface  => unfitted_cell_accessor_get_phys_coords_of_subface
    procedure, non_overridable :: get_ref_coords_of_subface   => unfitted_cell_accessor_get_ref_coords_of_subface
    
    procedure, non_overridable :: is_cut      => unfitted_cell_accessor_is_cut
    procedure, non_overridable :: is_interior => unfitted_cell_accessor_is_interior
    procedure, non_overridable :: is_exterior => unfitted_cell_accessor_is_exterior
    procedure, non_overridable :: is_interior_subcell => unfitted_cell_accessor_is_interior_subcell
    procedure, non_overridable :: is_exterior_subcell => unfitted_cell_accessor_is_exterior_subcell

    ! Private TBPs
    procedure, non_overridable, private :: get_number_of_subnodes         => unfitted_cell_accessor_get_number_of_subnodes
    procedure, non_overridable, private :: subcell_has_been_reoriented    => unfitted_cell_accessor_subcell_has_been_reoriented
    procedure, non_overridable, private :: subface_touches_interior_reoriented_subcell => unfitted_cell_accessor_subface_touches_reoriented_subcell

    
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

  type :: marching_cubes_t
    private

    ! The underlying triangulation
    class(base_static_triangulation_t), pointer :: triangulation => null()

    ! The level set funciton
    class(level_set_function_t), pointer :: level_set_function => null()

    ! Look up-tables (precomputed off-line, for each cell type)
    integer(ip)                :: mc_table_num_cases
    integer(ip)                :: mc_table_max_num_subcells
    integer(ip)                :: mc_table_max_num_subfaces
    integer(ip)                :: mc_table_max_num_cut_edges
    integer(ip)                :: mc_table_num_nodes_subcell
    integer(ip)                :: mc_table_num_nodes_subface
    integer(ip),   allocatable :: mc_table_num_subcells_per_case(:)
    integer(ip),   allocatable :: mc_table_num_subfaces_per_case(:)
    integer(ip),   allocatable :: mc_table_num_cut_edges_per_case(:)
    integer(ip),   allocatable :: mc_table_inout_subcells_per_case(:,:)
    integer(ip),   allocatable :: mc_table_subcell_node_ids_per_case(:,:,:)
    integer(ip),   allocatable :: mc_table_subface_node_ids_per_case(:,:,:) !TODO Provide also boundary sub-faces??
    logical :: mc_tables_init = .false.

    ! Info related to cut cells on this triangulation (this is computed at runtime)
    integer(ip),   allocatable :: mc_case_per_cell(:)
    integer(ip),   allocatable :: mc_ptr_to_intersections(:)
    type(point_t), allocatable :: mc_intersection_points(:)
    logical :: mc_runtime_init = .false.

    ! Info related to the sub-tessellation
    ! When located on a cell, and the sub-triangulation is updated,
    ! these memeber variables contain info about the current sub-tessalation
    integer(ip),        allocatable :: subcells_nodal_connectivities(:,:)
    integer(ip),        allocatable :: subfaces_nodal_connectivities(:,:)
    type(point_t),      allocatable :: subnodes_ref_coords(:)
    logical,            allocatable :: subcell_has_been_reoriented(:)
    
    ! Auxiliary work data
    type(quadrature_t), allocatable :: subnodes_nodal_quadratures(:)
    type(fe_map_t),     allocatable :: subnodes_fe_maps(:)

  contains

    ! Creation / deletion methods
    procedure                  :: create                        => marching_cubes_create
    procedure                  :: free                          => marching_cubes_free
    procedure, non_overridable :: create_unfitted_cell_iterator => marching_cubes_create_unfitted_cell_iterator

    ! Getters
    procedure, non_overridable :: get_num_cut_cells             => marching_cubes_get_num_cut_cells
    procedure, non_overridable :: get_num_interior_cells        => marching_cubes_get_num_interior_cells
    procedure, non_overridable :: get_num_exterior_cells        => marching_cubes_get_num_exterior_cells
    procedure, non_overridable :: get_max_num_subcells_in_cell  => marching_cubes_get_max_subcells_in_cell
    procedure, non_overridable :: get_max_num_nodes_in_subcell  => marching_cubes_get_max_nodes_in_subcell
    procedure, non_overridable :: get_total_num_of_subcells     => marching_cubes_get_total_num_of_subcells
    procedure, non_overridable :: get_max_num_subfaces_in_cell  => marching_cubes_get_max_subfaces_in_cell
    procedure, non_overridable :: get_max_num_nodes_in_subface  => marching_cubes_get_max_nodes_in_subface
    procedure, non_overridable :: get_total_num_of_subfaces     => marching_cubes_get_total_num_of_subfaces
    procedure, non_overridable :: get_max_num_subnodes_in_cell  => marching_cubes_get_max_subnodes_in_cell
    
    ! TODO this getters should be private in the future
    ! Now the fe space uses them.
    ! The goal is that the fe space does not assume that the tesselation algorithm is the mc algorithm
    procedure, non_overridable :: get_num_mc_cases              => marching_cubes_get_num_mc_cases
    procedure, non_overridable :: get_num_subcells_mc_case      => marching_cubes_get_num_subcells_mc_case
    procedure, non_overridable :: get_num_subfaces_mc_case      => marching_cubes_get_num_subfaces_mc_case
    
    ! Printers
    procedure :: print                     => marching_cubes_print
    ! TODO move to output hadler
    procedure :: print_to_vtk_file         => marching_cubes_print_to_vtk_file
    
    ! Private TBP
    procedure, non_overridable, private :: fulfills_assumptions           => marching_cubes_fulfills_assumptions
    procedure, non_overridable, private :: mc_tables_create               => marching_cubes_mc_tables_create
    procedure, non_overridable, private :: mc_tables_free                 => marching_cubes_mc_tables_free
    procedure, non_overridable, private :: mc_runtime_info_create         => marching_cubes_mc_runtime_info_create
    procedure, non_overridable, private :: mc_runtime_info_free           => marching_cubes_mc_runtime_info_free
    procedure, non_overridable, private :: subnodes_data_create           => marching_cubes_subnodes_data_create
    procedure, non_overridable, private :: subnodes_data_free             => marching_cubes_subnodes_data_free
    procedure, non_overridable, private :: print_vtk_subcells             => marching_cubes_print_vtk_subcells
    procedure, non_overridable, private :: print_vtk_subfaces             => marching_cubes_print_vtk_subfaces
    

  end type marching_cubes_t

  type, extends(serial_triangulation_t) :: serial_unfitted_triangulation_t
    private
      type(marching_cubes_t) :: marching_cubes
    contains

      ! Creation / deletion methods
      generic             :: create                       => nsut_create
      procedure           :: free                         => nsut_free
      procedure,  private :: serial_triangulation_create  => nsut_serial_triangulation_create
      procedure,  private :: nsut_create

      ! Generate iterator
      procedure, non_overridable :: create_unfitted_cell_iterator => nsut_create_unfitted_cell_iterator

      ! Getters needed by other classes
      ! TODO put here all the other getters of the aggregate

      ! By vtk writer
      procedure, non_overridable :: get_max_num_nodes_in_subcell  => nsut_get_max_num_nodes_in_subcell
      procedure, non_overridable :: get_total_num_of_subcells     => nsut_get_total_num_of_subcells

      ! By unfitted fe space
      procedure, non_overridable :: get_total_num_of_subfaces     => nsut_get_total_num_of_subfaces

      ! TODO this getters should be removed in the future
      ! Now the fe space uses them.
      ! The goal is that the fe space does not assume that the tesselation algorithm is the mc algorithm
      procedure, non_overridable :: get_num_mc_cases              => nsut_get_num_mc_cases
      procedure, non_overridable :: get_num_subcells_mc_case      => nsut_get_num_subcells_mc_case
      procedure, non_overridable :: get_num_subfaces_mc_case      => nsut_get_num_subfaces_mc_case

      ! Printers
      procedure :: print                     => nsut_print
      procedure :: print_to_vtk_file         => nsut_print_to_vtk_file

  end type serial_unfitted_triangulation_t

  ! Derived types
  public :: unfitted_cell_accessor_t
  public :: unfitted_cell_iterator_t
  public :: serial_unfitted_triangulation_t

contains

#include "sbm_unfitted_cell_accessor.i90"
#include "sbm_unfitted_cell_iterator.i90"
#include "sbm_serial_unfitted_triangulation.i90"

! TODO put these in a submodule (or not?)

!========================================================================================
subroutine nsut_create(this, parameters, level_set_function)
  implicit none
  class(serial_unfitted_triangulation_t), target, intent(inout) :: this
  type(ParameterList_t)                         , intent(inout) :: parameters
  class(level_set_function_t),            target, intent(in)    :: level_set_function
  call this%free()
  call this%serial_triangulation_t%create(parameters)
  call this%marching_cubes%create(this, level_set_function)
end subroutine nsut_create

!========================================================================================
subroutine nsut_free(this)
  implicit none
  class(serial_unfitted_triangulation_t), target, intent(inout) :: this
  call this%serial_triangulation_t%free()
  call this%marching_cubes%free()
end subroutine nsut_free

!========================================================================================
subroutine nsut_serial_triangulation_create( this, parameters)
  implicit none
  class(serial_unfitted_triangulation_t), target, intent(inout) :: this
  type(ParameterList_t)             , intent(inout) :: parameters
  check(.false.)
end subroutine nsut_serial_triangulation_create

!========================================================================================
function nsut_create_unfitted_cell_iterator ( this ) result (cell_iter)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  type(unfitted_cell_iterator_t) :: cell_iter
  cell_iter = this%marching_cubes%create_unfitted_cell_iterator()
end function nsut_create_unfitted_cell_iterator

!========================================================================================
function nsut_get_max_num_nodes_in_subcell(this) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: val
  val = this%marching_cubes%get_max_num_nodes_in_subcell()
end function nsut_get_max_num_nodes_in_subcell

!========================================================================================
function nsut_get_total_num_of_subcells(this) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: val
  val = this%marching_cubes%get_total_num_of_subcells()
end function nsut_get_total_num_of_subcells

!========================================================================================
function nsut_get_num_mc_cases(this) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: val
  val = this%marching_cubes%get_num_mc_cases()
end function nsut_get_num_mc_cases

!========================================================================================
function nsut_get_total_num_of_subfaces(this) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: val
  val = this%marching_cubes%get_total_num_of_subfaces()
end function nsut_get_total_num_of_subfaces

!========================================================================================
function nsut_get_num_subcells_mc_case(this,mc_case) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: mc_case
  integer(ip) :: val
  val = this%marching_cubes%get_num_subcells_mc_case(mc_case)
end function nsut_get_num_subcells_mc_case

!========================================================================================
function nsut_get_num_subfaces_mc_case(this,mc_case) result (val)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  integer(ip) :: mc_case
  integer(ip) :: val
  val = this%marching_cubes%get_num_subfaces_mc_case(mc_case)
end function nsut_get_num_subfaces_mc_case

!========================================================================================
subroutine nsut_print(this)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  call this%marching_cubes%print()
end subroutine nsut_print

!========================================================================================
subroutine nsut_print_to_vtk_file(this)
  implicit none
  class(serial_unfitted_triangulation_t), intent(in)    :: this
  call this%marching_cubes%print_to_vtk_file()
end subroutine nsut_print_to_vtk_file

end module serial_unfitted_triangulation_names

