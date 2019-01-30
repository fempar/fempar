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

module unfitted_boundary_cutter_names
  use types_names
  use field_names
  use triangulation_names
  use reference_fe_names
  
  implicit none
  private
#include "debug.i90"
  
  type, abstract :: unfitted_boundary_cutter_t
    private
    class(triangulation_t) , pointer :: triangulation => null()
  contains
    ! Regular TBPs
    procedure, non_overridable :: set_triangulation     => ubc_set_triangulation
    procedure, non_overridable :: get_triangulation     => ubc_get_triangulation
    procedure, non_overridable :: nullify_triangulation => ubc_nullify_triangulation
    
    ! We are not sure whether this set of 4x TBPs should actually be here taking into
    ! account how we conceptually conceive unfitted_boundary_cutter_t. There are
    ! some smells in the relationship among unfitted_cell_iterator_t and unfitted_boundary_cutter_t
    ! that have to be tackled and resolved
    procedure, non_overridable :: create_cell_iterator  => ubc_create_cell_iterator
    procedure, non_overridable :: free_cell_iterator    => ubc_free_cell_iterator
    procedure, non_overridable :: create_vef_iterator   => ubc_create_vef_iterator
    procedure, non_overridable :: free_vef_iterator     => ubc_free_vef_iterator
  
    ! Deferred TBPs
    procedure(ubc_free_interface),                            deferred :: free
    procedure(ubc_get_cell_gid_interface)                   , deferred :: get_cell_gid
    procedure(ubc_get_num_cut_cells_interface)              , deferred :: get_num_cut_cells
    procedure(ubc_get_num_interior_cells_interface)         , deferred :: get_num_interior_cells
    procedure(ubc_get_num_exterior_cells_interface)         , deferred :: get_num_exterior_cells
    procedure(ubc_get_max_num_nodes_in_subcell_interface)   , deferred :: get_max_num_nodes_in_subcell
    procedure(ubc_get_total_num_subcells_interface)         , deferred :: get_total_num_subcells
    procedure(ubc_get_max_num_nodes_in_subfacet_interface)  , deferred :: get_max_num_nodes_in_subfacet
    procedure(ubc_get_total_num_subfacets_interface)        , deferred :: get_total_num_subfacets
    procedure(ubc_get_total_num_fitted_sub_facets_interface), deferred :: get_total_num_fitted_sub_facets
    procedure(ubc_update_sub_triangulation_interface)       , deferred :: update_sub_triangulation
    procedure(ubc_get_num_subcells_interface)               , deferred :: get_num_subcells
    procedure(ubc_get_num_subcell_nodes_interface)          , deferred :: get_num_subcell_nodes
    procedure(ubc_get_phys_coords_of_subcell_interface)     , deferred :: get_phys_coords_of_subcell
    procedure(ubc_get_ref_coords_of_subcell_interface)      , deferred :: get_ref_coords_of_subcell
    procedure(ubc_get_num_subfacets_interface)              , deferred :: get_num_subfacets
    procedure(ubc_get_num_subfacet_nodes_interface)         , deferred :: get_num_subfacet_nodes
    procedure(ubc_get_phys_coords_of_subfacet_interface)    , deferred :: get_phys_coords_of_subfacet
    procedure(ubc_get_ref_coords_of_subfacet_interface)     , deferred :: get_ref_coords_of_subfacet
    procedure(ubc_is_cut_interface)                         , deferred :: is_cut
    procedure(ubc_is_interior_interface)                    , deferred :: is_interior
    procedure(ubc_is_exterior_interface)                    , deferred :: is_exterior
    procedure(ubc_is_interior_subcell_interface)            , deferred :: is_interior_subcell
    procedure(ubc_is_exterior_subcell_interface)            , deferred :: is_exterior_subcell
    procedure(ubc_get_num_fitted_subfacets_interface)       , deferred :: get_num_fitted_subfacets
    procedure(ubc_get_phys_coords_of_subvef_interface)      , deferred :: get_phys_coords_of_subvef
    procedure(ubc_get_ref_coords_of_subvef_interface)       , deferred :: get_ref_coords_of_subvef
    procedure(ubc_is_cut_facet_interface)                   , deferred :: is_cut_facet
    procedure(ubc_is_interior_facet_interface)              , deferred :: is_interior_facet
    procedure(ubc_is_exterior_facet_interface)              , deferred :: is_exterior_facet
    procedure(ubc_is_interior_subfacet_interface)           , deferred :: is_interior_subfacet
    procedure(ubc_is_exterior_subfacet_interface)           , deferred :: is_exterior_subfacet
    procedure(ubc_print_interface)                          , deferred :: print
  end type unfitted_boundary_cutter_t
  
  abstract interface
     subroutine ubc_free_interface ( this )
       import :: unfitted_boundary_cutter_t
       class(unfitted_boundary_cutter_t), target, intent(inout) :: this
     end subroutine ubc_free_interface

     function ubc_get_cell_gid_interface ( this ) result( gid )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: gid
     end function ubc_get_cell_gid_interface

     function ubc_get_num_cut_cells_interface ( this ) result( num_cut_cells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_cut_cells
     end function ubc_get_num_cut_cells_interface

     function ubc_get_num_interior_cells_interface ( this ) result( num_interior_cells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_interior_cells
     end function ubc_get_num_interior_cells_interface

     function ubc_get_num_exterior_cells_interface ( this ) result( num_exterior_cells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_exterior_cells
     end function ubc_get_num_exterior_cells_interface

     function ubc_get_max_num_nodes_in_subcell_interface ( this ) result( max_num_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: max_num_nodes
     end function ubc_get_max_num_nodes_in_subcell_interface

     function ubc_get_total_num_subcells_interface ( this ) result( total_num_subcells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: total_num_subcells
     end function ubc_get_total_num_subcells_interface

     function ubc_get_max_num_nodes_in_subfacet_interface ( this ) result( max_num_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: max_num_nodes
     end function ubc_get_max_num_nodes_in_subfacet_interface

     function ubc_get_total_num_subfacets_interface ( this ) result( total_num_subfacets )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: total_num_subfacets
     end function ubc_get_total_num_subfacets_interface

     function ubc_get_total_num_fitted_sub_facets_interface ( this ) result( total_num_fitted_sub_facets )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: total_num_fitted_sub_facets
     end function ubc_get_total_num_fitted_sub_facets_interface

     subroutine ubc_update_sub_triangulation_interface ( this, cell_iterator )
       import :: unfitted_boundary_cutter_t, cell_iterator_t, ip
       class(unfitted_boundary_cutter_t), target, intent(inout)    :: this
       class(cell_iterator_t),            intent(in)       :: cell_iterator
     end subroutine ubc_update_sub_triangulation_interface

     function ubc_get_num_subcells_interface ( this ) result( num_subcells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subcells
     end function ubc_get_num_subcells_interface

     function ubc_get_num_subcell_nodes_interface ( this ) result( num_subcell_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subcell_nodes
     end function ubc_get_num_subcell_nodes_interface

     subroutine ubc_get_phys_coords_of_subcell_interface ( this, subcell, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       integer(ip),                               intent(in)    :: subcell
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_phys_coords_of_subcell_interface

     subroutine ubc_get_ref_coords_of_subcell_interface ( this, subcell, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       integer(ip),                               intent(in)    :: subcell
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_ref_coords_of_subcell_interface

     function ubc_get_num_subfacets_interface ( this ) result( num_subfacets )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subfacets
     end function ubc_get_num_subfacets_interface

     function ubc_get_num_subfacet_nodes_interface ( this ) result( num_subfacet_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subfacet_nodes
     end function ubc_get_num_subfacet_nodes_interface

     subroutine ubc_get_phys_coords_of_subfacet_interface ( this, subfacet, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       integer(ip),                               intent(in)    :: subfacet
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_phys_coords_of_subfacet_interface

     subroutine ubc_get_ref_coords_of_subfacet_interface ( this, subfacet, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       integer(ip),                               intent(in)    :: subfacet
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_ref_coords_of_subfacet_interface

     function ubc_is_cut_interface ( this, gid ) result( is_cut )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       logical :: is_cut
     end function ubc_is_cut_interface

     function ubc_is_interior_interface ( this, gid ) result( is_interior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       logical :: is_interior
     end function ubc_is_interior_interface

     function ubc_is_exterior_interface ( this, gid ) result( is_exterior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       logical :: is_exterior
     end function ubc_is_exterior_interface

     function ubc_is_interior_subcell_interface ( this, subcell ) result( is_interior_subcell )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: subcell
       logical :: is_interior_subcell
     end function ubc_is_interior_subcell_interface

     function ubc_is_exterior_subcell_interface ( this, subcell ) result( is_exterior_subcell )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: subcell
       logical :: is_exterior_subcell
     end function ubc_is_exterior_subcell_interface

     function ubc_get_num_fitted_subfacets_interface ( this, gid, facet_lid ) result (num_fitted_subfacets)
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in) :: this
       integer(ip),                       intent(in) :: gid
       integer(ip),                       intent(in) :: facet_lid
       integer(ip) :: num_fitted_subfacets
     end function ubc_get_num_fitted_subfacets_interface

     subroutine ubc_get_phys_coords_of_subvef_interface ( this, gid, facet_lid, subvef, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       integer(ip),                               intent(in)    :: gid
       integer(ip),                               intent(in)    :: facet_lid
       integer(ip),                               intent(in)    :: subvef
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_phys_coords_of_subvef_interface

     subroutine ubc_get_ref_coords_of_subvef_interface ( this, reference_fe, gid, facet_lid, subvef, points )
       import :: unfitted_boundary_cutter_t, point_t, reference_fe_t, ip
       class(unfitted_boundary_cutter_t), target, intent(in)    :: this
       class(reference_fe_t),                     intent(in)    :: reference_fe
       integer(ip),                               intent(in)    :: gid
       integer(ip),                               intent(in)    :: facet_lid
       integer(ip),                               intent(in)    :: subvef
       type(point_t),                             intent(inout) :: points(:)
     end subroutine ubc_get_ref_coords_of_subvef_interface

     function ubc_is_cut_facet_interface ( this, gid, facet_lid ) result ( is_cut )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_cut 
     end function ubc_is_cut_facet_interface

     function ubc_is_interior_facet_interface ( this, gid, facet_lid ) result ( is_interior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_interior 
     end function ubc_is_interior_facet_interface

     function ubc_is_exterior_facet_interface ( this, gid, facet_lid ) result ( is_exterior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_exterior 
     end function ubc_is_exterior_facet_interface

     function ubc_is_interior_subfacet_interface ( this, gid, facet_lid, subvef ) result ( is_in )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       integer(ip),                       intent(in)    :: facet_lid
       integer(ip),                       intent(in)    :: subvef
       logical :: is_in
     end function ubc_is_interior_subfacet_interface

     function ubc_is_exterior_subfacet_interface ( this, gid, facet_lid, subvef ) result ( is_out )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: gid
       integer(ip),                       intent(in)    :: facet_lid
       integer(ip),                       intent(in)    :: subvef
       logical :: is_out
     end function ubc_is_exterior_subfacet_interface

     subroutine ubc_print_interface ( this )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
     end subroutine ubc_print_interface
 end interface

  ! Derived types
  public :: unfitted_boundary_cutter_t

contains 

  function ubc_get_triangulation ( this )
    implicit none
    class(unfitted_boundary_cutter_t), target, intent(in)    :: this
    class(triangulation_t), pointer :: ubc_get_triangulation
    ubc_get_triangulation => this%triangulation
  end function ubc_get_triangulation 
    
  subroutine ubc_set_triangulation ( this, triangulation )
    implicit none
    class(unfitted_boundary_cutter_t), intent(inout) :: this
    class(triangulation_t), target   , intent(in)    :: triangulation
    this%triangulation => triangulation
  end subroutine ubc_set_triangulation 
  
  subroutine ubc_nullify_triangulation ( this )
    implicit none
    class(unfitted_boundary_cutter_t), intent(inout) :: this
    nullify(this%triangulation)
  end subroutine ubc_nullify_triangulation
  
  subroutine ubc_create_cell_iterator(this, cell)
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    class(cell_iterator_t), allocatable, intent(inout) :: cell
    assert(associated(this%triangulation))
    call this%triangulation%create_cell_iterator(cell)
  end subroutine ubc_create_cell_iterator
  
  subroutine ubc_free_cell_iterator(this, cell)
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    class(cell_iterator_t), allocatable, intent(inout) :: cell
    assert(associated(this%triangulation))
    call this%triangulation%free_cell_iterator(cell)
  end subroutine ubc_free_cell_iterator
  
  subroutine ubc_create_vef_iterator(this, vef)
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    class(vef_iterator_t), allocatable, intent(inout) :: vef
    assert(associated(this%triangulation))
    call this%triangulation%create_vef_iterator(vef)
  end subroutine ubc_create_vef_iterator
  
  subroutine ubc_free_vef_iterator(this, vef)
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    class(vef_iterator_t), allocatable, intent(inout) :: vef
    assert(associated(this%triangulation))
    call this%triangulation%free_vef_iterator(vef)
  end subroutine ubc_free_vef_iterator
 
end module unfitted_boundary_cutter_names

