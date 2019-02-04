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
  
  ! In my view (@amartinhuertas), unfitted_boundary_cutter_t is conceptually seen as follows (perhaps we may
  ! think of a different/better conceptual definition). It is a mesh-like container
  ! able to describe the subtriangulation of cut cells. However, as its clients are cell iterators,
  ! it has to become aware of the current cell on which the iterator is positioned to resolve the queries
  ! related to the current cell.  I have made this concept explicit in the abstract class, i.e., the client
  ! has to signal explicitly via a call to "set_current_cell" where it is currently positioned.
  type, abstract :: unfitted_boundary_cutter_t
    private
    class(triangulation_t) , pointer :: triangulation => null()
    class(cell_iterator_t) , pointer :: current_cell  => null()
  contains
    ! Regular TBPs
    procedure, non_overridable :: set_triangulation     => ubc_set_triangulation
    procedure, non_overridable :: get_triangulation     => ubc_get_triangulation
    procedure, non_overridable :: nullify_triangulation => ubc_nullify_triangulation
    procedure                  :: set_current_cell      => ubc_set_current_cell
    procedure                  :: get_current_cell      => ubc_get_current_cell
    procedure                  :: nullify_current_cell  => ubc_nullify_current_cell
    
    ! We are not sure whether this set of 4x TBPs should actually be here taking into
    ! account how we conceptually conceive unfitted_boundary_cutter_t. In my view (@amartin)
    ! there are some smells in the relationship among unfitted_cell_iterator_t and unfitted_boundary_cutter_t
    ! that have to be tackled and resolved
    procedure, non_overridable :: create_cell_iterator  => ubc_create_cell_iterator
    procedure, non_overridable :: free_cell_iterator    => ubc_free_cell_iterator
    procedure, non_overridable :: create_vef_iterator   => ubc_create_vef_iterator
    procedure, non_overridable :: free_vef_iterator     => ubc_free_vef_iterator

    ! These procedures are overridable since the concrete implementations can provide
    ! more efficient code
    procedure :: get_num_cut_cells               => ubc_get_num_cut_cells
    procedure :: get_num_interior_cells          => ubc_get_num_interior_cells
    procedure :: get_num_exterior_cells          => ubc_get_num_exterior_cells
    procedure :: get_total_num_subcells          => ubc_get_total_num_subcells
    procedure :: get_total_num_subfacets         => ubc_get_total_num_subfacets
    procedure :: get_total_num_fitted_sub_facets => ubc_get_total_num_fitted_sub_facets
    procedure :: get_max_num_nodes_in_subcell    => ubc_get_max_num_nodes_in_subcell
    procedure :: get_max_num_nodes_in_subfacet   => ubc_get_max_num_nodes_in_subfacet
    procedure :: get_max_num_subcells_in_cell    => ubc_get_max_num_subcells_in_cell
    procedure :: get_max_num_subfacets_in_cell   => ubc_get_max_num_subfacets_in_cell
  
    ! Deferred BPs
    procedure(free_interface),                            deferred :: free
    procedure(get_num_subcells_interface)               , deferred :: get_num_subcells
    procedure(get_num_subcell_nodes_interface)          , deferred :: get_num_subcell_nodes
    procedure(get_phys_coords_of_subcell_interface)     , deferred :: get_phys_coords_of_subcell
    procedure(get_ref_coords_of_subcell_interface)      , deferred :: get_ref_coords_of_subcell
    procedure(get_num_subfacets_interface)              , deferred :: get_num_subfacets
    procedure(get_num_subfacet_nodes_interface)         , deferred :: get_num_subfacet_nodes
    procedure(get_phys_coords_of_subfacet_interface)    , deferred :: get_phys_coords_of_subfacet
    procedure(get_ref_coords_of_subfacet_interface)     , deferred :: get_ref_coords_of_subfacet
    procedure(is_cut_interface)                         , deferred :: is_cut
    procedure(is_interior_interface)                    , deferred :: is_interior
    procedure(is_exterior_interface)                    , deferred :: is_exterior
    procedure(is_interior_subcell_interface)            , deferred :: is_interior_subcell
    procedure(is_exterior_subcell_interface)            , deferred :: is_exterior_subcell
    procedure(get_num_fitted_subfacets_interface)       , deferred :: get_num_fitted_subfacets
    procedure(get_phys_coords_of_subvef_interface)      , deferred :: get_phys_coords_of_subvef
    procedure(get_ref_coords_of_subvef_interface)       , deferred :: get_ref_coords_of_subvef
    procedure(is_cut_facet_interface)                   , deferred :: is_cut_facet
    procedure(is_interior_facet_interface)              , deferred :: is_interior_facet
    procedure(is_exterior_facet_interface)              , deferred :: is_exterior_facet
    procedure(is_interior_subfacet_interface)           , deferred :: is_interior_subfacet
    procedure(is_exterior_subfacet_interface)           , deferred :: is_exterior_subfacet
    procedure(print_interface)                          , deferred :: print
  end type unfitted_boundary_cutter_t
  
  abstract interface
     subroutine free_interface ( this )
       import :: unfitted_boundary_cutter_t
       class(unfitted_boundary_cutter_t), intent(inout) :: this
     end subroutine free_interface

     function get_num_subcells_interface ( this ) result( num_subcells )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subcells
     end function get_num_subcells_interface

     function get_num_subcell_nodes_interface ( this ) result( num_subcell_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subcell_nodes
     end function get_num_subcell_nodes_interface

     subroutine get_phys_coords_of_subcell_interface ( this, subcell, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip)                      , intent(in)    :: subcell
       type(point_t)                    , intent(inout) :: points(:)
     end subroutine get_phys_coords_of_subcell_interface

     subroutine get_ref_coords_of_subcell_interface ( this, subcell, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip)                      , intent(in)    :: subcell
       type(point_t)                    , intent(inout) :: points(:)
     end subroutine get_ref_coords_of_subcell_interface

     function get_num_subfacets_interface ( this ) result( num_subfacets )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subfacets
     end function get_num_subfacets_interface

     function get_num_subfacet_nodes_interface ( this ) result( num_subfacet_nodes )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip) :: num_subfacet_nodes
     end function get_num_subfacet_nodes_interface

     subroutine get_phys_coords_of_subfacet_interface ( this, subfacet, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip)                      , intent(in)    :: subfacet
       type(point_t)                    , intent(inout) :: points(:)
     end subroutine get_phys_coords_of_subfacet_interface

     subroutine get_ref_coords_of_subfacet_interface ( this, subfacet, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip)                      , intent(in)    :: subfacet
       type(point_t)                    , intent(inout) :: points(:)
     end subroutine get_ref_coords_of_subfacet_interface

     function is_cut_interface ( this ) result( is_cut )
       import :: unfitted_boundary_cutter_t
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       logical :: is_cut
     end function is_cut_interface

     function is_interior_interface ( this ) result( is_interior )
       import :: unfitted_boundary_cutter_t
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       logical :: is_interior
     end function is_interior_interface

     function is_exterior_interface ( this ) result( is_exterior )
       import :: unfitted_boundary_cutter_t
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       logical :: is_exterior
     end function is_exterior_interface

     function is_interior_subcell_interface ( this, subcell ) result( is_interior_subcell )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: subcell
       logical :: is_interior_subcell
     end function is_interior_subcell_interface

     function is_exterior_subcell_interface ( this, subcell ) result( is_exterior_subcell )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: subcell
       logical :: is_exterior_subcell
     end function is_exterior_subcell_interface

     function get_num_fitted_subfacets_interface ( this, facet_lid ) result (num_fitted_subfacets)
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in) :: this
       integer(ip),                       intent(in) :: facet_lid
       integer(ip) :: num_fitted_subfacets
     end function get_num_fitted_subfacets_interface

     subroutine get_phys_coords_of_subvef_interface ( this, facet_lid, subvef, points )
       import :: unfitted_boundary_cutter_t, point_t, ip
       class(unfitted_boundary_cutter_t),  intent(in)    :: this
       integer(ip)                      ,  intent(in)    :: facet_lid
       integer(ip)                      ,  intent(in)    :: subvef
       type(point_t)                    ,  intent(inout) :: points(:)
     end subroutine get_phys_coords_of_subvef_interface

     subroutine get_ref_coords_of_subvef_interface ( this, reference_fe, facet_lid, subvef, points )
       import :: unfitted_boundary_cutter_t, point_t, reference_fe_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       class(reference_fe_t)            , intent(in)    :: reference_fe
       integer(ip)                      , intent(in)    :: facet_lid
       integer(ip)                      , intent(in)    :: subvef
       type(point_t)                    , intent(inout) :: points(:)
     end subroutine get_ref_coords_of_subvef_interface

     function is_cut_facet_interface ( this, facet_lid ) result ( is_cut )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_cut 
     end function is_cut_facet_interface

     function is_interior_facet_interface ( this, facet_lid ) result ( is_interior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_interior 
     end function is_interior_facet_interface

     function is_exterior_facet_interface ( this, facet_lid ) result ( is_exterior )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: facet_lid
       logical :: is_exterior 
     end function is_exterior_facet_interface

     function is_interior_subfacet_interface ( this, facet_lid, subvef ) result ( is_in )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: facet_lid
       integer(ip),                       intent(in)    :: subvef
       logical :: is_in
     end function is_interior_subfacet_interface

     function is_exterior_subfacet_interface ( this, facet_lid, subvef ) result ( is_out )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
       integer(ip),                       intent(in)    :: facet_lid
       integer(ip),                       intent(in)    :: subvef
       logical :: is_out
     end function is_exterior_subfacet_interface

     subroutine print_interface ( this )
       import :: unfitted_boundary_cutter_t, ip
       class(unfitted_boundary_cutter_t), intent(in)    :: this
     end subroutine print_interface
 end interface

  ! Derived types
  public :: unfitted_boundary_cutter_t
  public :: ubc_set_current_cell

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
  
  subroutine ubc_set_current_cell ( this, cell ) 
    implicit none
    class(unfitted_boundary_cutter_t), intent(inout) :: this
    class(cell_iterator_t), target   , intent(in)    :: cell
#ifdef DEBUG
    class(triangulation_t), pointer :: triangulation
    assert(associated(this%triangulation))
    triangulation => cell%get_triangulation()
    assert(associated(this%triangulation,triangulation))
#endif    
    this%current_cell => cell    
  end subroutine ubc_set_current_cell 
  
  function ubc_get_current_cell ( this ) 
    implicit none
    class(unfitted_boundary_cutter_t), target, intent(in) :: this
    class(cell_iterator_t), pointer :: ubc_get_current_cell
    ubc_get_current_cell => this%current_cell
  end function ubc_get_current_cell 
  
  subroutine ubc_nullify_current_cell (this) 
    implicit none
    class(unfitted_boundary_cutter_t), intent(inout) :: this
    nullify(this%current_cell)
  end subroutine ubc_nullify_current_cell
  
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

  function ubc_get_num_cut_cells( this ) result ( num_cut_cells )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: num_cut_cells
    class(cell_iterator_t), allocatable  :: cell
    num_cut_cells = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( cell%is_cut() ) then
        num_cut_cells = num_cut_cells + 1_ip
      end if
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_num_cut_cells

  function ubc_get_num_interior_cells( this ) result ( num_interior_cells )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: num_interior_cells
    class(cell_iterator_t), allocatable  :: cell
    num_interior_cells = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( cell%is_interior() ) then
        num_interior_cells = num_interior_cells + 1_ip
      end if
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_num_interior_cells
  
  function ubc_get_num_exterior_cells( this ) result ( num_exterior_cells )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: num_exterior_cells
    class(cell_iterator_t), allocatable  :: cell
    num_exterior_cells = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( cell%is_exterior() ) then
        num_exterior_cells = num_exterior_cells + 1_ip
      end if
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_num_exterior_cells

  function ubc_get_total_num_subcells( this ) result ( total_num )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: total_num
    class(cell_iterator_t), allocatable  :: cell
    total_num = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if (cell%is_ghost()) then
        call cell%next(); cycle
      end if
      total_num = total_num + cell%get_num_subcells()
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_total_num_subcells
  
  function ubc_get_total_num_subfacets( this ) result ( total_num )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: total_num
    class(cell_iterator_t), allocatable  :: cell
    total_num = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      total_num = total_num + cell%get_num_subfacets()
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_total_num_subfacets
  
  function ubc_get_total_num_fitted_sub_facets( this ) result ( total_num )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: total_num
    class(vef_iterator_t), allocatable  :: vef
    total_num = 0_ip
    call this%create_vef_iterator(vef)
    do while ( .not. vef%has_finished() )
      if ( vef%is_facet() ) then
        total_num = total_num + vef%get_num_subvefs()
      end if
      call vef%next()
    end do
    call this%free_vef_iterator(vef)
  end function ubc_get_total_num_fitted_sub_facets

  function ubc_get_max_num_nodes_in_subcell( this ) result ( max_nodes_in_subcell )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: max_nodes_in_subcell
    class(cell_iterator_t), allocatable  :: cell
    call this%create_cell_iterator(cell)
    max_nodes_in_subcell = 0_ip
    do while ( .not. cell%has_finished() )
      max_nodes_in_subcell = max( max_nodes_in_subcell, cell%get_num_subcell_nodes() )
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_max_num_nodes_in_subcell
  
  function ubc_get_max_num_nodes_in_subfacet( this ) result ( max_nodes_in_subfacet)
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: max_nodes_in_subfacet
    class(cell_iterator_t), allocatable  :: cell
    call this%create_cell_iterator(cell)
    max_nodes_in_subfacet = 0_ip
    do while ( .not. cell%has_finished() )
      max_nodes_in_subfacet = max( max_nodes_in_subfacet, cell%get_num_subfacet_nodes() )
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_max_num_nodes_in_subfacet

  function ubc_get_max_num_subcells_in_cell( this ) result ( max_subcells_in_cell )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: max_subcells_in_cell
    class(cell_iterator_t), allocatable  :: cell
    max_subcells_in_cell = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( cell%is_local() ) then
        if ( cell%is_cut() ) then
          max_subcells_in_cell = max(max_subcells_in_cell,cell%get_num_subcells())
        end if
      end if
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_max_num_subcells_in_cell
  
  function ubc_get_max_num_subfacets_in_cell( this ) result ( max_subfacets_in_cell )
    implicit none
    class(unfitted_boundary_cutter_t)  , intent(in)    :: this
    integer(ip) :: max_subfacets_in_cell
    class(cell_iterator_t), allocatable  :: cell
    max_subfacets_in_cell = 0_ip
    call this%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( cell%is_local() ) then
        if ( cell%is_cut() ) then
          max_subfacets_in_cell = max(max_subfacets_in_cell,cell%get_num_subfacets())
        end if
      end if
      call cell%next()
    end do
    call this%free_cell_iterator(cell)
  end function ubc_get_max_num_subfacets_in_cell

end module unfitted_boundary_cutter_names

