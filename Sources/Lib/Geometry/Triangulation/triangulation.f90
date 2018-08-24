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
module triangulation_names
  
  ! Modules required by triangulation_t
  use types_names
  use memor_names
  use reference_fe_names
  use cell_import_names
  use hash_table_names
  use list_types_names
  use field_names
  use environment_names
  use std_vector_integer_ip_names
  
  ! Modules required by triangulation_t implementors
  ! Serial modules
  use sort_names
  use reference_fe_names
  use mesh_names
  use mesh_distribution_names
  use par_io_names
  use stdio_names
  use FPL
  use uniform_hex_mesh_generator_names
  
  ! Geometry modules
  use sisl_names
  use geometry_names
  use gid_geometry_reader_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: cell_iterator_t
    private
    integer(ip)                     :: gid = -1
    class(triangulation_t), pointer :: triangulation => NULL()
  contains
    procedure                            :: create                  => cell_iterator_create
    procedure                            :: free                    => cell_iterator_free
    procedure                            :: first                   => cell_iterator_first
    procedure                            :: next                    => cell_iterator_next
    procedure                            :: has_finished            => cell_iterator_has_finished
    procedure, non_overridable           :: get_gid                 => cell_iterator_get_gid
    procedure                            :: set_gid                 => cell_iterator_set_gid
    procedure, non_overridable           :: get_triangulation       => cell_iterator_get_triangulation
    procedure, non_overridable           :: get_my_subpart          => cell_iterator_get_mysubpart
    procedure, non_overridable           :: get_my_subpart_lid      => cell_iterator_get_mysubpart_lid
        
    ! Cell-data related TBPs
    procedure(cell_get_ggid_interface)         , deferred :: get_ggid
    procedure(get_my_part_interface)           , deferred :: get_my_part
    procedure(cell_is_local_interface)         , deferred :: is_local
    procedure(cell_is_ghost_interface)         , deferred :: is_ghost

    ! Topology-data related TBPs 
    procedure(get_num_vefs_interface)          , deferred :: get_num_vefs
    procedure(get_vef_interface)               , deferred :: get_vef
    procedure(get_vef_gid_interface)           , deferred :: get_vef_gid
    procedure(get_vef_ggid_interface)          , deferred :: get_vef_ggid
    procedure(get_vefs_gid_interface)          , deferred :: get_vefs_gid
    procedure(get_vef_lid_from_gid_interface)  , deferred :: get_vef_lid_from_gid
    procedure(get_vef_lid_from_ggid_interface) , deferred :: get_vef_lid_from_ggid
    
    ! Transformation among n-face-wise node local numberings
    procedure(get_permutation_index_interface), deferred :: get_permutation_index
    
    ! Geometry-related data
    procedure(get_reference_fe_interface)          , deferred :: get_reference_fe
    procedure(get_reference_fe_id_interface)       , deferred :: get_reference_fe_id
    procedure(cell_get_num_nodes_interface)        , deferred :: get_num_nodes
    procedure(cell_get_nodes_coordinates_interface), deferred :: get_nodes_coordinates
    
    ! Set IDs-related TBPs
    procedure(cell_get_set_id_interface)             , deferred :: get_set_id
    procedure(cell_set_set_id_interface)             , deferred :: set_set_id
    ! Disconnected part cell set id 
    procedure(cell_get_disconnected_set_id_interface), deferred :: get_disconnected_set_id
 
    ! XFEM-related TBPs
    procedure(update_sub_triangulation_interface)   , deferred :: update_sub_triangulation
    procedure(get_num_subcells_interface)           , deferred :: get_num_subcells
    procedure(get_num_subcell_nodes_interface)      , deferred :: get_num_subcell_nodes
    procedure(get_phys_coords_of_subcell_interface) , deferred :: get_phys_coords_of_subcell
    procedure(get_ref_coords_of_subcell_interface)  , deferred :: get_ref_coords_of_subcell
    procedure(get_num_subfacets_interface)          , deferred :: get_num_subfacets
    procedure(get_num_subfacet_nodes_interface)     , deferred :: get_num_subfacet_nodes
    procedure(get_phys_coords_of_subfacet_interface), deferred :: get_phys_coords_of_subfacet
    procedure(get_ref_coords_of_subfacet_interface) , deferred :: get_ref_coords_of_subfacet
    procedure(is_cut_interface)                     , deferred :: is_cut
    procedure(is_exterior_interface)                , deferred :: is_exterior
    procedure(is_interior_interface)                , deferred :: is_interior
    procedure(is_exterior_subcell_interface)        , deferred :: is_exterior_subcell
    procedure(is_interior_subcell_interface)        , deferred :: is_interior_subcell

    ! h-adaptivity related TBPs
    procedure(get_level_interface)                  , deferred :: get_level
    procedure(set_for_refinement_interface)         , deferred :: set_for_refinement
    procedure(set_for_coarsening_interface)         , deferred :: set_for_coarsening
    procedure(set_for_do_nothing_interface)         , deferred :: set_for_do_nothing
    procedure(set_weight_interface)                 , deferred :: set_weight
    ! set_for_do_nothing ?
    ! get_transformation_flag
    ! ??? What else required
  end type cell_iterator_t
 
  type, abstract :: vef_iterator_t
    private
    integer(ip)                     :: gid = -1
    class(triangulation_t), pointer :: triangulation => NULL()
  contains
     procedure                           :: create                    => vef_iterator_create
     ! This regular TBP is required in extensions of vef_iterator_t that require to :
     !       (1) overwrite create; 
     !       (2) call the create TBP of the parent;
     ! Note that the Fortran200X standard does not allow to call the overrided regular method of a parent
     ! data type when the parent is abstract
     procedure, non_overridable          :: create_for_extensions     => vef_iterator_create_for_extensions
     procedure                           :: free                      => vef_iterator_free
     procedure, non_overridable          :: free_for_extensions       => vef_iterator_free_for_extensions
     procedure                           :: first                     => vef_iterator_first
     procedure                           :: next                      => vef_iterator_next
     procedure                           :: has_finished              => vef_iterator_has_finished
     procedure, non_overridable          :: get_gid                   => vef_iterator_get_gid
     procedure                           :: set_gid                   => vef_iterator_set_gid
     procedure, non_overridable          :: get_triangulation         => vef_iterator_get_triangulation
     procedure, non_overridable          :: is_facet                  => vef_iterator_is_facet
     
     ! Topology-data related TBPs
     procedure(get_num_cells_around_interface), deferred :: get_num_cells_around
     procedure(get_cell_around_interface)     , deferred :: get_cell_around
     
     ! Geometry related-data TBPs
     procedure(get_num_nodes_interface)        , deferred :: get_num_nodes
     procedure(get_nodes_coordinates_interface), deferred :: get_nodes_coordinates
     
     ! Misc TBPs
     procedure(get_ggid_interface)            , deferred :: get_ggid
     procedure(is_at_interior_interface)      , deferred :: is_at_interior
     procedure(is_at_boundary_interface)      , deferred :: is_at_boundary
     procedure(get_dim_interface)             , deferred :: get_dim
     procedure(get_set_id_interface)          , deferred :: get_set_id
     procedure(set_set_id_interface)          , deferred :: set_set_id
     procedure(is_ghost_interface)            , deferred :: is_ghost
     procedure(is_local_interface)            , deferred :: is_local 
     procedure(is_at_interface)               , deferred :: is_at_interface

     ! XFEM-related TBPs
     procedure               :: update_sub_triangulation  => vef_iterator_update_sub_triangulation
     procedure               :: get_num_subvefs           => vef_iterator_get_num_subvefs
     procedure               :: get_num_subvef_nodes      => vef_iterator_get_num_subvef_nodes
     procedure               :: get_phys_coords_of_subvef => vef_iterator_get_phys_coords_of_subvef
     procedure               :: get_ref_coords_of_subvef  => vef_iterator_get_ref_coords_of_subvef
     procedure               :: is_cut                    => vef_iterator_is_cut
     procedure               :: is_exterior               => vef_iterator_is_exterior
     procedure               :: is_interior               => vef_iterator_is_interior
     procedure               :: is_exterior_subvef        => vef_iterator_is_exterior_subvef
     procedure               :: is_interior_subvef        => vef_iterator_is_interior_subvef

     ! h-adaptive FEM
     procedure(is_proper_interface)                       , deferred :: is_proper
     procedure(get_num_improper_cells_around_interface)   , deferred :: get_num_improper_cells_around
     procedure(get_improper_cell_around_interface)        , deferred :: get_improper_cell_around
     procedure(get_improper_cell_around_ivef_interface)   , deferred :: get_improper_cell_around_ivef
     procedure(get_improper_cell_around_subvef_interface) , deferred :: get_improper_cell_around_subvef
     procedure(get_num_half_cells_around_interface)       , deferred :: get_num_half_cells_around
     procedure(get_half_cell_around_interface)            , deferred :: get_half_cell_around
  end type vef_iterator_t

  type, extends(vef_iterator_t) :: itfc_vef_iterator_t
    private
    integer(ip)                        :: itfc_gid = -1
    class(vef_iterator_t), allocatable :: vef
  contains
     procedure          :: create                 => itfc_vef_iterator_create
     procedure          :: free                   => itfc_vef_iterator_free
     procedure          :: first                  => itfc_vef_iterator_first
     procedure          :: next                   => itfc_vef_iterator_next
     procedure          :: has_finished           => itfc_vef_iterator_has_finished
     
     ! Topology-data related TBPs
     procedure          :: get_num_cells_around   => itfc_vef_iterator_get_num_cells_around
     procedure          :: get_cell_around        => itfc_vef_iterator_get_cell_around

     ! Geometry related-data TBPs
     procedure          :: get_num_nodes          => itfc_vef_iterator_get_num_nodes
     procedure          :: get_nodes_coordinates  => itfc_vef_iterator_get_nodes_coordinates

     ! Misc TBPs
     procedure          :: get_ggid               => itfc_vef_iterator_get_ggid
     procedure          :: is_at_interior         => itfc_vef_iterator_is_at_interior
     procedure          :: is_at_boundary         => itfc_vef_iterator_is_at_boundary
     procedure          :: get_dim                => itfc_vef_iterator_get_dim
     procedure          :: get_set_id             => itfc_vef_iterator_get_set_id
     procedure          :: set_set_id             => itfc_vef_iterator_set_set_id
     procedure          :: is_ghost               => itfc_vef_iterator_is_ghost
     procedure          :: is_local               => itfc_vef_iterator_is_local
     procedure          :: is_at_interface        => itfc_vef_iterator_is_at_interface
     procedure          :: is_cut                 => itfc_vef_iterator_is_cut
     
     ! h-adaptive FEM
     procedure          :: is_proper                        => itfc_vef_iterator_is_proper
     procedure          :: get_num_improper_cells_around    => itfc_vef_iterator_get_num_improper_cells_around
     procedure          :: get_improper_cell_around         => itfc_vef_iterator_get_improper_cell_around 
     procedure          :: get_improper_cell_around_ivef    => itfc_vef_iterator_get_improper_cell_around_ivef
     procedure          :: get_improper_cell_around_subvef  => itfc_vef_iterator_get_improper_cell_around_subvef
     procedure          :: get_num_half_cells_around        => itfc_vef_iterator_get_num_half_cells_around
     procedure          :: get_half_cell_around             => itfc_vef_iterator_get_half_cell_around 
  end type itfc_vef_iterator_t

  abstract interface
     function cell_get_ggid_interface ( this )
       import :: cell_iterator_t, igp
       class(cell_iterator_t), intent(in) :: this
       integer(igp) :: cell_get_ggid_interface 
     end function cell_get_ggid_interface 

     function get_my_part_interface ( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip) :: get_my_part_interface
     end function get_my_part_interface

     function cell_is_local_interface ( this )
       import :: cell_iterator_t
       class(cell_iterator_t), intent(in) :: this
       logical :: cell_is_local_interface
     end function cell_is_local_interface

     function cell_is_ghost_interface ( this )
       import :: cell_iterator_t
       class(cell_iterator_t), intent(in) :: this
       logical :: cell_is_ghost_interface
     end function cell_is_ghost_interface
     
     function get_num_vefs_interface ( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip) :: get_num_vefs_interface
     end function get_num_vefs_interface
     
     subroutine get_vef_interface ( this, vef_lid, vef)
       import :: cell_iterator_t, vef_iterator_t, ip
       class(cell_iterator_t), intent(in)    :: this
       integer(ip)           , intent(in)    :: vef_lid
       class(vef_iterator_t) , intent(inout) :: vef
     end subroutine get_vef_interface
     
     function get_vef_gid_interface ( this, vef_lid )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in)    :: this
       integer(ip)           , intent(in)    :: vef_lid
       integer(ip) :: get_vef_gid_interface
     end function get_vef_gid_interface
     
     function get_vef_ggid_interface ( this, vef_lid )
       import :: cell_iterator_t, ip, igp
       class(cell_iterator_t), intent(in)    :: this
       integer(ip)           , intent(in)    :: vef_lid
       integer(igp) :: get_vef_ggid_interface
     end function get_vef_ggid_interface
     
     function get_vefs_gid_interface ( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip), pointer  :: get_vefs_gid_interface(:)
     end function get_vefs_gid_interface
     
     function get_vef_lid_from_gid_interface ( this, vef_gid )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in)    :: this
       integer(ip)           , intent(in)    :: vef_gid
       integer(ip) :: get_vef_lid_from_gid_interface 
     end function get_vef_lid_from_gid_interface
     
     function get_vef_lid_from_ggid_interface ( this, vef_ggid )
       import :: cell_iterator_t, ip, igp
       class(cell_iterator_t), intent(in)    :: this
       integer(igp)          , intent(in)    :: vef_ggid
       integer(ip) :: get_vef_lid_from_ggid_interface 
     end function get_vef_lid_from_ggid_interface
     
     function get_permutation_index_interface(this, target_cell, source_vef_lid, target_vef_lid )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       class(cell_iterator_t), intent(in) :: target_cell
       integer(ip)           , intent(in) :: source_vef_lid
       integer(ip)           , intent(in) :: target_vef_lid
       integer(ip) :: get_permutation_index_interface
     end function get_permutation_index_interface
     
     function get_reference_fe_interface ( this )
       import :: cell_iterator_t, reference_fe_t 
       class(cell_iterator_t), intent(in) :: this 
       class(reference_fe_t), pointer :: get_reference_fe_interface
     end function get_reference_fe_interface 
     
     function get_reference_fe_id_interface ( this )
       import :: cell_iterator_t, ip 
       class(cell_iterator_t), intent(in) :: this 
       integer(ip) :: get_reference_fe_id_interface
     end function get_reference_fe_id_interface 
     
     function cell_get_num_nodes_interface ( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip) :: cell_get_num_nodes_interface
     end function cell_get_num_nodes_interface
     
     subroutine cell_get_nodes_coordinates_interface( this, nodes_coordinates )
       import :: cell_iterator_t, point_t
       class(cell_iterator_t), intent(in)    :: this
       type(point_t)         , intent(inout) :: nodes_coordinates(:)
     end subroutine cell_get_nodes_coordinates_interface
       
     function cell_get_set_id_interface ( this )
       import :: ip, cell_iterator_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip) :: cell_get_set_id_interface
     end function cell_get_set_id_interface 
     
     subroutine cell_set_set_id_interface ( this, set_id )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(inout)  :: this
       integer(ip)           , intent(in)     :: set_id
     end subroutine cell_set_set_id_interface
     
     function cell_get_disconnected_set_id_interface ( this )
       import :: ip, cell_iterator_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip) :: cell_get_disconnected_set_id_interface
     end function cell_get_disconnected_set_id_interface

     subroutine update_sub_triangulation_interface( this )
       import :: cell_iterator_t
       class(cell_iterator_t), intent(inout) :: this
     end subroutine update_sub_triangulation_interface
     
     function get_num_subcells_interface( this ) result ( num_subcells )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip) :: num_subcells
     end function get_num_subcells_interface

     function get_num_subcell_nodes_interface( this ) result ( num_nodes_subcell )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       integer(ip) :: num_nodes_subcell
     end function get_num_subcell_nodes_interface

     subroutine get_phys_coords_of_subcell_interface( this, subcell, points)
       import :: cell_iterator_t, ip, point_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip),                     intent(in)    :: subcell
       type(point_t),                   intent(inout) :: points(:)
     end subroutine get_phys_coords_of_subcell_interface

     subroutine get_ref_coords_of_subcell_interface( this, subcell, points)
       import :: cell_iterator_t, ip, point_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip),                     intent(in)    :: subcell
       type(point_t),                   intent(inout) :: points(:)
     end subroutine get_ref_coords_of_subcell_interface

     function get_num_subfacets_interface( this ) result ( num_subfacets )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in)    :: this
       integer(ip) :: num_subfacets
     end function get_num_subfacets_interface

     function get_num_subfacet_nodes_interface( this ) result ( num_nodes_subfacet )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in)    :: this
       integer(ip) :: num_nodes_subfacet
     end function get_num_subfacet_nodes_interface

     subroutine get_phys_coords_of_subfacet_interface( this, subfacet, points )
       import :: cell_iterator_t, ip, point_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip),                     intent(in)    :: subfacet
       type(point_t),                   intent(inout) :: points(:)
     end subroutine get_phys_coords_of_subfacet_interface

     subroutine get_ref_coords_of_subfacet_interface( this, subfacet, points )
       import :: cell_iterator_t, ip, point_t
       class(cell_iterator_t), intent(in)    :: this
       integer(ip),                     intent(in)    :: subfacet
       type(point_t),                   intent(inout) :: points(:)
     end subroutine get_ref_coords_of_subfacet_interface
     
     function is_cut_interface( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       logical :: is_cut_interface
     end function is_cut_interface
     
     function is_exterior_interface( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       logical :: is_exterior_interface 
     end function is_exterior_interface
     
     function is_interior_interface( this )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in) :: this
       logical :: is_interior_interface 
     end function is_interior_interface
     
     function is_interior_subcell_interface( this, subcell ) result ( is_in )
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(in)  :: this
       integer(ip), intent(in) :: subcell
       logical :: is_in
     end function is_interior_subcell_interface

     function is_exterior_subcell_interface( this, subcell ) result ( is_out )
        import :: cell_iterator_t, ip
        class(cell_iterator_t), intent(in)  :: this
        integer(ip), intent(in) :: subcell
        logical :: is_out
     end function is_exterior_subcell_interface
     
     function get_level_interface(this)
        import :: cell_iterator_t, ip
        class(cell_iterator_t), intent(in)  :: this
        integer(ip) :: get_level_interface
     end function get_level_interface   
     
     subroutine set_for_refinement_interface(this)
       import :: cell_iterator_t
       class(cell_iterator_t), intent(inout)  :: this
     end subroutine set_for_refinement_interface 
     
     subroutine set_for_coarsening_interface(this)
       import :: cell_iterator_t
       class(cell_iterator_t), intent(inout)  :: this
     end subroutine set_for_coarsening_interface
     
     subroutine set_for_do_nothing_interface(this)
       import :: cell_iterator_t
       class(cell_iterator_t), intent(inout)  :: this
     end subroutine set_for_do_nothing_interface
     
     subroutine set_weight_interface(this, weight)
       import :: cell_iterator_t, ip
       class(cell_iterator_t), intent(inout) :: this
       integer(ip)           , intent(in)    :: weight
     end subroutine set_weight_interface
     
  end interface
  
   abstract interface
     function get_num_cells_around_interface ( this )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       integer(ip) :: get_num_cells_around_interface
     end function get_num_cells_around_interface

     subroutine get_cell_around_interface ( this, icell_around, cell )
       import :: vef_iterator_t, cell_iterator_t, ip
       class(vef_iterator_t) , intent(in)    :: this
       integer(ip)           , intent(in)    :: icell_around
       class(cell_iterator_t), intent(inout) :: cell 
     end subroutine get_cell_around_interface

     function get_num_nodes_interface ( this )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       integer(ip) :: get_num_nodes_interface
     end function get_num_nodes_interface

     subroutine get_nodes_coordinates_interface ( this, nodes_coordinates )
       import :: vef_iterator_t, point_t
       class(vef_iterator_t), intent(in)    :: this
       type(point_t)        , intent(inout) :: nodes_coordinates(:)
     end subroutine get_nodes_coordinates_interface
     
     function get_ggid_interface ( this )
       import :: vef_iterator_t, igp
       class(vef_iterator_t), intent(in) :: this
       integer(igp) :: get_ggid_interface
     end function get_ggid_interface
     
     function is_at_interior_interface ( this )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       logical :: is_at_interior_interface 
     end function is_at_interior_interface
     
     function is_at_boundary_interface ( this )
       import :: vef_iterator_t
       class(vef_iterator_t), intent(in) :: this
       logical :: is_at_boundary_interface 
     end function is_at_boundary_interface
     
     function get_dim_interface ( this )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       integer(ip) :: get_dim_interface
     end function get_dim_interface
     
     function get_set_id_interface ( this )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in)    :: this
       integer(ip) :: get_set_id_interface
     end function get_set_id_interface
     
     subroutine set_set_id_interface ( this, set_id )
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(inout) :: this
       integer(ip)          , intent(in)    :: set_id
     end subroutine set_set_id_interface
     
     function is_ghost_interface ( this )
       import :: vef_iterator_t, ip, igp
       class(vef_iterator_t), intent(in)    :: this
       logical :: is_ghost_interface 
     end function is_ghost_interface
     
     function is_local_interface(this)
       import :: vef_iterator_t
       class(vef_iterator_t), intent(in) :: this
       logical :: is_local_interface
     end function is_local_interface
     
     function is_at_interface ( this )
       import :: vef_iterator_t
       class(vef_iterator_t), intent(in) :: this 
       logical :: is_at_interface 
     end function is_at_interface 

     function is_proper_interface ( this )
       import :: vef_iterator_t
       class(vef_iterator_t), intent(in) :: this
       logical :: is_proper_interface
     end function is_proper_interface

     function get_num_improper_cells_around_interface (this)
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       integer(ip) :: get_num_improper_cells_around_interface
     end function get_num_improper_cells_around_interface

     subroutine get_improper_cell_around_interface (this, icell_around, cell)
       import :: vef_iterator_t, cell_iterator_t, ip
       class(vef_iterator_t) , intent(in)    :: this
       integer(ip)           , intent(in)    :: icell_around
       class(cell_iterator_t), intent(inout) :: cell
       integer(ip)                          :: position_in_lst_cells_around
       integer(ip)                          :: icell
     end subroutine get_improper_cell_around_interface

     function get_improper_cell_around_ivef_interface(this, icell_around)
       import :: vef_iterator_t, ip
       class(vef_iterator_t) , intent(in)    :: this
       integer(ip)           , intent(in)    :: icell_around
       integer(ip) :: get_improper_cell_around_ivef_interface
     end function get_improper_cell_around_ivef_interface

     function get_improper_cell_around_subvef_interface(this, icell_around)
       import :: vef_iterator_t, ip
       class(vef_iterator_t) , intent(in)    :: this
       integer(ip)           , intent(in)    :: icell_around
       integer(ip) :: get_improper_cell_around_subvef_interface
     end function get_improper_cell_around_subvef_interface

     function get_num_half_cells_around_interface (this)
       import :: vef_iterator_t, ip
       class(vef_iterator_t), intent(in) :: this
       integer(ip) :: get_num_half_cells_around_interface
     end function get_num_half_cells_around_interface

     subroutine get_half_cell_around_interface (this, icell_around, cell)
       import :: vef_iterator_t, cell_iterator_t, ip
       class(vef_iterator_t) , intent(in)    :: this
       integer(ip)           , intent(in)    :: icell_around
       class(cell_iterator_t), intent(inout) :: cell
     end subroutine get_half_cell_around_interface
  end interface
  
  type object_iterator_t
    private
    integer(ip)                     :: gid = -1
    type(list_iterator_t)           :: vefs_object_iterator
    class(triangulation_t), pointer :: triangulation => NULL()
  contains
    procedure                           :: create                          => object_iterator_create
    procedure                           :: free                            => object_iterator_free
    final                               ::                                    object_iterator_free_final
    procedure, non_overridable, private :: update_vefs_object_iterator     => object_iterator_update_vefs_object_iterator
    procedure                           :: first                           => object_iterator_first
    procedure                           :: next                            => object_iterator_next
    procedure                           :: set_gid                         => object_iterator_set_gid
    procedure                           :: has_finished                    => object_iterator_has_finished
    procedure, non_overridable          :: get_gid                         => object_iterator_get_gid
    procedure, non_overridable          :: get_ggid                        => object_iterator_get_ggid
    procedure, non_overridable          :: get_dim                         => object_iterator_get_dim
    procedure, non_overridable          :: get_num_parts_around            => object_iterator_get_num_parts_around
    procedure, non_overridable          :: get_num_subparts_around         => object_iterator_get_num_subparts_around
    procedure, non_overridable          :: create_parts_around_iterator    => object_iterator_create_parts_around_iterator
    procedure, non_overridable          :: create_subparts_around_iterator => object_iterator_create_subparts_around_iterator
    procedure, non_overridable          :: get_num_vefs                    => object_iterator_get_num_vefs
    procedure, non_overridable          :: get_vef                         => object_iterator_get_vef
  end type object_iterator_t
       
  type, abstract :: triangulation_t 
     private
  
     logical                               :: single_octree_mesh = .false.
  
     ! num of space dimensions
     integer(ip)                           :: num_dims = 0

     ! Data structures to store cell related information
     integer(ip)                           :: num_local_cells = 0
     integer(ip)                           :: num_ghost_cells = 0
     
     ! Data structures to store vef related information
     integer(ip)                           :: num_vefs = 0
     type(std_vector_integer_ip_t)         :: lst_itfc_vefs

     ! Parallel environment describing MPI tasks among which the triangulation is distributed
     type(environment_t), pointer          :: environment => NULL()
     logical                               :: environment_allocated = .false.
     
     ! Data type describing the layout in distributed-memory of the dual graph
     ! (It is required, e.g., for nearest neighbour comms on this graph)
     type(cell_import_t)                   :: cell_import   
     
     ! Data structures to create objects (coarse cell info)
     logical                               :: coarse_triangulation_set_up = .false.
     integer(ip)                           :: num_global_objects = 0
     integer(ip)                           :: num_objects = 0
     integer(igp), allocatable             :: objects_ggids(:)
     integer(ip) , allocatable             :: objects_dim(:)
     type(list_t)                          :: vefs_object
     type(list_t)                          :: faces_object
     type(list_t)                          :: parts_object
     integer(ip)                           :: max_cell_set_id = 0        ! Max cell_set_id among parts   
     integer(ip)                           :: num_subparts = 0           ! num of subparts around part (including those subparts which are local)
     character(len=:), allocatable         :: subparts_coupling_criteria 
     type(list_t)                          :: subparts_object            ! num and list of subparts GIDs around each coarse n_face
     type(hash_table_ip_ip_t)              :: g2l_subparts               ! Translator among the GIDs of subparts and LIDs
     type(coarse_triangulation_t), pointer :: coarse_triangulation => NULL()
     
     ! Scratch data required for non-conforming triangulations, assuming that 
     ! there might be more than one subpart per local subdomain
     type(std_vector_integer_ip_t), allocatable :: lst_subparts_vefwise(:)     
     type(std_vector_integer_ip_t)              :: num_subparts_vefs_cell_wise
     type(std_vector_integer_ip_t)              :: snd_ptrs_lst_subparts_cell_wise
     type(std_vector_integer_ip_t)              :: lst_subparts_pack_idx_cell_wise
     type(std_vector_integer_ip_t)              :: rcv_ptrs_lst_subparts_cell_wise
     type(std_vector_integer_ip_t)              :: lst_subparts_vefs_cell_wise
     type(std_vector_integer_ip_t)              :: rcv_lst_subparts_vefs_cell_wise
     type(std_vector_integer_ip_t)              :: ptrs_to_rcv_lst_subparts_vefs_cell_wise
 contains  
     ! Will the triangulation_t be ALWAYS conforming? (e.g., no matter 
     ! whether it is transformed, refined, coarsened, etc.)
     procedure(is_conforming_interface)         , deferred :: is_conforming
  
     ! Getters
     ! Returns .true. if the triangulation is octree-like, and  it is composed of a 
     ! single octree. Single octree-meshes are such that, for all cells, the mapping 
     ! that transforms among the reference coordinate system and the real cell
     ! coordinate system is composed of translation and/or scalings. Therefore, e.g., 
     ! no rotations are permitted.
     procedure, non_overridable :: is_single_octree_mesh    => triangulation_is_single_octree_mesh
     procedure, non_overridable :: get_num_dims             => triangulation_get_num_dims
     procedure, non_overridable :: get_num_cells            => triangulation_get_num_cells
     procedure, non_overridable :: get_num_local_cells      => triangulation_get_num_local_cells
     procedure, non_overridable :: get_num_ghost_cells      => triangulation_get_num_ghost_cells
     procedure, non_overridable :: get_num_vefs             => triangulation_get_num_vefs
     
     procedure(get_num_proper_vefs_interface)   , deferred :: get_num_proper_vefs
     procedure(get_num_improper_vefs_interface) , deferred :: get_num_improper_vefs
     
     procedure, non_overridable :: set_single_octree_mesh   => triangulation_set_single_octree_mesh
     procedure, non_overridable :: set_num_dims             => triangulation_set_num_dims
     procedure, non_overridable :: set_num_local_cells      => triangulation_set_num_local_cells
     procedure, non_overridable :: set_num_ghost_cells      => triangulation_set_num_ghost_cells
     procedure, non_overridable :: set_num_vefs             => triangulation_set_num_vefs
     
     
     procedure, non_overridable :: get_environment                => triangulation_get_environment
     procedure, non_overridable :: get_cell_import                => triangulation_get_cell_import
     procedure, non_overridable :: get_num_objects                => triangulation_get_num_objects
     procedure, non_overridable :: get_coarse_triangulation       => triangulation_get_coarse_triangulation
     procedure, non_overridable :: get_num_itfc_vefs              => triangulation_get_num_itfc_vefs
     procedure, non_overridable :: coarse_triangulation_is_set_up => triangulation_coarse_triangulation_is_set_up
     
     procedure, non_overridable :: set_environment          => triangulation_set_environment
     !procedure, non_overridable :: allocate_environment     => triangulation_allocate_environment
     !procedure, non_overridable :: free_environment         => triangulation_free_environment
     
     procedure                  :: free                     => triangulation_free
     
     ! Cell traversals-related TBPs
     procedure(create_cell_iterator_interface)  , deferred :: create_cell_iterator
     procedure(free_cell_iterator_interface)    , deferred :: free_cell_iterator
     
     ! Vef traversals-related TBPs
     procedure(create_vef_iterator_interface)   , deferred  :: create_vef_iterator     
     procedure(free_vef_iterator_interface)     , deferred  :: free_vef_iterator
     procedure :: create_itfc_vef_iterator => triangulation_create_itfc_vef_iterator
    
     ! Get num Reference FEs
     procedure(get_num_reference_fes_interface)      , deferred :: get_num_reference_fes
     procedure(tria_get_reference_fe_interface)      , deferred :: get_reference_fe
     procedure(get_max_num_shape_functions_interface), deferred :: get_max_num_shape_functions
 
     ! Objects-related traversals
     procedure, non_overridable :: create_object_iterator  => triangulation_create_object_iterator
     procedure, non_overridable :: free_object_iterator    => triangulation_free_object_iterator
     
     ! Set up TBPs of lst_itfc_vefs
     procedure, non_overridable :: set_up_lst_itfc_vefs    => triangulation_set_up_lst_itfc_vefs
     procedure, non_overridable :: free_lst_itfc_vefs      => triangulation_free_lst_itfc_vefs

     ! Private methods to compute objects
     procedure, non_overridable          :: get_num_subparts                                   => triangulation_get_num_subparts
     procedure, non_overridable          :: get_subpart_lid                                    => triangulation_get_subpart_lid
     procedure, non_overridable, private :: compute_vefs_and_parts_object                      => triangulation_compute_vefs_and_parts_object
     procedure, non_overridable, private :: compute_vefs_and_parts_object_body                 => triangulation_compute_vefs_and_parts_object_body
     procedure, non_overridable, private :: compute_parts_itfc_vefs                            => triangulation_compute_parts_itfc_vefs
     procedure, non_overridable, private :: compute_subparts_itfc_vefs_conforming_mesh         => triangulation_compute_subparts_itfc_vefs_conforming_mesh
     procedure, non_overridable, private :: compute_subparts_itfc_vefs_non_conforming_mesh     => triangulation_compute_subparts_itfc_vefs_non_conforming_mesh
     procedure, non_overridable, private :: compute_parts_object_from_subparts_object          => triangulation_compute_parts_object_from_subparts_object
     procedure, non_overridable, private :: compute_part_id_from_subpart_gid                   => triangulation_compute_part_id_from_subpart_gid
     procedure, non_overridable, private :: compute_objects_dim                                => triangulation_compute_objects_dim
     procedure, non_overridable, private :: compute_objects_neighbours_exchange_data           => triangulation_compute_objects_neighbours_exchange_data
     procedure, non_overridable, private :: compute_num_global_objects_and_their_gids          => triangulation_compute_num_global_objs_and_their_gids
     procedure, non_overridable, private :: free_objects_ggids_and_dim                         => triangulation_free_objects_ggids_and_dim
     
     ! Private methods for coarser triangulation set-up
     procedure, non_overridable          :: setup_coarse_triangulation                 => triangulation_setup_coarse_triangulation
     procedure, non_overridable          :: free_coarse_triangulation_l1_data          => triangulation_free_coarse_triangulation_l1_data
     procedure, non_overridable          :: free_coarse_triangulation_lgt1_data        => triangulation_free_coarse_triangulation_lgt1_data
        
     procedure, non_overridable, private :: gather_coarse_cell_gids                    => triangulation_gather_coarse_cell_gids
     procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs   => triangulation_gather_coarse_vefs_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_vefs_gids                    => triangulation_gather_coarse_vefs_gids
     procedure, non_overridable, private :: gather_coarse_vefs_dim                     => triangulation_gather_coarse_vefs_dim
     procedure, non_overridable, private :: fetch_l2_part_id_neighbours                => triangulation_fetch_l2_part_id_neighbours
     procedure, non_overridable, private :: gather_coarse_dgraph_rcv_counts_and_displs => triangulation_gather_coarse_dgraph_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_dgraph_lextn_and_lextp       => triangulation_gather_coarse_dgraph_lextn_and_lextp
     procedure, non_overridable, private :: adapt_coarse_raw_arrays                    => triangulation_adapt_coarse_raw_arrays
     
     ! Private methods for set up of scratch data related to non-conforming triangulations
     procedure, non_overridable, private :: compute_local_lst_subparts_vefwise              => t_compute_local_lst_subparts_vefwise
     procedure, non_overridable, private :: exchange_lst_subparts_round                     => t_exchange_lst_subparts_round
     procedure, non_overridable, private :: fetch_num_subparts_vefs_cell_wise               => t_fetch_num_subparts_vefs_cell_wise
     procedure, non_overridable, private :: compute_near_neigh_ctrl_data_lst_subparts       => t_compute_near_neigh_ctrl_data_lst_subparts
     procedure, non_overridable, private :: fetch_lst_subparts_vefs_cell_wise               => t_fetch_lst_subparts_vefs_cell_wise
     procedure, non_overridable, private :: compute_ptrs_to_rcv_lst_subparts_vefs_cell_wise => t_compute_ptrs_to_rcv_lst_subparts_vefs_cell_wise
     procedure, non_overridable, private :: update_lst_subparts_vefwise_after_exchange      => t_update_lst_subparts_vefwise_after_exchange
     procedure, non_overridable, private :: free_non_conforming_scratch_data                => triangulation_free_non_conforming_scratch_data
     
     ! Methods for disconnected parts identification 
     procedure                           :: get_max_cell_set_id                         => triangulation_get_max_cell_set_id
     procedure, non_overridable          :: set_subparts_coupling_criteria              => triangulation_set_subparts_coupling_criteria 
     procedure, non_overridable, private :: allocate_and_fill_disconnected_cells_set_id => triangulation_allocate_and_fill_disconnected_cells_set_id 
     procedure, non_overridable, private :: compute_disconnected_cells_set_id           => triangulation_compute_disconnected_cells_set_id
     procedure, non_overridable, private :: generate_dual_graph                         => triangulation_generate_dual_graph  
     procedure(compute_max_cells_set_id_interface)     , deferred :: compute_max_cells_set_id 
     procedure(resize_disconnected_cells_set_interface), deferred :: resize_disconnected_cells_set
     procedure(fill_disconnected_cells_set_interface)  , deferred :: fill_disconnected_cells_set

  end type triangulation_t
  
  abstract interface
     function is_conforming_interface ( this )
       import :: triangulation_t
       class(triangulation_t) , intent(in) :: this
       logical :: is_conforming_interface 
     end function is_conforming_interface 
     
     function get_num_proper_vefs_interface ( this )
       import :: triangulation_t, ip
       class(triangulation_t) , intent(in) :: this
       integer(ip) :: get_num_proper_vefs_interface
     end function get_num_proper_vefs_interface

     function get_num_improper_vefs_interface ( this )
       import :: triangulation_t, ip
       class(triangulation_t) , intent(in) :: this
       integer(ip) :: get_num_improper_vefs_interface
     end function get_num_improper_vefs_interface

     subroutine create_cell_iterator_interface ( this, cell )
       import :: triangulation_t, cell_iterator_t
       class(triangulation_t) , intent(in) :: this
       class(cell_iterator_t), allocatable, intent(inout) :: cell
     end subroutine create_cell_iterator_interface

     subroutine free_cell_iterator_interface ( this, cell )
       import :: triangulation_t, cell_iterator_t
       class(triangulation_t) , intent(in) :: this
       class(cell_iterator_t), allocatable, intent(inout) :: cell
     end subroutine free_cell_iterator_interface

     subroutine create_vef_iterator_interface ( this, vef )
       import :: triangulation_t, vef_iterator_t
       class(triangulation_t) , intent(in) :: this
       class(vef_iterator_t) , allocatable,  intent(inout) :: vef
     end subroutine create_vef_iterator_interface

     subroutine free_vef_iterator_interface ( this, vef )
       import :: triangulation_t, vef_iterator_t
       class(triangulation_t) , intent(in) :: this
       class(vef_iterator_t), allocatable, intent(inout) :: vef
     end subroutine free_vef_iterator_interface

     function get_num_reference_fes_interface ( this )
       import :: triangulation_t, ip
       class(triangulation_t), intent(in) :: this
       integer(ip) :: get_num_reference_fes_interface
     end function get_num_reference_fes_interface
     
     function tria_get_reference_fe_interface ( this, ref_fe_geo_id )
       import :: triangulation_t, reference_fe_t, ip
       class(triangulation_t), target, intent(in) :: this 
       integer(ip)                   , intent(in) :: ref_fe_geo_id
       class(reference_fe_t), pointer :: tria_get_reference_fe_interface
     end function tria_get_reference_fe_interface 

     function compute_max_cells_set_id_interface ( this ) 
       import :: triangulation_t, ip
       class(triangulation_t), intent(inout)   :: this
       integer(ip) :: compute_max_cells_set_id_interface
     end function compute_max_cells_set_id_interface
     
     subroutine resize_disconnected_cells_set_interface ( this ) 
       import :: triangulation_t
       class(triangulation_t), intent(inout)   :: this
     end subroutine resize_disconnected_cells_set_interface 
     
     subroutine fill_disconnected_cells_set_interface ( this, disconnected_cells_set ) 
       import :: triangulation_t, ip 
       class(triangulation_t), intent(inout)   :: this
       integer(ip)           , intent(in)      :: disconnected_cells_set(:)
     end subroutine fill_disconnected_cells_set_interface 

     function get_max_num_shape_functions_interface ( this )
       import :: triangulation_t, ip
       class(triangulation_t), intent(in) :: this
       integer(ip) :: get_max_num_shape_functions_interface
     end function get_max_num_shape_functions_interface

  end interface
  
  ! Parameters to define vef_type. Observe that
  ! mod(vef_type,10)     = dimension
  ! mod(vef_type/10,10)  = interior (0) or boundary (1)
  ! vef_type/100         = local (0), interface(1) or ghost(2)
  integer(ip), parameter :: local_interior_dim0 =  0
  integer(ip), parameter :: local_interior_dim1 =  1
  integer(ip), parameter :: local_interior_dim2 =  2
  integer(ip), parameter :: local_boundary_dim0 = 10
  integer(ip), parameter :: local_boundary_dim1 = 11
  integer(ip), parameter :: local_boundary_dim2 = 12

  integer(ip), parameter :: interface_interior_dim0 = 100
  integer(ip), parameter :: interface_interior_dim1 = 101
  integer(ip), parameter :: interface_interior_dim2 = 102
  integer(ip), parameter :: interface_boundary_dim0 = 110
  integer(ip), parameter :: interface_boundary_dim1 = 111
  integer(ip), parameter :: interface_boundary_dim2 = 112

  integer(ip), parameter :: ghost_dim0 = 200
  integer(ip), parameter :: ghost_dim1 = 201
  integer(ip), parameter :: ghost_dim2 = 202
  
  integer(ip), parameter :: triangulation_generate_from_mesh  = 0
  integer(ip), parameter :: triangulation_generate_structured = 1
  public :: triangulation_generate_from_mesh
  public :: triangulation_generate_structured
  
  character(len=*), parameter :: geometry_interpolation_order_key    = 'geometry_interpolation_order'
  character(len=*), parameter :: triangulation_generate_key          = 'triangulation_generate'
  public :: geometry_interpolation_order_key
  public :: triangulation_generate_key
  
  character(len=*), parameter :: subparts_coupling_criteria_key = 'subparts_coupling_criteria'
  character(len=*), parameter :: all_coupled                    = 'all_coupled'
  character(len=*), parameter :: loose_coupling                 = 'loose_coupling' 
  character(len=*), parameter :: strong_coupling                = 'strong_coupling' 
  public :: subparts_coupling_criteria_key 
  public :: loose_coupling
  public :: strong_coupling 
  
  type, extends(cell_iterator_t) :: bst_cell_iterator_t
    private
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
    procedure                            :: create                  => bst_cell_iterator_create
    procedure                            :: free                    => bst_cell_iterator_free
    final                                ::                            bst_cell_iterator_free_final
    procedure, non_overridable, private  :: set_ggid                => bst_cell_iterator_set_ggid
    procedure, non_overridable, private  :: set_mypart              => bst_cell_iterator_set_mypart
    procedure                            :: get_reference_fe        => bst_cell_iterator_get_reference_fe
    procedure                            :: get_reference_fe_id     => bst_cell_iterator_get_reference_fe_id
    procedure                            :: get_nodes_coordinates   => bst_cell_iterator_get_nodes_coordinates
    procedure, non_overridable           :: set_nodes_coordinates   => bst_cell_iterator_set_nodes_coordinates
    procedure                            :: get_nodes_coordinates_ref_space => bst_cell_iterator_get_nodes_coordinates_ref_space
    procedure                            :: get_ggid                 => bst_cell_iterator_get_ggid
    procedure                            :: get_my_part             => bst_cell_iterator_get_mypart
    procedure                            :: set_set_id              => bst_cell_iterator_set_set_id
    procedure                            :: get_set_id              => bst_cell_iterator_get_set_id
    procedure                            :: get_disconnected_set_id => bst_cell_iterator_get_disconnected_set_id
    procedure                            :: get_num_vefs            => bst_cell_iterator_get_num_vefs
    procedure                            :: get_num_nodes           => bst_cell_iterator_get_num_nodes
    procedure                            :: get_node_gid            => bst_cell_iterator_get_node_gid
    procedure                            :: get_vef_gid             => bst_cell_iterator_get_vef_gid
    procedure                            :: get_vef_ggid            => bst_cell_iterator_get_vef_ggid
    procedure                            :: get_vefs_gid            => bst_cell_iterator_get_vefs_gid
    procedure                            :: get_vef_lid_from_gid    => bst_cell_iterator_get_vef_lid_from_gid
    procedure                            :: get_vef_lid_from_ggid   => bst_cell_iterator_get_vef_lid_from_ggid
    procedure                            :: get_vef                 => bst_cell_iterator_get_vef
    procedure                            :: is_local                => bst_cell_iterator_is_local
    procedure                            :: is_ghost                => bst_cell_iterator_is_ghost
    procedure                            :: get_permutation_index   => bst_cell_iterator_get_permutation_index

    ! Declare dummy procedures to be implemented in the corresponding derived classes 
    procedure :: update_sub_triangulation    => bst_cell_iterator_update_sub_triangulation
    procedure :: get_num_subcells      => bst_cell_iterator_get_num_subcells
    procedure :: get_num_subcell_nodes => bst_cell_iterator_get_num_subcell_nodes
    procedure :: get_phys_coords_of_subcell  => bst_cell_iterator_get_phys_coords_of_subcell
    procedure :: get_ref_coords_of_subcell   => bst_cell_iterator_get_ref_coords_of_subcell
    procedure :: get_num_subfacets      => bst_cell_iterator_get_num_subfacets
    procedure :: get_num_subfacet_nodes => bst_cell_iterator_get_num_subfacet_nodes
    procedure :: get_phys_coords_of_subfacet  => bst_cell_iterator_get_phys_coords_of_subfacet
    procedure :: get_ref_coords_of_subfacet   => bst_cell_iterator_get_ref_coords_of_subfacet
    procedure :: is_cut                      => bst_cell_iterator_is_cut
    procedure :: is_interior                 => bst_cell_iterator_is_interior
    procedure :: is_exterior                 => bst_cell_iterator_is_exterior
    procedure :: is_interior_subcell         => bst_cell_iterator_is_interior_subcell
    procedure :: is_exterior_subcell         => bst_cell_iterator_is_exterior_subcell
    
    procedure :: get_level                   => bst_cell_iterator_get_level
    procedure :: set_for_refinement          => bst_cell_iterator_set_for_refinement
    procedure :: set_for_coarsening          => bst_cell_iterator_set_for_coarsening
    procedure :: set_for_do_nothing          => bst_cell_iterator_set_for_do_nothing
    procedure :: set_weight                  => bst_cell_iterator_set_weight
    
    procedure, non_overridable, private  :: fill_nodes_on_vertices        => bst_cell_iterator_fill_nodes_on_vertices
    procedure, non_overridable, private  :: fill_nodes_on_vef_new         => bst_cell_iterator_fill_nodes_on_vef_new
    procedure, non_overridable, private  :: fill_nodes_on_vef_from_source => bst_cell_iterator_fill_nodes_on_vef_from_source
    procedure, non_overridable, private  :: fill_internal_nodes_new       => bst_cell_iterator_fill_internal_nodes_new
    
  end type bst_cell_iterator_t
  
  type, extends(vef_iterator_t) :: bst_vef_iterator_t
    private
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
     procedure                           :: create                    => bst_vef_iterator_create
     procedure                           :: free                      => bst_vef_iterator_free
     final                               ::                              bst_vef_iterator_free_final
     procedure                           :: get_num_nodes             => bst_vef_iterator_get_num_nodes
     procedure                           :: get_nodes_coordinates     => bst_vef_iterator_get_nodes_coordinates
     procedure                           :: get_ggid                  => bst_vef_iterator_get_ggid

     procedure, non_overridable          :: set_geom_id               => bst_vef_iterator_set_geom_id
     procedure                           :: set_set_id                => bst_vef_iterator_set_set_id
     procedure, non_overridable          :: get_geom_id               => bst_vef_iterator_get_geom_id
     procedure                           :: get_set_id                => bst_vef_iterator_get_set_id

     procedure, non_overridable          :: set_dim                   => bst_vef_iterator_set_dim
     procedure, non_overridable          :: set_it_at_boundary        => bst_vef_iterator_set_it_at_boundary
     procedure, non_overridable          :: set_it_as_local           => bst_vef_iterator_set_it_as_local
     procedure, non_overridable          :: set_it_as_ghost           => bst_vef_iterator_set_it_as_ghost
     procedure, non_overridable          :: set_it_at_interface       => bst_vef_iterator_set_it_at_interface

     procedure                           :: get_dim                   => bst_vef_iterator_get_dim
     procedure                           :: is_at_interior            => bst_vef_iterator_is_at_interior
     procedure                           :: is_at_boundary            => bst_vef_iterator_is_at_boundary
     procedure                           :: is_local                  => bst_vef_iterator_is_local
     procedure                           :: is_ghost                  => bst_vef_iterator_is_ghost
     procedure                           :: is_at_interface           => bst_vef_iterator_is_at_interface
     procedure                           :: is_cut                    => bst_vef_iterator_is_cut
     
     procedure                           :: get_num_cells_around      => bst_vef_iterator_get_num_cells_around
     procedure                           :: get_cell_around           => bst_vef_iterator_get_cell_around

     procedure                           :: is_proper                       => bst_vef_iterator_is_proper
     procedure                           :: get_num_improper_cells_around   => bst_vef_iterator_get_num_improper_cells_around
     procedure                           :: get_improper_cell_around        => bst_vef_iterator_get_improper_cell_around
     procedure                           :: get_improper_cell_around_ivef   => bst_vef_iterator_get_improper_cell_around_ivef
     procedure                           :: get_improper_cell_around_subvef => bst_vef_iterator_get_improper_cell_around_subvef
     procedure                           :: get_num_half_cells_around       => bst_vef_iterator_get_num_half_cells_around
     procedure                           :: get_half_cell_around            => bst_vef_iterator_get_half_cell_around
  end type bst_vef_iterator_t
      
  integer(ip), parameter      :: max_num_reference_fes_geo = 3

  type, extends(triangulation_t) :: base_static_triangulation_t ! Base class for coarse_triangulation_t and fine_triangulation_t
     private
    
     integer(ip)                           :: num_vertices
     
     integer(igp), allocatable             :: cells_ggid(:)                 ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_mypart(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_set(:)                  ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: disconnected_cells_set_ids(:) ! Num local cells + num ghost cells 
     integer(ip) , allocatable             :: ptr_vefs_x_cell(:)            ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable             :: lst_vefs_gids(:)
 
     
     ! Data structures to store vef related information
     integer(igp), allocatable             :: vefs_ggid(:)         ! num_vefs
     integer(ip) , allocatable             :: vefs_set(:)          ! num_vefs
     integer(ip) , allocatable             :: vefs_geometry(:)     ! num_vefs
     integer(ip) , allocatable             :: vefs_type(:)         ! num_vefs
     integer(ip) , allocatable             :: ptrs_cells_around(:) ! num_itfc_vefs+1
     integer(ip) , allocatable             :: lst_cells_around(:)  ! ptrs_cells_around(num_itfc_vefs+1)-1

     ! Data structures that should be defined in fine_triangulation_t (which requires extensive refactoring)     
     type(geometry_t)                      :: geometry
     type(p_lagrangian_reference_fe_t)     :: reference_fe_geo_list(max_num_reference_fes_geo)
     type(hash_table_ip_ip_t)              :: reference_fe_geo_index
     
     ! Geometry interpolation
     integer(ip)                           :: num_nodes
     integer(ip)  , allocatable            :: ptr_nodes_x_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip)  , allocatable            :: lst_nodes(:)
     type(point_t), allocatable            :: coordinates(:)
 contains  
     procedure                           :: get_num_proper_vefs              => bst_get_num_proper_vefs
     procedure                           :: get_num_improper_vefs            => bst_get_num_improper_vefs
     
     ! Cell traversals-related TBPs
     procedure                           :: create_cell_iterator             => bst_create_cell_iterator
     procedure                           :: free_cell_iterator               => bst_free_cell_iterator
     
     ! Vef traversals-related TBPs
     procedure                           :: create_vef_iterator              => bst_create_vef_iterator
     procedure                           :: free_vef_iterator                => bst_free_vef_iterator
     
     ! Other
     procedure                           :: print                            => bst_print    
     procedure                           :: free                             => bst_free
     
     ! Getters
     procedure                           :: get_num_reference_fes            => bst_get_num_reference_fes
     procedure                           :: get_reference_fe                 => bst_get_reference_fe
     procedure                           :: get_max_num_shape_functions      => bst_get_max_num_shape_functions
     procedure                           :: is_tet_mesh                      => bst_is_tet_mesh
     procedure                           :: is_hex_mesh                      => bst_is_hex_mesh
     procedure                           :: is_mix_mesh                      => bst_is_mix_mesh
     
     procedure                           :: is_conforming                    => bst_is_conforming
   
     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_ptr_vefs_x_cell   => bst_allocate_and_fill_ptr_vefs_x_cell
     procedure, non_overridable, private :: allocate_cells_ggid                 => bst_allocate_cells_ggid
     procedure, non_overridable, private :: fill_local_cells_ggid               => bst_fill_local_cells_ggid
     procedure, non_overridable, private :: allocate_cells_mypart               => bst_allocate_cells_mypart
     procedure, non_overridable, private :: fill_local_cells_mypart             => bst_fill_local_cells_mypart
     procedure, non_overridable, private :: allocate_cells_set                  => bst_allocate_cells_set
     procedure, non_overridable, private :: allocate_disconnected_cells_set     => bst_allocate_disconnected_cells_set
     procedure, non_overridable          :: fill_cells_set                      => bst_fill_cells_set
     procedure                           :: compute_max_cells_set_id            => bst_compute_max_cells_set_id
     procedure                           :: resize_disconnected_cells_set       => bst_resize_disconnected_cells_set
     procedure                           :: fill_disconnected_cells_set         => bst_fill_disconnected_cells_set
     procedure, non_overridable, private :: free_ptr_vefs_x_cell                => bst_free_ptr_vefs_x_cell
     procedure, non_overridable, private :: free_lst_vefs_gids                  => bst_free_lst_vefs_gids 
     procedure, non_overridable, private :: free_cells_ggid                     => bst_free_cells_ggid
     procedure, non_overridable, private :: free_cells_mypart                   => bst_free_cells_mypart
     procedure, non_overridable, private :: free_cells_set                      => bst_free_cells_set 
     procedure, non_overridable, private :: free_disconnected_cells_set         => bst_free_disconnected_cells_set
     procedure, non_overridable, private :: orient_tet_mesh                     => bst_orient_tet_mesh

     ! Private methods to perform nearest neighbor exchange
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_ggids
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_ggids_and_dim
     procedure, non_overridable, nopass, private :: bst_cell_pack_vef_ggids_and_coordinates
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_ggids
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_ggids_and_dim
     procedure, non_overridable, nopass, private :: bst_cell_unpack_vef_ggids_and_coordinates
     procedure, non_overridable, private :: fetch_ghost_cells_data         => bst_fetch_ghost_cells_data
     procedure, non_overridable, nopass, private :: cell_size              => bst_cell_size
     generic,                            private :: cell_pack              => bst_cell_pack_vef_ggids, &
                                                                              bst_cell_pack_vef_ggids_and_dim, &
                                                                              bst_cell_pack_vef_ggids_and_coordinates
     generic,                            private :: cell_unpack            => bst_cell_unpack_vef_ggids, &
                                                                              bst_cell_unpack_vef_ggids_and_dim, &
                                                                              bst_cell_unpack_vef_ggids_and_coordinates

     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: compute_num_vefs                    => bst_compute_num_vefs
     procedure, non_overridable, private :: allocate_and_fill_vefs_ggid         => bst_allocate_and_fill_vefs_ggid
     procedure, non_overridable, private :: free_vefs_ggid                      => bst_free_vefs_ggid
     procedure, non_overridable, private :: free_vefs_type                      => bst_triangulation_free_vefs_type
     procedure, non_overridable, private :: allocate_and_fill_cells_around      => bst_allocate_and_fill_cells_around
     procedure, non_overridable, private :: free_cells_around                   => bst_free_cells_around
     procedure, non_overridable, private :: find_local_ghost_vefs               => bst_find_local_ghost_vefs
  end type base_static_triangulation_t
  
  type, extends(base_static_triangulation_t) :: fine_triangulation_t
     private
  contains
     ! Private methods to create vefs (these functions make use of the reference fe and therefore are not bounded
     ! to the mother class)
     procedure, non_overridable, private :: fill_reference_fe_geo_list          => fine_triangulation_fill_reference_fe_geo_list
     procedure, non_overridable, private :: generate_vefs                       => fine_triangulation_generate_vefs
     procedure, non_overridable, private :: allocate_and_fill_geometry_and_set  => fine_triangulation_allocate_and_fill_geometry_and_set
     procedure, non_overridable, private :: free_geometry_and_set               => fine_triangulation_free_geometry_and_set
     procedure, non_overridable, private :: compute_vefs_dim                    => fine_triangulation_compute_vefs_dim
     procedure, non_overridable, private :: find_vefs_at_boundary               => fine_triangulation_find_vefs_at_boundary
     ! Geometry interpolation 
     procedure, non_overridable, private :: allocate_and_fill_nodes             => fine_triangulation_allocate_and_fill_nodes
     procedure, non_overridable, private :: free_nodes                          => fine_triangulation_free_nodes
     procedure, non_overridable          :: allocate_and_fill_coordinates       => fine_triangulation_allocate_and_fill_coordinates
     procedure, non_overridable          :: free_coordinates                    => fine_triangulation_free_coordinates
     procedure                           :: free                                => fine_triangulation_free     
  end type fine_triangulation_t

  type, extends(fine_triangulation_t) :: serial_triangulation_t
  contains
     procedure                 , private :: serial_triangulation_create
     generic                             :: create                              => serial_triangulation_create
  end type serial_triangulation_t
  
  type, extends(fine_triangulation_t) :: par_triangulation_t
  contains
     generic                             :: create                              => par_triangulation_create
     procedure, private                  :: par_triangulation_create
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_gids     => par_triangulation_allocate_and_fill_lst_vefs_gids 
  end type par_triangulation_t
  
  type, extends(base_static_triangulation_t) :: coarse_triangulation_t
  contains
     procedure, non_overridable          :: create                              => coarse_triangulation_create
     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_gids     => coarse_triangulation_allocate_and_fill_lst_vefs_gids 
     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: allocate_and_fill_vefs_dim          => coarse_triangulation_allocate_and_fill_vefs_dim
     procedure                           :: print                               => coarse_triangulation_print
     procedure                           :: get_ptr_vefs_x_cell                 => coarse_triangulation_get_ptr_vefs_x_cell
  end type coarse_triangulation_t
  
  public :: triangulation_t, serial_triangulation_t, par_triangulation_t, coarse_triangulation_t
  public :: triangulation_free
  public :: cell_iterator_t
  public :: vef_iterator_t
  public :: itfc_vef_iterator_t, object_iterator_t
  public :: cell_iterator_create
  public :: cell_iterator_free
  public :: cell_iterator_next
  public :: cell_iterator_first
  public :: cell_iterator_set_gid
  public :: vef_iterator_create
  public :: vef_iterator_free
  
  public :: base_static_triangulation_t
  public :: bst_cell_iterator_t, bst_vef_iterator_t
  
  
contains

#include "sbm_cell_iterator.i90"
#include "sbm_vef_iterator.i90"
#include "sbm_itfc_vef_iterator.i90"
#include "sbm_object_iterator.i90"
#include "sbm_triangulation.i90"

#include "sbm_bst_cell_iterator.i90"
#include "sbm_bst_vef_iterator.i90"
#include "sbm_base_static_triangulation.i90"
#include "sbm_fine_triangulation.i90"
#include "sbm_serial_triangulation.i90"
#include "sbm_par_triangulation.i90"
#include "sbm_coarse_triangulation.i90"

end module triangulation_names
