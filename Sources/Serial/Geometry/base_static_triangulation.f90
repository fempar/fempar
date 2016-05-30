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
module base_static_triangulation_names
  ! Serial modules
  use types_names
  use memor_names
  use sort_names
  use reference_fe_names
  use reference_fe_factory_names
  use element_import_names
  use hash_table_names
  use list_types_names
  use mesh_names

  ! Par modules
  use par_environment_names

  implicit none
# include "debug.i90"
  private
  
  type cell_accessor_t
    private
    integer(ip)                                 :: lid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation
  contains
    procedure, non_overridable, private  :: create               => cell_accessor_create
    procedure, non_overridable, private  :: free                 => cell_accessor_free
    procedure, non_overridable, private  :: next                 => cell_accessor_next
    procedure, non_overridable, private  :: set_lid              => cell_accessor_set_lid
    procedure, non_overridable, private  :: set_gid              => cell_accessor_set_gid
    procedure, non_overridable, private  :: set_mypart           => cell_accessor_set_mypart
    procedure, non_overridable, private  :: get_triangulation    => cell_accessor_get_triangulation
    procedure, non_overridable           :: past_the_end         => cell_accessor_past_the_end
    procedure, non_overridable           :: get_reference_fe_geo => cell_accessor_get_reference_fe_geo
    procedure, non_overridable           :: get_lid              => cell_accessor_get_lid
    procedure, non_overridable           :: get_gid              => cell_accessor_get_gid
    procedure, non_overridable           :: get_mypart           => cell_accessor_get_mypart
    procedure, non_overridable           :: get_num_vefs         => cell_accessor_get_num_vefs
    procedure, non_overridable           :: get_vef_lid          => cell_accessor_get_vef_lid
    procedure, non_overridable           :: get_vef_gid          => cell_accessor_get_vef_gid
    procedure, non_overridable           :: get_vef              => cell_accessor_get_vef
    procedure, non_overridable           :: is_local             => cell_accessor_is_local
    procedure, non_overridable           :: is_ghost             => cell_accessor_is_ghost
  end type cell_accessor_t
  
  type cell_iterator_t
    private
    type(cell_accessor_t) :: current_cell_accessor
  contains
     procedure, non_overridable, private :: create       => cell_iterator_create
     procedure, non_overridable          :: free         => cell_iterator_free
     procedure, non_overridable          :: init         => cell_iterator_init
     procedure, non_overridable          :: next         => cell_iterator_next
     procedure, non_overridable          :: has_finished => cell_iterator_has_finished
     procedure, non_overridable          :: current      => cell_iterator_current
  end type cell_iterator_t  
  
  type vef_accessor_t
    private
    integer(ip)                                 :: lid = -1
    class(base_static_triangulation_t), pointer :: base_static_triangulation => NULL()
  contains
     procedure, non_overridable, private :: create                   => vef_accessor_create
     procedure, non_overridable, private :: free                     => vef_accessor_free
     procedure, non_overridable, private :: next                     => vef_accessor_next
     procedure, non_overridable, private :: set_lid                  => vef_accessor_set_lid
     procedure, non_overridable, private :: past_the_end             => vef_accessor_past_the_end
     procedure, non_overridable, private :: get_triangulation        => vef_accessor_get_triangulation
     procedure, non_overridable          :: get_lid                  => vef_accessor_get_lid
     procedure, non_overridable          :: get_gid                  => vef_accessor_get_gid
     procedure, non_overridable          :: is_local                 => vef_accessor_is_local
     procedure, non_overridable          :: is_ghost                 => vef_accessor_is_ghost
     procedure, non_overridable          :: at_interface             => vef_accessor_at_interface
     procedure, non_overridable          :: get_dimension            => vef_accessor_get_dimension
     procedure, non_overridable          :: get_num_cells_around     => vef_accessor_get_num_cells_around
     procedure, non_overridable          :: get_cell_around          => vef_accessor_get_cell_around
  end type vef_accessor_t
  
  type vef_iterator_t
    private
    type(vef_accessor_t) :: current_vef_accessor
  contains
     procedure, non_overridable, private :: create       => vef_iterator_create
     procedure, non_overridable          :: free         => vef_iterator_free
     procedure, non_overridable          :: init         => vef_iterator_init
     procedure, non_overridable          :: next         => vef_iterator_next
     procedure, non_overridable          :: has_finished => vef_iterator_has_finished
     procedure, non_overridable          :: current      => vef_iterator_current
  end type vef_iterator_t  
  
  type :: itfc_vef_iterator_t
    private
    integer(ip)          :: itfc_lid = -1
    type(vef_accessor_t) :: current_vef_accessor
  contains
    procedure, non_overridable, private :: create       => itfc_vef_iterator_create
    procedure, non_overridable          :: free         => itfc_vef_iterator_free
    procedure, non_overridable          :: init         => itfc_vef_iterator_init
    procedure, non_overridable          :: next         => itfc_vef_iterator_next
    procedure, non_overridable          :: has_finished => itfc_vef_iterator_has_finished
    procedure, non_overridable          :: current      => itfc_vef_iterator_current
  end type itfc_vef_iterator_t
  
   ! JP-TODO: implement states: discuss with Alberto and Victor.
   !
   ! State transition diagram for type(base_static_triangulation_t). The
   ! creation of reference elements must be performed together
   ! with reading from mesh when going from created to elements_filled.
   !
   ! ------------------------------------------------...--------
   ! Input State      | Action                    | Output State 
   ! -----------------------------------------------------------
   ! not_created      | create                    | created 
   ! not_created      | free                      | not_created
   ! created          | free                      | not_created
   ! created          | external (read from mesh) | elements_filled
   ! elements_filled  | create_dual               | vefs_filled
   ! elements_filled  | free                      | not_created
   ! vefs_filled      | free                      | not_created

   !integer(ip), parameter :: base_static_triangulation_not_created     = 0 
   !integer(ip), parameter :: base_static_triangulation_created         = 1 
   !integer(ip), parameter :: base_static_triangulation_elements_filled = 2 
   !integer(ip), parameter :: base_static_triangulation_vefs_filled     = 3 

  integer(ip), parameter      :: max_num_elem_types = 3
  
  type base_static_triangulation_t ! Base class for serial_triangulation_t and par_base_static_triangulation_t
     private
     integer(ip)                           :: num_dimensions  = -1
     integer(ip)                           :: num_local_cells = -1
     integer(ip)                           :: num_ghost_cells = -1
     integer(ip)                           :: max_vefs_per_cell = -1
     
     integer(igp), allocatable             :: cells_gid(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: cells_mypart(:)            ! Num local cells + num ghost cells
     
     type(p_reference_fe_t)                :: reference_fe_geo_list(max_num_elem_types)
     type(position_hash_table_t)           :: reference_fe_geo_index(max_num_elem_types)
     ! The reference fe for the geometry of each element need not to be stored as it can
     ! be recovered from the number of vefs
     integer(ip) , allocatable             :: elems_reference_fe_geo(:)  ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: ptr_vefs_per_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable             :: lst_vefs_lids(:)
          
     integer(ip)                           :: num_vefs = -1        ! = num_local_vefs+num_ghost_vefs
     integer(ip)                           :: num_local_vefs = -1
     integer(ip)                           :: num_ghost_vefs = -1
     integer(igp), allocatable             :: vefs_gid(:)          ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_set(:)          ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_geometry(:)     ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_dimension(:)    ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_itfc_lid(:)     ! num_local_vefs + num_ghost_vefs

     integer(ip)                           :: num_itfc_vefs  = -1
     integer(ip), allocatable              :: lst_itfc_vefs(:)
     integer(ip), allocatable              :: ptrs_cells_around(:) ! num_itfc_vefs+1
     integer(ip), allocatable              :: lst_cells_around(:)  ! ptrs_cells_around(num_itfc_vefs+1)-1
     
     integer(ip)                           :: num_vertices
     real(rp)   , allocatable              :: coordinates(:,:)
  contains  
  
     ! Private methods for creating vef-related data
     procedure, non_overridable, private :: free_ptr_vefs_per_cell             => base_static_triangulation_free_ptr_vefs_per_cell
     procedure, non_overridable, private :: free_lst_vefs_lids                 => base_static_triangulation_free_lst_vefs_lids 

     procedure, non_overridable, private :: compute_num_local_vefs             => base_static_triangulation_compute_num_local_vefs
     procedure, non_overridable, private :: compute_num_ghost_vefs             => base_static_triangulation_compute_num_ghost_vefs
     procedure, non_overridable, private :: allocate_and_fill_vefs_gid         => base_static_triangulation_allocate_and_fill_vefs_gid
     procedure, non_overridable, private :: free_vefs_gid                      => base_static_triangulation_free_vefs_gid
     procedure, non_overridable, private :: allocate_and_fill_vefs_dimension   => base_static_triangulation_allocate_and_fill_vefs_dimension
     procedure, non_overridable, private :: free_vefs_dimension                => base_static_triangulation_free_vefs_dimension
     
     procedure, non_overridable, private :: compute_num_itfc_vefs              => base_static_triangulation_compute_num_itfc_vefs
     procedure, non_overridable, private :: allocate_and_fill_lst_itfc_vefs    => base_static_triangulation_allocate_and_fill_lst_itfc_vefs
     procedure, non_overridable, private :: free_lst_itfc_vefs                 => base_static_triangulation_free_lst_itfc_vefs
     procedure, non_overridable, private :: allocate_and_fill_vefs_itfc_lid    => base_static_triangulation_allocate_and_fill_vefs_itfc_lid
     procedure, non_overridable, private :: free_vefs_itfc_lid                 => base_static_triangulation_free_vefs_itfc_lid
     procedure, non_overridable, private :: allocate_and_fill_cells_around     => base_static_triangulation_allocate_and_fill_cells_around
     procedure, non_overridable, private :: free_cells_around                  => base_static_triangulation_free_cells_around

     procedure, non_overridable          :: generate_vefs               => base_static_triangulation_generate_vefs

     ! Cell traversals-related TBPs
     procedure, non_overridable            :: create_cell_iterator      => base_static_triangulation_create_cell_iterator
  
     ! Vef traversals-related TBPs
     procedure, non_overridable          :: create_vef_iterator                => base_static_triangulation_create_vef_iterator
     procedure, non_overridable          :: create_itfc_vef_iterator           => base_static_triangulation_create_itfc_vef_iterator
  end type base_static_triangulation_t
  
  type, extends(base_static_triangulation_t) :: par_base_static_triangulation_t
     private
     ! Parallel environment describing MPI tasks among which the triangulation is distributed
     type(par_environment_t),   pointer      :: p_env => NULL()
     
     ! Data type describing the layout in distributed-memory of the dual graph
     ! (It is required, e.g., for nearest neighbour comms on this graph)
     type(element_import_t)                  :: element_import   
     
     ! Perhaps the following three member variables should be packed within type(map_t) ?
     ! I didn't do that because type(map_t) has extra members that do not make sense anymore
     ! for the current situation with objects (i.e., interior, boundary, external) etc.
     integer(ip)                             :: number_global_objects = -1
     integer(ip)                             :: number_objects = -1
     integer(igp), allocatable               :: objects_gids(:)
     integer(ip), allocatable                :: objects_dimension(:)
     type(list_t)                            :: vefs_object
     type(list_t)                            :: parts_object
     type(coarse_triangulation_t), pointer   :: coarse_triangulation
  contains
     ! Private methods for creating coarse objects-related data
     procedure, non_overridable, private :: compute_parts_itfc_vefs                        => par_base_static_tria_compute_parts_itfc_vefs
     procedure, non_overridable, private :: compute_vefs_and_parts_object                  => par_base_static_tria_compute_vefs_and_parts_object
     procedure, non_overridable, private :: compute_objects_dimension                      => par_base_static_tria_compute_objects_dimension
     procedure, non_overridable, private :: compute_objects_neighbours_exchange_data       => par_base_static_tria_compute_objects_neighbours_exchange_data
     procedure, non_overridable, private :: compute_number_global_objects_and_their_gids   => par_base_static_tria_compute_num_global_objs_and_their_gids
     
     ! Private methods for coarser triangulation set-up
     procedure, non_overridable, private :: setup_coarse_triangulation                     => par_base_static_tria_setup_coarse_triangulation
     procedure, non_overridable, private :: gather_coarse_cell_gids                        => par_base_static_tria_gather_coarse_cell_gids
     procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs       => par_base_static_tria_gather_coarse_vefs_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_vefs_gids                        => par_base_static_tria_gather_coarse_vefs_gids
     procedure, non_overridable, private :: gather_coarse_vefs_dimension                   => par_base_static_tria_gather_coarse_vefs_dimension
     procedure, non_overridable, private :: fetch_l2_part_id_neighbours                    => par_base_static_tria_fetch_l2_part_id_neighbours
     procedure, non_overridable, private :: gather_coarse_dgraph_rcv_counts_and_displs     => par_base_static_tria_gather_coarse_dgraph_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_dgraph_lextn_and_lextp           => par_base_static_tria_gather_coarse_dgraph_lextn_and_lextp
     procedure, non_overridable, private :: adapt_coarse_raw_arrays                        => par_base_static_tria_adapt_coarse_raw_arrays
  end type par_base_static_triangulation_t
  
  type, extends(base_static_triangulation_t) :: serial_triangulation_t
  contains
     procedure, non_overridable          :: create                              => serial_triangulation_create
     procedure, non_overridable          :: free                                => serial_triangulation_free
     procedure, non_overridable          :: print                               => serial_triangulation_print
  end type serial_triangulation_t
  
  ! type, extends(par_base_static_triangulation_t) :: par_triangulation_t
  ! end type par_triangulation_t
  
  type, extends(par_base_static_triangulation_t) :: coarse_triangulation_t
  contains
     procedure, non_overridable          :: create                              => coarse_triangulation_create
     procedure, non_overridable          :: free                                => coarse_triangulation_free
     procedure, non_overridable          :: print                               => coarse_triangulation_print
     
     ! Private methods for creating cell-related data
     procedure, non_overridable, private :: allocate_and_fill_ptr_vefs_per_cell => coarse_triangulation_allocate_and_fill_ptr_vefs_per_cell
     procedure, non_overridable, private :: allocate_cells_gid                  => coarse_triangulation_allocate_cells_gid
     procedure, non_overridable, private :: free_cells_gid                      => coarse_triangulation_free_cells_gid
     procedure, non_overridable, private :: fill_local_cells_gid                => coarse_triangulation_fill_local_cells_gid
     procedure, non_overridable, private :: allocate_cells_mypart               => coarse_triangulation_allocate_cells_mypart
     procedure, non_overridable, private :: free_cells_mypart                   => coarse_triangulation_free_cells_mypart
     procedure, non_overridable, private :: fill_local_cells_mypart             => coarse_triangulation_fill_local_cells_mypart
     procedure, non_overridable, private :: fetch_ghost_cells_data              => coarse_triangulation_fetch_ghost_cells_data
     procedure, non_overridable, nopass, private :: cell_size                   => coarse_triangulation_cell_size
     procedure, non_overridable, nopass, private :: cell_pack                   => coarse_triangulation_cell_pack
     procedure, non_overridable, nopass, private :: cell_unpack                 => coarse_triangulation_cell_unpack
     procedure, non_overridable, private :: allocate_and_fill_lst_vefs_lids     => coarse_triangulation_allocate_and_fill_lst_vefs_lids 
  end type coarse_triangulation_t

  integer(ieep), parameter :: mold(1) = [0_ieep]
  integer(ip)  , parameter :: size_of_ip = size(transfer(1_ip, mold))
  integer(ip)  , parameter :: size_of_igp = size(transfer(1_igp ,mold))

  public :: serial_triangulation_t, coarse_triangulation_t 
  
contains

#include "sbm_cell_accessor.i90"
#include "sbm_cell_iterator.i90"
#include "sbm_vef_accessor.i90"
#include "sbm_vef_iterator.i90"
#include "sbm_base_static_triangulation.i90"
#include "sbm_par_base_static_triangulation.i90"
#include "sbm_serial_triangulation.i90"
#include "sbm_coarse_triangulation.i90"

end module base_static_triangulation_names
