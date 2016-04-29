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
  use element_import_names
  use hash_table_names
  use list_types_names

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
    procedure, non_overridable, private  :: get_triangulation    => cell_accessor_get_triangulation
    procedure, non_overridable           :: past_the_end         => cell_accessor_past_the_end
    procedure, non_overridable           :: get_reference_fe_geo => cell_accessor_get_reference_fe_geo
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
  
  type base_static_triangulation_t ! Base class for serial_triangulation_t and base_par_static_triangulation_t
     private
     integer(ip)                           :: num_dimensions  = -1
     integer(ip)                           :: num_local_cells = -1
     integer(ip)                           :: num_ghost_cells = -1
     
     integer(igp), allocatable             :: elems_gid(:)               ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: elems_part(:)              ! Num local cells + num ghost cells
     type(p_reference_fe_t)                :: reference_fe_geo_list(1)
     integer(ip) , allocatable             :: elems_reference_fe_geo(:)  ! Num local cells + num ghost cells
     integer(ip) , allocatable             :: ptr_vefs_per_cell(:)       ! Num local cells + num ghost cells + 1
     integer(ip) , allocatable             :: lst_vefs_lids(:)
          
     integer(ip)                           :: num_local_vefs = -1
     integer(ip)                           :: num_ghost_vefs = -1
     integer(igp), allocatable             :: vefs_gid(:)          ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_dimension(:)    ! num_local_vefs + num_ghost_vefs
     integer(ip) , allocatable             :: vefs_itfc_lid(:)     ! num_local_vefs + num_ghost_vefs
     
     integer(ip)                           :: num_itfc_vefs  = -1
     integer(ip), allocatable              :: lst_itfc_vefs(:)
     integer(ip), allocatable              :: ptrs_cells_around(:) ! num_itfc_vefs+1
     integer(ip), allocatable              :: lst_cells_around(:)  ! ptrs_cells_around(num_itfc_vefs+1)-1
  contains     
     ! Cell traversals-related TBPs
     procedure, non_overridable            :: create_cell_iterator      => base_static_triangulation_create_cell_iterator
  
     ! Vef traversals-related TBPs
     procedure, non_overridable            :: create_vef_iterator       => base_static_triangulation_create_vef_iterator
     procedure, non_overridable            :: create_itfc_vef_iterator  => base_static_triangulation_create_itfc_vef_iterator
  end type base_static_triangulation_t

  public :: base_static_triangulation_t
  
contains

  subroutine cell_accessor_create ( this, lid, base_static_triangulation )
    implicit none
    class(cell_accessor_t)                   , intent(inout) :: this
    integer(ip)                              , intent(in)    :: lid
    type(base_static_triangulation_t), target, intent(in)    :: base_static_triangulation
    call this%free()
    this%lid = lid
    this%base_static_triangulation => base_static_triangulation
  end subroutine cell_accessor_create 
  
  subroutine cell_accessor_free ( this)
    implicit none
    class(cell_accessor_t), intent(inout) :: this
    this%lid = -1
    nullify ( this%base_static_triangulation )
  end subroutine cell_accessor_free
  
  subroutine cell_accessor_next(this)
    implicit none
    class(cell_accessor_t), intent(inout) :: this
    this%lid = this%lid + 1
  end subroutine cell_accessor_next
  
  subroutine cell_accessor_set_lid(this, lid)
    implicit none
    class(cell_accessor_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: lid
    this%lid = lid
  end subroutine cell_accessor_set_lid
  
  function cell_accessor_past_the_end(this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    logical :: cell_accessor_past_the_end
    cell_accessor_past_the_end = (this%lid > this%base_static_triangulation%num_local_cells + &
                                             this%base_static_triangulation%num_ghost_cells)
  end function cell_accessor_past_the_end
  
  function cell_accessor_get_triangulation(this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    type(base_static_triangulation_t), pointer :: cell_accessor_get_triangulation
    cell_accessor_get_triangulation => this%base_static_triangulation
  end function cell_accessor_get_triangulation
  
  function cell_accessor_get_lid (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip) :: cell_accessor_get_lid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_lid = this%lid
  end function cell_accessor_get_lid
  
  function cell_accessor_get_gid (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(igp) :: cell_accessor_get_gid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_gid = this%base_static_triangulation%elems_gid(this%lid)
  end function cell_accessor_get_gid
  
  function cell_accessor_get_reference_fe_geo (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    class(reference_fe_t), pointer     :: cell_accessor_get_reference_fe_geo
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    assert ( allocated ( this%base_static_triangulation%elems_reference_fe_geo ) )
    cell_accessor_get_reference_fe_geo => &
      this%base_static_triangulation%reference_fe_geo_list(this%base_static_triangulation%elems_reference_fe_geo(this%lid))%p
  end function cell_accessor_get_reference_fe_geo
  
  function cell_accessor_get_mypart (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip) :: cell_accessor_get_mypart
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_mypart = this%base_static_triangulation%elems_part(this%lid)
  end function cell_accessor_get_mypart
  
  function cell_accessor_get_num_vefs (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip)                        :: cell_accessor_get_num_vefs
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_num_vefs = this%base_static_triangulation%ptr_vefs_per_cell(this%lid+1) - &
                                 this%base_static_triangulation%ptr_vefs_per_cell(this%lid)
  end function cell_accessor_get_num_vefs
  
  function cell_accessor_get_vef_lid (this, ivef)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip)                        :: ivef
    integer(ip)                        :: cell_accessor_get_vef_lid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_vef_lid = this%base_static_triangulation%lst_vefs_lids(this%base_static_triangulation%lst_vefs_lids(this%lid)+ivef-1)
  end function cell_accessor_get_vef_lid
  
  function cell_accessor_get_vef_gid (this, ivef)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip)                        :: ivef
    integer(igp)                       :: cell_accessor_get_vef_gid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_get_vef_gid = this%base_static_triangulation%vefs_gid(this%get_vef_lid(ivef))
  end function cell_accessor_get_vef_gid
  
  function cell_accessor_get_vef (this, ivef)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    integer(ip)                        :: ivef
    type(vef_iterator_t)               :: cell_accessor_get_vef
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    call cell_accessor_get_vef%create(this%get_vef_lid(ivef), this%base_static_triangulation)
  end function cell_accessor_get_vef
  
  function cell_accessor_is_local (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    logical                            :: cell_accessor_is_local
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_is_local = (this%lid <= this%base_static_triangulation%num_local_cells)
  end function cell_accessor_is_local
  
  function cell_accessor_is_ghost (this)
    implicit none
    class(cell_accessor_t), intent(in) :: this
    logical                            :: cell_accessor_is_ghost
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    cell_accessor_is_ghost = (this%lid > this%base_static_triangulation%num_local_cells)
  end function cell_accessor_is_ghost
  
    subroutine cell_iterator_create ( this, lid, base_static_triangulation ) 
    implicit none
    class(cell_iterator_t)            , intent(inout) :: this
    integer(ip)                      , intent(in)    :: lid
    type(base_static_triangulation_t), intent(in)    :: base_static_triangulation
    call this%free()
    call this%current_cell_accessor%create(lid=lid, base_static_triangulation=base_static_triangulation )
  end subroutine cell_iterator_create
  
  subroutine cell_iterator_free ( this ) 
    implicit none
    class(cell_iterator_t), intent(inout) :: this
    call this%current_cell_accessor%free()
  end subroutine cell_iterator_free
  
  subroutine cell_iterator_init ( this ) 
    implicit none
    class(cell_iterator_t), intent(inout) :: this
    call this%current_cell_accessor%set_lid(lid=1)
  end subroutine cell_iterator_init
  
  subroutine cell_iterator_next ( this ) 
    implicit none
    class(cell_iterator_t), intent(inout) :: this
    call this%current_cell_accessor%next()
  end subroutine cell_iterator_next
  
  function cell_iterator_has_finished ( this ) 
    implicit none
    class(cell_iterator_t), intent(in) :: this
    logical                                  :: cell_iterator_has_finished
    cell_iterator_has_finished = this%current_cell_accessor%past_the_end()
  end function cell_iterator_has_finished
  
  function cell_iterator_current ( this ) 
    implicit none
    class(cell_iterator_t), target, intent(in) :: this
    type(cell_accessor_t), pointer :: cell_iterator_current
    cell_iterator_current => this%current_cell_accessor
  end function cell_iterator_current
  
  subroutine vef_accessor_create ( this, lid, base_static_triangulation )
    implicit none
    class(vef_accessor_t)               , intent(inout) :: this
    integer(ip)                         , intent(in)    :: lid
    type(base_static_triangulation_t), target, intent(in)    :: base_static_triangulation
    call this%free()
    this%lid = lid
    this%base_static_triangulation => base_static_triangulation
  end subroutine vef_accessor_create 
  
  subroutine vef_accessor_free ( this)
    implicit none
    class(vef_accessor_t), intent(inout) :: this
    this%lid = -1
    nullify ( this%base_static_triangulation )
  end subroutine vef_accessor_free
  
  subroutine vef_accessor_next(this)
    implicit none
    class(vef_accessor_t), intent(inout) :: this
    this%lid = this%lid + 1
  end subroutine vef_accessor_next
  
  subroutine vef_accessor_set_lid(this, lid)
    implicit none
    class(vef_accessor_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: lid
    this%lid = lid
  end subroutine vef_accessor_set_lid
  
  function vef_accessor_past_the_end(this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    logical :: vef_accessor_past_the_end
    vef_accessor_past_the_end = (this%lid > this%base_static_triangulation%num_local_vefs + &
                                            this%base_static_triangulation%num_ghost_vefs)
  end function vef_accessor_past_the_end
  
  function vef_accessor_get_triangulation(this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    type(base_static_triangulation_t), pointer :: vef_accessor_get_triangulation
    vef_accessor_get_triangulation => this%base_static_triangulation
  end function vef_accessor_get_triangulation
  
  function vef_accessor_get_lid (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    integer(ip) :: vef_accessor_get_lid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_get_lid = this%lid
  end function vef_accessor_get_lid
  
  function vef_accessor_get_gid (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    integer(igp) :: vef_accessor_get_gid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_get_gid = this%base_static_triangulation%vefs_gid(this%lid)
  end function vef_accessor_get_gid
  
  function vef_accessor_at_interface (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    logical :: vef_accessor_at_interface 
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_at_interface  = (this%base_static_triangulation%vefs_itfc_lid(this%lid) /= -1 )
  end function vef_accessor_at_interface 
  
  function vef_accessor_get_dimension(this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    logical :: vef_accessor_get_dimension
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_get_dimension  = this%base_static_triangulation%vefs_dimension(this%lid)
 end function vef_accessor_get_dimension
  
  function vef_accessor_get_num_cells_around (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    integer(ip) :: vef_accessor_get_num_cells_around
    integer(ip) :: vef_itfc_lid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    assert ( this%at_interface() )
    vef_itfc_lid = this%base_static_triangulation%vefs_itfc_lid(this%lid)
    vef_accessor_get_num_cells_around = this%base_static_triangulation%ptrs_cells_around(vef_itfc_lid+1)- &
                                      this%base_static_triangulation%ptrs_cells_around(vef_itfc_lid)
  end function vef_accessor_get_num_cells_around
  
  function vef_accessor_get_cell_around (this, icell_around)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    integer(ip)          , intent(in) :: icell_around
    integer(ip)                     :: vef_accessor_get_cell_around
    integer(ip)                     :: position_in_lst_cells_around
    integer(ip)                     :: vef_itfc_lid
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    assert ( this%at_interface() )
    assert ( icell_around >= 1 .and. icell_around <= this%get_num_cells_around() )
    !vef_itfc_lid = this%base_static_triangulation%vefs_itfc_lid(this%lid)
    !position_in_lst_cells_around = this%base_static_triangulation%ptrs_cells_around(vef_itfc_lid) + icell_around-1
    !vef_accessor_get_cell_around => this%base_static_triangulation%cells(this%base_static_triangulation%lst_cells_around(position_in_lst_cells_around))
  end function vef_accessor_get_cell_around
  
  function vef_accessor_is_local (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    logical                           :: vef_accessor_is_local
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_is_local = (this%lid <= this%base_static_triangulation%num_local_vefs)
  end function vef_accessor_is_local
  
  function vef_accessor_is_ghost (this)
    implicit none
    class(vef_accessor_t), intent(in) :: this
    logical                            :: vef_accessor_is_ghost
    assert ( this%lid >= 1 .and. .not. this%past_the_end() )
    vef_accessor_is_ghost = (this%lid > this%base_static_triangulation%num_ghost_vefs)
  end function vef_accessor_is_ghost
  
  subroutine vef_iterator_create ( this, lid, base_static_triangulation ) 
    implicit none
    class(vef_iterator_t)            , intent(inout) :: this
    integer(ip)                      , intent(in)    :: lid
    type(base_static_triangulation_t), intent(in)    :: base_static_triangulation
    call this%free()
    call this%current_vef_accessor%create(lid=lid, base_static_triangulation=base_static_triangulation )
  end subroutine vef_iterator_create
  
  subroutine vef_iterator_free ( this ) 
    implicit none
    class(vef_iterator_t), intent(inout) :: this
    call this%current_vef_accessor%free()
  end subroutine vef_iterator_free
  
  subroutine vef_iterator_init ( this ) 
    implicit none
    class(vef_iterator_t), intent(inout) :: this
    call this%current_vef_accessor%set_lid(lid=1)
  end subroutine vef_iterator_init
  
  subroutine vef_iterator_next ( this ) 
    implicit none
    class(vef_iterator_t), intent(inout) :: this
    call this%current_vef_accessor%next()
  end subroutine vef_iterator_next
  
  function vef_iterator_has_finished ( this ) 
    implicit none
    class(vef_iterator_t), intent(in) :: this
    logical                                  :: vef_iterator_has_finished
    vef_iterator_has_finished = this%current_vef_accessor%past_the_end()
  end function vef_iterator_has_finished
  
  function vef_iterator_current ( this ) 
    implicit none
    class(vef_iterator_t), target, intent(in) :: this
    type(vef_accessor_t), pointer :: vef_iterator_current
    vef_iterator_current => this%current_vef_accessor
  end function vef_iterator_current
  
  subroutine itfc_vef_iterator_create ( this, base_static_triangulation ) 
    implicit none
    class(itfc_vef_iterator_t), intent(inout) :: this
    type(base_static_triangulation_t)     , intent(in)    :: base_static_triangulation
    call this%free()
    this%itfc_lid = 1
    if ( base_static_triangulation%num_itfc_vefs == 0 ) then
       call this%current_vef_accessor%create(lid=base_static_triangulation%num_local_vefs+base_static_triangulation%num_ghost_vefs+1, &
                                             base_static_triangulation=base_static_triangulation)
    else
       call this%current_vef_accessor%create(lid=base_static_triangulation%lst_itfc_vefs(this%itfc_lid), &
                                    base_static_triangulation=base_static_triangulation)
    end if
  end subroutine itfc_vef_iterator_create
  
  subroutine itfc_vef_iterator_free ( this ) 
    implicit none
    class(itfc_vef_iterator_t), intent(inout) :: this
    this%itfc_lid = -1
    call this%current_vef_accessor%free()
  end subroutine itfc_vef_iterator_free
  
  subroutine itfc_vef_iterator_init (this) 
    implicit none
    class(itfc_vef_iterator_t), intent(inout) :: this
    type(base_static_triangulation_t), pointer            :: base_static_triangulation
    base_static_triangulation => this%current_vef_accessor%get_triangulation()
    this%itfc_lid = 1
    if ( base_static_triangulation%num_itfc_vefs == 0 ) then
      call this%current_vef_accessor%set_lid(lid=base_static_triangulation%num_local_vefs+base_static_triangulation%num_ghost_vefs+1)
    else
      call this%current_vef_accessor%set_lid(lid=base_static_triangulation%lst_itfc_vefs(this%itfc_lid))
    end if
  end subroutine itfc_vef_iterator_init
  
  subroutine itfc_vef_iterator_next ( this ) 
    implicit none
    class(itfc_vef_iterator_t), intent(inout) :: this
    type(base_static_triangulation_t), pointer            :: base_static_triangulation
    base_static_triangulation => this%current_vef_accessor%get_triangulation()
    this%itfc_lid = this%itfc_lid + 1
    if ( this%itfc_lid > base_static_triangulation%num_itfc_vefs ) then
      call this%current_vef_accessor%set_lid(lid=base_static_triangulation%num_local_vefs+base_static_triangulation%num_ghost_vefs+1)
    else
      call this%current_vef_accessor%set_lid(lid=base_static_triangulation%lst_itfc_vefs(this%itfc_lid))
    end if
  end subroutine itfc_vef_iterator_next
  
  function itfc_vef_iterator_has_finished ( this ) 
    implicit none
    class(itfc_vef_iterator_t), intent(in) :: this
    logical                                  :: itfc_vef_iterator_has_finished
    itfc_vef_iterator_has_finished = this%current_vef_accessor%past_the_end()
  end function itfc_vef_iterator_has_finished
  
  function itfc_vef_iterator_current ( this ) 
    implicit none
    class(itfc_vef_iterator_t), target, intent(in) :: this
    type(vef_accessor_t), pointer            :: itfc_vef_iterator_current
    itfc_vef_iterator_current => this%current_vef_accessor
  end function itfc_vef_iterator_current
  
  function base_static_triangulation_create_vef_iterator ( this )
    implicit none
    class(base_static_triangulation_t), intent(in)    :: this
    type(vef_iterator_t) :: base_static_triangulation_create_vef_iterator
    call base_static_triangulation_create_vef_iterator%create(1, this)
  end function base_static_triangulation_create_vef_iterator
  
  function base_static_triangulation_create_itfc_vef_iterator ( this )
    implicit none
    class(base_static_triangulation_t), intent(in)    :: this
    type(itfc_vef_iterator_t) :: base_static_triangulation_create_itfc_vef_iterator
    call base_static_triangulation_create_itfc_vef_iterator%create(this)
  end function base_static_triangulation_create_itfc_vef_iterator
  
  function base_static_triangulation_create_cell_iterator ( this )
    implicit none
    class(base_static_triangulation_t), intent(in)    :: this
    type(cell_iterator_t) :: base_static_triangulation_create_cell_iterator
    call base_static_triangulation_create_cell_iterator%create(1, this)
  end function base_static_triangulation_create_cell_iterator
  
end module base_static_triangulation_names
