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
module coarse_triangulation_names
  ! Serial modules
  use types_names
  use memor_names
  use migratory_element_names
  use element_import_names
  use hash_table_names
  use list_types_names
  
  ! Parallel modules
  use par_environment_names
  use par_element_exchange_names

  implicit none
# include "debug.i90"
  private
  
  type, extends(migratory_element_t) :: coarse_cell_t
    private
    integer(ip)                           :: gid      = -1
    integer(ip)                           :: mypart   = -1
    integer(ip)                           :: num_vefs = -1
    integer(ip), allocatable              :: vefs_lids(:)
    integer(ip), allocatable              :: vefs_gids(:)
    ! A pointer to type(coarse_triangulation_t) is pending here to increase performance
  contains
    procedure, non_overridable :: create             => coarse_cell_create
    procedure, non_overridable :: free               => coarse_cell_free 
    procedure, non_overridable :: match_vefs_lids    => coarse_cell_match_vefs_lids
    procedure, non_overridable :: get_gid            => coarse_cell_get_gid
    procedure, non_overridable :: get_num_vefs       => coarse_cell_get_num_vefs
    procedure, non_overridable :: get_vef_lid        => coarse_cell_get_vef_lid
    procedure, non_overridable :: get_vef_gid        => coarse_cell_get_vef_gid
    procedure, non_overridable :: print              => coarse_cell_print
    procedure                  :: size               => coarse_cell_size
    procedure                  :: pack               => coarse_cell_pack
    procedure                  :: unpack             => coarse_cell_unpack
  end type coarse_cell_t
  
  type coarse_vef_t
     private
     ! An integer member variable with the dimension (vertex=0, edge=1, face=2) 
     ! of the vef will be needed in the future
     integer(ip)                           :: gid              = -1  
     integer(ip)                           :: num_cells_around = -1 
     integer(ip), allocatable              :: cells_around(:)
     type(coarse_triangulation_t), pointer :: coarse_triangulation
  contains
     procedure, non_overridable :: create               => coarse_vef_create
     procedure, non_overridable :: free                 => coarse_vef_free
     procedure, non_overridable :: print                => coarse_vef_print
     procedure, non_overridable :: get_gid              => coarse_vef_get_gid
     procedure, non_overridable :: get_num_cells_around => coarse_vef_get_num_cells_around
     procedure, non_overridable :: get_cell_around      => coarse_vef_get_cell_around
  end type coarse_vef_t
  
  type coarse_triangulation_t
     private
     integer(ip)                              :: num_dimensions  = -1 
     
     integer(ip)                              :: num_local_cells = -1
     integer(ip)                              :: num_ghost_cells = -1
     type(coarse_cell_t), allocatable         :: cells(:)
     
     integer(ip)                              :: num_local_vefs = -1
     type(coarse_vef_t) , allocatable         :: vefs(:)
     
     integer(ip)                              :: num_itfc_vefs  = -1
     integer(ip), allocatable                 :: lst_itfc_vefs(:)
 
     ! Parallel environment describing MPI tasks among which coarse_triangulation_t is distributed
     type(par_environment_t),   pointer       :: p_env => NULL()
     
     ! Data type describing the layout in distributed-memory of the dual graph
     ! (It is required, e.g., for nearest neighbour comms on this graph)
     type(element_import_t)                   :: element_import   
     
     ! Perhaps the following three member variables should be packed within type(map_t) ?
     ! I didn't do that because type(map_t) has extra members that do not make sense anymore
     ! for the current situation with objects (i.e., interior, boundary, external) etc.
     integer(ip)                             :: number_global_objects
     integer(ip)                             :: number_objects
     integer(ip), allocatable                :: objects_gids(:)
     
     type(list_t)                            :: vefs_object
     type(list_t)                            :: parts_object
     
     
     
     type(coarse_triangulation_t), pointer    :: next_level
  contains
     procedure, non_overridable          :: create                               => coarse_triangulation_create
     procedure, non_overridable          :: free                                 => coarse_triangulation_free
     procedure, non_overridable          :: print                                => coarse_triangulation_print
     procedure, non_overridable, private :: allocate_cell_array                  => coarse_triangulation_allocate_cell_array
     procedure, non_overridable, private :: free_cell_array                      => coarse_triangulation_free_cell_array
     procedure, non_overridable, private :: fill_local_cells                     => coarse_triangulation_fill_local_cells
     procedure, non_overridable, private :: fill_ghost_cells                     => coarse_triangulation_fill_ghost_cells
     procedure, non_overridable, private :: match_vefs_lids_of_ghost_cells       => coarse_triangulation_match_vefs_lids_of_ghost_cells
     procedure, non_overridable, private :: allocate_vef_array                   => coarse_triangulation_allocate_vef_array
     procedure, non_overridable, private :: free_vef_array                       => coarse_triangulation_free_vef_array
     procedure, non_overridable, private :: fill_vef_array                       => coarse_triangulation_fill_vef_array
     procedure, non_overridable, private :: compute_num_local_vefs               => coarse_triangulation_compute_num_local_vefs
     procedure, non_overridable, private :: get_global_vef_gid_array             => coarse_triangulation_get_global_vef_gid_array
     procedure, non_overridable, private :: get_ptrs_lst_cells_around            => coarse_triangulation_get_ptrs_lst_cells_around
     procedure, non_overridable, private :: compute_num_itfc_vefs                => coarse_triangulation_compute_num_itfc_vefs
     procedure, non_overridable, private :: allocate_lst_itfc_vefs               => coarse_triangulation_allocate_lst_itfc_vefs
     procedure, non_overridable, private :: free_lst_itfc_vefs                   => coarse_triangulation_free_lst_itfc_vefs
     procedure, non_overridable, private :: fill_lst_itfc_vefs                   => coarse_triangulation_fill_lst_itfc_vefs
     
     !procedure, non_overridable, private :: compute_parts_itfc_vefs                        => coarse_triangulation_compute_parts_itfc_vefs
     !procedure, non_overridable, private :: compute_vefs_and_parts_object                  => coarse_triangulation_compute_vefs_and_parts_object
     !procedure, non_overridable, private :: compute_objects_neighbours_exchange_data       => coarse_triangulation_compute_objects_neighbours_exchange_data
     !procedure, non_overridable, private :: compute_number_global_objects_and_their_gids   => coarse_triangulation_compute_number_global_objects_and_their_gids
     !procedure, non_overridable, private :: setup_coarse_triangulation                     => coarse_triangulation_setup_coarse_triangulation
     !procedure, non_overridable, private :: gather_coarse_cell_gids                        => coarse_triangulation_gather_coarse_cell_gidsm
     !procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs       => coarse_triangulation_gather_coarse_vefs_rcv_counts_and_displs
     !procedure, non_overridable, private :: gather_coarse_vefs_gids                        => coarse_triangulation_gather_coarse_vefs_gids
     !procedure, non_overridable, private :: fetch_l2_part_id_neighbours                    => coarse_triangulation_fetch_l2_part_id_neighbours
     !procedure, non_overridable, private :: gather_coarse_dgraph_rcv_counts_and_displs     => coarse_triangulation_gather_coarse_dgraph_rcv_counts_and_displs
     !procedure, non_overridable, private :: gather_coarse_dgraph_lextn_and_lextp           => coarse_triangulation_gather_coarse_dgraph_lextn_and_lextp
     !procedure, non_overridable, private :: adapt_coarse_raw_arrays                        => coarse_triangulation_adapt_coarse_raw_arrays
     
     
  end type coarse_triangulation_t

  public :: coarse_triangulation_t
  
contains

  subroutine coarse_cell_create ( this, gid, mypart, num_local_vefs, vefs_lids, vefs_gids )
    implicit none
    class(coarse_cell_t), intent(inout) :: this
    integer(ip)         , intent(in)    :: gid
    integer(ip)         , intent(in)    :: mypart
    integer(ip)         , intent(in)    :: num_local_vefs
    integer(ip)         , intent(in)    :: vefs_lids(num_local_vefs)
    integer(ip)         , intent(in)    :: vefs_gids(num_local_vefs)
    call this%free()
    this%gid = gid
    this%mypart = mypart
    this%num_vefs = num_local_vefs
    call memalloc ( num_local_vefs, this%vefs_lids, __FILE__, __LINE__ ) 
    call memalloc ( num_local_vefs, this%vefs_gids, __FILE__, __LINE__ )
    this%vefs_lids = vefs_lids
    this%vefs_gids = vefs_gids
  end subroutine coarse_cell_create
  
  subroutine coarse_cell_free ( this)
    implicit none
    class(coarse_cell_t), intent(inout) :: this
    this%gid = -1
    this%mypart = -1
    this%num_vefs = -1
    if ( allocated (this%vefs_lids) ) call memfree ( this%vefs_lids, __FILE__, __LINE__)
    if ( allocated (this%vefs_gids) ) call memfree ( this%vefs_gids, __FILE__, __LINE__)
  end subroutine coarse_cell_free
  
  subroutine coarse_cell_match_vefs_lids ( this, source_cell )
    implicit none
    class(coarse_cell_t), intent(inout) :: this
    type(coarse_cell_t) , intent(in)    :: source_cell
    integer(ip) :: i, j
    
    do i=1, this%num_vefs
      if ( this%vefs_lids(i) == -1 ) then
        do j=1, source_cell%num_vefs
          if ( this%vefs_gids(i) == source_cell%vefs_gids(j) ) then
            this%vefs_lids(i) = source_cell%vefs_lids(j)
          end if
        end do
      end if
    end do
 end subroutine coarse_cell_match_vefs_lids 
 
 function coarse_cell_get_gid( this )
   implicit none
   class(coarse_cell_t), intent(in) :: this
   integer(ip)                      :: coarse_cell_get_gid
   coarse_cell_get_gid = this%gid
 end function coarse_cell_get_gid
 
 function coarse_cell_get_mypart( this )
   implicit none
   class(coarse_cell_t), intent(in) :: this
   integer(ip)                      :: coarse_cell_get_mypart
   coarse_cell_get_mypart = this%mypart
 end function coarse_cell_get_mypart
 
 function coarse_cell_get_num_vefs( this )
   implicit none
   class(coarse_cell_t), intent(in) :: this
   integer(ip)                      :: coarse_cell_get_num_vefs
   coarse_cell_get_num_vefs = this%num_vefs
 end function coarse_cell_get_num_vefs
 
 function coarse_cell_get_vef_lid ( this, ivef )
   implicit none
   class(coarse_cell_t), intent(in) :: this
   integer(ip)         , intent(in) :: ivef
   integer(ip)                      :: coarse_cell_get_vef_lid
   assert ( ivef >= 1 .and. ivef <= this%num_vefs )
   coarse_cell_get_vef_lid = this%vefs_lids(ivef)
 end function coarse_cell_get_vef_lid
 
 function coarse_cell_get_vef_gid ( this, ivef )
   implicit none
   class(coarse_cell_t), intent(in) :: this
   integer(ip)         , intent(in) :: ivef
   integer(ip)                      :: coarse_cell_get_vef_gid
   assert ( ivef >= 1 .and. ivef <= this%num_vefs )
   coarse_cell_get_vef_gid = this%vefs_gids(ivef)
 end function coarse_cell_get_vef_gid
 
  subroutine coarse_cell_print ( this)
    implicit none
    class(coarse_cell_t), intent(in) :: this
    write (*,'(a)') '****print type(coarse_cell_t)'
    write (*,'(a,i10)'  ) 'gid      :', this%gid
    write (*,'(a,i10)'  ) 'mypart   :', this%mypart
    write (*,'(a,i10)'  ) 'num_vefs :', this%num_vefs
    write (*,'(a,10i10)') 'vefs_lids:', this%vefs_lids
    write (*,'(a,10i10)') 'vefs_gids:', this%vefs_gids
    write (*,'(a)') '****end print type(coarse_cell_t)'
  end subroutine coarse_cell_print
 
 subroutine coarse_cell_size (my, n)
    implicit none
    class(coarse_cell_t), intent(in)  :: my
    integer(ip)         , intent(out) :: n
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    size_of_ip = size(transfer(1_ip ,mold))
    n = size_of_ip*3 + size_of_ip*(my%num_vefs)
 end subroutine coarse_cell_size

  subroutine coarse_cell_pack (my, n, buffer)
    implicit none
    class(coarse_cell_t), intent(in)  :: my
    integer(ip)         , intent(in)  :: n
    integer(ieep)       , intent(out) :: buffer(n)
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    integer(ip)   :: start, end
    size_of_ip   = size(transfer(1_ip ,mold))
    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%gid,mold)
    start = end + 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%mypart,mold)
    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(my%num_vefs,mold)
    start = end + 1
    end   = start + my%num_vefs*size_of_ip - 1
    buffer(start:end) = transfer(my%vefs_gids,mold)
  end subroutine coarse_cell_pack

  subroutine coarse_cell_unpack (my, n, buffer)
    implicit none
    class(coarse_cell_t), intent(inout) :: my
    integer(ip)         , intent(in)    :: n
    integer(ieep)       , intent(in)    :: buffer(n)
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip
    integer(ip)   :: start, end
    call my%free()
    size_of_ip   = size(transfer(1_ip ,mold))
    start = 1
    end = start + size_of_ip -1
    my%gid = transfer(buffer(start:end), my%gid) 
    start = end + 1
    end = start + size_of_ip - 1
    my%mypart = transfer(buffer(start:end), my%mypart)
    start = end + 1
    end   = start + size_of_ip - 1
    my%num_vefs  = transfer(buffer(start:end), my%num_vefs)
    call memalloc( my%num_vefs, my%vefs_gids, __FILE__, __LINE__ )
    start = end + 1
    end   = start + my%num_vefs*size_of_ip - 1
    my%vefs_gids = transfer(buffer(start:end), my%vefs_gids)
    call memalloc( my%num_vefs, my%vefs_lids, __FILE__, __LINE__ )
    my%vefs_lids = -1 
  end subroutine coarse_cell_unpack
  
  subroutine coarse_vef_create ( this, gid, num_cells_around, cells_around, coarse_triangulation )
    implicit none
    class(coarse_vef_t)                 , intent(inout) :: this 
    integer(ip)                         , intent(in)    :: gid
    integer(ip)                         , intent(in)    :: num_cells_around
    integer(ip)                         , intent(in)    :: cells_around(num_cells_around)
    type(coarse_triangulation_t), target, intent(in)    :: coarse_triangulation 
    call this%free()
    this%gid              = gid
    this%num_cells_around = num_cells_around
    call memalloc ( num_cells_around, this%cells_around, __FILE__, __LINE__ )
    this%cells_around = cells_around
    this%coarse_triangulation => coarse_triangulation
  end subroutine coarse_vef_create
  
  subroutine coarse_vef_free ( this)
    implicit none
    class(coarse_vef_t), intent(inout) :: this
    this%gid = -1
    this%num_cells_around = -1
    if ( allocated (this%cells_around) ) call memfree ( this%cells_around, __FILE__, __LINE__)
    nullify ( this%coarse_triangulation )
  end subroutine coarse_vef_free
  
  subroutine coarse_vef_print ( this )
    implicit none
    class(coarse_vef_t), intent(inout) :: this
    integer(ip) :: i
    write (*,'(a)') '****print type(coarse_vef_t)****'
    write (*,'(a,i10)'  ) 'gid                :', this%gid
    write (*,'(a,i10)'  ) 'num_cells_around   :', this%num_cells_around
    write (*,'(a,10i10)') 'cells_around       :', this%cells_around
    write (*,'(a)') '****end print type(coarse_vef_t)****'
  end subroutine coarse_vef_print 
  
  function coarse_vef_get_gid (this)
    implicit none
    class(coarse_vef_t), intent(in) :: this
    integer(ip) :: coarse_vef_get_gid
    coarse_vef_get_gid = this%gid
  end function coarse_vef_get_gid
  
  function coarse_vef_get_num_cells_around (this)
    implicit none
    class(coarse_vef_t), intent(in) :: this
    integer(ip) :: coarse_vef_get_num_cells_around
    coarse_vef_get_num_cells_around = this%num_cells_around
  end function coarse_vef_get_num_cells_around
  
  function coarse_vef_get_cell_around (this, icell_around)
    implicit none
    class(coarse_vef_t), intent(in) :: this
    integer(ip)        , intent(in) :: icell_around
    type(coarse_cell_t), pointer    :: coarse_vef_get_cell_around
    coarse_vef_get_cell_around => this%coarse_triangulation%cells(this%cells_around(icell_around))
  end function coarse_vef_get_cell_around
  
  subroutine coarse_triangulation_create ( this, &
                                           par_environment, &
                                           num_dimensions, &
                                           num_local_cells, &
                                           cell_gids, &
                                           ptr_vefs_per_cell, &
                                           lst_vefs_gids, &
                                           num_itfc_cells, &
                                           lst_itfc_cells, &
                                           ptr_ext_neighs_per_itfc_cell, &
                                           lst_ext_neighs_gids, &
                                           lst_ext_neighs_part_ids)
    implicit none
    class(coarse_triangulation_t)      , intent(inout) :: this
    type(par_environment_t)     ,target, intent(in)    :: par_environment
    integer(ip)                        , intent(in)    :: num_dimensions
    integer(ip)                        , intent(in)    :: num_local_cells
    integer(ip)                        , intent(in)    :: cell_gids(*)
    integer(ip)                        , intent(in)    :: ptr_vefs_per_cell(*)
    integer(ip)                        , intent(in)    :: lst_vefs_gids(*)
    integer(ip)                        , intent(in)    :: num_itfc_cells
    integer(ip)                        , intent(in)    :: lst_itfc_cells(*)
    integer(ip)                        , intent(in)    :: ptr_ext_neighs_per_itfc_cell(*)
    integer(ip)                        , intent(in)    :: lst_ext_neighs_gids(*)
    integer(ip)                        , intent(in)    :: lst_ext_neighs_part_ids(*)
  
    call this%free()
    
    this%p_env => par_environment
    if(this%p_env%am_i_l1_task()) then
      ! We need to fill the element_import data structure first
      ! in order to determine the number of ghost elements. This
      ! in turn is required as a precondition for the allocate_cell_array
      ! TBP below.
      call this%element_import%create  ( this%p_env%get_l1_rank()+1, &
                                         this%p_env%get_l1_size(), &
                                         num_local_cells, &
                                         num_itfc_cells, &
                                         lst_itfc_cells(1:num_itfc_cells), & ! I was for to provide l/u bounds to let gfortran 5.3.0 compile
                                         ptr_ext_neighs_per_itfc_cell(1:num_itfc_cells+1), &
                                         lst_ext_neighs_gids(1:ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1), &
                                         lst_ext_neighs_part_ids(1:ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1))
      this%num_dimensions  = num_dimensions      
      this%num_local_cells = num_local_cells
      this%num_ghost_cells = this%element_import%get_number_ghost_elements()
      call this%allocate_cell_array()
      call this%fill_local_cells ( cell_gids, &
                                   ptr_vefs_per_cell, &
                                   lst_vefs_gids )
      call this%fill_ghost_cells(num_itfc_cells, &
                                 lst_itfc_cells, &
                                 ptr_ext_neighs_per_itfc_cell, &
                                 lst_ext_neighs_gids)
      call this%compute_num_local_vefs()
      call this%allocate_vef_array()
      call this%fill_vef_array()
      
      call this%compute_num_itfc_vefs()
      call this%allocate_lst_itfc_vefs()
      call this%fill_lst_itfc_vefs()
    end if
    call this%print()
  end subroutine coarse_triangulation_create
  
  subroutine coarse_triangulation_free ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: icell, ivef
    if ( associated(this%p_env) ) then
      if (this%p_env%am_i_l1_task()) then
        do icell=1, this%num_local_cells + this%num_ghost_cells
          call this%cells(icell)%free()
        end do
        call this%free_cell_array()
        call this%element_import%free()
        this%num_local_cells = -1
        this%num_ghost_cells = -1
        do ivef=1, this%num_local_vefs
         call this%vefs(ivef)%free()
        end do
        call this%free_vef_array()
        this%num_local_vefs = -1
        call this%free_lst_itfc_vefs()
        this%num_itfc_vefs = -1
        this%num_dimensions = -1
      end if
      nullify(this%p_env)
    end if 
  end subroutine coarse_triangulation_free
  
  subroutine coarse_triangulation_print ( this )
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: i
    assert( associated(this%p_env) )
    if ( this%p_env%am_i_l1_task() ) then
      write (*,'(a)') '****print type(coarse_triangulation_t)****'
      write (*,'(a,i10)'  ) 'num_dimensions:', this%num_dimensions
      do i = 1, this%num_local_cells
        write(*,'(a,i10,a)') '**** local cell: ',i,'****'
        call this%cells(i)%print()
      end do
      do i = this%num_local_cells+1, this%num_local_cells+this%num_ghost_cells
        write(*,'(a,i10,a)') '**** ghost cell: ',i,'****'
        call this%cells(i)%print()
      end do
      do i = 1, this%num_local_vefs
        write(*,'(a,i10,a)') '**** local vef: ',i,'****'
        call this%vefs(i)%print()
      end do
      write(*,'(a,i10)') '**** num_itfc_vefs: ', this%num_itfc_vefs
      write(*,'(a,10i10)') '**** lst_itfc_vefs: ', this%lst_itfc_vefs
      write (*,'(a)') '****end print type(coarse_triangulation_t)****'
    end if
  end subroutine coarse_triangulation_print 
  
  subroutine coarse_triangulation_allocate_cell_array ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: istat
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    call this%free_cell_array()
    allocate ( this%cells(this%num_local_cells+this%num_ghost_cells), stat=istat)
    check(istat == 0)
  end subroutine coarse_triangulation_allocate_cell_array 
  
  subroutine coarse_triangulation_free_cell_array ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: istat
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    if (allocated(this%cells)) then
      deallocate (this%cells, stat=istat)
      check(istat==0)
    end if
  end subroutine coarse_triangulation_free_cell_array 
  
  subroutine coarse_triangulation_fill_local_cells ( this, &
                                                                         cell_gids, &
                                                                         ptr_vefs_per_cell,&
                                                                         lst_vefs_gids)                                                     
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: cell_gids(this%num_local_cells)
    integer(ip)                  , intent(in)    :: ptr_vefs_per_cell(this%num_local_cells+1)
    integer(ip)                  , intent(in)    :: lst_vefs_gids(ptr_vefs_per_cell(this%num_local_cells+1)-1)
    
    type(position_hash_table_t) :: next_vef_lid_avail
    integer(ip), allocatable :: lst_vefs_lids(:)
    integer(ip)              :: icell, istat, j, init_pos, end_pos                    

    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    assert ( this%num_local_cells>=0 )

    call memalloc ( ptr_vefs_per_cell(this%num_local_cells+1)-1, lst_vefs_lids, __FILE__, __LINE__ )
    call next_vef_lid_avail%init ( max(int(real(ptr_vefs_per_cell(this%num_local_cells+1)-1,rp)*0.1_rp),5) )
    do icell=1, this%num_local_cells
      init_pos = ptr_vefs_per_cell(icell)
      end_pos  = ptr_vefs_per_cell(icell+1)-1
      do j=init_pos, end_pos
        call next_vef_lid_avail%get(key=lst_vefs_gids(j), val=lst_vefs_lids(j), stat=istat)
      end do
      call this%cells(icell)%create(cell_gids(icell), &
                                    this%p_env%get_l1_rank()+1, &
                                    end_pos-init_pos+1, &
                                    lst_vefs_lids(init_pos:end_pos), &
                                    lst_vefs_gids(init_pos:end_pos) )
    end do
    call next_vef_lid_avail%free()
    call memfree ( lst_vefs_lids, __FILE__, __LINE__ )
  end subroutine coarse_triangulation_fill_local_cells 
  
  subroutine coarse_triangulation_fill_ghost_cells ( this, &
                                                     num_itfc_cells, &
                                                     lst_itfc_cells, &
                                                     ptr_ext_neighs_per_itfc_cell, &
                                                     lst_ext_neighs_gids)
   implicit none
   class(coarse_triangulation_t), intent(inout) :: this
   integer(ip)                  , intent(in)    :: num_itfc_cells
   integer(ip)                  , intent(in)    :: lst_itfc_cells(num_itfc_cells)
   integer(ip)                  , intent(in)    :: ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)
   integer(ip)                  , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1)
    
    
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    call ghost_elements_exchange ( this%p_env, &
                                   this%element_import, &
                                   this%cells ) 
    
    call this%match_vefs_lids_of_ghost_cells ( num_itfc_cells, &
                                               lst_itfc_cells, &
                                               ptr_ext_neighs_per_itfc_cell, &
                                               lst_ext_neighs_gids )
  end subroutine coarse_triangulation_fill_ghost_cells
 
  subroutine coarse_triangulation_match_vefs_lids_of_ghost_cells ( this, &
                                                                   num_itfc_cells, &
                                                                   lst_itfc_cells, &
                                                                   ptr_ext_neighs_per_itfc_cell, &
                                                                   lst_ext_neighs_gids)
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: num_itfc_cells
    integer(ip)                  , intent(in)    :: lst_itfc_cells(num_itfc_cells)
    integer(ip)                  , intent(in)    :: ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)
    integer(ip)                  , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_per_itfc_cell(num_itfc_cells+1)-1)
    
    type(hash_table_ip_ip_t) :: ht_g2l_ghosts
    integer(ip)              :: icell, jcell, istat, icell_itfc, j
    
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    
    ! Hash table global to local for ghost elements
    call ht_g2l_ghosts%init(this%num_ghost_cells)
    do icell = this%num_local_cells+1, this%num_local_cells + this%num_ghost_cells
       call ht_g2l_ghosts%put( key = this%cells(icell)%gid, val = icell, stat=istat)
    end do
    
    do icell_itfc = 1, num_itfc_cells     
       icell = lst_itfc_cells(icell_itfc) 
       do j = ptr_ext_neighs_per_itfc_cell(icell), ptr_ext_neighs_per_itfc_cell(icell+1)-1
          call ht_g2l_ghosts%get(key = lst_ext_neighs_gids(j), val=jcell, stat=istat) 
          call this%cells(jcell)%match_vefs_lids(this%cells(icell))
       end do
    end do
    
    call ht_g2l_ghosts%free()
  end subroutine coarse_triangulation_match_vefs_lids_of_ghost_cells
  
  subroutine coarse_triangulation_compute_num_local_vefs ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    type(hash_table_ip_ip_t)                     :: visited_vefs
    integer(ip)                                  :: icell, ivef, vef_lid, istat
    
    call visited_vefs%init(max(5,int(real(this%num_local_cells,rp)*0.2_rp,ip))) 
    this%num_local_vefs = 0
    do icell=1, this%num_local_cells
       do ivef=1, this%cells(icell)%get_num_vefs()
          vef_lid = this%cells(icell)%get_vef_lid(ivef)
          call visited_vefs%put(key=vef_lid, val=1, stat=istat)
          if (istat == now_stored) this%num_local_vefs = this%num_local_vefs + 1
       end do
    end do
    call visited_vefs%free()
  end subroutine coarse_triangulation_compute_num_local_vefs
  
    subroutine coarse_triangulation_allocate_vef_array ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: istat
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    call this%free_vef_array()
    allocate ( this%vefs(this%num_local_vefs), stat=istat)
    check(istat == 0)
  end subroutine coarse_triangulation_allocate_vef_array 
  
  subroutine coarse_triangulation_free_vef_array ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip) :: istat
    assert ( associated ( this%p_env ) )
    assert ( this%p_env%am_i_l1_task() )
    if (allocated(this%vefs)) then
      deallocate (this%vefs, stat=istat)
      check(istat==0)
    end if
  end subroutine coarse_triangulation_free_vef_array 
  
  subroutine coarse_triangulation_fill_vef_array ( this ) 
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip), allocatable :: global_vef_gid_array(:)
    integer(ip), allocatable :: ptr_cells_around(:)
    integer(ip), allocatable :: lst_cells_around(:)
    integer(ip)              :: ivef, vef_lid, j, init_pos, end_pos

    
    call this%get_global_vef_gid_array(global_vef_gid_array)
    call this%get_ptrs_lst_cells_around(ptr_cells_around, lst_cells_around)
    
    do vef_lid=1, this%num_local_vefs
      init_pos = ptr_cells_around(vef_lid)
      end_pos  = ptr_cells_around(vef_lid+1)-1
      call this%vefs(vef_lid)%create(global_vef_gid_array(vef_lid), &
                                     end_pos-init_pos+1, &
                                     lst_cells_around(init_pos:end_pos),&
                                     this )
    end do
    call memfree(global_vef_gid_array,__FILE__,__LINE__)
    call memfree(ptr_cells_around,__FILE__,__LINE__)
    call memfree(lst_cells_around,__FILE__,__LINE__)
  end subroutine coarse_triangulation_fill_vef_array
  
  subroutine coarse_triangulation_get_global_vef_gid_array ( this, global_vef_gid_array )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip), allocatable     , intent(inout) :: global_vef_gid_array(:)
    integer(ip)                                  :: icell, ivef, vef_lid, vef_gid
    assert ( this%num_local_vefs >= 0 ) 
    
    if (allocated(global_vef_gid_array)) call memfree(global_vef_gid_array,__FILE__,__LINE__)
    call memalloc(this%num_local_vefs, global_vef_gid_array,__FILE__,__LINE__)
    do icell=1, this%num_local_cells
       do ivef=1, this%cells(icell)%get_num_vefs()
          vef_lid = this%cells(icell)%get_vef_lid(ivef)
          vef_gid = this%cells(icell)%get_vef_gid(ivef)
          global_vef_gid_array(vef_lid) = vef_gid
       end do
    end do
  end subroutine coarse_triangulation_get_global_vef_gid_array
  
  subroutine coarse_triangulation_get_ptrs_lst_cells_around ( this, ptr_cells_around, lst_cells_around )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip), allocatable     , intent(inout) :: ptr_cells_around(:)
    integer(ip), allocatable     , intent(inout) :: lst_cells_around(:)
    
    integer(ip)                                  :: icell, ivef, vef_lid
    
    if (allocated(ptr_cells_around)) call memfree(ptr_cells_around,__FILE__,__LINE__)
    if (allocated(lst_cells_around)) call memfree(lst_cells_around,__FILE__,__LINE__)
    
    call memalloc ( this%num_local_vefs+1, ptr_cells_around, __FILE__, __LINE__ )
    ptr_cells_around = 0
    
    ! Count elements around each vef
    do icell=1, this%num_local_cells + this%num_ghost_cells
       do ivef=1, this%cells(icell)%get_num_vefs()
          vef_lid = this%cells(icell)%get_vef_lid(ivef)
          if (vef_lid /= -1) then ! vef_lid == -1 then vef belongs to neighbouring processor
             ptr_cells_around(vef_lid+1) = ptr_cells_around(vef_lid+1) + 1
          end if
       end do
    end do
    ptr_cells_around(1) = 1
    do ivef=2, this%num_local_vefs+1
      ptr_cells_around(ivef) = ptr_cells_around(ivef) + ptr_cells_around(ivef-1)
    end do
    
    call memalloc ( ptr_cells_around(this%num_local_vefs+1)-1, lst_cells_around, __FILE__, __LINE__ )
    do icell=1, this%num_local_cells + this%num_ghost_cells
       do ivef=1, this%cells(icell)%get_num_vefs()
          vef_lid = this%cells(icell)%get_vef_lid(ivef)
          if (vef_lid /= -1) then ! vef_lid == -1 then vef belongs to neighbouring processor
             lst_cells_around(ptr_cells_around(vef_lid)) = icell
             ptr_cells_around(vef_lid) = ptr_cells_around(vef_lid) + 1
          end if
       end do
    end do
    
    do ivef=this%num_local_vefs+1,2,-1 
      ptr_cells_around(ivef) = ptr_cells_around(ivef-1)
    end do
    ptr_cells_around(1) = 1
  end subroutine coarse_triangulation_get_ptrs_lst_cells_around
  
  subroutine coarse_triangulation_compute_num_itfc_vefs ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)                                  :: icell, vef_lid, ivef, istat
    type(hash_table_ip_ip_t)                     :: visited_vefs
    
    call visited_vefs%init(max(5,int(real(this%num_local_vefs,rp)*0.2_rp,ip)))
    
    this%num_itfc_vefs = 0
    ! Traverse local ghost elements (all interface vefs are there)
    do icell=this%num_local_cells+1, this%num_local_cells+this%num_ghost_cells
       do ivef=1, this%cells(icell)%get_num_vefs()
          vef_lid = this%cells(icell)%get_vef_lid(ivef)
          if ( vef_lid /= -1 ) then
            call visited_vefs%put(key=vef_lid, val=1, stat=istat)
            if ( istat == now_stored ) then
              this%num_itfc_vefs = this%num_itfc_vefs + 1
            end if
          end if
       end do
    end do
    call visited_vefs%free()
  end subroutine coarse_triangulation_compute_num_itfc_vefs
  
  subroutine coarse_triangulation_allocate_lst_itfc_vefs ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    assert ( this%num_itfc_vefs >= 0 )
    call this%free_lst_itfc_vefs()
    call memalloc(this%num_itfc_vefs, this%lst_itfc_vefs,__FILE__,__LINE__)
  end subroutine coarse_triangulation_allocate_lst_itfc_vefs
  
  subroutine coarse_triangulation_free_lst_itfc_vefs( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    if (allocated(this%lst_itfc_vefs)) call memfree(this%lst_itfc_vefs,__FILE__,__LINE__)
  end subroutine coarse_triangulation_free_lst_itfc_vefs
  
  subroutine coarse_triangulation_fill_lst_itfc_vefs( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)                                  :: icell, vef_lid, ivef, i, istat
    type(hash_table_ip_ip_t)                     :: visited_vefs

    !call visited_vefs%init(max(5,int(real(this%num_local_vefs,rp)*0.2_rp,ip)))
    !i=0
    !! Traverse local ghost elements (all interface vefs are there)
    !do icell=this%num_local_cells+1, this%num_local_cells+this%num_ghost_cells
    !   do ivef=1, this%cells(icell)%get_num_vefs()
    !      vef_lid = this%cells(icell)%get_vef_lid(ivef)
    !      if ( vef_lid /= -1 ) then
    !        call visited_vefs%put(key=vef_lid, val=1, stat=istat)
    !        if ( istat == now_stored ) then
    !          i=i+1
    !          this%lst_itfc_vefs(i) = vef_lid
    !        end if
    !      end if
    !   end do
    !end do
    !call visited_vefs%free()
  end subroutine coarse_triangulation_fill_lst_itfc_vefs
  
    subroutine coarse_triangulation_compute_parts_itfc_vefs ( this, parts_itfc_vefs, perm_itfc_vefs )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip), allocatable  , intent(inout) :: parts_itfc_vefs(:,:)
    integer(ip), allocatable  , intent(inout) :: perm_itfc_vefs(:)
    
    integer(ip)                               :: num_neighbours
    logical, allocatable                      :: touched_neighbours(:)
    integer(ip)                               :: nparts_around, mypart_id, part_id, local_part_id
    integer(ip)                               :: ivef_itfc, ielem, vef_lid
    integer(ip)                               :: elem_lid
    integer(ip)                               :: num_rows_parts_itfc_vefs
    integer(ip), allocatable                  :: work1(:), work2(:)
    
    !assert ( this%p_env%am_i_l1_task() )
    
    !if (allocated(parts_itfc_vefs)) call memfree(parts_itfc_vefs,__FILE__,__LINE__)
    !if (allocated(perm_itfc_vefs)) call memfree(perm_itfc_vefs,__FILE__,__LINE__)
    
    !mypart_id = this%p_env%get_l1_rank() + 1 
   
    !num_neighbours = this%element_import%get_number_neighbours()    
    !call memalloc ( num_neighbours, touched_neighbours, __FILE__, __LINE__ )
    
    !! The two extra rows in parts_per_itfc_vef are required in order to: (1) hold the number of parts around an interface vef
    !!                                                                    (2) to hold mypart_id, which should be also listed among 
    !!                                                                        the parts around each vef
    !num_rows_parts_itfc_vefs = num_neighbours + 2
    !call memalloc ( num_rows_parts_itfc_vefs, this%num_itfc_vefs, parts_itfc_vefs, __FILE__, __LINE__ )
    !parts_itfc_vefs = 0
    
    !do ivef_itfc = 1, this%num_itfc_vefs
    !  vef_lid = this%lst_itfc_vefs(ivef_itfc)
    !  touched_neighbours = .false.
    !  
    !  nparts_around = 1 
    !  parts_itfc_vefs(nparts_around+1,ivef_itfc) = mypart_id
    !  
    !  do ielem=1, this%triangulation%vefs(vef_lid)%num_elems_around
    !    elem_lid = this%triangulation%vefs(vef_lid)%elems_around(ielem)
    !    part_id = this%elems(elem_lid)%mypart
    !    
    !    if ( part_id /= mypart_id ) then
    !     local_part_id = this%element_import%get_local_neighbour_id(part_id)
    !     if (.not. touched_neighbours (local_part_id)) then
    !       touched_neighbours (local_part_id) = .true.
    !       nparts_around = nparts_around + 1 
    !       parts_itfc_vefs(nparts_around+1,ivef_itfc) = part_id
    !     end if
    !    end if
    !  end do
    !  parts_itfc_vefs(1,ivef_itfc) = nparts_around
    !  ! Sort list of parts in increasing order by part identifiers
    !  ! This is required by the call to icomp subroutine below 
    !  call sort ( nparts_around, parts_itfc_vefs(2:nparts_around+1, ivef_itfc) )
    !end do
    
    !call memalloc ( this%num_itfc_vefs, perm_itfc_vefs, __FILE__, __LINE__ )
    !do ivef_itfc = 1, this%num_itfc_vefs
    !  perm_itfc_vefs(ivef_itfc) = ivef_itfc 
    !end do
    
    
    !! Re-number vefs in increasing order by the number of parts that share them, 
    !! and among vefs sharing the same list of parts, in increasing order by the list 
    !! of parts shared by the vef 
    !call memalloc ( num_rows_parts_itfc_vefs, work1, __FILE__,__LINE__ )
    !call memalloc ( num_rows_parts_itfc_vefs, work2, __FILE__,__LINE__ )
    !call sort_array_cols_by_row_section( num_rows_parts_itfc_vefs, & 
    !   &                                 num_rows_parts_itfc_vefs, & 
    !   &                                 this%num_itfc_vefs, & 
    !   &                                 parts_itfc_vefs, & 
    !   &                                 perm_itfc_vefs, &
    !   &                                 work1, &
    !   &                                 work2 ) 
    !call memfree ( work2, __FILE__,__LINE__ )
    !call memfree ( work1, __FILE__,__LINE__ )
    !call memfree ( touched_neighbours, __FILE__, __LINE__ )
    
    !do ivef_itfc=1,this%num_itfc_vefs
    !  write(6,'(10i10)') ivef_itfc, this%lst_itfc_vefs(perm_itfc_vefs(ivef_itfc)), parts_itfc_vefs(:, ivef_itfc) 
    !end do
  end subroutine coarse_triangulation_compute_parts_itfc_vefs
  
  subroutine coarse_triangulation_compute_vefs_and_parts_object(this)
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)                               :: nparts_around
    integer(ip)                               :: ivef_itfc, init_vef, end_vef
    integer(ip)                               :: iobj, ipart
    integer(ip)                               :: num_rows_parts_itfc_vefs
    integer(ip), allocatable                  :: parts_itfc_vefs (:,:)
    integer(ip), allocatable                  :: perm_itfc_vefs(:)
    type(list_iterator_t)                     :: vefs_object_iterator, parts_object_iterator

    !assert ( this%p_env%am_i_l1_task() )
    
    !call this%compute_parts_itfc_vefs(parts_itfc_vefs,perm_itfc_vefs)
    !num_rows_parts_itfc_vefs = size(parts_itfc_vefs,1)
    
    !! Count number_objects
    !ivef_itfc = 1
    !this%number_objects = 0
    !do while ( ivef_itfc <= this%num_itfc_vefs ) 
    !  if ( ivef_itfc < this%num_itfc_vefs ) then
    !    do while (icomp(num_rows_parts_itfc_vefs,parts_itfc_vefs(:,ivef_itfc),parts_itfc_vefs(:,ivef_itfc+1)) == 0)
    !      ivef_itfc = ivef_itfc + 1
    !      if ( ivef_itfc == this%num_itfc_vefs  ) exit
    !    end do
    !  end if  
    !  this%number_objects = this%number_objects + 1
    !  ivef_itfc = ivef_itfc + 1
    !end do
    !    
    !! Count number_vefs_per_object and number_parts_per_object
    !call this%vefs_object%create(n=this%number_objects)
    !call this%parts_object%create(n=this%number_objects)
    !ivef_itfc = 1
    !this%number_objects = 0
    !do while ( ivef_itfc <= this%num_itfc_vefs ) 
    !  init_vef = ivef_itfc
    !  if ( ivef_itfc < this%num_itfc_vefs ) then
    !    do while (icomp(num_rows_parts_itfc_vefs,parts_itfc_vefs(:,ivef_itfc),parts_itfc_vefs(:,ivef_itfc+1)) == 0)
    !      ivef_itfc = ivef_itfc + 1
    !      if ( ivef_itfc == this%num_itfc_vefs  ) exit
    !    end do
    !  end if  
    !  end_vef = ivef_itfc
    !  nparts_around = parts_itfc_vefs(1,end_vef)
    !  this%number_objects = this%number_objects + 1
    !  call this%parts_object%sum_to_pointer_index(this%number_objects, nparts_around)
    !  call this%vefs_object%sum_to_pointer_index(this%number_objects, end_vef-init_vef+1 )
    !  ivef_itfc = ivef_itfc + 1
    !end do
    
    !call this%vefs_object%calculate_header()
    !call this%parts_object%calculate_header()
    !call this%vefs_object%allocate_list_from_pointer()
    !call this%parts_object%allocate_list_from_pointer()
    
    !! List number_vefs_per_object and number_parts_per_object
    !ivef_itfc=1
    !do iobj=1, this%vefs_object%get_num_pointers()
    !   vefs_object_iterator = this%vefs_object%get_iterator(iobj)
    !   parts_object_iterator = this%parts_object%get_iterator(iobj)
    !   
    !   nparts_around = parts_itfc_vefs(1,ivef_itfc)
    !   do ipart=1, nparts_around
    !     call parts_object_iterator%set_current(parts_itfc_vefs(1+ipart,ivef_itfc))
    !     call parts_object_iterator%next()
    !   end do
    !   
    !   do while(.not. vefs_object_iterator%is_upper_bound())
    !    call vefs_object_iterator%set_current(this%lst_itfc_vefs(perm_itfc_vefs(ivef_itfc)))
    !    call vefs_object_iterator%next()
    !    ivef_itfc = ivef_itfc + 1
    !   end do
    !end do
    
    !call this%vefs_object%print(6)
    !call this%parts_object%print(6)
    
    !call memfree ( parts_itfc_vefs, __FILE__, __LINE__ )
    !call memfree ( perm_itfc_vefs, __FILE__, __LINE__ )
  end subroutine coarse_triangulation_compute_vefs_and_parts_object
  
  subroutine coarse_triangulation_compute_objects_neighbours_exchange_data ( this, &
                                                                          num_rcv,&
                                                                          list_rcv, &
                                                                          rcv_ptrs,&
                                                                          unpack_idx, &
                                                                          num_snd, &
                                                                          list_snd,&
                                                                          snd_ptrs,&
                                                                          pack_idx )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(out)   :: num_rcv
    integer(ip), allocatable  , intent(inout) :: list_rcv(:)    
    integer(ip), allocatable  , intent(inout) :: rcv_ptrs(:)
    integer(ip), allocatable  , intent(inout) :: unpack_idx(:)
    integer(ip)               , intent(out)   :: num_snd
    integer(ip), allocatable  , intent(inout) :: list_snd(:)    
    integer(ip), allocatable  , intent(inout) :: snd_ptrs(:)
    integer(ip), allocatable  , intent(inout) :: pack_idx(:)
    
    ! Locals
    integer(ip)                 :: part_id, my_part_id, num_neighbours
    integer(ip)                 :: i, iobj, istat
    type(list_iterator_t)       :: parts_object_iterator
    type(position_hash_table_t) :: position_parts_rcv
    integer(ip)                 :: current_position_parts_rcv
    type(position_hash_table_t) :: position_parts_snd
    integer(ip)                 :: current_position_parts_snd
    
    !assert ( this%p_env%am_i_l1_task() )
    
    !if (allocated(list_rcv)) call memfree(list_rcv,__FILE__,__LINE__)
    !if (allocated(rcv_ptrs)) call memfree(rcv_ptrs,__FILE__,__LINE__)
    !if (allocated(unpack_idx)) call memfree(unpack_idx,__FILE__,__LINE__)
    !if (allocated(list_snd)) call memfree(list_snd,__FILE__,__LINE__)
    !if (allocated(snd_ptrs)) call memfree(snd_ptrs,__FILE__,__LINE__)
    !if (allocated(pack_idx)) call memfree(pack_idx,__FILE__,__LINE__)
    
    !my_part_id     = this%p_env%get_l1_rank() + 1
    !num_neighbours = this%element_import%get_number_neighbours()  
    
    !call position_parts_rcv%init(num_neighbours)
    !call position_parts_snd%init(num_neighbours)
    !call memalloc ( num_neighbours  , list_rcv, __FILE__, __LINE__ )
    !call memalloc ( num_neighbours+1, rcv_ptrs, __FILE__, __LINE__ )
    !rcv_ptrs = 0 
    
    !call memalloc ( num_neighbours  , list_snd, __FILE__, __LINE__ )
    !call memalloc ( num_neighbours+1, snd_ptrs, __FILE__, __LINE__ )
    !snd_ptrs = 0
    
    !num_rcv = 0
    !num_snd = 0
    !do iobj=1, this%number_objects
    !   parts_object_iterator = this%parts_object%get_iterator(iobj)
    !   part_id = parts_object_iterator%get_current()
    !   if ( my_part_id == part_id ) then
    !     ! I am owner of the present object
    !     call parts_object_iterator%next()
    !     do while ( .not. parts_object_iterator%is_upper_bound() ) 
    !        part_id = parts_object_iterator%get_current()
    !        ! Insert part_id in the list of parts I have to send data
    !        ! Increment by +1 the amount of data I have to send to part_id
    !        call position_parts_snd%get(key=part_id, val=current_position_parts_snd, stat=istat)
    !        if ( istat == new_index ) then
    !          list_snd ( current_position_parts_snd ) = part_id
    !        end if
    !        snd_ptrs(current_position_parts_snd+1) = snd_ptrs(current_position_parts_snd+1)+1
    !        call parts_object_iterator%next()
    !     end do
    !   else
    !     ! I am non-owner of the present object
    !     call position_parts_rcv%get(key=part_id, val=current_position_parts_rcv, stat=istat)
    !     if ( istat == new_index ) then
    !       list_rcv ( current_position_parts_rcv ) = part_id
    !     end if
    !     rcv_ptrs(current_position_parts_rcv+1) = rcv_ptrs(current_position_parts_rcv+1)+1 
    !   end if
    !end do
   
    !num_rcv = position_parts_rcv%last()
    !num_snd = position_parts_snd%last() 
    !rcv_ptrs(1) = 1 
    !do i=1, num_rcv
    !  rcv_ptrs(i+1) = rcv_ptrs(i+1) + rcv_ptrs(i)
    !end do
    
    !snd_ptrs(1) = 1 
    !do i=1, num_snd
    !  snd_ptrs(i+1) = snd_ptrs(i+1) + snd_ptrs(i)
    !end do
    
    !call memrealloc ( num_snd+1, snd_ptrs, __FILE__, __LINE__ )
    !call memrealloc ( num_rcv+1, rcv_ptrs, __FILE__, __LINE__ )
    !call memrealloc ( num_snd, list_snd, __FILE__, __LINE__ )
    !call memrealloc ( num_rcv, list_rcv, __FILE__, __LINE__ )
    !call memalloc ( snd_ptrs(num_snd+1)-1, pack_idx, __FILE__, __LINE__ )
    !call memalloc ( rcv_ptrs(num_rcv+1)-1, unpack_idx, __FILE__, __LINE__ )
    
    !do iobj=1, this%number_objects
    !   parts_object_iterator = this%parts_object%get_iterator(iobj)
    !   part_id = parts_object_iterator%get_current()
    !   if ( my_part_id == part_id ) then
    !     ! I am owner of the present object
    !     call parts_object_iterator%next()
    !     do while ( .not. parts_object_iterator%is_upper_bound() ) 
    !       part_id = parts_object_iterator%get_current()
    !       call position_parts_snd%get(key=part_id, val=current_position_parts_snd, stat=istat)
    !       pack_idx (snd_ptrs(current_position_parts_snd)) = iobj
    !       snd_ptrs(current_position_parts_snd) = snd_ptrs(current_position_parts_snd)+1
    !       call parts_object_iterator%next()
    !     end do
    !   else
    !     ! I am non-owner of the present object
    !     call position_parts_rcv%get(key=part_id, val=current_position_parts_rcv, stat=istat)
    !     unpack_idx (rcv_ptrs(current_position_parts_rcv)) = iobj
    !     rcv_ptrs(current_position_parts_rcv) = rcv_ptrs(current_position_parts_rcv)+1 
    !   end if
    !end do
    
    !do i=num_snd, 2, -1
    !  snd_ptrs(i) = snd_ptrs(i-1) 
    !end do
    !snd_ptrs(1) = 1 
    
    !do i=num_rcv, 2, -1
    !  rcv_ptrs(i) = rcv_ptrs(i-1) 
    !end do
    !rcv_ptrs(1) = 1
    
    !call position_parts_rcv%free()
    !call position_parts_snd%free()
  end subroutine coarse_triangulation_compute_objects_neighbours_exchange_data 
  
  subroutine coarse_triangulation_compute_num_global_objects_and_their_gids ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this

    integer(ip)               :: num_rcv
    integer(ip), allocatable  :: list_rcv(:)    
    integer(ip), allocatable  :: rcv_ptrs(:)
    integer(ip), allocatable  :: unpack_idx(:)
    
    integer(ip)               :: num_snd
    integer(ip), allocatable  :: list_snd(:)    
    integer(ip), allocatable  :: snd_ptrs(:)
    integer(ip), allocatable  :: pack_idx(:)
   
    integer(ip)               :: number_local_objects_with_gid
    integer(ip), allocatable  :: local_objects_with_gid(:)
    integer(ip), allocatable  :: per_rank_objects_with_gid(:)
    integer(ip)               :: start_object_gid
    type(list_iterator_t)     :: parts_object_iterator
    integer(ip)               :: my_part_id, number_parts
    integer(ip)               :: i, iobj
    integer                   :: my_rank
    
    integer       , parameter :: root_pid = 0
    integer(ip)               :: dummy_integer_array(1)

    !assert ( this%p_env%am_i_l1_task() )
    !my_rank      = this%p_env%get_l1_rank() 
    !my_part_id   = my_rank + 1 
    !number_parts = this%p_env%get_l1_size()
    
    !! 1. Count/list how many local objects I am responsible to assign a global ID
    !call memalloc ( this%number_objects, local_objects_with_gid, __FILE__, __LINE__ )
    !number_local_objects_with_gid = 0
    !do iobj=1, this%number_objects
    !  parts_object_iterator = this%parts_object%get_iterator(iobj)
    !  if ( my_part_id == parts_object_iterator%get_current() ) then
    !    number_local_objects_with_gid = number_local_objects_with_gid + 1
    !    local_objects_with_gid (number_local_objects_with_gid) = iobj
    !  end if
    !end do
    
    !! 2. Gather + Scatter
    !if ( my_rank == root_pid ) then
    !  call memalloc( number_parts+1, per_rank_objects_with_gid, __FILE__,__LINE__ )
    !  call this%p_env%l1_gather (root=root_pid, &
    !                             input_data=number_local_objects_with_gid, &
    !                             output_data=per_rank_objects_with_gid(2:) ) 
    !   ! Transform length to header
    !   per_rank_objects_with_gid(1)=1 
    !   do i=1, number_parts
    !      per_rank_objects_with_gid(i+1) = per_rank_objects_with_gid(i) + per_rank_objects_with_gid(i+1) 
    !   end do
    !   this%number_global_objects = per_rank_objects_with_gid(number_parts+1)-1 
    !else
    !  call this%p_env%l1_gather (root=root_pid, &
    !                             input_data=number_local_objects_with_gid, &
    !                             output_data=dummy_integer_array ) 
    !end if
    
    !call this%p_env%l1_bcast (root=root_pid, data = this%number_global_objects )
    
    !if ( my_rank == root_pid ) then
    !  call this%p_env%l1_scatter (root=root_pid, &
    !                              input_data=per_rank_objects_with_gid, &
    !                              output_data=start_object_gid) 
    !  call memfree( per_rank_objects_with_gid, __FILE__,__LINE__ )
    !else
    !  call this%p_env%l1_scatter (root=root_pid, &
    !                              input_data=dummy_integer_array, &
    !                              output_data=start_object_gid) 
    !end if
    
    
    !call memalloc (this%number_objects, this%objects_gids)
    !do i=1, number_local_objects_with_gid
    !  this%objects_gids ( local_objects_with_gid(i) ) = start_object_gid
    !  start_object_gid = start_object_gid + 1 
    !end do
    
    !! Set-up objects nearest neighbour exchange data
    !! num_rcv, rcv_ptrs, lst_rcv, unpack_idx
    !! num_snd, snd_ptrs, lst_snd, pack_idx    
    !call this%compute_objects_neighbours_exchange_data ( num_rcv, &
    !                                                     list_rcv,&
    !                                                     rcv_ptrs,&
    !                                                     unpack_idx,&
    !                                                     num_snd,&
    !                                                     list_snd,&
    !                                                     snd_ptrs,&
    !                                                     pack_idx )
    
    !call this%p_env%l1_neighbours_exchange ( num_rcv, &
    !                                         list_rcv,&
    !                                         rcv_ptrs,&
    !                                         unpack_idx,&
    !                                         num_snd,&
    !                                         list_snd,&
    !                                         snd_ptrs,&
    !                                         pack_idx,&
    !                                         this%objects_gids )
    
    !call memfree ( list_rcv, __FILE__, __LINE__ )
    !call memfree ( rcv_ptrs, __FILE__, __LINE__ )
    !call memfree ( unpack_idx, __FILE__, __LINE__ )
    !call memfree ( list_snd, __FILE__, __LINE__ )
    !call memfree ( snd_ptrs, __FILE__, __LINE__ )
    !call memfree ( pack_idx, __FILE__, __LINE__ )
    !call memfree ( local_objects_with_gid, __FILE__, __LINE__ )
  end subroutine coarse_triangulation_compute_num_global_objects_and_their_gids
  
  subroutine coarse_triangulation_setup_coarse_triangulation ( this )
    implicit none
    class(coarse_triangulation_t), intent(inout) :: this
    integer(ip)               , allocatable   :: coarse_cell_gids(:)
    integer(ip)               , allocatable   :: coarse_vefs_recv_counts(:)
    integer(ip)               , allocatable   :: coarse_vefs_displs(:)
    integer(ip)               , allocatable   :: lst_coarse_vef_gids(:)
    integer(ip)               , allocatable   :: l2_part_id_neighbours(:)
    integer(ip)               , allocatable   :: coarse_dgraph_recv_counts(:)
    integer(ip)               , allocatable   :: coarse_dgraph_displs(:)
    integer(ip)               , allocatable   :: lextn(:)
    integer(ip)               , allocatable   :: lextp(:)
    
    integer(ip)                      :: i, istat
    integer(ip)                      :: num_dimensions
    integer(ip)                      :: num_local_coarse_cells
    integer(ip)                      :: num_itfc_coarse_cells
    
    !! All MPI tasks (even if they are not involved in the L2 from L1 gather) should also allocate the
    !! allocatable arrays due to the fact that non-allocated allocatable arrays cannot
    !! be passed as actual arguments of dummy arguments that do not have the allocatable attribute 
    !! (see e.g. coarse_triangulation%create() below). Otherwise, the code crashes with a segmentation fault. 
    !! Likewise, actual arguments which are used as input dummy arguments to size another array-type dummy arguments should also
    !! be initialized on all MPI tasks
    !num_local_coarse_cells = 0
    !num_itfc_coarse_cells  = 0
    !call memalloc (0, coarse_cell_gids, __FILE__, __LINE__)
    !call memalloc (0, coarse_vefs_recv_counts, __FILE__, __LINE__)
    !call memalloc (0, coarse_vefs_displs, __FILE__, __LINE__)
    !call memalloc (0, lst_coarse_vef_gids, __FILE__, __LINE__)
    !call memalloc (0, l2_part_id_neighbours, __FILE__, __LINE__)
    !call memalloc (0, coarse_dgraph_recv_counts, __FILE__, __LINE__)
    !call memalloc (0, coarse_dgraph_displs, __FILE__, __LINE__)
    !call memalloc (0, lextn, __FILE__, __LINE__)
    !call memalloc (0, lextp, __FILE__, __LINE__)
    !    
    !! L2 tasks gather from L1 tasks all raw data required to set-up the coarse triangulation on L2 tasks
    !if ( this%p_env%am_i_l1_to_l2_task() ) then
    !  num_dimensions = this%triangulation%num_dims
    !  call this%p_env%l1_to_l2_transfer ( num_dimensions ) 
    !  call this%gather_coarse_cell_gids (coarse_cell_gids)
    !  call this%gather_coarse_vefs_rcv_counts_and_displs (coarse_vefs_recv_counts, coarse_vefs_displs)
    !  call this%gather_coarse_vefs_gids (coarse_vefs_recv_counts, coarse_vefs_displs, lst_coarse_vef_gids)
    !  call this%fetch_l2_part_id_neighbours(l2_part_id_neighbours)
    !  call this%gather_coarse_dgraph_rcv_counts_and_displs ( l2_part_id_neighbours, &
    !                                                         coarse_dgraph_recv_counts, &
    !                                                         coarse_dgraph_displs )
    !  call this%gather_coarse_dgraph_lextn_and_lextp ( l2_part_id_neighbours, &
    !                                                   coarse_dgraph_recv_counts, &
    !                                                   coarse_dgraph_displs, &
    !                                                   lextn, &
    !                                                   lextp )
    !  ! Evaluate number of local coarse cells
    !  num_local_coarse_cells = this%p_env%get_l1_to_l2_size()-1
    !  
    !  ! Evaluate number of interface coarse cells
    !  ! Adapt and re-use coarse_vefs_displs/coarse_dgraph_recv_counts/coarse_dgraph_displs
    !  ! as required by this%coarse_triangulation%create below
    !  num_itfc_coarse_cells = this%adapt_coarse_raw_arrays (coarse_vefs_displs, &
    !                                                        coarse_dgraph_recv_counts, &
    !                                                        coarse_dgraph_displs )
    !end if
    
    !if ( this%p_env%am_i_lgt1_task() ) then
    !  ! lgt1 MPI tasks (recursively) build coarse triangulation
    !  allocate  ( this%coarse_triangulation, stat = istat )
    !  check( istat == 0 )
    !  call this%coarse_triangulation%create ( par_environment              = this%p_env%get_next_level(), &
    !                                          num_dimensions               = num_dimensions, &
    !                                          num_local_cells              = num_local_coarse_cells, &
    !                                          cell_gids                    = coarse_cell_gids, &
    !                                          ptr_vefs_per_cell            = coarse_vefs_displs, &
    !                                          lst_vefs_gids                = lst_coarse_vef_gids, &
    !                                          num_itfc_cells               = num_itfc_coarse_cells, &
    !                                          lst_itfc_cells               = coarse_dgraph_recv_counts, &
    !                                          ptr_ext_neighs_per_itfc_cell = coarse_dgraph_displs, &
    !                                          lst_ext_neighs_gids          = lextn, &
    !                                          lst_ext_neighs_part_ids      = lextp )
    !else
    !  ! L1 tasks do not hold any piece of the coarse triangulation
    !  nullify(this%coarse_triangulation)
    !end if
    
    !! All tasks free raw data (see actual reason on the top part of this subroutine)
    !call memfree (coarse_cell_gids, __FILE__, __LINE__)
    !call memfree (coarse_vefs_recv_counts, __FILE__, __LINE__)
    !call memfree (coarse_vefs_displs, __FILE__, __LINE__)
    !call memfree (lst_coarse_vef_gids, __FILE__, __LINE__)
    !call memfree (l2_part_id_neighbours, __FILE__, __LINE__)
    !call memfree (coarse_dgraph_recv_counts, __FILE__, __LINE__)
    !call memfree (coarse_dgraph_displs, __FILE__, __LINE__)
    !call memfree (lextn, __FILE__, __LINE__)
    !call memfree (lextp, __FILE__, __LINE__)
  end subroutine coarse_triangulation_setup_coarse_triangulation
  
  subroutine coarse_triangulation_gather_coarse_cell_gids( this, coarse_cell_gids)
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: coarse_cell_gids(:)
    
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)
    
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  l1_to_l2_size = this%p_env%get_l1_to_l2_size()
    !  if ( allocated (coarse_cell_gids) ) call memfree ( coarse_cell_gids, __FILE__, __LINE__ )
    !  call memalloc ( l1_to_l2_size, coarse_cell_gids, __FILE__, __LINE__ )
    !  call this%p_env%l2_from_l1_gather( input_data = 0, &
    !                                     output_data = coarse_cell_gids ) 
    !else
    !  call this%p_env%l2_from_l1_gather( input_data  = this%p_env%get_l1_rank()+1, &
    !                                     output_data = dummy_integer_array ) 
    !end if
  end subroutine coarse_triangulation_gather_coarse_cell_gids
  
  
  subroutine coarse_triangulation_gather_coarse_vefs_rcv_counts_and_displs( this, recv_counts, displs )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)

    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  l1_to_l2_size = this%p_env%get_l1_to_l2_size()
    !  if ( allocated (recv_counts) ) call memfree ( recv_counts, __FILE__, __LINE__ )
    !  if ( allocated (displs) ) call memfree ( displs, __FILE__, __LINE__ )
    !  call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
    !  call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    !  call this%p_env%l2_from_l1_gather( input_data = 0, &
    !                                     output_data = recv_counts ) 
    !  displs(1) = 0
    !  do i=2, l1_to_l2_size
    !    displs(i) = displs(i-1) + recv_counts(i-1)
    !  end do
    !else
    !  call this%p_env%l2_from_l1_gather( input_data  = this%number_objects, &
    !                                     output_data = dummy_integer_array ) 
    !end if
  end subroutine coarse_triangulation_gather_coarse_vefs_rcv_counts_and_displs
  
  subroutine coarse_triangulation_gather_coarse_vefs_gids ( this, recv_counts, displs, lst_gids )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: recv_counts(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(in)    :: displs(this%p_env%get_l1_to_l2_size())
    integer(ip), allocatable  , intent(inout) :: lst_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)
    
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  l1_to_l2_size = this%p_env%get_l1_to_l2_size()
    !  if (allocated(lst_gids)) call memfree ( lst_gids, __FILE__, __LINE__ )
    !  call memalloc ( displs(l1_to_l2_size), lst_gids, __FILE__, __LINE__ )
    !  call this%p_env%l2_from_l1_gather( input_data_size = 0, &
    !                                     input_data      = dummy_integer_array, &
    !                                     recv_counts     = recv_counts, &
    !                                     displs          = displs, &
    !                                     output_data     = lst_gids )
    !else
    !  call this%p_env%l2_from_l1_gather( input_data_size = this%number_objects, &
    !                                     input_data      = this%objects_gids, &
    !                                     recv_counts     = dummy_integer_array, &
    !                                     displs          = dummy_integer_array, &
    !                                     output_data     = dummy_integer_array )
    !end if    
  end subroutine coarse_triangulation_gather_coarse_vefs_gids
  
  subroutine coarse_triangulation_fetch_l2_part_id_neighbours ( this, l2_part_id_neighbours )    
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: l2_part_id_neighbours(:)
    integer(ip) :: my_l2_part_id
    integer(ip) :: num_neighbours
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if (this%p_env%am_i_l1_task()) then
    !  num_neighbours = this%element_import%get_number_neighbours()
    !  my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
    !  if (allocated(l2_part_id_neighbours)) call memfree ( l2_part_id_neighbours, __FILE__, __LINE__ )
    !  call memalloc ( num_neighbours, l2_part_id_neighbours, __FILE__, __LINE__ )
    !  call this%p_env%l1_neighbours_exchange ( num_neighbours  = num_neighbours, &
    !                                           list_neighbours = this%element_import%get_neighbours_ids(), &
    !                                           input_data      = my_l2_part_id,&
    !                                           output_data     = l2_part_id_neighbours)
    !end if
  end subroutine coarse_triangulation_fetch_l2_part_id_neighbours
  
  subroutine coarse_triangulation_gather_coarse_dgraph_rcv_counts_and_displs ( this, &
                                                                            l2_part_id_neighbours, &
                                                                            recv_counts, &
                                                                            displs )
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: l2_part_id_neighbours(this%element_import%get_number_neighbours())
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip) :: i
    integer(ip) :: l1_to_l2_size 
    integer(ip) :: my_l2_part_id
    integer(ip) :: num_neighbours
    integer(ip) :: num_external_l2_elements
    integer(ip) :: dummy_integer_array(0)
    
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  l1_to_l2_size = this%p_env%get_l1_to_l2_size()
    !  if (allocated(recv_counts)) call memfree ( recv_counts, __FILE__, __LINE__ )
    !  if (allocated(displs)) call memfree ( displs, __FILE__, __LINE__ )
    !  call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
    !  call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
    !  call this%p_env%l2_from_l1_gather( input_data = 0, &
    !                                     output_data = recv_counts ) 
    !  displs(1) = 0
    !  do i=2, l1_to_l2_size
    !    displs(i) = displs(i-1) + recv_counts(i-1)
    !  end do
    !else
    !  assert ( this%p_env%am_i_l1_task() )
    !  num_neighbours = this%element_import%get_number_neighbours()
    !  my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
    !  num_external_l2_elements = 0
    !  do i = 1, num_neighbours
    !    if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
    !      num_external_l2_elements = num_external_l2_elements + 1
    !    end if
    !  end do
    !  call this%p_env%l2_from_l1_gather( input_data = num_external_l2_elements, &
    !                                     output_data = dummy_integer_array ) 
    !end if
  end subroutine coarse_triangulation_gather_coarse_dgraph_rcv_counts_and_displs 
  
  subroutine coarse_triangulation_gather_coarse_dgraph_lextn_and_lextp( this,                  & 
                                                                     l2_part_id_neighbours, &
                                                                     recv_counts,           &
                                                                     displs,                &
                                                                     lextn,                 &
                                                                     lextp)
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: l2_part_id_neighbours(this%element_import%get_number_neighbours())
    integer(ip)               , intent(in)    :: recv_counts(this%p_env%get_l1_to_l2_size()) 
    integer(ip)               , intent(in)    :: displs(this%p_env%get_l1_to_l2_size())
    integer(ip), allocatable  , intent(inout) :: lextn(:)
    integer(ip), allocatable  , intent(inout) :: lextp(:)
    
    integer(ip)              :: i
    integer(ip)              :: l1_to_l2_size 
    integer(ip)              :: my_l2_part_id
    integer(ip)              :: num_neighbours
    integer(ip), pointer     :: neighbours_ids(:)
    integer(ip)              :: num_external_l2_elements
    integer(ip), allocatable :: lst_external_l2_element_gids(:)
    integer(ip), allocatable :: lst_external_l2_part_ids(:)
    integer(ip)              :: dummy_integer_array(0)
    
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  l1_to_l2_size = this%p_env%get_l1_to_l2_size()
    !  if (allocated(lextn)) call memfree ( lextn, __FILE__, __LINE__ )
    !  if (allocated(lextp)) call memfree ( lextp, __FILE__, __LINE__ )
    !  call memalloc ( displs(l1_to_l2_size), lextn, __FILE__, __LINE__ )
    !  call memalloc ( displs(l1_to_l2_size), lextp, __FILE__, __LINE__ )
    !  ! Gather lextn
    !  call this%p_env%l2_from_l1_gather( input_data_size = 0, &
    !                                     input_data      = dummy_integer_array, &
    !                                     recv_counts     = recv_counts, &
    !                                     displs          = displs, &
    !                                     output_data     = lextn )
    !  ! Gather lextp
    !  call this%p_env%l2_from_l1_gather( input_data_size = 0, &
    !                                     input_data      = dummy_integer_array, &
    !                                     recv_counts     = recv_counts, &
    !                                     displs          = displs, &
    !                                     output_data     = lextp )
    !else
    !  assert ( this%p_env%am_i_l1_task() )
    !  num_neighbours =  this%element_import%get_number_neighbours()
    !  neighbours_ids => this%element_import%get_neighbours_ids()
    !  my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
    !  num_external_l2_elements = 0
    !  do i = 1, num_neighbours
    !    if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
    !      num_external_l2_elements = num_external_l2_elements + 1
    !    end if
    !  end do
    !  
    !  call memalloc (num_external_l2_elements, lst_external_l2_part_ids, __FILE__, __LINE__)
    !  call memalloc (num_external_l2_elements, lst_external_l2_element_gids,__FILE__, __LINE__)
    !  num_external_l2_elements = 0
    !  neighbours_ids => this%element_import%get_neighbours_ids()
    !  do i = 1, num_neighbours
    !    if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
    !      num_external_l2_elements = num_external_l2_elements + 1
    !      lst_external_l2_element_gids(num_external_l2_elements) = neighbours_ids(i)
    !      lst_external_l2_part_ids(num_external_l2_elements) = l2_part_id_neighbours(i)
    !    end if
    !  end do
    !  call this%p_env%l2_from_l1_gather( input_data_size = num_external_l2_elements, &
    !                                     input_data      = lst_external_l2_element_gids, &
    !                                     recv_counts     = dummy_integer_array, &
    !                                     displs          = dummy_integer_array, &
    !                                     output_data     = dummy_integer_array )
    !  
    !  call this%p_env%l2_from_l1_gather( input_data_size = num_external_l2_elements, &
    !                                     input_data      = lst_external_l2_part_ids, &
    !                                     recv_counts     = dummy_integer_array, &
    !                                     displs          = dummy_integer_array, &
    !                                     output_data     = dummy_integer_array )
    !  
    !  call memfree (lst_external_l2_part_ids   , __FILE__, __LINE__)
    !  call memfree (lst_external_l2_element_gids,__FILE__, __LINE__)
    !end if
  end subroutine coarse_triangulation_gather_coarse_dgraph_lextn_and_lextp 
  
  function coarse_triangulation_adapt_coarse_raw_arrays( this, &
                                                      coarse_vefs_displs, &
                                                      coarse_dgraph_recv_counts, &
                                                      coarse_dgraph_displs ) result(num_itfc_coarse_cells)
    implicit none
    class(coarse_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(inout) :: coarse_vefs_displs(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(inout) :: coarse_dgraph_recv_counts(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(inout) :: coarse_dgraph_displs(this%p_env%get_l1_to_l2_size())
    integer(ip)                               :: num_itfc_coarse_cells
    
    integer(ip) :: i 
    !assert ( this%p_env%am_i_l1_to_l2_task() )
    !if ( this%p_env%am_i_l1_to_l2_root() ) then
    !  ! Re-use coarse_vefs_displs as ptr_vefs_gids
    !  do i=1, size(coarse_vefs_displs)
    !    coarse_vefs_displs(i)=coarse_vefs_displs(i)+1
    !  end do
   
    !  ! Re-use coarse_dgraph_displs as ptr_ext_neighs_per_itfc_cell
    !  num_itfc_coarse_cells = 0
    !  coarse_dgraph_displs(1) = 1 
    !  do i=1, size(coarse_dgraph_recv_counts)
    !    if (coarse_dgraph_recv_counts(i) /= 0) then
    !      num_itfc_coarse_cells = num_itfc_coarse_cells+1
    !      coarse_dgraph_displs(num_itfc_coarse_cells+1) = coarse_dgraph_displs(num_itfc_coarse_cells) + &
    !                                                      coarse_dgraph_recv_counts(i)                                           
    !    end if  
    !  end do
    !  
    !  ! Re-use coarse_dgraph_recv_counts as lst_itfc_cells
    !  num_itfc_coarse_cells = 0
    !  do i=1, size(coarse_dgraph_recv_counts)
    !    if (coarse_dgraph_recv_counts(i) /= 0) then
    !      num_itfc_coarse_cells = num_itfc_coarse_cells+1
    !      coarse_dgraph_recv_counts(num_itfc_coarse_cells) = i
    !    end if  
    !  end do
    !else
    !  ! L1 tasks do not hold any itfc_coarse_cells
    !  num_itfc_coarse_cells = 0
    !end if
  end function coarse_triangulation_adapt_coarse_raw_arrays
  
end module coarse_triangulation_names
