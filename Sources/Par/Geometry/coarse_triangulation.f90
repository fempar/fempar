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
                                         lst_itfc_cells, &
                                         ptr_ext_neighs_per_itfc_cell, &
                                         lst_ext_neighs_gids, &
                                         lst_ext_neighs_part_ids)
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
  
end module coarse_triangulation_names
