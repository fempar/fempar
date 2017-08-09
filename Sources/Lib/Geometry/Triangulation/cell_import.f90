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
module cell_import_names
  use types_names
  use memor_names
  use hash_table_names
  
# include "debug.i90"  
  implicit none
  private

  ! Element_Import
  ! Host for data needed to perform nearest neighbour communications for
  ! the distributed dual graph associated to the triangulation
  type cell_import_t
     private
     integer(ip)               :: part_id
     integer(ip)               :: num_parts
     integer(ip)               :: num_ghost_elements    
     integer(ip)               :: num_neighbours          
     integer(ip), allocatable  :: neighbours_ids(:)          
     integer(ip), allocatable  :: rcv_ptrs(:)                ! How many elements does the part receive from each neighbour ?
     integer(ip), allocatable  :: rcv_leids(:)               ! Local position in the elements array in which the elements to be received by this part  
                                                             ! are going to be placed
     integer(ip), allocatable  :: snd_ptrs(:)                ! How many elements does the part send to each neighbour?
     integer(ip), allocatable  :: snd_leids(:)               ! Local elements IDs of the elements to be sent to each neighbour ?
  contains
     
     procedure, private, non_overridable :: cell_import_create_ip
     procedure, private, non_overridable :: cell_import_create_igp
     generic                             :: create => cell_import_create_ip, &
                                                      cell_import_create_igp
     procedure, private, non_overridable :: compute_num_neighbours     => cell_import_compute_num_neighbours
     procedure, private, non_overridable :: compute_neighbours_ids        => cell_import_compute_neighbour_ids
     procedure, private, non_overridable :: compute_snd_rcv_ptrs          => cell_import_compute_snd_rcv_ptrs
     procedure, private, non_overridable :: compute_snd_rcv_leids         => cell_import_compute_snd_rcv_leids
     procedure, non_overridable :: free                                   => cell_import_free
     procedure, non_overridable :: print                                  => cell_import_print
     procedure, non_overridable :: get_num_parts                       => cell_import_get_num_parts
     procedure, non_overridable :: get_num_ghost_elements              => cell_import_get_num_ghost_elements
     procedure, non_overridable :: get_num_neighbours                  => cell_import_get_num_neighbours
     procedure, non_overridable :: get_neighbours_ids                     => cell_import_get_neighbours_ids
     procedure, non_overridable :: get_rcv_ptrs                           => cell_import_get_rcv_ptrs
     procedure, non_overridable :: get_rcv_leids                          => cell_import_get_rcv_leids
     procedure, non_overridable :: get_snd_ptrs                           => cell_import_get_snd_ptrs
     procedure, non_overridable :: get_snd_leids                          => cell_import_get_snd_leids
     procedure, non_overridable :: get_global_neighbour_id                => cell_import_get_global_neighbour_id
     procedure, non_overridable :: get_local_neighbour_id                 => cell_import_get_local_neighbour_id
  end type cell_import_t

  ! Types
  public :: cell_import_t

contains

  subroutine cell_import_create_igp (this,                                 & ! Passed-object dummy argument
                                        part_id,                              & ! My part identifier
                                        num_parts,                         & ! Number of parts
                                        num_elements,                      & ! Number of local elements
                                        num_itfc_elems,                       & ! Number of (local) itfc elements
                                        lst_itfc_elems,                       & ! Local IDs of itfc elements
                                        ptr_ext_neighs_x_itfc_elem,         & ! Ptrs to start-end of lst of external neighbours per itfc elem
                                        lst_ext_neighs_gids,    & ! List of GIDs of external neighbours
                                        lst_ext_neighs_part_ids ) ! Part IDs the external neighbours are mapped to
    implicit none
    class(cell_import_t)  , intent(inout) :: this
    integer(ip)              , intent(in)    :: part_id
    integer(ip)              , intent(in)    :: num_parts
    integer(ip)              , intent(in)    :: num_elements
    integer(ip)              , intent(in)    :: num_itfc_elems
    integer(ip)              , intent(in)    :: lst_itfc_elems(num_itfc_elems)
    integer(ip)              , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(igp)             , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    integer(ip)              , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    
    call this%free()    
    this%part_id = part_id
    this%num_parts  = num_parts
    call this%compute_num_neighbours( num_itfc_elems, &
                                         ptr_ext_neighs_x_itfc_elem, &
                                         lst_ext_neighs_part_ids )
    
    call this%compute_neighbours_ids ( num_itfc_elems, &
                                       ptr_ext_neighs_x_itfc_elem, &
                                       lst_ext_neighs_part_ids )
    
    call this%compute_snd_rcv_ptrs (num_itfc_elems, &
                                    lst_itfc_elems, &
                                    ptr_ext_neighs_x_itfc_elem, &
                                    lst_ext_neighs_gids, &
                                    lst_ext_neighs_part_ids )
    
    call this%compute_snd_rcv_leids (num_elements, &
                                     num_itfc_elems, &
                                     lst_itfc_elems, &
                                     ptr_ext_neighs_x_itfc_elem, &
                                     lst_ext_neighs_gids, &
                                     lst_ext_neighs_part_ids)
    
    this%num_ghost_elements = this%rcv_ptrs(this%num_neighbours+1)-1    
  end subroutine cell_import_create_igp
  
  subroutine cell_import_create_ip (this,                                 & ! Passed-object dummy argument
                                       part_id,                              & ! My part identifier
                                       num_parts,                         & ! Number of parts
                                       num_elements,                      & ! Number of local elements
                                       num_itfc_elems,                       & ! Number of (local) itfc elements
                                       lst_itfc_elems,                       & ! Local IDs of itfc elements
                                       ptr_ext_neighs_x_itfc_elem,         & ! Ptrs to start-end of lst of external neighbours per itfc elem
                                       lst_ext_neighs_gids,                  & ! List of GIDs of external neighbours
                                       lst_ext_neighs_part_ids )               ! Part IDs the external neighbours are mapped to
    implicit none
    class(cell_import_t)  , intent(inout) :: this
    integer(ip)              , intent(in)    :: part_id
    integer(ip)              , intent(in)    :: num_parts
    integer(ip)              , intent(in)    :: num_elements
    integer(ip)              , intent(in)    :: num_itfc_elems
    integer(ip)              , intent(in)    :: lst_itfc_elems(num_itfc_elems)
    integer(ip)              , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(ip)              , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    integer(ip)              , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    
    integer(igp), allocatable :: lst_ext_neighs_gids_igp(:)
    
    call memalloc ( ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1, &
                    lst_ext_neighs_gids_igp, &
                    __FILE__, __LINE__ )
    
    lst_ext_neighs_gids_igp = lst_ext_neighs_gids
    
    call this%cell_import_create_igp ( part_id, & 
                                          num_parts, & 
                                          num_elements, & 
                                          num_itfc_elems, & 
                                          lst_itfc_elems, & 
                                          ptr_ext_neighs_x_itfc_elem, & 
                                          lst_ext_neighs_gids_igp, & 
                                          lst_ext_neighs_part_ids )
    
    call memfree ( lst_ext_neighs_gids_igp, __FILE__, __LINE__)
  end subroutine cell_import_create_ip
  
  subroutine cell_import_compute_num_neighbours ( this, &
                                                        num_itfc_elems, &
                                                        ptr_ext_neighs_x_itfc_elem, &
                                                        lst_ext_neighs_part_ids)
    implicit none
    class(cell_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: num_itfc_elems
    integer(ip)            , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(ip)            , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)

    ! Locals
    type(hash_table_ip_ip_t)   :: parts_visited
    integer(ip)                :: i, istat
    
    ! *WARNING*: We need an estimation for the number of neighbours
    !            in place of a hard-coded, magic number!
    call parts_visited%init(20)

    this%num_neighbours = 0
    do i=1, ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1
       call parts_visited%put(key=lst_ext_neighs_part_ids(i),val=1,stat=istat)
       if(istat==now_stored) this%num_neighbours = this%num_neighbours + 1 
    end do
    call parts_visited%free()

  end subroutine cell_import_compute_num_neighbours

  subroutine cell_import_compute_neighbour_ids ( this, &
                                                    num_itfc_elems, &
                                                    ptr_ext_neighs_x_itfc_elem, &
                                                    lst_ext_neighs_part_ids)
    implicit none
    class(cell_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: num_itfc_elems
    integer(ip)            , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(ip)            , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)

    ! Locals
    type(hash_table_ip_ip_t)   :: parts_visited
    integer(ip)                :: i, j, istat
    
    call memalloc ( this%num_neighbours, this%neighbours_ids, __FILE__, __LINE__ )
    call parts_visited%init(this%num_neighbours)
    j = 1
    do i=1, ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1
       call parts_visited%put(key=lst_ext_neighs_part_ids(i),val=1,stat=istat)
       if(istat==now_stored) then
          this%neighbours_ids(j) = lst_ext_neighs_part_ids(i)
          j = j + 1
       end if
    end do
    call parts_visited%free()
  end subroutine cell_import_compute_neighbour_ids

  subroutine cell_import_compute_snd_rcv_ptrs (this, &
                                                  num_itfc_elems, &
                                                  lst_itfc_elems, &
                                                  ptr_ext_neighs_x_itfc_elem, &
                                                  lst_ext_neighs_gids, &
                                                  lst_ext_neighs_part_ids)
    implicit none
    class(cell_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: num_itfc_elems
    integer(ip)            , intent(in)    :: lst_itfc_elems(num_itfc_elems)
    integer(ip)            , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(igp)           , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    integer(ip)            , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)

    ! Locals
    integer(ip)                            :: i, j, istat, local_neighbour_id
    type(hash_table_ip_ip_t) , allocatable :: snd_lids_x_proc(:) 
    type(hash_table_igp_ip_t), allocatable :: rcv_gids_x_proc(:)

    allocate(snd_lids_x_proc(this%num_neighbours))
    allocate(rcv_gids_x_proc(this%num_neighbours))

    do local_neighbour_id=1, this%num_neighbours
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_x_proc(local_neighbour_id)%init(max( int( real(num_itfc_elems,rp)*0.1_rp, ip), 5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_x_proc(local_neighbour_id)%init( max ( int( real(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1),rp)*0.05_rp,ip), 5) )
    end do

    call memalloc ( this%num_neighbours+1, this%rcv_ptrs, __FILE__, __LINE__)
    call memalloc ( this%num_neighbours+1, this%snd_ptrs, __FILE__, __LINE__)
    this%snd_ptrs = 0
    this%rcv_ptrs = 0
    
    do i=1, num_itfc_elems
       do j=ptr_ext_neighs_x_itfc_elem(i), ptr_ext_neighs_x_itfc_elem(i+1)-1
          local_neighbour_id = this%get_local_neighbour_id(lst_ext_neighs_part_ids(j))
          call snd_lids_x_proc(local_neighbour_id)%put(key=lst_itfc_elems(i), val=1, stat=istat)
          if ( istat == now_stored ) then
             this%snd_ptrs(local_neighbour_id+1) =  this%snd_ptrs(local_neighbour_id+1) + 1
          end if 
          call rcv_gids_x_proc(local_neighbour_id)%put(key=lst_ext_neighs_gids(j), val=1, stat=istat) 
          if ( istat == now_stored ) then
             this%rcv_ptrs(local_neighbour_id+1) =  this%rcv_ptrs(local_neighbour_id+1) + 1
          end if
       end do
    end do

    this%snd_ptrs(1) = 1
    this%rcv_ptrs(1) = 1
    do local_neighbour_id=1, this%num_neighbours
       call snd_lids_x_proc(local_neighbour_id)%free()
       call rcv_gids_x_proc(local_neighbour_id)%free()
       this%snd_ptrs(local_neighbour_id+1) = this%snd_ptrs(local_neighbour_id) + this%snd_ptrs(local_neighbour_id+1)
       this%rcv_ptrs(local_neighbour_id+1) = this%rcv_ptrs(local_neighbour_id) + this%rcv_ptrs(local_neighbour_id+1)
    end do

  end subroutine cell_import_compute_snd_rcv_ptrs

  subroutine cell_import_compute_snd_rcv_leids (this, &
                                                   num_elements, &
                                                   num_itfc_elems, &
                                                   lst_itfc_elems, &
                                                   ptr_ext_neighs_x_itfc_elem, &
                                                   lst_ext_neighs_gids, &
                                                   lst_ext_neighs_part_ids)
    implicit none
    class(cell_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: num_elements
    integer(ip)            , intent(in)    :: num_itfc_elems
    integer(ip)            , intent(in)    :: lst_itfc_elems(num_itfc_elems)
    integer(ip)            , intent(in)    :: ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)
    integer(igp)           , intent(in)    :: lst_ext_neighs_gids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)
    integer(ip)            , intent(in)    :: lst_ext_neighs_part_ids(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1)-1)

    ! Locals
    integer(ip)                            :: i, j, istat, local_neighbour_id
    integer(ip)                            :: num_elements_to_be_received
    type(hash_table_ip_ip_t) , allocatable :: snd_lids_x_proc(:) 
    type(hash_table_igp_ip_t), allocatable :: rcv_gids_x_proc(:)
    
    allocate(snd_lids_x_proc(this%num_neighbours))
    allocate(rcv_gids_x_proc(this%num_neighbours))
    do i=1, this%num_neighbours
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_x_proc(i)%init(max(int(real(num_itfc_elems,rp)*0.1_rp,ip),5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_x_proc(i)%init(max(int(real(ptr_ext_neighs_x_itfc_elem(num_itfc_elems+1),rp)*0.05_rp,ip),5))
    end do

    call memalloc ( this%rcv_ptrs(this%num_neighbours+1)-1, this%rcv_leids, __FILE__, __LINE__)
    call memalloc ( this%snd_ptrs(this%num_neighbours+1)-1, this%snd_leids, __FILE__, __LINE__)
    
    num_elements_to_be_received = 0
    do i=1, num_itfc_elems
       do j=ptr_ext_neighs_x_itfc_elem(i), ptr_ext_neighs_x_itfc_elem(i+1)-1
          local_neighbour_id = this%get_local_neighbour_id(lst_ext_neighs_part_ids(j))
          call snd_lids_x_proc(local_neighbour_id)%put(key=lst_itfc_elems(i), val=1, stat=istat)
          if ( istat == now_stored ) then
             this%snd_leids(this%snd_ptrs(local_neighbour_id)) = lst_itfc_elems(i)
             this%snd_ptrs(local_neighbour_id) = this%snd_ptrs(local_neighbour_id) + 1
          end if 
          
          call rcv_gids_x_proc(local_neighbour_id)%put(key=lst_ext_neighs_gids(j), val=1, stat=istat) 
          if ( istat == now_stored ) then
             num_elements_to_be_received = num_elements_to_be_received + 1 
             this%rcv_leids(this%rcv_ptrs(local_neighbour_id)) = num_elements + num_elements_to_be_received
             this%rcv_ptrs(local_neighbour_id) = this%rcv_ptrs(local_neighbour_id) + 1
          end if
       end do
    end do

    do i=this%num_neighbours+1, 2, -1
       call snd_lids_x_proc(i-1)%free()
       call rcv_gids_x_proc(i-1)%free()
       this%snd_ptrs(i) = this%snd_ptrs(i-1)
       this%rcv_ptrs(i) = this%rcv_ptrs(i-1)
    end do
    this%snd_ptrs(1) = 1
    this%rcv_ptrs(1) = 1
  end subroutine cell_import_compute_snd_rcv_leids
  
  !=============================================================================
  subroutine cell_import_free ( this )
    implicit none
    class(cell_import_t), intent(inout) :: this
    if (allocated(this%neighbours_ids)) call memfree ( this%neighbours_ids,__FILE__,__LINE__)
    if (allocated(this%rcv_ptrs)) call memfree ( this%rcv_ptrs,__FILE__,__LINE__)
    if (allocated(this%rcv_leids)) call memfree ( this%rcv_leids,__FILE__,__LINE__)
    if (allocated(this%snd_ptrs)) call memfree ( this%snd_ptrs,__FILE__,__LINE__)
    if (allocated(this%snd_leids)) call memfree ( this%snd_leids,__FILE__,__LINE__)
    this%part_id           = -1
    this%num_parts      = -1
    this%num_ghost_elements     = -1
    this%num_neighbours = -1
  end subroutine cell_import_free

  !=============================================================================
  subroutine cell_import_print (this, lu_out)
    !-----------------------------------------------------------------------
    ! This routine prints an cell_import object
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    class(cell_import_t), intent(in)  :: this
    integer(ip)            , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then
       write(lu_out,'(a)') '*** begin cell_import_t data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
            &  this%num_parts

       write(lu_out,'(a,i10)') 'Number of neighbours:', &
            &  this%num_neighbours

       write(lu_out,'(a,i10)') 'Number of ghost elements:', &
            &  this%num_ghost_elements

       write(lu_out,'(a)') 'List of neighbours:'
       write(lu_out,'(10i10)') this%neighbours_ids(1:this%num_neighbours)

       write(lu_out,'(a)') 'Rcv_ptrs:'
       write(lu_out,'(10i10)') this%rcv_ptrs(1:this%num_neighbours+1)
       write(lu_out,'(a)') 'Rcv_leids:'
       write(lu_out,'(10i10)') this%rcv_leids
       

       write(lu_out,'(a)') 'Snd_ptrs:'
       write(lu_out,'(10i10)') this%snd_ptrs(1:this%num_neighbours+1)
       write(lu_out,'(a)') 'Snd_leids:'
       write(lu_out,'(10i10)') this%snd_leids
    end if
  end subroutine cell_import_print

  function cell_import_get_num_parts ( this )
    implicit  none
    class(cell_import_t), intent(in) :: this
    integer(ip) :: cell_import_get_num_parts
    cell_import_get_num_parts = this%num_parts
  end function cell_import_get_num_parts
  
  pure function cell_import_get_num_ghost_elements ( this )
    implicit none
    class(cell_import_t), intent(in) :: this
    integer(ip)                      :: cell_import_get_num_ghost_elements
    cell_import_get_num_ghost_elements = this%num_ghost_elements
  end function cell_import_get_num_ghost_elements
  
  pure function cell_import_get_num_neighbours ( this )
    implicit none
    class(cell_import_t), intent(in) :: this
    integer(ip)                      :: cell_import_get_num_neighbours
    cell_import_get_num_neighbours = this%num_neighbours
  end function cell_import_get_num_neighbours
  
  function cell_import_get_neighbours_ids ( this )
    implicit none
    class(cell_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: cell_import_get_neighbours_ids(:)
    cell_import_get_neighbours_ids => this%neighbours_ids
  end function cell_import_get_neighbours_ids
  
  function cell_import_get_rcv_ptrs ( this )
    implicit none
    class(cell_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: cell_import_get_rcv_ptrs(:)
    cell_import_get_rcv_ptrs => this%rcv_ptrs
  end function cell_import_get_rcv_ptrs
 
  function cell_import_get_rcv_leids ( this )
    implicit none
    class(cell_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: cell_import_get_rcv_leids(:)
    cell_import_get_rcv_leids => this%rcv_leids
  end function cell_import_get_rcv_leids
  
  function cell_import_get_snd_ptrs ( this )
    implicit none
    class(cell_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: cell_import_get_snd_ptrs(:)
    cell_import_get_snd_ptrs => this%snd_ptrs
  end function cell_import_get_snd_ptrs
  
  function cell_import_get_snd_leids ( this )
    implicit none
    class(cell_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: cell_import_get_snd_leids(:)
    cell_import_get_snd_leids => this%snd_leids
  end function cell_import_get_snd_leids
  
  function cell_import_get_global_neighbour_id ( this, local_neighbour_id )
    implicit none
    class(cell_import_t), intent(in) :: this
    integer(ip)         , intent(in) :: local_neighbour_id
    integer(ip)                      :: cell_import_get_global_neighbour_id
    assert ( local_neighbour_id >=1 .and. local_neighbour_id <= this%num_neighbours )
    cell_import_get_global_neighbour_id = this%neighbours_ids(local_neighbour_id)
  end function cell_import_get_global_neighbour_id
  
  function cell_import_get_local_neighbour_id ( this, global_neighbour_id )
    implicit none
    class(cell_import_t), intent(in) :: this
    integer(ip)            , intent(in) :: global_neighbour_id
    integer(ip)                      :: cell_import_get_local_neighbour_id
        
    do cell_import_get_local_neighbour_id = 1, this%num_neighbours
      if ( this%neighbours_ids(cell_import_get_local_neighbour_id) == global_neighbour_id ) return
    end do
    assert ( .not. (cell_import_get_local_neighbour_id == this%num_neighbours+1) )
  end function cell_import_get_local_neighbour_id
  
end module cell_import_names
