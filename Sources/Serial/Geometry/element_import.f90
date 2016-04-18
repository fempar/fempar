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
module element_import_names
  use types_names
  use memor_names
  use hash_table_names
  
# include "debug.i90"  
  implicit none
  private

  ! Element_Import
  ! Host for data needed to perform nearest neighbour communications for
  ! the distributed dual graph associated to the triangulation
  type element_import_t
     private
     integer(ip)               :: part_id
     integer(ip)               :: number_parts
     integer(ip)               :: number_ghost_elements    
     integer(ip)               :: number_neighbours          
     integer(ip), allocatable  :: neighbours_ids(:)          
     integer(ip), allocatable  :: rcv_ptrs(:)                ! How many elements does the part receive from each neighbour ?
     integer(ip), allocatable  :: rcv_leids(:)               ! Local position in the elements array in which the elements to be received by this part  
                                                             ! are going to be placed
     integer(ip), allocatable  :: snd_ptrs(:)                ! How many elements does the part send to each neighbour?
     integer(ip), allocatable  :: snd_leids(:)               ! Local elements IDs of the elements to be sent to each neighbour ?
  contains
     procedure, non_overridable :: create                                 => element_import_create
     procedure, private, non_overridable :: compute_number_neighbours     => element_import_compute_number_neighbours
     procedure, private, non_overridable :: compute_neighbours_ids        => element_import_compute_neighbour_ids
     procedure, private, non_overridable :: compute_snd_rcv_ptrs          => element_import_compute_snd_rcv_ptrs
     procedure, private, non_overridable :: compute_snd_rcv_leids         => element_import_compute_snd_rcv_leids
     procedure, non_overridable :: free                                   => element_import_free
     procedure, non_overridable :: print                                  => element_import_print
     procedure, non_overridable :: get_number_ghost_elements              => element_import_get_number_ghost_elements
     procedure, non_overridable :: get_number_neighbours                  => element_import_get_number_neighbours
     procedure, non_overridable :: get_neighbours_ids                     => element_import_get_neighbours_ids
     procedure, non_overridable :: get_rcv_ptrs                           => element_import_get_rcv_ptrs
     procedure, non_overridable :: get_rcv_leids                          => element_import_get_rcv_leids
     procedure, non_overridable :: get_snd_ptrs                           => element_import_get_snd_ptrs
     procedure, non_overridable :: get_snd_leids                          => element_import_get_snd_leids
     procedure, non_overridable :: get_global_neighbour_id                => element_import_get_global_neighbour_id
     procedure, non_overridable :: get_local_neighbour_id                 => element_import_get_local_neighbour_id
  end type element_import_t

  ! Types
  public :: element_import_t

  ! Functions
  public :: element_import_free, element_import_print

contains

  subroutine element_import_create (this,            & ! Passed-object dummy argument
                                    part_id,         & ! My part identifier
                                    number_parts,    & ! Number of parts
                                    number_elements, & ! Number of local elements
                                    nebou,           & ! Number of (local) itfc elements
                                    lebou,           & ! Local IDs of itfc elements
                                    pextn,           & ! Pointers to start-end of the list of the external neighbours of each itfc element
                                    lextn,           & ! List of GIDs of external neighbours
                                    lextp )            ! Part IDs the external neighbours are mapped to
    implicit none
    class(element_import_t)  , intent(inout) :: this
    integer(ip)              , intent(in)    :: part_id
    integer(ip)              , intent(in)    :: number_parts
    integer(ip)              , intent(in)    :: number_elements
    integer(ip)              , intent(in)    :: nebou
    integer(ip)              , intent(in)    :: lebou(nebou)
    integer(ip)              , intent(in)    :: pextn(nebou+1)
    integer(igp)             , intent(in)    :: lextn(pextn(nebou+1)-1)
    integer(ip)              , intent(in)    :: lextp(pextn(nebou+1)-1)
    
    call this%free()    
    this%part_id = part_id
    this%number_parts  = number_parts
    call this%compute_number_neighbours( nebou, pextn, lextp )
    call this%compute_neighbours_ids ( nebou, pextn, lextp )
    call this%compute_snd_rcv_ptrs (nebou, lebou, pextn, lextn, lextp )
    call this%compute_snd_rcv_leids (number_elements, nebou, lebou, pextn, lextn, lextp)
    this%number_ghost_elements = this%rcv_ptrs(this%number_neighbours+1)-1    
  end subroutine element_import_create

  subroutine element_import_compute_number_neighbours ( this, nebou, pextn, lextp)
    implicit none
    class(element_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: nebou
    integer(ip)            , intent(in)    :: pextn(nebou+1)
    integer(ip)            , intent(in)    :: lextp(pextn(nebou+1)-1)

    ! Locals
    type(hash_table_ip_ip_t)   :: parts_visited
    integer(ip)                :: i, istat
    
    ! *WARNING*: We need an estimation for the number of neighbours
    !            in place of a hard-coded, magic number!
    call parts_visited%init(20)

    this%number_neighbours = 0
    do i=1, pextn(nebou+1)-1
       call parts_visited%put(key=lextp(i),val=1,stat=istat)
       if(istat==now_stored) this%number_neighbours = this%number_neighbours + 1 
    end do
    call parts_visited%free()

  end subroutine element_import_compute_number_neighbours

  subroutine element_import_compute_neighbour_ids ( this, nebou, pextn, lextp)
    implicit none
    class(element_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: nebou
    integer(ip)            , intent(in)    :: pextn(nebou+1)
    integer(ip)            , intent(in)    :: lextp(pextn(nebou+1)-1)

    ! Locals
    type(hash_table_ip_ip_t)   :: parts_visited
    integer(ip)                :: i, j, istat
    
    call memalloc ( this%number_neighbours, this%neighbours_ids, __FILE__, __LINE__ )
    call parts_visited%init(this%number_neighbours)
    j = 1
    do i=1, pextn(nebou+1)-1
       call parts_visited%put(key=lextp(i),val=1,stat=istat)
       if(istat==now_stored) then
          this%neighbours_ids(j) = lextp(i)
          j = j + 1
       end if
    end do
    call parts_visited%free()
  end subroutine element_import_compute_neighbour_ids

  subroutine element_import_compute_snd_rcv_ptrs (this, nebou, lebou, pextn, lextn, lextp)
    implicit none
    class(element_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: nebou
    integer(ip)            , intent(in)    :: lebou(nebou)
    integer(ip)            , intent(in)    :: pextn(nebou+1)
    integer(igp)           , intent(in)    :: lextn(pextn(nebou+1)-1)
    integer(ip)            , intent(in)    :: lextp(pextn(nebou+1)-1)

    ! Locals
    integer(ip)                            :: i, j, istat, local_neighbour_id
    type(hash_table_ip_ip_t) , allocatable :: snd_lids_per_proc(:) 
    type(hash_table_igp_ip_t), allocatable :: rcv_gids_per_proc(:)

    allocate(snd_lids_per_proc(this%number_neighbours))
    allocate(rcv_gids_per_proc(this%number_neighbours))

    do local_neighbour_id=1, this%number_neighbours
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_per_proc(local_neighbour_id)%init(max( int( real(nebou,rp)*0.1_rp, ip), 5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_per_proc(local_neighbour_id)%init( max ( int( real(pextn(nebou+1),rp)*0.05_rp,ip), 5) )
    end do

    call memalloc ( this%number_neighbours+1, this%rcv_ptrs, __FILE__, __LINE__)
    call memalloc ( this%number_neighbours+1, this%snd_ptrs, __FILE__, __LINE__)
    this%snd_ptrs = 0
    this%rcv_ptrs = 0
    
    do i=1, nebou
       do j=pextn(i), pextn(i+1)-1
          local_neighbour_id = this%get_local_neighbour_id(lextp(j))
          call snd_lids_per_proc(local_neighbour_id)%put(key=lebou(i), val=1, stat=istat)
          if ( istat == now_stored ) then
             this%snd_ptrs(local_neighbour_id+1) =  this%snd_ptrs(local_neighbour_id+1) + 1
          end if 
          call rcv_gids_per_proc(local_neighbour_id)%put(key=lextn(j), val=1, stat=istat) 
          if ( istat == now_stored ) then
             this%rcv_ptrs(local_neighbour_id+1) =  this%rcv_ptrs(local_neighbour_id+1) + 1
          end if
       end do
    end do

    this%snd_ptrs(1) = 1
    this%rcv_ptrs(1) = 1
    do local_neighbour_id=1, this%number_neighbours
       call snd_lids_per_proc(local_neighbour_id)%free()
       call rcv_gids_per_proc(local_neighbour_id)%free()
       this%snd_ptrs(local_neighbour_id+1) = this%snd_ptrs(local_neighbour_id) + this%snd_ptrs(local_neighbour_id+1)
       this%rcv_ptrs(local_neighbour_id+1) = this%rcv_ptrs(local_neighbour_id) + this%rcv_ptrs(local_neighbour_id+1)
    end do

  end subroutine element_import_compute_snd_rcv_ptrs

  subroutine element_import_compute_snd_rcv_leids (this, number_elements, nebou, lebou, pextn, lextn, lextp)
    implicit none
    class(element_import_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: number_elements
    integer(ip)            , intent(in)    :: nebou
    integer(ip)            , intent(in)    :: lebou(nebou)
    integer(ip)            , intent(in)    :: pextn(nebou+1)
    integer(igp)           , intent(in)    :: lextn(pextn(nebou+1)-1)
    integer(ip)            , intent(in)    :: lextp(pextn(nebou+1)-1)

    ! Locals
    integer(ip)                            :: i, j, istat, local_neighbour_id
    integer(ip)                            :: number_elements_to_be_received
    type(hash_table_ip_ip_t) , allocatable :: snd_lids_per_proc(:) 
    type(hash_table_igp_ip_t), allocatable :: rcv_gids_per_proc(:)
    
    allocate(snd_lids_per_proc(this%number_neighbours))
    allocate(rcv_gids_per_proc(this%number_neighbours))
    do i=1, this%number_neighbours
       ! Assume that nebou elements will be sent to each processor, 
       ! and take its 10% for the size of the hash table
       call snd_lids_per_proc(i)%init(max(int(real(nebou,rp)*0.1_rp,ip),5))

       ! Assume that as many elements will be received from each processor
       ! as the number of remote edges communicating elements in this part
       ! and any other remote part, and take its 5%
       call rcv_gids_per_proc(i)%init(max(int(real(pextn(nebou+1),rp)*0.05_rp,ip),5))
    end do

    call memalloc ( this%rcv_ptrs(this%number_neighbours+1)-1, this%rcv_leids, __FILE__, __LINE__)
    call memalloc ( this%snd_ptrs(this%number_neighbours+1)-1, this%snd_leids, __FILE__, __LINE__)
    
    number_elements_to_be_received = 0
    do i=1, nebou
       do j=pextn(i), pextn(i+1)-1
          local_neighbour_id = this%get_local_neighbour_id(lextp(j))
          call snd_lids_per_proc(local_neighbour_id)%put(key=lebou(i), val=1, stat=istat)
          if ( istat == now_stored ) then
             this%snd_leids(this%snd_ptrs(local_neighbour_id)) = lebou(i)
             this%snd_ptrs(local_neighbour_id) = this%snd_ptrs(local_neighbour_id) + 1
          end if 
          
          call rcv_gids_per_proc(local_neighbour_id)%put(key=lextn(j), val=1, stat=istat) 
          if ( istat == now_stored ) then
             number_elements_to_be_received = number_elements_to_be_received + 1 
             this%rcv_leids(this%rcv_ptrs(local_neighbour_id)) = number_elements + number_elements_to_be_received
             this%rcv_ptrs(local_neighbour_id) = this%rcv_ptrs(local_neighbour_id) + 1
          end if
       end do
    end do

    do i=this%number_neighbours+1, 2, -1
       call snd_lids_per_proc(i-1)%free()
       call rcv_gids_per_proc(i-1)%free()
       this%snd_ptrs(i) = this%snd_ptrs(i-1)
       this%rcv_ptrs(i) = this%rcv_ptrs(i-1)
    end do
    this%snd_ptrs(1) = 1
    this%rcv_ptrs(1) = 1
  end subroutine element_import_compute_snd_rcv_leids
  
  !=============================================================================
  subroutine element_import_free ( this )
    implicit none
    class(element_import_t), intent(inout) :: this
    if (allocated(this%neighbours_ids)) call memfree ( this%neighbours_ids,__FILE__,__LINE__)
    if (allocated(this%rcv_ptrs)) call memfree ( this%rcv_ptrs,__FILE__,__LINE__)
    if (allocated(this%rcv_leids)) call memfree ( this%rcv_leids,__FILE__,__LINE__)
    if (allocated(this%snd_ptrs)) call memfree ( this%snd_ptrs,__FILE__,__LINE__)
    if (allocated(this%snd_leids)) call memfree ( this%snd_leids,__FILE__,__LINE__)
    this%part_id           = -1
    this%number_parts      = -1
    this%number_ghost_elements     = -1
    this%number_neighbours = -1
  end subroutine element_import_free

  !=============================================================================
  subroutine element_import_print (this, lu_out)
    !-----------------------------------------------------------------------
    ! This routine prints an element_import object
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    class(element_import_t), intent(in)  :: this
    integer(ip)            , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then
       write(lu_out,'(a)') '*** begin element_import_t data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
            &  this%number_parts

       write(lu_out,'(a,i10)') 'Number of neighbours:', &
            &  this%number_neighbours

       write(lu_out,'(a,i10)') 'Number of ghost elements:', &
            &  this%number_ghost_elements

       write(lu_out,'(a)') 'List of neighbours:'
       write(lu_out,'(10i10)') this%neighbours_ids(1:this%number_neighbours)

       write(lu_out,'(a)') 'Rcv_ptrs:'
       write(lu_out,'(10i10)') this%rcv_ptrs(1:this%number_neighbours+1)
       write(lu_out,'(a)') 'Rcv_leids:'
       write(lu_out,'(10i10)') this%rcv_leids
       

       write(lu_out,'(a)') 'Snd_ptrs:'
       write(lu_out,'(10i10)') this%snd_ptrs(1:this%number_neighbours+1)
       write(lu_out,'(a)') 'Snd_leids:'
       write(lu_out,'(10i10)') this%snd_leids
    end if
  end subroutine element_import_print

  pure function element_import_get_number_ghost_elements ( this )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)                      :: element_import_get_number_ghost_elements
    element_import_get_number_ghost_elements = this%number_ghost_elements
  end function element_import_get_number_ghost_elements
  
  pure function element_import_get_number_neighbours ( this )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)                      :: element_import_get_number_neighbours
    element_import_get_number_neighbours = this%number_neighbours
  end function element_import_get_number_neighbours
  
  function element_import_get_neighbours_ids ( this )
    implicit none
    class(element_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: element_import_get_neighbours_ids(:)
    element_import_get_neighbours_ids => this%neighbours_ids
  end function element_import_get_neighbours_ids
  
  function element_import_get_rcv_ptrs ( this )
    implicit none
    class(element_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: element_import_get_rcv_ptrs(:)
    element_import_get_rcv_ptrs => this%rcv_ptrs
  end function element_import_get_rcv_ptrs
 
  function element_import_get_rcv_leids ( this )
    implicit none
    class(element_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: element_import_get_rcv_leids(:)
    element_import_get_rcv_leids => this%rcv_leids
  end function element_import_get_rcv_leids
  
  function element_import_get_snd_ptrs ( this )
    implicit none
    class(element_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: element_import_get_snd_ptrs(:)
    element_import_get_snd_ptrs => this%snd_ptrs
  end function element_import_get_snd_ptrs
  
  function element_import_get_snd_leids ( this )
    implicit none
    class(element_import_t), target, intent(in) :: this
    integer(ip), pointer                        :: element_import_get_snd_leids(:)
    element_import_get_snd_leids => this%snd_leids
  end function element_import_get_snd_leids
  
  function element_import_get_global_neighbour_id ( this, local_neighbour_id )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)         , intent(in) :: local_neighbour_id
    integer(ip)                      :: element_import_get_global_neighbour_id
    assert ( local_neighbour_id >=1 .and. local_neighbour_id <= this%number_neighbours )
    element_import_get_global_neighbour_id = this%neighbours_ids(local_neighbour_id)
  end function element_import_get_global_neighbour_id
  
  function element_import_get_local_neighbour_id ( this, global_neighbour_id )
    implicit none
    class(element_import_t), intent(in) :: this
    integer(ip)            , intent(in) :: global_neighbour_id
    integer(ip)                      :: element_import_get_local_neighbour_id
        
    do element_import_get_local_neighbour_id = 1, this%number_neighbours
      if ( this%neighbours_ids(element_import_get_local_neighbour_id) == global_neighbour_id ) return
    end do
    assert ( .not. (element_import_get_local_neighbour_id == this%number_neighbours+1) )
  end function element_import_get_local_neighbour_id
  
end module element_import_names
