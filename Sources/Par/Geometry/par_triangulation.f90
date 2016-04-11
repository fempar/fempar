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
module par_triangulation_names
  ! Serial modules
  use types_names
  use list_types_names
  use memor_names
  use sort_names
  use migratory_element_names
  use triangulation_names
  use element_import_names
  use hash_table_names

  ! Parallel modules
  use par_context_names
  use par_environment_names

  implicit none
# include "debug.i90"
  private

  ! AFM: To think. Is it required to know at some point within FEMPAR the number of global vefs ?
  !                With the current design, it is known the number of local vefs
  !                of each subdomain, but not the number of global vefs.
  !                If it would be required, one possibility would be to store it at type(par_triangulation).

  !                On the other hand, would it be necessary/convenient/helpful to store 
  !                the number of interface vefs, and the list of interface vefs local IDs?
  !                At this point, I know that this would be convenient/helpful for par_dof_descriptor_partition.f90.
  !                (That requires fast location of interface vefs). Any other point?

  !                I already did store these data on num_itfc_objs and lst_itfc_objs members of type(par_triangulation)
  !                but at the current state I do not know whether this is the most appropriate place. 

  integer(ip), parameter :: par_triangulation_not_created  = 0 ! Initial state
  integer(ip), parameter :: par_triangulation_filled       = 1 ! Elems + Vefs arrays allocated and filled 


  type par_vef_topology_t
     integer(ip)  :: interface  = -1 ! Interface local id of this vef
     integer(igp) :: globalID   = -1 ! Global ID of this vef
                                     ! Local ID is the position in the array of vefs
  end type par_vef_topology_t

  type, extends(migratory_element_t) :: par_elem_topology_t
     integer(ip)  :: interface  = -1              ! The boundary number ieboun (if this element is a interface element)
     integer(ip)  :: mypart     = -1              ! To which part this element is mapped to ?
     integer(igp) :: globalID   = -1              ! Global ID of this element
                                                  ! Local ID is the position in the array of elements
     integer(ip)  :: num_vefs = -1             ! Number of vefs
     integer(igp), allocatable :: vefs_GIDs(:) ! List of the GIDs of the vefs that make up this element
     
   contains
     procedure :: size   => par_elem_topology_size
     procedure :: pack   => par_elem_topology_pack
     procedure :: unpack => par_elem_topology_unpack
  end type par_elem_topology_t

  type :: par_triangulation_t
     integer(ip)                              :: state = par_triangulation_not_created  
     type(triangulation_t)                    :: triangulation             ! Data common with a centralized (serial) triangulation
     integer(ip)                              :: num_elems   = -1    ! should it match f_trian%num_elems or not? 
     integer(ip)                              :: num_ghosts  = -1    ! number of ghost elements (remote neighbors)
     class(migratory_element_t), allocatable  :: mig_elems(:)        ! Migratory elements list_t
     type(par_elem_topology_t),   pointer     :: elems(:) => NULL()  ! Array of elements in the mesh
     type(par_vef_topology_t), allocatable    :: vefs(:)          ! array of vefs in the mesh
     integer(ip)                              :: num_itfc_vefs = -1  ! Number of vefs in the interface among subdomains
     integer(ip), allocatable                 :: lst_itfc_vefs(:)    ! List of vefs local IDs in the interface among subdomains 
     integer(ip)                              :: num_itfc_elems= -1  ! Number of elements in the interface
     integer(ip), allocatable                 :: lst_itfc_elems(:)   ! List of elements local IDs in the interface
     type(par_environment_t),   pointer       :: p_env => NULL()     ! Parallel environment describing MPI tasks among which par_triangulation is distributed
     type(element_import_t)                   :: element_import         ! Data type describing the layout in distributed-memory of the dual graph
                                                                    ! (It is required, e.g., for nearest neighbour comms on this graph)
     
     ! Perhaps the following three member variables should be packed within type(map_t) ?
     ! I didn't do that because type(map_t) has extra members that do not make sense anymore
     ! for the current situation with objects (i.e., interior, boundary, external) etc.
     integer(ip)                             :: number_global_objects
     integer(ip)                             :: number_objects
     integer(ip), allocatable                :: objects_gids(:)
     
     type(list_t)                            :: vefs_object
     type(list_t)                            :: parts_object
  contains
     procedure, private :: compute_parts_itfc_vefs                      => par_triangulation_compute_parts_itfc_vefs
     procedure, private :: compute_vefs_and_parts_object                => par_triangulation_compute_vefs_and_parts_object
     procedure, private :: compute_objects_neighbours_exchange_data     => par_triangulation_compute_objects_neighbours_exchange_data
     procedure, private :: compute_number_global_objects_and_their_gids => par_triangulation_compute_number_global_objects_and_their_gids
  end type par_triangulation_t

  ! Types
  public :: par_triangulation_t, par_elem_topology_t

  ! Functions
  public :: par_triangulation_free, par_triangulation_print, par_triangulation_to_dual

  ! Auxiliary Subroutines (should only be used by modules that have control over type(par_triangulation))
  public :: free_par_elem_topology, free_par_vef_topology, par_triangulation_free_elems_data, par_triangulation_free_objs_data 

  ! Constants (should only be used by modules that have control over type(triangulation_t))
  public :: par_triangulation_not_created, par_triangulation_filled

contains

  subroutine par_triangulation_print ( lunou,  p_trian ) ! (SBmod)
    implicit none
    ! Parameters
    integer(ip)              , intent(in) :: lunou
    type(par_triangulation_t), intent(in) :: p_trian
    if(p_trian%p_env%am_i_l1_task()) then
       call triangulation_print ( lunou, p_trian%triangulation, p_trian%num_elems + p_trian%num_ghosts ) 
    end if

  end subroutine par_triangulation_print


  !=============================================================================
  subroutine par_triangulation_free(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: iobj, istat, ielem

    if(.not. p_trian%p_env%am_i_l1_task()) return

    assert(p_trian%state == par_triangulation_filled) 
    
    call par_triangulation_free_objs_data (p_trian)
    call par_triangulation_free_elems_data(p_trian)
    
    ! Free local portion of triangulation
    call triangulation_free(p_trian%triangulation)

    p_trian%state = par_triangulation_not_created

  end subroutine par_triangulation_free

  !=============================================================================
  subroutine par_triangulation_free_objs_data(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: iobj, istat

    if ( p_trian%state == par_triangulation_filled ) then
       do iobj=1, p_trian%triangulation%num_vefs 
          call free_par_vef_topology(p_trian%vefs(iobj)) 
       end do
       
       p_trian%num_itfc_vefs = -1
       call memfree ( p_trian%lst_itfc_vefs, __FILE__, __LINE__ )
        
       p_trian%number_objects        = -1
       call p_trian%vefs_object%free()
       call p_trian%parts_object%free()
       
       p_trian%number_global_objects = -1
       call memfree ( p_trian%objects_gids, __FILE__, __LINE__ )
       
       ! Deallocate the vef structure array 
       deallocate(p_trian%vefs, stat=istat)
       check(istat==0)
    end if
  end subroutine par_triangulation_free_objs_data

  !=============================================================================
  subroutine par_triangulation_free_elems_data(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: ielem, istat

    if ( p_trian%state == par_triangulation_filled ) then
       do ielem=1, p_trian%triangulation%elem_array_len 
          call free_par_elem_topology(p_trian%elems(ielem)) 
       end do
       
       p_trian%num_itfc_elems = -1
       call memfree ( p_trian%lst_itfc_elems, __FILE__, __LINE__ )
       
       call element_import_free ( p_trian%element_import )
       
       ! Deallocate the element structure array */
       deallocate(p_trian%mig_elems, stat=istat)
       check(istat==0)
       nullify(p_trian%elems)
       
       p_trian%num_elems  = -1 
       p_trian%num_ghosts = -1
       nullify(p_trian%p_env)
    end if

  end subroutine par_triangulation_free_elems_data

  !=============================================================================
  subroutine par_triangulation_to_dual(p_trian)  
    implicit none
    ! Parameters
    type(par_triangulation_t), intent(inout) :: p_trian

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat

    call triangulation_to_dual ( p_trian%triangulation, p_trian%num_elems+p_trian%num_ghosts )

    ! Allocate p_trian%vefs array
    allocate(p_trian%vefs(p_trian%triangulation%num_vefs), stat=istat)
    check(istat==0)

    ! Initialize par_vef_topology vefs
    do iobj=1, p_trian%triangulation%num_vefs 
       call initialize_par_vef_topology(p_trian%vefs(iobj)) 
    end do

    ! Assign global ID for all local vefs
    do ielem=1, p_trian%triangulation%num_elems
       do iobj=1, p_trian%elems(ielem)%num_vefs
          jobj = p_trian%triangulation%elems(ielem)%vefs(iobj)
          if ( jobj /= -1 ) then
             ! write(*,*) 'ZZZ', p_trian%elems(ielem)%vefs_GIDs(iobj)
             p_trian%vefs(jobj)%globalID = p_trian%elems(ielem)%vefs_GIDs(iobj)
          end if
       end do
    end do

    ! Count number of interface vefs
    p_trian%num_itfc_vefs = 0
    ! Traverse local view of ghost elements (all the interface vefs are there!!!)
    do ielem=p_trian%num_elems+1, p_trian%num_ghosts+p_trian%num_elems
       do iobj=1, p_trian%triangulation%elems(ielem)%num_vefs
          jobj = p_trian%triangulation%elems(ielem)%vefs(iobj)
          if ( jobj /= -1 ) then
             if (p_trian%vefs(jobj)%interface == -1) then
                p_trian%num_itfc_vefs = p_trian%num_itfc_vefs + 1
                p_trian%vefs(jobj)%interface = p_trian%num_itfc_vefs
             end if
          end if
       end do
    end do

    ! List interface vefs
    call memalloc ( p_trian%num_itfc_vefs, p_trian%lst_itfc_vefs, __FILE__, __LINE__ )
    p_trian%num_itfc_vefs = 0 
    ! Traverse local view of ghost elements (all the interface vefs are there!!!)
    do ielem=p_trian%num_elems+1, p_trian%num_ghosts+p_trian%num_elems
       do iobj=1, p_trian%triangulation%elems(ielem)%num_vefs
          jobj = p_trian%triangulation%elems(ielem)%vefs(iobj)
          if ( jobj /= -1 ) then
             if (p_trian%vefs(jobj)%interface == (p_trian%num_itfc_vefs + 1)) then
                p_trian%num_itfc_vefs = p_trian%num_itfc_vefs + 1
                p_trian%lst_itfc_vefs(p_trian%num_itfc_vefs) = jobj 
             end if
          end if
       end do
    end do

    ! Compute p_trian%max_nparts, p_trian%nobjs, p_trian%lobjs
    ! Re-order (permute) p_trian%lst_itfc_vefs accordingly
    call p_trian%compute_vefs_and_parts_object()
    call p_trian%compute_number_global_objects_and_their_gids()
  end subroutine par_triangulation_to_dual
 
  subroutine par_triangulation_compute_parts_itfc_vefs ( this, parts_itfc_vefs, perm_itfc_vefs )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip), allocatable  , intent(inout) :: parts_itfc_vefs(:,:)
    integer(ip), allocatable  , intent(inout) :: perm_itfc_vefs(:)
    
    integer(ip)                               :: num_neighbours
    logical, allocatable                      :: touched_neighbours(:)
    integer(ip)                               :: nparts_around, mypart_id, part_id, local_part_id
    integer(ip)                               :: ivef_itfc, ielem, vef_lid
    integer(ip)                               :: elem_lid
    integer(ip)                               :: num_rows_parts_itfc_vefs
    integer(ip), allocatable                  :: work1(:), work2(:)
    
    assert ( this%p_env%am_i_l1_task() )
    
    if (allocated(parts_itfc_vefs)) call memfree(parts_itfc_vefs,__FILE__,__LINE__)
    if (allocated(perm_itfc_vefs)) call memfree(perm_itfc_vefs,__FILE__,__LINE__)
    
    mypart_id = this%p_env%get_l1_rank() + 1 
   
    num_neighbours = this%element_import%get_number_neighbours()    
    call memalloc ( num_neighbours, touched_neighbours, __FILE__, __LINE__ )
    
    ! The two extra rows in parts_per_itfc_vef are required in order to: (1) hold the number of parts around an interface vef
    !                                                                    (2) to hold mypart_id, which should be also listed among 
    !                                                                        the parts around each vef
    num_rows_parts_itfc_vefs = num_neighbours + 2
    call memalloc ( num_rows_parts_itfc_vefs, this%num_itfc_vefs, parts_itfc_vefs, __FILE__, __LINE__ )
    parts_itfc_vefs = 0
    
    do ivef_itfc = 1, this%num_itfc_vefs
      vef_lid = this%lst_itfc_vefs(ivef_itfc)
      touched_neighbours = .false.
      
      nparts_around = 1 
      parts_itfc_vefs(nparts_around+1,ivef_itfc) = mypart_id
      
      do ielem=1, this%triangulation%vefs(vef_lid)%num_elems_around
        elem_lid = this%triangulation%vefs(vef_lid)%elems_around(ielem)
        part_id = this%elems(elem_lid)%mypart
        
        if ( part_id /= mypart_id ) then
         local_part_id = this%element_import%get_local_neighbour_id(part_id)
         if (.not. touched_neighbours (local_part_id)) then
           touched_neighbours (local_part_id) = .true.
           nparts_around = nparts_around + 1 
           parts_itfc_vefs(nparts_around+1,ivef_itfc) = part_id
         end if
        end if
      end do
      parts_itfc_vefs(1,ivef_itfc) = nparts_around
      ! Sort list of parts in increasing order by part identifiers
      ! This is required by the call to icomp subroutine below 
      call sort ( nparts_around, parts_itfc_vefs(2:nparts_around+1, ivef_itfc) )
    end do
    
    call memalloc ( this%num_itfc_vefs, perm_itfc_vefs, __FILE__, __LINE__ )
    do ivef_itfc = 1, this%num_itfc_vefs
      perm_itfc_vefs(ivef_itfc) = ivef_itfc 
    end do
    
    
    ! Re-number vefs in increasing order by the number of parts that share them, 
    ! and among vefs sharing the same list of parts, in increasing order by the list 
    ! of parts shared by the vef 
    call memalloc ( num_rows_parts_itfc_vefs, work1, __FILE__,__LINE__ )
    call memalloc ( num_rows_parts_itfc_vefs, work2, __FILE__,__LINE__ )
    call sort_array_cols_by_row_section( num_rows_parts_itfc_vefs, & 
       &                                 num_rows_parts_itfc_vefs, & 
       &                                 this%num_itfc_vefs, & 
       &                                 parts_itfc_vefs, & 
       &                                 perm_itfc_vefs, &
       &                                 work1, &
       &                                 work2 ) 
    call memfree ( work2, __FILE__,__LINE__ )
    call memfree ( work1, __FILE__,__LINE__ )
    call memfree ( touched_neighbours, __FILE__, __LINE__ )
    
    do ivef_itfc=1,this%num_itfc_vefs
      write(6,'(10i10)') ivef_itfc, this%lst_itfc_vefs(perm_itfc_vefs(ivef_itfc)), parts_itfc_vefs(:, ivef_itfc) 
    end do
  end subroutine par_triangulation_compute_parts_itfc_vefs
  
  subroutine par_triangulation_compute_vefs_and_parts_object(this)
    implicit none
    class(par_triangulation_t), intent(inout) :: this
    integer(ip)                               :: nparts_around
    integer(ip)                               :: ivef_itfc, init_vef, end_vef
    integer(ip)                               :: iobj, ipart
    integer(ip)                               :: num_rows_parts_itfc_vefs
    integer(ip), allocatable                  :: parts_itfc_vefs (:,:)
    integer(ip), allocatable                  :: perm_itfc_vefs(:)
    type(list_iterator_t)                     :: vefs_object_iterator, parts_object_iterator

    assert ( this%p_env%am_i_l1_task() )
    
    call this%compute_parts_itfc_vefs(parts_itfc_vefs,perm_itfc_vefs)
    num_rows_parts_itfc_vefs = size(parts_itfc_vefs,1)
    
    ! Count number_objects
    ivef_itfc = 1
    this%number_objects = 0
    do while ( ivef_itfc <= this%num_itfc_vefs ) 
      if ( ivef_itfc < this%num_itfc_vefs ) then
        do while (icomp(num_rows_parts_itfc_vefs,parts_itfc_vefs(:,ivef_itfc),parts_itfc_vefs(:,ivef_itfc+1)) == 0)
          ivef_itfc = ivef_itfc + 1
          if ( ivef_itfc == this%num_itfc_vefs  ) exit
        end do
      end if  
      this%number_objects = this%number_objects + 1
      ivef_itfc = ivef_itfc + 1
    end do
        
    ! Count number_vefs_per_object and number_parts_per_object
    call this%vefs_object%create(n=this%number_objects)
    call this%parts_object%create(n=this%number_objects)
    ivef_itfc = 1
    this%number_objects = 0
    do while ( ivef_itfc <= this%num_itfc_vefs ) 
      init_vef = ivef_itfc
      if ( ivef_itfc < this%num_itfc_vefs ) then
        do while (icomp(num_rows_parts_itfc_vefs,parts_itfc_vefs(:,ivef_itfc),parts_itfc_vefs(:,ivef_itfc+1)) == 0)
          ivef_itfc = ivef_itfc + 1
          if ( ivef_itfc == this%num_itfc_vefs  ) exit
        end do
      end if  
      end_vef = ivef_itfc
      nparts_around = parts_itfc_vefs(1,end_vef)
      this%number_objects = this%number_objects + 1
      call this%parts_object%sum_to_pointer_index(this%number_objects, nparts_around)
      call this%vefs_object%sum_to_pointer_index(this%number_objects, end_vef-init_vef+1 )
      ivef_itfc = ivef_itfc + 1
    end do
    
    call this%vefs_object%calculate_header()
    call this%parts_object%calculate_header()
    call this%vefs_object%allocate_list_from_pointer()
    call this%parts_object%allocate_list_from_pointer()
    
    ! List number_vefs_per_object and number_parts_per_object
    ivef_itfc=1
    do iobj=1, this%vefs_object%get_num_pointers()
       vefs_object_iterator = this%vefs_object%get_iterator(iobj)
       parts_object_iterator = this%parts_object%get_iterator(iobj)
       
       nparts_around = parts_itfc_vefs(1,ivef_itfc)
       do ipart=1, nparts_around
         call parts_object_iterator%set_current(parts_itfc_vefs(1+ipart,ivef_itfc))
         call parts_object_iterator%next()
       end do
       
       do while(.not. vefs_object_iterator%is_upper_bound())
        call vefs_object_iterator%set_current(this%lst_itfc_vefs(perm_itfc_vefs(ivef_itfc)))
        call vefs_object_iterator%next()
        ivef_itfc = ivef_itfc + 1
       end do
    end do
    
    call this%vefs_object%print(6)
    call this%parts_object%print(6)
    
    call memfree ( parts_itfc_vefs, __FILE__, __LINE__ )
    call memfree ( perm_itfc_vefs, __FILE__, __LINE__ )
  end subroutine par_triangulation_compute_vefs_and_parts_object
  
  subroutine par_triangulation_compute_objects_neighbours_exchange_data ( this, &
                                                                          num_rcv,&
                                                                          list_rcv, &
                                                                          rcv_ptrs,&
                                                                          unpack_idx, &
                                                                          num_snd, &
                                                                          list_snd,&
                                                                          snd_ptrs,&
                                                                          pack_idx )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
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
    
    assert ( this%p_env%am_i_l1_task() )
    
    if (allocated(list_rcv)) call memfree(list_rcv,__FILE__,__LINE__)
    if (allocated(rcv_ptrs)) call memfree(rcv_ptrs,__FILE__,__LINE__)
    if (allocated(unpack_idx)) call memfree(unpack_idx,__FILE__,__LINE__)
    if (allocated(list_snd)) call memfree(list_snd,__FILE__,__LINE__)
    if (allocated(snd_ptrs)) call memfree(snd_ptrs,__FILE__,__LINE__)
    if (allocated(pack_idx)) call memfree(pack_idx,__FILE__,__LINE__)
    
    my_part_id     = this%p_env%get_l1_rank() + 1
    num_neighbours = this%element_import%get_number_neighbours()  
    
    call position_parts_rcv%init(num_neighbours)
    call position_parts_snd%init(num_neighbours)
    call memalloc ( num_neighbours  , list_rcv, __FILE__, __LINE__ )
    call memalloc ( num_neighbours+1, rcv_ptrs, __FILE__, __LINE__ )
    rcv_ptrs = 0 
    
    call memalloc ( num_neighbours  , list_snd, __FILE__, __LINE__ )
    call memalloc ( num_neighbours+1, snd_ptrs, __FILE__, __LINE__ )
    snd_ptrs = 0
    
    num_rcv = 0
    num_snd = 0
    do iobj=1, this%number_objects
       parts_object_iterator = this%parts_object%get_iterator(iobj)
       part_id = parts_object_iterator%get_current()
       if ( my_part_id == part_id ) then
         ! I am owner of the present object
         call parts_object_iterator%next()
         do while ( .not. parts_object_iterator%is_upper_bound() ) 
            part_id = parts_object_iterator%get_current()
            ! Insert part_id in the list of parts I have to send data
            ! Increment by +1 the amount of data I have to send to part_id
            call position_parts_snd%get(key=part_id, val=current_position_parts_snd, stat=istat)
            if ( istat == new_index ) then
              list_snd ( current_position_parts_snd ) = part_id
            end if
            snd_ptrs(current_position_parts_snd+1) = snd_ptrs(current_position_parts_snd+1)+1
            call parts_object_iterator%next()
         end do
       else
         ! I am non-owner of the present object
         call position_parts_rcv%get(key=part_id, val=current_position_parts_rcv, stat=istat)
         if ( istat == new_index ) then
           list_rcv ( current_position_parts_rcv ) = part_id
         end if
         rcv_ptrs(current_position_parts_rcv+1) = rcv_ptrs(current_position_parts_rcv+1)+1 
       end if
    end do
   
    num_rcv = position_parts_rcv%last()
    num_snd = position_parts_snd%last() 
    rcv_ptrs(1) = 1 
    do i=1, num_rcv
      rcv_ptrs(i+1) = rcv_ptrs(i+1) + rcv_ptrs(i)
    end do
    
    snd_ptrs(1) = 1 
    do i=1, num_snd
      snd_ptrs(i+1) = snd_ptrs(i+1) + snd_ptrs(i)
    end do
    
    call memrealloc ( num_snd+1, snd_ptrs, __FILE__, __LINE__ )
    call memrealloc ( num_rcv+1, rcv_ptrs, __FILE__, __LINE__ )
    call memrealloc ( num_snd, list_snd, __FILE__, __LINE__ )
    call memrealloc ( num_rcv, list_rcv, __FILE__, __LINE__ )
    call memalloc ( snd_ptrs(num_snd+1)-1, pack_idx, __FILE__, __LINE__ )
    call memalloc ( rcv_ptrs(num_rcv+1)-1, unpack_idx, __FILE__, __LINE__ )
    
    do iobj=1, this%number_objects
       parts_object_iterator = this%parts_object%get_iterator(iobj)
       part_id = parts_object_iterator%get_current()
       if ( my_part_id == part_id ) then
         ! I am owner of the present object
         call parts_object_iterator%next()
         do while ( .not. parts_object_iterator%is_upper_bound() ) 
           part_id = parts_object_iterator%get_current()
           call position_parts_snd%get(key=part_id, val=current_position_parts_snd, stat=istat)
           pack_idx (snd_ptrs(current_position_parts_snd)) = iobj
           snd_ptrs(current_position_parts_snd) = snd_ptrs(current_position_parts_snd)+1
           call parts_object_iterator%next()
         end do
       else
         ! I am non-owner of the present object
         call position_parts_rcv%get(key=part_id, val=current_position_parts_rcv, stat=istat)
         unpack_idx (rcv_ptrs(current_position_parts_rcv)) = iobj
         rcv_ptrs(current_position_parts_rcv) = rcv_ptrs(current_position_parts_rcv)+1 
       end if
    end do
    
    do i=num_snd, 2, -1
      snd_ptrs(i) = snd_ptrs(i-1) 
    end do
    snd_ptrs(1) = 1 
    
    do i=num_rcv, 2, -1
      rcv_ptrs(i) = rcv_ptrs(i-1) 
    end do
    rcv_ptrs(1) = 1
    
    call position_parts_rcv%free()
    call position_parts_snd%free()
  end subroutine par_triangulation_compute_objects_neighbours_exchange_data 
  
  subroutine par_triangulation_compute_number_global_objects_and_their_gids ( this )
    implicit none
    class(par_triangulation_t), intent(inout) :: this

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

    assert ( this%p_env%am_i_l1_task() )
    my_rank      = this%p_env%get_l1_rank() 
    my_part_id   = my_rank + 1 
    number_parts = this%p_env%get_l1_size()
    
    ! 1. Count/list how many local objects I am responsible to assign a global ID
    call memalloc ( this%number_objects, local_objects_with_gid, __FILE__, __LINE__ )
    number_local_objects_with_gid = 0
    do iobj=1, this%number_objects
      parts_object_iterator = this%parts_object%get_iterator(iobj)
      if ( my_part_id == parts_object_iterator%get_current() ) then
        number_local_objects_with_gid = number_local_objects_with_gid + 1
        local_objects_with_gid (number_local_objects_with_gid) = iobj
      end if
    end do
    
    ! 2. Gather + Scatter
    if ( my_rank == root_pid ) then
      call memalloc( number_parts+1, per_rank_objects_with_gid, __FILE__,__LINE__ )
      call this%p_env%l1_gather (root=root_pid, &
                                 input_data=number_local_objects_with_gid, &
                                 output_data=per_rank_objects_with_gid(2:) ) 
       ! Transform length to header
       per_rank_objects_with_gid(1)=1 
       do i=1, number_parts
          per_rank_objects_with_gid(i+1) = per_rank_objects_with_gid(i) + per_rank_objects_with_gid(i+1) 
       end do
       this%number_global_objects = per_rank_objects_with_gid(number_parts+1)-1 
    else
      call this%p_env%l1_gather (root=root_pid, &
                                 input_data=number_local_objects_with_gid, &
                                 output_data=dummy_integer_array ) 
    end if
    
    call this%p_env%l1_bcast (root=root_pid, data = this%number_global_objects )
    
    if ( my_rank == root_pid ) then
      call this%p_env%l1_scatter (root=root_pid, &
                                  input_data=per_rank_objects_with_gid, &
                                  output_data=start_object_gid) 
      call memfree( per_rank_objects_with_gid, __FILE__,__LINE__ )
    else
      call this%p_env%l1_scatter (root=root_pid, &
                                  input_data=dummy_integer_array, &
                                  output_data=start_object_gid) 
    end if
    
    
    call memalloc (this%number_objects, this%objects_gids)
    do i=1, number_local_objects_with_gid
      this%objects_gids ( local_objects_with_gid(i) ) = start_object_gid
      start_object_gid = start_object_gid + 1 
    end do
    
    ! Set-up objects nearest neighbour exchange data
    ! num_rcv, rcv_ptrs, lst_rcv, unpack_idx
    ! num_snd, snd_ptrs, lst_snd, pack_idx    
    call this%compute_objects_neighbours_exchange_data ( num_rcv, &
                                                         list_rcv,&
                                                         rcv_ptrs,&
                                                         unpack_idx,&
                                                         num_snd,&
                                                         list_snd,&
                                                         snd_ptrs,&
                                                         pack_idx )
    
    call this%p_env%l1_neighbours_exchange ( num_rcv, &
                                             list_rcv,&
                                             rcv_ptrs,&
                                             unpack_idx,&
                                             num_snd,&
                                             list_snd,&
                                             snd_ptrs,&
                                             pack_idx,&
                                             this%objects_gids )
    
    call memfree ( list_rcv, __FILE__, __LINE__ )
    call memfree ( rcv_ptrs, __FILE__, __LINE__ )
    call memfree ( unpack_idx, __FILE__, __LINE__ )
    call memfree ( list_snd, __FILE__, __LINE__ )
    call memfree ( snd_ptrs, __FILE__, __LINE__ )
    call memfree ( pack_idx, __FILE__, __LINE__ )
    call memfree ( local_objects_with_gid, __FILE__, __LINE__ )
  end subroutine par_triangulation_compute_number_global_objects_and_their_gids

  !=============================================================================
  ! Auxiliary subroutines
  subroutine initialize_par_vef_topology (vef)
    implicit none
    type(par_vef_topology_t), intent(inout) :: vef
    
    vef%interface   = -1
    vef%globalID = -1 
  end subroutine initialize_par_vef_topology
  
  subroutine free_par_vef_topology (vef)
    implicit none
    type(par_vef_topology_t), intent(inout) :: vef
    call initialize_par_vef_topology(vef)
  end subroutine free_par_vef_topology

  subroutine free_par_elem_topology(element)
    implicit none
    type(par_elem_topology_t), intent(inout) :: element
    
    if (allocated(element%vefs_GIDs)) then
       call memfree(element%vefs_GIDs, __FILE__, __LINE__)
    end if
    element%interface = -1
    element%mypart = -1 
  end subroutine free_par_elem_topology
  
  subroutine initialize_par_elem_topology(element)
    implicit none
    type(par_elem_topology_t), intent(inout) :: element

    assert(allocated(element%vefs_GIDs))
    element%interface      = -1
    element%mypart      = -1
    element%globalID    = -1
    element%num_vefs = -1
  end subroutine initialize_par_elem_topology

  subroutine par_elem_topology_size (my, n)
    implicit none
    class(par_elem_topology_t), intent(in)  :: my
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    n = size_of_ip*3 + size_of_igp*(my%num_vefs+1)

  end subroutine par_elem_topology_size

  subroutine par_elem_topology_pack (my, n, buffer)
    implicit none
    class(par_elem_topology_t), intent(in)  :: my
    integer(ip)            , intent(in)   :: n
    integer(ieep)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%mypart,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(my%interface,mold)

    start = end + 1
    end   = start + size_of_igp - 1 
    buffer(start:end) = transfer(my%globalID,mold)

    start = end + 1
    end   = start + size_of_ip - 1
    buffer(start:end) = transfer(my%num_vefs,mold)

    start = end + 1
    end   = start + my%num_vefs*size_of_igp - 1
    buffer(start:end) = transfer(my%vefs_GIDs,mold)

  end subroutine par_elem_topology_pack

  subroutine par_elem_topology_unpack(my, n, buffer)
    implicit none
    class(par_elem_topology_t), intent(inout) :: my
    integer(ip)            , intent(in)     :: n
    integer(ieep)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))
    
    start = 1
    end   = start + size_of_ip -1
    my%mypart  = transfer(buffer(start:end), my%mypart)

    start = end + 1
    end   = start + size_of_ip - 1
    my%interface  = transfer(buffer(start:end), my%interface)

    start = end + 1
    end   = start + size_of_igp - 1 
    my%globalID = transfer(buffer(start:end), my%globalID)

    start = end + 1
    end   = start + size_of_ip - 1
    my%num_vefs = transfer(buffer(start:end), my%num_vefs)

    call memalloc( my%num_vefs, my%vefs_GIDs, __FILE__, __LINE__ )

    start = end + 1
    end   = start + my%num_vefs*size_of_igp - 1
    my%vefs_GIDs = transfer(buffer(start:end), my%vefs_GIDs)

  end subroutine par_elem_topology_unpack


end module par_triangulation_names
