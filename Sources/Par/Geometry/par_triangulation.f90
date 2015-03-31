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
  use types
  use memor
  use sort_names
  use migratory_element_names
  use fem_triangulation_names
  use fem_element_import_names
  use hash_table_names

  ! Parallel modules
  use par_context_names
  use par_partition_names

  implicit none
# include "debug.i90"
  private

  ! AFM: To think. Is it required to know at some point within FEMPAR the number of global objects ?
  !                With the current design, it is known the number of local objects
  !                of each subdomain, but not the number of global objects.
  !                If it would be required, one possibility would be to store it at type(par_triangulation).

  !                On the other hand, would it be necessary/convenient/helpful to store 
  !                the number of interface objects, and the list of interface objects local IDs?
  !                At this point, I know that this would be convenient/helpful for par_dof_handler_partition.f90.
  !                (That requires fast location of interface objects). Any other point?

  !                I already did store these data on num_itfc_objs and lst_itfc_objs members of type(par_triangulation)
  !                but at the current state I do not know whether this is the most appropriate place. 

  type par_object_topology
     integer(ip)  :: interface  = -1 ! Interface local id of this object
     integer(igp) :: globalID   = -1 ! Global ID of this object
                                     ! Local ID is the position in the array of objects
  end type par_object_topology

  type, extends(migratory_element) :: par_elem_topology
     integer(ip)  :: interface     = -1              ! The boundary number ieboun (if this element is a interface element)
     integer(ip)  :: mypart     = -1              ! To which part this element is mapped to ?
     integer(igp) :: globalID   = -1              ! Global ID of this element
                                                  ! Local ID is the position in the array of elements
     integer(ip)  :: num_objects = -1             ! Number of objects
     integer(igp), allocatable :: objects_GIDs(:) ! List of the GIDs of the objects that make up this element
     
   contains
     procedure :: size   => par_elem_topology_size
     procedure :: pack   => par_elem_topology_pack
     procedure :: unpack => par_elem_topology_unpack

  end type par_elem_topology

  type :: par_triangulation
     type(fem_triangulation)                 :: f_trian             ! Data common with a centralized (serial) triangulation
     integer(ip)                             :: num_elems   = -1    ! should it match f_trian%num_elems or not? 
     integer(ip)                             :: num_ghosts  = -1    ! number of ghost elements (remote neighbors)
     class(migratory_element), allocatable   :: mig_elems(:)        ! Migratory elements list
     type(par_elem_topology),   pointer      :: elems(:) => NULL()  ! Array of elements in the mesh
     type(par_object_topology), allocatable  :: objects(:)          ! array of objects in the mesh
     integer(ip)                             :: num_itfc_objs = -1  ! Number of objects in the interface among subdomains
     integer(ip), allocatable                :: lst_itfc_objs(:)    ! List of objects local IDs in the interface among subdomains 
     integer(ip)                             :: num_itfc_elems= -1  ! Number of elements in the interface
     integer(ip), allocatable                :: lst_itfc_elems(:)   ! List of elements local IDs in the interface
     type(par_context),   pointer            :: p_context => NULL() ! Parallel context describing MPI tasks among which par_triangulation is distributed
     type(fem_element_import)                :: f_el_import         ! Object describing the layout in distributed-memory of the dual graph
                                                                    ! (It is required for nearest neighbour comms on this graph)
     integer(ip)                             :: max_nparts          ! Maximum number of parts around any vef communication object
     integer(ip)                             :: nobjs               ! Number of local vef communication objects
     integer(ip), allocatable                :: lobjs(:,:)          ! List of local vef communication objects
  end type par_triangulation

  ! Types
  public :: par_triangulation, par_elem_topology

  ! Functions
  public :: par_triangulation_free, par_triangulation_to_dual

  ! Auxiliary Subroutines (should only be used by modules that have control over type(par_triangulation))
  public :: free_par_elem_topology, free_par_object_topology, par_triangulation_free_elems_data, par_triangulation_free_objs_data 

contains

  !=============================================================================
  subroutine par_triangulation_free(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: iobj, istat, ielem

    assert(.not. p_trian%f_trian%state == triangulation_not_created) 
    assert(.not. p_trian%f_trian%state == triangulation_created)
    
    call par_triangulation_free_objs_data (p_trian)
    call par_triangulation_free_elems_data(p_trian)
    
    ! Free local portion of triangulation
    call fem_triangulation_free(p_trian%f_trian)
  end subroutine par_triangulation_free

  !=============================================================================
  subroutine par_triangulation_free_objs_data(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: iobj, istat

    if ( p_trian%f_trian%state == triangulation_elems_objects_filled ) then
       do iobj=1, p_trian%f_trian%num_objects 
          call free_par_object_topology(p_trian%objects(iobj)) 
       end do
  
       p_trian%num_itfc_objs = -1
       call memfree ( p_trian%lst_itfc_objs, __FILE__, __LINE__ )
       
       p_trian%max_nparts = -1 
       p_trian%nobjs      = -1 
       call memfree ( p_trian%lobjs, __FILE__, __LINE__ )

       ! Deallocate the object structure array 
       deallocate(p_trian%objects, stat=istat)
       check(istat==0)
    end if
  end subroutine par_triangulation_free_objs_data

  !=============================================================================
  subroutine par_triangulation_free_elems_data(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: ielem, istat
    
    if ( p_trian%f_trian%state == triangulation_elems_objects_filled .or.  &
         p_trian%f_trian%state == triangulation_elems_filled ) then
       
       do ielem=1, p_trian%f_trian%elem_array_len 
          call free_par_elem_topology(p_trian%elems(ielem)) 
       end do
       
       p_trian%num_itfc_elems = -1
       call memfree ( p_trian%lst_itfc_elems, __FILE__, __LINE__ )
       
       call fem_element_import_free ( p_trian%f_el_import )
       
       ! Deallocate the element structure array */
       deallocate(p_trian%mig_elems, stat=istat)
       check(istat==0)
       nullify(p_trian%elems)
       
       p_trian%num_elems  = -1 
       p_trian%num_ghosts = -1
       nullify(p_trian%p_context)
    end if

  end subroutine par_triangulation_free_elems_data

  !=============================================================================
  subroutine par_triangulation_to_dual(p_trian)  
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat

    assert(p_trian%f_trian%state == triangulation_elems_filled .or. p_trian%f_trian%state == triangulation_elems_objects_filled) 

    call par_triangulation_free_objs_data(p_trian)

    call fem_triangulation_to_dual ( p_trian%f_trian, p_trian%num_ghosts )

    ! Allocate p_trian%objects array
    allocate(p_trian%objects(p_trian%f_trian%num_objects), stat=istat)
    check(istat==0)

    ! Initialize par_object_topology objects
    do iobj=1, p_trian%f_trian%num_objects 
       call initialize_par_object_topology(p_trian%objects(iobj)) 
    end do

    ! Assign global ID for all local objects
    do ielem=1, p_trian%f_trian%num_elems
       do iobj=1, p_trian%elems(ielem)%num_objects
          jobj = p_trian%f_trian%elems(ielem)%objects(iobj)
          if ( jobj /= -1 ) then
             p_trian%objects(jobj)%globalID = p_trian%elems(ielem)%objects_GIDs(iobj)
          end if
       end do
    end do

    ! Count number of interface objects
    p_trian%num_itfc_objs = 0
    ! Traverse local view of ghost elements (all the interface objects are there!!!)
    do ielem=p_trian%num_elems+1, p_trian%num_ghosts+p_trian%num_elems
       do iobj=1, p_trian%f_trian%elems(ielem)%num_objects
          jobj = p_trian%f_trian%elems(ielem)%objects(iobj)
          if ( jobj /= -1 ) then
             if (p_trian%objects(jobj)%interface == -1) then
                p_trian%num_itfc_objs = p_trian%num_itfc_objs + 1
                p_trian%objects(jobj)%interface = p_trian%num_itfc_objs
             end if
          end if
       end do
    end do

    ! List interface objects
    call memalloc ( p_trian%num_itfc_objs, p_trian%lst_itfc_objs, __FILE__, __LINE__ )
    p_trian%num_itfc_objs = 0 
    ! Traverse local view of ghost elements (all the interface objects are there!!!)
    do ielem=p_trian%num_elems+1, p_trian%num_ghosts+p_trian%num_elems
       do iobj=1, p_trian%f_trian%elems(ielem)%num_objects
          jobj = p_trian%f_trian%elems(ielem)%objects(iobj)
          if ( jobj /= -1 ) then
             if (p_trian%objects(jobj)%interface == (p_trian%num_itfc_objs + 1)) then
                p_trian%num_itfc_objs = p_trian%num_itfc_objs + 1
                p_trian%lst_itfc_objs(p_trian%num_itfc_objs) = jobj 
             end if
          end if
       end do
    end do

    ! Compute p_trian%max_nparts, p_trian%nobjs, p_trian%lobjs
    ! Re-order (permute) p_trian%lst_itfc_objs accordingly
    call par_triangulation_create_lobjs(p_trian)

  end subroutine par_triangulation_to_dual

  !=============================================================================
  subroutine par_triangulation_create_lobjs(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip)               :: i, j, k, iobj, jelem, est_max_nparts
    integer(ip)               :: ipart, istat, count 
    integer(igp), allocatable :: lst_parts_per_itfc_obj (:,:)
    integer(igp), allocatable :: ws_lobjs_temp (:,:)
    integer(igp), allocatable :: sort_parts_per_itfc_obj_l1 (:)
    integer(igp), allocatable :: sort_parts_per_itfc_obj_l2 (:)
    type(hash_table_ip_ip)    :: ws_parts_visited
    integer(ip), parameter    :: tbl_length = 100

   
    ipart = p_trian%p_context%iam + 1

    ! Compute an estimation (upper bound) of the maximum number of parts around any local interface vef.
    ! This estimation assumes that all elements around all local interface vefs are associated to different parts.
    est_max_nparts = 0
    do i=1, p_trian%num_itfc_objs
       iobj = p_trian%lst_itfc_objs(i)
       est_max_nparts = max(p_trian%f_trian%objects(iobj)%num_elems_around, est_max_nparts)
    end do

    call memalloc ( est_max_nparts+2, p_trian%num_itfc_objs, lst_parts_per_itfc_obj, __FILE__,__LINE__ )

    p_trian%max_nparts = 0
    lst_parts_per_itfc_obj = 0
    do i=1, p_trian%num_itfc_objs
       call ws_parts_visited%init(tbl_length)
       call ws_parts_visited%put(key=ipart,val=1,stat=istat)

       lst_parts_per_itfc_obj (2,i) = ipart
       count = 1

       iobj = p_trian%lst_itfc_objs(i)

       ! Count/list parts around iobj 
       do j=1, p_trian%f_trian%objects(iobj)%num_elems_around 
          jelem=p_trian%f_trian%objects(iobj)%elems_around(j)
          call ws_parts_visited%put(key=p_trian%elems(jelem)%mypart,val=1,stat=istat)
          if ( istat == now_stored ) then
               count = count + 1
               lst_parts_per_itfc_obj (count+1,i) = p_trian%elems(jelem)%mypart 
          end if
       end do

       ! Finish a new column by setting up first and last entries
       lst_parts_per_itfc_obj(1,i) = count
       lst_parts_per_itfc_obj(est_max_nparts+2,i) = p_trian%objects(iobj)%globalID 

       ! Sort list of parts in increasing order by part identifiers
       ! This is required by the call to icomp subroutine below 
       call sort ( count, lst_parts_per_itfc_obj( 2:(count+1), i) )

       p_trian%max_nparts = max(p_trian%max_nparts, count)
       call ws_parts_visited%free
    end do


    ! Re-number vefs in increasing order by the number of parts that share them, 
    ! and among vefs sharing the same list of parts, in increasing order by the list 
    ! of parts shared by the vef 
    call memalloc ( est_max_nparts+2, sort_parts_per_itfc_obj_l1, __FILE__,__LINE__ )
    call memalloc ( est_max_nparts+2, sort_parts_per_itfc_obj_l2, __FILE__,__LINE__ )
    call sort_array_cols_by_row_section( est_max_nparts+2,           & ! Rows of lst_parts_per_itfc_obj 
       &                                 est_max_nparts+2,           & ! LD of lst_parts_per_itfc_obj
       &                                 p_trian%num_itfc_objs,      & ! Cols of lst_parts_per_itfc_obj
       &                                 lst_parts_per_itfc_obj,     & 
       &                                 p_trian%lst_itfc_objs,      &
       &                                 sort_parts_per_itfc_obj_l1, &
       &                                 sort_parts_per_itfc_obj_l2 )
    call memfree ( sort_parts_per_itfc_obj_l2, __FILE__,__LINE__ )
    call memfree ( sort_parts_per_itfc_obj_l1, __FILE__,__LINE__ )

    ! Re-compute p_trian%objects(:)%interface to reflect the new status of p_trian%lst_itfc_objs
    do i=1, p_trian%num_itfc_objs
       iobj = p_trian%lst_itfc_objs(i)
       assert ( p_trian%objects(iobj)%interface /= -1 )
       p_trian%objects(iobj)%interface = i
    end do

    ! Identify communication objects 
    call memalloc ( p_trian%max_nparts+4, p_trian%num_itfc_objs, ws_lobjs_temp, __FILE__,__LINE__ )
    
    p_trian%nobjs=1

    ! Prepare first object
    ws_lobjs_temp(1, p_trian%nobjs)= -1  ! Unused entry (maintained for historical reasons) 
    ws_lobjs_temp(2, p_trian%nobjs) = 1  ! Begin obj

    k  = 1 
    do i=1,p_trian%num_itfc_objs-1
       if(icomp(p_trian%max_nparts+1,lst_parts_per_itfc_obj(:,i),lst_parts_per_itfc_obj(:,i+1)) /= 0) then
          ! Complete current object
          ws_lobjs_temp(3,p_trian%nobjs)=k ! end obj
          ws_lobjs_temp(4,p_trian%nobjs)=lst_parts_per_itfc_obj(1,i)
          ws_lobjs_temp(5:4+lst_parts_per_itfc_obj(1,i),p_trian%nobjs) = &
               & lst_parts_per_itfc_obj(2:1+lst_parts_per_itfc_obj(1,i),i)
          ws_lobjs_temp(4+lst_parts_per_itfc_obj(1,i)+1:,p_trian%nobjs) = 0 

          ! Prepare next object
          p_trian%nobjs=p_trian%nobjs+1
          ws_lobjs_temp(1,p_trian%nobjs)=-1     ! Unused entry (maintained for historical reasons) 
          ws_lobjs_temp(2,p_trian%nobjs)=k+1    ! begin obj
       end if
       k=k+1
    end do

    ! Complete last object
    ws_lobjs_temp(3,p_trian%nobjs)=k        ! end obj
    ws_lobjs_temp(4,p_trian%nobjs)=lst_parts_per_itfc_obj(1,i)
    ws_lobjs_temp(5:4+lst_parts_per_itfc_obj(1,i),p_trian%nobjs) = &
         & lst_parts_per_itfc_obj(2:1+lst_parts_per_itfc_obj(1,i),i)
    ws_lobjs_temp(4+lst_parts_per_itfc_obj(1,i)+1:,p_trian%nobjs) = 0 

    ! Reallocate lobjs and add internal object first
    call memalloc (p_trian%max_nparts+4, p_trian%nobjs, p_trian%lobjs, __FILE__, __LINE__)

    ! Copy ws_lobjs_temp to lobjs ...
    do i=1,p_trian%nobjs
       p_trian%lobjs(:,i)=ws_lobjs_temp(:,i)
    end do

    call memfree ( ws_lobjs_temp, __FILE__,__LINE__ )
    call memfree ( lst_parts_per_itfc_obj, __FILE__,__LINE__ )


    write(*,'(a,i10)') 'Number of local objects:', &
       &  p_trian%nobjs

    write(*,'(a)') 'List of local objects:'
    do i=1,p_trian%nobjs 
       write(*,'(10i10)') i, p_trian%lobjs(:,i)
    end do

  end subroutine par_triangulation_create_lobjs
   

  !=============================================================================
  ! Auxiliary subroutines
  subroutine initialize_par_object_topology (object)
    implicit none
    type(par_object_topology), intent(inout) :: object
    
    object%interface   = -1
    object%globalID = -1 
  end subroutine initialize_par_object_topology
  
  subroutine free_par_object_topology (object)
    implicit none
    type(par_object_topology), intent(inout) :: object
    call initialize_par_object_topology(object)
  end subroutine free_par_object_topology

  subroutine free_par_elem_topology(element)
    implicit none
    type(par_elem_topology), intent(inout) :: element
    
    if (allocated(element%objects_GIDs)) then
       call memfree(element%objects_GIDs, __FILE__, __LINE__)
    end if
    element%interface = -1
    element%mypart = -1 
  end subroutine free_par_elem_topology
  
  subroutine initialize_par_elem_topology(element)
    implicit none
    type(par_elem_topology), intent(inout) :: element

    assert(allocated(element%objects_GIDs))
    element%interface      = -1
    element%mypart      = -1
    element%globalID    = -1
    element%num_objects = -1
  end subroutine initialize_par_elem_topology

  subroutine par_elem_topology_size (my, n)
    implicit none
    class(par_elem_topology), intent(in)  :: my
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    n = size_of_ip*3 + size_of_igp*(my%num_objects+1)

  end subroutine par_elem_topology_size

  subroutine par_elem_topology_pack (my, n, buffer)
    implicit none
    class(par_elem_topology), intent(in)  :: my
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
    buffer(start:end) = transfer(my%num_objects,mold)

    start = end + 1
    end   = start + my%num_objects*size_of_igp - 1
    buffer(start:end) = transfer(my%objects_GIDs,mold)

  end subroutine par_elem_topology_pack

  subroutine par_elem_topology_unpack(my, n, buffer)
    implicit none
    class(par_elem_topology), intent(inout) :: my
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
    my%num_objects = transfer(buffer(start:end), my%num_objects)

    call memalloc( my%num_objects, my%objects_GIDs, __FILE__, __LINE__ )

    start = end + 1
    end   = start + my%num_objects*size_of_igp - 1
    my%objects_GIDs = transfer(buffer(start:end), my%objects_GIDs)

  end subroutine par_elem_topology_unpack


end module par_triangulation_names
