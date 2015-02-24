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
module par_triangulation_class
  ! Serial modules
  use types
  use memor
  use migratory_element_class
  use fem_triangulation_class

  ! Parallel modules
  use par_context_class
  use par_partition_class

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
     integer(ip)  :: border     = -1 ! The boundary number ioboun (if this object is a border object)
     integer(igp) :: globalID   = -1 ! Global ID of this object
                                     ! Local ID is the position in the array of objects
  end type par_object_topology

  type, extends(migratory_element) :: par_elem_topology
     integer(ip)  :: border     = -1              ! The boundary number ieboun (if this element is a border element)
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
     type(par_partition),   pointer          :: p_part => NULL()    ! How are the elements and objects distributed among processors ?
  end type par_triangulation

  ! Types
  public :: par_triangulation, par_elem_topology

  ! Functions
  public :: par_triangulation_free, par_triangulation_to_dual

  ! Auxiliary Subroutines (should only be used by modules that have control over type(par_triangulation))
  public :: free_par_elem_topology, free_par_object_topology  

contains

  !=============================================================================
  subroutine par_triangulation_free(p_trian)
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip) :: iobj, istat, ielem

    assert(.not. p_trian%f_trian%state == triangulation_not_created) 
    
    if ( p_trian%f_trian%state == triangulation_elems_objects_filled ) then
       do iobj=1, p_trian%f_trian%num_objects 
          call free_par_object_topology(p_trian%objects(iobj)) 
       end do
  
       p_trian%num_itfc_objs = -1
       call memfree ( p_trian%lst_itfc_objs, __FILE__, __LINE__ )

       ! Deallocate the object structure array 
       deallocate(p_trian%objects, stat=istat)
       check(istat==0)
    end if

    if ( p_trian%f_trian%state == triangulation_elems_objects_filled .or. & 
         p_trian%f_trian%state == triangulation_elems_filled ) then
       do ielem=1, p_trian%f_trian%elem_array_len 
          call free_par_elem_topology(p_trian%elems(ielem)) 
       end do
    end if

    ! Deallocate the element structure array */
    deallocate(p_trian%mig_elems, stat=istat)
    check(istat==0)
    nullify(p_trian%elems)

    p_trian%num_elems  = -1 
    p_trian%num_ghosts = -1
    nullify(p_trian%p_part)
    
    ! Free local portion of triangulation
    call fem_triangulation_free(p_trian%f_trian)
    
  end subroutine par_triangulation_free

  subroutine par_triangulation_to_dual(p_trian)  
    implicit none
    ! Parameters
    type(par_triangulation), intent(inout) :: p_trian

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat

    assert(p_trian%f_trian%state == triangulation_elems_filled .or. p_trian%f_trian%state == triangulation_elems_objects_filled) 


    if ( p_trian%f_trian%state == triangulation_elems_objects_filled ) then
       do iobj=1, p_trian%f_trian%num_objects 
          call free_par_object_topology(p_trian%objects(iobj)) 
       end do

       p_trian%num_itfc_objs = -1
       call memfree ( p_trian%lst_itfc_objs, __FILE__, __LINE__ )

       ! Deallocate the object structure array 
       deallocate(p_trian%objects, stat=istat)
       check(istat==0)
    end if

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
             if (p_trian%objects(jobj)%border == -1) then
                p_trian%num_itfc_objs = p_trian%num_itfc_objs + 1
                p_trian%objects(jobj)%border = p_trian%num_itfc_objs
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
             if (p_trian%objects(jobj)%border == (p_trian%num_itfc_objs + 1)) then
                p_trian%num_itfc_objs = p_trian%num_itfc_objs + 1
                p_trian%lst_itfc_objs(p_trian%num_itfc_objs) = jobj 
             end if
          end if
       end do
    end do
  end subroutine par_triangulation_to_dual
   

  !=============================================================================

  ! Auxiliary subroutines
  subroutine initialize_par_object_topology (object)
    implicit none
    type(par_object_topology), intent(inout) :: object
    
    object%border   = -1
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
    element%border = -1
  end subroutine free_par_elem_topology
  
  subroutine initialize_par_elem_topology(element)
    implicit none
    type(par_elem_topology), intent(inout) :: element

    assert(allocated(element%objects_GIDs))
    element%border      = -1
    element%globalID    = -1
    element%num_objects = -1
  end subroutine initialize_par_elem_topology

  subroutine par_elem_topology_size (my, n)
    implicit none
    class(par_elem_topology), intent(in)  :: my
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    n = size_of_ip*2 + size_of_igp*(my%num_objects+1)

  end subroutine par_elem_topology_size

  subroutine par_elem_topology_pack (my, n, buffer)
    implicit none
    class(par_elem_topology), intent(in)  :: my
    integer(ip)            , intent(in)   :: n
    integer(ip)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(my%border,mold)

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
    integer(ip)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ip) :: mold(1)
    integer(ip) :: size_of_ip, size_of_igp
    integer(ip) :: start, end
    
    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_igp  = size(transfer(1_igp,mold))
    
    start = 1
    end   = start + size_of_ip -1
    my%border  = transfer(buffer(start:end), my%border)

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


end module par_triangulation_class
