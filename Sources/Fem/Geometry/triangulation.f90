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
module fem_triangulation_class
  use types
  use memor
  use fem_space_types
  use hash_table_class
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: triangulation_not_created          = 0 ! Initial state
  integer(ip), parameter :: triangulation_created              = 1 ! Elems array allocated 
  integer(ip), parameter :: triangulation_elems_filled         = 2 ! Elems array filled with data 
  integer(ip), parameter :: triangulation_elems_objects_filled = 3 ! Objects array allocated and filled 


  ! State graph for type(fem_triangulation) instances
  ! case (state) 
  !  triangulation_not_created:
  !     case (action)
  !        triangulation_create:
  !            create_triangulation() 
  !            state = triangulation_created
  !        else:
  !           ERROR
  !  triangulation_created:
  !     case (action):
  !        mesh_to_triangulation:
  !           mesh_to_triangulation()
  !           state = triangulation_elems_filled
  !        triangulation_free:
  !           free_triangulation()
  !           state = triangulation_not_created
  !        else:
  !          ERROR
  ! triangulation_elems_filled:
  !     case (action):
  !        triangulation_to_dual:
  !           triangulation_to_dual()
  !           state = triangulation_elems_objects_filled
  !        mesh_to_triangulation:
  !          free_element_topology for all elements
  !          free_element_array
  !          mesh_to_triangulation()
  !        triangulation_free:
  !           free_triangulation()
  !           state = triangulation_not_created
  !        else:
  !          ERROR
  ! triangulation_elems_objects_filled:
  !     case(action):
  !        triangulation_to_dual:
  !          free_object_topology for all objects
  !          free_object_array
  !          triangulation_to_dual()
  !        mesh_to_triangulation: 
  !          free_object_topology for all objects
  !          free_object_array
  !          free_element_topology for all elements
  !          free_element_array
  !          mesh_to_triangulation()
  !          state = triangulation_elems_filled
  !        triangulation_free:
  !           free_triangulation()
  !           state = triangulation_not_created

  type elem_topology
     integer(ip)               :: num_objects = -1    ! Number of objects
     integer(ip) , allocatable :: objects(:)          ! List of Local IDs of the objects (vertices, edges, faces) that make up this element
     type(fem_fixed_info), pointer :: topology => NULL() ! Topological info of the geometry (SBmod)
  end type elem_topology

  type object_topology
     integer(ip)  :: border     = -1 ! Border local id of this object
     integer(ip)               :: dimension             ! Object dimension (SBmod)
     integer(ip)               :: num_elems_around = -1 ! Number of elements around object 
     integer(ip), allocatable  :: elems_around(:)       ! List of elements around object 
  end type object_topology

  type fem_triangulation
     integer(ip) :: state          =  triangulation_not_created  
     integer(ip) :: num_objects    = -1  ! number of objects 
     integer(ip) :: num_elems      = -1  ! number of elements
     integer(ip) :: num_dims       = -1  ! number of dimensions
     integer(ip) :: elem_array_len = -1  ! length that the elements array is allocated for. 
     type(elem_topology), allocatable :: elems(:) ! array of elements in the mesh.
     type(object_topology) , allocatable :: objects(:)    ! array of objects in the mesh.
     type (hash_table_ip_ip)             :: ht_elem_info  ! Topological info hash table (SBmod)
     integer(ip)                         :: cur_elinf = 0 !  "
     type (fem_fixed_info), allocatable  :: lelem_info(:) ! List of topological info's
     !integer(ip), allocatable            :: lst_boundary_objs(:)    ! List of objects local IDs in the boundary 
  end type fem_triangulation

  ! Types
  public :: fem_triangulation

  ! Main Subroutines 
  public :: fem_triangulation_create, fem_triangulation_free, fem_triangulation_to_dual

  ! Auxiliary Subroutines (should only be used by modules that have control over type(fem_triangulation))
  public :: free_elem_topology, free_object_topology, put_topology_element_triangulation, local_id_from_vertices

  ! Constants (should only be used by modules that have control over type(fem_triangulation))
  public :: triangulation_not_created, triangulation_created, triangulation_elems_filled, triangulation_elems_objects_filled

contains

  !=============================================================================
  subroutine fem_triangulation_create(len,trian)
    implicit none
    integer(ip)            , intent(in)    :: len
    type(fem_triangulation), intent(inout) :: trian
    integer(ip) :: istat,ielem

    assert(trian%state == triangulation_not_created)

    trian%elem_array_len = len
    trian%num_objects = -1
    trian%num_elems = -1
    trian%num_dims = -1 

    ! Allocate the element structure array 
    allocate(trian%elems(trian%elem_array_len), stat=istat)
    check(istat==0)

    ! Initialize all of the element structs
    do ielem = 1, trian%elem_array_len
       call initialize_elem_topology(trian%elems(ielem))
    end do

    trian%state = triangulation_created

    ! Initialization of element fixed info parameters (SBmod)
    call trian%ht_elem_info%init(ht_length)
    trian%cur_elinf = 1
    allocate(trian%lelem_info(max_elinf), stat=istat)
    if(istat/=0) then 
       write(6,'(a)') '[Fempar Fatal Error] ***Memory (de)allocation failed.'
       write(6,'(a,i20)') 'Error code: ', istat
       write(6,'(a)') 'Caller routine: triangulation::fem_triangulation_create::lelem_info'
    end if
  end subroutine fem_triangulation_create

  !=============================================================================
  subroutine fem_triangulation_free(trian)
    implicit none
    type(fem_triangulation), intent(inout) :: trian
    integer(ip) :: istat,ielem, iobj

    assert(.not. trian%state == triangulation_not_created) 

    if ( trian%state == triangulation_elems_objects_filled ) then
       do iobj=1, trian%num_objects 
          call free_object_topology(trian%objects(iobj)) 
       end do
       ! Deallocate the object structure array 
       deallocate(trian%objects, stat=istat)
       check(istat==0)
    end if

    if ( trian%state == triangulation_elems_objects_filled .or. & 
         trian%state == triangulation_elems_filled ) then
       do ielem=1, trian%elem_array_len 
          call free_elem_topology(trian%elems(ielem)) 
       end do
    end if

    ! Deallocate the element structure array */
    deallocate(trian%elems, stat=istat)
    check(istat==0)

    ! Deallocate fixed info
    do iobj = 1,trian%cur_elinf-1
       call fem_element_fixed_info_free (trian%lelem_info(iobj))
    end do
    deallocate (trian%lelem_info)
    trian%cur_elinf = 0
    call trian%ht_elem_info%free

    trian%elem_array_len = -1 
    trian%num_objects = -1
    trian%num_elems = -1
    trian%num_dims = -1 
  end subroutine fem_triangulation_free

  ! Auxiliary subroutines
  subroutine initialize_object_topology (object)
    implicit none
    type(object_topology), intent(inout) :: object

    assert(.not.allocated(object%elems_around))
    object%num_elems_around = 0 
  end subroutine initialize_object_topology

  subroutine free_object_topology (object)
    implicit none
    type(object_topology), intent(inout) :: object

    if (allocated(object%elems_around)) then
       call memfree(object%elems_around, __FILE__, __LINE__)
    end if
    object%num_elems_around = -1
  end subroutine free_object_topology

  subroutine free_elem_topology(element)
    implicit none
    type(elem_topology), intent(inout) :: element

    if (allocated(element%objects)) then
       call memfree(element%objects, __FILE__, __LINE__)
    end if
    element%num_objects = -1
    nullify( element%topology )
  end subroutine free_elem_topology

  subroutine initialize_elem_topology(element)
    implicit none
    type(elem_topology), intent(inout) :: element

    assert(.not. allocated(element%objects))
    element%num_objects = -1
  end subroutine initialize_elem_topology

  subroutine fem_triangulation_to_dual(trian, num_ghosts)  
    implicit none
    ! Parameters
    type(fem_triangulation), intent(inout) :: trian
    integer(ip), optional, intent(in)      :: num_ghosts

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat, idime, num_elems
    type(hash_table_ip_ip)   :: visited
    integer(ip), allocatable :: elems_around_pos(:)

    assert(trian%state == triangulation_elems_filled .or. trian%state == triangulation_elems_objects_filled) 

    if ( trian%state == triangulation_elems_objects_filled ) then
       do iobj=1, trian%num_objects 
          call free_object_topology(trian%objects(iobj)) 
       end do
       ! Deallocate the object structure array 
       deallocate(trian%objects, stat=istat)
       check(istat==0)
    end if

    num_elems = trian%num_elems
    if ( present(num_ghosts) ) then
       num_elems = num_elems + num_ghosts
    end if

    ! Count objects
    call visited%init(max(5,int(real(num_elems,rp)*0.2_rp,ip))) 
    trian%num_objects = 0
    do ielem=1, num_elems
       do iobj=1, trian%elems(ielem)%num_objects
          jobj = trian%elems(ielem)%objects(iobj)
          if (jobj /= -1) then ! jobj == -1 if object belongs to neighbouring processor
             call visited%put(key=jobj, val=1, stat=istat)
             if (istat == now_stored) trian%num_objects = trian%num_objects + 1
          end if
       end do
    end do
    call visited%free

    ! Allocate the object structure array 
    allocate(trian%objects(trian%num_objects), stat=istat)
    check(istat==0)
    do iobj=1, trian%num_objects
       call initialize_object_topology(trian%objects(iobj))
    end do

    ! Count elements around each object
    do ielem=1, num_elems
       do iobj=1, trian%elems(ielem)%num_objects
          jobj = trian%elems(ielem)%objects(iobj)
          if (jobj /= -1) then ! jobj == -1 if object belongs to neighbouring processor
             trian%objects(jobj)%num_elems_around = trian%objects(jobj)%num_elems_around + 1 
          end if
       end do
    end do

    call memalloc ( trian%num_objects, elems_around_pos, __FILE__, __LINE__ )
    elems_around_pos = 1

    ! List elements and add object dimensio
    do ielem=1, num_elems
       do idime =1, trian%num_dims    ! (SBmod)
          do iobj = trian%elems(ielem)%topology%nobje_dim(idime), &
               trian%elems(ielem)%topology%nobje_dim(idime+1)-1 
             !do iobj=1, trian%elems(ielem)%num_objects
             jobj = trian%elems(ielem)%objects(iobj)
             if (jobj /= -1) then ! jobj == -1 if object belongs to neighbouring processor
                trian%objects(jobj)%dimension = idime
                if (elems_around_pos(jobj) == 1) then
                   call memalloc( trian%objects(jobj)%num_elems_around, trian%objects(jobj)%elems_around, __FILE__, __LINE__ )
                end if
                trian%objects(jobj)%elems_around(elems_around_pos(jobj)) = ielem
                elems_around_pos(jobj) = elems_around_pos(jobj) + 1 
             end if
          end do
       end do
    end do

    call memfree ( elems_around_pos, __FILE__, __LINE__ )
    trian%state = triangulation_elems_objects_filled

    ! Put dimension of the object
  end subroutine fem_triangulation_to_dual

  subroutine put_topology_element_triangulation( ielem, trian ) !(SBmod)
    implicit none
    type(fem_triangulation), intent(inout), target :: trian
    integer(ip),             intent(in)            :: ielem
    ! Locals
    integer(ip)              :: nobje, v_key, ndime, etype, pos_elinf, istat
    logical(ip) :: created

    nobje = trian%elems(ielem)%num_objects 
    ndime = trian%num_dims
    ! Variable values depending of the element ndime
    if(ndime == 2) then        ! 2D
       if(nobje == 6) then     ! Linear triangles (P1)
          etype = P_type_id
       elseif(nobje == 8) then ! Linear quads (Q1)
          etype = Q_type_id
       end if
    elseif(ndime == 3) then    ! 3D
       if(nobje == 14) then     ! Linear tetrahedra (P1)
          etype = P_type_id
       elseif(nobje == 26) then ! Linear hexahedra (Q1)
          etype = Q_type_id
       end if
    end if
    ! Assign pointer to topological information
    v_key = ndime + (max_ndime+1)*etype + (max_ndime+1)*(max_FE_types+1)
    call trian%ht_elem_info%put(key=v_key,val=trian%cur_elinf,stat=istat)
    if ( istat == now_stored) then 
       ! Create fixed info if not constructed
       call fem_element_fixed_info_create(trian%lelem_info(trian%cur_elinf),etype,  &
            &                             1,ndime,created)
       assert(created)
       pos_elinf = trian%cur_elinf
       trian%cur_elinf = trian%cur_elinf + 1
    else if ( istat == was_stored ) then
       call trian%ht_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
       assert ( istat == key_found )
    end if
    trian%elems(ielem)%topology => trian%lelem_info(pos_elinf)

  end subroutine put_topology_element_triangulation

  subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
    implicit none
    type(elem_topology), intent(in) :: e
    integer(ip), intent(in)  :: nd, list(:), no
    integer(ip), intent(out) :: lid
    ! Locals
    integer(ip)              :: first, last, io, iv, jv, ivl, c
    lid = -1

    do io = e%topology%nobje_dim(nd), e%topology%nobje_dim(nd+1)-1
       first =  e%topology%crxob%p(io)
       last = e%topology%crxob%p(io+1) -1
       if ( last - first + 1  == no ) then 
          do iv = first,last
             ivl = e%objects(e%topology%crxob%l(iv)) ! LID of vertices of the ef
             c = 0
             do jv = 1,no
                if ( ivl ==  list(jv) ) then
                   c  = 1 ! vertex in the external ef
                   exit
                end if
             end do
             if (c == 0) exit
          end do
          if (c == 1) then ! object in the external element
             lid = e%objects(io)
             exit
          end if
       end if
    end do
  end subroutine local_id_from_vertices

end module fem_triangulation_class
