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
module serial_triangulation_names
  use types_names
  use memor_names
  use fe_space_types_names
  use hash_table_names  
  implicit none
  private
# include "debug.i90"

  integer(ip), parameter :: triangulation_not_created  = 0 ! Initial state
  integer(ip), parameter :: triangulation_filled       = 1 ! Elems + Vefs arrays allocated and filled 

  ! In the future this class will probably be abstract and the function deferred.
  ! I'm starting with this one to implement iterators without recode everything.
  !type, abstract :: element_topology_t
  type :: element_topology_t
     integer(ip)               :: local_id = 0
     integer(ip)               :: num_vefs = -1    ! Number of vefs
     integer(ip), allocatable  :: vefs(:)          ! List of Local IDs of the vefs (vertices, edges, faces) that make up this element
     type(reference_element_t), pointer :: geo_reference_element => NULL() ! Topological info of the geometry (SBmod)
     real(rp), allocatable     :: coordinates(:,:)
     integer(ip)               :: order
     !contains
     !procedure(get_local_id_interface), deferred  :: get_local_id
  end type element_topology_t
  !abstract interface
  !   function get_local_id_interface (this) result(lid)
  !     import :: element_topology_t
  !     implicit none
  !     class(element_topology_t), intent(in)  :: this
  !     integer(ip) :: lid
  !   end function get_local_id_interface
  !end interface

  type, abstract :: element_iterator_t
     private
   contains
     procedure(next_interface)    , deferred :: next
     procedure(has_next_interface), deferred :: has_next
  end type element_iterator_t
  abstract interface
     function next_interface (this) result(p)
       import :: element_iterator_t, element_topology_t
       implicit none
       class(element_iterator_t), intent(inout)  :: this
       class(element_topology_t), pointer        :: p
     end function next_interface
     function has_next_interface(this) result(res)
       import :: element_iterator_t
       implicit none
       class(element_iterator_t), intent(inout) :: this
       logical :: res
     end function has_next_interface	
  end interface

  type vef_topology_t
     integer(ip)               :: border     = -1       ! Border local id of this vef (only for faces)
     integer(ip)               :: dimension             ! Vef dimension (SBmod)
     integer(ip)               :: num_elems_around = -1 ! Number of elements around vef 
     integer(ip), allocatable  :: elems_around(:)       ! List of elements around vef 
  end type vef_topology_t

  type, abstract :: serial_triangulation_t
     integer(ip)                            :: state     =  triangulation_not_created  
     integer(ip)                            :: num_vefs  = -1  ! number of vefs (vertices, edges, and faces) 
     integer(ip)                            :: num_elems = -1  ! number of elements
     integer(ip)                            :: num_dims  = -1  ! number of dimensions
     class(element_iterator_t), allocatable :: element_iterator
     type (vef_topology_t)    , allocatable :: vefs(:) ! array of vefs in the mesh.
     type (position_hash_table_t)           :: pos_elem_info  ! Topological info hash table (SBmod)
     type (reference_element_t)             :: reference_elements(max_elinf) ! List of topological info's
     integer(ip)                            :: num_boundary_faces ! Number of faces in the boundary 
     integer(ip), allocatable               :: lst_boundary_faces(:) ! List of faces LIDs in the boundary
   contains
     procedure(create_element_iterator_interface), deferred :: create_element_iterator
     procedure(free_element_iterator_interface)  , deferred :: free_element_iterator
     procedure :: create => serial_triangulation_create
     procedure :: free   => serial_triangulation_free
  end type serial_triangulation_t
  
  abstract interface
      subroutine create_element_iterator_interface(this,iterator)
       import :: serial_triangulation_t, element_iterator_t
       implicit none
       class(serial_triangulation_t), target, intent(in)  :: this
       class(element_iterator_t)    , intent(out) :: iterator
     end subroutine create_element_iterator_interface
      subroutine free_element_iterator_interface(this,iterator)
       import :: serial_triangulation_t, element_iterator_t
       implicit none
       class(serial_triangulation_t), intent(in)    :: this
       class(element_iterator_t)    , intent(inout) :: iterator
     end subroutine free_element_iterator_interface
  end interface
   
  ! Types

  public :: serial_triangulation_t, element_topology_t, element_iterator_t
  
  ! Main Subroutines 
  public :: triangulation_to_dual, triangulation_print

  ! Auxiliary Subroutines (should only be used by modules that have control over type(triangulation_t))
  public :: put_topology_element_triangulation, local_id_from_vertices

  ! Constants (should only be used by modules that have control over type(triangulation_t))
  public :: triangulation_not_created, triangulation_filled
  !public :: triangulation_free_elems_data, triangulation_free_objs_data

contains

  !=============================================================================
  subroutine serial_triangulation_create(trian,len)
    implicit none
    integer(ip)                  , intent(in)    :: len
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip) :: istat,ielem

    trian%num_vefs = -1
    trian%num_elems = -1
    trian%num_dims = -1 

    ! Create iterator using the virtual function
    call trian%create_element_iterator(trian%element_iterator)

    ! Initialize all of the element structs (using the iterator)
    do while(trian%element_iterator%has_next())
       call initialize_elem_topology(trian%element_iterator%next())
    end do

    ! Initialize all of the element structs (hard coded, old stuff)
    !do ielem = 1, trian%elem_array_len
    !   call initialize_elem_topology(trian%elems(ielem))
    !end do

    ! Initialization of element fixed info parameters (SBmod)
    call trian%pos_elem_info%init(ht_length)

  end subroutine serial_triangulation_create

  !=============================================================================
  subroutine serial_triangulation_free(trian)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip) :: iobj

    assert(trian%state == triangulation_filled) 

    do while(trian%element_iterator%has_next())
       call free_elem_topology(trian%element_iterator%next()) 
    end do
    ! Free iterator using the virtual function
    call trian%free_element_iterator(trian%element_iterator)

    call triangulation_free_objs_data(trian)
    call memfree ( trian%lst_boundary_faces, __FILE__, __LINE__ )

    ! Deallocate fixed info
    do iobj = 1,trian%pos_elem_info%last()
       call reference_element_free (trian%reference_elements(iobj))
    end do
    call trian%pos_elem_info%free

    trian%num_vefs = -1
    trian%num_elems = -1
    trian%num_dims = -1 

    trian%state = triangulation_not_created

  end subroutine serial_triangulation_free
  !=============================================================================

  subroutine triangulation_free_objs_data(trian)
    implicit none
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip) :: istat,ielem, iobj
    
    if ( trian%state == triangulation_filled ) then
       do iobj=1, trian%num_vefs 
          call free_vef_topology(trian%vefs(iobj)) 
       end do
       ! Deallocate the vef structure array 
       deallocate(trian%vefs, stat=istat)
       check(istat==0)
    end if
  end subroutine triangulation_free_objs_data

  !=============================================================================
  ! Auxiliary subroutines
  subroutine initialize_vef_topology (vef)
    implicit none
    type(vef_topology_t), intent(inout) :: vef

    assert(.not.allocated(vef%elems_around))
    vef%num_elems_around = 0 
  end subroutine initialize_vef_topology

  subroutine free_vef_topology (vef)
    implicit none
    type(vef_topology_t), intent(inout) :: vef

    if (allocated(vef%elems_around)) then
       call memfree(vef%elems_around, __FILE__, __LINE__)
    end if
    vef%num_elems_around = -1
  end subroutine free_vef_topology

  subroutine free_elem_topology(element)
    implicit none
    type(element_topology_t), intent(inout) :: element

    if (allocated(element%vefs)) then
       call memfree(element%vefs, __FILE__, __LINE__)
    end if

    if (allocated(element%coordinates)) then
       call memfree(element%coordinates, __FILE__, __LINE__)
    end if

    element%num_vefs = -1
    nullify( element%geo_reference_element )
  end subroutine free_elem_topology

  subroutine initialize_elem_topology(element)
    implicit none
    type(element_topology_t), intent(inout) :: element

    assert(.not. allocated(element%vefs))
    element%num_vefs = -1
  end subroutine initialize_elem_topology

  !=============================================================================
  subroutine triangulation_to_dual(trian, length_trian)  
    implicit none
    ! Parameters
    class(serial_triangulation_t), intent(inout) :: trian
    integer(ip), optional, intent(in)      :: length_trian

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat, idime, touch,  length_trian_
    type(hash_table_ip_ip_t) :: visited
    integer(ip), allocatable :: elems_around_pos(:)
    class(element_topology_t), pointer :: elem

    if (present(length_trian)) then
       length_trian_ = length_trian
    else
       length_trian_ = trian%num_elems 
    endif

    ! Count vefs
    call visited%init(max(5,int(real(length_trian_,rp)*0.2_rp,ip))) 
    trian%num_vefs = 0
    touch = 1
    do while(trian%element_iterator%has_next())
       elem => trian%element_iterator%next()
       do iobj=1, elem%num_vefs
          jobj = elem%vefs(iobj)
          if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
             !call visited%put(key=jobj, val=1, stat=istat)
             call visited%put(key=jobj, val=touch, stat=istat)
             if (istat == now_stored) trian%num_vefs = trian%num_vefs + 1
          end if
       end do
    end do
    ! do ielem=1, length_trian_
    !    do iobj=1, trian%elems(ielem)%num_vefs
    !       jobj = trian%elems(ielem)%vefs(iobj)
    !       if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
    !          !call visited%put(key=jobj, val=1, stat=istat)
    !          call visited%put(key=jobj, val=touch, stat=istat)
    !          if (istat == now_stored) trian%num_vefs = trian%num_vefs + 1
    !       end if
    !    end do
    ! end do
    call visited%free

    ! Allocate the vef structure array 
    allocate(trian%vefs(trian%num_vefs), stat=istat)
    check(istat==0)
    do iobj=1, trian%num_vefs
       call initialize_vef_topology(trian%vefs(iobj))
    end do

    ! Count elements around each vef
    do while(trian%element_iterator%has_next())
       elem => trian%element_iterator%next()
       do iobj=1, elem%num_vefs
          jobj = elem%vefs(iobj)
          if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
             trian%vefs(jobj)%num_elems_around = trian%vefs(jobj)%num_elems_around + 1 
          end if
       end do
    end do
    ! do ielem=1, length_trian_
    !    do iobj=1, trian%elems(ielem)%num_vefs
    !       jobj = trian%elems(ielem)%vefs(iobj)
    !       if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
    !          trian%vefs(jobj)%num_elems_around = trian%vefs(jobj)%num_elems_around + 1 
    !       end if
    !    end do
    ! end do
    call memalloc ( trian%num_vefs, elems_around_pos, __FILE__, __LINE__ )
    elems_around_pos = 1

    !call triangulation_print( 6, trian, length_trian_ )

    ! List elements and add vef dimension
    do while(trian%element_iterator%has_next())
       elem => trian%element_iterator%next()
       do idime =1, trian%num_dims
          do iobj = elem%geo_reference_element%nvef_dim(idime), &
                    elem%geo_reference_element%nvef_dim(idime+1)-1 
             jobj = elem%vefs(iobj)
             if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
                trian%vefs(jobj)%dimension = idime-1
                if (elems_around_pos(jobj) == 1) then
                   call memalloc( trian%vefs(jobj)%num_elems_around, trian%vefs(jobj)%elems_around, __FILE__, __LINE__ )
                end if
                trian%vefs(jobj)%elems_around(elems_around_pos(jobj)) = elem%local_id ! ielem
                elems_around_pos(jobj) = elems_around_pos(jobj) + 1 
             end if
          end do
       end do
    end do
    ! do ielem=1, length_trian_
    !    do idime =1, trian%num_dims    ! (SBmod)
    !       do iobj = trian%elems(ielem)%geo_reference_element%nvef_dim(idime), &
    !            trian%elems(ielem)%geo_reference_element%nvef_dim(idime+1)-1 
    !          !do iobj=1, trian%elems(ielem)%num_vefs
    !          jobj = trian%elems(ielem)%vefs(iobj)
    !          if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
    !             trian%vefs(jobj)%dimension = idime-1
    !             if (elems_around_pos(jobj) == 1) then
    !                call memalloc( trian%vefs(jobj)%num_elems_around, trian%vefs(jobj)%elems_around, __FILE__, __LINE__ )
    !             end if
    !             trian%vefs(jobj)%elems_around(elems_around_pos(jobj)) = ielem
    !             elems_around_pos(jobj) = elems_around_pos(jobj) + 1 
    !          end if
    !       end do
    !    end do
    ! end do

    ! Assign border and count boundary faces
    trian%num_boundary_faces = 0
    do iobj = 1, trian%num_vefs
       if ( trian%vefs(iobj)%dimension == trian%num_dims -1 ) then
          if ( trian%vefs(iobj)%num_elems_around == 1 ) then 
             trian%num_boundary_faces = trian%num_boundary_faces + 1
             trian%vefs(iobj)%border = trian%num_boundary_faces
          end if
       end if
    end do

    ! List boundary faces
    call memalloc (  trian%num_boundary_faces, trian%lst_boundary_faces,  __FILE__, __LINE__ )
    do iobj = 1, trian%num_vefs
       if ( trian%vefs(iobj)%dimension == trian%num_dims -1 ) then
          if ( trian%vefs(iobj)%num_elems_around == 1 ) then 
             trian%lst_boundary_faces(trian%vefs(iobj)%border) = iobj
          end if
       end if
    end do

    call memfree ( elems_around_pos, __FILE__, __LINE__ )

    trian%state = triangulation_filled

  end subroutine triangulation_to_dual

  subroutine put_topology_element_triangulation( elem, trian )
    implicit none
    type(element_topology_t)    , intent(inout)         :: elem
    class(serial_triangulation_t), intent(inout), target :: trian
    ! Locals
    integer(ip) :: nvef, v_key, ndime, etype, pos_elinf, istat
    logical :: created
    integer(ip) :: aux_val
    nvef  = elem%num_vefs 
    ndime = trian%num_dims
    ! Variable values depending of the element ndime
    etype = 0
    if(ndime == 2) then       ! 2D
       if(nvef == 6) then     ! Linear triangles (P1)
          etype = P_type_id
       elseif(nvef == 8) then ! Linear quads (Q1)
          etype = Q_type_id
       end if
    elseif(ndime == 3) then    ! 3D
       if(nvef == 14) then     ! Linear tetrahedra (P1)
          etype = P_type_id
       elseif(nvef == 26) then ! Linear hexahedra (Q1)
          etype = Q_type_id
       end if
    end if
    assert( etype /= 0 )

    ! Assign pointer to topological information
    v_key = ndime + (max_ndime+1)*etype + (max_ndime+1)*(max_FE_types+1)
    call trian%pos_elem_info%get(key=v_key,val=pos_elinf,stat=istat)
    if ( istat == new_index) then
       ! Create fixed info if not constructed
       call reference_element_create(trian%reference_elements(pos_elinf),etype,1,ndime)
    end if
    elem%geo_reference_element => trian%reference_elements(pos_elinf)
  end subroutine put_topology_element_triangulation

  subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
    implicit none
    type(element_topology_t), intent(in) :: e
    integer(ip), intent(in)  :: nd, list(:), no
    integer(ip), intent(out) :: lid
    ! Locals
    integer(ip)              :: first, last, io, iv, jv, ivl, c
    lid = -1

    do io = e%geo_reference_element%nvef_dim(nd), e%geo_reference_element%nvef_dim(nd+1)-1
       first =  e%geo_reference_element%crxob%p(io)
       last = e%geo_reference_element%crxob%p(io+1) -1
       if ( last - first + 1  == no ) then 
          do iv = first,last
             ivl = e%vefs(e%geo_reference_element%crxob%l(iv)) ! LID of vertices of the ef
             c = 0
             do jv = 1,no
                if ( ivl ==  list(jv) ) then
                   c  = 1 ! vertex in the external ef
                   exit
                end if
             end do
             if (c == 0) exit
          end do
          if (c == 1) then ! vef in the external element
             lid = e%vefs(io)
             exit
          end if
       end if
    end do
  end subroutine local_id_from_vertices

  subroutine triangulation_print ( lunou,  trian ) ! (SBmod)
    implicit none
    ! Parameters
    integer(ip)                  , intent(in)    :: lunou
    class(serial_triangulation_t), intent(inout) :: trian
    ! Locals
    class(element_topology_t)    , pointer       :: elem
    integer(ip) :: iobje

    assert(trian%state == triangulation_filled) 

    write (lunou,*) '****PRINT TOPOLOGY****'
    write (lunou,*) 'state:'         , trian%state
    write (lunou,*) 'num_vefs:'      , trian%num_vefs
    write (lunou,*) 'num_elems:'     , trian%num_elems
    write (lunou,*) 'num_dims:'      , trian%num_dims

    do while(trian%element_iterator%has_next())
       elem => trian%element_iterator%next()
       write (lunou,*) '****PRINT ELEMENT ',elem%local_id,' INFO****'
       write (lunou,*) 'num_vefs:'   , elem%num_vefs
       write (lunou,*) 'vefs:'       , elem%vefs
       write (lunou,*) 'coordinates:', elem%coordinates
       write (lunou,*) 'order:'      , elem%order
       !call reference_element_write ( elem%geo_reference_element )
       write (lunou,*) '****END PRINT ELEMENT ',elem%local_id,' INFO****'
    end do
    ! do ielem = 1, length_trian_
    !    write (lunou,*) '****PRINT ELEMENT ',ielem,' INFO****'
    !    write (lunou,*) 'num_vefs:', trian%elems(ielem)%num_vefs
    !    write (lunou,*) 'vefs:', trian%elems(ielem)%vefs
    !    write (lunou,*) 'coordinates:', trian%elems(ielem)%coordinates
    !    write (lunou,*) 'order:', trian%elems(ielem)%order
    !    !call reference_element_write ( trian%elems(ielem)%geo_reference_element )
    !    write (lunou,*) '****END PRINT ELEMENT ',ielem,' INFO****'
    ! end do

    do iobje = 1, trian%num_vefs
       write (lunou,*) '****PRINT VEF ',iobje,' INFO****'
       write (lunou,*) 'border', trian%vefs(iobje)%border
       write (lunou,*) 'dimension', trian%vefs(iobje)%dimension
       write (lunou,*) 'num_elems_around', trian%vefs(iobje)%num_elems_around
       write (lunou,*) 'elems_around', trian%vefs(iobje)%elems_around
       write (lunou,*) '****END PRINT VEF ',iobje,' INFO****'
    end do

    write (lunou,*) '****END PRINT TOPOLOGY****'
  end subroutine triangulation_print

end module serial_triangulation_names
