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
module JP_triangulation_names
  use types_names
  use memor_names
  use fe_space_types_names
  use element_id_names
  use migratory_element_names
  use hash_table_names  
  implicit none
  private
# include "debug.i90"

  integer(ip), parameter :: JP_triangulation_not_created  = 0 ! Initial state
  integer(ip), parameter :: JP_triangulation_filled       = 1 ! Elems + Vefs arrays allocated and filled 

  type, extends(migratory_element_t) :: JP_element_topology_t
     integer(ip)               :: num_vefs = -1    ! Number of vefs
     integer(ip), allocatable  :: vefs(:)          ! List of Local IDs of the vefs (vertices, edges, faces) that make up this element
     type(reference_element_t), pointer :: geo_reference_element => NULL() ! Topological info of the geometry (SBmod)
     real(rp), allocatable     :: coordinates(:,:)
     integer(ip)               :: order
   contains
     procedure :: size   => element_topology_size
     procedure :: pack   => element_topology_pack
     procedure :: unpack => element_topology_unpack
     procedure :: free   => element_topology_free_unpacked
     procedure :: assign => element_topology_assign
  end type JP_element_topology_t

  type vef_topology_t
     integer(ip)               :: border     = -1       ! Border local id of this vef (only for faces)
     integer(ip)               :: dimension             ! Vef dimension (SBmod)
     integer(ip)               :: num_elems_around = -1 ! Number of elements around vef 
     !integer(ip), allocatable  :: elems_around(:)              ! List of elements around vef 
     class(element_id_t), allocatable  :: elems_around(:)       ! List of elements around vef 
     !class(migratory_element_pointer), allocatable :: elems_around(:) ! We can get any type of element around!?
  end type vef_topology_t

  type, abstract :: JP_triangulation_t
     integer(ip)                            :: state     =  JP_triangulation_not_created  
     integer(ip)                            :: num_vefs  = -1  ! number of vefs (vertices, edges, and faces) 
     integer(ip)                            :: num_elems = -1  ! number of elements
     integer(ip)                            :: num_dims  = -1  ! number of dimensions
     class(migratory_element_set_t)     , allocatable :: element_set
     class(migratory_element_iterator_t), allocatable :: element_iterator
     type (vef_topology_t)              , allocatable :: vefs(:) ! array of vefs in the mesh.
     type (position_hash_table_t)           :: pos_elem_info     ! Topological info hash table (SBmod)
     type (reference_element_t)             :: reference_elements(max_elinf) ! List of topological info's
     integer(ip)                            :: num_boundary_faces ! Number of faces in the boundary 
     integer(ip), allocatable               :: lst_boundary_faces(:) ! List of faces LIDs in the boundary
   contains
     procedure(create_element_iterator_interface)      , deferred :: create_element_iterator
     procedure(free_element_iterator_interface)        , deferred :: free_element_iterator
     !procedure(get_element_topology_pointer_interface) , deferred :: get_element_topology_pointer
  end type JP_triangulation_t
  
  abstract interface
     subroutine create_element_iterator_interface(this,iterator)
       import :: JP_triangulation_t, migratory_element_iterator_t
       implicit none
       class(JP_triangulation_t), target, intent(in)  :: this
       class(migratory_element_iterator_t)    , intent(out) :: iterator
     end subroutine create_element_iterator_interface
     subroutine free_element_iterator_interface(this,iterator)
       import :: JP_triangulation_t, migratory_element_iterator_t
       implicit none
       class(JP_triangulation_t), intent(in)    :: this
       class(migratory_element_iterator_t)    , intent(inout) :: iterator
     end subroutine free_element_iterator_interface
     ! function get_element_topology_pointer_interface(this,id) result(p)
     !   import :: JP_triangulation_t, element_id_t, JP_element_topology_t
     !   implicit none
     !   class(JP_triangulation_t), target, intent(in) :: this
     !   class(element_id_t)              , intent(in) :: id
     !   class(JP_element_topology_t)        , pointer    :: p
     ! end function get_element_topology_pointer_interface
  end interface
   
  ! Types
  public :: JP_triangulation_t , JP_element_topology_t !, element_iterator_t
  
  ! Main Subroutines 
  public :: JP_triangulation_to_dual,  JP_triangulation_create, JP_triangulation_free, JP_triangulation_print

  ! Auxiliary Subroutines (should only be used by modules that have control over type(JP_triangulation_t))
  public :: downcast_to_element_topology, put_topology_element_JP_triangulation, local_id_from_vertices

  ! Constants (should only be used by modules that have control over type(JP_triangulation_t))
  public :: JP_triangulation_not_created, JP_triangulation_filled
  !public :: JP_triangulation_free_elems_data, JP_triangulation_free_objs_data

contains

  !=============================================================================
  subroutine JP_triangulation_create(trian,size,id_mold)
    implicit none
    integer(ip)              , intent(in)    :: size
    class(JP_triangulation_t), intent(inout) :: trian
    class(element_id_t)       :: id_mold

    integer(ip) :: istat,ielem
    class(JP_element_topology_t), pointer :: elem

    trian%num_vefs = -1
    trian%num_elems = -1
    trian%num_dims = -1 

    ! Create iterator using the virtual function
    call trian%create_element_iterator(trian%element_iterator)

    ! Initialize all of the element structs (using the iterator)
    do while(trian%element_iterator%has_next())
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       call elem%create_id(id_mold)
       call initialize_elem_topology( elem )
    end do

    ! Initialize all of the element structs (hard coded, old stuff)
    !do ielem = 1, trian%elem_array_len
    !   call initialize_elem_topology(trian%elems(ielem))
    !end do

    ! Initialization of element fixed info parameters (SBmod)
    call trian%pos_elem_info%init(ht_length)

  end subroutine JP_triangulation_create

  !=============================================================================
  subroutine JP_triangulation_free(trian)
    implicit none
    class(JP_triangulation_t), intent(inout) :: trian
    integer(ip) :: iobj
    class(JP_element_topology_t), pointer :: elem

    assert(trian%state == JP_triangulation_filled) 

    do while(trian%element_iterator%has_next())
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       call free_elem_topology( elem )
    end do
    ! Free iterator using the virtual function
    call trian%free_element_iterator(trian%element_iterator)

    call JP_triangulation_free_objs_data(trian)
    call memfree ( trian%lst_boundary_faces, __FILE__, __LINE__ )

    ! Deallocate fixed info
    do iobj = 1,trian%pos_elem_info%last()
       call reference_element_free (trian%reference_elements(iobj))
    end do
    call trian%pos_elem_info%free

    trian%num_vefs = -1
    trian%num_elems = -1
    trian%num_dims = -1 

    trian%state = JP_triangulation_not_created

  end subroutine JP_triangulation_free
  !=============================================================================

  subroutine JP_triangulation_free_objs_data(trian)
    implicit none
    class(JP_triangulation_t), intent(inout) :: trian
    integer(ip) :: istat,ielem, iobj
    
    if ( trian%state == JP_triangulation_filled ) then
       do iobj=1, trian%num_vefs 
          call free_vef_topology(trian%vefs(iobj)) 
       end do
       ! Deallocate the vef structure array 
       deallocate(trian%vefs, stat=istat)
       check(istat==0)
    end if
  end subroutine JP_triangulation_free_objs_data

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
       !call memfree(vef%elems_around, __FILE__, __LINE__)
       deallocate(vef%elems_around)
    end if
    vef%num_elems_around = -1
  end subroutine free_vef_topology

  subroutine free_elem_topology(element)
    implicit none
    type(JP_element_topology_t), intent(inout) :: element

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
    type(JP_element_topology_t), intent(inout) :: element

    assert(.not. allocated(element%vefs))
    element%num_vefs = -1
  end subroutine initialize_elem_topology

  !=============================================================================
  subroutine JP_triangulation_to_dual(trian, length_trian)  
    implicit none
    ! Parameters
    class(JP_triangulation_t), intent(inout) :: trian
    integer(ip), optional, intent(in)      :: length_trian

    ! Locals
    integer(ip)              :: ielem, iobj, jobj, istat, idime, touch,  length_trian_
    type(hash_table_ip_ip_t) :: visited
    integer(ip), allocatable :: elems_around_pos(:)
    class(JP_element_topology_t), pointer :: elem
    !type(element_id_t)                 :: elem_id

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
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       !elem_id = trian%element_iterator%next()
       !elem => trian%get_element_topology_pointer( elem_id )
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
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       !elem_id = trian%element_iterator%next()
       !elem => trian%get_element_topology_pointer( elem_id )
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

    !call JP_triangulation_print( 6, trian, length_trian_ )

    ! List elements and add vef dimension
    do while(trian%element_iterator%has_next())
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       !elem_id = trian%element_iterator%next()
       !elem => trian%get_element_topology_pointer( elem_id )
       do idime =1, trian%num_dims
          do iobj = elem%geo_reference_element%nvef_dim(idime), &
                    elem%geo_reference_element%nvef_dim(idime+1)-1 
             jobj = elem%vefs(iobj)
             if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
                trian%vefs(jobj)%dimension = idime-1
                if (elems_around_pos(jobj) == 1) then
                   allocate( trian%vefs(jobj)%elems_around (trian%vefs(jobj)%num_elems_around), mold=elem%get_id() )
                   !call memalloc( trian%vefs(jobj)%num_elems_around, trian%vefs(jobj)%elems_around, __FILE__, __LINE__ )
                end if
                trian%vefs(jobj)%elems_around(elems_around_pos(jobj)) = elem%get_id() ! ielem
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

    trian%state = JP_triangulation_filled

  end subroutine JP_triangulation_to_dual

  subroutine put_topology_element_JP_triangulation( elem, trian )
    implicit none
    type(JP_element_topology_t)    , intent(inout)         :: elem
    class(JP_triangulation_t), intent(inout), target :: trian
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
  end subroutine put_topology_element_JP_triangulation

  subroutine local_id_from_vertices( e, nd, list, no, lid ) ! (SBmod)
    implicit none
    type(JP_element_topology_t), intent(in) :: e
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

  subroutine JP_triangulation_print ( lunou,  trian ) ! (SBmod)
    implicit none
    ! Parameters
    integer(ip)                , intent(in)    :: lunou
    class(JP_triangulation_t)  , intent(inout) :: trian
    ! Locals
    class(JP_element_topology_t), pointer      :: elem
    class(element_id_t)         , pointer      :: elem_id
    integer(ip)                                :: iobje, ielem

    assert(trian%state == JP_triangulation_filled) 

    write (lunou,*) '****PRINT TOPOLOGY****'
    write (lunou,*) 'state:'         , trian%state
    write (lunou,*) 'num_vefs:'      , trian%num_vefs
    write (lunou,*) 'num_elems:'     , trian%num_elems
    write (lunou,*) 'num_dims:'      , trian%num_dims

    do while(trian%element_iterator%has_next())
       elem => downcast_to_element_topology( trian%element_iterator%next() )
       !elem_id = trian%element_iterator%next()
       !elem => trian%get_element_topology_pointer( elem_id )
       write (lunou,*) '****PRINT ELEMENT INFO****'
       elem_id => elem%get_id()
       call elem_id%print(lunou)
       write (lunou,*) 'num_vefs:'   , elem%num_vefs
       write (lunou,*) 'vefs:'       , elem%vefs
       write (lunou,*) 'coordinates:', elem%coordinates
       write (lunou,*) 'order:'      , elem%order
       !call reference_element_write ( elem%geo_reference_element )
       write (lunou,*) '****END PRINT ELEMENT INFO****'
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
       !write (lunou,*) 'elems_around', trian%vefs(iobje)%elems_around
       do ielem = 1,trian%vefs(iobje)%num_elems_around
          call  trian%vefs(iobje)%elems_around(ielem)%print(lunou)
       end do
       write (lunou,*) '****END PRINT VEF ',iobje,' INFO****'
    end do

    write (lunou,*) '****END PRINT TOPOLOGY****'
  end subroutine JP_triangulation_print


 !=============================================================================
 !=============================================================================
 !=============================================================================

  subroutine element_topology_size (this, n)
    implicit none
    class(JP_element_topology_t), intent(in)  :: this
    integer(ip)            , intent(out) :: n
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip, size_of_rp

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_rp   = size(transfer(1.0_rp ,mold))
    n = 2*size_of_ip + size(this%coordinates)*size_of_rp

  end subroutine element_topology_size

  subroutine element_topology_pack (this, n, buffer)
    implicit none
    class(JP_element_topology_t), intent(in)  :: this
    integer(ip)              , intent(in)   :: n
    integer(ieep)            , intent(out)  :: buffer(n)
    
    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip,size_of_rp
    integer(ip)   :: start, end

    size_of_ip   = size(transfer(1_ip ,mold))
    size_of_rp   = size(transfer(1.0_rp ,mold))

    start = 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(size(this%coordinates,1),mold)

    start = end + 1
    end   = start + size_of_ip -1
    buffer(start:end) = transfer(this%num_vefs,mold)

    !start = end + 1
    !end   = start + size_of_rp*size(this%coordinates) - 1
    !buffer(start:end) = transfer(this%coordinates,mold)

  end subroutine element_topology_pack

  subroutine element_topology_unpack(this, n, buffer)
    implicit none
    class(JP_element_topology_t), intent(inout)  :: this
    integer(ip)              , intent(in)     :: n
    integer(ieep)            , intent(in)     :: buffer(n)

    ! Locals
    integer(ieep) :: mold(1)
    integer(ip)   :: size_of_ip,size_of_rp
    integer(ip)   :: start, end
    integer(ip)   :: num_dims
    
    size_of_ip = size(transfer(1_ip ,mold))
    size_of_rp = size(transfer(1_rp,mold))

    start = 1
    end   = start + size_of_ip -1
    this%num_vefs  = transfer(buffer(start:end), this%num_vefs)

    start = end + 1
    end   = start + size_of_ip - 1
    num_dims  = transfer(buffer(start:end), num_dims)

   ! call memalloc( num_dims, this%num_vefs, this%coordinates, __FILE__, __LINE__ )
   ! start = end + 1
   ! end   = start + size_of_rp*size(this%coordinates) - 1
   ! this%coordinates  = transfer(buffer(start:end), this%coordinates)
    
  end subroutine element_topology_unpack

  function downcast_to_element_topology(parent) result(this)
    implicit none
    class(migratory_element_t), pointer, intent(in) :: parent
    class(JP_element_topology_t) , pointer             :: this
    select type(parent)
    class is(JP_element_topology_t)
       this => parent
    class default
       write(*,*) 'Cannot downcast to element_topology'
       check(.false.)
    end select
  end function downcast_to_element_topology

  subroutine element_topology_free_unpacked(this)
    implicit none
    class(JP_element_topology_t), intent(inout) :: this
    call memfree( this%vefs, __FILE__, __LINE__ )
    call memfree( this%coordinates, __FILE__, __LINE__ )
  end subroutine element_topology_free_unpacked

  subroutine element_topology_assign(this, that)
    implicit none
    class(JP_element_topology_t), intent(inout) :: this
    class(migratory_element_t), intent(in)    :: that

    select type(that)
    class is(JP_element_topology_t)
       this = that
    class default
       write(*,*) 'Error calling JP_element_topology_t assignment'
       write(*,*) 'cannot assign object of another class'
       check(.false.)
    end select

  end subroutine element_topology_assign

 !=============================================================================
 !=============================================================================
 !=============================================================================

end module JP_triangulation_names
