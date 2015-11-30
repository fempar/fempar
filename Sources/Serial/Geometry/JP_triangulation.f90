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
  use JP_element_topology_names
  use migratory_element_names
  use hash_table_names  
  implicit none
  private
# include "debug.i90"

  integer(ip), parameter :: JP_triangulation_not_created  = 0 ! Initial state
  integer(ip), parameter :: JP_triangulation_filled       = 1 ! Elems + Vefs arrays allocated and filled 

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
  public :: JP_triangulation_t , JP_element_topology_t
  
  ! Main Subroutines 
  public :: JP_triangulation_to_dual,  JP_triangulation_create, JP_triangulation_free, JP_triangulation_print
  public :: create_reference_elements

  ! Constants
  public :: JP_triangulation_not_created, JP_triangulation_filled

contains

  !=============================================================================
  subroutine JP_triangulation_create(trian,size,id_mold)
    implicit none
    integer(ip)              , intent(in)    :: size
    class(JP_triangulation_t), intent(inout) :: trian
    class(element_id_t)      , intent(in)    :: id_mold

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
       call elem%create ! (id_mold)
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
       call elem%free()
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
    class(element_id_t)         , allocatable  :: elem_id
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
       !elem => trian%element_iterator%next_element_topology()
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
       call trian%element_iterator%current_id(elem_id)
       do idime =1, trian%num_dims
          do iobj = elem%geo_reference_element%nvef_dim(idime), &
                    elem%geo_reference_element%nvef_dim(idime+1)-1 
             jobj = elem%vefs(iobj)
             if (jobj /= -1) then ! jobj == -1 if vef belongs to neighbouring processor
                trian%vefs(jobj)%dimension = idime-1
                if (elems_around_pos(jobj) == 1) then
                   allocate( trian%vefs(jobj)%elems_around (trian%vefs(jobj)%num_elems_around), mold=elem_id )
                   !call memalloc( trian%vefs(jobj)%num_elems_around, trian%vefs(jobj)%elems_around, __FILE__, __LINE__ )
                end if
                trian%vefs(jobj)%elems_around(elems_around_pos(jobj)) = elem_id ! elem%get_id() ! ielem
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

  ! This function replaces put_topology_element_triangulation
  ! and it performs the element loop internally (not to be called
  ! inside an element loop as the previous one.
  subroutine create_reference_elements(trian)
    implicit none
    class(JP_triangulation_t), intent(inout), target :: trian
    ! Locals
    class(JP_element_topology_t), pointer :: elem
    integer(ip) :: nvef, v_key, ndime, etype, pos_elinf, istat
    integer(ip) :: aux_val

    do while(trian%element_iterator%has_next())
       elem => downcast_to_element_topology( trian%element_iterator%next() )

       ! Legacy code
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
    end do

  end subroutine create_reference_elements

  subroutine JP_triangulation_print ( lunou,  trian ) ! (SBmod)
    implicit none
    ! Parameters
    integer(ip)                , intent(in)    :: lunou
    class(JP_triangulation_t)  , intent(inout) :: trian
    ! Locals
    class(JP_element_topology_t), pointer      :: elem
    class(element_id_t)         , allocatable  :: elem_id
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
       !elem_id => elem%get_id() 
       call trian%element_iterator%current_id(elem_id)
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

end module JP_triangulation_names
