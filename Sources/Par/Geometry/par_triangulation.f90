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
  use mesh_to_triangulation_names
  use base_static_triangulation_names 
 
  ! Parallel modules
  use par_context_names
  use par_environment_names
  use psi_penv_mod_names
  use par_mesh_names
  use par_element_exchange_names
  use par_conditions_names
 

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
     integer(ip)  :: itfc       = -1 ! Interface local id of this vef
     integer(igp) :: globalID   = -1 ! Global ID of this vef
                                     ! Local ID is the position in the array of vefs
  end type par_vef_topology_t

  type, extends(migratory_element_t) :: par_elem_topology_t
     integer(ip)  :: itfc       = -1           ! The boundary number ieboun (if this element is a itfc element)
     integer(ip)  :: mypart     = -1           ! To which part this element is mapped to ?
     integer(igp) :: globalID   = -1           ! Global ID of this element
                                               ! Local ID is the position in the array of elements
     integer(ip)  :: num_vefs   = -1           ! Number of vefs
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
     integer(ip)                             :: number_global_objects = -1
     integer(ip)                             :: number_objects        = -1
     integer(igp), allocatable               :: objects_gids(:)
     integer(ip), allocatable                :: objects_dimension(:)
     
     type(list_t)                            :: vefs_object
     type(list_t)                            :: parts_object
     
     type(coarse_triangulation_t), pointer   :: coarse_triangulation
  contains
     procedure, non_overridable          :: get_coarse_triangulation                       => par_triangulation_get_coarse_triangulation
  
     procedure, non_overridable, private :: compute_parts_itfc_vefs                        => par_triangulation_compute_parts_itfc_vefs
     procedure, non_overridable, private :: compute_vefs_and_parts_object                  => par_triangulation_compute_vefs_and_parts_object
     procedure, non_overridable, private :: compute_objects_dimension                      => par_triangulation_compute_objects_dimension
     procedure, non_overridable, private :: compute_objects_neighbours_exchange_data       => par_triangulation_compute_objects_neighbours_exchange_data
     procedure, non_overridable, private :: compute_number_global_objects_and_their_gids   => par_triangulation_compute_number_global_objects_and_their_gids
     procedure, non_overridable, private :: setup_coarse_triangulation                     => par_triangulation_setup_coarse_triangulation
     procedure, non_overridable, private :: gather_coarse_cell_gids                        => par_triangulation_gather_coarse_cell_gids
     procedure, non_overridable, private :: gather_coarse_vefs_rcv_counts_and_displs       => par_triangulation_gather_coarse_vefs_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_vefs_dimension                   => par_triangulation_gather_coarse_vefs_dimension
     procedure, non_overridable, private :: gather_coarse_vefs_gids                        => par_triangulation_gather_coarse_vefs_gids
     procedure, non_overridable, private :: fetch_l2_part_id_neighbours                    => par_triangulation_fetch_l2_part_id_neighbours
     procedure, non_overridable, private :: gather_coarse_dgraph_rcv_counts_and_displs     => par_triangulation_gather_coarse_dgraph_rcv_counts_and_displs
     procedure, non_overridable, private :: gather_coarse_dgraph_lextn_and_lextp           => par_triangulation_gather_coarse_dgraph_lextn_and_lextp
     procedure, non_overridable, private :: adapt_coarse_raw_arrays                        => par_triangulation_adapt_coarse_raw_arrays
  end type par_triangulation_t
  
  ! Types
  public :: par_triangulation_t, par_elem_topology_t

  ! Functions
  public :: par_mesh_to_triangulation, par_triangulation_free, par_triangulation_print, par_triangulation_to_dual

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
    integer(ip) :: iobj, istat, icell

    if(p_trian%p_env%am_i_lgt1_task()) then
     call p_trian%coarse_triangulation%free()
     deallocate (p_trian%coarse_triangulation, stat=istat)
     check (istat==0)
    else 
     assert(p_trian%state == par_triangulation_filled) 
     call par_triangulation_free_objs_data (p_trian)
     call par_triangulation_free_elems_data(p_trian)
     call triangulation_free(p_trian%triangulation)
    end if
    nullify(p_trian%coarse_triangulation)
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
        
       p_trian%number_objects = -1
       call p_trian%vefs_object%free()
       call p_trian%parts_object%free()
       
       p_trian%number_global_objects = -1
       call memfree ( p_trian%objects_gids, __FILE__, __LINE__ )
       call memfree ( p_trian%objects_dimension, __FILE__, __LINE__ )
       
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
       
       call p_trian%element_import%free()
       
       ! Deallocate the element structure array */
       deallocate(p_trian%mig_elems, stat=istat)
       check(istat==0)
       nullify(p_trian%elems)
       
       p_trian%num_elems  = -1 
       p_trian%num_ghosts = -1
       nullify(p_trian%p_env)
    end if

  end subroutine par_triangulation_free_elems_data
  
    subroutine par_mesh_to_triangulation (p_gmesh, p_trian, p_cond)
    implicit none
    ! Parameters
    type(par_mesh_t)         , target  , intent(in)    :: p_gmesh ! Geometry mesh
    type(par_triangulation_t), target  , intent(inout) :: p_trian 
    type(par_conditions_t)   , optional, intent(inout) :: p_cond

    ! Locals
    integer(ip) :: istat, ielem, iobj, jobj, state
    integer(ip) :: num_elems, num_ghosts, num_verts
    type (hash_table_igp_ip_t) :: hash ! Topological info hash table (SBmod)
    integer(ip) :: ilele, nvert, jelem, jlele, idime, count, ivere 
    integer(igp):: iobjg
    integer(igp), allocatable :: aux_igp(:)
    integer(ip), allocatable :: aux(:)
    integer(ip) :: aux_val
    type(list_t), pointer :: vertices_vef
    type(par_context_t), pointer :: l1_context
    
    ! Set a reference to the type(par_environment_t) instance describing the set of MPI tasks
    ! among which this type(par_triangulation) instance is going to be distributed 
    p_trian%p_env => p_gmesh%p_env
    if(p_trian%p_env%am_i_l1_task()) then
       l1_context    => p_trian%p_env%get_l1_context()
       state = p_trian%state
       assert(state == par_triangulation_not_created .or. state == par_triangulation_filled )

       if ( state == par_triangulation_filled ) call par_triangulation_free(p_trian)

       call p_trian%element_import%create  ( p_trian%p_env%get_l1_rank()+1, &
                                             p_trian%p_env%get_l1_size(), &
                                             p_gmesh%f_mesh%nelem, &
                                             p_gmesh%f_mesh_dist%nebou, &
                                             p_gmesh%f_mesh_dist%lebou, &
                                             p_gmesh%f_mesh_dist%pextn, &
                                             p_gmesh%f_mesh_dist%lextn, &
                                             p_gmesh%f_mesh_dist%lextp)

       ! Now we are sure that the local portion of p_trian, i.e., p_trian%f_trian is in triangulation_not_created state
       ! Let's create and fill it
       num_elems  = p_gmesh%f_mesh%nelem
       num_ghosts = p_trian%element_import%get_number_ghost_elements()

       p_trian%num_ghosts = num_ghosts
       p_trian%num_elems  = num_elems

       ! Fill local portion with local data
       call mesh_to_triangulation_fill_elements ( p_gmesh%f_mesh, p_trian%triangulation, num_elems + num_ghosts, p_cond%f_conditions )

       ! p_trian%f_trian%num_elems = num_elems+num_ghosts
       ! **IMPORTANT NOTE**: the code that comes assumes that edges and faces in p_trian%f_trian are numbered after vertices.
       ! This requirement is currently fulfilled by mesh_to_triangulation (in particular, by geom2topo within) but we should
       ! keep this in mind all the way through. If this were not assumed, we should find a way to identify a corner within
       ! p_trian%f_trian, and map from a corner local ID in p_trian%f_trian to a vertex local ID in p_gmesh%f_mesh.
       num_verts = p_gmesh%f_mesh%npoin

       ! Create array of elements with room for ghost elements
       allocate( par_elem_topology_t :: p_trian%mig_elems(num_elems + num_ghosts), stat=istat)
       check(istat==0)

       select type( this => p_trian%mig_elems )
       type is(par_elem_topology_t)
          p_trian%elems => this
       end select

       p_trian%elems(:)%itfc = -1 
       do ielem=1, p_gmesh%f_mesh_dist%nebou
          p_trian%elems(p_gmesh%f_mesh_dist%lebou(ielem))%itfc = ielem
       end do

       p_trian%num_itfc_elems = p_gmesh%f_mesh_dist%nebou
       call memalloc( p_trian%num_itfc_elems, p_trian%lst_itfc_elems, __FILE__, __LINE__ )
       p_trian%lst_itfc_elems = p_gmesh%f_mesh_dist%lebou

       ! Fill array of elements (local ones)
       do ielem=1, num_elems
          p_trian%elems(ielem)%mypart      = l1_context%get_rank() + 1
          p_trian%elems(ielem)%globalID    = p_gmesh%f_mesh_dist%emap%l2g(ielem)
          p_trian%elems(ielem)%num_vefs = p_trian%triangulation%elems(ielem)%num_vefs
          call memalloc( p_trian%elems(ielem)%num_vefs, p_trian%elems(ielem)%vefs_GIDs, __FILE__, __LINE__ )
          do iobj=1, p_trian%elems(ielem)%num_vefs
             jobj = p_trian%triangulation%elems(ielem)%vefs(iobj)
             if ( jobj <= num_verts ) then ! It is a corner => re-use global ID
                p_trian%elems(ielem)%vefs_GIDs(iobj) = p_gmesh%f_mesh_dist%nmap%l2g(jobj)
             else ! It is an edge or face => generate new local-global ID (non-consistent, non-consecutive)
                ! The ISHFT(1,50) is used to start numbering efs after vertices, assuming nvert < 2**60
                p_trian%elems(ielem)%vefs_GIDs(iobj) = ISHFT(int(p_gmesh%f_mesh_dist%ipart,igp),int(32,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(60,igp))
                !p_trian%elems(ielem)%vefs_GIDs(iobj) = ISHFT(int(p_gmesh%f_mesh_dist%ipart,igp),int(6,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(6,igp))
             end if
          end do
       end do

       ! Get vefs_GIDs from ghost elements
       call ghost_elements_exchange ( p_trian%p_env, p_trian%element_import, p_trian%elems )

       ! Allocate elem_topology in triangulation for ghost elements  (SBmod)
       do ielem = num_elems+1, num_elems+num_ghosts       
          p_trian%triangulation%elems(ielem)%num_vefs = p_trian%elems(ielem)%num_vefs
          call memalloc(p_trian%triangulation%elems(ielem)%num_vefs, p_trian%triangulation%elems(ielem)%vefs, &
               & __FILE__, __LINE__)
          p_trian%triangulation%elems(ielem)%vefs = -1
       end do

       ! Put the topology info in the ghost elements
       do ielem= num_elems+1, num_elems+num_ghosts
          call put_topology_element_triangulation ( ielem, p_trian%triangulation )
       end do

       ! Hash table global to local for ghost elements
       call hash%init(num_ghosts)
       do ielem=num_elems+1, num_elems+num_ghosts
          aux_val = ielem
          call hash%put( key = p_trian%elems(ielem)%globalID, val = aux_val, stat=istat)
       end do

       do ielem = 1, p_gmesh%f_mesh_dist%nebou     ! Loop interface elements 
          ! Step 1: Put LID of vertices in the ghost elements (f_mesh_dist)
          ilele = p_gmesh%f_mesh_dist%lebou(ielem) ! local ID element
          ! aux : array of ilele (LID) vertices in GID
          nvert  = p_trian%triangulation%elems(ilele)%reference_fe_geo%get_number_vertices()
          call memalloc( nvert, aux_igp, __FILE__, __LINE__  )
          do iobj = 1, nvert                        ! vertices only
             aux_igp(iobj) = p_trian%elems(ilele)%vefs_GIDs(iobj) ! extract GIDs vertices
          end do
          do jelem = p_gmesh%f_mesh_dist%pextn(ielem), & 
               p_gmesh%f_mesh_dist%pextn(ielem+1)-1  ! external neighbor elements
             call hash%get(key = p_gmesh%f_mesh_dist%lextn(jelem), val=jlele, stat=istat) ! LID external element
             do jobj = 1, p_trian%triangulation%elems(jlele)%reference_fe_geo%get_number_vertices() ! vertices external 
                if ( p_trian%triangulation%elems(jlele)%vefs(jobj) == -1) then
                   do iobj = 1, nvert
                      if ( aux_igp(iobj) == p_trian%elems(jlele)%vefs_GIDs(jobj) ) then
                         ! Put LID of vertices for ghost_elements
                         p_trian%triangulation%elems(jlele)%vefs(jobj) =  p_trian%triangulation%elems(ilele)%vefs(iobj)
                      end if
                   end do
                end if
             end do
          end do
          call memfree(aux_igp, __FILE__, __LINE__) 
       end do

       do ielem = 1, p_gmesh%f_mesh_dist%nebou     ! Loop interface elements 
          ! Step 2: Put LID of efs in the ghost elements (f_mesh_dist) 
          ilele = p_gmesh%f_mesh_dist%lebou(ielem) ! local ID element
          do jelem = p_gmesh%f_mesh_dist%pextn(ielem), &
               p_gmesh%f_mesh_dist%pextn(ielem+1)-1  ! external neighbor elements
             call hash%get(key = p_gmesh%f_mesh_dist%lextn(jelem), val=jlele, stat=istat) ! LID external element
             vertices_vef => p_trian%triangulation%elems(jlele)%reference_fe_geo%get_vertices_vef()
             ! loop over all efs of external elements
             do idime =2,p_trian%triangulation%num_dims
                do iobj = p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime-1), &
                     p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime)-1 
                   if ( p_trian%triangulation%elems(jlele)%vefs(iobj) == -1) then ! efs not assigned yet
                      count = 1
                      ! loop over vertices of every ef
                      do jobj = vertices_vef%p(iobj), vertices_vef%p(iobj+1)-1    
                         ivere = vertices_vef%l(jobj)
                         if (p_trian%triangulation%elems(jlele)%vefs(ivere) == -1) then
                            count = 0 ! not an vef of the local triangulation
                            exit
                         end if
                      end do
                      if (count == 1) then
                         nvert = vertices_vef%p(iobj+1)-vertices_vef%p(iobj)
                         call memalloc( nvert, aux, __FILE__, __LINE__)
                         count = 1
                         do jobj = vertices_vef%p(iobj), vertices_vef%p(iobj+1)-1
                            ivere = vertices_vef%l(jobj)
                            aux(count) = p_trian%triangulation%elems(jlele)%vefs(ivere)
                            count = count+1
                         end do
                         call local_id_from_vertices( p_trian%triangulation%elems(ilele), idime, aux, nvert, &
                              p_trian%triangulation%elems(jlele)%vefs(iobj) )
                         call memfree(aux, __FILE__, __LINE__) 
                      end if
                   end if
                end do
             end do
          end do
       end do

       call hash%free

       call par_triangulation_to_dual ( p_trian )

       !pause

       ! Check results
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do ielem = 1,num_elems+num_ghosts
       !       write (*,*) '****ielem:',ielem          
       !       write (*,*) '****LID_vefs ghost:',p_trian%f_trian%elems(ielem)%vefs
       !       write (*,*) '****GID_vefs ghost:',p_trian%elems(ielem)%vefs_GIDs
       !    end do
       ! end if
       ! pause

       ! write (*,*) '*********************************'
       ! write (*,*) '*********************************' 
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do iobj = 1, p_trian%f_trian%num_vefs 
       !       write (*,*) 'is interface vef',iobj, ' :', p_trian%vefs(iobj)%interface 
       !       write (*,*) 'is interface vef',iobj, ' :', p_trian%f_trian%vefs(iobj)%dimension
       !    end do
       ! end if
       ! write (*,*) '*********************************'
       ! write (*,*) '*********************************'

       ! Step 3: Make GID consistent among processors (p_part%elems%vefs_GIDs)
       do iobj = 1, p_trian%triangulation%num_vefs 
          if ( (p_trian%vefs(iobj)%itfc .ne. -1) .and. &
               (p_trian%triangulation%vefs(iobj)%dime >= 1) ) then
             idime = p_trian%triangulation%vefs(iobj)%dime+1
             iobjg = -1

             do jelem = 1,p_trian%triangulation%vefs(iobj)%num_elems_around
                jlele = p_trian%triangulation%vefs(iobj)%elems_around(jelem)

                do jobj = p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime-1), &
                     & p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime)-1 ! efs of neighbor els
                   if ( p_trian%triangulation%elems(jlele)%vefs(jobj) == iobj ) then
                      if ( iobjg == -1 ) then 
                         iobjg  = p_trian%elems(jlele)%vefs_GIDs(jobj)
                      else
                         iobjg  = min(iobjg,p_trian%elems(jlele)%vefs_GIDs(jobj))
                         exit
                      end if
                   end if
                end do

             end do


             do jelem = 1,p_trian%triangulation%vefs(iobj)%num_elems_around
                jlele = p_trian%triangulation%vefs(iobj)%elems_around(jelem)
                do jobj = p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime-1), &
                     & p_trian%triangulation%elems(jlele)%reference_fe_geo%get_first_vef_id_of_dimension(idime)-1 ! efs of neighbor els
                   if ( p_trian%triangulation%elems(jlele)%vefs(jobj) == iobj) then
                      p_trian%elems(jlele)%vefs_GIDs(jobj) = iobjg
                      exit
                   end if
                end do
             end do
             p_trian%vefs(iobj)%globalID = iobjg          
          end if
       end do

       p_trian%state = par_triangulation_filled
       ! Check results
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do ielem = 1,num_elems+num_ghosts
       !       write (*,*) '****ielem:',ielem          
       !       write (*,*) '****LID_vefs ghost:',p_trian%f_trian%elems(ielem)%vefs
       !       write (*,*) '****GID_vefs ghost:',p_trian%elems(ielem)%vefs_GIDs
       !    end do
       ! end if
       ! pause
    else
       ! AFM: TODO: Partially broadcast p_trian%f_trian from 1st level tasks to 2nd level tasks (e.g., num_dims)
    end if
    call p_trian%setup_coarse_triangulation()

  end subroutine par_mesh_to_triangulation

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
             if (p_trian%vefs(jobj)%itfc == -1) then
                p_trian%num_itfc_vefs = p_trian%num_itfc_vefs + 1
                p_trian%vefs(jobj)%itfc = p_trian%num_itfc_vefs
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
             if (p_trian%vefs(jobj)%itfc == (p_trian%num_itfc_vefs + 1)) then
                p_trian%num_itfc_vefs = p_trian%num_itfc_vefs + 1
                p_trian%lst_itfc_vefs(p_trian%num_itfc_vefs) = jobj 
             end if
          end if
       end do
    end do

    ! Compute p_trian%max_nparts, p_trian%nobjs, p_trian%lobjs
    ! Re-order (permute) p_trian%lst_itfc_vefs accordingly
    call p_trian%compute_vefs_and_parts_object()
    call p_trian%compute_objects_dimension()
    call p_trian%compute_number_global_objects_and_their_gids()
  end subroutine par_triangulation_to_dual
  
  ! Required (at least) by type(mlbddc_t)
  function par_triangulation_get_coarse_triangulation(this)
    implicit none
    class(par_triangulation_t), target, intent(in)  :: this
    type(coarse_triangulation_t), pointer :: par_triangulation_get_coarse_triangulation
    par_triangulation_get_coarse_triangulation => this%coarse_triangulation
  end function par_triangulation_get_coarse_triangulation
 
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
       vefs_object_iterator = this%vefs_object%create_iterator(iobj)
       parts_object_iterator = this%parts_object%create_iterator(iobj)
       
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
  
  subroutine par_triangulation_compute_objects_dimension(this)
    implicit none
    class(par_triangulation_t), intent(inout) :: this
    integer(ip)                               :: iobj, vef_lid
    type(list_iterator_t)                     :: vefs_object_iterator
    
    call memalloc ( this%number_objects, this%objects_dimension, __FILE__, __LINE__ )
    do iobj=1, this%number_objects
      vefs_object_iterator = this%vefs_object%create_iterator(iobj)
      this%objects_dimension(iobj) = 0
      do while(.not. vefs_object_iterator%is_upper_bound())
        vef_lid = vefs_object_iterator%get_current()
        this%objects_dimension(iobj) = max ( this%objects_dimension(iobj), this%triangulation%vefs(vef_lid)%dime )  
        call vefs_object_iterator%next()
       end do
    end do
  end subroutine par_triangulation_compute_objects_dimension
  
  
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
       parts_object_iterator = this%parts_object%create_iterator(iobj)
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
       parts_object_iterator = this%parts_object%create_iterator(iobj)
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
      parts_object_iterator = this%parts_object%create_iterator(iobj)
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
  
  subroutine par_triangulation_setup_coarse_triangulation ( this )
    implicit none
    class(par_triangulation_t), intent(inout) :: this
    integer(igp)              , allocatable   :: coarse_cell_gids(:)
    integer(ip)               , allocatable   :: coarse_vefs_recv_counts(:)
    integer(ip)               , allocatable   :: coarse_vefs_displs(:)
    integer(igp)              , allocatable   :: lst_coarse_vef_gids(:)
    integer(ip)               , allocatable   :: lst_coarse_vef_dimension(:)
    integer(ip)               , allocatable   :: l2_part_id_neighbours(:)
    integer(ip)               , allocatable   :: coarse_dgraph_recv_counts(:)
    integer(ip)               , allocatable   :: coarse_dgraph_displs(:)
    integer(ip)               , allocatable   :: lextn(:)
    integer(ip)               , allocatable   :: lextp(:)
    
    integer(ip)                      :: i, istat
    integer(ip)                      :: num_dimensions
    integer(ip)                      :: num_local_coarse_cells
    integer(ip)                      :: num_itfc_coarse_cells
    
    ! All MPI tasks (even if they are not involved in the L2 from L1 gather) should also allocate the
    ! allocatable arrays due to the fact that non-allocated allocatable arrays cannot
    ! be passed as actual arguments of dummy arguments that do not have the allocatable attribute 
    ! (see e.g. coarse_triangulation%create() below). Otherwise, the code crashes with a segmentation fault. 
    ! Likewise, actual arguments which are used as input dummy arguments to size another array-type dummy arguments should also
    ! be initialized on all MPI tasks
    num_local_coarse_cells = 0
    num_itfc_coarse_cells  = 0
    call memalloc (0, coarse_cell_gids, __FILE__, __LINE__)
    call memalloc (0, coarse_vefs_recv_counts, __FILE__, __LINE__)
    call memalloc (0, coarse_vefs_displs, __FILE__, __LINE__)
    call memalloc (0, lst_coarse_vef_gids, __FILE__, __LINE__)
    call memalloc (0, lst_coarse_vef_dimension, __FILE__, __LINE__)
    call memalloc (0, l2_part_id_neighbours, __FILE__, __LINE__)
    call memalloc (0, coarse_dgraph_recv_counts, __FILE__, __LINE__)
    call memalloc (0, coarse_dgraph_displs, __FILE__, __LINE__)
    call memalloc (0, lextn, __FILE__, __LINE__)
    call memalloc (0, lextp, __FILE__, __LINE__)
        
    ! L2 tasks gather from L1 tasks all raw data required to set-up the coarse triangulation on L2 tasks
    if ( this%p_env%am_i_l1_to_l2_task() ) then
      call this%p_env%l1_to_l2_transfer ( this%triangulation%num_dims, num_dimensions ) 
      call this%gather_coarse_cell_gids (coarse_cell_gids)
      call this%gather_coarse_vefs_rcv_counts_and_displs (coarse_vefs_recv_counts, coarse_vefs_displs)
      call this%gather_coarse_vefs_gids (coarse_vefs_recv_counts, coarse_vefs_displs, lst_coarse_vef_gids)
      call this%gather_coarse_vefs_dimension (coarse_vefs_recv_counts, coarse_vefs_displs, lst_coarse_vef_dimension)
      call this%fetch_l2_part_id_neighbours(l2_part_id_neighbours)
      call this%gather_coarse_dgraph_rcv_counts_and_displs ( l2_part_id_neighbours, &
                                                             coarse_dgraph_recv_counts, &
                                                             coarse_dgraph_displs )
      call this%gather_coarse_dgraph_lextn_and_lextp ( l2_part_id_neighbours, &
                                                       coarse_dgraph_recv_counts, &
                                                       coarse_dgraph_displs, &
                                                       lextn, &
                                                       lextp )
      ! Evaluate number of local coarse cells
      num_local_coarse_cells = this%p_env%get_l1_to_l2_size()-1
      
      ! Evaluate number of interface coarse cells
      ! Adapt and re-use coarse_vefs_displs/coarse_dgraph_recv_counts/coarse_dgraph_displs
      ! as required by this%coarse_triangulation%create below
      num_itfc_coarse_cells = this%adapt_coarse_raw_arrays (coarse_vefs_displs, &
                                                            coarse_dgraph_recv_counts, &
                                                            coarse_dgraph_displs )
    end if
    
    if ( this%p_env%am_i_lgt1_task() ) then
      ! lgt1 MPI tasks (recursively) build coarse triangulation
      allocate  ( this%coarse_triangulation, stat = istat )
      check( istat == 0 )
      call this%coarse_triangulation%create ( par_environment              = this%p_env%get_next_level(), &
                                              num_dimensions               = num_dimensions, &
                                              num_local_cells              = num_local_coarse_cells, &
                                              cells_gid                    = coarse_cell_gids, &
                                              ptr_vefs_per_cell            = coarse_vefs_displs, &
                                              lst_vefs_gids                = lst_coarse_vef_gids, &
                                              lst_vefs_dimension           = lst_coarse_vef_dimension, &
                                              num_itfc_cells               = num_itfc_coarse_cells, &
                                              lst_itfc_cells               = coarse_dgraph_recv_counts, &
                                              ptr_ext_neighs_per_itfc_cell = coarse_dgraph_displs, &
                                              lst_ext_neighs_gids          = lextn, &
                                              lst_ext_neighs_part_ids      = lextp )
    else
      ! L1 tasks do not hold any piece of the coarse triangulation
      nullify(this%coarse_triangulation)
    end if
    
    ! All tasks free raw data (see actual reason on the top part of this subroutine)
    call memfree (coarse_cell_gids, __FILE__, __LINE__)
    call memfree (coarse_vefs_recv_counts, __FILE__, __LINE__)
    call memfree (coarse_vefs_displs, __FILE__, __LINE__)
    call memfree (lst_coarse_vef_gids, __FILE__, __LINE__)
    call memfree (lst_coarse_vef_dimension, __FILE__, __LINE__)
    call memfree (l2_part_id_neighbours, __FILE__, __LINE__)
    call memfree (coarse_dgraph_recv_counts, __FILE__, __LINE__)
    call memfree (coarse_dgraph_displs, __FILE__, __LINE__)
    call memfree (lextn, __FILE__, __LINE__)
    call memfree (lextp, __FILE__, __LINE__)
  end subroutine par_triangulation_setup_coarse_triangulation
  
  subroutine par_triangulation_gather_coarse_cell_gids( this, coarse_cell_gids)
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(igp) , allocatable , intent(inout) :: coarse_cell_gids(:)
    
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(igp)                              :: dummy_integer_array(0)
    
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if ( allocated (coarse_cell_gids) ) call memfree ( coarse_cell_gids, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, coarse_cell_gids, __FILE__, __LINE__ )
      call this%p_env%l2_from_l1_gather( input_data  = 0_igp, &
                                         output_data = coarse_cell_gids ) 
    else
      call this%p_env%l2_from_l1_gather( input_data  = int(this%p_env%get_l1_rank()+1,igp), &
                                         output_data = dummy_integer_array ) 
    end if
  end subroutine par_triangulation_gather_coarse_cell_gids
  
  
  subroutine par_triangulation_gather_coarse_vefs_rcv_counts_and_displs( this, recv_counts, displs )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip)                               :: i
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array(0)

    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if ( allocated (recv_counts) ) call memfree ( recv_counts, __FILE__, __LINE__ )
      if ( allocated (displs) ) call memfree ( displs, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
      call this%p_env%l2_from_l1_gather( input_data = 0, &
                                         output_data = recv_counts ) 
      displs(1) = 0
      do i=2, l1_to_l2_size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
    else
      call this%p_env%l2_from_l1_gather( input_data  = this%number_objects, &
                                         output_data = dummy_integer_array ) 
    end if
  end subroutine par_triangulation_gather_coarse_vefs_rcv_counts_and_displs
  
  subroutine par_triangulation_gather_coarse_vefs_gids ( this, recv_counts, displs, lst_gids )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: recv_counts(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(in)    :: displs(this%p_env%get_l1_to_l2_size())
    integer(igp), allocatable , intent(inout) :: lst_gids(:)
    integer(ip)                               :: l1_to_l2_size
    integer(ip)                               :: dummy_integer_array_ip(0)
    integer(igp)                              :: dummy_integer_array_igp(0)
    
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if (allocated(lst_gids)) call memfree ( lst_gids, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lst_gids, __FILE__, __LINE__ )
      call this%p_env%l2_from_l1_gather( input_data_size = 0, &
                                         input_data      = dummy_integer_array_igp, &
                                         recv_counts     = recv_counts, &
                                         displs          = displs, &
                                         output_data     = lst_gids )
    else
      call this%p_env%l2_from_l1_gather( input_data_size = this%number_objects, &
                                         input_data      = this%objects_gids, &
                                         recv_counts     = dummy_integer_array_ip, &
                                         displs          = dummy_integer_array_ip, &
                                         output_data     = dummy_integer_array_igp )
    end if    
  end subroutine par_triangulation_gather_coarse_vefs_gids
  
  subroutine par_triangulation_gather_coarse_vefs_dimension ( this, recv_counts, displs, lst_vefs_dimension )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip)                  , intent(in)    :: recv_counts(this%p_env%get_l1_to_l2_size())
    integer(ip)                  , intent(in)    :: displs(this%p_env%get_l1_to_l2_size())
    integer(ip), allocatable     , intent(inout) :: lst_vefs_dimension(:)
    integer(ip)                                  :: l1_to_l2_size
    integer(ip)                                  :: dummy_integer_array(0)
    
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if (allocated(lst_vefs_dimension)) call memfree ( lst_vefs_dimension, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lst_vefs_dimension, __FILE__, __LINE__ )
      call this%p_env%l2_from_l1_gather( input_data_size = 0, &
                                         input_data      = dummy_integer_array, &
                                         recv_counts     = recv_counts, &
                                         displs          = displs, &
                                         output_data     = lst_vefs_dimension )
    else
      call this%p_env%l2_from_l1_gather( input_data_size = this%number_objects, &
                                         input_data      = this%objects_dimension, &
                                         recv_counts     = dummy_integer_array, &
                                         displs          = dummy_integer_array, &
                                         output_data     = dummy_integer_array )
    end if    
  end subroutine par_triangulation_gather_coarse_vefs_dimension
  
  
  subroutine par_triangulation_fetch_l2_part_id_neighbours ( this, l2_part_id_neighbours )    
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip) , allocatable , intent(inout) :: l2_part_id_neighbours(:)
    integer(ip) :: my_l2_part_id
    integer(ip) :: num_neighbours
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if (this%p_env%am_i_l1_task()) then
      num_neighbours = this%element_import%get_number_neighbours()
      my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
      if (allocated(l2_part_id_neighbours)) call memfree ( l2_part_id_neighbours, __FILE__, __LINE__ )
      call memalloc ( num_neighbours, l2_part_id_neighbours, __FILE__, __LINE__ )
      call this%p_env%l1_neighbours_exchange ( num_neighbours  = num_neighbours, &
                                               list_neighbours = this%element_import%get_neighbours_ids(), &
                                               input_data      = my_l2_part_id,&
                                               output_data     = l2_part_id_neighbours)
    end if
  end subroutine par_triangulation_fetch_l2_part_id_neighbours
  
  subroutine par_triangulation_gather_coarse_dgraph_rcv_counts_and_displs ( this, &
                                                                            l2_part_id_neighbours, &
                                                                            recv_counts, &
                                                                            displs )
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: l2_part_id_neighbours(this%element_import%get_number_neighbours())
    integer(ip) , allocatable , intent(inout) :: recv_counts(:) 
    integer(ip) , allocatable , intent(inout) :: displs(:)
    integer(ip) :: i
    integer(ip) :: l1_to_l2_size 
    integer(ip) :: my_l2_part_id
    integer(ip) :: num_neighbours
    integer(ip) :: num_external_l2_elements
    integer(ip) :: dummy_integer_array(0)

    
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if (allocated(recv_counts)) call memfree ( recv_counts, __FILE__, __LINE__ )
      if (allocated(displs)) call memfree ( displs, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, recv_counts, __FILE__, __LINE__ )
      call memalloc ( l1_to_l2_size, displs, __FILE__, __LINE__ )
      call this%p_env%l2_from_l1_gather( input_data = 0, &
                                         output_data = recv_counts ) 
      displs(1) = 0
      do i=2, l1_to_l2_size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
    else
      assert ( this%p_env%am_i_l1_task() )
      num_neighbours = this%element_import%get_number_neighbours()
      my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
      num_external_l2_elements = 0
      do i = 1, num_neighbours
        if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
          num_external_l2_elements = num_external_l2_elements + 1
        end if
      end do
      call this%p_env%l2_from_l1_gather( input_data = num_external_l2_elements, &
                                         output_data = dummy_integer_array ) 
    end if
  end subroutine par_triangulation_gather_coarse_dgraph_rcv_counts_and_displs 
  
  subroutine par_triangulation_gather_coarse_dgraph_lextn_and_lextp( this,                  & 
                                                                     l2_part_id_neighbours, &
                                                                     recv_counts,           &
                                                                     displs,                &
                                                                     lextn,                 &
                                                                     lextp)
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(in)    :: l2_part_id_neighbours(this%element_import%get_number_neighbours())
    integer(ip)               , intent(in)    :: recv_counts(this%p_env%get_l1_to_l2_size()) 
    integer(ip)               , intent(in)    :: displs(this%p_env%get_l1_to_l2_size())
    integer(ip), allocatable  , intent(inout) :: lextn(:)
    integer(ip), allocatable  , intent(inout) :: lextp(:)
    
    integer(ip)              :: i
    integer(ip)              :: l1_to_l2_size 
    integer(ip)              :: my_l2_part_id
    integer(ip)              :: num_neighbours
    integer(ip), pointer     :: neighbours_ids(:)
    integer(ip)              :: num_external_l2_elements
    integer(ip), allocatable :: lst_external_l2_element_gids(:)
    integer(ip), allocatable :: lst_external_l2_part_ids(:)
    integer(ip)              :: dummy_integer_array(0)
    
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      l1_to_l2_size = this%p_env%get_l1_to_l2_size()
      if (allocated(lextn)) call memfree ( lextn, __FILE__, __LINE__ )
      if (allocated(lextp)) call memfree ( lextp, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lextn, __FILE__, __LINE__ )
      call memalloc ( displs(l1_to_l2_size), lextp, __FILE__, __LINE__ )
      ! Gather lextn
      call this%p_env%l2_from_l1_gather( input_data_size = 0, &
                                         input_data      = dummy_integer_array, &
                                         recv_counts     = recv_counts, &
                                         displs          = displs, &
                                         output_data     = lextn )
      ! Gather lextp
      call this%p_env%l2_from_l1_gather( input_data_size = 0, &
                                         input_data      = dummy_integer_array, &
                                         recv_counts     = recv_counts, &
                                         displs          = displs, &
                                         output_data     = lextp )
    else
      assert ( this%p_env%am_i_l1_task() )
      num_neighbours =  this%element_import%get_number_neighbours()
      neighbours_ids => this%element_import%get_neighbours_ids()
      my_l2_part_id  = this%p_env%get_l2_part_id_l1_task_is_mapped_to()
      num_external_l2_elements = 0
      do i = 1, num_neighbours
        if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
          num_external_l2_elements = num_external_l2_elements + 1
        end if
      end do
      
      call memalloc (num_external_l2_elements, lst_external_l2_part_ids, __FILE__, __LINE__)
      call memalloc (num_external_l2_elements, lst_external_l2_element_gids,__FILE__, __LINE__)
      num_external_l2_elements = 0
      neighbours_ids => this%element_import%get_neighbours_ids()
      do i = 1, num_neighbours
        if ( my_l2_part_id /= l2_part_id_neighbours(i) ) then
          num_external_l2_elements = num_external_l2_elements + 1
          lst_external_l2_element_gids(num_external_l2_elements) = neighbours_ids(i)
          lst_external_l2_part_ids(num_external_l2_elements) = l2_part_id_neighbours(i)
        end if
      end do
      call this%p_env%l2_from_l1_gather( input_data_size = num_external_l2_elements, &
                                         input_data      = lst_external_l2_element_gids, &
                                         recv_counts     = dummy_integer_array, &
                                         displs          = dummy_integer_array, &
                                         output_data     = dummy_integer_array )
      
      call this%p_env%l2_from_l1_gather( input_data_size = num_external_l2_elements, &
                                         input_data      = lst_external_l2_part_ids, &
                                         recv_counts     = dummy_integer_array, &
                                         displs          = dummy_integer_array, &
                                         output_data     = dummy_integer_array )
      
      call memfree (lst_external_l2_part_ids   , __FILE__, __LINE__)
      call memfree (lst_external_l2_element_gids,__FILE__, __LINE__)
    end if
  end subroutine par_triangulation_gather_coarse_dgraph_lextn_and_lextp 
  
  function par_triangulation_adapt_coarse_raw_arrays( this, &
                                                      coarse_vefs_displs, &
                                                      coarse_dgraph_recv_counts, &
                                                      coarse_dgraph_displs ) result(num_itfc_coarse_cells)
    implicit none
    class(par_triangulation_t), intent(in)    :: this
    integer(ip)               , intent(inout) :: coarse_vefs_displs(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(inout) :: coarse_dgraph_recv_counts(this%p_env%get_l1_to_l2_size())
    integer(ip)               , intent(inout) :: coarse_dgraph_displs(this%p_env%get_l1_to_l2_size())
    integer(ip)                               :: num_itfc_coarse_cells
    
    integer(ip) :: i 
    assert ( this%p_env%am_i_l1_to_l2_task() )
    if ( this%p_env%am_i_l1_to_l2_root() ) then
      ! Re-use coarse_vefs_displs as ptr_vefs_gids
      do i=1, size(coarse_vefs_displs)
        coarse_vefs_displs(i)=coarse_vefs_displs(i)+1
      end do
   
      ! Re-use coarse_dgraph_displs as ptr_ext_neighs_per_itfc_cell
      num_itfc_coarse_cells = 0
      coarse_dgraph_displs(1) = 1 
      do i=1, size(coarse_dgraph_recv_counts)
        if (coarse_dgraph_recv_counts(i) /= 0) then
          num_itfc_coarse_cells = num_itfc_coarse_cells+1
          coarse_dgraph_displs(num_itfc_coarse_cells+1) = coarse_dgraph_displs(num_itfc_coarse_cells) + &
                                                          coarse_dgraph_recv_counts(i)                                           
        end if  
      end do
      
      ! Re-use coarse_dgraph_recv_counts as lst_itfc_cells
      num_itfc_coarse_cells = 0
      do i=1, size(coarse_dgraph_recv_counts)
        if (coarse_dgraph_recv_counts(i) /= 0) then
          num_itfc_coarse_cells = num_itfc_coarse_cells+1
          coarse_dgraph_recv_counts(num_itfc_coarse_cells) = i
        end if  
      end do
    else
      ! L1 tasks do not hold any itfc_coarse_cells
      num_itfc_coarse_cells = 0
    end if
  end function par_triangulation_adapt_coarse_raw_arrays
  
  !=============================================================================
  ! Auxiliary subroutines
  subroutine initialize_par_vef_topology (vef)
    implicit none
    type(par_vef_topology_t), intent(inout) :: vef
    
    vef%itfc   = -1
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
    element%itfc = -1
    element%mypart = -1 
  end subroutine free_par_elem_topology
  
  subroutine initialize_par_elem_topology(element)
    implicit none
    type(par_elem_topology_t), intent(inout) :: element

    assert(allocated(element%vefs_GIDs))
    element%itfc      = -1
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
    buffer(start:end) = transfer(my%itfc,mold)

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
    my%itfc  = transfer(buffer(start:end), my%itfc)

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
