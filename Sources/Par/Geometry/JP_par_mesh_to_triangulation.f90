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
module JP_par_mesh_to_triangulation_names
  ! Serial modules
  use types_names
  use memor_names
  use JP_triangulation_names
  use element_id_names
  use migratory_element_names
  use element_import_names
  use element_import_create_names
  use hash_table_names
  use JP_element_topology_names
  use JP_mesh_to_triangulation_names
  use psi_penv_mod_names

  ! Parallel modules
  use par_element_topology_names
  use JP_par_triangulation_names
  use par_mesh_names
  use par_element_exchange_names
  use par_conditions_names

  implicit none
# include "debug.i90"
  private

  public :: JP_par_mesh_to_triangulation

contains

  !*********************************************************************************
  ! This subroutine takes as input a 'plain' parallel mesh and creates a parallel 
  ! triangulation. The code only works for JP_par_triangulation_t because the acces
  ! to the data is based on an index ielem using the functions set_index and get_index
  !*********************************************************************************
  subroutine JP_par_mesh_to_triangulation (p_gmesh, p_trian)
    implicit none
    ! Parameters
    type(par_mesh_t)             , target  , intent(in)   :: p_gmesh
    type(JP_par_triangulation_t), target  , intent(inout) :: p_trian 

    ! Locals
    integer(ip) :: istat, ielem, iobj, jobj
    integer(ip) :: num_elems, num_ghost, num_verts
    type (hash_table_igp_ip_t) :: hash ! Topological info hash table (SBmod)
    integer(ip) :: ilele, nvert, jelem, jlele, idime, count, ivere 
    integer(igp):: iobjg
    integer(igp), allocatable :: aux_igp(:)
    integer(ip) , allocatable :: aux(:)
    integeR(ip) :: aux_val

    class(migratory_element_t)  , pointer :: migratory_elements(:) => NULL()
    type(par_element_topology_t), pointer :: elements(:) => NULL()
    type(par_vef_topology_t)    , pointer :: vefs(:)     => NULL()

    !class(par_element_topology_t), pointer     :: elem, neighbor
    !class(element_id_t)          , allocatable :: elem_id

    ! Set a reference to the type(par_environment_t) instance describing the set of MPI tasks
    ! among which this type(par_triangulation) instance is going to be distributed 
    p_trian%p_env => p_gmesh%p_env

    if(p_trian%p_env%p_context%iam >= 0) then

       if ( p_trian%state /= JP_triangulation_not_created ) call p_trian%free

       ! Create element_import from geometry mesh partition data
       ! AFM: CURRENTLY element_import_create is the only way to create a type(element_import) instance.
       !      I have stored it inside type(par_triangulation) as I do not have a better guess.
       !      In the future, we should get rid of finite_element_t_import_create, and provide a new
       !      subroutine which allows to create this instance using the dual graph and the gluing
       !      data describing its distributed-memory layout. Both the dual graph and associated gluing
       !      data are to be stored in type(par_neighborhood) according to Javier's UML diagram (i.e., fempar.dia). 
       !      From this point of view, type(par_neighborhood) should somehow aggregate/reference an instance of 
       !      type(element_import) and manage its creation. However, this is not currently reflected in fempar.dia, 
       !      which defines a type(triangulation_partition) which in turn includes an instance of type(element_import) inside. 
       !      Assuming we agree in the first option, how type(par_triangulation) is going to access type(par_neighbourhood) ??? 
       !      This is related with a parallel discussion about the possibility of enriching type(par_triangulation) with the dual graph 
       !      and associated gluing data. Does it make sense? If yes, does type(par_neighborhood) still makes any sense?
       ! JP:  Apart from the previous discussion, an important point is that the element_import is CURRENTLY concrete, i.e.
       !      it can only be used with a par_triangulation based on a static set. We need to provide a new implementation
       !      based on abstract element_id_t instead of global ids to be able to use it with adaptive triangultions.
       call element_import_create ( p_gmesh%f_mesh_dist, p_trian%f_el_import )

       ! Create triangulation and fill it with local data
       num_elems = p_trian%f_el_import%nelem
       num_ghost = p_trian%f_el_import%nghost

       call p_trian%create(num_elems + num_ghost)
       call JP_mesh_to_triangulation( p_gmesh%f_mesh, p_trian )

       p_trian%num_ghost = num_ghost
       p_trian%num_elems = num_elems

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       ! From now onwards, the code is valid for a par_triangulation based on a plain element_set. It is the legacy
       ! code with minor modification (renaming of elems)
       !
       select type( set => p_trian%element_set )
       class is(plain_migratory_element_set_t)
          migratory_elements => set%get_elements()
          select type( migratory_elements )
          class is(par_element_topology_t)
             elements => migratory_elements
          class default
             write(*,*) 'The procedure par_mesh_to_triangulation can only work with elements of type par_element_topology'
             check(.false.)
          end select
       class default
          write(*,*) 'The procedure par_mesh_to_triangulation can only work with a par_triangulation defined'
          write(*,*) 'using a plain element set'
          check(.false.)
       end select
       select type( item => p_trian%vefs)
       class is(par_vef_topology_t)
          vefs => item
       class default
          write(*,*) 'The procedure par_mesh_to_triangulation can only work with vefs of type par_vef_topology'
          check(.false.)
       end select


       elements(:)%interface = -1 
       do ielem=1, p_gmesh%f_mesh_dist%nebou
          elements(p_gmesh%f_mesh_dist%lebou(ielem))%interface = ielem
       end do

       p_trian%num_itfc_elems = p_gmesh%f_mesh_dist%nebou
       call memalloc( p_trian%num_itfc_elems, p_trian%lst_itfc_elems, __FILE__, __LINE__ )
       p_trian%lst_itfc_elems = p_gmesh%f_mesh_dist%lebou


       ! **IMPORTANT NOTE**: the code that comes assumes that edges and faces in p_trian are numbered after vertices.
       ! This requirement is currently fulfilled by mesh_to_triangulation (in particular, by geom2topo within) but we should
       ! keep this in mind all the way through. If this were not assumed, we should find a way to identify a corner within
       ! p_trian, and map from a corner local ID in p_trian to a vertex local ID in p_gmesh%f_mesh.
       num_verts = p_gmesh%f_mesh%npoin

       ! Fill array of elements (local ones)
       do ielem=1, num_elems
          elements(ielem)%mypart      = p_trian%p_env%p_context%iam + 1
          elements(ielem)%globalID    = p_gmesh%f_mesh_dist%emap%l2g(ielem)
          elements(ielem)%num_vefs = elements(ielem)%num_vefs
          call memalloc( elements(ielem)%num_vefs, elements(ielem)%vefs_GIDs, __FILE__, __LINE__ )
          do iobj=1, elements(ielem)%num_vefs
             jobj = elements(ielem)%vefs(iobj)
             if ( jobj <= num_verts ) then ! It is a corner => re-use global ID
                elements(ielem)%vefs_GIDs(iobj) = p_gmesh%f_mesh_dist%nmap%l2g(jobj)
             else ! It is an edge or face => generate new local-global ID (non-consistent, non-consecutive)
                ! The ISHFT(1,50) is used to start numbering efs after vertices, assuming nvert < 2**60
                elements(ielem)%vefs_GIDs(iobj) = ISHFT(int(p_gmesh%f_mesh_dist%ipart,igp),int(32,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(60,igp))
                !elements(ielem)%vefs_GIDs(iobj) = ISHFT(int(p_gmesh%f_mesh_dist%ipart,igp),int(6,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(6,igp))
             end if
          end do
       end do

       ! Get vefs_GIDs from ghost elements and create reference elements
       call ghost_elements_exchange ( p_trian%p_env%p_context%icontxt, p_trian%f_el_import, elements )

       ! Create reference elements, both in the local and ghost elements
       call create_reference_elements( p_trian )

       ! Hash table global to local for ghost elements
       call hash%init(num_ghost)
       do ielem=num_elems+1, num_elems+num_ghost
          aux_val = ielem
          call hash%put( key = elements(ielem)%globalID, val = aux_val, stat=istat)
       end do

       do ielem = 1, p_gmesh%f_mesh_dist%nebou     ! Loop interface elements 
          ! Step 1: Put LID of vertices in the ghost elements (f_mesh_dist)
          ilele = p_gmesh%f_mesh_dist%lebou(ielem) ! local ID element
          ! aux : array of ilele (LID) vertices in GID
          nvert  = elements(ilele)%geo_reference_element%nvef_dim(2)-1
          call memalloc( nvert, aux_igp, __FILE__, __LINE__  )
          do iobj = 1, nvert                        ! vertices only
             aux_igp(iobj) = elements(ilele)%vefs_GIDs(iobj) ! extract GIDs vertices
          end do
          do jelem = p_gmesh%f_mesh_dist%pextn(ielem), & 
               p_gmesh%f_mesh_dist%pextn(ielem+1)-1  ! external neighbor elements
             call hash%get(key = p_gmesh%f_mesh_dist%lextn(jelem), val=jlele, stat=istat) ! LID external element
             do jobj = 1, elements(jlele)%geo_reference_element%nvef_dim(2)-1 ! vertices external 
                if ( elements(jlele)%vefs(jobj) == -1) then
                   do iobj = 1, nvert
                      if ( aux_igp(iobj) == elements(jlele)%vefs_GIDs(jobj) ) then
                         ! Put LID of vertices for ghost_elements
                         elements(jlele)%vefs(jobj) =  elements(ilele)%vefs(iobj)
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
             ! loop over all efs of external elements
             do idime =2,p_trian%num_dims
                do iobj = elements(jlele)%geo_reference_element%nvef_dim(idime), &
                     elements(jlele)%geo_reference_element%nvef_dim(idime+1)-1 
                   if ( elements(jlele)%vefs(iobj) == -1) then ! efs not assigned yet
                      count = 1
                      ! loop over vertices of every ef
                      do jobj = elements(jlele)%geo_reference_element%crxob%p(iobj), &
                           elements(jlele)%geo_reference_element%crxob%p(iobj+1)-1    
                         ivere = elements(jlele)%geo_reference_element%crxob%l(jobj)
                         if (elements(jlele)%vefs(ivere) == -1) then
                            count = 0 ! not an vef of the local triangulation
                            exit
                         end if
                      end do
                      if (count == 1) then
                         nvert = elements(jlele)%geo_reference_element%crxob%p(iobj+1)- &
                              elements(jlele)%geo_reference_element%crxob%p(iobj)
                         call memalloc( nvert, aux, __FILE__, __LINE__)
                         count = 1
                         do jobj = elements(jlele)%geo_reference_element%crxob%p(iobj), &
                              elements(jlele)%geo_reference_element%crxob%p(iobj+1)-1 
                            ivere = elements(jlele)%geo_reference_element%crxob%l(jobj)
                            aux(count) = elements(jlele)%vefs(ivere)
                            count = count+1
                         end do
                         call local_id_from_vertices( elements(ilele), idime, aux, nvert, &
                              elements(jlele)%vefs(iobj) )
                         call memfree(aux, __FILE__, __LINE__) 
                      end if
                   end if
                end do
             end do
          end do
       end do

       call hash%free

       !call par_triangulation_to_dual ( p_trian )

       !pause

       ! Check results
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do ielem = 1,num_elems+num_ghost
       !       write (*,*) '****ielem:',ielem          
       !       write (*,*) '****LID_vefs ghost:',elements(ielem)%vefs
       !       write (*,*) '****GID_vefs ghost:',elements(ielem)%vefs_GIDs
       !    end do
       ! end if
       ! pause

       ! write (*,*) '*********************************'
       ! write (*,*) '*********************************' 
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do iobj = 1, p_trian%num_vefs 
       !       write (*,*) 'is interface vef',iobj, ' :', vefs(iobj)%interface 
       !       write (*,*) 'is interface vef',iobj, ' :', vefs(iobj)%dimension
       !    end do
       ! end if
       ! write (*,*) '*********************************'
       ! write (*,*) '*********************************'

       ! Step 3: Make GID consistent among processors (p_part%elems%vefs_GIDs)
       do iobj = 1, p_trian%num_vefs 
          if ( (vefs(iobj)%interface .ne. -1) .and. &
               (vefs(iobj)%dimension >= 1) ) then
             idime = vefs(iobj)%dimension+1
             iobjg = -1

             do jelem = 1,vefs(iobj)%num_elems_around
                jlele = vefs(iobj)%elems_around(jelem)%get_index()

                do jobj = elements(jlele)%geo_reference_element%nvef_dim(idime), &
                     & elements(jlele)%geo_reference_element%nvef_dim(idime+1)-1 ! efs of neighbor els
                   if ( elements(jlele)%vefs(jobj) == iobj ) then
                      if ( iobjg == -1 ) then 
                         iobjg  = elements(jlele)%vefs_GIDs(jobj)
                      else
                         iobjg  = min(iobjg,elements(jlele)%vefs_GIDs(jobj))
                         exit
                      end if
                   end if
                end do

             end do


             do jelem = 1,vefs(iobj)%num_elems_around
                jlele = vefs(iobj)%elems_around(jelem)%get_index()
                do jobj = elements(jlele)%geo_reference_element%nvef_dim(idime), &
                     & elements(jlele)%geo_reference_element%nvef_dim(idime+1)-1 ! efs of neighbor els
                   if ( elements(jlele)%vefs(jobj) == iobj) then
                      elements(jlele)%vefs_GIDs(jobj) = iobjg
                      exit
                   end if
                end do
             end do
             vefs(iobj)%globalID = iobjg          
          end if
       end do

       p_trian%state = JP_triangulation_vefs_filled

       ! Check results
       ! if ( p_trian%p_env%p_context%iam == 0) then
       !    do ielem = 1,num_elems+num_ghost
       !       write (*,*) '****ielem:',ielem          
       !       write (*,*) '****LID_vefs ghost:',elements(ielem)%vefs
       !       write (*,*) '****GID_vefs ghost:',elements(ielem)%vefs_GIDs
       !    end do
       ! end if
       ! pause
    else
       ! AFM: TODO: Partially broadcast p_trian from 1st level tasks to 2nd level tasks (e.g., num_dims)
    end if

  end subroutine JP_par_mesh_to_triangulation

end module JP_par_mesh_to_triangulation_names
