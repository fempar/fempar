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
module par_mesh_triangulation
  ! Serial modules
  use types
  use memor
  use fem_triangulation_class
  use fem_partition_class
  use hash_table_class
  use mesh_triangulation
  use psi_penv_mod

  ! Parallel modules
  use par_triangulation_class
  use par_mesh_class
  use par_partition_class
  use par_element_exchange

# include "debug.i90"
  implicit none
  private

  public :: par_mesh_to_triangulation

contains

  subroutine par_mesh_to_triangulation (p_gmesh, p_trian)
    implicit none
    ! Parameters
    type(par_mesh)         , target, intent(inout) :: p_gmesh ! Geometry mesh
    type(par_triangulation), target, intent(inout) :: p_trian 

    ! Locals
    integer(ip) :: istat, ielem, iobj, jobj, state
    integer(ip) :: num_elems, num_ghosts, num_verts
    type (hash_table_igp_ip)            :: hash ! Topological info hash table (SBmod)
    integer(ip) :: ilele, nvert, jelem, jlele, idime, count, ivere 
    integer(igp):: iobjg
    integer(igp), allocatable :: aux_igp(:)
    integer(ip), allocatable :: aux(:)
    !     integer(ip), allocatable :: is_eboun(:)

    state = p_trian%f_trian%state

    assert(state == triangulation_not_created .or. state == triangulation_elems_filled .or. state == triangulation_elems_objects_filled)
    assert(p_gmesh%p_part%f_part%pinfo == extended_adjacencies )

    if ( state == triangulation_elems_objects_filled ) then
       do iobj=1, p_trian%f_trian%num_objects 
          call free_par_object_topology(p_trian%objects(iobj)) 
       end do

       p_trian%num_itfc_objs = -1
       call memfree ( p_trian%lst_itfc_objs, __FILE__, __LINE__ )

       ! Deallocate the object structure array 
       deallocate(p_trian%objects, stat=istat)
       check(istat==0)
    end if

    if ( state == triangulation_elems_objects_filled .or. & 
         state == triangulation_elems_filled ) then
       do ielem=1, p_trian%f_trian%elem_array_len 
          call free_par_elem_topology(p_trian%elems(iobj)) 
       end do
       ! Deallocate the element structure array */
       deallocate(p_trian%mig_elems, stat=istat)
       check(istat==0)
       nullify(p_trian%elems)

       call fem_triangulation_free ( p_trian%f_trian )
    end if

    ! Now we are sure that the local portion of p_trian, i.e., p_trian%f_trian is in triangulation_not_created state
    ! Let's create and fill it
    num_elems  = p_gmesh%p_part%f_el_import%nelem
    num_ghosts = p_gmesh%p_part%f_el_import%nghost

    p_trian%num_ghosts = num_ghosts
    p_trian%num_elems  = num_elems

    call fem_triangulation_create ( num_elems + num_ghosts, p_trian%f_trian )

    ! Fill local portion with local data
    call mesh_to_triangulation ( p_gmesh%f_mesh, p_trian%f_trian )
    ! p_trian%f_trian%num_elems = num_elems+num_ghosts
    ! **IMPORTANT NOTE**: the code that comes assumes that edges and faces in p_trian%f_trian are numbered after vertices.
    ! This requirement is currently fulfilled by mesh_to_triangulation (in particular, by geom2topo within) but we should
    ! keep this in mind all the way through. If this were not assumed, we should find a way to identify a corner within
    ! p_trian%f_trian, and map from a corner local ID in p_trian%f_trian to a vertex local ID in p_gmesh%f_mesh.
    num_verts = p_gmesh%p_part%f_part%npoin

    ! Create array of elements with room for ghost elements
    allocate( par_elem_topology :: p_trian%mig_elems(num_elems + num_ghosts), stat=istat)
    check(istat==0)

    select type( this => p_trian%mig_elems )
    type is(par_elem_topology)
       p_trian%elems => this
    end select

    p_trian%elems(:)%border = -1 
    do ielem=1, p_gmesh%p_part%f_part%nebou
       p_trian%elems(p_gmesh%p_part%f_part%lebou(ielem))%border = ielem
    end do

    ! Fill array of elements (local ones)
    do ielem=1, num_elems
       p_trian%elems(ielem)%globalID    = p_gmesh%p_part%f_part%l2ge(ielem)
       p_trian%elems(ielem)%num_objects = p_trian%f_trian%elems(ielem)%num_objects
       call memalloc( p_trian%elems(ielem)%num_objects, p_trian%elems(ielem)%objects_GIDs, __FILE__, __LINE__ )
       do iobj=1, p_trian%elems(ielem)%num_objects
          jobj = p_trian%f_trian%elems(ielem)%objects(iobj)
          if ( jobj <= num_verts ) then ! It is a corner => re-use global ID
             p_trian%elems(ielem)%objects_GIDs(iobj) = p_gmesh%p_part%f_part%l2gn(jobj)
          else ! It is an edge or face => generate new local-global ID (non-consistent, non-consecutive)
             ! The ISHFT(1,50) is used to start numbering efs after vertices, assuming nvert < 2**60
             !p_trian%elems(ielem)%objects_GIDs(iobj) = ISHFT(int(p_gmesh%p_part%f_part%ipart,igp),int(32,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(60,igp))
             p_trian%elems(ielem)%objects_GIDs(iobj) = ISHFT(int(p_gmesh%p_part%f_part%ipart,igp),int(6,igp)) + int(jobj, igp) + ISHFT(int(1,igp),int(6,igp))
          end if
       end do
    end do

    ! Get objects_GIDs from ghost elements
    call ghost_elements_exchange ( p_gmesh%p_part%p_context%icontxt, p_gmesh%p_part%f_el_import, p_trian%elems )

    ! Allocate elem_topology in triangulation for ghost elements  (SBmod)
    do ielem = num_elems+1, num_elems+num_ghosts       
       p_trian%f_trian%elems(ielem)%num_objects = p_trian%elems(ielem)%num_objects
       call memalloc(p_trian%f_trian%elems(ielem)%num_objects, p_trian%f_trian%elems(ielem)%objects, &
            & __FILE__, __LINE__)
       p_trian%f_trian%elems(ielem)%objects = -1
    end do

    ! Put the topology info in the ghost elements
    do ielem= num_elems+1, num_elems+num_ghosts
       call put_topology_element_triangulation ( ielem, p_trian%f_trian )
    end do

    ! Hash table global to local for ghost elements
    call hash%init(num_ghosts)
    do ielem=num_elems+1, num_elems+num_ghosts

       call hash%put( key = p_trian%elems(ielem)%globalID, val = ielem, stat=istat)
    end do

    do ielem = 1, p_gmesh%p_part%f_part%nebou     ! Loop interface elements 
       ! Step 1: Put LID of vertices in the ghost elements (p_part%f_part)
       ilele = p_gmesh%p_part%f_part%lebou(ielem) ! local ID element
       ! aux : array of ilele (LID) vertices in GID
       nvert  = p_trian%f_trian%elems(ilele)%topology%nobje_dim(2)-1
       call memalloc( nvert, aux_igp, __FILE__, __LINE__  )
       do iobj = 1, nvert                        ! vertices only
          aux_igp(iobj) = p_trian%elems(ilele)%objects_GIDs(iobj) ! extract GIDs vertices
       end do
       do jelem = p_gmesh%p_part%f_part%pextn(ielem), & 
            p_gmesh%p_part%f_part%pextn(ielem+1)-1  ! external neighbor elements
          call hash%get(key = p_gmesh%p_part%f_part%lextn(jelem), val=jlele, stat=istat) ! LID external element
          do jobj = 1, p_trian%f_trian%elems(jlele)%topology%nobje_dim(2)-1 ! vertices external 
             if ( p_trian%f_trian%elems(jlele)%objects(jobj) == -1) then
                do iobj = 1, nvert
                   if ( aux_igp(iobj) == p_trian%elems(jlele)%objects_GIDs(jobj) ) then
                      ! Put LID of vertices for ghost_elements
                      p_trian%f_trian%elems(jlele)%objects(jobj) =  p_trian%f_trian%elems(ilele)%objects(iobj)
                   end if
                end do
             end if
          end do
       end do
       call memfree(aux_igp, __FILE__, __LINE__) 
    end do

    do ielem = 1, p_gmesh%p_part%f_part%nebou     ! Loop interface elements 
       ! Step 2: Put LID of efs in the ghost elements (p_part%f_part) 
       ilele = p_gmesh%p_part%f_part%lebou(ielem) ! local ID element
       do jelem = p_gmesh%p_part%f_part%pextn(ielem), &
            p_gmesh%p_part%f_part%pextn(ielem+1)-1  ! external neighbor elements
          call hash%get(key = p_gmesh%p_part%f_part%lextn(jelem), val=jlele, stat=istat) ! LID external element
          ! loop over all efs of external elements
          do idime =2,p_trian%f_trian%num_dims
             do iobj = p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime), &
                  p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime+1)-1 
                if ( p_trian%f_trian%elems(jlele)%objects(iobj) == -1) then ! efs not assigned yet
                   count = 1
                   ! loop over vertices of every ef
                   do jobj = p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj), &
                        p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj+1)-1    
                      ivere = p_trian%f_trian%elems(jlele)%topology%crxob_j(jobj)
                      if (p_trian%f_trian%elems(jlele)%objects(ivere) == -1) then
                         count = 0 ! not an object of the local triangulation
                         exit
                      end if
                   end do
                   if (count == 1) then
                      nvert = p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj+1)- &
                           p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj)
                      call memalloc( nvert, aux, __FILE__, __LINE__)
                      count = 1
                      do jobj = p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj), &
                           p_trian%f_trian%elems(jlele)%topology%crxob_i(iobj+1)-1 
                         ivere = p_trian%f_trian%elems(jlele)%topology%crxob_j(jobj)
                         aux(count) = p_trian%f_trian%elems(jlele)%objects(ivere)
                         count = count+1
                      end do
                      call local_id_from_vertices( p_trian%f_trian%elems(ilele), idime, aux, nvert, &
                           p_trian%f_trian%elems(jlele)%objects(iobj) )
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

    ! Step 3: Make GID consistent among processors (p_part%elems%objects_GIDs)
    do iobj = 1, p_trian%f_trian%num_objects 
       if ( (p_trian%objects(iobj)%border .ne. -1) .and. &
            (p_trian%f_trian%objects(iobj)%dimension > 1) ) then
          idime = p_trian%f_trian%objects(iobj)%dimension
          iobjg = -1
          do jelem = 1,p_trian%f_trian%objects(iobj)%num_elems_around
             jlele = p_trian%f_trian%objects(iobj)%elems_around(jelem)
             do jobj = p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime), &
                  & p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime+1)-1 ! efs of neighbor els
                if ( p_trian%f_trian%elems(jlele)%objects(jobj) == iobj ) then
                   if ( iobjg == -1 ) then 
                      iobjg  = p_trian%elems(jlele)%objects_GIDs(jobj)
                   else
                      iobjg  = min(iobjg,p_trian%elems(jlele)%objects_GIDs(jobj))
                      exit
                   end if
                end if
             end do
          end do
          do jelem = 1,p_trian%f_trian%objects(iobj)%num_elems_around
             jlele = p_trian%f_trian%objects(iobj)%elems_around(jelem)
             do jobj = p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime), &
                  & p_trian%f_trian%elems(jlele)%topology%nobje_dim(idime+1)-1 ! efs of neighbor els
                if ( p_trian%f_trian%elems(jlele)%objects(jobj) == iobj) then
                   p_trian%elems(jlele)%objects_GIDs(jobj) = iobjg
                   exit
                end if
             end do
          end do
          p_trian%objects(iobj)%globalID = iobjg          
       end if
    end do


    ! Check results
    do ielem = 1,num_elems+num_ghosts
       write (*,*) '****ielem:',ielem          
       write (*,*) '****LID_objects ghost:',p_trian%f_trian%elems(ielem)%objects
       write (*,*) '****GID_objects ghost:',p_trian%elems(ielem)%objects_GIDs
    end do
    pause

    p_trian%p_part => p_gmesh%p_part 
  end subroutine par_mesh_to_triangulation

end module par_mesh_triangulation
