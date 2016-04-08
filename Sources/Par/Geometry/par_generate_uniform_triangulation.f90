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
module par_generate_uniform_triangulation_names
  use types_names
  use list_types_names
  use memor_names
  use triangulation_names
  use generate_uniform_triangulation_names
  use mesh_distribution_names
  use element_import_create_names
  use hash_table_names
  use reference_fe_names
  
  use par_environment_names
  use par_context_names
  use par_triangulation_names
  use par_conditions_names
  use par_element_exchange_names
  
  implicit none
# include "debug.i90"
  private

  public :: par_generate_uniform_triangulation
 
contains
  
  !==================================================================================================
  subroutine par_generate_uniform_triangulation(p_env,gdata,bdata,geo_reference_element,p_trian,p_cond,material)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates a parallel triangulation.                                         !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(par_environment_t)  , target, intent(in)  :: p_env
    type(uniform_mesh_descriptor_t)                , intent(in)  :: gdata
    type(uniform_conditions_descriptor_t)               , intent(in)  :: bdata
    class(reference_fe_t)           , intent(in)  :: geo_reference_element
    type(par_triangulation_t), target, intent(out) :: p_trian 
    type(par_conditions_t)           , intent(out) :: p_cond
    integer(ip), allocatable       , intent(out) :: material(:)
    ! Locals
    type(mesh_distribution_t)    :: mdist
    type (hash_table_igp_ip_t)   :: hash
    integer(ip)                  :: num_elems, num_ghosts, aux_val
    integer(ip)                  :: istat, ielem, iobj, jobj, state
    integer(ip)                  :: ilele, nvert, jelem, jlele, idime, count, ivere 
    integer(igp), allocatable    :: aux_igp(:)
    integer(ip) , allocatable    :: aux(:)
    type(list_t), pointer        :: vertices_vef
    type(par_context_t), pointer :: l1_context
    
    ! Assign environment
    p_trian%p_env => p_env
    p_cond%p_env  => p_env
    l1_context    => p_env%get_l1_context()

    if(p_trian%p_env%am_i_l1_task()) then
       
       ! Checks
       state = p_trian%state
       check(state == par_triangulation_not_created)

       ! Create Local triangulation and mesh_distribution
       call generate_uniform_triangulation(l1_context%get_rank()+1,gdata,bdata,geo_reference_element,p_trian%triangulation, &
            &                              p_cond%f_conditions,material,mdist)

       ! Create element_import
       call element_import_create(mdist,p_trian%element_import)
       num_elems  = p_trian%element_import%nelem
       num_ghosts = p_trian%element_import%nghost
       p_trian%num_ghosts = num_ghosts
       p_trian%num_elems  = num_elems

       ! Resize elem_array_len
       p_trian%triangulation%elem_array_len = num_elems+num_ghosts

       ! Create array of elements with room for ghost elements
       allocate( par_elem_topology_t :: p_trian%mig_elems(num_elems + num_ghosts), stat=istat)
       check(istat==0)

       select type( this => p_trian%mig_elems )
       type is(par_elem_topology_t)
          p_trian%elems => this
       end select

       p_trian%elems(:)%interface = -1 
       do ielem=1, mdist%nebou
          p_trian%elems(mdist%lebou(ielem))%interface = ielem
       end do

       p_trian%num_itfc_elems = mdist%nebou
       call memalloc( p_trian%num_itfc_elems, p_trian%lst_itfc_elems, __FILE__, __LINE__ )
       p_trian%lst_itfc_elems = mdist%lebou

       ! Fill array of elements (local ones)
       do ielem=1, num_elems
          p_trian%elems(ielem)%mypart      = l1_context%get_rank() + 1
          p_trian%elems(ielem)%globalID    = mdist%emap%l2g(ielem)
          p_trian%elems(ielem)%num_vefs = p_trian%triangulation%elems(ielem)%num_vefs
          call memalloc(p_trian%elems(ielem)%num_vefs, p_trian%elems(ielem)%vefs_GIDs, __FILE__, __LINE__ )
          do iobj=1, p_trian%elems(ielem)%num_vefs
             jobj = p_trian%triangulation%elems(ielem)%vefs(iobj)
             p_trian%elems(ielem)%vefs_GIDs(iobj) = mdist%nmap%l2g(jobj)
          end do
       end do

       ! Get vefs_GIDs from ghost elements
       call ghost_elements_exchange ( l1_context%get_icontxt(), p_trian%element_import, p_trian%elems )

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

       do ielem = 1, mdist%nebou     ! Loop interface elements 
          ! Step 1: Put LID of vertices in the ghost elements (f_mesh_dist)
          ilele = mdist%lebou(ielem) ! local ID element
          ! aux : array of ilele (LID) vertices in GID
          nvert  = p_trian%triangulation%elems(ilele)%reference_fe_geo%get_number_vertices()
          call memalloc( nvert, aux_igp, __FILE__, __LINE__  )
          do iobj = 1, nvert                        ! vertices only
             aux_igp(iobj) = p_trian%elems(ilele)%vefs_GIDs(iobj) ! extract GIDs vertices
          end do
          do jelem = mdist%pextn(ielem), mdist%pextn(ielem+1)-1  ! external neighbor elements
             call hash%get(key = mdist%lextn(jelem), val=jlele, stat=istat) ! LID external element
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

       do ielem = 1, mdist%nebou     ! Loop interface elements 
          ! Step 2: Put LID of efs in the ghost elements (f_mesh_dist) 
          ilele = mdist%lebou(ielem) ! local ID element
          do jelem = mdist%pextn(ielem), mdist%pextn(ielem+1)-1  ! external neighbor elements
             call hash%get(key = mdist%lextn(jelem), val=jlele, stat=istat) ! LID external element
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

       ! Deallocate
       call mesh_distribution_free(mdist)

       p_trian%state = par_triangulation_filled

    else
       call memalloc(0,material,__FILE__,__LINE__)
       ! AFM: TODO: Partially broadcast p_trian%f_trian from 1st level tasks to 2nd level tasks (e.g., num_dims)
    end if
    
  end subroutine par_generate_uniform_triangulation

end module par_generate_uniform_triangulation_names
