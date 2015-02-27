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
module fem_space_faces
  use types
  use memor
  use hash_table_names
  use fem_mesh_names
  use integration_names
!  use fem_mesh_faces
  use fem_space_names
  use fem_space_types
  use fem_partition_names
  use array_names

# include "debug.i90"
  implicit none

  private
  integer(ip), parameter :: max_subfaces = 4
  ! Functions
  public :: fem_space_faces_list_create

  interface fill_fem_faces_list
     module procedure fill_fem_faces_list_nostruct, fill_fem_faces_list_struct
  end interface fill_fem_faces_list

contains

 !===================================================================================================
  subroutine fem_space_faces_list_create( fspac, ndime,prob_nunk,prob_list_nunk, fpart )
    implicit none
    ! Parameters
    type(fem_space)      , intent(inout),target  :: fspac
    integer(ip)          , intent(in)            :: ndime
    integer(ip), optional, intent(in)            :: prob_nunk
    integer(ip), optional, intent(in)            :: prob_list_nunk(:)
    type(fem_partition), optional, intent(in)    :: fpart

    ! Local variables
    integer(ip) :: nunkn, fadof, aux, fadof_aux(2),nodfac, iobje
    integer(ip) :: iface, jelem, kelem,  istat, nnode, num_interpolations, iint
    integer(ip) :: pos_famat, pos_aux_famat(2), pos_favec, pos_faint,iunk
    integer(ip) :: v_key, gtype(2), utype(2), g_ord(2), u_ord(2),subface(2)!, nface(2)
    type(fem_fixed_info_pointer) :: gfinf(2),ufinf(2)

    ! To be copy-pasted in the future from FEMPAR

  end subroutine fem_space_faces_list_create


 !====================================================================================
  subroutine fill_fem_faces_list_nostruct( lface, tmesh )

    implicit none
    ! Parameters
    type(fem_mesh), intent(in)    :: tmesh
    type(fem_face), intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)

    ! Local variables
    integer(ip)    :: iface, ifacb, jface, kface, jelem, kelem
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh


    call mesh_to_dual(tmesh,d_msh)
    faceint = 0
    iface = 0
    ifacb = 0
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                !BOUNDARY FACE
                ifacb = ifacb +1    
                lface(tmesh%tface + ifacb)%nei_elem(1) = jelem
                lface(tmesh%tface + ifacb)%nei_elem(2) = 0
                lface(tmesh%tface + ifacb)%pos_elem(1) = jface
                lface(tmesh%tface + ifacb)%pos_elem(2) = 0
                lface(tmesh%tface + ifacb)%subface    = 0
                lface(tmesh%tface + ifacb)%refinement_level = 0 
             elseif (nelfa==2) then
                !INTERIOR FACE
                do i = 0,1
                   kelem = d_msh%lnods(d_msh%pnods(gface)+i)
                   if (kelem .ne. jelem) exit
                end do
                assert((i==1) .or. (d_msh%lnods(d_msh%pnods(gface)+1) == jelem))
                mface = tmesh%pnods(kelem+1) - tmesh%p_face(kelem)
                do kface = 1, mface
                   if (tmesh%lnods(tmesh%p_face(kelem)+kface-1) == gface) exit
                end do
                assert(kface .le. mface)
                iface = iface +1           
                faceint(kelem,kface) = iface
                lface(iface)%nei_elem(1) = jelem
                lface(iface)%nei_elem(2) = kelem
                lface(iface)%pos_elem(1) = jface
                lface(iface)%pos_elem(2) = kface  
                lface(iface)%subface     = 0
                lface(iface)%refinement_level = 0  
             else
                write(*,*) __FILE__,__LINE__,'ERROR! More than 2 elements per face'
             end if            
          end if
       end do
    end do
    call fem_mesh_free(d_msh)
    assert(iface == tmesh%tface)
    assert(ifacb == tmesh%bou_tface)

  end subroutine fill_fem_faces_list_nostruct

 !====================================================================================
  subroutine fill_fem_faces_list_struct( lface, tmesh, fpart)
    implicit none
    ! Parameters
    type(fem_mesh)     , intent(in)    :: tmesh
    type(fem_partition), intent(in)    :: fpart
    type(fem_face)     , intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)

    ! Local variables
    integer(ip)    :: iface, ifacbe, ifacbi, jface, kface, jelem, kelem, bface
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh


    call mesh_to_dual(tmesh,d_msh)
    faceint = 0
    iface  = 0
    ifacbe = 0
    ifacbi = 0
    bface = tmesh%bou_tface - tmesh%shaface
    do jelem = 1, tmesh%nelem
       nface = tmesh%pnods(jelem+1) - tmesh%p_face(jelem)
       do jface = 1, nface
          if (faceint(jelem,jface)==0) then
             gface = tmesh%lnods(tmesh%p_face(jelem)+jface-1)
             nelfa = d_msh%pnods(gface+1) - d_msh%pnods(gface)
             if (nelfa == 1) then
                if (gface <= fpart%nmap%ni) then
                   ! EXTERNAL BOUNDARY FACE
                   ifacbe = ifacbe + 1    
                   lface(tmesh%tface + ifacbe)%nei_elem(1) = jelem
                   lface(tmesh%tface + ifacbe)%nei_elem(2) = 0
                   lface(tmesh%tface + ifacbe)%pos_elem(1) = jface
                   lface(tmesh%tface + ifacbe)%pos_elem(2) = 0
                   lface(tmesh%tface + ifacbe)%subface     = 0
                   lface(tmesh%tface + ifacbe)%refinement_level = 0 
                else
                   ! INTERFACE BOUNDARY FACE
                   ifacbi = ifacbi + 1    
                   lface(tmesh%tface + bface + ifacbi)%nei_elem(1) = jelem
                   lface(tmesh%tface + bface + ifacbi)%nei_elem(2) = 0
                   lface(tmesh%tface + bface + ifacbi)%pos_elem(1) = jface
                   lface(tmesh%tface + bface + ifacbi)%pos_elem(2) = 0
                   lface(tmesh%tface + bface + ifacbi)%subface     = 0
                   lface(tmesh%tface + bface + ifacbi)%refinement_level = 0 
                end if
             elseif (nelfa==2) then
                !INTERIOR FACE
                do i = 0,1
                   kelem = d_msh%lnods(d_msh%pnods(gface)+i)
                   if (kelem .ne. jelem) exit
                end do
                assert((i==1) .or. (d_msh%lnods(d_msh%pnods(gface)+1) == jelem))
                mface = tmesh%pnods(kelem+1) - tmesh%p_face(kelem)
                do kface = 1, mface
                   if (tmesh%lnods(tmesh%p_face(kelem)+kface-1) == gface) exit
                end do
                assert(kface .le. mface)
                iface = iface +1           
                faceint(kelem,kface) = iface
                lface(iface)%nei_elem(1) = jelem
                lface(iface)%nei_elem(2) = kelem
                lface(iface)%pos_elem(1) = jface
                lface(iface)%pos_elem(2) = kface   
                lface(iface)%subface     = 0
                lface(iface)%refinement_level = 0 
             else
                write(*,*) __FILE__,__LINE__,'ERROR! More than 2 elements per face'
             end if            
          end if
       end do
    end do
    call fem_mesh_free(d_msh)

    assert(iface == tmesh%tface)
    assert((ifacbi+ifacbe) == tmesh%bou_tface)
  end subroutine fill_fem_faces_list_struct

 !====================================================================================
  subroutine fill_fem_faces_list_nostruct_adapt( lface,lelem, tmesh )
    implicit none
    ! Parameters
    type(fem_mesh)   , intent(in)    :: tmesh
    type(fem_face)   , intent(inout) :: lface(tmesh%tface+tmesh%bou_tface)
    type(fem_element), intent(in)    :: lelem(:)

    ! Local variables
    integer(ip)    :: iface, ifacb, jface, kface, jelem, kelem
    integer(ip)    :: i, nface, mface, gface, nelfa
    integer(ip)    :: faceint(tmesh%nelem,tmesh%nface)
    type(fem_mesh) :: d_msh

    integer(ip)                   :: l_obje, nnodexface, num_constraints, p_node, l_node, g_node
    integer(ip)                   :: first_node, p_elem, neigh_elem, neigh_obje, ndime, clist_id
    logical                       :: is_the_face
    type(fem_fixed_info), pointer :: elem_info,neigh_elem_info

    ! To be copy-pasted in the future from fempar

  end subroutine fill_fem_faces_list_nostruct_adapt

  !====================================================================================
  subroutine check_face (num_nodes,list_nodes, ielem, iobje, elem_info, pnods, lnods, is_the_face)
    implicit none
    ! Parameters
    integer(ip)         , intent(in)  :: num_nodes, list_nodes(num_nodes)
    integer(ip)         , intent(in)  :: ielem, iobje
    type(fem_fixed_info), intent(in)  :: elem_info
    integer(ip)         , intent(in)  :: pnods(:),lnods(:)
    logical             , intent(out) :: is_the_face

    ! Local variables
    integer(ip) :: p_first_node, p_l_node, l_node, g_node, inode

    p_first_node = pnods(ielem) -1
    
    ! Loop over all the nodes in the iobje of ielem
    ploop: do p_l_node = elem_info%crxob_i(iobje),elem_info%crxob_i(iobje+1)-1
       ! local node id
       l_node = elem_info%crxob_j(p_l_node)
       ! Global node id
       g_node = lnods(p_first_node + l_node)
       ! Loop over all the nodes of the face; check if the node is in the list
       iloop: do inode = 1, num_nodes
          ! If in the list, check next node
          if(list_nodes(inode) == g_node) cycle ploop
       end do iloop
       ! If not in the list, this is not the face we were looking for
       is_the_face = .false.
       return
    end do ploop
    is_the_face = .true.

    
  end subroutine check_face

 !====================================================================================
  ! This subroutines checkes for neigh_elem wich subfaces of neigh_obje corresponds to the obje of 
  ! elem. It is assuming that there ist at most one level of refinement and thus, they have one 
  ! common node.
  subroutine check_subface (elem, obje, neigh_elem, neigh_obje, elem_info, neigh_elem_info, pnods,   &
       &                    lnods, subface)
    implicit none
    ! Parameters
    integer(ip)         , intent(in)  :: elem, obje
    integer(ip)         , intent(in)  :: neigh_elem, neigh_obje
    type(fem_fixed_info), intent(in)  :: elem_info, neigh_elem_info
    integer(ip)         , intent(in)  :: pnods(:),lnods(:)
    integer(ip)         , intent(out) :: subface

    ! Local variables
    integer(ip) :: p_first_node, p_l_node, l_node, g_node
    integer(ip) :: p_neigh_first_node, p_l_neigh_node, l_neigh_node, g_neigh_node

    p_first_node       = pnods(elem) - 1
    p_neigh_first_node = pnods(neigh_elem) - 1

    subface = 1
    ! Loop over all the nodes in the neigh_obje of neigh_elem 
    ploop: do p_l_neigh_node = neigh_elem_info%crxob_i(neigh_obje),                                 &
         &                     neigh_elem_info%crxob_i(neigh_obje+1)-1
       ! local node id
       l_neigh_node = elem_info%crxob_j(p_l_neigh_node)
       ! Global node id
       g_neigh_node = lnods(p_neigh_first_node + l_neigh_node)
       ! Loop over all the nodes of the obje of elem; check if the node is in the list
       do p_l_node = elem_info%crxob_i(obje),elem_info%crxob_i(obje+1)-1
          ! local node id
          l_node = elem_info%crxob_j(p_l_node)
          ! Global node id
          g_node = lnods(p_first_node + l_node)
          ! Loop over all the nodes of the face; check if the node is in the list
          if(g_neigh_node == g_node) exit ploop
       end do
       subface = subface + 1
    end do ploop

    assert(subface <= neigh_elem_info%crxob_i(neigh_obje+1)- neigh_elem_info%crxob_i(neigh_obje))
    
  end subroutine check_subface

end module fem_space_faces
