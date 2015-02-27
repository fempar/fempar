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
module fem_mesh_dual
  use types
  use memor
  use fem_mesh_names
  use fem_space_names 
 
# include "debug.i90"
  implicit none
  private

  public :: mesh_to_dual_w_faces

contains
  !=============================================================================
  subroutine mesh_to_dual_w_faces(primal_mesh,dual_mesh, femsp)
    !-----------------------------------------------------------------------
    ! This routine generates the dual mesh (list of elements around
    ! points) of a given primal mesh considering the interior faces as elements.
    ! The dual mesh allways has nelty/=1
    ! and, as it is not used, we set it to 0.
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh
    type(fem_mesh) , intent(out) :: dual_mesh
    type(fem_space), intent(in)  :: femsp
    
    ! Set constants
    dual_mesh%nelty=0
    dual_mesh%ndime=primal_mesh%ndime
    dual_mesh%npoin=primal_mesh%nelem + 2*primal_mesh%tface
    dual_mesh%nelem=primal_mesh%npoin ! + primal_mesh%nelem
    dual_mesh%nelpo=0
    dual_mesh%nboun=0

    ! Count elements and faces around points and generate pointers to store them 
    ! Create and fill dual_mesh%pnods
    call memalloc (dual_mesh%nelem+1, dual_mesh%pnods, __FILE__,__LINE__)    
    call count_elements_and_faces_around_points( primal_mesh, dual_mesh,femsp)
    
    ! List elements and faces around points 
    ! Create and fill dual_mesh%lnods
    call memalloc (dual_mesh%pnods(dual_mesh%nelem+1), dual_mesh%lnods, __FILE__,__LINE__)
    call list_elements_and_faces_around_points(primal_mesh,dual_mesh,femsp)
    
  end subroutine mesh_to_dual_w_faces

 !=============================================================================
  subroutine count_elements_and_faces_around_points( primal_mesh, dual_mesh,femsp)
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh
    type(fem_mesh) , intent(inout) :: dual_mesh
    type(fem_space), intent(in)    :: femsp

    ! Local variables
    integer(ip) :: ielem, inode, ipoin, iface, one_or_two, gelem, neigh_elem
    integer(ip) :: iobje, jpoin, l_obje, jnode
    logical     :: should_be_added

    ! Initialize
    dual_mesh%pnods = 0

    ! Loop over elements to count how many elements around each object
    do ielem = 1, primal_mesh%nelem
       ! Loop over objects in ielem
       do inode = primal_mesh%pnods(ielem), primal_mesh%pnods(ielem+1)-1
          ! Global id of the object
          ipoin = primal_mesh%lnods(inode)
          ! Add another element to the counter
          dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
       end do
!!$       ! Add now the id of interior objects, not considered in the topological mesh
!!$       ! The id of the interior of ielem will be npoin + ielem
!!$       ipoin = primal_mesh%npoin + ielem
!!$       dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
    end do

    ! Loop over face to count how many of them around each object
    do iface = 1, primal_mesh%tface
       ! For each face we construct two entities
       do one_or_two = 1,2
          ! One of the objects around iface
          gelem = femsp%lface(iface)%nei_elem(one_or_two)
          ! Loop over the objects of gelem
          do inode= primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
             ! Global id of the object
             ipoin = primal_mesh%lnods(inode)
             ! Add another face to the counter
             dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
          end do
          
          ! The other object (we will only take in account its objects on the face)
          ! TODO: For 3D, we are not considering the edges in the face
          neigh_elem = femsp%lface(iface)%nei_elem(3-one_or_two)
          ! The local object corresponding to the face
          iobje = femsp%lelem(neigh_elem)%p_geo_info%nobje_dim(primal_mesh%ndime) +                 &
               &  femsp%lface(iface)%pos_elem(3-one_or_two) - 1
          ! Loop over the local corners in iobje
          do inode = femsp%lelem(neigh_elem)%p_geo_info%crxob_i(iobje),                             &
               &     femsp%lelem(neigh_elem)%p_geo_info%crxob_i(iobje+1)-1
             l_obje = femsp%lelem(neigh_elem)%p_geo_info%crxob_j(inode)
             ! global id of the object
             ipoin = primal_mesh%lnods(primal_mesh%pnods(neigh_elem) + l_obje - 1)

             ! Check the object is not in gelem already
             should_be_added = .true.
             do jnode = primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
                jpoin = primal_mesh%lnods(jnode)
                if (jpoin == ipoin) then
                   should_be_added = .false.
                   exit
                end if
             end do

             ! If it is not in gelem, add the face to the counter
             if (should_be_added) then
                dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
             end if
          end do
          
          ! We repeat the process for the face object itself
          ! Notice that for 3D the edges are still missing (ahierro 01/2015)
          ipoin = primal_mesh%lnods(primal_mesh%pnods(neigh_elem) + iobje - 1)

          ! Check the object is not in gelem already
          should_be_added = .true.
          do jnode = primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
             jpoin = primal_mesh%lnods(jnode)
             if (jpoin == ipoin) then
                should_be_added = .false.
                exit
             end if
          end do

          ! If it is not in gelem, add the face to the counter
          if (should_be_added) then
             dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
          end if

!!$          ! Add now the id of interior objects, not considered in the topological mesh
!!$          ! The id of the interior of ielem will be npoin + ielem
!!$          ipoin = primal_mesh%npoin + gelem
!!$          dual_mesh%pnods(ipoin+1) = dual_mesh%pnods(ipoin+1) + 1
       end do
    end do

    ! Find the maximum number of entities around a object
    dual_mesh%nnode=0
    do ipoin = 1,dual_mesh%nelem
       dual_mesh%nnode = max(dual_mesh%nnode,dual_mesh%pnods(ipoin+1))
    end do
    
    ! Compute pointers to the starting position list
    dual_mesh%pnods(1) = 1
    do ipoin = 1, dual_mesh%nelem
       dual_mesh%pnods(ipoin+1) = dual_mesh%pnods (ipoin+1) + dual_mesh%pnods (ipoin)
    end do
  end subroutine count_elements_and_faces_around_points

!=============================================================================
  subroutine list_elements_and_faces_around_points( primal_mesh, dual_mesh,femsp)
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh
    type(fem_mesh) , intent(inout) :: dual_mesh
    type(fem_space), intent(in)    :: femsp

    ! Local variables
    integer(ip) :: ielem, inode, ipoin, iface, one_or_two, gelem, neigh_elem
    integer(ip) :: iobje, jpoin, l_obje, jnode
    logical     :: should_be_added

    ! Initialize
    dual_mesh%lnods = 0

    ! Loop over elements to count how many elements around each object
    do ielem = 1, primal_mesh%nelem
       ! Loop over objects in ielem
       do inode = primal_mesh%pnods(ielem), primal_mesh%pnods(ielem+1)-1
          ! Global id of the object
          ipoin = primal_mesh%lnods(inode)
          ! Add another element to the list
          dual_mesh%lnods(dual_mesh%pnods(ipoin)) = ielem
          dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
       end do
!!$       ! Add now the id of interior objects, not considered in the topological mesh
!!$       ! The id of the interior of ielem will be npoin + ielem
!!$       ipoin = primal_mesh%npoin + ielem
!!$       dual_mesh%lnods(dual_mesh%pnods(ipoin)) = ielem
!!$       dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
    end do

    ! Loop over face to count how many of them around each object
    do iface = 1, primal_mesh%tface
       ! For each face we construct two entities
       do one_or_two = 1,2
          ! One of the objects around iface
          gelem = femsp%lface(iface)%nei_elem(one_or_two)
          ! Loop over the objects of gelem
          do inode= primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
             ! Global id of the object
             ipoin = primal_mesh%lnods(inode)
             ! Add another face to the list
             dual_mesh%lnods(dual_mesh%pnods(ipoin)) = primal_mesh%nelem + 2*(iface-1) + one_or_two
             dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
          end do
          
          ! The other object (we will only take in account its objects on the face)
          ! TODO: For 3D, we are not considering the edges in the face
          neigh_elem = femsp%lface(iface)%nei_elem(3-one_or_two)
          ! The local object corresponding to the face
          iobje = femsp%lelem(neigh_elem)%p_geo_info%nobje_dim(primal_mesh%ndime) +                 &
               &  femsp%lface(iface)%pos_elem(3-one_or_two) - 1
          ! Loop over the local corners in iobje
          do inode = femsp%lelem(neigh_elem)%p_geo_info%crxob_i(iobje),                             &
               &     femsp%lelem(neigh_elem)%p_geo_info%crxob_i(iobje+1)-1
             l_obje = femsp%lelem(neigh_elem)%p_geo_info%crxob_j(inode)
             ! global id of the object
             ipoin = primal_mesh%lnods(primal_mesh%pnods(neigh_elem) + l_obje - 1)

             ! Check the object is not in gelem already
             should_be_added = .true.
             do jnode = primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
                jpoin = primal_mesh%lnods(jnode)
                if (jpoin == ipoin) then
                   should_be_added = .false.
                   exit
                end if
             end do

             ! If it is not in gelem, add the face to the list
             if (should_be_added) then
                dual_mesh%lnods(dual_mesh%pnods(ipoin)) = primal_mesh%nelem + 2*(iface-1)+one_or_two
                dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
             end if
          end do
          
          ! We repeat the process for the face object itself
          ! Notice that for 3D the edges are still missing (ahierro 01/2015)
          ipoin = primal_mesh%lnods(primal_mesh%pnods(neigh_elem) + iobje - 1)

          ! Check the object is not in gelem already
          should_be_added = .true.
          do jnode = primal_mesh%pnods(gelem),primal_mesh%pnods(gelem+1)-1
             jpoin = primal_mesh%lnods(jnode)
             if (jpoin == ipoin) then
                should_be_added = .false.
                exit
             end if
          end do

          ! If it is not in gelem, add the face to the list
          if (should_be_added) then
             dual_mesh%lnods(dual_mesh%pnods(ipoin)) = primal_mesh%nelem + 2*(iface-1)+one_or_two
             dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
          end if

!!$          ! Add now the id of interior objects, not considered in the topological mesh
!!$          ! The id of the interior of ielem will be npoin + ielem
!!$          ipoin = primal_mesh%npoin + gelem
!!$          dual_mesh%lnods(dual_mesh%pnods(ipoin)) = primal_mesh%nelem + 2*(iface-1)+one_or_two
!!$          dual_mesh%pnods(ipoin) = dual_mesh%pnods(ipoin) + 1
       end do
    end do
    
    ! Recover pnods
    do ipoin =  dual_mesh%nelem+1,2,-1
       dual_mesh%pnods(ipoin) = dual_mesh%pnods (ipoin-1)
    end do
    dual_mesh%pnods(1) = 1

  end subroutine list_elements_and_faces_around_points
end module fem_mesh_dual
  
