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
module fem_mesh_refine_names
use types_names
use memor_names
  use fem_mesh_names
  use fem_conditions_names
  use fem_materials_names
# include "debug.i90"
  implicit none
  private

  ! Functions
  !public :: fem_mesh_refine_tetrahedra, fem_mesh_lin2quad_tetrahedra, fem_mesh_hexa9
  public :: fem_mesh_refine_tetrahedra

contains
  !================================================================================================
  !
  ! Public methods:
  ! - fem_mesh_refine_tetrahedra
  !
  ! Auxiliar routines
  ! - fem_mesh_lin2quad_tetrahedra
  !
  !================================================================================================
  subroutine fem_mesh_refine_tetrahedra(msh,nodes,mat)
    !------------------------------------------------------------------------
    !
    ! This routine generates a uniform refinement of a tetrahedral mesh.
    !
    !------------------------------------------------------------------------
    implicit none
    type(fem_mesh),               intent(inout) :: msh
    type(fem_conditions)        , intent(inout) :: nodes
    type(fem_materials),optional, intent(inout) :: mat
    integer(ip), allocatable                    :: lnods_aux(:),tetr_perm(:,:),list_aux(:)
    integer(ip)                                 :: ielem,inode,ispos,ispos_old,isele,itetr,nnode,ntetr

    assert(msh%ndime==3)

    if(msh%nnode==4) then
       call fem_mesh_lin2quad_tetrahedra(msh,nodes)
    else if(msh%nnode==8) then
       call fem_mesh_hexa9(msh,nodes)
    end if

    assert(msh%nnode==10.or.msh%nnode==9)

    ! Local numeration of 8 or 12 lineal sub-tetrahedra
    nnode = 4
    if(msh%nnode==10) then
       ntetr = 8
    else if(msh%nnode==9) then
       ntetr =12
    end if
    call memalloc (nnode,ntetr,tetr_perm, __FILE__,__LINE__)
    if(msh%nnode==10) then
       tetr_perm = reshape((/1,5,7,8,5,2,6,9,7,6,3,10,8,9,10,4,5,8,9,7,5,9,6,7,10,7,6,9,10,8,7,9/),(/nnode,ntetr/))
    else if(msh%nnode==9) then
       tetr_perm = reshape((/3,4,2,9,1,2,4,9,7,6,8,9,5,8,6,9,4,8,1,9,5,1,8,9,3,2,7,9,6,7,2,9,3,7,4,9, &
            &                8,4,7,9,2,1,6,9,5,6,1,9/),(/nnode,ntetr/))
    end if

    ! Copy original lnods to lnods_aux
    call memalloc (msh%nnode*msh%nelem,lnods_aux, __FILE__,__LINE__)
    lnods_aux = msh%lnods

    ! Copy original mat%list to list_aux
    if(present(mat)) then
       call memalloc (msh%nelem,list_aux, __FILE__,__LINE__)
       list_aux = mat%list
    end if

    ! Now, we generate the new (8 or 12)*nelem tetrahedra
    call memfree (msh%lnods,__FILE__,__LINE__)
    call memalloc (nnode*msh%nelem*ntetr,msh%lnods, __FILE__,__LINE__)
    if(present(mat)) then
       call memfree (mat%list,__FILE__,__LINE__)
       call memalloc (msh%nelem*ntetr,mat%list, __FILE__,__LINE__)
    end if

    ispos = 1
    do ielem=1,msh%nelem
       ispos_old = (ielem-1)*msh%nnode
       do itetr=1,ntetr
          do inode=1,nnode
             msh%lnods(ispos) = lnods_aux(ispos_old + tetr_perm(inode,itetr))
             ispos = ispos + 1
          end do
       end do
    end do
    if(present(mat)) then
       isele = 1
       do ielem=1,msh%nelem
          do itetr=1,ntetr
             mat%list(isele) = list_aux(ielem)
             isele = isele + 1
          end do
       end do
    end if
    msh%nelem = ntetr*msh%nelem
    msh%nnode = nnode
    if(present(mat)) mat%nelem = msh%nelem

    ! Deallocate auxiliar arrays
    call memfree (lnods_aux,__FILE__,__LINE__)
    call memfree (tetr_perm,__FILE__,__LINE__)
    if(present(mat)) call memfree (list_aux,__FILE__,__LINE__)

    ! Write data from refined mesh
    write(*,*) 'Refined lineal tetrahedra mesh generated'
    write(*,*) 'nnode:',msh%nnode
    write(*,*) 'npoin:',msh%npoin
    write(*,*) 'nelem:',msh%nelem

  end subroutine fem_mesh_refine_tetrahedra

  !=============================================================================
  subroutine fem_mesh_lin2quad_tetrahedra(msh,nodes)
    !------------------------------------------------------------------------
    !
    ! This routine generates a quadratic tetrahedral mesh from a lineal one
    !
    !------------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(inout)       :: msh
    type(fem_conditions), intent(inout) :: nodes
    integer(ip)                         :: nedge,nodfac,nnode_new
    integer(ip)                         :: iedge,ielem,jelem,inode,jnode,knode,ispos,ispos_new,ipoin
    integer(ip) , allocatable           :: nelpo_aux(:),aux1(:),lelpo_aux(:),lnods_aux(:)
    integer(ip) , allocatable           :: edgeint(:,:),perm(:,:),code_aux(:,:),edge_elem(:)
    real(rp) , allocatable              :: coord_aux(:,:),valu_aux(:,:)
    integer(ip)                         :: iiaux1,iiaux2,iiaux3,jnods,jpoin,nod_i(2),nod_j(2),idime
    integer(ip)                         :: counter,rnode,snode,already_counted,iedgb,ncode,nvalu
    real(rp)                            :: val

    assert(msh%ndime==3)
    assert(msh%nnode==4)  ! Lineal tetrahedra

    nedge=6               ! Edges in tetrahedra
    nodfac=2              ! Nodes in edge
    nnode_new=10          ! Nodes in new elements
    ncode = nodes%ncode
    nvalu = nodes%nvalu

    ! Local numeration of 6 edges
    call memalloc (nodfac, nedge, perm, __FILE__,__LINE__)
    perm = reshape((/1,2,2,3,1,3,1,4,2,4,3,4/),(/nodfac,nedge/))

    ! Allocate arrays that store mesh info
    call memalloc(msh%npoin+1, nelpo_aux, __FILE__,__LINE__)
    call memalloc(msh%npoin, aux1, __FILE__,__LINE__)
    call memalloc(msh%nelem*msh%nnode, lelpo_aux, __FILE__,__LINE__)
    call memalloc(msh%nelem*nnode_new, lnods_aux, __FILE__,__LINE__)
    call memalloc(msh%nelem, nedge, edgeint, __FILE__,__LINE__)
    call memalloc(msh%nelem*nedge, edge_elem, __FILE__,__LINE__)

    ! Initialization
    nelpo_aux = 0      ! Number of elements where every node is in (nelpo_aux(i+1)-nelpo_aux(i))
    lnods_aux = 0      ! Array of connectivities of new quadratic elements
    lelpo_aux = 0
    aux1 = 0
    edgeint = 0
    edge_elem = 0

    ! Auxiliar arrays
    do ielem = 1,msh%nelem
       ispos = msh%nnode*(ielem-1)
       ispos_new = nnode_new*(ielem-1)
       do inode = 1,msh%nnode
          nelpo_aux(msh%lnods(ispos+inode)) = nelpo_aux(msh%lnods(ispos+inode)) + 1
       end do
       do inode = 1,msh%nnode
          lnods_aux(ispos_new+inode) = msh%lnods(ispos+inode)
       end do
    end do
    nelpo_aux(msh%npoin+1) = msh%nelem*msh%nnode+1

    do inode=msh%npoin,1,-1
       nelpo_aux(inode) = nelpo_aux(inode+1) - nelpo_aux(inode)
    end do

    do ielem=1,msh%nelem
       ispos = msh%nnode*(ielem-1)
       do inode=1,msh%nnode
          ipoin = msh%lnods(ispos+inode)
          lelpo_aux(nelpo_aux(ipoin)+aux1(ipoin)) = ielem
          aux1(ipoin) = aux1(ipoin) + 1
       end do
    end do
    call memfree(aux1,__FILE__,__LINE__)

    ! First, we generate the common interior edges
    iedge = 0
    do ielem=1,msh%nelem
       do inode=1,nedge
          if(edgeint(ielem,inode)==0) then
             already_counted=0
             ispos = msh%nnode*(ielem-1)
             ipoin = msh%lnods(ispos+perm(1,inode))
             jpoin = msh%lnods(ispos+perm(2,inode))
             ! Number of elements for ipoin
             iiaux1 = nelpo_aux(ipoin)
             iiaux2 = nelpo_aux(ipoin+1)-1

             do iiaux3 = iiaux1,iiaux2
                jelem = lelpo_aux(iiaux3)

                if (jelem.gt.ielem) then
                   ispos = msh%nnode*(jelem-1)

                   nod_i(1) = ipoin
                   nod_i(2) = jpoin

                   do knode=1,nedge
                      nod_j(1) = msh%lnods(ispos+perm(1,knode))
                      nod_j(2) = msh%lnods(ispos+perm(2,knode))

                      counter=0
                      do rnode=1,nodfac
                         do snode=1,nodfac
                            if (nod_i(snode) == nod_j(rnode)) then
                               counter = counter + 1
                            end if
                         end do
                      end do

                      if (counter == nodfac) then
                         if(already_counted==0) iedge = iedge + 1
                         edgeint(ielem,inode) = iedge
                         lnods_aux(nnode_new*(ielem-1)+msh%nnode+inode)=msh%npoin+iedge
                         edgeint(jelem,knode) = iedge
                         lnods_aux(nnode_new*(jelem-1)+msh%nnode+knode)=msh%npoin+iedge
                         edge_elem(iedge)=ielem
                         already_counted=1
                      end if
                   end do
                end if
             end do
          end if
       end do
    end do
    call memfree(nelpo_aux,__FILE__,__LINE__)
    call memfree(lelpo_aux,__FILE__,__LINE__)

    ! Next, we generate the edges on the boundary (nodes left to fill)
    iedgb = 0
    do inode=1,nnode_new*msh%nelem
       if(lnods_aux(inode)==0) then
          iedgb = iedgb + 1
          lnods_aux(inode) = msh%npoin+iedge+iedgb
          !edge_elem(iedge+iedgb) = inode/nnode_new
          edge_elem(iedge+iedgb) = CEILING(REAL(inode)/REAL(nnode_new))-1
       end if
    end do

    ! Allocate array for new points & conditions
    call memalloc(msh%ndime, msh%npoin+iedge+iedgb, coord_aux, __FILE__,__LINE__)
    call memalloc(nodes%ncode, msh%npoin+iedge+iedgb, code_aux, __FILE__,__LINE__)
    code_aux=0
    call memalloc(nodes%nvalu, msh%npoin+iedge+iedgb, valu_aux, __FILE__,__LINE__)
    valu_aux=0.0_rp

    ! First, we copy the original nodes & conditions
    do inode=1,msh%npoin
       do idime=1,msh%ndime
          coord_aux(idime,inode) = msh%coord(idime,inode)
       end do
       do jnode=1,ncode
          code_aux(jnode,inode) = nodes%code(jnode,inode)
       end do
       do jnode=1,nvalu
          valu_aux(jnode,inode) = nodes%valu(jnode,inode)
       end do
    end do

    ! Second, we add the new nodes on the common edges & bc's
    do inode=1,iedge
       ielem = edge_elem(inode)
       do jelem=1,nedge
          if (edgeint(ielem,jelem)==inode) then
             ispos = msh%nnode*(ielem-1)
             nod_i(1) = msh%lnods(ispos+perm(1,jelem))
             nod_i(2) = msh%lnods(ispos+perm(2,jelem))

             do idime=1,msh%ndime
                val = 0.0_rp
                do jnode=1,nodfac
                   val = val + msh%coord(idime,nod_i(jnode))
                end do
                coord_aux(idime,msh%npoin+inode) = val/real(nodfac)
             end do
             do jnode=1,ncode
                if(code_aux(jnode,nod_i(1))==code_aux(jnode,nod_i(2))) then
                   code_aux(jnode,msh%npoin+inode) = code_aux(jnode,nod_i(1))
                end if
             end do
             do jnode=1,nvalu
                valu_aux(jnode,msh%npoin+inode) = 0.5_rp*(valu_aux(jnode,nod_i(1))+valu_aux(jnode,nod_i(2)))
             end do
             exit
          end if
       end do
    end do
    call memfree(edgeint,__FILE__,__LINE__)

    ! Third, we add the nodes on the boundary edges & create new bc's
    do inode=1,iedgb
       ielem = edge_elem(iedge+inode)
       ispos = (ielem)*msh%nnode
       ispos_new = (ielem)*nnode_new
       do ipoin=1,nedge
          if(lnods_aux(ispos_new+msh%nnode+ipoin)==msh%npoin+iedge+inode) then
             nod_i(1) = msh%lnods(ispos+perm(1,ipoin))
             nod_i(2) = msh%lnods(ispos+perm(2,ipoin))
             do idime=1,msh%ndime
                val = 0.0_rp
                do jnode=1,nodfac
                   val = val + msh%coord(idime,nod_i(jnode))
                end do
                coord_aux(idime,msh%npoin+iedge+inode) = val/real(nodfac)
             end do
             do jnode=1,ncode
                if(code_aux(jnode,nod_i(1))==code_aux(jnode,nod_i(2))) then
                   code_aux(jnode,msh%npoin+iedge+inode) = code_aux(jnode,nod_i(1))
                end if
             end do
             do jnode=1,nvalu
                valu_aux(jnode,msh%npoin+iedge+inode) = 0.5_rp*(valu_aux(jnode,nod_i(1))+valu_aux(jnode,nod_i(2)))
             end do
             exit
          end if
       end do
    end do
    call memfree(perm,__FILE__,__LINE__)
    call memfree(edge_elem,__FILE__,__LINE__)

    ! Finally, create the mesh & conditions fields
    call memfree(msh%lnods,__FILE__,__LINE__)
    call memfree(msh%coord,__FILE__,__LINE__)
    call memfree(nodes%code,__FILE__,__LINE__)
    call memfree(nodes%valu,__FILE__,__LINE__)
    msh%npoin = msh%npoin + iedge + iedgb
    msh%nnode = nnode_new
    call memalloc(msh%nnode*msh%nelem, msh%lnods, __FILE__,__LINE__)
    call memalloc(msh%ndime, msh%npoin, msh%coord, __FILE__,__LINE__)
    call memalloc(ncode, msh%npoin, nodes%code, __FILE__,__LINE__)
    call memalloc(nvalu, msh%npoin, nodes%valu, __FILE__,__LINE__)
    msh%lnods = lnods_aux
    msh%coord = coord_aux
    nodes%code = code_aux
    nodes%valu = valu_aux
    nodes%ncond = msh%npoin
    call memfree(lnods_aux,__FILE__,__LINE__)
    call memfree(coord_aux,__FILE__,__LINE__)
    call memfree(code_aux,__FILE__,__LINE__)
    call memfree(valu_aux,__FILE__,__LINE__)

    ! Write data from quadratic mesh
    write(*,*) 'Quadratic tetrahedra mesh generated'
    write(*,*) 'nnode:',msh%nnode
    write(*,*) 'npoin:',msh%npoin
    write(*,*) 'nelem:',msh%nelem

  end subroutine fem_mesh_lin2quad_tetrahedra

  !=============================================================================
  subroutine fem_mesh_hexa9(msh,nodes)
    !------------------------------------------------------------------------
    !
    ! This routine generates an hexahedra with an interior node
    !
    !------------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(inout)       :: msh
    type(fem_conditions), intent(inout) :: nodes

    integer(ip)                         :: nnode_new,ielem,inode,jnode,idime,ipoin,inipo,endpo,ncode,nvalu
    integer(ip) , allocatable           :: lnods_aux(:),code_aux(:,:)
    real(rp)    , allocatable           :: coord_aux(:,:),valu_aux(:,:)
    real(rp)                            :: val

    assert(msh%ndime==3)
    assert(msh%nnode==8)  ! Lineal hexahedra

    nnode_new=9          ! Nodes in new elements
    ncode = nodes%ncode
    nvalu = nodes%nvalu

    ! Allocate arrays that store mesh and conditions info
    call memalloc(msh%nelem*nnode_new, lnods_aux, __FILE__,__LINE__)
    call memalloc(msh%ndime, msh%npoin+msh%nelem, coord_aux, __FILE__,__LINE__)
    call memalloc(nodes%ncode, msh%npoin+msh%nelem, code_aux, __FILE__,__LINE__)
    call memalloc(nodes%nvalu, msh%npoin+msh%nelem, valu_aux, __FILE__,__LINE__)

    ! Initialization
    lnods_aux = 0      ! Array of connectivities of new elements
    code_aux=0
    valu_aux=0.0_rp

    ! First, we copy the original nodes & conditions
    do inode=1,msh%npoin
       do idime=1,msh%ndime
          coord_aux(idime,inode) = msh%coord(idime,inode)
       end do
       do jnode=1,ncode
          code_aux(jnode,inode) = nodes%code(jnode,inode)
       end do
       do jnode=1,nvalu
          valu_aux(jnode,inode) = nodes%valu(jnode,inode)
       end do
    end do

    ! Second, we add the new node on the interior of the hexahedra
    do ielem=1,msh%nelem
       do idime=1,msh%ndime
          val = 0.0_rp
          do inode=1,msh%nnode
             ipoin = msh%lnods((ielem-1)*msh%nnode+inode)
             val = val + msh%coord(idime,ipoin)
          end do
          coord_aux(idime,msh%npoin+ielem) = val/real(msh%nnode)
       end do
       inipo = (ielem-1)*nnode_new+1
       endpo = (ielem-1)*nnode_new+msh%nnode
       lnods_aux(inipo:endpo) = msh%lnods((ielem-1)*msh%nnode+1:ielem*msh%nnode)
       lnods_aux(endpo+1) = msh%npoin+ielem
    end do
   
    ! Finally, create the mesh & conditions fields
    call memfree(msh%lnods,__FILE__,__LINE__)
    call memfree(msh%coord,__FILE__,__LINE__)
    call memfree(nodes%code,__FILE__,__LINE__)
    call memfree(nodes%valu,__FILE__,__LINE__)
    msh%npoin = msh%npoin + msh%nelem
    msh%nnode = nnode_new
    call memalloc(msh%nnode*msh%nelem, msh%lnods, __FILE__,__LINE__)
    call memalloc(msh%ndime, msh%npoin, msh%coord, __FILE__,__LINE__)
    call memalloc(ncode, msh%npoin, nodes%code, __FILE__,__LINE__)
    call memalloc(nvalu, msh%npoin, nodes%valu, __FILE__,__LINE__)
    msh%lnods = lnods_aux
    msh%coord = coord_aux
    nodes%code = code_aux
    nodes%valu = valu_aux
    nodes%ncond = msh%npoin
    call memfree(lnods_aux,__FILE__,__LINE__)
    call memfree(coord_aux,__FILE__,__LINE__)
    call memfree(code_aux,__FILE__,__LINE__)
    call memfree(valu_aux,__FILE__,__LINE__)

    ! Write data from quadratic mesh
    write(*,*) 'Interior node on hexahedra generated'
    write(*,*) 'nnode:',msh%nnode
    write(*,*) 'npoin:',msh%npoin
    write(*,*) 'nelem:',msh%nelem

  end subroutine fem_mesh_hexa9

end module fem_mesh_refine_names

