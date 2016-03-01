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
module generate_vefs_mesh_conditions_names
  use types_names
  use list_types_names
  use memor_names
  use mesh_names
  use conditions_names
  use reference_fe_names
  
  implicit none
# include "debug.i90"
  private

  ! Functions
  public ::  generate_vefs_mesh_conditions

contains

  !==================================================================================================
  subroutine generate_vefs_mesh_conditions(gmsh,omsh,reference_fe,gcnd,ocnd)
    implicit none
    !------------------------------------------------------------------------------------------------
    !
    ! This routine generates the topology (nodes, edges, faces) from a mesh of linear elements.
    !
    !------------------------------------------------------------------------------------------------
    ! Parameters
    type(mesh_t)                , intent(in)    :: gmsh
    type(mesh_t)                , intent(out)   :: omsh
    class(reference_fe_t)       , intent(in)    :: reference_fe
    type(conditions_t), optional, intent(in)    :: gcnd
    type(conditions_t), optional, intent(out)   :: ocnd

    ! Local variables
    integer(ip)              :: etype,nodim(3),nndim(3)
    integer(ip)              :: i,j,k,r,s, t
    integer(ip)              :: nvef,gp,op,nd_i(2),nd_j(2)
    integer(ip)              :: iedge,iedgb, iface, ifacb,ielem, jelem
    integer(ip)              :: counter, already_counted, npoin_aux, icode
    integer(ip), allocatable :: nd_jf(:), fnode(:)
    integer(ip), allocatable :: nelpo_aux(:), lelpo_aux(:), aux1(:)
    integer(ip), allocatable :: edgeint(:,:), faceint(:,:)
    logical                  :: kfl_bc
    
    integer(ip)              :: vef_lid
    type(list_t), pointer    :: vertices_vef

    ! Create ocnd flag
    kfl_bc = (present(gcnd) .and. present(ocnd))
    ! Set scalar parameters
    omsh%nelem     = gmsh%nelem     ! Number of elements
    omsh%ndime     = gmsh%ndime     ! Number of space dimensions

    vertices_vef => reference_fe%get_vertices_vef()
    ! Construct the array of #objects(nodim) and #nodesxobject(nndim) for each dimension
    nodim = 0
    do i = 1,gmsh%ndime
       nodim(i) = reference_fe%get_number_vefs_of_dimension(i-1)
       vef_lid  = reference_fe%get_first_vef_id_of_dimension(i-1) 
       nndim(i) = reference_fe%get_number_vertices_vef(vef_lid)
    end do
    nvef=gmsh%nnode+nodim(2)+nodim(3) ! Total number of objects per element
   
    ! Allocation of auxiliar arrays
    call memalloc(gmsh%npoin+1,          nelpo_aux, __FILE__, __LINE__ )
    call memalloc(  gmsh%npoin,               aux1, __FILE__, __LINE__ )
    call memalloc(gmsh%nelem*gmsh%nnode, lelpo_aux, __FILE__, __LINE__ )
    call memalloc(gmsh%nelem,  nodim(2),   edgeint,__FILE__, __LINE__ )
    if(gmsh%ndime==3) call memalloc(gmsh%nelem, nodim(3), faceint, __FILE__, __LINE__ )

    ! Allocate omsh vectors
    call memalloc(       gmsh%nelem+1,  omsh%pnods, __FILE__, __LINE__ )
    call memalloc(gmsh%nelem*nvef, omsh%lnods , __FILE__, __LINE__ )

    ! Initialization
    omsh%npoin = gmsh%npoin ! Number of objects (nodes + edges + faces)
    nelpo_aux  = 0          ! Number of elements where every node is in (nelpo_aux(i+1)-nelpo_aux(i))
    omsh%pnods = 0          ! Array of pointers to omsh%lnods
    omsh%lnods = 0          ! Array of topologic connectivities
    lelpo_aux  = 0
    aux1       = 0
    edgeint    = 0

    ! Compute nelpo and omsh%pnods
    do ielem = 1,gmsh%nelem
       gp = gmsh%nnode*(ielem-1)
       op = nvef*(ielem-1)
       omsh%pnods(ielem) = op+1

       ! nelpo(i) computes how many elements surround i
       do i = 1,gmsh%nnode
          nelpo_aux(gmsh%lnods(gp+i)) = nelpo_aux(gmsh%lnods(gp+i)) + 1
       end do

       ! Fill the positions of omsh%lnods corresponding to the corners (objects of dimension 0
       do i = 1,gmsh%nnode
          omsh%lnods(op+i) = gmsh%lnods(gp+i)
       end do
    end do
    omsh%pnods(gmsh%nelem+1) = nvef*gmsh%nelem +1
    nelpo_aux(gmsh%npoin+1)  = gmsh%nelem*gmsh%nnode+1

    ! nelpo_aux will be the pointer to lelpo_aux, the array of elems around corners
    do i=gmsh%npoin,1,-1
       nelpo_aux(i) = nelpo_aux(i+1) - nelpo_aux(i)
    end do

    ! Fill lelpo_aux
    do ielem=1,gmsh%nelem
       gp = gmsh%nnode*(ielem-1)
       do i=1,gmsh%nnode
          j = gmsh%lnods(gp+i)
          lelpo_aux(nelpo_aux(j)+aux1(j)) = ielem
          aux1(j) = aux1(j) + 1
       end do
    end do
    call memfree(aux1,__FILE__,__LINE__)

    ! -------------------------------------------- EDGES -------------------------------------------- 
    ! First:: generate the common edges
    vef_lid = reference_fe%get_first_vef_id_of_dimension(1)
    iedge = 0
    do ielem=1,gmsh%nelem
       do i=1,nodim(2)
          if(edgeint(ielem,i)==0) then
             already_counted=0
             gp = gmsh%nnode*(ielem-1)
             nd_i(1) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)))
             nd_i(2) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)+1))

             ! Loop over the elements around corner j
             do t = nelpo_aux(nd_i(1)),nelpo_aux(nd_i(1)+1)-1
                jelem = lelpo_aux(t)

                if (jelem.gt.ielem) then
                   gp = gmsh%nnode*(jelem-1)

                   do k=1,nodim(2)
                      nd_j(1) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+k-1)))
                      nd_j(2) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+k-1)+1))

                      counter=0
                      do r=1,nndim(2)
                         do s=1,nndim(2)
                            if (nd_i(s) == nd_j(r)) then
                               counter = counter + 1
                            end if
                         end do
                      end do

                      if (counter == nndim(2)) then
                         if(already_counted==0) iedge = iedge + 1
                         edgeint(ielem,i) = iedge
                         omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)=omsh%npoin+iedge
                         edgeint(jelem,k) = iedge
                         omsh%lnods(nvef*(jelem-1)+gmsh%nnode+k)=omsh%npoin+iedge
                         already_counted=1
                      end if
                   end do
                end if
             end do
          end if
       end do
    end do
    omsh%npoin = omsh%npoin + iedge

    !Count Edges Alone
    iedgb = 0
    do i=1,gmsh%nelem
       do j=gmsh%nnode+1,gmsh%nnode+nodim(2)
          if(omsh%lnods((i-1)*nvef+j)==0) then
             iedgb = iedgb + 1
          end if
       end do
    end do
    omsh%npoin = omsh%npoin + iedgb

    ! ------------------------------------------ FACES ----------------------------------------------
    if(gmsh%ndime==3) then
       ! First:: generate the common faces
       vef_lid = reference_fe%get_first_vef_id_of_dimension(1)
       faceint = 0
       iface = 0
       call memalloc(nndim(3), nd_jf, __FILE__, __LINE__ )
       call memalloc(nndim(3), fnode, __FILE__, __LINE__ )
       do ielem=1,gmsh%nelem
          do i=1,nodim(3)
             if(faceint(ielem,i)==0) then
                already_counted=0
                gp = gmsh%nnode*(ielem-1)
                do j=1,nndim(3)
                   fnode(j) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)+j-1))
                end do

                do t = nelpo_aux(fnode(1)),nelpo_aux(fnode(1)+1)-1
                   jelem = lelpo_aux(t)

                   if (jelem.gt.ielem) then
                      gp = gmsh%nnode*(jelem-1)

                      do k=1,nodim(3)
                         do j=1,nndim(3)
                            nd_jf(j) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+k-1)+j-1))
                         end do

                         counter=0
                         do r=1,nndim(3)
                            do s=1,nndim(3)
                               if (fnode(s) == nd_jf(r)) then
                                  counter = counter + 1
                               end if
                            end do
                         end do

                         if (counter == nndim(3)) then
                            if(already_counted==0) iface = iface + 1
                            faceint(ielem,i) = iface
                            omsh%lnods(nvef*(ielem-1)+gmsh%nnode+nodim(2)+i)=omsh%npoin+iface
                            faceint(jelem,k) = iface
                            omsh%lnods(nvef*(jelem-1)+gmsh%nnode+nodim(2)+k)=omsh%npoin+iface
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end do
       call memfree(nd_jf,__FILE__,__LINE__)
       call memfree(fnode,__FILE__,__LINE__)
       omsh%npoin = omsh%npoin + iface

       ! Count faces alone
       ifacb = 0
       do i=1,gmsh%nelem
          do j=gmsh%nnode+nodim(2)+1,gmsh%nnode+nodim(2)+nodim(3)
             if(omsh%lnods((i-1)*nvef+j)==0) then
                ifacb = ifacb + 1
             end if
          end do
       end do
       omsh%npoin = omsh%npoin + ifacb
    end if

    !------------------------------------------------------------------------------------------------
    !---------------- BOUNDARY OBJECTS (AND CONSTRUCTION OF TOPOLOGICAL BC's ) ----------------------

    if (kfl_bc) then
       ! Initialize topological conditions
       call conditions_create(gcnd%ncode, gcnd%nvalu, omsh%npoin, ocnd)
       
       ! Assign Conditions on boundary vertices
       do i = 1, gmsh%npoin 
          ocnd%code(:,i) = gcnd%code(:,i)
          ocnd%valu(:,i) = gcnd%valu(:,i)
       end do
    end if

    ! Next:: we generate the unique edges (edges left to fill)
    vef_lid = reference_fe%get_first_vef_id_of_dimension(1)
    npoin_aux  = gmsh%npoin + iedge
    iedgb = 0
    do ielem=1,gmsh%nelem
       do i=1,nodim(2)
          if(omsh%lnods((ielem-1)*nvef+gmsh%nnode+i)==0) then
             iedgb = iedgb + 1
             omsh%lnods((ielem-1)*nvef+gmsh%nnode+i) = npoin_aux + iedgb
             if (kfl_bc) then
                gp = gmsh%nnode*(ielem-1)
                nd_i(1) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)))
                nd_i(2) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)+1))
                do icode = 1, ocnd%ncode
                   if (ocnd%code(icode,nd_i(1)) ==  ocnd%code(icode,nd_i(2))) then
                      ocnd%code(icode,omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)) =                   &
                           &          ocnd%code(icode,nd_i(1))
                   else
                      !write(*,*) "mesh_topology:: WARNING!! Arbitrary BC's on object",          &
                      !     &        omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)
                      ocnd%code(icode,omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)) =  0
                   end if
                end do
                ocnd%valu(:,omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)) =                             &
                     &                (ocnd%valu(:,nd_i(1)) + ocnd%valu(:,nd_i(2)))/2.0_rp
             end if
          end if
       end do
    end do

    if(gmsh%ndime==3) then
       vef_lid = reference_fe%get_first_vef_id_of_dimension(1)
       if (kfl_bc) then
          edgeint = 0
          ! Find edges on the boundary
          do ielem=1,gmsh%nelem
             do i=1,nodim(2)
                if(edgeint(ielem,i)==0) then
                   already_counted=1
                   gp = gmsh%nnode*(ielem-1)
                   nd_i(1)=gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)))
                   nd_i(2)=gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+i-1)+1))

                   ! Number of elements for j
                   do t = nelpo_aux(nd_i(1)), nelpo_aux(nd_i(1)+1)-1
                      jelem = lelpo_aux(t)

                      if (jelem.gt.ielem) then
                         gp = gmsh%nnode*(jelem-1)

                         do k=1,nodim(2)
                            nd_j(1)=gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+k-1)))
                            nd_j(2)=gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+k-1)+1))

                            counter=0
                            do r=1,nndim(2)
                               do s=1,nndim(2)
                                  if (nd_i(s) == nd_j(r)) then
                                     counter = counter + 1
                                  end if
                               end do
                            end do

                            if (counter == nndim(2)) then
                               edgeint(ielem,i) = 1
                               edgeint(jelem,k) = 1
                               already_counted= already_counted + 1
                            end if
                         end do
                      end if
                   end do
                   do icode = 1, ocnd%ncode
                      if (ocnd%code(icode,nd_i(1)) == 1 .and.  ocnd%code(icode,nd_i(2))==1) then
                         ocnd%code(icode,omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)) =  1
                      end if
                   end do
                   ocnd%valu(:,omsh%lnods(nvef*(ielem-1)+gmsh%nnode+i)) = &
                        & (ocnd%valu(:,nd_i(1)) + ocnd%valu(:,nd_i(2)))/2.0_rp
                end if
             end do
          end do
       end if

       ! Next, we generate the unique faces (faces left to fill)
       npoin_aux  = gmsh%npoin + iedge + iedgb + iface
       vef_lid = reference_fe%get_first_vef_id_of_dimension(2)
       ifacb = 0
       do ielem=1,gmsh%nelem
          j = 0
          do i=gmsh%nnode+nodim(2)+1,gmsh%nnode+nodim(2)+nodim(3)
             j = j + 1
             if(omsh%lnods((ielem-1)*nvef+i)==0) then
                ifacb = ifacb + 1
                omsh%lnods((ielem-1)*nvef+i) =  npoin_aux+ifacb
                if (kfl_bc) then
                   gp = gmsh%nnode*(ielem-1)
                   op = nvef*(ielem-1)
                   nd_i(1) = gmsh%lnods(gp+vertices_vef%l(vertices_vef%p(vef_lid+j-1)))
                   do icode = 1, ocnd%ncode
                      ocnd%valu(icode,omsh%lnods(op+i)) =                                           &
                           &            ocnd%valu(icode,omsh%lnods(op+i))                           &
                           &          + ocnd%valu(icode,nd_i(1))
                      s = 1
                      do k = 2,nndim(3)
                         nd_i(2) = gmsh%lnods(gp+vertices_vef%l                                      &
                              &    (vertices_vef%p(vef_lid+j-1)+k-1))
                         if (ocnd%code(icode,nd_i(1)) .ne.  ocnd%code(icode,nd_i(2))) s = 0
                         ocnd%valu(icode,omsh%lnods(op+i)) =                                        &
                              &     ocnd%valu(icode,omsh%lnods(op+i))                               &
                              &   + ocnd%valu(icode,nd_i(2))
                      end do
                      if (s == 1) then
                         ocnd%code(icode,omsh%lnods(op+i)) =  ocnd%code(icode,nd_i(1)) 
                      else
                         !write(*,*) 'mesh_topology:: WARNING!! Arbitrary BCs on object',        &
                         !     &        omsh%lnods(op+i)
                         ocnd%code(icode,omsh%lnods(op+i)) =  0 
                      end if
                      ocnd%valu(icode,omsh%lnods(op+i)) =                                           &
                           &    ocnd%valu(icode,omsh%lnods(op+i))/int(nndim(3))
                   end do
                end if
             end if
          end do
       end do
    end if

    ! Deallocate auxiliar arrays
    call memfree(nelpo_aux,__FILE__,__LINE__)
    call memfree(lelpo_aux,__FILE__,__LINE__)

    ! Finally, create the mesh & conditions fields
    omsh%nnode = gmsh%nnode + nodim(2) + nodim(3)

    ! Deallocate auxiliar arrays
    call memfree(edgeint,__FILE__,__LINE__)
    if(gmsh%ndime==3) call memfree(faceint,__FILE__,__LINE__)
  end subroutine generate_vefs_mesh_conditions

end module generate_vefs_mesh_conditions_names
