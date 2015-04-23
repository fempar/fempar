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
module fem_mesh_gen_partition
  use types
  use memor
  use maps_names
  use mesh_graph
  use fem_mesh_names
  ! use fem_partition_names
  use fem_conditions_names
  use fem_materials_names
  use fem_mesh_gen
  use sort_names
# include "debug.i90"
  implicit none
  private

  type geom_size 
     integer(ip) :: nnode      ! Number of nodes on each element
     integer(ip) :: nedir(3)   ! Number of elements on each direction
     integer(ip) :: npdir(3)   ! Number of parts on each direction
     integer(ip) :: nsckt(3)   ! Number of parts on each socket and direction
     integer(ip) :: npsoc(3)   ! Number of sockets on each direction
     integer(ip) :: nedom(3)   ! Number of elements on each part and direction
     integer(ip) :: npdom(3)   ! Number of points on each part and direction
     integer(ip) :: ncorn(3)   ! Number of partition corners on each direction
     integer(ip) :: nedge(3,3) ! Number of partition edges on each direction
     integer(ip) :: nface(3,3) ! Number of partition faces on each direction
     integer(ip) :: npdomt     ! Number of points on each domain
     integer(ip) :: nedomt     ! Number of elements on each domain
     integer(ip) :: ncornt     ! Number of partition corners on each domain
     integer(ip) :: nedget(3)  ! Number of edges on each direction
     integer(ip) :: nfacet(3)  ! Number of faces on each direction
     integer(ip) :: nedgett    ! Number of edges on each domain
     integer(ip) :: nfacett    ! Number of faces on each domain
  end type geom_size

  type topo_size
     integer(ip) :: notot      ! Total amount of elemental objects of the partition
     integer(ip) :: nctot      ! Total amount of elemental corners of the partition
     integer(ip) :: ncglb      ! Total amount of elemental corners of the domain
     integer(ip) :: ndtot      ! Total amount of elemental edges of the partition
     integer(ip) :: nftot      ! Total amount of elemental faces of the partition
     integer(ip) :: nddir(3)   ! Total amount of elemental edges of the partition for each direction
     integer(ip) :: nfdir(3)   ! Total amount of elemental faces of the partition for each direction
     integer(ip) :: nddom(3,3) ! # edges for each direction given the edge direction (local)
     integer(ip) :: nfdom(3,3) ! # faces for each direction given the face normal direction (local)
     integer(ip) :: ndglb(3,3) ! # edges for each direction given the edge direction (global)
     integer(ip) :: nfglb(3,3) ! # faces for each direction given the face normal direction (global)
     integer(ip) :: ndsum(3)   ! Total amount of elemental edges of the domain for each direction
     integer(ip) :: nfsum(3)   ! Total amount of elemental faces of the domain for each direction
  end type topo_size

  !interface globalid
  !   module procedure globalid_ip, globalid_igp
  !end interface


  ! Declare constants
  integer(ip), parameter :: geom =0, topo=1 
  integer(ip), parameter :: inter=0, bound=1
  integer(ip), parameter :: do_count=0, do_list=1

!!$  ! Functions
!!$  public :: fem_mesh_gen_partition_create,fem_mesh_gen_partition_set,fem_mesh_gen_partition_bcs, &
!!$       &    global_to_ijk

  public :: geom, topo

contains
!!$  !================================================================================================
!!$  subroutine fem_mesh_gen_partition_create(lpart,ndime,isper,nedir,npdir,nsckt,msize,poin,line,surf, &
!!$       &                                   msh,parts,nodes,mater,mtype_)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),                   intent(in)    :: isper(3),nedir(3),npdir(3),nsckt(3)
!!$    integer(ip),                   intent(in)    :: ndime,lpart
!!$    type(mesh_size),               intent(in)    :: msize
!!$    type(fem_conditions),          intent(in)    :: poin,line,surf
!!$    type(fem_partition),           intent(out)   :: parts
!!$    type(fem_mesh),                intent(out)   :: msh
!!$    type(fem_conditions),          intent(out)   :: nodes
!!$    type(fem_materials), optional, intent(inout) :: mater
!!$    integer(ip),         optional, intent(in)    :: mtype_
!!$
!!$    ! Local variables
!!$    integer(ip)     :: nparts,npadj,mtype,i, tface, bou_tface, intface, bouface, shaface
!!$    integer(ip)     :: ijkpart(3)
!!$    integer(ip)     :: nobjs,neigh_list(3**ndime-1)
!!$    type(map)       :: omap
!!$    type(geom_size) :: gsize
!!$    type(topo_size) :: tsize
!!$    
!!$    ! Allocatables
!!$    integer(ip) , allocatable   :: lobjs(:,:),p(:),l(:)
!!$    integer(ip) , allocatable   :: pextn(:),lextp(:),lexte(:),lextm(:)
!!$    integer(igp), allocatable   :: lextn(:)
!!$
!!$    ! Set geometrical/topological mesh generation (0: geo, 1: topo)
!!$    if(present(mtype_)) then
!!$       mtype = mtype_
!!$    else
!!$       mtype = geom
!!$    end if
!!$ 
!!$    ! Total number of parts
!!$    nparts = npdir(1)*npdir(2)*npdir(3)
!!$
!!$    ! Geometrical sizes
!!$    call fem_structured_geom_size_create(npdir,nedir,nsckt,isper,ndime,msize,gsize)
!!$
!!$    ! Topological sizes
!!$    call fem_structured_topo_size_create(ndime,gsize,tsize)
!!$
!!$    ! Transform lpart to ijk numeration
!!$    call global_to_ijk(lpart,nsckt,gsize%npsoc,ndime,ijkpart) 
!!$         
!!$    ! Create mesh
!!$    if(mtype==topo) call fem_structured_topo_mesh_alloc(ndime,gsize,tsize,msh)
!!$    if(mtype==geom) call fem_structured_geom_mesh_alloc(ndime,gsize%nnode,gsize%nedom,msh)
!!$
!!$    ! Allocate local to global maps:
!!$    if(mtype==topo) call map_alloc(tsize%notot ,msh%npoin,parts%nmap)
!!$    if(mtype==geom) call map_alloc(gsize%npdomt,msh%npoin,parts%nmap)
!!$
!!$    call map_alloc(gsize%nedomt,int(msh%nelem,igp),parts%emap)
!!$
!!$    ! Create a void belem map 
!!$    call map_alloc(0, 0, parts%bmap)
!!$
!!$    ! Mesh, nodes map, elements map and objects map generation
!!$    call structured_mesh_gen(ijkpart,nparts,ndime,gsize,tsize,msize,isper,poin,line,surf,mtype, &
!!$         &                   msh%coord,msh%lnods,parts%nmap,parts%emap,parts%omap,pextn,lextn,lextp,  &
!!$         &                   lexte,lextm,nobjs,lobjs,p,l,nodes,intface,bouface,shaface,mater)
!!$               
!!$       ! Fill topologic face info
!!$    if(mtype==topo) then
!!$       msh%tface     = intface           ! Total amount of faces shared by 2 elements
!!$       msh%bou_tface = bouface + shaface ! Total amount of faces only in 1 element
!!$       msh%shaface   = shaface           ! Total amoutn of faces in the shared boundary of a subdomain
!!$    end if
!!$
!!$    ! Fill partition data
!!$    parts%pinfo      = interfaces 
!!$    parts%ptype      = element_based
!!$    parts%ipart      = lpart
!!$    parts%nparts     = nparts
!!$    parts%max_nparts = 2**ndime
!!$    parts%nobjs      = nobjs
!!$
!!$    ! List of neighbor parts
!!$    call neighbor_part(ijkpart,npdir,nsckt,gsize%npsoc,ndime,isper,neigh_list,npadj)
!!$    parts%npadj = npadj
!!$    call memalloc(npadj,parts%lpadj,__FILE__,__LINE__)
!!$    do i=1,npadj
!!$       parts%lpadj(i) = neigh_list(i)
!!$    end do
!!$
!!$    ! Allocate memory for the pointers to the lists of objects on each edge
!!$    call memalloc(parts%npadj+1,parts%int_objs%p,__FILE__,__LINE__)
!!$    
!!$    ! Set the pointers to the list of local objects on each edge
!!$    parts%int_objs%n = parts%npadj
!!$    parts%int_objs%p = p
!!$    call memalloc(parts%int_objs%p(parts%int_objs%n+1)-1,parts%int_objs%l,__FILE__,__LINE__)
!!$    parts%int_objs%l = l(1:(p(npadj+1)-1))
!!$
!!$    ! Allocate memory for the list of local objects
!!$    call memalloc(parts%max_nparts+4,parts%nobjs,parts%lobjs,__FILE__,__LINE__)
!!$    parts%lobjs = lobjs
!!$
!!$    ! Fill adjacencies
!!$    call memalloc(parts%emap%nb+1,parts%pextn,__FILE__,__LINE__)
!!$    call memalloc(pextn(parts%emap%nb+1)-1,parts%lextn,__FILE__,__LINE__)
!!$    call memalloc(pextn(parts%emap%nb+1)-1,parts%lextp,__FILE__,__LINE__)
!!$    call memalloc(pextn(parts%emap%nb+1)-1,parts%lexte,__FILE__,__LINE__)
!!$    call memalloc(pextn(parts%emap%nb+1)-1,parts%lextm,__FILE__,__LINE__)
!!$    parts%pextn = pextn
!!$    parts%lextn = lextn
!!$    parts%lextp = lextp
!!$    parts%lexte = lexte
!!$    parts%lextm = lextm
!!$
!!$    ! Deallocate objects
!!$    call memfree(p,__FILE__,__LINE__)
!!$    call memfree(l,__FILE__,__LINE__)
!!$    call memfree(lobjs,__FILE__,__LINE__)
!!$
!!$    ! Deallocate auxiliar adjacencies
!!$    call memfree(pextn,__FILE__,__LINE__)
!!$    call memfree(lextn,__FILE__,__LINE__)
!!$    call memfree(lextp,__FILE__,__LINE__)
!!$    call memfree(lexte,__FILE__,__LINE__)
!!$    call memfree(lextm,__FILE__,__LINE__)
!!$   
!!$  end subroutine fem_mesh_gen_partition_create
!!$
!!$  !================================================================================================
!!$  subroutine structured_mesh_gen(ijkpart,nparts,ndime,gsize,tsize,msize,isper,poin,line,surf,mtype, &
!!$       &                         coord,lnods,nmap,emap,omap,pextn,lextn,lextp,lexte,lextm,nobjs,    &
!!$       &                         lobjs,p,l,nodes,intface,bouface,shaface, mater)
!!$
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),              intent(in)    :: ijkpart(3),isper(3)
!!$    integer(ip),              intent(in)    :: nparts,ndime,mtype
!!$    type(geom_size),          intent(in)    :: gsize
!!$    type(topo_size),          intent(in)    :: tsize
!!$    type(mesh_size),          intent(in)    :: msize
!!$    type(fem_conditions),     intent(in)    :: poin,line,surf
!!$    integer(ip),              intent(inout) :: lnods(:)
!!$    real(rp),                 intent(inout) :: coord(:,:)
!!$    type(map),                intent(inout) :: nmap
!!$    type(map_igp),            intent(inout) :: emap
!!$    type(map),                intent(out)   :: omap
!!$    integer(ip),              intent(out)   :: nobjs
!!$    type(fem_conditions),     intent(out)   :: nodes
!!$    integer(ip),              intent(out)   :: intface,bouface,shaface
!!$    integer(ip), allocatable, intent(out)   :: pextn(:),lextp(:),lexte(:),lextm(:)
!!$    integer(igp), allocatable, intent(out)  :: lextn(:)
!!$    integer(ip), allocatable, intent(out)   :: lobjs(:,:),p(:),l(:)
!!$    type(fem_materials), optional, intent(inout) :: mater
!!$
!!$    ! Local variables
!!$    integer(ip)              :: cnt(2),subgl(ndime),case,ni,nb,i,j,k,lpart,pni,eni,pnb,enb,pdime
!!$    integer(ip)              :: lobjs_int(4+2**ndime)
!!$    integer(ip)              :: lobjs_face(4+2**ndime,6),l2g_face(2*ndime)
!!$    integer(ip)              :: lobjs_edge(4+2**ndime,2*(ndime-1)*ndime),l2g_edge(2*(ndime-1)*ndime)
!!$    integer(ip)              :: lobjs_corn(4+2**ndime,4*(ndime-1)),l2g_corn(4*(ndime-1))
!!$    integer(ip)              :: nobjs_face,nobjs_edge,nobjs_corn,nfacet_tot,nedget_tot,nobjs_gl
!!$    integer(ip)              :: touch_face(3**ndime),touch_edge(3**ndime),touch_corn(3**ndime)
!!$    integer(ip)              :: neigh_obj_face(2*ndime,3**ndime)
!!$    integer(ip)              :: neigh_obj_edge(2*(ndime-1)*ndime,3**ndime)
!!$    integer(ip)              :: neigh_obj_corn(4*(ndime-1),3**ndime)
!!$    integer(ip)              :: neigh_obj_list(3**ndime,3**ndime),touch(3**ndime)
!!$    integer(ip)              :: neigh_obj(3**ndime,3**ndime)
!!$    integer(ip)              :: p_aux(3**ndime)
!!$    integer(ip)              :: pextn_cnt(2),nelbl(3)
!!$
!!$    ! Local allocatables 
!!$    integer(ip)              :: idime,aux_isp
!!$    integer(ip) , allocatable :: npnum(:),nenum(:),l2gp(:)    
!!$    integer(igp), allocatable :: l2ge(:)    
!!$
!!$    ! Allocate 
!!$    if(mtype==topo) then
!!$       call memalloc(tsize%notot,l2gp,__FILE__,__LINE__)
!!$       call memalloc(tsize%notot,npnum,__FILE__,__LINE__)
!!$    elseif(mtype==geom) then
!!$       call memalloc(gsize%npdomt,l2gp,__FILE__,__LINE__)
!!$       call memalloc(gsize%npdomt,npnum,__FILE__,__LINE__)
!!$    end if
!!$    call memalloc(gsize%nedomt,l2ge,__FILE__,__LINE__)
!!$    call memalloc(gsize%nedomt,nenum,__FILE__,__LINE__)
!!$
!!$    ! Initialize counter
!!$    cnt(1) = 1; cnt(2) = 1
!!$
!!$    intface   = 0  ! Total amount of faces in the interior of a subdomain
!!$    bouface   = 0  ! Total amount of faces in the non-shared boundary of a subdomain
!!$    shaface   = 0  ! Total amoutn of faces in the shared boundary of a subdomain
!!$
!!$    ! Part id
!!$    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)
!!$
!!$    ! Create bounary conditions
!!$    if(mtype==geom) call fem_conditions_create(poin%ncode,poin%nvalu,gsize%npdomt,nodes)
!!$    if(mtype==topo) call fem_conditions_create(poin%ncode,poin%nvalu,tsize%notot ,nodes)
!!$   
!!$    ! Interior nodes
!!$
!!$    call volu_loop(ijkpart,ndime,mtype,gsize,tsize,msize,npnum,nenum,coord,cnt, intface)
!!$
!!$    call face_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,surf,npnum,nenum,coord,cnt,inter, &
!!$         &         lobjs_face,nobjs_face,l2g_face,neigh_obj_face,touch_face,bouface, nodes)
!!$
!!$    call edge_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,line,surf,npnum,nenum,coord,cnt,  &
!!$         &         inter,lobjs_edge,nobjs_edge,l2g_edge,neigh_obj_edge,touch_edge, bouface, nodes)
!!$
!!$
!!$    call corn_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,poin,line,surf,npnum,nenum,coord, &
!!$         &         cnt,inter,lobjs_corn,nobjs_corn,l2g_corn,neigh_obj_corn,touch_corn,nodes)
!!$    pni = cnt(1)-1
!!$    eni = cnt(2)-1
!!$
!!$    ! Interior objects
!!$    call inte_object(ijkpart,gsize%nsckt,gsize%npsoc,ndime,cnt,lobjs_int)
!!$    ! Boundary nodes/objects
!!$    call face_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,surf,npnum,nenum,coord,cnt,bound, &
!!$         &         lobjs_face,nobjs_face,l2g_face,neigh_obj_face,touch_face, shaface) 
!!$
!!$    call edge_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,line,surf,npnum,nenum,coord,cnt,  &
!!$         &         bound,lobjs_edge,nobjs_edge,l2g_edge,neigh_obj_edge,touch_edge, shaface, nodes)
!!$
!!$    call corn_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,poin,line,surf,npnum,nenum,coord, &
!!$         &         cnt,bound,lobjs_corn,nobjs_corn,l2g_corn,neigh_obj_corn,touch_corn, nodes)
!!$    pnb = cnt(1) - pni - 1
!!$    enb = cnt(2) - eni - 1
!!$
!!$    ! Global part numeration
!!$    do i=1,ndime
!!$       subgl(i) = (ijkpart(i)-1)*gsize%nedom(i)
!!$    end do
!!$
!!$    ! Initialize materia
!!$    if(msize%mater>0 .and. present(mater)) then
!!$       call fem_materials_create(gsize%nedomt,mater)
!!$    end if
!!$
!!$    ! Calculate lnods and l2g vector for emap and nmap
!!$    call generic_l2g(subgl,npnum,nenum,ndime,mtype,gsize,tsize,msize,l2ge,l2gp,lnods,mater)
!!$
!!$    ! Fill nmap
!!$    nmap%ni  = pni
!!$    nmap%nb  = pnb
!!$    nmap%ne  = 0
!!$    nmap%l2g = l2gp
!!$
!!$    ! Fill emap
!!$    emap%ni  = eni
!!$    emap%nb  = enb
!!$    emap%ne  = 0
!!$    emap%l2g = l2ge
!!$    
!!$    ! Compute adjacencies
!!$    call memalloc(enb+1,pextn,__FILE__,__LINE__)
!!$    nelbl(1) = msize%neblx
!!$    nelbl(2) = msize%nebly
!!$    nelbl(3) = msize%neblz
!!$    pextn_cnt = 0
!!$    pextn(1) = 1
!!$    call face_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,0)
!!$    call edge_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,0)
!!$    call corn_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,0)
!!$    pextn_cnt = 0
!!$    call memalloc(pextn(enb+1)-1,lextn,__FILE__,__LINE__)
!!$    call memalloc(pextn(enb+1)-1,lextp,__FILE__,__LINE__)
!!$    call memalloc(pextn(enb+1)-1,lexte,__FILE__,__LINE__)
!!$    call memalloc(pextn(enb+1)-1,lextm,__FILE__,__LINE__)
!!$    lextn = 0; lextp = 0; lexte = 0; lextm = 0
!!$    call face_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,1, &
!!$         &             lextn,lextp,lexte,lextm,msize%mater,nelbl)
!!$    call edge_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,1, &
!!$         &             lextn,lextp,lexte,lextm,msize%mater,nelbl)
!!$    call corn_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,emap%nb,pextn_cnt,pextn,1, &
!!$         &             lextn,lextp,lexte,lextm,msize%mater,nelbl)
!!$
!!$    ! Total number of objects
!!$    nobjs = 1 + nobjs_face + nobjs_edge + nobjs_corn
!!$    nobjs_gl = nparts + gsize%nfacett + gsize%nedgett + gsize%ncornt
!!$
!!$    ! Fill omap
!!$    call map_alloc(nobjs,nobjs_gl,omap)
!!$    omap%ni  = 1
!!$    omap%nb  = nobjs - 1
!!$    omap%ne  = 0
!!$
!!$    ! Auxiliar integer
!!$    aux_isp = 1
!!$    do idime=1,ndime
!!$       aux_isp = aux_isp*isper(idime)
!!$    end do
!!$
!!$    ! List of objects and l2g vector for omap
!!$    call memalloc((4+2**ndime),nobjs,lobjs,__FILE__,__LINE__)
!!$    lobjs_int(6:4+2**ndime) = 0
!!$    neigh_obj = 0
!!$    lobjs(:,1) = lobjs_int(:)
!!$    omap%l2g(1) = lpart
!!$    do i=1,nobjs_face
!!$       lobjs(:,i+1) = lobjs_face(:,i)
!!$       omap%l2g(i+1) = l2g_face(i) + nparts
!!$       neigh_obj(i+1,:) = neigh_obj_face(i,:)
!!$    end do
!!$    do i=1,nobjs_edge
!!$       lobjs(:,i+nobjs_face+1) = lobjs_edge(:,i)
!!$       omap%l2g(i+nobjs_face+1) = l2g_edge(i) + nparts + gsize%nfacett
!!$       neigh_obj(i+nobjs_face+1,:) = neigh_obj_edge(i,:)
!!$    end do
!!$    do i=1,nobjs_corn
!!$       lobjs(:,i+nobjs_edge+nobjs_face+1) = lobjs_corn(:,i)
!!$       omap%l2g(i+nobjs_edge+nobjs_face+1) = l2g_corn(i) + nparts + gsize%nfacett + gsize%nedgett
!!$       neigh_obj(i+nobjs_edge+nobjs_face+1,:) = neigh_obj_corn(i,:)
!!$       if(l2g_corn(i)==1.and.aux_isp==1.and.nodes%ncode==(ndime+1)) then
!!$          ! Fix pressure at one point when periodic bc's (Only for NSI)
!!$          nodes%code(ndime+1,lobjs_corn(2,i)) = 1
!!$          nodes%valu(ndime+1,lobjs_corn(2,i)) = 0.0_rp
!!$       end if
!!$    end do
!!$
!!$    ! List of local objects
!!$    ! Pointers
!!$    p_aux=0
!!$    touch = touch_face + touch_edge + touch_corn
!!$    p_aux(1) = 1
!!$    k = 2
!!$    do i=1,(3**ndime-1)/2
!!$       if(touch(i)/=0) then
!!$          p_aux(k) = p_aux(k-1) + touch(i)
!!$          k = k+1
!!$       end if
!!$    end do
!!$    do i=(3**ndime+1)/2+1,3**ndime
!!$       if(touch(i)/=0) then
!!$          p_aux(k) = p_aux(k-1) + touch(i)
!!$          k = k+1
!!$       end if
!!$    end do
!!$
!!$    ! Select nonzero pointers
!!$    k=0
!!$    do i=1,3**ndime
!!$       if(p_aux(i)/=0) then
!!$          k=k+1
!!$       end if
!!$    end do
!!$    call memalloc(k,p,__FILE__,__LINE__)
!!$    p = p_aux(1:k)
!!$
!!$    ! Objects list
!!$    call memalloc(p(size(p))-1,l,__FILE__,__LINE__)
!!$    do i=1,3**ndime
!!$       k=1
!!$       do j=1,nobjs
!!$          if(neigh_obj(j,i)==1) then
!!$             neigh_obj_list(i,k)=j
!!$             k = k+1
!!$          end if
!!$       end do
!!$    end do
!!$    k=1
!!$    do i=1,(3**ndime-1)/2
!!$       do j=1,touch(i)
!!$          if(touch(i)/=0) then
!!$             l(k) = neigh_obj_list(i,j) 
!!$             k=k+1
!!$          end if
!!$       end do
!!$    end do
!!$    do i=(3**ndime+1)/2+1,3**ndime
!!$       do j=1,touch(i)
!!$          if(touch(i)/=0) then
!!$             l(k) = neigh_obj_list(i,j) 
!!$             k = k+1
!!$          end if
!!$       end do
!!$    end do
!!$
!!$    ! Deallocate auxiliar vectors
!!$    call memfree( l2gp,__FILE__,__LINE__)
!!$    call memfree( l2ge,__FILE__,__LINE__)
!!$    call memfree(npnum,__FILE__,__LINE__)
!!$    call memfree(nenum,__FILE__,__LINE__)
!!$
!!$  end subroutine structured_mesh_gen
!!$
!!$  !================================================================================================
!!$  subroutine neighbor_part(ijkpart,npdir,nsckt,npsoc,ndime,isper,neigh_list,npadj)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),npdir(3),isper(3),nsckt(3),npsoc(3),ndime
!!$    integer(ip), intent(out) :: neigh_list(3**ndime-1),npadj
!!$
!!$    integer(ip)              :: cnt,i,j,k,ijkpartn(3),neigh,k_aux
!!$   
!!$    cnt = 1
!!$    npadj = 0
!!$    neigh_list = 0
!!$    k_aux=1
!!$
!!$    do i=-1,1
!!$       do j=-1,1
!!$          do k=-1,1
!!$
!!$             if(ndime==2) k_aux=0
!!$
!!$             ! Identifier
!!$             ijkpartn = ijkpart + (/i,j,k/)
!!$
!!$             ! Check boundary
!!$             call check_part_boundary(ijkpartn,npdir,ndime,isper)
!!$        
!!$             ! Global numbering
!!$             call globalid_2l(ijkpartn,nsckt,npsoc,ndime,neigh)
!!$                        
!!$             if(neigh>0.and.((i/=0).or.(j/=0).or.((k*k_aux)/=0))) then
!!$                neigh_list(cnt) = neigh
!!$                npadj = npadj + 1
!!$                ! Update counter
!!$                cnt = cnt + 1
!!$             end if
!!$
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine neighbor_part
!!$
!!$  !================================================================================================
!!$
!!$  subroutine volu_loop(ijkpart,ndime,mtype,gsize,tsize,msize,npnum,nenum,coord,cnt, nface)
!!$
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),     intent(in)    :: ijkpart(3),ndime,mtype
!!$    type(geom_size), intent(in)    :: gsize
!!$    type(topo_size), intent(in)    :: tsize
!!$    type(mesh_size), intent(in)    :: msize
!!$    integer(ip),     intent(inout) :: npnum(:),nenum(:),cnt(2), nface
!!$
!!$    real(rp),        intent(inout) :: coord(:,:)
!!$
!!$    integer(ip) :: i,j,k,ijkpoin(3),ijkelem(3),glnum,npdomk,nedomk,nddomk,aux_cnt
!!$    integer(ip) :: pdime,ijkface(3),ijkedge(3),jdime
!!$    integer(ip), allocatable :: auxv(:,:)
!!$
!!$    ! Auxiliar
!!$    if(ndime==2) then
!!$       npdomk=2
!!$       nedomk=2
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$    else if(ndime==3) then
!!$       npdomk = gsize%npdom(3) - 1
!!$       nedomk = gsize%nedom(3) - 1
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$    end if
!!$    
!!$    ! Interior nodes loop
!!$    do i=2,gsize%npdom(1)-1
!!$       do j=2,gsize%npdom(2)-1
!!$          do k=2,npdomk
!!$
!!$             ! Point identifier (inside the part)
!!$             ijkpoin = (/i,j,k/)
!!$             call globalid(ijkpoin,gsize%npdom,ndime,glnum)
!!$
!!$             ! Global to local numeration (inside the part)
!!$             npnum(glnum) = cnt(1)
!!$
!!$             ! Coordinates
!!$             if(mtype==geom) call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!!$                  &                         ndime,msize,coord(:,cnt(1)))
!!$
!!$             ! Update counter
!!$             cnt(1) = cnt(1) + 1
!!$
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    if(mtype==topo) then
!!$       ! Interior edges loop
!!$       aux_cnt = tsize%nctot
!!$       do pdime=ndime,1,-1
!!$          if(ndime==2) then
!!$             nddomk=2
!!$          elseif(ndime==3) then
!!$             nddomk=tsize%nddom(auxv(pdime,2),pdime)-1
!!$          end if
!!$          do i=1,tsize%nddom(pdime,pdime)
!!$             do j=2,tsize%nddom(auxv(pdime,1),pdime)-1
!!$                do k=2,nddomk
!!$
!!$                   ! Face identifier (inside the part)
!!$                   ijkedge(pdime) = i
!!$                   if(ndime==2) then
!!$                      ijkedge(auxv(pdime,:)) =  j
!!$                   else
!!$                      ijkedge(auxv(pdime,:)) =  (/j,k/)
!!$                   end if
!!$                   call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)
!!$
!!$                   ! Global to local numeration (inside the part)
!!$                   npnum(aux_cnt+glnum) = cnt(1)
!!$
!!$                   ! Update counter
!!$                   cnt(1) = cnt(1) + 1
!!$                   
!!$                   ! Update face counter
!!$                   if(ndime==2) nface = nface+1
!!$
!!$                end do
!!$             end do
!!$          end do
!!$          aux_cnt = aux_cnt + tsize%nddir(pdime)
!!$       end do
!!$
!!$       ! Interior faces loop
!!$       aux_cnt = tsize%nctot + tsize%ndtot
!!$       if(ndime==3) then
!!$          do pdime=ndime,1,-1
!!$             do i=2,tsize%nfdom(pdime,pdime)-1
!!$                do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!!$                   do k=1,tsize%nfdom(auxv(pdime,2),pdime)
!!$
!!$                      ! Face identifier (inside the part)
!!$                      ijkface(pdime) = i
!!$                      ijkface(auxv(pdime,:)) = (/j,k/)
!!$                      call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)
!!$
!!$                      ! Global to local numeration (inside the part)
!!$                      npnum(aux_cnt+glnum) = cnt(1)
!!$
!!$                      ! Update counter
!!$                      cnt(1) = cnt(1) + 1
!!$
!!$                      ! Update face counter
!!$                      nface = nface+1
!!$
!!$                   end do
!!$                end do
!!$             end do
!!$             aux_cnt = aux_cnt + tsize%nfdir(pdime)
!!$          end do
!!$       end if
!!$    end if
!!$    
!!$    ! Interior elements loop
!!$    do i=2,gsize%nedom(1)-1
!!$       do j=2,gsize%nedom(2)-1
!!$          do k=2,nedomk
!!$
!!$             ! Element identifier (inside the part)
!!$             ijkelem = (/i,j,k/)
!!$             call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$             ! Global to local numeration (inside the part)
!!$             nenum(glnum) = cnt(2)
!!$
!!$             ! Update counter
!!$             cnt(2) = cnt(2) + 1
!!$
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    ! Deallocate auxiliar vector
!!$    deallocate(auxv)
!!$
!!$  end subroutine volu_loop
!!$
!!$  !================================================================================================
!!$  subroutine face_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,surf,npnum,nenum,coord,cnt,case, &
!!$       &               lobjs_face,obj,l2g_face,neigh_obj,touch, nface, nodes)
!!$
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),          intent(in)    :: ijkpart(3),isper(3),ndime,case,mtype
!!$    type(geom_size),      intent(in)    :: gsize
!!$    type(topo_size),      intent(in)    :: tsize
!!$    type(mesh_size),      intent(in)    :: msize
!!$    type(fem_conditions), intent(in)    :: surf
!!$    integer(ip),          intent(inout) :: npnum(:),nenum(:),cnt(2), nface
!!$    real(rp),             intent(inout) :: coord(:,:)
!!$    integer(ip),          intent(out)   :: lobjs_face(4+2**ndime,6),obj,l2g_face(2*ndime)
!!$    integer(ip),          intent(out)   :: neigh_obj(2*ndime,3**ndime),touch(3**ndime)
!!$    type(fem_conditions), optional, intent(inout) :: nodes
!!$    
!!$    integer(ip) :: auxv(3,2),i,j,k,pdime,iface,ijkpoin(3),ijkelem(3),ijkface(3),glnum,flag,neigh(2)
!!$    integer(ip) :: inipo,lpart,glface,lface(2),ijklocal(3),glocal,surf_cnt,nface_aux(3),aux_vec(3)
!!$    integer(ip) :: idir,jdir,idime,jdime,aux_cnt_f,aux_cnt_d,ijkedge(3)
!!$
!!$    ! Auxiliar vector of components
!!$    auxv(1,:) = (/2,3/)
!!$    auxv(2,:) = (/1,3/)
!!$    auxv(3,:) = (/1,2/)
!!$
!!$    ! Objects counter
!!$    obj = 0
!!$    lobjs_face = 0
!!$    neigh_obj = 0
!!$    touch = 0
!!$
!!$    ! Face nodes loop
!!$    surf_cnt = 1
!!$    aux_cnt_f = tsize%nctot + tsize%ndtot
!!$    if(ndime==3) then
!!$       do pdime = ndime,1,-1
!!$          do iface = 0,1
!!$
!!$             ! Check if it is a boundary
!!$             call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!!$                  &                   isper,auxv,flag,neigh,lface)
!!$
!!$             ! Check case
!!$             if(flag==case) then
!!$
!!$                ! Initial object point
!!$                inipo = cnt(1)
!!$
!!$                ! Point identifier (inside the part)
!!$                ijkpoin(pdime) = (gsize%npdom(pdime)-1)*iface+1
!!$
!!$                do i=2,gsize%npdom(auxv(pdime,1))-1
!!$                   do j=2,gsize%npdom(auxv(pdime,2))-1
!!$
!!$                      ! Point identifier (inside the part)
!!$                      ijkpoin(auxv(pdime,:)) = (/i,j/)
!!$                      call globalid(ijkpoin,gsize%npdom,ndime,glnum)
!!$
!!$                      ! Global to local numeration (inside the part)
!!$                      npnum(glnum) = cnt(1)
!!$
!!$                      ! Coordinates
!!$                      if(mtype==geom) call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!!$                           &                         ndime,msize,coord(:,cnt(1)))
!!$
!!$                      ! Boundary conditions
!!$                      if(case==inter) then
!!$                         nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!!$                         nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!!$                      end if
!!$
!!$                      ! Update point counter
!!$                      cnt(1) = cnt(1) + 1
!!$
!!$                   end do
!!$                end do
!!$
!!$                if(mtype==topo) then
!!$                   ! Elemental edges loop
!!$                   do idir = 1,2
!!$                      jdir = 3-idir
!!$                      idime = auxv(pdime,idir)
!!$                      i = (tsize%nddom(pdime,idime)-1)*iface+1
!!$                      aux_cnt_d = tsize%nctot
!!$                      do jdime=1,3-idime
!!$                         aux_cnt_d = aux_cnt_d + tsize%nddir(4-jdime)
!!$                      end do
!!$                      do j=1,tsize%nddom(idime,idime)
!!$                         do k=2,tsize%nddom(auxv(pdime,jdir),idime)-1
!!$
!!$                            ! Face identifier (inside the part)
!!$                            ijkedge(pdime) = i
!!$                            ijkedge(auxv(pdime,idir)) = j
!!$                            ijkedge(auxv(pdime,jdir)) = k
!!$                            call globalid(ijkedge,tsize%nddom(:,idime),ndime,glnum)
!!$
!!$                            ! Global to local numeration (inside the part)
!!$                            npnum(aux_cnt_d+glnum) = cnt(1)
!!$
!!$                            ! Boundary conditions
!!$                            if(case==inter) then
!!$                               nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!!$                               nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!!$                            end if
!!$
!!$                            ! Update counter
!!$                            cnt(1) = cnt(1) + 1
!!$
!!$                         end do
!!$                      end do
!!$                   end do
!!$
!!$                   ! Elemental faces loop
!!$                   i = (tsize%nfdom(pdime,pdime)-1)*iface+1
!!$                   do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!!$                      do k=1,tsize%nfdom(auxv(pdime,2),pdime)
!!$
!!$                         ! Face identifier (inside the part)
!!$                         ijkface(pdime) = i
!!$                         ijkface(auxv(pdime,:)) = (/j,k/)
!!$                         call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)
!!$
!!$                         ! Global to local numeration (inside the part)
!!$                         npnum(aux_cnt_f+glnum) = cnt(1)
!!$
!!$                         ! Boundary conditions
!!$                         if(case==inter) then
!!$                            nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
!!$                            nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
!!$                         end if
!!$
!!$                         ! Update counter
!!$                         cnt(1) = cnt(1) + 1
!!$
!!$                         ! Update face counter
!!$                         nface = nface+1
!!$
!!$                      end do
!!$                   end do
!!$                end if
!!$
!!$                if(case==bound) then
!!$
!!$                   ! Update counter
!!$                   obj = obj + 1
!!$
!!$                   ! Object list
!!$                   lobjs_face(1,obj) = 0              ! separator id, useless in the structured case
!!$                   lobjs_face(2,obj) = inipo          ! Starting node in local numbering 
!!$                   lobjs_face(3,obj) = cnt(1)-1       ! End node in local numbering
!!$                   lobjs_face(4,obj) = 2              ! number of subdomains sharing the object
!!$                   lobjs_face(5,obj) = neigh(1)       ! 1st subdomain sharing the object
!!$                   lobjs_face(6,obj) = neigh(2)       ! 2nd subdomain sharing the object
!!$!                   call psb_hsort(lobjs_face(5:6,obj))
!!$                   call sort( 1, lobjs_face(5:6,obj) )
!!$
!!$                   ! Global face
!!$                   !if(ijkpart(pdime)==npdir(pdime).and.isper(pdime)==1.and.iface==1) then
!!$                   !   ijkface(pdime) = 1
!!$                   !else
!!$                   !   ijkface(pdime) = ijkpart(pdime) + iface
!!$                   !end if
!!$                   ijkface(pdime) = ijkpart(pdime) + iface
!!$                   if(ijkface(pdime)>gsize%nface(pdime,pdime)) ijkface(pdime) = 1
!!$                   ijkface(auxv(pdime,:)) = ijkpart(auxv(pdime,:))
!!$                   do i=1,ndime-1
!!$                      if(ijkface(auxv(pdime,i))>gsize%nface(pdime,auxv(pdime,i))) then
!!$                         ijkface(auxv(pdime,i)) = 1
!!$                      end if
!!$                   end do
!!$                   nface_aux = gsize%nface(pdime,:)
!!$                   call globalid(ijkface,nface_aux,ndime,glface)
!!$                   if(pdime==ndime) then
!!$                      l2g_face(obj) = glface
!!$                   else if(pdime==ndime-1) then
!!$                      l2g_face(obj) = gsize%nfacet(pdime+1) + glface
!!$                   else if(pdime==ndime-2) then
!!$                      l2g_face(obj) = gsize%nfacet(pdime+1) + gsize%nfacet(pdime+2) + glface
!!$                   end if
!!$
!!$                   ! Local objects
!!$                   do i=1,2
!!$                      if(lface(i)==1) then
!!$                         ijklocal(pdime) = i + iface
!!$                         ijklocal(auxv(pdime,:)) = 2
!!$                         aux_vec = (/3,3,3/)
!!$                         call globalid(ijklocal,aux_vec,ndime,glocal)
!!$                         touch(glocal) = touch(glocal) + 1
!!$                         neigh_obj(obj,glocal) = 1
!!$                      end if
!!$                   end do
!!$                end if
!!$
!!$                if(iface.gt.min(1,abs(gsize%nedom(pdime)-1))) cycle
!!$
!!$                ! Element identifier (inside the part)
!!$                ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1
!!$
!!$                do i=2,gsize%nedom(auxv(pdime,1))-1
!!$                   do j=2,gsize%nedom(auxv(pdime,2))-1
!!$
!!$                      ! Point identifier (inside the part)
!!$                      ijkelem(auxv(pdime,:)) = (/i,j/)
!!$                      call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                      ! Global to local numeration (inside the part)
!!$                      nenum(glnum) = cnt(2)
!!$
!!$                      ! Update element counter
!!$                      cnt(2) = cnt(2) + 1
!!$
!!$                   end do
!!$                end do
!!$             end if
!!$
!!$             ! Update surface counter
!!$             surf_cnt = surf_cnt + 1
!!$
!!$          end do
!!$          
!!$          ! Update elemental faces counter
!!$          aux_cnt_f = aux_cnt_f + tsize%nfdir(pdime)
!!$       end do
!!$    end if
!!$
!!$  end subroutine face_loop
!!$
!!$ !================================================================================================
!!$  subroutine face_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,nb,pextn_cnt,pextn,case, &
!!$       &                   lextn,lextp,lexte,lextm,mcase,nelbl)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),           intent(in)    :: ijkpart(3),isper(3),ndime,case,nb,mtype
!!$    type(geom_size),       intent(in)    :: gsize
!!$    type(topo_size),       intent(in)    :: tsize
!!$    integer(ip),           intent(inout) :: pextn_cnt(2),pextn(nb+1)
!!$    integer(ip), optional, intent(inout) :: lextp(:),lexte(:),lextm(:)
!!$    integer(igp), optional, intent(inout):: lextn(:)
!!$    integer(ip), optional, intent(in)    :: mcase,nelbl(3)
!!$    
!!$    integer(ip)  :: auxv(3,2),ijkelem(3),neigh(2),ijkneigh(3),subgl(3),subgl_ijk(3),ijkpart_neigh(3)
!!$    integer(ip)  :: i,j,pdime,iface,ielem,jelem,glnum,flag,gpart,lface(2),k,kedge,kdime,idime
!!$    integer(igp) :: gneigh 
!!$    integer(ip)  :: knode,mnode,auxva(2**ndime),auxk(3),auxkv(3,2),nnode(3),ijknode(3),gnode,ipos
!!$    integer(ip)  :: ijelem(2),auxvd(3,4),auxvf(3,2),auxdf(3,3)
!!$
!!$    if(ndime==3) then
!!$
!!$    ! Auxiliar vector of components
!!$    auxv(1,:) = (/2,3/)
!!$    auxv(2,:) = (/1,3/)
!!$    auxv(3,:) = (/1,2/)
!!$    ! auxva = (/1,5,4,8,2,6,3,7/)
!!$    auxva = (/1,5,3,7,2,6,4,8/)
!!$    auxk = (/1,2,1/)
!!$    auxkv(1,:) = (/0,0/)
!!$    auxkv(2,:) = (/0,1/)    
!!$    auxkv(3,:) = (/1,0/)
!!$    nnode = (/2,2,2/)
!!$    ! auxvd(1,:) = (/9,17,11,19/)!(/9,13,11,15/)
!!$    ! auxvd(2,:) = (/12,20,10,18/)!(/12,16,10,14/)
!!$    ! auxvd(3,:) = (/13,16,14,15/)!(/17,20,18,19/)
!!$    auxvd(1,:) = (/9,11,10,12/)!(/9,13,11,15/)
!!$    auxvd(2,:) = (/13,15,14,16/)!(/12,16,10,14/)
!!$    auxvd(3,:) = (/17,19,18,20/)!(/17,20,18,19/)
!!$    ! auxvf(1,:) = (/25,23/)
!!$    ! auxvf(2,:) = (/22,24/)
!!$    ! auxvf(3,:) = (/21,26/)
!!$    auxvf(1,:) = (/25,26/)
!!$    auxvf(2,:) = (/23,24/)
!!$    auxvf(3,:) = (/21,22/)
!!$    auxdf = reshape((/0,2,2,1,0,2,1,1,0/),(/3,3/))
!!$
!!$    ! Face nodes loop
!!$    do pdime = ndime,1,-1
!!$       do iface = 0,1
!!$
!!$          ! Check if it is a boundary
!!$          call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!!$               &                   isper,auxv,flag,neigh,lface)
!!$
!!$          ! Only boundary elements
!!$          if(flag==bound) then
!!$
!!$             if(case==do_count) then
!!$
!!$                do i=2,gsize%nedom(auxv(pdime,1))-1
!!$                   do j=2,gsize%nedom(auxv(pdime,2))-1
!!$
!!$                      ! Update boundary element counter
!!$                      pextn_cnt(1) = pextn_cnt(1) + 1
!!$
!!$                      do ielem = -1,1
!!$                         do jelem = -1,1
!!$
!!$                            ! Update boundary neighbour counter
!!$                            pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                         end do
!!$                      end do
!!$
!!$                      ! Fill pextn
!!$                      pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1
!!$
!!$                   end do
!!$                end do
!!$
!!$             elseif(case==do_list) then
!!$
!!$                ! Element identifier (inside the part)
!!$                ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1
!!$
!!$                ! Neighbour element identifier
!!$                ijkneigh(pdime) = (gsize%nedom(pdime)-1)*(1-iface)+1
!!$
!!$                ! Neighbour part
!!$                ijkpart_neigh = ijkpart
!!$                ijkpart_neigh(pdime) = ijkpart(pdime) + (2*iface-1)
!!$                call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!!$                call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!!$
!!$                ! Global part numeration
!!$                do i=1,ndime
!!$                   subgl(i) = (ijkpart_neigh(i)-1)*gsize%nedom(i)
!!$                end do
!!$
!!$                do i=2,gsize%nedom(auxv(pdime,1))-1
!!$                   do j=2,gsize%nedom(auxv(pdime,2))-1
!!$
!!$                      ! Update boundary element counter
!!$                      pextn_cnt(1) = pextn_cnt(1) + 1
!!$
!!$                      ! Point identifier (inside the part)
!!$                      ijkelem(auxv(pdime,:)) = (/i,j/)
!!$                      call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                      do ielem = -1,1
!!$                         do jelem = -1,1
!!$
!!$                            ! Update boundary neighbour counter
!!$                            pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                            ! Neighbour element identifier (local)
!!$                            ijkneigh(auxv(pdime,:)) = (/i+ielem,j+jelem/)
!!$                            subgl_ijk = subgl + ijkneigh
!!$                            call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)
!!$
!!$                            ! Fill adjacencies info
!!$                            lextn(pextn_cnt(2)) = gneigh
!!$                            lextp(pextn_cnt(2)) = gpart
!!$                            lextm(pextn_cnt(2)) = 1
!!$                            call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,lextm(pextn_cnt(2)))
!!$                            lexte(pextn_cnt(2)) = 0
!!$
!!$                            ! Construct bit maps
!!$                            ! Nodes
!!$                            do knode=1,auxk(ielem+2)
!!$                               do mnode=1,auxk(jelem+2)
!!$
!!$                                  ! ijk node identification
!!$                                  ijknode(pdime) = iface + 1
!!$                                  ijknode(auxv(pdime,:)) = (/auxkv(ielem+2,knode),auxkv(jelem+2,mnode)/) + (/1,1/)
!!$
!!$                                  ! Global element numeration
!!$                                  call globalid(ijknode,nnode,ndime,gnode)
!!$
!!$                                  ! Fill bit map
!!$                                  lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxva(gnode)-1)
!!$
!!$                               end do
!!$                            end do
!!$
!!$                            if(mtype==topo) then
!!$                               ! Edges
!!$                               ijelem = (/jelem,ielem/)
!!$                               if(ielem*jelem==0) then
!!$                                  do kdime=1,2
!!$                                     if(ijelem(3-kdime)==0) then
!!$                                        do kedge=1,auxk(ijelem(kdime)+2)
!!$
!!$                                           ! ijk edge identification
!!$                                           ijknode(auxdf(pdime,auxv(pdime,3-kdime))) = iface + 1
!!$                                           ijknode(3-auxdf(pdime,auxv(pdime,3-kdime))) = auxkv(ijelem(kdime)+2,kedge) + 1
!!$
!!$                                           ! Global element numeration
!!$                                           call globalid(ijknode(1:2),nnode,2,gnode)
!!$
!!$                                           ! Fill bit map
!!$                                           lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvd(auxv(pdime,kdime),gnode)-1)
!!$
!!$                                        end do
!!$                                     end if
!!$                                  end do
!!$                               end if
!!$                               ! Faces
!!$                               if(ielem==0.and.jelem==0) then
!!$
!!$                                  ! Fill bit map
!!$                                  lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvf(pdime,iface+1)-1)
!!$
!!$                               end if
!!$                            end if
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end do
!!$             end if
!!$          end if
!!$       end do
!!$    end do
!!$
!!$  end if
!!$
!!$  end subroutine face_loop_adj
!!$
!!$  !================================================================================================
!!$  subroutine edge_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,line,surf,npnum,nenum,coord,cnt, &
!!$       &               case,lobjs_edge,obj,l2g_edge,neigh_obj,touch, nface, nodes)
!!$
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),          intent(in)    :: ijkpart(3),isper(3),ndime,case,mtype
!!$    type(geom_size),      intent(in)    :: gsize
!!$    type(topo_size),      intent(in)    :: tsize
!!$    type(mesh_size),      intent(in)    :: msize
!!$    type(fem_conditions), intent(in)    :: line,surf
!!$    integer(ip),          intent(inout) :: npnum(:),nenum(:),cnt(2), nface
!!$    real(rp),             intent(inout) :: coord(:,:)
!!$    integer(ip),          intent(out)   :: lobjs_edge(4+2**ndime,2*(ndime-1)*ndime),obj
!!$    integer(ip),          intent(out)   :: l2g_edge(2*(ndime-1)*ndime)
!!$    integer(ip),          intent(out)   :: neigh_obj(2*(ndime-1)*ndime,3**ndime),touch(3**ndime)
!!$    type(fem_conditions), intent(inout) :: nodes
!!$ 
!!$    integer(ip) :: i,j,k,pdime,iedge,jedge,ijedge(2),ijkpoin(3),ijkelem(3),ijkvert(3),ijkedge(3)
!!$    integer(ip) :: neigh(2*(ndime-1)),nedge_aux(ndime),aux_vec(3),jdime
!!$    integer(ip) :: glnum,flag,kron(ndime),gsurf,inode,aux_surf(4,2),gsurf_aux(4),aux_cnt,aux
!!$    integer(ip) :: inipo,lpart,work(2*(ndime-1)),gledge,ledge(2,2),ijklocal(3),glocal,line_cnt
!!$    integer(ip), allocatable :: auxv(:,:)
!!$
!!$    ! Auxiliar vector of components
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$    end if
!!$
!!$    ! Objects counter
!!$    obj = 0
!!$    lobjs_edge = 0
!!$    neigh_obj = 0
!!$    touch = 0
!!$    line_cnt = 0
!!$    aux_cnt = tsize%nctot
!!$
!!$    ! Edge nodes loop
!!$    do pdime = ndime,1,-1
!!$       do iedge = 0,1!min(1,nedom(auxv(pdime,1))-1)
!!$          do jedge =0,1!min(1,nedom(auxv(pdime,2))-1)
!!$
!!$             ! Check if it is a boundary
!!$             ijedge = (/iedge,jedge/)
!!$             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
!!$                  &                   isper,auxv,flag,neigh,ledge)
!!$
!!$             ! Update edge counter
!!$             line_cnt = line_cnt + 1
!!$
!!$             ! Check case
!!$             if(flag==case) then
!!$
!!$                ! Initial object point
!!$                inipo = cnt(1)
!!$
!!$                ! Point identifier (inside the part)
!!$                ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
!!$                ijkvert(pdime) = 0
!!$                do i=1,ndime
!!$                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
!!$                end do
!!$
!!$                do i=1,gsize%npdom(pdime)-2
!!$
!!$                   ! Point identifier (inside the part)
!!$                   call kronec(pdime,ndime,kron)
!!$                   ijkpoin(1:ndime) = ijkpoin(1:ndime) + kron
!!$                   call globalid(ijkpoin,gsize%npdom,ndime,glnum)
!!$
!!$                   ! Global to local numeration (inside the part)
!!$                   npnum(glnum) = cnt(1)
!!$               
!!$                   ! Coordinates
!!$                   if(mtype==geom) call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!!$                        &                         ndime,msize,coord(:,cnt(1)))
!!$                   
!!$                   ! Boundary conditions
!!$                   if(case==inter) then
!!$                      nodes%code(:,cnt(1)) = line%code(:,line_cnt)
!!$                      nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
!!$                   end if
!!$
!!$                   ! Update point counter
!!$                   cnt(1) = cnt(1) + 1
!!$                end do
!!$
!!$                if(mtype==topo) then
!!$                   ! Elemental edges loop
!!$                   j = iedge*(tsize%nddom(auxv(pdime,1),pdime)-1) + 1
!!$                   if(ndime==3) k = jedge*(tsize%nddom(auxv(pdime,2),pdime)-1) + 1
!!$                   do i = 1,tsize%nddom(pdime,pdime)
!!$
!!$                      ! Face identifier (inside the part)
!!$                      ijkedge(pdime) = i
!!$                      ijkedge(auxv(pdime,1)) = j
!!$                      if(ndime==3) ijkedge(auxv(pdime,2)) = k
!!$                      call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)
!!$
!!$                      ! Global to local numeration (inside the part)
!!$                      npnum(aux_cnt+glnum) = cnt(1)
!!$
!!$                      ! Boundary conditions
!!$                      if(case==inter) then
!!$                         nodes%code(:,cnt(1)) = line%code(:,line_cnt)
!!$                         nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
!!$                      end if
!!$
!!$                      ! Update counter
!!$                      cnt(1) = cnt(1) + 1
!!$
!!$                      ! Update face counter
!!$                      if(ndime==2) nface = nface + 1
!!$
!!$                   end do
!!$                end if
!!$
!!$                if(case==bound) then
!!$
!!$                   ! Update object counter
!!$                   obj = obj + 1
!!$
!!$                   ! Object list
!!$                   lobjs_edge(1,obj) = 0              ! separator id, useless in the structured case
!!$                   lobjs_edge(2,obj) = inipo          ! Starting node in local numbering 
!!$                   lobjs_edge(3,obj) = cnt(1)-1       ! End node in local numbering
!!$                   lobjs_edge(4,obj) = 2*(ndime-1)    ! number of subdomains sharing the object
!!$                   lobjs_edge(5,obj) = neigh(1)       ! 1st subdomain sharing the object
!!$                   lobjs_edge(6,obj) = neigh(2)       ! 2nd subdomain sharing the object
!!$                   if(ndime==3) then
!!$                      lobjs_edge(7,obj) = neigh(3)    ! 3rd subdomain sharing the object
!!$                      lobjs_edge(8,obj) = neigh(4)    ! 4th subdomain sharing the object
!!$ !                     call psb_hsort(lobjs_edge(5:8,obj))
!!$                      call sort( 3, lobjs_edge(5:8,obj) )
!!$                   else
!!$!                      call psb_hsort(lobjs_edge(5:6,obj))
!!$                      call sort( 1, lobjs_edge(5:6,obj) )
!!$                   end if
!!$
!!$                   ! Search non-zero sharing objects and compact
!!$                   work = 0
!!$                   k = 0
!!$                   do i=5,4+2*(ndime-1)
!!$                      if(lobjs_edge(i,obj)>0) then
!!$                         k=k+1
!!$                         work(k) = lobjs_edge(i,obj)
!!$                      end if
!!$                   end do
!!$                   lobjs_edge(4,obj) = k
!!$                   lobjs_edge(5:(4+2**ndime),obj) = 0
!!$                   do j=1,k
!!$                      lobjs_edge(4+j,obj) = work(j)
!!$                   end do
!!$
!!$                   ! Boundary conditions for non-global edges
!!$                   if(ndime==3.and.k==2) then
!!$                      call aux_surf_num(pdime,aux_surf,gsurf_aux)
!!$                      do i=1,4
!!$                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0) then
!!$                            gsurf = gsurf_aux(i)
!!$                         end if
!!$                      end do
!!$                      do inode=inipo,cnt(1)-1
!!$                         nodes%code(:,inode) = surf%code(:,gsurf)
!!$                         nodes%valu(:,inode) = surf%valu(:,gsurf)
!!$                      end do
!!$                   end if
!!$                  
!!$                   ! Global edge
!!$                   ijkedge(pdime) = ijkpart(pdime)
!!$                   ijkedge(auxv(pdime,1)) = ijkpart(auxv(pdime,1)) + iedge
!!$                   if(ndime==3) then
!!$                         ijkedge(auxv(pdime,2)) = ijkpart(auxv(pdime,2)) + jedge
!!$                   end if
!!$                   do i=1,ndime-1
!!$                      if(ijkedge(auxv(pdime,i))>gsize%nedge(pdime,auxv(pdime,i))) then
!!$                         ijkedge(auxv(pdime,i)) = 1
!!$                      end if
!!$                   end do
!!$                   nedge_aux=gsize%nedge(pdime,1:ndime)
!!$                   call globalid(ijkedge,nedge_aux,ndime,gledge)
!!$                   if(pdime==ndime) then
!!$                      l2g_edge(obj) = gledge
!!$                   else if(pdime==ndime-1) then
!!$                      l2g_edge(obj) = gsize%nedget(pdime+1) + gledge
!!$                   else if(pdime==ndime-2) then
!!$                      l2g_edge(obj) = gsize%nedget(pdime+1) + gsize%nedget(pdime+2) + gledge
!!$                   end if
!!$                  
!!$                   ! Local objects
!!$                   do i=1,2
!!$                      do j=1,2
!!$                         if(ledge(i,j)==1) then
!!$                            ijklocal(pdime) = 2
!!$                            ijklocal(auxv(pdime,1)) = i + iedge
!!$                            if(ndime==3) then
!!$                               ijklocal(auxv(pdime,2)) = j + jedge
!!$                            end if
!!$                            aux_vec = (/3,3,3/)
!!$                            call globalid(ijklocal,aux_vec,ndime,glocal)
!!$                            touch(glocal) = touch(glocal) + 1
!!$                            neigh_obj(obj,glocal) = 1
!!$                         end if
!!$                         if(ndime==2) exit
!!$                      end do
!!$                   end do
!!$                end if
!!$
!!$                if(iedge.gt.min(1,abs(gsize%nedom(auxv(pdime,1))-1))) cycle
!!$                if(ndime==3) then
!!$                   if(jedge.gt.min(1,abs(gsize%nedom(auxv(pdime,2))-1))) cycle
!!$                end if
!!$
!!$                ! Element identifier (inside the part)
!!$                do i=1,ndime
!!$                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!!$                end do
!!$
!!$                do i=2,gsize%nedom(pdime)-1
!!$
!!$                   ! Element identifier (inside the part)
!!$                   call kronec(pdime,ndime,kron)
!!$                   ijkelem(1:ndime) = ijkelem(1:ndime) + kron
!!$                   call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                   ! Global to local numeration (inside the part)
!!$                   nenum(glnum) = cnt(2)
!!$
!!$                   ! Update element counter
!!$                   cnt(2) = cnt(2) + 1
!!$                end do
!!$
!!$             end if
!!$             
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$       
!!$       ! Update elemental edge counter
!!$       aux_cnt = aux_cnt + tsize%nddir(pdime)
!!$    end do
!!$
!!$    ! Deallocate auxiliar vector
!!$    deallocate(auxv)
!!$
!!$  end subroutine edge_loop
!!$
!!$ !================================================================================================
!!$  subroutine edge_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn, &
!!$       &                   lextp,lexte,lextm,mcase,nelbl)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),           intent(in)    :: ijkpart(3),isper(3),ndime,case,nb,mtype
!!$    type(geom_size),       intent(in)    :: gsize
!!$    type(topo_size),       intent(in)    :: tsize
!!$    integer(ip),           intent(inout) :: pextn_cnt(2),pextn(nb+1)
!!$    integer(ip) , optional, intent(inout) :: lextp(:),lexte(:),lextm(:)
!!$    integer(igp), optional, intent(inout) :: lextn(:)
!!$    integer(ip), optional, intent(in)    :: mcase,nelbl(3)
!!$    
!!$    integer(ip) :: ijkelem(3),neigh(2*(ndime-1)),ijedge(2),ijkvert(3),ijkneigh(3),ijkvert_aux(3),ijknode(3)
!!$    integer(ip) :: subgl(3),subgl_ijk(3),ijkpart_neigh(3),ledge(2,2),kron(ndime),lpart
!!$    integer(ip) :: i,j,pdime,iedge,jedge,kedge,medge,ielem,jelem,glnum,flag,gpart,idime,melem
!!$    integer(ip) :: inode,jnode,auxva(2**ndime),auxk(3),auxkv(3,2),kelem,nnode(3),gnode,ipos,iaux,jaux,iaux2,jaux2
!!$    integer(ip) :: ijk(3),kdime,idir,jdir,ijkedge(3),auxvd(ndime,2*(ndime-1)),auxvf(3,2)
!!$    integer(igp) :: gneigh
!!$    integer(ip), allocatable :: auxv(:,:)
!!$
!!$    ! Auxiliar vector of components
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$       ! auxva = (/1,4,2,3/)
!!$       auxva = (/1,3,2,4/)
!!$       ! auxvd(1,:) = (/5,7/)
!!$       ! auxvd(2,:) = (/8,6/)
!!$       auxvd(1,:) = (/5,6/)
!!$       auxvd(2,:) = (/7,8/)
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$       ! auxva = (/1,5,4,8,2,6,3,7/)
!!$       auxva = (/1,5,3,7,2,6,4,8/)
!!$       ! auxvd(1,:) = (/9,17,11,19/)!(/9,13,11,15/)
!!$       ! auxvd(2,:) = (/12,20,10,18/)!(/12,16,10,14/)
!!$       ! auxvd(3,:) = (/13,16,14,15/)!(/17,20,18,19/)
!!$       auxvd(1,:) = (/9,11,10,12/)!(/9,13,11,15/)
!!$       auxvd(2,:) = (/13,15,14,16/)!(/12,16,10,14/)
!!$       auxvd(3,:) = (/17,19,18,20/)!(/17,20,18,19/)
!!$    end if
!!$    auxk = (/1,2,1/)
!!$    auxkv(1,:) = (/0,0/)
!!$    auxkv(2,:) = (/0,1/)    
!!$    auxkv(3,:) = (/1,0/)
!!$    nnode = (/2,2,2/)
!!$    ! auxvf(1,:) = (/25,23/)
!!$    ! auxvf(2,:) = (/22,24/)
!!$    ! auxvf(3,:) = (/21,26/)
!!$    auxvf(1,:) = (/25,26/)
!!$    auxvf(2,:) = (/23,24/)
!!$    auxvf(3,:) = (/21,22/)
!!$    ijk = 0
!!$
!!$    ! Global numbering
!!$    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)
!!$
!!$    ! Edge nodes loop
!!$    do pdime = ndime,1,-1
!!$       do iedge = 0,1
!!$          iaux = 2*iedge-1
!!$          do jedge =0,1
!!$             jaux = 2*jedge-1
!!$
!!$             ! Check if it is a boundary
!!$             ijedge = (/iedge,jedge/)
!!$             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!!$                  &                   ndime,isper,auxv,flag,neigh,ledge)
!!$
!!$             ! Only boundary elements
!!$             if(flag==bound) then
!!$
!!$                if(case==do_count) then
!!$                
!!$                   do i=2,gsize%nedom(pdime)-1
!!$
!!$                      ! Update boundary element counter
!!$                      pextn_cnt(1) = pextn_cnt(1) + 1
!!$
!!$                      do ielem = -1,1
!!$                         do kedge = -1,1
!!$                            do medge = -1,1
!!$                               if(kedge==iaux.or.medge==jaux) then
!!$
!!$                                  iaux2 = kedge
!!$                                  jaux2 = medge
!!$                                  if(kedge.ne.iaux) iaux2 = 0
!!$                                  if(medge.ne.jaux) jaux2 = 0
!!$
!!$                                  ! Neighbour part
!!$                                  ijkpart_neigh(pdime) = ijkpart(pdime)
!!$                                  if(ndime==2) then
!!$                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
!!$                                  else
!!$                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
!!$                                  end if
!!$                                  call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!!$                                  call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!!$
!!$                                  if(gpart>0.and.gpart.ne.lpart) then
!!$
!!$                                     ! Update boundary neighbour counter
!!$                                     pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                                  end if
!!$                               end if
!!$                               if(ndime==2) exit
!!$                            end do
!!$                         end do
!!$                      end do
!!$
!!$                      ! Fill pextn
!!$                      pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1             
!!$
!!$                   end do
!!$
!!$                elseif(case==do_list) then
!!$
!!$                   ! Element identifier (inside the part)
!!$                   ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
!!$                   ijkvert(pdime) = 0
!!$                   do i=1,ndime
!!$                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!!$                   end do
!!$
!!$                   do i=2,gsize%nedom(pdime)-1
!!$
!!$                      ! Element identifier (inside the part)
!!$                      call kronec(pdime,ndime,kron)
!!$                      ijkelem(1:ndime) = ijkelem(1:ndime) + kron
!!$                      call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                      do kedge = -1,1
!!$                         do medge = -1,1
!!$                            if(kedge==iaux.or.medge==jaux) then
!!$                               
!!$                               iaux2 = kedge
!!$                               jaux2 = medge
!!$                               if(kedge.ne.iaux) iaux2 = 0
!!$                               if(medge.ne.jaux) jaux2 = 0
!!$
!!$                               ! Neighbour part
!!$                               ijkpart_neigh(pdime) = ijkpart(pdime)
!!$                               if(ndime==2) then
!!$                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
!!$                               else
!!$                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
!!$                               end if
!!$
!!$                               ! Global part numeration
!!$                               do idime=1,ndime
!!$                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
!!$                                  if(ijkpart_neigh(idime)<1) then
!!$                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
!!$                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
!!$                                     subgl(idime) = -gsize%nedom(idime)
!!$                                  end if
!!$                               end do
!!$
!!$                               ! Check part boundary
!!$                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!!$                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!!$
!!$                               if(gpart>0.and.gpart.ne.lpart) then
!!$
!!$                                  do ielem = -1,1
!!$
!!$                                     ! Update boundary neighbour counter
!!$                                     pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                                     ! Neighbour element identifier (local)
!!$                                     ijkneigh(pdime) = ijkelem(pdime) + ielem
!!$                                     if(ndime==2) then
!!$                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) + kedge
!!$                                     else
!!$                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) +(/kedge,medge/)
!!$                                     end if
!!$                                     subgl_ijk = subgl + ijkneigh
!!$                                     call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)
!!$
!!$                                     ! Fill adjacencies info
!!$                                     lextn(pextn_cnt(2)) = gneigh
!!$                                     lextp(pextn_cnt(2)) = gpart
!!$                                     lextm(pextn_cnt(2)) = 1
!!$                                     call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,lextm(pextn_cnt(2)))
!!$                                     lexte(pextn_cnt(2)) = 0
!!$
!!$                                     ! Construct bit maps
!!$                                     ! Nodes
!!$                                     do jelem=1,auxk(ielem+2)
!!$                                        do kelem=1,auxk(kedge+2)
!!$                                           do melem=1,auxk(medge+2)
!!$
!!$                                              ! ijk node identification
!!$                                              ijknode(pdime) = auxkv(ielem+2,jelem) + 1
!!$                                              if(ndime==2) then
!!$                                                 ijknode(auxv(pdime,:)) = auxkv(kedge+2,kelem) + 1
!!$                                              else
!!$                                                 ijknode(auxv(pdime,:)) = (/auxkv(kedge+2,kelem),auxkv(medge+2,melem)/) &
!!$                                                      &                   + (/1,1/)
!!$                                              end if
!!$
!!$                                              ! Global element numeration
!!$                                              call globalid(ijknode,nnode,ndime,gnode)
!!$
!!$                                              ! Fill bit map
!!$                                              lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxva(gnode)-1)
!!$
!!$                                           end do
!!$                                        end do
!!$                                     end do
!!$
!!$                                     if(mtype==topo) then
!!$                                        ! ijk element identification
!!$                                        ijk(pdime) = ielem
!!$                                        if(ndime==2) then
!!$                                           ijk(auxv(pdime,:)) = kedge
!!$                                        else
!!$                                           ijk(auxv(pdime,:)) = (/kedge,medge/)
!!$                                        end if
!!$                                        
!!$                                        ! Elements sharing 4 edges
!!$                                        if(abs(ijk(1))+abs(ijk(2))+abs(ijk(3))==ndime-2) then
!!$
!!$                                           ! Identify the sharing direction
!!$                                           do idime=1,ndime
!!$                                              if(ijk(idime).ne.0) kdime=idime
!!$                                           end do
!!$                                           
!!$                                           ! Loop over edges direction in the face
!!$                                           do idime=1,ndime-1
!!$                                              idir=auxv(kdime,idime)
!!$                                              jdir=auxv(kdime,3-idime)
!!$                                              
!!$                                              ! ijk edge identification
!!$                                              ijkedge(kdime) = (3+ijk(kdime))/2
!!$                                              do ipos=1,ndime-1
!!$                                                 ijkedge(jdir)=ipos
!!$                                                 ijknode(1:ndime-1) = ijkedge(auxv(idir,1:ndime-1))
!!$
!!$                                                 ! Global element numeration
!!$                                                 call globalid(ijknode(1:ndime-1),nnode,ndime-1,gnode)
!!$
!!$                                                 ! Fill bit map
!!$                                                 lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvd(idir,gnode)-1)
!!$
!!$                                              end do
!!$                                           end do
!!$
!!$                                           ! Faces
!!$                                           ! Fill bit map
!!$                                           lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvf(kdime,ijkedge(kdime))-1)
!!$
!!$                                        end if
!!$                                        
!!$                                        ! Elements sharing 2 edges
!!$                                        if(abs(ijk(1))+abs(ijk(2))+abs(ijk(3))==ndime-1) then
!!$                                           
!!$                                           ! Identify the sharing direction
!!$                                           do idime=1,ndime
!!$                                              if(ijk(idime)==0) kdime=idime
!!$                                           end do
!!$
!!$                                           ! ijk edge identification
!!$                                           ijkedge(auxv(kdime,1)) = (3+ijk(auxv(kdime,1)))/2
!!$                                           if(ndime==3) ijkedge(auxv(kdime,2)) = (3+ijk(auxv(kdime,2)))/2
!!$                                           ijknode(1:ndime-1) = ijkedge(auxv(kdime,1:ndime-1))
!!$
!!$                                           ! Global element numeration
!!$                                           call globalid(ijknode(1:ndime-1),nnode,ndime-1,gnode)
!!$
!!$                                           ! Fill bit map
!!$                                           lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvd(kdime,gnode)-1)
!!$
!!$                                         end if
!!$                                     end if
!!$                                  end do
!!$                               end if
!!$                            end if
!!$                            if(ndime==2) exit
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end if
!!$             end if
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    ! Deallocate auxiliar vector
!!$    deallocate(auxv)
!!$
!!$  end subroutine edge_loop_adj
!!$
!!$  !================================================================================================
!!$  subroutine corn_loop(ijkpart,ndime,mtype,gsize,tsize,msize,isper,poin,line,surf,npnum,nenum,coord, & 
!!$       &               cnt,case,lobjs_corn,obj,l2g_corn,neigh_obj,touch,nodes)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),          intent(in)    :: ijkpart(3),isper(3),ndime,case,mtype
!!$    type(geom_size),      intent(in)    :: gsize
!!$    type(topo_size),      intent(in)    :: tsize
!!$    type(mesh_size),      intent(in)    :: msize
!!$    type(fem_conditions), intent(in)    :: poin,line,surf
!!$    integer(ip),          intent(inout) :: npnum(:),nenum(:),cnt(2)
!!$    real(rp),             intent(inout) :: coord(:,:)
!!$    integer(ip),          intent(out)   :: lobjs_corn(4+2**ndime,4*(ndime-1)),obj,l2g_corn(4*(ndime-1))
!!$    integer(ip),          intent(out)   :: neigh_obj(4*(ndime-1),3**ndime),touch(3**ndime)
!!$    type(fem_conditions), intent(inout) :: nodes
!!$    
!!$    integer(ip) :: ivert,jvert,kvert,ijkpoin(3),ijkvert(3),glnum,flag,kron(ndime),neigh(2**ndime)
!!$    integer(ip) :: i,j,k,work(2**ndime),glcorn,lcorn(2,2,2),ijkelem(3),ijklocal(3),glocal,lpart
!!$    integer(ip) :: ijkcorn(3),poin_cnt,gline,gsurf,aux_vec(3)
!!$    integer(ip) :: aux_line(2*(ndime-1)*ndime,2),aux_surf(2*ndime,4)
!!$    integer(ip), allocatable :: auxv(:,:)
!!$
!!$    ! Auxiliar vector of components
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$       aux_line(:,1) = (/3,1,2,1/)
!!$       aux_line(:,2) = (/4,2,4,3/)
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$       aux_line(:,1) = (/7,5,3,1,6,5,2,1,4,3,2,1/)
!!$       aux_line(:,2) = (/8,6,4,2,8,7,4,3,8,7,6,5/)
!!$       aux_surf(:,1) = (/2,1,3,1,5,1/)
!!$       aux_surf(:,2) = (/4,3,4,2,6,2/)
!!$       aux_surf(:,3) = (/6,5,7,5,7,3/)
!!$       aux_surf(:,4) = (/8,7,8,6,8,4/)
!!$    end if
!!$
!!$    ! Objects counter
!!$    obj = 0
!!$    lobjs_corn = 0
!!$    neigh_obj = 0
!!$    touch = 0
!!$    poin_cnt = 1
!!$
!!$    ! Corner nodes loop
!!$    do ivert = 0,1
!!$       do jvert = 0,1
!!$          do kvert =0,1
!!$
!!$             ! Check if it is a boundary
!!$             ijkvert = (/ivert,jvert,kvert/)
!!$             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!!$                  &                   ndime,isper,auxv,flag,neigh,lcorn)
!!$
!!$             ! Check case
!!$             if(flag==case) then
!!$
!!$                ! Point identifier (inside the part)
!!$                do i=1,ndime
!!$                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
!!$                end do
!!$                call globalid(ijkpoin,gsize%npdom,ndime,glnum)
!!$
!!$                ! Global to local numeration (inside the part)
!!$                npnum(glnum) = cnt(1)
!!$              
!!$                ! Coordinates
!!$                if(mtype==geom) call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
!!$                     &                         ndime,msize,coord(:,cnt(1)))
!!$
!!$                ! Boundary conditions
!!$                if(case==inter) then
!!$                   nodes%code(:,cnt(1)) = poin%code(:,poin_cnt)
!!$                   nodes%valu(:,cnt(1)) = poin%valu(:,poin_cnt)
!!$                end if
!!$
!!$                ! Update point counter
!!$                cnt(1) = cnt(1) + 1
!!$
!!$                if(case==bound) then
!!$
!!$                   ! Update object counter
!!$                   obj = obj + 1
!!$
!!$                   ! Object list
!!$                   lobjs_corn(1,obj) = 0              ! separator id, useless in the structured case
!!$                   lobjs_corn(2,obj) = cnt(1)-1       ! Starting node in local numbering 
!!$                   lobjs_corn(3,obj) = cnt(1)-1       ! End node in local numbering
!!$                   lobjs_corn(4,obj) = 2**ndime       ! number of subdomains sharing the object
!!$                   lobjs_corn(5,obj) = neigh(1)       ! 1st subdomain sharing the object
!!$                   lobjs_corn(6,obj) = neigh(2)       ! 2nd subdomain sharing the object
!!$                   lobjs_corn(7,obj) = neigh(3)       ! 3rd subdomain sharing the object
!!$                   lobjs_corn(8,obj) = neigh(4)       ! 4th subdomain sharing the object
!!$                   if(ndime==3) then
!!$                      lobjs_corn(9,obj) = neigh(5)    ! 5th subdomain sharing the object
!!$                      lobjs_corn(10,obj) = neigh(6)   ! 6th subdomain sharing the object
!!$                      lobjs_corn(11,obj) = neigh(7)   ! 7th subdomain sharing the object
!!$                      lobjs_corn(12,obj) = neigh(8)   ! 8th subdomain sharing the object
!!$                      !call psb_hsort(lobjs_corn(5:12,obj))
!!$                      call sort( 7, lobjs_corn(5:12,obj) )
!!$                   else
!!$                      !call psb_hsort(lobjs_corn(5:8,obj))
!!$                      call sort( 3, lobjs_corn(5:8,obj) )
!!$                   end if
!!$
!!$                   ! Search non-zero sharing objects and compact
!!$                   work = 0
!!$                   k = 0
!!$                   do i=5,4+2**ndime 
!!$                      if(lobjs_corn(i,obj)>0) then
!!$                         k=k+1
!!$                         work(k) = lobjs_corn(i,obj)
!!$                      end if
!!$                   end do
!!$                   lobjs_corn(4,obj) = k
!!$                   lobjs_corn(5:(4+2**ndime),obj) = 0
!!$                   do j=1,k
!!$                      lobjs_corn(4+j,obj) = work(j)
!!$                   end do
!!$
!!$                   ! Boundary conditions for non-global corners
!!$                   if(k==2) then
!!$                      do i=1,2*(ndime-1)*ndime
!!$                         if(neigh(aux_line(i,1))>0.and.neigh(aux_line(i,2))>0) then
!!$                            gline = i
!!$                         end if
!!$                      end do
!!$                      nodes%code(:,cnt(1)-1) = line%code(:,gline)
!!$                      nodes%valu(:,cnt(1)-1) = line%valu(:,gline)
!!$                   else if(k>2.and.k<2**ndime.and.ndime==3) then
!!$                      do i=1,2*ndime
!!$                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0.and. &
!!$                             neigh(aux_surf(i,3))>0.and.neigh(aux_surf(i,4))>0) then
!!$                            gsurf = i
!!$                         end if
!!$                      end do
!!$                      nodes%code(:,cnt(1)-1) = surf%code(:,gsurf)
!!$                      nodes%valu(:,cnt(1)-1) = surf%valu(:,gsurf)
!!$                   end if
!!$
!!$                   ! Fix pressure at one point when periodic bc's
!!$                   if(ijkpart(1)==1.and.ijkpart(2)==1.and.ijkpart(3)==1.and. &
!!$                        & ijkvert(1)==1.and.ijkvert(2)==1.and.ijkvert(3)==1)  then
!!$                      if(ndime==2) then
!!$                         if(isper(1)==1.and.isper(2)==1) then
!!$                            !nodes%code(ndime+1,cnt(1)-1) = poin%code(ndime+1,1)
!!$                            !nodes%valu(ndime+1,cnt(1)-1) = poin%valu(ndime+1,1)
!!$                            !nodes%code(ndime+1,1) = 1!poin%code(ndime+1,1)
!!$                            !nodes%valu(ndime+1,1) = 0.0_rp!poin%valu(ndime+1,1)
!!$                         end if
!!$                      else if(ndime==3) then
!!$                         if(isper(1)==1.and.isper(2)==1.and.isper(3)==1) then
!!$                            !nodes%code(ndime+1,1) = 1!poin%code(ndime+1,1)
!!$                            !nodes%valu(ndime+1,1) = 0.0_rp!poin%valu(ndime+1,1)
!!$                         end if
!!$                      end if
!!$                   end if
!!$                                   
!!$                   ! Global corner
!!$                   ijkcorn(1) = ijkpart(1) + ivert
!!$                   ijkcorn(2) = ijkpart(2) + jvert
!!$                   ijkcorn(3) = ijkpart(3) + kvert 
!!$                   do i=1,ndime
!!$                      if(ijkcorn(i)>gsize%ncorn(i)) then
!!$                         ijkcorn(i) =1
!!$                      end if
!!$                   end do
!!$                   call globalid(ijkcorn,gsize%ncorn,ndime,glcorn)
!!$                   l2g_corn(obj) = glcorn
!!$
!!$                   ! Local objects
!!$                   do i=1,2
!!$                      do j=1,2
!!$                         do k=1,2
!!$                            if(lcorn(i,j,k)==1) then
!!$                               ijklocal = ijkvert + (/i,j,k/)
!!$                               aux_vec = (/3,3,3/)
!!$                               call globalid(ijklocal,aux_vec,ndime,glocal)
!!$                               touch(glocal) = touch(glocal) + 1
!!$                               neigh_obj(obj,glocal) = 1
!!$                            end if
!!$                            if(ndime==2) exit
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end if
!!$
!!$                if(ivert.gt.min(1,abs(gsize%nedom(1)-1))) cycle
!!$                if(jvert.gt.min(1,abs(gsize%nedom(2)-1))) cycle
!!$                if(kvert.gt.min(1,abs(gsize%nedom(3)-1))) cycle
!!$
!!$                ! Element identifier (inside the part)
!!$                do i=1,ndime
!!$                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!!$                end do
!!$                call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                ! Global to local numeration (inside the part)
!!$                nenum(glnum) = cnt(2)
!!$
!!$                ! Update element counter
!!$                cnt(2) = cnt(2) + 1
!!$
!!$             end if
!!$
!!$             ! Update corner counter
!!$             poin_cnt = poin_cnt + 1
!!$
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    ! Deallocate auxiliar vector
!!$    deallocate(auxv)
!!$
!!$  end subroutine corn_loop
!!$
!!$ !================================================================================================
!!$  subroutine corn_loop_adj(ijkpart,ndime,mtype,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn, &
!!$       &                   lextp,lexte,lextm,mcase,nelbl)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),           intent(in)     :: ijkpart(3),isper(3),ndime,case,nb,mtype
!!$    type(geom_size),       intent(in)     :: gsize
!!$    type(topo_size),       intent(in)     :: tsize
!!$    integer(ip),           intent(inout)  :: pextn_cnt(2),pextn(nb+1)
!!$    integer(ip), optional, intent(inout)  :: lextp(:),lexte(:),lextm(:)
!!$    integer(igp), optional, intent(inout) :: lextn(:)
!!$
!!$    integer(ip), optional, intent(in)    :: mcase,nelbl(3)
!!$    
!!$    integer(ip) :: ijkelem(3),neigh(2**ndime),ijkvert(3),ijkneigh(3),subgl(3),subgl_ijk(3),ijkpart_neigh(3)
!!$    integer(ip) :: i,j,k,ivert,jvert,kvert,ielem,jelem,kelem,glnum,flag,gpart,lcorn(2,2,2),idime
!!$    integer(ip) :: iaux,jaux,kaux,iaux2,jaux2,kaux2,lpart,auxk(3),auxkv(3,2),nnode(3),ipos,gnode
!!$    integer(ip) :: auxva(2**ndime),ijknode(3),ijk(3),kdime,idir,jdir,ijkedge(3),auxvd(ndime,2*(ndime-1)),auxvf(3,2)
!!$    integer(ip), allocatable :: auxv(:,:)
!!$    integer(igp) :: gneigh
!!$
!!$    ! Auxiliar vector of components
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$       ! auxva = (/1,4,2,3/)
!!$       auxva = (/1,3,2,4/)
!!$       ! auxvd(1,:) = (/5,7/)
!!$       ! auxvd(2,:) = (/8,6/)
!!$       auxvd(1,:) = (/5,6/)
!!$       auxvd(2,:) = (/7,8/)
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$       ! auxva = (/1,5,4,8,2,6,3,7/) 
!!$       auxva = (/1,5,3,7,2,6,4,8/)
!!$       ! auxvd(1,:) = (/9,17,11,19/)!(/9,13,11,15/)
!!$       ! auxvd(2,:) = (/12,20,10,18/)!(/12,16,10,14/)
!!$       ! auxvd(3,:) = (/13,16,14,15/)!(/17,20,18,19/)
!!$       auxvd(1,:) = (/9,11,10,12/)!(/9,13,11,15/)
!!$       auxvd(2,:) = (/13,15,14,16/)!(/12,16,10,14/)
!!$       auxvd(3,:) = (/17,19,18,20/)!(/17,20,18,19/)
!!$    end if       
!!$    auxk = (/1,2,1/)
!!$    auxkv(1,:) = (/0,0/)
!!$    auxkv(2,:) = (/0,1/)    
!!$    auxkv(3,:) = (/1,0/)
!!$    nnode = (/2,2,2/)
!!$    ! auxvf(1,:) = (/25,23/)
!!$    ! auxvf(2,:) = (/22,24/)
!!$    ! auxvf(3,:) = (/21,26/)
!!$    auxvf(1,:) = (/25,26/)
!!$    auxvf(2,:) = (/23,24/)
!!$    auxvf(3,:) = (/21,22/)
!!$
!!$    ! Global numbering
!!$    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)
!!$
!!$    !Corner nodes loop
!!$    do ivert = 0,1
!!$       iaux = 2*ivert-1
!!$       do jvert = 0,1
!!$          jaux = 2*jvert-1
!!$          do kvert =0,1
!!$             kaux = 2*kvert-1
!!$    
!!$             ! Check if it is a boundary
!!$             ijkvert = (/ivert,jvert,kvert/)
!!$             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
!!$                  &                   ndime,isper,auxv,flag,neigh,lcorn)
!!$
!!$             ! Only boundary elements
!!$             if(flag==bound) then
!!$
!!$                if(case==do_count) then
!!$                   
!!$                   ! Update boundary element counter
!!$                   pextn_cnt(1) = pextn_cnt(1) + 1
!!$
!!$                   do i = -1,1
!!$                      do j = -1,1
!!$                         do k = -1,1
!!$                            if(i==iaux.or.j==jaux.or.k==kaux) then
!!$
!!$                               iaux2 = i
!!$                               jaux2 = j
!!$                               kaux2 = k
!!$                               if(i.ne.iaux) iaux2 = 0
!!$                               if(j.ne.jaux) jaux2 = 0
!!$                               if(k.ne.kaux) kaux2 = 0
!!$                               
!!$                               ! Neighbour part
!!$                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)
!!$                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!!$                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!!$                               
!!$                               if(gpart>0.and.gpart.ne.lpart) then
!!$                                  
!!$                                  ! Update boundary neighbour counter
!!$                                  pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                               end if
!!$
!!$                            end if
!!$                            if(ndime==2) exit
!!$                         end do
!!$                      end do
!!$                   end do
!!$
!!$                   ! Fill pextn
!!$                   pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1
!!$
!!$                elseif(case==do_list) then
!!$
!!$                   ! Element identifier (inside the part)
!!$                   do i=1,ndime
!!$                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
!!$                   end do
!!$
!!$                   ! Element identifier (inside the part)
!!$                   call globalid(ijkelem,gsize%nedom,ndime,glnum)
!!$
!!$                   do i = -1,1
!!$                      do j = -1,1
!!$                         do k = -1,1
!!$                            if(i==iaux.or.j==jaux.or.k==kaux) then
!!$
!!$                               iaux2 = i
!!$                               jaux2 = j
!!$                               kaux2 = k
!!$                               if(i.ne.iaux) iaux2 = 0
!!$                               if(j.ne.jaux) jaux2 = 0
!!$                               if(k.ne.kaux) kaux2 = 0
!!$                               
!!$                               ! Neighbour part
!!$                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)
!!$
!!$                               ! Global part numeration
!!$                               do idime=1,ndime
!!$                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
!!$                                  if(ijkpart_neigh(idime)<1) then
!!$                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
!!$                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
!!$                                     subgl(idime) = -gsize%nedom(idime)
!!$                                  end if
!!$                               end do
!!$
!!$                               ! Check boundary partition
!!$                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
!!$                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
!!$                               
!!$                               if(gpart>0.and.gpart.ne.lpart) then
!!$
                                  ! Global part numeration
                                  ! do idime=1,ndime
                                  !    subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
                                  ! end do
!!$
!!$                                  ! Update boundary neighbour counter
!!$                                  pextn_cnt(2) = pextn_cnt(2) + 1
!!$
!!$                                  ! Neighbour element identifier (local)
!!$                                  ijkneigh = ijkelem + (/i,j,k/)
!!$                                  subgl_ijk = subgl + ijkneigh
!!$                                  call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)
!!$
!!$                                  ! Fill adjacencies info
!!$                                  lextn(pextn_cnt(2)) = gneigh
!!$                                  lextp(pextn_cnt(2)) = gpart
!!$                                  lextm(pextn_cnt(2)) = 1
!!$                                  call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,lextm(pextn_cnt(2)))
!!$                                  lexte(pextn_cnt(2)) = 0
!!$
!!$                                  ! Construct bit maps
!!$                                  ! Nodes
!!$                                  do ielem=1,auxk(i+2)
!!$                                     do jelem=1,auxk(j+2)
!!$                                        do kelem=1,auxk(k+2)
!!$
!!$                                           ! ijk node identification
!!$                                           ijknode = (/auxkv(i+2,ielem),auxkv(j+2,jelem),auxkv(k+2,kelem)/) + (/1,1,1/)
!!$
!!$                                           ! Global element numeration
!!$                                           call globalid(ijknode,nnode,ndime,gnode)
!!$
!!$                                           ! Fill bit map
!!$                                           lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxva(gnode)-1)
!!$                                           
!!$                                        end do
!!$                                     end do
!!$                                  end do
!!$                                  
!!$                                  if(mtype==topo) then
!!$                                     ! ijk element identification
!!$                                     ijk = (/i,j,k/)
!!$                                     
!!$                                     ! Elements sharing 4 edges
!!$                                     if(abs(ijk(1))+abs(ijk(2))+abs(ijk(3))==1) then
!!$
!!$                                        ! Identify the sharing direction
!!$                                        do idime=1,ndime
!!$                                           if(ijk(idime).ne.0) kdime=idime
!!$                                        end do
!!$
!!$                                        ! Loop over edges direction in the face
!!$                                        do idime=1,ndime-1
!!$                                           idir=auxv(kdime,idime)
!!$                                           jdir=auxv(kdime,3-idime)
!!$
!!$                                           ! ijk edge identification
!!$                                           ijkedge(kdime) = (3+ijk(kdime))/2
!!$                                           do ipos=1,ndime-1
!!$                                              ijkedge(jdir)=ipos
!!$                                              ijknode(1:ndime-1) = ijkedge(auxv(idir,1:ndime-1))
!!$
!!$                                              ! Global element numeration
!!$                                              call globalid(ijknode(1:ndime-1),nnode,ndime-1,gnode)
!!$
!!$                                              ! Fill bit map
!!$                                              lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvd(idir,gnode)-1)
!!$
!!$                                           end do
!!$                                        end do
!!$
!!$                                        ! Faces
!!$                                        ! Fill bit map
!!$                                        lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvf(kdime,ijkedge(kdime))-1)
!!$
!!$                                     end if
!!$
!!$                                     ! Elements sharing 2 edges
!!$                                     if(abs(ijk(1))+abs(ijk(2))+abs(ijk(3))==2) then
!!$
!!$                                        ! Identify the sharing direction
!!$                                        do idime=1,ndime
!!$                                           if(ijk(idime)==0) kdime=idime
!!$                                        end do
!!$
!!$                                        ! ijk edge identification
!!$                                        ijkedge(auxv(kdime,1)) = (3+ijk(auxv(kdime,1)))/2
!!$                                        if(ndime==3) ijkedge(auxv(kdime,2)) = (3+ijk(auxv(kdime,2)))/2
!!$                                        ijknode(1:ndime-1) = ijkedge(auxv(kdime,1:ndime-1))
!!$
!!$                                        ! Global element numeration
!!$                                        call globalid(ijknode(1:ndime-1),nnode,ndime-1,gnode)
!!$
!!$                                        ! Fill bit map
!!$                                        lexte(pextn_cnt(2)) = ibset(lexte(pextn_cnt(2)),auxvd(kdime,gnode)-1)
!!$
!!$                                     end if
!!$                                  end if
!!$                               end if
!!$                            end if
!!$                            if(ndime==2) exit
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end if
!!$             end if
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    ! Deallocate auxiliar vector
!!$    deallocate(auxv)
!!$
!!$  end subroutine corn_loop_adj
!!$
!!$  !================================================================================================
!!$  subroutine inte_object(ijkpart,nsckt,npsoc,ndime,cnt,lobjs_int)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)    :: ijkpart(3),nsckt(3),npsoc(3),ndime,cnt(2)
!!$    integer(ip), intent(out)   :: lobjs_int(5)
!!$
!!$    integer(ip) :: lpart
!!$
!!$    ! Part id
!!$    call globalid_2l(ijkpart,nsckt,npsoc,ndime,lpart)
!!$
!!$    ! Interior objects
!!$    lobjs_int(1) = 0        ! separator id, useless in the structured case
!!$    lobjs_int(2) = 1        ! Starting node in local numbering 
!!$    lobjs_int(3) = cnt(1)-1 ! End node in local numbering
!!$    lobjs_int(4) = 1        ! number of subdomains sharing the object
!!$    lobjs_int(5) = lpart    ! 1st subdomain sharing the object
!!$
!!$  end subroutine inte_object
!!$
!!$  !================================================================================================
!!$  subroutine coord_ijk(ijkpoin,ijkpart,npdir,nedom,nedir,ndime,msize,coord)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),     intent(in)    :: ijkpoin(3),ijkpart(3),npdir(3),nedom(3),nedir(3)
!!$    integer(ip),     intent(in)    :: ndime
!!$    type(mesh_size), intent(in)    :: msize
!!$    real(rp),        intent(inout) :: coord(:)
!!$
!!$    integer(ip) :: i,ntdis(3),pdegr,nebl(3)
!!$    real(rp)    :: leng(3),coor0(3),stret(3),old_stret,istret,lengbl(3)
!!$    integer(ip) :: newijkpoin,newnedir,newnedom,info
!!$    real(rp)    :: newleng
!!$
!!$    ! Unpack mesh_size
!!$    pdegr   = msize%pdegr
!!$    leng(1) = msize%xleng
!!$    leng(2) = msize%yleng
!!$    leng(3) = msize%zleng
!!$    ntdis(1) = msize%ntdix
!!$    ntdis(2) = msize%ntdiy
!!$    ntdis(3) = msize%ntdiz
!!$    stret(1) = msize%xstret
!!$    stret(2) = msize%ystret
!!$    stret(3) = msize%zstret
!!$    old_stret = msize%stret
!!$    nebl(1) = msize%neblx
!!$    nebl(2) = msize%nebly
!!$    nebl(3) = msize%neblz
!!$    lengbl(1) = msize%xlengbl
!!$    lengbl(2) = msize%ylengbl
!!$    lengbl(3) = msize%zlengbl
!!$
!!$    ! Set origin of coordinates
!!$    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0
!!$
!!$    ! Set origin of coordinates
!!$    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0
!!$
!!$    do i=1,ndime 
!!$       assert ( ntdis(i).ne.1 ) ! This case is not implemented
!!$       if(ntdis(i)==0) then
!!$          coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)
!!$       elseif(ntdis(i)==2) then
!!$          if(old_stret>0.0_rp) then
!!$             istret=old_stret
!!$          else
!!$             istret=stret(i)
!!$          end if
!!$          coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)))/tanh(istret)
!!$       elseif(ntdis(i)==3.or.ntdis(i)==4) then ! boundary layers (Hunt's case: solid uniform + fluid unif==3 or tanh==4)
!!$          if(old_stret>0.0_rp) then
!!$             istret=old_stret
!!$          else
!!$             istret=stret(i)
!!$          end if
!!$          
          !TO DO
          !if((nebl(i)>0).and.(nedom(i)<nebl(i))) then
          !   write(*,*) 'boundary layer cannot be splitted over several sbds, dir ',i,'nedom',nedom,'nebl',nebl
          !   call mpi_barrier (1, info)
          !   call mpi_finalize( info )
          !   stop
          !end if
!!$
!!$          !first boundary layer
!!$          if((ijkpart(i)==1).and.(ijkpoin(i)<=nebl(i))) then 
!!$             ! coord(i) = -leng(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
!!$             if(ntdis(i)==3) then
!!$                coord(i) = coor0(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
!!$             elseif(ntdis(i)==4) then
!!$                coord(i) = coor0(i) - leng(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
!!$             end if
!!$          !latest boundary layer
!!$          elseif((ijkpart(i)==npdir(i)).and.(ijkpoin(i)>(nedom(i)-nebl(i)+1))) then 
!!$             newijkpoin=ijkpoin(i)-(nedom(i)-nebl(i))
!!$             ! coord(i) = (leng(i)-lengbl(i)) + lengbl(i)*(newijkpoin-1)/nebl(i)
!!$             coord(i) = coor0(i) + leng(i) + lengbl(i)*(newijkpoin-1)/nebl(i)
!!$          else !core region

             !newleng=leng(i)-lengbl(i)
             !newijkpoin=ijkpoin(i)+(ijkpart(i)-1)*nedom(i)-nebl(i)
             !newnedir=nedir(i)-2*nebl(i)
             !newnedom=nedom(i)
             !coord(i) = coor0(i)-newleng*tanh(istret*(1.0_rp-2.0_rp*(newijkpoin-1)/newnedir))/tanh(istret)

!!$             if(ntdis(i)==3) then
!!$                coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))
!!$             elseif(ntdis(i)==4) then
!!$                coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+&
!!$                     (ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))))/tanh(istret)
!!$             end if
!!$          end if
!!$       end if
!!$    end do
!!$   
!!$  end subroutine coord_ijk
!!$
!!$  !================================================================================================
!!$  subroutine check_part_boundary(ijk,npdir,ndime,isper)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)    :: npdir(3),ndime,isper(3)
!!$    integer(ip), intent(inout) :: ijk(3)
!!$
!!$    integer(ip) :: i
!!$    
    !if(ndime==2) then
    !   if(isper(1)==0.and.isper(2)==0) then
    !      ! Not periodic mesh
    !      do i=1,ndime
    !         if(ijk(i)<1.or.ijk(i)>npdir(i)) then
    !            ijk = 0*ijk
    !         end if
    !      end do
          
    !   else if(isper(1)==1.and.isper(2)==1) then
    !      ! Periodic mesh with respect to all directions
    !      do i=1,ndime
    !         if(ijk(i)<1) then
    !            ijk(i) = npdir(i)
    !         else if(ijk(i)>npdir(i)) then
    !            ijk(i) = 1!npdir(i)
    !         end if
    !      end do
    !   end if
    !else if(ndime==3) then
!!$
!!$       ! Periodic mesh with respect to all directions
!!$       do i=1,ndime
!!$          if(isper(i)==1) then
!!$             if(ijk(i)<1) then
!!$                ijk(i) = npdir(i)
!!$             else if(ijk(i)>npdir(i)) then
!!$                ijk(i) = 1!npdir(i)
!!$             end if
!!$          else
!!$             if(ijk(i)<1.or.ijk(i)>npdir(i)) then
!!$                ijk = 0*ijk
!!$             end if
!!$          end if
!!$       end do
!!$
       !if(isper(1)==0.and.isper(2)==0.and.isper(3)==0) then
       !   ! Not periodic mesh
       !   do i=1,ndime
       !      if(ijk(i)<1.or.ijk(i)>npdir(i)) then
       !         ijk = 0*ijk
       !      end if
       !   end do  
       !else if(isper(1)==1.and.isper(2)==1.and.isper(3)==1) then
       !   ! Periodic mesh with respect to all directions
       !   do i=1,ndime
       !      if(ijk(i)<1) then
       !         ijk(i) = npdir(i)
       !      else if(ijk(i)>npdir(i)) then
       !         ijk(i) = 1!npdir(i)
       !      end if
       !   end do
       !else if(isper(1)==1.and.isper(2)==0.and.isper(3)==1) then
       !   ! Non periodic condition on y direction
       !   if(ijk(1)<1) then
       !      ijk(1) = npdir(1)
       !   else if(ijk(1)>npdir(1)) then
       !      ijk(1) = 1!npdir(i)
       !   end if
       !   if(ijk(3)<1) then
       !      ijk(3) = npdir(3)
       !   else if(ijk(3)>npdir(3)) then
       !      ijk(3) = 1!npdir(i)
       !   end if
       !   if(ijk(2)<1.or.ijk(2)>npdir(2)) then
       !      ijk = 0*ijk
       !   end if
       !else if(isper(1)==1.and.isper(2)==0.and.isper(3)==0) then 
       !   ! Periodic mesh only with respect to x direction
       !   if(ijk(1)<1) then
       !      ijk(1) = npdir(1)
       !   else if(ijk(1)>npdir(1)) then
       !      ijk(1) = 1!npdir(i)
       !   end if
       !   do i=2,ndime
       !      if(ijk(i)<1.or.ijk(i)>npdir(i)) then
       !         ijk = 0*ijk 
       !      end if
       !   end do
       !end if
    !end if
!!$
!!$  end subroutine check_part_boundary
!!$
!!$  !================================================================================================
!!$  subroutine check_face_boundary(ijkpart,pdime,iface,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
!!$       &                         neigh,lface)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),pdime,iface,npdir(3),ndime,isper(3),auxv(3,2)
!!$    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!!$    integer(ip), intent(out) :: flag,neigh(2),lface(2)
!!$
!!$    integer(ip)              :: i,ijkneig(3),newpo,cnt
!!$
!!$    cnt=1
!!$    flag = -1
!!$    lface = 0
!!$
!!$    do i=-1,0
!!$
!!$       ! Auxiliar variables
!!$       newpo = iface + i
!!$
!!$       ! Check neighbors
!!$       call neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)
!!$
!!$       ! Check part boundary
!!$       call check_part_boundary(ijkneig,npdir,ndime,isper)
!!$
!!$       ! Global numbering
!!$       call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))
!!$       
!!$       ! Set flag
!!$       if(neigh(cnt)>0) then
!!$          flag = min(flag+1,1)
!!$          lface(i+2) = 1 
!!$       end if
!!$
!!$       ! Update counter
!!$       cnt = cnt + 1
!!$
!!$    end do
!!$
!!$  end subroutine check_face_boundary
!!$
!!$  !================================================================================================
!!$  subroutine check_edge_boundary(ijkpart,pdime,ijedge,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
!!$       &                         neigh,ledge)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),pdime,ijedge(2),npdir(3),ndime,isper(3),auxv(:,:)
!!$    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!!$    integer(ip), intent(out) :: flag,neigh(2*(ndime-1)),ledge(2,2)
!!$
!!$    integer(ip)              :: i,j,ijkneig(3),newpo(2),cnt
!!$
!!$    cnt=1
!!$    flag = -1
!!$    ledge = 0
!!$
!!$    do i=-1,0
!!$       do j=-1,0
!!$
!!$          ! Auxiliar variables
!!$          newpo = ijedge + (/i,j/)
!!$
!!$          ! Check neighbors
!!$          call neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)
!!$
!!$          ! Check part boundary
!!$          call check_part_boundary(ijkneig,npdir,ndime,isper)
!!$
!!$          ! Global numbering
!!$          call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))
!!$
!!$          ! Set flag
!!$          if(neigh(cnt)>0) then
!!$             flag = min(flag+1,1)
!!$             ledge(i+2,j+2) = 1
!!$          end if
!!$
!!$          ! Update counter
!!$          cnt = cnt + 1
!!$
!!$          if(ndime==2) exit
!!$       end do
!!$    end do
!!$
!!$  end subroutine check_edge_boundary
!!$
!!$  !================================================================================================
!!$  subroutine check_corn_boundary(ijkpart,ijkvert,npdir,nsckt,npsoc,ndime,isper,auxv,flag,neigh, &
!!$       &                         lcorn)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),ijkvert(3),npdir(3),ndime,isper(3),auxv(:,:)
!!$    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
!!$    integer(ip), intent(out) :: flag,neigh(2**ndime),lcorn(2,2,2)
!!$
!!$    integer(ip)              :: i,j,k,ijkneig(3),newpo(3),cnt
!!$
!!$    cnt = 1
!!$    flag = -1
!!$    lcorn = 0
!!$
!!$    do i=-1,0
!!$       do j=-1,0
!!$          do k=-1,0
!!$
!!$             ! Auxiliar variables
!!$             newpo = ijkvert + (/i,j,k/)
!!$
!!$             ! Check neighbors
!!$             call neigh_corn(newpo,ijkpart,ijkneig)
!!$
!!$             ! Check part boundary
!!$             call check_part_boundary(ijkneig,npdir,ndime,isper)
!!$
!!$             ! Global numbering
!!$             call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))
!!$
!!$             ! Set flag
!!$             if(neigh(cnt)>0) then
!!$                flag = min(flag+1,1)
!!$                lcorn(i+2,j+2,k+2) = 1
!!$             end if
!!$             ! Update counter
!!$             cnt = cnt + 1
!!$
!!$             if(ndime==2) exit
!!$          end do
!!$       end do
!!$    end do
!!$
!!$  end subroutine check_corn_boundary
!!$
!!$  !================================================================================================
!!$  subroutine neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),newpo,pdime,auxv(3,2)
!!$    integer(ip), intent(out) :: ijkneig(3)
!!$
!!$    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:))
!!$    ijkneig(pdime) = ijkpart(pdime) + newpo
!!$
!!$  end subroutine neigh_face
!!$
!!$  !================================================================================================
!!$  subroutine neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),newpo(2),pdime,auxv(:,:)
!!$    integer(ip), intent(out) :: ijkneig(3)
!!$
!!$    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + newpo(1:size(auxv(pdime,:)))
!!$    ijkneig(pdime) = ijkpart(pdime)
!!$
!!$  end subroutine neigh_edge
!!$
!!$  !================================================================================================
!!$  subroutine neigh_corn(newpo,ijkpart,ijkneig)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijkpart(3),newpo(3)
!!$    integer(ip), intent(out) :: ijkneig(3)
!!$
!!$    ijkneig = ijkpart + newpo
!!$
!!$  end subroutine neigh_corn
!!$
!!$  !================================================================================================
!!$  subroutine generic_l2g(subgl,npnum,nenum,ndime,mtype,gsize,tsize,msize,l2ge,l2gp,lnods,mater)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip),     intent(in)  :: subgl(ndime),ndime,mtype
!!$    integer(ip),     intent(in)  :: npnum(:),nenum(:)
!!$    type(geom_size), intent(in)  :: gsize
!!$    type(topo_size), intent(in)  :: tsize
!!$    type(mesh_size), intent(in)  :: msize
!!$    integer(ip),     intent(out) :: l2gp(:)
!!$    integer(igp),     intent(out) :: l2ge(:)
!!$    integer(ip),     intent(out) :: lnods(:)
!!$    type(fem_materials), optional, intent(inout) :: mater
!!$
!!$    integer(ip)               :: i,j,k,m,n,l,ijk(3),ijklnods(3),num,count,auxva((msize%pdegr+1)**ndime)
!!$    integer(ip)               :: ne_aux(3),np_aux(3),pdegr,mcase,auxnum,nedir_l2g(ndime),nobje,aux_glb
!!$    integer(ip), allocatable  :: po_l2g(:),auxv(:,:)
!!$    integer(igp), allocatable :: el_l2g(:)
!!$    integer(ip)              :: subgl_aux(3),subgl_ijk(3),nelbl(3)
!!$    integer(ip)              :: cnt,ipoin,jpoin,ielem,jelem,pdime,iedge,jedge,iface,nddomk
!!$    integer(ip)              :: auxvc(2**ndime),auxvf(2*ndime),auxvd(2*(ndime-1)*ndime),aux_cnt
!!$
!!$    ! Allocate
!!$    call memalloc(gsize%nedomt,el_l2g,__FILE__,__LINE__)
!!$    if(mtype==geom) call memalloc(gsize%npdomt,po_l2g,__FILE__,__LINE__)
!!$    if(mtype==topo) call memalloc(tsize%notot ,po_l2g,__FILE__,__LINE__)
!!$
!!$    ! Auxiliar vectro
!!$    pdegr = msize%pdegr
!!$    check( pdegr == 1 )
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$       if(pdegr==1) then
!!$          !auxva = (/1,4,2,3/)
!!$          auxva = (/1,3,2,4/)
!!$       else if(pdegr==2) then
!!$          auxva = (/1,8,4,5,9,7,2,6,3/)
!!$       else if(pdegr==3) then
!!$          auxva = (/1,12,11,4,5,13,16,10,6,14,15,9,2,7,8,3/)
!!$       end if
!!$       ne_aux(1)=1
!!$       ne_aux(2)=gsize%nedom(1)
!!$       ne_aux(3)=gsize%nedom(2)
!!$       np_aux(1)=1
!!$       np_aux(2)=gsize%npdom(1)
!!$       np_aux(3)=gsize%npdom(2)
!!$       subgl_aux(1:2) = subgl
!!$       subgl_aux(3) = 0
!!$       !auxvc = (/1,4,2,3/)
!!$       auxvc = (/1,3,2,4/)
!!$       !auxvd = (/8,6,5,7/)
!!$       auxvd = (/7,8,5,6/)
!!$       nobje = 8
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$       if(pdegr==1) then
!!$          !auxva = (/4,1,8,5,3,2,7,6/)
!!$          ! auxva = (/1,5,4,8,2,6,3,7/)
!!$          auxva = (/1,5,3,7,2,6,4,8/)
!!$       else if(pdegr==2) then
!!$          !auxva = (/5,20,8,13,25,16,1,12,4,17,26,19,22,27,24,9,21,11,6,18,7,14,23,15,2,10,3/)
!!$          auxva = (/1,13,5,12,25,20,4,16,8,9,22,17,21,27,26,11,24,19,2,14,6,10,23,18,3,15,7/)
!!$          !auxva = (/7,15,3,18,23,10,6,14,2,19,24,11,26,27,21,17,22,9,8,16,4,20,25,12,5,13,1/)
!!$       else if(pdegr==3) then
!!$          auxva = (/1,17,21,5,16,44,52,32,15,43,51,31,4,20,24,8,9,37,45,25,33,57,61,53,36,60,64,56,14,42, &
!!$               &    50,30,10,38,46,26,34,58,62,54,35,59,63,55,13,41,49,29,2,18,22,6,11,39,47,27,12,40,48, &
!!$               &    28,3,19,23,7/)
!!$       end if
!!$       ne_aux=gsize%nedom
!!$       np_aux=gsize%npdom
!!$       subgl_aux = subgl
!!$       !auxvc = (/1,5,4,8,2,6,3,7/)
!!$       auxvc = (/1,5,3,7,2,6,4,8/)
!!$       !auxvd = (/13,16,14,15,12,20,10,18,9,17,11,19/)!(/17,20,18,19,12,16,10,14,9,13,11,15/)
!!$       auxvd = (/17,19,18,20,13,15,14,16,9,11,10,12/)
!!$       !auxvf = (/21,26,22,24,25,23/)
!!$       auxvf = (/21,22,23,24,25,26/)
!!$       nobje = 26
!!$    end if
!!$    mcase = msize%mater
!!$    nelbl(1) = msize%neblx
!!$    nelbl(2) = msize%nebly
!!$    nelbl(3) = msize%neblz
!!$
!!$    ! Elements
!!$    do i=1,ne_aux(1)
!!$       do j=1,ne_aux(2)
!!$          do k=1,ne_aux(3)
!!$
!!$             ! Identifier
!!$             if(ndime==2) then
!!$                ijk=(/j,k,i/)
!!$             else if(ndime==3) then
!!$                ijk=(/i,j,k/)
!!$             end if
!!$             call globalid(ijk,gsize%nedom,ndime,num)
!!$             
!!$             ! Local to global vector
!!$             subgl_ijk = subgl_aux+ijk
!!$             call globalid(subgl_ijk,gsize%nedir,ndime,el_l2g(num))
!!$
!!$             ! Fill material
!!$             if(mcase>0 .and. present(mater)) then
!!$                call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,mater%list(nenum(num)))
!!$             end if
!!$           
!!$             count = 1
!!$             if(mtype==geom) then
!!$                do m=0,pdegr
!!$                   do n=0,pdegr
!!$                      do l=0,pdegr
!!$
!!$                         ! Identifier (point in the element)
!!$                         if(ndime==2) then
!!$                            ijklnods=(/((j-1)*pdegr+1),((k-1)*pdegr+1),((i-1)*pdegr+1)/) + (/m,n,l/) !ijk + (/m,n,l/)
!!$                         else if(ndime==3) then
!!$                            ijklnods=(/((i-1)*pdegr+1),((j-1)*pdegr+1),((k-1)*pdegr+1)/) + (/m,n,l/) !ijk + (/m,n,l/)
!!$                         end if
!!$                         call globalid(ijklnods,gsize%npdom,ndime,auxnum)
!!$
!!$                         ! Generate lnods
!!$                         lnods((pdegr+1)**ndime*(nenum(num)-1)+auxva(count)) = npnum(auxnum)
!!$
!!$                         ! Update counter
!!$                         count = count +1
!!$
!!$                         if(ndime==2) exit
!!$                      end do
!!$                   end do
!!$                end do
!!$             elseif(mtype==topo) then
!!$                ! Elemental corners
!!$                do m=0,1
!!$                   do n=0,1
!!$                      do l=0,1
!!$
!!$                         ! Corner identifier in the element
!!$                         if(ndime==2) then
!!$                            ijklnods= ijk + (/m,n,l/)
!!$                         else if(ndime==3) then
!!$                            ijklnods= ijk + (/m,n,l/)
!!$                         end if
!!$                         call globalid(ijklnods,gsize%npdom,ndime,auxnum)
!!$
!!$                         ! Generate lnods
!!$                         lnods(nobje*(nenum(num)-1)+auxvc(count)) = npnum(auxnum)
!!$
!!$                         ! Update counter
!!$                         count = count +1
!!$
!!$                         if(ndime==2) exit
!!$                      end do
!!$                   end do
!!$                end do
!!$                
!!$                ! Elemental edges
!!$                count = 1
!!$                aux_cnt = tsize%nctot
!!$                do pdime=ndime,1,-1
!!$                   do iedge=0,1
!!$                      do jedge=0,1
!!$
!!$                         ! Edge identifier in the element
!!$                         ijklnods(pdime) = ijk(pdime)
!!$                         ijklnods(auxv(pdime,1)) = ijk(auxv(pdime,1)) + iedge
!!$                         if(ndime==3) ijklnods(auxv(pdime,2)) = ijk(auxv(pdime,2)) + jedge
!!$                         call globalid(ijklnods,tsize%nddom(:,pdime),ndime,auxnum)
!!$
!!$                         ! Generate lnods
!!$                         lnods(nobje*(nenum(num)-1)+auxvd(count)) = npnum(aux_cnt + auxnum)
!!$
!!$                         ! Update counter
!!$                         count = count +1
!!$
!!$                         if(ndime==2) exit
!!$                      end do
!!$                   end do
!!$                   aux_cnt = aux_cnt + tsize%nddir(pdime)
!!$                end do
!!$
!!$                ! Elemental faces
!!$                count = 1
!!$                if(ndime==3) then
!!$                   aux_cnt = tsize%nctot + tsize%ndtot
!!$                   do pdime=ndime,1,-1
!!$                      do iface = 0,1
!!$
!!$                         ! Face identifier in the element
!!$                         ijklnods(pdime) = ijk(pdime) + iface
!!$                         ijklnods(auxv(pdime,:)) = ijk(auxv(pdime,:))
!!$                         call globalid(ijklnods,tsize%nfdom(:,pdime),ndime,auxnum)
!!$
!!$                         ! Generate lnods
!!$                         lnods(nobje*(nenum(num)-1)+auxvf(count)) = npnum(aux_cnt + auxnum)
!!$                         
!!$                         ! Update counter
!!$                         count = count + 1
!!$
!!$                      end do
!!$                      aux_cnt = aux_cnt + tsize%nfdir(pdime)
!!$                   end do
!!$                end if
!!$             end if
!!$
!!$          end do
!!$       end do
!!$       if(ndime==2) exit
!!$    end do
!!$
!!$    ! Construct l2g emap (ordered by objects)
!!$    do ielem=1,gsize%nedomt
!!$       l2ge(nenum(ielem)) = el_l2g(ielem)
!!$    end do
!!$
!!$    ! Nodes
!!$    do i=1,np_aux(1)
!!$       do j=1,np_aux(2)
!!$          do k=1,np_aux(3)
!!$
!!$             ! Identifier
!!$             if(ndime==2) then
!!$                ijk=(/j,k,i/)
!!$             else if(ndime==3) then
!!$                ijk=(/i,j,k/)
!!$             end if
!!$             call globalid(ijk,gsize%npdom,ndime,num)
!!$
!!$             ! Local to global vector
!!$             nedir_l2g = gsize%nedir(1:ndime)*pdegr + 1
!!$             subgl_ijk = subgl_aux*pdegr+ijk
!!$             call globalid(subgl_ijk,nedir_l2g,ndime,po_l2g(num))
!!$    
!!$          end do
!!$       end do
!!$       if(ndime==2) exit
!!$    end do
!!$
!!$    if(mtype==topo) then
!!$       ! Edges
!!$       aux_cnt = tsize%nctot
!!$       aux_glb = tsize%ncglb
!!$       do pdime=ndime,1,-1
!!$          if(ndime==2) then
!!$             nddomk=1
!!$          elseif(ndime==3) then
!!$             nddomk=tsize%nddom(auxv(pdime,2),pdime)
!!$          end if
!!$          do i=1,tsize%nddom(pdime,pdime)
!!$             do j=1,tsize%nddom(auxv(pdime,1),pdime)
!!$                do k=1,nddomk
!!$
!!$                   ! Identifier
!!$                   ijk(pdime) = i
!!$                   ijk(auxv(pdime,1)) = j
!!$                   if(ndime==3) ijk(auxv(pdime,2)) = k
!!$                   call globalid(ijk,tsize%nddom(:,pdime),ndime,num)
!!$
!!$                   ! Local to global vector
!!$                   subgl_ijk = subgl_aux + ijk
!!$                   call globalid(subgl_ijk,tsize%ndglb(:,pdime),ndime,po_l2g(num+aux_cnt))
!!$                   po_l2g(num+aux_cnt) = po_l2g(num+aux_cnt) + aux_glb
!!$
!!$                end do
!!$             end do
!!$          end do
!!$          aux_cnt = aux_cnt + tsize%nddir(pdime)
!!$          aux_glb = aux_glb + tsize%ndsum(pdime)
!!$       end do
!!$
!!$       ! Faces
!!$       if(ndime==3) then
!!$          aux_cnt = tsize%nctot + tsize%ndtot
!!$          do pdime=ndime,1,-1
!!$             do i=1,tsize%nfdom(pdime,pdime)
!!$                do j=1,tsize%nfdom(auxv(pdime,1),pdime)
!!$                   do k=1,tsize%nfdom(auxv(pdime,2),pdime)
!!$
!!$                      ! Identifier
!!$                      ijk(pdime) = i
!!$                      ijk(auxv(pdime,:)) = (/j,k/)
!!$                      call globalid(ijk,tsize%nfdom(:,pdime),ndime,num)
!!$
!!$                      ! Local to global vector
!!$                      subgl_ijk = subgl_aux + ijk
!!$                      call globalid(subgl_ijk,tsize%nfglb(:,pdime),ndime,po_l2g(num+aux_cnt))
!!$                      po_l2g(num+aux_cnt) = po_l2g(num+aux_cnt) + aux_glb
!!$
!!$                   end do
!!$                end do
!!$             end do
!!$             aux_cnt = aux_cnt + tsize%nfdir(pdime)
!!$             aux_glb = aux_glb + tsize%nfsum(pdime)
!!$          end do
!!$       end if
!!$    end if
!!$
!!$    ! Construct l2g nmap (ordered by objects)
!!$    ! Elemental corners
!!$    do ipoin=1,gsize%npdomt
!!$       l2gp(npnum(ipoin)) = po_l2g(ipoin)
!!$    end do
!!$    if(mtype==topo) then
!!$       ! Elemental edges
!!$       aux_cnt = tsize%nctot
!!$       do pdime=ndime,1,-1
!!$          do i=1,tsize%nddir(pdime)
!!$             l2gp(npnum(i+aux_cnt)) = po_l2g(i+aux_cnt)
!!$          end do
!!$          aux_cnt = aux_cnt + tsize%nddir(pdime)
!!$       end do
!!$       ! Elemental faces
!!$       aux_cnt = tsize%nctot + tsize%ndtot
!!$       do pdime=ndime,1,-1
!!$          do i=1,tsize%nfdir(pdime)
!!$             l2gp(npnum(i+aux_cnt)) = po_l2g(i+aux_cnt)
!!$          end do
!!$          aux_cnt = aux_cnt + tsize%nfdir(pdime)
!!$       end do
!!$    end if        
!!$
!!$    ! Deallocate
!!$    call memfree(el_l2g,__FILE__,__LINE__)
!!$    call memfree(po_l2g,__FILE__,__LINE__)
!!$    deallocate(auxv)
!!$
!!$  end subroutine generic_l2g
!!$
!!$  !================================================================================================
!!$  subroutine kronec(i,n,kron)
!!$    !-----------------------------------------------------------------------
!!$    ! 
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)  :: i,n
!!$    integer(ip), intent(out) :: kron(n)
!!$
!!$    kron = 0
!!$    kron(i) = 1
!!$
!!$  end subroutine kronec
!!$
!!$  !================================================================================================
!!$  subroutine globalid_ip(ijk,nd,ndime,gl)
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
!!$    integer(ip), intent(out) :: gl
!!$
!!$    if(ndime==1) then
!!$       gl = ijk(1)
!!$    elseif(ndime==2) then
!!$       gl = (ijk(1)-1)*nd(2)+ijk(2)
!!$    else if(ndime==3) then
!!$       gl = (ijk(1)-1)*nd(2)*nd(3)+(ijk(2)-1)*nd(3)+ijk(3)
!!$    end if
!!$
!!$  end subroutine globalid_ip
!!$
!!$  !================================================================================================
!!$  subroutine globalid_igp(ijk,nd,ndime,gl)
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
!!$    integer(igp), intent(out) :: gl
!!$    
!!$    if(ndime==1) then
!!$       gl = int(ijk(1),igp)
!!$    elseif(ndime==2) then
!!$       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)+int(ijk(2),igp)
!!$    else if(ndime==3) then
!!$       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)*int(nd(3),igp)+(int(ijk(2),igp)-1)*int(nd(3),igp)+int(ijk(3),igp)
!!$    end if
!!$    
!!$  end subroutine globalid_igp
!!$  
!!$  !================================================================================================
!!$  subroutine globalid_2l(ijk,nsckt,npsoc,ndime,gl_2l)
!!$    implicit none
!!$    integer(ip), intent(in)  :: ijk(3),nsckt(3),npsoc(3),ndime
!!$    integer(ip), intent(out) :: gl_2l
!!$
!!$    integer(ip) :: aux1(ndime),aux2(ndime),aux3,aux4,i,gl,local_ijk(ndime)
!!$    real(rp)    :: work
!!$
!!$    ! Auxiliar variables
!!$    aux3=1
!!$    do i=1,ndime
!!$       work = (ijk(i)-1)/npsoc(i)
!!$       aux1(i) = floor(work) + 1;
!!$       aux2(i) = floor(work)*npsoc(i)
!!$       aux3 = aux3*npsoc(i)
!!$    end do
!!$
!!$    ! 2-level global numeration
!!$    call globalid(aux1,nsckt,ndime,aux4)
!!$    aux4 = (aux4-1)*aux3
!!$    local_ijk = ijk(1:ndime) - aux2
!!$    call globalid(local_ijk,npsoc,ndime,gl)
!!$    gl_2l = aux4 + gl
!!$    
!!$  end subroutine globalid_2l
!!$
!!$  !=============================================================================
!!$  subroutine fem_structured_geom_mesh_alloc (ndime,nnode,nedir,msh)
!!$    !-----------------------------------------------------------------------
!!$    ! This routine allocates a partitioned structured mesh
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip)   , intent(in)  :: ndime,nnode,nedir(3)
!!$    type(fem_mesh), intent(out) :: msh
!!$
!!$    integer(ip) :: npoix,npoiy,npoiz,i
!!$
!!$    msh%mtype=structured
!!$    msh%nelty=1
!!$    msh%ndime=ndime
!!$    msh%nnode=nnode
!!$    msh%nedir=nedir
!!$
!!$    if(msh%ndime==2) then
!!$
!!$       if(msh%nnode==3.or.msh%nnode==4) then ! P1 or Q1
!!$          npoix=msh%nedir(1)+1
!!$          npoiy=msh%nedir(2)+1
!!$       else if(msh%nnode==6.or.msh%nnode==9) then ! P2 or Q2
!!$          npoix=2*msh%nedir(1)+1
!!$          npoiy=2*msh%nedir(2)+1
!!$       else if(msh%nnode==10.or.msh%nnode==16) then ! P3 or Q3
!!$          npoix=3*msh%nedir(1)+1
!!$          npoiy=3*msh%nedir(2)+1
!!$       end if
!!$       msh%npoin=npoix*npoiy
!!$
!!$       ! Number of elements
!!$       if(msh%nnode==3.or.msh%nnode==6.or.msh%nnode==10) then ! P1, P2 or P3
!!$          msh%nelem=2*nedir(1)*nedir(2)
!!$       else if(msh%nnode==4.or.msh%nnode==9.or.msh%nnode==16) then ! Q1, Q2 or Q3
!!$          msh%nelem=nedir(1)*nedir(2)
!!$       end if
!!$
!!$    else if(msh%ndime==3) then
!!$       if(msh%nnode==4.or.msh%nnode==8) then ! P1 or Q1
!!$          npoix=msh%nedir(1)+1
!!$          npoiy=msh%nedir(2)+1
!!$          npoiz=msh%nedir(3)+1
!!$       else if(msh%nnode==10.or.msh%nnode==27) then ! P2 or Q2
!!$          npoix=2*msh%nedir(1)+1
!!$          npoiy=2*msh%nedir(2)+1
!!$          npoiz=2*msh%nedir(3)+1
!!$       else if(msh%nnode==20.or.msh%nnode==64) then ! P3 or Q3
!!$          npoix=3*msh%nedir(1)+1
!!$          npoiy=3*msh%nedir(2)+1
!!$          npoiz=3*msh%nedir(3)+1
!!$       end if
!!$       msh%npoin=npoix*npoiy*npoiz
!!$
!!$       ! Number of elements
!!$       if(msh%nnode==4.or.msh%nnode==10.or.msh%nnode==20) then ! P1, P2 or P3
!!$          msh%nelem=4*nedir(1)*nedir(2)*nedir(3)
!!$       else if(msh%nnode==8.or.msh%nnode==27.or.msh%nnode==64) then ! Q1, Q2 or Q3
!!$          msh%nelem=nedir(1)*nedir(2)*nedir(3)
!!$       end if
!!$    end if
!!$
!!$    !call memalloc (1,msh%pnods, __FILE__,__LINE__)
!!$    call memalloc (msh%nelem+1,msh%pnods, __FILE__,__LINE__)
!!$    call memalloc (msh%nnode*msh%nelem,msh%lnods, __FILE__,__LINE__)
!!$    call memalloc (msh%ndime,msh%npoin,msh%coord, __FILE__,__LINE__)
!!$    !msh%pnods(1) = msh%nnode
!!$    msh%pnods(1) = 1
!!$    do i = 1,msh%nelem
!!$       msh%pnods(i+1) = msh%pnods(i) + msh%nnode
!!$    end do
!!$
!!$    return
!!$
!!$  end subroutine fem_structured_geom_mesh_alloc
!!$
!!$  !=============================================================================
!!$  subroutine fem_structured_topo_mesh_alloc (ndime,gsize,tsize,msh)
!!$    !-----------------------------------------------------------------------
!!$    ! This routine allocates a partitioned structured mesh
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip)    , intent(in)  :: ndime
!!$    type(geom_size), intent(in)  :: gsize
!!$    type(topo_size), intent(in)  :: tsize
!!$    type(fem_mesh),  intent(out) :: msh
!!$    integer(ip) :: npoix,npoiy,npoiz,i, aux!,nobje
!!$
!!$
!!$    msh%mtype = structured
!!$    msh%nelty = 1
!!$    msh%ndime = ndime
!!$    msh%nedir = gsize%nedir
!!$    msh%npoin = tsize%notot
!!$    msh%nelem = gsize%nedomt
!!$    !msh%nnode = gsize%nnode
!!$    
!!$    if(msh%ndime==2) then
!!$       !nobje = 8
!!$       msh%nnode = 8
!!$       msh%nface = 4
!!$       msh%nodfac = 2
!!$    else if(msh%ndime==3) then
!!$       !nobje = 26
!!$       msh%nnode = 26
!!$       msh%nface = 6
!!$       msh%nodfac = 4
!!$    end if
!!$
!!$   aux = msh%nnode - msh%nface
!!$
!!$
!!$    msh%tface = 0       ! TO DO
!!$    msh%bou_tface = 0   ! TO DO
!!$
!!$    call memalloc (msh%nelem+1,msh%pnods, __FILE__,__LINE__)
!!$    call memalloc (msh%nnode*msh%nelem,msh%lnods, __FILE__,__LINE__)
!!$    call memalloc (msh%nelem, msh%p_face, __FILE__, __LINE__ )
!!$
!!$    !call memalloc (nobje*msh%nelem,msh%lnods, __FILE__,__LINE__)
!!$    !call memalloc (msh%ndime,msh%npoin,msh%coord, __FILE__,__LINE__)
!!$    
!!$   msh%pnods(1) = 1
!!$    do i = 1,msh%nelem
!!$
!!$       msh%pnods(i+1) = msh%pnods(i) + msh%nnode 
!!$       msh%p_face(i) = msh%pnods(i) + aux
!!$
!!$       !msh%pnods(i+1) = msh%pnods(i) + nobje
!!$    end do
!!$
!!$    return
!!$
!!$  end subroutine fem_structured_topo_mesh_alloc
!!$
!!$  !================================================================================================
!!$  subroutine aux_surf_num(pdime,aux_surf,gsurf_aux)
!!$    implicit none
!!$    integer(ip), intent(in)  :: pdime
!!$    integer(ip), intent(out) :: aux_surf(4,2),gsurf_aux(4)
!!$
!!$    if(pdime==3) then
!!$       gsurf_aux = (/3,4,5,6/)
!!$    else if(pdime==2) then
!!$       gsurf_aux = (/1,2,5,6/)
!!$    else if(pdime==1) then
!!$       gsurf_aux = (/1,2,3,4/)
!!$    end if
!!$    aux_surf(1,:) = (/2,4/)
!!$    aux_surf(2,:) = (/1,3/)
!!$    aux_surf(3,:) = (/3,4/)
!!$    aux_surf(4,:) = (/1,2/)    
!!$
!!$  end subroutine aux_surf_num
!!$
!!$  !================================================================================================
!!$  subroutine global_to_ijk(lpart,nsckt,npsoc,ndime,ijk)
!!$    implicit none
!!$    integer(ip), intent(in)  :: lpart,ndime,nsckt(ndime),npsoc(ndime)
!!$    integer(ip), intent(out) :: ijk(3)
!!$
!!$    integer(ip) :: isckt,jsckt,ksckt,lsckt
!!$    integer(ip) :: ipart_aux,jpart_aux,kpart_aux,lpart_aux
!!$    integer(ip) :: ipart,jpart,kpart
!!$    real(rp) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7
!!$
!!$    if(ndime==2) then
!!$       ! First level (Socket level)
!!$       aux1 = (lpart-1)/(npsoc(1)*npsoc(2))
!!$       lsckt = floor(aux1) + 1
!!$       aux2 = (lsckt-1)/nsckt(2)
!!$       jsckt = lsckt - floor(aux2)*nsckt(2)
!!$       isckt = floor(aux2) + 1
!!$
!!$       ! Second level (Part inside the socket)
!!$       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)
!!$       aux3 = (lpart_aux-1)/npsoc(2)
!!$       jpart_aux = lpart_aux - floor(aux3)*npsoc(2)
!!$       ipart_aux = floor(aux3) + 1
!!$
!!$       ! ijk numeration
!!$       ipart = ipart_aux + (isckt-1)*npsoc(1)
!!$       jpart = jpart_aux + (jsckt-1)*npsoc(2)
!!$       kpart = 1
!!$       ijk = (/ipart,jpart,kpart/)
!!$    else if(ndime==3) then
!!$       ! First level (Socket level)
!!$       aux1 = (lpart-1)/(npsoc(1)*npsoc(2)*npsoc(3))
!!$       lsckt = floor(aux1) + 1
!!$       aux2 = (lsckt-1)/(nsckt(2)*nsckt(3))
!!$       isckt = floor(aux2) + 1
!!$       aux3 = (lsckt - (isckt-1)*nsckt(2)*nsckt(3) - 1)/nsckt(3)
!!$       jsckt = floor(aux3) + 1
!!$       aux4 = lsckt - (isckt-1)*nsckt(2)*nsckt(3) - (jsckt-1)*nsckt(3) - 1
!!$       ksckt = floor(aux4) + 1
!!$
!!$       ! Second level (Part inside the socket)
!!$       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)*npsoc(3)
!!$       aux5 = (lpart_aux-1)/(npsoc(2)*npsoc(3))
!!$       ipart_aux = floor(aux5) + 1
!!$       aux6 = (lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - 1)/npsoc(3)
!!$       jpart_aux = floor(aux6) + 1
!!$       aux7 = lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - (jpart_aux-1)*npsoc(3) - 1
!!$       kpart_aux = floor(aux7) + 1
!!$
!!$       ! ijk numeration
!!$       ipart = ipart_aux + (isckt-1)*npsoc(1)
!!$       jpart = jpart_aux + (jsckt-1)*npsoc(2)
!!$       kpart = kpart_aux + (ksckt-1)*npsoc(3)
!!$       ijk = (/ipart,jpart,kpart/)       
!!$    end if
!!$  
!!$  end subroutine global_to_ijk
!!$
!!$  !================================================================================================
!!$  subroutine materialid(mcase,ijkelem,nedir,nelbl,mater)
!!$    implicit none
!!$    integer(ip), intent(in)    :: mcase,ijkelem(3),nedir(3),nelbl(3)
!!$    integer(ip), intent(inout) :: mater
!!$    
!!$    if(mcase==1) then
!!$       ! (x<leng(x)/2 --> mat == 1)
!!$       ! (x>leng(x)/2 --> mat == 2)
!!$       assert(mod(nedir(1),2)==0)
!!$       if(ijkelem(1).le.nedir(1)/2) then
!!$          mater = 1
!!$       else
!!$          mater = 2
!!$       end if
!!$    elseif(mcase==2) then
!!$       ! (y<leng(y)/2 --> mat == 1)
!!$       ! (y>leng(y)/2 --> mat == 2)
!!$       assert(mod(nedir(2),2)==0)
!!$       if(ijkelem(2).le.nedir(2)/2) then
!!$          mater = 1
!!$       else
!!$          mater = 2
!!$       end if
!!$    elseif(mcase==3) then   ! Hunt's case
!!$       ! (z >= -1.0 and z <= 1.0 --> mat == 1) IMH
!!$       ! (z < -1.0 or z > 1.0 --> mat == 2) DCY
!!$       if(ijkelem(3)>nelbl(3) .and. ijkelem(3)<=nedir(3)-nelbl(3)) then  ! We have nelbl elements in the layer (DCY-solid)
!!$          mater = 1
!!$       else
!!$          mater = 2
!!$       end if
!!$    elseif(mcase==4) then   ! Backward-facing step
!!$       ! (x >= 0.0 or y >= 0.0 --> mat == 1) NSI
!!$       ! (x < 0.0 and y < 0.0 --> mat == 2) Boundary
!!$       if(ijkelem(1)>nelbl(1) .or. ijkelem(2)>nelbl(2)) then 
!!$          mater = 1
!!$       else
!!$          mater = 2
!!$       end if
!!$    end if
!!$    
!!$  end subroutine materialid
!!$
!!$  !================================================================================================
!!$  subroutine fem_structured_geom_size_create(npdir,nedir,nsckt,isper,ndime,msize,gsize)
!!$    implicit none
!!$    integer(ip),     intent(in)  :: npdir(3),nedir(3),nsckt(3),isper(3),ndime
!!$    type(mesh_size), intent(in)  :: msize
!!$    type(geom_size), intent(out) :: gsize
!!$
!!$    ! Local variables
!!$    integer(ip) :: pdegr,idime,jdime
!!$
!!$    pdegr = msize%pdegr 
!!$
!!$    ! Directional sizes
!!$    gsize%npdir = npdir
!!$    gsize%nedir = nedir
!!$    gsize%nsckt = nsckt
!!$    gsize%nnode = (pdegr+1)**ndime  
!!$    gsize%npsoc = npdir/nsckt   
!!$    gsize%nedom = nedir/npdir   
!!$    gsize%npdom = pdegr*gsize%nedom + 1  
!!$    gsize%ncorn = npdir + (1-isper) 
!!$    do idime=1,ndime
!!$       gsize%nedge(idime,:) = npdir + 1 - isper
!!$       gsize%nface(idime,:) = npdir
!!$       if(isper(idime)==0) then
!!$          gsize%nface(idime,idime) = gsize%nface(idime,idime) + 1 
!!$       end if
!!$    end do
!!$    if(ndime==2) gsize%nface = 0
!!$
!!$    ! Total sizes
!!$    gsize%npdomt  = 1
!!$    gsize%nedomt  = 1
!!$    gsize%ncornt  = 1
!!$    gsize%nedget  = 1
!!$    gsize%nfacet  = 1
!!$    gsize%nedgett = 0
!!$    gsize%nfacett = 0
!!$    do idime=1,ndime
!!$       gsize%npdomt = gsize%npdomt*gsize%npdom(idime)
!!$       gsize%nedomt = gsize%nedomt*gsize%nedom(idime)
!!$       gsize%ncornt = gsize%ncornt*gsize%ncorn(idime)
!!$       do jdime=1,ndime
!!$          gsize%nedget(idime) = gsize%nedget(idime)*gsize%nedge(idime,jdime)
!!$          gsize%nfacet(idime) = gsize%nfacet(idime)*gsize%nface(idime,jdime) 
!!$       end do
!!$       gsize%nedgett = gsize%nedgett + gsize%nedget(idime)
!!$       gsize%nfacett = gsize%nfacett + gsize%nfacet(idime)
!!$    end do
!!$    
!!$  end subroutine fem_structured_geom_size_create
!!$
!!$  !================================================================================================
!!$  subroutine fem_structured_topo_size_create(ndime,gsize,tsize)
!!$    implicit none
!!$    integer(ip),     intent(in)  :: ndime
!!$    type(geom_size), intent(in)  :: gsize
!!$    type(topo_size), intent(out) :: tsize
!!$
!!$    ! Local variables
!!$    integer(ip) :: pdime,i,nfaux,ndaux
!!$    integer(ip), allocatable :: auxv(:,:)
!!$    
!!$    ! Auxiliar vector of components
!!$    if(ndime==2) then
!!$       allocate(auxv(2,1))
!!$       auxv(1,1) = 2
!!$       auxv(2,1) = 1
!!$    else if(ndime==3) then
!!$       allocate(auxv(3,2))
!!$       auxv(1,:) = (/2,3/)
!!$       auxv(2,:) = (/1,3/)
!!$       auxv(3,:) = (/1,2/)
!!$    end if
!!$
!!$    tsize%notot = 0
!!$    tsize%nctot = 1
!!$    tsize%ncglb = 1
!!$    tsize%ndtot = 0
!!$    tsize%nftot = 0
!!$    do pdime=ndime,1,-1
!!$       nfaux = 1
!!$       ndaux = 1
!!$       do i=1,ndime-1
!!$          nfaux = nfaux*gsize%nedom(auxv(pdime,i))
!!$          ndaux = ndaux*(gsize%nedom(auxv(pdime,i))+1)
!!$          tsize%nfdom(auxv(pdime,i),pdime) = gsize%nedom(auxv(pdime,i))
!!$          tsize%nfglb(auxv(pdime,i),pdime) = gsize%nedir(auxv(pdime,i))          
!!$          tsize%nddom(auxv(pdime,i),pdime) = (gsize%nedom(auxv(pdime,i))+1)
!!$          tsize%ndglb(auxv(pdime,i),pdime) = (gsize%nedir(auxv(pdime,i))+1)        
!!$       end do
!!$       tsize%nfdir(pdime) = (gsize%nedom(pdime)+1)*nfaux
!!$       tsize%nddir(pdime) = gsize%nedom(pdime)*ndaux
!!$       tsize%nfdom(pdime,pdime) = (gsize%nedom(pdime)+1)
!!$       tsize%nfglb(pdime,pdime) = (gsize%nedir(pdime)+1)
!!$       tsize%nddom(pdime,pdime) = gsize%nedom(pdime)
!!$       tsize%ndglb(pdime,pdime) = gsize%nedir(pdime)
!!$       tsize%ncglb = tsize%ncglb*(gsize%nedir(pdime)+1)
!!$       tsize%nctot = tsize%nctot*gsize%npdom(pdime)
!!$       tsize%ndtot = tsize%ndtot + tsize%nddir(pdime)
!!$       tsize%nftot = tsize%nftot + tsize%nfdir(pdime)
!!$    end do 
!!$    do pdime=ndime,1,-1
!!$       nfaux = 1
!!$       ndaux = 1
!!$       do i=1,ndime-1
!!$          nfaux = nfaux*gsize%nedir(auxv(pdime,i))
!!$          ndaux = ndaux*(gsize%nedir(auxv(pdime,i))+1)
!!$       end do
!!$       tsize%ndsum(pdime) = gsize%nedir(pdime)*ndaux
!!$       tsize%nfsum(pdime) = (gsize%nedir(pdime)+1)*nfaux
!!$    end do
!!$    if(ndime==2) then
!!$       tsize%nftot = 0
!!$       tsize%nfdir = 0
!!$       tsize%nfdom = 0
!!$       tsize%nfglb = 0
!!$    end if
!!$    tsize%notot = tsize%nftot + tsize%ndtot + tsize%nctot
!!$
!!$    ! Deallocate auxiliar vector of components
!!$    deallocate(auxv)
!!$    
!!$  end subroutine fem_structured_topo_size_create
!!$
!!$  !================================================================================================
!!$  subroutine fem_mesh_gen_partition_set(npartx,nparty,npartz,nx,ny,nz,lx,ly,lz,npsocx,npsocy,npsocz, &
!!$       &                                isper,nedir,npdir,nsckt,nparts,ndime,msize)
!!$    implicit none
!!$    integer(ip),     intent(in)  :: npartx,nparty,npartz,nx,ny,nz,npsocx,npsocy,npsocz
!!$    real(rp),        intent(in)  :: lx,ly,lz
!!$    integer(ip),     intent(out) :: ndime,nparts
!!$    integer(ip),     intent(out) :: isper(3),nedir(3),npdir(3),nsckt(3)
!!$    type(mesh_size), intent(out) :: msize
!!$
!!$    ! Initialize data
!!$    isper = 0
!!$    nedir = 0
!!$    npdir = 1
!!$    nsckt = 1
!!$
!!$    ! Mesh generation data
!!$    if(npartz==0) then
!!$       ndime = 2
!!$    else
!!$       ndime=3
!!$    end if
!!$    isper(1)=0                   ! If direction x is periodic
!!$    isper(2)=0                   ! If direction y is periodic
!!$    nedir(1)=nx                  ! Number of elements on x direction
!!$    nedir(2)=ny                  ! Number of elements on y direction
!!$    npdir(1)=npartx              ! Number of parts on x direction
!!$    npdir(2)=nparty              ! Number of parts on y direction
!!$    nsckt(1)=npartx/npsocx       ! Number of sockets on x direction
!!$    nsckt(2)=nparty/npsocy       ! Number of sockets on y direction
!!$    nparts  =npartx*nparty
!!$    if(ndime==3) then
!!$       isper(3)=0                   ! If direction z is periodic
!!$       nedir(3)=nz                  ! Number of elements on z direction
!!$       npdir(3)=npartz              ! Number of parts on z direction
!!$       nsckt(3)=npartz/npsocz       ! Number of sockets on z direction
!!$       nparts = nparts*npartz
!!$    end if
!!$
!!$    ! Mesh size
!!$    msize%ntdix = 0        ! Type of discretization in x (0=uniform, 1=refined near walls)
!!$    msize%ntdiy = 0        ! Type of discretization in y (0=uniform, 1=refined near walls)
!!$    msize%xleng = lx       ! Size of the domain in x
!!$    msize%yleng = ly       ! Size of the domain in y
!!$    msize%zleng = lz       ! Size of the domain in z
!!$    msize%zx1   = 0.1_rp   ! size of the elements at x=0   (left)  
!!$    msize%zx2   = 0.1_rp   ! size of the elements at x=a/2 (center)
!!$    msize%zy1   = 0.1_rp   ! size of the elements at y=0   (bottom)
!!$    msize%zy2   = 0.1_rp   ! size of the elements at y=b/2 (center)
!!$
!!$  end subroutine fem_mesh_gen_partition_set
!!$
!!$  !================================================================================================
!!$  subroutine fem_mesh_gen_partition_bcs(ndime,prob,poin,line,surf)
!!$    implicit none
!!$    integer(ip),          intent(in)  :: ndime,prob
!!$    type(fem_conditions), intent(out) :: poin,line,surf
!!$
!!$    integer(ip)          :: ndofn,idof,ipoin,iline,isurf
!!$
!!$    ! Degrees of freedom
!!$    if(prob==1) then       ! CDR
!!$       ndofn=1
!!$    elseif(prob==2) then   ! Stokes
!!$       ndofn=ndime+1
!!$    elseif(prob==3) then   ! MHD
!!$       ndofn=2*(ndime+1)
!!$    end if
!!$    call fem_conditions_create(ndofn,ndofn,2**ndime,poin)
!!$    call fem_conditions_create(ndofn,ndofn,2*(ndime-1)*ndime,line)
!!$    if(ndime==2) then
!!$       call fem_conditions_create(ndofn,ndofn,1,surf)
!!$    else if(ndime==3) then
!!$       call fem_conditions_create(ndofn,ndofn,2*ndime,surf)
!!$    end if
!!$
!!$    ! 2D:
!!$    ! Points code: 1=bottom left; 2=bottom right; 3=top right; 4=top left
!!$    ! Lines code:  1=left line; 2=right line;   3=bottom line;  4=top line 
!!$
!!$    ! 3D:
!!$    ! Points code:  1=(0,0,0), 2=(0,0,1), 3=(0,1,0), 4=(0,1,1), 5=(1,0,0),
!!$    !               6=(1,0,1), 7=(1,1,0), 8=(1,1,1)
!!$    ! Lines code:   1=(0,0,z), 2=(0,1,z), 3=(1,0,z), 4=(1,1,z), 5=(0,y,0),
!!$    !               6=(0,y,1), 7=(1,y,0), 8=(1,y,1), 9=(x,0,0), 10=(x,0,1),
!!$    !               11=(x,1,0), 12=(x,1,1)
!!$    ! Surface code: 1=(x,y,0), 2=(x,y,1), 3=(x,0,z), 4=(x,1,z), 5=(0,y,z),
!!$    !               6=(1,y,z)
!!$
!!$    if(prob==1) then            ! CDR
!!$       ! Point conditions
!!$       do ipoin=1,2**ndime
!!$          poin%code(ndofn,ipoin)=1; poin%valu(:,ipoin)=0.0_rp
!!$       end do
!!$
!!$       ! Line conditions
!!$       do iline=1,2*(ndime-1)*ndime
!!$          line%code(ndofn,iline)=1; line%valu(:,iline)=0.0_rp
!!$       end do
!!$
!!$       ! Suface conditions
!!$       if(ndime==2) then
!!$          surf%code(ndofn,1)=0; surf%valu(ndofn,1)=0.0_rp
!!$       else if(ndime==3) then
!!$          do isurf=1,2*ndime
!!$             surf%code(ndofn,isurf)=0; surf%valu(ndofn,isurf)=0.0_rp
!!$          end do
!!$       end if
!!$
!!$    elseif(prob==2) then          ! Stokes
!!$       ! Point conditions
!!$       do ipoin=1,2**ndime
!!$          do idof=1, ndime
!!$             poin%code(idof,ipoin)=1
!!$          end do
!!$          poin%code(ndofn,ipoin)=0
!!$       end do
!!$       poin%code(ndofn,1)=1  ! hay que fijar la presion en un punto en flujo cerrado
!!$
!!$       ! Point values
!!$       do ipoin=1,2**ndime
!!$          poin%valu(:,ipoin)=0.0_rp
!!$       end do
!!$
!!$       ! Line conditions
!!$       do iline=1,2*(ndime-1)*ndime
!!$          do idof=1, ndime
!!$             line%code(idof,iline)=1
!!$          end do
!!$          line%code(ndofn,iline)=0
!!$       end do
!!$
!!$       ! Line values
!!$       do iline=1,2*(ndime-1)*ndime
!!$          line%valu(:,iline)=0.0_rp
!!$       end do
!!$
!!$       ! Suface conditions
!!$       if(ndime==2) then
!!$          do idof=1, ndime
!!$             surf%code(idof,1)=0; surf%valu(idof,1)=0.0_rp
!!$          end do
!!$          surf%code(ndofn,1)=0; surf%valu(ndofn,1)=0.0_rp
!!$       else if(ndime==3) then
!!$          do isurf=1,2*ndime
!!$             do idof=1, ndime
!!$                surf%code(idof,isurf)=1; surf%valu(idof,isurf)=0.0_rp
!!$             end do
!!$             surf%code(ndofn,1)=0; surf%valu(ndofn,1)=0.0_rp
!!$          end do
!!$       end if
!!$
!!$       ! Stokes driven cavity flow
!!$       !if(ndime==2) then
!!$       !    line%code(1,3)=1
!!$       !    line%valu(1,3)=1.0_rp
!!$       !else if(ndime==3) then
!!$       !    surf%code(1,2)=1
!!$       !    surf%valu(1,2)=1.0_rp
!!$       !end if
!!$
!!$    elseif(prob==3) then          ! MHD
!!$       ! Point conditions 
!!$       do ipoin=1,2**ndime    
!!$          do idof=1, ndime        
!!$             poin%code(idof,ipoin)=1          ! Velocity
!!$             poin%code(ndime+1+idof,ipoin)=1  ! Induction
!!$          end do
!!$          poin%code(ndime+1,ipoin)=0          ! Pressure
!!$          poin%code(ndofn,ipoin)=1            ! Magnetic pressure
!!$       end do
!!$       poin%code(ndime+1,1)=1  ! Pressure (hay que fijar la presion en un punto en flujo cerrado)
!!$
!!$       ! Point values
!!$       do ipoin=1,2**ndime 
!!$          poin%valu(:,ipoin)=0.0_rp
!!$       end do
!!$
!!$       ! Line conditions
!!$       do iline=1,2*(ndime-1)*ndime
!!$          do idof=1, ndime        
!!$             line%code(idof,iline)=1  ! Velocity
!!$             if(ndime==3) then
!!$                line%code(ndime+1+idof,iline)=1
!!$             end if
!!$          end do
!!$          line%code(ndime+1,iline)=0  ! Pressure
!!$          line%code(ndofn,iline)=1    ! Magnetic pressure
!!$       end do
!!$       if(ndime==2) then
!!$          line%code(ndime+2,1)=1; line%code(ndime+3,1)=0  ! Induction (bx; by) we impose the tangential comp.
!!$          line%code(ndime+2,2)=0; line%code(ndime+3,2)=1
!!$          line%code(ndime+2,3)=1; line%code(ndime+3,3)=0
!!$          line%code(ndime+2,4)=0; line%code(ndime+3,4)=1
!!$       end if
!!$
!!$       ! Line values
!!$       do iline=1,2*(ndime-1)*ndime
!!$          line%valu(:,iline)=0.0_rp
!!$       end do
!!$
!!$       ! Suface conditions
!!$       if(ndime==2) then
!!$          do idof=1, ndime
!!$             surf%code(idof,1)=0; surf%valu(idof,1)=0.0_rp
!!$          end do
!!$          surf%code(ndofn,1)=0; surf%valu(ndofn,1)=0.0_rp
!!$       else if(ndime==3) then
!!$          do isurf=1,2*ndime
!!$             do idof=1, ndime
!!$                surf%code(idof,isurf)=0; surf%valu(idof,isurf)=0.0_rp
!!$             end do
!!$             surf%code(ndofn,1)=0; surf%valu(ndofn,1)=0.0_rp
!!$          end do
!!$          surf%code(ndime+2,1)=1; surf%code(ndime+3,1)=1; surf%code(ndime+4,1)=0;
!!$          surf%code(ndime+2,2)=1; surf%code(ndime+3,2)=1; surf%code(ndime+4,2)=0;
!!$          surf%code(ndime+2,3)=1; surf%code(ndime+3,3)=0; surf%code(ndime+4,3)=1;
!!$          surf%code(ndime+2,4)=1; surf%code(ndime+3,4)=0; surf%code(ndime+4,4)=1;
!!$          surf%code(ndime+2,5)=0; surf%code(ndime+3,5)=1; surf%code(ndime+4,5)=1;
!!$          surf%code(ndime+2,6)=0; surf%code(ndime+3,6)=1; surf%code(ndime+4,6)=1;
!!$          do isurf=1,2*ndime
!!$             surf%code(ndofn,1)=0; 
!!$          end do
!!$          ! Surface values
!!$          do isurf=1,2*ndime
!!$             do idof=1, ndime
!!$                surf%valu(ndofn+1+idof,1)=0.0_rp; 
!!$             end do
!!$          end do
!!$       end if
!!$    end if
!!$
!!$  end subroutine fem_mesh_gen_partition_bcs

end module fem_mesh_gen_partition
