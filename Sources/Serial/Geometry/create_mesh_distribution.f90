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
module create_mesh_distribution_names
  use types_names
  use list_types_names
  use memor_names
  use mesh_distribution_names
  use metis_interface_names
  use rcm_renumbering_names
  use mesh_names
  use hash_table_names
  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: create_mesh_distribution

contains

  subroutine create_mesh_distribution( femesh, prt_pars, distr, lmesh)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    class(mesh_t)                   , intent(inout)      :: femesh
    type(mesh_distribution_params_t), intent(in)         :: prt_pars
    type(mesh_distribution_t) , allocatable, intent(out) :: distr(:) ! Mesh distribution instances
    type(mesh_t)              , allocatable, intent(out) :: lmesh(:) ! Local mesh instances

    ! Local variables
    type(list_t)                 :: fe_graph    ! Dual graph (to be partitioned)
    integer(ip)   , allocatable  :: ldome(:)    ! Part of each element
    integer(ip)                  :: ipart
    integer                      :: istat

    ! Generate dual mesh (i.e., list of elements around points)
    call femesh%to_dual()

    ! Create dual (i.e. list of elements around elements)
    call create_dual_graph(femesh,fe_graph)

    ! Partition dual graph to assign a domain to each element (in ldome)
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)
    call graph_pt_renumbering(prt_pars,fe_graph,ldome)

    ! Now free fe_graph, not needed anymore?
    allocate(distr(prt_pars%nparts), stat=istat)
    check(istat==0)
    allocate(lmesh(prt_pars%nparts), stat=istat)
    check(istat==0) 

    do ipart=1, prt_pars%nparts
       distr(ipart)%ipart  = ipart
       distr(ipart)%nparts = prt_pars%nparts
    end do
    call build_maps(prt_pars%nparts, ldome, femesh, distr)

    ! Build local meshes and their duals and generate partition adjacency
    do ipart=1,prt_pars%nparts

       ! Generate Local mesh
       call mesh_g2l(distr(ipart)%num_local_vertices,  &
                     distr(ipart)%l2g_vertices,        &
                     distr(ipart)%num_local_cells,     &
                     distr(ipart)%l2g_cells,           &
                     femesh,                           &
                     lmesh(ipart))

       call build_adjacency_new (femesh, ldome,             &
            &                    ipart,                     &
            &                    lmesh(ipart),              &
            &                    distr(ipart)%l2g_vertices, &
            &                    distr(ipart)%l2g_cells,    &
            &                    distr(ipart)%nebou,        &
            &                    distr(ipart)%nnbou,        &
            &                    distr(ipart)%lebou,        &
            &                    distr(ipart)%lnbou,        &
            &                    distr(ipart)%pextn,        &
            &                    distr(ipart)%lextn,        &
            &                    distr(ipart)%lextp )
    end do
    call fe_graph%free()
    call memfree(ldome,__FILE__,__LINE__)
  end subroutine create_mesh_distribution

  !================================================================================================
  subroutine create_dual_graph(mesh,graph)
    ! Parameters
    type(mesh_t) , intent(in)  :: mesh
    type(list_t),  intent(out) :: graph
    ! Locals
    integer(ip), allocatable :: lelem(:)
    integer(ip), allocatable :: keadj(:)

    call memalloc(           mesh%nelem,lelem,__FILE__,__LINE__)
    call memalloc(mesh%nelpo*mesh%nnode,keadj,__FILE__,__LINE__)
    lelem=0
    call graph%create(mesh%nelem)

    call count_elemental_graph(mesh%ndime,mesh%npoin,mesh%nelem, &
         &                     mesh%pnods,mesh%lnods,mesh%nnode,mesh%nelpo, &
         &                     mesh%pelpo,mesh%lelpo,lelem,keadj,graph%p)

    call graph%allocate_list_from_pointer()

    call list_elemental_graph(mesh%ndime,mesh%npoin,mesh%nelem, &
         &                    mesh%pnods,mesh%lnods,mesh%nnode,mesh%nelpo, &
         &                    mesh%pelpo,mesh%lelpo,lelem,keadj,graph%p,graph%l)

    call memfree(lelem,__FILE__,__LINE__)
    call memfree(keadj,__FILE__,__LINE__)

  end subroutine create_dual_graph

  !-----------------------------------------------------------------------
  subroutine count_elemental_graph(ncomm,npoin,nelem,pnods,lnods,nnode, &
       &                           nelpo,pelpo,lelpo,lelem,keadj,ieadj)
    implicit none
    integer(ip), intent(in)  :: ncomm,npoin,nelem
    integer(ip), intent(in)  :: pnods(nelem+1),lnods(pnods(nelem+1))
    integer(ip), intent(in)  :: nnode             ! Number of nodes of each element (max.)
    integer(ip), intent(in)  :: nelpo             ! Number of elements around points (max.)
    integer(ip), intent(in)  :: pelpo(npoin+1)    ! Number of elements around points
    integer(ip), intent(in)  :: lelpo(pelpo(npoin+1))    ! List of elements around points
    integer(ip), intent(out) :: ieadj(nelem+1)           ! Number of edges on each element
    integer(ip), intent(out) :: lelem(nelem)             ! Auxiliar array
    integer(ip), intent(out) :: keadj(nelpo*nnode)       ! Auxiliar array
    integer(ip)              :: ielem,jelem,inode,knode  ! Indices
    integer(ip)              :: ipoin,ielpo,index        ! Indices
    integer(ip)              :: neadj,ielel,jelel,nelel  ! Indices

    lelem=0
    neadj=1
    knode=nnode
    do ielem=1,nelem
       ieadj(ielem)=neadj
       ! Loop over nodes and their surrounding elements and count
       ! how many times they are repeated as neighbors of ielem
       nelel=0
       keadj=0
       index=pnods(ielem)-1
       knode=pnods(ielem+1)-pnods(ielem)
       do inode=1,knode
          ipoin=lnods(index+inode)
          do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
             jelem=lelpo(ielpo)
             if(lelem(jelem)==0) then
                nelel=nelel+1
                keadj(nelel)=jelem
             end if
             lelem(jelem)=lelem(jelem)+1
          end do
       end do

       ! Now we loop over the elements around ielem and define neighbors
       ! as those sharing ncomm nodes. The smaller this number is the
       ! higher the connectivity of elemental graph is. Note that prisms,
       ! for example, could share 3 or 4 nodes depending on the face, so
       ! ndime (the number of space dimensions) is a good choice.
       jelel=0
       do ielel=1,nelel
          jelem=keadj(ielel)
          if(lelem(jelem)>=ncomm) jelel=jelel+1
       end do
       neadj=neadj+jelel

       ! Reset lelem
       do ielel=1,nelel
          jelem=keadj(ielel)
          lelem(jelem)=0
       end do
    end do

    ieadj(nelem+1)=neadj

  end subroutine count_elemental_graph
  !-----------------------------------------------------------------------
  subroutine list_elemental_graph(ncomm,npoin,nelem,pnods,lnods,nnode, &
       &                          nelpo,pelpo,lelpo,lelem,keadj,ieadj,jeadj)
    implicit none
    integer(ip), intent(in)  :: ncomm,npoin,nelem
    integer(ip), intent(in)  :: pnods(nelem+1),lnods(pnods(nelem+1))
    integer(ip), intent(in)  :: nnode             ! Number of nodes of each element (max.)
    integer(ip), intent(in)  :: nelpo             ! Number of elements around points (max.)
    integer(ip), intent(in)  :: pelpo(npoin+1)    ! Number of elements around points
    integer(ip), intent(in)  :: lelpo(pelpo(npoin+1))    ! List of elements around points
    integer(ip), intent(in)  :: ieadj(nelem+1)           ! Number of edges on each element
    integer(ip), intent(out) :: jeadj(ieadj(nelem+1))    ! List of edges on each element
    integer(ip), intent(out) :: lelem(nelem)             ! Auxiliar array
    integer(ip), intent(out) :: keadj(nelpo*nnode)       ! Auxiliar array
    integer(ip)              :: ielem,jelem,inode,knode  ! Indices
    integer(ip)              :: ipoin,ielpo,index        ! Indices
    integer(ip)              :: neadj,ielel,jelel,nelel  ! Indices

    lelem=0
    knode=nnode
    do ielem=1,nelem
       ! Loop over nodes and their surrounding elements and count
       ! how many times they are repeated as neighbors of ielem
       nelel=0
       keadj=0
       index=pnods(ielem)-1
       knode=pnods(ielem+1)-pnods(ielem)
       do inode=1,knode
          ipoin=lnods(index+inode)
          do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
             jelem=lelpo(ielpo)
             if(lelem(jelem)==0) then
                nelel=nelel+1
                keadj(nelel)=jelem
             end if
             lelem(jelem)=lelem(jelem)+1
          end do
       end do

       ! Now we loop over the elements around ielem and define neighbors
       jelel=ieadj(ielem)-1
       do ielel=1,nelel
          jelem=keadj(ielel)
          if(lelem(jelem)>=ncomm) then
             jelel=jelel+1
             jeadj(jelel)=jelem
          end if
       end do

       ! Reset lelem
       do ielel=1,nelel
          jelem=keadj(ielel)
          lelem(jelem)=0
       end do
    end do

  end subroutine list_elemental_graph

  !================================================================================================
  subroutine build_adjacency_new ( gmesh, ldome, my_part, lmesh, l2gn, l2ge, &
       &                           nebou, nnbou, lebou, lnbou, pextn, lextn, lextp)
    implicit none
    integer(ip)   , intent(in)  :: my_part
    type(mesh_t)  , intent(in)  :: gmesh,lmesh
    !type(mesh_t)  , intent(in)  :: dual_lmesh
    integer(ip)   , intent(in)  :: ldome(gmesh%nelem)
    integer(igp)  , intent(in)  :: l2gn(lmesh%npoin)
    integer(igp)  , intent(in)  :: l2ge(lmesh%nelem)
    !integer(ip)   , intent(in)  :: dual_parts( dual_lmesh%pnods(dual_lmesh%nelem+1)-1)
    integer(ip)   , intent(out) :: nebou
    integer(ip)   , intent(out) :: nnbou
    integer(ip)   , allocatable, intent(out) ::  lebou(:)    ! List of boundary elements
    integer(ip)   , allocatable, intent(out) ::  lnbou(:)    ! List of boundary nodes
    integer(ip)   , allocatable, intent(out) ::  pextn(:)    ! Pointers to the lextn
    integer(igp)  , allocatable, intent(out) ::  lextn(:)    ! List of (GID of) external neighbors
    integer(ip)   , allocatable, intent(out) ::  lextp(:)    ! List of parts of external neighbors

    integer(ip) :: lelem, ielem, jelem, pelem, pnode, inode1, inode2, ipoin, lpoin, jpart, iebou, istat, touch
    integer(ip) :: nextn, nexte, nepos
    integer(ip), allocatable :: local_visited(:)
    type(hash_table_ip_ip_t)   :: external_visited

    ! Count boundary nodes
    nnbou = 0 
    do lpoin=1, lmesh%npoin
       ipoin = l2gn(lpoin)
       do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
          ielem = gmesh%lelpo(pelem)
          jpart = ldome(ielem)
          if ( jpart /= my_part ) then 
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! List boundary nodes
    call memalloc ( nnbou, lnbou, __FILE__, __LINE__ ) 
    nnbou = 0
    do lpoin=1, lmesh%npoin
       ipoin = l2gn(lpoin)
       do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
          ielem = gmesh%lelpo(pelem)
          jpart = ldome(ielem)
          if ( jpart /= my_part ) then 
             lnbou(nnbou+1) = ipoin
             nnbou = nnbou +1
             exit
          end if
       end do
    end do

    ! As the dual mesh is given with global IDs we need a hash table to do the touch.
    call memalloc(lmesh%nelem, local_visited,__FILE__,__LINE__)
    local_visited = 0
    call external_visited%init(20)

    ! 1) Count boundary elements and external edges
    touch = 1
    nebou = 0 ! number of boundary elements
    nextn = 0 ! number of external edges
    do lelem = 1, lmesh%nelem
       nexte = 0   ! number of external neighbours of this element
       ielem = l2ge(lelem)
       inode1 = gmesh%pnods(lelem)
       inode2 = gmesh%pnods(lelem+1)-1
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   if(local_visited(lelem) == 0 ) nebou = nebou +1        ! Count it
                   !call external_visited%put(key=jelem,val=1, stat=istat) ! Touch jelem as external neighbor of lelem.
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) nexte = nexte + 1                ! Count external neighbours of lelem
                   local_visited(lelem) = nexte                           ! Touch lelem also storing the number
                end if                                                    ! of external neighbours it has
             end if
          end do
       end do
       nextn = nextn + nexte
       ! Clean hash table
       if(local_visited(lelem) /= 0 ) then 
          do pnode = inode1, inode2
             ipoin = gmesh%lnods(pnode)
             do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
                jelem = gmesh%lelpo(pelem)
                if(jelem/=ielem) then
                   jpart = ldome(jelem)
                   if(jpart/=my_part) then
                      call external_visited%del(key=jelem, stat=istat)
                   end if
                end if
             end do
          end do
       end if
       call external_visited%print
    end do

    ! 2) Allocate arrays and store list and pointers to externals
    call memalloc(nebou  , lebou,__FILE__,__LINE__)
    call memalloc(nebou+1, pextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextn,__FILE__,__LINE__)
    call memalloc(nextn  , lextp,__FILE__,__LINE__)

    iebou = 0
    pextn(1) = 1
    do lelem = 1, lmesh%nelem
       if(local_visited(lelem) /= 0 ) then
          iebou = iebou +1
          lebou(iebou) = lelem
          pextn(iebou+1) = local_visited(lelem) + pextn(iebou)
       end if
    end do

    ! 3) Store boundary elements and external edges
    !do lelem = 1, lmesh%nelem
    do iebou = 1, nebou
       lelem = lebou(iebou)
       ielem = l2ge(lelem)
       nexte = 0   ! number of external neighbours of this element
       inode1 = gmesh%pnods(lelem)
       inode2 = gmesh%pnods(lelem+1)-1
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   call external_visited%put(key=jelem,val=touch, stat=istat) ! Touch jelem as external neighbor of lelem.
                   if(istat==now_stored) then
                      lextn(pextn(iebou)+nexte) = jelem
                      lextp(pextn(iebou)+nexte) = jpart
                      nexte = nexte + 1
                   end if
                end if
             end if
          end do
       end do
       ! Clean hash table
       do pnode = inode1, inode2
          ipoin = gmesh%lnods(pnode)
          do pelem = gmesh%pelpo(ipoin), gmesh%pelpo(ipoin+1) - 1
             jelem = gmesh%lelpo(pelem)
             if(jelem/=ielem) then
                jpart = ldome(jelem)
                if(jpart/=my_part) then                                   ! This is an external element
                   call external_visited%del(key=jelem, stat=istat)
                end if
             end if
          end do
       end do
    end do

    call external_visited%free
    call memfree(local_visited,__FILE__,__LINE__)

  end subroutine build_adjacency_new

  !================================================================================================
  subroutine build_maps( nparts, ldome, femesh, distr )
    ! This routine builds (node and element) partition maps without using the objects
    ! and (unlike parts_sizes, parts_maps, etc.) does not generate a new global numbering.
    implicit none
    integer(ip)                , intent(in)    :: nparts
    type(mesh_t)             , intent(in)    :: femesh
    integer(ip)                , intent(in)    :: ldome(femesh%nelem)
    type(mesh_distribution_t), intent(inout) :: distr(nparts)

    integer(ip)   , allocatable  :: nedom(:) ! Number of points per part (here is not header!)
    integer(ip)   , allocatable  :: npdom(:) ! Number of elements per part (here is not header!)
    integer(ip)   , allocatable  :: work1(:)
    integer(ip)   , allocatable  :: work2(:)
    integer(ip) :: ielem, ipart, inode, iboun

    ! Number of elements of each part and global to local element map (is one to one)
    call memalloc (nparts, nedom,__FILE__,__LINE__)
    nedom=0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
    end do
    ! Allocate local to global maps
    do ipart=1,nparts
       distr(ipart)%num_local_cells  = nedom(ipart)
       distr(ipart)%num_global_cells = int(femesh%nelem,igp)
       call memalloc(distr(ipart)%num_local_cells, distr(ipart)%l2g_cells, __FILE__, __LINE__)
    end do
    nedom = 0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
       distr(ipart)%l2g_cells(nedom(ipart)) = ielem
    end do

    call memfree ( nedom,__FILE__,__LINE__)

    ! Number of nodes of each part and global to local node map (is NOT one to one)
    call memalloc ( nparts, npdom,__FILE__,__LINE__)
    call memalloc ( femesh%npoin, work1,__FILE__,__LINE__)
    call memalloc ( femesh%npoin, work2,__FILE__,__LINE__)
    npdom=0
    do ipart = 1, nparts
       work1 = 0
       work2 = 0
       do ielem=1,femesh%nelem
          if(ldome(ielem)==ipart) then
             do inode = femesh%pnods(ielem), femesh%pnods(ielem+1) - 1 
                if(work1(femesh%lnods(inode)) == 0 ) then
                   npdom(ipart) = npdom(ipart)+1
                   work1(femesh%lnods(inode)) = 1
                   work2(npdom(ipart)) = femesh%lnods(inode)
                end if
             end do
          end if
       end do
       distr(ipart)%num_local_vertices  = npdom(ipart)
       distr(ipart)%num_global_vertices = int(femesh%npoin,igp)
       call memalloc(distr(ipart)%num_local_vertices, distr(ipart)%l2g_vertices, __FILE__, __LINE__)
       distr(ipart)%l2g_vertices = work2(1:npdom(ipart))
    end do
    call memfree ( work1,__FILE__,__LINE__)
    call memfree ( work2,__FILE__,__LINE__)
    call memfree ( npdom,__FILE__,__LINE__)
  end subroutine build_maps

  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search.
  ! Given a mesh (m) and its dual graph (g), it computes the list 
  ! of nodes (lconn) of each connected component in m. Can be very
  ! useful as a tool to determine whether the mesh partitioning process
  ! leads to disconnected subdomains or not.
  subroutine mesh_graph_compute_connected_components (m, g, lconn)
    implicit none

    ! Parameters
    type(mesh_t) , intent(in)   :: m   
    type(list_t),  intent(in)   :: g
    type(list_t),  intent(out)  :: lconn

    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
         j, l, k, inods1d, inods2d, p_ipoin, ipoin, graph_num_rows
    type(list_iterator_t)    :: graph_column_iterator

    graph_num_rows = g%get_num_pointers()
    call memalloc ( graph_num_rows   , auxe     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , auxv     , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , q        , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%npoin          , vmarked  , __FILE__,__LINE__)
    call memalloc ( graph_num_rows   ,  e       , __FILE__,__LINE__)

    lconn%n  = 0
    emarked  = 0
    current  = 1 

    do i=1, graph_num_rows
       if (emarked(i) == 0) then
          ! New connected component
          lconn%n = lconn%n +1
          esize   = 0
          vsize   = 0
          vmarked = 0 
!!$1  procedure BFS(G,v):
!!$2      create a queue Q
          head=1
          tail=1
!!$3      enqueue v onto Q
          q(tail)=i
          tail=tail+1
!!$4      mark v
          emarked(i)=1
          e(current)=i
          esize  = esize + 1
          current = current + 1  

!!$5      while Q is not empty:
          do while (head/=tail)
!!$6         t ← Q.dequeue()
             j=q(head)
             head = head + 1

             ! Traverse the nodes of the element number j
             inods1d = m%pnods(j)
             inods2d = m%pnods(j+1)-1

             do p_ipoin = inods1d, inods2d
                ipoin = m%lnods(p_ipoin)
                if (vmarked(ipoin)==0) then
                   vmarked(ipoin)=1
                   vsize = vsize+1
                end if
             end do

!!$9         for all edges e in G.adjacentEdges(t) do
             graph_column_iterator = g%create_iterator(j)
             do while(.not. graph_column_iterator%is_upper_bound())
!!$12           u ← G.adjacentVertex(t,e)
                l=graph_column_iterator%get_current()
!!$13           if u is not emarked:
                if (emarked(l)==0) then
!!$14              mark u
                   emarked(l)=1
                   e(current)=l
                   esize  = esize + 1
                   current = current + 1  

!!$15              enqueue u onto Q
                   q(tail)=l
                   tail=tail+1
                end if
                call graph_column_iterator%next()
             end do
          end do
          auxe(lconn%n) = esize
          auxv(lconn%n) = vsize
       end if
    end do

    call lconn%create(lconn%n)

    lconn%p(1) = 1
    do i=1, lconn%n
       lconn%p(i+1) = lconn%p(i) + auxv(i)
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call lconn%allocate_list_from_pointer()

    current=1
    l=1
    do i=1, lconn%n
       vmarked = 0
       ! Traverse elements of current connected component  
       do current=current,current+auxe(i)-1
          j=e(current)

          ! Traverse the nodes of the element number j
          inods1d = m%pnods(j)
          inods2d = m%pnods(j+1)-1

          do p_ipoin = inods1d, inods2d
             ipoin = m%lnods(p_ipoin)
             if (vmarked(ipoin)==0) then
                vmarked(ipoin)=1
                lconn%l(l)=ipoin
                l=l+1
             end if
          end do

       end do
    end do

    ! write(*,*) 'ZZ', g%nv, m%npoin, l
    ! write(*,*) 'XX', lconn%n
    ! write(*,*) 'YY', lconn%p
    ! write(*,*) 'PP', lconn%l

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)

  end subroutine mesh_graph_compute_connected_components

  !=================================================================================================
  subroutine graph_nd_renumbering(prt_parts, gp, iperm, lperm)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_distribution_params_t), intent(in)         :: prt_parts
    type(list_t)                    , target, intent(in) :: gp
    integer(ip)                     , target, intent(out):: iperm(gp%n)
    integer(ip)                     , target, intent(out):: lperm(gp%n)
    
    if ( gp%get_num_pointers() == 1 ) then
       lperm(1) = 1
       iperm(1) = 1
    else
#ifdef ENABLE_METIS
       ierr = metis_setdefaultoptions(c_loc(options))
       assert(ierr == METIS_OK) 
       
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ierr = metis_nodend ( c_loc(gp%n) ,c_loc(gp%p) ,c_loc(gp%l),C_NULL_PTR,c_loc(options), &
            &                c_loc(iperm),c_loc(lperm))
       
       assert(ierr == METIS_OK)
#else
       call enable_metis_error_message
#endif
    end if
  end subroutine graph_nd_renumbering

  !=================================================================================================
  subroutine graph_pt_renumbering(prt_parts,gp,ldomn)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of the input graph gp
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_distribution_params_t), target, intent(in)    :: prt_parts
    type(list_t)                    , target, intent(inout) :: gp
    integer(ip)                     , target, intent(out)   :: ldomn(gp%n)

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)

   
#ifdef ENABLE_METIS
    ierr = metis_setdefaultoptions(c_loc(options))
    assert(ierr == METIS_OK) 

!!$      From METIS 5.0 manual:
!!$
!!$      The following options are valid for METIS PartGraphRecursive:
!!$      
!!$      METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
!!$      METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
!!$      METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL
!!$     
!!$      The following options are valid for METIS PartGraphKway:
!!$ 
!!$      METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
!!$      METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
!!$      METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
!!$      METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL

    if ( prt_parts%strat == part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = prt_parts%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the subdomain graph
       options(METIS_OPTION_MINCONN)   = prt_parts%metis_option_minconn
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor
       
       ncon = 1 
       ierr = metis_partgraphkway( c_loc(gp%n) , c_loc(ncon), c_loc(gp%p)   , c_loc(gp%l) , & 
                                   C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                   C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       
       assert(ierr == METIS_OK) 
       
    else if ( prt_parts%strat == part_recursive ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor

       ncon = 1 
       ierr = metis_partgraphrecursive( c_loc(gp%n) , c_loc(ncon), c_loc(gp%p)   , c_loc(gp%l) , & 
                                        C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
                                        C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
    end if    
#else
    call enable_metis_error_message
#endif

    if ( prt_parts%strat == part_strip ) then
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( prt_parts%strat == part_rcm_strip ) then
       call memalloc ( gp%get_num_pointers(), iperm, __FILE__,__LINE__ )
       call genrcm ( gp%get_num_pointers(), gp%p(gp%get_num_pointers()+1)-1, gp%p, gp%l, iperm )
       j = gp%get_num_pointers()
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(iperm(m+i)) = ipart
          end do
          m = m + k
          j = j - k
       end do
       call memfree ( iperm,__FILE__,__LINE__)
    end if

  end subroutine graph_pt_renumbering


  !================================================================================================
  subroutine mesh_g2l(num_local_vertices, l2g_vertices, num_local_cells, l2g_cells, gmesh, lmesh)
    implicit none
    integer(ip),     intent(in)    :: num_local_vertices
    integer(igp),    intent(in)    :: l2g_vertices(num_local_vertices)
    integer(ip),     intent(in)    :: num_local_cells
    integer(igp),    intent(in)    :: l2g_cells(num_local_cells)
    type(mesh_t)   , intent(in)    :: gmesh
    type(mesh_t)   , intent(inout) :: lmesh
    type(hash_table_igp_ip_t)      :: ws_inmap
    type(hash_table_igp_ip_t)      :: el_inmap
    integer(ip)    , allocatable   :: node_list(:)
    integer(ip)                    :: aux, ipoin,inode,inodb,knode,knodb,lnodb_size,istat
    integer(ip)                    :: ielem_lmesh,ielem_gmesh,iboun_lmesh,iboun_gmesh
    integer(ip)                    :: p_ielem_gmesh,p_ipoin_lmesh,p_ipoin_gmesh
    logical :: count_it


    lmesh%order=gmesh%order
    lmesh%nelty=gmesh%nelty
    lmesh%ndime=gmesh%ndime
    lmesh%npoin=num_local_vertices
    lmesh%nelem=num_local_cells

    call ws_inmap%init(max(int(num_local_vertices*0.25,ip),10))
    do ipoin=1,num_local_vertices
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call ws_inmap%put(key=l2g_vertices(ipoin),val=aux,stat=istat) 
    end do

    call el_inmap%init(max(int(num_local_cells*0.25,ip),10))
    do ipoin=1,num_local_cells
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call el_inmap%put(key=l2g_cells(ipoin),val=aux,stat=istat) 
    end do

    ! Elements
    call memalloc(lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    call memalloc(lmesh%nelem  , lmesh%legeo, __FILE__,__LINE__)
    call memalloc(lmesh%nelem  , lmesh%leset, __FILE__,__LINE__)
    lmesh%nnode=0
    lmesh%pnods=0
    lmesh%pnods(1)=1
    do ielem_lmesh=1,lmesh%nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
       lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
       lmesh%nnode=max(lmesh%nnode,knode)
       lmesh%legeo(ielem_lmesh)=gmesh%legeo(ielem_gmesh)
       lmesh%leset(ielem_lmesh)=gmesh%leset(ielem_gmesh)
    end do
    call memalloc (lmesh%pnods(lmesh%nelem+1), lmesh%lnods, __FILE__,__LINE__)
    do ielem_lmesh=1,lmesh%nelem
       ielem_gmesh = l2g_cells(ielem_lmesh)
       p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
       p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
       knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
       do inode=1,knode
          call ws_inmap%get(key=int(gmesh%lnods(p_ipoin_gmesh+inode),igp),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
       end do
    end do

    ! Boundary elements
    iboun_lmesh=0
    lmesh%nnodb=0
    lnodb_size=0
    do iboun_gmesh=1,gmesh%nboun
       p_ipoin_gmesh = gmesh%pnodb(iboun_gmesh)-1
       knodb = gmesh%pnodb(iboun_gmesh+1)-gmesh%pnodb(iboun_gmesh)
       count_it=.true.
       do inode=1,knodb
          call ws_inmap%get(key=int(gmesh%lnodb(p_ipoin_gmesh+inode),igp),val=knode,stat=istat)
          if(istat==key_not_found) then
             count_it=.false.
             exit
          end if
       end do
       if(count_it) then
          lnodb_size=lnodb_size+knodb
          lmesh%nnodb=max(lmesh%nnodb,knodb)
          iboun_lmesh=iboun_lmesh+1
       end if
    end do

    if(iboun_lmesh>0) then
       lmesh%nboun=iboun_lmesh
       call memalloc (  lmesh%nnodb,   node_list, __FILE__,__LINE__)
       call memalloc (lmesh%nboun+1, lmesh%pnodb, __FILE__,__LINE__)
       call memalloc (   lnodb_size, lmesh%lnodb, __FILE__,__LINE__)
       call memalloc(   lmesh%nboun, lmesh%lbgeo, __FILE__,__LINE__)
       call memalloc(   lmesh%nboun, lmesh%lbset, __FILE__,__LINE__)

       lmesh%pnodb=0
       lmesh%pnodb(1)=1
       iboun_lmesh=0
       do iboun_gmesh=1,gmesh%nboun
          p_ipoin_gmesh = gmesh%pnodb(iboun_gmesh)-1
          knodb = gmesh%pnodb(iboun_gmesh+1)-gmesh%pnodb(iboun_gmesh)
          count_it=.true.
          do inode=1,knodb
             call ws_inmap%get(key=int(gmesh%lnodb(p_ipoin_gmesh+inode),igp),val=node_list(inode),stat=istat)
             if(istat==key_not_found) then
                count_it=.false.
                exit
             end if
          end do
          if(count_it) then
             iboun_lmesh=iboun_lmesh+1
             lmesh%pnodb(iboun_lmesh+1)=lmesh%pnodb(iboun_lmesh)+knodb
             p_ipoin_lmesh = lmesh%pnodb(iboun_lmesh)-1
             lmesh%lnodb(p_ipoin_lmesh+1:p_ipoin_lmesh+knodb)=node_list(1:knodb)
             lmesh%lbgeo(iboun_lmesh)=gmesh%lbgeo(iboun_gmesh)
             lmesh%lbset(iboun_lmesh)=gmesh%lbset(iboun_gmesh)
          end if
       end do
       call memfree (node_list, __FILE__,__LINE__)
    end if
    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    !call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord)
    do ipoin=1,num_local_vertices
       lmesh%coord(:,ipoin)=gmesh%coord(:,l2g_vertices(ipoin))
    end do

  end subroutine mesh_g2l

end module create_mesh_distribution_names
