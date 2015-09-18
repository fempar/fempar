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
  use memor_names
  use mesh_distribution_names
  use map_names
  use map_apply_names
  use graph_names
  use graph_renumbering_names
  use mesh_names
  use partitioning_params_names
  use hash_table_names
# include "debug.i90"
  implicit none
  private

   ! Functions
  public :: create_mesh_distribution, create_graph_from_mesh

contains

  subroutine create_mesh_distribution( prt_pars, femesh, distr, lmesh)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(partitioning_params_t)                        , intent(in)  :: prt_pars
    type(mesh_t)                           , intent(in)  :: femesh
    type(mesh_distribution_t) , allocatable, intent(out) :: distr(:) ! Mesh distribution instances
    type(mesh_t)              , allocatable, intent(out) :: lmesh(:) ! Local mesh instances

    ! Local variables
    type(mesh_t)               :: dual_femesh, dual_lmesh
    type(graph_t)              :: fe_graph    ! Dual graph (to be partitioned)
    integer(ip)   , allocatable  :: ldome(:)    ! Part of each element
    integer(ip)   , allocatable  :: dual_parts(:)
    integer(ip)                  :: ipart
    integer                      :: istat
    ! AFM. Dummy map declared and used for backward compatibility.
    ! It should be re-considered when we decide how to implement
    ! the boundary mesh.
    type(map_t)                    :: dummy_bmap


    ! Generate dual mesh (i.e., list of elements around points)
    call mesh_to_dual(femesh, dual_femesh)
    
    ! Allocate working arrays
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)

    ! dual graph (elements around elements)
    call create_graph_from_mesh (dual_femesh, femesh, dual_femesh%ndime, fe_graph)
    !out0(call graph_print(6, fe_graph))
          
    ! write (*,*) 'fe_graph%nv', fe_graph%nv, 'fe_graph%nnz', fe_graph%ia(fe_graph%nv+1) 
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
       call mesh_g2l(distr(ipart)%nmap, distr(ipart)%emap, dummy_bmap, femesh, lmesh(ipart))

       ! Mesh to dual with global element numbers in the dual (just get numbers applying g2l map to dual_femesh)
       ! Store dual_part too.
       call dual_mesh_g2l(distr(ipart)%nmap, dual_femesh, ldome, lmesh(ipart), dual_lmesh, dual_parts)

       call build_adjacency (ipart,                 &
            &                lmesh(ipart),          &
            &                distr(ipart)%emap%l2g, &
            &                dual_lmesh,            &
            &                dual_parts,            &
            &                distr(ipart)%nebou,    &
            &                distr(ipart)%nnbou,    &
            &                distr(ipart)%lebou,    &
            &                distr(ipart)%lnbou,    &
            &                distr(ipart)%pextn,    &
            &                distr(ipart)%lextn,    &
            &                distr(ipart)%lextp )

       call mesh_free(dual_lmesh)
       call memfree(dual_parts,__FILE__,__LINE__)
    end do

    call mesh_free(dual_femesh)
    call fe_graph%free()
    call memfree(ldome,__FILE__,__LINE__)
  end subroutine create_mesh_distribution

  !================================================================================================
   subroutine build_adjacency ( my_part, lmesh, l2ge, dual_lmesh, dual_parts, &
        &                       nebou, nnbou, lebou, lnbou, pextn, lextn, lextp)
     implicit none
     integer(ip)   , intent(in)  :: my_part
     type(mesh_t), intent(in)  :: lmesh
     type(mesh_t), intent(in)  :: dual_lmesh
     integer(igp)  , intent(in)  :: l2ge(lmesh%nelem)
     integer(ip)   , intent(in)  :: dual_parts( dual_lmesh%pnods(dual_lmesh%nelem+1)-1)
     integer(ip)   , intent(out) :: nebou
     integer(ip)   , intent(out) :: nnbou
     integer(ip)   , allocatable, intent(out) ::  lebou(:)    ! List of boundary elements
     integer(ip)   , allocatable, intent(out) ::  lnbou(:)    ! List of boundary nodes
     integer(ip)   , allocatable, intent(out) ::  pextn(:)    ! Pointers to the lextn
     integer(igp)  , allocatable, intent(out) ::  lextn(:)    ! List of (GID of) external neighbors
     integer(ip)   , allocatable, intent(out) ::  lextp(:)    ! List of parts of external neighbors

     integer(ip) :: lelem, ielem, jelem, pelem, pnode, inode1, inode2, ipoin, jpart, iebou, istat, touch
     integer(ip) :: nextn, nexte, nepos
     integer(ip), allocatable :: local_visited(:)
     type(hash_table_ip_ip_t)   :: external_visited

     ! Count boundary nodes
     nnbou = 0 
     do ipoin=1, lmesh%npoin
        do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
           jpart = dual_parts(pelem)
           if ( jpart /= my_part ) then 
              nnbou = nnbou +1
              exit
           end if
        end do
     end do

     ! List boundary nodes
     call memalloc ( nnbou, lnbou, __FILE__, __LINE__ ) 
     nnbou = 0
     do ipoin=1, lmesh%npoin
        do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
           jpart = dual_parts(pelem)
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
        inode1 = lmesh%pnods(lelem)
        inode2 = lmesh%pnods(lelem+1)-1
        do pnode = inode1, inode2
           ipoin = lmesh%lnods(pnode)
           do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
              jelem = dual_lmesh%lnods(pelem)
              if(jelem/=ielem) then
                 jpart = dual_parts(pelem)
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
              ipoin = lmesh%lnods(pnode)
              do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
                 jelem = dual_lmesh%lnods(pelem)
                 if(jelem/=ielem) then
                    jpart = dual_parts(pelem)
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


     ! 3) Store Count boundary elements and external edges
     !do lelem = 1, lmesh%nelem
     do iebou = 1, nebou
        lelem = lebou(iebou)
        ielem = l2ge(lelem)
        nexte = 0   ! number of external neighbours of this element
        inode1 = lmesh%pnods(lelem)
        inode2 = lmesh%pnods(lelem+1)-1
        do pnode = inode1, inode2
           ipoin = lmesh%lnods(pnode)
           do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
              jelem = dual_lmesh%lnods(pelem)
              if(jelem/=ielem) then
                 jpart = dual_parts(pelem)
                 if(jpart/=my_part) then                            ! This is an external element
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
           ipoin = lmesh%lnods(pnode)
           do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
              jelem = dual_lmesh%lnods(pelem)
              if(jelem/=ielem) then
                 jpart = dual_parts(pelem)
                 if(jpart/=my_part) then
                    call external_visited%del(key=jelem, stat=istat)
                 end if
              end if
           end do
        end do
     end do

     call external_visited%free
     call memfree(local_visited,__FILE__,__LINE__)
   end subroutine build_adjacency

  subroutine dual_mesh_g2l(nmap, dual_mesh, ldome, lmesh, dual_lmesh, dual_parts)
    implicit none
    type(map_igp_t) , intent(in)  :: nmap
    type(mesh_t), intent(in)  :: dual_mesh
    integer(ip)   , intent(in)  :: ldome(dual_mesh%npoin)
    type(mesh_t), intent(in)  :: lmesh
    type(mesh_t), intent(inout) :: dual_lmesh
    integer(ip)   , allocatable, intent(inout)  :: dual_parts(:)

    integer(ip) :: ipart,lelem,ielem, pnode,i

    dual_lmesh%nelem = lmesh%npoin
    dual_lmesh%npoin = lmesh%nelem
    call memalloc (dual_lmesh%nelem+1, dual_lmesh%pnods, __FILE__,__LINE__)
    dual_lmesh%pnods(1) = 1
    do lelem = 1, dual_lmesh%nelem
       ielem = nmap%l2g(lelem)
       dual_lmesh%pnods(lelem+1) = dual_mesh%pnods(ielem+1) - dual_mesh%pnods(ielem) + dual_lmesh%pnods(lelem)
    end do
    call memalloc (dual_lmesh%pnods(dual_lmesh%nelem+1), dual_lmesh%lnods, __FILE__,__LINE__)
    call memalloc (dual_lmesh%pnods(dual_lmesh%nelem+1), dual_parts, __FILE__,__LINE__)
    do lelem = 1, dual_lmesh%nelem
       ielem = nmap%l2g(lelem)
       pnode = dual_lmesh%pnods(lelem+1) - dual_lmesh%pnods(lelem)
       ! assert( pnode == dual_mesh%pnods(ielem+1) - dual_mesh%pnods(ielem))
       do i = 0, pnode-1
          dual_lmesh%lnods( dual_lmesh%pnods(lelem) + i ) =  dual_mesh%lnods( dual_mesh%pnods(ielem) + i )
          dual_parts( dual_lmesh%pnods(lelem) + i ) =  ldome(dual_mesh%lnods( dual_mesh%pnods(ielem) + i ))
       end do
    end do
    
  end subroutine dual_mesh_g2l

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
       call map_alloc(nedom(ipart),int(femesh%nelem,igp),distr(ipart)%emap)
    end do
    nedom = 0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
       distr(ipart)%emap%l2g(nedom(ipart)) = ielem
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
       call map_alloc(npdom(ipart),int(femesh%npoin,igp),distr(ipart)%nmap)
       distr(ipart)%nmap%l2g = work2(1:npdom(ipart))
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
    type(graph_t), intent(in)   :: g
    type(list_t)     , intent(out)  :: lconn
 
    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
                                j, l, k, inods1d, inods2d, p_ipoin, ipoin

    call memalloc ( g%nv   , auxe     , __FILE__,__LINE__)
    call memalloc ( g%nv   , auxv     , __FILE__,__LINE__)
    call memalloc ( g%nv   , q        , __FILE__,__LINE__)
    call memalloc ( g%nv   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%npoin, vmarked  , __FILE__,__LINE__)
    call memalloc ( g%nv  ,  e        , __FILE__,__LINE__)

    lconn%n  = 0
    emarked  = 0
    current  = 1 

    do i=1, g%nv
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
             do k=g%ia(j), g%ia(j+1)-1
!!$12           u ← G.adjacentVertex(t,e)
                l=g%ja(k)
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
             end do
          end do
          auxe(lconn%n) = esize
          auxv(lconn%n) = vsize
       end if
    end do
    
    call memalloc( lconn%n+1, lconn%p, __FILE__,__LINE__)

    lconn%p(1) = 1
    do i=1, lconn%n
       lconn%p(i+1) = lconn%p(i) + auxv(i)
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call memalloc( lconn%p(lconn%n+1)-1, lconn%l, __FILE__,__LINE__)

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

  !================================================================================================
  subroutine create_graph_from_mesh (primal_mesh, dual_mesh, min_freq_neig, primal_graph)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! The routine exploits the duality among primal/dual meshes. In order to 
    ! work properly pnods has to be an array of pointers to address lnods
    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
    !
    ! min_freq_neig is an input integer which determines which are
    ! the neighbours of a particular vertex of the graph. In particular, 
    ! neighbours are those vertices which are listed at least min_freq_neig times 
    ! when traversing the list of primal elements around a primal point
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, create_graph_from_mesh_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(mesh_t) , intent(in)  :: primal_mesh, dual_mesh
    integer(ip), intent(in)      :: min_freq_neig
    type(graph_t), intent(out) :: primal_graph

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(4)            ! Pointers to work space

    ! Allocate space for ia on the primal graph
    primal_graph%nv                = primal_mesh%npoin
    primal_graph%nv2               = primal_mesh%npoin
	primal_graph%symmetric_storage = .false.
    call memalloc (primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_primal_graph and list_primal_graph routines
    ! (TOTAL WS SIZE = primal mesh npoin + 2*maximum number of neighbours of any primal
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + primal_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    pwork(4) = pwork(3) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(4), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_primal_graph' ! DBG:
    ! Calculate ia array
    call count_primal_graph ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
         &                           iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)), &
         &                           iwork(pwork(3):pwork(4)) ) 
    ! write(*,*) 'end count_primal_graph'   ! DBG:

    ! Allocate space for ja on the primal graph 
    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,           __FILE__,__LINE__)

    ! write(*,*) 'begin list_primal_graph' ! DBG:
    ! Calculate ja array
    call list_primal_graph  (primal_mesh, dual_mesh, min_freq_neig, primal_graph, & 
         &                   iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)),  &
         &                   iwork(pwork(3):pwork(4)) )
    ! write(*,*) 'end list_primal_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine create_graph_from_mesh

  !============================================================================================
  subroutine  count_primal_graph ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                  ws_position, ws_freq, ws_neighbors)
    implicit none

    ! Parameters
    type(mesh_t) , intent(in)     :: primal_mesh, dual_mesh
    integer(ip), intent(in)         :: min_freq_neig
    type(graph_t), intent(inout)  :: primal_graph
    integer(ip), intent(out)        :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)        :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                 :: ineigh  
    integer(ip)                 :: ipoindm ! Dual   mesh point identifier
    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
    integer(ip)                 :: ipoinpg ! Primal graph point identifier 
    integer(ip)                 :: p_ipoindm   
    integer(ip)                 :: p_ipoinpm   
    integer(ip)                 :: inods1 , inods2   
    integer(ip)                 :: inods1d, inods2d

    integer(ip)                 :: first_free_pos
    integer(ip)                 :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_mesh%npoin
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0
       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)

       inods1d=dual_mesh%pnods(ipoinpg)
       inods2d=dual_mesh%pnods(ipoinpg+1)-1
       do p_ipoindm = inods1d, inods2d
       !do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          inods1=primal_mesh%pnods(ipoindm)
          inods2=primal_mesh%pnods(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm) 

             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq (ws_position(ipoinpm)) = ws_freq( ws_position(ipoinpm) ) + 1
             end if

          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if (.not.ws_neighbors(ineigh)==ipoinpg) then ! Exclude self-edges
                primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

  end subroutine count_primal_graph

  !============================================================================================
  subroutine  list_primal_graph ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                 ws_position, ws_freq, ws_neighbors )
    implicit none
    ! Parameters
    type(mesh_t) , intent(in)    :: primal_mesh, dual_mesh
    type(graph_t), intent(inout) :: primal_graph
    integer(ip), intent(in)        :: min_freq_neig
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_mesh%npoin       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       inods1d=dual_mesh%pnods(ipoinpg)
       inods2d=dual_mesh%pnods(ipoinpg+1)-1

       do p_ipoindm = inods1d, inods2d
       ! do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm
          inods1=primal_mesh%pnods(ipoindm)
          inods2=primal_mesh%pnods(ipoindm+1)-1
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq ( ws_position(ipoinpm) ) = ws_freq( ws_position(ipoinpm) ) + 1
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if ( .not.ws_neighbors(ineigh)==ipoinpg ) then ! Exclude self-edges
                primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
                first_free_ja_pos = first_free_ja_pos + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_primal_graph


end module create_mesh_distribution_names
