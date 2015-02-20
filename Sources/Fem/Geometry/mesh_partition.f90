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
module fem_mesh_partition
  use types
  use memor
  use sort_class
  use fem_partition_class
  use sep_tree_class
  use renum_class
  use maps_class
  use map_apply
  use fem_graph_class
  use graph_renum
  use fem_mesh_class
  use mesh_graph
  use cartesian
  use post
  use fem_materials_class
  use fem_mesh_partition_base
  use hash_table_class
# include "debug.i90"
  implicit none
  private

  interface fem_mesh_partition_create
     module procedure fem_mesh_partition_create_new, fem_mesh_partition_create_old
  end interface fem_mesh_partition_create

  ! Functions
  public :: fem_mesh_partition_create, fem_mesh_partition_from_periodic, build_partition_adjacency

contains
  !================================================================================================
  ! Methods defined for fem_mesh are:
  !
  ! - fem_mesh_partition_create
  ! - fem_mesh_get_parts
  ! - fem_mesh_partition_write
  !
  !================================================================================================

  subroutine fem_mesh_partition_create_new( prt_pars, femesh, parts, lmesh)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(part_params)   , intent(in)  :: prt_pars
    type(fem_mesh)      , intent(in)  :: femesh
    type(fem_partition) , intent(out) :: parts(prt_pars%nparts) ! Partition
    type(fem_mesh)      , intent(out) :: lmesh(prt_pars%nparts) ! Local meshes

    ! Local variables
    type(fem_mesh)               :: dual_femesh, dual_lmesh
    type(fem_graph)              :: fe_graph    ! Dual graph (to be partitioned)
    type(fem_graph)              :: parts_graph
    type(renum)                  :: eren        ! Element renumbering
    integer(ip)   , allocatable  :: ldome(:)    ! Part of each element
    integer(ip)   , allocatable  :: dual_parts(:)
    integer(ip)                  :: ipart

    assert ( prt_pars%ptype==element_based )
    assert ( prt_pars%use_graph==dual )

    ! Generate dual mesh (i.e., list of elements around points)
    call mesh_to_dual(femesh, dual_femesh)
    
    ! Allocate working arrays
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)

    ! dual graph (elements around elements)
    call mesh_to_graph (dual_femesh, femesh, dual_femesh%ndime, fe_graph)
    !out0(call fem_graph_print(6, fe_graph))
          
    ! write (*,*) 'fe_graph%nv', fe_graph%nv, 'fe_graph%nnz', fe_graph%ia(fe_graph%nv+1) 
    call graph_pt_renumbering(prt_pars,fe_graph,eren,ldome)

    call build_parts_graph (prt_pars%nparts, ldome, fe_graph, parts_graph)

    ! Now free fe_graph, not needed anymore?

    ! In the multilevel setting, simply call the last two procedures recursively
    ! reallocating data
    ! fe_graph <- parts_graph
    ! call graph_pt_renumbering(prt_pars,fe_graph,eren,ldome)
    ! parts(ipart)%id_parts(ilevel) = ldome(ipart)
    ! call build_parts_graph (prt_pars%nparts, ldome, fe_graph, parts_graph)

    do ipart=1,prt_pars%nparts
       parts(ipart)%pinfo  = extended_adjacencies
       parts(ipart)%ptype  = element_based
       parts(ipart)%ipart  = ipart
       parts(ipart)%nparts = prt_pars%nparts
    end do

    call build_partition_maps(prt_pars%nparts, ldome, femesh, parts)

    ! Build local meshes and their duals and generate partition adjacency
    do ipart=1,prt_pars%nparts

       parts(ipart)%pinfo=extended_adjacencies

       ! Local mesh
       call fem_mesh_g2l(parts(ipart)%nmap,parts(ipart)%emap, parts(ipart)%bmap, femesh, lmesh(ipart))

       ! Mesh to dual with global element numbers in the dual (just get numbers applying g2l map to dual_femesh)
       ! Store dual_part too.
       call dual_mesh_g2l(parts(ipart)%nmap, dual_femesh, ldome, lmesh(ipart), dual_lmesh, dual_parts)

       call build_partition_adjacency (ipart, &
            &                          lmesh(ipart), &
            &                          parts(ipart)%emap%l2g, &
            &                          dual_lmesh, &
            &                          dual_parts, &
            &                          parts(ipart)%nebou, &
            &                          parts(ipart)%nnbou, &
            &                          parts(ipart)%lebou, &
            &                          parts(ipart)%lnbou, &
            &                          parts(ipart)%pextn, &
            &                          parts(ipart)%lextn, &
            &                          parts(ipart)%lextp, &
            &                          parts(ipart)%lexte, &
            &                          parts(ipart)%lextm)

       call fem_mesh_free(dual_lmesh)

    end do

    call fem_mesh_free(dual_femesh)
    call fem_graph_free(fe_graph)
    call fem_graph_free(parts_graph)
    call memfree(ldome,__FILE__,__LINE__)
    call memfree(dual_parts,__FILE__,__LINE__)

  end subroutine fem_mesh_partition_create_new

  !================================================================================================
   subroutine build_partition_adjacency ( my_part, lmesh, l2ge, dual_lmesh, dual_parts, &
        &                                 nebou, nnbou, lebou, lnbou, pextn, lextn, lextp, lexte, lextm)
     implicit none
     integer(ip)   , intent(in)  :: my_part
     type(fem_mesh), intent(in)  :: lmesh
     type(fem_mesh), intent(in)  :: dual_lmesh
     integer(igp)  , intent(in)  :: l2ge(lmesh%npoin)
     integer(ip)   , intent(in)  :: dual_parts( dual_lmesh%pnods(dual_lmesh%nelem+1) )
     integer(ip)   , intent(out) :: nebou
     integer(ip)   , intent(out) :: nnbou
     integer(ip)   , allocatable, intent(out) ::  lebou(:)    ! List of boundary elements
     integer(ip)   , allocatable, intent(out) ::  lnbou(:)    ! List of boundary nodes
     integer(ip)   , allocatable, intent(out) ::  pextn(:)    ! Pointers to the lextn
     integer(igp)  , allocatable, intent(out) ::  lextn(:)    ! List of (GID of) external neighbors
     integer(ip)   , allocatable, intent(out) ::  lextp(:)    ! List of parts of external neighbors
     integer(ip)   , allocatable, intent(out) ::  lexte(:)    ! Edge information of external neighbors
     integer(ip)   , allocatable, optional, intent(out) ::  lextm(:)    ! Edge information of external neighbors

     integer(ip) :: lelem, ielem, jelem, pelem, pnode, inode1, inode2, ipoin, jpart, iebou, istat
     integer(ip) :: nextn, nexte, nepos
     integer(ip), allocatable :: local_visited(:)
     type(hash_table_ip_ip)   :: external_visited

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
     nebou = 0 ! number of boundary elements
     nextn = 0 ! number of external edges
     do lelem = 1, lmesh%nelem
        nexte = 0   ! number of external neighbours of this element
        ielem = l2ge(lelem)
        if(lmesh%nelty==1) then
           inode1 = (lelem-1)*lmesh%nnode + 1
           inode2 = lelem*lmesh%nnode
        else
           inode1 = lmesh%pnods(lelem)
           inode2 = lmesh%pnods(lelem+1)-1
        end if
        do pnode = inode1, inode2
           ipoin = lmesh%lnods(pnode)
           do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
              jelem = dual_lmesh%lnods(pelem)
              if(jelem/=ielem) then
                 jpart = dual_parts(pelem)
                 if(jpart/=my_part) then                                   ! This is an external element
                    if(local_visited(lelem) == 0 ) nebou = nebou +1        ! Count it
                    call external_visited%put(key=jelem,val=1, stat=istat) ! Touch jelem as external neighbor of lelem.
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
     end do

     ! 2) Allocate arrays and store list and pointers to externals
     call memalloc(nebou  , lebou,__FILE__,__LINE__)
     call memalloc(nebou+1, pextn,__FILE__,__LINE__)
     call memalloc(nextn  , lextn,__FILE__,__LINE__)
     call memalloc(nextn  , lextp,__FILE__,__LINE__)
     call memalloc(nextn  , lexte,__FILE__,__LINE__)
     lexte = 0
     if(present(lextm)) then ! TODO: material
        call memalloc(nextn  , lextm,__FILE__,__LINE__)
        lextm = 0
     end if

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
        if(lmesh%nelty==1) then
           inode1 = (lelem-1)*lmesh%nnode + 1
           inode2 = lelem*lmesh%nnode
        else
           inode1 = lmesh%pnods(lelem)
           inode2 = lmesh%pnods(lelem+1)-1
        end if
        do pnode = inode1, inode2
           ipoin = lmesh%lnods(pnode)
           do pelem = dual_lmesh%pnods(ipoin), dual_lmesh%pnods(ipoin+1) - 1
              jelem = dual_lmesh%lnods(pelem)
              if(jelem/=ielem) then
                 jpart = dual_parts(pelem)
                 if(jpart/=my_part) then                            ! This is an external element
                    call external_visited%put(key=jelem,val=nexte, stat=istat) ! Touch jelem as external neighbor of lelem.
                    if(istat==now_stored) then
                       lextn(pextn(iebou)+nexte) = jelem
                       lextp(pextn(iebou)+nexte) = jpart
                       lexte(pextn(iebou)+nexte) =  ibset( lexte(pextn(iebou)+nexte), pnode-inode1 )
                       nexte = nexte + 1
                    else
                       assert(istat==was_stored)
                       call external_visited%get(key=jelem,val=nepos, stat=istat) ! Touch jelem as external neighbor of lelem.
                       lexte(pextn(iebou)+nepos) =  ibset( lexte(pextn(iebou)+nepos), pnode-inode1 )
                    end if
                 end if
              end if
           end do
        end do
        nextn = nextn + nexte
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

   end subroutine build_partition_adjacency

  subroutine dual_mesh_g2l(nmap, dual_mesh, ldome, lmesh, dual_lmesh, dual_parts)
    implicit none
    type(map)     , intent(in)  :: nmap
    type(fem_mesh), intent(in)  :: dual_mesh
    integer(ip)   , intent(in)  :: ldome(dual_mesh%npoin)
    type(fem_mesh), intent(in)  :: lmesh
    type(fem_mesh), intent(out) :: dual_lmesh
    integer(ip)   , allocatable, intent(out)  :: dual_parts(:)

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

  subroutine build_parts_graph (nparts, ldome, fe_graph, parts_graph)
    ! This procedure is order nparts**2 both in memory and complexity. This number could be
    ! reduced using hash_tables but it is not easy: we need to use npart hash_tables to store
    ! the touched parts in a sparse structure. Each table should be of size max_nparts, which
    ! could be bounded by the maximum number of elements per node. 
    !
    ! However, this procedure will be executed serially, and therefore, we do not expect
    ! nparts to be huge.
    ! 
    implicit none
    integer(ip)                  , intent(in)  :: nparts
    type(fem_graph)              , intent(in)  :: fe_graph
    integer(ip)                  , intent(in)  :: ldome(fe_graph%nv)
    type(fem_graph)              , intent(out) :: parts_graph
    integer(ip), allocatable :: work(:,:)
    integer(ip)              :: ielem,jelem,pelem,iz, ipart,jpart

    parts_graph%nv = nparts
    call memalloc(parts_graph%nv+1, parts_graph%ia, __FILE__,__LINE__)
    parts_graph%ia = 0

    !write(*,*) fe_graph%nv, fe_graph%ia(fe_graph%nv+1)

    call memalloc(nparts, nparts, work, __FILE__,__LINE__)
    work = 0
    parts_graph%nzt = 0
    do ielem = 1, fe_graph%nv
       ipart = ldome(ielem)
       do pelem = fe_graph%ia(ielem), fe_graph%ia(ielem+1) - 1
          jelem = fe_graph%ja(pelem)
          jpart = ldome(jelem)
          if(work(jpart,ipart) == 0) then
             work(jpart,ipart) = 1
             parts_graph%ia(ipart+1) = parts_graph%ia(ipart+1) +1
             parts_graph%nzt = parts_graph%nzt + 1
          end if
       end do
    end do

    call memalloc(parts_graph%nzt, parts_graph%ja, __FILE__,__LINE__)
    parts_graph%nzt = parts_graph%nzt + 1

    ! This loop could be replaced by a loop over the fe_graph vertices, as the previous one,
    ! whose complexity is of order fe_graph%nv compared to the parts_graph%nv**2.
    ! However, that would require
    ! work = 0
    ! which also of order parts_graph%nv**2.
    iz = 0
    parts_graph%ia(1) = 1
    do ipart = 1, parts_graph%nv
       ! Transform ia from length to header
       parts_graph%ia(ipart+1) = parts_graph%ia(ipart+1) + parts_graph%ia(ipart)
       do jpart = 1, parts_graph%nv
          if(work(jpart,ipart) == 1) then
             iz = iz + 1
             parts_graph%ja(iz) = jpart
          end if
       end do
    end do

    call memfree(work, __FILE__,__LINE__)

  end subroutine build_parts_graph

  !================================================================================================
  subroutine build_partition_maps(nparts, ldome, femesh, parts)
    ! This routine builds (node and element) partition maps without using the objects
    ! and (unlike parts_sizes, parts_maps, etc.) does not generate a new global numbering.
    implicit none
    integer(ip)                  , intent(in)  :: nparts
    type(fem_mesh)               , intent(in)  :: femesh
    integer(ip)                  , intent(in)  :: ldome(femesh%nelem)
    type(fem_partition)          , intent(inout) :: parts(nparts)
    integer(ip)   , allocatable  :: nedom(:)    ! Number of points per part (here is not header!)
    integer(ip)   , allocatable  :: npdom(:)    ! Number of elements per part (here is not header!)
    integer(ip)   , allocatable  :: nbdom(:)    ! Number of boundary elements per part
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
       call map_alloc(nedom(ipart),int(femesh%nelem,igp),parts(ipart)%emap)
    end do
    nedom = 0
    do ielem=1,femesh%nelem
       ipart = ldome(ielem)
       nedom(ipart)=nedom(ipart)+1
       parts(ipart)%emap%l2g(nedom(ipart)) = ielem
    end do

    ! Copy emap to l2ge
    do ipart = 1, nparts
       if (parts(ipart)%pinfo==extended_adjacencies) then
          call memalloc (nedom(ipart),parts(ipart)%l2ge,__FILE__,__LINE__)
          parts(ipart)%nelem = nedom(ipart)
          parts(ipart)%l2ge = parts(ipart)%emap%l2g
       end if
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
       call map_alloc(npdom(ipart),femesh%npoin,parts(ipart)%nmap)
       parts(ipart)%nmap%l2g = work2(1:npdom(ipart))
    end do
    call memfree ( work1,__FILE__,__LINE__)
    call memfree ( work2,__FILE__,__LINE__)

    ! Copy emap to l2gn
    do ipart = 1, nparts
       if (parts(ipart)%pinfo==extended_adjacencies) then
          call memalloc (npdom(ipart),parts(ipart)%l2gn,__FILE__,__LINE__)
          parts(ipart)%npoin = npdom(ipart)
          parts(ipart)%l2gn = parts(ipart)%nmap%l2g
       end if
    end do
    call memfree ( npdom,__FILE__,__LINE__)

    ! Create map for boundary elements
    if(femesh%nboun>0) then
       call memalloc (nparts, nbdom, __FILE__,__LINE__)
       nbdom=0
       do ielem=1,femesh%nboun
          nbdom(ldome(femesh%lboel(femesh%nnodb+1,ielem))) = &
               nbdom(ldome(femesh%lboel(femesh%nnodb+1,ielem))) + 1
       end do
       nbdom=0
       do iboun=1,femesh%nboun
          ipart = ldome(femesh%lboel(femesh%nnodb+1,iboun))
          nbdom(ipart) = nbdom(ipart) + 1
          parts(ipart)%bmap%l2g(nbdom(ipart)) = iboun
       end do
       call memfree (nbdom, __FILE__,__LINE__)
    end if

  end subroutine build_partition_maps

  !================================================================================================
  !================================================================================================
  !================================================================================================
  !================================================================================================

  subroutine fem_mesh_partition_create_old( prt_pars, femesh, nren, eren, parts, material, flag)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(part_params)            , intent(in)  :: prt_pars
    type(fem_mesh)               , intent(in)  :: femesh
    type(renum)                  , intent(out) :: nren        ! Node renumbering
    type(renum)                  , intent(out) :: eren        ! Element renumbering
    type(fem_partition)          , intent(out) :: parts(prt_pars%nparts)
    type(fem_materials), optional, intent(in)  :: material    ! Material associated to elements
    integer(ip),         optional, intent(in)  :: flag        ! Flag for new mesh_to_graph subroutine

    ! Local variables
    type(fem_mesh)               :: dual_femesh
    type(fem_graph)              :: fe_graph    ! Either primal or dual (to be partitioned)
    type(sep_tree)               :: stree       ! Separator_tree
    integer(ip)   , allocatable  :: ldome(:)    ! Part (or sep_tree_node) of each element
    integer(ip)   , allocatable  :: ldomn(:)    ! Part (or sep_tree_node) of each vertex
    type(fem_graph)              :: parts_graph
    integer(ip)                  :: max_nparts  ! Maximum number of parts around a vertex
    integer(ip)                  :: nobjs       ! number of objects
    integer(ip)   , allocatable  :: lobjs(:,:)  ! List of objects
    type(list)                   :: part_objs   ! Currently pdobj,ldobj
    type(list)                   :: int_objs    ! Currently piobj,liobj
    integer(ip)   , allocatable  :: edges(:)    ! Edge information
    integer(ip)   , allocatable  :: nedom(:)    ! Number of points per part
    integer(ip)   , allocatable  :: npdom(:)    ! Number of elements per part
    integer(ip)   , allocatable  :: nbdom(:)    ! Number of boundary elements per part
    type(map)     , allocatable  :: nmaps(:)    ! Nodal maps array
    type(map)     , allocatable  :: emaps(:)    ! Elemental maps array
    type(map)     , allocatable  :: bmaps(:)    ! Boundary elemental maps array
    integer(ip)                  :: ipart,flag_

    if(present(flag)) then
       flag_ = flag
    else
       flag_ = 0
    end if

    assert(flag_==0.or.flag_==1)

    ! Create separator tree
    call sep_tree_create(prt_pars%nparts, stree)

    ! Generate dual mesh (i.e., list of elements around points)
    call mesh_to_dual(femesh, dual_femesh)
    
    ! Allocate working arrays
    call memalloc (femesh%npoin, ldomn, __FILE__,__LINE__)
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)

    ! Allocate renumberings
    call renum_alloc(femesh%npoin, nren)
    call renum_alloc(femesh%nelem, eren)

    if(prt_pars%ptype==element_based) then

       if(prt_pars%use_graph==primal) then

          if(stree%nlevel==0) then
             write(*,*) 'Graph partition by nested disection only possible when'
             write(*,*) 'nparts=2^k for some integer k.'
             stop
          end if

          !primal graph (points around points)
          call mesh_to_graph(femesh,dual_femesh,1,fe_graph)
          ! call fem_graph_print(6, fe_graph) ! DBG:

          call graph_nd_renumbering(prt_pars,fe_graph,stree%nlevel,nren,ldomn) 

          ! call cartesian_nd_renumbering(femesh%ndime,femesh%npoin,nparts,femesh%coord,ldomn)
          call renum_by_sets(stree%nnode,ldomn,nren)

          call sep_tree_set_pointers(femesh%npoin,ldomn,stree)

          call sep_to_part(stree,femesh,dual_femesh,ldomn,ldome)

          call renum_by_sets(prt_pars%nparts,ldome,eren)

       else if(prt_pars%use_graph==dual) then

          ! dual graph (elements around elements)
          if(flag_==0) then
             call mesh_to_graph (dual_femesh, femesh, dual_femesh%ndime, fe_graph)
          elseif(flag_==1) then
             call mesh_to_graph( femesh, dual_femesh, fe_graph, flag_ )
          end if
          ! call graph_write(6, fe_graph) DBG:
          
          ! write (*,*) 'fe_graph%nv', fe_graph%nv, 'fe_graph%nnz', fe_graph%ia(fe_graph%nv+1) 

          call graph_pt_renumbering(prt_pars,fe_graph,eren,ldome)

          call part_to_sep(stree,dual_femesh,ldome,ldomn)

          call sep_tree_set_pointers(femesh%npoin,ldomn,stree)

          call renum_by_sets(femesh%npoin,ldomn,nren)

       end if

       ! Create graph of parts
       call parts_graph_create(dual_femesh,stree,ldome,nren,parts_graph)

       ! Create objects on the interfaces
       call objects_create (dual_femesh, ldomn, ldome, nren, max_nparts, &
          &                 nobjs, lobjs, stree)

       ! Create list of objects of each part 
       call part_objects_create(prt_pars%nparts,nobjs,lobjs,part_objs)

       ! Create list of objects on each edge of the graph of parts 
       call int_objects_create(nobjs,lobjs,parts_graph,int_objs)

       call memalloc (prt_pars%nparts+1, nedom,__FILE__,__LINE__)
       call memalloc (prt_pars%nparts+1, npdom,__FILE__,__LINE__)
       call parts_sizes(zero,dual_femesh,nren,ldome,prt_pars%nparts,max_nparts, &
          &             nobjs,lobjs,part_objs,npdom,nedom)

       ! Allocate and compute local to global maps
       allocate(nmaps(prt_pars%nparts))
       allocate(emaps(prt_pars%nparts))
       do ipart=1,prt_pars%nparts
          call map_alloc(npdom(ipart+1),femesh%npoin,parts(ipart)%nmap)
          call map_alloc(nedom(ipart+1),int(femesh%nelem,igp),parts(ipart)%emap)
          !call map_alloc(npdom(ipart+1),femesh%npoin,nmaps(ipart))
          !call map_alloc(nedom(ipart+1),femesh%nelem,emaps(ipart))
       end do
       call parts_maps(zero,dual_femesh,nren,eren,ldome,prt_pars%nparts,max_nparts, &
          &            nobjs,lobjs,part_objs,npdom,nedom,nmaps,emaps)
       do ipart=1,prt_pars%nparts
          parts(ipart)%nmap = nmaps(ipart)
          ! AFM parts(ipart)%emap = emaps(ipart)
          call map_free(nmaps(ipart))
          call map_free(emaps(ipart))
       end do

       ! Create map for boundary elements
       if(femesh%nboun>0) then
          call memalloc (prt_pars%nparts+1, nbdom, __FILE__,__LINE__)
          call parts_sizes_boundary(zero,femesh,ldome,ldomn,prt_pars%nparts,nbdom)

          ! Allocate and compute local to global map
          allocate(bmaps(prt_pars%nparts))
          do ipart=1,prt_pars%nparts
             if(nbdom(ipart+1).gt.0) then
                call map_alloc(nbdom(ipart+1),femesh%nboun,parts(ipart)%bmap)
                call map_alloc(nbdom(ipart+1),femesh%nboun,bmaps(ipart))
             end if
          end do
          call parts_maps_boundary(zero,femesh,ldome,ldomn,prt_pars%nparts,nbdom,bmaps)
          do ipart=1,prt_pars%nparts
             parts(ipart)%bmap = bmaps(ipart)
             call map_free(bmaps(ipart))
          end do
       end if

       if(prt_pars%use_graph==dual) then
          ! Computation of external elements information.
          ! To this end we need to regenerate the dual graph but now with
          ! ncommon=1 (which is implicitly the case when calling mesh_to_graph
          ! with this interface, which includes edges).
          call fem_graph_free (fe_graph)
          call mesh_to_graph (dual_femesh, femesh, fe_graph, edges)
          do ipart=1,prt_pars%nparts
!!$            call neighbors_graph_extract(parts(ipart)%emap,femesh%nelem,ldome,eren,fe_graph,edges, &
!!$                 &                       parts(ipart)%pextn,parts(ipart)%lextn,parts(ipart)%lextp, &
!!$                 &                       parts(ipart)%lexte,parts(ipart)%lextm,material)
          end do
          call memfree (edges,__FILE__,__LINE__)
       end if

    else if(prt_pars%ptype==vertex_based) then

       if(prt_pars%use_graph==primal) then

          !primal graph (points around points)
          call mesh_to_graph (femesh, dual_femesh, 1, fe_graph)

          call graph_pt_renumbering(prt_pars,fe_graph, nren, ldomn)

          call part_to_sep(stree, femesh, ldomn, ldome)

          call sep_tree_set_pointers(femesh%nelem,ldome,stree)

          call renum_by_sets(femesh%nelem,ldome,eren)

       else if(prt_pars%use_graph==dual) then

          if(stree%nlevel==0) then
             write(*,*) 'Graph partition by nested disection only possible when'
             write(*,*) 'nparts=2^k for some integer k.'
             stop
          end if

          ! dual graph (elements around elements)
          call mesh_to_graph (dual_femesh, femesh, 1, fe_graph)
          ! call graph_write(6, fe_graph) ! DBG: 

          call graph_nd_renumbering(prt_pars,fe_graph, stree%nlevel, eren, ldome)

          call sep_tree_set_pointers(femesh%nelem, ldome, stree)

          call sep_to_part(stree,dual_femesh,femesh,ldome,ldomn)

          call renum_by_sets(femesh%npoin,ldomn,nren)

       end if

       ! Create graph of parts
       call parts_graph_create(femesh, stree, ldomn, eren, parts_graph)

       ! Create objects on the interfaces
       call objects_create (femesh,ldome,ldomn,eren,max_nparts, &
          &                 nobjs,lobjs,stree)

       ! Create list of objects of each part 
       call part_objects_create (prt_pars%nparts,nobjs,lobjs,part_objs)

       ! Create list of objects on each edge of the graph of parts 
       call int_objects_create  (nobjs,lobjs,parts_graph,int_objs)

       call memalloc (prt_pars%nparts+1, nedom,__FILE__,__LINE__)
       call memalloc (prt_pars%nparts+1, npdom,__FILE__,__LINE__)
       call parts_sizes(one,femesh,eren,ldomn,prt_pars%nparts,max_nparts, &
          &             nobjs,lobjs,part_objs,nedom,npdom)

       ! Allocate and compute local to global maps
       allocate(nmaps(prt_pars%nparts))
       allocate(emaps(prt_pars%nparts))
       do ipart=1,prt_pars%nparts
          call map_alloc(npdom(ipart+1),femesh%npoin,parts(ipart)%nmap)
          call map_alloc(nedom(ipart+1),int(femesh%nelem,igp),parts(ipart)%emap)
          !call map_alloc(npdom(ipart+1),femesh%npoin,nmaps(ipart))
          !call map_alloc(nedom(ipart+1),femesh%nelem,emaps(ipart))
       end do
       call parts_maps(one,femesh,eren,nren,ldomn,prt_pars%nparts,max_nparts, &
          &            nobjs,lobjs,part_objs,nedom,npdom,emaps,nmaps)
       do ipart=1,prt_pars%nparts
          parts(ipart)%nmap = nmaps(ipart)
          ! AFM parts(ipart)%emap = emaps(ipart)
          call map_free(nmaps(ipart))
          call map_free(emaps(ipart))
       end do

       ! Create map for boundary elements
       if(femesh%nboun>0) then
          call memalloc (prt_pars%nparts+1, nbdom, __FILE__,__LINE__)
          call parts_sizes_boundary(one,femesh,ldome,ldomn,prt_pars%nparts,nbdom)

          ! Allocate and compute local to global map
          allocate(bmaps(prt_pars%nparts))
          do ipart=1,prt_pars%nparts
             if(nbdom(ipart+1).gt.0) then
                call map_alloc(nbdom(ipart+1),femesh%nboun,parts(ipart)%bmap)
                call map_alloc(nbdom(ipart+1),femesh%nboun,bmaps(ipart))
             end if
          end do
          call parts_maps_boundary(one,femesh,ldome,ldomn,prt_pars%nparts,nbdom,bmaps)
          do ipart=1,prt_pars%nparts
             parts(ipart)%bmap = bmaps(ipart)
             call map_free(bmaps(ipart))
          end do
       end if

    end if

    !call sep_tree_write(6, stree)

    if(prt_pars%debug>0) then 
       call fem_mesh_partition_outinfo(6,prt_pars%nparts,npdom,nedom,stree,parts_graph, &
            &                          nobjs,lobjs,part_objs,int_objs)

!!$       call fem_mesh_partition_outinfo_bddc_matlab_code (prt_pars%nparts,npdom,nedom,stree,parts_graph, &
!!$    &                                      nobjs,lobjs,part_objs,int_objs)

!!$       call fem_mesh_partition_outobjects(nren,eren,femesh,max_nparts,nobjs,lobjs)

    end if

    call fem_partition_g2l(prt_pars%ptype,prt_pars%nparts,max_nparts,nobjs,parts_graph, &
       &                   lobjs,part_objs,int_objs,parts)
    do ipart=1,prt_pars%nparts
       parts(ipart)%pinfo = interfaces
    end do

    call memfree(ldomn,__FILE__,__LINE__)
    call memfree(ldome,__FILE__,__LINE__)

    call memfree(nedom,__FILE__,__LINE__)
    call memfree(npdom,__FILE__,__LINE__)
    call memfree(lobjs,__FILE__,__LINE__)

    call memfree(part_objs%p,__FILE__,__LINE__)
    call memfree(part_objs%l,__FILE__,__LINE__)
    call memfree(int_objs%p,__FILE__,__LINE__)
    call memfree(int_objs%l,__FILE__,__LINE__)

    call sep_tree_free  (stree) 
    call fem_graph_free (parts_graph)

    call fem_mesh_free  (dual_femesh)
    call fem_graph_free (fe_graph)

  end subroutine fem_mesh_partition_create_old

  !================================================================================================
  function global_to_local(ipoin,ipart,nparts,nobjs,max_nparts,pdobj,ldobj,lobjs)
    implicit none
    integer(ip), intent(in)  :: ipoin,ipart,nparts,nobjs,max_nparts
    integer(ip), intent(in)  :: pdobj(nparts+1)
    integer(ip), intent(in)  :: ldobj(pdobj(nparts+1))
    integer(ip), intent(in)  :: lobjs(max_nparts+4,nobjs)
    integer(ip)              :: global_to_local
    integer(ip)              :: j,jpoin

    jpoin=0
    j=pdobj(ipart)
    do while( ipoin>lobjs(3,ldobj(j)) )
       jpoin=jpoin+lobjs(3,ldobj(j))-lobjs(2,ldobj(j))+1
       j=j+1
    end do
    global_to_local=jpoin+ipoin-lobjs(2,ldobj(j))+1

  end function global_to_local

  !================================================================================================
  subroutine fem_partition_outobjects(nren,eren,femesh,max_nparts,nobjs,lobjs)
    ! To write objects in GiD
    use stdio
    use fem_mesh_io
    use post
    implicit none
    ! Parameters
    integer(ip)        , intent(in) :: max_nparts, nobjs
    type(renum)        , intent(in) :: nren
    type(renum)        , intent(in) :: eren
    type(fem_mesh)     , intent(in) :: femesh
    integer(ip)        , intent(in) :: lobjs(max_nparts+4,nobjs)
    ! Locals
    character(len=256)           :: name,preord
    integer(ip)                  :: lunio,iobj
    integer(ip)   , allocatable  :: ws_lobj_field(:)
    type(post_file)              :: lupos

    call memalloc ( femesh%npoin, ws_lobj_field,__FILE__,__LINE__ )

    ! Define the object each node belongs
    do iobj=1,nobjs
       ws_lobj_field(lobjs(2,iobj):lobjs(3,iobj))=iobj
    end do

    ! Write a mesh
    preord = 'objects'
    call postpro_compose_mesh_name( preord, name ) 
    lunio = io_open(name)
    call fem_mesh_write(lunio,femesh,nren,eren)
    call io_close(lunio)

    ! Write objects as a postprocessed field
    call postpro_compose_result_name( preord, name ) 
    call postpro_open_file(1,name,lupos)
    call postpro(lupos,ws_lobj_field,'OBJS',1,1.0)
    call postpro_close_file(lupos)

    ! Free ws
    call memfree ( ws_lobj_field,__FILE__,__LINE__)

  end subroutine fem_partition_outobjects

  !================================================================================================
  subroutine fem_partition_g2l(ptype,nparts,max_nparts,nobjs,parts_graph, &
       &                       lobjs,part_objs,int_objs,parts)
    implicit none
    ! Parameters
    integer(ip)        , intent(in)    :: nparts, ptype, max_nparts, nobjs
    type(fem_graph)    , intent(in)    :: parts_graph
    integer(ip)        , intent(in)    :: lobjs(max_nparts+4,nobjs)
    type(list)         , intent(in)    :: part_objs
    type(list)         , intent(in)    :: int_objs
    type(fem_partition), intent(inout) :: parts(nparts)

    ! Locals
    integer(ip)                  :: ipart
    integer(ip)                  :: pi_obj, igobj, ilobj, igedge, iledge
    integer(ip)   , allocatable  :: ws_omap_g2l(:)

    call memalloc ( nobjs, ws_omap_g2l,__FILE__,__LINE__ )

    ! Distribute global information of the 
    ! mesh partitioning among each part
    do ipart=1,nparts
       parts(ipart)%ptype      = ptype
       parts(ipart)%ipart      = ipart
       parts(ipart)%nparts     = nparts
       parts(ipart)%max_nparts = max_nparts
       parts(ipart)%nobjs      = part_objs%p(ipart+1)-part_objs%p(ipart)

       ! Extract the number of neighbours
       parts(ipart)%npadj = parts_graph%ia(ipart+1)-parts_graph%ia(ipart)

       ! Allocate memory for the list of neighbours
       call memalloc ( parts(ipart)%npadj, parts(ipart)%lpadj,   __FILE__,__LINE__ )

       ! Extract the list of neighbours
       parts(ipart)%lpadj(1:parts(ipart)%npadj) = &
            & parts_graph%ja(parts_graph%ia(ipart):parts_graph%ia(ipart+1)-1)

       ! write(*,*) max_nparts+4, part_objs%p(ipart+1)-part_objs%p(ipart)
       ! Allocate memory for the list of local objects
       call memalloc ( max_nparts+4, part_objs%p(ipart+1)-part_objs%p(ipart),   parts(ipart)%lobjs, __FILE__,__LINE__)

       ! Allocate memory for the local to global mapping of the objects
       call map_alloc( part_objs%p(ipart+1)-part_objs%p(ipart), nobjs, parts(ipart)%omap ) 

       ! Allocate memory for the pointers to the lists of objects on each edge
       call memalloc ( parts(ipart)%npadj+1, parts(ipart)%int_objs%p,  __FILE__,__LINE__ )


       ! Set the pointers to the list of local objects on each edge
       parts(ipart)%int_objs%n = parts(ipart)%npadj
       parts(ipart)%int_objs%p(1) = 1
       iledge = 2
       do igedge = parts_graph%ia(ipart), parts_graph%ia(ipart+1)-1
          parts(ipart)%int_objs%p(iledge) =  parts(ipart)%int_objs%p(iledge-1) + &
               &  int_objs%p (igedge+1) - int_objs%p (igedge)
          iledge = iledge + 1
       end do

       ! Allocate memory for the lists of objects on each edge
       call memalloc ( parts(ipart)%int_objs%p(parts(ipart)%int_objs%n+1)-1,        parts(ipart)%int_objs%l,  __FILE__,__LINE__ )

       ! Copy the list of objects on each edge (objects id's have 
       ! still to be translated from global to local)
       iledge = 1
       do igedge = parts_graph%ia(ipart), parts_graph%ia(ipart+1)-1
          ! Extract the list of objects on igedge
          parts(ipart)%int_objs%l( parts(ipart)%int_objs%p(iledge): &
               &                   parts(ipart)%int_objs%p(iledge+1)-1) = &
               &       int_objs%l(int_objs%p(igedge):int_objs%p(igedge+1)-1)
          iledge = iledge + 1
       end do

       ! Set number of internal, boundary and external objects
       parts(ipart)%omap%ni = 1
       parts(ipart)%omap%nb = parts(ipart)%omap%nl - parts(ipart)%omap%ni
       parts(ipart)%omap%ne = 0      


       ! Traverse objects on part ipart
       ilobj = 1
       do pi_obj = part_objs%p(ipart), part_objs%p(ipart+1)-1
          igobj = part_objs%l(pi_obj)

          ! Set local to global object mapping
          parts(ipart)%omap%l2g ( ilobj ) = igobj
          ws_omap_g2l ( igobj ) = ilobj


          ! Copy object from the global list of objects to the local one
          parts(ipart)%lobjs(:, ilobj) = lobjs(:, igobj)

          ! Transform start/end pointers from global to local space
          parts(ipart)%lobjs(2, ilobj) =  &
               &   global_to_local( lobjs(2, igobj), ipart, nparts, nobjs, &
               &                    max_nparts, part_objs%p, part_objs%l, lobjs)

          ! Transform start/end pointers from global to local space
          parts(ipart)%lobjs(3, ilobj) = &
               &   global_to_local( lobjs(3, igobj), ipart, nparts, nobjs, &
               &                    max_nparts, part_objs%p, part_objs%l, lobjs)

          ilobj = ilobj + 1
       end do

       ! Translate object id's from global to local in
       ! the list of objects on each edge
       do iledge = 1, parts(ipart)%int_objs%p(parts(ipart)%npadj+1)-1
          ! Extract the list of objects on igedge
          parts(ipart)%int_objs%l( iledge ) =  ws_omap_g2l ( parts(ipart)%int_objs%l( iledge ) )
       end do

    end do

    call memfree ( ws_omap_g2l,__FILE__,__LINE__)

  end subroutine fem_partition_g2l

!!$  !================================================================================================
!!$  subroutine objects_hierarchy_create (lunou,ptype,nparts,max_nparts,nobjs,npdom,nedom,part_objs, &
!!$       &                               parts_graph,int_objs,lobjs,stree,parts,num_nested,list_nested)
!!$    implicit none
!!$    ! Parameters
!!$    integer(ip)        , intent(in)    :: lunou,nparts, ptype, max_nparts, nobjs
!!$    integer(ip)        , intent(in)    :: nedom(:)    ! Number of points per part
!!$    integer(ip)        , intent(in)    :: npdom(:)    ! Number of elements per part
!!$    type(fem_graph)    , intent(inout) :: parts_graph
!!$    integer(ip)        , intent(in)    :: lobjs(max_nparts+4,nobjs)
!!$    type(list)         , intent(inout) :: part_objs
!!$    type(list)         , intent(inout) :: int_objs
!!$    type(sep_tree)     , intent(in)    :: stree
!!$    type(fem_partition), intent(inout) :: parts(nparts,num_nested)
!!$    integer(ip)        , intent(in)    :: num_nested
!!$    integer(ip)        , intent(in)    :: list_nested(num_nested+1)
!!$
!!$    ! Locals
!!$    integer(ip)   , allocatable  :: parts_groups(:),parts_subgroups(:)
!!$    integer(ip)   , allocatable  :: groups_id(:),subgroups_id(:)
!!$    integer(ip)   , allocatable  :: lobjs_subgroups(:,:)
!!$    integer(ip)   , allocatable  :: lgrouped_objs(:,:)
!!$    integer(ip)   , allocatable  :: lgrouped_objs_num(:)
!!$    integer(ip)   , allocatable  :: ws_objs_num(:)
!!$    integer(ip)   , allocatable  :: ws_sort_l1(:)
!!$    integer(ip)   , allocatable  :: ws_sort_l2(:)
!!$    type(fem_graph)              :: parts_graph_grouped
!!$    type(list)                   :: part_objs_grouped
!!$    type(list)                   :: int_objs_grouped
!!$
!!$    integer(ip)                  :: ipart,iobj,inest,inode,igroup,i
!!$    integer(ip)                  :: ngroups,nsubgroups,isubgroup,jsubgroup
!!$
!!$    integer(ip)                  :: pi_obj, igobj, ilobj, igedge, iledge
!!$    integer(ip)   , allocatable  :: ws_omap_g2l(:)
!!$
!!$    ! An old stuff to change....some day...
!!$    !integer  icomp
!!$    !external icomp
!!$
!!$    call memalloc ( nparts, parts_groups    ,__FILE__,__LINE__ )
!!$    call memalloc ( nparts, parts_subgroups ,__FILE__,__LINE__ )
!!$    ! Always ngroups<=nparts (and nsubgroups<=nparts, of course)
!!$    call memalloc ( nparts, groups_id    ,__FILE__,__LINE__ )
!!$    call memalloc ( nparts, subgroups_id ,__FILE__,__LINE__ )
!!$
!!$    call memalloc ( max_nparts+1, nobjs, lobjs_subgroups,__FILE__,__LINE__ )
!!$    call memalloc ( max_nparts+4, nobjs, lgrouped_objs,__FILE__,__LINE__)
!!$    call memalloc ( nobjs, lgrouped_objs_num,__FILE__,__LINE__)
!!$    call memalloc ( nobjs, ws_objs_num,__FILE__,__LINE__)
!!$    call memalloc ( max_nparts+1, ws_sort_l1,__FILE__,__LINE__)
!!$    call memalloc ( max_nparts+1, ws_sort_l2,__FILE__,__LINE__)
!!$
!!$    inest=1
!!$    do inest=1,num_nested
!!$       ngroups=0
!!$       nsubgroups=0
!!$       groups_id=0
!!$       subgroups_id=0
!!$       parts_groups=0
!!$       parts_subgroups=0
!!$       do inode=1,stree%nnode
!!$          if(stree%nodes(inode)%level==list_nested(inest)) then
!!$             ! Set groups
!!$             ngroups=ngroups+1
!!$             igroup=stree%nodes(inode)%part
!!$             groups_id(ngroups)=igroup
!!$             call sep_tree_set_group(nparts,stree%nnode,igroup,inode,parts_groups,stree%nodes)
!!$          else if(stree%nodes(inode)%level==list_nested(inest+1)) then
!!$             ! Set subgroups
!!$             nsubgroups=nsubgroups+1
!!$             igroup=stree%nodes(inode)%part
!!$             subgroups_id(nsubgroups)=igroup
!!$             call sep_tree_set_group(nparts,stree%nnode,igroup,inode,parts_subgroups,stree%nodes)
!!$          end if
!!$       end do
!!$       ! DBG:
!!$       write(*,'(a,10i10)') 'Nested level', inest
!!$       write(*,'(a,10i10)') 'Number of groups', ngroups
!!$       write(*,'(a,10i10)') 'Number of subgroups', nsubgroups
!!$          write(*,'(a,10i10)') 'Groups id ', groups_id
!!$          write(*,'(a,10i10)') 'parts_groups',parts_groups
!!$          write(*,'(a,10i10)') 'Subgroups id ', subgroups_id
!!$          write(*,'(a,10i10)') 'parts_subgroups',parts_subgroups
!!$
!!$       ! For each group a partition needs to be generated
!!$       ! taking (only) subgroups communication into account
!!$       do igroup=1,ngroups
!!$          write(*,*) 'Group', groups_id(igroup)
!!$          write(*,*) 'Part  subgroup'
!!$          do ipart=1,nparts
!!$             if(parts_groups(ipart) == groups_id(igroup)) write(*,*) ipart, parts_subgroups(ipart)
!!$          end do
!!$
!!$          write(*,*) 'Objects'
!!$          lobjs_subgroups=0
!!$          lgrouped_objs=0
!!$          do iobj=1,nobjs
!!$             i=0
!!$             isubgroup = 0
!!$             ! Copy parts connection in lobjs (last max_part+1 columns)
!!$             ! into lobjs_subgroups replacing parts by subgroups
!!$             ! and eliminating connections outside the current 
!!$             ! group (whose id is groups_id(igroup)). Also copy lobjs into
!!$             ! lgrouped_objs eliminating connections outside the current
!!$             ! group (but keeping all the rest).
!!$             lgrouped_objs(1:3,iobj) = lobjs(1:3,iobj)
!!$             do ipart=1,lobjs(4,iobj)
!!$                if( parts_groups( lobjs(4+ipart,iobj) ) == groups_id(igroup)) then
!!$                   i=i+1
!!$                   lobjs_subgroups(1+i,iobj) = parts_subgroups( lobjs(4+ipart,iobj) )
!!$                   lgrouped_objs(4+i,iobj) = lobjs(4+ipart,iobj)
!!$                end if
!!$             end do
!!$             lobjs_subgroups(1,iobj)=i
!!$             lgrouped_objs(4,iobj)=i
!!$             ! Now eliminate repeated subgroups in lobjs_subgroups (which come
!!$             ! from the communication of parts inside the same subgroup).
!!$             if(i>1) then
!!$                call sort_eliminate_repeated(i,lobjs_subgroups(1,iobj),lobjs_subgroups(2:1+i,iobj))
!!$                lobjs_subgroups(2+lobjs_subgroups(1,iobj):max_nparts+1,iobj)=0
!!$             end if
!!$             ! DBG:
!!$             !write(*,'(10i10)') iobj,lobjs_subgroups(:,iobj)
!!$             ! END DBG:
!!$          end do
!!$
!!$          ! Set the id numbering and call sort procedure
!!$          do iobj=1,nobjs
!!$             ws_objs_num(iobj)=iobj
!!$          end do
!!$          call intsort(max_nparts+1,max_nparts+1,nobjs,lobjs_subgroups,ws_objs_num,ws_sort_l1, ws_sort_l2)
!!$          ! DBG:
!!$          !write(*,*) 'Sort result'
!!$          !do iobj=1,nobjs
!!$          !   write(*,'(10i10)') iobj,ws_objs_num(iobj),lobjs_subgroups(:,iobj)
!!$          !end do
!!$          ! END_DBG:
!!$
!!$          i=1
!!$          lgrouped_objs_num=0
!!$          if(lobjs_subgroups(1,nobjs)>1) then
!!$             lgrouped_objs_num(ws_objs_num(nobjs))=i
!!$          else
!!$             lgrouped_objs_num(ws_objs_num(nobjs))=0
!!$          end if
!!$          do iobj=nobjs-1,1,-1
!!$             if(lobjs_subgroups(1,iobj)>1) then
!!$                ! This if is the one generating the same number for "coarse" objects
!!$                if(icomp(max_nparts+1,lobjs_subgroups(:,iobj+1),lobjs_subgroups(:,iobj))/=0) then
!!$                   i=i+1
!!$                end if
!!$                lgrouped_objs_num(ws_objs_num(iobj))=i
!!$             else
!!$                ! Set number of subdomains to 0
!!$                lgrouped_objs(4,ws_objs_num(iobj))=0
!!$             end if
!!$          end do
!!$          ! DBG:
!!$          !write(*,*) 'Final numbering and subgroups'
!!$          !do iobj=1,nobjs
!!$          !   write(*,'(10i10)') iobj,lgrouped_objs_num(ws_objs_num(iobj)),lobjs_subgroups(:,iobj)
!!$          !end do
!!$          ! END_DBG:
!!$
!!$          ! Print objects with global id
!!$          ! DBG:
!!$          !write(*,*) 'Final numbering of objects'
!!$          !do iobj=1,nobjs
!!$          !   write(*,'(10i10)') iobj,lgrouped_objs_num(iobj),lobjs(:,iobj)
!!$          !end do
!!$          write(*,*) 'Final numbering of grouped objects'
!!$          do iobj=1,nobjs
!!$             write(*,'(10i10)') iobj,lgrouped_objs_num(iobj),lgrouped_objs(:,iobj)
!!$          end do
!!$          write(*,*) 'End objects'
!!$          ! END_DBG:
!!$
!!$          ! Create (reduced) graph of parts 
!!$          call parts_graph_create_from_objects (nparts,max_nparts,nobjs,lgrouped_objs, &
!!$               &                                parts_graph_grouped)
!!$          ! Create list of objects of each part  (using the same space)
!!$          call part_objects_create (nparts,nobjs,lgrouped_objs,part_objs_grouped)
!!$          ! Create list of objects on each edge of the graph of parts (using the same space)
!!$          call int_objects_create  (nobjs,lgrouped_objs,parts_graph_grouped,int_objs_grouped)
!!$          if(lunou>0) &
!!$               &  call fem_mesh_partition_outinfo(lunou,nparts,npdom,nedom,stree,parts_graph_grouped, &
!!$               &                                  nobjs,lgrouped_objs,part_objs_grouped,int_objs_grouped)
!!$
!!$          call fem_graph_free (parts_graph_grouped)
!!$          call memfree(part_objs_grouped%p,__FILE__,__LINE__)
!!$          call memfree(part_objs_grouped%l,__FILE__,__LINE__)
!!$          call memfree(int_objs_grouped%p,__FILE__,__LINE__)
!!$          call memfree(int_objs_grouped%l,__FILE__,__LINE__)
!!$
!!$       end do
!!$
!!$    end do
!!$
!!$    call memfree ( parts_groups   ,__FILE__,__LINE__)
!!$    call memfree ( parts_subgroups,__FILE__,__LINE__)
!!$
!!$    call memfree ( groups_id   ,__FILE__,__LINE__)
!!$    call memfree ( subgroups_id,__FILE__,__LINE__)
!!$
!!$    call memfree ( lobjs_subgroups,__FILE__,__LINE__)
!!$    call memfree ( lgrouped_objs,__FILE__,__LINE__)
!!$    call memfree ( lgrouped_objs_num,__FILE__,__LINE__)
!!$    call memfree ( ws_sort_l1,__FILE__,__LINE__)
!!$    call memfree ( ws_sort_l2,__FILE__,__LINE__)
!!$
!!$  end subroutine objects_hierarchy_create

  !================================================================================================
  subroutine parts_graph_create_from_objects (nparts,max_nparts,nobjs,lobjs,p_graph)
    implicit none
    integer(ip)    , intent(in)    :: nparts, max_nparts, nobjs
    type(fem_graph), intent(inout) :: p_graph
    integer(ip)    , intent(in)    :: lobjs(max_nparts+4,nobjs)
    integer(ip)    , allocatable   :: iwork(:),jwork(:)
    integer(ip) :: max_size,iobj

    ! An upper bound to the number of edges
    max_size=0
    do iobj=1,nobjs
       max_size=max_size+lobjs(4,iobj)*(lobjs(4,iobj)-1)
    end do

    call memalloc (nparts+1, iwork, __FILE__,__LINE__)
    call memalloc (max_size, jwork, __FILE__,__LINE__)

    call parts_graph_create_from_objects_aux (nparts,max_nparts,nobjs,lobjs,max_size,iwork,jwork)

    ! Alloc and fill values
    call memalloc (nparts+1, p_graph%ia, __FILE__,__LINE__)
    call memalloc(iwork(nparts+1)-1, p_graph%ja,__FILE__,__LINE__)
    p_graph%nv=nparts
    p_graph%ia=iwork
    p_graph%ja=jwork(iwork(1):iwork(nparts+1)-1)

    call memfree (iwork,__FILE__,__LINE__)
    call memfree (jwork,__FILE__,__LINE__)

    return

  end subroutine parts_graph_create_from_objects

  subroutine parts_graph_create_from_objects_aux (nparts,max_nparts,nobjs,lobjs, &
       &                                          max_size,iwork,jwork)
    implicit none
    integer(ip)        , intent(in)    :: nparts, max_nparts, nobjs, max_size
    integer(ip)        , intent(in)    :: lobjs(max_nparts+4,nobjs)
    integer(ip)        , intent(inout) :: iwork(nparts+1)
    integer(ip)        , intent(inout) :: jwork(max_size)
    integer(ip)                        :: ipart,iobj,i,j,k,l,m,n,next,nz

    ! Count (repeteaded) edges 
    iwork=0
    do iobj=1,nobjs
       do i=1,lobjs(4,iobj)
          iwork(lobjs(4+i,iobj)+1)=iwork(lobjs(4+i,iobj)+1)+lobjs(4,iobj)-1
       end do
    end do

    ! Compress iwork (pointing to repeated edges)
    iwork(1)=1
    do ipart=1,nparts
       iwork(ipart+1)=iwork(ipart+1)+iwork(ipart)
    end do

    ! Store (repeteaded) edges in jwork
    do iobj=1,nobjs
       do i=1,lobjs(4,iobj)
          do j=1,i-1
             jwork(iwork(lobjs(4+i,iobj)))=lobjs(4+j,iobj)
             iwork(lobjs(4+i,iobj))=iwork(lobjs(4+i,iobj))+1
          end do
          do j=i+1,lobjs(4,iobj)
             jwork(iwork(lobjs(4+i,iobj)))=lobjs(4+j,iobj)
             iwork(lobjs(4+i,iobj))=iwork(lobjs(4+i,iobj))+1
          end do
       end do
    end do

    ! Recover iwork
    do ipart=nparts+1, 2, -1
       iwork(ipart) = iwork(ipart-1)
    end do
    iwork(ipart) = 1

    ! Eliminate repeated edges
    next=1
    do ipart=1,nparts
       k=iwork(ipart)
       n=iwork(ipart+1)-iwork(ipart)
       iwork(ipart)=next
       if(n>0) then
          call sort_eliminate_repeated(n,m,jwork(k:k+n-1))
          next=next+m
          jwork(iwork(ipart):next-1)=jwork(k:k+n-1)
       end if
    end do
    iwork(nparts+1)=next

  end subroutine parts_graph_create_from_objects_aux

  !================================================================================================
  subroutine fem_mesh_partition_outinfo(lu_out,nparts,npdom,nedom,stree,parts_graph, &
     &                                  nobjs,lobjs,part_objs,int_objs)
    implicit none
    integer(ip)    , intent(in) :: lu_out
    integer(ip)    , intent(in) :: nparts
    type(sep_tree) , intent(in) :: stree       ! Separator_tree
    type(fem_graph), intent(in) :: parts_graph ! Parts graph
    integer(ip)    , intent(in) :: nobjs       ! number of objects
    integer(ip)    , intent(in) :: lobjs(:,:)  ! List of objects
    type(list)                  :: part_objs   ! Currently pdobj,ldobj
    type(list)                  :: int_objs    ! Currently piobj,liobj
    integer(ip)   , intent(in)  :: nedom(:)    ! Number of points per part
    integer(ip)   , intent(in)  :: npdom(:)    ! Number of elements per part
    integer(ip) :: i,j

    if(lu_out>0) then

       write(lu_out,*)
       write(lu_out,'(a)') 'Nodes and elements distribution:'
       write(lu_out,'(a)') '      Part   1st_node  lst_node  1st_elem  lst_elem'
       do i=1,nparts
          write(lu_out,'(5(i8,2x))') i,npdom(i),npdom(i+1), &
             &                         nedom(i),nedom(i+1)
       end do
       write(lu_out,*)
       write(lu_out,'(a)') 'Vertices on each node of the separator tree:'
       write(lu_out,'(10i10)') stree%node_ptrs
       write(lu_out,*)

       write(lu_out,'(a15)')     'Parts graph:'
       write(lu_out,'(a18,i10)') 'Number of edges:', &
          &  parts_graph%ia(parts_graph%nv+1)-1
       do i=1,parts_graph%nv
          write(lu_out,'(10i10)') i,parts_graph%ja( &
             &  parts_graph%ia(i):parts_graph%ia(i+1)-1)
       end do

       write(lu_out,'(a18,i10)') 'Number of elements on the boundary:', &
          &  parts_graph%ia(parts_graph%nv+1)-1
       do i=1,parts_graph%nv
          write(lu_out,'(10i10)') i,parts_graph%ja( &
             &  parts_graph%ia(i):parts_graph%ia(i+1)-1)
       end do



       write(lu_out,'(a)') 'List of objects:'
       do i=1,nobjs
          write(lu_out,'(10i10)') i,lobjs(:,i)
       end do

       write(lu_out,'(a)') 'List of part objects:'
       do i=1,nparts
          write(lu_out,'(10i10)') i, &
             & (part_objs%l(j),j=part_objs%p(i),part_objs%p(i+1)-1)
       end do

       !write(lu_out,'(a)') 'Pointers to interface objects:'
       !do i=1,parts_graph%ia(parts_graph%nv+1)
       !   write(lu_out,*) int_objs%p(i)
       !end do
       write(lu_out,'(a)') 'List of interface objects:'
       do i=1,parts_graph%ia(parts_graph%nv+1)-1
          write(lu_out,'(10i10)') i, &
             & (int_objs%l(j),j=int_objs%p(i),int_objs%p(i+1)-1)
       end do

    end if
    
  end subroutine fem_mesh_partition_outinfo
   
  !================================================================================================
  subroutine fem_mesh_partition_from_periodic (nparts, uparts, ueren, gmesh, gnren, gnmaps) 
    !-----------------------------------------------------------------------
    ! This routine computes the minimum data required (i.e., gnmaps+gneren) 
    ! in order to build the local parts of the original mesh (i.e., gmesh)
    ! from a partition/renumeration of the periodic mesh (i.e., umesh) given 
    ! on uparts+ueren
    !-----------------------------------------------------------------------
    implicit none  
    ! Parameters
    type(fem_mesh)        , intent(in)  :: gmesh
    integer(ip)           , intent(in)  :: nparts 
    type(fem_partition)   , intent(in)  :: uparts(nparts) ! partition of the periodic mesh
    type(renum)           , intent(in)  :: ueren
    type(renum)           , intent(out) :: gnren         ! node re-numbering for gmesh
    type(map), allocatable, intent(out) :: gnmaps(:)      ! node maps for gmesh

    ! Locals
    integer(ip)    , allocatable :: work(:) 
    integer(ip)                  :: ielpo,ipart, lielem, gielem, inods1, inods2,ipoin, k

    allocate(gnmaps(nparts))
    call memalloc (gmesh%npoin, work, __FILE__,__LINE__)

    ! Allocates and initializes permutation to identity
    ! This suffices in the current context
    call renum_alloc (gmesh%npoin, gnren)

    do ipart=1,nparts
       work = 0
       k    = 0
       do lielem=1, uparts(ipart)%emap%nl
          gielem = uparts(ipart)%emap%l2g(lielem)
          gielem = ueren%iperm(gielem)
          if(gmesh%nelty==1) then
             inods1=(gielem-1)*gmesh%nnode+1
             inods2=gielem*gmesh%nnode
          else
             inods1=gmesh%pnods(gielem)
             inods2=gmesh%pnods(gielem+1)-1
          end if
          do ielpo = inods1,inods2
             ipoin=gmesh%lnods(ielpo)
             if (work(ipoin)==0) then
                k = k+1
                work(ipoin) = 1
             end if
          end do
       end do
       call map_alloc(k,gmesh%npoin,gnmaps(ipart))
    end do

    do ipart=1,nparts
       work = 0
       k    = 0
       do lielem=1, uparts(ipart)%emap%nl
          gielem = uparts(ipart)%emap%l2g(lielem)
          gielem = ueren%iperm(gielem)
          if(gmesh%nelty==1) then
             inods1=(gielem-1)*gmesh%nnode+1
             inods2=gielem*gmesh%nnode
          else
             inods1=gmesh%pnods(gielem)
             inods2=gmesh%pnods(gielem+1)-1
          end if
          do ielpo = inods1,inods2
             ipoin=gmesh%lnods(ielpo)
             if (work(ipoin)==0) then
                k = k+1
                gnmaps(ipart)%l2g(k) = ipoin
                work(ipoin) = 1
             end if
          end do
       end do
       ! write (*,*) 'XXX', ipart, gnmaps(ipart)%l2g DBG:
    end do
    
    call memfree (work,__FILE__,__LINE__)

  end subroutine fem_mesh_partition_from_periodic

  !================================================================================================
  subroutine fem_mesh_partition_outinfo_bddc_matlab_code (nparts,npdom,nedom,stree,parts_graph, &
     &                                                     nobjs,lobjs,part_objs,int_objs)
    use stdio
    implicit none
    ! Parameters
    integer(ip)    , intent(in) :: nparts
    type(sep_tree) , intent(in) :: stree       ! Separator_tree
    type(fem_graph), intent(in) :: parts_graph ! Parts graph
    integer(ip)    , intent(in) :: nobjs       ! number of objects
    integer(ip)    , intent(in) :: lobjs(:,:)  ! List of objects
    type(list)                  :: part_objs   ! Currently pdobj,ldobj
    type(list)                  :: int_objs    ! Currently piobj,liobj
    integer(ip)   , intent(in)  :: nedom(:)    ! Number of points per part
    integer(ip)   , intent(in)  :: npdom(:)    ! Number of elements per part

    ! Locals
    integer      :: fileo
    integer (ip) :: npsep, nodesep, idomn, i, j

    fileo = io_open ('partition_matlab.m', 'write')

    ! Compute npsep
    npsep=0
    do nodesep=1, stree%nnode
       write (*,*) stree%nodes(nodesep)%level, stree%nlevel, stree%node_ptrs(nodesep+1) - stree%node_ptrs(nodesep) 
       if (stree%nodes(nodesep)%lson /= 0 .or. stree%nodes(nodesep)%rson /= 0) then
          npsep = npsep + (stree%node_ptrs(nodesep+1) - stree%node_ptrs(nodesep))
       end if
    end do

    write(fileo,'(a)') 'function partition_matlab'
    write(fileo,*) 'npsep = ', npsep, ';'
    write(fileo,*) 'gmesh = [', stree%node_ptrs, '];'
    write(fileo,*)
    !
    write(fileo,'(a15,i10)') 'ndadj = ', parts_graph%ia(parts_graph%nv+1)-1 ,';'
    do idomn=1,nparts
       write(fileo,'(a6,i5,a2)',advance='no') 'jdadj{',idomn,'}='
       write(fileo,'(a1,10i10)') '[', idomn, &
            &  parts_graph%ja(parts_graph%ia(idomn):parts_graph%ia(idomn+1)-1)
       write(fileo,'(a2)') '];'
    end do
    !
    write(fileo,'(a)',advance='no') 'lobjs=['
    do i=1,nobjs
       write(fileo,'(a15,10i10)',advance='no') '[',i,lobjs(:,i)
       write(fileo,'(a15)') '];'
    end do
    write(fileo,'(a)') '];'
    !
    do i=1,nparts
       write(fileo,'(a6,i5,a2)',advance='no') 'ldobj{',i,'}='
       write(fileo,'(a1,20i10)',advance='no') '[', i, &
            & (part_objs%l(j),j=part_objs%p(i),part_objs%p(i+1)-1)
       write (fileo,'(a1)') ' ';
       write(fileo,'(a2)') '];'
    end do

    do i=1,parts_graph%ia(parts_graph%nv+1)-1
       write(fileo,'(a6,i5,a2)',advance='no') 'liobj{',i,'}='
       write(fileo,'(a1,10i10)',advance='no') '[', i, &
            & (int_objs%l(j),j=int_objs%p(i),int_objs%p(i+1)-1)
       write(fileo,'(a2)') '];'
    end do

    write(fileo,'(a)') 'npdom = ['
    do i=1,nparts
       write(fileo,'(2i10)') npdom(i+1)-1
    end do
    write(fileo,'(a)') '];'
    write(fileo,'(a)') 'nedom = ['
    do i=1,nparts
       write(fileo,'(2i10)') nedom(i+1)-1
    end do
    write(fileo,'(a)') '];'
    write(fileo,'(a)') 'save partition_matlab.mat'
    write(fileo,'(a)') 'end'
    call io_close (fileo)
  end subroutine fem_mesh_partition_outinfo_bddc_matlab_code

  !=================================================================================================
  !================================================================================================
  subroutine sep_to_part(stree, primal_mesh, dual_mesh, ldomn, ldome)
    !-----------------------------------------------------------------------
    ! Assign part to dual vertices based on separator of primal vertices,i.e.
    !
    ! Input:  ldomn
    ! Output: ldome
    !-----------------------------------------------------------------------
    implicit none
    type(sep_tree), intent(in)  :: stree
    type(fem_mesh), intent(in)  :: primal_mesh
    type(fem_mesh), intent(in)  :: dual_mesh
    integer(ip)   , intent(in)  :: ldomn(primal_mesh%npoin)
    integer(ip)   , intent(out) :: ldome(primal_mesh%nelem)
    integer(ip)                 :: inods1,inods2
    integer(ip)                 :: inods1d, inods2d  


    ! Local variables 
    integer(ip)                       :: ipoin,ielpo,ielem,untouched,jpoin, jelem, max_lev,inode

    ! I do not agree with the following
    ! assertion. dual_mesh can actually be
    ! a primal finite element mesh. This last
    ! situation was not taken into account in 
    ! the code below. I have modified the code 
    ! below so that it is also able to compute dual graphs.
    
    ! !!!! IMPORTANT NOTE: INCORRECT ASSERTION !!!
    ! assert(dual_mesh%nelty/=1)

    ! CORRECT ASSERTION: Either primal_mesh or dual_mesh has
    !                    to be a dual_mesh
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 

    ! Initialize ldome
    ldome = 0

    do ipoin=1,primal_mesh%npoin
       !if (stree%nodes(ldomn(ipoin))%level==stree%nlevel) then
       if (stree%nodes(ldomn(ipoin))%rson==0) then
          if(dual_mesh%nelty==1) then
             inods1d=(ipoin-1)*dual_mesh%nnode+1
             inods2d=ipoin*dual_mesh%nnode
          else
             inods1d=dual_mesh%pnods(ipoin)
             inods2d=dual_mesh%pnods(ipoin+1)-1
          end if
          do ielpo=inods1d,inods2d
          ! do ielpo=dual_mesh%pnods(ipoin),dual_mesh%pnods(ipoin+1)-1
             ielem=dual_mesh%lnods(ielpo)
             ldome(ielem) = stree%nodes(ldomn(ipoin))%part
          end do
       end if
    end do

    do ielem=1,primal_mesh%nelem
       if ( ldome(ielem) == 0 ) then
          ! TO-DO: This message should only be shown if a debug option
          ! is activated by the user within the params data structure
          write(*,*) '** [Fempar Warning] ** sep_to_part: element ', ielem, &
           &          ' untouched!!!! Starting second pass  ...'  ! DBG: 
          ! The solution required to define the following strategy for those
          ! elements which are not assigned by the previous loop

          ! Traverse nodes of ielem. Determine those which are closer 
          ! to the leaves in the separator tree. Among those which are 
          ! closer to the leaves, select one arbitrarily
          max_lev=0
          if(primal_mesh%nelty==1) then
             inods1=(ielem-1)*primal_mesh%nnode+1
             inods2=ielem*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ielem)
             inods2=primal_mesh%pnods(ielem+1)-1
          end if
          do inode = inods1,inods2
          !do inode=primal_mesh%pnods(ielem), primal_mesh%pnods(ielem+1)-1
             ipoin=primal_mesh%lnods(inode)
             ! write(*,*) ipoin, ielem, jpoin, max_lev, stree%nodes(ldomn(ipoin))%level ! DBG:
             if ( stree%nodes(ldomn(ipoin))%level > max_lev ) then
                max_lev = stree%nodes(ldomn(ipoin))%level
                jpoin   = ipoin
             end if
          end do

          ! write(*,*) ielem, jpoin, max_lev ! DBG: 

          ! Traverse elements around jpoin, and pick the first one  
          ! already mapped to a part
          do ielpo=dual_mesh%pnods(jpoin),dual_mesh%pnods(jpoin+1)-1
             jelem=dual_mesh%lnods(ielpo)
             if ( ldome(jelem) /= 0 ) then
                ldome(ielem) = ldome(jelem)
                ! TO-DO: This message should only be shown if a debug option
                ! is activated by the user within the params data structure
                write(*,*) '** [Fempar Warning] ** sep_to_part: element ', ielem, &
                  &         ' mapped to part #', ldome(jelem)  ! DBG: 
                exit
             end if
          end do
          ! write(*,*) ielem, ldome(ielem) DBG:

       end if
    end do

    ! TO-DO: This ERROR check should only be performed if a debug option
    ! is activated by the user within the params data structure
    do ielem=1,primal_mesh%nelem
       if ( ldome(ielem) == 0 ) then
          write(*,*) '** [Fempar Warning] ** sep_to_part: element ', ielem, &
             &       ' untouched!!!! Unexpected behaviour from here ...'  
       end if
    end do

  end subroutine sep_to_part

  !================================================================================================
  subroutine part_to_sep(stree,msh,ldome,ldomn)
    !-----------------------------------------------------------------------
    ! Assign separator to primal vertices based on part to dual vertices,i.e.
    !
    ! Input:  ldome
    ! Output: ldomn
    !
    ! Here mesh is dual (elements are nodes and viceversa).
    !-----------------------------------------------------------------------
    implicit none
    type(sep_tree), intent(in)  :: stree
    type(fem_mesh), intent(in)  :: msh
    integer(ip)   , intent(in)  :: ldome(msh%npoin)
    integer(ip)   , intent(out) :: ldomn(msh%nelem)
    integer(ip)                 :: ipoin,ielpo,ielem,idmin,idmax
    integer(ip)                 :: inods1,inods2

    do ielem=1,msh%nelem
       ! Find minimum and maximum tree node
       idmin=stree%nnode+1
       idmax=0

       if(msh%nelty==1) then
          inods1=(ielem-1)*msh%nnode+1
          inods2=ielem*msh%nnode
       else
          inods1=msh%pnods(ielem)
          inods2=msh%pnods(ielem+1)-1
       end if
       do ielpo = inods1,inods2
      !do ielpo=msh%pnods(ielem),msh%pnods(ielem+1)-1
          ipoin=msh%lnods(ielpo)
          idmin=min(idmin,stree%part_to_node(ldome(ipoin)))
          idmax=max(idmax,stree%part_to_node(ldome(ipoin)))
       end do
       ! Go up to the same label (incomplete trees)
       if(stree%nodes(idmin)%level<stree%nodes(idmax)%level) idmax=stree%nodes(idmax)%parent
       if(stree%nodes(idmin)%level>stree%nodes(idmax)%level) idmin=stree%nodes(idmin)%parent

       !do while(idmin/=idmax)
       !   if(stree%nodes(idmin)%level==stree%nodes(idmax)%level) idmin=stree%nodes(idmin)%parent
       !   idmax=stree%nodes(idmax)%parent
       !end do

!!$       ! Search separator tree to find a common ancestor
!!$       do while(idmin/=idmax)
!!$          idmin=stree%nodes(idmin)%parent
!!$          idmax=stree%nodes(idmax)%parent
!!$       end do

       ! DEBUG: Search separator tree to find a common ancestor
       do while(idmin/=idmax)
          if(idmin==0) then
             write(*,*) '--------------'
             idmin=stree%nnode+1
             idmax=0
             do ielpo = inods1,inods2
                ipoin=msh%lnods(ielpo)
                idmin=min(idmin,stree%part_to_node(ldome(ipoin)))
                idmax=max(idmax,stree%part_to_node(ldome(ipoin)))
             end do
             write(*,*) idmin,idmax
             stop
          end if
          idmin=stree%nodes(idmin)%parent
          idmax=stree%nodes(idmax)%parent
       end do
       ! END DEBUG

       ldomn(ielem)=idmin
    end do

  end subroutine part_to_sep

  !================================================================================================
  subroutine parts_graph_create( dual_mesh, stree, ldome, ren, parts_graph )
    !-----------------------------------------------------------------------
    ! Generate a fem_graph of parts
    !-----------------------------------------------------------------------
    implicit none
    type(fem_mesh) , intent(in)  :: dual_mesh
    type(sep_tree) , intent(in)  :: stree
    type(renum)    , intent(in)  :: ren
    integer(ip)    , intent(in)  :: ldome(dual_mesh%npoin)
    type(fem_graph), intent(out) :: parts_graph

    ! Local variables 
    integer(ip)     , allocatable      :: iwork(:)    ! Integer ip working array
    integer(ip)                        :: pwork(4)

    !if(ren%n/=dual_mesh%nelem) stop DBG: 

    parts_graph%nv = stree%nparts
    call memalloc (parts_graph%nv+1, parts_graph%ia, __FILE__,__LINE__)

    ! Allocate work space
    pwork(1) = 1
    pwork(2) = pwork(1) + dual_mesh%nnode               ! List of parts around a given point
    pwork(3) = pwork(2) + parts_graph%nv                ! 0/1 vector of visited parts
    pwork(4) = pwork(3) + parts_graph%nv*parts_graph%nv ! (i,j) position holds how many vertices 
                                                        ! reside on the interface of the i-th and 
                                                        ! j-th part

    call memalloc (pwork(4), iwork, __FILE__,__LINE__)

    call parts_graph_create_aux (dual_mesh, stree, ldome, ren%iperm, parts_graph, & 
       &                         iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)), &
       &                         iwork(pwork(3):pwork(4)))

    call memfree (iwork,__FILE__,__LINE__)

    return

  end subroutine parts_graph_create

  !------------------------------------------------------------------------------------
  ! This routine creates a graph of parts (parts adjacencies). 
  ! As a by product it also counts the number of points on each interface 
  ! (stored temporarily in ws_itfc_verts) which must be used to sort this graph,
  ! numbering edges with more nodes first (TO DO).
  ! A coloring of the edges is also needed for an efficient communication
  ! in an MPI context (TO DO).
  ! In contrast to the old routine domains_graph, this routine does not longer
  ! store self-edges for the graph of parts.
  !------------------------------------------------------------------------------------
  !  BEGIN OLD INTERFACES
  !  subroutine domains_graph(ndoms,nleaf,npoin,nelem,nelpo,pelpo,lelpo,ltree,ptree, &
  !     &                     ldome,pleaf,iperm,ldomp,ldaux,ndadj,idadj,jdadj)
  !  subroutine domains_graph(ndoms,nleaf,npoin,nelem,nelpo,pelpo,lelpo,ltree,ptree, &
  !    &                    ldome,pleaf,iperm,ldomp,ldaux,idadj,jdadj)
  !  END OLD INTERFACES
  !------------------------------------------------------------------------------------

  subroutine parts_graph_create_aux (dual_mesh , stree, ldome, ipren, parts_graph, &
     &                               ws_parts_list, ws_parts_visited, ws_itfc_verts)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)    :: dual_mesh
    type(sep_tree) , intent(in)    :: stree
    integer(ip)    , intent(in)    :: ldome(dual_mesh%npoin)
    integer(ip)    , intent(in)    :: ipren(dual_mesh%nelem)
    type(fem_graph), intent(inout) :: parts_graph
    integer(ip)                    :: inods1d, inods2d  


    ! Work space:
    ! List of parts around a given point
    integer(ip), intent(inout)           :: ws_parts_list(dual_mesh%nnode)
    ! 0/1 vector of visited parts
    integer(ip), intent(inout)           :: ws_parts_visited(parts_graph%nv)
    ! (i,j) position holds how many vertices reside on the interface of the i-th and j-th part
    integer(ip), intent(inout)           :: ws_itfc_verts(parts_graph%nv,parts_graph%nv)

    ! Local variables
    integer(ip)  :: ipart, jpart, inode, ipoinpm, jpoinpm, posdm, ipoindm, pospartr, pospartc
    integer(ip)  :: nparts, nedges

    ! Initializations
    ws_itfc_verts = 0

    ! Loop on the tree nodes
    do inode=1, stree%nnode
       !if ( stree%nodes(inode)%level /= stree%nlevel  ) then 
       if ( stree%nodes(inode)%rson /= 0 ) then 
         ! Loop over the primal mesh points of the current tree node
          do ipoinpm = stree%node_ptrs(inode), stree%node_ptrs(inode+1)-1
             ws_parts_visited = 0    
             jpoinpm       = ipren(ipoinpm)
             nparts        = 0

             ! Loop over the dual mesh points (primal elements) of the current primal mesh point
              if( dual_mesh%nelty==1 ) then
                inods1d=(jpoinpm-1)*dual_mesh%nnode+1
                inods2d=jpoinpm*dual_mesh%nnode
              else
                inods1d=dual_mesh%pnods(jpoinpm)
                inods2d=dual_mesh%pnods(jpoinpm+1)-1
             end if
             do posdm = inods1d, inods2d
             ! do posdm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                ipoindm = dual_mesh%lnods( posdm )
                jpart   = ldome( ipoindm )

                if ( ws_parts_visited(jpart) == 0 ) then
                   nparts = nparts + 1
                   ws_parts_list (nparts) = jpart
                   ws_parts_visited(jpart) = 1 
                end if
             end do ! posdm

             ! Increment one unit the vertices residing on the
             ! interface of each pair (ws_parts_list(:),ws_parts_list(:)) 
             do pospartc=1, nparts
                do pospartr=1, nparts
                   ws_itfc_verts(ws_parts_list(pospartr), ws_parts_list(pospartc)) = &
                      & ws_itfc_verts(ws_parts_list(pospartr), ws_parts_list(pospartc)) + 1
                end do
             end do ! nparts

          end do ! ipoinpm 
       end if
    end do  ! inode

    ! write(*,*) 'Number of vertices on each interface:' ! DBG: 
    ! do ipart=1,nparts                                  ! DBG:
    !    write(*,'(8i10)') ws_itfc_verts(:,ipart)        ! DBG:
    ! end do

    ! Count edges in parts_graph and set ia pointers
    nedges = 1
    do ipart=1, stree%nparts
       parts_graph%ia (ipart) = nedges 
       do jpart=1, stree%nparts
          if( ( ipart /= jpart ) .and. &   ! Exclude self-edges
             & ( ws_itfc_verts(jpart, ipart) > 0) )  nedges = nedges + 1
          !if(  &                             ! Include self-edges
          !& ( ws_itfc_verts(jpart, ipart) > 0) )  nedges = nedges + 1
       end do
    end do

    parts_graph%ia( parts_graph%nv+1 ) = nedges
    call memalloc(nedges-1, parts_graph%ja, __FILE__,__LINE__)

    ! write (*,*) 'parts_graph%ia='                      ! DBG:
    ! write (*,*) parts_graph%ia (1:(parts_graph%nv+1))  ! DBG:

    ! List edges in parts_graph, i.e., set ja pointers
    nedges = 1
    do ipart=1, stree%nparts
       do jpart=1, stree%nparts
          !if( &                                ! Include self-edges
          !& ( ws_itfc_verts(jpart, ipart) > 0) ) then
          !  parts_graph%ja (nedges)=jpart
          !  nedges = nedges + 1
          !end if 
          if( ( ipart /= jpart )  .and. &   ! Exclude self-edges
             & ( ws_itfc_verts(jpart, ipart) > 0) ) then
             parts_graph%ja (nedges)=jpart
             nedges = nedges + 1
          end if
       end do
    end do

  end subroutine parts_graph_create_aux

  !================================================================================================
  subroutine parts_maps(add_ext,dual_mesh,nren,eren,ldome,nparts,max_nparts,nobjs, &
     &                  lobjs,part_objs,npdom,nedom,nmaps,emaps)
    implicit none
    integer(ip)   , intent(in)    :: add_ext,nparts,max_nparts,nobjs
    type(fem_mesh), intent(in)    :: dual_mesh
    type(renum)   , intent(in)    :: nren,eren
    integer(ip)   , intent(inout) :: ldome(dual_mesh%npoin)
    integer(ip)   , intent(in)    :: lobjs(max_nparts+4,nobjs)
    type(list)    , intent(in)    :: part_objs
    integer(ip)   , intent(out)   :: npdom(nparts+1),nedom(nparts+1)
    type(map)     , intent(out)   :: nmaps(nparts),emaps(nparts)
    ! Locals
    integer(ip), allocatable :: ws_elm_touched(:)
    integer(ip)              :: ipart,p_iobj,iobj,ipoin,kpoin,ielem,p_ielem,jelem
    integer(ip)              :: inods1d, inods2d  


    do ipart=1,nparts
       call map_alloc(npdom(ipart+1),dual_mesh%nelem,nmaps(ipart))
       call map_alloc(nedom(ipart+1),dual_mesh%npoin,emaps(ipart))
    end do

    call memalloc(dual_mesh%npoin,ws_elm_touched,__FILE__,__LINE__)

    ! List nodes assigned to each part (in the order given by the objects)
    ! Reuse npdom to count
    npdom=0
    npdom(1)=1
    do ipart=1,nparts
       ! Loop on part objects
       p_iobj = part_objs%p(ipart)
       iobj = part_objs%l(p_iobj)
       ! Loop on object nodes
       do ipoin=lobjs(2,iobj),lobjs(3,iobj)
          npdom(ipart+1) = npdom(ipart+1) + 1
          nmaps(ipart)%l2g(npdom(ipart+1)) = ipoin
       end do
       nmaps(ipart)%ni = npdom(ipart+1)
       do p_iobj = part_objs%p(ipart)+1,part_objs%p(ipart+1)-1
          iobj = part_objs%l(p_iobj)
          ! Loop on object nodes
          do ipoin=lobjs(2,iobj),lobjs(3,iobj)
             npdom(ipart+1) = npdom(ipart+1) + 1
             nmaps(ipart)%l2g(npdom(ipart+1)) = ipoin
          end do
       end do
       nmaps(ipart)%nb = npdom(ipart+1) - nmaps(ipart)%ni
    end do
    do ipart=1,nparts
       npdom(ipart+1)=npdom(ipart+1)+npdom(ipart)
    end do

    ! List elements assigned to each part
    nedom=0
    nedom(1)=1
    ! Loop on parts
    do ipart=1,nparts
       ws_elm_touched=0
       if(add_ext==one) then       ! List external elements first
          ! Loop on part (boundary) objects
          do p_iobj=part_objs%p(ipart)+1,part_objs%p(ipart+1)-1
             iobj=part_objs%l(p_iobj)
             ! Loop on object nodes
             do kpoin=lobjs(2,iobj),lobjs(3,iobj)
                ipoin=nren%iperm(kpoin)
                ! List untouched external elements around ipoin
                 if(dual_mesh%nelty==1) then
                   inods1d=(ipoin-1)*dual_mesh%nnode+1
                   inods2d=ipoin*dual_mesh%nnode
                else
                   inods1d=dual_mesh%pnods(ipoin)
                   inods2d=dual_mesh%pnods(ipoin+1)-1
                end if
                do p_ielem=inods1d, inods2d
                ! do p_ielem=dual_mesh%pnods(ipoin),dual_mesh%pnods(ipoin+1)-1
                   ielem=dual_mesh%lnods(p_ielem)
                   if(ldome(ielem)/=ipart.and.ldome(ielem)>0) then
                      ldome(ielem)=-ldome(ielem)
                      nedom(ipart+1)=nedom(ipart+1)+1
                      ws_elm_touched(nedom(ipart+1))=ielem
                      emaps(ipart)%l2g(nedom(ipart+1))=eren%lperm(ielem)
                   end if
                end do
             end do
          end do
          emaps(ipart)%ne = nedom(ipart+1)
       else
          ! In case of an element-based partition (i.e., add_ext==zero)
          ! external elements are not part of to the local mesh, so that
          ! we set to zero the number of external elements
          emaps(ipart)%ne = 0
       end if
       ! List elements adjacent to external nodes first
       ! Loop on part (boundary) objects
       do p_iobj=part_objs%p(ipart)+1,part_objs%p(ipart+1)-1
          iobj=part_objs%l(p_iobj)
          ! Loop on object nodes
          do kpoin=lobjs(2,iobj),lobjs(3,iobj)
             ipoin=nren%iperm(kpoin)
             ! List untouched internal elements around ipoin
             if(dual_mesh%nelty==1) then
               inods1d=(ipoin-1)*dual_mesh%nnode+1
               inods2d=ipoin*dual_mesh%nnode
             else
               inods1d=dual_mesh%pnods(ipoin)
               inods2d=dual_mesh%pnods(ipoin+1)-1
             end if
             do p_ielem=inods1d, inods2d
             ! do p_ielem=dual_mesh%pnods(ipoin),dual_mesh%pnods(ipoin+1)-1
                ielem=dual_mesh%lnods(p_ielem)
                if(ldome(ielem)==ipart.and.ldome(ielem)>0) then
                   ldome(ielem)=-ldome(ielem)
                   nedom(ipart+1)=nedom(ipart+1)+1
                   ws_elm_touched(nedom(ipart+1))=ielem
                   emaps(ipart)%l2g(nedom(ipart+1))=eren%lperm(ielem)
                end if
             end do
          end do
       end do
       emaps(ipart)%nb = nedom(ipart+1) - emaps(ipart)%ne
       ! List elements adjacent to external nodes last
       ! Loop on part (internal) objects
       p_iobj=part_objs%p(ipart)
       iobj=part_objs%l(p_iobj)
       ! Loop on object nodes
       do kpoin=lobjs(2,iobj),lobjs(3,iobj)
          ipoin=nren%iperm(kpoin)
          ! List untouched internal elements around ipoin
          if(dual_mesh%nelty==1) then
             inods1d=(ipoin-1)*dual_mesh%nnode+1
             inods2d=ipoin*dual_mesh%nnode
          else
             inods1d=dual_mesh%pnods(ipoin)
             inods2d=dual_mesh%pnods(ipoin+1)-1
          end if
          do p_ielem=inods1d, inods2d
          ! do p_ielem=dual_mesh%pnods(ipoin),dual_mesh%pnods(ipoin+1)-1
             ielem=dual_mesh%lnods(p_ielem)
             if(ldome(ielem)==ipart.and.ldome(ielem)>0) then
                ldome(ielem)=-ldome(ielem)
                nedom(ipart+1)=nedom(ipart+1)+1
                ws_elm_touched(nedom(ipart+1))=ielem
                emaps(ipart)%l2g(nedom(ipart+1))=eren%lperm(ielem)
             end if
          end do
       end do
       emaps(ipart)%ni = nedom(ipart+1) - emaps(ipart)%nb - emaps(ipart)%ne

       ! Reverse numbering and untouch ws
       do ielem=1,nedom(ipart+1)/2
          jelem=emaps(ipart)%l2g(ielem)
          emaps(ipart)%l2g(ielem)=emaps(ipart)%l2g(nedom(ipart+1)+1-ielem)
          emaps(ipart)%l2g(nedom(ipart+1)+1-ielem)=jelem
          ldome(ws_elm_touched(ielem))=-ldome(ws_elm_touched(ielem))
       end do
       do ielem=nedom(ipart+1)/2+1,nedom(ipart+1)
          ldome(ws_elm_touched(ielem))=-ldome(ws_elm_touched(ielem))
       end do

    end do
    do ipart=1,nparts
       nedom(ipart+1)=nedom(ipart+1)+nedom(ipart)
    end do

    call memfree(ws_elm_touched,__FILE__,__LINE__)

  end subroutine parts_maps

  !================================================================================================
  subroutine parts_sizes(add_ext,dual_mesh,ren,ldome,nparts,max_nparts,nobjs, &
     &                   lobjs,part_objs,npdom,nedom)
    implicit none
    integer(ip), intent(in)    :: add_ext,nparts,max_nparts,nobjs
    type(fem_mesh), intent(in) :: dual_mesh
    type(renum), intent(in)    :: ren
    integer(ip), intent(in)    :: lobjs(max_nparts+4,nobjs)
    type(list) , intent(in)    :: part_objs
    integer(ip), intent(inout) :: ldome(dual_mesh%npoin)
    integer(ip), intent(out)   :: npdom(nparts+1),nedom(nparts+1)
    ! Locals
    integer(ip), allocatable   :: ws_elm_touched(:)
    integer(ip)                :: ipart,p_iobj,iobj,ipoin,kpoin,ielem,p_ielem,ielem_ext,nelem_ext
    integer(ip)                :: inods1d, inods2d  

    call memalloc(dual_mesh%npoin,ws_elm_touched,__FILE__,__LINE__)

    ! Count nodes assigned to each part
    npdom=0
    npdom(1)=1
    do ipart=1,nparts
       do p_iobj = part_objs%p(ipart),part_objs%p(ipart+1)-1
          iobj = part_objs%l(p_iobj)
          npdom(ipart+1) = npdom(ipart+1) + lobjs(3,iobj)-lobjs(2,iobj)+1
       end do
    end do
!!$    do ipart=1,nparts
!!$       npdom(ipart+1)=npdom(ipart+1)+npdom(ipart)
!!$    end do

    ! Count elements assigned to each part
    nedom(1:nparts+1)=0
    nedom(1)=1
    do ielem=1,dual_mesh%npoin
       nedom(ldome(ielem)+1)=nedom(ldome(ielem)+1)+1
    end do
    if(add_ext==one) then
       ! Count external elements
       ! Loop on parts
       do ipart=1,nparts
          nelem_ext=0
          ws_elm_touched=0
          ! Loop on part objects
          do p_iobj=part_objs%p(ipart)+1,part_objs%p(ipart+1)-1
             iobj=part_objs%l(p_iobj)
             ! Loop on object nodes
             do kpoin=lobjs(2,iobj),lobjs(3,iobj)
                ipoin=ren%iperm(kpoin)

                if(dual_mesh%nelty==1) then
                  inods1d=(ipoin-1)*dual_mesh%nnode+1
                  inods2d=ipoin*dual_mesh%nnode
                else
                  inods1d=dual_mesh%pnods(ipoin)
                  inods2d=dual_mesh%pnods(ipoin+1)-1
                end if
                do p_ielem=inods1d, inods2d
                !do p_ielem=dual_mesh%pnods(ipoin),dual_mesh%pnods(ipoin+1)-1
                   ielem=dual_mesh%lnods(p_ielem)
                   if(ldome(ielem)/=ipart.and.ldome(ielem)>0) then
                      nelem_ext=nelem_ext+1
                      ws_elm_touched(nelem_ext)=ielem
                      ldome(ielem)=-ldome(ielem)
                      !nedom(ipart+1)=nedom(ipart+1)+1
                   end if
                end do
             end do
          end do
          nedom(ipart+1)=nedom(ipart+1)+nelem_ext
          do ielem_ext=1,nelem_ext
             ldome(ws_elm_touched(ielem_ext))=-ldome(ws_elm_touched(ielem_ext))
          end do
       end do
    end if
!!$    do ipart=1,nparts
!!$       nedom(ipart+1)=nedom(ipart+1)+nedom(ipart)
!!$    end do

    call memfree(ws_elm_touched,__FILE__,__LINE__)

  end subroutine parts_sizes

  !================================================================================================
  subroutine parts_sizes_boundary(add_ext,mesh,ldome,ldomn,nparts,nbdom)
    implicit none
    integer(ip), intent(in)    :: add_ext,nparts
    type(fem_mesh), intent(in) :: mesh
    integer(ip), intent(inout) :: ldome(mesh%nelem),ldomn(mesh%npoin)
    integer(ip), intent(out)   :: nbdom(nparts+1)
    ! Locals
    integer(ip)                :: ielem,inode,jnode,touched(mesh%nnodb),gnode,count

    ! Count boundary elements assigned to each part
    if(add_ext == zero) then              ! Element based partition
       nbdom=0
       nbdom(1)=1
       do ielem=1,mesh%nboun
          nbdom(ldome(mesh%lboel(mesh%nnodb+1,ielem))+1)=&
               nbdom(ldome(mesh%lboel(mesh%nnodb+1,ielem))+1) + 1
       end do
    elseif(add_ext == one) then           ! Vertex based partition
       nbdom=0
       nbdom(1)=1
       do ielem=1,mesh%nboun
          touched=0
          do inode=1,mesh%nnodb
             gnode = mesh%lboun((ielem-1)*mesh%nnodb+inode)
             ! Avoid repeating the addition of a boundary element to a part because
             ! the element has several nodes belonging to the same part.
             count = 0
             do jnode=1,inode-1
                if (ldomn(gnode)==touched(jnode)) count=count+1
             end do
             if (count==0) then
                nbdom(ldomn(gnode)+1) = nbdom(ldomn(gnode)+1) + 1
                touched(inode) = ldomn(gnode)
             end if
          end do
       end do
    end if

  end subroutine parts_sizes_boundary

  !================================================================================================
  subroutine parts_maps_boundary(add_ext,mesh,ldome,ldomn,nparts,nbdom,bmaps)
    implicit none
    integer(ip), intent(in)    :: add_ext,nparts
    type(fem_mesh), intent(in) :: mesh
    integer(ip), intent(in)    :: ldome(mesh%nelem),ldomn(mesh%npoin)
    integer(ip), intent(in)    :: nbdom(nparts+1)
    type(map)  , intent(inout) :: bmaps(nparts)
    ! Locals
    integer(ip)                :: iboun,ipart,icnt(nparts)
    integer(ip)                :: inode,jnode,touched(mesh%nnodb),gnode,count

    if(add_ext == zero) then        ! Element based
       icnt = 1
       ! Loop on boundary elements
       do iboun=1,mesh%nboun
          ipart = ldome(mesh%lboel(mesh%nnodb+1,iboun))
          bmaps(ipart)%l2g(icnt(ipart)) = iboun
          bmaps(ipart)%ni = bmaps(ipart)%ni + 1
          icnt(ipart) = icnt(ipart) + 1
       end do

    elseif(add_ext == one) then     ! Vertex based
       ! Interior boundary elements
       icnt = 1
       ! Loop on boundary elements
       do iboun=1,mesh%nboun
          ! We check that all nodes belong to the same part
          count = 0
          gnode = mesh%lboun((iboun-1)*mesh%nnodb+1)
          ipart = ldomn(gnode)
          do inode=2,mesh%nnodb
             if(ldomn(mesh%lboun((iboun-1)*mesh%nnodb+inode)).ne.ipart) count = count + 1
          end do
          if(count==0) then
             bmaps(ipart)%l2g(icnt(ipart)) = iboun
             bmaps(ipart)%ni = bmaps(ipart)%ni + 1
             icnt(ipart) = icnt(ipart) + 1
          end if
       end do

       ! Next, we add the boundary elements on the partition boundaries
       do iboun=1,mesh%nboun
          ! We check that there are nodes belonging to different parts
          count = 0
          gnode = mesh%lboun((iboun-1)*mesh%nnodb+1)
          ipart = ldomn(gnode)
          do inode=2,mesh%nnodb
             if(ldomn(mesh%lboun((iboun-1)*mesh%nnodb+inode)).ne.ipart) count = count + 1
          end do
          if(count.gt.0) then
             touched=0
             do inode=1,mesh%nnodb
                gnode = mesh%lboun((iboun-1)*mesh%nnodb+inode)
                ! Avoid repeating the addition of a boundary element to a part because
                ! the element has several nodes belonging to the same part.
                count = 0
                do jnode=1,inode-1
                   if (ldomn(gnode)==touched(jnode)) count=count+1
                end do
                if (count==0) then
                   bmaps(ldomn(gnode))%l2g(icnt(ldomn(gnode))) = iboun
                   bmaps(ldomn(gnode))%nb = bmaps(ldomn(gnode))%nb + 1
                   icnt(ldomn(gnode)) = icnt(ldomn(gnode)) + 1
                   touched(inode) = ldomn(gnode)
                end if
             end do
          end if
       end do

    end if

  end subroutine parts_maps_boundary


  !================================================================================================
  subroutine neighbors_graph_extract(emap,nelem,ldome,eren,graph,edges,pextn,lextn,lextp,lexte,lextm,material)
    implicit none
    type(map)      , intent(in)   :: emap
    integer(ip)    , intent(in)   :: nelem
    integer(ip)    , intent(in)   :: ldome(nelem)
    type(renum)    , intent(in)   :: eren
    type(fem_graph), intent(in)   :: graph
    integer(ip)    , intent(in)   :: edges(graph%ia(nelem+1)-1)
    type(fem_materials), optional, intent(in)  :: material ! Material associated to elements

    integer(ip)    , allocatable,  intent(out)  :: pextn(:)
    integer(ip)    , allocatable,  intent(out)  :: lextn(:)
    integer(ip)    , allocatable,  intent(out)  :: lextp(:)
    integer(ip)    , allocatable,  intent(out)  :: lexte(:)
    integer(ip)    , allocatable,  intent(out)  :: lextm(:)
    ! Locals
    integer(ip)                   :: ipart, ieloc, ielem, index, neighbor, next

    ! The first element is always local and that gives current part.
    ipart = ldome(eren%iperm(emap%l2g(1)))
    call memalloc (emap%nb+1, pextn, __FILE__,__LINE__)
    !pextn = 0
    pextn(1) = 1
    do ieloc = 1,emap%nb
       pextn(ieloc+1) = 0
       ielem =  eren%iperm(emap%l2g(emap%ni+ieloc))
       do index = graph%ia(ielem), graph%ia(ielem+1)-1
          neighbor = graph%ja(index)
          if(ldome(neighbor)/=ipart) pextn(ieloc+1) = pextn(ieloc+1) + 1
       end do
    end do
    do ieloc = 1,emap%nb
       pextn(ieloc+1) = pextn(ieloc+1) + pextn(ieloc)
    end do
    !write(*,*) emap%nb,pextn(emap%nb+1)

    call memalloc (pextn(emap%nb+1)-1, lextn, __FILE__,__LINE__)
    call memalloc (pextn(emap%nb+1)-1, lextp, __FILE__,__LINE__)
    call memalloc (pextn(emap%nb+1)-1, lexte, __FILE__,__LINE__)

    if(present(material) ) then

       assert(material%nelem==nelem)
       call memalloc (pextn(emap%nb+1)-1, lextm, __FILE__,__LINE__)
       next=1
       do ieloc = 1,emap%nb
          ielem =  eren%iperm(emap%l2g(emap%ni+ieloc))
          do index = graph%ia(ielem), graph%ia(ielem+1)-1
             neighbor = graph%ja(index)
             if(ldome(neighbor)/=ipart) then
                lextn(next) = eren%lperm(neighbor)
                lextp(next) = ldome(neighbor)
                lexte(next) = edges(index)
                lextm(next) = material%list(neighbor)
                next = next + 1
             end if
          end do
       end do

    else

       next=1
       do ieloc = 1,emap%nb
          ielem =  eren%iperm(emap%l2g(emap%ni+ieloc))
          do index = graph%ia(ielem), graph%ia(ielem+1)-1
             neighbor = graph%ja(index)
             if(ldome(neighbor)/=ipart) then
                lextn(next) = eren%lperm(neighbor)
                lextp(next) = ldome(neighbor)
                lexte(next) = edges(index)
                next = next + 1
             end if
          end do
       end do
       call memalloc (pextn(emap%nb+1)-1, lextm, __FILE__,__LINE__,0)

    end if

  end subroutine neighbors_graph_extract
  !================================================================================================
  subroutine objects_create (dual_mesh, ldomn, ldome, ren, max_nparts, nobjs, lobjs, stree)
    !-----------------------------------------------------------------------
    ! Generates objects
    ! General comment: "+1", "+4", "+3" are hard-coded magic numbers that
    ! should be avoided. It would increase the cleanliness/quality of the code to have
    ! parameters/macros defining these values.
    !-----------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(in)    :: dual_mesh
    integer(ip)   , intent(in)    :: ldomn(dual_mesh%nelem)
    integer(ip)   , intent(in)    :: ldome(dual_mesh%npoin)
    type(renum)   , intent(inout) :: ren
    integer(ip)   , intent(out)   :: max_nparts, nobjs 
    integer(ip)   , allocatable   :: lobjs(:,:)
    type(sep_tree), intent(in)    :: stree

    !type(mesh_partition), intent(inout) :: mesh_part

    ! Local variables 
    integer(ip)     , allocatable :: iwork(:)    ! Adjustable size working space
    integer(ip)                   :: pwork(8)
    integer(ip)                   :: npsep 

    !if(ren%n/=dual_mesh%nelem) stop DBG: 

    ! npsep = total number of nodes of the mesh that belong to the separators
    !       = sum(number of vertices of inner stree nodes, i.e., non-leaf nodes)
    ! npsep (estimate) >> nobjs (actually computed in objects_create_aux)
    call sep_tree_sep_size(stree, npsep) 

    pwork(1) = 1
    pwork(2) = pwork(1) + stree%nparts              ! 0/1 vector of visited parts
    pwork(3) = pwork(2) + (dual_mesh%nnode+4)*npsep    ! Temporary list of objects
    pwork(4) = pwork(3) + (dual_mesh%nnode+1)*npsep    ! List of parts around points on the 
                                                       ! current separator
    pwork(5) = pwork(4) + npsep                        ! Permutation array for the list of objects
    pwork(6) = pwork(5) + dual_mesh%nnode+1            ! Sort workspace
    pwork(7) = pwork(6) + dual_mesh%nnode+1            ! Sort workspace

    ! Allocate work space
    call memalloc (pwork(7), iwork, __FILE__,__LINE__)

    call objects_create_aux (dual_mesh, ldomn, ldome, ren%lperm, ren%iperm, &
       &                     max_nparts, nobjs, lobjs, stree, npsep, &
       &                     iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)),  &
       &                     iwork(pwork(3):pwork(4)), iwork(pwork(4):pwork(5)),  &
       &                     iwork(pwork(5):pwork(6)), iwork(pwork(6):pwork(7)) )

    ! Deallocate work space
    call memfree (iwork,__FILE__,__LINE__)

    return

  end subroutine objects_create

  !================================================================================================
  subroutine objects_create_aux ( dual_mesh, ldomn, ldome, lpren, ipren,  &
     &                            max_nparts, nobjs, lobjs, stree, nobjs_estimate,    &
     &                            ws_parts_visited, ws_lobjs_temp, ws_parts_list_sep, &
     &                            ws_lobjs_perm, ws_sort_l1, ws_sort_l2)

    !-----------------------------------------------------------------------
    ! Generates objects, auxiliary routine
    ! General comment: expresions like "1", "2", "+1", "+4", "+3" are hard-coded magic numbers that
    ! should be avoided. It would increase the cleanliness/quality of the code to have
    ! parameters/macros defining these values (TODO).
    !-----------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(in)    :: dual_mesh
    integer(ip)   , intent(in)    :: ldomn(dual_mesh%nelem)
    integer(ip)   , intent(in)    :: ldome(dual_mesh%npoin)
    integer(ip)   , intent(inout) :: lpren(dual_mesh%nelem)
    integer(ip)   , intent(inout) :: ipren(dual_mesh%nelem)
    !type(mesh_partition), intent(inout) :: mesh_part
    integer(ip)   , intent(out)   :: max_nparts, nobjs
    integer(ip)   , allocatable   :: lobjs(:,:)
    type(sep_tree), intent(in)    :: stree

    ! Work space
    integer(ip), intent(in)      :: nobjs_estimate
    ! 0/1 vector of visited parts
    integer(ip), intent(inout)   :: ws_parts_visited  (stree%nparts)
    ! Temporary list of objects
    integer(ip), intent(inout)   :: ws_lobjs_temp     (dual_mesh%nnode+4,nobjs_estimate)
    ! List of parts around points of the current separator
    integer(ip), intent(inout)   :: ws_parts_list_sep (dual_mesh%nnode+1,nobjs_estimate)
    ! Permutation of the list of objects   
    integer(ip), intent(inout)   :: ws_lobjs_perm     (nobjs_estimate)
    ! Local work space for sorting algorithms and comparison function
    integer(ip), intent(inout)   :: ws_sort_l1(dual_mesh%nnode+1)
    integer(ip), intent(inout)   :: ws_sort_l2(dual_mesh%nnode+1)

    ! Local variables 
    integer(ip)                  :: inode, nparts, max_nparts_cur_sep !, max_nparts
    integer(ip)                  :: ipoinpm,lipoinpm,ipoindm,jpoinpm,posdm,ipart,jpart
    integer(ip)                  :: i,j
    integer(ip)                  :: inods1d, inods2d  


    ! A. Martin: This is dirty and deprecated. We should declare an interface 
    ! to icomp somewhere !!! I guess that the same has to
    ! be applied for intsort and sortix sort subroutines (TODO)
    !integer icomp
    !external icomp    

    max_nparts = 0
    nobjs      = 0

    ! Initialize temporary list of objects 
    ! (not actually required)
    ws_lobjs_temp = 0

    ! Loop over separator tree nodes ...
    do inode=1, stree%nnode
       ! If inode is a separator ...
       !if ( stree%nodes(inode)%level /= stree%nlevel  ) then
       if ( stree%nodes(inode)%rson /= 0 ) then
          max_nparts_cur_sep = 0
          ws_parts_list_sep  = 0
          do ipoinpm = stree%node_ptrs(inode), stree%node_ptrs(inode+1)-1
             ws_parts_visited = 0    
             jpoinpm          = ipren(ipoinpm)
             nparts           = 0
             lipoinpm         = ipoinpm+1-stree%node_ptrs(inode)

             if(dual_mesh%nelty==1) then
               inods1d=(jpoinpm-1)*dual_mesh%nnode+1
               inods2d=jpoinpm*dual_mesh%nnode
             else
               inods1d=dual_mesh%pnods(jpoinpm)
               inods2d=dual_mesh%pnods(jpoinpm+1)-1
              end if

             ! Loop over the dual mesh points (primal elements) of the current primal mesh point
             do posdm=inods1d, inods2d
             ! do posdm = dual_mesh%pnods(jpoinpm), dual_mesh%pnods(jpoinpm+1)-1
                ipoindm = dual_mesh%lnods( posdm )
                jpart   = ldome( ipoindm )            
                if ( ws_parts_visited(jpart) == 0 ) then
                   nparts = nparts + 1
                   ws_parts_list_sep (nparts+1,lipoinpm) = jpart
                   ws_parts_visited(jpart) = 1 
                end if
             end do ! posdm
             ws_parts_list_sep(1,lipoinpm) = nparts

             ! Sort list of parts in increasing order by part identifiers
             ! This is required by the call to icomp subroutine below 
             call simple_sort ( nparts, ws_parts_list_sep( 2:(nparts+1), lipoinpm) )

             ! write(*,*) 'List of parts around node ', ipoinpm, ' is ', &
             !   &         ws_parts_list_sep (1:(nparts+1),lipoinpm) ! DBG:

             max_nparts_cur_sep = max(max_nparts_cur_sep, nparts)

          end do ! ipoinpm      

          ! write(*,*) 'Max. number of parts around separator ', inode, ' is ',  &
          ! &           max_nparts_cur_sep                  ! DBG:  

          max_nparts = max(max_nparts_cur_sep, max_nparts)

          ! Re-number nodes of the current separator in increasing order by the 
          ! number of parts they belong and, among vertices sharing the same parts,
          ! in increasing order by the list of parts shared by each vertex

          call intsort(max_nparts_cur_sep+1, & ! Rows of ws_parts_list_sep
             &         dual_mesh%nnode+1,    & ! Leading dimension of ws_parts_list_sep
             &         stree%node_ptrs(inode+1)- stree%node_ptrs(inode), &
                                               ! Cols of ws_parts_list_sep
             &         ws_parts_list_sep,    &
             &         ipren(stree%node_ptrs(inode):  &
             &                        (stree%node_ptrs(inode+1)-1)), &
             &         ws_sort_l1, ws_sort_l2)

          ! Identify interface objects 
          ! (faces, edges and vertices)
          nobjs=nobjs+1

          ! Prepare first object
          ws_lobjs_temp(1,nobjs)=inode                               ! separator node id
          ws_lobjs_temp(2,nobjs)=stree%node_ptrs(inode)           ! begin obj
          ws_lobjs_temp(4,nobjs)=ws_parts_list_sep(1,1)
          ws_lobjs_temp(5:4+ws_parts_list_sep(1,1),nobjs) = &
             &  ws_parts_list_sep(2:1+ws_parts_list_sep(1,1),1)

          ! Loop over vertices of the current separator (i.e., inode)
          do i=1,stree%node_ptrs(inode+1)-stree%node_ptrs(inode)-1
             if(icomp(max_nparts_cur_sep+1,ws_parts_list_sep(:,i),ws_parts_list_sep(:,i+1))/=0) then
                ! Complete current object
                ws_lobjs_temp(3,nobjs)=stree%node_ptrs(inode)+i-1 ! end obj
                ws_lobjs_temp(4,nobjs)=ws_parts_list_sep(1,i)
                ws_lobjs_temp(5:4+ws_parts_list_sep(1,i),nobjs) = &
                   & ws_parts_list_sep(2:1+ws_parts_list_sep(1,i),i)

                ! Prepare next object
                nobjs=nobjs+1
                ws_lobjs_temp(1,nobjs)=inode                          ! separator node id
                ws_lobjs_temp(2,nobjs)=stree%node_ptrs(inode)+i   ! begin obj
             end if
          end do
          ! Complete last object
          ws_lobjs_temp(3,nobjs)=stree%node_ptrs(inode)+i-1       ! end obj
          ws_lobjs_temp(4,nobjs)=ws_parts_list_sep(1,i)
          ws_lobjs_temp(5:4+ws_parts_list_sep(1,i),nobjs) = &
             & ws_parts_list_sep(2:1+ws_parts_list_sep(1,i),i)

          ! Finally re-compute permutation from its inverse
          do ipoinpm=stree%node_ptrs(inode),stree%node_ptrs(inode+1)-1
             lpren(ipren(ipoinpm))=ipoinpm
          end do
       end if
    end do ! inode

    ! write(*,*) 'Total number of objects ', nobjs               ! DBG:  


    ! Add one object per part to ws_lobjs_temp
    ! in order to accomodate interior nodes of the mesh
    do ipart=1,stree%nparts
       inode=stree%part_to_node(ipart)
       ws_lobjs_temp(1,nobjs+1) = inode
       ws_lobjs_temp(2,nobjs+1) = stree%node_ptrs(inode)
       ws_lobjs_temp(3,nobjs+1) = stree%node_ptrs(inode+1)-1
       ws_lobjs_temp(4,nobjs+1) = 1
       ws_lobjs_temp(5,nobjs+1) = stree%nodes(inode)%part
       ws_lobjs_temp(6:max_nparts+4,nobjs+1)=0
       nobjs=nobjs+1
    end do

    ! write(*,*) 'Total number of objects ', nobjs               ! DBG:  

    do i=1, nobjs_estimate
       ws_lobjs_perm(i) = i
    end do

    ! Sort objects in ws_lobjs_temp accordingly to their 
    ! starting position in ipren (i.e.,  ws_lobjs_temp(2,:)).
    ! Note that here is the point of the code where we associate each object with 
    ! an identifier. After sorting, interior objects will be positioned first 
    ! in the list of objects of each part, so that these objects can be easily located 
    ! without searching
    !call sortix(2, dual_mesh%nnode+4, nobjs, ws_lobjs_temp, ws_lobjs_perm)
    call sort_array_cols_by_row_element(2, dual_mesh%nnode+4, nobjs, ws_lobjs_temp, ws_lobjs_perm)

    ! Allocate space for the list of objects
    call memalloc (max_nparts+4, nobjs, lobjs,__FILE__,__LINE__)

    ! Copy permuted ws_lobjs_temp to mesh_part ...
    do i=1,nobjs
       lobjs(:,i)=ws_lobjs_temp(1:(max_nparts+4), ws_lobjs_perm(i))
    end do

    ! write(*,'(a)') 'List of objects:'               ! DBG:
    ! do i=1,nobjs                          ! DBG:
    !    write(*,'(10i10)') i, lobjs(:,i)   ! DBG:
    ! end do                                          ! DBG:

  end subroutine objects_create_aux

  !================================================================================================
  subroutine part_objects_create (nparts,nobjs,lobjs,part_objs) !mesh_part)
    !-----------------------------------------------------------------------
    ! Generates list of objects of each part
    ! General comment: "4", "+4", etc. are hard-coded magic numbers that
    ! should be avoided. It would increase the cleanliness/quality of the code to have
    ! parameters/macros defining these values.
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nparts, nobjs
    integer(ip), intent(in)  :: lobjs(:,:)
    type(list) , intent(out) :: part_objs

    ! Parameters 
    !type(mesh_partition), intent(inout) :: mesh_part

    ! Local variables
    integer(ip)                                :: i, j, ipart

    part_objs%n = nparts
    call memalloc (part_objs%n+1, part_objs%p, __FILE__,__LINE__)

    ! Count objects of each part
    part_objs%p = 0
    do i=1,nobjs
       do j=1,lobjs(4,i)
          part_objs%p(lobjs(4+j,i)+1) = &
             & part_objs%p(lobjs(4+j,i)+1)+1
       end do
    end do

    ! Generate pointers to the list of objects of each part
    part_objs%p(1) = 1
    do i=1, nparts
       part_objs%p(i+1)=part_objs%p(i+1)+part_objs%p(i)
    end do

    ! Generate list of objects of each part
    call memalloc (part_objs%p(nparts+1)-1, part_objs%l,     __FILE__,__LINE__)
    do i=1,nobjs
       do j=1,lobjs(4,i)
          part_objs%l(part_objs%p(lobjs(4+j,i)))=i
          part_objs%p(lobjs(4+j,i)) = &
             & part_objs%p(lobjs(4+j,i))+1
       end do
    end do

    ! Recover part_objs%p
    do ipart=nparts+1, 2, -1
       part_objs%p(ipart) = part_objs%p(ipart-1)
    end do
    part_objs%p(ipart) = 1


    !write(*,'(a)') 'List of part objects:'                                   ! DBG:
    !do i=1,nparts                                                  ! DBG:
    !    write(*,'(10i10)') i, (part_objs%l(j), &                   ! DBG:
    !    & j=part_objs%p(i),part_objs%p(i+1)-1)           ! DBG:
    !end do

    return
  end subroutine part_objects_create

  subroutine int_objects_create (nobjs,lobjs,parts_graph,int_objs) ! (mesh_part)
    !-----------------------------------------------------------------------
    ! Generates list of objects on each edge of the graph of parts
    ! General comment: "4", "+4", etc. are hard-coded magic numbers that
    ! should be avoided. It would increase the cleanliness/quality of the code to have
    ! parameters/macros defining these values.
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)    , intent(in)  :: nobjs
    integer(ip)    , intent(in)  :: lobjs(:,:)
    type(fem_graph), intent(in)  :: parts_graph
    type(list)     , intent(out) :: int_objs

    ! Parameters 
    !type(mesh_partition), intent(inout) :: mesh_part

    ! Local variables
    integer(ip)                                :: iobj, i, j, ipart, jpart, iedge

    int_objs%n = parts_graph%ia(parts_graph%nv+1)-1    

    ! write (*,*) int_objs%n  ! DBG:

    call memalloc (int_objs%n+1, int_objs%p,         __FILE__,__LINE__)

    ! write (*,*) int_objs%n  ! DBG:


    ! Count objects on each edge of the graph of parts
    int_objs%p = 0

    ! write (*,*)  int_objs%p ! DBG:

    do iobj=1,nobjs
       do i=1,lobjs(4,iobj)
          do j=1,lobjs(4,iobj)
             ipart = lobjs(4+i,iobj)
             jpart = lobjs(4+j,iobj)
             if (ipart /= jpart) then  ! Exclude self-edges 
                ! Locate edge identifier of ipart => jpart on the list of edges of ipart 
                iedge = parts_graph%ia(ipart)
                do while ((parts_graph%ja(iedge) /= jpart).and. &
                   &      (iedge<parts_graph%ia(ipart+1)))
                   iedge = iedge + 1
                end do

                ! Increment by one the number of objects on the edge ipart => jpart
                int_objs%p(iedge+1) = int_objs%p(iedge+1) + 1 
             end if
          end do
       end do
    end do

    ! Generate pointers to the list of objects on each edge
    int_objs%p(1) = 1
    do iedge=1, int_objs%n
       int_objs%p(iedge+1)=int_objs%p(iedge+1)+int_objs%p(iedge)
    end do

    ! Generate list of objects on each edge
    call memalloc (int_objs%p(int_objs%n+1)-1, int_objs%l,            __FILE__,__LINE__)

    do iobj=1,nobjs
       do i=1,lobjs(4,iobj)
          do j=1,lobjs(4,iobj)
             ipart = lobjs(4+i,iobj)
             jpart = lobjs(4+j,iobj)
             if (ipart /= jpart) then ! Exclude self-edges
                ! Locate edge identifier of ipart => jpart on the list of edges of ipart
                iedge = parts_graph%ia(ipart)
                do while ((parts_graph%ja(iedge) /= jpart).and. &
                   &      (iedge<parts_graph%ia(ipart+1)))
                   iedge = iedge + 1 
                end do

                ! Insert current object on the edge ipart => jpart
                int_objs%l(int_objs%p(iedge)) = iobj
                int_objs%p(iedge) = int_objs%p(iedge) + 1 
             end if
          end do
       end do
    end do

    ! Recover int_objs%p
    do iedge=int_objs%n+1, 2, -1
       int_objs%p(iedge) = int_objs%p(iedge-1)
    end do
    int_objs%p(iedge) = 1

    return
  end subroutine int_objects_create

!!$  !============================================================================================
!!$  subroutine sort_subdomains(primal_mesh,dual_mesh,gph,ldomn,ldome,nren,eren,mesh_part)
!!$    !------------------------------------------------------------------------------------
!!$    ! This routine performs CM renumbering on each subdomain starting from separators.
!!$    ! Dual objects are renumbered accordingly, i.e. in an element based partition starting
!!$    ! from nodes on the separator adjacent elements
!!$    ! It renumbers both nodes and elements and uses the nodal graph.
!!$    !------------------------------------------------------------------------------------
!!$    implicit none
!!$    type(mesh)          , intent(in)    :: primal_mesh,dual_mesh
!!$    type(graph)         , intent(inout) :: gph
!!$    integer(ip)         , intent(inout) :: ldomn(primal_mesh%npoin)
!!$    integer(ip)         , intent(inout) :: ldome(primal_mesh%nelem)
!!$    type(renum)         , intent(inout) :: nren
!!$    type(renum)         , intent(inout) :: eren
!!$    type(mesh_partition), intent(inout) :: mesh_part
!!$    integer(ip)         , allocatable   :: iwork(:)    ! Integer ip working array
!!$
!!$    allocate(iwork(primal_mesh%npoin))
!!$    if(mesh_part%use_graph==primal) then
!!$       call sort_subdomains_primal( &
!!$          &      0,mesh_part%nparts,primal_mesh%nelty,primal_mesh%npoin,primal_mesh%nelem, &
!!$          &      primal_mesh%nnode,dual_mesh%pnods,dual_mesh%lnods,mesh_part%stree%part_to_node, &
!!$          &      gph%ia,gph%ja,primal_mesh%pnods,primal_mesh%lnods,ldomn,ldome,mesh_part%max_nparts, &
!!$          &      mesh_part%nobjs,mesh_part%lobjs,mesh_part%part_objs%p,mesh_part%part_objs%l, &
!!$          &      nren%iperm,nren%lperm,eren%iperm,iwork)
!!$    else if(mesh_part%use_graph==dual) then
!!$       call sort_subdomains_dual( &
!!$          &      mesh_part%nparts,primal_mesh%nelty,primal_mesh%npoin,primal_mesh%nelem, &
!!$          &      primal_mesh%nnode,dual_mesh%pnods,dual_mesh%lnods,mesh_part%stree%part_to_node, &
!!$          &      gph%ia,gph%ja,primal_mesh%pnods,primal_mesh%lnods,ldomn,ldome,mesh_part%max_nparts, &
!!$          &      mesh_part%nobjs,mesh_part%lobjs,mesh_part%part_objs%p,mesh_part%part_objs%l, &
!!$          &      nren%iperm,nren%lperm,eren%iperm,iwork)
!!$    end if
!!$    deallocate(iwork)
!!$    return
!!$
!!$  end subroutine sort_subdomains
!!$
!!$  !------------------------------------------------------------------------------------
!!$  ! This routine performs CM renumbering on each subdomain starting from separators.
!!$  ! It renumbers both nodes and elements and uses the element's graph.
!!$  !------------------------------------------------------------------------------------
!!$  subroutine sort_subdomains_dual(ndoms,nelty,npoin,nelem,nnode,pelpo,lelpo, &
!!$     &                                ltree,ia,ja,pnods,lnods,ldomn,ldome,mdomp, &
!!$     &                                nobjs,lobjs,pdobj,ldobj,iperm,lperm,ieren,level)
!!$    implicit none
!!$    integer(ip), intent(in)    :: ndoms,nelty,npoin,nelem,nnode
!!$    integer(ip), intent(in)    :: pelpo(npoin+1),lelpo(pelpo(npoin+1))
!!$    integer(ip), intent(in)    :: ltree(ndoms)
!!$    integer(ip), intent(in)    :: pnods(*),lnods(*)
!!$
!!$    integer(ip), intent(inout) :: ia(nelem+1)
!!$    integer(ip), intent(inout) :: ja(*)
!!$
!!$    integer(ip), intent(in)    :: ldomn(npoin)
!!$    integer(ip), intent(inout) :: ldome(nelem)
!!$    integer(ip), intent(in)    :: mdomp,nobjs
!!$    integer(ip), intent(in)    :: lobjs(mdomp+4,nobjs)
!!$    integer(ip), intent(in)    :: pdobj(ndoms+1)
!!$    integer(ip), intent(in)    :: ldobj(pdobj(ndoms+1))
!!$    integer(ip), intent(inout) :: iperm(npoin)
!!$    integer(ip), intent(inout) :: lperm(npoin)
!!$    integer(ip), intent(inout) :: ieren(nelem)
!!$    integer(ip), intent(inout) :: level(*)
!!$
!!$    integer(ip)                :: idoms,ipoin,jpoin,ielem,ienew,jenew,ielpo
!!$    integer(ip)                :: inode,knode,inods,jnods,i,j,k,iz,ileve,nleve
!!$
!!$    ieren(1:nelem)=0
!!$
!!$    ienew=1
!!$    do idoms=1,ndoms
!!$       ! Loop over the elements around nodes on the separator to define first level
!!$       jenew=ienew
!!$       do i=pdobj(idoms)+1,pdobj(idoms+1)-1
!!$          do j=lobjs(2,ldobj(i)),lobjs(3,ldobj(i))
!!$             ipoin=iperm(j)   ! This is the original node number
!!$             lperm(ipoin)=0   ! Fix its number
!!$             if(i==pdobj(idoms)+1) then ! Only first object as starting level
!!$                do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
!!$                   ielem=lelpo(ielpo)
!!$                   if(ldome(ielem)==idoms.and.ia(ielem)>0) then
!!$                      ieren(jenew)= ielem
!!$                      ia(ielem)=-ia(ielem) ! This is required by rcm_renum (as a touch)
!!$                      jenew=jenew+1
!!$                   end if
!!$                end do
!!$             end if
!!$          end do
!!$       end do
!!$       level(1)=ienew
!!$       level(2)=jenew
!!$       ! do rcm of the elements graph to obtain the element renumbering
!!$       ! using a mask ldome(ielem)==idoms
!!$       call rcm_renum_masked(nelem,ia(nelem+1),ia,ja,ldome,idoms,ieren,nleve,level)
!!$       ! Loop on levels and renumber interior nodes backwards
!!$       i=pdobj(idoms)
!!$       j=lobjs(3,ldobj(i))
!!$       do ileve=1,nleve
!!$          do ienew=level(ileve),level(ileve+1)-1
!!$             ielem=ieren(ienew)
!!$             ia(ielem)=-ia(ielem)
!!$             ! Loop over its nodes and number them
!!$             knode=nnode
!!$             inods=(ielem-1)*knode
!!$             if(nelty>1) then
!!$                knode=pnods(ielem+1)-pnods(ielem)
!!$                inods=pnods(ielem)
!!$             end if
!!$             do inode=1,knode
!!$                ipoin=lnods(inods+inode)
!!$                if(lperm(ipoin)/=0) then
!!$                   lperm(ipoin)=0
!!$                   iperm(j)=ipoin
!!$                   j=j-1
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$       ienew=level(nleve+1)
!!$       ! Finally recompute permutation from its inverse
!!$       do i=pdobj(idoms),pdobj(idoms+1)-1
!!$          do j=lobjs(2,ldobj(i)),lobjs(3,ldobj(i))
!!$             lperm(iperm(j))=j
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$  end subroutine sort_subdomains_dual
!!$
!!$  !------------------------------------------------------------------------------------
!!$  ! This routine performs CM renumbering on each subdomain starting from separators.
!!$  ! It renumbers both nodes and elements and uses the nodal graph.
!!$  !
!!$  ! On each part:
!!$  !
!!$  ! 1) Loop over separators (vertices in objects on the interfaces) and touch elements
!!$  !    defining a first level
!!$  ! 2) 
!!$  !
!!$  !------------------------------------------------------------------------------------
!!$  subroutine sort_subdomains_primal(krenu,ndoms,nelty,npoin,nelem,nnode,pelpo, &
!!$     &                              lelpo,ltree,ia,ja,pnods,lnods,ldomn,ldome, &
!!$     &                              mdomp,nobjs,lobjs,pdobj,ldobj,iperm,lperm,ieren,level)
!!$    implicit none
!!$    integer(ip), intent(in)    :: krenu,ndoms,nelty,npoin,nelem,nnode
!!$    integer(ip), intent(in)    :: pelpo(npoin+1),lelpo(pelpo(npoin+1))
!!$    integer(ip), intent(in)    :: ltree(ndoms)
!!$    integer(ip), intent(in)    :: pnods(*),lnods(*)
!!$
!!$    integer(ip), intent(inout) :: ia(npoin+1)
!!$    integer(ip), intent(inout) :: ja(*)
!!$
!!$    integer(ip), intent(in)    :: ldomn(npoin)
!!$    integer(ip), intent(inout) :: ldome(nelem)
!!$    integer(ip), intent(in)    :: mdomp,nobjs
!!$    integer(ip), intent(in)    :: lobjs(mdomp+4,nobjs)
!!$    integer(ip), intent(in)    :: pdobj(ndoms+1)
!!$    integer(ip), intent(in)    :: ldobj(pdobj(ndoms+1))
!!$    integer(ip), intent(inout) :: iperm(npoin)
!!$    integer(ip), intent(inout) :: lperm(npoin)
!!$    integer(ip), intent(inout) :: ieren(nelem)
!!$    integer(ip), intent(inout) :: level(*)
!!$
!!$    integer(ip)                :: idoms,ipoin,jpoin,ielem,ienew,jenew,ielpo
!!$    integer(ip)                :: inode,knode,inods,jnods,i,j,k,iz,ileve,nleve
!!$
!!$    ieren(1:nelem)=0
!!$
!!$    if(krenu==0) then
!!$
!!$       ! rcm renum
!!$       ienew=1
!!$       do idoms=1,ndoms
!!$          i=pdobj(idoms)
!!$          k=lobjs(2,ldobj(i))
!!$          level(1)=k
!!$          ! Loop over the nodes on the separator to define first level
!!$          ! TO DO: only first object as starting level
!!$          do i=pdobj(idoms)+1,pdobj(idoms+1)-1
!!$             !i=pdobj(idoms)+1
!!$             do j=lobjs(2,ldobj(i)),lobjs(3,ldobj(i))
!!$                ipoin=iperm(j)       ! This is the global node number
!!$                do iz=ia(ipoin),ia(ipoin+1)-1
!!$                   jpoin=ja(iz)
!!$                   if(ldomn(jpoin)==ltree(idoms).and.ia(jpoin)>0) then
!!$                      ia(jpoin)=-ia(jpoin) ! This is required by rcm_renum (as a touch)
!!$                      iperm(k)=jpoin
!!$                      k=k+1
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$          level(2)=k
!!$          ! do rcm using the nodal graph to obtain the nodal renumbering
!!$          ! using a mask ldomn(ipoin)=idoms
!!$          call rcm_renum_masked(npoin,ia(npoin+1),ia,ja,ldomn,ltree(idoms),iperm,nleve,level)
!!$          ! Reverse node numbering
!!$          do i=1,(level(nleve+1)-level(1))/2
!!$             j=iperm(level(1)+i-1)
!!$             iperm(level(1)+i-1)=iperm(level(nleve+1)-i)
!!$             iperm(level(nleve+1)-i)=j
!!$          end do
!!$
!!$          ! Finally recompute permutation from its inverse (on interior nodes only)
!!$          i=pdobj(idoms)
!!$          do j=lobjs(2,ldobj(i)),lobjs(3,ldobj(i))
!!$             ia(iperm(j))=-ia(iperm(j))
!!$             lperm(iperm(j))=j
!!$          end do
!!$
!!$          ! Loop on the levels and number the elements touching nodes on each level.
!!$          jenew=ienew
!!$          do ileve=1,nleve
!!$             do j=level(ileve),level(ileve+1)-1
!!$                ipoin=iperm(j)   ! This is the original node number
!!$                do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
!!$                   ielem=lelpo(ielpo)
!!$                   if(ldome(ielem)==idoms) then
!!$                      ieren(jenew)= ielem
!!$                      ldome(ielem)=0       ! Touch this element
!!$                      jenew=jenew+1
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$          ! Untouch elements
!!$          do i=ienew,jenew-1
!!$             ielem=ieren(i)
!!$             ldome(ielem)=idoms
!!$          end do
!!$          ienew=jenew
!!$       end do
!!$
!!$    else if(krenu==1) then
!!$
!!$       ! keep original (nd) renum of nodes but renumber elements accordingly
!!$       ienew=1
!!$       do idoms=1,ndoms
!!$          ! Loop over interior nodes and touch elements
!!$          jenew=ienew
!!$          i=pdobj(idoms)
!!$          do j=lobjs(2,ldobj(i)),lobjs(3,ldobj(i))
!!$             ipoin=iperm(j)            ! This is the global node number
!!$             do ielpo=pelpo(ipoin),pelpo(ipoin+1)-1
!!$                ielem=lelpo(ielpo)
!!$                if(ldome(ielem)==idoms) then
!!$                   ieren(jenew)= ielem
!!$                   ldome(ielem)=0      ! Touch this element
!!$                   jenew=jenew+1
!!$                end if
!!$             end do
!!$          end do
!!$          ! Untouch elements
!!$          do i=ienew,jenew-1
!!$             ielem=ieren(i)
!!$             ldome(ielem)=idoms
!!$          end do
!!$          ienew=jenew
!!$       end do
!!$
!!$    end if
!!$
!!$
!!$  end subroutine sort_subdomains_primal
!!$

!============================================================================================
!
! The following two routines implement stupid n^2 sort algorithms. They 
! should be used for n small (<10, say).
!
!============================================================================================
  !-----------------------------------------------------------------------
  ! This routine sorts the array l in place
  !-----------------------------------------------------------------------
  subroutine simple_sort (n, l)
    implicit none
    integer(ip), intent(in)    :: n
    integer(ip), intent(inout) :: l(n)
    integer(ip) :: i,j,k

    do i=1,n-1
       do j=i+1,n
          if (l(i)>l(j)) then
             k    = l(i)
             l(i) = l(j) 
             l(j) = k
          end if
       end do
    end do
  end subroutine simple_sort


!-----------------------------------------------------------------------
! This routine sorts the array l in place and eliminates repeated 
! elements
!-----------------------------------------------------------------------
  subroutine sort_eliminate_repeated(n,m,l)
    implicit none
    integer(ip), intent(in)    :: n
    integer(ip), intent(out)   :: m
    integer(ip), intent(inout) :: l(n)
    integer(ip) :: i,j,k

    do i=1,n-1
       do j=i+1,n
          if (l(i)>l(j)) then
             k    = l(i)
             l(i) = l(j) 
             l(j) = k
          end if
       end do
    end do

    i=1
    m=n
    do while(i<m)
       if(l(i+1)==l(i)) then
          do j=i+1,m-1
             l(j)=l(j+1)
          end do
          m=m-1
       else
          i=i+1
       end if
    end do

  end subroutine sort_eliminate_repeated



end module fem_mesh_partition
