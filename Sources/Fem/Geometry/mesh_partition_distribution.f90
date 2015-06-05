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
module fem_mesh_partition_distribution
  use types
  use memor
  use fem_mesh_distribution_names
  use maps_names
  use map_apply
  use fem_graph_names
  use graph_renum
  use fem_mesh_names
  use mesh_graph
  use fem_mesh_partition_base
  use hash_table_names
# include "debug.i90"
  implicit none
  private

   ! Functions
  public :: fem_mesh_distribution_create !, build_partition_adjacency

contains

  subroutine fem_mesh_distribution_create( prt_pars, femesh, distr, lmesh)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(part_params)                        , intent(in)  :: prt_pars
    type(fem_mesh)                           , intent(in)  :: femesh
    type(fem_mesh_distribution) , allocatable, intent(out) :: distr(:) ! Mesh distribution instances
    type(fem_mesh)              , allocatable, intent(out) :: lmesh(:) ! Local mesh instances

    ! Local variables
    type(fem_mesh)               :: dual_femesh, dual_lmesh
    type(fem_graph)              :: fe_graph    ! Dual graph (to be partitioned)
    integer(ip)   , allocatable  :: ldome(:)    ! Part of each element
    integer(ip)   , allocatable  :: dual_parts(:)
    integer(ip)                  :: ipart
    integer                      :: istat
    ! AFM. Dummy map declared and used for backward compatibility.
    ! It should be re-considered when we decide how to implement
    ! the boundary mesh.
    type(map)                    :: dummy_bmap


    ! Generate dual mesh (i.e., list of elements around points)
    call mesh_to_dual(femesh, dual_femesh)
    
    ! Allocate working arrays
    call memalloc (femesh%nelem, ldome, __FILE__,__LINE__)

    ! dual graph (elements around elements)
    call mesh_to_graph (dual_femesh, femesh, dual_femesh%ndime, fe_graph)
    !out0(call fem_graph_print(6, fe_graph))
          
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
       call fem_mesh_g2l(distr(ipart)%nmap, distr(ipart)%emap, dummy_bmap, femesh, lmesh(ipart))

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

       call fem_mesh_free(dual_lmesh)
       call memfree(dual_parts,__FILE__,__LINE__)
    end do

    call fem_mesh_free(dual_femesh)
    call fem_graph_free(fe_graph)
    call memfree(ldome,__FILE__,__LINE__)
  end subroutine fem_mesh_distribution_create

  !================================================================================================
   subroutine build_adjacency ( my_part, lmesh, l2ge, dual_lmesh, dual_parts, &
        &                       nebou, nnbou, lebou, lnbou, pextn, lextn, lextp)
     implicit none
     integer(ip)   , intent(in)  :: my_part
     type(fem_mesh), intent(in)  :: lmesh
     type(fem_mesh), intent(in)  :: dual_lmesh
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
    type(map_igp) , intent(in)  :: nmap
    type(fem_mesh), intent(in)  :: dual_mesh
    integer(ip)   , intent(in)  :: ldome(dual_mesh%npoin)
    type(fem_mesh), intent(in)  :: lmesh
    type(fem_mesh), intent(inout) :: dual_lmesh
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
    type(fem_mesh)             , intent(in)    :: femesh
    integer(ip)                , intent(in)    :: ldome(femesh%nelem)
    type(fem_mesh_distribution), intent(inout) :: distr(nparts)
    
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

end module fem_mesh_partition_distribution
