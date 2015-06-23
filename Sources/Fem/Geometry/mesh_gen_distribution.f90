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
module fem_mesh_gen_distribution_names
  use types
  use fem_conditions_names
  use fem_triangulation_names
  use fem_mesh_distribution_names
  use fem_materials_names
  use fem_space_types
  use maps_names
  implicit none
# include "debug.i90"
  private

  type geom_data
     integer(ip)            :: &
          ntdix=0,             &         ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
          ntdiy=0,             &         ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
          ntdiz=0,             &         ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
          mater=0,             &         ! Material case
          neblx=0,             &         ! Number of elements in the x boundary layer
          nebly=0,             &         ! Number of elements in the y boundary layer
          neblz=0,             &         ! Number of elements in the z boundary layer
          ndime,               &         ! Dimension 
          nparts=1,            &         ! Number of partitions
          nedir(3),            &         ! Number of elements in each direction
          npdir(3)=(/1,1,1/),  &         ! Number of parts on each direction           
          nsckt(3)=(/1,1,1/),  &         ! Number of parts on each socket and direction
          isper(3)=(/0,0,0/)             ! Flag for periodic boundary conditions on each direction

     real(rp)               :: &
          xleng   = 1.0_rp,    &         ! Size of the domain in x
          yleng   = 1.0_rp,    &         ! Size of the domain in y
          zleng   = 1.0_rp,    &         ! Size of the domain in z
          zx1     = 0.1_rp,    &         ! size of the elements at x=0   (left)  
          zx2     = 0.1_rp,    &         ! size of the elements at x=a/2 (center)
          zy1     = 0.1_rp,    &         ! size of the elements at y=0   (bottom)
          zy2     = 0.1_rp,    &         ! size of the elements at y=b/2 (center)
          x0      = 0.0_rp,    &         ! Origin x-coordinate
          y0      = 0.0_rp,    &         ! Origin y-coordinate
          z0      = 0.0_rp,    &         ! Origin z-coordinate
          xstret  = 2.75_rp,   &         ! Stretching parameter
          ystret  = 2.75_rp,   &         ! Stretching parameter
          zstret  = 2.75_rp,   &         ! Stretching parameter
          xlengbl = 0.0_rp,    &         ! Size of the boundary layer in x
          ylengbl = 0.0_rp,    &         ! Size of the boundary layer in y
          zlengbl = 0.0_rp               ! Size of the boundary layer in z
  end type geom_data     

  type bound_data
     type(fem_conditions)   :: &
          poin,                &         ! Boundary conditions on geometry conrners
          line,                &         ! Boundary conditions on geometry edges
          surf                           ! Boundary conditions on geometry faces
  end type bound_data
  
  type geom_size 
     integer(ip)            :: &
          nnode,               &         ! Number of nodes on each element
          nedir(3),            &         ! Number of elements in each direction
          npdir(3),            &         ! Number of parts on each direction           
          nsckt(3),            &         ! Number of parts on each socket and direction
          npsoc(3),            &         ! Number of sockets on each direction
          nedom(3),            &         ! Number of elements on each part and direction
          npdom(3),            &         ! Number of points on each part and direction
          ncorn(3),            &         ! Number of partition corners on each direction
          nedge(3,3),          &         ! Number of partition edges on each direction
          nface(3,3),          &         ! Number of partition faces on each direction
          npdomt,              &         ! Number of points on each domain
          nedomt,              &         ! Number of elements on each domain
          ncornt,              &         ! Number of partition corners on each domain
          nedget(3),           &         ! Number of edges on each direction
          nfacet(3),           &         ! Number of faces on each direction
          nedgett,             &         ! Number of edges on each domain
          nfacett,             &         ! Number of faces on each domain
          neghost                        ! Maximum of elements on each domain counting ghosts
  end type geom_size

  type topo_size
     integer(ip)            :: &
          notot,               &         ! Total amount of elemental objects of the partition
          nctot,               &         ! Total amount of elemental corners of the partition
          ncglb,               &         ! Total amount of elemental corners of the domain
          noglb,               &         ! Total amount of elemental corners of the domain
          neglb,               &         ! Total amount of elements of the domain
          ndtot,               &         ! Total amount of elemental edges of the partition
          nftot,               &         ! Total amount of elemental faces of the partition
          nddir(3),            &         ! Total amount of elemental edges of the partition for each direction
          nfdir(3),            &         ! Total amount of elemental faces of the partition for each direction
          nddom(3,3),          &         ! # edges for each direction given the edge direction (local)
          nfdom(3,3),          &         ! # faces for each direction given the face normal direction (local)
          ndglb(3,3),          &         ! # edges for each direction given the edge direction (global)
          nfglb(3,3),          &         ! # faces for each direction given the face normal direction (global)
          ndsum(3),            &         ! Total amount of elemental edges of the domain for each direction
          nfsum(3)                       ! Total amount of elemental faces of the domain for each direction
  end type topo_size

  ! Declare constants
  integer(ip), parameter :: geom =0, topo=1 
  integer(ip), parameter :: inter=0, bound=1
  integer(ip), parameter :: do_count=0, do_list=1

  interface globalid
     module procedure globalid_ip, globalid_igp
  end interface

  ! Types
  public :: geom_data, bound_data

  ! Functions
  public :: geom_data_create, bound_data_create, gen_triangulation
  public :: bound_data_free

contains

  !==================================================================================================
  subroutine geom_data_create(gdata,nex,ney,nez,npx,npy,npz,nsx,nsy,nsz,ntdix,ntdiy,ntdiz,mater,  &
       &                      neblx,nebly,neblz,perix,periy,periz,lx,ly,lz,x0,y0,z0,zx1,zx2,zy1, &
       &                      zy2,xstret,ystret,zstret,xlengbl,ylengbl,zlengbl)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates geometry data to construct a structured mesh                      !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(geom_data)      , intent(out) :: gdata
    integer(ip)          , intent(in)  :: nex,ney,nez
    integer(ip), optional, intent(in)  :: npx,npy,npz,nsx,nsy,nsz
    integer(ip), optional, intent(in)  :: ntdix,ntdiy,ntdiz
    integer(ip), optional, intent(in)  :: mater
    integer(ip), optional, intent(in)  :: neblx,nebly,neblz
    integer(ip), optional, intent(in)  :: perix,periy,periz
    real(rp)   , optional, intent(in)  :: lx,ly,lz
    real(rp)   , optional, intent(in)  :: x0,y0,z0  
    real(rp)   , optional, intent(in)  :: zx1,zx2
    real(rp)   , optional, intent(in)  :: zy1,zy2  
    real(rp)   , optional, intent(in)  :: xstret,ystret,zstret  
    real(rp)   , optional, intent(in)  :: xlengbl,ylengbl,zlengbl 

    ! Fill non default geom_data
    if(present(ntdix))   gdata%ntdix   = ntdix   ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    if(present(ntdiy))   gdata%ntdiy   = ntdiy   ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    if(present(ntdiz))   gdata%ntdiz   = ntdiz   ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
    if(present(mater))   gdata%mater   = mater   ! Material case
    if(present(neblx))   gdata%neblx   = neblx   ! Number of elements in the x boundary layer
    if(present(nebly))   gdata%nebly   = nebly   ! Number of elements in the y boundary layer
    if(present(neblz))   gdata%neblz   = neblz   ! Number of elements in the z boundary layer
    if(present(lx))      gdata%xleng   = lx      ! Size of the domain in x
    if(present(ly))      gdata%yleng   = ly      ! Size of the domain in y
    if(present(lz))      gdata%zleng   = lz      ! Size of the domain in z
    if(present(zx1))     gdata%zx1     = zx1     ! size of the elements at x=0   (left)  
    if(present(zx2))     gdata%zx2     = zx2     ! size of the elements at x=a/2 (center)
    if(present(zy1))     gdata%zy1     = zy1     ! size of the elements at y=0   (bottom)
    if(present(zy2))     gdata%zy2     = zy2     ! size of the elements at y=b/2 (center)
    if(present(x0))      gdata%x0      = x0      ! Origin x-coordinate
    if(present(y0))      gdata%y0      = y0      ! Origin y-coordinate
    if(present(z0))      gdata%z0      = z0      ! Origin z-coordinate
    if(present(xstret))  gdata%xstret  = xstret  ! Stretching parameter
    if(present(ystret))  gdata%ystret  = ystret  ! Stretching parameter
    if(present(zstret))  gdata%zstret  = zstret  ! Stretching parameter
    if(present(xlengbl)) gdata%xlengbl = xlengbl ! Size of the boundary layer in x
    if(present(ylengbl)) gdata%ylengbl = ylengbl ! Size of the boundary layer in y
    if(present(zlengbl)) gdata%zlengbl = zlengbl ! Size of the boundary layer in z
  
    ! Dimension
    if(nez==0) then
       gdata%ndime=2                
    else
       gdata%ndime=3
    end if 

    ! Partition info
    gdata%nedir(1) = nex
    gdata%nedir(2) = ney
    gdata%nedir(3) = nez
    if(present(npx)) gdata%npdir(1) = npx
    if(present(npy)) gdata%npdir(2) = npy
    if(present(npz)) gdata%npdir(3) = npz
    if(present(nsx)) then
       check(nsx>0)
       gdata%nsckt(1) = gdata%npdir(1)/nsx
    else
       gdata%nsckt(1) = gdata%npdir(1)
    end if
    if(present(nsy)) then
       check(nsy>0)
       gdata%nsckt(2) = gdata%npdir(2)/nsy
    else
       gdata%nsckt(2) = gdata%npdir(2)
    end if
    gdata%nparts = gdata%npdir(1)*gdata%npdir(2)
    if(gdata%ndime==3) then
       if(present(nsz)) then
          check(nsz>0)
          gdata%nsckt(3) = gdata%npdir(3)/nsz
       else
          gdata%nsckt(3) = gdata%npdir(3)
       end if
       gdata%nparts = gdata%nparts*gdata%npdir(3)
    else
       gdata%nsckt(3) = 1
       gdata%npdir(3) = 1
    end if
    if(present(perix)) gdata%isper(1) = perix
    if(present(periy)) gdata%isper(2) = periy
    if(present(periz)) gdata%isper(3) = periz    
    
  end subroutine geom_data_create

  !==================================================================================================
  subroutine bound_data_create(ncode,nvalu,ndime,bdata)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates boundary conditions for the geometrycal entities of a structured  !
    !   domain. By default, homogeneous Dirichlet BCs are condidered in all boundary                !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)     , intent(in)  :: ncode,nvalu,ndime
    type(bound_data), intent(out) :: bdata

    ! 2D:
    ! Points code: 1=bottom left; 2=top left;   3=bottom right; 4=top right
    ! Lines code:  1=left line;   2=right line; 3=bottom line;  4=top line 

    ! 3D:
    ! Points code:  1=(0,0,0), 2=(0,0,1), 3=(0,1,0), 4=(0,1,1), 5=(1,0,0),
    !               6=(1,0,1), 7=(1,1,0), 8=(1,1,1)
    ! Lines code:   1=(0,0,z), 2=(0,1,z), 3=(1,0,z), 4=(1,1,z), 5=(0,y,0),
    !               6=(0,y,1), 7=(1,y,0), 8=(1,y,1), 9=(x,0,0), 10=(x,0,1),
    !               11=(x,1,0), 12=(x,1,1)
    ! Surface code: 1=(x,y,0), 2=(x,y,1), 3=(x,0,z), 4=(x,1,z), 5=(0,y,z),
    !               6=(1,y,z)
    
    ! Create fem conditions
    call fem_conditions_create(ncode,nvalu,2**ndime,bdata%poin)
    call fem_conditions_create(ncode,nvalu,2*(ndime-1)*ndime,bdata%line)
    if(ndime==2) then
       call fem_conditions_create(ncode,nvalu,1,bdata%surf)
    else
       call fem_conditions_create(ncode,nvalu,2*ndime,bdata%surf)
    end if
    
    ! Initialize code and value
    bdata%poin%code(:,:) = 1
    bdata%line%code(:,:) = 1
    bdata%surf%code(:,:) = 1
    bdata%poin%valu(:,:) = 0.0_rp
    bdata%line%valu(:,:) = 0.0_rp
    bdata%surf%valu(:,:) = 0.0_rp   

  end subroutine bound_data_create  

  !==================================================================================================
  subroutine bound_data_free(bdata)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates boundary conditions for the geometrycal entities of a           !
    !   structured domain                                                                           !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(bound_data), intent(inout) :: bdata
    
    ! Free fem conditions
    call fem_conditions_free(bdata%poin)
    call fem_conditions_free(bdata%line)
    call fem_conditions_free(bdata%surf)

  end subroutine bound_data_free

  !==================================================================================================
  subroutine structured_geom_size_create(gdata,ginfo,gsize)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates geom_size type from geom_data and fem_fixed_info types            !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(geom_data)     , intent(in)  :: gdata
    type(fem_fixed_info), intent(in)  :: ginfo
    type(geom_size)     , intent(out) :: gsize

    ! Local variables
    integer(ip) :: pdegr,idime,jdime

    pdegr = ginfo%order

    ! Directional sizes
    gsize%npdir = gdata%npdir
    gsize%nedir = gdata%nedir
    gsize%nsckt = gdata%nsckt
    gsize%nnode = (pdegr+1)**gdata%ndime  
    gsize%npsoc = gdata%npdir/gdata%nsckt   
    gsize%nedom = gdata%nedir/gdata%npdir   
    gsize%npdom = pdegr*gsize%nedom + 1  
    gsize%ncorn = gdata%npdir + (1-gdata%isper) 
    do idime=1,gdata%ndime
       gsize%nedge(idime,:) = gdata%npdir + 1 - gdata%isper
       gsize%nface(idime,:) = gdata%npdir
       if(gdata%isper(idime)==0) then
          gsize%nface(idime,idime) = gsize%nface(idime,idime) + 1 
       end if
    end do
    if(gdata%ndime==2) gsize%nface = 0

    ! Total sizes
    gsize%npdomt  = 1
    gsize%nedomt  = 1
    gsize%ncornt  = 1
    gsize%nedget  = 1
    gsize%nfacet  = 1
    gsize%nedgett = 0
    gsize%nfacett = 0
    gsize%neghost = 1
    do idime=1,gdata%ndime
       gsize%npdomt = gsize%npdomt*gsize%npdom(idime)
       gsize%nedomt = gsize%nedomt*gsize%nedom(idime)
       gsize%ncornt = gsize%ncornt*gsize%ncorn(idime)
       gsize%neghost = gsize%neghost*(gsize%nedom(idime)+2)
       do jdime=1,gdata%ndime
          gsize%nedget(idime) = gsize%nedget(idime)*gsize%nedge(idime,jdime)
          gsize%nfacet(idime) = gsize%nfacet(idime)*gsize%nface(idime,jdime) 
       end do
       gsize%nedgett = gsize%nedgett + gsize%nedget(idime)
       gsize%nfacett = gsize%nfacett + gsize%nfacet(idime)
    end do
    
  end subroutine structured_geom_size_create

  !==================================================================================================
  subroutine structured_topo_size_create(gdata,gsize,tsize)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates topo_size type from geom_data and geom_size types            !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    type(geom_data), intent(in)  :: gdata
    type(geom_size), intent(in)  :: gsize
    type(topo_size), intent(out) :: tsize

    ! Local variables
    integer(ip) :: pdime,i,nfaux,ndaux
    integer(ip), allocatable :: auxv(:,:)
    
    ! Auxiliar vector of components
    if(gdata%ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
    else if(gdata%ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
    end if

    tsize%notot = 0
    tsize%nctot = 1
    tsize%ncglb = 1
    tsize%neglb = 1
    tsize%ndtot = 0
    tsize%nftot = 0
    do pdime=gdata%ndime,1,-1
       nfaux = 1
       ndaux = 1
       do i=1,gdata%ndime-1
          nfaux = nfaux*gsize%nedom(auxv(pdime,i))
          ndaux = ndaux*(gsize%nedom(auxv(pdime,i))+1)
          tsize%nfdom(auxv(pdime,i),pdime) = gsize%nedom(auxv(pdime,i))
          tsize%nfglb(auxv(pdime,i),pdime) = gsize%nedir(auxv(pdime,i))          
          tsize%nddom(auxv(pdime,i),pdime) = (gsize%nedom(auxv(pdime,i))+1)
          tsize%ndglb(auxv(pdime,i),pdime) = (gsize%nedir(auxv(pdime,i))+1)        
       end do
       tsize%nfdir(pdime) = (gsize%nedom(pdime)+1)*nfaux
       tsize%nddir(pdime) = gsize%nedom(pdime)*ndaux
       tsize%nfdom(pdime,pdime) = (gsize%nedom(pdime)+1)
       tsize%nfglb(pdime,pdime) = (gsize%nedir(pdime)+1)
       tsize%nddom(pdime,pdime) = gsize%nedom(pdime)
       tsize%ndglb(pdime,pdime) = gsize%nedir(pdime)
       tsize%ncglb = tsize%ncglb*(gsize%nedir(pdime)+1)
       tsize%neglb = tsize%neglb*gsize%nedir(pdime)
       tsize%nctot = tsize%nctot*gsize%npdom(pdime)
       tsize%ndtot = tsize%ndtot + tsize%nddir(pdime)
       tsize%nftot = tsize%nftot + tsize%nfdir(pdime)
    end do 
    do pdime=gdata%ndime,1,-1
       nfaux = 1
       ndaux = 1
       do i=1,gdata%ndime-1
          nfaux = nfaux*gsize%nedir(auxv(pdime,i))
          ndaux = ndaux*(gsize%nedir(auxv(pdime,i))+1)
       end do
       tsize%ndsum(pdime) = gsize%nedir(pdime)*ndaux
       tsize%nfsum(pdime) = (gsize%nedir(pdime)+1)*nfaux
    end do
    if(gdata%ndime==2) then
       tsize%nftot = 0
       tsize%nfdir = 0
       tsize%nfdom = 0
       tsize%nfglb = 0
    end if
    tsize%notot = tsize%nftot + tsize%ndtot + tsize%nctot
    tsize%noglb = tsize%notot*gdata%nparts

    ! Deallocate auxiliar vector of components
    call memfree(auxv,__FILE__,__LINE__)
    
  end subroutine structured_topo_size_create
  
  !==================================================================================================
  subroutine gen_triangulation(lpart,gdata,bdata,ginfo,trian,bcond,mater,mdist)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine generates a triangulation, boundary conditions and mesh_distribution for a  !
    !   structured mesh                                                                             !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)                          , intent(in)  :: lpart
    type(geom_data)                      , intent(in)  :: gdata
    type(bound_data)                     , intent(in)  :: bdata
    type(fem_fixed_info)                 , intent(in)  :: ginfo
    type(fem_triangulation)              , intent(out) :: trian
    type(fem_conditions)                 , intent(out) :: bcond
    integer(ip), allocatable             , intent(out) :: mater(:)
    type(fem_mesh_distribution), optional, intent(out) :: mdist
    
    ! Locals
    type(geom_size) :: gsize
    type(topo_size) :: tsize
    integer(ip)     :: ijkpart(3),ielem,inode

    ! Allocatables
    integer(ip) , allocatable   :: pextn(:),lextp(:),lexte(:),lextm(:)
    integer(igp), allocatable   :: lextn(:)

    ! Checks
    check(ginfo%ftype == Q_type_id)

    ! Geometrical sizes
    call structured_geom_size_create(gdata,ginfo,gsize)

    ! Topological sizes
    call structured_topo_size_create(gdata,gsize,tsize)

    ! Transform lpart to ijk numeration
    call global_to_ijk(lpart,gdata%nsckt(1:gdata%ndime),gsize%npsoc(1:gdata%ndime),gdata%ndime,ijkpart)

    ! Create triangulation
    if(present(mdist)) then
       call fem_triangulation_create(gsize%neghost,trian)
    else
       call fem_triangulation_create(gsize%nedomt,trian)
    end if
    trian%num_elems = gsize%nedomt
    trian%num_dims  = gdata%ndime

    ! Triangulation, bcond and distribution generation
    if(present(mdist)) then
       call structured_mesh_gen(ijkpart,gdata,gsize,tsize,ginfo,bdata%poin,bdata%line,bdata%surf,    &
            &                   trian,bcond,mater,nmap=mdist%nmap,emap=mdist%emap,pextn=mdist%pextn, &
            &                   lextn=mdist%lextn,lextp=mdist%lextp)
    else
       call structured_mesh_gen(ijkpart,gdata,gsize,tsize,ginfo,bdata%poin,bdata%line,bdata%surf, &
            &                   trian,bcond,mater)
    end if

    ! Dual triangulation
    if(.not.present(mdist)) call fem_triangulation_to_dual(trian)

    ! Create mesh_distribution
    if(present(mdist)) then
       mdist%ipart  = lpart
       mdist%nparts = gdata%nparts
       mdist%nebou  = mdist%emap%nb
       mdist%nnbou  = mdist%nmap%nb
       call memalloc(mdist%nebou,mdist%lebou,__FILE__,__LINE__)
       do ielem = mdist%emap%ni+1,mdist%emap%ni+mdist%emap%nb
          mdist%lebou(ielem-mdist%emap%ni) = ielem
       end do
       call memalloc(mdist%nnbou,mdist%lnbou,__FILE__,__LINE__)
       do inode = mdist%nmap%ni+1,mdist%nmap%ni+mdist%nmap%nb
          mdist%lnbou(inode-mdist%nmap%ni) = inode
       end do
    end if

  end subroutine gen_triangulation

  !================================================================================================
  subroutine structured_mesh_gen(ijkpart,gdata,gsize,tsize,ginfo,poin,line,surf,trian,nodes,mater, &
       &                         nmap,emap,pextn,lextn,lextp)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)                        , intent(in)    :: ijkpart(3)
    type(geom_data)                    , intent(in)    :: gdata
    type(geom_size)                    , intent(in)    :: gsize
    type(topo_size)                    , intent(in)    :: tsize
    type(fem_fixed_info)               , intent(in)    :: ginfo
    type(fem_triangulation)            , intent(inout) :: trian
    type(fem_conditions)               , intent(in)    :: poin,line,surf
    type(fem_conditions)               , intent(out)   :: nodes
    integer(ip), allocatable           , intent(out)   :: mater(:)
    type(map_igp)            , optional, intent(inout) :: nmap,emap
    integer(ip) , allocatable, optional, intent(out)   :: pextn(:),lextp(:)
    integer(igp), allocatable, optional, intent(out)   :: lextn(:)

    ! Local variables
    integer(ip)              :: nparts,ndime,isper(3)
    integer(ip)              :: cnt(3),subgl(gdata%ndime),idime,lpart,pni,eni,pnb,enb
    integer(ip)              :: pextn_cnt(2),nelbl(3)

    ! Local allocatables 
    integer(ip) , allocatable :: npnumg(:),npnumt(:),nenum(:)    
    integer(igp), allocatable :: l2ge(:),l2gp(:)   
    real(rp)    , allocatable :: coord(:,:)

    ! Unpack variables
    nparts = gdata%nparts
    ndime  = gdata%ndime
    isper  = gdata%isper

    ! Allocate 
    !call memalloc(gsize%npdomt,l2gp,__FILE__,__LINE__)
    call memalloc(tsize%notot,l2gp,__FILE__,__LINE__)
    call memalloc(gsize%nedomt,l2ge,__FILE__,__LINE__)
    call memalloc(gsize%nedomt,nenum,__FILE__,__LINE__)
    call memalloc(gsize%npdomt,npnumg,__FILE__,__LINE__)
    call memalloc(tsize%notot,npnumt,__FILE__,__LINE__)
    call memalloc(ndime,gsize%npdomt,coord,__FILE__,__LINE__)
    call memalloc(gsize%nedomt,mater,__FILE__,__LINE__)
    
    ! Initialize material
    mater = 1

    ! Initialize counter
    cnt(1) = 1; cnt(2) = 1; cnt(3) = 1

    ! Part id
    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

    ! Create bounary conditions
    call fem_conditions_create(poin%ncode,poin%nvalu,tsize%notot,nodes)
   
    ! Interior nodes/elements
    call volu_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,npnumg,npnumt,nenum,coord,cnt)
    call face_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,surf,npnumg,npnumt,nenum,coord, &
         &         cnt,inter,nodes)
    call edge_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,line,surf,npnumg,npnumt,nenum,  &
         &         coord,cnt,inter,nodes)
    call corn_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,poin,line,surf,npnumg,npnumt,   &
         &         nenum,coord,cnt,inter,nodes)
    !pni = cnt(3)-1
    pni = cnt(1)-1
    eni = cnt(2)-1

    ! Boundary nodes/elements
    call face_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,surf,npnumg,npnumt,nenum,coord, &
         &         cnt,bound)
    call edge_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,line,surf,npnumg,npnumt,nenum,  &
         &         coord,cnt,bound,nodes)
    call corn_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,gdata%isper,poin,line,surf,npnumg,npnumt,   &
         &         nenum,coord,cnt,bound,nodes)
    !pnb = cnt(3) - pni - 1
    pnb = cnt(1) - pni - 1
    enb = cnt(2) - eni - 1

    ! Global part numeration
    do idime=1,ndime
       subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
    end do

    ! Calculate lnods, coord and l2g vector for emap and nmap
    call generic_l2g(subgl,npnumg,npnumt,nenum,ndime,ginfo%order,gsize,tsize,gdata,l2ge,l2gp,trian, &
         &           coord,mater)

    ! Fill nmap
    if(present(nmap)) then
       !call map_alloc(gsize%npdomt,int(tsize%ncglb,igp),nmap)
       call map_alloc(tsize%notot,int(tsize%noglb,igp),nmap)
       nmap%ni  = pni
       nmap%nb  = pnb
       nmap%ne  = 0
       nmap%l2g = l2gp
    end if

    ! Fill emap
    if(present(emap)) then
       call map_alloc(gsize%nedomt,int(tsize%neglb,igp),emap)
       emap%ni  = eni
       emap%nb  = enb
       emap%ne  = 0
       emap%l2g = l2ge
    end if
    
    ! Compute adjacencies
    if(present(pextn)) then
       check(present(lextp))
       check(present(lextn))
       call memalloc(enb+1,pextn,__FILE__,__LINE__)
       nelbl(1) = gdata%neblx
       nelbl(2) = gdata%nebly
       nelbl(3) = gdata%neblz
       pextn_cnt = 0
       pextn(1) = 1
       call face_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
       call edge_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
       call corn_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,0)
       pextn_cnt = 0
       call memalloc(pextn(enb+1)-1,lextn,__FILE__,__LINE__)
       call memalloc(pextn(enb+1)-1,lextp,__FILE__,__LINE__)
       lextn = 0; lextp = 0
       call face_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
            &             lextn,lextp,gdata%mater,nelbl)
       call edge_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
            &             lextn,lextp,gdata%mater,nelbl)
       call corn_loop_adj(ijkpart,ndime,gsize,tsize,gdata%isper,emap%nb,pextn_cnt,pextn,1, &
            &             lextn,lextp,gdata%mater,nelbl)
    end if

    ! Deallocate auxiliar vectors
    call memfree( l2gp,__FILE__,__LINE__)
    call memfree( l2ge,__FILE__,__LINE__)
    call memfree(npnumt,__FILE__,__LINE__)
    call memfree(npnumg,__FILE__,__LINE__)
    call memfree(nenum,__FILE__,__LINE__)
    call memfree(coord,__FILE__,__LINE__)

  end subroutine structured_mesh_gen

  !==================================================================================================
  subroutine volu_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,npnumg,npnumt,nenum,coord,cnt)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in)    :: ijkpart(3),ndime
    type(geom_size)     , intent(in)    :: gsize
    type(topo_size)     , intent(in)    :: tsize
    type(geom_data)     , intent(in)    :: gdata
    type(fem_fixed_info), intent(in)    :: ginfo
    integer(ip)         , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
    real(rp)            , intent(inout) :: coord(:,:)

    integer(ip) :: i,j,k,ijkpoin(3),ijkelem(3),glnum,npdomk,nedomk,nddomk,aux_cnt
    integer(ip) :: pdime,ijkface(3),ijkedge(3),jdime
    integer(ip), allocatable :: auxv(:,:)

    ! Auxiliar
    if(ndime==2) then
       npdomk=2
       nedomk=2
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
    else if(ndime==3) then
       npdomk = gsize%npdom(3) - 1
       nedomk = gsize%nedom(3) - 1
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
    end if
    
    ! Interior nodes loop
    do i=2,gsize%npdom(1)-1
       do j=2,gsize%npdom(2)-1
          do k=2,npdomk

             ! Point identifier (inside the part)
             ijkpoin = (/i,j,k/)
             call globalid(ijkpoin,gsize%npdom,ndime,glnum)

             ! Global to local numeration (inside the part)
             npnumg(glnum) = cnt(3)
             npnumt(glnum) = cnt(1)

             ! Coordinates
             call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
                  &         ndime,gdata,ginfo,coord(:,cnt(3)))

             ! Update counter
             cnt(1) = cnt(1) + 1
             cnt(3) = cnt(3) + 1

          end do
       end do
    end do

    ! Interior edges loop
    aux_cnt = tsize%nctot
    do pdime=ndime,1,-1
       if(ndime==2) then
          nddomk=2
       elseif(ndime==3) then
          nddomk=tsize%nddom(auxv(pdime,2),pdime)-1
       end if
       do i=1,tsize%nddom(pdime,pdime)
          do j=2,tsize%nddom(auxv(pdime,1),pdime)-1
             do k=2,nddomk

                ! Face identifier (inside the part)
                ijkedge(pdime) = i
                if(ndime==2) then
                   ijkedge(auxv(pdime,:)) =  j
                else
                   ijkedge(auxv(pdime,:)) =  (/j,k/)
                end if
                call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)

                ! Global to local numeration (inside the part)
                npnumt(aux_cnt+glnum) = cnt(1)

                ! Update counter
                cnt(1) = cnt(1) + 1

             end do
          end do
       end do
       aux_cnt = aux_cnt + tsize%nddir(pdime)
    end do

    ! Interior faces loop
    aux_cnt = tsize%nctot + tsize%ndtot
    if(ndime==3) then
       do pdime=ndime,1,-1
          do i=2,tsize%nfdom(pdime,pdime)-1
             do j=1,tsize%nfdom(auxv(pdime,1),pdime)
                do k=1,tsize%nfdom(auxv(pdime,2),pdime)

                   ! Face identifier (inside the part)
                   ijkface(pdime) = i
                   ijkface(auxv(pdime,:)) = (/j,k/)
                   call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)

                   ! Global to local numeration (inside the part)
                   npnumt(aux_cnt+glnum) = cnt(1)

                   ! Update counter
                   cnt(1) = cnt(1) + 1

                end do
             end do
          end do
          aux_cnt = aux_cnt + tsize%nfdir(pdime)
       end do
    end if
    
    ! Interior elements loop
    do i=2,gsize%nedom(1)-1
       do j=2,gsize%nedom(2)-1
          do k=2,nedomk

             ! Element identifier (inside the part)
             ijkelem = (/i,j,k/)
             call globalid(ijkelem,gsize%nedom,ndime,glnum)

             ! Global to local numeration (inside the part)
             nenum(glnum) = cnt(2)

             ! Update counter
             cnt(2) = cnt(2) + 1

          end do
       end do
    end do

    ! Deallocate auxiliar vector
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine volu_loop

  !==================================================================================================
  subroutine face_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,isper,surf,npnumg,npnumt,nenum,coord, &
       &               cnt,case,nodes)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)                   , intent(in)    :: ijkpart(3),isper(3),ndime,case
    type(geom_size)               , intent(in)    :: gsize
    type(topo_size)               , intent(in)    :: tsize
    type(geom_data)               , intent(in)    :: gdata
    type(fem_fixed_info)          , intent(in)    :: ginfo
    type(fem_conditions)          , intent(in)    :: surf
    integer(ip)                   , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
    real(rp)                      , intent(inout) :: coord(:,:)
    type(fem_conditions), optional, intent(inout) :: nodes
    
    integer(ip) :: auxv(3,2),i,j,k,pdime,iface,ijkpoin(3),ijkelem(3),ijkface(3),glnum,flag,neigh(2)
    integer(ip) :: lface(2),surf_cnt
    integer(ip) :: idir,jdir,idime,jdime,aux_cnt_f,aux_cnt_d,ijkedge(3)

    ! Auxiliar vector of components
    auxv(1,:) = (/2,3/)
    auxv(2,:) = (/1,3/)
    auxv(3,:) = (/1,2/)

    ! Face nodes loop
    surf_cnt = 1
    aux_cnt_f = tsize%nctot + tsize%ndtot
    if(ndime==3) then
       do pdime = ndime,1,-1
          do iface = 0,1

             ! Check if it is a boundary
             call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
                  &                   isper,auxv,flag,neigh,lface)

             ! Check case
             if(flag==case) then

                ! Point identifier (inside the part)
                ijkpoin(pdime) = (gsize%npdom(pdime)-1)*iface+1

                do i=2,gsize%npdom(auxv(pdime,1))-1
                   do j=2,gsize%npdom(auxv(pdime,2))-1

                      ! Point identifier (inside the part)
                      ijkpoin(auxv(pdime,:)) = (/i,j/)
                      call globalid(ijkpoin,gsize%npdom,ndime,glnum)

                      ! Global to local numeration (inside the part)
                      npnumg(glnum) = cnt(3)
                      npnumt(glnum) = cnt(1)

                      ! Coordinates
                      call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
                           &         ndime,gdata,ginfo,coord(:,cnt(3)))

                      ! Boundary conditions
                      if(case==inter) then
                         nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
                         nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
                      end if

                      ! Update point counter
                      cnt(1) = cnt(1) + 1
                      cnt(3) = cnt(3) + 1

                   end do
                end do

                ! Elemental edges loop
                do idir = 1,2
                   jdir = 3-idir
                   idime = auxv(pdime,idir)
                   i = (tsize%nddom(pdime,idime)-1)*iface+1
                   aux_cnt_d = tsize%nctot
                   do jdime=1,3-idime
                      aux_cnt_d = aux_cnt_d + tsize%nddir(4-jdime)
                   end do
                   do j=1,tsize%nddom(idime,idime)
                      do k=2,tsize%nddom(auxv(pdime,jdir),idime)-1

                         ! Face identifier (inside the part)
                         ijkedge(pdime) = i
                         ijkedge(auxv(pdime,idir)) = j
                         ijkedge(auxv(pdime,jdir)) = k
                         call globalid(ijkedge,tsize%nddom(:,idime),ndime,glnum)

                         ! Global to local numeration (inside the part)
                         npnumt(aux_cnt_d+glnum) = cnt(1)

                         ! Boundary conditions
                         if(case==inter) then
                            nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
                            nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
                         end if

                         ! Update counter
                         cnt(1) = cnt(1) + 1

                      end do
                   end do
                end do

                ! Elemental faces loop
                i = (tsize%nfdom(pdime,pdime)-1)*iface+1
                do j=1,tsize%nfdom(auxv(pdime,1),pdime)
                   do k=1,tsize%nfdom(auxv(pdime,2),pdime)

                      ! Face identifier (inside the part)
                      ijkface(pdime) = i
                      ijkface(auxv(pdime,:)) = (/j,k/)
                      call globalid(ijkface,tsize%nfdom(:,pdime),ndime,glnum)

                      ! Global to local numeration (inside the part)
                      npnumt(aux_cnt_f+glnum) = cnt(1)

                      ! Boundary conditions
                      if(case==inter) then
                         nodes%code(:,cnt(1)) = surf%code(:,surf_cnt)
                         nodes%valu(:,cnt(1)) = surf%valu(:,surf_cnt)
                      end if

                      ! Update counter
                      cnt(1) = cnt(1) + 1

                   end do
                end do

                if(iface.gt.min(1,abs(gsize%nedom(pdime)-1))) cycle

                ! Element identifier (inside the part)
                ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1

                do i=2,gsize%nedom(auxv(pdime,1))-1
                   do j=2,gsize%nedom(auxv(pdime,2))-1

                      ! Point identifier (inside the part)
                      ijkelem(auxv(pdime,:)) = (/i,j/)
                      call globalid(ijkelem,gsize%nedom,ndime,glnum)

                      ! Global to local numeration (inside the part)
                      nenum(glnum) = cnt(2)

                      ! Update element counter
                      cnt(2) = cnt(2) + 1

                   end do
                end do
             end if

             ! Update surface counter
             surf_cnt = surf_cnt + 1

          end do
          
          ! Update elemental faces counter
          aux_cnt_f = aux_cnt_f + tsize%nfdir(pdime)
       end do
    end if

  end subroutine face_loop

 !===================================================================================================
  subroutine face_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
       &                   mcase,nelbl)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
    type(geom_size)       , intent(in)    :: gsize
    type(topo_size)       , intent(in)    :: tsize
    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
    integer(ip) , optional, intent(inout) :: lextp(:)
    integer(igp), optional, intent(inout) :: lextn(:)
    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
    
    integer(ip)  :: auxv(3,2),ijkelem(3),ijkneigh(3),subgl(3),subgl_ijk(3),ijkpart_neigh(3)
    integer(ip)  :: i,j,pdime,iface,ielem,jelem,glnum,flag,gpart,lface(2),k,kedge,idime,neigh(2)
    integer(igp) :: gneigh 

    if(ndime==3) then

       ! Auxiliar vector of components
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)

       ! Face nodes loop
       do pdime = ndime,1,-1
          do iface = 0,1

             ! Check if it is a boundary
             call check_face_boundary(ijkpart,pdime,iface,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
                  &                   isper,auxv,flag,neigh,lface)

             ! Only boundary elements
             if(flag==bound) then

                if(case==do_count) then

                   do i=2,gsize%nedom(auxv(pdime,1))-1
                      do j=2,gsize%nedom(auxv(pdime,2))-1

                         ! Update boundary element counter
                         pextn_cnt(1) = pextn_cnt(1) + 1

                         do ielem = -1,1
                            do jelem = -1,1

                               ! Update boundary neighbour counter
                               pextn_cnt(2) = pextn_cnt(2) + 1

                            end do
                         end do

                         ! Fill pextn
                         pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1

                      end do
                   end do

                elseif(case==do_list) then

                   ! Element identifier (inside the part)
                   ijkelem(pdime) = (gsize%nedom(pdime)-1)*iface+1

                   ! Neighbour element identifier
                   ijkneigh(pdime) = (gsize%nedom(pdime)-1)*(1-iface)+1

                   ! Neighbour part
                   ijkpart_neigh = ijkpart
                   ijkpart_neigh(pdime) = ijkpart(pdime) + (2*iface-1)
                   call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
                   call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

                   ! Global part numeration
                   do i=1,ndime
                      subgl(i) = (ijkpart_neigh(i)-1)*gsize%nedom(i)
                   end do

                   do i=2,gsize%nedom(auxv(pdime,1))-1
                      do j=2,gsize%nedom(auxv(pdime,2))-1

                         ! Update boundary element counter
                         pextn_cnt(1) = pextn_cnt(1) + 1

                         ! Point identifier (inside the part)
                         ijkelem(auxv(pdime,:)) = (/i,j/)
                         call globalid(ijkelem,gsize%nedom,ndime,glnum)

                         do ielem = -1,1
                            do jelem = -1,1

                               ! Update boundary neighbour counter
                               pextn_cnt(2) = pextn_cnt(2) + 1

                               ! Neighbour element identifier (local)
                               ijkneigh(auxv(pdime,:)) = (/i+ielem,j+jelem/)
                               subgl_ijk = subgl + ijkneigh
                               call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

                               ! Fill adjacencies info
                               lextn(pextn_cnt(2)) = gneigh
                               lextp(pextn_cnt(2)) = gpart

                            end do
                         end do

                      end do
                   end do
                end if
             end if
          end do
       end do

    end if

  end subroutine face_loop_adj

  !==================================================================================================
  subroutine edge_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,isper,line,surf,npnumg,npnumt,nenum, &
       &               coord,cnt,case,nodes)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip),          intent(in)    :: ijkpart(3),isper(3),ndime,case
    type(geom_size),      intent(in)    :: gsize
    type(topo_size),      intent(in)    :: tsize
    type(geom_data),      intent(in)    :: gdata
    type(fem_fixed_info), intent(in)    :: ginfo
    type(fem_conditions), intent(in)    :: line,surf
    integer(ip),          intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
    real(rp),             intent(inout) :: coord(:,:)
    type(fem_conditions), intent(inout) :: nodes
 
    integer(ip) :: i,j,k,pdime,iedge,jedge,ijedge(2),ijkpoin(3),ijkelem(3),ijkvert(3),ijkedge(3)
    integer(ip) :: inipo,ledge(2,2),line_cnt,neigh(2*(ndime-1))
    integer(ip) :: glnum,flag,kron(ndime),gsurf,inode,aux_surf(4,2),gsurf_aux(4),aux_cnt
    integer(ip), allocatable :: auxv(:,:)

    ! Auxiliar vector of components 
    if(ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
    else if(ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
    end if

    ! Edge nodes loop
    line_cnt = 0
    aux_cnt = tsize%nctot
    do pdime = ndime,1,-1
       do iedge = 0,1!min(1,nedom(auxv(pdime,1))-1)
          do jedge =0,1!min(1,nedom(auxv(pdime,2))-1)

             ! Check if it is a boundary
             ijedge = (/iedge,jedge/)
             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc,ndime, &
                  &                   isper,auxv,flag,neigh,ledge)

             ! Update edge counter
             line_cnt = line_cnt + 1

             ! Check case
             if(flag==case) then

                ! Initial object point
                inipo = cnt(1)

                ! Point identifier (inside the part)
                ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
                ijkvert(pdime) = 0
                do i=1,ndime
                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
                end do

                do i=1,gsize%npdom(pdime)-2

                   ! Point identifier (inside the part)
                   call kronec(pdime,ndime,kron)
                   ijkpoin(1:ndime) = ijkpoin(1:ndime) + kron
                   call globalid(ijkpoin,gsize%npdom,ndime,glnum)

                   ! Global to local numeration (inside the part)
                   npnumg(glnum) = cnt(3)
                   npnumt(glnum) = cnt(1)

                   ! Coordinates
                   call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
                        &         ndime,gdata,ginfo,coord(:,cnt(3)))

                   ! Boundary conditions
                   if(case==inter) then
                      nodes%code(:,cnt(1)) = line%code(:,line_cnt)
                      nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
                   end if

                   ! Update point counter
                   cnt(1) = cnt(1) + 1
                   cnt(3) = cnt(3) + 1
                end do

                ! Elemental edges loop
                j = iedge*(tsize%nddom(auxv(pdime,1),pdime)-1) + 1
                if(ndime==3) k = jedge*(tsize%nddom(auxv(pdime,2),pdime)-1) + 1
                do i = 1,tsize%nddom(pdime,pdime)

                   ! Face identifier (inside the part)
                   ijkedge(pdime) = i
                   ijkedge(auxv(pdime,1)) = j
                   if(ndime==3) ijkedge(auxv(pdime,2)) = k
                   call globalid(ijkedge,tsize%nddom(:,pdime),ndime,glnum)

                   ! Global to local numeration (inside the part)
                   npnumt(aux_cnt+glnum) = cnt(1)

                   ! Boundary conditions
                   if(case==inter) then
                      nodes%code(:,cnt(1)) = line%code(:,line_cnt)
                      nodes%valu(:,cnt(1)) = line%valu(:,line_cnt)
                   end if

                   ! Update counter
                   cnt(1) = cnt(1) + 1

                end do

                if(case==bound) then

                   ! Boundary conditions for non-global edges
                   k = 0
                   do i=1,2*(ndime-1)
                      if(neigh(i)>0) then
                         k=k+1
                      end if
                   end do
                   if(ndime==3.and.k==2) then
                      call aux_surf_num(pdime,aux_surf,gsurf_aux)
                      do i=1,4
                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0) then
                            gsurf = gsurf_aux(i)
                         end if
                      end do
                      do inode=inipo,cnt(1)-1
                         nodes%code(:,inode) = surf%code(:,gsurf)
                         nodes%valu(:,inode) = surf%valu(:,gsurf)
                      end do
                   end if

                end if

                if(iedge.gt.min(1,abs(gsize%nedom(auxv(pdime,1))-1))) cycle
                if(ndime==3) then
                   if(jedge.gt.min(1,abs(gsize%nedom(auxv(pdime,2))-1))) cycle
                end if

                ! Element identifier (inside the part)
                do i=1,ndime
                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
                end do

                do i=2,gsize%nedom(pdime)-1

                   ! Element identifier (inside the part)
                   call kronec(pdime,ndime,kron)
                   ijkelem(1:ndime) = ijkelem(1:ndime) + kron
                   call globalid(ijkelem,gsize%nedom,ndime,glnum)

                   ! Global to local numeration (inside the part)
                   nenum(glnum) = cnt(2)

                   ! Update element counter
                   cnt(2) = cnt(2) + 1
                end do

             end if

             if(ndime==2) exit
          end do
       end do

       ! Update elemental edge counter
       aux_cnt = aux_cnt + tsize%nddir(pdime)
    end do

    ! Deallocate auxiliar vector
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine edge_loop

  !===================================================================================================
  subroutine edge_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
       &                   mcase,nelbl)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
    type(geom_size)       , intent(in)    :: gsize
    type(topo_size)       , intent(in)    :: tsize
    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
    integer(ip) , optional, intent(inout) :: lextp(:)
    integer(igp), optional, intent(inout) :: lextn(:)
    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
    
    integer(ip) :: ijkelem(3),neigh(2*(ndime-1)),ijedge(2),ijkvert(3),ijkneigh(3)
    integer(ip) :: subgl(3),subgl_ijk(3),ijkpart_neigh(3),ledge(2,2),kron(ndime),lpart
    integer(ip) :: i,j,pdime,iedge,jedge,kedge,medge,ielem,glnum,flag,gpart,idime
    integer(ip) :: iaux,jaux,iaux2,jaux2
    integer(igp) :: gneigh
    integer(ip), allocatable :: auxv(:,:)

    ! Auxiliar vector of components
    if(ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
    else if(ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
    end if

    ! Global numbering
    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

    ! Edge nodes loop
    do pdime = ndime,1,-1
       do iedge = 0,1
          iaux = 2*iedge-1
          do jedge =0,1
             jaux = 2*jedge-1

             ! Check if it is a boundary
             ijedge = (/iedge,jedge/)
             call check_edge_boundary(ijkpart,pdime,ijedge,gsize%npdir,gsize%nsckt,gsize%npsoc, &
                  &                   ndime,isper,auxv,flag,neigh,ledge)

             ! Only boundary elements
             if(flag==bound) then

                if(case==do_count) then

                   do i=2,gsize%nedom(pdime)-1

                      ! Update boundary element counter
                      pextn_cnt(1) = pextn_cnt(1) + 1

                      do ielem = -1,1
                         do kedge = -1,1
                            do medge = -1,1
                               if(kedge==iaux.or.medge==jaux) then

                                  iaux2 = kedge
                                  jaux2 = medge
                                  if(kedge.ne.iaux) iaux2 = 0
                                  if(medge.ne.jaux) jaux2 = 0

                                  ! Neighbour part
                                  ijkpart_neigh(pdime) = ijkpart(pdime)
                                  if(ndime==2) then
                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
                                  else
                                     ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
                                  end if
                                  call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
                                  call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

                                  if(gpart>0.and.gpart.ne.lpart) then

                                     ! Update boundary neighbour counter
                                     pextn_cnt(2) = pextn_cnt(2) + 1

                                  end if
                               end if
                               if(ndime==2) exit
                            end do
                         end do
                      end do

                      ! Fill pextn
                      pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1             

                   end do

                elseif(case==do_list) then

                   ! Element identifier (inside the part)
                   ijkvert(auxv(pdime,:)) = ijedge(1:ndime-1)
                   ijkvert(pdime) = 0
                   do i=1,ndime
                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
                   end do

                   do i=2,gsize%nedom(pdime)-1

                      ! Element identifier (inside the part)
                      call kronec(pdime,ndime,kron)
                      ijkelem(1:ndime) = ijkelem(1:ndime) + kron
                      call globalid(ijkelem,gsize%nedom,ndime,glnum)

                      do kedge = -1,1
                         do medge = -1,1
                            if(kedge==iaux.or.medge==jaux) then
                               
                               iaux2 = kedge
                               jaux2 = medge
                               if(kedge.ne.iaux) iaux2 = 0
                               if(medge.ne.jaux) jaux2 = 0

                               ! Neighbour part
                               ijkpart_neigh(pdime) = ijkpart(pdime)
                               if(ndime==2) then
                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + iaux2
                               else
                                  ijkpart_neigh(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + (/iaux2,jaux2/)
                               end if

                               ! Global part numeration
                               do idime=1,ndime
                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
                                  if(ijkpart_neigh(idime)<1) then
                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
                                     subgl(idime) = -gsize%nedom(idime)
                                  end if
                               end do

                               ! Check part boundary
                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)

                               if(gpart>0.and.gpart.ne.lpart) then

                                  do ielem = -1,1

                                     ! Update boundary neighbour counter
                                     pextn_cnt(2) = pextn_cnt(2) + 1

                                     ! Neighbour element identifier (local)
                                     ijkneigh(pdime) = ijkelem(pdime) + ielem
                                     if(ndime==2) then
                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) + kedge
                                     else
                                        ijkneigh(auxv(pdime,:)) = ijkelem(auxv(pdime,:)) +(/kedge,medge/)
                                     end if
                                     subgl_ijk = subgl + ijkneigh
                                     call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

                                     ! Fill adjacencies info
                                     lextn(pextn_cnt(2)) = gneigh
                                     lextp(pextn_cnt(2)) = gpart

                                  end do
                               end if
                            end if
                            if(ndime==2) exit
                         end do
                      end do
                   end do
                end if
             end if
             if(ndime==2) exit
          end do
       end do
    end do

    ! Deallocate auxiliar vector
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine edge_loop_adj

  !================================================================================================
  subroutine corn_loop(ijkpart,ndime,gsize,tsize,gdata,ginfo,isper,poin,line,surf,npnumg,npnumt, & 
       &               nenum,coord,cnt,case,nodes)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in)    :: ijkpart(3),isper(3),ndime,case
    type(geom_size)     , intent(in)    :: gsize
    type(topo_size)     , intent(in)    :: tsize
    type(geom_data)     , intent(in)    :: gdata
    type(fem_fixed_info), intent(in)    :: ginfo
    type(fem_conditions), intent(in)    :: poin,line,surf
    integer(ip)         , intent(inout) :: npnumg(:),npnumt(:),nenum(:),cnt(3)
    real(rp)            , intent(inout) :: coord(:,:)
    type(fem_conditions), intent(inout) :: nodes
    
    integer(ip) :: ivert,jvert,kvert,ijkpoin(3),ijkvert(3),glnum,flag,neigh(2**ndime)
    integer(ip) :: i,j,k,lcorn(2,2,2),ijkelem(3)
    integer(ip) :: poin_cnt,gline,gsurf
    integer(ip) :: aux_line(2*(ndime-1)*ndime,2),aux_surf(2*ndime,4)
    integer(ip), allocatable :: auxv(:,:)

    ! Auxiliar vector of components
    if(ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
       aux_line(:,1) = (/3,1,2,1/)
       aux_line(:,2) = (/4,2,4,3/)
    else if(ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
       aux_line(:,1) = (/7,5,3,1,6,5,2,1,4,3,2,1/)
       aux_line(:,2) = (/8,6,4,2,8,7,4,3,8,7,6,5/)
       aux_surf(:,1) = (/2,1,3,1,5,1/)
       aux_surf(:,2) = (/4,3,4,2,6,2/)
       aux_surf(:,3) = (/6,5,7,5,7,3/)
       aux_surf(:,4) = (/8,7,8,6,8,4/)
    end if

    ! Corner nodes loop
    poin_cnt = 1
    do ivert = 0,1
       do jvert = 0,1
          do kvert =0,1

             ! Check if it is a boundary
             ijkvert = (/ivert,jvert,kvert/)
             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
                  &                   ndime,isper,auxv,flag,neigh,lcorn)

             ! Check case
             if(flag==case) then

                ! Point identifier (inside the part)
                do i=1,ndime
                   ijkpoin(i) = ijkvert(i)*(gsize%npdom(i)-1) + 1
                end do
                call globalid(ijkpoin,gsize%npdom,ndime,glnum)

                ! Global to local numeration (inside the part)
                npnumg(glnum) = cnt(3)
                npnumt(glnum) = cnt(1)
              
                ! Coordinates
                call coord_ijk(ijkpoin,ijkpart,gsize%npdir,gsize%nedom,gsize%nedir, &
                     &         ndime,gdata,ginfo,coord(:,cnt(3)))

                ! Boundary conditions
                if(case==inter) then
                   nodes%code(:,cnt(1)) = poin%code(:,poin_cnt)
                   nodes%valu(:,cnt(1)) = poin%valu(:,poin_cnt)
                end if

                ! Update point counter
                cnt(1) = cnt(1) + 1
                cnt(3) = cnt(3) + 1

                if(case==bound) then

                   ! Boundary conditions for non-global corners
                   k = 0
                   do i=1,2**ndime 
                      if(neigh(i)>0) then
                         k=k+1
                      end if
                   end do
                   if(k==2) then
                      do i=1,2*(ndime-1)*ndime
                         if(neigh(aux_line(i,1))>0.and.neigh(aux_line(i,2))>0) then
                            gline = i
                         end if
                      end do
                      nodes%code(:,cnt(1)-1) = line%code(:,gline)
                      nodes%valu(:,cnt(1)-1) = line%valu(:,gline)
                   else if(k>2.and.k<2**ndime.and.ndime==3) then
                      do i=1,2*ndime
                         if(neigh(aux_surf(i,1))>0.and.neigh(aux_surf(i,2))>0.and. &
                              neigh(aux_surf(i,3))>0.and.neigh(aux_surf(i,4))>0) then
                            gsurf = i
                         end if
                      end do
                      nodes%code(:,cnt(1)-1) = surf%code(:,gsurf)
                      nodes%valu(:,cnt(1)-1) = surf%valu(:,gsurf)
                   end if

                end if

                if(ivert.gt.min(1,abs(gsize%nedom(1)-1))) cycle
                if(jvert.gt.min(1,abs(gsize%nedom(2)-1))) cycle
                if(kvert.gt.min(1,abs(gsize%nedom(3)-1))) cycle

                ! Element identifier (inside the part)
                do i=1,ndime
                   ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
                end do
                call globalid(ijkelem,gsize%nedom,ndime,glnum)

                ! Global to local numeration (inside the part)
                nenum(glnum) = cnt(2)

                ! Update element counter
                cnt(2) = cnt(2) + 1

             end if

             ! Update corner counter
             poin_cnt = poin_cnt + 1

             if(ndime==2) exit
          end do
       end do
    end do

    ! Deallocate auxiliar vector
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine corn_loop

  !==================================================================================================
  subroutine corn_loop_adj(ijkpart,ndime,gsize,tsize,isper,nb,pextn_cnt,pextn,case,lextn,lextp, &
       &                   mcase,nelbl)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)           , intent(in)    :: ijkpart(3),isper(3),ndime,case,nb
    type(geom_size)       , intent(in)    :: gsize
    type(topo_size)       , intent(in)    :: tsize
    integer(ip)           , intent(inout) :: pextn_cnt(2),pextn(nb+1)
    integer(ip) , optional, intent(inout) :: lextp(:)
    integer(igp), optional, intent(inout) :: lextn(:)
    integer(ip) , optional, intent(in)    :: mcase,nelbl(3)
    
    integer(ip) :: ijkelem(3),neigh(2**ndime),ijkvert(3),ijkneigh(3),subgl(3),subgl_ijk(3)
    integer(ip) :: i,j,k,ivert,jvert,kvert,glnum,flag,gpart,lcorn(2,2,2),idime,ijkpart_neigh(3)
    integer(ip) :: iaux,jaux,kaux,iaux2,jaux2,kaux2,lpart
    integer(ip), allocatable :: auxv(:,:)
    integer(igp) :: gneigh

    ! Auxiliar vector of components
    if(ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
    else if(ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
    end if       

    ! Global numbering
    call globalid_2l(ijkpart,gsize%nsckt,gsize%npsoc,ndime,lpart)

    !Corner nodes loop
    do ivert = 0,1
       iaux = 2*ivert-1
       do jvert = 0,1
          jaux = 2*jvert-1
          do kvert =0,1
             kaux = 2*kvert-1
    
             ! Check if it is a boundary
             ijkvert = (/ivert,jvert,kvert/)
             call check_corn_boundary(ijkpart,ijkvert,gsize%npdir,gsize%nsckt,gsize%npsoc, &
                  &                   ndime,isper,auxv,flag,neigh,lcorn)

             ! Only boundary elements
             if(flag==bound) then

                if(case==do_count) then
                   
                   ! Update boundary element counter
                   pextn_cnt(1) = pextn_cnt(1) + 1

                   do i = -1,1
                      do j = -1,1
                         do k = -1,1
                            if(i==iaux.or.j==jaux.or.k==kaux) then

                               iaux2 = i
                               jaux2 = j
                               kaux2 = k
                               if(i.ne.iaux) iaux2 = 0
                               if(j.ne.jaux) jaux2 = 0
                               if(k.ne.kaux) kaux2 = 0
                               
                               ! Neighbour part
                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)
                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
                               
                               if(gpart>0.and.gpart.ne.lpart) then
                                  
                                  ! Update boundary neighbour counter
                                  pextn_cnt(2) = pextn_cnt(2) + 1

                               end if

                            end if
                            if(ndime==2) exit
                         end do
                      end do
                   end do

                   ! Fill pextn
                   pextn(pextn_cnt(1)+1) = pextn_cnt(2) + 1

                elseif(case==do_list) then

                   ! Element identifier (inside the part)
                   do i=1,ndime
                      ijkelem(i) = ijkvert(i)*(gsize%nedom(i)-1) + 1
                   end do

                   ! Element identifier (inside the part)
                   call globalid(ijkelem,gsize%nedom,ndime,glnum)

                   do i = -1,1
                      do j = -1,1
                         do k = -1,1
                            if(i==iaux.or.j==jaux.or.k==kaux) then

                               iaux2 = i
                               jaux2 = j
                               kaux2 = k
                               if(i.ne.iaux) iaux2 = 0
                               if(j.ne.jaux) jaux2 = 0
                               if(k.ne.kaux) kaux2 = 0
                               
                               ! Neighbour part
                               ijkpart_neigh = ijkpart + (/iaux2,jaux2,kaux2/)

                               ! Global part numeration
                               do idime=1,ndime
                                  subgl(idime) = (ijkpart(idime)-1)*gsize%nedom(idime)
                                  if(ijkpart_neigh(idime)<1) then
                                     subgl(idime) = gsize%npdir(idime)*gsize%nedom(idime)
                                  elseif(ijkpart_neigh(idime)>gsize%npdir(idime)) then
                                     subgl(idime) = -gsize%nedom(idime)
                                  end if
                               end do

                               ! Check boundary partition
                               call check_part_boundary(ijkpart_neigh,gsize%npdir,ndime,isper)
                               call globalid_2l(ijkpart_neigh,gsize%nsckt,gsize%npsoc,ndime,gpart)
                               
                               if(gpart>0.and.gpart.ne.lpart) then

                                  ! Update boundary neighbour counter
                                  pextn_cnt(2) = pextn_cnt(2) + 1

                                  ! Neighbour element identifier (local)
                                  ijkneigh = ijkelem + (/i,j,k/)
                                  subgl_ijk = subgl + ijkneigh
                                  call globalid(subgl_ijk,gsize%nedir,ndime,gneigh)

                                  ! Fill adjacencies info
                                  lextn(pextn_cnt(2)) = gneigh
                                  lextp(pextn_cnt(2)) = gpart

                               end if
                            end if
                            if(ndime==2) exit
                         end do
                      end do
                   end do
                end if
             end if
             if(ndime==2) exit
          end do
       end do
    end do

    ! Deallocate auxiliar vector
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine corn_loop_adj

  !================================================================================================
  subroutine coord_ijk(ijkpoin,ijkpart,npdir,nedom,nedir,ndime,msize,ginfo,coord)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)         , intent(in)    :: ijkpoin(3),ijkpart(3),npdir(3),nedom(3),nedir(3)
    integer(ip)         , intent(in)    :: ndime
    type(geom_data)     , intent(in)    :: msize
    type(fem_fixed_info), intent(in)    :: ginfo
    real(rp)            , intent(inout) :: coord(:)

    integer(ip) :: i,ntdis(3),pdegr,nebl(3)
    real(rp)    :: leng(3),coor0(3),stret(3),istret,lengbl(3)
    integer(ip) :: newijkpoin,newnedir,newnedom,info
    real(rp)    :: newleng

    ! Unpack geom_data
    pdegr   = ginfo%order
    leng(1) = msize%xleng
    leng(2) = msize%yleng
    leng(3) = msize%zleng
    ntdis(1) = msize%ntdix
    ntdis(2) = msize%ntdiy
    ntdis(3) = msize%ntdiz
    stret(1) = msize%xstret
    stret(2) = msize%ystret
    stret(3) = msize%zstret
    nebl(1) = msize%neblx
    nebl(2) = msize%nebly
    nebl(3) = msize%neblz
    lengbl(1) = msize%xlengbl
    lengbl(2) = msize%ylengbl
    lengbl(3) = msize%zlengbl

    ! Set origin of coordinates
    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0

    ! Set origin of coordinates
    coor0(1) = msize%x0; coor0(2) = msize%y0; coor0(3) = msize%z0

    do i=1,ndime 
       assert ( ntdis(i).ne.1 ) ! This case is not implemented
       if(ntdis(i)==0) then
          coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)
       elseif(ntdis(i)==2) then
          istret=stret(i)
          coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-1)/real(pdegr))/nedir(i)))/tanh(istret)
       elseif(ntdis(i)==3.or.ntdis(i)==4) then ! boundary layers (Hunt's case: solid uniform + fluid unif==3 or tanh==4)
          istret=stret(i)
          
!!$          !TO DO
!!$          if((nebl(i)>0).and.(nedom(i)<nebl(i))) then
!!$             write(*,*) 'boundary layer cannot be splitted over several sbds, dir ',i,'nedom',nedom,'nebl',nebl
!!$             call mpi_barrier (1, info)
!!$             call mpi_finalize( info )
!!$             stop
!!$          end if

          !first boundary layer
          if((ijkpart(i)==1).and.(ijkpoin(i)<=nebl(i))) then 
             if(ntdis(i)==3) then
                coord(i) = coor0(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
             elseif(ntdis(i)==4) then
                coord(i) = coor0(i) - leng(i) - lengbl(i) + lengbl(i)*(ijkpoin(i)-1)/nebl(i)
             end if
          !latest boundary layer
          elseif((ijkpart(i)==npdir(i)).and.(ijkpoin(i)>(nedom(i)-nebl(i)+1))) then 
             newijkpoin=ijkpoin(i)-(nedom(i)-nebl(i))
             coord(i) = coor0(i) + leng(i) + lengbl(i)*(newijkpoin-1)/nebl(i)
          else !core region
             if(ntdis(i)==3) then
                coord(i) = coor0(i) + leng(i)*((ijkpart(i)-1)*nedom(i)+(ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))
             elseif(ntdis(i)==4) then
                coord(i) = coor0(i) - leng(i)*tanh(istret*(1.0_rp-2.0_rp*((ijkpart(i)-1)*nedom(i)+&
                     (ijkpoin(i)-nebl(i)-1)/real(pdegr))/(nedir(i)-2*nebl(i))))/tanh(istret)
             end if
          end if
       end if
    end do
   
  end subroutine coord_ijk

  !================================================================================================
  subroutine check_part_boundary(ijk,npdir,ndime,isper)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: npdir(3),ndime,isper(3)
    integer(ip), intent(inout) :: ijk(3)

    integer(ip) :: i

    ! Periodic mesh with respect to all directions
    do i=1,ndime
       if(isper(i)==1) then
          if(ijk(i)<1) then
             ijk(i) = npdir(i)
          else if(ijk(i)>npdir(i)) then
             ijk(i) = 1!npdir(i)
          end if
       else
          if(ijk(i)<1.or.ijk(i)>npdir(i)) then
             ijk = 0*ijk
          end if
       end if
    end do

  end subroutine check_part_boundary

  !================================================================================================
  subroutine check_face_boundary(ijkpart,pdime,iface,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
       &                         neigh,lface)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),pdime,iface,npdir(3),ndime,isper(3),auxv(3,2)
    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
    integer(ip), intent(out) :: flag,neigh(2),lface(2)

    integer(ip)              :: i,ijkneig(3),newpo,cnt

    cnt=1
    flag = -1
    lface = 0

    do i=-1,0

       ! Auxiliar variables
       newpo = iface + i

       ! Check neighbors
       call neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)

       ! Check part boundary
       call check_part_boundary(ijkneig,npdir,ndime,isper)

       ! Global numbering
       call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))
       
       ! Set flag
       if(neigh(cnt)>0) then
          flag = min(flag+1,1)
          lface(i+2) = 1 
       end if

       ! Update counter
       cnt = cnt + 1

    end do

  end subroutine check_face_boundary

  !================================================================================================
  subroutine check_edge_boundary(ijkpart,pdime,ijedge,npdir,nsckt,npsoc,ndime,isper,auxv,flag, &
       &                         neigh,ledge)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),pdime,ijedge(2),npdir(3),ndime,isper(3),auxv(:,:)
    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
    integer(ip), intent(out) :: flag,neigh(2*(ndime-1)),ledge(2,2)

    integer(ip)              :: i,j,ijkneig(3),newpo(2),cnt

    cnt=1
    flag = -1
    ledge = 0

    do i=-1,0
       do j=-1,0

          ! Auxiliar variables
          newpo = ijedge + (/i,j/)

          ! Check neighbors
          call neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)

          ! Check part boundary
          call check_part_boundary(ijkneig,npdir,ndime,isper)

          ! Global numbering
          call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))

          ! Set flag
          if(neigh(cnt)>0) then
             flag = min(flag+1,1)
             ledge(i+2,j+2) = 1
          end if

          ! Update counter
          cnt = cnt + 1

          if(ndime==2) exit
       end do
    end do

  end subroutine check_edge_boundary

  !================================================================================================
  subroutine check_corn_boundary(ijkpart,ijkvert,npdir,nsckt,npsoc,ndime,isper,auxv,flag,neigh, &
       &                         lcorn)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),ijkvert(3),npdir(3),ndime,isper(3),auxv(:,:)
    integer(ip), intent(in)  :: nsckt(3),npsoc(3)
    integer(ip), intent(out) :: flag,neigh(2**ndime),lcorn(2,2,2)

    integer(ip)              :: i,j,k,ijkneig(3),newpo(3),cnt

    cnt = 1
    flag = -1
    lcorn = 0

    do i=-1,0
       do j=-1,0
          do k=-1,0

             ! Auxiliar variables
             newpo = ijkvert + (/i,j,k/)

             ! Check neighbors
             call neigh_corn(newpo,ijkpart,ijkneig)

             ! Check part boundary
             call check_part_boundary(ijkneig,npdir,ndime,isper)

             ! Global numbering
             call globalid_2l(ijkneig,nsckt,npsoc,ndime,neigh(cnt))

             ! Set flag
             if(neigh(cnt)>0) then
                flag = min(flag+1,1)
                lcorn(i+2,j+2,k+2) = 1
             end if
             ! Update counter
             cnt = cnt + 1

             if(ndime==2) exit
          end do
       end do
    end do

  end subroutine check_corn_boundary

  !================================================================================================
  subroutine neigh_face(newpo,pdime,auxv,ijkpart,ijkneig)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),newpo,pdime,auxv(3,2)
    integer(ip), intent(out) :: ijkneig(3)

    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:))
    ijkneig(pdime) = ijkpart(pdime) + newpo

  end subroutine neigh_face

  !================================================================================================
  subroutine neigh_edge(newpo,pdime,auxv,ijkpart,ijkneig)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),newpo(2),pdime,auxv(:,:)
    integer(ip), intent(out) :: ijkneig(3)

    ijkneig(auxv(pdime,:)) = ijkpart(auxv(pdime,:)) + newpo(1:size(auxv(pdime,:)))
    ijkneig(pdime) = ijkpart(pdime)

  end subroutine neigh_edge

  !================================================================================================
  subroutine neigh_corn(newpo,ijkpart,ijkneig)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ijkpart(3),newpo(3)
    integer(ip), intent(out) :: ijkneig(3)

    ijkneig = ijkpart + newpo

  end subroutine neigh_corn

  !==================================================================================================
  subroutine generic_l2g(subgl,npnumg,npnumt,nenum,ndime,pdegr,gsize,tsize,gdata,l2ge,l2gp,trian, &
       &                 coord,mater)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)             , intent(in)    :: subgl(ndime),ndime,pdegr
    integer(ip)             , intent(in)    :: npnumg(:),npnumt(:),nenum(:)
    type(geom_size)         , intent(in)    :: gsize
    type(topo_size)         , intent(in)    :: tsize
    type(geom_data)         , intent(in)    :: gdata
    integer(igp)            , intent(out)   :: l2ge(:),l2gp(:)
    type(fem_triangulation) , intent(inout) :: trian
    real(rp)                , intent(in)    :: coord(:,:)
    integer(ip)             , intent(inout) :: mater(:)

    integer(ip)               :: i,j,k,m,n,l,ijk(3),ijklnods(3),num,count,auxva((pdegr+1)**ndime)
    integer(ip)               :: ne_aux(3),np_aux(3),mcase,auxnum,nedir_l2g(ndime),nobje,aux_glb
    integer(ip)               :: subgl_aux(3),subgl_ijk(3),nelbl(3),ncorn
    integer(ip)               :: cnt,ipoin,jpoin,ielem,pdime,iedge,jedge,iface,nddomk
    integer(ip)               :: auxvc(2**ndime),auxvf(2*ndime),auxvd(2*(ndime-1)*ndime),aux_cnt
    integer(ip) , allocatable :: auxv(:,:)
    integer(igp), allocatable :: po_l2g(:),el_l2g(:)

    ! Allocate
    call memalloc(gsize%nedomt,el_l2g,__FILE__,__LINE__)
    !call memalloc(gsize%npdomt,po_l2g,__FILE__,__LINE__)
    call memalloc(tsize%notot ,po_l2g,__FILE__,__LINE__)

    ! Auxiliar vectro
    check( pdegr == 1 )  ! Only working for linear geometrical elements
    if(ndime==2) then
       call memalloc(2,1,auxv,__FILE__,__LINE__)
       auxv(1,1) = 2
       auxv(2,1) = 1
       if(pdegr==1) then
          auxva = (/1,3,2,4/)
       else if(pdegr==2) then
          auxva = (/1,8,4,5,9,7,2,6,3/)
       else if(pdegr==3) then
          auxva = (/1,12,11,4,5,13,16,10,6,14,15,9,2,7,8,3/)
       end if
       ne_aux(1)=1
       ne_aux(2)=gsize%nedom(1)
       ne_aux(3)=gsize%nedom(2)
       np_aux(1)=1
       np_aux(2)=gsize%npdom(1)
       np_aux(3)=gsize%npdom(2)
       subgl_aux(1:2) = subgl
       subgl_aux(3) = 0
       auxvc = (/1,3,2,4/)
       auxvd = (/7,8,5,6/)
       nobje = 8
       ncorn = 4
    else if(ndime==3) then
       call memalloc(3,2,auxv,__FILE__,__LINE__)
       auxv(1,:) = (/2,3/)
       auxv(2,:) = (/1,3/)
       auxv(3,:) = (/1,2/)
       if(pdegr==1) then
          auxva = (/1,5,3,7,2,6,4,8/)
       else if(pdegr==2) then
          auxva = (/1,13,5,12,25,20,4,16,8,9,22,17,21,27,26,11,24,19,2,14,6,10,23,18,3,15,7/)
       else if(pdegr==3) then
          auxva = (/1,17,21,5,16,44,52,32,15,43,51,31,4,20,24,8,9,37,45,25,33,57,61,53,36,60,64,56,14,42, &
               &    50,30,10,38,46,26,34,58,62,54,35,59,63,55,13,41,49,29,2,18,22,6,11,39,47,27,12,40,48, &
               &    28,3,19,23,7/)
       end if
       ne_aux=gsize%nedom
       np_aux=gsize%npdom
       subgl_aux = subgl
       auxvc = (/1,5,3,7,2,6,4,8/)
       auxvd = (/17,19,18,20,13,15,14,16,9,11,10,12/)
       auxvf = (/21,22,23,24,25,26/)
       nobje = 26
       ncorn = 8
    end if
    mcase = gdata%mater
    nelbl(1) = gdata%neblx
    nelbl(2) = gdata%nebly
    nelbl(3) = gdata%neblz

    ! Elements
    do i=1,ne_aux(1)
       do j=1,ne_aux(2)
          do k=1,ne_aux(3)

             ! Identifier
             if(ndime==2) then
                ijk=(/j,k,i/)
             else if(ndime==3) then
                ijk=(/i,j,k/)
             end if
             call globalid(ijk,gsize%nedom,ndime,num)
             
             ! Local to global vector
             subgl_ijk = subgl_aux+ijk
             call globalid(subgl_ijk,gsize%nedir,ndime,el_l2g(num))

             ! Fill material
             if(mcase>0) then
                call materialid(mcase,subgl_ijk,gsize%nedir,nelbl,mater(nenum(num)))
             end if

             ! Allocate triangulation elemental objects
             trian%elems(nenum(num))%num_objects = nobje
             call memalloc(trian%elems(nenum(num))%num_objects,trian%elems(nenum(num))%objects, &
                  &        __FILE__,__LINE__)
             call memalloc(trian%num_dims,ncorn,trian%elems(nenum(num))%coordinates, __FILE__, __LINE__ )
             call put_topology_element_triangulation(nenum(num),trian)
             trian%elems(nenum(num))%order = get_order(trian%elems(nenum(num))%topology%ftype,ncorn,trian%num_dims)
           
             count = 1
             ! Elemental corners
             do m=0,1
                do n=0,1
                   do l=0,1

                      ! Corner identifier in the element
                      ijklnods= ijk + (/m,n,l/)
                      call globalid(ijklnods,gsize%npdom,ndime,auxnum)

                      ! Generate lnods
                      trian%elems(nenum(num))%objects(auxvc(count)) = npnumt(auxnum)

                      ! Generate coords
                      trian%elems(nenum(num))%coordinates(:,auxva(count)) = coord(:,npnumg(auxnum))

                      ! Update counter
                      count = count +1

                      if(ndime==2) exit
                   end do
                end do
             end do
                
             ! Elemental edges
             count = 1
             aux_cnt = tsize%nctot
             do pdime=ndime,1,-1
                do iedge=0,1
                   do jedge=0,1

                      ! Edge identifier in the element
                      ijklnods(pdime) = ijk(pdime)
                      ijklnods(auxv(pdime,1)) = ijk(auxv(pdime,1)) + iedge
                      if(ndime==3) ijklnods(auxv(pdime,2)) = ijk(auxv(pdime,2)) + jedge
                      call globalid(ijklnods,tsize%nddom(:,pdime),ndime,auxnum)

                      ! Generate lnods
                      trian%elems(nenum(num))%objects(auxvd(count)) = npnumt(aux_cnt + auxnum)

                      ! Update counter
                      count = count +1

                      if(ndime==2) exit
                   end do
                end do
                aux_cnt = aux_cnt + tsize%nddir(pdime)
             end do

             ! Elemental faces
             count = 1
             if(ndime==3) then
                aux_cnt = tsize%nctot + tsize%ndtot
                do pdime=ndime,1,-1
                   do iface = 0,1

                      ! Face identifier in the element
                      ijklnods(pdime) = ijk(pdime) + iface
                      ijklnods(auxv(pdime,:)) = ijk(auxv(pdime,:))
                      call globalid(ijklnods,tsize%nfdom(:,pdime),ndime,auxnum)

                      ! Generate lnods
                      trian%elems(nenum(num))%objects(auxvf(count)) = npnumt(aux_cnt + auxnum)

                      ! Update counter
                      count = count + 1

                   end do
                   aux_cnt = aux_cnt + tsize%nfdir(pdime)
                end do
             end if

          end do
       end do
       if(ndime==2) exit
    end do

    ! Construct l2g emap (ordered by objects)
    do ielem=1,gsize%nedomt
       l2ge(nenum(ielem)) = el_l2g(ielem)
    end do

    ! Nodes
    do i=1,np_aux(1)
       do j=1,np_aux(2)
          do k=1,np_aux(3)

             ! Identifier
             if(ndime==2) then
                ijk=(/j,k,i/)
             else if(ndime==3) then
                ijk=(/i,j,k/)
             end if
             call globalid(ijk,gsize%npdom,ndime,num)

             ! Local to global vector
             nedir_l2g = gsize%nedir(1:ndime)*pdegr + 1
             subgl_ijk = subgl_aux*pdegr+ijk
             call globalid(subgl_ijk,nedir_l2g,ndime,po_l2g(num))
    
          end do
       end do
       if(ndime==2) exit
    end do

    ! Construct l2g nmap (ordered by objects)
    ! Elemental corners
    do ipoin=1,gsize%npdomt
       !l2gp(npnumg(ipoin)) = po_l2g(ipoin
       l2gp(npnumt(ipoin)) = po_l2g(ipoin)
    end do
    ! Elemental edges
    aux_cnt = tsize%nctot
    do pdime=ndime,1,-1
       do i=1,tsize%nddir(pdime)
          l2gp(npnumt(i+aux_cnt)) = po_l2g(i+aux_cnt)
       end do
       aux_cnt = aux_cnt + tsize%nddir(pdime)
    end do
    ! Elemental faces
    aux_cnt = tsize%nctot + tsize%ndtot
    do pdime=ndime,1,-1
       do i=1,tsize%nfdir(pdime)
          l2gp(npnumt(i+aux_cnt)) = po_l2g(i+aux_cnt)
       end do
       aux_cnt = aux_cnt + tsize%nfdir(pdime)
    end do

    ! Deallocate
    call memfree(el_l2g,__FILE__,__LINE__)
    call memfree(po_l2g,__FILE__,__LINE__)
    call memfree(auxv,__FILE__,__LINE__)

  end subroutine generic_l2g

  !==================================================================================================
  subroutine global_to_ijk(lpart,nsckt,npsoc,ndime,ijk)
    implicit none
    integer(ip), intent(in)  :: lpart,ndime,nsckt(ndime),npsoc(ndime)
    integer(ip), intent(out) :: ijk(3)

    integer(ip) :: isckt,jsckt,ksckt,lsckt
    integer(ip) :: ipart_aux,jpart_aux,kpart_aux,lpart_aux
    integer(ip) :: ipart,jpart,kpart
    real(rp) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7

    if(ndime==2) then
       ! First level (Socket level)
       aux1 = (lpart-1)/(npsoc(1)*npsoc(2))
       lsckt = floor(aux1) + 1
       aux2 = (lsckt-1)/nsckt(2)
       jsckt = lsckt - floor(aux2)*nsckt(2)
       isckt = floor(aux2) + 1

       ! Second level (Part inside the socket)
       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)
       aux3 = (lpart_aux-1)/npsoc(2)
       jpart_aux = lpart_aux - floor(aux3)*npsoc(2)
       ipart_aux = floor(aux3) + 1

       ! ijk numeration
       ipart = ipart_aux + (isckt-1)*npsoc(1)
       jpart = jpart_aux + (jsckt-1)*npsoc(2)
       kpart = 1
       ijk = (/ipart,jpart,kpart/)
    else if(ndime==3) then
       ! First level (Socket level)
       aux1 = (lpart-1)/(npsoc(1)*npsoc(2)*npsoc(3))
       lsckt = floor(aux1) + 1
       aux2 = (lsckt-1)/(nsckt(2)*nsckt(3))
       isckt = floor(aux2) + 1
       aux3 = (lsckt - (isckt-1)*nsckt(2)*nsckt(3) - 1)/nsckt(3)
       jsckt = floor(aux3) + 1
       aux4 = lsckt - (isckt-1)*nsckt(2)*nsckt(3) - (jsckt-1)*nsckt(3) - 1
       ksckt = floor(aux4) + 1

       ! Second level (Part inside the socket)
       lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)*npsoc(3)
       aux5 = (lpart_aux-1)/(npsoc(2)*npsoc(3))
       ipart_aux = floor(aux5) + 1
       aux6 = (lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - 1)/npsoc(3)
       jpart_aux = floor(aux6) + 1
       aux7 = lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - (jpart_aux-1)*npsoc(3) - 1
       kpart_aux = floor(aux7) + 1

       ! ijk numeration
       ipart = ipart_aux + (isckt-1)*npsoc(1)
       jpart = jpart_aux + (jsckt-1)*npsoc(2)
       kpart = kpart_aux + (ksckt-1)*npsoc(3)
       ijk = (/ipart,jpart,kpart/)       
    end if
  
  end subroutine global_to_ijk

  !================================================================================================
  subroutine aux_surf_num(pdime,aux_surf,gsurf_aux)
    implicit none
    integer(ip), intent(in)  :: pdime
    integer(ip), intent(out) :: aux_surf(4,2),gsurf_aux(4)

    if(pdime==3) then
       gsurf_aux = (/3,4,5,6/)
    else if(pdime==2) then
       gsurf_aux = (/1,2,5,6/)
    else if(pdime==1) then
       gsurf_aux = (/1,2,3,4/)
    end if
    aux_surf(1,:) = (/2,4/)
    aux_surf(2,:) = (/1,3/)
    aux_surf(3,:) = (/3,4/)
    aux_surf(4,:) = (/1,2/)    

  end subroutine aux_surf_num

  !==================================================================================================
  subroutine kronec(i,n,kron)
    !-----------------------------------------------------------------------
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: i,n
    integer(ip), intent(out) :: kron(n)

    kron = 0
    kron(i) = 1

  end subroutine kronec

  !==================================================================================================
  subroutine globalid_ip(ijk,nd,ndime,gl)
    implicit none
    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
    integer(ip), intent(out) :: gl

    if(ndime==1) then
       gl = ijk(1)
    elseif(ndime==2) then
       gl = (ijk(1)-1)*nd(2)+ijk(2)
    else if(ndime==3) then
       gl = (ijk(1)-1)*nd(2)*nd(3)+(ijk(2)-1)*nd(3)+ijk(3)
    end if

  end subroutine globalid_ip

  !==================================================================================================
  subroutine globalid_igp(ijk,nd,ndime,gl)
    implicit none
    integer(ip), intent(in)  :: ijk(ndime),nd(ndime),ndime
    integer(igp), intent(out) :: gl
    
    if(ndime==1) then
       gl = int(ijk(1),igp)
    elseif(ndime==2) then
       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)+int(ijk(2),igp)
    else if(ndime==3) then
       gl = (int(ijk(1),igp)-1)*int(nd(2),igp)*int(nd(3),igp)+(int(ijk(2),igp)-1)*int(nd(3),igp)+int(ijk(3),igp)
    end if
    
  end subroutine globalid_igp
  
  !==================================================================================================
  subroutine globalid_2l(ijk,nsckt,npsoc,ndime,gl_2l)
    implicit none
    integer(ip), intent(in)  :: ijk(3),nsckt(3),npsoc(3),ndime
    integer(ip), intent(out) :: gl_2l

    integer(ip) :: aux1(ndime),aux2(ndime),aux3,aux4,i,gl,local_ijk(ndime)
    real(rp)    :: work

    ! Auxiliar variables
    aux3=1
    do i=1,ndime
       work = (ijk(i)-1)/npsoc(i)
       aux1(i) = floor(work) + 1;
       aux2(i) = floor(work)*npsoc(i)
       aux3 = aux3*npsoc(i)
    end do

    ! 2-level global numeration
    call globalid(aux1,nsckt,ndime,aux4)
    aux4 = (aux4-1)*aux3
    local_ijk = ijk(1:ndime) - aux2
    call globalid(local_ijk,npsoc,ndime,gl)
    gl_2l = aux4 + gl
    
  end subroutine globalid_2l

  !==================================================================================================
  subroutine materialid(mcase,ijkelem,nedir,nelbl,mater)
    implicit none
    integer(ip), intent(in)    :: mcase,ijkelem(3),nedir(3),nelbl(3)
    integer(ip), intent(inout) :: mater
    
    if(mcase==1) then
       ! (x<leng(x)/2 --> mat == 1)
       ! (x>leng(x)/2 --> mat == 2)
       assert(mod(nedir(1),2)==0)
       if(ijkelem(1).le.nedir(1)/2) then
          mater = 1
       else
          mater = 2
       end if
    elseif(mcase==2) then
       ! (y<leng(y)/2 --> mat == 1)
       ! (y>leng(y)/2 --> mat == 2)
       assert(mod(nedir(2),2)==0)
       if(ijkelem(2).le.nedir(2)/2) then
          mater = 1
       else
          mater = 2
       end if
    elseif(mcase==3) then   ! Hunt's case
       ! (z >= -1.0 and z <= 1.0 --> mat == 1) IMH
       ! (z < -1.0 or z > 1.0 --> mat == 2) DCY
       if(ijkelem(3)>nelbl(3) .and. ijkelem(3)<=nedir(3)-nelbl(3)) then  ! We have nelbl elements in the layer (DCY-solid)
          mater = 1
       else
          mater = 2
       end if
    elseif(mcase==4) then   ! Backward-facing step
       ! (x >= 0.0 or y >= 0.0 --> mat == 1) NSI
       ! (x < 0.0 and y < 0.0 --> mat == 2) Boundary
       if(ijkelem(1)>nelbl(1) .or. ijkelem(2)>nelbl(2)) then 
          mater = 1
       else
          mater = 2
       end if
    end if
    
  end subroutine materialid
                
end module fem_mesh_gen_distribution_names
