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
module mesh_names
  use types_names
  use memor_names
  use list_types_names
  use hash_table_names
  !use materials_names
  use stdio_names
  use reference_fe_names
  use reference_fe_factory_names
  !use conditions_names
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: c_order = 0
  integer(ip), parameter :: z_order = 1

  integer(ip), target :: permu_2DP1(3) = (/ 1, 2, 3/)
  integer(ip), target :: permu_2DQ1(4) = (/ 1, 2, 4, 3/)
  integer(ip), target :: permu_3DP1(4) = (/ 1, 2, 3, 4/)
  integer(ip), target :: permu_3DPR(6) = (/ 1, 2, 3, 4, 5, 6/)
  integer(ip), target :: permu_3DQ1(8) = (/ 1, 2, 4, 3, 5, 6, 8, 7/)
  integer(ip), target :: permu_id  (8) = (/ 1, 2, 3, 4, 5, 6, 7, 8/)

  type mesh_t
     ! Sizes
     integer(ip)                :: &
          order=c_order,           &         ! GiD element order (c)
          nelty=1,                 &         ! Number of element types
          ndime,                   &         ! Number of space dimensions
          npoin,                   &         ! Number of nodes (vertices)
          nvefs,                   &         ! Number of vefs
          nelem,                   &         ! Number of elements
          nnode,                   &         ! Maximum number of nodes per element
          nboun,                   &         ! Number of boundary elements
          nnodb                              ! Maximum number of nodes per boundary element

     ! Elements
     integer(ip), allocatable ::  &
          pnods(:),               &         ! pointers to the lnods
          lnods(:),               &         ! list of vefs of each element
          legeo(:),               &         ! List of geometry (volume) each element lies in
          leset(:),               &         ! List of sets associated to each element
          pvefs(:),               &         ! pointers to the lvefs
          lvefs(:),               &         ! list of vefs of each element
          lvef_geo(:),            &         ! List of geometric entities (volume, surface, point) each vef lies in
          lvef_set(:)                       ! List of sets associated to each vef

     ! Booundary (understood as a subset of vefs on which the definition of a set and a geometry is relevant,
     ! e.g. an internal interface separating materials over which a force has to be computed)
     integer(ip), allocatable ::  &
          pnodb(:),               &         ! pointers to the lnodb
          lnodb(:),               &         ! list of vertices of each boundary element (edges and faces)
          lbgeo(:),               &         ! List of geometric entities (volume, surface, point) each boundary lies in
          lbset(:)                          ! List of sets associated to each boundary

     ! Dual mesh (elements around vertices)
     integer(ip)              ::  &
          nelpo = 0                         ! Nonzero when created
     integer(ip), allocatable ::  &
          pelpo(:),               &
          lelpo(:)

     real(rp), allocatable ::     &
          coord(:,:)                         ! Node coordinates
    contains
     ! JP-TODO: program get and set for variables.
     procedure, non_overridable :: to_dual              => mesh_to_dual_new
     procedure, non_overridable :: generate_vefs        => mesh_generate_vefs
     procedure, non_overridable :: free                 => mesh_free
     procedure, non_overridable :: read                 => mesh_read
     procedure, non_overridable :: read_file            => mesh_read_file
     procedure, non_overridable, nopass :: compose_name => mesh_compose_name
     procedure, non_overridable :: write_file_for_postprocess => mesh_write_file_for_postprocess
  end type mesh_t

  ! Types
  public :: mesh_t

  ! Constants
  public :: c_order, z_order

  ! Functions
  public :: mesh_write_file, mesh_write_files, mesh_write_files_for_postprocess

contains
!=============================================================================
  subroutine mesh_to_dual_new(this)
    class(mesh_t), intent(inout)     :: this
    
    ! Local variables
    integer(ip)              :: inode, ipoin, ielem, size_lnods

    if(this%nelpo>0) return

    call memalloc (this%npoin+1, this%pelpo, __FILE__,__LINE__)
    size_lnods = this%pnods(this%nelem+1)-1 
    
    ! Compute the number of elements around each point
    this%pelpo=0
    do inode=1, size_lnods
       ipoin=this%lnods(inode)
       this%pelpo(ipoin+1)=this%pelpo(ipoin+1)+1
    end do
    
    ! Find the maximum number of elements around a point
    this%nelpo=0
    do ipoin=1,this%npoin
       this%nelpo=max(this%nelpo,this%pelpo(ipoin+1))
    end do
    
    ! Compute pointers to the starting position of the list
    ! of elements around each point
    this%pelpo(1)=1
    do ipoin=1,this%npoin
       this%pelpo(ipoin+1)=this%pelpo(ipoin+1)+this%pelpo(ipoin)
    end do

    ! Allocate lelpo and fill it
    call memalloc (this%pelpo(this%npoin+1), this%lelpo, __FILE__,__LINE__)

    ! Compute the list of elements around each point.
    ! pelpo is used instead of auxiliary work space.
    do ielem=1,this%nelem 
       do inode=this%pnods(ielem),this%pnods(ielem+1)-1 
          ipoin=this%lnods(inode)
          this%lelpo(this%pelpo(ipoin))=ielem
          this%pelpo(ipoin)=this%pelpo(ipoin)+1
       end do
    end do
    
    ! Recover pelpo
    do ipoin=this%npoin+1, 2, -1
       this%pelpo(ipoin)=this%pelpo(ipoin-1)
    end do
    this%pelpo(1) = 1

  end subroutine mesh_to_dual_new

  !=============================================================================
  subroutine mesh_free (msh)
    !-----------------------------------------------------------------------
    ! This routine generates deallocates a mesh
    !-----------------------------------------------------------------------
    implicit none
    class(mesh_t), intent(inout)  :: msh

    call memfree (msh%pnods,__FILE__,__LINE__)
    call memfree (msh%lnods,__FILE__,__LINE__)

    ! This data structure can hold (depending on the context):
    ! (1) Permanent geometrical mesh (e.g., mesh read from GiD)
    ! (2) Temporary topological mesh (resulting from generate_vefs_mesh_conditions)
    ! (2) Temporary Dual mesh
    ! (3) Coarse-grid mesh in MLBDDC
    ! In the case of (2), (3), and (4) coord is not allocated!!!
    if (allocated(msh%coord)) call memfree (msh%coord,__FILE__,__LINE__)

    msh%ndime=0
    msh%npoin=0
    msh%nelem=0
    msh%nnode=0

  end subroutine mesh_free

  !===============================================================================================
  subroutine mesh_copy(msh_old,msh_new)
    implicit none
    type(mesh_t), intent(in)    :: msh_old
    type(mesh_t), intent(inout) :: msh_new

    msh_new%ndime = msh_old%ndime
    msh_new%npoin = msh_old%npoin
    msh_new%nelem = msh_old%nelem
    msh_new%nnode = msh_old%nnode
    
    call memalloc(msh_new%nelem+1,msh_new%pnods,__FILE__,__LINE__)
    msh_new%pnods = msh_old%pnods

    call memalloc(msh_new%pnods(msh_new%nelem+1),msh_new%lnods,__FILE__,__LINE__)
    msh_new%lnods = msh_old%lnods

    if (allocated(msh_old%coord)) then
       call memalloc(msh_new%ndime,msh_new%npoin,msh_new%coord,__FILE__,__LINE__)
       msh_new%coord = msh_old%coord
    end if

  end subroutine mesh_copy
  
  !=============================================================================
  subroutine mesh_read_file(msh,lunio,permute_c2z)
    !------------------------------------------------------------------------
    !
    ! This routine reads a mesh writen by GiD according to fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: lunio
    class(mesh_t)     , intent(out) :: msh
    logical, optional, intent(in)  :: permute_c2z
    integer(ip)     :: idime,ipoin,inode,ielem,iboun,nnode,nnodb ! ,istat,jpoin,jelem,iboun
    character(14)   :: dum1
    character(7)    :: dum2
    character(10)   :: dum3
    character(7)    :: dum4
    character(12)   :: dum5
    character(1000) :: tel
    integer(ip), allocatable :: lnods_aux(:)
    integer(ip), pointer     :: permu(:)
    logical                  :: permute_c2z_

    if(present(permute_c2z)) then
       permute_c2z_ = permute_c2z
    else
       permute_c2z_ = .true.
    end if

    ! Read first line: "MESH dimension  3  types  1  elements        352  nodes        100  boundaries        144"
    read(lunio,'(a14,1x,i2,a7,1x,i2,a10,1x,i10,a7,1x,i10,a12,1x,i10)') &
         & dum1,msh%ndime,dum2,msh%nelty,dum3,msh%nelem,dum4,msh%npoin,dum5,msh%nboun

    write(*,*) msh%ndime,msh%order,msh%nelty,msh%nelem,msh%npoin,msh%nboun

    ! Read nodes
    call memalloc(msh%ndime,msh%npoin,msh%coord,__FILE__,__LINE__)
    do while(tel(1:5).ne.'coord')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end c')
       read(tel,*) ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
       read(lunio,'(a)') tel
    end do

    ! Read elements' size (pnods)
    call memalloc(msh%nelem+1,msh%pnods,__FILE__,__LINE__)
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       read(tel,*) ielem,msh%pnods(ielem+1)
       read(lunio,'(a)') tel
    end do
    ! Transform length to header and get mesh%nnode
    msh%pnods(1) = 1
    msh%nnode    = 0
    do ielem = 2, msh%nelem+1
       msh%nnode = max(msh%nnode,msh%pnods(ielem))
       msh%pnods(ielem) = msh%pnods(ielem)+msh%pnods(ielem-1)
    end do

    ! Read elements
    call memalloc(msh%pnods(msh%nelem+1),msh%lnods,__FILE__,__LINE__)
    call memalloc(msh%nelem,msh%legeo,__FILE__,__LINE__)
    call memalloc(msh%nelem,msh%leset,__FILE__,__LINE__)
    call io_rewind(lunio)
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       read(tel,*) ielem,nnode,(msh%lnods(msh%pnods(ielem)-1+inode),inode=1,nnode),msh%leset(ielem),msh%legeo(ielem)
       read(lunio,'(a)') tel
    end do

    ! Read boundary elements' size (pnodb)
    call memalloc(msh%nboun+1,msh%pnodb,__FILE__,__LINE__)
    do while(tel(1:5).ne.'bound')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end b')
       read(tel,*) iboun,msh%pnodb(iboun+1)
       read(lunio,'(a)') tel
    end do
    ! Transform length to header and get mesh%nnodb
    msh%pnodb(1) = 1
    msh%nnodb    = 0
    do iboun = 2, msh%nboun+1
       msh%nnodb = max(msh%nnodb,msh%pnodb(iboun))
       msh%pnodb(iboun) = msh%pnodb(iboun)+msh%pnodb(iboun-1)
    end do

    ! Read boundary elements
    call memalloc(msh%pnodb(msh%nboun+1),msh%lnodb,__FILE__,__LINE__)
    call memalloc(msh%nboun,msh%lbgeo,__FILE__,__LINE__)
    call memalloc(msh%nboun,msh%lbset,__FILE__,__LINE__)
    call io_rewind(lunio)
    do while(tel(1:5).ne.'bound')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end b')
       read(tel,*) iboun,nnodb,(msh%lnodb(msh%pnodb(iboun)-1+inode),inode=1,nnodb),msh%lbset(iboun),msh%lbgeo(iboun)
       read(lunio,'(a)') tel
    end do

    ! Reordering (c to z) the nodes of the mesh, if needed
    !if(permute_c2z_) then
    if(msh%order==c_order) then
       msh%order= z_order
       call memalloc(msh%nnode, lnods_aux, __FILE__, __LINE__)
       do ielem = 1,msh%nelem
          nnode = msh%pnods(ielem+1) - msh%pnods(ielem)
          lnods_aux(1:nnode) = msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1)
          if(msh%ndime == 2) then
             if(nnode == 3)  then    ! Linear triangles (2DP1)
                permu => permu_2DP1
             elseif(nnode == 4) then ! Linear quadrilaterals(2DQ1)
                permu => permu_2DQ1
             end if
          elseif(msh%ndime == 3) then
             if(nnode == 4) then     ! Linear tetrahedra (3DP1)
                permu => permu_3DP1
             elseif(nnode == 8) then ! Linear hexahedra (3DQ1)
                permu => permu_3DQ1
             end if
          end if
          do inode = 1, nnode
             msh%lnods(msh%pnods(ielem)+inode-1) = lnods_aux(permu(inode))
          end do
       end do
       call memfree(lnods_aux,__FILE__,__LINE__)
    end if

  end subroutine mesh_read_file

  !=============================================================================
  subroutine mesh_write_file (msh,lunio,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in the format defined by GiD fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)  , intent(in)           :: lunio
    class(mesh_t) , intent(in)           :: msh
    character(*) , intent(in), optional :: title

    integer(ip)                    :: ielem, idime, ipoin, inode, iboun

    ! Read first line: "MESH dimension  2  order  0  types  2  elements         86  nodes         63  boundaries         28"
    write(lunio,'(a14,1x,i2,a7,1x,i2,a7,1x,i2,a10,1x,i10,a7,1x,i10,a12,1x,i10)') &
         & 'MESH dimension',msh%ndime,'  order',msh%order,'  types',msh%nelty,'  elements', &
         & msh%nelem,'  nodes',msh%npoin,'  boundaries',msh%nboun

    ! Coordinates
    write(lunio,'(a)')'coordinates'
    assert(allocated(msh%coord))
    do ipoin=1,msh%npoin
       write(lunio,'(i10,3(1x,e16.8e3))') ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
    end do
    write(lunio,'(a)')'end coordinates'

    ! Elements
    write(lunio,'(a)')'elements'
    do ielem=1,msh%nelem
       write(lunio,'(i10,65(1x,i10))') ielem, msh%pnods(ielem+1)-msh%pnods(ielem),&
            &  msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1),msh%leset(ielem),msh%legeo(ielem)
    end do
    write(lunio,'(a)') 'end elements'

    ! Boundary elements
    write(lunio,'(a)')'boundaries'
    do iboun=1,msh%nboun
       write(lunio,'(i10,65(1x,i10))') iboun, msh%pnodb(iboun+1)-msh%pnodb(iboun), &
            &  msh%lnodb(msh%pnodb(iboun):msh%pnodb(iboun+1)-1),msh%lbset(iboun),msh%lbgeo(iboun)
    end do
    write(lunio,'(a)') 'end boundaries'

  end subroutine mesh_write_file

  !=============================================================================
  subroutine mesh_write_post_file (msh,lunio,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in GiD format (only works for linear elements).
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)           :: lunio
    class(mesh_t)     , intent(in)           :: msh
    character(*)     , intent(in), optional :: title

    integer(ip)                    :: ielem, idime, ipoin, inode, nnode
    character(13)                  :: elemt
    character(len=:), allocatable  :: title_

    integer(ip)     , pointer      :: permu(:)

    permu => permu_id
    if(msh%ndime==2) then
       if(msh%nnode==3) then
          elemt='Triangle' 
          if(msh%order==z_order) permu => permu_2DP1
       else
          elemt='Quadrilateral'
          if(msh%order==z_order) permu => permu_2DQ1
       end if
    else
       if(msh%nnode==4) then 
          elemt='Tetrahedra'
          if(msh%order==z_order) permu => permu_3DP1
       else if(msh%nnode==6) then 
          elemt='Prism'
          if(msh%order==z_order) permu => permu_3DPR
       else if(msh%nnode==8) then 
          elemt='Hexahedra'
          if(msh%order==z_order) permu => permu_3DQ1
       end if
    end if

    ! Header
    title_ = 'TITLE'
    if(present(title)) title_=title
    write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),msh%nnode

    ! Coordinates
    write(lunio,2)'coordinates'
    if (allocated(msh%coord)) then
       do ipoin=1,msh%npoin
          write(lunio,3) ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
       end do
    end if
    write(lunio,2)'end coordinates'

    ! Connectivity
    if(msh%nelty==1) then
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          nnode = msh%pnods(ielem+1)-msh%pnods(ielem)
          write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
    else
       ! Write hexahedra or prismas (3D) or quads(2)
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          nnode = msh%pnods(ielem+1)-msh%pnods(ielem)
          if(nnode == msh%nnode) &
               write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
       ! Now write tetrahedra (3D) or triangles (2D)
       if(msh%ndime==2) then
          nnode = 3
          elemt = 'Triangle' 
          if(msh%order==z_order) permu => permu_2DP1
       else if(msh%ndime==3) then
          nnode = 4
          elemt='Tetrahedra'
          if(msh%order==z_order) permu => permu_3DP1
       end if
       write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),nnode
       write(lunio,2)'coordinates'
       write(lunio,2)'end coordinates'
       write(lunio,2)'elements'
       do ielem=1,msh%nelem
          if(msh%pnods(ielem+1)-msh%pnods(ielem) == nnode) &
               write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
       end do
       write(lunio,2) 'end elements'
       ! Eventually write prismas (3D)
       if(msh%ndime==3.and.msh%nnode==8) then
          nnode = 4
          elemt='Prism'
          if(msh%order==z_order) permu => permu_3DPR
          write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),nnode
          write(lunio,2)'coordinates'
          write(lunio,2)'end coordinates'
          write(lunio,2)'elements'
          do ielem=1,msh%nelem
             if(msh%pnods(ielem+1)-msh%pnods(ielem) == nnode) &
                  write(lunio,4) ielem, (msh%lnods(msh%pnods(ielem)-1+permu(inode)),inode=1,nnode),1
          end do
          write(lunio,2) 'end elements'
       end if
    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i10, 3(1x,e16.8e3))
4   format(i10,65(1x,i10))
5   format('BOUNDARY ',a,' Nnodb ',i2)
6   format(i6,10(1x,i6))

  end subroutine mesh_write_post_file

  !=============================================================================
  subroutine mesh_compose_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.mesh'
  end subroutine mesh_compose_name
  !=============================================================================
  subroutine mesh_compose_post_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.post.msh'
  end subroutine mesh_compose_post_name

  !=============================================================================
  subroutine mesh_write_files ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*)   , intent(in)  :: dir_path 
     character(*)   , intent(in)  :: prefix
     integer(ip)     , intent(in)  :: nparts
     type(mesh_t)  , intent(in)  :: lmesh (nparts)

     character(len=:), allocatable :: name, rename ! Deferred-length allocatable character arrays

     ! Locals 
     integer (ip) :: i,lunio

     call mesh_compose_name ( prefix, name )

     do i=nparts, 1, -1  
        rename=name
        call numbered_filename_compose(i,nparts,rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
        call mesh_write_file(lmesh(i),lunio)
        call io_close(lunio)
     end do
     
     ! name, and rename should be automatically deallocated by the compiler when they
     ! go out of scope. Should we deallocate them explicitly for safety reasons?
   end subroutine mesh_write_files

  !=============================================================================
  subroutine mesh_write_files_for_postprocess ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*)    , intent(in)  :: dir_path 
     character(*)    , intent(in)  :: prefix
     integer(ip)     , intent(in)  :: nparts
     type(mesh_t)    , intent(in)  :: lmesh (nparts)
     character(len=:), allocatable :: name, rename ! Deferred-length allocatable character arrays

     ! Locals 
     integer (ip) :: i,lunio

     do i=nparts, 1, -1  
        name=prefix
        call numbered_filename_compose(i,nparts,name)
        call mesh_compose_post_name (name, rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
        call mesh_write_post_file(lmesh(i),lunio)
        call io_close(lunio)
     end do
     
   end subroutine mesh_write_files_for_postprocess

  !=============================================================================
   subroutine mesh_read_files ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*), intent(in)       :: dir_path 
     character(*), intent(in)       :: prefix
     integer(ip)   , intent(in)     :: nparts
     type(mesh_t), intent(out)    :: lmesh (nparts)
     character(len=:), allocatable  :: name, rename ! Deferred-length allocatable character arrays
     
     ! Locals 
     integer (ip)                     :: i,lunio
     
     call mesh_compose_name ( prefix, name )
     do i=nparts, 1, -1  
        rename=name
        call numbered_filename_compose(i,nparts,rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'read' )
        call lmesh(i)%read_file(lunio)
        call io_close(lunio)
     end do
     
   end subroutine mesh_read_files

  !=============================================================================
   subroutine mesh_read (f_mesh,  dir_path, prefix, permute_c2z )
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     class(mesh_t)                , intent(out) :: f_mesh
     logical, optional, intent(in)  :: permute_c2z
     ! Locals
     integer                        :: iam, num_procs
     integer(ip)                    :: j, ndigs_iam, ndigs_num_procs, lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call mesh_compose_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call f_mesh%read_file(lunio, permute_c2z)
     call io_close(lunio)
   end subroutine mesh_read

   !=============================================================================
   subroutine mesh_write_file_for_postprocess ( f_mesh, dir_path, prefix)
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     class(mesh_t)                , intent(in)  :: f_mesh

     ! Locals
     integer(ip)                    :: lunio
     character(len=:), allocatable  :: name

     call mesh_compose_post_name ( prefix, name )
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' )
     call mesh_write_post_file(f_mesh,lunio)
     call io_close(lunio)

   end subroutine mesh_write_file_for_postprocess

  !==================================================================================================
  subroutine mesh_generate_vefs(mesh)
    implicit none
    ! Parameters
    class(mesh_t), intent(inout) :: mesh

    type(list_t), pointer    :: vertices_ivef
    type(list_t), pointer    :: vertices_jvef
    logical     :: equal
    integer(ip) :: istat, count, iboun, ivert, jvert, nnodb
    integer(ip) :: ielem, ivef, ielem_type, ielem_num_nodes
    integer(ip) :: ielem_num_vefs, ielem_first_vef_id, ielem_num_vef_verts
    integer(ip) :: vertex_of_ivef(4)
    integer(ip) :: jelpo, jelem, jvef, jelem_type, jelem_num_nodes
    integer(ip) :: jelem_num_vefs, jelem_first_vef_id, jelem_num_vef_verts
    integer(ip) :: vertex_of_jvef(4)

    ! We only require linear continuous elements so we could only have
    ! tetrahedra, hexahedra and prisms combined.
    integer(ip), parameter      :: max_num_elem_types = 3
    type(position_hash_table_t) :: pos_reference_fe
    type(p_reference_fe_t)      :: reference_fe(max_num_elem_types)

    ! Build reference_fe list
    call pos_reference_fe%init(max_num_elem_types)
    do ielem=1,mesh%nelem
       ielem_num_nodes=mesh%pnods(ielem+1)-mesh%pnods(ielem)
       call pos_reference_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       if ( istat == new_index ) then
          if(mesh%ndime==2) then
             if(ielem_num_nodes==3) then ! Triangle
                assert(.false.)
             else ! Quadrilateral
                reference_fe(ielem_type) = &
                     &    make_reference_fe ( topology = topology_quad, fe_type = fe_type_lagrangian, &
                     &                        number_dimensions = mesh%ndime, order = 1,          &
                     &                        field_type = field_type_vector, continuity = .true. )
             end if
          else
             if(ielem_num_nodes==4) then ! Tetrahedra
                assert(.false.)
             else if(ielem_num_nodes==6) then ! Prism
                assert(.false.)
             else if(ielem_num_nodes==8) then ! Hexahedra
                assert(.false.)
             end if
          end if
       end if
    end do

    ! Compute mesh%pvefs and allocate mesh%lvefs
    call memalloc(mesh%nelem+1, mesh%pvefs, __FILE__, __LINE__ )
    mesh%pvefs(1)=1
    do ielem=1,mesh%nelem
       ielem_num_nodes=mesh%pnods(ielem+1)-mesh%pnods(ielem)
       call pos_reference_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       assert(istat==old_index)
       ielem_num_vefs = ielem_num_nodes &
            &   + reference_fe(ielem_type)%p%get_number_vefs_of_dimension(1) 
       if(mesh%ndime==3) ielem_num_vefs = ielem_num_vefs &
            &   + reference_fe(ielem_type)%p%get_number_vefs_of_dimension(2)
       mesh%pvefs(ielem+1)=mesh%pvefs(ielem)+ielem_num_vefs
    end do
    call memalloc(mesh%pvefs(mesh%nelem+1), mesh%lvefs , __FILE__, __LINE__ )
    mesh%lvefs=0
    !write(*,*) mesh%pvefs

    ! Fill vefs
    mesh%nvefs = mesh%npoin
    do ielem=1,mesh%nelem
       ielem_num_nodes = mesh%pnods(ielem+1)-mesh%pnods(ielem)
       ! Fill vertices
       mesh%lvefs(mesh%pvefs(ielem):mesh%pvefs(ielem)+ielem_num_nodes-1)=mesh%lnods(mesh%pnods(ielem):mesh%pnods(ielem+1)-1)
       call pos_reference_fe%get(key=ielem_num_nodes,val=ielem_type,stat=istat)
       assert(istat==old_index)
       ! Fill edges
       ielem_num_vefs     = reference_fe(ielem_type)%p%get_number_vefs_of_dimension(1)
       ielem_first_vef_id = reference_fe(ielem_type)%p%get_first_vef_id_of_dimension(1)
       vertices_ivef => reference_fe(ielem_type)%p%get_vertices_vef()
       do ivef=1,ielem_num_vefs
          if(mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)==0) then ! Not filled yet
             mesh%nvefs=mesh%nvefs+1                                       ! Count it
             mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)=mesh%nvefs ! Fill it
             vertex_of_ivef(1) = mesh%lnods(mesh%pnods(ielem)-1+vertices_ivef%l(vertices_ivef%p(ielem_first_vef_id+ivef-1)))
             vertex_of_ivef(2) = mesh%lnods(mesh%pnods(ielem)-1+vertices_ivef%l(vertices_ivef%p(ielem_first_vef_id+ivef-1)+1))
             do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
                jelem=mesh%lelpo(jelpo)
                if(jelem>ielem) then
                   jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
                   call pos_reference_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
                   assert(istat==old_index)
                   jelem_num_vefs     = reference_fe(jelem_type)%p%get_number_vefs_of_dimension(1)
                   jelem_first_vef_id = reference_fe(jelem_type)%p%get_first_vef_id_of_dimension(1)
                   vertices_jvef => reference_fe(jelem_type)%p%get_vertices_vef()
                   do jvef=1,jelem_num_vefs
                      vertex_of_jvef(1) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+ivef-1)))
                      vertex_of_jvef(2) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+ivef-1)+1))
                      ! Compare, here we are using that edges have two vertices, hard coded
                      equal = (vertex_of_ivef(1)==vertex_of_jvef(1).and.vertex_of_ivef(2)==vertex_of_jvef(2)).or. &
                           &  (vertex_of_ivef(1)==vertex_of_jvef(2).and.vertex_of_ivef(2)==vertex_of_jvef(1))
                      if(equal) then ! Fill it
                         mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)=mesh%nvefs
                         exit
                      end if
                   end do
                end if
             end do
          end if
       end do
       ! Fill faces (similar code except for the number of vertices of each face that is variable)
       if(mesh%ndime==3) then
          ielem_num_vefs      = reference_fe(ielem_type)%p%get_number_vefs_of_dimension(2)
          ielem_first_vef_id  = reference_fe(ielem_type)%p%get_first_vef_id_of_dimension(2)
          do ivef=1,ielem_num_vefs
             if(mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)==0) then ! Not filled yet
                mesh%nvefs=mesh%nvefs+1                                       ! Count it
                mesh%lvefs(mesh%pvefs(ielem)-1+ielem_first_vef_id-1+ivef)=mesh%nvefs ! Fill it
                ielem_num_vef_verts = reference_fe(ielem_type)%p%get_number_vertices_vef(ielem_first_vef_id+ivef-1)
                vertices_ivef => reference_fe(ielem_type)%p%get_vertices_vef()
                vertex_of_ivef = 0
                do ivert=1,ielem_num_vef_verts
                   vertex_of_ivef(ivert)=mesh%lnods( mesh%pnods(ielem)-1+vertices_ivef%l(vertices_ivef%p(ielem_first_vef_id+ivef-1)+ivert-1))
                end do
                do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
                   jelem=mesh%lelpo(jelpo)
                   if(jelem>ielem) then
                      jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
                      call pos_reference_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
                      assert(istat==old_index)
                      jelem_num_vefs      = reference_fe(jelem_type)%p%get_number_vefs_of_dimension(2)
                      jelem_first_vef_id  = reference_fe(jelem_type)%p%get_first_vef_id_of_dimension(2)
                      vertices_jvef => reference_fe(jelem_type)%p%get_vertices_vef()
                      do jvef=1,jelem_num_vefs
                         jelem_num_vef_verts = reference_fe(jelem_type)%p%get_number_vertices_vef(jelem_first_vef_id+jvef-1)
                         if(jelem_num_vef_verts==ielem_num_vef_verts) then
                            vertex_of_jvef = 0
                            do jvert=1,jelem_num_vef_verts
                               vertex_of_jvef(jvert)=mesh%lnods( mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+jvef-1)+jvert-1))
                            end do
                            count=0
                            do ivert=1,ielem_num_vef_verts
                               do jvert=1,jelem_num_vef_verts
                                  if(vertex_of_ivef(ivert)==vertex_of_jvef(jvert)) then
                                     count=count+1
                                     exit
                                  end if
                               end do
                            end do
                            equal=(count==ielem_num_vef_verts)
                            if(equal) then ! Fill it
                               mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)=mesh%nvefs
                               exit
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    ! Identify boundary faces and assign set and geometry to vefs
    call memalloc(mesh%nvefs, mesh%lvef_geo, __FILE__, __LINE__ )
    call memalloc(mesh%nvefs, mesh%lvef_set, __FILE__, __LINE__ )
    do iboun=1,mesh%nboun
       nnodb=mesh%pnodb(iboun+1)-mesh%pnodb(iboun)
       if(nnodb==1) then      ! Vertex
          ivert=mesh%lnodb(mesh%pnodb(iboun))
          mesh%lvef_geo(ivert)=mesh%lbgeo(iboun)
          mesh%lvef_set(ivert)=mesh%lbset(iboun)
       else if(nnodb==2) then ! Edge
          vertex_of_ivef(1) = mesh%lnodb(mesh%pnodb(iboun)) 
          vertex_of_ivef(2) = mesh%lnodb(mesh%pnodb(iboun)+1)
          elems1: do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
             jelem=mesh%lelpo(jelpo)
             jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
             call pos_reference_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
             assert(istat==old_index)
             jelem_num_vefs     = reference_fe(jelem_type)%p%get_number_vefs_of_dimension(1)
             jelem_first_vef_id = reference_fe(jelem_type)%p%get_first_vef_id_of_dimension(1)
             vertices_jvef => reference_fe(jelem_type)%p%get_vertices_vef()
             do jvef=1,jelem_num_vefs
                vertex_of_jvef(1) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+jvef-1)))
                vertex_of_jvef(2) = mesh%lnods(mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+jvef-1)+1))
                ! Compare, here we are using that edges have two vertices, hard coded
                equal = (vertex_of_ivef(1)==vertex_of_jvef(1).and.vertex_of_ivef(2)==vertex_of_jvef(2)).or. &
                     &  (vertex_of_ivef(1)==vertex_of_jvef(2).and.vertex_of_ivef(2)==vertex_of_jvef(1))
                if(equal) then ! Fill it
                   mesh%lvef_geo( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbgeo(iboun)
                   mesh%lvef_set( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbset(iboun)
                   exit elems1
                end if
             end do
          end do elems1
       else                   ! Face
          ivert=mesh%lnodb(mesh%pnodb(iboun))
          vertex_of_ivef = 0
          do ivert=1,nnodb
             vertex_of_ivef(ivert)= mesh%lnodb(mesh%pnodb(iboun)-1+ivert)
          end do
          elems2: do jelpo=mesh%pelpo(vertex_of_ivef(1)),mesh%pelpo(vertex_of_ivef(1)+1)-1
             jelem=mesh%lelpo(jelpo)
             jelem_num_nodes=mesh%pnods(jelem+1)-mesh%pnods(jelem)
             call pos_reference_fe%get(key=jelem_num_nodes,val=jelem_type,stat=istat)
             assert(istat==old_index)
             jelem_num_vefs     = reference_fe(jelem_type)%p%get_number_vefs_of_dimension(2)
             jelem_first_vef_id = reference_fe(jelem_type)%p%get_first_vef_id_of_dimension(2)
             vertices_jvef => reference_fe(jelem_type)%p%get_vertices_vef()
             do jvef=1,jelem_num_vefs
                jelem_num_vef_verts = reference_fe(jelem_type)%p%get_number_vertices_vef(jelem_first_vef_id+jvef-1)
                if(jelem_num_vef_verts==nnodb) then
                   vertex_of_jvef = 0
                   do jvert=1,jelem_num_vef_verts
                      vertex_of_jvef(jvert)=mesh%lnods( mesh%pnods(jelem)-1+vertices_jvef%l(vertices_jvef%p(jelem_first_vef_id+jvef-1)+jvert-1))
                   end do
                   count=0
                   do ivert=1,ielem_num_vef_verts
                      do jvert=1,jelem_num_vef_verts
                         if(vertex_of_ivef(ivert)==vertex_of_jvef(jvert)) then
                            count=count+1
                            exit
                         end if
                      end do
                   end do
                   equal=(count==ielem_num_vef_verts)
                   if(equal) then ! Fill it
                      mesh%lvef_geo( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbgeo(iboun)
                      mesh%lvef_set( mesh%lvefs(mesh%pvefs(jelem)-1+jelem_first_vef_id-1+jvef)  ) = mesh%lbset(iboun)
                      exit elems2
                   end if
                end if
             end do
          end do elems2
       end if
    end do

  end subroutine mesh_generate_vefs

end module mesh_names
