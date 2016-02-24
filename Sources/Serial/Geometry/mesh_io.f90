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
module mesh_io_names
  use types_names
  use stdio_names
  use memor_names
  use mesh_names
  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: mesh_read, mesh_write_gid, &
          & mesh_compose_name, mesh_write_file, mesh_write_gid_files, &
          & mesh_write_files, mesh_read_file, mesh_read_files

  !integer(ip), parameter :: permu_2DQ1(4) = (/ 1, 2, 4, 3/)
  !integer(ip), parameter :: permu_3DQ1(8) = (/ 1, 2, 4, 3, 5, 6, 8, 7/)
  integer(ip), target :: permu_2DQ1(4) = (/ 1, 2, 4, 3/)
  integer(ip), target :: permu_3DQ1(8) = (/ 1, 2, 4, 3, 5, 6, 8, 7/)
  ! permu_2D_nnode_8 = (/ 1, 2, 4, 3, 5, 7, 8, 6/)
  ! permu_3D_nnode_26 = (/ 1, 9, 2, 12, 21, 10, 4, 11, 3, 13, 22, 14, 25, 27, 23, 16, 24, 15, 5, 17, 6, 20, 26, 18, 8, 19, 7/)

contains

  !=============================================================================
  ! TODO:
  !  * use stdio functions instead of line to read from files (thus allowing
  !    the posibility of using comments, checking errors, etc.)
  !  * implement other formats (when needed)
  !
  !=============================================================================
  subroutine mesh_read_file(lunio,msh,permute_c2z)
    !------------------------------------------------------------------------
    !
    ! This routine reads a mesh writen by GiD according to fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: lunio
    type(mesh_t)     , intent(out) :: msh
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
    logical                  :: permute_c2z_, apply_perm

    if(present(permute_c2z)) then
       permute_c2z_ = permute_c2z
    else
       permute_c2z_ = .true.
    end if

    ! Read first line: "MESH dimension  3  types  1  elements        352  nodes        100  boundaries        144"
    read(lunio,'(a14,1x,i2,a7,1x,i2,a10,1x,i10,a7,1x,i10,a12,1x,i10)') &
         & dum1,msh%ndime,dum2,msh%nelty,dum3,msh%nelem,dum4,msh%npoin,dum5,msh%nboun

    write(*,*) msh%ndime,msh%nelty,msh%nelem,msh%npoin,msh%nboun

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
    ! Transform length to header and get mesh%nnode
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
    if(msh%ndime == 2) then
       if(msh%nnode == 3) then     ! Linear triangles (2DP1)
          permute_c2z_ = .false. 
       elseif(msh%nnode == 4) then ! Linear quadrilaterals(2DQ1)
          permu => permu_2DQ1
       end if
    elseif(msh%ndime == 3) then
       if(msh%nnode == 4) then     ! Linear tetrahedra (3DP1)
          permute_c2z_ = .false. 
       elseif(msh%nnode == 8) then ! Linear hexahedra (3DQ1)
          permu => permu_3DQ1
       end if
    end if
    if(permute_c2z_) then
       call memalloc(msh%nnode, lnods_aux, __FILE__, __LINE__)
       do ielem = 1,msh%nelem
          lnods_aux = msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1)
          do inode = 1, msh%pnods(ielem+1) - msh%pnods(ielem)
             msh%lnods(msh%pnods(ielem)+inode-1) = lnods_aux(permu(inode))
          end do
       end do
       call memfree(lnods_aux,__FILE__,__LINE__)
    end if
 
  end subroutine mesh_read_file

  !=============================================================================
  subroutine mesh_write_file (lunio,msh,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in the format defined by GiD fempar problem type.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)  , intent(in)           :: lunio
    type(mesh_t) , intent(in)           :: msh
    character(*) , intent(in), optional :: title

    integer(ip)                    :: ielem, idime, ipoin, inode, iboun

    ! Read first line: "MESH dimension  3  types  1  elements        352  nodes        100  boundaries        144"
    write(lunio,'(a14,1x,i2,a7,1x,i2,a10,1x,i10,a7,1x,i10,a12,1x,i10)') &
         & 'MESH dimension',msh%ndime,'  types',msh%nelty,'  elements', &
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
  subroutine mesh_write_gid_file (lunio,msh,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a mesh in GiD format. Only works for nelty=1
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)           :: lunio
    type(mesh_t)   , intent(in)           :: msh
    character(*)     , intent(in), optional :: title

    integer(ip)                    :: ielem, idime, ipoin, inode
    character(13)                  :: elemt
    character(len=:), allocatable  :: title_

    assert(msh%nelty==1)

    if(msh%ndime==2) then
       if(msh%nnode==3.or.msh%nnode==6.or.msh%nnode==7) then
          elemt='Triangle' 
       else
          elemt='Quadrilateral'
       end if
    else
       if(msh%nnode==4.or.msh%nnode==10) then 
          elemt='Tetrahedra'
       else if(msh%nnode==8.or.msh%nnode==20.or.msh%nnode==27.or.msh%nnode==64) then 
          elemt='Hexahedra'
       else if(msh%nnode==6.or.msh%nnode==15) then 
          elemt='Prism'
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
    write(lunio,2)'elements'
    do ielem=1,msh%nelem
       write(lunio,4) ielem, &
            &  (msh%lnods(inode+(ielem-1)*msh%nnode),inode=1,msh%nnode),1
    end do
    write(lunio,2) 'end elements'

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i10, 3(1x,e16.8e3))
4   format(i10,65(1x,i10))
5   format('BOUNDARY ',a,' Nnodb ',i2)
6   format(i6,10(1x,i6))

  end subroutine mesh_write_gid_file

  subroutine mesh_compose_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.mesh'
  end subroutine mesh_compose_name
  subroutine mesh_compose_gid_name ( prefix, name ) 
    implicit none
    character(len=*)             , intent(in)    :: prefix 
    character(len=:), allocatable, intent(inout) :: name
    name = trim(prefix) // '.post.msh'
  end subroutine mesh_compose_gid_name

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
        call mesh_write_file(lunio,lmesh(i))
        call io_close(lunio)
     end do
     
     ! name, and rename should be automatically deallocated by the compiler when they
     ! go out of scope. Should we deallocate them explicitly for safety reasons?
   end subroutine mesh_write_files

  subroutine mesh_write_gid_files ( dir_path, prefix, nparts, lmesh )
     implicit none
     ! Parameters 
     character(*)   , intent(in)  :: dir_path 
     character(*)   , intent(in)  :: prefix
     integer(ip)     , intent(in)  :: nparts
     type(mesh_t)  , intent(in)  :: lmesh (nparts)

     character(len=:), allocatable :: name, rename ! Deferred-length allocatable character arrays

     ! Locals 
     integer (ip) :: i,lunio

     call mesh_compose_gid_name ( prefix, name )

     do i=nparts, 1, -1  
        rename=name
        call numbered_filename_compose(i,nparts,rename)
        lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
        call mesh_write_gid_file(lunio,lmesh(i))
        call io_close(lunio)
     end do
     
   end subroutine mesh_write_gid_files


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
        call mesh_read_file(lunio,lmesh(i))
        call io_close(lunio)
     end do
     
   end subroutine mesh_read_files

  !=============================================================================
   subroutine mesh_read ( dir_path, prefix, f_mesh, permute_c2z )
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     type(mesh_t)               , intent(out) :: f_mesh
     logical, optional, intent(in)  :: permute_c2z

     ! Locals
     integer                        :: iam, num_procs
     integer(ip)                    :: j, ndigs_iam, ndigs_num_procs, lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call mesh_compose_name ( prefix, name )

     lunio = io_open( trim(dir_path)//'/'//trim(name), 'read', status='old' )
     call mesh_read_file(lunio, f_mesh, permute_c2z)
     call io_close(lunio)

   end subroutine mesh_read

   !=============================================================================
   subroutine mesh_write_gid ( dir_path, prefix, f_mesh )
     implicit none 
     ! Parameters
     character (*)                , intent(in)  :: dir_path
     character (*)                , intent(in)  :: prefix
     type(mesh_t)               , intent(in)  :: f_mesh

     ! Locals
     integer(ip)                    :: lunio
     character(len=:), allocatable  :: name

     ! Read mesh
     call mesh_compose_gid_name ( prefix, name )
     
     lunio = io_open( trim(dir_path)//'/'//trim(name), 'write' )
     call mesh_write_gid_file(lunio, f_mesh)
     call io_close(lunio)

   end subroutine mesh_write_gid
   
end module mesh_io_names
