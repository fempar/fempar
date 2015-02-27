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
module fem_mesh_io
  use types
  use stdio
  use memor
  use renum_names
  use fem_mesh_names
  use fem_space_names
  implicit none
  private

  ! Functions
  public :: fem_mesh_read, fem_mesh_write, &
          & fem_mesh_compose_name, fem_mesh_write_files, &
          & fem_mesh_read_files, &
          & fem_mesh_write_vtk_w_dof_handler

contains

  !=============================================================================
  ! TODO:
  !  * use stdio functions instead of line to read from files (thus allowing
  !    the posibility of using comments, checking errors, etc.)
  !  * implement other formats (when needed)
  !
  !=============================================================================
  subroutine fem_mesh_read(lunio,msh,flag)
    !------------------------------------------------------------------------
    !
    ! This routine reads a fem_mesh in GiD format.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)      , intent(in)  :: lunio
    type(fem_mesh)   , intent(out) :: msh
    logical, optional, intent(in)  :: flag
    integer(ip)       :: i,inode,idime,ipoin,jpoin,ielem,jelem,istat,iboun
    character(1000)     :: tel
    integer(ip), allocatable :: aux(:),permu(:)
    logical                  :: flag_

    if(present(flag)) then
       flag_ = flag
    else
       flag_ = .false.
    end if

    ! Read first line to get ndime and nnode 
    read(lunio,'(a)') tel

    do i=1,75
       if(tel(i:i+8)=='dimension') then
          read(tel(i+9:i+10),*) msh%ndime
       end if
       if(tel(i:i+4)=='Nnode') then
          read(tel(i+5:80),*) msh%nnode
       end if
    end do
    msh%nelty=1
   
    ! Count nodes
    do while(tel(1:5).ne.'coord')
       read(lunio,'(a)') tel
    end do
    ipoin = 0
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end c')
       ipoin = ipoin + 1
       read(lunio,'(a)') tel
    end do
    msh%npoin=ipoin
    if(msh%npoin>0) call memalloc(msh%ndime,msh%npoin,msh%coord,__FILE__,__LINE__)
    
    call io_rewind(lunio)
    ! Read nodes
    do while(tel(1:5).ne.'coord')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end c')
       read(tel,*) jpoin,(msh%coord(idime,jpoin),idime=1,msh%ndime)
       read(lunio,'(a)') tel
    end do

    ! Count elements
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    ielem = 0
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       ielem = ielem + 1
       read(lunio,'(a)') tel
    end do
    msh%nelem = ielem
    call memalloc(msh%nelem+1,msh%pnods,__FILE__,__LINE__)
    do ielem = 1, msh%nelem+1
       msh%pnods(ielem) = (ielem-1)*msh%nnode+1
    end do
    call memalloc(msh%nnode*msh%nelem,msh%lnods,__FILE__,__LINE__)

    call io_rewind(lunio)
    ! Read elements
    do while(tel(1:5).ne.'eleme')
       read(lunio,'(a)') tel
    end do
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end e')
       read(tel,*) jelem,(msh%lnods(inode+(jelem-1)*msh%nnode),inode=1,msh%nnode)
       read(lunio,'(a)') tel
    end do

    ! Read nodes (if there aren't coordinates)
    ipoin = 0
    call memalloc(msh%nnode*msh%nelem,aux,__FILE__,__LINE__)
    aux = 0
    do ielem = 1,msh%nelem
       do inode = 1,msh%nnode
          if(aux(msh%lnods(inode+(ielem-1)*msh%nnode)) == 0) then
             aux(msh%lnods(inode+(ielem-1)*msh%nnode)) = 1
             ipoin = ipoin + 1
          end if
       end do
    end do
    msh%npoin = ipoin       
    call memfree(aux,__FILE__,__LINE__)

    ! Reordering the nodes of the mesh
    call memalloc(msh%nnode,permu,__FILE__,__LINE__)
    if(msh%ndime == 2) then        ! 2D
       if(msh%nnode == 3) then     ! Linear triangles (P1)
          permu = (/ 1, 2, 3/)
       elseif(msh%nnode == 4) then ! Linear quads (Q1)
          if(flag_) then
             permu = (/ 1, 2, 3, 4/) ! No permuted
          else
             permu = (/ 1, 2, 4, 3/)
          end if
       elseif(msh%nnode == 8) then ! Quadratic quads (Q2)
          if(flag_) then
             permu = (/ 1, 2, 3, 4, 5, 6, 7, 8/) ! No permuted
          else
             permu = (/ 1, 2, 4, 3, 5, 7, 8, 6/)
          end if
       end if
    elseif(msh%ndime == 3) then    ! 3D
       if(msh%nnode == 4) then     ! Linear tetrahedra (P1)
          permu = (/ 1, 2, 3, 4/)
       elseif(msh%nnode == 8) then ! Linear hexahedra (Q1)
          if(flag_) then
             permu = (/ 1, 2, 3, 4, 5, 6, 7, 8/) ! No permuted
          else
             permu = (/ 1, 2, 4, 3, 5, 6, 8, 7/)
          end if
       elseif(msh%nnode == 26) then ! Quadratic hexahedra (Q2)
          if(flag_) then
             permu = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, &
                  &     22, 23, 24, 25, 26/) ! No permuted
          else
             permu = (/ 1, 9, 2, 12, 21, 10, 4, 11, 3, 13, 22, 14, 25, 27, 23, 16, 24, 15, 5, 17, 6, &
                  &     20, 26, 18, 8, 19, 7/)
          end if
       end if
    end if
    call memalloc(msh%nnode,aux,__FILE__,__LINE__)
    do ielem = 1,msh%nelem
       aux = msh%lnods(msh%pnods(ielem):msh%pnods(ielem+1)-1)
       do i = 1,msh%nnode
          msh%lnods(msh%pnods(ielem)+i-1) = aux(permu(i))
       end do
    end do
    call memfree(aux,__FILE__,__LINE__)
    call memfree(permu,__FILE__,__LINE__)
 
    ! Count boundary elements
    do while(tel(1:5).ne.'bound')
       read(lunio,'(a)',IOSTAT=istat) tel
       if(istat<0) return                      ! Check End-of-file for old files
    end do

    iboun = 0
    read(lunio,'(a)') tel
    do while(tel(1:5).ne.'end b')
       iboun = iboun + 1
       read(lunio,'(a)') tel
    end do
    msh%nboun = iboun

    if(msh%nboun.gt.0) then
       call io_rewind(lunio)

       ! Read nnodb
       do while(tel(1:5).ne.'end e')
          read(lunio,'(a)') tel
       end do
        read(lunio,'(a)') tel
       do i=1,50
          if(tel(i:i+4)=='Nnodb') then
             read(tel(i+5:50),*) msh%nnodb
          end if
       end do
       call memalloc(1,msh%pboun,__FILE__,__LINE__)
       call memalloc(msh%nnodb*msh%nboun,msh%lboun,__FILE__,__LINE__)
       call memalloc(msh%nnodb+1,msh%nboun,msh%lboel,__FILE__,__LINE__)

       ! Read boundary elements
       read(lunio,'(a)') tel
       read(lunio,'(a)') tel   ! We need 2 reads to get to the boundary elements
       do while(tel(1:5).ne.'end b')
          read(tel,*) iboun,(msh%lboun(inode+(iboun-1)*msh%nnodb),inode=1,msh%nnodb),&
               &            (msh%lboel(inode,iboun),inode=1,msh%nnodb), msh%lboel(msh%nnodb+1,iboun)
          read(lunio,'(a)') tel
       end do
    else
      msh%nnodb=0
    end if

    msh%nelpo = 0

    
    return

  end subroutine fem_mesh_read

  !=============================================================================
  subroutine fem_mesh_write(lunio,msh,nren,eren,title)
    !------------------------------------------------------------------------
    !
    ! This routine writes a fem_mesh in GiD format.
    !
    !------------------------------------------------------------------------
    implicit none
    integer(ip)  , intent(in) :: lunio
    type(fem_mesh)   , intent(in) :: msh
    type(renum)  , intent(in), optional :: nren,eren
    character(*) , intent(in), optional :: title
    integer(ip)    :: ielem,idime,ipoin,inode
    character(13)  :: elemt
    character(80)  :: title_

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
write(*,*) lunio
    write(lunio,1) adjustl(trim(title_)),msh%ndime,adjustl(trim(elemt)),msh%nnode

    ! Coordinates
    write(lunio,2)'coordinates'
    if (allocated(msh%coord)) then
       if(present(nren)) then
          do ipoin=1,msh%npoin
             write(lunio,3) nren%lperm(ipoin),(msh%coord(idime,ipoin),idime=1,msh%ndime)
          end do
       else
          do ipoin=1,msh%npoin
             write(lunio,3) ipoin,(msh%coord(idime,ipoin),idime=1,msh%ndime)
          end do
       end if
    end if
    write(lunio,2)'end coordinates'

    ! Connectivity
    write(lunio,2)'elements'
    if(present(nren).and.present(eren)) then
       do ielem=1,msh%nelem
          write(lunio,4) eren%lperm(ielem), &
             &  (nren%lperm(msh%lnods(inode+(ielem-1)*msh%nnode)),inode=1,msh%nnode),1
       end do
    else if(present(nren)) then
       do ielem=1,msh%nelem
          write(lunio,4) ielem, &
             &  (nren%lperm(msh%lnods(inode+(ielem-1)*msh%nnode)),inode=1,msh%nnode),1
       end do
    else if(present(eren)) then
       do ielem=1,msh%nelem
          write(lunio,4) eren%lperm(ielem), &
             &  (msh%lnods(inode+(ielem-1)*msh%nnode),inode=1,msh%nnode),1
       end do
    else
       do ielem=1,msh%nelem
          write(lunio,4) ielem, &
             &  (msh%lnods(inode+(ielem-1)*msh%nnode),inode=1,msh%nnode),1
       end do
    end if
    write(lunio,2)'end elements'

    ! Boundary elements
    if(msh%nboun>0) then
       write(lunio,5) adjustl(trim(title_)),msh%nnodb
       write(lunio,2)'boundaries'
       do ielem=1,msh%nboun
          write(lunio,6) ielem, (msh%lboun(inode+(ielem-1)*msh%nnodb),inode=1,msh%nnodb), &
               &         (msh%lboel(inode,ielem),inode=1,msh%nnodb), msh%lboel(msh%nnodb+1,ielem)
       end do
       write(lunio,2)'end boundaries'
    else
       write(lunio,2)'boundaries'
       write(lunio,2)'end boundaries'
    end if

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i10, 3(1x,e16.8e3))
4   format(i10,65(1x,i10))
5   format('BOUNDARY ',a,' Nnodb ',i2)
6   format(i6,10(1x,i6))

  end subroutine fem_mesh_write

  subroutine fem_mesh_compose_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.msh'
  end subroutine 

  subroutine fem_mesh_write_files ( dir_path, prefix, nparts, lmesh, nren, eren )
     implicit none
     ! Parameters 
     character *(*), intent(in)       :: dir_path 
     character *(*), intent(in)       :: prefix
     integer(ip)   , intent(in)       :: nparts
     type(fem_mesh), intent(in)       :: lmesh (nparts)
     type(renum)  , intent(in), optional :: nren(nparts),eren(nparts)
     character(256)                   :: name, rename

     ! Locals 
     integer (ip)                     :: i,lunio

     call fem_mesh_compose_name ( prefix, name )

     do i=nparts, 1, -1  
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'write' )
       call fem_mesh_write(lunio,lmesh(i),nren(i),eren(i))
       call io_close(lunio)
     end do

  end subroutine fem_mesh_write_files 

  subroutine fem_mesh_read_files ( dir_path, prefix, nparts, lmesh )
    implicit none
    ! Parameters 
    character *(*), intent(in)       :: dir_path 
    character *(*), intent(in)       :: prefix
    integer(ip)   , intent(in)       :: nparts
    type(fem_mesh), intent(out)      :: lmesh (nparts)
    character(256)                   :: name, rename
    
    ! Locals 
    integer (ip)                     :: i,lunio
    
    call fem_mesh_compose_name ( prefix, name )
    
    do i=nparts, 1, -1  
       rename=name
       call numbered_filename_compose(i,nparts,rename)
       lunio = io_open( trim(dir_path) // '/' // trim(rename), 'read' )
       call fem_mesh_read(lunio,lmesh(i))
       call io_close(lunio)
    end do
    
  end subroutine fem_mesh_read_files

  !=============================================================================
  subroutine fem_mesh_write_vtk_w_dof_handler ( femsp, X, Y, Z, T, offset, Ttype )
     implicit none
     ! Parameters 
     type(fem_space),          intent(in)    :: femsp
     real(rp), allocatable,    intent(inout) :: X(:),Y(:),Z(:)
     integer(ip), allocatable, intent(inout) :: T(:), offset(:)
     integer(1), allocatable,  intent(inout) :: Ttype(:)
     ! Locals 
     integer (ip) :: ipoin, ielem, iobje, npoin, count

     ! Count nodes
     npoin = 0
     do ielem = 1,femsp%g_mesh%nelem
        npoin = npoin + femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1
     end do

     ! Allocate mesh
     call memalloc (npoin, X, __FILE__,__LINE__)
     call memalloc (npoin, Y, __FILE__,__LINE__)
     call memalloc (npoin, Z, __FILE__,__LINE__)
     ! Allocate connectivities
     call memalloc (npoin, T, __FILE__,__LINE__)
     call memalloc (femsp%g_mesh%nelem, offset, __FILE__,__LINE__)

     ! Write connectivites
     !!!! IMPORTANT: Paraview counts the points starting from 0 !!!!!
     do ipoin=1,npoin
        T(ipoin) = ipoin-1
     end do
     
     ! Generate a list of coordinates for each node per element (consecutively)
     X=0.0_rp; Y=0.0_rp; Z=0.0_rp
     offset = 0
     count = 1
     do ielem = 1,femsp%g_mesh%nelem
        do iobje = 1,femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1
           ipoin = femsp%g_mesh%lnods(femsp%g_mesh%pnods(ielem)+iobje-1)
           ! Write mesh
           X(count) = femsp%g_mesh%coord(1,ipoin)
           Y(count) = femsp%g_mesh%coord(2,ipoin)
           if(femsp%g_mesh%ndime==3) Z(count) = femsp%g_mesh%coord(3,ipoin)
           count = count + 1
        end do
        offset(ielem) = count - 1
     end do
     
     ! Cell types ( Linear hexahedron = 12; Linear tetrahedra = 10 )
     !call memalloc (gelem, Ttype, __FILE__,__LINE__)
     allocate(Ttype(femsp%g_mesh%nelem))
     do ielem = 1,femsp%g_mesh%nelem
        if(femsp%g_mesh%ndime==2) then
           if(femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1==4) then
              Ttype(ielem) = 8
           elseif(femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1==3) then
              Ttype(ielem) = 5
           end if
        elseif(femsp%g_mesh%ndime==3) then
           if(femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1==8) then
              Ttype(ielem) = 11
           elseif(femsp%lelem(ielem)%f_inf(1)%p%nobje_dim(2)-1==4) then
              Ttype(ielem) = 10
           end if
        end if
     end do   
  
   end subroutine fem_mesh_write_vtk_w_dof_handler

end module fem_mesh_io
