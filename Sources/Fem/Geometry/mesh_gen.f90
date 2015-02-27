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
module fem_mesh_gen
  use types
  use memor
  use fem_mesh_names
!  use fem_conditions_names
# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: periodic = -1

  type mesh_size
     integer(ip)                :: &
        lverb=0,                   &         ! Verbosity level
        ntdix=0,                   &         ! Type of discretization in x (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        ntdiy=0,                   &         ! Type of discretization in y (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        ntdiz=0,                   &         ! Type of discretization in z (0=uniform, 1=cubic, 2=tanh, 3=imh+unif, 4:imh+tanh)
        pdegr=1,                   &         ! Interpolation order p (1=liniar, 2=quadratic, 3=cubic) 
        mater=0,                   &         ! Material case
        ierrc=0,                   &         ! Error code
        neblx=0,                   &         ! Number of elements in the x boundary layer
        nebly=0,                   &         ! Number of elements in the y boundary layer
        neblz=0                              ! Number of elements in the z boundary layer

     real(rp)                   :: &
        xleng   = 1.0_rp,          &         ! Size of the domain in x
        yleng   = 1.0_rp,          &         ! Size of the domain in y
        zleng   = 1.0_rp,          &         ! Size of the domain in z
        zx1     = 0.1_rp,          &         ! size of the elements at x=0   (left)  
        zx2     = 0.1_rp,          &         ! size of the elements at x=a/2 (center)
        zy1     = 0.1_rp,          &         ! size of the elements at y=0   (bottom)
        zy2     = 0.1_rp,          &         ! size of the elements at y=b/2 (center)
        x0      = 0.0_rp,          &         ! Origin x-coordinate
        y0      = 0.0_rp,          &         ! Origin y-coordinate
        z0      = 0.0_rp,          &         ! Origin z-coordinate
        atol    = 1.0e-8_rp,       &         ! Absolute tolerance
        stret   = 0.0_rp,          &         ! Stretching parameter (old)
        xstret  = 2.75_rp,         &         ! Stretching parameter
        ystret  = 2.75_rp,         &         ! Stretching parameter
        zstret  = 2.75_rp,         &         ! Stretching parameter
        xlengbl = 0.0_rp,          &         ! Size of the boundary layer in x
        ylengbl = 0.0_rp,          &         ! Size of the boundary layer in y
        zlengbl = 0.0_rp                     ! Size of the boundary layer in z
  end type mesh_size

  ! Types
  public :: mesh_size

  ! Functions
  public :: fem_mesh_box, fem_mesh_fold

contains
  !================================================================================================
  !
  ! Public methods:
  ! - fem_mesh_box
  ! - fem_mesh_fold
  !
  ! Auxiliar routines
  ! mesh_rectangular_elements
  ! mesh_rectangular_nodes
  !
  !================================================================================================
  subroutine fem_mesh_fold(orig,isper,fold,list)
    !------------------------------------------------------------------------
    !
    ! This routine generates a periodic mesh by folding a standard one,
    ! and a list of the nodes correspondence.
    !
    !------------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(in)  :: orig
    integer(ip)   , intent(in)  :: isper(3)
    type(fem_mesh), intent(out) :: fold
    integer(ip)   , intent(out) :: list(orig%npoin)

    integer(ip) :: npoix,npoiy,npoiz,nelex,neley,nelez,nedir(3)
    integer(ip) :: ic,i,j,k

    assert(orig%mtype==structured)

    ! Allocate fold mesh
    !nedir=orig%nedir-isper
    call fem_mesh_alloc(orig%ndime,orig%nnode,orig%nedir,isper,fold)

    if(orig%ndime==2) then

       nelex=orig%nedir(1)
       neley=orig%nedir(2)
       if(orig%nnode==16.or.orig%nnode==10) then
          npoiy=3*neley+1
          npoix=3*nelex+1
       else if(orig%nnode==9.or.orig%nnode==6) then
          npoiy=2*neley+1
          npoix=2*nelex+1
       else if(orig%nnode==4.or.orig%nnode==3) then
          npoiy=neley+1
          npoix=nelex+1
       end if

       ! Periodic only in x direction
       if(fold%isper(1)==1 .and. fold%isper(2)==0) then
          ic=0
          do i=1,npoix-1
             do j=1,npoiy
                ic=ic+1
                !ic=(i-1)*npoiy+j
                list(ic)=(i-1)*npoiy+j
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             end do
          end do
          do j=1,npoiy
             ic=ic+1
             list(ic)=j
          end do

       else
          ic=0
          do i=1,npoix-1
             do j=1,npoiy-1
                ic=ic+1
                !ic=(i-1)*npoiy+j
                list(ic)=(i-1)*(npoiy-1)+j
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             end do
             j=npoiy
             ic=ic+1
             !ic=(i-1)*npoiy+j
             if(fold%isper(2)==0) then
                list(ic)=(i-1)*(npoiy-1)+j
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(2)==1) then
                list(ic)=(i-1)*(npoiy-1)+1
             end if
          end do
          i=npoix
          if(fold%isper(1)==0) then
             do j=1,npoiy-1
                ic=ic+1
                !ic=(i-1)*npoiy+j
                list(ic)=(i-1)*(npoiy-1)+j
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             end do
             j=npoiy
             ic=ic+1
             !ic=(i-1)*npoiy+j
             if(fold%isper(2)==0) then
                list(ic)=(i-1)*(npoiy-1)+j
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(2)==1) then
                list(ic)=(i-1)*(npoiy-1)+1
             end if
          else if(fold%isper(1)==1) then
             do j=1,npoiy-1
                ic=ic+1
                !ic=(i-1)*npoiy+j
                list(ic)=j
             end do
             j=npoiy
             ic=ic+1
             !ic=(i-1)*npoiy+j
             if(fold%isper(2)==0) then
                list(ic)=j
             else if(fold%isper(2)==1) then
                list(ic)=1
             end if
          end if
       end if

    else if(orig%ndime==3) then

       nelex=orig%nedir(1)
       neley=orig%nedir(2)
       nelez=orig%nedir(3)
       if(orig%nnode==8) then  !Q1
          npoiy=neley+1
          npoix=nelex+1
          npoiz=nelez+1
       end if

       ic=0
       do i=1,npoix-1
          do j=1,npoiy-1
             do k=1,npoiz-1
                ic=ic+1
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             end do
             k=npoiz
             ic=ic+1
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+1
             end if
          end do
          j=npoiy
          do k=1,npoiz-1
             ic=ic+1
             if(fold%isper(2)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(2)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+k
             end if
          end do
          k=npoiz
          ic=ic+1
          if(fold%isper(2)==0) then
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+1
             end if
          else if(fold%isper(2)==1) then
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+k
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+1
             end if
          end if
       end do
       i=npoix
       if(fold%isper(1)==0) then
          do j=1,npoiy-1
             do k=1,npoiz-1
                ic=ic+1
                !ic=(i-1)*npoiy+j
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             end do
             k=npoiz
             ic=ic+1
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+1
             end if
          end do
          j=npoiy
          do k=1,npoiz-1
             ic=ic+1
             if(fold%isper(2)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(2)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+k
             end if
          end do
          k=npoiz
          ic=ic+1
          if(fold%isper(2)==0) then
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+k
                if(allocated(fold%coord)) then
                   fold%coord(:,list(ic))=orig%coord(:,ic)
                end if
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+(j-1)*(npoiz-1)+1
             end if
          else if(fold%isper(2)==1) then
             if(fold%isper(3)==0) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+k
             else if(fold%isper(3)==1) then
                list(ic)=(i-1)*(npoiy-1)*(npoiz-1)+1
             end if
          end if
       else if(fold%isper(1)==1) then
          do j=1,npoiy-1
             do k=1,npoiz-1
                ic=ic+1
                list(ic)=(j-1)*(npoiz-1)+k
             end do
             k=npoiz
             ic=ic+1
             if(fold%isper(3)==0) then
                list(ic)=(j-1)*(npoiz-1)+k
             else if(fold%isper(3)==1) then
                list(ic)=(j-1)*(npoiz-1)+1
             end if
          end do
          j=npoiy
          do k=1,npoiz-1
             ic=ic+1
             if(fold%isper(2)==0) then
                list(ic)=(j-1)*(npoiz-1)+k
             else if(fold%isper(2)==1) then
                list(ic)=k
             end if
          end do
          k=npoiz
          ic=ic+1
          if(fold%isper(2)==0) then
             if(fold%isper(3)==0) then
                list(ic)=(j-1)*(npoiz-1)+k
             else if(fold%isper(3)==1) then
                list(ic)=(j-1)*(npoiz-1)+1
             end if
          else if(fold%isper(2)==1) then
             if(fold%isper(3)==0) then
                list(ic)=k
             else if(fold%isper(3)==1) then
                list(ic)=1
             end if
          end if
       end if
    end if

    ! Generate lnods, coord
    fold%lnods = list(orig%lnods)

  end subroutine fem_mesh_fold

  !================================================================================================
  subroutine fem_mesh_box(size,msh)
    !------------------------------------------------------------------------
    !
    ! This routine generates a structured mesh in a box.
    !
    !------------------------------------------------------------------------
    implicit none
    type(mesh_size)      , intent(in)    :: size
    type(fem_mesh)       , intent(inout) :: msh
    integer(ip) :: npoix,npoiy,npoiz,nelex,neley,nelez

    assert(msh%mtype==structured)

    if(msh%ndime==2) then

       nelex=msh%nedir(1)
       neley=msh%nedir(2)

       if(msh%nnode==16.or.msh%nnode==10) then
          npoiy=3*neley+1
          npoix=3*nelex+1
       else if(msh%nnode==9.or.msh%nnode==6) then
          npoiy=2*neley+1
          npoix=2*nelex+1
       else if(msh%nnode==4.or.msh%nnode==3) then
          npoiy=neley+1
          npoix=nelex+1
       end if
       ! Generate elements
       call mesh_rectangular_elements_2d(nelex,neley,msh%nnode,msh%nelem,msh%lnods)

       if(allocated(msh%coord)) &    ! Generate nodes
            & call mesh_rectangular_nodes_2d(npoix,npoiy,msh%ndime,msh%npoin,size%ntdix, &
            &                             size%ntdiy,size%xleng,size%yleng,size%zx1, &
            &                             size%zx2,size%zy1,size%zy2,msh%coord)

    else if(msh%ndime==3) then
       
       nelex=msh%nedir(1)
       neley=msh%nedir(2)
       nelez=msh%nedir(3)

       if(msh%nnode==8) then
          npoiy=neley+1
          npoix=nelex+1
          npoiz=nelez+1
       end if
       ! Generate elements
       call mesh_rectangular_elements_3d(nelex,neley,nelez,msh%nnode,msh%nelem,msh%lnods)

       if(allocated(msh%coord)) &    ! Generate nodes
            & call mesh_rectangular_nodes_3d(npoix,npoiy,npoiz,msh%ndime,msh%npoin,size%ntdix, &
            &                             size%ntdiy,size%ntdiz,size%xleng,size%yleng,size%zleng, &
            &                             size%zx1,size%zx2,size%zy1,size%zy2,msh%coord)

    end if

    return

  end subroutine fem_mesh_box

  !================================================================================================
  subroutine mesh_rectangular_elements_2d (nelex,neley,nnode,nelem,lnods)
    !subroutine mesh_rectangular_elements (nelex,neley,nnode,nelem,ln)
    implicit none
    integer(ip), intent(in)  :: nelex,neley,nnode,nelem
    integer(ip), intent(out) :: lnods(nnode,nelem)
    ! If there is some compiler problem, try
    !integer(ip) ln(nnode*nelem)
    !#define lnods(i,j) ln(((j)-1)*(nnode)+(i))

    integer(ip)              :: i,j,k,ic,ib,npoiy

    ic=0

    if(nnode==4) then       ! Q1 elements

       do i=1,nelex
          do j=1,neley
             ic=ic+1
             !ic=(i-1)*neley+j
             lnods(1,ic)=(i-1)*(neley+1)+j
             lnods(2,ic)=i*(neley+1)+j
             lnods(4,ic)=i*(neley+1)+j+1
             lnods(3,ic)=(i-1)*(neley+1)+j+1
          end do
       end do

! This is a stupid idea. It is much better to generate
! standard mesh and then fold it (so no recode of higher
! order elements is needed and code is far more clear).
! Now isper is deleted from the argument list.

!!$       do i=1,nelex-1
!!$          do j=1,neley-1
!!$             ic=ic+1
!!$             !ic=(i-1)*neley+j
!!$             lnods(1,ic)=(i-1)*(neley+1)+j
!!$             lnods(2,ic)=i*(neley+1)+j
!!$             lnods(3,ic)=i*(neley+1)+j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          end do
!!$          y=neley
!!$          ic=ic+1
!!$          !ic=(i-1)*neley+j
!!$          lnods(1,ic)=(i-1)*(neley+1)+j
!!$          lnods(2,ic)=i*(neley+1)+j
!!$          if(isper(2)==0) then
!!$             lnods(3,ic)=i*(neley+1)+j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          else if(isper(2)==1) then
!!$             lnods(3,ic)=i*(neley+1)+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+1
!!$          end if
!!$       end do
!!$
!!$       i=nelex
!!$       if(isper(1)==0) then
!!$          do j=1,neley-1
!!$             ic=ic+1
!!$             !ic=(i-1)*neley+j
!!$             lnods(1,ic)=(i-1)*(neley+1)+j
!!$             lnods(2,ic)=i*(neley+1)+j
!!$             lnods(3,ic)=i*(neley+1)+j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          end do
!!$          y=neley
!!$          ic=ic+1
!!$          !ic=(i-1)*neley+j
!!$          lnods(1,ic)=(i-1)*(neley+1)+j
!!$          lnods(2,ic)=i*(neley+1)+j
!!$          if(isper(2)==0) then
!!$             lnods(3,ic)=i*(neley+1)+j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          else if(isper(2)==1) then
!!$             lnods(3,ic)=i*(neley+1)+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+1
!!$          end if
!!$       else if(isper(1)==1) then
!!$          do j=1,neley-1
!!$             ic=ic+1
!!$             !ic=(i-1)*neley+j
!!$             lnods(1,ic)=(i-1)*(neley+1)+j
!!$             lnods(2,ic)=j
!!$             lnods(3,ic)=j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          end do
!!$          y=neley
!!$          ic=ic+1
!!$          !ic=(i-1)*neley+j
!!$          lnods(1,ic)=(i-1)*(neley+1)+j
!!$          lnods(2,ic)=j
!!$          if(isper(2)==0) then
!!$             lnods(3,ic)=j+1
!!$             lnods(4,ic)=(i-1)*(neley+1)+j+1
!!$          else if(isper(2)==1) then
!!$             lnods(3,ic)=1
!!$             lnods(4,ic)=(i-1)*(neley+1)+1
!!$          end if
!!$       end if

    else if(nnode.eq.3) then       ! P1 from left to rigth 

       do i=1,nelex
          j=1
          ic=ic+1
          !ic=(i-1)*2*neley+j
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j
          lnods(3,ic)=i*(neley+1)+j+1
          ic=ic+1
          !ic=(i-1)*2*neley+j+1
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j+1
          lnods(3,ic)=(i-1)*(neley+1)+j+1
          do j=2,neley-1
             ic=ic+1
             !ic=(i-1)*2*neley+j
             lnods(1,ic)=(i-1)*(neley+1)+j
             lnods(2,ic)=i*(neley+1)+j
             lnods(3,ic)=i*(neley+1)+j+1
             ic=ic+1
             !ic=(i-1)*2*neley+j+1
             lnods(1,ic)=(i-1)*(neley+1)+j
             lnods(2,ic)=i*(neley+1)+j+1
             lnods(3,ic)=(i-1)*(neley+1)+j+1
          end do
          j=neley
          ic=ic+1
          !ic=(i-1)*2*neley+j
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j
          lnods(3,ic)=i*(neley+1)+j+1
          ic=ic+1
          !ic=(i-1)*2*neley+j+1
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j+1
          lnods(3,ic)=(i-1)*(neley+1)+j+1
       end do

    else if(nnode.eq.-3) then      ! P1 from rigth to left

       do i=1,nelex
          j=1
          ic=ic+1
          !ic=(i-1)*2*neley+j
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j
          lnods(3,ic)=(i-1)*(neley+1)+j+1
          ic=ic+1
          !ic=(i-1)*2*neley+j+1
          lnods(1,ic)=i*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j+1
          lnods(3,ic)=(i-1)*(neley+1)+j+1
          do j=2,neley-1
             ic=ic+1
             !ic=(i-1)*2*neley+j
             lnods(1,ic)=(i-1)*(neley+1)+j
             lnods(2,ic)=i*(neley+1)+j
             lnods(3,ic)=(i-1)*(neley+1)+j+1
             ic=ic+1
             !ic=(i-1)*2*neley+j+1
             lnods(1,ic)=i*(neley+1)+j
             lnods(2,ic)=i*(neley+1)+j+1
             lnods(3,ic)=(i-1)*(neley+1)+j+1
          end do
          j=neley
          ic=ic+1
          !ic=(i-1)*2*neley+j
          lnods(1,ic)=(i-1)*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j
          lnods(3,ic)=(i-1)*(neley+1)+j+1
          ic=ic+1
          !ic=(i-1)*2*neley+j+1
          lnods(1,ic)=i*(neley+1)+j
          lnods(2,ic)=i*(neley+1)+j+1
          lnods(3,ic)=(i-1)*(neley+1)+j+1
       end do

    else if(nnode.eq.9)then        ! Q2 elements
       npoiy=2*neley+1
       do i=1,nelex
          do j=1,neley
             ic=ic+1
             !ic=(i-1)*neley+j
             lnods(1,ic)=(i-1)*2*npoiy        +2*(j-1)+1
             lnods(2,ic)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
             lnods(3,ic)=(i-1)*2*npoiy+2*npoiy+2*j+1
             lnods(4,ic)=(i-1)*2*npoiy        +2*j+1
             lnods(5,ic)=(i-1)*2*npoiy+  npoiy+2*(j-1)+1
             lnods(6,ic)=(i-1)*2*npoiy+2*npoiy+2*j
             lnods(7,ic)=(i-1)*2*npoiy+  npoiy+2*j+1
             lnods(8,ic)=(i-1)*2*npoiy        +2*j
             lnods(9,ic)=(i-1)*2*npoiy+  npoiy+2*j
          end do
       end do

    else if(nnode.eq.6)then        ! P2 elements
       npoiy=2*neley+1
       do i=1,nelex
          do j=1,neley
             ic=ic+1
             !ic=(i-1)*2*neley+j
             lnods(1,ic)=(i-1)*2*npoiy        +2*(j-1)+1
             lnods(2,ic)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
             lnods(3,ic)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+3
             lnods(4,ic)=(i-1)*2*npoiy+  npoiy+2*(j-1)+1
             lnods(5,ic)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+2
             lnods(6,ic)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2
             ic=ic+1
             !ic=(i-1)*2*neley+j+1
             lnods(1,ic)=(i-1)*2*npoiy        +2*(j-1)+1
             lnods(2,ic)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+3
             lnods(3,ic)=(i-1)*2*npoiy        +2*(j-1)+3
             lnods(4,ic)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2
             lnods(5,ic)=(i-1)*2*npoiy+  npoiy+2*(j-1)+3
             lnods(6,ic)=(i-1)*2*npoiy        +2*(j-1)+2
          end do
       end do

       ! Correct corner elements that do not satisfy BB condition
       i=1
       j=neley
       ib=2*(i-1)*neley+2*(j-1)+1
       lnods(1,ib)=(i-1)*2*npoiy        +2*(j-1)+1
       lnods(2,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
       lnods(3,ib)=(i-1)*2*npoiy        +2*(j-1)+3
       lnods(4,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+1
       lnods(5,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2
       lnods(6,ib)=(i-1)*2*npoiy        +2*(j-1)+2
       ib=ib+1
       lnods(1,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
       lnods(2,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+3
       lnods(3,ib)=(i-1)*2*npoiy        +2*(j-1)+3
       lnods(4,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+2
       lnods(5,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+3
       lnods(6,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2
       i=nelex
       j=1
       ib=2*(i-1)*neley+2*(j-1)+1
       lnods(1,ib)=(i-1)*2*npoiy        +2*(j-1)+1
       lnods(2,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
       lnods(3,ib)=(i-1)*2*npoiy        +2*(j-1)+3
       lnods(4,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+1
       lnods(5,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2
       lnods(6,ib)=(i-1)*2*npoiy        +2*(j-1)+2
       ib=ib+1
       lnods(1,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+1
       lnods(2,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+3
       lnods(3,ib)=(i-1)*2*npoiy        +2*(j-1)+3
       lnods(4,ib)=(i-1)*2*npoiy+2*npoiy+2*(j-1)+2
       lnods(5,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+3
       lnods(6,ib)=(i-1)*2*npoiy+  npoiy+2*(j-1)+2

    else if(nnode.eq.16)then            ! Q3 elements
       npoiy=3*neley+1
       do i=1,nelex
          do j=1,neley
             ic=ic+1
             !ic=(i-1)*neley+j
             lnods( 1,ic)=(i-1)*3*npoiy+3*(j-1)+1
             lnods( 2,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
             lnods( 3,ic)=(i-1)*3*npoiy+3*npoiy+3*j+1
             lnods( 4,ic)=(i-1)*3*npoiy+3*j+1
             lnods( 5,ic)=(i-1)*3*npoiy+npoiy+3*(j-1)+1
             lnods( 6,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+1
             lnods( 7,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+2
             lnods( 8,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+3
             lnods( 9,ic)=(i-1)*3*npoiy+2*npoiy+3*j+1
             lnods(10,ic)=(i-1)*3*npoiy+npoiy+3*j+1
             lnods(11,ic)=(i-1)*3*npoiy+3*(j-1)+3
             lnods(12,ic)=(i-1)*3*npoiy+3*(j-1)+2
             lnods(13,ic)=(i-1)*3*npoiy+npoiy+3*(j-1)+2
             lnods(14,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
             lnods(15,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+3
             lnods(16,ic)=(i-1)*3*npoiy+npoiy+3*(j-1)+3
          end do
       end do

    else if(nnode.eq.10)then       ! P3 elements
       npoiy=3*neley+1
       do i=1,nelex
          do j=1,neley
             ic=ic+1
             !ic=(i-1)*2*neley+j
             lnods( 1,ic)=(i-1)*3*npoiy        +3*(j-1)+1
             lnods( 2,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
             lnods( 3,ic)=(i-1)*3*npoiy+3*npoiy+3*j+1
             lnods( 4,ic)=(i-1)*3*npoiy+  npoiy+3*(j-1)+1
             lnods( 5,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+1
             lnods( 6,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+2
             lnods( 7,ic)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+3
             lnods( 8,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+3
             lnods( 9,ic)=(i-1)*3*npoiy+  npoiy+3*(j-1)+2
             lnods(10,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
             ic=ic+1
             !ic=(i-1)*2*neley+j+1
             lnods( 1,ic)=(i-1)*3*npoiy+        3*(j-1)+1
             lnods( 2,ic)=(i-1)*3*npoiy+3*npoiy+3*j+1
             lnods( 3,ic)=(i-1)*3*npoiy+        3*j+1
             lnods( 4,ic)=(i-1)*3*npoiy+  npoiy+3*(j-1)+2
             lnods( 5,ic)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+3
             lnods( 6,ic)=(i-1)*3*npoiy+2*npoiy+3*j+1
             lnods( 7,ic)=(i-1)*3*npoiy+  npoiy+3*j+1
             lnods( 8,ic)=(i-1)*3*npoiy+        3*(j-1)+3
             lnods( 9,ic)=(i-1)*3*npoiy+        3*(j-1)+2
             lnods(10,ic)=(i-1)*3*npoiy+  npoiy+3*(j-1)+3
          end do
       end do

       ! Correct corner elements that do not satisfy BB condition
       i=1
       j=neley
       ib=2*(i-1)*neley+2*(j-1)+1
       lnods( 1,ib)=(i-1)*3*npoiy+        3*(j-1)+1
       lnods( 2,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
       lnods( 3,ib)=(i-1)*3*npoiy+        3*j+1
       lnods( 4,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+1
       lnods( 5,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+1
       lnods( 6,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
       lnods( 7,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+3
       lnods( 8,ib)=(i-1)*3*npoiy+        3*(j-1)+3
       lnods( 9,ib)=(i-1)*3*npoiy+        3*(j-1)+2
       lnods(10,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+2
       ib=ib+1
       lnods( 1,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
       lnods( 2,ib)=(i-1)*3*npoiy+3*npoiy+3*j+1
       lnods( 3,ib)=(i-1)*3*npoiy+        3*j+1
       lnods( 4,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+2
       lnods( 5,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+3
       lnods( 6,ib)=(i-1)*3*npoiy+2*npoiy+3*j+1
       lnods( 7,ib)=(i-1)*3*npoiy+  npoiy+3*j+1
       lnods( 8,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+3
       lnods( 9,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
       lnods(10,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+3
       i=nelex
       j=1
       ib=2*(i-1)*neley+2*(j-1)+1
       lnods( 1,ib)=(i-1)*3*npoiy+        3*(j-1)+1
       lnods( 2,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
       lnods( 3,ib)=(i-1)*3*npoiy+        3*j+1
       lnods( 4,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+1
       lnods( 5,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+1
       lnods( 6,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
       lnods( 7,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+3
       lnods( 8,ib)=(i-1)*3*npoiy+        3*(j-1)+3
       lnods( 9,ib)=(i-1)*3*npoiy+        3*(j-1)+2
       lnods(10,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+2
       ib=ib+1
       lnods( 1,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+1
       lnods( 2,ib)=(i-1)*3*npoiy+3*npoiy+3*j+1
       lnods( 3,ib)=(i-1)*3*npoiy+        3*j+1
       lnods( 4,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+2
       lnods( 5,ib)=(i-1)*3*npoiy+3*npoiy+3*(j-1)+3
       lnods( 6,ib)=(i-1)*3*npoiy+2*npoiy+3*j+1
       lnods( 7,ib)=(i-1)*3*npoiy+  npoiy+3*j+1
       lnods( 8,ib)=(i-1)*3*npoiy+  npoiy+3*(j-1)+3
       lnods( 9,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+2
       lnods(10,ib)=(i-1)*3*npoiy+2*npoiy+3*(j-1)+3

    end if

    if(ic.ne.nelem) write(*,*) 'Elements: Wrong number of elements!'
    return

  end subroutine mesh_rectangular_elements_2d

!================================================================================================
  subroutine mesh_rectangular_elements_3d (nelex,neley,nelez,nnode,nelem,lnods)
    !subroutine mesh_rectangular_elements (nelex,neley,nnode,nelem,ln)
    implicit none
    integer(ip), intent(in)  :: nelex,neley,nelez,nnode,nelem
    integer(ip), intent(out) :: lnods(nnode,nelem)
    ! If there is some compiler problem, try
    !integer(ip) ln(nnode*nelem)
    !#define lnods(i,j) ln(((j)-1)*(nnode)+(i))

    integer(ip)              :: i,j,k,ic,ib,npoiy

    ic=0

    if(nnode==8) then       ! Q1 elements

       do i=1,nelex
          do j=1,neley
             do k=1,nelez
                ic=ic+1
                !ic=(i-1)*neley+j
                lnods(1,ic)=(i-1)*(neley+1)*(nelez+1)+(j-1)*(nelez+1)+k
                lnods(2,ic)=(i-1)*(neley+1)*(nelez+1)+(j-1)*(nelez+1)+k+1
                lnods(3,ic)=(i-1)*(neley+1)*(nelez+1)+j*(nelez+1)+k+1
                lnods(4,ic)=(i-1)*(neley+1)*(nelez+1)+j*(nelez+1)+k
                lnods(5,ic)=i*(neley+1)*(nelez+1)+(j-1)*(nelez+1)+k
                lnods(6,ic)=i*(neley+1)*(nelez+1)+(j-1)*(nelez+1)+k+1
                lnods(7,ic)=i*(neley+1)*(nelez+1)+j*(nelez+1)+k+1    
                lnods(8,ic)=i*(neley+1)*(nelez+1)+j*(nelez+1)+k
             end do
          end do
       end do
    end if

    if(ic.ne.nelem) write(*,*) 'Elements: Wrong number of elements!'
    return
    
  end subroutine mesh_rectangular_elements_3d

  !================================================================================================
  subroutine mesh_rectangular_nodes_2d(npoix,npoiy,ndime,npoin,ntdix,ntdiy, &
     &                              aleng,bleng,zx1,zx2,zy1,zy2,coord)
    implicit none
    integer(ip), intent(in)  :: ndime,npoin,ntdix,ntdiy,npoix,npoiy
    real(rp)   , intent(in)  :: aleng,bleng,zx1,zx2,zy1,zy2
    real(rp)   , intent(out) :: coord(ndime,npoin)
    integer(ip) :: ic,j,k
    real(rp)    :: a1,a2,a3,a4,b1,b2,b3,b4
    real(rp)    :: deltx,delty

    ic=0

    if(ntdix.eq.0.and.ntdiy.eq.0) then

       deltx=aleng/dble(npoix-1)
       delty=bleng/dble(npoiy-1)
       do k=1,npoix
          do j=1,npoiy
             ic=ic+1
             coord(1,ic)=(k-1)*deltx
             coord(2,ic)=(j-1)*delty
          end do
       end do

    else if(ntdix.eq.0.and.ntdiy.eq.1) then

       deltx=aleng/dble(npoix-1)
       call mesh_rectangular_nodes_cubic(0.0_rp,bleng/2.0_rp,zy1,zy2,npoiy/2,1.4_rp,b1,b2,b3,b4)
       do k=1,npoix
          ic=ic+1
          coord(1,ic)=(k-1)*deltx
          coord(2,ic)=0.0
          do j=1,npoiy/2
             ic=ic+1
             coord(1,ic)=(k-1)*deltx
             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
          end do
          do j=npoiy/2-1,1,-1
             ic=ic+1
             coord(1,ic)=(k-1)*deltx
             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
          end do
          ic=ic+1
          coord(1,ic)=(k-1)*deltx
          coord(2,ic)=bleng
       end do

    else if(ntdix.eq.1.and.ntdiy.eq.0) then

       call mesh_rectangular_nodes_cubic(0.0_rp,aleng/2.0_rp,zx1,zx2,npoix/2,1.4_rp,a1,a2,a3,a4)
       delty=bleng/dble(npoiy-1)
       do j=1,npoiy
          ic=ic+1
          coord(1,ic)=0.0
          coord(2,ic)=(j-1)*delty
       end do
       do k=1,npoix/2
          do j=1,npoiy
             ic=ic+1
             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
             coord(2,ic)=(j-1)*delty
          end do
       end do
       do k=npoix/2-1,1,-1
          do j=1,npoiy
             ic=ic+1
             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
             coord(2,ic)=(j-1)*delty
          end do
       end do
       do j=1,npoiy
          ic=ic+1
          coord(1,ic)=aleng
          coord(2,ic)=(j-1)*delty
       end do

    else if(ntdix.eq.1.and.ntdiy.eq.1) then

       call mesh_rectangular_nodes_cubic(0.0_rp,aleng/2.0_rp,zx1,zx2,npoix/2,1.4_rp,a1,a2,a3,a4)
       call mesh_rectangular_nodes_cubic(0.0_rp,bleng/2.0_rp,zy1,zy2,npoiy/2,1.4_rp,b1,b2,b3,b4)
       ic=ic+1
       coord(1,ic)=0.0
       coord(2,ic)=0.0
       do j=1,npoiy/2
          ic=ic+1
          coord(1,ic)=0.0
          coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
       end do
       do j=npoiy/2-1,1,-1
          ic=ic+1
          coord(1,ic)=0.0
          coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
       end do
       ic=ic+1
       coord(1,ic)=0.0
       coord(2,ic)=bleng
       do k=1,npoix/2
          ic=ic+1
          coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
          coord(2,ic)=0.0
          do j=1,npoiy/2
             ic=ic+1
             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
          end do
          do j=npoiy/2-1,1,-1
             ic=ic+1
             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
          end do
          ic=ic+1
          coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
          coord(2,ic)=bleng
       end do
       do k=npoix/2-1,1,-1
          ic=ic+1
          coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
          coord(2,ic)=0.0
          do j=1,npoiy/2
             ic=ic+1
             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
          end do
          do j=npoiy/2-1,1,-1
             ic=ic+1
             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
          end do
          ic=ic+1
          coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
          coord(2,ic)=bleng
       end do
       ic=ic+1
       coord(1,ic)=aleng
       coord(2,ic)=0.0
       do j=1,npoiy/2
          ic=ic+1
          coord(1,ic)=aleng
          coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
       end do
       do j=npoiy/2-1,1,-1
          ic=ic+1
          coord(1,ic)=aleng
          coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
       end do
       ic=ic+1
       coord(1,ic)=aleng
       coord(2,ic)=bleng

    end if

    return 

  end subroutine mesh_rectangular_nodes_2d

  !================================================================================================
  subroutine mesh_rectangular_nodes_3d(npoix,npoiy,npoiz,ndime,npoin,ntdix,ntdiy,ntdiz, &
     &                              aleng,bleng,cleng,zx1,zx2,zy1,zy2,coord)
    implicit none
    integer(ip), intent(in)  :: ndime,npoin,ntdix,ntdiy,ntdiz,npoix,npoiy,npoiz
    real(rp)   , intent(in)  :: aleng,bleng,cleng,zx1,zx2,zy1,zy2
    real(rp)   , intent(out) :: coord(ndime,npoin)
    integer(ip) :: ic,j,k,i
    real(rp)    :: a1,a2,a3,a4,b1,b2,b3,b4
    real(rp)    :: deltx,delty,deltz

    ic=0

    if(ntdix.eq.0.and.ntdiy.eq.0.and.ntdiz.eq.0) then

       deltx=aleng/dble(npoix-1)
       delty=bleng/dble(npoiy-1)
       deltz=cleng/dble(npoiz-1)
       do k=1,npoix
          do j=1,npoiy
             do i=1,npoiz
                ic=ic+1
                coord(1,ic)=(k-1)*deltx
                coord(2,ic)=(j-1)*delty
                coord(3,ic)=(i-1)*deltz
             end do
          end do
       end do

    else if(ntdix.eq.0.and.ntdiy.eq.1) then

!!$       deltx=aleng/dble(npoix-1)
!!$       call mesh_rectangular_nodes_cubic(0.0_rp,bleng/2.0_rp,zy1,zy2,npoiy/2,1.4_rp,b1,b2,b3,b4)
!!$       do k=1,npoix
!!$          ic=ic+1
!!$          coord(1,ic)=(k-1)*deltx
!!$          coord(2,ic)=0.0
!!$          do j=1,npoiy/2
!!$             ic=ic+1
!!$             coord(1,ic)=(k-1)*deltx
!!$             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
!!$          end do
!!$          do j=npoiy/2-1,1,-1
!!$             ic=ic+1
!!$             coord(1,ic)=(k-1)*deltx
!!$             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
!!$          end do
!!$          ic=ic+1
!!$          coord(1,ic)=(k-1)*deltx
!!$          coord(2,ic)=bleng
!!$       end do
!!$
!!$    else if(ntdix.eq.1.and.ntdiy.eq.0) then
!!$
!!$       call mesh_rectangular_nodes_cubic(0.0_rp,aleng/2.0_rp,zx1,zx2,npoix/2,1.4_rp,a1,a2,a3,a4)
!!$       delty=bleng/dble(npoiy-1)
!!$       do j=1,npoiy
!!$          ic=ic+1
!!$          coord(1,ic)=0.0
!!$          coord(2,ic)=(j-1)*delty
!!$       end do
!!$       do k=1,npoix/2
!!$          do j=1,npoiy
!!$             ic=ic+1
!!$             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
!!$             coord(2,ic)=(j-1)*delty
!!$          end do
!!$       end do
!!$       do k=npoix/2-1,1,-1
!!$          do j=1,npoiy
!!$             ic=ic+1
!!$             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
!!$             coord(2,ic)=(j-1)*delty
!!$          end do
!!$       end do
!!$       do j=1,npoiy
!!$          ic=ic+1
!!$          coord(1,ic)=aleng
!!$          coord(2,ic)=(j-1)*delty
!!$       end do
!!$
!!$    else if(ntdix.eq.1.and.ntdiy.eq.1) then
!!$
!!$       call mesh_rectangular_nodes_cubic(0.0_rp,aleng/2.0_rp,zx1,zx2,npoix/2,1.4_rp,a1,a2,a3,a4)
!!$       call mesh_rectangular_nodes_cubic(0.0_rp,bleng/2.0_rp,zy1,zy2,npoiy/2,1.4_rp,b1,b2,b3,b4)
!!$       ic=ic+1
!!$       coord(1,ic)=0.0
!!$       coord(2,ic)=0.0
!!$       do j=1,npoiy/2
!!$          ic=ic+1
!!$          coord(1,ic)=0.0
!!$          coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
!!$       end do
!!$       do j=npoiy/2-1,1,-1
!!$          ic=ic+1
!!$          coord(1,ic)=0.0
!!$          coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
!!$       end do
!!$       ic=ic+1
!!$       coord(1,ic)=0.0
!!$       coord(2,ic)=bleng
!!$       do k=1,npoix/2
!!$          ic=ic+1
!!$          coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
!!$          coord(2,ic)=0.0
!!$          do j=1,npoiy/2
!!$             ic=ic+1
!!$             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
!!$             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
!!$          end do
!!$          do j=npoiy/2-1,1,-1
!!$             ic=ic+1
!!$             coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
!!$             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
!!$          end do
!!$          ic=ic+1
!!$          coord(1,ic)=a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4
!!$          coord(2,ic)=bleng
!!$       end do
!!$       do k=npoix/2-1,1,-1
!!$          ic=ic+1
!!$          coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
!!$          coord(2,ic)=0.0
!!$          do j=1,npoiy/2
!!$             ic=ic+1
!!$             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
!!$             coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
!!$          end do
!!$          do j=npoiy/2-1,1,-1
!!$             ic=ic+1
!!$             coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
!!$             coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
!!$          end do
!!$          ic=ic+1
!!$          coord(1,ic)=aleng-(a1*(k-1)**3+a2*(k-1)**2+a3*(k-1)+a4)
!!$          coord(2,ic)=bleng
!!$       end do
!!$       ic=ic+1
!!$       coord(1,ic)=aleng
!!$       coord(2,ic)=0.0
!!$       do j=1,npoiy/2
!!$          ic=ic+1
!!$          coord(1,ic)=aleng
!!$          coord(2,ic)=b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4
!!$       end do
!!$       do j=npoiy/2-1,1,-1
!!$          ic=ic+1
!!$          coord(1,ic)=aleng
!!$          coord(2,ic)=bleng-(b1*(j-1)**3+b2*(j-1)**2+b3*(j-1)+b4)
!!$       end do
!!$       ic=ic+1
!!$       coord(1,ic)=aleng
!!$       coord(2,ic)=bleng

    end if

    return 

  end subroutine mesh_rectangular_nodes_3d

  !================================================================================================
  subroutine mesh_rectangular_nodes_cubic(h1,h2,z1,z2,N,p,a,b,c,d)
    !-------------------------------------------------------------------------
    !-------------------------------------------------------------------------
    implicit none
    integer(ip) :: n
    real(rp)    :: h1,h2,z1,z2,p,a,b,c,d
    real(rp)    :: r,q

    c=z1
    d=h1+z1
    r=dble(N-1)/dble(2*(N-1)*(N-p)-3*(N-p)**2)
    q=(c*dble(N-1)+d-h2)/dble((N-1)**3)
    b=(z2-c+3*dble((N-p)**2)*q)*r
    a=(-c*(N-1)-d-b*dble((N-1)**2)+h2)/dble((N-1)**3)

    return
  end subroutine mesh_rectangular_nodes_cubic


end module fem_mesh_gen

