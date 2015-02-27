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
module fem_conditions_names
  use types
  use memor
  use stdio
  use fem_mesh_names
  use renum_names
# include "debug.i90"
  implicit none
  private

  type fem_conditions
     integer(ip)                :: &
          ncode=1,                 &         ! Number of codes  (=nvars)
          nvalu=1,                 &         ! Number of values (=nvars)
          ncond=1                            ! Number of conds  (=npoin)
     integer(ip), allocatable   :: &
          code(:,:)                          ! Codes
     real(rp), allocatable      :: &
          valu(:,:)                          ! Values
  end type fem_conditions

  ! Types
  public :: fem_conditions

  ! Methods
  public :: fem_conditions_create, fem_conditions_box, fem_conditions_free, fem_conditions_copy, fem_conditions_apply_renum

contains

  !===============================================================================================
  subroutine fem_conditions_create(ncode,nvalu,ncond,cnd)
    implicit none
    integer(ip)         , intent(in)  :: ncode,nvalu,ncond
    type(fem_conditions), intent(out) :: cnd

    cnd%ncode=ncode
    cnd%nvalu=nvalu
    cnd%ncond=ncond

    call memalloc (cnd%ncode,cnd%ncond,cnd%code, __FILE__,__LINE__)
    call memalloc (cnd%nvalu,cnd%ncond,cnd%valu, __FILE__,__LINE__)
    cnd%code=0
    cnd%valu=0.0_rp

  end subroutine fem_conditions_create

  !===============================================================================================
  subroutine fem_conditions_copy(cnd_old,cnd_new)
    implicit none
    type(fem_conditions), intent(in)  :: cnd_old
    type(fem_conditions), intent(out) :: cnd_new

    call fem_conditions_create( cnd_old%ncode, cnd_old%nvalu, cnd_old%ncond, cnd_new)
    cnd_new%code=cnd_old%code
    cnd_new%valu=cnd_old%code

  end subroutine fem_conditions_copy

  !===============================================================================================
  subroutine fem_conditions_apply_renum(cnd, renumeration)
    implicit none
    type(fem_conditions), intent(inout)  :: cnd
    type(renum),              intent(in) :: renumeration
    
    integer(ip)                          :: tmp_int1(cnd%ncond),tmp_int2(cnd%ncond)
    real(rp)                             :: tmp_real1(cnd%ncond), tmp_real2(cnd%ncond)
    integer(ip)                          :: idx

    ! todo: do it without 2 temporary arrays.. if I try it, there is a warning that they were created it anyway
    do idx = 1, cnd%ncode
       tmp_int1 = cnd%code(idx, :)
       call renum_apply(renumeration, tmp_int1, tmp_int2)
       cnd%code(idx, :) = tmp_int2
    end do

    do idx = 1, cnd%nvalu
       tmp_real1 = cnd%valu(idx, :)
       call renum_apply(renumeration, tmp_real1, tmp_real2)
       cnd%valu(idx, :) = tmp_real2
    end do
  end subroutine fem_conditions_apply_renum

  !===============================================================================================
  subroutine fem_conditions_box(msh,poin,line,surf,node)
    !------------------------------------------------------------------------
    !
    ! This routine applies boundary conditions to a structured mesh in a box.
    !
    !------------------------------------------------------------------------
    implicit none
    type(fem_mesh)       , intent(in)  :: msh
    type(fem_conditions) , intent(in)  :: poin,line,surf
    type(fem_conditions) , intent(out) :: node
    integer(ip) :: npoix,npoiy,npoiz,nelex,neley,nelez

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
       if(msh%isper(1)==1) npoix=npoix-1
       if(msh%isper(2)==1) npoiy=npoiy-1

       assert(poin%ncode==line%ncode)
       assert(poin%nvalu==line%nvalu)
       assert(poin%ncode==surf%ncode)
       assert(poin%nvalu==surf%nvalu)
       assert(poin%ncond==4)
       assert(line%ncond==4)
       assert(surf%ncond==1)
       call fem_conditions_create(poin%ncode,poin%nvalu,msh%npoin,node)

       call mesh_rectangular_conditions_2d(npoix,npoiy,poin%ncode,poin%nvalu,poin%code,poin%valu,  &
            &                           line%code,line%valu,surf%code,surf%valu,node%code, &
            &                           node%valu)

    else if(msh%ndime==3) then
       
       nelex=msh%nedir(1)
       neley=msh%nedir(2)
       nelez=msh%nedir(3)

       if(msh%nnode==8) then
          npoiy=neley+1
          npoix=nelex+1
          npoiz=nelez+1
       end if

       if(msh%isper(1)==1) npoix=npoix-1
       if(msh%isper(2)==1) npoiy=npoiy-1
       if(msh%isper(3)==1) npoiz=npoiz-1

       assert(poin%ncode==line%ncode)
       assert(poin%nvalu==line%nvalu)
       assert(poin%ncode==surf%ncode)
       assert(poin%nvalu==surf%nvalu)
       assert(poin%ncond==8)
       assert(line%ncond==12)
       assert(surf%ncond==6)
       call fem_conditions_create(poin%ncode,poin%nvalu,msh%npoin,node)

       call mesh_rectangular_conditions_3d(npoix,npoiy,npoiz,poin%ncode,poin%nvalu,poin%code,  &
            &                           poin%valu,line%code,line%valu,surf%code,surf%valu,node%code, &
            &                           node%valu)

    end if

  end subroutine fem_conditions_box

  !===============================================================================================
  subroutine mesh_rectangular_conditions_2d(npoix,npoiy,ncode,nvalu,codep,valup, &
       &                                 codel,valul,codes,valus,coden,valun)
    implicit none
    integer(ip), intent(in)  :: npoix,npoiy
    integer(ip), intent(in)  :: ncode,nvalu
    integer(ip), intent(in)  :: codep(ncode,4)
    integer(ip), intent(in)  :: codel(ncode,4)
    integer(ip), intent(in)  :: codes(ncode,1)
    integer(ip), intent(out) :: coden(ncode,npoix*npoiy)
    real(rp)   , intent(in)  :: valup(nvalu,4)
    real(rp)   , intent(in)  :: valul(nvalu,4)
    real(rp)   , intent(in)  :: valus(nvalu,1)
    real(rp)   , intent(out) :: valun(nvalu,npoix*npoiy)
    integer(ip) :: ic,j,i


    ! Boundary conditions
    ic=0

    ! Left bottom corner
    i=1
    j=1
    ic=ic+1
    coden(:,ic)=codep(:,1)
    valun(:,ic)=valup(:,1)

    ! Left line
    do j=2,npoiy-1
       ic=ic+1
       coden(:,ic)=codel(:,4)
       valun(:,ic)=valul(:,4)
    end do

    ! Left top corner
    j=npoiy
    ic=ic+1
    coden(:,ic)=codep(:,4)
    valun(:,ic)=valup(:,4)

    do i=2,npoix-1
       ! Bottom line
       j=1
       ic=ic+1
       coden(:,ic)=codel(:,1)
       valun(:,ic)=valul(:,1)
       ! Interior
       do j=2,npoiy-1
          ic=ic+1
          coden(:,ic)=codes(:,1)
          valun(:,ic)=valus(:,1)
       end do
       ! Top line
       j=npoiy
       ic=ic+1
       coden(:,ic)=codel(:,3)
       valun(:,ic)=valul(:,3)
    end do

    ! Right bottom corner
    i=npoix
    j=1
    ic=ic+1
    coden(:,ic)=codep(:,2)
    valun(:,ic)=valup(:,2)
    ! Right line
    do j=2,npoiy-1
       ic=ic+1
       coden(:,ic)=codel(:,2)
       valun(:,ic)=valul(:,2)
    end do
    ! Right top corner
    j=npoiy
    ic=ic+1
    coden(:,ic)=codep(:,3)
    valun(:,ic)=valup(:,3)

  end subroutine mesh_rectangular_conditions_2d

  !===============================================================================================
  subroutine mesh_rectangular_conditions_3d(npoix,npoiy,npoiz,ncode,nvalu,codep,valup, &
       &                                 codel,valul,codes,valus,coden,valun)
    implicit none
    integer(ip), intent(in)  :: npoix,npoiy,npoiz
    integer(ip), intent(in)  :: ncode,nvalu
    integer(ip), intent(in)  :: codep(ncode,8)
    integer(ip), intent(in)  :: codel(ncode,12)
    integer(ip), intent(in)  :: codes(ncode,6)
    integer(ip), intent(out) :: coden(ncode,npoix*npoiy*npoiz)
    real(rp)   , intent(in)  :: valup(nvalu,8)
    real(rp)   , intent(in)  :: valul(nvalu,12)
    real(rp)   , intent(in)  :: valus(nvalu,6)
    real(rp)   , intent(out) :: valun(nvalu,npoix*npoiy*npoiz)
    integer(ip) :: ic,j,i,k

    ! 3D:
    ! Points code:  1=(0,0,0), 2=(0,0,1), 3=(0,1,0), 4=(0,1,1), 5=(1,0,0),
    !               6=(1,0,1), 7=(1,1,0), 8=(1,1,1)
    ! Lines code:   1=(0,0,z), 2=(0,1,z), 3=(1,0,z), 4=(1,1,z), 5=(0,y,0),
    !               6=(0,y,1), 7=(1,y,0), 8=(1,y,1), 9=(x,0,0), 10=(x,0,1),
    !               11=(x,1,0), 12=(x,1,1)
    ! Surface code: 1=(x,y,0), 2=(x,y,1), 3=(x,0,z), 4=(x,1,z), 5=(0,y,z),
    !               6=(1,y,z)

    ! Boundary conditions
    ic=0

    i=1
    j=1
    ! Point 1
    k=1
    ic=ic+1
    coden(:,ic)=codep(:,1)
    valun(:,ic)=valup(:,1)
    ! Line 1
    do k=2,npoiz-1
       ic=ic+1
       coden(:,ic)=codel(:,1)
       valun(:,ic)=valul(:,1)
    end do
    ! Point 2
    k=npoiz
    ic=ic+1
    coden(:,ic)=codep(:,2)
    valun(:,ic)=valup(:,2)
       
    do j=2,npoiy-1
       ! Line 5
       k=1
       ic=ic+1
       coden(:,ic)=codel(:,5)
       valun(:,ic)=valul(:,5)
       ! Surface 5
       do k=2,npoiz-1
          ic=ic+1
          coden(:,ic)=codes(:,5)
          valun(:,ic)=valus(:,5)
       end do   
       ! Line 6
       k=npoiz
       ic=ic+1
       coden(:,ic)=codel(:,6)
       valun(:,ic)=valul(:,6)
    end do
    
    j=npoiy
    ! Point 3
    k=1
    ic=ic+1
    coden(:,ic)=codep(:,3)
    valun(:,ic)=valup(:,3)
    ! Line 2
    do k=2,npoiz-1
       ic=ic+1
       coden(:,ic)=codel(:,2)
       valun(:,ic)=valul(:,2)
    end do
    ! Point 4
    k=npoiz
    ic=ic+1
    coden(:,ic)=codep(:,4)
    valun(:,ic)=valup(:,4)

    do i=2,npoix-1
       j=1
       ! Line 9
       k=1
       ic=ic+1
       coden(:,ic)=codel(:,9)
       valun(:,ic)=valul(:,9)
       ! Surface 3
       do k=2,npoiz-1
          ic=ic+1
          coden(:,ic)=codes(:,3)
          valun(:,ic)=valus(:,3)
       end do
       ! Line 10
       k=npoiz
       ic=ic+1
       coden(:,ic)=codel(:,10)
       valun(:,ic)=valul(:,10)
       
       do j=2,npoiy-1
          ! Surface 1
          k=1
          ic=ic+1
          coden(:,ic)=codes(:,1)
          valun(:,ic)=valus(:,1)
          ! Interior nodes
          do k=2,npoiz-1
             ic=ic+1
          end do
          ! Surface 2
          k=npoiz
          ic=ic+1
          coden(:,ic)=codes(:,2)
          valun(:,ic)=valus(:,2)
       end do
       
       j=npoiy
       ! Line 11
       k=1
       ic=ic+1
       coden(:,ic)=codel(:,11)
       valun(:,ic)=valul(:,11)
       ! Surface 4
       do k=2,npoiz-1
          ic=ic+1
          coden(:,ic)=codes(:,4)
          valun(:,ic)=valus(:,4)
       end do
       ! Line 12
       k=npoiz
       ic=ic+1
       coden(:,ic)=codel(:,12)
       valun(:,ic)=valul(:,12)
    end do

    i=npoix
    j=1
    ! Point 5
    k=1
    ic=ic+1
    coden(:,ic)=codep(:,5)
    valun(:,ic)=valup(:,5)
    ! Line 3
    do k=2,npoiz-1
       ic=ic+1
       coden(:,ic)=codel(:,3)
       valun(:,ic)=valul(:,3)
    end do
    ! Point 6
    k=npoiz
    ic=ic+1
    coden(:,ic)=codep(:,6)
    valun(:,ic)=valup(:,6)

    do j=2,npoiy-1
       ! Line 7
       k=1
       ic=ic+1
       coden(:,ic)=codel(:,7)
       valun(:,ic)=valul(:,7)
       ! Surface 6
       do k=2,npoiz-1
          ic=ic+1
          coden(:,ic)=codes(:,6)
          valun(:,ic)=valus(:,6)
       end do
       ! Line 8
       k=npoiz
       ic=ic+1
       coden(:,ic)=codel(:,8)
       valun(:,ic)=valul(:,8)
    end do

    j=npoiy
    ! Point 7
    k=1
    ic=ic+1
    coden(:,ic)=codep(:,7)
    valun(:,ic)=valup(:,7)
    ! Line 4
    do k=2,npoiz-1
       ic=ic+1
       coden(:,ic)=codel(:,4)
       valun(:,ic)=valul(:,4)
    end do
    ! Point 8
    k=npoiz
    ic=ic+1
    coden(:,ic)=codep(:,8)
    valun(:,ic)=valup(:,8)
       
  end subroutine mesh_rectangular_conditions_3d

  !===============================================================================================
  subroutine fem_conditions_free(cnd)
    implicit none
    type(fem_conditions), intent(inout) :: cnd

    cnd%ncode=-1
    cnd%nvalu=-1
    cnd%ncond=-1

    call memfree (cnd%code,__FILE__,__LINE__)
    call memfree (cnd%valu,__FILE__,__LINE__)

    return

  end subroutine fem_conditions_free

end module fem_conditions_names
