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
module bomap_interp
  use types
  use memor
  use bomap_names
  use interpolation_names
  implicit none
# include "debug.i90"
  private

  ! Functions
  public :: bomap_from_interp

contains

  !==============================================================================
  !
  !==============================================================================
  subroutine bomap_from_interp(int,bocod,map)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    type(interpolation), intent(in)    :: int
    real(rp)           , intent(in)    :: bocod(:,:)
    type(bomap)        , intent(inout) :: map
    ! Locals
    integer(ip) :: ndime,ndimb,nnodb,nlocs
    integer(ip) :: ilocs, inode, idime

    ! Check and get data from int
    assert(int%kder==1)
    ndimb=size(int%deriv,dim=1)
    nnodb=size(int%deriv,dim=2)
    nlocs=size(int%deriv,dim=3)
    ndime=ndimb+1

    ! Check elcod
    assert(ndime==size(bocod,dim=1))
    assert(nnodb==size(bocod,dim=2))

    ! Check map allocation (baloc)
    assert(ndime==size(map%baloc,dim=1))
    assert(nlocs==size(map%baloc,dim=3))

    do ilocs = 1,nlocs
       map%clocs(:,ilocs)=0.0_rp
       do inode=1,nnodb
          do idime=1,ndime
             map%clocs(idime,ilocs) = map%clocs(idime,ilocs) &
                  + bocod(idime,inode)*int%shape(inode,ilocs)
          end do
       end do
       ! Evaluate tangent vectors
       call mbmabt(map%baloc(:,:,ilocs),bocod,int%deriv,ndime,ndimb,nnodb)
       ! map%baloc(:,1:ndimb,ilocs)=matmul(bocod(:,:),transpose(derib(:,:,ilocs)))

       ! Compute normal vector (last one in baloc)
       if(ndime==2) then
          map%baloc(1,2,ilocs) =  map%baloc(2,1,ilocs)
          map%baloc(2,2,ilocs) = -map%baloc(1,1,ilocs)
       else
          call vecpro(map%baloc(1,1,ilocs),map%baloc(1,2,ilocs),map%baloc(1,3,ilocs),3)
       end if

       ! Compute Euclidian norm of the normal vector
       call vecno2(map%baloc(1,ndime,ilocs),ndime,map%detjm(ilocs))

    end do

  end subroutine bomap_from_interp

  !-----------------------------------------------------------------------
  subroutine mbmabt(a,b,c,n1,n2,n3)

    !-----------------------------------------------------------------------
    !
    ! This routine evaluates the matrix product A = B Ct, where
    ! A -> Mat(n1,n2), B -> Mat(n1,n3), C -> Mat(n2,n3)
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip) :: n1,n2,n3,i,j,k
    real(rp)    :: a(n1,n2), b(n1,n3), c(n2,n3)

    do i=1,n1
       do j=1,n2
          a(i,j)=0.0_rp
          do k=1,n3
             a(i,j)=a(i,j)+b(i,k)*c(j,k)
          end do
       end do
    end do

  end subroutine mbmabt

  subroutine vecpro(v1,v2,v3,n)
    !-----------------------------------------------------------------------
    !
    ! Two and three-dimensional vectorial product of two vectors  v3 = v1 x v2.
    ! The same pointer as for v1 or v2 may be used for v3. If N = 2, it is
    !  assumed that v1 = (0,0,v1_3) and v2 = (v2_1,v2_2,0).      
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v2(n),v1(3)
    real(rp),    intent(out) :: v3(n)
    real(rp)                 :: c1,c2,c3

    if(n==2) then
       c1=-v1(3)*v2(2)
       c2= v1(3)*v2(1)
       v3(1)=c1
       v3(2)=c2
    else if(n==3) then
       c1=v1(2)*v2(3)-v1(3)*v2(2)
       c2=v1(3)*v2(1)-v1(1)*v2(3)
       c3=v1(1)*v2(2)-v1(2)*v2(1)
       v3(1)=c1
       v3(2)=c2
       v3(3)=c3
    end if

  end subroutine vecpro

  subroutine vecno2(v,n,vnor)
    implicit none
    integer(ip), intent(in)  :: n
    real(rp),    intent(in)  :: v(n)
    real(rp),    intent(out) :: vnor
    integer(ip)              :: i

    vnor=0.0_rp
    do i=1,n
       vnor=vnor+v(i)*v(i)
    end do
    vnor=sqrt(vnor)

  end subroutine vecno2

end module bomap_interp
