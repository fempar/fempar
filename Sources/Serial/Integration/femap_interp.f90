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
module femap_interp_names
use types_names
use memor_names
  use femap_names
  use interpolation_names
  use face_interpolation_names
  implicit none
# include "debug.i90"
  private
  real(rp), parameter   :: tol = 1.e-8_rp

  ! Functions
  public :: femap_from_interp, femap_apply_to_interp,femap_face_to_interp, femap_from_face_interp

contains

  !==============================================================================
  !
  !==============================================================================
  subroutine femap_apply_to_interp(map,ref,phy)
    implicit none
    type(interpolation_t), intent(in)    :: ref
    type(femap_t)        , intent(in)    :: map
    type(interpolation_t), intent(inout) :: phy
    real(rp), allocatable :: wmat1(:,:,:)
    real(rp), allocatable :: wmat2(:,:,:)
    integer(ip) :: ndime,nnode,nlocs,ntens
    integer(ip) :: ilocs,idime,jdime,kdime,ldime,inode

    assert(ref%nlocs==phy%nlocs)
    assert(ref%nlocs==map%nlocs)

    ! Shape functions do not change
    phy%shape = ref%shape

    ! First derivatives do
    if(phy%kder==1) then
       phy%deriv=0.0_rp
       do ilocs=1,phy%nlocs
          do inode=1,phy%nnode
             do idime=1,phy%ndime
                do jdime=1,ref%ndime
                   phy%deriv(idime,inode,ilocs) = phy%deriv(idime,inode,ilocs) &
                      + map%jainv(jdime,idime,ilocs)*ref%deriv(jdime,inode,ilocs)
                end do
             end do
          end do
       end do
    end if

    ! Second derivatives are
    !
    !    d^2 N / d x_i d x_j
    !       = (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j)
    !       + (d N / d s_k) (d^2 s_k / d x_i d x_j) 
    !
    if(phy%khes==1) then

       call memalloc(phy%ndime,phy%ndime,phy%nnode,wmat1,__FILE__,__LINE__)
       call memalloc(phy%ndime,phy%ndime,phy%nnode,wmat2,__FILE__,__LINE__)

       do ilocs=1,phy%nlocs

          if(ref%khes==1) then
             ! Transforms the array HESSI to a symmetric matrix WMAT1
             do inode=1,ref%nnode
                call vetoma(ref%hessi(1,inode,ilocs),wmat1(1,1,inode),ref%ndime,ref%ntens)
             end do
             ! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for each node
             do inode=1,ref%nnode
                call btdbma(wmat2(1,1,inode),wmat1(1,1,inode),map%jainv(:,:,ilocs), &
                   &        map%ndime,map%ndime)
             end do
          end if

          if(map%khes==1) then
             ! Add second cartesian derivatives of the map times 
             ! first derivatives of shape functions
             do inode=1,ref%nnode
                do idime=1,map%ndime
                   do jdime=1,map%ndime
                      do kdime=1,map%ndime
                         wmat2(idime,jdime,inode)=wmat2(idime,jdime,inode) &
                            & + ref%deriv(kdime,inode,ilocs) &
                            &   * map%d2sdx(kdime,idime,jdime,ilocs)
                      end do
                   end do
                end do
             end do
          end if

          ! Writes the Hessian matrix as an array
          do inode=1,phy%nnode
             call matove(wmat2(1,1,inode),phy%hessi(1,inode,ilocs),phy%ndime,phy%ntens)
          end do

       end do

       call memfree(wmat1,__FILE__,__LINE__)
       call memfree(wmat2,__FILE__,__LINE__)

    end if

  end subroutine femap_apply_to_interp


 !==============================================================================
  subroutine femap_face_to_interp(map,ref,fac,phy)
    implicit none
    type(face_interpolation_t), intent(in)    :: ref
    type(femap_t)             , intent(in)    :: map
    integer(ip)             , intent(in)    :: fac
    type(face_interpolation_t), intent(inout) :: phy
 
    real(rp)    :: norm2
    integer(ip) :: ndime,nnode,nlocs
    integer(ip) :: ilocs,idime,jdime,inode

    assert(ref%ngaus==phy%ngaus)

    ! Shape functions do not change
    phy%shapf(:,1,:) = ref%shapf(:,fac,:)

    ! First derivatives on the boundary
    phy%derif=0.0_rp
    do ilocs=1,phy%ngaus
       do inode=1,phy%nnode
          do idime=1,phy%ndime
             do jdime=1,ref%ndime
                phy%derif(idime,inode,1,ilocs) = phy%derif(idime,inode,1,ilocs) &
                     + map%jainv(jdime,idime,ilocs)*ref%derif(jdime,inode,fac,ilocs)
             end do
          end do
       end do
    end do

    ! Outside Normal
    phy%outno=0.0_rp
    do ilocs=1,phy%ngaus
       norm2 = 0.0_rp
       do idime=1,phy%ndime
          do jdime=1,ref%ndime
             phy%outno(idime,1,ilocs) = phy%outno(idime,1,ilocs) &
                    + map%jainv(jdime,idime,ilocs)*ref%outno(jdime,fac,ilocs)
          end do
          norm2 = norm2 + phy%outno(idime,1,ilocs)**2
       end do
       phy%outno(:,:,ilocs) = phy%outno(:,:,ilocs)/sqrt(norm2)
    end do  

  end subroutine femap_face_to_interp

  !==============================================================================
  subroutine femap_from_interp(int,elcod,map)
    !-----------------------------------------------------------------------
    ! A map obtained from the (usually isoparametric) interpolation of the geometry
    !-----------------------------------------------------------------------
    implicit none
    type(interpolation_t), intent(in)    :: int
!    real(rp)           , intent(inout) :: hnatu
    real(rp)           , intent(in)    :: elcod(:,:)
    type(femap_t)        , intent(inout) :: map
    ! Locals
    real(rp), allocatable :: wmat1(:,:,:)
    real(rp), allocatable :: wmat2(:,:,:)
    real(rp)    :: hnatu
    real(rp)    :: enor0,h_tem
    integer(ip) :: ndime,nnode,nlocs,ntens
    integer(ip) :: ilocs,idime,jdime,kdime,ldime,inode

    ! Check and get data from int
    assert(int%kder==1)
    ndime=size(int%deriv,dim=1)
    nnode=size(int%deriv,dim=2)
    nlocs=size(int%deriv,dim=3)

    ! Check elcod
    assert(ndime==size(elcod,dim=1))
    assert(nnode==size(elcod,dim=2))

    ! Check map allocation (jainv, detjm and hleng
    ! are assumed to be allocated if jacob is)
    assert(ndime==size(map%jacob,dim=1))
    assert(nlocs==size(map%jacob,dim=3))

    ! Jacobian and its inverse
    do ilocs=1,nlocs
       ! Matmul is not thread safe
       !map%jacob(:,:,ilocs)=matmul(elcod,transpose(int%deriv(:,:,ilocs)))
       map%jacob(:,:,ilocs)=0.0_rp
       do inode=1,nnode
          do jdime=1,ndime
             do idime=1,ndime
                map%jacob(idime,jdime,ilocs) = map%jacob(idime,jdime,ilocs) &
                     + elcod(idime,inode)*int%deriv(jdime,inode,ilocs)
             end do
          end do
       end do
       ! J^(-t)
       call invmtx(map%jacob(:,:,ilocs),map%jainv(:,:,ilocs),map%detjm(ilocs),ndime)
    end do

    ! Evaluation (Gauss) point coordinates
    do ilocs=1,nlocs
       map%clocs(:,ilocs)=0.0_rp
       do inode=1,nnode
          do idime=1,ndime
                map%clocs(idime,ilocs) = map%clocs(idime,ilocs) &
                     + elcod(idime,inode)*int%shape(inode,ilocs)
          end do
       end do
    end do

    ! Natural element length (de hecho hay que dividir por p, no?)
    if(ndime==2) then
       if(nnode==3) then ! P1
          hnatu=1.0_rp
       else if(nnode==4) then ! Q1
          hnatu=2.0_rp
       else if(nnode==6) then ! P2
          hnatu=1.0_rp
       else if(nnode==9) then ! Q2
          hnatu=2.0_rp
       else if(nnode==10) then ! P3
          hnatu=1.0_rp
       else if(nnode==16) then ! Q3
          hnatu=2.0_rp
       end if
    else
       if(nnode==4) then ! P1
          hnatu=1.0_rp
       else if(nnode==8) then ! Q1
          hnatu=2.0_rp
       else if(nnode==10) then ! P2
          hnatu=1.0_rp
       else if(nnode==27) then ! Q2
          hnatu=2.0_rp
       else if(nnode==20) then ! P3
          hnatu=1.0_rp
       else if(nnode==64) then ! Q3
          hnatu=2.0_rp
       end if
    end if

    do ilocs=1,nlocs
       ! Element length
       do idime=1,ndime
          enor0=0.0_rp
          do jdime=1,ndime
             enor0=enor0+map%jainv(idime,jdime,ilocs)**2
          end do
          enor0=sqrt(enor0)
          map%hleng(idime,ilocs)=hnatu/enor0
       end do

       ! Sorth hleng: hleng(1)=max; hleng(ndime)=min
       do idime=1,ndime-1                   
          do jdime=idime+1,ndime
             if(map%hleng(jdime,ilocs)>map%hleng(idime,ilocs)) then
                h_tem       =map%hleng(jdime,ilocs)
                map%hleng(jdime,ilocs)=map%hleng(idime,ilocs)
                map%hleng(idime,ilocs)=h_tem
             end if
          end do
       end do
    end do

    ! Second derivatives of the map
    if(map%khes>0) then

       assert(int%khes==1)
       ntens=size(int%hessi,dim=1)
       ! Check that second derivativesof the map have been allocated.
       assert(ndime==size(map%d2sdx,dim=1))
       assert(nlocs==size(map%d2sdx,dim=4))

       call memalloc(ndime,ndime,nnode,wmat1,__FILE__,__LINE__)
       call memalloc(ndime,ndime,nnode,wmat2,__FILE__,__LINE__)

       do ilocs=1,nlocs

          ! Transforms the array HESSI to a symmetric matrix WMAT1
          do inode=1,nnode
             call vetoma(int%hessi(1,inode,ilocs),wmat1(1,1,inode),ndime,ntens)
          end do

          ! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for
          ! each node
          do inode=1,nnode
             call btdbma(wmat2(1,1,inode),wmat1(1,1,inode), &
                  &        map%jainv(:,:,ilocs),ndime,ndime)
          end do

          ! Obtains (d^2 s_k / d x_i d x_j) as the solution of the system
          ! (d x_l / d s_k) (d^2 s_k / d x_i d x_j) 
          !     = - (d^2 x_l / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j), 
          ! for l,i,j = 1,...,NDIME
          do kdime=1,ndime
             do idime=1,ndime
                do jdime=1,ndime
                   map%d2sdx(kdime,idime,jdime,ilocs)=0.0_rp
                   do ldime=1,ndime
                      do inode=1,nnode
                         map%d2sdx(kdime,idime,jdime,ilocs) =    &
                              & map%d2sdx(kdime,idime,jdime,ilocs) &
                              & - map%jainv(kdime,ldime,ilocs)     &
                              &   * wmat2(idime,jdime,inode) * elcod(ldime,inode)
                      end do
                   end do
                end do
             end do
          end do

       end do

       call memfree(wmat1,__FILE__,__LINE__)
       call memfree(wmat2,__FILE__,__LINE__)

    end if


  end subroutine femap_from_interp

 !==============================================================================
  subroutine femap_from_face_interp(int,elcod,iface,map)
    !-----------------------------------------------------------------------
    ! A map obtained from the (usually isoparametric) interpolation of the geometry
    !-----------------------------------------------------------------------
    implicit none
    type(face_interpolation_t), intent(in)    :: int
    real(rp)                , intent(in)    :: elcod(:,:)
    integer(ip)             , intent(in)    :: iface
    type(femap_t)             , intent(inout) :: map

    ! Locals
    real(rp)    :: hnatu
    real(rp)    :: enor0,h_tem
    integer(ip) :: ndime,nnode,nlocs,ntens
    integer(ip) :: ilocs,idime,jdime,kdime,ldime,inode

    ! Check and get data from int
    ndime=size(int%derif,dim=1)
    nnode=size(int%derif,dim=2)
    nlocs=size(int%derif,dim=4)

    ! Check elcod
    assert(ndime==size(elcod,dim=1))
    assert(nnode==size(elcod,dim=2))

    ! Check map allocation (jainv, detjm and hleng
    ! are assumed to be allocated if jacob is)
    assert(ndime==size(map%jacob,dim=1))
    assert(nlocs==size(map%jacob,dim=3))

    ! Jacobian and its inverse
    do ilocs=1,nlocs
       ! Matmul is not thread safe
       !map%jacob(:,:,ilocs)=matmul(elcod,transpose(int%deriv(:,:,ilocs)))
       map%jacob(:,:,ilocs)=0.0_rp
       do inode=1,nnode
          do jdime=1,ndime
             do idime=1,ndime
                map%jacob(idime,jdime,ilocs) = map%jacob(idime,jdime,ilocs) &
                     + elcod(idime,inode)*int%derif(jdime,inode,iface,ilocs)
             end do
          end do
       end do
       ! J^(-t)
       call invmtx(map%jacob(:,:,ilocs),map%jainv(:,:,ilocs),map%detjm(ilocs),ndime)
    end do

    ! Evaluation (Gauss) point coordinates
    do ilocs=1,nlocs
       map%clocs(:,ilocs)=0.0_rp
       do inode=1,nnode
          do idime=1,ndime
                map%clocs(idime,ilocs) = map%clocs(idime,ilocs) &
                     + elcod(idime,inode)*int%shapf(inode,iface,ilocs)
          end do
       end do
    end do

    ! Natural element length (de hecho hay que dividir por p, no?)
    if(ndime==2) then
       if(nnode==3) then ! P1
          hnatu=1.0_rp
       else if(nnode==4) then ! Q1
          hnatu=2.0_rp
       else if(nnode==6) then ! P2
          hnatu=1.0_rp
       else if(nnode==9) then ! Q2
          hnatu=2.0_rp
       else if(nnode==10) then ! P3
          hnatu=1.0_rp
       else if(nnode==16) then ! Q3
          hnatu=2.0_rp
       end if
    else
       if(nnode==4) then ! P1
          hnatu=1.0_rp
       else if(nnode==8) then ! Q1
          hnatu=2.0_rp
       else if(nnode==10) then ! P2
          hnatu=1.0_rp
       else if(nnode==27) then ! Q2
          hnatu=2.0_rp
       else if(nnode==20) then ! P3
          hnatu=1.0_rp
       else if(nnode==64) then ! Q3
          hnatu=2.0_rp
       end if
    end if

    do ilocs=1,nlocs
       ! Element length
       do idime=1,ndime
          enor0=0.0_rp
          do jdime=1,ndime
             enor0=enor0+map%jainv(idime,jdime,ilocs)**2
          end do
          enor0=sqrt(enor0)
          map%hleng(idime,ilocs)=hnatu/enor0
       end do

       ! Sorth hleng: hleng(1)=max; hleng(ndime)=min
       do idime=1,ndime-1                   
          do jdime=idime+1,ndime
             if(map%hleng(jdime,ilocs)>map%hleng(idime,ilocs)) then
                h_tem       =map%hleng(jdime,ilocs)
                map%hleng(jdime,ilocs)=map%hleng(idime,ilocs)
                map%hleng(idime,ilocs)=h_tem
             end if
          end do
       end do
    end do

!!$    ! Second derivatives of the map
!!$    if(map%khes>0) then
!!$
!!$       assert(int%khes==1)
!!$       ntens=size(int%hessi,dim=1)
!!$       ! Check that second derivativesof the map have been allocated.
!!$       assert(ndime==size(map%d2sdx,dim=1))
!!$       assert(nlocs==size(map%d2sdx,dim=4))
!!$
!!$       call memalloc(ndime,ndime,nnode,wmat1,__FILE__,__LINE__)
!!$       call memalloc(ndime,ndime,nnode,wmat2,__FILE__,__LINE__)
!!$
!!$       do ilocs=1,nlocs
!!$
!!$          ! Transforms the array HESSI to a symmetric matrix WMAT1
!!$          do inode=1,nnode
!!$             call vetoma(int%hessi(1,inode,ilocs),wmat1(1,1,inode),ndime,ntens)
!!$          end do
!!$
!!$          ! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for
!!$          ! each node
!!$          do inode=1,nnode
!!$             call btdbma(wmat2(1,1,inode),wmat1(1,1,inode), &
!!$                  &        map%jainv(:,:,ilocs),ndime,ndime)
!!$          end do
!!$
!!$          ! Obtains (d^2 s_k / d x_i d x_j) as the solution of the system
!!$          ! (d x_l / d s_k) (d^2 s_k / d x_i d x_j) 
!!$          !     = - (d^2 x_l / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j), 
!!$          ! for l,i,j = 1,...,NDIME
!!$          do kdime=1,ndime
!!$             do idime=1,ndime
!!$                do jdime=1,ndime
!!$                   map%d2sdx(kdime,idime,jdime,ilocs)=0.0_rp
!!$                   do ldime=1,ndime
!!$                      do inode=1,nnode
!!$                         map%d2sdx(kdime,idime,jdime,ilocs) =    &
!!$                              & map%d2sdx(kdime,idime,jdime,ilocs) &
!!$                              & - map%jainv(kdime,ldime,ilocs)     &
!!$                              &   * wmat2(idime,jdime,inode) * elcod(ldime,inode)
!!$                      end do
!!$                   end do
!!$                end do
!!$             end do
!!$          end do
!!$
!!$       end do
!!$
!!$       call memfree(wmat1,__FILE__,__LINE__)
!!$       call memfree(wmat2,__FILE__,__LINE__)
!!$
!!$    end if
!!$
!!$
  end subroutine femap_from_face_interp

  !-----------------------------------------------------------------------
  subroutine matove(xmatr,vecto,ndime,ntens)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine stores a symmetric matrix XMATR into a vector VECTO
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ntens
    real(rp)   , intent(in)  :: xmatr(ndime,ndime)
    real(rp)   , intent(out) :: vecto(ntens)

    if(ndime.eq.2) then
       vecto(1)=xmatr(1,1)
       vecto(3)=xmatr(1,2)
       vecto(2)=xmatr(2,2)
    else
       vecto(1)=xmatr(1,1)
       vecto(4)=xmatr(1,2)
       vecto(2)=xmatr(2,2)
       vecto(5)=xmatr(1,3)
       vecto(6)=xmatr(2,3)
       vecto(3)=xmatr(3,3)
    end if

  end subroutine matove

  !-----------------------------------------------------------------------
  subroutine vetoma(vecto,xmatr,ndime,ntens)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine stores a vector VECTO as a symmetric matrix XMATR
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ntens
    real(rp)   , intent(in)  :: vecto(ntens)
    real(rp)   , intent(out) :: xmatr(ndime,ndime)

    if(ndime.eq.2) then
       xmatr(1,1)=vecto(1)
       xmatr(1,2)=vecto(3)
       xmatr(2,1)=vecto(3)
       xmatr(2,2)=vecto(2)
    else
       xmatr(1,1)=vecto(1)
       xmatr(1,2)=vecto(4)
       xmatr(1,3)=vecto(5)
       xmatr(2,1)=vecto(4)
       xmatr(2,2)=vecto(2)
       xmatr(2,3)=vecto(6)
       xmatr(3,1)=vecto(5)
       xmatr(3,2)=vecto(6)
       xmatr(3,3)=vecto(3)
    end if

  end subroutine vetoma

  !-----------------------------------------------------------------------
  subroutine btdbma(aglob,aloca,bmatr,n1,n2)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine computes Ag = Bt Al B  when Ag and Al are stored as full
    ! matrices (Ag := aglob, Al := aloca, B := bmatr). The dimensions are
    ! Al -> Mat(n1,n1), Ag -> Mat(n2,n2), B -> Mat(n2,n1) 
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: n1,n2
    real(rp)   , intent(in)  :: aloca(n1,n1), bmatr(n1,n2)
    real(rp)   , intent(out) :: aglob(n2,n2)
    integer(ip)              :: i,j,k,l

    do i=1,n2
       do j=1,n2
          aglob(i,j)=0.0
          do k=1,n1
             do l=1,n1
                aglob(i,j)=aglob(i,j)+bmatr(k,i)*aloca(k,l)*bmatr(l,j)
             end do
          end do
       end do
    end do

  end subroutine btdbma


  subroutine invmtx(a,b,deter,nsize)
    !-----------------------------------------------------------------------
    !
    ! This routine inverts a square matrix A -> Mat(nsize,nsize). The
    ! inverse is stored in B. Its determinant is DETER
    !    
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nsize
    real(rp),    intent(in)  :: a(nsize,nsize)
    real(rp),    intent(out) :: b(nsize,nsize),deter
    integer(ip)              :: isize,jsize
    real(rp)                 :: denom,t1,t2,t3,t4

    select case (nsize)

    case(1)
       deter=a(1,1)
       if(deter==0.0_rp) return
       b(1,1) = 1.0_rp/a(1,1)

    case(2)
       deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
       if(deter/=0.0_rp) then
          denom=1.0_rp/deter
          b(1,1) = a(2,2)*denom
          b(2,2) = a(1,1)*denom
          b(2,1) =-a(2,1)*denom
          b(1,2) =-a(1,2)*denom 
       end if

    case(3)
       t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
       t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
       t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
       deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
       if(deter==0.0_rp) return
       denom = 1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
       b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
       b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
       b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
       b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
       b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

    case(4)
       t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
       t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
       t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
       t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
       deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
       if(deter==0.0_rp) return
       denom=1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(4,1) = t4*denom
       b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
          - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
          + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
       b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
          + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
          - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
       b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
          - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
          + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
       b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
          + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
          - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
       b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
          + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
          - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
       b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
          - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
          + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
       b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
          + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
          - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
       b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
          - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
          + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
       b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
          - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
          + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
       b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
          + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
          - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
       b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
          - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
          + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
       b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
          + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
          - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

    case default
       b=a
       call invert(b,nsize,nsize)

    end select

  end subroutine invmtx

  !-----------------------------------------------------------------------
  subroutine invert(a,nmax,ndm)
    !-----------------------------------------------------------------------
    !
    ! This routine performs the inversion of a ndm*ndm square matrix 
    ! or just part of it (nmax*nmax)
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: ndm,nmax
    real(rp),    intent(inout) :: a(ndm,ndm)
    real(rp)                   :: d
    integer(ip)                :: n,j,i

    do n = 1,nmax
       d = a(n,n)
       do j = 1,nmax
          a(n,j) = -a(n,j)/d
       end do
       do i = 1,nmax
          if(n/=i) then
             do j = 1,nmax
                if(n/=j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
             end do
          end if
          a(i,n) = a(i,n)/d
       end do
       a(n,n) = 1.0_rp/d
    end do

  end subroutine invert

end module femap_interp_names
