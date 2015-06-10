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
module eltrm_gen_names
  use types
  implicit none

contains
  ! ---------------------------------------------------
  !         ELMVIS_GAL
  ! ---------------------------------------------------
  subroutine elmvis_gal(dvolt0,visac,deriv,ndime,nnode,elmat,work)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: visac,dvolt0
    real(rp),    intent(inout) :: elmat(nnode,nnode)
    real(rp),    intent(inout) :: work(:)
    ! Locals
    integer(ip)                :: inode,jnode

    work(1) = dvolt0*visac
    do inode=1,nnode
       do jnode=1,nnode
          elmat(inode,jnode) = work(1)*dot_product(deriv(:,inode),deriv(:,jnode)) &
               + elmat(inode,jnode)
       end do
    end do

  end subroutine elmvis_gal
  ! ---------------------------------------------------
  !         ELMVIS_GAL_SYM
  ! ---------------------------------------------------
  subroutine elmvis_gal_sym(dvolt0,visac,deriv,ndime,nnode,elmat,work)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: visac,dvolt0
    real(rp),    intent(inout) :: elmat(ndime,ndime,nnode,nnode)
    real(rp),    intent(inout) :: work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime,jdime

    ! Add only cross terms
    work(1) = dvolt0*visac
    do inode=1,nnode
       do jnode=1,nnode
          do idime=1,ndime
             do jdime=1,ndime
                elmat(idime,jdime,inode,jnode) = work(1)*deriv(jdime,inode)*deriv(idime,jnode) &
                     + elmat(idime,jdime,inode,jnode)
             end do
          end do
       end do
    end do

  end subroutine elmvis_gal_sym
  ! ---------------------------------------------------
  !         ELMVIS_DIV_GAL
  ! ---------------------------------------------------
  subroutine elmvis_div_gal(dvolu,visac,deriv,ndime,nnode,elmat,work)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: visac,dvolu
    real(rp),    intent(inout) :: elmat(ndime,ndime,nnode,nnode)
    real(rp),    intent(inout) :: work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime,jdime

    work(1) = dvolu*visac
    do inode=1,nnode
       do jnode=1,nnode
          do idime=1,ndime
             do jdime=1,ndime
                elmat(idime,jdime,inode,jnode) = work(1)*deriv(idime,jnode)*deriv(jdime,inode) &
               + elmat(idime,jdime,inode,jnode)
             end do
          end do
       end do
    end do

  end subroutine elmvis_div_gal
  ! ---------------------------------------------------
  !         ELMMSS_GAL
  ! ---------------------------------------------------
  subroutine elmmss_gal(dvolu,dtinv,shape,nnode,elmat,work)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: dtinv,dvolu
    real(rp),    intent(inout) :: elmat(nnode,nnode)
    real(rp),    intent(inout) :: work(:)
    ! Locals
    integer(ip)                :: inode,jnode

    work(1) = dvolu*dtinv
    do inode=1,nnode
       do jnode=1,nnode
          elmat(inode,jnode) = work(1)*shape(inode)*shape(jnode) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmmss_gal
  ! ---------------------------------------------------
  !         ELMDIV_STAB
  ! ---------------------------------------------------
  subroutine elmdiv_stab(tidiv,dvolu,deriv,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for ASGS & OSS for block U,V 
    !     tau*(div v, div u)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: deriv(ndime,nnode),tidiv,dvolu
    real(rp),    intent(inout) :: elmat(ndime,ndime,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime,jdime

    work(1) = tidiv*dvolu
    do jnode=1,nnode
       do inode=1,nnode
          do jdime=1,ndime
             do idime=1,ndime
                elmat(idime,jdime,inode,jnode) =  elmat(idime,jdime,inode,jnode) + &
                     work(1)*deriv(idime,inode)*deriv(jdime,jnode)
             end do
          end do
       end do
    end do

  end subroutine elmdiv_stab
  ! ---------------------------------------------------
  !         ELMBUV_GAL
  ! ---------------------------------------------------
  subroutine elmbuv_gal(dvolu,porac,dtinv,shape,agran,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,V 
    !   (v, a·grad u) + s*(v,u) + (v, u/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: shape(nnode),agran(nnode)
    real(rp),    intent(in)    :: porac,dvolu,dtinv
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode
    real(rp)                   :: elmat_aux(nnode,nnode)

    do jnode=1,nnode
       work(1) = dvolu*(shape(jnode)*(porac + dtinv) + agran(jnode))
       do inode=1,nnode
          elmat(inode,jnode) = shape(inode)*work(1) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmbuv_gal

! ---------------------------------------------------
  !         ELMBUV_GAL
  ! ---------------------------------------------------
  subroutine elmbvu_gal(dvolu,shape,agran,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for block U,V 
    !   (u, a·grad v)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: shape(nnode),agran(nnode)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode
    real(rp)                   :: elmat_aux(nnode,nnode)

    do jnode=1,nnode
       work(1) = dvolu*shape(jnode)
       do inode=1,nnode
          elmat(inode,jnode) = agran(inode)*work(1) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmbvu_gal

  ! ---------------------------------------------------
  !         ELMBUV_GAL_SKEWSYMMETRIC(1)
  ! ---------------------------------------------------
  subroutine elmbuv_gal_skew1(dvolu,porac,dtinv,shape,agran,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,V 
    !   1/2(v, a·grad u) - 1/2(u,a·grad v) + s*(v,u) + (v, u/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: shape(nnode),agran(nnode)
    real(rp),    intent(in)    :: porac,dvolu,dtinv
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode
    real(rp)                   :: elmat_aux(nnode,nnode)

    do jnode=1,nnode
       work(1) = dvolu*shape(jnode)*(porac + dtinv)
       work(2) = dvolu*agran(jnode)*0.5
       do inode=1,nnode
          elmat_aux(inode,jnode) = shape(inode)*work(2)
          elmat(inode,jnode) = shape(inode)*work(1) + elmat_aux(inode,jnode) + elmat(inode,jnode)
       end do
    end do
    do jnode=1,nnode
       do inode=1,nnode
          elmat(inode,jnode) = - elmat_aux(jnode,inode) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmbuv_gal_skew1
  ! ---------------------------------------------------
  !         ELMBUV_GAL_SKEWYMMETRIC(2)
  ! ---------------------------------------------------
  subroutine elmbuv_gal_skew2(dvolu,shape,grvel,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the antisymmetric part of the convective term 
    !   1/2(v·u, div u)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: dvolu,grvel(ndime,ndime)
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime
    real(rp)                   :: divel

    divel = 0.0_rp
    do idime = 1,ndime
       divel = divel + grvel(idime,idime)
    end do

    work(1)=0.5*dvolu*divel
    do jnode=1,nnode
       work(2) = work(1)*shape(jnode)
       do inode=1,nnode
          elmat(inode,jnode) = work(2)*shape(inode) + &
               elmat(inode,jnode)
       end do
    end do

  end subroutine elmbuv_gal_skew2
  ! ---------------------------------------------------
  !         ELMBPV_GAL_DIV
  ! ---------------------------------------------------
  subroutine elmbpv_gal_div(dvolu,shape,deriv,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block P,V 
    !     - (div v, p)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: dvolu,deriv(ndime,nnode)
    real(rp),    intent(inout) :: elmat(ndime,1,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnode
       work(1) = shape(jnode)*dvolu
       do inode=1,nnode
          do idime=1,ndime
             elmat(idime,1,inode,jnode) = -deriv(idime,inode)*work(1) + &
                  elmat(idime,1,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbpv_gal_div
  ! ---------------------------------------------------
  !         ELMBPV_GAL_DIV_INF-SUP_STABLE
  ! ---------------------------------------------------
  subroutine elmbpv_gal_div_iss(dvolu,shape,deriv,ndime,nnodu,nnodp,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block P,V 
    !     - (div v, p)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnodu,nnodp,ndime
    real(rp),    intent(in)    :: shape(nnodp)
    real(rp),    intent(in)    :: dvolu,deriv(ndime,nnodu)
    real(rp),    intent(inout) :: elmat(ndime,1,nnodu,nnodp),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnodp
       work(1) = shape(jnode)*dvolu
       do inode=1,nnodu
          do idime=1,ndime
             elmat(idime,1,inode,jnode) = -deriv(idime,inode)*work(1) + &
                  elmat(idime,1,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbpv_gal_div_iss
  ! ---------------------------------------------------
  !         ELMBPV_GAL_GRD
  ! ---------------------------------------------------
  subroutine elmbpv_gal_grd(dvolu,shape,deriv,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block P,V 
    !     + (v, grad p)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: dvolu,deriv(ndime,nnode)
    real(rp),    intent(inout) :: elmat(ndime,1,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnode
       do inode=1,nnode
          work(1) = shape(inode)*dvolu
          do idime=1,ndime
             elmat(idime,1,inode,jnode) = deriv(idime,jnode)*work(1) + &
                  elmat(idime,1,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbpv_gal_grd
  ! ---------------------------------------------------
  !         ELMBPV_GAL_GRD_INF-SUP_STABLE
  ! ---------------------------------------------------
  subroutine elmbpv_gal_grd_iss(dvolu,shape,deriv,ndime,nnodu,nnodp,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block P,V 
    !     + (v, grad p)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnodu,nnodp,ndime
    real(rp),    intent(in)    :: shape(nnodu)
    real(rp),    intent(in)    :: dvolu,deriv(ndime,nnodp)
    real(rp),    intent(inout) :: elmat(ndime,1,nnodu,nnodp),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnodp
       do inode=1,nnodu
          work(1) = shape(inode)*dvolu
          do idime=1,ndime
             elmat(idime,1,inode,jnode) = deriv(idime,jnode)*work(1) + &
                  elmat(idime,1,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbpv_gal_grd_iss
  ! ---------------------------------------------------
  !         ELMBUQ_GAL_DIV
  ! ---------------------------------------------------
  subroutine elmbuq_gal_div(dvolu,shape,deriv,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,Q
    !     (div u, q)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(1,ndime,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    
    do inode=1,nnode
       work(1) = shape(inode)*dvolu
       do jnode=1,nnode
          do idime=1,ndime
             elmat(1,idime,inode,jnode) = deriv(idime,jnode)*work(1) + elmat(1,idime,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbuq_gal_div
  ! ---------------------------------------------------
  !         ELMBUQ_GAL_DIV_INF-SUP_STABLE
  ! ---------------------------------------------------
  subroutine elmbuq_gal_div_iss(dvolu,shape,deriv,ndime,nnodu,nnodp,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,Q
    !     (div u, q)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnodu,nnodp,ndime
    real(rp),    intent(in)    :: shape(nnodp)
    real(rp),    intent(in)    :: deriv(ndime,nnodu)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(1,ndime,nnodp,nnodu),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do inode=1,nnodp
       work(1) = shape(inode)*dvolu
       do jnode=1,nnodu
          do idime=1,ndime
             elmat(1,idime,inode,jnode) = deriv(idime,jnode)*work(1) + elmat(1,idime,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbuq_gal_div_iss
  ! ---------------------------------------------------
  !         ELMBUQ_GAL_GRD
  ! ---------------------------------------------------
  subroutine elmbuq_gal_grd(dvolu,shape,deriv,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,Q
    !     - (u, grad q)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(1,ndime,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnode
       work(1) = shape(jnode)*dvolu
       do inode=1,nnode
          do idime=1,ndime
             elmat(1,idime,inode,jnode) = - deriv(idime,inode)*work(1) + elmat(1,idime,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbuq_gal_grd
  ! ---------------------------------------------------
  !         ELMBUQ_GAL_GRD_INF-SUP_STABLE
  ! ---------------------------------------------------
  subroutine elmbuq_gal_grd_iss(dvolu,shape,deriv,ndime,nnodu,nnodp,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for Galerkin for block U,Q
    !     - (u, grad q)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnodu,nnodp,ndime
    real(rp),    intent(in)    :: shape(nnodu)
    real(rp),    intent(in)    :: deriv(ndime,nnodp)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(1,ndime,nnodp,nnodu),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnodu
       work(1) = shape(jnode)*dvolu
       do inode=1,nnodp
          do idime=1,ndime
             elmat(1,idime,inode,jnode) = - deriv(idime,inode)*work(1) + elmat(1,idime,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbuq_gal_grd_iss
  ! ---------------------------------------------------
  !         ELMRHU_GAL
  ! ---------------------------------------------------
  subroutine elmrhu_gal(dvolu,dtinv,shape,gpveln,elext,nnode,ndime,elrhs,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for Galerkin
    !    (v, f) + (v, u_n/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: elext(ndime),gpveln(ndime)
    real(rp),    intent(in)    :: dvolu,dtinv
    real(rp),    intent(inout) :: elrhs(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime

    do inode=1,nnode
       work(1) = shape(inode)*dvolu
       do idime=1,ndime
          elrhs(idime,inode) = (elext(idime) + gpveln(idime)*dtinv)*work(1) + &
               elrhs(idime,inode)
       end do
    end do

  end subroutine elmrhu_gal
  ! ---------------------------------------------------
  !         ELMRHS_FCE
  ! ---------------------------------------------------
  subroutine elmrhs_fce(dvolu,shape,elext,nnode,ndime,elrhs,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for Galerkin
    !    (v, f)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode),elext(ndime)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elrhs(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime

    do inode=1,nnode
       work(1) = shape(inode)*dvolu
       do idime=1,ndime
          elrhs(idime,inode) = elext(idime)*work(1) + elrhs(idime,inode)
       end do
    end do

  end subroutine elmrhs_fce
  ! ---------------------------------------------------
  !         OSS_CONVU
  ! ---------------------------------------------------
  subroutine oss_convu_arg_chk(dvolu,vegau,grave,shape,ndime,nnode,elrhs,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs term for projection OSS 
    !     (u · grad u, v) 
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode),dvolu
    real(rp),    intent(in)    :: vegau(ndime),grave(ndime,ndime)
    real(rp),    intent(inout) :: elrhs(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime,jdime

    work=0.0_rp
    do idime = 1,ndime
       do jdime = 1,ndime
          work(idime) = work(idime) + vegau(jdime)*grave(jdime,idime) 
       end do
    end do

    do inode=1,nnode
       work(4) = shape(inode)*dvolu
       do idime=1,ndime
          elrhs(idime,inode) = work(4)*work(idime) + elrhs(idime,inode)
       end do
    end do

  end subroutine oss_convu_arg_chk
  ! ---------------------------------------------------
  !         OSS_GRADP
  ! ---------------------------------------------------
  subroutine oss_gradp_arg_chk(dvolu,grapr,shape,ndime,nnode,elrhs,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs term for projection OSS 
    !     (grad p, v) 
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode),dvolu
    real(rp),    intent(in)    :: grapr(ndime)
    real(rp),    intent(inout) :: elrhs(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime

    do inode=1,nnode
       work(1) = shape(inode)*dvolu
       do idime=1,ndime
          elrhs(idime,inode) = work(1)*grapr(idime) + elrhs(idime,inode)
       end do
    end do

  end subroutine oss_gradp_arg_chk
  ! ---------------------------------------------------
  !         RESIDUAL
  ! ---------------------------------------------------
  subroutine gal_resid_arg_chk(dvolu,shape,ndime,ndofn,nnode,acvis,acpor,dtnsi,veadv,vegau,grapr, &
       &                       grave,elext,resim,work,ntens,hsvel)
    !-----------------------------------------------------------------------
    !
    !    This routine computes the residual of the finite element component
    !    
    !    IMPORTANT NOTE: R = (Lu - f), => positive sign on the RHS terms
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip),        intent(in) :: ndime,ndofn,nnode,ntens
    real(rp),           intent(in) :: dtnsi,acvis,acpor,dvolu
    real(rp),           intent(in) :: veadv(ndime),vegau(ndime,2)
    real(rp),           intent(in) :: grapr(ndime)
    real(rp),           intent(in) :: grave(ndime,ndime)
    real(rp),           intent(in) :: elext(ndime)
    real(rp),           intent(in) :: shape(nnode)
    real(rp),          intent(out) :: resim(ndofn)
    real(rp),        intent(inout) :: work(:)
    real(rp), optional, intent(in) :: hsvel(ntens,ndime)

    integer(ip)                :: idime,jdime,inode,idofn,itens

    resim = 0.0_rp

    resim(1:ndime) = resim(1:ndime) + grapr     ! Contribution from pressure term

    do idime=1,ndime                            ! Contribution from temporal derivative
       resim(idime) = resim(idime) + dtnsi*(vegau(idime,1)-vegau(idime,2))
    end do

    if(present(hsvel)) then                     ! Contribution from viscous term
        do itens=1,ndime
            do idime=1,ndime
                resim(idime) = resim(idime) - acvis*hsvel(itens,idime)
            end do
        end do
    end if

    do idime = 1,ndime                          ! Contribution from the convective term
       do jdime = 1,ndime
          resim(idime) = resim(idime) + veadv(jdime)*grave(jdime,idime) 
       end do
    end do

    resim(1:ndime) = resim(1:ndime) + acpor*vegau(:,1) ! Contribution from the porosity term

    resim(1:ndime) = resim(1:ndime) - elext         ! Contribution from the external forces

    do idime = 1,ndime                              ! Contribution from the divergence term
       resim(ndime+1) = resim(ndime+1) + grave(idime,idime)
    end do

  end subroutine gal_resid_arg_chk
  ! ---------------------------------------------------
  !         ELMRHS_RES
  ! ---------------------------------------------------
  subroutine elmrhs_res(ndofn,nnode,dvolu,resim,shape,elrhs,work)
    !-----------------------------------------------------------------------
    !
    !    This routine computes the RHS of the residual of the finite 
    !    element component
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: nnode,ndofn
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(in)    :: resim(ndofn),shape(nnode)
    real(rp),    intent(inout) :: elrhs(ndofn,nnode),work(:)

    integer(ip)                :: inode,idofn

    ! Compute only rhsid
    do inode=1,nnode
       work(1) = shape(inode)*dvolu
       do idofn=1,ndofn
          elrhs(idofn,inode) = elrhs(idofn,inode) + work(1)*resim(idofn)
       end do
    end do

  end subroutine elmrhs_res
  ! ---------------------------------------------------
  !         ELMBUV_ASGS
  ! ---------------------------------------------------
  subroutine elmbuv_asgs(dvolu,dtinv,shape,testf,operu,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for ASGS for block U,V 
    !   - tau1*(L*v, Lu) - tau1*(L*v, u/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: shape(nnode)
    real(rp),    intent(in)    :: operu(nnode),testf(nnode)
    real(rp),    intent(in)    :: dvolu,dtinv
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode

    do jnode=1,nnode
       work(1) = (operu(jnode) + shape(jnode)*dtinv)*dvolu
       do inode=1,nnode
          elmat(inode,jnode) = testf(inode)*work(1) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmbuv_asgs
  ! ---------------------------------------------------
  !         ELMBUV_OSS
  ! ---------------------------------------------------
  subroutine elmbuv_oss(dvolu,testf,operu,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for OSS for block U,V 
    !   - tau1*(L*v, Lu)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode
    real(rp),    intent(in)    :: operu(nnode),testf(nnode)
    real(rp),    intent(in)    :: dvolu
    real(rp),    intent(inout) :: elmat(nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode

    do jnode=1,nnode
       work(1) = operu(jnode)*dvolu
       do inode=1,nnode
          elmat(inode,jnode) = testf(inode)*work(1) + elmat(inode,jnode)
       end do
    end do

  end subroutine elmbuv_oss
  ! ---------------------------------------------------
  !         ELMBUQ_ASGS
  ! ---------------------------------------------------
  subroutine elmbuq_asgs(timom,dvolu,dtinv,shape,deriv,operu,ndime,nnode,elmat,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs terms for ASGS for block U,Q
    !     tau1*(grad q, Lu) + tau1*(grad q, u/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: shape(nnode),operu(nnode)
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: dvolu,timom,dtinv
    real(rp),    intent(inout) :: elmat(1,ndime,nnode,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,jnode,idime

    do jnode=1,nnode
       work(1) = (operu(jnode)+dtinv*shape(jnode))*timom*dvolu
       do inode=1,nnode
          do idime=1,ndime
             elmat(1,idime,inode,jnode) = &
                  deriv(idime,inode)*work(1) + elmat(1,idime,inode,jnode)
          end do
       end do
    end do

  end subroutine elmbuq_asgs
  ! ---------------------------------------------------
  !         ELMRHP_ASGS
  ! ---------------------------------------------------
  subroutine elmrhp_asgs(timom,dvolu,dtinv,deriv,gpveln,elext,nnode,ndime,elrhs,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the rhs terms for ASGS
    !    tau*(f, grad q) + tau*(grad q, u_n/dt)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: elext(ndime),gpveln(ndime)
    real(rp),    intent(in)    :: timom,dvolu,dtinv
    real(rp),    intent(inout) :: elrhs(nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime

    work(1) = dvolu*timom
    work(2) = work(1)*dtinv
    do inode=1,nnode
       do idime=1,ndime
          elrhs(inode) = deriv(idime,inode)*(elext(idime)*work(1) + &
               gpveln(idime)*work(2)) + elrhs(inode)
       end do
    end do

  end subroutine elmrhp_asgs
  ! ---------------------------------------------------
  !         ELMRHU_GRADU_GRADV
  ! ---------------------------------------------------
  subroutine elmrhu_gradgrad(dvolu,diffu,grvel,deriv,ndime,nnode,elvec,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs term grad-grad for block V
    !     nu*(grad u, grad v)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: grvel(ndime,ndime)
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: dvolu,diffu
    real(rp),    intent(inout) :: elvec(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime,jdime

    work(1) = dvolu*diffu
    do inode=1,nnode
       do idime=1,ndime
          do jdime=1,ndime
             elvec(idime,inode) =  grvel(idime,jdime)*deriv(jdime,inode)*work(1) + elvec(idime,inode)
          end do
       end do
    end do

  end subroutine elmrhu_gradgrad
  ! ---------------------------------------------------
  !         ELMRHU_GRADU_GRADV_SYM
  ! ---------------------------------------------------
  subroutine elmrhu_gradgrad_sym(dvolu,diffu,grvel,deriv,ndime,nnode,elvec,work)
    !-----------------------------------------------------------------------
    !
    ! This routine computes the lhs term grad-grad for block V
    !     nu*(grad u, grad v)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: grvel(ndime,ndime)
    real(rp),    intent(in)    :: deriv(ndime,nnode)
    real(rp),    intent(in)    :: dvolu,diffu
    real(rp),    intent(inout) :: elvec(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime,jdime

    work(1) = dvolu*diffu
    do inode=1,nnode
       do idime=1,ndime
          do jdime=1,ndime
             elvec(idime,inode) =  grvel(jdime,idime)*deriv(jdime,inode)*work(1) + elvec(idime,inode)
          end do
       end do
    end do

  end subroutine elmrhu_gradgrad_sym
  ! ---------------------------------------------------
  !         ELMRHU_DIVV_P
  ! ---------------------------------------------------
  subroutine elmrhu_pdivv(dvolu,gppre,deriv,ndime,nnode,elvec,work)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the lhs term p-div(v) block V 
    !     - (div v, p)
    !
    !-----------------------------------------------------------------------
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nnode,ndime
    real(rp),    intent(in)    :: gppre
    real(rp),    intent(in)    :: dvolu,deriv(ndime,nnode)
    real(rp),    intent(inout) :: elvec(ndime,nnode),work(:)
    ! Locals
    integer(ip)                :: inode,idime

    work(1) = gppre*dvolu
    do inode=1,nnode
       do idime=1,ndime
          elvec(idime,inode) = - deriv(idime,inode)*work(1) + &
               elvec(idime,inode)
       end do
    end do

  end subroutine elmrhu_pdivv

end module eltrm_gen_names
