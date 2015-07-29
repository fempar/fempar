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
module element_gather_tools_names
use types_names
  implicit none
  private

  interface gather
     module procedure gather_scalar, gather_vector
  end interface gather

  ! Functions
  ! public :: gather,interp_gauss,interp_grad,interp_hess,elem2gauss,elem2gauss_grad, &
  !      &    elem2gauss_hess
   public :: interpolate,interpolate_grad,interpolate_hess,elem2var, &
        &    elmuv2elmat,elmuv2elmat_lump

contains
  !==============================================================================
  ! 
  !==============================================================================
  subroutine gather_scalar (nnode,lnods,glval,elval)
    !-----------------------------------------------------------------------
    ! This routine performs a gather of an array
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nnode,lnods(nnode)
    real(rp),    intent(in)  :: glval(*)
    real(rp),    intent(out) :: elval(nnode)
    integer(ip)              :: inode

    do inode=1,nnode
       elval(inode) = glval(lnods(inode))
    end do

  end subroutine gather_scalar

  subroutine gather_vector (ncomp,nnode,lnods,glval,elval)
    !-----------------------------------------------------------------------
    ! This routine performs a gather of an array
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,lnods(nnode)
    real(rp),    intent(in)  :: glval(ncomp,*)
    real(rp),    intent(out) :: elval(ncomp,nnode)
    integer(ip)              :: inode

    do inode=1,nnode
       elval(1:ncomp,inode) = glval(1:ncomp,lnods(inode))
    end do

  end subroutine gather_vector

  subroutine interpolate (ncomp,nnode,ngaus,shape,elval,gpval)
    !-----------------------------------------------------------------------
    ! This routine computes the interpolation at gauss points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus
    real(rp),    intent(in)  :: shape(nnode,ngaus)
    real(rp),    intent(in)  :: elval(ncomp,nnode)
    real(rp),    intent(out) :: gpval(ncomp,ngaus)
    integer(ip)              :: inode,igaus

    gpval=0.0_rp
    do igaus=1,ngaus
       do inode=1,nnode
          gpval(1:ncomp,igaus) = gpval(1:ncomp,igaus) &
               &               + shape(inode,igaus)*elval(1:ncomp,inode)
       end do
    end do

  end subroutine interpolate

  subroutine interpolate_grad (ncomp,nnode,ngaus,cartd,elval,gpval)
    !-----------------------------------------------------------------------
    ! This routine computes the interpolation of gradients at gauss points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus
    real(rp),    intent(in)  :: cartd(ncomp,nnode,ngaus)
    real(rp),    intent(in)  :: elval(ncomp,nnode)
    real(rp),    intent(out) :: gpval(ncomp,ncomp,ngaus)
    integer(ip)              :: inode,igaus,idime,jdime

    do igaus=1,ngaus
       do idime=1,ncomp
          do jdime=1,ncomp
             gpval(idime,jdime,igaus)=0.0_rp
             do inode=1,nnode
                gpval(idime,jdime,igaus) = gpval(idime,jdime,igaus) + &
                     &                     cartd(idime,inode,igaus)*elval(jdime,inode)
             end do
          end do
       end do
    end do
    
  end subroutine interpolate_grad

  subroutine interpolate_hess (ncomp,nnode,ntens,ngaus,hessi,elval,gphes)
    !-----------------------------------------------------------------------
    ! This routine computes the interpolation of the Hessian at gauss points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus,ntens
    real(rp),    intent(in)  :: hessi(ntens,nnode,ngaus)
    real(rp),    intent(in)  :: elval(ncomp,ngaus)
    real(rp),    intent(out) :: gphes(ntens,ncomp,ngaus)
    integer(ip)              :: inode,igaus,idime,itens

    do igaus=1,ngaus
       do idime=1,ncomp
          do itens=1,ntens
             gphes(itens,idime,igaus)=0.0_rp
             do inode=1,nnode
                gphes(itens,idime,igaus) = gphes(itens,idime,igaus) + &
                     &                     hessi(itens,inode,igaus)*elval(idime,inode)
             end do
          end do
       end do
    end do
    
  end subroutine interpolate_hess

  subroutine interp_gauss (ncomp,nnode,shape,elval,gpval)
    !-----------------------------------------------------------------------
    ! This routine computes the interpolation in a given gauss point
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode
    real(rp),    intent(in)  :: shape(nnode)
    real(rp),    intent(in)  :: elval(ncomp,nnode)
    real(rp),    intent(out) :: gpval(ncomp)
    integer(ip)              :: inode

    gpval=0.0_rp
    do inode=1,nnode
       gpval(1:ncomp) = gpval(1:ncomp) + shape(inode)*elval(1:ncomp,inode)
    end do

  end subroutine interp_gauss

  subroutine interp_grad(nnode,ndime,cartd,elvar,gpvar)
    !------------------------------------------------------------------------
    ! This routine computes the interpolation of gradients
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in) :: nnode,ndime
    real(rp), intent(in)    :: cartd(ndime,nnode),elvar(ndime,nnode)
    real(rp), intent(inout) :: gpvar(ndime,ndime)
    integer(ip)             :: inode,idime,jdime

    do idime=1,ndime
       do jdime=1,ndime
          gpvar(idime,jdime)=0.0_rp
          do inode=1,nnode
             gpvar(idime,jdime) = gpvar(idime,jdime) + cartd(idime,inode)*elvar(jdime,inode)
          end do
       end do
    end do

  end subroutine interp_grad

  subroutine interp_hess(nnode,ndime,ntens,hessi,elvar,gpvar)
    !------------------------------------------------------------------------
    ! This routine computes the interpolation of gradients
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in) :: nnode,ndime,ntens
    real(rp),    intent(in) :: hessi(ntens,nnode),elvar(ndime,nnode)
    real(rp), intent(inout) :: gpvar(ntens,ndime)
    integer(ip)             :: inode,idime,jdime,itens

    do idime=1,ndime
       do itens=1,ntens
          gpvar(itens,idime)=0.0_rp
          do inode=1,nnode
             gpvar(itens,idime) = gpvar(itens,idime) + hessi(itens,inode)*elvar(idime,inode)
          end do
       end do
    end do

  end subroutine interp_hess

  subroutine elem2gauss(ncomp,nnode,ngaus,shape,lnods,glval,gpval)
    !------------------------------------------------------------------------
    ! This routine performs a gather and computes the interpolation at gauss points
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus,lnods(nnode)
    real(rp),    intent(in)  :: glval(ncomp,*)
    real(rp),    intent(in)  :: shape(nnode,ngaus)
    real(rp),    intent(out) :: gpval(ncomp,ngaus)

    real(rp)    :: elval(ncomp,nnode)
    integer(ip) :: igaus

    call gather(ncomp,nnode,lnods,glval,elval)
    !call interpolate(ncomp,nnode,ngaus,shape,elval,gpval)
    do igaus=1,ngaus
       call interp_gauss(ncomp,nnode,shape(:,igaus),elval,gpval(:,igaus))
    end do

  end subroutine elem2gauss

  subroutine elem2gauss_grad(ncomp,nnode,ngaus,cartd,lnods,glval,gpval)
    !------------------------------------------------------------------------
    ! This routine performs a gather and computes the interpolation of gradients
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus,lnods(nnode)
    real(rp),    intent(in)  :: glval(ncomp,*)
    real(rp),    intent(in)  :: cartd(ncomp,nnode,ngaus)
    real(rp),    intent(out) :: gpval(ncomp,ncomp,ngaus)

    integer(ip) :: igaus
    real(rp)    :: elval(ncomp,nnode)

    call gather(ncomp,nnode,lnods,glval,elval)
    do igaus=1,ngaus
       call interp_grad(nnode,ncomp,cartd(:,:,igaus),elval,gpval(:,:,igaus))
    end do

  end subroutine elem2gauss_grad

  subroutine elem2gauss_hess(ncomp,nnode,ntens,ngaus,hessi,lnods,glval,gphes)
    !------------------------------------------------------------------------
    ! This routine performs a gather and computes the velocity hessian at gauss points
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,ngaus,ntens,lnods(nnode)
    real(rp),    intent(in)  :: glval(ncomp,*)
    real(rp),    intent(in)  :: hessi(ntens,nnode,ngaus)
    real(rp),    intent(out) :: gphes(ntens,ncomp,ngaus)

    integer(ip) :: igaus
    real(rp)    :: elval(ncomp,nnode)

    call gather(ncomp,nnode,lnods,glval,elval)
    do igaus=1,ngaus
       call interp_hess(nnode,ncomp,ntens,hessi(:,:,igaus),elval,gphes(:,:,igaus))
    end do

  end subroutine elem2gauss_hess
  
  subroutine elem2var(ncomp,nnode,inipo,ndime,elvar,elval)
    !------------------------------------------------------------------------
    ! This routine extracts a given components of an elemental array
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ncomp,nnode,inipo,ndime
    real(rp),    intent(in)  :: elvar(nnode,ncomp)
    real(rp),    intent(out) :: elval(ndime,nnode)
    
    integer(ip) :: idime,kdime
    
    kdime = 0
    do idime=inipo,inipo+ndime-1
       kdime = kdime + 1
       elval(kdime,:) = elvar(:,idime)
    end do

  end subroutine elem2var

  subroutine elmuv2elmat(ncomp,nnode,elmuv,elmat)
    !------------------------------------------------------------------------
    ! This routine assemble elmuv matrix to elemental matrix emat
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: ncomp,nnode
    real(rp),    intent(in)    :: elmuv(nnode,nnode)
    real(rp),    intent(inout) :: elmat(ncomp,ncomp,nnode,nnode)

    integer(ip) :: jdime,inode,jnode
    do jnode=1,nnode
       do inode=1,nnode
          do jdime=1,ncomp
             elmat(jdime,jdime,inode,jnode) = elmuv(inode,jnode) + elmat(jdime,jdime,inode,jnode)       
          end do
       end do
    end do

  end subroutine elmuv2elmat
  
  subroutine elmuv2elmat_lump(ncomp,nnode,elmuv,elmat)
    !------------------------------------------------------------------------
    ! This routine assemble elmuv matrix to elemental lumped matrix emat
    !------------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: ncomp,nnode
    real(rp),    intent(in)    :: elmuv(nnode,nnode)
    real(rp),    intent(inout) :: elmat(ncomp,ncomp,nnode,nnode)

    integer(ip) :: jdime,inode,jnode
    do jnode=1,nnode
       do inode=1,nnode
          do jdime=1,ncomp
             elmat(jdime,jdime,inode,inode) = elmuv(inode,jnode) + elmat(jdime,jdime,inode,inode)       
          end do
       end do
    end do

  end subroutine elmuv2elmat_lump

end module element_gather_tools_names
