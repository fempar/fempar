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
!=============================================================================
! Auxiliar modules to handle ndofs by a template like technique
!=============================================================================
module matvec_ndof11
 use types
 integer(ip), parameter :: ndof1 = 1, ndof2 = 1
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof12
 use types
 integer(ip), parameter :: ndof1 = 1, ndof2 = 2
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof13
 use types
 integer(ip), parameter :: ndof1 = 1, ndof2 = 3
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof14
 use types
 integer(ip), parameter :: ndof1 = 1, ndof2 = 4
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof21
 use types
 integer(ip), parameter :: ndof1 = 2, ndof2 = 1
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof22
 use types
 integer(ip), parameter :: ndof1 = 2, ndof2 = 2
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof23
 use types
 integer(ip), parameter :: ndof1 = 2, ndof2 = 3
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof24
 use types
 integer(ip), parameter :: ndof1 = 2, ndof2 = 4
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof31
 use types
 integer(ip), parameter :: ndof1 = 3, ndof2 = 1
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof32
 use types
 integer(ip), parameter :: ndof1 = 3, ndof2 = 2
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof33
 use types
 integer(ip), parameter :: ndof1 = 3, ndof2 = 3
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof34
 use types
 integer(ip), parameter :: ndof1 = 3, ndof2 = 4
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof41
 use types
 integer(ip), parameter :: ndof1 = 4, ndof2 = 1
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof42
 use types
 integer(ip), parameter :: ndof1 = 4, ndof2 = 2
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof43
 use types
 integer(ip), parameter :: ndof1 = 4, ndof2 = 3
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------
module matvec_ndof44
 use types
 integer(ip), parameter :: ndof1 = 4, ndof2 = 4
 contains
# include "matvec_dof.i90"
end module
!!$!-----------------------------------------------------------------------
module matvec_ndofgg
 use types
 integer(ip), save :: ndof1, ndof2
 contains
# include "matvec_dof.i90"
end module
!-----------------------------------------------------------------------

module matvec_dof
  use types
  use matvec_ndof11, only : matvec_css11 => matvec_css, matvec_csr11 => matvec_csr, matvec_csr_symm11 => matvec_csr_symm,             & 
                          & matvec_csr_trans11 => matvec_csr_trans, matvec_csc11 => matvec_csc, matvec_css11_scal => matvec_css_scal, & 
                          & matvec_csr11_scal => matvec_csr_scal, matvec_csr_symm11_scal => matvec_csr_symm_scal,                     & 
                          & matvec_csr_trans11_scal => matvec_csr_trans_scal, matvec_csc11_scal => matvec_csc_scal,                   & 
                          & matvec_csr_symm_trans11_scal => matvec_csr_symm_trans_scal

  use matvec_ndof12, only : matvec_css12 => matvec_css, matvec_csr12 => matvec_csr, matvec_csr_symm12 => matvec_csr_symm,             & 
                          & matvec_csr_trans12 => matvec_csr_trans,matvec_csc12 => matvec_csc, matvec_css12_scal => matvec_css_scal,  &
                          & matvec_csr12_scal => matvec_csr_scal, matvec_csr_symm12_scal => matvec_csr_symm_scal,                     & 
                          & matvec_csr_trans12_scal => matvec_csr_trans_scal, matvec_csc12_scal => matvec_csc_scal,                   & 
                          & matvec_csr_symm_trans12_scal => matvec_csr_symm_trans_scal
!!$ 

  use matvec_ndof13, only : matvec_css13 => matvec_css, matvec_csr13 => matvec_csr, matvec_csr_symm13 => matvec_csr_symm,             & 
                          & matvec_csr_trans13 => matvec_csr_trans, matvec_csc13 => matvec_csc, matvec_css13_scal => matvec_css_scal, &
                          & matvec_csr13_scal => matvec_csr_scal, matvec_csr_symm13_scal => matvec_csr_symm_scal,                     & 
                          & matvec_csr_trans13_scal => matvec_csr_trans_scal, matvec_csc13_scal => matvec_csc_scal,                   &
                          & matvec_csr_symm_trans13_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof14, only : matvec_css14 => matvec_css, matvec_csr14 => matvec_csr, matvec_csr_symm14 => matvec_csr_symm,             & 
                          & matvec_csr_trans14 => matvec_csr_trans, matvec_csc14 => matvec_csc, matvec_css14_scal => matvec_css_scal, &
                          & matvec_csr14_scal => matvec_csr_scal, matvec_csr_symm14_scal => matvec_csr_symm_scal,                     & 
                          & matvec_csr_trans14_scal => matvec_csr_trans_scal, matvec_csc14_scal => matvec_csc_scal,                   &
                          & matvec_csr_symm_trans14_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof21, only : matvec_css21 => matvec_css, matvec_csr21 => matvec_csr, matvec_csr_symm21 => matvec_csr_symm,              & 
                          & matvec_csr_trans21 => matvec_csr_trans, matvec_csc21 => matvec_csc,  matvec_css21_scal => matvec_css_scal, & 
                          & matvec_csr21_scal => matvec_csr_scal, matvec_csr_symm21_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans21_scal => matvec_csr_trans_scal, matvec_csc21_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans21_scal => matvec_csr_symm_trans_scal


  use matvec_ndof22, only : matvec_css22 => matvec_css, matvec_csr22 => matvec_csr, matvec_csr_symm22 => matvec_csr_symm,              & 
                          & matvec_csr_trans22 => matvec_csr_trans, matvec_csc22 => matvec_csc, matvec_css22_scal => matvec_css_scal,  &
                          & matvec_csr22_scal => matvec_csr_scal, matvec_csr_symm22_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans22_scal => matvec_csr_trans_scal, matvec_csc22_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans22_scal => matvec_csr_symm_trans_scal

  use matvec_ndof23, only : matvec_css23 => matvec_css, matvec_csr23 => matvec_csr, matvec_csr_symm23 => matvec_csr_symm,              & 
                          & matvec_csr_trans23 => matvec_csr_trans, matvec_csc23 => matvec_csc, matvec_css23_scal => matvec_css_scal,  & 
                          & matvec_csr23_scal => matvec_csr_scal, matvec_csr_symm23_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans23_scal => matvec_csr_trans_scal, matvec_csc23_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans23_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof24, only : matvec_css24 => matvec_css, matvec_csr24 => matvec_csr, matvec_csr_symm24 => matvec_csr_symm,              & 
                          & matvec_csr_trans24 => matvec_csr_trans, matvec_csc24 => matvec_csc, matvec_css24_scal => matvec_css_scal,  & 
                          & matvec_csr24_scal => matvec_csr_scal, matvec_csr_symm24_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans24_scal => matvec_csr_trans_scal, matvec_csc24_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans24_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof31, only : matvec_css31 => matvec_css, matvec_csr31 => matvec_csr, matvec_csr_symm31 => matvec_csr_symm,              & 
                          & matvec_csr_trans31 => matvec_csr_trans, matvec_csc31 => matvec_csc, matvec_css31_scal => matvec_css_scal,  &
                          & matvec_csr31_scal => matvec_csr_scal, matvec_csr_symm31_scal => matvec_csr_symm_scal,                      &  
                          & matvec_csr_trans31_scal => matvec_csr_trans_scal, matvec_csc31_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans31_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof32, only : matvec_css32 => matvec_css, matvec_csr32 => matvec_csr, matvec_csr_symm32 => matvec_csr_symm,              & 
                          & matvec_csr_trans32 => matvec_csr_trans, matvec_csc32 => matvec_csc, matvec_css32_scal => matvec_css_scal,  &
                          & matvec_csr32_scal => matvec_csr_scal, matvec_csr_symm32_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans32_scal => matvec_csr_trans_scal, matvec_csc32_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans32_scal => matvec_csr_symm_trans_scal


  use matvec_ndof33, only : matvec_css33 => matvec_css, matvec_csr33 => matvec_csr, matvec_csr_symm33 => matvec_csr_symm,              & 
                          & matvec_csr_trans33 => matvec_csr_trans, matvec_csc33 => matvec_csc, matvec_css33_scal => matvec_css_scal,  &
                          & matvec_csr33_scal => matvec_csr_scal, matvec_csr_symm33_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans33_scal => matvec_csr_trans_scal, matvec_csc33_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans33_scal => matvec_csr_symm_trans_scal

  use matvec_ndof34, only : matvec_css34 => matvec_css, matvec_csr34 => matvec_csr, matvec_csr_symm34 => matvec_csr_symm,              & 
                          & matvec_csr_trans34 => matvec_csr_trans, matvec_csc34 => matvec_csc, matvec_css34_scal => matvec_css_scal,  & 
                          & matvec_csr34_scal => matvec_csr_scal, matvec_csr_symm34_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans34_scal => matvec_csr_trans_scal, matvec_csc34_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans34_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof41, only : matvec_css41 => matvec_css, matvec_csr41 => matvec_csr, matvec_csr_symm41 => matvec_csr_symm,              & 
                          & matvec_csr_trans41 => matvec_csr_trans, matvec_csc41 => matvec_csc, matvec_css41_scal => matvec_css_scal,  &
                          & matvec_csr41_scal => matvec_csr_scal, matvec_csr_symm41_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans41_scal => matvec_csr_trans_scal, matvec_csc41_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans41_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof42, only : matvec_css42 => matvec_css, matvec_csr42 => matvec_csr, matvec_csr_symm42 => matvec_csr_symm,              & 
                          & matvec_csr_trans42 => matvec_csr_trans, matvec_csc42 => matvec_csc,  matvec_css42_scal => matvec_css_scal, & 
                          & matvec_csr42_scal => matvec_csr_scal, matvec_csr_symm42_scal => matvec_csr_symm_scal,                      &  
                          & matvec_csr_trans42_scal => matvec_csr_trans_scal, matvec_csc42_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans42_scal => matvec_csr_symm_trans_scal
!!$

  use matvec_ndof43, only : matvec_css43 => matvec_css, matvec_csr43 => matvec_csr, matvec_csr_symm43 => matvec_csr_symm,              & 
                          & matvec_csr_trans43 => matvec_csr_trans, matvec_csc43 => matvec_csc,  matvec_css43_scal => matvec_css_scal, & 
                          & matvec_csr43_scal => matvec_csr_scal, matvec_csr_symm43_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans43_scal => matvec_csr_trans_scal, matvec_csc43_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans43_scal => matvec_csr_symm_trans_scal


  use matvec_ndof44, only : matvec_css44 => matvec_css, matvec_csr44 => matvec_csr, matvec_csr_symm44 => matvec_csr_symm,              & 
                          & matvec_csr_trans44 => matvec_csr_trans, matvec_csc44 => matvec_csc, matvec_css44_scal => matvec_css_scal,  & 
                          & matvec_csr44_scal => matvec_csr_scal, matvec_csr_symm44_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_trans44_scal => matvec_csr_trans_scal, matvec_csc44_scal => matvec_csc_scal,                    &
                          & matvec_csr_symm_trans44_scal => matvec_csr_symm_trans_scal

  use matvec_ndofgg, only : matvec_cssgg => matvec_css, matvec_csrgg => matvec_csr, matvec_csr_symmgg => matvec_csr_symm,              & 
                          & matvec_csr_transgg => matvec_csr_trans, matvec_cscgg => matvec_csc, matvec_cssgg_scal => matvec_css_scal,  & 
                          & matvec_csrgg_scal => matvec_csr_scal, matvec_csr_symmgg_scal => matvec_csr_symm_scal,                      & 
                          & matvec_csr_transgg_scal => matvec_csr_trans_scal, matvec_cscgg_scal => matvec_csc_scal, ndof1, ndof2,      &
                          & matvec_csr_symm_transgg_scal => matvec_csr_symm_trans_scal
  private

  public :: matvec_css, matvec_css_scal, matvec_csr, matvec_csr_scal, & 
          & matvec_csr_symm, matvec_csr_symm_scal,  matvec_csr_trans, & 
          & matvec_csr_trans_scal, matvec_csc, matvec_csc_scal,       &
          & matvec_csr_symm_trans_scal

contains

  subroutine matvec_css(nd1,nd2,nv,ia,is,ja,da,la,ua,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,ia(nv+1),is(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: x(nd2,nv),da(nd1,nd2,nv)
    real(rp)   , intent(in)  :: la(nd1,nd2,is(nv+1)-1),ua(nd2,nd1,is(nv+1)-1)
    real(rp)   , intent(out) :: y(nd1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_css11(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_css12(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_css13(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_css14(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_css21(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_css22(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_css23(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_css24(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_css31(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_css32(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_css33(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_css34(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_css41(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_css42(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_css43(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_css44(nv,ia,is,ja,da,la,ua,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_cssgg(nv,ia,is,ja,da,la,ua,x,y)
    end if

  end subroutine matvec_css

  subroutine matvec_css_scal(nd1,nd2,nv,ia,is,ja,da,la,ua,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,ia(nv+1),is(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: x(1,nv),da(1,1,nv)
    real(rp)   , intent(in)  :: la(1,1,is(nv+1)-1),ua(1,1,is(nv+1)-1)
    real(rp)   , intent(out) :: y(1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_css11_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_css12_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_css13_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_css14_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_css21_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_css22_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_css23_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_css24_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_css31_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_css32_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_css33_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_css34_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_css41_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_css42_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_css43_scal(nv,ia,is,ja,da,la,ua,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_css44_scal(nv,ia,is,ja,da,la,ua,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_cssgg_scal(nv,ia,is,ja,da,la,ua,x,y)
    end if

  end subroutine matvec_css_scal

  subroutine matvec_csr(nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1),x(nd2,nv2)
    real(rp)   , intent(out) :: y(nd1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr11(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr12(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr13(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr14(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr21(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr22(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr23(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr24(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr31(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr32(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr33(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr34(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr41(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr42(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr43(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr44(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csrgg(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr

  subroutine matvec_csr_scal(nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1),x(1,nv2)
    real(rp)   , intent(out) :: y(1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr11_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr12_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr13_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr14_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr21_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr22_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr23_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr24_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr31_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr32_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr33_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr34_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr41_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr42_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr43_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr44_scal(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csrgg_scal(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_scal

  subroutine matvec_csr_symm (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1),x(nd2,nv2)
    real(rp)   , intent(out) :: y(nd1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr_symm11(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr_symm12(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr_symm13(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr_symm14(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr_symm21(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr_symm22(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr_symm23(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr_symm24(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr_symm31(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr_symm32(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr_symm33(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr_symm34(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr_symm41(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr_symm42(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr_symm43(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr_symm44(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csr_symmgg(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_symm

  subroutine matvec_csr_symm_scal (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1),x(1,nv2)
    real(rp)   , intent(out) :: y(1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr_symm11_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr_symm12_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr_symm13_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr_symm14_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr_symm21_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr_symm22_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr_symm23_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr_symm24_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr_symm31_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr_symm32_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr_symm33_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr_symm34_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr_symm41_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr_symm42_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr_symm43_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr_symm44_scal(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csr_symmgg_scal(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_symm_scal


  subroutine matvec_csr_trans (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1),x(nd1,nv)
    real(rp)   , intent(out) :: y(nd2,nv2)


    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr_trans11(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr_trans12(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr_trans13(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr_trans14(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr_trans21(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr_trans22(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr_trans23(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr_trans24(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr_trans31(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr_trans32(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr_trans33(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr_trans34(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr_trans41(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr_trans42(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr_trans43(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr_trans44(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csr_transgg(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_trans

  subroutine matvec_csr_trans_scal (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1),x(1,nv)
    real(rp)   , intent(out) :: y(1,nv2)


    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr_trans11_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr_trans12_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr_trans13_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr_trans14_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr_trans21_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr_trans22_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr_trans23_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr_trans24_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr_trans31_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr_trans32_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr_trans33_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr_trans34_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr_trans41_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr_trans42_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr_trans43_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr_trans44_scal(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csr_transgg_scal(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_trans_scal

  subroutine matvec_csr_symm_trans_scal (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1),x(1,nv)
    real(rp)   , intent(out) :: y(1,nv2)


    if((nd1==1) .and. (nd2==1)) then
       call matvec_csr_symm_trans11_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csr_trans12_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csr_trans13_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csr_trans14_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csr_trans21_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csr_symm_trans22_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csr_trans23_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csr_trans24_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csr_trans31_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csr_trans32_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csr_symm_trans33_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csr_trans34_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csr_trans41_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csr_trans42_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csr_trans43_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csr_symm_trans44_scal(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_csr_symm_transgg_scal(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csr_symm_trans_scal

  subroutine matvec_csc(nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(nd2,nd1,ia(nv+1)-1),x(nd2,nv2)
    real(rp)   , intent(out) :: y(nd1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csc11(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csc12(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csc13(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csc14(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csc21(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csc22(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csc23(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csc24(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csc31(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csc32(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csc33(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csc34(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csc41(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csc42(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csc43(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csc44(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_cscgg(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csc

  subroutine matvec_csc_scal (nd1,nd2,nv,nv2,ia,ja,a,x,y)
    implicit none
    integer(ip), intent(in)  :: nd1,nd2,nv,nv2,ia(nv+1),ja(ia(nv+1)-1)
    real(rp)   , intent(in)  :: a(1,1,ia(nv+1)-1),x(1,nv2)
    real(rp)   , intent(out) :: y(1,nv)

    if((nd1==1) .and. (nd2==1)) then
       call matvec_csc11_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==2)) then
       call matvec_csc12_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==3)) then
       call matvec_csc13_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==1) .and. (nd2==4)) then
       call matvec_csc14_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==1)) then
       call matvec_csc21_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==2)) then
       call matvec_csc22_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==3)) then
       call matvec_csc23_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==2) .and. (nd2==4)) then
       call matvec_csc24_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==1)) then
       call matvec_csc31_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==2)) then
       call matvec_csc32_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==3)) then
       call matvec_csc33_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==3) .and. (nd2==4)) then
       call matvec_csc34_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==1)) then
       call matvec_csc41_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==2)) then
       call matvec_csc42_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==3)) then
       call matvec_csc43_scal(nv,nv2,ia,ja,a,x,y)
    else if((nd1==4) .and. (nd2==4)) then
       call matvec_csc44_scal(nv,nv2,ia,ja,a,x,y)
    else
       ndof1=nd1
       ndof2=nd2
       call matvec_cscgg(nv,nv2,ia,ja,a,x,y)
    end if

  end subroutine matvec_csc_scal

end module matvec_dof




       
