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
  public :: fem_conditions_create, fem_conditions_free, fem_conditions_copy, fem_conditions_apply_renum

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
  subroutine fem_conditions_free(cnd)
    implicit none
    type(fem_conditions), intent(inout) :: cnd

    cnd%ncode=-1
    cnd%nvalu=-1
    cnd%ncond=-1

    call memfree (cnd%code,__FILE__,__LINE__)
    call memfree (cnd%valu,__FILE__,__LINE__)

  end subroutine fem_conditions_free

end module fem_conditions_names
