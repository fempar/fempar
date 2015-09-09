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
module conditions_names
use types_names
use memor_names
use stdio_names
  use renumbering_names
# include "debug.i90"
  implicit none
  private

  type conditions_t
     integer(ip)                :: &
          ncode=1,                 &         ! Number of codes  (=nvars)
          nvalu=1,                 &         ! Number of values (=nvars)
          ncond=1                            ! Number of conds  (=npoin)
     integer(ip), allocatable   :: &
          code(:,:)                          ! Codes
     real(rp), allocatable      :: &
          valu(:,:)                          ! Values
  end type conditions_t

  ! Types
  public :: conditions_t

  ! Methods
  public :: conditions_create, conditions_free, conditions_copy, conditions_apply_renumbering

contains

  !===============================================================================================
  subroutine conditions_create(ncode,nvalu,ncond,cnd)
    implicit none
    integer(ip)         , intent(in)    :: ncode,nvalu,ncond
    type(conditions_t), intent(inout) :: cnd

    cnd%ncode=ncode
    cnd%nvalu=nvalu
    cnd%ncond=ncond

    call memalloc (cnd%ncode,cnd%ncond,cnd%code, __FILE__,__LINE__)
    call memalloc (cnd%nvalu,cnd%ncond,cnd%valu, __FILE__,__LINE__)
    cnd%code=0
    cnd%valu=0.0_rp

  end subroutine conditions_create

  !===============================================================================================
  subroutine conditions_copy(cnd_old,cnd_new)
    implicit none
    type(conditions_t), intent(in)    :: cnd_old
    type(conditions_t), intent(inout) :: cnd_new

    call conditions_create( cnd_old%ncode, cnd_old%nvalu, cnd_old%ncond, cnd_new)
    cnd_new%code=cnd_old%code
    cnd_new%valu=cnd_old%valu

  end subroutine conditions_copy

  !===============================================================================================
  subroutine conditions_apply_renumbering(renumbering, cnd)
    implicit none
    type(renumbering_t)         , intent(in)    :: renumbering
    type(conditions_t), intent(inout) :: cnd

    integer(ip), allocatable :: tmp_int(:,:)
    real(rp)   , allocatable :: tmp_real(:,:)

    ! Renum code 2D array
    call memalloc ( cnd%ncode, cnd%ncond, tmp_int, __FILE__, __LINE__ )
    tmp_int = cnd%code
    call renumbering_apply(cnd%ncode, renumbering, tmp_int, cnd%code)
    
    ! Renum valu 2D array
    call memalloc ( cnd%nvalu, cnd%ncond, tmp_real, __FILE__, __LINE__ )
    tmp_real = cnd%valu
    call renumbering_apply(cnd%nvalu, renumbering, tmp_real, cnd%valu)

    ! Free work space
    call memfree ( tmp_real, __FILE__, __LINE__ )
    call memfree ( tmp_int, __FILE__, __LINE__ )

  end subroutine conditions_apply_renumbering

  !===============================================================================================
  subroutine conditions_free(cnd)
    implicit none
    type(conditions_t), intent(inout) :: cnd

    cnd%ncode=-1
    cnd%nvalu=-1
    cnd%ncond=-1

    call memfree (cnd%code,__FILE__,__LINE__)
    call memfree (cnd%valu,__FILE__,__LINE__)

  end subroutine conditions_free

end module conditions_names
