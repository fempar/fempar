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
module types_names
  !-----------------------------------------------------------------------
  !    This module contains kind and type definitions.
  !-----------------------------------------------------------------------
  use, intrinsic :: iso_fortran_env, only: INT8, INT32, INT64, REAL64
  implicit none

  integer, parameter       :: ieep = INT8   ! Integer precision for buffers in element exchanges
  integer, parameter       :: ip   = INT32  ! Integer precision
  integer, parameter       :: rp   = REAL64 ! Real precision
  !integer, parameter       :: lg   = 1     ! Logical precision

  integer(ip)  , parameter :: imp  = INT64  ! Integer precision, 
  ! memory consumption

  ! integer(ip), parameter :: igp = 8    ! Integer precision, 
  !                                      ! for global ids

  ! This is a 8-byte integer, and typically different from default integer(ip).
  ! (A.F.M.: fem has to be responsible for the following two definitions)
  ! psb_long_int_k_ kind parameter is required by some psb modules. I extracted
  ! the definition of this parameter from psb_const_mod PSBLAS 2.4.0 module
  integer, parameter  :: longndig=12
  integer, parameter  :: psb_long_int_k_ = selected_int_kind(longndig)

  integer(ip)  , parameter :: igp = psb_long_int_k_ ! Integer precision, 
  ! for global ids

  ! Error constant psb_success_ (extracted from PSBLAS 2.4.0 
  ! psb_const_mod module)
  integer, parameter  :: psb_success_=0
  
  type i1p_t
     integer(ip), pointer :: p(:) => NULL()
  end type i1p_t
  type i2p_t
     integer(ip), pointer :: p(:,:) => NULL()
  end type i2p_t
  type i3p_t
     integer(ip), pointer :: p(:,:,:) => NULL()
  end type i3p_t
  type r1p_t
     real(rp),    pointer :: p(:) => NULL()
  end type r1p_t
  type r2p_t
     real(rp),    pointer :: p(:,:) => NULL()
  end type r2p_t
  type r3p_t
     real(rp),    pointer :: p(:,:,:) => NULL()
  end type r3p_t

  ! Frequently used mathematical constants:
  real(rp),    parameter :: pi    = 3.141592653589793238462643383279502884197_rp
  real(rp),    parameter :: pio2  = 1.570796326794896619231321691639751442099_rp
  real(rp),    parameter :: twopi = 6.283185307179586476925286766559005768394_rp
  real(rp),    parameter :: sqrt2 = 1.414213562373095048801688724209698078570_rp
  real(rp),    parameter :: euler = 0.577215664901532860606512090082402431042_rp

  ! Actions related to free routines
  integer (ip), parameter  :: free_numerical_setup = 7
  integer (ip), parameter  :: free_symbolic_setup = 8
  integer (ip), parameter  :: free_clean  = 9
  
  ! Number of space dimensions for statically allocated data types (see, e.g., field.f90)
  integer(ip), parameter :: SPACE_DIM = 3
  
  !integer(ip), parameter :: MAX_NUM_LEVELS = 5

  integer(ieep), parameter :: mold(1) = [0_ieep]
  integer(ip)  , parameter :: size_of_ip = size(transfer(1_ip, mold))
  integer(ip)  , parameter :: size_of_igp = size(transfer(1_igp ,mold))

  
  character(len=*), parameter :: dir_path_key     = 'dir_path'
  character(len=*), parameter :: prefix_key       = 'prefix'
  character(len=*), parameter :: dir_path_out_key = 'dir_path_out'

  interface
     subroutine runend
     end subroutine runend
  end interface

end module types_names
		
