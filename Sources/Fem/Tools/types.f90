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
module types
  !-----------------------------------------------------------------------
  !    This module contains kind and type definitions.
  !-----------------------------------------------------------------------
  implicit none

  integer, parameter       :: ieep = 1    ! Integer precision for buffers in element exchanges
  integer, parameter       :: ip   = 4    ! Integer precision
  integer, parameter       :: rp   = 8    ! Real precision
  integer, parameter       :: lg   = 1    ! Logical precision

  integer(ip)  , parameter :: imp = 8    ! Integer precision, 
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

  integer(ip), parameter :: nsi_code=1, ela_code=2, cdr_code=3, adr_code=4 ! Temp here, probably better in prob
  integer(ip), parameter :: imh_code=5, dcy_code=6, mss_code=7, lap_code=8 ! Temp here, probably better in prob

  ! Integration rules
  integer(ip), parameter :: ruqope_id = 1
  integer(ip), parameter :: ruqclo_id = 2
  integer(ip), parameter :: rutope_id = 3
  integer(ip), parameter :: rutclo_id = 4
  integer(ip), parameter :: rupope_id = 5
  integer(ip), parameter :: rupclo_id = 6
  integer(ip), parameter :: ruqope_tp_id = 7

!!$  ! List of convergence criteria available for Krylov subspace 
!!$  ! iterative solvers (#1 is the one implemented in FELAP)
!!$  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
!!$  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
!!$  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
!!$  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
!!$  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
!!$  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
!!$  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
!!$  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
!!$  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
!!$  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
!!$                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
!!$                                                            ! r(i) is the residual at the i-th iteration


!!$  integer (ip), parameter :: mgs  = 1 ! mgs : Modified Gram-Schmidt (appropriate for serial GMRES)
!!$  integer (ip), parameter :: icgs = 2 ! icgs: Iterative Classical Gram-Schmidt
!!$                                      ! (appropriate for distributed GMRES)
!!$
!!$  integer (ip), parameter :: default_luout    = 6
!!$  integer (ip), parameter :: default_stopc    = res_nrmgiven_res_nrmgiven
!!$  integer (ip), parameter :: default_trace    = 0
!!$  integer (ip), parameter :: default_itmax    = 1000
!!$  real    (rp), parameter :: default_rtol     = 1.0e-08_rp
!!$  integer (ip), parameter :: default_dkrymax  = 30
!!$  integer (ip), parameter :: default_orto     = icgs
!!$  real    (rp), parameter :: default_relax    = 1.0_rp
!!$
!!$  ! List of Krylov subspace methods available
!!$  integer (ip), parameter :: cg               = 0
!!$  integer (ip), parameter :: lgmres           = 1
!!$  integer (ip), parameter :: rgmres           = 2
!!$  integer (ip), parameter :: fgmres           = 3
!!$  integer (ip), parameter :: richard          = 4  ! Richardson method ( fixed-point iteration,
!!$                                                   ! not actually a Krylov subspace method)
!!$  integer (ip), parameter :: default_kry_meth = lgmres 

  type i1p
     integer(ip), pointer :: l(:) => NULL()
  end type i1p
  type i2p
     integer(ip), pointer :: l(:,:) => NULL()
  end type i2p
  type i3p
     integer(ip), pointer :: l(:,:,:) => NULL()
  end type i3p
  type r1p
     real(rp),    pointer :: a(:) => NULL()
  end type r1p
  type r2p
     real(rp),    pointer :: a(:,:) => NULL()
  end type r2p
  type r3p
     real(rp),    pointer :: a(:,:,:) => NULL()
  end type r3p

  type list
     integer(ip) :: n
     integer(ip), allocatable :: p(:) 
     integer(ip), allocatable :: l(:) 
  end type list

  type list_2d
     integer(ip) :: n1
     integer(ip) :: n2
     integer(ip), allocatable :: p(:) 
     integer(ip), allocatable :: l(:,:) 
  end type list_2d

  type list_pointer
     type(list)          , pointer :: p => NULL()
  end type list_pointer

  ! Frequently used mathematical constants:
  real(rp),    parameter :: pi    = 3.141592653589793238462643383279502884197_rp
  real(rp),    parameter :: pio2  = 1.570796326794896619231321691639751442099_rp
  real(rp),    parameter :: twopi = 6.283185307179586476925286766559005768394_rp
  real(rp),    parameter :: sqrt2 = 1.414213562373095048801688724209698078570_rp
  real(rp),    parameter :: euler = 0.577215664901532860606512090082402431042_rp
  real(rp),    parameter :: zero_rp=epsilon(0.0_rp)
  integer(ip), parameter :: mone=-1,zero=0,one=1,two=2,three=3,four=4
  integer(ip), parameter :: five=5,six=6,seven=7,eight=8,nine=9

  ! Actions related to free routines
  integer (ip), parameter  :: free_only_values = 7
  integer (ip), parameter  :: free_only_struct = 8
  integer (ip), parameter  :: free_clean       = 9

  interface
     subroutine runend
     end subroutine runend
  end interface

  ! Functions
  public :: print_list_2d

contains

  subroutine print_list_2d( lunou, list2 )
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(list_2d)    , intent(in) :: list2

    integer(ip) :: i,j

    write(lunou,*) '****PRINT LIST_2D****'
    write(lunou,*) 'size total list:',list2%n1
    write(lunou,*) 'size components:',list2%n2
    do j = 1,list2%n2
       do i = 1,list2%n1
          write(lunou,*) 'l(',i,',',j,')',list2%l(list2%p(i):list2%p(i+1)-1,j)
       end do
    end do
    write(lunou,*) '****END PRINT LIST_2D****'

  end subroutine print_list_2d

  end module types
