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
module serial_environment_names
use types_names
  use environment_names
  implicit none
  
  private
  ! Abstract environment
  type, extends(environment_t) :: serial_environment_t
   contains
     procedure :: info                        => serial_environment_info
     procedure :: am_i_l1_task                => serial_environment_am_i_l1_task
     procedure :: l1_lgt1_bcast               => serial_environment_l1_lgt1_bcast
     procedure :: l1_barrier                  => serial_environment_l1_barrier
     procedure :: l1_sum_scalar_rp            => serial_environment_l1_sum_scalar_rp
     procedure :: l1_sum_vector_rp            => serial_environment_l1_sum_vector_rp
  end type serial_environment_t
  
  public :: serial_environment_t
  
contains

  subroutine serial_environment_l1_barrier(this) 
    implicit none
    class(serial_environment_t),intent(in)  :: this
  end subroutine serial_environment_l1_barrier

  subroutine serial_environment_info(this,me,np) 
    implicit none
    class(serial_environment_t),intent(in)  :: this
    integer(ip)              ,intent(out) :: me
    integer(ip)              ,intent(out) :: np
    me = 0
    np = 1 
  end subroutine serial_environment_info
  
  function serial_environment_am_i_l1_task(this) 
    implicit none
    class(serial_environment_t) ,intent(in)  :: this
    logical                                  :: serial_environment_am_i_l1_task 
    serial_environment_am_i_l1_task = .true.
  end function serial_environment_am_i_l1_task
  
  subroutine serial_environment_l1_lgt1_bcast (this, condition)
    implicit none
    class(serial_environment_t) ,intent(in)    :: this
    logical                     ,intent(inout) :: condition
  end subroutine serial_environment_l1_lgt1_bcast
  
  subroutine serial_environment_l1_sum_scalar_rp (this,alpha)
    implicit none
    class(serial_environment_t) , intent(in)    :: this
    real(rp)                    , intent(inout) :: alpha
  end subroutine serial_environment_l1_sum_scalar_rp
     
 subroutine serial_environment_l1_sum_vector_rp(this,alpha)
    implicit none
    class(serial_environment_t) , intent(in)    :: this
    real(rp)                    , intent(inout) :: alpha(:) 
 end subroutine serial_environment_l1_sum_vector_rp
  
end module serial_environment_names
