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
module matrix_names
  use types_names
  use operator_names

  implicit none
# include "debug.i90"
  private

  type, abstract, extends(operator_t) :: matrix_t
    contains
	  procedure (allocate_interface)      , deferred :: allocate
	  procedure (free_in_stages_interface), deferred :: free_in_stages
	  ! This subroutine is an instance of the Template Method pattern with
	  ! free_in_stages being the primitive method. According to this pattern,
	  ! template methods cannot be overrided by subclasses
  	  procedure, non_overridable :: free => matrix_free_template_method
  end type
  
  abstract interface
     ! Allocates the entries of the matrix once it has been created and symbolically set-up
     subroutine allocate_interface(this) 
	   import :: matrix_t
	   implicit none
       class(matrix_t)       , intent(inout) :: this
     end subroutine allocate_interface
	 ! Progressively free a matrix_t in three stages: action={free_numeric,free_symbolic,free_clean}
	 subroutine free_in_stages_interface(this,action) 
       import :: matrix_t, ip
       implicit none
       class(matrix_t)       , intent(inout) :: this
	   integer(ip)           , intent(in)    :: action
     end subroutine free_in_stages_interface
  end interface
	 
  ! Data types
  public :: matrix_t

contains  
  
   subroutine matrix_free_template_method ( this )
     implicit none
	 class(matrix_t), intent(inout) :: this
	 call this%free_in_stages(free_values)
	 call this%free_in_stages(free_struct)
	 call this%free_in_stages(free_clean)
   end subroutine matrix_free_template_method  
  
end module matrix_names
