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
module scalar_names
  use types_names
  use assembler_names
  use dof_descriptor_names
  use finite_element_names
  implicit none
  private

  type, abstract, extends(assembler_t) :: scalar_t
  contains
	procedure (init_interface)     , deferred  :: init
	procedure (sum_interface)      , deferred  :: sum
    procedure (reduce_interface)   , deferred  :: reduce
	procedure (get_value_interface), deferred  :: get_value
	procedure (free_interface)     , deferred  :: free
	procedure :: assembly => scalar_assembly
  end type scalar_t
  
  ! Types
  public :: scalar_t

  abstract interface
     subroutine init_interface (this,b)
	   import :: scalar_t, rp
       implicit none
       class(scalar_t), intent(inout) :: this
       real(rp)       ,  intent(in)   :: b
     end subroutine init_interface

     subroutine sum_interface (this,b)
	   import :: scalar_t, rp
       implicit none
       class(scalar_t)       , intent(inout) :: this
       real(rp)              , intent(in)    :: b
     end subroutine sum_interface
  
     subroutine reduce_interface (this)
	   import :: scalar_t
       implicit none
       class(scalar_t), intent(inout) :: this
     end subroutine reduce_interface

     function get_value_interface (this)
	   import :: scalar_t, rp
       implicit none
       class(scalar_t), intent(inout) :: this
       real(rp)                       :: get_value_interface
     end function get_value_interface

     subroutine free_interface (this)
	   import :: scalar_t
       implicit none
       class(scalar_t), intent(inout) :: this
     end subroutine free_interface
  end interface

contains

   subroutine scalar_assembly(this, dof_descriptor, finite_element)
     implicit none
	 class(scalar_t), intent(inout) :: this
	 type(dof_descriptor_t)   , intent(in)    :: dof_descriptor
     type(finite_element_t)   , intent(in)    :: finite_element
	 call this%sum(finite_element%scalar)
   end subroutine scalar_assembly

end module scalar_names
  
