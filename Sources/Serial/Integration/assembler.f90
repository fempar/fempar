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
module assembler_names
  use types_names
  use dof_descriptor_names
  use finite_element_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: assembler_t
    contains
	  procedure (assembly_interface), deferred :: assembly
  end type
  
  abstract interface
     subroutine assembly_interface(this,dof_descriptor,finite_element) 
       import :: assembler_t, dof_descriptor_t, finite_element_t
       implicit none
       class(assembler_t)       , intent(inout) :: this
       type(dof_descriptor_t)   , intent(in)    :: dof_descriptor
       type(finite_element_t)   , intent(in)    :: finite_element 
     end subroutine assembly_interface
  end interface
	 
  ! Data types
  public :: assembler_t

end module assembler_names
