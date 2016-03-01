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
module discrete_integration_names
  use field_names
  use reference_fe_names
  use types_names
  use assembler_names
  use serial_fe_space_names
  use memor_names

  implicit none
# include "debug.i90"
  private

  type, abstract :: discrete_integration_t
   contains
     procedure (integrate_interface), deferred :: integrate
     procedure                                 :: impose_strong_dirichlet_data
  end type discrete_integration_t

  type p_discrete_integration_t
     class(discrete_integration_t), pointer :: p => NULL()  
  end type p_discrete_integration_t

  public :: discrete_integration_t, p_discrete_integration_t

  abstract interface
     subroutine integrate_interface ( this, fe_space, assembler  )
       import :: discrete_integration_t, serial_fe_space_t, assembler_t
       implicit none
       class(discrete_integration_t), intent(in)    :: this
       class(serial_fe_space_t)     , intent(inout) :: fe_space
       class(assembler_t)           , intent(inout) :: assembler
     end subroutine integrate_interface
  end interface

contains

subroutine impose_strong_dirichlet_data ( this, elmat, elvec, code, value, nnode, num_fe_spaces )
 implicit none
 class(discrete_integration_t) :: this
 real(rp), intent(in) :: elmat(:,:)
 real(rp), intent(inout) :: elvec(:)  
 type(i1p_t), intent(in) :: code(:)
 type(r1p_t), intent(in) :: value(:)
 integer(ip), intent(in) :: nnode(:), num_fe_spaces
 integer(ip) :: i, c, ifes

 c = 0
 do ifes = 1, num_fe_spaces
    do i = 1, nnode(ifes)
       c = c + 1
       if ( code(ifes)%p(i) /= 0 ) then
          elvec = elvec - elmat(:,c)*value(ifes)%p(i)
       end if
    end do
 end do

end subroutine impose_strong_dirichlet_data

end module discrete_integration_names
