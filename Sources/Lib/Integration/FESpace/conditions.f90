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
  use function_names
  implicit none
# include "debug.i90"
  private
  
  ! Data type which is pretended to be extended by the user s.t.
  ! he might customize which function to be imposed on each boundary
  ! indicator + component combination
  type, abstract :: conditions_t
   contains
     procedure(get_number_components_interface), deferred :: get_number_components
     procedure(get_components_code_interface)  , deferred :: get_components_code
     procedure(get_function_interface)         , deferred :: get_function
  end type conditions_t

  ! Types
  public :: conditions_t
  
  abstract interface
     function get_number_components_interface(this)
       import :: conditions_t, ip
       implicit none
       class(conditions_t), intent(in) :: this
       integer(ip) :: get_number_components_interface
     end function get_number_components_interface

     subroutine get_components_code_interface(this, boundary_id, components_code)
       import :: conditions_t, ip
       implicit none
       class(conditions_t), intent(in)  :: this
       integer(ip)            , intent(in)  :: boundary_id
       logical                , intent(out) :: components_code(:)
     end subroutine get_components_code_interface
  
     subroutine get_function_interface ( this, boundary_id, component_id, function )
       import :: conditions_t, ip, scalar_function_t
       implicit none
       class(conditions_t), target  , intent(in)  :: this
       integer(ip)                      , intent(in)  :: boundary_id
       integer(ip)                      , intent(in)  :: component_id
       class(scalar_function_t), pointer, intent(out) :: function
     end subroutine get_function_interface 
  end interface

end module conditions_names
