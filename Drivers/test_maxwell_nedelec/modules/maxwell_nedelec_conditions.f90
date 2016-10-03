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
module maxwell_nedelec_conditions_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(conditions_t) :: maxwell_nedelec_conditions_t
     private
     integer(ip)                       :: num_dimensions
	 class(scalar_function_t), pointer :: boundary_function_Hx
	 class(scalar_function_t), pointer :: boundary_function_Hy 
   contains
     procedure :: set_num_dimensions          => maxwell_nedelec_conditions_set_num_dimensions
	 procedure :: set_boundary_function_Hx    => maxwell_nedelec_conditions_set_boundary_function_Hx
	 procedure :: set_boundary_function_Hy    => maxwell_nedelec_conditions_set_boundary_function_Hy
     procedure :: get_number_components       => maxwell_nedelec_conditions_get_number_components  
     procedure :: get_components_code         => maxwell_nedelec_conditions_get_components_code
     procedure :: get_function                => maxwell_nedelec_conditions_get_function
  end type maxwell_nedelec_conditions_t
  
  public :: maxwell_nedelec_conditions_t
  
contains

  subroutine maxwell_nedelec_conditions_set_num_dimensions (this, num_dimensions)
    implicit none
    class(maxwell_nedelec_conditions_t), intent(inout) :: this
    integer(ip)                           , intent(in)    :: num_dimensions
    this%num_dimensions = num_dimensions
  end subroutine maxwell_nedelec_conditions_set_num_dimensions 
  
    subroutine maxwell_nedelec_conditions_set_boundary_function_Hx (this, scalar_function)
    implicit none
    class(maxwell_nedelec_conditions_t), intent(inout)      :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_Hx => scalar_function
  end subroutine maxwell_nedelec_conditions_set_boundary_function_Hx
  
    subroutine maxwell_nedelec_conditions_set_boundary_function_Hy (this, scalar_function)
    implicit none
    class(maxwell_nedelec_conditions_t), intent(inout)      :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_Hy => scalar_function
  end subroutine maxwell_nedelec_conditions_set_boundary_function_Hy

  function maxwell_nedelec_conditions_get_number_components(this)
    implicit none
    class(maxwell_nedelec_conditions_t), intent(in) :: this
    integer(ip) :: maxwell_nedelec_conditions_get_number_components
    assert ( this%num_dimensions == 2 .or. this%num_dimensions == 3 ) 
    maxwell_nedelec_conditions_get_number_components = this%num_dimensions
  end function maxwell_nedelec_conditions_get_number_components

  subroutine maxwell_nedelec_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(maxwell_nedelec_conditions_t), intent(in)  :: this
    integer(ip)                       , intent(in)  :: boundary_id
    logical                           , intent(out) :: components_code(:)
    assert ( size(components_code) == 2 .or. size(components_code) == 3 )
    components_code(1:size(components_code)) = .false.
    if ( boundary_id == 1 ) then
      components_code(1:size(components_code)) = .true.
    end if
  end subroutine maxwell_nedelec_conditions_get_components_code
  
  subroutine maxwell_nedelec_conditions_get_function ( this, boundary_id, component_id, function )
    implicit none
    class(maxwell_nedelec_conditions_t), target     , intent(in) :: this
    integer(ip)                                    , intent(in)  :: boundary_id
    integer(ip)                                    , intent(in)  :: component_id
    class(scalar_function_t)          , pointer    , intent(out) :: function
    assert ( component_id == 1 .or. component_id == 2 .or. component_id == 3  )

	if ( component_id == 1) then 
	function => this%boundary_function_Hx
 	else if ( component_id == 2 ) then 
	function => this%boundary_function_Hy
	end if 

  end subroutine maxwell_nedelec_conditions_get_function 

end module maxwell_nedelec_conditions_names
