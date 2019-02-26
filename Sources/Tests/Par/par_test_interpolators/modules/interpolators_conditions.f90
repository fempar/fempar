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
module interpolators_conditions_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(conditions_t) :: interpolators_conditions_t
     private
     integer(ip)                       :: num_components
     class(scalar_function_t), pointer :: boundary_function_Hx
     class(scalar_function_t), pointer :: boundary_function_Hy 
     class(scalar_function_t), pointer :: boundary_function_Hz 
     class(scalar_function_t), pointer :: boundary_function_pressure
   contains
     procedure :: set_num_components          => interpolators_conditions_set_num_components
     procedure :: set_boundary_function_Hx    => interpolators_conditions_set_boundary_function_Hx
     procedure :: set_boundary_function_Hy    => interpolators_conditions_set_boundary_function_Hy
     procedure :: set_boundary_function_Hz    => interpolators_conditions_set_boundary_function_Hz
     procedure :: set_boundary_function_pressure => interpolators_conditions_set_boundary_function_pressure
     procedure :: get_num_components       => interpolators_conditions_get_num_components  
     procedure :: get_components_code         => interpolators_conditions_get_components_code
     procedure :: get_function                => interpolators_conditions_get_function
  end type interpolators_conditions_t
  
  public :: interpolators_conditions_t
  
contains

  subroutine interpolators_conditions_set_num_components (this, num_components)
    implicit none
    class(interpolators_conditions_t), intent(inout) :: this
    integer(ip)                      , intent(in)    :: num_components
    this%num_components = num_components
  end subroutine interpolators_conditions_set_num_components 
  
    subroutine interpolators_conditions_set_boundary_function_Hx (this, scalar_function)
    implicit none
    class(interpolators_conditions_t), intent(inout)      :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_Hx => scalar_function
  end subroutine interpolators_conditions_set_boundary_function_Hx
  
    subroutine interpolators_conditions_set_boundary_function_Hy (this, scalar_function)
    implicit none
    class(interpolators_conditions_t), intent(inout)      :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_Hy => scalar_function
  end subroutine interpolators_conditions_set_boundary_function_Hy
  
    subroutine interpolators_conditions_set_boundary_function_Hz (this, scalar_function)
    implicit none
    class(interpolators_conditions_t), intent(inout)      :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_Hz => scalar_function
  end subroutine interpolators_conditions_set_boundary_function_Hz
  
  subroutine interpolators_conditions_set_boundary_function_pressure (this, scalar_function)
    implicit none
    class(interpolators_conditions_t), intent(inout)          :: this
    class(scalar_function_t)           , target, intent(in) :: scalar_function
    this%boundary_function_pressure => scalar_function
  end subroutine interpolators_conditions_set_boundary_function_pressure

  function interpolators_conditions_get_num_components(this)
    implicit none
    class(interpolators_conditions_t), intent(in) :: this
    integer(ip) :: interpolators_conditions_get_num_components
    assert ( this%num_components == 2 .or. this%num_components == 3 .or. this%num_components == 4) 
    interpolators_conditions_get_num_components = this%num_components
  end function interpolators_conditions_get_num_components

  subroutine interpolators_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(interpolators_conditions_t), intent(in)  :: this
    integer(ip)                       , intent(in)  :: boundary_id
    logical                           , intent(out) :: components_code(:)
    assert ( size(components_code) == 2 .or. size(components_code) == 3 .or. size(components_code)==4 )
    components_code(1:size(components_code)) = .false.
    if ( boundary_id == 1 ) then
      components_code(1:size(components_code)) = .true.
    end if
  end subroutine interpolators_conditions_get_components_code
  
  subroutine interpolators_conditions_get_function ( this, boundary_id, component_id, function )
    implicit none
    class(interpolators_conditions_t), target     , intent(in) :: this
    integer(ip)                                    , intent(in)  :: boundary_id
    integer(ip)                                    , intent(in)  :: component_id
    class(scalar_function_t)          , pointer    , intent(out) :: function
    assert ( component_id == 1 .or. component_id == 2 .or. component_id == 3 .or. component_id == 4 )
 
  if (this%num_components == 3) then 
       if ( component_id == 1) then 
          function => this%boundary_function_Hx
       else if ( component_id == 2 ) then 
          function => this%boundary_function_Hy
       else if ( component_id == 3 ) then 
          function => this%boundary_function_pressure 
       end if

    else if (this%num_components == 4) then 
       if ( component_id == 1) then 
          function => this%boundary_function_Hx
       else if ( component_id == 2 ) then 
          function => this%boundary_function_Hy
       else if ( component_id == 3 ) then 
          function => this%boundary_function_Hz   
       else if ( component_id == 4 ) then 
          function => this%boundary_function_pressure
       end if
    end if

  end subroutine interpolators_conditions_get_function 

end module interpolators_conditions_names
