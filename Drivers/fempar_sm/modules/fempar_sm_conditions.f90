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
module fempar_sm_conditions_names
  use fempar_names
  use fempar_sm_analytical_functions_names

  implicit none
# include "debug.i90"
  private

  type, extends(conditions_t) :: fempar_sm_conditions_t
     private
     integer(ip)                       :: number_components
     integer(ip)                       :: number_dimensions
     class(vector_function_t), pointer :: boundary_function   ! This is indeed not needed
     type(vector_component_function_t) :: boundary_function_x
     type(vector_component_function_t) :: boundary_function_y
     type(vector_component_function_t) :: boundary_function_z
  contains
     procedure :: set_boundary_function       => fempar_sm_conditions_set_boundary_function
     procedure :: set_number_components       => fempar_sm_conditions_set_number_components
     procedure :: set_number_dimensions       => fempar_sm_conditions_set_number_dimensions
     procedure :: get_number_components       => fempar_sm_conditions_get_number_components  
     procedure :: get_components_code         => fempar_sm_conditions_get_components_code
     procedure :: get_function                => fempar_sm_conditions_get_function
  end type fempar_sm_conditions_t
  
  public :: fempar_sm_conditions_t
  
contains

  function fempar_sm_conditions_get_number_components(this)
    implicit none
    class(fempar_sm_conditions_t), intent(in) :: this
    integer(ip) :: fempar_sm_conditions_get_number_components
    fempar_sm_conditions_get_number_components = this%number_components
  end function fempar_sm_conditions_get_number_components
  
  subroutine fempar_sm_conditions_set_number_components (this, number_components)
    implicit none
    class(fempar_sm_conditions_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: number_components
    this%number_components = number_components
  end subroutine fempar_sm_conditions_set_number_components

  subroutine fempar_sm_conditions_set_number_dimensions (this, number_dimensions)
    implicit none
    class(fempar_sm_conditions_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: number_dimensions
    this%number_dimensions = number_dimensions
  end subroutine fempar_sm_conditions_set_number_dimensions
  
  subroutine fempar_sm_conditions_set_boundary_function(this,boundary_function)
    implicit none
    class(fempar_sm_conditions_t), intent(inout) :: this
    class(vector_function_t), target, intent(in) :: boundary_function
    this%boundary_function => boundary_function
    call this%boundary_function_x%set(boundary_function,1)
    call this%boundary_function_y%set(boundary_function,2)
    call this%boundary_function_z%set(boundary_function,SPACE_DIM)
  end subroutine fempar_sm_conditions_set_boundary_function

  subroutine fempar_sm_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(fempar_sm_conditions_t), intent(in)  :: this
    integer(ip)                  , intent(in)  :: boundary_id
    logical                      , intent(out) :: components_code(:)
    assert ( size(components_code) >= this%number_components )
    components_code(1:this%number_components) = .false.
    if(boundary_id >= 1 ) components_code(1:this%number_dimensions) = .true.
    !if ( boundary_id == 1 ) then
    !  components_code(1:this%number_dimensions) = .true.
    !else  if ( boundary_id == 2 ) then
    !  components_code(1:this%number_dimensions) = .true.
    !end if
  end subroutine fempar_sm_conditions_get_components_code
  
  subroutine fempar_sm_conditions_get_function ( this, boundary_id, component_id, function )
    implicit none
    class(fempar_sm_conditions_t), target, intent(in)  :: this
    integer(ip)                        , intent(in)  :: boundary_id
    integer(ip)                        , intent(in)  :: component_id
    class(scalar_function_t), pointer  , intent(out) :: function

    nullify(function)
    if ( boundary_id >= 1 ) then
       ! Dirichlet conditions using the exact solution
       if(this%number_dimensions==2) then
          select case ( component_id )
          case (1)
             function => this%boundary_function_x
          case (2)
             function => this%boundary_function_y
          case default
             check(.false.)
          end select
       else if(this%number_dimensions==3) then
          select case ( component_id )
          case (1)
             function => this%boundary_function_x
          case (2)
             function => this%boundary_function_y
          case (3)
             function => this%boundary_function_z
          case default
             check(.false.)
          end select
       else
          ! Think about how to code other Dirichlet conditions
          check(.false.)
       end if
    end if
     
   end subroutine fempar_sm_conditions_get_function 

end module fempar_sm_conditions_names
