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
module par_nsi_conditions_names
  use fempar_names
  use par_nsi_analytical_functions_names

  implicit none
# include "debug.i90"
  private

  type, extends(conditions_t) :: par_nsi_conditions_t
     private
     integer(ip)                       :: number_components
     integer(ip)                       :: number_dimensions
     class(scalar_function_t), pointer :: boundary_function_p
     type(vector_component_function_t) :: boundary_function_x
     type(vector_component_function_t) :: boundary_function_y
     type(vector_component_function_t) :: boundary_function_z
  contains
     procedure :: set_boundary_function       => par_nsi_conditions_set_boundary_function
     procedure :: set_number_components       => par_nsi_conditions_set_number_components
     procedure :: set_number_dimensions       => par_nsi_conditions_set_number_dimensions
     procedure :: get_number_components       => par_nsi_conditions_get_number_components  
     procedure :: get_num_components          => par_nsi_conditions_get_num_components
     procedure :: get_components_code         => par_nsi_conditions_get_components_code
     procedure :: get_function                => par_nsi_conditions_get_function
  end type par_nsi_conditions_t
  
  public :: par_nsi_conditions_t
  
contains

  function par_nsi_conditions_get_number_components(this)
    implicit none
    class(par_nsi_conditions_t), intent(in) :: this
    integer(ip) :: par_nsi_conditions_get_number_components
    par_nsi_conditions_get_number_components = this%number_components
  end function par_nsi_conditions_get_number_components
  
  subroutine par_nsi_conditions_set_number_components (this, number_components)
    implicit none
    class(par_nsi_conditions_t), intent(inout) :: this
    integer(ip)            , intent(in)    :: number_components
    this%number_components = number_components
  end subroutine par_nsi_conditions_set_number_components

  subroutine par_nsi_conditions_set_number_dimensions (this, number_dimensions)
    implicit none
    class(par_nsi_conditions_t), intent(inout) :: this
    integer(ip)                  , intent(in)    :: number_dimensions
    this%number_dimensions = number_dimensions
  end subroutine par_nsi_conditions_set_number_dimensions
  
  subroutine par_nsi_conditions_set_boundary_function(this,boundary_function_u,boundary_function_p)
    implicit none
    class(par_nsi_conditions_t), intent(inout) :: this
    class(vector_function_t), target, intent(in) :: boundary_function_u
    class(scalar_function_t), target, intent(in) :: boundary_function_p
    this%boundary_function_p => boundary_function_p
    call this%boundary_function_x%set(boundary_function_u,1)
    call this%boundary_function_y%set(boundary_function_u,2)
    call this%boundary_function_z%set(boundary_function_u,SPACE_DIM)
  end subroutine par_nsi_conditions_set_boundary_function

  subroutine par_nsi_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(par_nsi_conditions_t), intent(in)  :: this
    integer(ip)                  , intent(in)  :: boundary_id
    logical                      , intent(out) :: components_code(:)
    assert ( size(components_code) >= this%number_components )
    components_code(1:this%number_dimensions) = .true.
    components_code(this%number_dimensions+1) = .false.
    if ( boundary_id == 2 ) then
       components_code(this%number_dimensions+1) = .true.
    end if
  end subroutine par_nsi_conditions_get_components_code
  
  function par_nsi_conditions_get_num_components(this)
    implicit none
    class(par_nsi_conditions_t), intent(in)  :: this
    integer(ip) :: par_nsi_conditions_get_num_components
    par_nsi_conditions_get_num_components = this%number_components 
  end function par_nsi_conditions_get_num_components
  
  subroutine par_nsi_conditions_get_function ( this, boundary_id, component_id, function )
    implicit none
    class(par_nsi_conditions_t), target, intent(in)  :: this
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
          case (3)
             function => this%boundary_function_p
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
          case (4)
             function => this%boundary_function_p
          case default
             check(.false.)
          end select
       else
          ! Think about how to code other Dirichlet conditions
          check(.false.)
       end if
    end if
     
   end subroutine par_nsi_conditions_get_function 

end module par_nsi_conditions_names
