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

  ! Here we could use parametrized data types!
  type, extends(scalar_function_t) :: boundary_function_t
     class(vector_function_t), pointer :: solution_function_u
    private
   contains
     procedure :: set_solution_function => boundary_function_set_solution_function
  end type boundary_function_t
  
  type, extends(boundary_function_t) :: boundary_function_ux_t
    private
   contains
     procedure :: get_value_space => boundary_function_ux_get_value_space
  end type boundary_function_ux_t
  
  type, extends(boundary_function_t) :: boundary_function_uy_t
    private
   contains
     procedure :: get_value_space => boundary_function_uy_get_value_space
  end type boundary_function_uy_t
  
  type, extends(boundary_function_t) :: boundary_function_uz_t
    private
   contains
     procedure :: get_value_space => boundary_function_uz_get_value_space
  end type boundary_function_uz_t
  
  type, extends(scalar_function_t) :: boundary_function_p_t
    private
   contains
     procedure :: get_value_space => boundary_function_p_get_value_space
  end type boundary_function_p_t
  
  !===============================================================================================
  
  type, extends(conditions_t) :: fempar_sm_conditions_t
     private
     integer(ip)                       :: number_components
     integer(ip)                       :: number_dimensions
     type(fempar_sm_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(boundary_function_ux_t) :: boundary_function_ux
     type(boundary_function_uy_t) :: boundary_function_uy
     type(boundary_function_uz_t) :: boundary_function_uz
     type(boundary_function_p_t)  :: boundary_function_p
  contains
     procedure :: set_analytical_functions    => fempar_sm_conditions_set_analytical_functions
     procedure :: set_number_components       => fempar_sm_conditions_set_number_components
     procedure :: set_number_dimensions       => fempar_sm_conditions_set_number_dimensions
     procedure :: get_number_components       => fempar_sm_conditions_get_number_components  
     procedure :: get_components_code         => fempar_sm_conditions_get_components_code
     procedure :: get_function                => fempar_sm_conditions_get_function
  end type fempar_sm_conditions_t
  
  public :: fempar_sm_conditions_t
  
contains

  subroutine fempar_sm_conditions_set_analytical_functions ( this, analytical_functions )
     implicit none
     class(fempar_sm_conditions_t)         , intent(inout)      :: this
     type(fempar_sm_analytical_functions_t), target, intent(in) :: analytical_functions
     this%analytical_functions => analytical_functions
     call this%boundary_function_ux%set_solution_function(this%analytical_functions%get_solution_function_u())
     call this%boundary_function_uy%set_solution_function(this%analytical_functions%get_solution_function_u())
     call this%boundary_function_uz%set_solution_function(this%analytical_functions%get_solution_function_u())
  end subroutine fempar_sm_conditions_set_analytical_functions

  subroutine boundary_function_set_solution_function ( this, solution_function)
    implicit none
    class(boundary_function_t), intent(inout)    :: this
    class(vector_function_t), target, intent(in) :: solution_function
    this%solution_function_u  => solution_function
  end subroutine boundary_function_set_solution_function

  subroutine boundary_function_ux_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_ux_t), intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    real(rp)                     , intent(inout) :: result
    type(vector_field_t)                         :: u_field
    call this%solution_function_u%get_value(point,u_field)
    result = u_field%get(1)
  end subroutine boundary_function_ux_get_value_space

  subroutine boundary_function_uy_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_uy_t), intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    real(rp)                     , intent(inout) :: result
    type(vector_field_t)                         :: u_field
    call this%solution_function_u%get_value(point,u_field)
    result = u_field%get(2)
  end subroutine boundary_function_uy_get_value_space

  subroutine boundary_function_uz_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_uz_t), intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    real(rp)                     , intent(inout) :: result
    type(vector_field_t)                        :: u_field
    call this%solution_function_u%get_value(point,u_field)
    result = u_field%get(3)
  end subroutine boundary_function_uz_get_value_space

  subroutine boundary_function_p_get_value_space ( this, point, result )
    implicit none
    class(boundary_function_p_t), intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    class(vector_function_t)    , pointer       :: p_function
    result = 0.0_rp
    !p_function => this%analytical_functions%get_solution_function_p()    
    !call p_function%get_value(point,result)
  end subroutine boundary_function_p_get_value_space

  !===============================================================================================

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
  
  subroutine fempar_sm_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(fempar_sm_conditions_t), intent(in)  :: this
    integer(ip)                  , intent(in)  :: boundary_id
    logical                      , intent(out) :: components_code(:)
    assert ( size(components_code) >= this%number_components )
    components_code(1:this%number_components) = .false.
    if ( boundary_id == 1 ) then
      components_code(1:this%number_components-1) = .true.
    else  if ( boundary_id == 2 ) then
      components_code(1:this%number_components-1) = .true.
    end if
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
             function => this%boundary_function_ux
          case (2)
             function => this%boundary_function_uy
          case (3)
             function => this%boundary_function_p
          case default
             check(.false.)
          end select
       else if(this%number_dimensions==3) then
          select case ( component_id )
          case (1)
             function => this%boundary_function_ux
          case (2)
             function => this%boundary_function_uy
          case (3)
             function => this%boundary_function_uz
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
     
   end subroutine fempar_sm_conditions_get_function 

end module fempar_sm_conditions_names
