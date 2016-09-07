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

module mixed_laplacian_rt_analytical_functions_names
  use serial_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: pressure_source_term_t
   contains
     procedure :: get_value_space => pressure_source_term_get_value_space
  end type pressure_source_term_t

  type, extends(scalar_function_t) :: pressure_boundary_function_t
   contains
     procedure :: get_value_space => pressure_boundary_function_get_value_space
  end type pressure_boundary_function_t

  type mixed_laplacian_rt_analytical_functions_t
     private
     type(pressure_source_term_t)         :: pressure_source_term
     type(pressure_boundary_function_t)   :: pressure_boundary_function
   contains
     procedure :: get_pressure_source_term         => mlrtaf_get_pressure_source_term
     procedure :: get_pressure_boundary_function   => mlrtaf_get_pressure_boundary_function
  end type mixed_laplacian_rt_analytical_functions_t

  public :: mixed_laplacian_rt_analytical_functions_t

contains  

  !===============================================================================================
  subroutine pressure_source_term_get_value_space ( this, point, result )
    implicit none
    class(pressure_source_term_t), intent(in)    :: this
    type(point_t)                , intent(in)    :: point
    real(rp)                     , intent(inout) :: result
    result = 0.0_rp
  end subroutine pressure_source_term_get_value_space

  !===============================================================================================
  subroutine pressure_boundary_function_get_value_space ( this, point, result )
    implicit none
    class(pressure_boundary_function_t), intent(in)    :: this
    type(point_t)           , intent(in)    :: point
    real(rp)                , intent(inout) :: result
    result = point%get(1) + point%get(2)
  end subroutine pressure_boundary_function_get_value_space
  
  !===============================================================================================
  function mlrtaf_get_pressure_source_term ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mlrtaf_get_pressure_source_term
    mlrtaf_get_pressure_source_term => this%pressure_source_term
  end function mlrtaf_get_pressure_source_term
  
  !===============================================================================================
  function mlrtaf_get_pressure_boundary_function ( this )
    implicit none
    class(mixed_laplacian_rt_analytical_functions_t), target, intent(in)    :: this
    class(scalar_function_t), pointer :: mlrtaf_get_pressure_boundary_function
    mlrtaf_get_pressure_boundary_function => this%pressure_boundary_function
  end function mlrtaf_get_pressure_boundary_function
  
end module mixed_laplacian_rt_analytical_functions_names
!***************************************************************************************************


