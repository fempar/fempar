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
module poisson_unfitted_dG_discrete_integration_names
  use fempar_names
  use poisson_unfitted_analytical_functions_names
  use poisson_unfitted_conditions_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_unfitted_dG_discrete_integration_t
     type(poisson_unfitted_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(poisson_unfitted_conditions_t)          , pointer :: poisson_unfitted_conditions   => NULL()
   contains
     procedure :: set_analytical_functions
     procedure :: set_poisson_unfitted_conditions
     procedure :: integrate_galerkin
  end type poisson_unfitted_dG_discrete_integration_t
  
  public :: poisson_unfitted_dG_discrete_integration_t
  
contains

  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_unfitted_dG_discrete_integration_t)        , intent(inout) :: this
     type(poisson_unfitted_analytical_functions_t)    , target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions
  
  subroutine set_poisson_unfitted_conditions ( this, poisson_unfitted_conditions )
     implicit none
     class(poisson_unfitted_dG_discrete_integration_t)        , intent(inout) :: this
     type(poisson_unfitted_conditions_t)              , target, intent(in)    :: poisson_unfitted_conditions
     this%poisson_unfitted_conditions => poisson_unfitted_conditions
  end subroutine set_poisson_unfitted_conditions
  
  
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(poisson_unfitted_dG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)         , intent(inout) :: assembler

    ! TO BE DONE
    
  end subroutine integrate_galerkin
  
end module poisson_unfitted_dG_discrete_integration_names
