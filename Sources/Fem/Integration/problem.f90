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
module problem_names
  use types_names
  use memor_names
  use array_names
  use integration_tools_names
  use fem_element_names
  implicit none
  private

  type :: physical_problem
     integer(ip)        ::             &
          nvars,                       &       ! Number of variables
          nunks,                       &       ! Number of unknowns (groups of variables)
          ndime,                       &       ! Number of space dimensions
          ntens                                ! Number of tensor components
     integer(ip), allocatable ::       &
          vars_of_unk(:)                       ! Number of variables of each unknown (size nunks)
     character(len=:), allocatable ::  &
          unkno_names(:)                       ! Names for the gauss_properties (nunks)
   ! contains
   !    procedure(create_interface), deferred :: create
   !    procedure(matvec_interface), deferred :: matvec
   !    procedure(free_interface), deferred :: free
  end type physical_problem

  type :: p_physical_problem
     type(physical_problem), pointer :: p 
  end type p_physical_problem

  type, abstract :: discrete_problem
     integer(ip)        ::             &
          nvars                                ! Number of discrete variables
     integer(ip), allocatable ::       &
          l2g_var(:)                           ! Order chosen for variables (size nvars)
   contains
     procedure(create_problem_interface), deferred :: create
     procedure :: free => discrete_problem_free
  end type discrete_problem

   type :: discrete_problem_pointer
      class(discrete_problem), pointer :: p
   end type discrete_problem_pointer

  type, abstract :: discrete_integration
    contains
      procedure(create_integration_interface) , deferred :: create 
      procedure(compute_integration_interface), deferred :: compute
      procedure(free_integration_interface)   , deferred :: free 
   end type discrete_integration

   type :: discrete_integration_pointer
      class(discrete_integration), pointer :: p
   end type discrete_integration_pointer

  abstract interface
     subroutine create_problem_interface(discret,physics,l2g)
       import :: physical_problem, discrete_problem, ip
       implicit none
       class(discrete_problem), intent(out) :: discret
       class(physical_problem), intent(in)  :: physics
       integer(ip), optional  , intent(in)  :: l2g(:)
     end subroutine create_problem_interface
     subroutine create_integration_interface( approx, physics, discret )
       import :: discrete_integration, discrete_problem, physical_problem
       implicit none
       class(discrete_integration)   , intent(inout) :: approx
       class(physical_problem), target, intent(in)    :: physics
       class(discrete_problem), target, intent(in)    :: discret
     end subroutine create_integration_interface
     subroutine compute_integration_interface(approx,start,elem)
       import :: discrete_integration, fem_element_t, ip
       implicit none
       class(discrete_integration), intent(inout) :: approx
       integer(ip)                , intent(in)    :: start(:)
       type(fem_element_t)          , intent(inout) :: elem
     end subroutine compute_integration_interface
     subroutine free_integration_interface(approx)
       import :: discrete_integration
       implicit none
       class(discrete_integration), intent(inout) :: approx
     end subroutine free_integration_interface
  end interface

  public :: physical_problem, p_physical_problem, discrete_problem, &
            discrete_problem_pointer, discrete_problem_free,        &
            discrete_integration, discrete_integration_pointer

contains 

  subroutine discrete_problem_free( prob  )
    implicit none
    class(discrete_problem), intent(inout) :: prob
    if(allocated(prob%l2g_var)) call memfree( prob%l2g_var, __FILE__, __LINE__ )
  end subroutine discrete_problem_free

end module problem_names
