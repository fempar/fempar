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
  use finite_element_names
  implicit none
  private

  type :: physical_problem_t
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
  end type physical_problem_t

  type :: p_physical_problem_t
     type(physical_problem_t), pointer :: p 
  end type p_physical_problem_t

  type, abstract :: discrete_problem_t
     integer(ip)        ::             &
          nvars                                ! Number of discrete variables
     integer(ip), allocatable ::       &
          l2g_var(:)                           ! Order chosen for variables (size nvars)
   contains
     procedure(create_problem_interface), deferred :: create
     procedure :: free => discrete_problem_free
  end type discrete_problem_t

   type :: discrete_problem_pointer_t
      class(discrete_problem_t), pointer :: p
   end type discrete_problem_pointer_t

  type, abstract :: discrete_integration_t
     integer(ip) :: integration_stage = update_nonlinear
     integer(ip), allocatable :: working_vars(:)
    contains
      procedure(create_integration_interface) , deferred :: create 
      procedure(compute_integration_interface), deferred :: compute
      procedure(free_integration_interface)   , deferred :: free 
   end type discrete_integration_t

   type :: discrete_integration_pointer_t
      class(discrete_integration_t), pointer :: p
   end type discrete_integration_pointer_t

  abstract interface
     subroutine create_problem_interface(discret,physics,l2g)
       import :: physical_problem_t, discrete_problem_t, ip
       implicit none
       class(discrete_problem_t), intent(out) :: discret
       class(physical_problem_t), intent(in)  :: physics
       integer(ip), optional  , intent(in)  :: l2g(:)
     end subroutine create_problem_interface
     subroutine create_integration_interface( approx, physics, discret )
       import :: discrete_integration_t, discrete_problem_t, physical_problem_t
       implicit none
       class(discrete_integration_t)   , intent(inout) :: approx
       class(physical_problem_t), target, intent(in)    :: physics
       class(discrete_problem_t), target, intent(in)    :: discret
     end subroutine create_integration_interface
     subroutine compute_integration_interface(approx,finite_element)
       import :: discrete_integration_t, finite_element_t
       implicit none
       class(discrete_integration_t), intent(inout) :: approx
       type(finite_element_t)       , intent(inout) :: finite_element
     end subroutine compute_integration_interface
     subroutine free_integration_interface(approx)
       import :: discrete_integration_t
       implicit none
       class(discrete_integration_t), intent(inout) :: approx
     end subroutine free_integration_interface
  end interface

  public :: physical_problem_t, p_physical_problem_t, discrete_problem_t, &
            discrete_problem_pointer_t, discrete_problem_free,            &
            discrete_integration_t, discrete_integration_pointer_t

contains 

  subroutine discrete_problem_free( prob  )
    implicit none
    class(discrete_problem_t), intent(inout) :: prob
    if(allocated(prob%l2g_var)) call memfree( prob%l2g_var, __FILE__, __LINE__ )
  end subroutine discrete_problem_free

end module problem_names
