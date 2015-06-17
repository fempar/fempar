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
  use types
  use memor
  use array_names
  use fem_space_types
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
      procedure(create_interface) , deferred :: create
      procedure(compute_interface), deferred :: compute
      procedure :: free => discrete_problem_free
   end type discrete_problem

   type :: discrete_problem_pointer
      class(discrete_problem), pointer :: p
   end type discrete_problem_pointer

   type, abstract :: discrete_data
      contains
        procedure(create_data_interface), deferred :: create
   end type discrete_data

  abstract interface
     subroutine create_interface(approx,prob,data,l2g)
       import :: physical_problem, discrete_problem, discrete_data, ip
       implicit none
       class(discrete_problem)        , intent(out) :: approx
       class(physical_problem), target, intent(in)  :: prob
       class(discrete_data)   , target, intent(in)  :: data
       integer(ip), intent(in), optional :: l2g(:)
     end subroutine create_interface
     subroutine compute_interface(approx,start,elem)
       import :: discrete_problem, fem_element, ip
       implicit none
       class(discrete_problem), intent(inout) :: approx
       integer(ip)            , intent(in)    :: start(:)
       type(fem_element)      , intent(inout) :: elem
     end subroutine compute_interface
     subroutine free_interface(approx)
       import :: discrete_problem
       implicit none
       class(discrete_problem), intent(inout) :: approx
     end subroutine free_interface
     subroutine create_data_interface(data)
       import :: discrete_data
       implicit none
       class(discrete_data), intent(out) :: data
     end subroutine create_data_interface
  end interface

  public :: physical_problem, p_physical_problem, discrete_problem, &
            discrete_problem_pointer, discrete_problem_free, discrete_data

contains 

  subroutine discrete_problem_free( prob  )
    implicit none
    class(discrete_problem), intent(inout) :: prob
    if(allocated(prob%l2g_var)) call memfree( prob%l2g_var, __FILE__, __LINE__ )
  end subroutine discrete_problem_free

end module problem_names
