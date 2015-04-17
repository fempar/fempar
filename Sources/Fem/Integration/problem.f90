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
  implicit none
  private

  type :: physical_problem
     integer(ip)        ::             &
          nvars,                       &       ! Number of different variables
          nunks,                       &       ! Number of unknowns (groups of variables)
          ndime,                       &       ! Number of space dimensions
          ntens                                ! Number of tensor components
     integer(ip), allocatable ::       &
          vars_of_unk(:),              &       ! Number of variables of each unknown (size nunks)
          l2g_var(:)                           ! Order chosen for variables (size nvars)
     !integer(ip)        ::             &
     !     problem_code                         ! An internal code that defines a problem in FEMPAR
     character(len=:), allocatable ::  &
          unkno_names(:)                       ! Names for the gauss_properties (nunks)
  end type physical_problem


  type, abstract :: discrete_problem
   contains
      procedure(create_interface), deferred :: create
      procedure(matvec_interface), deferred :: matvec
   end type discrete_problem

   type :: discrete_problem_pointer
      class(discrete_problem), pointer :: p
   end type discrete_problem_pointer

  abstract interface
     subroutine create_interface(aprox,prob)
       import :: physical_problem, discrete_problem
       implicit none
       class(discrete_problem)        , intent(out) :: aprox
       class(physical_problem), target, intent(in)  :: prob
     end subroutine create_interface
     subroutine matvec_interface(aprox,integ,unkno,start,mat,vec)
       import :: discrete_problem, volume_integrator_pointer, array_rp2, array_rp1, rp, ip
       implicit none
       class(discrete_problem)        , intent(in) :: aprox
       type(volume_integrator_pointer), intent(in) :: integ(:)
       real(rp)                       , intent(in) :: unkno(:,:,:)
       integer(ip)                    , intent(in) :: start(:)
       type(array_rp2), intent(inout) :: mat
       type(array_rp1), intent(inout) :: vec
     end subroutine matvec_interface
  end interface

  public :: physical_problem, discrete_problem, discrete_problem_pointer

end module problem_names
