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

!****************************************************************************************************

!> The program of the *Test Transient Poisson* docummented [here](../page/tests/test_transient_poisson.html)
!> which validates the *time stepping*

!* The `[[test_transient_poisson]]`,
program test_transient_poisson
!* uses the `fempar_names` and `test_transient_poisson_driver_names`:
  use fempar_names
  use test_transient_poisson_driver_names
  !* First, declare the `test_driver` and the `world_contest`  
  implicit none
  type(test_transient_poisson_driver_t) :: test_driver
  type(serial_context_t)      :: world_context
  !* Then, initialite
  call world_context%create()
  call fempar_init()
  !* and, setup according to the input parameters
  call test_driver%parse_command_line_parameters()
  call test_driver%setup_environment(world_context)
  !* Perform the core computations of the test, placed into [`run_simulation`](#run-simulation)
  call test_driver%run_simulation()
  !* Finally, free and finalize the used variables: 
  call test_driver%free_environment()
  call fempar_finalize()
  call world_context%free(finalize=.true.)
end program test_transient_poisson
