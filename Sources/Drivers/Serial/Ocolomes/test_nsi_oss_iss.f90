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

!***************************************************************************************************!
! COMAND LINE PARAMETERS:                                                                           !
!***************************************************************************************************!
module command_line_parameters_names
  use types_names
  use Data_Type_Command_Line_Interface
# include "debug.i90"
  implicit none
  private

  type test_nsi_oss_iss_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
   contains
     procedure :: set_default    => set_default_params
     !procedure :: set_analytical => set_default_analytical_params
     !procedure :: add_to_cli     => add_params_to_cli
  end type test_nsi_oss_iss_params_t

  ! Types
  public :: test_nsi_oss_iss_params_t

contains

  !==================================================================================================
  subroutine set_default_params(this)
    implicit none
    class(test_nsi_oss_iss_params_t), intent(inout) :: this

    ! IO parameters
    this%default_dir_path     = 'data'
    this%default_prefix       = 'test'
    this%default_dir_path_out = 'output'

    !...

  end subroutine set_default_params
     
end module command_line_parameters_names

!***************************************************************************************************!
! DISCRETE INTEGRATION: NSI_ISS_OSS                                                                 ! 
! Navier-Stokes with Inf-Sup stable discretization using Orthogonal Subscales stabilization.        !
!***************************************************************************************************!
module nsi_iss_oss_discrete_integration_names
  use serial_names
# include "debug.i90"
  implicit none
  private

  type, extends(discrete_integration_t) :: nsi_iss_oss_discrete_integration_t
     real(rp)                 :: viscosity
     real(rp)                 :: c1
     real(rp)                 :: c2
     real(rp)                 :: cc
     integer(ip)              :: elemental_length_flag
     logical                  :: convection_activated
     class(vector_t), pointer :: dof_values => NULL() 
   contains
     procedure                  :: integrate
     !procedure, non_overridable :: initialize_from_cli
     procedure, non_overridable :: compute_stabilization_parameters
     procedure, non_overridable :: compute_characteristic_length
  end type nsi_iss_oss_discrete_integration_t

  integer(ip), parameter :: characteristic_elemental_length = 0
  integer(ip), parameter :: minimum_elemental_length = 1
  integer(ip), parameter :: maximum_elemental_length = 2

contains

# include "sbm_nsi_iss_oss_integrate.i90"
  
end module nsi_iss_oss_discrete_integration_names
