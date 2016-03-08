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

  type test_nsi_iss_oss_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     character(len=:), allocatable :: default_structured_mesh
     character(len=:), allocatable :: default_velocity_order 
     character(len=:), allocatable :: default_pressure_order 
     character(len=:), allocatable :: default_viscosity            
     character(len=:), allocatable :: default_c1                   
     character(len=:), allocatable :: default_c2                   
     character(len=:), allocatable :: default_cc                   
     character(len=:), allocatable :: default_elemental_length_flag
     character(len=:), allocatable :: default_convection_activated 
     character(len=:), allocatable :: default_is_analytical
     character(len=:), allocatable :: default_is_initial
     character(len=:), allocatable :: default_is_temporal
     character(len=:), allocatable :: default_analytical_velocity_id
     character(len=:), allocatable :: default_analytical_pressure_id
     character(len=:), allocatable :: default_initial_time
   contains
     procedure, non_overridable :: set_default    => set_default_params
     !procedure :: set_analytical => set_default_analytical_params
     procedure, non_overridable :: add_to_cli     => add_params_to_cli
  end type test_nsi_iss_oss_params_t

  ! Types
  public :: test_nsi_iss_oss_params_t

contains

# include "sbm_nsi_iss_oss_cli.i90"
     
end module command_line_parameters_names

!***************************************************************************************************!
! ANALYTICAL FUNCTIONS                                                                              ! 
! Definition of the analytical functions for the NSI_ISS_OSS problem.                               !
!***************************************************************************************************!
module nsi_iss_oss_analytical_functions_names
  use serial_names
# include "debug.i90"
  implicit none
  private

  type, extends(vector_function_t) :: velocity_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => velocity_get_value_space_time
  end type velocity_function_t

  type, extends(vector_function_t) :: dt_velocity_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => dt_velocity_get_value_space_time
  end type dt_velocity_function_t

  type, extends(tensor_function_t) :: velocity_gradient_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => velocity_gradient_get_value_space_time
  end type velocity_gradient_function_t

  type, extends(vector_function_t) :: velocity_grad_div_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => velocity_grad_div_get_value_space_time
  end type velocity_grad_div_function_t

  type, extends(scalar_function_t) :: pressure_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => pressure_get_value_space_time
  end type pressure_function_t

  type, extends(vector_function_t) :: pressure_gradient_function_t
     integer(ip) :: swich
   contains
     procedure, non_overridable :: get_value_space_time => pressure_gradient_get_value_space_time
  end type pressure_gradient_function_t

  type nsi_iss_oss_analytical_functions_t
     type(velocity_function_t)          :: velocity
     type(dt_velocity_function_t)       :: dt_velocity
     type(velocity_gradient_function_t) :: velocity_gradient
     type(velocity_grad_div_function_t) :: velocity_grad_div
     type(pressure_function_t)          :: pressure
     type(pressure_gradient_function_t) :: pressure_gradient
     integer(ip)                        :: velocity_function_id
     integer(ip)                        :: pressure_function_id
   contains
     procedure, non_overridable :: initialize_from_cli
  end type nsi_iss_oss_analytical_functions_t

  ! Types
  public :: nsi_iss_oss_analytical_functions_t

contains

# include "sbm_nsi_iss_oss_analytical.i90"
  
end module nsi_iss_oss_analytical_functions_names

!***************************************************************************************************!
! DISCRETE INTEGRATION: NSI_ISS_OSS                                                                 ! 
! Navier-Stokes with Inf-Sup stable discretization using Orthogonal Subscales stabilization.        !
!***************************************************************************************************!
module nsi_iss_oss_discrete_integration_names
  use serial_names
  use command_line_parameters_names
  use nsi_iss_oss_analytical_functions_names
# include "debug.i90"
  implicit none
  private

  type, extends(discrete_integration_t) :: nsi_iss_oss_discrete_integration_t
     real(rp)                 :: viscosity
     real(rp)                 :: c1
     real(rp)                 :: c2
     real(rp)                 :: cc
     real(rp)                 :: current_time
     integer(ip)              :: elemental_length_flag
     logical                  :: convection_activated
     logical                  :: is_analytical_solution
     logical                  :: is_initial_solution
     logical                  :: is_temporal_solution
     class(vector_t), pointer :: dof_values => NULL() 
     type(nsi_iss_oss_analytical_functions_t) :: analytical_functions
   contains
     procedure                  :: integrate
     procedure, non_overridable :: initialize_from_cli
     procedure, non_overridable :: compute_stabilization_parameters
     procedure, non_overridable :: compute_characteristic_length
     procedure, non_overridable :: compute_mean_elemental_velocity
     procedure, non_overridable :: compute_analytical_force
     procedure, non_overridable :: update_boundary_conditions_analytical
  end type nsi_iss_oss_discrete_integration_t

  integer(ip), parameter :: characteristic_elemental_length = 0
  integer(ip), parameter :: minimum_elemental_length = 1
  integer(ip), parameter :: maximum_elemental_length = 2

  ! Types
  public :: nsi_iss_oss_discrete_integration_t

contains

# include "sbm_nsi_iss_oss_integrate.i90"
  
end module nsi_iss_oss_discrete_integration_names

!***************************************************************************************************!
! SERIAL DRIVER:                                                                                    !
!***************************************************************************************************!
program test_nsi_iss_oss
  use serial_names
  use command_line_parameters_names
  use nsi_iss_oss_discrete_integration_names
  implicit none
#include "debug.i90"

  ! Geometry
  type(mesh_t)          :: f_mesh
  type(conditions_t)    :: f_cond
  type(triangulation_t) :: f_trian

  ! Problem
  type(nsi_iss_oss_discrete_integration_t) :: nsi_iss_oss_integration
  
  ! Finite Elment
  type(serial_fe_space_t)              :: fe_space
  type(p_reference_fe_t)               :: reference_fe_array(3)
  type(fe_affine_operator_t)           :: fe_affine_operator
  type(vector_space_t), pointer        :: fe_affine_operator_range_vector_space
  class(vector_t), allocatable, target :: dof_values
  class(vector_t), allocatable, target :: residual

  ! Solver
  type(linear_solver_t)      :: linear_solver
  type(serial_environment_t) :: senv

  ! Arguments
  type(Type_Command_Line_Interface) :: cli
  character(len=:), allocatable     :: group
  character(len=256)                :: dir_path
  character(len=256)                :: prefix
  character(len=256)                :: dir_path_out
  logical                           :: is_structured_mesh
  integer(ip)                       :: velocity_order
  integer(ip)                       :: pressure_order

  ! Locals
  integer(ip) :: number_dimensions
  integer(ip) :: istat
  integer(ip) :: max_nonlinear_iterations
  integer(ip) :: counter
  real(rp)    :: nonlinear_tolerance
  real(rp)    :: residual_norm  
  
# include "sbm_nsi_iss_oss_driver.i90"
 
contains
  
  !==================================================================================================
  subroutine read_flap_cli_nsi_iss_oss(cli)
    implicit none
    type(Type_Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(test_nsi_iss_oss_params_t) :: test_params
    logical                         :: authors_print
    integer(ip)                     :: error

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'nsi_iss_oss_main',                                                 &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         & 'Serial FEMPAR driver to solve the Navier-Stokes Incompressible problem with Inf-Sup stable discretization using Orthogonal Subscales stabilization.',   &
         &        examples    = ['nsi_iss_oss_main            -h ',                                 &
                                 'nsi_iss_oss_main analytical -h ' ]) 

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical solution')

    ! Set Command Line Arguments for each group
    call test_params%set_default()
    call test_params%add_to_cli(cli,'analytical')

  end subroutine read_flap_cli_nsi_iss_oss

end program test_nsi_iss_oss
