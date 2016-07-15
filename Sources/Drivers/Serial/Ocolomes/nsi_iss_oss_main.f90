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
  use serial_names
# include "debug.i90"
  implicit none
  private

  type test_nsi_iss_oss_params_t
     private
     ! Defaults
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     ! Mesh
     character(len=:), allocatable :: default_structured_mesh        
     character(len=:), allocatable :: default_number_elements        
     character(len=:), allocatable :: default_number_parts           
     character(len=:), allocatable :: default_number_sockets         
     character(len=:), allocatable :: default_discretization_type    
     character(len=:), allocatable :: default_domain_length          
     character(len=:), allocatable :: default_periodic_boundaries    
     character(len=:), allocatable :: default_origin                 
     character(len=:), allocatable :: default_stretching_parameter   
     character(len=:), allocatable :: default_elements_boundary_layer
     character(len=:), allocatable :: default_size_boundary_layer    
     character(len=:), allocatable :: default_material_case          
     ! FE space
     character(len=:), allocatable :: default_velocity_order 
     character(len=:), allocatable :: default_pressure_order 
     ! Problem
     character(len=:), allocatable :: default_viscosity            
     character(len=:), allocatable :: default_c1                   
     character(len=:), allocatable :: default_c2                   
     character(len=:), allocatable :: default_cc                   
     character(len=:), allocatable :: default_elemental_length_flag
     character(len=:), allocatable :: default_convection_activated 
     ! Solution
     character(len=:), allocatable :: default_is_analytical
     character(len=:), allocatable :: default_is_initial
     character(len=:), allocatable :: default_is_temporal
     character(len=:), allocatable :: default_analytical_function_name
     ! Time integration
     character(len=:), allocatable :: default_initial_time
     ! CLI
     type(Command_Line_Interface)  :: cli
     ! Parameters
     ! IO parameters
     character(len=2014)      :: dir_path
     character(len=2014)      :: prefix
     character(len=2014)      :: dir_path_out
     ! Mesh
     logical                  :: is_structured_mesh
     integer(ip), allocatable :: number_elements(:)
     integer(ip), allocatable :: number_parts(:)
     integer(ip), allocatable :: number_sockets(:)
     integer(ip), allocatable :: discretization_type(:)
     integer(ip), allocatable :: periodic_boundaries(:)
     integer(ip), allocatable :: number_elements_boundary(:)
     integer(ip)              :: material_case
     real(rp)   , allocatable :: domain_length(:)
     real(rp)   , allocatable :: origin(:)
     real(rp)   , allocatable :: stretching_parameter(:)
     real(rp)   , allocatable :: size_boundary(:) 
     ! FE space
     integer(ip)              :: velocity_order
     integer(ip)              :: pressure_order   
     ! Problem 
     real(rp)                 :: viscosity
     real(rp)                 :: c1
     real(rp)                 :: c2
     real(rp)                 :: cc
     integer(ip)              :: elemental_length_flag
     logical                  :: convection_activated
     ! Solution
     logical                  :: is_analytical_solution
     logical                  :: is_initial_solution
     logical                  :: is_temporal_solution
     character(len=2014)      :: analytical_function_name
     ! Time integration
     real(rp)                 :: initial_time
   contains
     procedure, non_overridable          :: create
     procedure, non_overridable, private :: set_default
     procedure, non_overridable, private :: apply_default
     procedure, non_overridable, private :: add_to_cli 
     procedure, non_overridable          :: parse
     procedure, non_overridable          :: free  
     procedure, non_overridable          :: create_parameter_list_for_structured_mesh
     ! Getters
     procedure, non_overridable :: get_dir_path
     procedure, non_overridable :: get_prefix
     procedure, non_overridable :: get_dir_path_out
     procedure, non_overridable :: get_is_structured_mesh
     procedure, non_overridable :: get_number_elements
     procedure, non_overridable :: get_number_parts
     procedure, non_overridable :: get_number_sockets
     procedure, non_overridable :: get_discretization_type
     procedure, non_overridable :: get_periodic_boundaries
     procedure, non_overridable :: get_number_elements_boundary
     procedure, non_overridable :: get_material_case
     procedure, non_overridable :: get_domain_length
     procedure, non_overridable :: get_origin
     procedure, non_overridable :: get_stretching_parameter
     procedure, non_overridable :: get_size_boundary
     procedure, non_overridable :: get_velocity_order
     procedure, non_overridable :: get_pressure_order    
     procedure, non_overridable :: get_viscosity
     procedure, non_overridable :: get_c1
     procedure, non_overridable :: get_c2
     procedure, non_overridable :: get_cc
     procedure, non_overridable :: get_elemental_length_flag
     procedure, non_overridable :: get_convection_activated
     procedure, non_overridable :: get_is_analytical_solution
     procedure, non_overridable :: get_is_initial_solution
     procedure, non_overridable :: get_is_temporal_solution
     procedure, non_overridable :: get_analytical_function_name
     procedure, non_overridable :: get_initial_time     
  end type test_nsi_iss_oss_params_t

  ! Types
  public :: test_nsi_iss_oss_params_t

contains

# include "sbm_nsi_iss_oss_cli.i90"
     
end module command_line_parameters_names

!***************************************************************************************************!
! ANALYTICAL LINEAR-STEADY FUNCTION                                                                 !
! Definition of the linear-steady analytical function for the NSI_ISS_OSS problem.                  !
!***************************************************************************************************!
module nsi_iss_oss_analytical_linear_steady_functions_names
  use serial_names
# include "debug.i90"
  implicit none
  private

  ! u = (x,-y,0)
  ! p = (x+y)

  type, extends(vector_function_t) :: linear_steady_velocity_function_t
     integer(ip) :: random = 3
   contains
     procedure :: get_value_space_time => linear_steady_velocity_get_value_space_time
  end type linear_steady_velocity_function_t

  type, extends(vector_function_t) :: linear_steady_dt_velocity_function_t
   contains
     procedure :: get_value_space_time => linear_steady_dt_velocity_get_value_space_time
  end type linear_steady_dt_velocity_function_t

  type, extends(tensor_function_t) :: linear_steady_velocity_gradient_function_t
   contains
     procedure :: get_value_space_time => linear_steady_velocity_gradient_get_value_space_time
  end type linear_steady_velocity_gradient_function_t

  type, extends(vector_function_t) :: linear_steady_velocity_grad_div_function_t
   contains
     procedure :: get_value_space_time => linear_steady_velocity_grad_div_get_value_space_time
  end type linear_steady_velocity_grad_div_function_t

  type, extends(scalar_function_t) :: linear_steady_pressure_function_t
   contains
     procedure :: get_value_space_time => linear_steady_pressure_get_value_space_time
  end type linear_steady_pressure_function_t

  type, extends(vector_function_t) :: linear_steady_pressure_gradient_function_t
   contains
     procedure :: get_value_space_time => linear_steady_pressure_gradient_get_value_space_time
  end type linear_steady_pressure_gradient_function_t

  public :: linear_steady_velocity_function_t
  public :: linear_steady_dt_velocity_function_t
  public :: linear_steady_velocity_gradient_function_t
  public :: linear_steady_velocity_grad_div_function_t
  public :: linear_steady_pressure_gradient_function_t
  public :: linear_steady_pressure_function_t

contains

# include "sbm_nsi_iss_oss_analytical_linear_steady.i90"

end module nsi_iss_oss_analytical_linear_steady_functions_names

!***************************************************************************************************!
! ANALYTICAL FUNCTIONS                                                                              ! 
! Definition of the analytical functions for the NSI_ISS_OSS problem.                               !
!***************************************************************************************************!
module nsi_iss_oss_analytical_functions_names
  use serial_names
  use command_line_parameters_names
  use nsi_iss_oss_analytical_linear_steady_functions_names
# include "debug.i90"
  implicit none
  private

  type nsi_iss_oss_analytical_functions_t
     private
     class(vector_function_t), allocatable :: velocity
     class(vector_function_t), allocatable :: dt_velocity
     class(tensor_function_t), allocatable :: velocity_gradient
     class(vector_function_t), allocatable :: velocity_grad_div
     class(scalar_function_t), allocatable :: pressure
     class(vector_function_t), allocatable :: pressure_gradient
     character(len=:)        , allocatable :: function_name
   contains
     procedure, non_overridable :: create                      => nsi_iss_oss_analytical_functions_create
     procedure, non_overridable :: free                        => nsi_iss_oss_analytical_functions_free
     procedure, non_overridable :: get_velocity_function       => nsi_iss_oss_analytical_functions_get_velocity_function
     procedure, non_overridable :: get_pressure_function       => nsi_iss_oss_analytical_functions_get_pressure_function
     procedure, non_overridable :: get_value_velocity          => nsi_iss_oss_analytical_functions_get_value_velocity
     procedure, non_overridable :: get_value_dt_velocity       => nsi_iss_oss_analytical_functions_get_value_dt_velocity
     procedure, non_overridable :: get_value_velocity_gradient => nsi_iss_oss_analytical_functions_get_value_velocity_gradient
     procedure, non_overridable :: get_value_velocity_grad_div => nsi_iss_oss_analytical_functions_get_value_velocity_grad_div
     procedure, non_overridable :: get_value_pressure_gradient => nsi_iss_oss_analytical_functions_get_value_pressure_gradient
     procedure, non_overridable :: get_value_pressure          => nsi_iss_oss_analytical_functions_get_value_pressure
  end type nsi_iss_oss_analytical_functions_t

  character(len=*), parameter :: nsi_linear_steady = 'NSI-LINEAR-STEADY'

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
     real(rp)                      :: viscosity
     real(rp)                      :: c1
     real(rp)                      :: c2
     real(rp)                      :: cc
     real(rp)                      :: current_time
     integer(ip)                   :: elemental_length_flag
     logical                       :: convection_activated
     logical                       :: is_analytical_solution
     logical                       :: is_initial_solution
     logical                       :: is_temporal_solution
     class(fe_function_t), pointer :: fe_values => NULL() 
     type(nsi_iss_oss_analytical_functions_t) :: analytical_functions
   contains
     procedure                  :: integrate
     procedure, non_overridable :: create
     procedure, non_overridable :: compute_stabilization_parameters
     procedure, non_overridable :: compute_characteristic_length
     procedure, non_overridable :: compute_mean_elemental_velocity
     procedure, non_overridable :: compute_analytical_force
     procedure, non_overridable :: update_boundary_conditions_analytical
     procedure, non_overridable :: interpolate_fe_function_analytical
     procedure, non_overridable :: compute_error_norms
     procedure, non_overridable :: free
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
  type(mesh_t)                          :: f_mesh
  type(conditions_t)                    :: f_cond
  type(triangulation_t)                 :: f_trian
  type(uniform_mesh_descriptor_t)       :: geometry_data
  type(uniform_conditions_descriptor_t) :: boundary_data

  ! Problem
  type(nsi_iss_oss_discrete_integration_t) :: nsi_iss_oss_integration
  
  ! Finite Elment
  type(serial_fe_space_t)              :: fe_space
  type(p_reference_fe_t)               :: reference_fe_array(3)
  type(fe_affine_operator_t)           :: fe_affine_operator
  type(fe_function_t)         , target :: fe_values
  class(vector_t), allocatable, target :: residual
  class(vector_t)            , pointer :: dof_values => null()

  ! Solver
  type(iterative_linear_solver_t) :: iterative_linear_solver
  type(serial_environment_t)      :: senv
  type(parameterlist_t)           :: solver_parameters

  ! Arguments
  type(test_nsi_iss_oss_params_t) :: params
  type(parameterlist_t)           :: mesh_parameters

  ! Locals
  integer(ip) :: number_dimensions,number_unknowns
  integer(ip) :: istat
  integer(ip) :: max_nonlinear_iterations
  integer(ip) :: counter
  real(rp)    :: nonlinear_tolerance
  real(rp)    :: residual_norm  

# include "sbm_nsi_iss_oss_driver.i90"
 
contains
  


  !==================================================================================================
  subroutine set_structured_conditions_nsi_iss_oss(geometry_data,boundary_data)
    implicit none
    type(uniform_mesh_descriptor_t)      , intent(in)    :: geometry_data
    type(uniform_conditions_descriptor_t), intent(inout) :: boundary_data
    integer(ip) :: number_dimensions
    
    number_dimensions = geometry_data%ndime
    
    ! Pressure codes
    boundary_data%poin%code(number_dimensions+1,1:2**number_dimensions-1) = 0
    boundary_data%poin%code(number_dimensions+1,2**number_dimensions)     = 1
    boundary_data%line%code(number_dimensions+1,:)                        = 0
    boundary_data%surf%code(number_dimensions+1,:)                        = 0
    ! Projection codes
    boundary_data%poin%code(number_dimensions+2:2*number_dimensions+1,:)  = 0
    boundary_data%line%code(number_dimensions+2:2*number_dimensions+1,:)  = 0
    boundary_data%surf%code(number_dimensions+2:2*number_dimensions+1,:)  = 0
    ! Velocity values
    boundary_data%poin%valu(1:number_dimensions,:) = 1.0_rp
    boundary_data%line%valu(1:number_dimensions,:) = 1.0_rp
    ! Pressure values
    boundary_data%poin%valu(number_dimensions+1,2**number_dimensions) = 0.0_rp

  end subroutine set_structured_conditions_nsi_iss_oss  

end program test_nsi_iss_oss
