!> date:   15-November-2016                                                                              
!> author: Oriol Colomes
!> summary:  NSI driver parmeters module
!>
!> This module contains the data types and subroutines to define the test_nsi parameters.
!>
!> The command line parameter parsing is based on FLAP, a library linked to FEMPAR.

module test_nsi_parameters_names
  use fempar_names
  use time_integrator_names
  use theta_method_names
  use nsi_keys_names
  implicit none
#include "debug.i90"
  private

  ! Types
  !> summary: NSI driver parameters
  !>
  !> This type contains the parameters needed by the test_nsi. It extends the class [[parameter_generator_t]]
  !> and defines the default values depending on the CLI execution group.
  !>
  !> Currently, the following CLI group are implemented:
  !>   - analytical (Default values are defined assuming 'analytical' group)
  type, extends(parameter_handler_t), public :: test_nsi_parameters_t
     private
   contains 
     ! Setters
     procedure                  :: define_parameters                 !< Set parameter generator lists
     procedure, non_overridable :: set_common_lists                  !< Set default parameter lists
     procedure, non_overridable :: set_output_group_lists            !< Set output parameter lists
     procedure, non_overridable :: set_mesh_group_lists              !< Set mesh parameter lists
     procedure, non_overridable :: set_problem_group_lists           !< Set physical problem parameter lists
     procedure, non_overridable :: set_space_group_lists             !< Set FE space parameter lists
     procedure, non_overridable :: set_formulation_group_lists       !< Set FE formulation parameter lists
     procedure, non_overridable :: set_solution_group_lists          !< Set FE solution parameter lists
     procedure, non_overridable :: set_time_integration_group_lists  !< Set time integration parameter lists
     procedure, non_overridable :: set_linear_solver_group_lists     !< Set linear solver parameter lists
     procedure, non_overridable :: set_nonlinear_solver_group_lists  !< Set nonlinear solver parameter lists
     ! Getters                                              
     procedure, non_overridable :: get_dir_path              !< Get INPUT directory path                     
     procedure, non_overridable :: get_triangulation_type    !< Get type of triangulation generator (from mesh or from structured)      
     procedure, non_overridable :: get_velocity_order        !< Get interpolation order for the velocity        
     procedure, non_overridable :: get_pressure_order        !< Get interpolation order for the pressure      
     procedure, non_overridable :: get_vms_method_type       !< Get VMS method type                    
     procedure, non_overridable :: get_time_integrator_type  !< Get type of the time integrator method        
     procedure, non_overridable :: get_solver_type           !< Get linear solver type                                           
     procedure, non_overridable :: get_solver_rtol           !< Get linear solver residual tolerance                
     procedure, non_overridable :: get_solver_trace          !< Get linear solver output trace                   
     procedure, non_overridable :: get_nonlinear_solver_type !< Get nonlinear solver type                         
  end type test_nsi_parameters_t

contains

  !==================================================================================================
  !> summary: Set parameter generator for the NSI problem
  !>
  !> This subroutine sets the default values of the ParameterLists_t's needed for the parameter
  !> generator:
  !>   - values      (list with parameter values)
  !>   - switches    (list with parameter names)
  !>   - switches_ab (list with parameter shortcut names)
  !>   - helpers     (list with help messages for each parameter)
  !>   - required    (list of logicals for the required/non-required parameters)
  !>
  !> The parameters are grouped in different sublists, depending their purpose. Only the input data
  !> (input path and prefix) is stored in the first level lists. The groups implemented in this
  !> subroutine are:
  !>   - Ouput parameters
  !>   - Mesh generation parameters
  !>   - Physical problem parameters
  !>   - FE space parameters
  !>   - FE formulation parameters
  !>   - FE solution parameters
  !>   - Time integration parameters
  !>   - Linear solver parameters
  !>   - Nonlinear solver parameters
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine define_parameters(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this

    ! Set default lists
    call this%set_common_lists()

    ! Set parameter lists for each CLI group
    call this%set_output_group_lists()          
    call this%set_mesh_group_lists()            
    call this%set_problem_group_lists()         
    call this%set_space_group_lists()           
    call this%set_formulation_group_lists()     
    call this%set_solution_group_lists()        
    call this%set_time_integration_group_lists()
    call this%set_linear_solver_group_lists()   
    call this%set_nonlinear_solver_group_lists()

  end subroutine define_parameters

  !==================================================================================================
  !> summary: Set default common parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_common_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    integer(ip)                    :: error = 0

    ! Point lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()
    
    ! Values
    ! ------
    error = values%set(key = dir_path_key, value = 'data'  ) ; check(error==0)
    error = values%set(key = prefix_key  , value = 'test'  ) ; check(error==0)
    
    ! Switches
    ! --------
    error = switches%set(key = dir_path_key, value = '--dir_path'); check(error==0)
    error = switches%set(key = prefix_key  , value = '--prefix')  ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab%set(key = dir_path_key, value = '-d'); check(error==0)
    error = switches_ab%set(key = prefix_key  , value = '-p'); check(error==0)

    ! Helpers
    ! -------
    error = helpers%set(key = dir_path_key, value = 'Directory of the source files'); check(error==0)
    error = helpers%set(key = prefix_key  , value = 'Prefix')                       ; check(error==0)

    ! Required
    ! --------
    error = required%set(key = dir_path_key, value = .false.); check(error==0)
    error = required%set(key = prefix_key  , value = .false.); check(error==0)

  end subroutine set_common_lists

  !==================================================================================================
  !> summary: Set default ouput parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_output_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(output_group_key)
    switches_ab_group => switches_ab%NewSubList(output_group_key)
    helpers_group     => helpers%NewSubList(output_group_key)
    required_group    => required%NewSubList(output_group_key)
    values_group      => values%NewSubList(output_group_key)

    ! Values
    ! ------
    error = values_group%set(key = dir_path_out_key    , value = 'output') ; check(error==0)
    error = values_group%set(key = write_solution_key  , value = .false. ) ; check(error==0)
    error = values_group%set(key = write_error_norm_key, value = .false. ) ; check(error==0)
    error = values_group%set(key = error_norm_type_key , value = l2_norm ) ; check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = dir_path_out_key    , value = '--dir_path_out')   ; check(error==0)
    error = switches_group%set(key = write_solution_key  , value = '--write-solution') ; check(error==0)
    error = switches_group%set(key = write_error_norm_key, value = '--write-error')    ; check(error==0)
    error = switches_group%set(key = error_norm_type_key , value = '--error-norm')     ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = dir_path_out_key    , value = '-o')        ; check(error==0)
    error = switches_ab_group%set(key = write_solution_key  , value = '-wsolution'); check(error==0)
    error = switches_ab_group%set(key = write_error_norm_key, value = '-werror')   ; check(error==0)
    error = switches_ab_group%set(key = error_norm_type_key , value = '-enorm')    ; check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = dir_path_out_key    , value = 'Output directory'); check(error==0)
    error = helpers_group%set(key = write_solution_key  , value = 'Write solution')  ; check(error==0)
    error = helpers_group%set(key = write_error_norm_key, value = 'Write error norm'); check(error==0)
    error = helpers_group%set(key = error_norm_type_key , value = 'Error norm type: mean_norm, l1_norm, l2_norm, lp_norm, linfty_norm, h1_seminorm, hdiv_seminorm, hcurl_seminorm, h1_norm, w1p_seminorm, w1p_norm, w1infty_seminorm, w1infty_norm') ; check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = dir_path_out_key    , value = .false.); check(error==0)
    error = required_group%set(key = write_solution_key  , value = .false.); check(error==0)
    error = required_group%set(key = write_error_norm_key, value = .false.); check(error==0)
    error = required_group%set(key = error_norm_type_key , value = .false.); check(error==0)

  end subroutine set_output_group_lists

  !==================================================================================================
  !> summary: Set default mesh parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_mesh_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0
    character(len=:), allocatable  :: msg

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(mesh_group_key)
    switches_ab_group => switches_ab%NewSubList(mesh_group_key)
    helpers_group     => helpers%NewSubList(mesh_group_key)
    required_group    => required%NewSubList(mesh_group_key)
    values_group      => values%NewSubList(mesh_group_key)

    ! Values
    ! ------
    ! Structured/non-structured
    error = values_group%set(key = triangulation_generate_key, value = triangulation_generate_structured); check(error==0)
    ! Domain definition
    error = values_group%set(key = domain_length_key   , value = [1.0_rp,1.0_rp,1.0_rp]);                   check(error==0)
    error = values_group%set(key = domain_origin_key   , value = [0.0_rp,0.0_rp,0.0_rp]);                   check(error==0)
    error = values_group%set(key = is_dir_periodic_key , value = [non_periodic,non_periodic,non_periodic]); check(error==0)
    ! Mesh characterization
    error = values_group%set(key = reference_fe_geo_order_key   , value = 1);                         check(error==0)
    error = values_group%set(key = num_dims_key     , value = 2);                         check(error==0)
    error = values_group%set(key = num_cells_per_dir_key  , value = [2,2,2]);                   check(error==0)
    error = values_group%set(key = discretization_type_key      , value = [uniform,uniform,uniform]); check(error==0)
    error = values_group%set(key = mesh_stretching_parameter_key, value = [2.75_rp,2.75_rp,2.75_rp]); check(error==0)
    ! Mesh parallelization
    error = values_group%set(key = num_levels_key       , value = 1);       check(error==0)
    error = values_group%set(key = num_parts_per_dir_key, value = [1,1,0]); check(error==0)
    
    ! Switches
    ! --------
    ! Structured/non-structured
    error = switches_group%set(key = triangulation_generate_key, value = '--triangulation-type'); check(error==0)
    ! Domain definition
    error = switches_group%set(key = domain_length_key   , value = '--domain_length'); check(error==0)
    error = switches_group%set(key = domain_origin_key   , value = '--domain_origin'); check(error==0)
    error = switches_group%set(key = is_dir_periodic_key , value = '--periodicity')  ; check(error==0)
    ! Mesh characterization
    error = switches_group%set(key = reference_fe_geo_order_key   , value = '--reference-fe-geo-order'); check(error==0)
    error = switches_group%set(key = num_dims_key     , value = '--dim')                   ; check(error==0)
    error = switches_group%set(key = num_cells_per_dir_key  , value = '--num_cells')       ; check(error==0)
    error = switches_group%set(key = discretization_type_key      , value = '--discretization_type')   ; check(error==0)
    error = switches_group%set(key = mesh_stretching_parameter_key, value = '--stretching')            ; check(error==0)
    ! Mesh parallelization
    error = switches_group%set(key = num_levels_key       , value = '--num_levels'); check(error==0)
    error = switches_group%set(key = num_parts_per_dir_key, value = '--num_parts') ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    ! Structured/non-structured
    error = switches_ab_group%set(key = triangulation_generate_key, value = '-tt'); check(error==0)
    ! Domain definition
    error = switches_ab_group%set(key = domain_length_key   , value = '-dl')    ; check(error==0)
    error = switches_ab_group%set(key = domain_origin_key   , value = '-do')    ; check(error==0)
    error = switches_ab_group%set(key = is_dir_periodic_key , value = '-period'); check(error==0)
    ! Mesh characterization
    error = switches_ab_group%set(key = reference_fe_geo_order_key   , value = '-gorder'); check(error==0)
    error = switches_ab_group%set(key = num_dims_key     , value = '-dm')    ; check(error==0)
    error = switches_ab_group%set(key = num_cells_per_dir_key  , value = '-n')     ; check(error==0)
    error = switches_ab_group%set(key = discretization_type_key      , value = '-dt')    ; check(error==0)
    error = switches_ab_group%set(key = mesh_stretching_parameter_key, value = '-st')    ; check(error==0)
    ! Mesh parallelization
    error = switches_ab_group%set(key = num_levels_key       , value = '-l') ; check(error==0)
    error = switches_ab_group%set(key = num_parts_per_dir_key, value = '-np'); check(error==0)

    ! Helpers
    ! -------
    ! Structured/non-structured
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers_group%set(key = triangulation_generate_key, value = msg); check(error==0)
    ! Domain definition
    error = helpers_group%set(key = domain_length_key  , value = 'Domain length per direction')     ; check(error==0)
    error = helpers_group%set(key = domain_origin_key  , value = 'Domain origin coordinates')       ; check(error==0)
    error = helpers_group%set(key = is_dir_periodic_key, value = 'Domain periodicity per direction'); check(error==0)
    ! Mesh characterization
    error = helpers_group%set(key = reference_fe_geo_order_key , value = 'Order of the triangulation reference Finite Element'); check(error==0)
    error = helpers_group%set(key = num_dims_key   , value = 'Number of space dimensions')                         ; check(error==0)
    error = helpers_group%set(key = num_cells_per_dir_key, value = 'Number of cells per direction')                      ; check(error==0)
    msg = 'Discretization type: (*) uniform, (*) cubic, (*) tanh'
    write(msg(23:23),'(i1)') uniform
    write(msg(36:36),'(i1)') cubic
    write(msg(47:47),'(i1)') tanh
    error = helpers_group%set(key = discretization_type_key      , value = msg)                                 ; check(error==0)
    error = helpers_group%set(key = mesh_stretching_parameter_key, value = 'Stretching parameter per direction'); check(error==0)
    ! Mesh parallelization
    error = helpers_group%set(key = num_levels_key       , value = 'Number of parallel solver coarsening levels'); check(error==0)
    error = helpers_group%set(key = num_parts_per_dir_key, value = 'Number of parts per direction and per level'); check(error==0)

    ! Required
    ! --------
    ! Structured/non-structured
    error = required_group%set(key = triangulation_generate_key, value = .false.); check(error==0)
    ! Domain definition
    error = required_group%set(key = domain_length_key   , value = .false.); check(error==0)
    error = required_group%set(key = domain_origin_key   , value = .false.); check(error==0)
    error = required_group%set(key = is_dir_periodic_key , value = .false.); check(error==0)
    ! Mesh characterization
    error = required_group%set(key = reference_fe_geo_order_key   , value = .false.); check(error==0)
    error = required_group%set(key = num_dims_key     , value = .false.); check(error==0)
    error = required_group%set(key = num_cells_per_dir_key  , value = .false.); check(error==0)
    error = required_group%set(key = discretization_type_key      , value = .false.); check(error==0)
    error = required_group%set(key = mesh_stretching_parameter_key, value = .false.); check(error==0)
    ! Mesh parallelization
    error = required_group%set(key = num_levels_key       , value = .false.); check(error==0)
    error = required_group%set(key = num_parts_per_dir_key, value = .false.); check(error==0)

  end subroutine set_mesh_group_lists

  !==================================================================================================
  !> summary: Set default physical problem parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_problem_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(problem_group_key)
    switches_ab_group => switches_ab%NewSubList(problem_group_key)
    helpers_group     => helpers%NewSubList(problem_group_key)
    required_group    => required%NewSubList(problem_group_key)
    values_group      => values%NewSubList(problem_group_key)

    ! Values
    ! ------
    error = values_group%set(key = viscosity_key, value = 1.0_rp);               check(error==0)
    error = values_group%set(key = is_convection_activated_key, value = .true.); check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = viscosity_key,               value = '--viscosity');               check(error==0)
    error = switches_group%set(key = is_convection_activated_key, value = '--is-convection-activated'); check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = viscosity_key              , value = '-nu')     ; check(error==0)
    error = switches_ab_group%set(key = is_convection_activated_key, value = '-is-conv'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = viscosity_key              , value = 'Viscosity')                        ; check(error==0)
    error = helpers_group%set(key = is_convection_activated_key, value = 'Is the convective term activated?'); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = viscosity_key              , value = .false.); check(error==0)
    error = required_group%set(key = is_convection_activated_key, value = .false.); check(error==0)

  end subroutine set_problem_group_lists

  !==================================================================================================
  !> summary: Set default FE space parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_space_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(space_group_key)
    switches_ab_group => switches_ab%NewSubList(space_group_key)
    helpers_group     => helpers%NewSubList(space_group_key)
    required_group    => required%NewSubList(space_group_key)
    values_group      => values%NewSubList(space_group_key)

    ! Values
    ! ------
    error = values_group%set(key = velocity_reference_fe_order_key, value = 2); check(error==0)
    error = values_group%set(key = pressure_reference_fe_order_key, value = 1); check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = velocity_reference_fe_order_key, value = '--velocity_order'); check(error==0)
    error = switches_group%set(key = pressure_reference_fe_order_key, value = '--pressure_order'); check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = velocity_reference_fe_order_key, value = '-vorder'); check(error==0)
    error = switches_ab_group%set(key = pressure_reference_fe_order_key, value = '-porder'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = velocity_reference_fe_order_key, value = 'Velocity Finite Element space order'); check(error==0)
    error = helpers_group%set(key = pressure_reference_fe_order_key, value = 'Pressure Finite Element space order'); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = velocity_reference_fe_order_key, value = .false.); check(error==0)
    error = required_group%set(key = pressure_reference_fe_order_key, value = .false.); check(error==0)

  end subroutine set_space_group_lists

  !==================================================================================================
  !> summary: Set default FE formulation parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_formulation_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0
    character(len=:), allocatable  :: msg

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(formulation_group_key)
    switches_ab_group => switches_ab%NewSubList(formulation_group_key)
    helpers_group     => helpers%NewSubList(formulation_group_key)
    required_group    => required%NewSubList(formulation_group_key)
    values_group      => values%NewSubList(formulation_group_key)

    ! Values
    ! ------
    error = values_group%set(key = vms_method_key ,                  value = none_name);                 check(error==0)
    error = values_group%set(key = c1_constant_key,                  value = 12.0_rp);                   check(error==0)
    error = values_group%set(key = c2_constant_key,                  value = 2.0_rp);                    check(error==0)
    error = values_group%set(key = cc_constant_key,                  value = 1.0_rp);                    check(error==0)
    error = values_group%set(key = stabilization_parameter_type_key, value = codina_name);               check(error==0)
    error = values_group%set(key = elemental_length_type_key,        value = standard_elemental_length); check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = vms_method_key ,                  value = '--vms_method')             ; check(error==0)
    error = switches_group%set(key = c1_constant_key,                  value = '--c1')                     ; check(error==0)
    error = switches_group%set(key = c2_constant_key,                  value = '--c2')                     ; check(error==0)
    error = switches_group%set(key = cc_constant_key,                  value = '--cc')                     ; check(error==0)
    error = switches_group%set(key = stabilization_parameter_type_key, value = '--stabilization_parameter'); check(error==0)
    error = switches_group%set(key = elemental_length_type_key,        value = '--elemental_length_type')  ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = vms_method_key ,                  value = '-vms')    ; check(error==0)
    error = switches_ab_group%set(key = c1_constant_key,                  value = '-c1')     ; check(error==0)
    error = switches_ab_group%set(key = c2_constant_key,                  value = '-c2')     ; check(error==0)
    error = switches_ab_group%set(key = cc_constant_key,                  value = '-cc')     ; check(error==0)
    error = switches_ab_group%set(key = stabilization_parameter_type_key, value = '-stab')   ; check(error==0)
    error = switches_ab_group%set(key = elemental_length_type_key,        value = '-elength'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = vms_method_key ,                  value = 'VMS method: NONE, ASGS, OSS, COP')                 ; check(error==0)
    error = helpers_group%set(key = c1_constant_key,                  value = 'Algorithmic stabilization parameter constant 1')   ; check(error==0)
    error = helpers_group%set(key = c2_constant_key,                  value = 'Algorithmic stabilization parameter constant 2')   ; check(error==0)
    error = helpers_group%set(key = cc_constant_key,                  value = 'Algorithmic stabilization parameter constant cc')  ; check(error==0)
    error = helpers_group%set(key = stabilization_parameter_type_key, value = 'Stabilization parameter definition: CODINA, HAUKE'); check(error==0)
    msg = 'Elemental length type: (*) standard, (*) maximum, (*) minimum'
    write(msg(25:25),'(i1)') standard_elemental_length
    write(msg(39:39),'(i1)') maximum_elemental_length
    write(msg(52:52),'(i1)') minimum_elemental_length
    error = helpers_group%set(key = elemental_length_type_key, value = msg); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = vms_method_key ,                  value = .false.); check(error==0)
    error = required_group%set(key = c1_constant_key,                  value = .false.); check(error==0)
    error = required_group%set(key = c2_constant_key,                  value = .false.); check(error==0)
    error = required_group%set(key = cc_constant_key,                  value = .false.); check(error==0)
    error = required_group%set(key = stabilization_parameter_type_key, value = .false.); check(error==0)
    error = required_group%set(key = elemental_length_type_key       , value = .false.); check(error==0)

  end subroutine set_formulation_group_lists

  !==================================================================================================
  !> summary: Set default FE solution parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_solution_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(solution_group_key)
    switches_ab_group => switches_ab%NewSubList(solution_group_key)
    helpers_group     => helpers%NewSubList(solution_group_key)
    required_group    => required%NewSubList(solution_group_key)
    values_group      => values%NewSubList(solution_group_key)

    ! Values
    ! ------
    error = values_group%set(key = is_analytical_solution_key  , value = .true.);        check(error==0)
    error = values_group%set(key = is_initial_solution_key     , value = .true.);        check(error==0)
    error = values_group%set(key = is_temporal_solution_key    , value = .false.);       check(error==0)
    error = values_group%set(key = analytical_function_type_key, value = nsi_linear_name); check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = is_analytical_solution_key  , value = '--is-analytical_solution'); check(error==0)
    error = switches_group%set(key = is_initial_solution_key     , value = '--is-initial_solution')   ; check(error==0)
    error = switches_group%set(key = is_temporal_solution_key    , value = '--is-temporal_solution')  ; check(error==0)
    error = switches_group%set(key = analytical_function_type_key, value = '--analytical_function'); check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = is_analytical_solution_key  , value = '-is-asol'); check(error==0)
    error = switches_ab_group%set(key = is_initial_solution_key     , value = '-is-isol'); check(error==0)
    error = switches_ab_group%set(key = is_temporal_solution_key    , value = '-is-tsol'); check(error==0)
    error = switches_ab_group%set(key = analytical_function_type_key, value = '-afun'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = is_analytical_solution_key  , value = 'Is it an analytical solution?'); check(error==0)
    error = helpers_group%set(key = is_initial_solution_key     , value = 'Is it an initial solution?'); check(error==0)
    error = helpers_group%set(key = is_temporal_solution_key    , value = 'Is it a temporal solution?'); check(error==0)
    error = helpers_group%set(key = analytical_function_type_key, value = 'Analytical function name'); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = is_analytical_solution_key  , value = .false.); check(error==0)
    error = required_group%set(key = is_initial_solution_key     , value = .false.); check(error==0)
    error = required_group%set(key = is_temporal_solution_key    , value = .false.); check(error==0)
    error = required_group%set(key = analytical_function_type_key, value = .false.); check(error==0)

  end subroutine set_solution_group_lists

  !==================================================================================================
  !> summary: Set default time integration parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_time_integration_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(time_integration_group_key)
    switches_ab_group => switches_ab%NewSubList(time_integration_group_key)
    helpers_group     => helpers%NewSubList(time_integration_group_key)
    required_group    => required%NewSubList(time_integration_group_key)
    values_group      => values%NewSubList(time_integration_group_key)

    ! Values
    ! ------
    error = values_group%set(key = time_integrator_type_key, value = theta_method_name); check(error==0)
    error = values_group%set(key = initial_time_key        , value = 0.0_rp);         check(error==0)
    error = values_group%set(key = final_time_key          , value = 1.0_rp);         check(error==0) 
    error = values_group%set(key = time_step_key           , value = 1.0_rp);         check(error==0)   
    error = values_group%set(key = max_steps_key           , value = 1000);           check(error==0)  
    error = values_group%set(key = theta_key               , value = 1.0_rp);         check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = time_integrator_type_key, value = '--time_integrator'); check(error==0)
    error = switches_group%set(key = initial_time_key        , value = '--initial_time')   ; check(error==0)
    error = switches_group%set(key = final_time_key          , value = '--final_time')     ; check(error==0) 
    error = switches_group%set(key = time_step_key           , value = '--time_step')      ; check(error==0)   
    error = switches_group%set(key = max_steps_key           , value = '--max_steps')      ; check(error==0)  
    error = switches_group%set(key = theta_key               , value = '--theta')          ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = time_integrator_type_key, value = '-tinteg'); check(error==0)
    error = switches_ab_group%set(key = initial_time_key        , value = '-intime'); check(error==0)
    error = switches_ab_group%set(key = final_time_key          , value = '-ftime') ; check(error==0) 
    error = switches_ab_group%set(key = time_step_key           , value = '-tstep') ; check(error==0)   
    error = switches_ab_group%set(key = max_steps_key           , value = '-msteps'); check(error==0)  
    error = switches_ab_group%set(key = theta_key               , value = '-theta') ; check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = time_integrator_type_key, value = 'Time integrator name: THETA-METHOD'); check(error==0)
    error = helpers_group%set(key = initial_time_key        , value = 'Initial time'); check(error==0)
    error = helpers_group%set(key = final_time_key          , value = 'Final time') ; check(error==0) 
    error = helpers_group%set(key = time_step_key           , value = 'Time step size') ; check(error==0)   
    error = helpers_group%set(key = max_steps_key           , value = 'Maximum number of time steps'); check(error==0)  
    error = helpers_group%set(key = theta_key               , value = 'Theta value for the theta method') ; check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = time_integrator_type_key, value = .false.); check(error==0)
    error = required_group%set(key = initial_time_key        , value = .false.); check(error==0)
    error = required_group%set(key = final_time_key          , value = .false.); check(error==0) 
    error = required_group%set(key = time_step_key           , value = .false.); check(error==0)   
    error = required_group%set(key = max_steps_key           , value = .false.); check(error==0)  
    error = required_group%set(key = theta_key               , value = .false.); check(error==0)

  end subroutine set_time_integration_group_lists

  !==================================================================================================
  !> summary: Set default linear solver parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_linear_solver_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(linear_solver_group_key)
    switches_ab_group => switches_ab%NewSubList(linear_solver_group_key)
    helpers_group     => helpers%NewSubList(linear_solver_group_key)
    required_group    => required%NewSubList(linear_solver_group_key)
    values_group      => values%NewSubList(linear_solver_group_key)

    ! Values
    ! ------
    error = values_group%set(key = solver_type_key , value = direct_solver_type); check(error==0)
    error = values_group%set(key = solver_rtol_key , value = 1.0e-12_rp);         check(error==0)
    error = values_group%set(key = solver_trace_key, value = 100);                check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = solver_type_key , value = '--solver_type') ; check(error==0)
    error = switches_group%set(key = solver_rtol_key , value = '--solver_rtol') ; check(error==0)
    error = switches_group%set(key = solver_trace_key, value = '--solver_trace'); check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = solver_type_key , value = '-stype') ; check(error==0)
    error = switches_ab_group%set(key = solver_rtol_key , value = '-srtol') ; check(error==0)
    error = switches_ab_group%set(key = solver_trace_key, value = '-strace'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = solver_type_key , value = 'Solver type: direct_solver_type, CG, FGMRES, ICG, LFOM, LGMRES, MINRES, RGMRES, RICHARDSON'); check(error==0)
    error = helpers_group%set(key = solver_rtol_key , value = 'Solver tolerance')   ; check(error==0)
    error = helpers_group%set(key = solver_trace_key, value = 'Solver output trace'); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = solver_type_key , value = .false.); check(error==0)
    error = required_group%set(key = solver_rtol_key , value = .false.); check(error==0)
    error = required_group%set(key = solver_trace_key, value = .false.); check(error==0)

  end subroutine set_linear_solver_group_lists

  !==================================================================================================
  !> summary: Set default nonlinear solver parameter values, switches, switches abbreviation, helpers and
  !>          required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_nonlinear_solver_group_lists(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    type(ParameterList_t), pointer :: switches_ab
    type(ParameterList_t), pointer :: helpers
    type(ParameterList_t), pointer :: required
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: switches_group
    type(ParameterList_t), pointer :: switches_ab_group
    type(ParameterList_t), pointer :: helpers_group
    type(ParameterList_t), pointer :: required_group
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error = 0
    character(len=:), allocatable  :: msg

    ! Point to lists
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()
    values      => this%get_values()

    ! Create new sublists
    switches_group    => switches%NewSubList(nonlinear_solver_group_key)
    switches_ab_group => switches_ab%NewSubList(nonlinear_solver_group_key)
    helpers_group     => helpers%NewSubList(nonlinear_solver_group_key)
    required_group    => required%NewSubList(nonlinear_solver_group_key)
    values_group      => values%NewSubList(nonlinear_solver_group_key)

    ! Values
    ! ------
    error = values_group%set(key = nonlinear_solver_type_key          , value = picard_name) ; check(error==0)
    error = values_group%set(key = nonlinear_absolute_tolerance_key   , value = 1.0e-8_rp)   ; check(error==0)
    error = values_group%set(key = nonlinear_relative_tolerance_key   , value = 1.0e-10_rp)  ; check(error==0)
    error = values_group%set(key = nonlinear_maximum_iterations_key   , value = 50)          ; check(error==0)
    error = values_group%set(key = nonlinear_convergence_criteria_key , value = abs_res_norm); check(error==0)
    
    ! Switches
    ! --------
    error = switches_group%set(key = nonlinear_solver_type_key          , value = '--nonlinear_stype')         ; check(error==0)
    error = switches_group%set(key = nonlinear_absolute_tolerance_key   , value = '--nonlinear_atol')          ; check(error==0)
    error = switches_group%set(key = nonlinear_relative_tolerance_key   , value = '--nonlinear_rtol')          ; check(error==0)
    error = switches_group%set(key = nonlinear_maximum_iterations_key   , value = '--nonlinear_max_iterations'); check(error==0)
    error = switches_group%set(key = nonlinear_convergence_criteria_key, value = '--nonlinear_conv_criteria')  ; check(error==0)
    
    ! Switches abbreviation
    ! ---------------------
    error = switches_ab_group%set(key = nonlinear_solver_type_key         , value = '-nltype'); check(error==0)
    error = switches_ab_group%set(key = nonlinear_absolute_tolerance_key  , value = '-nlatol'); check(error==0)
    error = switches_ab_group%set(key = nonlinear_relative_tolerance_key  , value = '-nlrtol'); check(error==0)
    error = switches_ab_group%set(key = nonlinear_maximum_iterations_key  , value = '-nlmaxi'); check(error==0)
    error = switches_ab_group%set(key = nonlinear_convergence_criteria_key, value = '-nlcrit'); check(error==0)

    ! Helpers
    ! -------
    error = helpers_group%set(key = nonlinear_solver_type_key         , value = 'Nonlinear solver type: PICARD')      ; check(error==0)
    error = helpers_group%set(key = nonlinear_absolute_tolerance_key  , value = 'Nonlinear solver absolute tolerance'); check(error==0)
    error = helpers_group%set(key = nonlinear_relative_tolerance_key  , value = 'Nonlinear solver relative tolerance'); check(error==0)
    error = helpers_group%set(key = nonlinear_maximum_iterations_key  , value = 'Maximum number of iterations for the nonlinear solver'); check(error==0)
    msg = 'Nonlinear convergence criteria: (*) absolute, (*) relative w.r.t. initial, (*) relative w.r.t. RHS'
    write(msg(34:34),'(i1)') abs_res_norm
    write(msg(48:48),'(i1)') rel_r0_res_norm
    write(msg(77:77),'(i1)') rel_rhs_res_norm
    error = helpers_group%set(key = nonlinear_convergence_criteria_key, value = msg); check(error==0)

    ! Required
    ! --------
    error = required_group%set(key = nonlinear_solver_type_key         , value = .false.); check(error==0)
    error = required_group%set(key = nonlinear_absolute_tolerance_key  , value = .false.); check(error==0)
    error = required_group%set(key = nonlinear_relative_tolerance_key  , value = .false.); check(error==0)
    error = required_group%set(key = nonlinear_maximum_iterations_key  , value = .false.); check(error==0)
    error = required_group%set(key = nonlinear_convergence_criteria_key, value = .false.); check(error==0)

  end subroutine set_nonlinear_solver_group_lists

  !==================================================================================================
  !> summary: Get INPUT directory path
  !> @param [in] this NSI driver parameters
  !> @return INPUT directory path
  !==================================================================================================
  function get_dir_path(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=:), allocatable  :: get_dir_path
    type(ParameterList_t), pointer :: values
    integer(ip)                    :: error
    values  => this%get_values()
    assert(values%isAssignable(dir_path_key, 'string'))
    error = values%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  !> summary: Get triangulation type
  !> @param [in] this NSI driver parameters
  !> @return integer triangulation type
  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_triangulation_type
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values  => this%get_values()
    error = values%getsublist(mesh_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = values_group%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  !> summary: Get interpolation order for the velocity
  !> @param [in] this NSI driver parameters
  !> @return interpolation order for the velocity
  !==================================================================================================
  function get_velocity_order(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_velocity_order
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(space_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(velocity_reference_fe_order_key, get_velocity_order))
    error = values_group%get(key = velocity_reference_fe_order_key, value = get_velocity_order)
    assert(error==0)
  end function get_velocity_order

  !==================================================================================================
  !> summary: Get interpolation order for the pressure
  !> @param [in] this NSI driver parameters
  !> @return interpolation order for the pressure
  !==================================================================================================
  function get_pressure_order(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_pressure_order
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(space_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(pressure_reference_fe_order_key, get_pressure_order))
    error = values_group%get(key = pressure_reference_fe_order_key, value = get_pressure_order)
    assert(error==0)
  end function get_pressure_order

  !==================================================================================================
  !> summary: Get name of the VMS method
  !> @param [in] this NSI driver parameters
  !> @return name of the VMS method
  !==================================================================================================
  function get_vms_method_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_vms_method_type
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(formulation_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(vms_method_key, get_vms_method_type))
    error = values_group%get(key = vms_method_key, value = get_vms_method_type)
    assert(error==0)
  end function get_vms_method_type

  !==================================================================================================
  !> summary: Get name of the time integrator method
  !> @param [in] this NSI driver parameters
  !> @return  name of the time integrator method
  !==================================================================================================
  function get_time_integrator_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_time_integrator_type
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(time_integration_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(time_integrator_type_key, get_time_integrator_type))
    error = values_group%get(key = time_integrator_type_key, value = get_time_integrator_type)
    assert(error==0)
  end function get_time_integrator_type

  !==================================================================================================
  !> summary: Get the solver type
  !> @param [in] this NSI driver parameters
  !> @return  name of the solver type
  !==================================================================================================
  function get_solver_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_solver_type
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(linear_solver_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(solver_type_key, get_solver_type))
    error = values_group%get(key = solver_type_key, value = get_solver_type)
    assert(error==0)
  end function get_solver_type
  
  !==================================================================================================
  !> summary: Get solver residual tolerance
  !> @param [in] this NSI driver parameters
  !> @return solver residual tolerance
  !==================================================================================================
  function get_solver_rtol(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    real(rp)                       :: get_solver_rtol
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(linear_solver_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(solver_rtol_key, get_solver_rtol))
    error = values_group%get(key = solver_rtol_key, value = get_solver_rtol)
    assert(error==0)
  end function get_solver_rtol
  
  !==================================================================================================
  !> summary: Get solver output trace
  !> @param [in] this NSI driver parameters
  !> @return solver output trace
  !==================================================================================================
  function get_solver_trace(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_solver_trace
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(linear_solver_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(solver_trace_key, get_solver_trace))
    error = values_group%get(key = solver_trace_key, value = get_solver_trace)
    assert(error==0)
  end function get_solver_trace

  !==================================================================================================
  !> summary: Get the nonlinear solver type
  !> @param [in] this NSI driver parameters
  !> @return  name of the nonlinear solver type
  !==================================================================================================
  function get_nonlinear_solver_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_nonlinear_solver_type
    type(ParameterList_t), pointer :: values
    type(ParameterList_t), pointer :: values_group
    integer(ip)                    :: error
    values => this%get_values()
    error = values%getsublist(nonlinear_solver_group_key, values_group)
    assert(error==0)
    assert(values_group%isAssignable(nonlinear_solver_type_key, get_nonlinear_solver_type))
    error = values_group%get(key = nonlinear_solver_type_key, value = get_nonlinear_solver_type)
    assert(error==0)
  end function get_nonlinear_solver_type

end module test_nsi_parameters_names
