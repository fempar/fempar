!> @file   test_nsi/parameters.f90                                                            
!> @date   15-November-2016                                                                              
!> @author Oriol Colomes
!>
!> @brief  NSI driver parmeters module
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
  !> @class test_nsi_parameters_t
  !> @brief NSI driver parameters
  !>
  !> This type contains the parameters needed by the test_nsi. It extends the class parameter_generator_t
  !> and defines the default values depending on the CLI execution group.
  !>
  !> Currently, the following CLI group are implemented:
  !>   - analytical (Default values are defined assuming 'analytical' group)
  type, extends(parameter_generator_t), public :: test_nsi_parameters_t
     private
   contains 
     procedure                  :: create                    !< Create CLI
     ! Setters
     procedure                  :: set_default               !< Set parameter generator lists
     procedure, non_overridable :: set_default_list          !< Set default parameter values
     procedure, non_overridable :: set_switches              !< Set parameter names
     procedure, non_overridable :: set_switches_ab           !< Set parameter shortcut names
     procedure, non_overridable :: set_helpers               !< Set parameter helper messages
     procedure, non_overridable :: set_required              !< Set parameter required logicals
     ! Getters                                              
     procedure, non_overridable :: get_dir_path              !< Get INPUT directory path                     
     procedure, non_overridable :: get_triangulation_type    !< Get type of triangulation generator (from mesh or from structured)      
     procedure, non_overridable :: get_velocity_order        !< Get interpolation order for the velocity        
     procedure, non_overridable :: get_pressure_order        !< Get interpolation order for the pressure      
     procedure, non_overridable :: get_vms_method_type       !< Get VMS method type                    
     procedure, non_overridable :: get_time_integrator_type  !< Get type of the time integrator method        
     procedure, non_overridable :: get_solver_type           !< Get solver type                                           
     procedure, non_overridable :: get_solver_rtol           !< Get solver residual tolerance                
     procedure, non_overridable :: get_solver_trace          !< Get solver output trace                   
     procedure, non_overridable :: get_group                 !< Get test group name                     
     procedure, non_overridable :: get_nonlinear_solver_type !< Get nonlinear solver type                         
  end type test_nsi_parameters_t

contains

  !==================================================================================================
  !> @brief Create the test parameter list
  !>
  !> This subroutine creates the CLI inside the test_nsi_parameters_t
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine create(this)
    implicit none
    class(test_nsi_parameters_t) , intent(inout) :: this
    ! Locals
    type(Command_Line_Interface), pointer :: cli
    type(test_nsi_parameters_t)           :: test_params
    integer(ip)                           :: error

    call this%free()

    ! Initialize Command Line Interface
    cli => this%get_cli()
    call cli%init(progname    = 'test_nsi',                                      &
         &        version     = '',                                              &
         &        authors     = 'Oriol Colomes',                                 & 
         &        license     = '',                                              &
         &        description = 'Serial FEMPAR driver to solve the Navier-Stokes &
                                &incompressible problem.',                       &
         &        examples    = ['test_nsi            -h ',                      &
         &                       'test_nsi analytical -h ' ]) 

    ! Initialize parameter lists
    call this%initialize_lists()

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical solution')

    ! Set Command Line Arguments for each group
    call this%set_default()
    call this%add_to_cli_group('analytical')
    call this%parse_group('analytical')

    ! If we had more than one group, instead of the previous three calls, we should have something
    ! like what follows:
    !
    ! call this%set_default_group1()
    ! call this%add_to_cli_group('group1')
    ! call this%set_default_group2()
    ! call this%add_to_cli_group('group2')
    ! ...
    ! call this%set_default_groupN()
    ! call this%add_to_cli_group('groupN')

    ! And then parse only the group we are executing:
    !
    ! if(cli%run_command('group1') call this%parse_group('group1')
    ! if(cli%run_command('group2') call this%parse_group('group2')
    ! ...
    ! if(cli%run_command('groupN') call this%parse_group('groupN')
    

  end subroutine create

  !==================================================================================================
  !> @brief Set parameter generator for the NSI problem
  !>
  !> This subroutine sets the default values of the ParameterLists_t's needed for the parameter
  !> generator:
  !>   - list        (list with default parameter values)
  !>   - switches    (list with parameter names)
  !>   - switches_ab (list with parameter shortcut names)
  !>   - helpers     (list with help messages for each parameter)
  !>   - required    (list of logicals for the required/non-required parameters)
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_default(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this

    call this%set_default_list()
    call this%set_switches()
    call this%set_switches_ab()
    call this%set_helpers()
    call this%set_required()

  end subroutine set_default

  !==================================================================================================
  !> @brief Set default parameter values
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_default_list(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error = 0

    ! Point lists
    list => this%get_values()

    ! IO parameters
    ! -------------
    error = list%set(key = dir_path_key      , value = 'data'  ) ; check(error==0)
    error = list%set(key = prefix_key        , value = 'test'  ) ; check(error==0)
    error = list%set(key = dir_path_out_key  , value = 'output') ; check(error==0)
    error = list%set(key = write_solution_key, value = .false. ) ; check(error==0)

    ! Mesh
    ! ----
    ! Structured/non-structured
    error = list%set(key = triangulation_generate_key, value = triangulation_generate_structured); check(error==0)
    ! Domain definition
    error = list%set(key = domain_length_key   , value = [1.0_rp,1.0_rp,1.0_rp]);                   check(error==0)
    error = list%set(key = domain_origin_key   , value = [0.0_rp,0.0_rp,0.0_rp]);                   check(error==0)
    error = list%set(key = is_dir_periodic_key , value = [non_periodic,non_periodic,non_periodic]); check(error==0)
    ! Mesh characterization
    error = list%set(key = reference_fe_geo_order_key   , value = 1);                         check(error==0)
    error = list%set(key = number_of_dimensions_key     , value = 2);                         check(error==0)
    error = list%set(key = number_of_cells_per_dir_key  , value = [2,2,2]);                   check(error==0)
    error = list%set(key = discretization_type_key      , value = [uniform,uniform,uniform]); check(error==0)
    error = list%set(key = mesh_stretching_parameter_key, value = [2.75_rp,2.75_rp,2.75_rp]); check(error==0)
    ! Mesh parallelization
    error = list%set(key = number_of_levels_key       , value = 1);       check(error==0)
    error = list%set(key = number_of_parts_per_dir_key, value = [1,1,0]); check(error==0)

    ! Physical problem parameters
    ! ---------------------------
    error = list%set(key = viscosity_key, value = 1.0_rp);               check(error==0)
    error = list%set(key = is_convection_activated_key, value = .true.); check(error==0)

    ! FE space
    ! --------
    error = list%set(key = velocity_reference_fe_order_key, value = 2); check(error==0)
    error = list%set(key = pressure_reference_fe_order_key, value = 1); check(error==0)

    ! FE formulation
    ! --------------
    error = list%set(key = vms_method_key ,                  value = none_name);                 check(error==0)
    error = list%set(key = c1_constant_key,                  value = 12.0_rp);                   check(error==0)
    error = list%set(key = c2_constant_key,                  value = 2.0_rp);                    check(error==0)
    error = list%set(key = cc_constant_key,                  value = 1.0_rp);                    check(error==0)
    error = list%set(key = stabilization_parameter_type_key, value = codina_name);               check(error==0)
    error = list%set(key = elemental_length_type_key,        value = standard_elemental_length); check(error==0)

    ! FE solution
    ! -----------
    error = list%set(key = is_analytical_solution_key  , value = .true.);        check(error==0)
    error = list%set(key = is_initial_solution_key     , value = .true.);        check(error==0)
    error = list%set(key = is_temporal_solution_key    , value = .false.);       check(error==0)
    error = list%set(key = analytical_function_type_key, value = nsi_linear_name); check(error==0)

    ! Time integration
    ! ----------------
    error = list%set(key = time_integrator_type_key, value = theta_method_name); check(error==0)
    error = list%set(key = initial_time_key        , value = 0.0_rp);         check(error==0)
    error = list%set(key = final_time_key          , value = 1.0_rp);         check(error==0) 
    error = list%set(key = time_step_key           , value = 1.0_rp);         check(error==0)   
    error = list%set(key = max_steps_key           , value = 1000);           check(error==0)  
    error = list%set(key = theta_key               , value = 1.0_rp);         check(error==0)

    ! Linear solver
    ! -------------
    error = list%set(key = solver_type_key , value = rgmres_name);   check(error==0)
    error = list%set(key = solver_rtol_key , value = 1.0e-12_rp); check(error==0)
    error = list%set(key = solver_trace_key, value = 100);        check(error==0)

    ! Nonlinear solver
    ! ----------------
    error = list%set(key = nonlinear_solver_type_key          , value = picard_name) ;   check(error==0)
    error = list%set(key = nonlinear_absolute_tolerance_key   , value = 1.0e-8_rp)   ; check(error==0)
    error = list%set(key = nonlinear_relative_tolerance_key   , value = 1.0e-10_rp)  ; check(error==0)
    error = list%set(key = nonlinear_maximum_iterations_key   , value = 50)          ; check(error==0)
    error = list%set(key = nonlinear_convergence_criteria_key , value = abs_res_norm); check(error==0)

  end subroutine set_default_list

  !==================================================================================================
  !> @brief Set switches
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_switches(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches
    integer(ip)                    :: error = 0

    ! Set pointer
    switches => this%get_switches()

    ! IO parameters
    ! -------------
    error = switches%set(key = dir_path_key      , value = '--dir_path')       ; check(error==0)
    error = switches%set(key = prefix_key        , value = '--prefix')         ; check(error==0)
    error = switches%set(key = dir_path_out_key  , value = '--dir_path_out')   ; check(error==0)
    error = switches%set(key = write_solution_key, value = '--write-solution') ; check(error==0)

    ! Mesh
    ! ----
    ! Structured/non-structured
    error = switches%set(key = triangulation_generate_key, value = '--triangulation-type'); check(error==0)
    ! Domain definition
    error = switches%set(key = domain_length_key   , value = '--domain_length'); check(error==0)
    error = switches%set(key = domain_origin_key   , value = '--domain_origin'); check(error==0)
    error = switches%set(key = is_dir_periodic_key , value = '--periodicity')  ; check(error==0)
    ! Mesh characterization
    error = switches%set(key = reference_fe_geo_order_key   , value = '--reference-fe-geo-order'); check(error==0)
    error = switches%set(key = number_of_dimensions_key     , value = '--dim')                   ; check(error==0)
    error = switches%set(key = number_of_cells_per_dir_key  , value = '--number_of_cells')       ; check(error==0)
    error = switches%set(key = discretization_type_key      , value = '--discretization_type')   ; check(error==0)
    error = switches%set(key = mesh_stretching_parameter_key, value = '--stretching')            ; check(error==0)
    ! Mesh parallelization
    error = switches%set(key = number_of_levels_key       , value = '--number_of_levels'); check(error==0)
    error = switches%set(key = number_of_parts_per_dir_key, value = '--number_of_parts') ; check(error==0)

    ! Physical problem parameters
    ! ---------------------------
    error = switches%set(key = viscosity_key,               value = '--viscosity');               check(error==0)
    error = switches%set(key = is_convection_activated_key, value = '--is-convection-activated'); check(error==0)

    ! FE space
    ! --------
    error = switches%set(key = velocity_reference_fe_order_key, value = '--velocity_order'); check(error==0)
    error = switches%set(key = pressure_reference_fe_order_key, value = '--pressure_order'); check(error==0)

    ! FE formulation
    ! --------------
    error = switches%set(key = vms_method_key ,                  value = '--vms_method')             ; check(error==0)
    error = switches%set(key = c1_constant_key,                  value = '--c1')                     ; check(error==0)
    error = switches%set(key = c2_constant_key,                  value = '--c2')                     ; check(error==0)
    error = switches%set(key = cc_constant_key,                  value = '--cc')                     ; check(error==0)
    error = switches%set(key = stabilization_parameter_type_key, value = '--stabilization_parameter'); check(error==0)
    error = switches%set(key = elemental_length_type_key,        value = '--elemental_length_type')  ; check(error==0)

    ! FE solution
    ! -----------
    error = switches%set(key = is_analytical_solution_key  , value = '--is-analytical_solution'); check(error==0)
    error = switches%set(key = is_initial_solution_key     , value = '--is-initial_solution')   ; check(error==0)
    error = switches%set(key = is_temporal_solution_key    , value = '--is-temporal_solution')  ; check(error==0)
    error = switches%set(key = analytical_function_type_key, value = '--analytical_function'); check(error==0)

    ! Time integration
    ! ----------------
    error = switches%set(key = time_integrator_type_key, value = '--time_integrator'); check(error==0)
    error = switches%set(key = initial_time_key        , value = '--initial_time')   ; check(error==0)
    error = switches%set(key = final_time_key          , value = '--final_time')     ; check(error==0) 
    error = switches%set(key = time_step_key           , value = '--time_step')      ; check(error==0)   
    error = switches%set(key = max_steps_key           , value = '--max_steps')      ; check(error==0)  
    error = switches%set(key = theta_key               , value = '--theta')          ; check(error==0)

    ! Linear solver
    ! -------------
    error = switches%set(key = solver_type_key , value = '--solver_type') ; check(error==0)
    error = switches%set(key = solver_rtol_key , value = '--solver_rtol') ; check(error==0)
    error = switches%set(key = solver_trace_key, value = '--solver_trace'); check(error==0)

    ! Nonlinear solver
    ! ----------------
    error = switches%set(key = nonlinear_solver_type_key          , value = '--nonlinear_stype')         ; check(error==0)
    error = switches%set(key = nonlinear_absolute_tolerance_key   , value = '--nonlinear_atol')          ; check(error==0)
    error = switches%set(key = nonlinear_relative_tolerance_key   , value = '--nonlinear_rtol')          ; check(error==0)
    error = switches%set(key = nonlinear_maximum_iterations_key   , value = '--nonlinear_max_iterations'); check(error==0)
    error = switches%set(key = nonlinear_convergence_criteria_key, value = '--nonlinear_conv_criteria')  ; check(error==0)
   
  end subroutine set_switches

  !==================================================================================================
  !> @brief Set switches shortcuts
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_switches_ab(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: switches_ab
    integer(ip)                    :: error = 0

    ! Set pointer
    switches_ab => this%get_switches_ab()

    ! IO parameters
    ! -------------
    error = switches_ab%set(key = dir_path_key      , value = '-d')        ; check(error==0)
    error = switches_ab%set(key = prefix_key        , value = '-p')        ; check(error==0)
    error = switches_ab%set(key = dir_path_out_key  , value = '-o')        ; check(error==0)
    error = switches_ab%set(key = write_solution_key, value = '-wsolution'); check(error==0)

    ! Mesh
    ! ----
    ! Structured/non-structured
    error = switches_ab%set(key = triangulation_generate_key, value = '-tt'); check(error==0)
    ! Domain definition
    error = switches_ab%set(key = domain_length_key   , value = '-dl')    ; check(error==0)
    error = switches_ab%set(key = domain_origin_key   , value = '-do')    ; check(error==0)
    error = switches_ab%set(key = is_dir_periodic_key , value = '-period'); check(error==0)
    ! Mesh characterization
    error = switches_ab%set(key = reference_fe_geo_order_key   , value = '-gorder'); check(error==0)
    error = switches_ab%set(key = number_of_dimensions_key     , value = '-dm')    ; check(error==0)
    error = switches_ab%set(key = number_of_cells_per_dir_key  , value = '-n')     ; check(error==0)
    error = switches_ab%set(key = discretization_type_key      , value = '-dt')    ; check(error==0)
    error = switches_ab%set(key = mesh_stretching_parameter_key, value = '-st')    ; check(error==0)
    ! Mesh parallelization
    error = switches_ab%set(key = number_of_levels_key       , value = '-l') ; check(error==0)
    error = switches_ab%set(key = number_of_parts_per_dir_key, value = '-np'); check(error==0)

    ! Physical problem parameters
    ! ---------------------------
    error = switches_ab%set(key = viscosity_key              , value = '-nu')     ; check(error==0)
    error = switches_ab%set(key = is_convection_activated_key, value = '-is-conv'); check(error==0)

    ! FE space
    ! --------
    error = switches_ab%set(key = velocity_reference_fe_order_key, value = '-vorder'); check(error==0)
    error = switches_ab%set(key = pressure_reference_fe_order_key, value = '-porder'); check(error==0)

    ! FE formulation
    ! --------------
    error = switches_ab%set(key = vms_method_key ,                  value = '-vms')    ; check(error==0)
    error = switches_ab%set(key = c1_constant_key,                  value = '-c1')     ; check(error==0)
    error = switches_ab%set(key = c2_constant_key,                  value = '-c2')     ; check(error==0)
    error = switches_ab%set(key = cc_constant_key,                  value = '-cc')     ; check(error==0)
    error = switches_ab%set(key = stabilization_parameter_type_key, value = '-stab')   ; check(error==0)
    error = switches_ab%set(key = elemental_length_type_key,        value = '-elength'); check(error==0)

    ! FE solution
    ! -----------
    error = switches_ab%set(key = is_analytical_solution_key  , value = '-is-asol'); check(error==0)
    error = switches_ab%set(key = is_initial_solution_key     , value = '-is-isol'); check(error==0)
    error = switches_ab%set(key = is_temporal_solution_key    , value = '-is-tsol'); check(error==0)
    error = switches_ab%set(key = analytical_function_type_key, value = '-afun'); check(error==0)

    ! Time integration
    ! ----------------
    error = switches_ab%set(key = time_integrator_type_key, value = '-tinteg'); check(error==0)
    error = switches_ab%set(key = initial_time_key        , value = '-intime'); check(error==0)
    error = switches_ab%set(key = final_time_key          , value = '-ftime') ; check(error==0) 
    error = switches_ab%set(key = time_step_key           , value = '-tstep') ; check(error==0)   
    error = switches_ab%set(key = max_steps_key           , value = '-msteps'); check(error==0)  
    error = switches_ab%set(key = theta_key               , value = '-theta') ; check(error==0)

    ! Linear solver
    ! -------------
    error = switches_ab%set(key = solver_type_key , value = '-stype') ; check(error==0)
    error = switches_ab%set(key = solver_rtol_key , value = '-srtol') ; check(error==0)
    error = switches_ab%set(key = solver_trace_key, value = '-strace'); check(error==0)

    ! Nonlinear solver
    ! ----------------
    error = switches_ab%set(key = nonlinear_solver_type_key         , value = '-nltype'); check(error==0)
    error = switches_ab%set(key = nonlinear_absolute_tolerance_key  , value = '-nlatol'); check(error==0)
    error = switches_ab%set(key = nonlinear_relative_tolerance_key  , value = '-nlrtol'); check(error==0)
    error = switches_ab%set(key = nonlinear_maximum_iterations_key  , value = '-nlmaxi'); check(error==0)
    error = switches_ab%set(key = nonlinear_convergence_criteria_key, value = '-nlcrit'); check(error==0)
   
  end subroutine set_switches_ab

  !==================================================================================================
  !> @brief Set parameter help messages
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_helpers(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: helpers
    integer(ip)                    :: error = 0
    character(len=:), allocatable  :: msg     

    ! Set pointer
    helpers => this%get_helpers()

    ! IO parameters
    ! -------------
    error = helpers%set(key = dir_path_key      , value = 'Directory of the source files''-d')        ; check(error==0)
    error = helpers%set(key = prefix_key        , value = 'Prefix''-p')        ; check(error==0)
    error = helpers%set(key = dir_path_out_key  , value = 'Output directory''-o')             ; check(error==0)
    error = helpers%set(key = write_solution_key, value = 'Write solution''-wsolution')       ; check(error==0)

    ! Mesh
    ! ----
    ! Structured/non-structured
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triangulation_generate_key, value = msg); check(error==0)
    ! Domain definition
    error = helpers%set(key = domain_length_key  , value = 'Domain length per direction')     ; check(error==0)
    error = helpers%set(key = domain_origin_key  , value = 'Domain origin coordinates')       ; check(error==0)
    error = helpers%set(key = is_dir_periodic_key, value = 'Domain periodicity per direction'); check(error==0)
    ! Mesh characterization
    error = helpers%set(key = reference_fe_geo_order_key , value = 'Order of the triangulation reference Finite Element'); check(error==0)
    error = helpers%set(key = number_of_dimensions_key   , value = 'Number of space dimensions')                         ; check(error==0)
    error = helpers%set(key = number_of_cells_per_dir_key, value = 'Number of cells per direction')                      ; check(error==0)
    msg = 'Discretization type: (*) uniform, (*) cubic, (*) tanh'
    write(msg(23:23),'(i1)') uniform
    write(msg(36:36),'(i1)') cubic
    write(msg(47:47),'(i1)') tanh
    error = helpers%set(key = discretization_type_key      , value = msg)                                 ; check(error==0)
    error = helpers%set(key = mesh_stretching_parameter_key, value = 'Stretching parameter per direction'); check(error==0)
    ! Mesh parallelization
    error = helpers%set(key = number_of_levels_key       , value = 'Number of parallel solver coarsening levels'); check(error==0)
    error = helpers%set(key = number_of_parts_per_dir_key, value = 'Number of parts per direction and per level'); check(error==0)

    ! Physical problem parameters
    ! ---------------------------
    error = helpers%set(key = viscosity_key              , value = 'Viscosity')                        ; check(error==0)
    error = helpers%set(key = is_convection_activated_key, value = 'Is the convective term activated?'); check(error==0)

    ! FE space
    ! --------
    error = helpers%set(key = velocity_reference_fe_order_key, value = 'Velocity Finite Element space order'); check(error==0)
    error = helpers%set(key = pressure_reference_fe_order_key, value = 'Pressure Finite Element space order'); check(error==0)

    ! FE formulation
    ! --------------
    error = helpers%set(key = vms_method_key ,                  value = 'VMS method: NONE, ASGS, OSS, COP')                 ; check(error==0)
    error = helpers%set(key = c1_constant_key,                  value = 'Algorithmic stabilization parameter constant 1')   ; check(error==0)
    error = helpers%set(key = c2_constant_key,                  value = 'Algorithmic stabilization parameter constant 2')   ; check(error==0)
    error = helpers%set(key = cc_constant_key,                  value = 'Algorithmic stabilization parameter constant cc')  ; check(error==0)
    error = helpers%set(key = stabilization_parameter_type_key, value = 'Stabilization parameter definition: CODINA, HAUKE'); check(error==0)
    msg = 'Elemental length type: (*) standard, (*) maximum, (*) minimum'
    write(msg(25:25),'(i1)') standard_elemental_length
    write(msg(39:39),'(i1)') maximum_elemental_length
    write(msg(52:52),'(i1)') minimum_elemental_length
    error = helpers%set(key = elemental_length_type_key, value = msg); check(error==0)

    ! FE solution
    ! -----------
    error = helpers%set(key = is_analytical_solution_key  , value = 'Is it an analytical solution?'); check(error==0)
    error = helpers%set(key = is_initial_solution_key     , value = 'Is it an initial solution?'); check(error==0)
    error = helpers%set(key = is_temporal_solution_key    , value = 'Is it a temporal solution?'); check(error==0)
    error = helpers%set(key = analytical_function_type_key, value = 'Analytical function name'); check(error==0)

    ! Time integration
    ! ----------------
    error = helpers%set(key = time_integrator_type_key, value = 'Time integrator name: THETA-METHOD'); check(error==0)
    error = helpers%set(key = initial_time_key        , value = 'Initial time'); check(error==0)
    error = helpers%set(key = final_time_key          , value = 'Final time') ; check(error==0) 
    error = helpers%set(key = time_step_key           , value = 'Time step size') ; check(error==0)   
    error = helpers%set(key = max_steps_key           , value = 'Maximum number of time steps'); check(error==0)  
    error = helpers%set(key = theta_key               , value = 'Theta value for the theta method') ; check(error==0)

    ! Solver
    ! ------
    error = helpers%set(key = solver_type_key , value = 'Solver type: direct_solver_type, CG, FGMRES, ICG, LFOM, LGMRES, MINRES, RGMRES, RICHARDSON'); check(error==0)
    error = helpers%set(key = solver_rtol_key , value = 'Solver tolerance')   ; check(error==0)
    error = helpers%set(key = solver_trace_key, value = 'Solver output trace'); check(error==0)

    ! Nonlinear solver
    ! ----------------
    error = helpers%set(key = nonlinear_solver_type_key         , value = 'Nonlinear solver type: PICARD')      ; check(error==0)
    error = helpers%set(key = nonlinear_absolute_tolerance_key  , value = 'Nonlinear solver absolute tolerance'); check(error==0)
    error = helpers%set(key = nonlinear_relative_tolerance_key  , value = 'Nonlinear solver relative tolerance'); check(error==0)
    error = helpers%set(key = nonlinear_maximum_iterations_key  , value = 'Maximum number of iterations for the nonlinear solver'); check(error==0)
    msg = 'Nonlinear convergence criteria: (*) absolute, (*) relative w.r.t. initial, (*) relative w.r.t. RHS'
    write(msg(34:34),'(i1)') abs_res_norm
    write(msg(48:48),'(i1)') rel_r0_res_norm
    write(msg(77:77),'(i1)') rel_rhs_res_norm
    error = helpers%set(key = nonlinear_convergence_criteria_key, value = msg); check(error==0)
   
  end subroutine set_helpers
  
  !==================================================================================================
  !> @brief Set required flags
  !> @param [in,out] this NSI driver parameters
  !==================================================================================================
  subroutine set_required(this)
    implicit none
    class(test_nsi_parameters_t), intent(inout) :: this
    ! Locals
    type(ParameterList_t), pointer :: required
    integer(ip)                    :: error = 0

    ! Set pointer
    required => this%get_required()

    ! IO parameters
    ! -------------
    error = required%set(key = dir_path_key      , value = .false.); check(error==0)
    error = required%set(key = prefix_key        , value = .false.); check(error==0)
    error = required%set(key = dir_path_out_key  , value = .false.); check(error==0)
    error = required%set(key = write_solution_key, value = .false.); check(error==0)

    ! Mesh
    ! ----
    ! Structured/non-structured
    error = required%set(key = triangulation_generate_key, value = .false.); check(error==0)
    ! Domain definition
    error = required%set(key = domain_length_key   , value = .false.); check(error==0)
    error = required%set(key = domain_origin_key   , value = .false.); check(error==0)
    error = required%set(key = is_dir_periodic_key , value = .false.); check(error==0)
    ! Mesh characterization
    error = required%set(key = reference_fe_geo_order_key   , value = .false.); check(error==0)
    error = required%set(key = number_of_dimensions_key     , value = .false.); check(error==0)
    error = required%set(key = number_of_cells_per_dir_key  , value = .false.); check(error==0)
    error = required%set(key = discretization_type_key      , value = .false.); check(error==0)
    error = required%set(key = mesh_stretching_parameter_key, value = .false.); check(error==0)
    ! Mesh parallelization
    error = required%set(key = number_of_levels_key       , value = .false.); check(error==0)
    error = required%set(key = number_of_parts_per_dir_key, value = .false.); check(error==0)

    ! Physical problem parameters
    ! ---------------------------
    error = required%set(key = viscosity_key              , value = .false.); check(error==0)
    error = required%set(key = is_convection_activated_key, value = .false.); check(error==0)

    ! FE space
    ! --------
    error = required%set(key = velocity_reference_fe_order_key, value = .false.); check(error==0)
    error = required%set(key = pressure_reference_fe_order_key, value = .false.); check(error==0)

    ! FE formulation
    ! --------------
    error = required%set(key = vms_method_key ,                  value = .false.); check(error==0)
    error = required%set(key = c1_constant_key,                  value = .false.); check(error==0)
    error = required%set(key = c2_constant_key,                  value = .false.); check(error==0)
    error = required%set(key = cc_constant_key,                  value = .false.); check(error==0)
    error = required%set(key = stabilization_parameter_type_key, value = .false.); check(error==0)
    error = required%set(key = elemental_length_type_key       , value = .false.); check(error==0)

    ! FE solution
    ! -----------
    error = required%set(key = is_analytical_solution_key  , value = .false.); check(error==0)
    error = required%set(key = is_initial_solution_key     , value = .false.); check(error==0)
    error = required%set(key = is_temporal_solution_key    , value = .false.); check(error==0)
    error = required%set(key = analytical_function_type_key, value = .false.); check(error==0)

    ! Time integration
    ! ----------------
    error = required%set(key = time_integrator_type_key, value = .false.); check(error==0)
    error = required%set(key = initial_time_key        , value = .false.); check(error==0)
    error = required%set(key = final_time_key          , value = .false.); check(error==0) 
    error = required%set(key = time_step_key           , value = .false.); check(error==0)   
    error = required%set(key = max_steps_key           , value = .false.); check(error==0)  
    error = required%set(key = theta_key               , value = .false.); check(error==0)

    ! Linear solver
    ! -------------
    error = required%set(key = solver_type_key , value = .false.); check(error==0)
    error = required%set(key = solver_rtol_key , value = .false.); check(error==0)
    error = required%set(key = solver_trace_key, value = .false.); check(error==0)

    ! Nonlinear solver
    ! ----------------
    error = required%set(key = nonlinear_solver_type_key         , value = .false.); check(error==0)
    error = required%set(key = nonlinear_absolute_tolerance_key  , value = .false.); check(error==0)
    error = required%set(key = nonlinear_relative_tolerance_key  , value = .false.); check(error==0)
    error = required%set(key = nonlinear_maximum_iterations_key  , value = .false.); check(error==0)
    error = required%set(key = nonlinear_convergence_criteria_key, value = .false.); check(error==0)
   
  end subroutine set_required

  !==================================================================================================
  !> @brief Get INPUT directory path
  !> @param [in] this NSI driver parameters
  !> @return INPUT directory path
  !==================================================================================================
  function get_dir_path(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=:), allocatable  :: get_dir_path
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list  => this%get_values()
    assert(list%isAssignable(dir_path_key, 'string'))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  !> @brief Get triangulation type
  !> @param [in] this NSI driver parameters
  !> @return integer triangulation type
  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_triangulation_type
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list  => this%get_values()
    assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key = triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type 

  !==================================================================================================
  !> @brief Get interpolation order for the velocity
  !> @param [in] this NSI driver parameters
  !> @return interpolation order for the velocity
  !==================================================================================================
  function get_velocity_order(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_velocity_order
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(velocity_reference_fe_order_key, get_velocity_order))
    error = list%get(key = velocity_reference_fe_order_key, value = get_velocity_order)
    assert(error==0)
  end function get_velocity_order

  !==================================================================================================
  !> @brief Get interpolation order for the pressure
  !> @param [in] this NSI driver parameters
  !> @return interpolation order for the pressure
  !==================================================================================================
  function get_pressure_order(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_pressure_order
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(pressure_reference_fe_order_key, get_pressure_order))
    error = list%get(key = pressure_reference_fe_order_key, value = get_pressure_order)
    assert(error==0)
  end function get_pressure_order

  !==================================================================================================
  !> @brief Get name of the VMS method
  !> @param [in] this NSI driver parameters
  !> @return name of the VMS method
  !==================================================================================================
  function get_vms_method_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_vms_method_type
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(vms_method_key, get_vms_method_type))
    error = list%get(key = vms_method_key, value = get_vms_method_type)
    assert(error==0)
  end function get_vms_method_type

  !==================================================================================================
  !> @brief Get name of the time integrator method
  !> @param [in] this NSI driver parameters
  !> @return  name of the time integrator method
  !==================================================================================================
  function get_time_integrator_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_time_integrator_type
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(time_integrator_type_key, get_time_integrator_type))
    error = list%get(key = time_integrator_type_key, value = get_time_integrator_type)
    assert(error==0)
  end function get_time_integrator_type

  !==================================================================================================
  !> @brief Get the solver type
  !> @param [in] this NSI driver parameters
  !> @return  name of the solver type
  !==================================================================================================
  function get_solver_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_solver_type
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(solver_type_key, get_solver_type))
    error = list%get(key = solver_type_key, value = get_solver_type)
    assert(error==0)
  end function get_solver_type
  
  !==================================================================================================
  !> @brief Get solver residual tolerance
  !> @param [in] this NSI driver parameters
  !> @return solver residual tolerance
  !==================================================================================================
  function get_solver_rtol(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    real(rp)                       :: get_solver_rtol
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(solver_rtol_key, get_solver_rtol))
    error = list%get(key = solver_rtol_key, value = get_solver_rtol)
    assert(error==0)
  end function get_solver_rtol
  
  !==================================================================================================
  !> @brief Get solver output trace
  !> @param [in] this NSI driver parameters
  !> @return solver output trace
  !==================================================================================================
  function get_solver_trace(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    integer(ip)                    :: get_solver_trace
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(solver_trace_key, get_solver_trace))
    error = list%get(key = solver_trace_key, value = get_solver_trace)
    assert(error==0)
  end function get_solver_trace

  !==================================================================================================
  !> @brief Get the nonlinear solver type
  !> @param [in] this NSI driver parameters
  !> @return  name of the nonlinear solver type
  !==================================================================================================
  function get_nonlinear_solver_type(this)
    implicit none
    class(test_nsi_parameters_t), intent(in) :: this
    character(len=2048)            :: get_nonlinear_solver_type
    type(ParameterList_t), pointer :: list
    integer(ip)                    :: error
    list => this%get_values()
    assert(list%isAssignable(nonlinear_solver_type_key, get_nonlinear_solver_type))
    error = list%get(key = nonlinear_solver_type_key, value = get_nonlinear_solver_type)
    assert(error==0)
  end function get_nonlinear_solver_type
  
  !==================================================================================================
  !> @brief Get test group name
  !> @param [in] this NSI driver parameters
  !> @return  test group name
  !==================================================================================================
  function get_group(this)
    implicit none
    class(test_nsi_parameters_t) , intent(in) :: this
    character(len=:), allocatable :: get_group
!!$    if(this%cli%run_command('analytical')) then
!!$       get_group = 'analytical'
!!$    else
!!$       write(*,*) 'There is any matching for the CLI group in ',__FILE__,__LINE__
!!$       check(.false.)
!!$    end if
  end function get_group

end module test_nsi_parameters_names
