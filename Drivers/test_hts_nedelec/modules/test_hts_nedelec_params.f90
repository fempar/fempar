module hts_nedelec_params_names
  use fempar_names
  implicit none
#include "debug.i90" 
  private

  ! IO parameters 
    character(len=*), parameter :: write_solution_key       = 'write_solution'   

    ! MESHING parameters 
    character(len=*), parameter :: domain_length_key             = 'domain_length'
    character(len=*), parameter :: hts_domain_length_key         = 'hts_domain_length'

    ! FE ELEMENT parameters 
    character(len=*), parameter :: magnetic_field_reference_fe_order_key      = 'magnetic_field_reference_fe_order'
    character(len=*), parameter :: magnetic_pressure_reference_fe_order_key   = 'magnetic_pressure_reference_fe_order'

    ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    character(len=*), parameter :: external_magnetic_field_amplitude_key      = 'external_magnetic_field_amplitude' 
    character(len=*), parameter :: external_magnetic_field_frequency_key      = 'external_magnetic_field_frequency' 
    character(len=*), parameter :: external_current_amplitude_key             = 'external_current_amplitude' 
    character(len=*), parameter :: external_current_frequency_key             = 'external_current_frequency' 
    character(len=*), parameter :: apply_current_density_constraint_key       = 'apply_current_density_constraint' 

    ! PHYSICAL PROPERTIES parameters
    character(len=*), parameter :: air_permeability_key           = 'air_permeability'
    character(len=*), parameter :: air_resistivity_key            = 'air_resistivity' 
    character(len=*), parameter :: hts_permeability_key           = 'hts_permeability' 
    character(len=*), parameter :: hts_resistivity_key            = 'hts_resistivity' 
    character(len=*), parameter :: critical_current_key           = 'critical_current' 
    character(len=*), parameter :: critical_electric_field_key    = 'critical_electric_field' 
    character(len=*), parameter :: nonlinear_exponent_key         = 'nonlinear_exponent' 

    ! TIME INTEGRATION parameters 
    character(len=*), parameter :: theta_value_key                = 'theta_value' 
    character(len=*), parameter :: initial_time_key               = 'initial_time' 
    character(len=*), parameter :: final_time_key                 = 'final_time' 
    character(len=*), parameter :: number_time_steps_key          = 'number_time_steps' 
    character(len=*), parameter :: is_adaptive_time_stepping_key  = 'is_adaptive_time_stepping' 
    character(len=*), parameter :: stepping_parameter_key         = 'stepping_parameter' 
    character(len=*), parameter :: max_time_step_key              = 'max_time_step' 
    character(len=*), parameter :: min_time_step_key              = 'min_time_step' 
    character(len=*), parameter :: save_solution_n_steps_key      = 'save_solution_n_steps' 

    ! NONLINEAR SOLVER parameters 
    character(len=*), parameter :: nonlinear_convergence_criteria_key  = 'nonlinear_convergence_criteria' 
    character(len=*), parameter :: absolute_nonlinear_tolerance_key    = 'absolute_nonlinear_tolerance' 
    character(len=*), parameter :: relative_nonlinear_tolerance_key    = 'relative_nonlinear_tolerance' 
    character(len=*), parameter :: max_nonlinear_iterations_key        = 'max_nonlinear_iterations' 
       
  type, extends(parameter_generator_t) :: hts_nedelec_params_t
     private
     contains
       procedure                              :: set_default  => hts_nedelec_params_set_default
       procedure, non_overridable             :: get_dir_path
       procedure, non_overridable             :: get_prefix
       procedure, non_overridable             :: get_dir_path_out
       procedure, non_overridable             :: get_write_solution
       procedure, non_overridable             :: get_triangulation_type
       procedure, non_overridable             :: get_domain_length 
       procedure, non_overridable             :: get_hts_domain_length
       procedure, non_overridable             :: get_magnetic_field_reference_fe_order
       procedure, non_overridable             :: get_magnetic_pressure_reference_fe_order
       procedure, non_overridable             :: get_external_magnetic_field_amplitude
       procedure, non_overridable             :: get_external_magnetic_field_frequency
       procedure, non_overridable             :: get_external_current_amplitude
       procedure, non_overridable             :: get_external_current_frequency
       procedure, non_overridable             :: get_apply_current_density_constraint
       procedure, non_overridable             :: get_air_permeability
       procedure, non_overridable             :: get_air_resistivity
       procedure, non_overridable             :: get_hts_permeability
       procedure, non_overridable             :: get_hts_resistivity
       procedure, non_overridable 	           :: get_critical_current           
       procedure, non_overridable             :: get_critical_electric_field   
       procedure, non_overridable             :: get_nonlinear_exponent
       procedure, non_overridable             :: get_theta_value 
       procedure, non_overridable             :: get_initial_time 
       procedure, non_overridable             :: get_final_time 
       procedure, non_overridable             :: get_number_time_steps
       procedure, non_overridable             :: get_is_adaptive_time_stepping 
       procedure, non_overridable             :: get_stepping_parameter
       procedure, non_overridable             :: get_max_time_step 
       procedure, non_overridable             :: get_min_time_step 
       procedure, non_overridable             :: get_save_solution_n_steps 
       procedure, non_overridable             :: get_nonlinear_convergence_criteria 
       procedure, non_overridable             :: get_absolute_nonlinear_tolerance
       procedure, non_overridable             :: get_relative_nonlinear_tolerance 
       procedure, non_overridable             :: get_max_nonlinear_iterations 
    end type hts_nedelec_params_t

  ! Types
  public :: hts_nedelec_params_t 

contains

  !==================================================================================================
  subroutine hts_nedelec_params_set_default(this)
    implicit none
    class(hts_nedelec_params_t), intent(inout) :: this
    type(ParameterList_t), pointer :: list, switches, switches_ab, helpers, required
    integer(ip)    :: error
    character(len=:), allocatable            :: msg

    list        => this%get_parameters()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()
    required    => this%get_required()

    ! SETTING DEFAULT PARAMETERS **************************************************************

    ! IO parameters 
    error = list%set(key = dir_path_key       , value = '.')      ; check(error==0)
    error = list%set(key = prefix_key         , value = 'square') ; check(error==0)
    error = list%set(key = dir_path_out_key   , value = '.')      ; check(error==0)
    error = list%set(key = write_solution_key , value =  .false.) ; check(error==0)

    ! MESHING parameters 
    error = list%set(key = number_of_dimensions_key   , value =  2)                                 ; check(error==0)
    error = list%set(key = number_of_cells_per_dir_key, value =  [2,2,2])                           ; check(error==0)
    error = list%set(key = is_dir_periodic_key        , value =  [0,0,0])                           ; check(error==0)
    error = list%set(key = triangulation_generate_key , value =  triangulation_generate_structured) ; check(error==0)
    error = list%set(key = domain_length_key          , value =  [1.0,1.0,1.0])                     ; check(error==0)
    error = list%set(key = hts_domain_length_key      , value =  [0.5,0.5,0.5])                     ; check(error==0)

    ! FE ELEMENT parameters 
    error = list%set(key = magnetic_field_reference_fe_order_key    , value = 1)   ; check(error==0)
    error = list%set(key = magnetic_pressure_reference_fe_order_key , value = 1)   ; check(error==0)

    ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    error = list%set(key = external_magnetic_field_amplitude_key , value = [0.0,0.0,0.0]) ; check(error==0) 
    error = list%set(key = external_magnetic_field_frequency_key , value = 50.0)      ; check(error==0)
    error = list%set(key = external_current_amplitude_key , value = [0.0,0.0,0.0])        ; check(error==0) 
    error = list%set(key = external_current_frequency_key , value = 50.0)             ; check(error==0)
    error = list%set(key = apply_current_density_constraint_key , value = .false.)  ; check(error==0)

    ! PHYSICAL PROPERTIES parameters
    error = list%set(key = air_permeability_key        , value = 1.0)      ; check(error==0)
    error = list%set(key = air_resistivity_key         , value = 1.0)      ; check(error==0)
    error = list%set(key = hts_permeability_key        , value = 1.0)      ; check(error==0)
    error = list%set(key = hts_resistivity_key         , value = 1.0)      ; check(error==0)
    error = list%set(key = critical_current_key        , value = 1.0)      ; check(error==0)
    error = list%set(key = critical_electric_field_key , value = 1.0)      ; check(error==0)
    error = list%set(key = nonlinear_exponent_key      , value = 1.0)      ; check(error==0)

    ! TIME INTEGRATION parameters 
    error = list%set(key = theta_value_key              , value = 1.0)      ; check(error==0)
    error = list%set(key = initial_time_key             , value = 0.0)      ; check(error==0)
    error = list%set(key = final_time_key               , value = 1.0)      ; check(error==0)
    error = list%set(key = number_time_steps_key        , value = 10)       ; check(error==0)
    error = list%set(key = is_adaptive_time_stepping_key, value = .true.)   ; check(error==0)
    error = list%set(key = stepping_parameter_key       , value = 10)       ; check(error==0)
    error = list%set(key = max_time_step_key            , value = 0.1)      ; check(error==0)
    error = list%set(key = min_time_step_key            , value = 0.1)      ; check(error==0)
    error = list%set(key = save_solution_n_steps_key    , value = 10)       ; check(error==0)

    ! NONLINEAR SOLVER parameters 
    error = list%set(key = nonlinear_convergence_criteria_key , value = 'rel_rhs_res_norm')   ; check(error==0)
    error = list%set(key = absolute_nonlinear_tolerance_key   , value = 1.0e-2)               ; check(error==0)
    error = list%set(key = relative_nonlinear_tolerance_key   , value = 1.0e-12)              ; check(error==0)
    error = list%set(key = max_nonlinear_iterations_key       , value = 50)                   ; check(error==0)
   
 
    ! CLI commands **************************************************************************************
    
    ! IO parameters 
    error = switches%set(key = dir_path_key                  , value = '--dir-path')         ; check(error==0)
    error = switches%set(key = prefix_key                    , value = '--prefix')           ; check(error==0)
    error = switches%set(key = dir_path_out_key              , value = '--dir-path-out')     ; check(error==0)
    error = switches%set(key = write_solution_key            , value = '--write_solution')   ; check(error==0)
    
    ! MESHING parameters 
    error = switches%set(key = number_of_dimensions_key      , value = '--number_dimensions') ; check(error==0)
    error = switches%set(key = number_of_cells_per_dir_key   , value = '--number_of_cells')   ; check(error==0)
    error = switches%set(key = is_dir_periodic_key           , value = '--is_dir_periodic')   ; check(error==0)
    error = switches%set(key = triangulation_generate_key    , value = '--triangulation-type'); check(error==0)
    error = switches%set(key = domain_length_key             , value = '--domain_length')     ; check(error==0)
    error = switches%set(key = hts_domain_length_key         , value = '--hts_domain_length') ; check(error==0)
    
    ! FE ELEMENT parameters 
    error = switches%set(key = magnetic_field_reference_fe_order_key     , value = '--magnetic_field_reference_fe_order')   ; check(error==0)
    error = switches%set(key = magnetic_pressure_reference_fe_order_key  , value = '--magnetic_pressure_reference_fe_order'); check(error==0)  
    
    ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    error = switches%set(key = external_magnetic_field_amplitude_key , value = '--external_magnetic_field_amplitude') ; check(error==0) 
    error = switches%set(key = external_magnetic_field_frequency_key , value = '--external_magnetic_field_frequency' ) ; check(error==0)
    error = switches%set(key = external_current_amplitude_key , value = '--external_current_amplitude') ; check(error==0) 
    error = switches%set(key = external_current_frequency_key , value = '--external_current_frequency'); check(error==0)
    error = switches%set(key = apply_current_density_constraint_key , value = '--apply_current_density_constraint')  ; check(error==0)
    
    ! PHYSICAL PROPERTIES parameters
    error = switches%set(key = air_permeability_key        , value = '--air_permeability')        ; check(error==0)
    error = switches%set(key = air_resistivity_key         , value = '--air_resistivity')         ; check(error==0)
    error = switches%set(key = hts_permeability_key        , value = '--hts_permeability')        ; check(error==0)
    error = switches%set(key = hts_resistivity_key         , value = '--hts_resistivity')         ; check(error==0)
    error = switches%set(key = critical_current_key        , value = '--critical_current')        ; check(error==0)
    error = switches%set(key = critical_electric_field_key , value = '--critical_electric_field') ; check(error==0)
    error = switches%set(key = nonlinear_exponent_key      , value = '--nonlinear_exponent')      ; check(error==0)
    
    ! TIME INTEGRATION parameters
    error = switches%set(key = theta_value_key              , value = '--theta_value')              ; check(error==0)
    error = switches%set(key = initial_time_key             , value = '--initial_time')             ; check(error==0)
    error = switches%set(key = final_time_key               , value = '--final_time')               ; check(error==0)
    error = switches%set(key = number_time_steps_key        , value = '--number_time_steps')        ; check(error==0)
    error = switches%set(key = is_adaptive_time_stepping_key, value ='--is_adaptive_time_stepping') ; check(error==0)
    error = switches%set(key = stepping_parameter_key       , value = '--stepping_parameter')       ; check(error==0)
    error = switches%set(key = max_time_step_key            , value = '--max_time_step')            ; check(error==0)
    error = switches%set(key = min_time_step_key            , value = '--min_time_step')            ; check(error==0)
    error = switches%set(key = save_solution_n_steps_key    , value = '--save_solution_interval')   ; check(error==0)
    
    ! NONLINEAR SOLVER parameters
    error = switches%set(key = nonlinear_convergence_criteria_key , value = '--convergence_criteria')    ;check(error==0)
    error = switches%set(key = absolute_nonlinear_tolerance_key   , value = '--absolute_residual')       ;check(error==0)
    error = switches%set(key = relative_nonlinear_tolerance_key   , value = '--relative_residual')       ; check(error==0)
    error = switches%set(key = max_nonlinear_iterations_key       , value = '--max_nonlinear_iterations');check(error==0)
    
    
    
    ! SHORT abreviations ***********************************************************************
    
    ! IO parameters
    error = switches_ab%set(key = dir_path_key       , value = '-dir-path')     ; check(error==0)
    error = switches_ab%set(key = prefix_key         , value = '-prf')          ; check(error==0)
    error = switches_ab%set(key = dir_path_out_key   , value = '-out')          ; check(error==0)
    
    ! MESHING parameters 
    error = switches_ab%set(key = write_solution_key         , value = '-wsolution') ; check(error==0)
    error = switches_ab%set(key = number_of_dimensions_key   , value = '-dim')       ; check(error==0)
    error = switches_ab%set(key = number_of_cells_per_dir_key, value = '-nelem')     ; check(error==0)
    error = switches_ab%set(key = is_dir_periodic_key        , value = '-idp')       ; check(error==0)
    error = switches_ab%set(key = triangulation_generate_key , value = '-tt')        ; check(error==0)
    error = switches_ab%set(key = domain_length_key          , value = '-dl')        ; check(error==0)
    error = switches_ab%set(key = hts_domain_length_key      , value = '-hts_dl')    ; check(error==0)
    
    ! FE ELEMENT parameters 
    error = switches_ab%set(key = magnetic_field_reference_fe_order_key    , value = '-mf_order')  ; check(error==0)
    error = switches_ab%set(key = magnetic_pressure_reference_fe_order_key , value = '-mp_order')  ; check(error==0)   
    
    ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    error = switches_ab%set(key = external_magnetic_field_amplitude_key , value = '-H')    ; check(error==0) 
    error = switches_ab%set(key = external_magnetic_field_frequency_key , value = '-w_H' ) ; check(error==0)
    error = switches_ab%set(key = external_current_amplitude_key        , value = '-J')    ; check(error==0) 
    error = switches_ab%set(key = external_current_frequency_key        , value = '-w_J')  ; check(error==0)
    error = switches_ab%set(key = apply_current_density_constraint_key  , value = '-cdc')  ; check(error==0)
    
    ! PHYSICAL PROPERTIES parameters
    error = switches_ab%set(key = air_permeability_key        , value = '-mu_air')     ; check(error==0)
    error = switches_ab%set(key = air_resistivity_key         , value = '-rho_air')    ; check(error==0)
    error = switches_ab%set(key = hts_permeability_key        , value = '-mu_hts')     ; check(error==0)
    error = switches_ab%set(key = hts_resistivity_key         , value = '-rho_hts')    ; check(error==0)
    error = switches_ab%set(key = critical_current_key        , value = '-Jc    ')     ; check(error==0)
    error = switches_ab%set(key = critical_electric_field_key , value = '-Ec')         ; check(error==0)
    error = switches_ab%set(key = nonlinear_exponent_key      , value = '-nl_exp')     ; check(error==0)
    
    ! TIME INTEGRATION parameters
    error = switches_ab%set(key = theta_value_key              , value = '-theta')  ; check(error==0)
    error = switches_ab%set(key = initial_time_key             , value = '-t0')     ; check(error==0)
    error = switches_ab%set(key = final_time_key               , value = '-tf')     ; check(error==0)
    error = switches_ab%set(key = number_time_steps_key        , value = '-nsteps') ; check(error==0)
    error = switches_ab%set(key = is_adaptive_time_stepping_key, value = '-iats')   ; check(error==0)
    error = switches_ab%set(key = stepping_parameter_key       , value = '-sp')     ; check(error==0)
    error = switches_ab%set(key = max_time_step_key            , value = '-max_ts') ; check(error==0)
    error = switches_ab%set(key = min_time_step_key            , value = '-min_ts') ; check(error==0)
    error = switches_ab%set(key = save_solution_n_steps_key    , value = '-ssi')    ; check(error==0)
    
    ! NONLINEAR SOLVER parameters
    error = switches_ab%set(key = nonlinear_convergence_criteria_key , value = '-conv_crit')  ; check(error==0)
    error = switches_ab%set(key = absolute_nonlinear_tolerance_key   , value = '-abs_tol')    ; check(error==0)
    error = switches_ab%set(key = relative_nonlinear_tolerance_key   , value = '-rel_tol')    ; check(error==0)
    error = switches_ab%set(key = max_nonlinear_iterations_key       , value = '-max_nl_its') ; check(error==0)


    ! HELPERS in some parameters ************************************************************************ 
    
    ! IO parameters
    error = helpers%set(key = dir_path_key       , value = 'Directory where input files are found')     ; check(error==0)
    error = helpers%set(key = prefix_key         , value = 'Name of the input GID files')     ; check(error==0)
    error = helpers%set(key = dir_path_out_key   , value = 'Output directory')                ; check(error==0)
    error = helpers%set(key = write_solution_key , value = 'Print solution?')    ; check(error==0)
    
    ! MESHING parameters
    error = helpers%set(key = number_of_dimensions_key   , value = 'Number of dimensions')  ; check(error==0)
    error = helpers%set(key = number_of_cells_per_dir_key, value = 'Number of elements in each direction'); check(error==0)
    error = helpers%set(key = is_dir_periodic_key        , value = 'Is periodic? in each direction'); check(error==0)
  
    msg = 'structured (*) or unstructured (*) triangulation?'
    write(msg(13:13),'(i1)') triangulation_generate_structured
    write(msg(33:33),'(i1)') triangulation_generate_from_mesh
    error = helpers%set(key = triangulation_generate_key     , value = msg)  ; check(error==0)
    error = helpers%set(key = domain_length_key              , value = 'Domain length for each direction')   ; check(error==0)
    error = helpers%set(key = hts_domain_length_key          , value = 'High Temperature Superconductor Device length ( concentric with the domain) ') ; check(error==0)
    
    ! FE ELEMENT parameters
    error = helpers%set(key = magnetic_field_reference_fe_order_key  , value = 'Magnetic field reference Finite Element Order '); check(error==0)
    error = helpers%set(key = magnetic_pressure_reference_fe_order_key , value = 'Magnetic Pressure reference Finite Element Order'); check(error==0)    
    
    ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    error = helpers%set(key = external_magnetic_field_amplitude_key , value = 'External magnetic field amplitude per direction') ; check(error==0) 
    error = helpers%set(key = external_magnetic_field_frequency_key , value = 'External magnetic field frequency' ) ; check(error==0)
    error = helpers%set(key = external_current_amplitude_key , value = 'External current amplitude') ; check(error==0) 
    error = helpers%set(key = external_current_frequency_key , value = 'External current frequency'); check(error==0)
    error = helpers%set(key = apply_current_density_constraint_key , value = 'Apply current constraint?')  ; check(error==0)
    
    ! PHYSICAL PROPERTIES parameters
    error = helpers%set(key = air_permeability_key , value = 'Air permeability')     ; check(error==0)
    error = helpers%set(key = air_resistivity_key  , value = 'Air resistivity')    ; check(error==0)
    error = helpers%set(key = hts_permeability_key , value = 'High Temperature Superconductor permeability')     ; check(error==0)
    error = helpers%set(key = hts_resistivity_key  , value = 'High Temperature Superconductor resistivity')    ; check(error==0)
    error = helpers%set(key = critical_current_key , value = 'Critical Current [Jc] ')     ; check(error==0)
    error = helpers%set(key = critical_electric_field_key , value = 'Critical Electric Field [Ec]'); check(error==0)
    error = helpers%set(key = nonlinear_exponent_key , value = 'Nonlinear exponent (E-J law)')  ; check(error==0)
    
    ! TIME INTEGRATION parameters
    error = helpers%set(key = theta_value_key        , value = 'Theta value')        ; check(error==0)
    error = helpers%set(key = initial_time_key       , value = 'Initial time')       ; check(error==0)
    error = helpers%set(key = final_time_key         , value = 'Final time')         ; check(error==0)
    error = helpers%set(key = number_time_steps_key  , value = 'Number of steps')  ; check(error==0)
    error = helpers%set(key = is_adaptive_time_stepping_key, value ='Is adaptive time stepping?'); check(error==0)
    error = helpers%set(key = stepping_parameter_key       , value = 'Ideal number of Newton-Raphson nonlinear iterations to converge') ; check(error==0)
    error = helpers%set(key = max_time_step_key            , value = 'Maximum time step size')  ; check(error==0)
    error = helpers%set(key = min_time_step_key            , value = 'Minimum time step size')  ; check(error==0)
    error = helpers%set(key = save_solution_n_steps_key    , value = 'Save solution in N steps ( time interval is divided into N steps to store solution) ');check(error==0)
    
    ! NONLINEAR SOLVER parameters
    error = helpers%set(key = nonlinear_convergence_criteria_key , value = 'Choose one of the convergence criteria: [abs_res_norm,rel_r0_res_norm,rel_rhs_res_norm]' );check(error==0)
    error = helpers%set(key = absolute_nonlinear_tolerance_key   , value = 'Absolute residual ');check(error==0)
    error = helpers%set(key = relative_nonlinear_tolerance_key   , value = 'Relative residual '); check(error==0)
    error = helpers%set(key = max_nonlinear_iterations_key       , value = 'Maximum number of nonlinear iterations allowed');check(error==0)

    
    
    ! REQUIRED attribute ********************************************************************************
    
    ! IO parameters
    error = required%set(key = dir_path_key       , value = .false.)    ; check(error==0)
    error = required%set(key = prefix_key         , value = .false.)    ; check(error==0)
    error = required%set(key = dir_path_out_key   , value = .false.)    ; check(error==0)
    error = required%set(key = write_solution_key , value = .false.)    ; check(error==0)
    
    ! MESHING parameters
    error = required%set(key = number_of_dimensions_key   , value = .false.)  ; check(error==0)
    error = required%set(key = number_of_cells_per_dir_key, value = .false.)  ; check(error==0)
    error = required%set(key = is_dir_periodic_key        , value = .false.)  ; check(error==0)
    error = required%set(key = triangulation_generate_key , value = .false.)  ; check(error==0)
    error = required%set(key = domain_length_key          , value = .false.)  ; check(error==0)
    error = required%set(key = hts_domain_length_key      , value = .false.)  ; check(error==0)
    
    ! FE ELEMENT parameters
    error = required%set(key = magnetic_field_reference_fe_order_key    , value = .false.)  ; check(error==0)
    error = required%set(key = magnetic_pressure_reference_fe_order_key , value = .false.)  ; check(error==0)    
    
     ! PARAMETRIZED ANALYTICAL FUNCTIONS parameters
    error = required%set(key = external_magnetic_field_amplitude_key , value = .false.)  ; check(error==0) 
    error = required%set(key = external_magnetic_field_frequency_key , value = .false.)  ; check(error==0)
    error = required%set(key = external_current_amplitude_key        , value = .false.)  ; check(error==0) 
    error = required%set(key = external_current_frequency_key        , value = .false.)  ; check(error==0)
    error = required%set(key = apply_current_density_constraint_key  , value = .false.)  ; check(error==0)
    
    ! PHYSICAL PROPERTIES parameters
    error = required%set(key = air_permeability_key        , value = .false.)    ; check(error==0)
    error = required%set(key = air_resistivity_key         , value = .false.)    ; check(error==0)
    error = required%set(key = hts_permeability_key        , value = .false.)    ; check(error==0)
    error = required%set(key = hts_resistivity_key         , value = .false.)    ; check(error==0)
    error = required%set(key = critical_current_key        , value = .false.)    ; check(error==0)
    error = required%set(key = critical_electric_field_key , value = .false.)    ; check(error==0)
    error = required%set(key = nonlinear_exponent_key      , value = .false.)    ; check(error==0)
    
    ! TIME INTEGRATION parameters
    error = required%set(key = theta_value_key              , value = .false.)   ; check(error==0)
    error = required%set(key = initial_time_key             , value = .false.)   ; check(error==0)
    error = required%set(key = final_time_key               , value = .false.)   ; check(error==0)
    error = required%set(key = number_time_steps_key        , value = .false.)   ; check(error==0)
    error = required%set(key = is_adaptive_time_stepping_key, value =.false.)    ; check(error==0)
    error = required%set(key = stepping_parameter_key       , value = .false.)   ; check(error==0)
    error = required%set(key = max_time_step_key            , value = .false.)   ; check(error==0)
    error = required%set(key = min_time_step_key            , value = .false.)   ; check(error==0)
    error = required%set(key = save_solution_n_steps_key    , value = .false.)   ;check(error==0)
    
    ! NONLINEAR SOLVER parameters
    error = required%set(key = nonlinear_convergence_criteria_key , value = .false.)  ;check(error==0)
    error = required%set(key = absolute_nonlinear_tolerance_key   , value = .false.)  ;check(error==0)
    error = required%set(key = relative_nonlinear_tolerance_key   , value = .false.)  ;check(error==0)
    error = required%set(key = max_nonlinear_iterations_key       , value = .false.)  ;check(error==0)

  end subroutine hts_nedelec_params_set_default

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
   ! assert(list%isAssignable(dir_path_key, get_dir_path))
    error = list%GetAsString(key = dir_path_key, string = get_dir_path)
    assert(error==0)
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_prefix
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
   ! assert(list%isAssignable(prefix_key, get_prefix))
    error = list%GetAsString(key = prefix_key, string = get_prefix)
    assert(error==0)
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=:),      allocatable            :: get_dir_path_out
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(dir_path_out_key, get_dir_path_out))
    error = list%GetAsString(key = dir_path_out_key, string = get_dir_path_out)
    assert(error==0)
  end function get_dir_path_out

  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical                                       :: get_write_solution
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(write_solution_key, get_write_solution))
    error = list%Get(key=write_solution_key, Value = get_write_solution)
    check(error==0)
  end function get_write_solution

  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_triangulation_type
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(triangulation_generate_key, get_triangulation_type))
    error = list%Get(key=triangulation_generate_key, Value = get_triangulation_type)
    assert(error==0)
  end function get_triangulation_type

   !==================================================================================================
  function get_domain_length(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                      :: get_domain_length(0:SPACE_DIM-1)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(domain_length_key, get_domain_length))
    error = list%Get(key = domain_length_key, Value = get_domain_length)
    assert(error==0)
  end function get_domain_length

   !==================================================================================================
  function get_hts_domain_length(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                      :: get_hts_domain_length(0:SPACE_DIM-1)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(hts_domain_length_key, get_hts_domain_length))
    error = list%Get(key = hts_domain_length_key, Value = get_hts_domain_length)
    assert(error==0)
  end function get_hts_domain_length

    !==================================================================================================
  function get_magnetic_field_reference_fe_order(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_magnetic_field_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(magnetic_field_reference_fe_order_key, get_magnetic_field_reference_fe_order))
    error = list%Get(key =magnetic_field_reference_fe_order_key, Value = get_magnetic_field_reference_fe_order)
    assert(error==0)
  end function get_magnetic_field_reference_fe_order

    !==================================================================================================
  function get_magnetic_pressure_reference_fe_order(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_magnetic_pressure_reference_fe_order
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(magnetic_pressure_reference_fe_order_key, get_magnetic_pressure_reference_fe_order))
    error = list%Get(key=magnetic_pressure_reference_fe_order_key, Value = get_magnetic_pressure_reference_fe_order)
    assert(error==0)
  end function get_magnetic_pressure_reference_fe_order

   !==================================================================================================
  function get_external_magnetic_field_amplitude(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                      :: get_external_magnetic_field_amplitude(0:SPACE_DIM-1)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(external_magnetic_field_amplitude_key, get_external_magnetic_field_amplitude))
    error = list%Get(key=external_magnetic_field_amplitude_key, Value = get_external_magnetic_field_amplitude)
    assert(error==0)
  end function get_external_magnetic_field_amplitude

   !==================================================================================================
  function get_external_magnetic_field_frequency(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_external_magnetic_field_frequency
    type(ParameterList_t), pointer             :: list
    integer(ip)                                :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(external_magnetic_field_frequency_key, get_external_magnetic_field_frequency))
    error = list%Get(key=external_magnetic_field_frequency_key, Value = get_external_magnetic_field_frequency)
    assert(error==0)
  end function get_external_magnetic_field_frequency

  !==================================================================================================
  function get_external_current_amplitude(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                      :: get_external_current_amplitude(0:SPACE_DIM-1)
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(external_current_amplitude_key, get_external_current_amplitude))
    error = list%Get(key=external_current_amplitude_key, Value = get_external_current_amplitude)
    assert(error==0)
  end function get_external_current_amplitude

   !==================================================================================================
  function get_external_current_frequency(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                      :: get_external_current_frequency 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(external_current_frequency_key, get_external_current_frequency))
    error = list%Get(key=external_current_frequency_key, Value = get_external_current_frequency)
    assert(error==0)
  end function get_external_current_frequency

 !==================================================================================================
  function get_apply_current_density_constraint(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical                                   :: get_apply_current_density_constraint
    type(ParameterList_t), pointer            :: list
    integer(ip)                               :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(apply_current_density_constraint_key, get_apply_current_density_constraint))
    error = list%Get(key=apply_current_density_constraint_key, Value = get_apply_current_density_constraint)
    assert(error==0)
  end function get_apply_current_density_constraint

  !==================================================================================================
  function get_air_permeability(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_air_permeability 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(air_permeability_key, get_air_permeability))
    error = list%Get(key=air_permeability_key, Value = get_air_permeability)
    assert(error==0)
  end function get_air_permeability

  !==================================================================================================
  function get_air_resistivity(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_air_resistivity 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(air_resistivity_key, get_air_resistivity))
    error = list%Get(key=air_resistivity_key, Value = get_air_resistivity)
    assert(error==0)
  end function get_air_resistivity


  !==================================================================================================
  function get_hts_permeability(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_hts_permeability 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(hts_permeability_key, get_hts_permeability))
    error = list%Get(key=hts_permeability_key, Value = get_hts_permeability)
    assert(error==0)
  end function get_hts_permeability

  !==================================================================================================
  function get_hts_resistivity(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_hts_resistivity 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(hts_resistivity_key, get_hts_resistivity))
    error = list%Get(key=hts_resistivity_key, Value = get_hts_resistivity)
    assert(error==0)
  end function get_hts_resistivity

  !==================================================================================================
  function get_critical_current(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_critical_current 
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(critical_current_key, get_critical_current))
    error = list%Get(key=critical_current_key, Value = get_critical_current)
    assert(error==0)
  end function get_critical_current 

  !==================================================================================================
  function get_critical_electric_field (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_critical_electric_field  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(critical_electric_field_key, get_critical_electric_field ))
    error = list%Get(key=critical_electric_field_key, Value = get_critical_electric_field )
    assert(error==0)
  end function get_critical_electric_field   

 !==================================================================================================
  function get_nonlinear_exponent  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_nonlinear_exponent   
    type(ParameterList_t), pointer             :: list
    integer(ip)                                :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(nonlinear_exponent_key, get_nonlinear_exponent  ))
    error = list%Get(key=nonlinear_exponent_key, Value = get_nonlinear_exponent  )
    assert(error==0)
  end function get_nonlinear_exponent 

 !==================================================================================================
  function get_theta_value  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                 :: get_theta_value   
    type(ParameterList_t), pointer           :: list
    integer(ip)                              :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(theta_value_key, get_theta_value  ))
    error = list%Get(key=theta_value_key, Value = get_theta_value  ) 
    assert(error==0)
  end function get_theta_value     

 !==================================================================================================
  function get_initial_time  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_initial_time   
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(initial_time_key, get_initial_time  ))
    error = list%Get(key=initial_time_key, Value = get_initial_time  )
    assert(error==0)
  end function get_initial_time   

!==================================================================================================
  function get_final_time  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                   :: get_final_time   
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(final_time_key, get_final_time  ))
    error = list%Get(key=final_time_key, Value = get_final_time  )
    assert(error==0)
  end function get_final_time   

!==================================================================================================
  function get_number_time_steps  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_number_time_steps   
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(number_time_steps_key, get_number_time_steps  ))
    error = list%Get(key=number_time_steps_key, Value = get_number_time_steps  )
    assert(error==0)
  end function get_number_time_steps   

!==================================================================================================
  function get_is_adaptive_time_stepping  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical                                  :: get_is_adaptive_time_stepping   
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(is_adaptive_time_stepping_key, get_is_adaptive_time_stepping  ))
    error = list%Get(key=is_adaptive_time_stepping_key, Value = get_is_adaptive_time_stepping  )
    assert(error==0)
  end function get_is_adaptive_time_stepping  

!==================================================================================================
  function get_stepping_parameter (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_stepping_parameter  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(stepping_parameter_key, get_stepping_parameter ))
    error = list%Get(key=stepping_parameter_key, Value = get_stepping_parameter )
    assert(error==0)
  end function get_stepping_parameter

!==================================================================================================
  function get_max_time_step (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                  :: get_max_time_step  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(max_time_step_key, get_max_time_step ))
    error = list%Get(key=max_time_step_key, Value = get_max_time_step )
    assert(error==0)
  end function get_max_time_step 

!==================================================================================================
  function get_min_time_step (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                  :: get_min_time_step  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(min_time_step_key, get_min_time_step ))
    error = list%Get(key=min_time_step_key, Value = get_min_time_step )
    assert(error==0)
  end function get_min_time_step

!==================================================================================================
  function get_save_solution_n_steps (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                 :: get_save_solution_n_steps  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(save_solution_n_steps_key, get_save_solution_n_steps ))
    error = list%Get(key=save_solution_n_steps_key, Value = get_save_solution_n_steps )
    assert(error==0)
  end function get_save_solution_n_steps 

!==================================================================================================
  function get_nonlinear_convergence_criteria (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=:)     , allocatable       :: get_nonlinear_convergence_criteria  
    type(ParameterList_t), pointer           :: list
    integer(ip)                              :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(nonlinear_convergence_criteria_key, get_nonlinear_convergence_criteria ))
    error = list%GetAsString(key=nonlinear_convergence_criteria_key, string = get_nonlinear_convergence_criteria )
    assert(error==0)
  end function get_nonlinear_convergence_criteria  

!==================================================================================================
  function get_absolute_nonlinear_tolerance (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                  :: get_absolute_nonlinear_tolerance  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(absolute_nonlinear_tolerance_key, get_absolute_nonlinear_tolerance ))
    error = list%Get(key=absolute_nonlinear_tolerance_key, Value = get_absolute_nonlinear_tolerance )
    assert(error==0)
  end function get_absolute_nonlinear_tolerance

!==================================================================================================
  function get_relative_nonlinear_tolerance (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp)                                  :: get_relative_nonlinear_tolerance  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(relative_nonlinear_tolerance_key, get_relative_nonlinear_tolerance ))
    error = list%Get(key=relative_nonlinear_tolerance_key, Value = get_relative_nonlinear_tolerance )
    assert(error==0)
  end function get_relative_nonlinear_tolerance

!==================================================================================================
  function get_max_nonlinear_iterations (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip)                                   :: get_max_nonlinear_iterations  
    type(ParameterList_t), pointer                :: list
    integer(ip)                                   :: error
    list  => this%get_parameters()
    !assert(list%isAssignable(max_nonlinear_iterations_key, get_max_nonlinear_iterations ))
    error = list%Get(key=max_nonlinear_iterations_key, Value = get_max_nonlinear_iterations )
    assert(error==0)
  end function get_max_nonlinear_iterations 
  

end module hts_nedelec_params_names
