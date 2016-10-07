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
module hts_nedelec_params_names
  use fempar_names
# include "debug.i90"

  implicit none
  private

  type hts_nedelec_params_t 
     private 
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     ! FE element
     character(len=:), allocatable :: default_reference_fe_geo_order
     character(len=:), allocatable :: default_reference_fe_order
     ! Structured meshing 
     character(len=:), allocatable :: default_triangulation_type
     character(len=:), allocatable :: default_num_dimensions
     character(len=:), allocatable :: default_nx
     character(len=:), allocatable :: default_ny
     character(len=:), allocatable :: default_nz
     character(len=:), allocatable :: default_is_periodic_in_x
     character(len=:), allocatable :: default_is_periodic_in_y
     character(len=:), allocatable :: default_is_periodic_in_z
     character(len=:), allocatable :: default_write_solution
     character(len=:), allocatable :: default_domain_length_lx
     character(len=:), allocatable :: default_domain_length_ly
     character(len=:), allocatable :: default_domain_length_lz
     character(len=:), allocatable :: default_hts_domain_length_lx
     character(len=:), allocatable :: default_hts_domain_length_ly
     character(len=:), allocatable :: default_hts_domain_length_lz
     ! Customized Problem conditions and source term  
     character(len=:), allocatable :: default_external_magnetic_field_frequency 
     character(len=:), allocatable :: default_external_current_frequency
     character(len=:), allocatable :: default_external_magnetic_field_Hx 
     character(len=:), allocatable :: default_external_magnetic_field_Hy
     character(len=:), allocatable :: default_external_magnetic_field_Hz
     character(len=:), allocatable :: default_external_current_Jx
     character(len=:), allocatable :: default_external_current_Jy
     character(len=:), allocatable :: default_external_current_Jz
     character(len=:), allocatable :: default_apply_current_density_constraint 
     ! Physical properties 
     character(len=:), allocatable :: default_air_permeability
     character(len=:), allocatable :: default_air_resistivity
     character(len=:), allocatable :: default_hts_permeability
     character(len=:), allocatable :: default_hts_resistivity
     character(len=:), allocatable :: default_critical_current 
     character(len=:), allocatable :: default_critical_electric_field
     character(len=:), allocatable :: default_nonlinear_exponent
     ! Time integration 
     character(len=:), allocatable :: default_theta_value 
     character(len=:), allocatable :: default_initial_time 
     character(len=:), allocatable :: default_final_time 
     character(len=:), allocatable :: default_number_of_steps
     character(len=:), allocatable :: default_is_adaptive_time_stepping 
     character(len=:), allocatable :: default_stepping_parameter
     character(len=:), allocatable :: default_max_time_step 
     character(len=:), allocatable :: default_min_time_step 
     ! Nonlinear solver tolerance 
     character(len=:), allocatable :: default_absolute_nonlinear_tolerance
     character(len=:), allocatable :: default_relative_nonlinear_tolerance 
     character(len=:), allocatable :: default_max_nonlinear_iterations 

     type(Command_Line_Interface):: cli 

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     character(len=256)            :: dir_path_out
     integer(ip)                   :: reference_fe_geo_order
     integer(ip)                   :: reference_fe_order
     character(len=256)            :: triangulation_type
     integer(ip)                   :: num_dimensions     
     integer(ip)                   :: number_of_cells_per_dir(0:SPACE_DIM-1)
     integer(ip)                   :: is_dir_periodic(0:SPACE_DIM-1)
     logical                       :: write_solution
     real(rp)                      :: domain_length(0:SPACE_DIM-1)
     real(rp)                      :: hts_domain_length(0:SPACE_DIM-1) 
     ! Customized Problem conditions and source term  
     integer(ip)                   :: external_magnetic_field_frequency 
     integer(ip)                   :: external_current_frequency 
     real(rp)                      :: external_magnetic_field(3) 
     real(rp)                      :: external_current(3) 
     logical                       :: apply_current_density_constraint 
     ! Physical properties 
     real(rp)                      :: air_permeability           ! [ H/m ]
     real(rp)                      :: air_resistivity            ! [ Ohm·m ]
     real(rp)                      :: hts_permeability           ! [ H/m ]
     real(rp)                      :: hts_resistivity            ! [ Ohm·m ]
     real(rp)                      :: critical_current           ! [ A/m² ]
     real(rp)                      :: critical_electric_field    ! [ V/m ]
     real(rp)                      :: nonlinear_exponent
     ! Time integration 
     real(rp)                      :: theta_value 
     real(rp)                      :: initial_time 
     real(rp)                      :: final_time 
     integer(ip)                   :: number_of_steps
     logical                       :: is_adaptive_time_stepping 
     integer(ip)                   :: stepping_parameter
     real(rp)                      :: max_time_step 
     real(rp)                      :: min_time_step 
     ! Nonlinear solver tolerance 
     real(rp)                      :: absolute_nonlinear_tolerance
     real(rp)                      :: relative_nonlinear_tolerance 
     real(rp)                      :: max_nonlinear_iterations 

   contains
     procedure, non_overridable             :: create       => hts_nedelec_create
     procedure, non_overridable, private    :: set_default  => hts_nedelec_set_default
     procedure, non_overridable, private    :: add_to_cli   => hts_nedelec_add_to_cli
     procedure, non_overridable             :: parse        => hts_nedelec_parse 
     procedure, non_overridable             :: free         => hts_nedelec_free
     procedure, non_overridable             :: get_dir_path
     procedure, non_overridable             :: get_prefix
     procedure, non_overridable             :: get_dir_path_out
     procedure, non_overridable             :: get_reference_fe_geo_order
     procedure, non_overridable             :: get_reference_fe_order
     procedure, non_overridable             :: get_triangulation_type
     procedure, non_overridable             :: get_write_solution
     procedure, non_overridable             :: get_domain_length 
     procedure, non_overridable             :: get_hts_domain_length
     procedure, non_overridable             :: get_external_magnetic_field_frequency
     procedure, non_overridable             :: get_external_current_frequency
     procedure, non_overridable             :: get_external_magnetic_field
     procedure, non_overridable             :: get_external_current
     procedure, non_overridable             :: get_apply_current_density_constraint
     procedure, non_overridable             :: get_air_permeability
     procedure, non_overridable             :: get_air_resistivity
     procedure, non_overridable             :: get_hts_permeability
     procedure, non_overridable             :: get_hts_resistivity
     procedure, non_overridable 	    :: get_critical_current           
     procedure, non_overridable             :: get_critical_electric_field   
     procedure, non_overridable             :: get_nonlinear_exponent
     procedure, non_overridable             :: get_theta_value 
     procedure, non_overridable             :: get_initial_time 
     procedure, non_overridable             :: get_final_time 
     procedure, non_overridable             :: get_number_of_steps
     procedure, non_overridable             :: get_is_adaptive_time_stepping 
     procedure, non_overridable             :: get_stepping_parameter
     procedure, non_overridable             :: get_max_time_step 
     procedure, non_overridable             :: get_min_time_step 
     procedure, non_overridable             :: get_absolute_nonlinear_tolerance
     procedure, non_overridable             :: get_relative_nonlinear_tolerance 
     procedure, non_overridable             :: get_max_nonlinear_iterations 
  end type hts_nedelec_params_t

  ! Types
  public :: hts_nedelec_params_t

contains

  subroutine hts_nedelec_create(this)
    implicit none
    class(hts_nedelec_params_t), intent(inout) :: this
    
    call this%free()
    
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'HTS_nedelec',                                                   &
         &        version     = '',                                                                   &
         &        authors     = '',                                                                   &
         &        license     = '',                                                                   &
         &        description =  'FEMPAR test to solve the High Temperature Superconductivity problem &
                                  proposed in the FORTISSIMO EXPERIMENT.',                            &
         &        examples    = ['test_hts_nedelec -h  ', 'test_hts_nedelec -h  ' ])
    
    call this%set_default()
    call this%add_to_cli()
  end subroutine hts_nedelec_create
  
  subroutine hts_nedelec_set_default(this)
    implicit none
    class(hts_nedelec_params_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = 'data/'
    this%default_prefix         = 'square'
    this%default_dir_path_out   = 'output/'
    this%default_reference_fe_geo_order = '1'
    this%default_reference_fe_order = '1'
    this%default_triangulation_type = 'structured'
    this%default_num_dimensions = '2'
    this%default_nx = '1'
    this%default_ny = '1'
    this%default_nz = '1'
    this%default_is_periodic_in_x = '0'
    this%default_is_periodic_in_y = '0'
    this%default_is_periodic_in_z = '0'
    this%default_write_solution = '.false.'
    this%default_domain_length_lx = '1.0'
    this%default_domain_length_ly = '1.0'
    this%default_domain_length_lz = '1.0'
    this%default_hts_domain_length_lx = '0.5'
    this%default_hts_domain_length_ly = '0.2'
    this%default_hts_domain_length_lz = '1.0'
    ! Customized Problem conditions and source term  
    this%default_external_magnetic_field_frequency = '50'
    this%default_external_current_frequency        = '50'
    this%default_external_magnetic_field_Hx        = '0'
    this%default_external_magnetic_field_Hy        = '1e6'
    this%default_external_magnetic_field_Hz        = '0'
    this%default_external_current_Jx                  = '0'
    this%default_external_current_Jy                  = '0'
    this%default_external_current_Jz                  = '0'
    this%default_apply_current_density_constraint  = '.false.'
    ! Physical properties 
    this%default_air_permeability        = '1.257e-6'  ! [ H/m ]
    this%default_air_resistivity         = '1e4'       ! [ Ohm·m ]
    this%default_hts_permeability        = '1.0'       ! [ H/m ]
    this%default_hts_resistivity         = '1.0'       ! [ Ohm·m ]
    this%default_critical_current        = '4.08e8'    ! [ A/m² ]
    this%default_critical_electric_field = '1e-4'      ! [ V/m ]
    this%default_nonlinear_exponent      = '0'
    ! Time integration 
    this%default_theta_value                = '1.0'
    this%default_initial_time               = '0.0'
    this%default_final_time                 = '1.0'
    this%default_number_of_steps            = '10'
    this%default_is_adaptive_time_stepping  = '.false.'
    this%default_stepping_parameter         = '20'
    this%default_max_time_step              = '0.1'
    this%default_min_time_step              = '1e-3'
    ! Nonlinear solver tolerance 
    this%default_absolute_nonlinear_tolerance       = '1e-3'
    this%default_relative_nonlinear_tolerance       = '1e-3'
    this%default_max_nonlinear_iterations           = '500'
  end subroutine hts_nedelec_set_default
  
  !==================================================================================================
  subroutine hts_nedelec_add_to_cli(this)
    implicit none
    class(hts_nedelec_params_t) , intent(inout) :: this

    ! Locals
    integer(ip) :: error

    ! IO parameters
    call this%cli%add(switch='--dir-path',switch_ab='-d',                                   &
         &            help='Directory of the source files',required=.false., act='store',   &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Name of the GiD files',  &
         &            required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--dir-path-out',switch_ab='-o',help='Output Directory',&
         &            required=.false.,act='store',def=trim(this%default_dir_path_out),error=error)
    check(error==0)
    call this%cli%add(switch='--reference-fe-geo-order',switch_ab='-gorder',help='Order of the triangulation reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_geo_order),error=error)
    check(error==0)  
    call this%cli%add(switch='--reference-fe-order',switch_ab='-order',help='Order of the fe space reference fe',&
         &            required=.false.,act='store',def=trim(this%default_reference_fe_order),error=error) 
    check(error==0) 
    call this%cli%add(switch='--triangulation-type',switch_ab='-tt',help='Structured or unstructured (GiD) triangulation?',&
         &            required=.false.,act='store',def=trim(this%default_triangulation_type),choices='structured,unstructured',error=error) 
    check(error==0) 
    call this%cli%add(switch='--number_of_dimensions',switch_ab='-dim',help='Number of space dimensions',&
         &            required=.false.,act='store',def=trim(this%default_num_dimensions),error=error) 
    check(error==0) 
    call this%cli%add(switch='--number_of_cells_in_x',switch_ab='-nx',help='Number of cells in x',&
         &            required=.false.,act='store',def=trim(this%default_nx),error=error) 
    check(error==0) 
    call this%cli%add(switch='--number_of_cells_in_y',switch_ab='-ny',help='Number of cells in y',&
         &            required=.false.,act='store',def=trim(this%default_ny),error=error) 
    check(error==0) 
    call this%cli%add(switch='--number_of_cells_in_z',switch_ab='-nz',help='Number of cells in z',&
         &            required=.false.,act='store',def=trim(this%default_nz),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_x',switch_ab='-px',help='Is the mesh periodic in x',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_x),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_y',switch_ab='-py',help='Is the mesh periodic in y',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_y),error=error) 
    check(error==0) 
    call this%cli%add(switch='--periodic_in_z',switch_ab='-pz',help='Is the mesh periodic in z',&
         &            required=.false.,act='store',def=trim(this%default_is_periodic_in_z),error=error) 
    check(error==0)
    call this%cli%add(switch='--write-solution',switch_ab='-wsolution',help='Write solution in VTK format',&
         &            required=.false.,act='store',def=trim(this%default_write_solution),error=error) 
    check(error==0)  
    call this%cli%add(switch='--domain_length_in_x',switch_ab='-dlx',help='Length of domain in X direction',&
         &            required=.false.,act='store',def=trim(this%default_domain_length_lx),error=error) 
    check(error==0)
    call this%cli%add(switch='--domain_length_in_y',switch_ab='-dly',help='Length of domain in Y direction',&
         &            required=.false.,act='store',def=trim(this%default_domain_length_ly),error=error) 
    check(error==0)
    call this%cli%add(switch='--domain_length_in_z',switch_ab='-dlz',help='Length of domain in Z direction',&
         &            required=.false.,act='store',def=trim(this%default_domain_length_lz),error=error) 
    check(error==0)
    call this%cli%add(switch='--hts_domain_length_in_x',switch_ab='-hts_dlx',help='Length of the HTS domain in X direction',&
         &            required=.false.,act='store',def=trim(this%default_hts_domain_length_lx),error=error) 
    check(error==0)
    call this%cli%add(switch='--hts_domain_length_in_y',switch_ab='-hts_dly',help='Length of the HTS domain in Y direction',&
         &            required=.false.,act='store',def=trim(this%default_hts_domain_length_ly),error=error) 
    check(error==0)
    call this%cli%add(switch='--hts_domain_length_in_z',switch_ab='-hts_dlz',help='Length of the HTS domain in Z direction',&
         &            required=.false.,act='store',def=trim(this%default_hts_domain_length_lz),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_magnetic_field_w',switch_ab='-w_H',help='Frequency of the applied external H',&
         &            required=.false.,act='store',def=trim(this%default_external_magnetic_field_frequency),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_current_w',switch_ab='-w_J',help='Frequency of the applied external J',&
         &            required=.false.,act='store',def=trim(this%default_external_current_frequency),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_magnetic_field_x_component',switch_ab='-Hx',help='Hx Applied external magnetic field',&
         &            required=.false.,act='store',def=trim(this%default_external_magnetic_field_Hx),error=error)
    check(error==0)
    call this%cli%add(switch='--external_magnetic_field_y_component',switch_ab='-Hy',help='Hy Applied external magnetic field',&
         &            required=.false.,act='store',def=trim(this%default_external_magnetic_field_Hy),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_magnetic_field_z_compoentn',switch_ab='-Hz',help='Hz Applied external magnetic field',&
         &            required=.false.,act='store',def=trim(this%default_external_magnetic_field_Hz),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_current_in_x',switch_ab='-Jx',help='Applied external current in X direction',&
         &            required=.false.,act='store',def=trim(this%default_external_current_Jx),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_current_in_y',switch_ab='-Jy',help='Applied external current in Y direction',&
         &            required=.false.,act='store',def=trim(this%default_external_current_Jy),error=error) 
    check(error==0)
    call this%cli%add(switch='--external_current_in_z',switch_ab='-Jz',help='Applied external current in Z direction',&
         &            required=.false.,act='store',def=trim(this%default_external_current_Jz),error=error) 
    check(error==0)
    call this%cli%add(switch='--apply_current_density_constraint',switch_ab='-cdc',help='Apply current constraint over HTS?',&
         &            required=.false.,act='store',def=trim(this%default_apply_current_density_constraint ),error=error) 
    check(error==0)
    call this%cli%add(switch='--air_permeability',switch_ab='-mu_air',help='Air permeability value',&
         &            required=.false.,act='store',def=trim(this%default_air_permeability ),error=error) 
    check(error==0)
    call this%cli%add(switch='--air_resistivity',switch_ab='-rho_air',help='Air resistivity value',&
         &            required=.false.,act='store',def=trim(this%default_air_resistivity ),error=error) 
    check(error==0)
    call this%cli%add(switch='--hts_permeability',switch_ab='-mu_hts',help='HTS permeability value',&
         &            required=.false.,act='store',def=trim(this%default_hts_permeability ),error=error) 
    check(error==0)
    call this%cli%add(switch='--hts_resistivity',switch_ab='-rho_hts',help='HTS resistivity value',&
         &            required=.false.,act='store',def=trim(this%default_hts_resistivity ),error=error) 
    check(error==0)
    call this%cli%add(switch='--critical_current',switch_ab='-Jc',help='Critical current Jc',&
         &            required=.false.,act='store',def=trim(this%default_critical_current ),error=error) 
    check(error==0)
    call this%cli%add(switch='--critical_electric_field',switch_ab='-Ec',help='Critical Electric field Ec',&
         &            required=.false.,act='store',def=trim(this%default_critical_electric_field ),error=error) 
    check(error==0)
    call this%cli%add(switch='--nonlinear_exponent',switch_ab='-nl_exp',help='Nonlinear exponent in the E-J law for HTS',&
         &            required=.false.,act='store',def=trim(this%default_nonlinear_exponent),error=error) 
    check(error==0)
    call this%cli%add(switch='--theta_value',switch_ab='-theta',help='Theta value in the Theta-method time integration',&
         &            required=.false.,act='store',def=trim(this%default_theta_value),error=error) 
    check(error==0)
    call this%cli%add(switch='--initial_time',switch_ab='-t0',help='Initial time for the simulation',&
         &            required=.false.,act='store',def=trim(this%default_initial_time),error=error) 
    check(error==0)
    call this%cli%add(switch='--final_time',switch_ab='-tf',help='Final time for the simulation',&
         &            required=.false.,act='store',def=trim(this%default_final_time),error=error) 
    check(error==0)
    call this%cli%add(switch='--number_of_steps',switch_ab='-nsteps',help='Number of steps in a regular partition of the time interval',&
         &            required=.false.,act='store',def=trim(this%default_number_of_steps),error=error) 
    check(error==0)
    call this%cli%add(switch='--is_adaptive_time_stepping',switch_ab='-iats',help='Enable adaptive time stepping technique?',&
         &            required=.false.,act='store',def=trim(this%default_is_adaptive_time_stepping),error=error) 
    check(error==0)
    call this%cli%add(switch='--time_stepping_parameter',switch_ab='-tsp',help='Number of NL iterations for which the time step remains unaltered ',&
         &            required=.false.,act='store',def=trim(this%default_stepping_parameter ),error=error) 
    check(error==0)
    call this%cli%add(switch='--max_time_step',switch_ab='-max_ts',help='Maximum time step length allowed ',&
         &            required=.false.,act='store',def=trim(this%default_max_time_step),error=error) 
    check(error==0)
    call this%cli%add(switch='--min_time_step',switch_ab='-min_ts',help='Minimum time step lenght allowed ',&
         &            required=.false.,act='store',def=trim(this%default_min_time_step),error=error) 
    check(error==0)
    call this%cli%add(switch='--abs_nl_tolerance',switch_ab='-abs_tol',help=' Absolute tolerance for the NL solver ',&
         &            required=.false.,act='store',def=trim(this%default_absolute_nonlinear_tolerance),error=error) 
    check(error==0)
    call this%cli%add(switch='--rel_nl_tolerance',switch_ab='-rel_tol',help=' Relative tolerance for the NL solver ',&
         &            required=.false.,act='store',def=trim(this%default_relative_nonlinear_tolerance),error=error) 
    check(error==0)
    call this%cli%add(switch='--max_NL_iterations',switch_ab='-max_nl_its',help=' Maximum nonlinear iterations allowed ',&
         &            required=.false.,act='store',def=trim(this%default_max_nonlinear_iterations),error=error) 
    check(error==0)

  end subroutine hts_nedelec_add_to_cli
  
  subroutine hts_nedelec_parse(this,parameter_list)
    implicit none
    class(hts_nedelec_params_t), intent(inout) :: this
    type(ParameterList_t)       , intent(inout) :: parameter_list
    integer(ip) :: istat
    
    call this%cli%parse(error=istat); check(istat==0)
    
    ! IO parameters
    call this%cli%get(switch='-d',val=this%dir_path    ,error=istat); check(istat==0)
    call this%cli%get(switch='-p',val=this%prefix      ,error=istat); check(istat==0)
    call this%cli%get(switch='-o',val=this%dir_path_out,error=istat); check(istat==0)
    call this%cli%get(switch='-gorder',val=this%reference_fe_geo_order,error=istat); check(istat==0)
    call this%cli%get(switch='-order',val=this%reference_fe_order,error=istat); check(istat==0)
    call this%cli%get(switch='-tt',val=this%triangulation_type,error=istat); check(istat==0)
    call this%cli%get(switch='-dim',val=this%num_dimensions,error=istat); check(istat==0)
    call this%cli%get(switch='-nx',val=this%number_of_cells_per_dir(0),error=istat); check(istat==0)
    call this%cli%get(switch='-ny',val=this%number_of_cells_per_dir(1),error=istat); check(istat==0)
    call this%cli%get(switch='-nz',val=this%number_of_cells_per_dir(2),error=istat); check(istat==0)
    call this%cli%get(switch='-px',val=this%is_dir_periodic(0),error=istat); check(istat==0)
    call this%cli%get(switch='-py',val=this%is_dir_periodic(1),error=istat); check(istat==0)
    call this%cli%get(switch='-pz',val=this%is_dir_periodic(2),error=istat); check(istat==0)
    call this%cli%get(switch='-wsolution',val=this%write_solution,error=istat); check(istat==0)
    call this%cli%get(switch='-dlx',val=this%domain_length(0),error=istat); check(istat==0)
    call this%cli%get(switch='-dly',val=this%domain_length(1),error=istat); check(istat==0)
    call this%cli%get(switch='-dlz',val=this%domain_length(2),error=istat); check(istat==0)
    call this%cli%get(switch='-hts_dlx',val=this%hts_domain_length(0),error=istat); check(istat==0)
    call this%cli%get(switch='-hts_dly',val=this%hts_domain_length(1),error=istat); check(istat==0)
    call this%cli%get(switch='-hts_dlz',val=this%hts_domain_length(2),error=istat); check(istat==0)
    call this%cli%get(switch='-w_H',val=this%external_magnetic_field_frequency,error=istat); check(istat==0)
    call this%cli%get(switch='-w_J',val=this%external_current_frequency,error=istat); check(istat==0)
    call this%cli%get(switch='-Hx',val=this%external_magnetic_field(1),error=istat); check(istat==0)
    call this%cli%get(switch='-Hy',val=this%external_magnetic_field(2),error=istat); check(istat==0)
    call this%cli%get(switch='-Hz',val=this%external_magnetic_field(3),error=istat); check(istat==0)
    call this%cli%get(switch='-Jx',val=this%external_current(1),error=istat); check(istat==0)
    call this%cli%get(switch='-Jy',val=this%external_current(2),error=istat); check(istat==0)
    call this%cli%get(switch='-Jz',val=this%external_current(3),error=istat); check(istat==0)
    call this%cli%get(switch='-cdc',val=this%apply_current_density_constraint,error=istat); check(istat==0)
    call this%cli%get(switch='-mu_air',val=this%air_permeability,error=istat); check(istat==0)
    call this%cli%get(switch='-rho_air',val=this%air_resistivity,error=istat); check(istat==0)
    call this%cli%get(switch='-mu_hts',val=this%hts_permeability,error=istat); check(istat==0)
    call this%cli%get(switch='-rho_hts',val=this%hts_resistivity,error=istat); check(istat==0)
    call this%cli%get(switch='-Jc',val=this%critical_current,error=istat); check(istat==0)
    call this%cli%get(switch='-Ec',val=this%critical_electric_field,error=istat); check(istat==0)
    call this%cli%get(switch='-nl_exp',val=this%nonlinear_exponent,error=istat); check(istat==0)
    call this%cli%get(switch='-theta',val=this%theta_value,error=istat); check(istat==0)
    call this%cli%get(switch='-t0',val=this%initial_time,error=istat); check(istat==0)
    call this%cli%get(switch='-tf',val=this%final_time,error=istat); check(istat==0)
    call this%cli%get(switch='-nsteps',val=this%number_of_steps,error=istat); check(istat==0)
    call this%cli%get(switch='-iats',val=this%is_adaptive_time_stepping,error=istat); check(istat==0)
    call this%cli%get(switch='-tsp',val=this%stepping_parameter,error=istat); check(istat==0)
    call this%cli%get(switch='-max_ts',val=this%max_time_step,error=istat); check(istat==0)
    call this%cli%get(switch='-min_ts',val=this%min_time_step,error=istat); check(istat==0)
    call this%cli%get(switch='-abs_tol',val=this%absolute_nonlinear_tolerance,error=istat); check(istat==0)
    call this%cli%get(switch='-rel_tol',val=this%relative_nonlinear_tolerance,error=istat); check(istat==0)
    call this%cli%get(switch='-max_nl_its',val=this%max_nonlinear_iterations,error=istat); check(istat==0)

    call parameter_list%init()
    istat = 0
    istat = istat + parameter_list%set(key = dir_path_key, value = this%dir_path)
    istat = istat + parameter_list%set(key = prefix_key  , value = this%prefix)
    istat = istat + parameter_list%set(key = geometry_interpolation_order_key  , value = this%reference_fe_geo_order)
    check(istat==0)
    
    if(trim(this%triangulation_type)=='unstructured') then
       istat = parameter_list%set(key = triangulation_generate_key, value = triangulation_generate_from_mesh)
    else if(trim(this%triangulation_type)=='structured') then
       istat = parameter_list%set(key = triangulation_generate_key         , value = triangulation_generate_structured)
       istat = istat + parameter_list%set(key = number_of_dimensions_key   , value = this%num_dimensions)
       istat = istat + parameter_list%set(key = number_of_cells_per_dir_key, value = this%number_of_cells_per_dir)
       istat = istat + parameter_list%set(key = is_dir_periodic_key        , value = this%is_dir_periodic)
    end if
    check(istat==0)
    
  end subroutine hts_nedelec_parse  

  subroutine hts_nedelec_free(this)
    implicit none
    class(hts_nedelec_params_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_dir_path_out)) deallocate(this%default_dir_path_out)
    if(allocated(this%default_reference_fe_geo_order)) deallocate(this%default_reference_fe_geo_order)
    if(allocated(this%default_reference_fe_order)) deallocate(this%default_reference_fe_order)
    call this%cli%free()
  end subroutine hts_nedelec_free

  ! GETTERS *****************************************************************************************
  function get_dir_path(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path
    get_dir_path = this%dir_path
  end function get_dir_path

  !==================================================================================================
  function get_prefix(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=256) :: get_prefix
    get_prefix = this%prefix
  end function get_prefix

  !==================================================================================================
  function get_dir_path_out(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    character(len=256) :: get_dir_path_out
    get_dir_path_out = this%dir_path_out
  end function get_dir_path_out
  
  !==================================================================================================
  function get_reference_fe_geo_order(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_geo_order
    get_reference_fe_geo_order = this%reference_fe_geo_order
  end function get_reference_fe_geo_order
  
  !==================================================================================================
  function get_reference_fe_order(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip) :: get_reference_fe_order
    get_reference_fe_order = this%reference_fe_order
  end function get_reference_fe_order

  !==================================================================================================
  function get_write_solution(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical :: get_write_solution
    get_write_solution = this%write_solution
  end function get_write_solution
  
  !==================================================================================================
  function get_triangulation_type(this)
    implicit none
    class(hts_nedelec_params_t) , target, intent(in) :: this
    character(:), pointer :: get_triangulation_type
    get_triangulation_type => this%triangulation_type
  end function get_triangulation_type
  
    !==================================================================================================
  function get_domain_length(this)
    implicit none
    class(hts_nedelec_params_t) , target, intent(in) :: this
    real(rp), pointer :: get_domain_length(:)
    get_domain_length => this%domain_length
  end function get_domain_length
  
      !==================================================================================================
  function get_hts_domain_length(this)
    implicit none
    class(hts_nedelec_params_t) , target, intent(in) :: this
    real(rp), pointer :: get_hts_domain_length(:)
    get_hts_domain_length => this%hts_domain_length
  end function get_hts_domain_length
  
    !==================================================================================================
  function get_external_magnetic_field_frequency(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_external_magnetic_field_frequency
    get_external_magnetic_field_frequency = this%external_magnetic_field_frequency
  end function get_external_magnetic_field_frequency
  
   !==================================================================================================
  function get_external_current_frequency(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_external_current_frequency
    get_external_current_frequency = this%external_current_frequency
  end function get_external_current_frequency
  
     !==================================================================================================
  function get_external_magnetic_field(this)
    implicit none
    class(hts_nedelec_params_t) , target, intent(in) :: this
    real(rp), pointer :: get_external_magnetic_field(:)
    get_external_magnetic_field => this%external_magnetic_field
  end function get_external_magnetic_field
  
       !==================================================================================================
  function get_external_current(this)
    implicit none
    class(hts_nedelec_params_t) , target, intent(in) :: this
    real(rp), pointer :: get_external_current(:)
    get_external_current => this%external_current
  end function get_external_current
  
        !==================================================================================================
  function get_apply_current_density_constraint(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical :: get_apply_current_density_constraint
    get_apply_current_density_constraint = this%apply_current_density_constraint
  end function get_apply_current_density_constraint
  
      !==================================================================================================
  function get_air_permeability(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_air_permeability
    get_air_permeability = this%air_permeability
  end function get_air_permeability
  
       !==================================================================================================
  function get_air_resistivity(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_air_resistivity
    get_air_resistivity = this%air_resistivity
  end function get_air_resistivity
  
      !==================================================================================================
  function get_hts_permeability(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_hts_permeability
    get_hts_permeability = this%hts_permeability
  end function get_hts_permeability
  
      !==================================================================================================
  function get_hts_resistivity(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_hts_resistivity
    get_hts_resistivity = this%hts_resistivity
  end function get_hts_resistivity
  
      !==================================================================================================
  function get_critical_current(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_critical_current
    get_critical_current = this%critical_current
  end function get_critical_current
  
      !==================================================================================================
  function get_critical_electric_field(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_critical_electric_field
    get_critical_electric_field = this%critical_electric_field
  end function get_critical_electric_field
  
      !==================================================================================================
  function get_nonlinear_exponent(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_nonlinear_exponent
    get_nonlinear_exponent = this%nonlinear_exponent
  end function get_nonlinear_exponent
  
      !==================================================================================================
  function get_theta_value(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_theta_value
    get_theta_value = this%theta_value
  end function get_theta_value
  
      !==================================================================================================
  function get_initial_time(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_initial_time 
    get_initial_time  = this%initial_time 
  end function get_initial_time 
  
      !==================================================================================================
  function get_final_time(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_final_time 
    get_final_time  = this%final_time 
  end function get_final_time 
  
      !==================================================================================================
  function get_number_of_steps(this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip) :: get_number_of_steps 
    get_number_of_steps  = this%number_of_steps
  end function get_number_of_steps 
  
     !==================================================================================================
  function get_is_adaptive_time_stepping  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    logical :: get_is_adaptive_time_stepping 
    get_is_adaptive_time_stepping   = this%is_adaptive_time_stepping 
  end function get_is_adaptive_time_stepping 
  
      !==================================================================================================
  function get_stepping_parameter  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip) :: get_stepping_parameter 
    get_stepping_parameter  = this%stepping_parameter
  end function get_stepping_parameter 
  
      !==================================================================================================
  function get_max_time_step  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_max_time_step 
    get_max_time_step  = this%max_time_step
  end function get_max_time_step
  
      !==================================================================================================
  function get_min_time_step  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_min_time_step 
    get_min_time_step  = this%min_time_step
  end function get_min_time_step
  
      !==================================================================================================
  function get_absolute_nonlinear_tolerance  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_absolute_nonlinear_tolerance 
    get_absolute_nonlinear_tolerance  = this%absolute_nonlinear_tolerance
  end function get_absolute_nonlinear_tolerance
  
       !==================================================================================================
  function get_relative_nonlinear_tolerance  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    real(rp) :: get_relative_nonlinear_tolerance
    get_relative_nonlinear_tolerance  = this%relative_nonlinear_tolerance
  end function get_relative_nonlinear_tolerance
  
        !==================================================================================================
  function get_max_nonlinear_iterations  (this)
    implicit none
    class(hts_nedelec_params_t) , intent(in) :: this
    integer(ip) :: get_max_nonlinear_iterations 
    get_max_nonlinear_iterations  = this%max_nonlinear_iterations
  end function get_max_nonlinear_iterations

end module hts_nedelec_params_names
