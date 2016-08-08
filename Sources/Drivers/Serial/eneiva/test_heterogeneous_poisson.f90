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

! This driver solves a SCALAR x VECTOR LAPLACIAN with constant STRONG BOUNDARY  
! conditions and HETEROGENEOUS material parameters, using a fixed-point (Picard) 
! iterative scheme. The input data is read from the output files generated 
! by COMET GiD problem types. COMET is the in-house code developed by Prof.  
! Michele Chiumenti to solve thermo-mechanical problems with finite elements.

module command_line_parameters_names
  use types_names
  use FLAP
# include "debug.i90"

  implicit none
  private

  type test_heterogeneous_poisson_params_t
     
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     
     ! Nonlinear loop parameters
     character(len=:), allocatable :: default_tolerance
     character(len=:), allocatable :: default_max_number_of_iterations
					
     ! Material parameters
     character(len=:), allocatable :: default_number_of_values
     character(len=:), allocatable :: default_temperature_points
     character(len=:), allocatable :: default_conductivity_values
     
  end type test_heterogeneous_poisson_params_t

  ! Types
  public :: test_heterogeneous_poisson_params_t

  ! Functions
  public :: set_default_params, cli_add_params

contains

  subroutine set_default_params(params)
    implicit none
    type(test_heterogeneous_poisson_params_t), intent(inout) :: params

    ! IO parameters
    params%default_dir_path                 = 'data'
    params%default_prefix                   = 'UnitSquare_Heterogeneous' ! Name of GiD project
    params%default_dir_path_out             = 'output'
    
    ! Nonlinear solver parameters
    params%default_tolerance                = '1.0e-03'
    params%default_max_number_of_iterations = '20'

    ! Material parameters
    params%default_number_of_values         = '2'
    params%default_temperature_points       = '0.0 100.0'
    params%default_conductivity_values      = '1.0 2.0'
				
  end subroutine set_default_params
  !==================================================================================================

  subroutine cli_add_params(cli,params,group)
    implicit none
    type(Command_Line_Interface)             , intent(inout) :: cli
    type(test_heterogeneous_poisson_params_t), intent(in)    :: params
    character(*)                             , intent(in)    :: group
    ! Locals
    integer(ip) :: error

    ! Set Command Line Arguments
    ! IO parameters
    call cli%add(group=trim(group),switch='--dir_path',switch_ab='-d',                              &
         &       help='Directory of the source files',required=.false., act='store',                &
         &       def=trim(params%default_dir_path),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--prefix',switch_ab='-pr',help='Name of the GiD files',  &
         &       required=.false.,act='store',def=trim(params%default_prefix),error=error) 
    check(error==0)
    call cli%add(group=trim(group),switch='--dir_path_out',switch_ab='-out',help='Output Directory',&
         &       required=.false.,act='store',def=trim(params%default_dir_path_out),error=error)
    check(error==0)
    
    ! Nonlinear solver parameters
    call cli%add(group=trim(group),switch='--tolerance',switch_ab='-tol',                           &
         &       help='Residual tolerance',required=.false.,act='store',                            &
         &       def=trim(params%default_tolerance),error=error)
    call cli%add(group=trim(group),switch='--max_number_of_iterations',switch_ab='-maxiter',        &
         &       help='Maximum number of iterations',required=.false.,act='store',                  &
         &       def=trim(params%default_max_number_of_iterations),error=error)
    
    ! Material parameters
    call cli%add(group=trim(group),switch='--number_of_values',switch_ab='-nvalu',                  &
         &       help='Number of values',required=.false.,act='store',                              &
         &       def=trim(params%default_number_of_values),error=error)
    call cli%add(group=trim(group),switch='--temperature_points',switch_ab='-temps',                &
         &       help='Temperature points',required=.false.,act='store',nargs='*',                  &
         &       def=trim(params%default_temperature_points),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--conductivity_values',switch_ab='-conds',               &
         &       help='Conductivity values',required=.false.,act='store',nargs='*',                 &
         &       def=trim(params%default_conductivity_values),error=error)
    check(error==0)
				
  end subroutine cli_add_params
  !==================================================================================================

end module command_line_parameters_names

!****************************************************************************************************

module heterogeneous_poisson_analytical_functions_names
  use serial_names
  implicit none
# include "debug.i90"
  private
  
  type, extends(scalar_function_t) :: temperature_function_t
     
   contains
     procedure :: get_value_space => temperature_function_get_value_space
  end type temperature_function_t
  
  type heterogeneous_poisson_analytical_functions_t
     type(temperature_function_t) :: temperature_constraint
   contains

  end type heterogeneous_poisson_analytical_functions_t

  public :: heterogeneous_poisson_analytical_functions_t
  
contains  

  !===============================================================================================
  subroutine temperature_function_get_value_space ( this, point, result )
     implicit none
     class(temperature_function_t), intent(in) :: this
     type(point_t), intent(in)                 :: point
     real(rp), intent(inout)                   :: result
     
     !   0ºC on the lower edge
     if ( point%get(2) == 0.0_rp ) result = 0.0_rp
     
     ! 100ºC on the upper edge
     if ( point%get(2) == 1.0_rp ) result = 100.0_rp
     
  end subroutine temperature_function_get_value_space

end module heterogeneous_poisson_analytical_functions_names

!****************************************************************************************************

module property_table_names
  use types_names
  use memor_names
  
  implicit none
# include "debug.i90"
  private
  
  type :: property_table_t
     ! User-defined property tables
     integer(ip)                  :: nvalu
     real(rp)       , allocatable :: data_points(:)
     real(rp)       , allocatable :: property_values(:)
   contains
     !procedure, non_overridable   :: create
     procedure, non_overridable   :: interpolate
     procedure, non_overridable   :: search_interval
     !procedure, non_overridable   :: free
  end type property_table_t

  public :: property_table_t
  
contains  

  !===============================================================================================
  function interpolate ( this, point ) result ( value )
     implicit none
     class(property_table_t), intent(in) :: this
     real(rp)               , intent(in) :: point
     real(rp)                            :: value

     integer(ip)                         :: ileft, iright
     real(rp)                            :: slope
     
     ileft = 1; iright = this%nvalu
     if (point <= this%data_points(ileft)) then
        value = this%property_values(ileft)
     else if (point >= this%property_values(iright)) then
        value = this%property_values(iright)
     else
        call this%search_interval( point, ileft, iright )
     endif
     
     slope = (this%property_values(iright)-this%property_values(ileft))/ &
             &  (this%data_points(iright)-this%data_points(ileft))
     
     value = slope*(point-this%data_points(ileft)) + this%property_values(ileft)
     
  end function interpolate
  
  !===============================================================================================
  recursive subroutine search_interval ( this, point, ileft, iright )
     implicit none
     class(property_table_t), intent(in)    :: this
     real(rp)               , intent(in)    :: point
     integer(ip)            , intent(inout) :: ileft, iright

     integer(ip)                            :: imid
     
     if (iright - ileft > 1) then
        imid = (ileft + iright)/2
        if (this%data_points(imid) < point) then
           ileft = imid
           call this%search_interval( point, ileft, iright )
        else
           iright = imid
           call this%search_interval( point, ileft, iright )
        endif 
     endif
     
  end subroutine search_interval

end module property_table_names

!****************************************************************************************************

module heterogeneous_poisson_discrete_integration_names
  use serial_names
  use property_table_names
  
  implicit none
# include "debug.i90"
  private

  type, extends(discrete_integration_t) :: heterogeneous_poisson_discrete_integration_t
     class(fe_function_t)   , pointer :: fe_function        => NULL()
     class(property_table_t), pointer :: conductivity_table => NULL()  
   contains
     procedure                        :: integrate
  end type heterogeneous_poisson_discrete_integration_t
  
  public :: heterogeneous_poisson_discrete_integration_t

contains
  
  !===============================================================================================
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(heterogeneous_poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                           , intent(inout) :: fe_space
    class(assembler_t)                                 , intent(inout) :: assembler
    ! Locals
    type(finite_element_t), pointer    :: fe
    type(volume_integrator_t), pointer :: vol_int_first_fe, vol_int_second_fe
    real(rp), allocatable              :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer            :: fe_map
    type(quadrature_t), pointer        :: quad

    integer(ip)                        :: qpoint,inode,jnode,ioffset,joffset,ngaus
    real(rp)                           :: factor

    type(vector_field_t)               :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t)               :: grad_test_vector, grad_trial_vector

    integer(ip)                        :: number_fe_spaces

    integer(ip), pointer               :: field_blocks(:)
    logical, pointer                   :: field_coupling(:,:)

    integer(ip)                        :: ielem, number_nodes
    type(i1p_t), pointer               :: elem2dof(:)
    integer(ip), allocatable           :: number_nodes_per_field(:)
    
    type(fe_function_scalar_t)         :: scalar_unknown
    type(fe_function_vector_t)         :: vector_unknown
    
    real(rp)                           :: qpoint_value_scalar
    type(vector_field_t)               :: qpoint_value_vector
    
    real(rp)                           :: conductivity
    type(tensor_field_t)               :: conductivity_matrix
				
    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks    => fe_space%get_field_blocks()
    field_coupling  => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )
    
    call fe_space%initialize_integration()
    
    quad => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    
    call fe_space%create_fe_function(1,scalar_unknown)
    call fe_space%create_fe_function(2,vector_unknown)
    
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map            => fe%get_fe_map()
       vol_int_first_fe  => fe%get_volume_integrator(1)
       vol_int_second_fe => fe%get_volume_integrator(2)
       elem2dof          => fe%get_elem2dof()
						
       call fe%update_values(scalar_unknown, this%fe_function)
       call fe%update_values(vector_unknown, this%fe_function)
       
       do qpoint = 1,ngaus
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
          call scalar_unknown%get_value(qpoint, qpoint_value_scalar)
          conductivity = this%conductivity_table%interpolate(qpoint_value_scalar)
          
          do inode = 1, number_nodes_per_field(1)
             call vol_int_first_fe%get_gradient(inode,qpoint,grad_trial_scalar)
             do jnode = 1, number_nodes_per_field(1)
                call vol_int_first_fe%get_gradient(jnode,qpoint,grad_test_scalar)
                elmat(inode,jnode) = elmat(inode,jnode) + & 
                   & factor * conductivity * grad_test_scalar * grad_trial_scalar
             end do
          end do
          
          call vector_unknown%get_value(qpoint, qpoint_value_vector)
          call conductivity_matrix%init(0.0_rp)
          call conductivity_matrix%set(1,1,this%conductivity_table%interpolate(qpoint_value_vector%get(1)))
          call conductivity_matrix%set(2,2,this%conductivity_table%interpolate(qpoint_value_vector%get(2)))
          
          do inode = 1, number_nodes_per_field(2)
             ioffset = number_nodes_per_field(1) + inode
             call vol_int_second_fe%get_gradient(inode,qpoint,grad_trial_vector)
             do jnode = 1, number_nodes_per_field(2)
                joffset = number_nodes_per_field(1) + jnode
                call vol_int_second_fe%get_gradient(jnode,qpoint,grad_test_vector)
                elmat(ioffset,joffset) = elmat(ioffset,joffset) + & 
                   & factor * double_contract(grad_test_vector,conductivity_matrix*grad_trial_vector)
             end do
          end do
          
       end do
       
       ! Apply boundary conditions
       call fe%impose_strong_dirichlet_bcs( elmat, elvec )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )
    end do
    
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    
    call scalar_unknown%free()
    call vector_unknown%free()
 
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    
  end subroutine integrate
  
end module heterogeneous_poisson_discrete_integration_names

!****************************************************************************************************

program test_heterogeneous_poisson
  use serial_names
  use command_line_parameters_names
  use heterogeneous_poisson_discrete_integration_names
  use heterogeneous_poisson_analytical_functions_names
  use property_table_names
  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                         :: mesh
  type(conditions_t)                   :: scalar_conds
  type(conditions_t)                   :: composed_conds
  type(triangulation_t)                :: triangulation
		
  class(vector_t), allocatable, target :: residual
  class(vector_t), pointer             :: dof_values

  type(iterative_linear_solver_t)      :: iterative_linear_solver
  type(ParameterList_t)                :: parameter_list
  type(ParameterList_t), pointer       :: iterative_linear_solver_params
  type(serial_environment_t)           :: serial_env

  character(len=256)                   :: dir_path, dir_path_out
  character(len=256)                   :: prefix, filename

  logical                              :: diagonal_blocks_symmetric_storage(1)
  logical                              :: diagonal_blocks_symmetric(1)
  integer(ip)                          :: diagonal_blocks_sign(1)

  integer(ip)                          :: lunio, istat

  type(Command_Line_Interface)         :: cli 
  character(len=:), allocatable        :: group

  type(serial_fe_space_t)              :: fe_space
  type(fe_affine_operator_t)           :: fe_affine_operator
  type(p_reference_fe_t)               :: reference_fe_array(2)
  
  type(heterogeneous_poisson_discrete_integration_t) :: heterogeneous_poisson_integration
  type(heterogeneous_poisson_analytical_functions_t) :: analytical_functions
  type(fe_function_t), target          :: fe_function
  type(property_table_t), target       :: conductivity_table
  
  integer(ip)                          :: count, max_number_iterations
  real(rp)                             :: tolerance, residual_nrm2

  call meminit
  
  ! Read IO parameters
  call read_flap_cli_test_cdr(cli)
  call cli%parse(error=istat)
  if(cli%run_command('default')) then
     group = 'default'
  else
     group = 'default'
  end if
		
  call cli%get(group=trim(group),switch='-d',val=dir_path,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-pr',val=prefix,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-out',val=dir_path_out,error=istat); check(istat==0)

  ! Read mesh
  call mesh_read (dir_path, prefix, mesh, permute_c2z=.true.)

  ! Read scalar conditions 
  call conditions_read (dir_path, prefix, mesh%npoin, scalar_conds)
  
  ! Create conditions. Currently, they are generated in a rather dirty way,
  ! by replicating the scalar conditions to get the composite space conditions.
  ! However, this allows to test the code with composite finite element spaces.
  call conditions_create(3,3,scalar_conds%ncond,composed_conds)
  composed_conds%code(1,1:scalar_conds%ncond) = scalar_conds%code(1,1:scalar_conds%ncond)
  composed_conds%code(2,1:scalar_conds%ncond) = scalar_conds%code(1,1:scalar_conds%ncond)
  composed_conds%code(3,1:scalar_conds%ncond) = scalar_conds%code(1,1:scalar_conds%ncond)
  composed_conds%valu(1,1:scalar_conds%ncond) = scalar_conds%valu(1,1:scalar_conds%ncond)
  composed_conds%valu(2,1:scalar_conds%ncond) = scalar_conds%valu(1,1:scalar_conds%ncond)
  composed_conds%valu(3,1:scalar_conds%ncond) = scalar_conds%valu(1,1:scalar_conds%ncond)

  ! Read default loop parameters
  call cli%get(group=trim(group), &
               switch='-tol', &
               val=tolerance, &
               error=istat); check(istat==0)
  call cli%get(group=trim(group), &
               switch='-maxiter', &
               val=max_number_iterations, &
               error=istat); check(istat==0)
  
  ! Read temperature - thermal conductivity mapping
  call cli%get(group=trim(group), &
               switch='-nvalu', &
               val=conductivity_table%nvalu, &
               error=istat); check(istat==0)
  call cli%get_varying(group=trim(group), &
                       switch='-temps', &
                       val=conductivity_table%data_points, &
                       error=istat); check(istat==0)
  call cli%get_varying(group=trim(group), &
                       switch='-conds', &
                       val=conductivity_table%property_values, &
                       error=istat); check(istat==0) 
  
  heterogeneous_poisson_integration%conductivity_table => conductivity_table

  ! Construct triangulation
  call mesh_to_triangulation ( mesh, triangulation, gcond = composed_conds )
  
  ! Create and fill scalar x vector space
  reference_fe_array(1) =  make_reference_fe ( topology = topology_tet, &
                                               fe_type = fe_type_lagrangian, &
                                               number_dimensions = 2, &
                                               order = 1, &
                                               field_type = field_type_scalar, &
                                               continuity = .true. )

  reference_fe_array(2) =  make_reference_fe ( topology = topology_tet, &
                                               fe_type = fe_type_lagrangian, &
                                               number_dimensions = 2, &
                                               order = 1, &
                                               field_type = field_type_vector, &
                                               continuity = .true. )
     
  call fe_space%create( triangulation = triangulation, &
                        boundary_conditions = composed_conds, &
                        reference_fe_phy = reference_fe_array, &
                        field_blocks = (/1,2/), &
                        field_coupling = reshape((/.true.,.false.,.false.,.true./),(/2,2/)) )

  call fe_space%fill_dof_info() 
  
  ! Impose strong Dirichlet values
  call fe_space%update_bc_value(analytical_functions%temperature_constraint,1,1)
  call fe_space%update_bc_value(analytical_functions%temperature_constraint,1,2,unknown_component=1)
  call fe_space%update_bc_value(analytical_functions%temperature_constraint,1,2,unknown_component=2)
  
  ! Create affine operator
  call fe_affine_operator%create ( sparse_matrix_storage_format= 'CSR', &
                                   diagonal_blocks_symmetric_storage=(/.true./), &
                                   diagonal_blocks_symmetric=(/.true.,.true./), &
                                   diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                   environment=serial_env, &
                                   fe_space=fe_space, &
                                   discrete_integration=heterogeneous_poisson_integration )
  
  ! Create and initialize fe_function_t and associate to discrete integration
  call fe_space%create_global_fe_function(fe_function)
  dof_values => fe_function%get_dof_values()
  call dof_values%init(0.0_rp)
  heterogeneous_poisson_integration%fe_function => fe_function
  
  ! Fill affine operator
  call fe_affine_operator%symbolic_setup()
  call fe_affine_operator%numerical_setup()
  
  ! Initialize default solver parameters
  count = 1
  call fe_affine_operator%create_range_vector(residual)
  residual_nrm2 = 1.0_rp
                                               
  call the_iterative_linear_solver_creational_methods_dictionary%init()
 
  do  while (residual_nrm2 > tolerance .and. count <= max_number_iterations)

     ! Create iterative linear solver, set operators and solve linear system
     call iterative_linear_solver%create(serial_env)                                                    
     call iterative_linear_solver%set_type_from_string(cg_name)                
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%set_initial_solution(dof_values)
     call iterative_linear_solver%solve(fe_affine_operator%get_translation(), dof_values)          
     call iterative_linear_solver%free()                                                            
     
     ! Update affine operator with solution at current iteration     
     call fe_affine_operator%free_in_stages(free_numerical_setup)
     call fe_affine_operator%numerical_setup()
     
     ! Evaluate norm of the residual
     call fe_affine_operator%apply(dof_values,residual)
     residual_nrm2 = residual%nrm2()
     
     ! Print current norm of the residual
     write(*,*) 'Norm of the residual at iteration',count,'=',residual_nrm2
     
     ! Update iteration counter
     count = count + 1
     
  end do
  
  ! Write convergence report and solution
  write(*,*) 'Fixed-point algorithm converged in',count,'iterations'
  select type(dof_values)
     class is(serial_scalar_array_t)
     call dof_values%print(6)
     class is(serial_block_array_t)
     call dof_values%print(6)
     class default
     check(.false.) 
  end select
		
  ! Release memory
  call iterative_linear_solver%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call residual%free()
  call dof_values%free()
  call fe_function%free()
  call reference_fe_array(1)%free()
  call reference_fe_array(2)%free()
  call triangulation_free(triangulation)
  call conditions_free(scalar_conds)
  call conditions_free(composed_conds)
  call mesh_free(mesh)

  call memstatus 

contains
  
  !==================================================================================================
  subroutine read_flap_cli_test_cdr(cli)
    implicit none
    type(Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(test_heterogeneous_poisson_params_t) :: test_params
    logical                                   :: authors_print
    integer(ip)                               :: error

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'test_heterogeneous_poisson',                                       &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         & 'Serial FEMPAR driver for heterogeneous poisson problems using Continuous-Galerkin .',   &
         &        examples    = ['test_heterogeneous_poisson            -h ',                       &
                                 'test_heterogeneous_poisson default    -h ' ]) ! ALERT: spaces

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='default',description='Solve an heterogeneous poisson problem using CG')

    ! Set Command Line Arguments for each group
    call set_default_params(test_params)
    call cli_add_params(cli,test_params,'default')

  end subroutine read_flap_cli_test_cdr
		
end program test_heterogeneous_poisson
