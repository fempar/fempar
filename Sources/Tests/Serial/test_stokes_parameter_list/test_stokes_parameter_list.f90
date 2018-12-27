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

!* The `[[test_stokes_parameter_list]]`,
program test_stokes_parameter_list
!* uses the `fempar_names` and `test_stokes_parameter_list_driver_names`:
  use fempar_names
  use stokes_discrete_integration_names  
  !* First, declare the `test_driver` and the `world_context`  
  implicit none  
# include "debug.i90"
  type(serial_context_t)               :: world_context
  type(environment_t)                  :: serial_environment
  !* The parameter_handler is an object that provides default values for all the keys that FEMPAR uses to run. They can be be
  !* 1) default values of the library, 2) the ones provided by the user through the command line (using the keys in 
  !* fempar_parameter_handler_t, 3) or overwritten by the user. Option 3) overwrites 2) which overwrites 1). In this tutorial
  !* we will explicitly provide the values in the code (option 3) but they could be provided by the command line argument instead.
  type(fempar_parameter_handler_t)     :: parameter_handler
  !* This is the object in parameter_handler_t that provides the list of parameters
  type(ParameterList_t), pointer       :: parameter_list
  !* The triangulation_t object provides the mesh. In this case, we consider a serial triangulation, i.e., not partitioned.
  type(serial_triangulation_t)         :: triangulation
  !* This is a pointer to the finite element used to map the reference cell to the physical space. Usually, it is first order.
  !* The FE space used for this purpose will be extracted from the triangulation.
  class(reference_fe_t), pointer       :: reference_fe_geo
  !* This is a pointer to the reference finite element used to approximate the stokes problem.
  type(p_reference_fe_t)               :: reference_fe
  type(p_reference_fe_t), allocatable  :: reference_fes(:) 
  !* The triangulation can provide a reference_fe... too complicated
  !* SB: Provide the topology in a different way!!!!
  class(cell_iterator_t), allocatable  :: cell
  !* The fe_space_t is the global finite element space to be used.
  type(serial_fe_space_t)              :: fe_space
  !* It is an extension of conditions_t that defines the Dirichlet boundary conditions using analytical functions.
  type(strong_boundary_conditions_t)   :: stokes_conditions
  !* An analytical_function_t with the expression of the desired source term, which is just 0 in this case.
  type(vector_field_t) :: zero_vector
  type(constant_vector_function_t)            :: source_term
  !type(source_term_t)                  :: source_term
  !* An analytical_function_t with the expression of the desired Dirichlet boundary condition.
  type(constant_scalar_function_t)            :: zero_function
  type(constant_scalar_function_t)            :: one_function
  !* An analytical_function_t with the expression of exact solution of the problem, to check the code.
  !type(solution_function_t)            :: exact_solution
  !* A fe_function_t belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                  :: solution
  !* stokes_discrete_integration_t provides the definition of the blilinear form and right-hand side of the problem at hand.
  type(stokes_discrete_integration_t) :: stokes_integration
  !* A fe_affine_operator_t that represents the affine operator the solution of which is the one we want, i.e., B = Ax-f.
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)           :: fe_affine_operator
  !* The problem will be solved with an iterative linear solver, to be defined later.
  type(iterative_linear_solver_t)      :: iterative_linear_solver
  !* The following object automatically compute error norms given a fe_function_t and the analytical solution.
  !type(error_norms_scalar_t) :: error_norm
  !* The output handler type is used to print the results
  type(output_handler_t) :: output_handler

  !* Local variables
  integer(ip) :: fe_order, istat, error, i, boundary_ids
  class(vector_t), pointer       :: dof_values
  class(vector_t), allocatable :: rhs
  real(rp) :: l2
  
  !* Initialize properly the FEMPAR library 
  call fempar_init()
  !* Initialize the FEMPAR context 
  call world_context%create()
  !* Initialize the list of parameters with all the options that are provided by FEMPAR
  !* It involves to create a parameter handler that is usually used to extract the values
  !* provided by the user through the command line. In this test, we assume that we are 
  !* not going to make use of the command line, and we are going to set the desired values
  !* in the driver instead.
  call parameter_handler%process_parameters(test_stokes_parameter_list_define_fempar_parameters)
  parameter_list => parameter_handler%get_values()

  !* Determine a serial execution mode (default case)
  call serial_environment%create(world_context,parameter_list)

  !* Create triangulation (we are not changing the default parameters from FEMPAR. Thus, we are solving in a serial
  !* triangulation of hex, structured, 10x10, 2D.
  call triangulation%create(serial_environment,parameter_list)

  !* Set the boundary Dirichlet data, using a user-defined analytical function (see below) with the expression we want.
  call zero_function%create(0.0_rp)
  call one_function%create(1.0_rp)
  !call old_stokes_conditions%set_boundary_function(boundary_function)

  !* Next, we build the global FE space. It only requires to know the triangulation, the Dirichlet data, and the reference FE to be
  !* used.
  ! The structured mesh generator has 8 (corner+edge) boundary objects (w/ different set id) in the triangulation, following the
  ! same numbering as the 2-cube in the FEMPAR article. Analogously, 26 (corner+edge+face) boundary objects in 3D.
  boundary_ids = 8
  if ( triangulation%get_num_dims() == 3) boundary_ids = 26
  call stokes_conditions%create()
  do i = 1, boundary_ids
    call stokes_conditions%insert_boundary_condition(boundary_id=i, field_id=1, &
                                                     cond_type=component_2, boundary_function=zero_function)
  end do
  ! 
  call stokes_conditions%insert_boundary_condition(boundary_id=1, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=2, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=5, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=8, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=7, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=3, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=4, field_id=1, &
                                                   cond_type=component_1, boundary_function=zero_function)
  call stokes_conditions%insert_boundary_condition(boundary_id=6, field_id=1, &
                                                   cond_type=component_1, boundary_function=one_function)
  error = 0
  error = error + parameter_list%set(key = fes_num_fields_key, value = 2)
  error = error + parameter_list%set(key = fes_num_ref_fes_key, value = 2)
  error = error + parameter_list%set(key = fes_ref_fe_types_key, value = 'Lagrangian Lagrangian' )
  error = error + parameter_list%set(key = fes_ref_fe_orders_key, value = [2, 1])
  error = error + parameter_list%set(key = fes_ref_fe_conformities_key, value =  [.true., .true.])
  error = error + parameter_list%set(key = fes_ref_fe_continuities_key, value = [.true., .false.])
  error = error + parameter_list%set(key = fes_field_types_key, value =  field_type_vector // " " // field_type_scalar)
  error = error + parameter_list%set(key = fes_field_blocks_key, value = [1, 1])
  mcheck(error==0,'Failed parameter set')
  
  call parameter_list%print()
  
  !
  call fe_space%create( triangulation            = triangulation,      &
                        conditions               = stokes_conditions, &
                        parameters               = parameter_list )
  ! We must explicitly say that we want to use integration arrays, e.g., quadratures, maps, etc. 
  call fe_space%set_up_cell_integration()

  ! Now, we define the source term with the function we have created in our module.
  !call source_term%set_num_dims(triangulation%get_num_dims())
  call zero_vector%init(0.0)
  call source_term%create(zero_vector)
  call stokes_integration%set_source_term(source_term)
  call stokes_integration%set_viscosity(1.0_rp)

  !* Now, we create the affine operator, i.e., b - Ax, providing the info for the matrix (storage, symmetric, etc.), and the form to
  !* be used to fill it, e.g., the bilinear form related that represents the weak form of the stokes problem and the right hand
  !* side.  
  call fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                   diagonal_blocks_symmetric_storage = [ .false. ], &
                                   diagonal_blocks_symmetric         = [ .false. ], &
                                   diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_INDEFINITE ], &
                                   fe_space                          = fe_space, &
                                   discrete_integration              = stokes_integration )
  
  !* In the next lines, we create a FE function of the FE space defined above, load the Dirichlet values defined above, and put it
  !* in the discrete integration (note that for linear operators, the Dirichlet values are used to compute the RHS).
  call solution%create(fe_space) 
  call fe_space%interpolate_dirichlet_values(solution)
  call stokes_integration%set_fe_function(solution)

  !* Now, the discrete integration has all the information needed. We can fill the affine operator b - Ax.
  call fe_affine_operator%compute()

  !* Next, we want to get the root of the operator, i.e., solve Ax = b. We are going to use an iterative solver. We overwrite the
  !* default values in the parameter list for tolerance, type of sover (Conjugate Gradient) and max number of iterations.
  error = 0
  error = error + parameter_list%set(key = ils_rtol_key, value = 1.0e-12_rp)
  error = error + parameter_list%set(key = ils_max_num_iterations_key, value = 5000)
  error = error + parameter_list%set(key = ils_type_key, value = fgmres_name )
  assert(error == 0)
  
  !* Now, we create a serial iterative solver with the values in the parameter list.
  call iterative_linear_solver%create(fe_space%get_environment())
  call iterative_linear_solver%set_type_and_parameters_from_pl(parameter_list)

  !* Next, we set the matrix in our system from the fe_affine operator (i.e., its tangent)
  call iterative_linear_solver%set_operators(fe_affine_operator%get_tangent(), .identity. fe_affine_operator) 
  !* We extract a pointer to the free nodal values of our FE function, which is the plae in which we will store the result.
  dof_values => solution%get_free_dof_values()
  !* We solve the problem with the matrix already associated, the RHS obtained from the fe_affine_operator using get_translation, and 
  !* putting the result in dof_values.
  call iterative_linear_solver%apply(-fe_affine_operator%get_translation(), &
                                          dof_values)   
  
  error = error + parameter_list%set(key = dir_path_out_key      , Value= 'RESULTS_STOKES_DRIVER')
  call output_handler%create()
  call output_handler%attach_fe_space(fe_space)
  call output_handler%add_fe_function(solution, 1, 'velocity')
  call output_handler%add_fe_function(solution, 2, 'pressure')
  call output_handler%open(parameter_handler%get_dir_path_out(), parameter_handler%get_prefix())
  call output_handler%write()
  call output_handler%close()
  call output_handler%free()
  
  
  
  !* Now, we want to compute error wrt an exact solution. It needs the FE space above to be constructed.                                       
  !call error_norm%create(fe_space,1)  
  !call exact_solution%set_num_dims(triangulation%get_num_dims())
  !* We compute the L2 norm of the difference between the exact and computed solution.
  !l2 = error_norm%compute(exact_solution, solution, l2_norm) 

  !* We finally write the result and check we have solved the problem "exactly", since the exact solution belongs to the FE space.
  !write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
!*
!* Free all the created objects
!*  call fempar_parameter_handler%free()
!*  call triangulation%free()
!*
!*   In the following code, we will use
  !call error_norm%free()
  call solution%free()
  call iterative_linear_solver%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call stokes_conditions%free()
  call triangulation%free()
  call serial_environment%free()
  call parameter_handler%free()
  call world_context%free(.true.)

  call fempar_finalize()
  contains
    subroutine test_stokes_parameter_list_define_fempar_parameters(this)
    implicit none
    class(fempar_parameter_handler_t), intent(inout) :: this
    type(parameterlist_t), pointer :: values, switches, switches_ab, helpers, required 
    integer(ip) :: error

    values      => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()

    error = 0
    error = error + helpers%set(key = 'viscosity'     , Value= 'Value of the viscosity')
    error = error + switches%set(key = 'viscosity'    , Value= '--VISCOSITY')
    error = error + values%set(key = 'viscosity'      , Value= '1.0')
    end subroutine  test_stokes_parameter_list_define_fempar_parameters
end program test_stokes_parameter_list


















