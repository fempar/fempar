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

!* The `[[test_steady_poisson]]`,
program test_steady_poisson
!* uses the `fempar_names` and `test_steady_poisson_driver_names`:
  use fempar_names
  use poisson_analytical_functions_names
  use poisson_discrete_integration_names  
  !* First, declare the `test_driver` and the `world_context`.
  implicit none  
# include "debug.i90"
  type(serial_context_t)               :: world_context
  type(environment_t)                  :: serial_environment
  !* The parameter_handler is an object that provides default values for all the keys that FEMPAR uses to run. They can be:
  !*   1) Default values of the library.
  !*   2) Values provided by the user through the command line (using the keys in fempar_parameter_handler_t3).
  !*   3) overwrites 2) which overwrites 1). 
  !* In this tutorial we will explicitly provide the values in the code (option 3), but they could be provided by the command line argument instead.
  type(fempar_parameter_handler_t)     :: parameter_handler
  !* This is the object in parameter_handler_t that provides the list of parameters.
  type(ParameterList_t), pointer       :: parameter_list
  !* The triangulation_t object provides the mesh. In this case, we consider a serial triangulation, i.e., not partitioned.
  type(serial_triangulation_t)         :: triangulation
  !* This is a pointer to the finite element used to map the reference cell to the physical space. Usually, it is first order.
  !* The FE space used for this purpose will be extracted from the triangulation.
  class(reference_fe_t), pointer       :: reference_fe_geo
  !* This is a pointer to the reference finite element used to approximate the Poisson problem.
  type(p_reference_fe_t)               :: reference_fe
  type(p_reference_fe_t), allocatable  :: reference_fes(:) 
  !* The triangulation can provide a reference_fe... too complicated.
  !* SB: Provide the topology in a different way!!!!
  class(cell_iterator_t), allocatable  :: cell
  !* The fe_space_t is the global finite element space to be used.
  type(serial_fe_space_t)              :: fe_space
  !* It is an extension of conditions_t that defines the Dirichlet boundary conditions using analytical functions.
  type(strong_boundary_conditions_t)   :: poisson_conditions
  !* An analytical_function_t with the expression of the desired source term, which is just 0 in this case.
  type(source_term_t)                  :: source_term
  !* An analytical_function_t with the expression of the desired Dirichlet boundary condition.
  type(boundary_function_t)            :: boundary_function
  !* An analytical_function_t with the expression of exact solution of the problem, to check the code.
  type(solution_function_t)            :: exact_solution
  !* A fe_function_t belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                  :: solution
  !* poisson_discrete_integration_t provides the definition of the blilinear form and right-hand side of the problem at hand.
  type(poisson_discrete_integration_t) :: poisson_integration
  !* A fe_affine_operator_t that represents the affine operator the solution of which is the one we want, i.e., B = Ax-f.
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)           :: fe_affine_operator
  !* The problem ca be solved using a direct solver (Pardiso) or an iterative linear solver (fempar cg). The solvers and their 
  !* parameters are defined later.
#ifdef ENABLE_MKL    
  !*  - Direct solver
  type(direct_solver_t)                :: direct_solver
  type(ParameterList_t), pointer       :: direct_solver_params
  integer                              :: iparm(64) = 0  
#else 
  !*  - Iterative solver
  type(iterative_linear_solver_t)      :: iterative_linear_solver
  type(ParameterList_t), pointer       :: iterative_solver_params
#endif    
  !* The output handler type is used to print the results
  type(output_handler_t) :: output_handler
  !* The following object automatically compute error norms given a fe_function_t and the analytical solution.
  type(error_norms_scalar_t) :: error_norm
  !*
  !* Local variables
  integer(ip) :: fe_order, istat, error, i, boundary_ids
  class(vector_t), pointer       :: dof_values
  class(vector_t), allocatable :: rhs
  real(rp) :: l2
  !*
  !* Initialize properly the FEMPAR library 
  call fempar_init()
  !* Initialize the FEMPAR context 
  call world_context%create()
  !* Initialize the list of parameters with all the options that are provided by FEMPAR
  !* It involves to create a parameter handler that is usually used to extract the values
  !* provided by the user through the command line. In this test, we assume that we are 
  !* not going to make use of the command line, and we are going to set the desired values
  !* in the driver instead.
  call parameter_handler%process_parameters()
  parameter_list => parameter_handler%get_values()
  !* Set parameters for the integration
  error = 0
  error = error + parameter_list%set(key = struct_hex_triang_num_dims_key, value = 2 )              ! Number of space dimensions
  error = error + parameter_list%set(key = struct_hex_triang_num_cells_dir, value = [10,10,10] ) ! Number of cells per each dimension
  error = error + parameter_list%set(key = fes_ref_fe_orders_key, value = [ 1 ] )                   ! Reference finite element orders
  error = error + parameter_list%set(key = fes_ref_fe_types_key, value = 'Lagrangian' )             ! Reference finite element types
  mcheck(error==0,'Failed parameter set')
  !* Print the list of parameters
  !    call parameter_list%print()
  !* Determine a serial execution mode (default case)
  call serial_environment%create(world_context, parameter_list)
  !* Create triangulation.  
  !* The structure mesh generator creates a domain that has:
  !*   8 boundary objects (corners+edges) in 2D
  !*   26 boundary objects (face+corners+edges) in 3D
  !* and 1 interior object 
  !*   1 face in 2D 
  !*   1 cell in 3D
  !* A graphical example for a 2D case is presented next
  !*           6
  !*      4 - - - - 3
  !*      -         -         1,2,3,4,5,6,7,8 are id for the boundary of the domain
  !*    7 -    9    - 8
  !*      -         -         9 is the id for the interior
  !*      1 - - - - 2
  !*           5     
  call triangulation%create(serial_environment, parameter_list)
  !* Set the boundary Dirichlet data, using a user-defined analytical function (see below) with the expression we want.
  call boundary_function%set_num_dims(triangulation%get_num_dims())
  boundary_ids = 8 
  if ( triangulation%get_num_dims() == 3) boundary_ids = 26
  call poisson_conditions%create()
  do i = 1, boundary_ids
    call poisson_conditions%insert_boundary_condition(boundary_id=i, field_id=1, cond_type=component_1, boundary_function=boundary_function)
  end do
  !* Next, we build the global FE space. It only requires to know:
  !*   The triangulation
  !*   The Dirichlet data
  !*   The reference FE to be used.
  call fe_space%create( triangulation            = triangulation,      &
                        conditions               = poisson_conditions, &
                        parameters               = parameter_list )
  ! We must explicitly say that we want to use integration arrays, e.g., quadratures, maps, etc. 
  call fe_space%set_up_cell_integration()
  ! Now, we define the source term with the function we have created in our module.
  call source_term%set_num_dims(triangulation%get_num_dims())
  call poisson_integration%set_source_term(source_term)
  !* Now, we create the affine operator, i.e., b - Ax, providing the info for the matrix (storage, symmetric, etc.), and the form to
  !* be used to fill it, e.g., the bilinear form related that represents the weak form of the Poisson problem and the right hand
  !* side.  
  call fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                   diagonal_blocks_symmetric_storage = [ .true. ], &
                                   diagonal_blocks_symmetric         = [ .true. ], &
                                   diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                   fe_space                          = fe_space, &
                                   discrete_integration              = poisson_integration )
  !* In the next lines, we create a FE function of the FE space defined above, load the Dirichlet values defined above, and put it
  !* in the discrete integration (note that for linear operators, the Dirichlet values are used to compute the RHS).
  call solution%create(fe_space) 
  call fe_space%interpolate_dirichlet_values(solution)
  call poisson_integration%set_fe_function(solution)
  !* Now, the discrete integration has all the information needed. We can fill the affine operator b - Ax.
  call fe_affine_operator%compute()
#ifdef ENABLE_MKL
  !* PARDISO MKL
  direct_solver_params => parameter_list%NewSubList(Key=pardiso_mkl)
  !* Set parameters of the direct solver
  error = 0
  error = error + direct_solver_params%set(key = dls_type_key,        value = pardiso_mkl)
  error = error + direct_solver_params%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_uss)
  error = error + direct_solver_params%set(key = pardiso_mkl_message_level, value = 0)
  error = error + direct_solver_params%set(key = pardiso_mkl_iparm,         value = iparm)
  check(error == 0)
  call direct_solver%set_type_from_pl(direct_solver_params)
  !* Next, we set the matrix in our system from the fe_affine operator
  call direct_solver%set_matrix(fe_affine_operator%get_matrix())
  !* We extract a pointer to the free nodal values of our FE function, which is the plae in which we will store the result.
  dof_values => solution%get_free_dof_values()
  !* We solve the problem with the matrix already associated, the RHS obtained from the fe_affine_operator using get_translation, and 
  !* putting the result in dof_values.
  call direct_solver%solve(fe_affine_operator%get_translation(),dof_values)
  call direct_solver%log_info()
#else   
  !* ITERATIVE SOLVER
  iterative_solver_params => parameter_list%NewSubList(Key=cg_solver)
  !* Set parameters of the direct solver
  error = 0
  error = error + iterative_solver_params%set(key = ils_rtol_key, value = 1.0e-12_rp)         ! Relative tolerance
  error = error + iterative_solver_params%set(key = ils_max_num_iterations_key, value = 5000) ! Max number of iterations
  error = error + iterative_solver_params%set(key = ils_type_key, value = cg_name )           ! Conjugate Gradient iterative solver
  assert(error == 0)
  !* Now, we create a serial iterative solver with the values in the parameter list.
  call iterative_linear_solver%create(fe_space%get_environment())
  call iterative_linear_solver%set_type_and_parameters_from_pl(iterative_solver_params)
  !* Next, we set the matrix in our system from the fe_affine operator (i.e., its tangent)
  call iterative_linear_solver%set_operators(fe_affine_operator%get_tangent(), .identity. fe_affine_operator) 
  !* We extract a pointer to the free nodal values of our FE function, which is the plae in which we will store the result.
  dof_values => solution%get_free_dof_values()
  !* We solve the problem with the matrix already associated, the RHS obtained from the fe_affine_operator using get_translation, and 
  !* putting the result in dof_values.
  call iterative_linear_solver%apply(fe_affine_operator%get_translation(),dof_values)   
#endif  
  !* Now, we want to compute error wrt an exact solution. It needs the FE space above to be constructed.                                       
  call error_norm%create(fe_space,1)  
  call exact_solution%set_num_dims(triangulation%get_num_dims())
  !* We compute the L2 norm of the difference between the exact and computed solution.
  l2 = error_norm%compute(exact_solution, solution, l2_norm) 
  !* We finally write the result and check we have solved the problem "exactly", since the exact solution belongs to the FE space.
  write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
  !* Last, we plot the results.
  !* We can set the direction where the results files are created using
  !*   error = 0
  !*   error = error + direct_solver_params%set(key = dir_path_out_key, value = '.')
  !*   check(error == 0)
  !* But in this example we use the default path.
  error = error + parameter_list%set(key = dir_path_out_key      , Value= 'RESULTS_PROBLEM')
  call output_handler%create()
  call output_handler%attach_fe_space(fe_space)
  call output_handler%add_fe_function(solution, 1, 'phi')
  call output_handler%open(parameter_handler%get_dir_path_out(), parameter_handler%get_prefix())
  call output_handler%write()
  call output_handler%close()
  call output_handler%free()
!*
!* Free all the created objects
!*  call fempar_parameter_handler%free()
!*  call triangulation%free()
  call error_norm%free()
  call solution%free()
#ifdef ENABLE_MKL        
    call direct_solver%free()
#else
    call iterative_linear_solver%free()
#endif
  call fe_affine_operator%free()
  call fe_space%free()
  call poisson_conditions%free()
  call triangulation%free()
  call serial_environment%free()
  call parameter_handler%free()
  call world_context%free(.true.)

  call fempar_finalize()
end program test_steady_poisson


















