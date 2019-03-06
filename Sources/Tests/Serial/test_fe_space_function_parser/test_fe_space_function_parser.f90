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
  !* First, declare the `test_driver` and the `world_context`  
  implicit none  
# include "debug.i90"

  type(scalar_function_parser_t)              :: source_term_parser
  type(scalar_function_parser_t)              :: boundary_function_parser
  type(scalar_function_parser_t)              :: solution_gradient_parser_comp1
  type(scalar_function_parser_t)              :: solution_gradient_parser_comp2
  type(scalar_function_parser_t)              :: solution_gradient_parser_comp3
  type(vector_function_parser_t)              :: solution_gradient_parser
  type(scalar_function_parser_t)              :: solution_function_parser
  type(scalar_function_and_gradient_parser_t) :: solution_parser

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
  !* This is a pointer to the reference finite element used to approximate the Poisson problem.
  type(p_reference_fe_t)               :: reference_fe
  type(p_reference_fe_t), allocatable  :: reference_fes(:) 
  !* The triangulation can provide a reference_fe... too complicated
  !* SB: Provide the topology in a different way!!!!
  class(cell_iterator_t), allocatable  :: cell
  !* The fe_space_t is the global finite element space to be used.
  type(serial_fe_space_t)              :: fe_space
  !* It is an extension of conditions_t that defines the Dirichlet boundary conditions using analytical functions.
  type(strong_boundary_conditions_t)   :: poisson_conditions
  !* A fe_function_t belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                  :: solution
  !* poisson_discrete_integration_t provides the definition of the blilinear form and right-hand side of the problem at hand.
  type(poisson_discrete_integration_t) :: poisson_integration
  !* A fe_affine_operator_t that represents the affine operator the solution of which is the one we want, i.e., B = Ax-f.
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)           :: fe_affine_operator
  !* The problem will be solved with an iterative linear solver, to be defined later.
  type(iterative_linear_solver_t)      :: iterative_linear_solver
  !* The following object automatically compute error norms given a fe_function_t and the analytical solution.
  type(error_norms_scalar_t) :: error_norm

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
  call parameter_handler%process_parameters()
  parameter_list => parameter_handler%get_values()

  !* Determine a serial execution mode (default case)
  call serial_environment%create(world_context, parameter_list)

  !* Create triangulation (we are not changing the default parameters from FEMPAR. Thus, we are solving in a serial
  !* triangulation of hex, structured, 10x10, 2D.
  call triangulation%create(serial_environment, parameter_list)


  if ( triangulation%get_num_dims() == 2) then  
    call source_term_parser%create(expression="-( (6*x*y^3) + (6*x^3*y) )", num_dims=2) 
    call boundary_function_parser%create(expression="x^3*y^3", num_dims=2)
    call solution_function_parser%create(expression="x^3*y^3", num_dims=2)
    call solution_gradient_parser_comp1%create(expression="3*x^2*y^3", num_dims=2)
    call solution_gradient_parser_comp2%create(expression="3*x^3*y^2", num_dims=2)
    call solution_gradient_parser%create( solution_gradient_parser_comp1, &
                                          solution_gradient_parser_comp2)
  else if ( triangulation%get_num_dims() == 3 ) then ! 
    call source_term_parser%create(expression="-( (6*x*y^3*z^3) + (6*x^3*y*z^3) + (6*x^3*y^3*z) )", num_dims=3)
    call boundary_function_parser%create(expression="x^3*y^3*z^3", num_dims=3) 
    call solution_function_parser%create(expression="x^3*y^3*z^3", num_dims=3) 
    call solution_gradient_parser_comp1%create(expression="3*x^2*y^3*z^3", num_dims=3)
    call solution_gradient_parser_comp2%create(expression="3*x^3*y^2*z^3", num_dims=3)
    call solution_gradient_parser_comp3%create(expression="3*x^3*y^3*z^2", num_dims=3)
    call solution_gradient_parser%create( solution_gradient_parser_comp1, &
                                          solution_gradient_parser_comp2, &
                                          solution_gradient_parser_comp3)
  end if
  call solution_parser%create(solution_function_parser, solution_gradient_parser)

  !* Set the boundary Dirichlet data, using a user-defined analytical function (see below) with the expression we want.
  call boundary_function_parser%set_num_dims(triangulation%get_num_dims())
  !call old_poisson_conditions%set_boundary_function(boundary_function)

  !* Next, we build the global FE space. It only requires to know the triangulation, the Dirichlet data, and the reference FE to be
  !* used.
  ! The structured mesh generator has 8 (corner+edge) boundary objects (w/ different set id) in the triangulation, following the
  ! same numbering as the 2-cube in the FEMPAR article. Analogously, 26 (corner+edge+face) boundary objects in 3D.
  boundary_ids = 8
  if ( triangulation%get_num_dims() == 3) boundary_ids = 26
  call poisson_conditions%create()
  do i = 1, boundary_ids
    call poisson_conditions%insert_boundary_condition(boundary_id=i, field_id=1, cond_type=component_1, boundary_function=boundary_function_parser)
  end do
  call fe_space%create( triangulation            = triangulation,      &
                        conditions               = poisson_conditions, &
                        parameters               = parameter_list )
  ! We must explicitly say that we want to use integration arrays, e.g., quadratures, maps, etc. 
  call fe_space%set_up_cell_integration()

  ! Now, we define the source term with the function we have created in our module.
  call source_term_parser%set_num_dims(triangulation%get_num_dims())
  call poisson_integration%set_source_term(source_term_parser)

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

  !* Next, we want to get the root of the operator, i.e., solve Ax = b. We are going to use an iterative solver. We overwrite the
  !* default values in the parameter list for tolerance, type of sover (Conjugate Gradient) and max number of iterations.
  error = 0
  error = error + parameter_list%set(key = ils_rtol_key, value = 1.0e-12_rp)
  error = error + parameter_list%set(key = ils_max_num_iterations_key, value = 5000)
  error = error + parameter_list%set(key = ils_type_key, value = cg_name )
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
  call iterative_linear_solver%apply(fe_affine_operator%get_translation(), &
                                          dof_values)   
  !* Now, we want to compute error wrt an exact solution. It needs the FE space above to be constructed.                                       
  call error_norm%create(fe_space,1)  
  
  call solution_parser%set_num_dims(triangulation%get_num_dims())
  
  !* We compute the L2 norm of the difference between the exact and computed solution.
  l2 = error_norm%compute(solution_parser, solution, l2_norm) 

  !* We finally write the result and check we have solved the problem "exactly", since the exact solution belongs to the FE space.
  write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
!*
!* Free all the created objects
!*  call fempar_parameter_handler%free()
!*  call triangulation%free()
!*
!*   In the following code, we will use
  call error_norm%free()
  call solution%free()
  call iterative_linear_solver%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call poisson_conditions%free()
  call triangulation%free()
  call serial_environment%free()
  call parameter_handler%free()
  call world_context%free(.true.)

  call fempar_finalize()

end program test_steady_poisson

