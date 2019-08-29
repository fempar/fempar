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

!* ### Module loading
!* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr]],
program tutorial_03_poisson_sharp_circular_wave_parallel_amr

!* uses the [[fempar_names]] to import all **`FEMPAR`** library symbols (i.e., derived types, parameter constants, system-wide variables, etc.),
!* and a set of tutorial-specific module units, 
!*
!* - [[tutorial_03_functions_names]]
!* - [[tutorial_03_error_estimator_names]]
!* - [[tutorial_03_discrete_integration_names]]
!*
!* These modules are not part of the **`FEMPAR`** library, but developed specifically for the problem at hand. Each of these modules defines a tutorial-specific data type and its TBPs.
  use fempar_names
  use tutorial_03_functions_names
  use tutorial_03_error_estimator_names
  use tutorial_03_discrete_integration_names

!* ### Variable definition
!* This tutorial declares:
!*
!* - a set of parameter constants, typically the tutorial name, authors, problem description, etc. (to be output on screen on demand by the user),
!* - the tutorial data type instances in charge of the FE simulation, such as the triangulation (mesh) of the computational domain or the FE space from which the approximate solution of the PDE is sought, 
!* - and a set of variables to hold the values of the Command-Line-Arguments (CLAs) which are specific to the tutorial.
  implicit none
#include "debug.i90"
  
  character(*), parameter :: tutorial_name = "tutorial_03_poisson_sharp_circular_wave_parallel_amr"
  character(*), parameter :: tutorial_version = "v1.0.0"
  character(*), parameter :: tutorial_description = "MPI-parallel program that solves a 2D/3D Poisson problem on the unit &
                                                     square/cube, where the solution has a sharp circular & 
                                                     wave front parametrizable via command-line arguments. &
                                                     The discretization of the domain is performed using an MPI-parallel p4est &
                                                     triangulation, i.e., a mesh which can by h-adapted and re-partitioned &
                                                     (re-distributed) among the MPI tasks during the simulation. For the linear &
                                                     solver step, it uses a Krylov subspace solver preconditioned with a highly &
                                                     scalable implementation of the 2-level BDDC preconditioner available at FEMPAR."
  character(*), parameter :: tutorial_authors = "Santiago Badia and Alberto F. Martín"

  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:world_context]] is a software abstraction 
  !* for a group of parallel tasks (a.k.a. processes) and the communication layer that orchestrates their concurrent execution.
  !* In this case is declared of type [[mpi_context_t]] as it is designed to work in serial computing environments.
  type(mpi_context_t)                         :: world_context
  !* the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:environment]], of type [[environment_t]], organizes the tasks of the context
  !* from which it is set up into subgroups of tasks, referred to as levels, and builds up additional communication mechanisms 
  !* to communicate tasks belonging to different levels.
  type(environment_t)                         :: environment
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:triangulation]] object, of type [[p4est_par_triangulation_t]], provides the mesh. 
  !* In this case, we consider a parallel p4est triangulation, i.e., 
  !* a triangulation that is distributed among processors, and that can be h-adapted during the simulation.
  type(p4est_par_triangulation_t)             :: triangulation
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:coarse_fe_handlers]] array, of type [[p_l1_coarse_fe_handler_t]], 
  !* holds polymorphic pointers to data type extensions of [[coarse_fe_handler_t]], as many as fields in the system of PDEs at hand.
  type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:coarse_fe_handler]] object, of type [[h_adaptive_algebraic_l1_coarse_fe_handler_t]], 
  !* a type extension of [[coarse_fe_handler_t]] suitable for grad-conforming FE spaces on h-adaptive meshes
  type(h_adaptive_algebraic_l1_coarse_fe_handler_t), target :: coarse_fe_handler
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:fe_space]], of type [[par_fe_space_t]] is the global finite element space to be used.
  type(par_fe_space_t)                        :: fe_space
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:strong_boundary_conditions]], of type [[strong_boundary_conditions_t]], is the object that provides the strong boundary conditions definition.
  type(strong_boundary_conditions_t)          :: strong_boundary_conditions 
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:source_term]], of type [[tutorial_03_functions_names:sharp_circular_wave_source_term_t]], is a scalar-valued function with the right hand side of the PDE problem at hand.
  type(sharp_circular_wave_source_term_t)     :: source_term
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:exact_solution]], of type [[tutorial_03_functions_names:sharp_circular_wave_solution_t]], is a scalar-valued function with the exact (analytical) solution of the PDE problem at hand.
  type(sharp_circular_wave_solution_t)        :: exact_solution
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:fe_affine_operator]], of type [[fe_affine_operator_t]],
  !* represents the affine operator the solution of which is the one we want, i.e., \(B = Ax-f\).
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)                  :: fe_affine_operator
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:discrete_solution]], of type [[fe_function_t]], 
  !* is the object belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                         :: discrete_solution
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:cg_discrete_integration]], 
  !* of type [[tutorial_03_discrete_integration_names:cg_discrete_integration_t]],
  !* provides the definition of the bilinear and linear forms for the problem at hand.
  type(cg_discrete_integration_t)             :: cg_discrete_integration
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:iterative_linear_solver]], of type [[iterative_linear_solver_t]], 
  !* provides native implementations of a bunch of preconditioner Krylov subspace iterative solvers
  type(iterative_linear_solver_t)             :: iterative_linear_solver
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc_parameters]], 
  !* of type [parameterlist_t](http://victorsndvg.github.io/FPL/type/parameterlist_t.html), 
  !* to store MLBDDC related parameters and configure the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc]] instance.
  type(parameterlist_t)                       :: mlbddc_parameters
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc]] object, of type [[mlbddc_t]], 
  !* is in charge of building and applying the BDDC pre-conditioner at each iterative solver iteration
  type(mlbddc_t)                              :: mlbddc
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler]] object, of type [[output_handler_t]], is used to generate the simulation data files for later visualization using, e.g., ParaView.
  type(output_handler_t)                      :: output_handler
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:error_estimator]] object, of type [[tutorial_01_error_estimator_names:poisson_error_estimator_t]], is used to compute the (square) of the error energy norm for all cells, i.e., \(||u-u_h||_E,K\).
  type(poisson_error_estimator_t)             :: error_estimator
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:refinement_strategy]] object, of type [[fixed_fraction_refinement_strategy_t]], 
  !* implements the strategy used for mesh refinement and coarsening.
  type(fixed_fraction_refinement_strategy_t)  :: refinement_strategy
  !* The [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:my_rank_cell_array]] object, of type [[std_vector_real_rp_t]], 
  !* is a dynamic array which is adapted along with the triangulation
  !* at each AMR loop iteration. It has as many entries as local cells in each parallel task, and
  !* for all of these entries, it holds the parallel task identifier. It is written into output data files,
  !* so that the user may visualize how the adaptive mesh is distributed among the parallel tasks.
  type(std_vector_real_rp_t)                  :: my_rank_cell_array
  
  !* Variables to store the values of the command-line arguments which are particular to this tutorial.
  !* Execute the tutorial with the "-h" command-line flag for more details
  real(rp)                  :: alpha
  real(rp)                  :: circle_radius
  real(rp), allocatable     :: circle_center(:)
  logical                   :: write_postprocess_data
  integer(ip)               :: num_uniform_refinement_steps
  integer(ip)               :: num_amr_steps
  integer(ip)               :: current_amr_step
  
  !* ### Main code  
  !* Here you can see the main code of this tutorial representing the entire flow
  !* of actions/computations performed by this tutorial.
  !*
  !* This code follows a similar structure as [tutorial 02 main code](../tutorial_02_poisson_sharp_circular_wave_amr/index.html#main-code),
  !* despite it is being executed in a significantly more complex, non-standard parallel execution environment.
  !* 
  !* In particular, all parallel tasks in world context execute the
  !* bulk of code lines of this tutorial, despite these are split into two levels by [[environment_t]]
  !* and assigned different duties (and data) at different levels. 
  !* We have devoted daunting efforts in order to hide as much as possible this
  !* complex execution environment to library users. The vast majority of TBPs associated to the
  !* library data types can be called safely from any task in world context.
  !* 
  !* The reader may also observe some new helper procedures in this main code.
  !* This tutorial declares extra procedures which are not necessary in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !*
  !* On one hand, This tutorial tackles a single-field PDE, this array is set up in the
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_coarse_fe_handler]] 
  !* helper subroutine to be a size-one array pointing to
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:coarse_fe_handler]].
  !* 
  !* On the other hand, [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_preconditioner()]]
  !* helper procedure customize the solver parameter required for each of the subproblems solved
  !* by the BDDC preconditioner at each level.

  call fempar_init()
  call setup_parameter_handler()
  call setup_context_and_environment()
  call get_tutorial_cla_values()
  
  current_amr_step = 0
  call setup_triangulation()
  call setup_problem_functions()
  call setup_strong_boundary_conditions()
  call setup_coarse_fe_handler()
  call setup_fe_space()
  call setup_discrete_solution()
  call setup_and_assemble_fe_affine_operator()
  call setup_preconditioner()
  call solve_system()  
  call compute_error()
  call setup_refinement_strategy()
  call output_handler_initialize()
  do
    call output_handler_write_current_amr_step()
    if ( current_amr_step == num_amr_steps ) exit
    current_amr_step = current_amr_step + 1
    call refinement_strategy%update_refinement_flags(triangulation)
    call setup_triangulation()
    call setup_fe_space()
    call redistribute_triangulation_and_fe_space()
    call setup_discrete_solution()
    call setup_and_assemble_fe_affine_operator()
    call solve_system()  
    call compute_error()
  end do
  call output_handler_finalize()
  
  call free_all_objects()
  call fempar_finalize()


  !* ### Helper procedures
contains

  !* #### Parameter handling
  !* The [[fempar_parameter_handler_names:parameter_handler]] is connected with the tutorial's command-line interface. 
  !* It parses the arguments provided to the command-line, and stores their values into a polymorphic dictionary of `<key,value>` pairs 
  !* (i.e., of type [parameterlist_t](http://victorsndvg.github.io/FPL/type/parameterlist_t.html)); a pointer to such dictionary can be obtained calling [[fempar_parameter_handler_t:get_values]], 
  !* as it will be required later on to create most of the objects provided by **`FEMPAR`**. 
  !*
  !* In [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_parameter_handler]] subroutine 
  !* we set up the system-wide variable referred to as [[fempar_parameter_handler_names:parameter_handler]].
  !* 
  !* **`FEMPAR`** users may define additional command-line arguments by means of a user-provided subroutine. 
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:define_tutorial_clas]] 
  !* subroutine (defined right below) is indeed the one that defines the command-line arguments for this particular 
  !* tutorial
  subroutine setup_parameter_handler()
    call parameter_handler%process_parameters(& 
         define_user_parameters_procedure = define_tutorial_clas, & 
         progname = tutorial_name, &
         version  = tutorial_version, &
         description = tutorial_description, &
         authors     = tutorial_authors)
  end subroutine setup_parameter_handler

  !* In [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:define_tutorial_clas]] subroutine 
  !* we use [[fempar_parameter_handler_names:parameter_handler]] in order to add and describe the command-line arguments 
  !* which are particular to this tutorial program. 
  !*
  !* Describing a command-line argument involves defining:
  !* 
  !* - a [[FPL]] dictionary _key_ (e.g., `ALPHA`),
  !* - a command-line argument _name_ (e.g., `--ALPHA`), 
  !* - a _default value_ for the command-line argument (in case it is not passed )(e.g., `200.0`), 
  !* - a _help message_, 
  !* - and (optionally) a set of admissible _choices_ for the command-line argument.
  !*
  !* The dictionary key can be used later on in order to get the value of the corresponding command-line argument
  !* or to override (update) it with a fixed value, thus ignoring the value provided to the command-line argument
  subroutine define_tutorial_clas()
    call parameter_handler%add("ALPHA",   &
                               "--ALPHA", & 
                               200.0_rp,  & 
                               "Scalar parameter which determines the sharpness of the solution's &
                               circular wave front")
    call parameter_handler%add("CIRCLE_RADIUS",   & 
                               "--CIRCLE_RADIUS", & 
                               0.7_rp, & 
                               "Radius of the circle/sphere that defines the solution's & 
                               circular wave front")
    call parameter_handler%add("CIRCLE_CENTER",     & 
                               "--CIRCLE_CENTER",   & 
                               [-0.05_rp,-0.05_rp,-0.05_rp], & 
                               "Cartesian coordinates of the center of the circle/sphere that & 
                               defines the solution's circular wave front")
    call parameter_handler%add("WRITE_POSTPROCESS_DATA",   & 
                               "--WRITE_POSTPROCESS_DATA", & 
                               .false., & 
                               "Enable/Disable generation of data files for later visualization & 
                               using, e.g., ParaView")
    call parameter_handler%add("NUM_UNIFORM_REFINEMENT_STEPS", &
                               "--NUM_UNIFORM_REFINEMENT_STEPS", & 
                               2, & 
                               "Number of uniform refinement steps right before starting the AMR loop")
    call parameter_handler%add("NUM_AMR_STEPS", &
                               "--NUM_AMR_STEPS", & 
                               5, & 
                               "Number of AMR loop steps")
  end subroutine define_tutorial_clas

  !* In [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:get_tutorial_cla_values]] subroutine 
  !* we use the [[fempar_parameter_handler_names:parameter_handler]] to obtain the values of the command-line arguments 
  !* which are particular to this tutorial.
  subroutine get_tutorial_cla_values()
    call parameter_handler%get("ALPHA", alpha)
    call parameter_handler%get("CIRCLE_RADIUS", circle_radius)
    call parameter_handler%getasarray("CIRCLE_CENTER", circle_center)
    call parameter_handler%get("NUM_UNIFORM_REFINEMENT_STEPS", num_uniform_refinement_steps)
    call parameter_handler%get("NUM_AMR_STEPS", num_amr_steps)
    call parameter_handler%get("WRITE_POSTPROCESS_DATA", write_postprocess_data)
    if ( environment%am_i_l1_task() ) then
      if ( environment%am_i_l1_root() ) then
        write(*,'(a70)')  repeat('=', 70)
        write(*,'(a70)') tutorial_name // ' CLA values'
        write(*,'(a70)')  repeat('=', 70)
        write(*,'(a46,f24.4)') 'ALPHA:'// repeat(' ', 80), alpha
        write(*,'(a46,f24.4)') 'CIRCLE_RADIUS:' // repeat(' ', 80), circle_radius
        if (size(circle_center) == 2) write(*,'(a46,3f12.4)') 'CIRCLE_CENTER:'// repeat(' ', 80), circle_center
        if (size(circle_center) == 3) write(*,'(a46,3f8.4)') 'CIRCLE_CENTER:' // repeat(' ', 80), circle_center
        write(*,'(a46,l24)') 'WRITE_POSTPROCESS_DATA:' // repeat(' ', 80), write_postprocess_data
        write(*,'(a46,i24)') 'NUM_UNIFORM_REFINEMENT_STEPS:' // repeat(' ', 80), num_uniform_refinement_steps
        write(*,'(a46,i24)') 'NUM_AMR_STEPS:' // repeat(' ', 80), num_amr_steps
        write(*,'(a)')  repeat('=', 70)
      end if
    end if   
  end subroutine get_tutorial_cla_values

  !* #### Environment initialization
  !* In [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_context_and_environment]] subroutine 
  !* a group of tasks (as many as specified to mpirun) is created.
  !* 
  !* Environment related parameters are also updated to preserve the coherence with this context
  !* (as [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:world_context]] is of type [[mpi_context_t]]).
  !* It split world_context into two subgroups of tasks (levels) composed by `world_context%get_num_tasks()-1` and a single task, resp.
  !* 
  !* When parameters are already update, the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:environment]] is created.
  subroutine setup_context_and_environment()
    call world_context%create()
    call parameter_handler%update(environment_num_levels_key, 2)
    call parameter_handler%update(environment_num_tasks_x_level_key, [world_context%get_num_tasks()-1,1])
    call environment%create(world_context, parameter_handler%get_values())
  end subroutine setup_context_and_environment

  !* #### Tringulation
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_triangulation]] helper procedure
  !* very much resemble [[tutorial_02_poisson_sharp_circular_wave_amr:setup_triangulation]] procedure 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !* In particular, when the value of current amr step is zero, setup triangulation calls the 
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_coarse_triangulation]] TBP of
  !* [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:triangulation]], 
  !* right after the latter triangulation is built. This procedure build a [[coarse_triangulation_t]] object. 
  !* This object is kept inside [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:triangulation]],
  !* although the user may have access to it via triangulation getters
  
  subroutine setup_triangulation()
     integer(ip) :: i, num_dims
     if ( current_amr_step == 0 ) then
       call parameter_handler%get(p4est_triang_num_dims_key, num_dims)
       if ( num_dims == 2 ) then
         call parameter_handler%update(p4est_triang_domain_limits_key, &
                                       [0.0_rp,1.0_rp,0.0_rp,1.0_rp])
       else
         call parameter_handler%update(p4est_triang_domain_limits_key, &
                                       [0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0,1.0_rp])
       end if
       call triangulation%create(environment, parameter_handler%get_values())
       do i=1, num_uniform_refinement_steps
         call flag_all_cells_for_refinement()
         call triangulation%refine_and_coarsen()
         call triangulation%redistribute()
       end do
       call triangulation%setup_coarse_triangulation()
     else
       call triangulation%refine_and_coarsen()
     end if 
     if ( write_postprocess_data ) then
       call setup_my_rank_cell_array()
     end if 
  end subroutine setup_triangulation

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:flag_all_cells_for_refinement]] 
  !* walks through over the mesh cells, and for each cell, 
  !* sets a per-cell flag that tells the triangulation to refine the cell in order to obtain a uniformly refined mesh.   
  subroutine flag_all_cells_for_refinement()
    class(cell_iterator_t), allocatable :: cell
    call triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() ) 
      if ( cell%is_local() ) then
        call cell%set_for_refinement()
      end if  
      call cell%next()
    end do
    call triangulation%free_cell_iterator(cell)
  end subroutine flag_all_cells_for_refinement

  !* #### Analytical expressions
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_problem_functions]] 
  !* looks like [[tutorial_02_poisson_sharp_circular_wave_amr:setup_problem_functions]] 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !* It sets up the exact solution and source term objects. 
  subroutine setup_problem_functions()
      call source_term%create(triangulation%get_num_dims(),& 
                              alpha,circle_radius,circle_center)
      call exact_solution%create(triangulation%get_num_dims(),&
                                 alpha,circle_radius,circle_center)
  end subroutine setup_problem_functions


  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_strong_boundary_conditions]] 
  !* looks like [[tutorial_02_poisson_sharp_circular_wave_amr:setup_strong_boundary_conditions]] 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !* It sets up the strong boundary conditions
  !* instance conformally with how the VEFs of the triangulation laying on the boundary are flagged
  !* with set identifiers. In the loop, it defines a strong boundary condition to be
  !* imposed for each of the regions that span the boundary of the unit box domain, and the same
  !* function, i.e., exact solution, to be imposed on all these regions, as required by the model
  !* problem of this tutorial.
  subroutine setup_strong_boundary_conditions()
    integer(ip) :: i, boundary_ids
    boundary_ids = merge(8, 26, triangulation%get_num_dims() == 2) 
    call strong_boundary_conditions%create()
    do i = 1, boundary_ids
       call strong_boundary_conditions%insert_boundary_condition(boundary_id=i, &
                                                                 field_id=1, &
                                                                 cond_type=component_1, &
                                                                 boundary_function=exact_solution)
    end do
  end subroutine setup_strong_boundary_conditions

  !* #### System setup
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_strong_boundary_conditions]] 
  !* implements the `CG` strategy as in [[tutorial_02_poisson_sharp_circular_wave_amr:setup_strong_boundary_conditions]] 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !*
  !* When the value of current AMR step is zero, it sets up the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:fe_space]], 
  !* instance using exactly the same sequence of calls as [[tutorial_01_poisson_sharp_circular_wave:setup_fe_space]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]].
  !* 
  !* When the value of current AMR step is larger than zero, it calls the [[serial_fe_space_t:refine_and_coarsen]] TBP
  !* of [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:fe_space]] to build a new global FE space \(V_h\)
  !* from the newly generated triangulation. Optionally, this TBP can by supplied with a FE function \(u_h\) 
  !* (or, more generally, an arbitrary number of the
  subroutine setup_fe_space()
    type(string) :: fes_field_types(1), fes_ref_fe_types(1)
    if ( current_amr_step == 0 ) then
       call parameter_handler%update(fes_num_fields_key, 1 )
       call parameter_handler%update(fes_same_ref_fes_all_cells_key, .true.)
       call parameter_handler%update(fes_ref_fe_types_key, fes_ref_fe_types )
       fes_field_types(1) = String(field_type_scalar)
       call parameter_handler%update(fes_field_types_key,fes_field_types)  
       fes_ref_fe_types(1) = String(fe_type_lagrangian)
       call parameter_handler%update(fes_ref_fe_types_key, fes_ref_fe_types )
       call parameter_handler%update(fes_ref_fe_conformities_key, [.true.] )
       call fe_space%create( triangulation  = triangulation, &
                             conditions     = strong_boundary_conditions, &
                             parameters     = parameter_handler%get_values())
       call fe_space%set_up_cell_integration()
       call fe_space%setup_coarse_fe_space(coarse_fe_handlers)
    else
       call fe_space%refine_and_coarsen()
    end if
 end subroutine setup_fe_space

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_strong_boundary_conditions]] 
  !* calls the [[triangulation_t:redistribute]] TBP of [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:triangulation]], 
  !* which dynamically balances the computational load by redistributing the adaptive mesh among the parallel tasks. 
  !* The default criteria is to balance the number of cells in each task. 
  !* Alternatively, the user might associate to each cell a partition weight. 
  !* In this case, the primitive balances the sums of the cell partition weights among processors. 
  !* The data that the user might have attached to the mesh objects (i.e., cells and VEFs set identifiers) is also migrated. 
  !* On the other hand, it also calls the [[par_fe_space_t::redistribute]] TBP of [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:fe_space]],
  !* which migrates the data that fe space holds conformally to how the triangulation has been
  !* redistributed. Optionally, this TBP can by supplied with a FE function \(u_h\) 
  !* (or, more generally,an arbitrary number of them)
 subroutine redistribute_triangulation_and_fe_space()
   call triangulation%redistribute()
   call fe_space%redistribute()
   if ( write_postprocess_data ) then
     call setup_my_rank_cell_array()
   end if
 end subroutine redistribute_triangulation_and_fe_space

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_discrete_solution]] 
  !* looks like [[tutorial_02_poisson_sharp_circular_wave_amr:setup_discrete_solution]] 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !*
  !* It sets up the discrete solution object, of type [[fe_function_t]]. This FEMPAR
  !* data type represents an element of \(V_h\) , the FE solution \(u_h\) in the case of this tutorial. 
 subroutine setup_discrete_solution()
   call discrete_solution%create(fe_space) 
   call fe_space%interpolate_dirichlet_values(discrete_solution)
 end subroutine setup_discrete_solution

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_and_assemble_fe_affine_operator]] 
  !* looks like [[tutorial_02_poisson_sharp_circular_wave_amr:setup_and_assemble_fe_affine_operator]] 
  !* in [tutorial 02](../tutorial_02_poisson_sharp_circular_wave_amr/index.html).
  !*
  !* When the value of current AMR step is zero, it sets up the [[discrete_integration_t]], 
  !* instance using exactly the same sequence of calls as [[tutorial_01_poisson_sharp_circular_wave:setup_and_assemble_fe_affine_operator]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]].
  !* 
  !* When the value of current AMR step is larger than zero, it calls the [[fe_affine_operator_t:reallocate_after_remesh]] TBP
  !* to adapt it to the new scenario after refining/coarsening the mesh.
 subroutine setup_and_assemble_fe_affine_operator()
   if ( current_amr_step == 0 ) then
      call cg_discrete_integration%set_source_term(source_term)
      call cg_discrete_integration%set_boundary_function(discrete_solution)
      call fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                       diagonal_blocks_symmetric_storage = [ .true. ], &
                                       diagonal_blocks_symmetric         = [ .true. ], &
                                       diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                       fe_space                          = fe_space, &
                                       discrete_integration              = cg_discrete_integration )
   else
      call fe_affine_operator%reallocate_after_remesh()
   end if
   call fe_affine_operator%compute()
 end subroutine setup_and_assemble_fe_affine_operator

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_and_assemble_fe_affine_operator]] 
  !* tackles a single-field PDE, this array is set up to be a size-one array pointing to coarse fe handler.
  subroutine setup_coarse_fe_handler()
    implicit none
    integer(ip) :: istat
    allocate(coarse_fe_handlers(1), stat=istat)
    check(istat==0)
    call coarse_fe_handler%create(parameter_handler%get_values()) 
    coarse_fe_handlers(1)%p => coarse_fe_handler
  end subroutine setup_coarse_fe_handler

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_preconditioner]]
  !* first builds the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc_parameters]] dictionary.
  !* This call takes the solver-related CLA values provided by parameter handler,
  !* and populates mlbddc parameters such that the same parameter values
  !* are used for the solvers of all subproblems that [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc_parameters]]
  !* handles internally (e.g., the Dirichlet and constrained Neumann subproblems). 
  !*
  !* The actual set up of the preconditioner occurs in [[mlbddc_t:create]] procedure. 
  !* We note that [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:mlbddc]] directly receives the FE affine operator, 
  !* instead of the coefficient matrix that it holds inside. This lets the BDDC framework to access to the 
  !* application FE discretization-related data, so that this information can be exploited when building
  !* an optimal preconditioner for the PDE problem at hand.
  subroutine setup_preconditioner()
    call setup_mlbddc_parameters_from_reference_parameters(environment, &
                                                           parameter_handler%get_values(), &
                                                           mlbddc_parameters)
    call mlbddc%create( fe_affine_operator, mlbddc_parameters )
  end subroutine setup_preconditioner

  !* #### System solving
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:solve_system]] 
  !* forces the Conjugate Gradients solver, as this is the most suitable iterative solver for the Poisson PDE. 
  !* 
  !* @note
  !* The rest of CLA values linked to iterative linear solver t are not forced, so that the user may choose,
  !* e.g., among several convergence criteria and related solver tolerances, or whether to print on
  !* screen or not the convergence history of the solver. An iterative solver needs a matrix and a
  !* preconditioner to solve the system.
  !* endnote
  !*
  !* Next, the coefficient matrix and preconditioner to be used by iterative linear solver are set.
  !* We extract a pointer to the free nodal values of our FE function, which is the place in which we will store the result.
  !* We solve the problem with the matrix already associated, the RHS obtained from the fe_affine_operator using get_translation, and 
  !* putting the result in dof_values.
  !* Once we have obtained the values corresponding to DoFs which are actually free, we have to update the values corresponding to DoFs
  !* which are hanging (constrained). This is required as, in this tutorial, we are using AMR on non-conforming meshes.
 subroutine solve_system()
   class(vector_t), pointer :: dof_values
   if ( current_amr_step == 0 ) then
      call iterative_linear_solver%create(environment)
      call parameter_handler%update(ils_type_key, cg_name)
      call iterative_linear_solver%set_type_from_pl(parameter_handler%get_values())
      call iterative_linear_solver%set_parameters_from_pl(parameter_handler%get_values())
      call iterative_linear_solver%set_operators(fe_affine_operator%get_matrix(), mlbddc)
   else 
      call iterative_linear_solver%reallocate_after_remesh()
   end if
   dof_values => discrete_solution%get_free_dof_values()
   call iterative_linear_solver%solve(fe_affine_operator%get_translation(),dof_values)
   call fe_space%update_hanging_dof_values(discrete_solution)
 end subroutine solve_system

  !* #### Error computation
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:compute_error]] 
  !* computes \(e^2_k\) for each \(K \in \mathcal{T}_h\), and the global error \(e\).
 subroutine compute_error()
   real(rp) :: global_error_energy_norm
   type(coarse_triangulation_t), pointer :: coarse_triangulation
   type(coarse_fe_space_t), pointer :: coarse_fe_space
   call error_estimator%create(fe_space,parameter_handler%get_values())
   call error_estimator%set_exact_solution(exact_solution)
   call error_estimator%set_discrete_solution(discrete_solution)
   if ( environment%am_i_l1_task() ) then
     if ( environment%am_i_l1_root() ) then
        write(*,'(a70)')  repeat('=', 70)
        write(*,'(a70)') tutorial_name // ' results'
        write(*,'(a70)')  repeat('=', 70)
      end if  
   end if
   call error_estimator%compute_local_true_errors()
   call error_estimator%compute_local_estimates()
   global_error_energy_norm = sum(error_estimator%get_sq_local_true_error_entries())
   call environment%w_sum(global_error_energy_norm)
   global_error_energy_norm = sqrt(global_error_energy_norm)
   if ( environment%am_i_l1_task() ) then
      if ( environment%am_i_l1_root() ) then
        write(*,'(a,i3)') 'AMR step:', current_amr_step
        write(*,'(a46,i24)') repeat(' ', 4) // 'NUM_CELLS:' // repeat(' ', 80), triangulation%get_num_global_cells()
        write(*,'(a46,i24)') repeat(' ', 4) // 'NUM_DOFS:' // repeat(' ', 80), fe_space%get_num_global_dofs()
        write(*,'(a46,e24.10)') repeat(' ', 4) // 'GLOBAL ERROR ENERGY NORM:'// repeat(' ', 80), global_error_energy_norm
      end if 
   end if
   call world_context%barrier() ! Synchronize all tasks
   if (environment%am_i_lgt1_task()) then
     coarse_triangulation => triangulation%get_coarse_triangulation()
     coarse_fe_space => fe_space%get_coarse_fe_space()
     write(*,'(a46,i24)') repeat(' ', 4) // 'NUM_COARSE_CELLS: ' // repeat(' ', 80), coarse_triangulation%get_num_cells()
     write(*,'(a46,i24)') repeat(' ', 4) // 'NUM_COARSE_DOFS:'// repeat(' ', 80)   , coarse_fe_space%get_total_num_dofs()
   end if
   call world_context%barrier() ! Synchronize all tasks
   if ( environment%am_i_l1_task() ) then
      if ( environment%am_i_l1_root() ) then
        if (current_amr_step == num_amr_steps) then
          write(*,'(a70)')  repeat('=', 70)
        end if
      end if 
   end if
  end subroutine compute_error

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_refinement_strategy]] 
  !* forces refinement_strategy to use parameter values \(\alpha_r = 0.1\) and \(\alpha_r = 0.05\)
  subroutine setup_refinement_strategy()
    call parameter_handler%update(key = ffrs_refinement_fraction_key, value = 0.10_rp)
    call parameter_handler%update(key = ffrs_coarsening_fraction_key, value = 0.05_rp)
    call parameter_handler%update(key = ffrs_max_num_mesh_iterations_key, value = num_amr_steps )
    call refinement_strategy%create(error_estimator,parameter_handler%get_values())
  end subroutine setup_refinement_strategy

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:setup_my_rank_cell_array]] 
  !* fills [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:my_rank_cell_array]] instance 
  !* with as many entries as local cells in each parallel task, and
  !* for all of these entries, it holds the parallel task identifier.
  subroutine setup_my_rank_cell_array()
    call my_rank_cell_array%resize(triangulation%get_num_local_cells())
    call my_rank_cell_array%init(real(environment%get_l1_rank()+1,rp))
  end subroutine setup_my_rank_cell_array

  !* #### Results plotting
  !* The handling the output writers is distributed in three different subroutines:
  !*
  !* - [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_initialize]] 
  !* - [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_write_current_amr_step]] 
  !* - [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_finalize]] 
  !* 
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_initialize]] 
  !* create the [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler]] 
  !* and add the select objects to be plotted.    
  subroutine output_handler_initialize()
    if (write_postprocess_data) then
      call parameter_handler%update(output_handler_static_grid_key, value=.false.)
      call output_handler%create(parameter_handler%get_values())
      call output_handler%attach_fe_space(fe_space)
      call output_handler%add_fe_function(discrete_solution, 1, 'solution')
      call output_handler%add_cell_vector(error_estimator%get_sq_local_estimates(), 'cell_error_energy_norm_squared')
      call output_handler%add_cell_vector(my_rank_cell_array, 'subdomain_partition')
      call output_handler%open()
    end if   
  end subroutine output_handler_initialize

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_write_current_amr_step]] 
  !* is in charge of plotting every AMR step.
  subroutine output_handler_write_current_amr_step()
    if (write_postprocess_data) then
      call output_handler%append_time_step(real(current_amr_step,rp))
      call output_handler%write()
    end if  
  end subroutine output_handler_write_current_amr_step

  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:output_handler_finalize]] 
  !* finalize the writing procedure.
  subroutine output_handler_finalize()
    if (write_postprocess_data) then
      call output_handler%close()
    end if  
  end subroutine output_handler_finalize

  !* #### Environment finalization
  !* Subroutine [[tutorial_03_poisson_sharp_circular_wave_parallel_amr:free_all_objects]] 
  !* free all the objects created in the course of the program.
  subroutine free_all_objects()
    call output_handler%free()
    call error_estimator%free()
    call iterative_linear_solver%free()
    call mlbddc%free()
    call mlbddc_parameters%free()
    call fe_affine_operator%free()
    call discrete_solution%free()
    call fe_space%free()
    call strong_boundary_conditions%free()
    call source_term%free()
    call exact_solution%free()
    call triangulation%free()
    call my_rank_cell_array%free()
    call environment%free()
    call world_context%free(.true.)  
  end subroutine free_all_objects
  
end program tutorial_03_poisson_sharp_circular_wave_parallel_amr
