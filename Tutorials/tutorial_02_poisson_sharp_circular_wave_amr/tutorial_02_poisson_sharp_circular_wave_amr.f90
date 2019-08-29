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
!* The [[tutorial_02_poisson_sharp_circular_wave_amr]],
program tutorial_02_poisson_sharp_circular_wave_amr
!* uses the [[fempar_names]] to import all **`FEMPAR`** library symbols (i.e., derived types, parameter constants, system-wide variables, etc.),
!* and a set of tutorial-specific module units, 
!*
!* - [[tutorial_01_functions_names]]
!* - [[tutorial_01_error_estimator_names]]
!* - [[tutorial_01_discrete_integration_names]]
!*
!* These modules are not part of the **`FEMPAR`** library, but developed specifically for the problem at hand. Each of these modules defines a tutorial-specific data type and its TBPs.
  use fempar_names
  use tutorial_01_functions_names
  use tutorial_01_error_estimator_names
  use tutorial_01_discrete_integration_names

!* ### Variable definition
!* This tutorial declares:
!*
!* - a set of parameter constants, typically the tutorial name, authors, problem description, etc. (to be output on screen on demand by the user),
!* - the tutorial data type instances in charge of the FE simulation, such as the triangulation (mesh) of the computational domain or the FE space from which the approximate solution of the PDE is sought, 
!* - and a set of variables to hold the values of the Command-Line-Arguments (CLAs) which are specific to the tutorial.
  implicit none
#include "debug.i90"
  
  character(*), parameter :: tutorial_name = "tutorial_02_poisson_sharp_circular_wave_amr"
  character(*), parameter :: tutorial_version = "v1.0.0"
  character(*), parameter :: tutorial_description = "Solves a 2D/3D Poisson problem on the unit &
                                                     square/cube, where the solution has a sharp circular & 
                                                     wave front parametrizable via command-line arguments. &
                                                     The discretization of the domain is performed using a p4est &
                                                     triangulation, i.e., a mesh which can by h-adapted  & 
                                                     during the simulation."
  character(*), parameter :: tutorial_authors = "Santiago Badia and Alberto F. Martín"
  
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:world_context]] is a software abstraction 
  !* for a group of parallel tasks (a.k.a. processes) and the communication layer that orchestrates their concurrent execution.
  !* In this case is declared of type [[serial_context_t]] as it is designed to work in serial computing environments.
  type(serial_context_t)                      :: world_context
  !* the [[tutorial_02_poisson_sharp_circular_wave_amr:environment]], of type [[environment_t]], organizes the tasks of the context
  !* from which it is set up into subgroups of tasks, referred to as levels, and builds up additional communication mechanisms 
  !* to communicate tasks belonging to different levels.
  type(environment_t)                         :: environment
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:triangulation]] object, of type [[p4est_serial_triangulation_t]], provides the mesh. 
  !* In this case, we consider a serial p4est triangulation, i.e., 
  !* a triangulation that is not distributed among processors, and that can be h-adapted during the simulation.
  type(p4est_serial_triangulation_t)          :: triangulation
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:fe_space]], of type [[serial_fe_space_t]] is the global finite element space to be used.
  type(serial_fe_space_t)                     :: fe_space
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:strong_boundary_conditions]], of type [[strong_boundary_conditions_t]], is the object that provides the strong boundary conditions definition.
  type(strong_boundary_conditions_t)          :: strong_boundary_conditions 
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:source_term]], of type [[tutorial_01_functions_names:sharp_circular_wave_source_term_t]], is a scalar-valued function with the right hand side of the PDE problem at hand.
  type(sharp_circular_wave_source_term_t)     :: source_term
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:exact_solution]], of type [[tutorial_01_functions_names:sharp_circular_wave_solution_t]], is a scalar-valued function with the exact (analytical) solution of the PDE problem at hand.
  type(sharp_circular_wave_solution_t)        :: exact_solution
  !* [[tutorial_02_poisson_sharp_circular_wave_amr:fe_affine_operator]], of type [[fe_affine_operator_t]],
  !* represents the affine operator the solution of which is the one we want, i.e., \(B = Ax-f\).
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)                  :: fe_affine_operator
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:discrete_solution]], of type [[fe_function_t]], 
  !* is the object belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                         :: discrete_solution
  !* [[tutorial_01_discrete_integration_names:base_discrete_integration_t]] provides the base definition of the bilinear and linear forms for the problem at hand.
  type(cg_discrete_integration_t), target     :: cg_discrete_integration
  type(dg_discrete_integration_t), target     :: dg_discrete_integration
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:direct_solver]], of type [[direct_solver_t]], provides an interface to several external sparse direct solver packages (PARDISO, UMFPACK).
  type(direct_solver_t)                       :: direct_solver
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler]] object, of type [[output_handler_t]], is used to generate the simulation data files for later visualization using, e.g., ParaView.
  type(output_handler_t)                      :: output_handler
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:error_estimator]] object, of type [[tutorial_01_error_estimator_names:poisson_error_estimator_t]], is used to compute the (square) of the error energy norm for all cells, i.e., \(||u-u_h||_E,K\).
  type(poisson_error_estimator_t)             :: error_estimator
  !* The [[tutorial_02_poisson_sharp_circular_wave_amr:refinement_strategy]] object, of type [[fixed_fraction_refinement_strategy_t]], 
  !* implements the strategy used for mesh refinement and coarsening.
  type(fixed_fraction_refinement_strategy_t)  :: refinement_strategy
  
  !* Variables to store the values of the command-line arguments which are particular to this tutorial.
  !* Execute the tutorial with the "-h" command-line flag for more details
  real(rp)                  :: alpha
  real(rp)                  :: circle_radius
  real(rp), allocatable     :: circle_center(:)
  character(:), allocatable :: fe_formulation
  logical                   :: write_postprocess_data
  integer(ip)               :: num_uniform_refinement_steps
  integer(ip)               :: num_amr_steps
  integer(ip)               :: current_amr_step
  
  !* ### Main code  
  !* Here you can see the main code of this tutorial representing the entire flow
  !* of actions/computations performed by this tutorial.
  !*
  !* As there is almost a one-to-one mapping among the data type instances and the
  !* helper procedures called in this tutorial, we will introduce them in the sequel step-by-step along
  !* with code snippets of the corresponding [helper procedures](#helper-procedures).
  !* 
  !* In order to illustrate the AMR capabilities in **`FEMPAR`**, while generating a suitable mesh for Poisson problem, 
  !* this tutorial performs an AMR loop comprising the following steps:
  !* 
  !* 1. Generate a conforming mesh \(\mathcal{T}_h\) by uniformly refining, a number user-defined steps, a single-cell coarse mesh \(\mathcal{C}_h\) representing the unit box domain (i.e., \(\Omega\) in Poisson problem.
  !* 2. Compute an approximate FE solution \(u_h\) using the current mesh \(\mathcal{T}_h\).
  !* 3. Compute \(e_K^2\) for all cells \(K \in \mathcal{T}_h\).
  !* 4. Given user-defined refinement and coarsening fractions, denoted by \(\alpha_r\) and \(\alpha_c\), resp., find thresholds \(\theta_r\) and \(\theta_c\) such that the number of cells with \(e_K >\theta_r\)  (resp., \(e_K < \theta_c\))  is (approximately) a fraction \(\alpha_r\)  (resp., \(\alpha_c\)) of the number of cells in \(\mathcal{T}_h\).
  !* 5. Refine and coarsen the mesh cells, i.e., generate a new mesh \(\mathcal{T}_h\),  accordingly to the input provided by the previous step.
  !* 6. Repeat steps (2)-(5) a number of user-defined steps.
  call fempar_init()
  call setup_parameter_handler()
  call get_tutorial_cla_values()
  call setup_context_and_environment()
  
  current_amr_step = 0
  call setup_triangulation()
  call setup_problem_functions()
  call setup_strong_boundary_conditions()
  call setup_fe_space()
  call setup_discrete_solution()
  call setup_and_assemble_fe_affine_operator()
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
  !* In [[tutorial_02_poisson_sharp_circular_wave_amr:setup_parameter_handler]] subroutine 
  !* we set up the system-wide variable referred to as [[fempar_parameter_handler_names:parameter_handler]].
  !* 
  !* **`FEMPAR`** users may define additional command-line arguments by means of a user-provided subroutine. 
  !* [[tutorial_02_poisson_sharp_circular_wave_amr:define_tutorial_clas]] 
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

  !* In [[tutorial_02_poisson_sharp_circular_wave_amr:define_tutorial_clas]] subroutine 
  !* we use [[fempar_parameter_handler_names:parameter_handler]] in order to add and describe the command-line arguments 
  !* which are particular to this tutorial program. 
  !*
  !* Describing a command-line argument involves defining:
  !* 
  !* - a [[FPL]] dictionary _key_ (e.g., `FE_FORMULATION`),
  !* - a command-line argument _name_ (e.g., `--FE_FORMULATION`), 
  !* - a _default value_ for the command-line argument (in case it is not passed )(e.g., `CG`), 
  !* - a _help message_, 
  !* - and (optionally) a set of admissible _choices_ for the command-line argument.
  !*
  !* The dictionary key can be used later on in order to get the value of the corresponding command-line argument
  !* or to override (update) it with a fixed value, thus ignoring the value provided to the command-line argument
  subroutine define_tutorial_clas()
    call parameter_handler%add("FE_FORMULATION",   &
                               "--FE_FORMULATION", & 
                               "CG",  & 
                               "Select Finite Element formulation for the problem at hand; & 
                               either Continuous (CG) or Discontinuous Galerkin (DG)", &
                               choices="CG,DG")
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

  !* In [[tutorial_02_poisson_sharp_circular_wave_amr:get_tutorial_cla_values]] subroutine 
  !* we use the [[fempar_parameter_handler_names:parameter_handler]] to obtain the values of the command-line arguments 
  !* which are particular to this tutorial.
  subroutine get_tutorial_cla_values()
    call parameter_handler%getasstring("FE_FORMULATION", fe_formulation)
    call parameter_handler%get("ALPHA", alpha)
    call parameter_handler%get("CIRCLE_RADIUS", circle_radius)
    call parameter_handler%getasarray("CIRCLE_CENTER", circle_center)
    call parameter_handler%get("NUM_UNIFORM_REFINEMENT_STEPS", num_uniform_refinement_steps)
    call parameter_handler%get("NUM_AMR_STEPS", num_amr_steps)
    call parameter_handler%get("WRITE_POSTPROCESS_DATA", write_postprocess_data)
    write(*,'(a54)')  repeat('=', 54)
    write(*,'(a54)') tutorial_name // ' CLA values'
    write(*,'(a54)')  repeat('=', 54)
    write(*,'(a30,a24)') 'FE_FORMULATION:' // repeat(' ', 80), fe_formulation
    write(*,'(a30,f24.4)') 'ALPHA:'// repeat(' ', 80), alpha
    write(*,'(a30,f24.4)') 'CIRCLE_RADIUS:' // repeat(' ', 80), circle_radius
    if (size(circle_center) == 2) write(*,'(a30,3f12.4)') 'CIRCLE_CENTER:'// repeat(' ', 80), circle_center
    if (size(circle_center) == 3) write(*,'(a30,3f8.4)') 'CIRCLE_CENTER:' // repeat(' ', 80), circle_center
    write(*,'(a30,l24)') 'WRITE_POSTPROCESS_DATA:' // repeat(' ', 80), write_postprocess_data
    write(*,'(a30,i24)') 'NUM_UNIFORM_REFINEMENT_STEPS:' // repeat(' ', 80), num_uniform_refinement_steps
    write(*,'(a30,i24)') 'NUM_AMR_STEPS:' // repeat(' ', 80), num_amr_steps
    write(*,'(a54)')  repeat('=', 54)
  end subroutine get_tutorial_cla_values

  !* #### Environment initialization
  !* In [[tutorial_02_poisson_sharp_circular_wave_amr:setup_context_and_environment]] subroutine 
  !* a single-task group of tasks is created.
  !* 
  !* Environment related parameters are also updated to preserve the coherence with this context
  !* (as [[tutorial_02_poisson_sharp_circular_wave_amr:world_context]] is of type [[serial_context_t]]).
  !* 
  !* When parameters are already update, the [[tutorial_02_poisson_sharp_circular_wave_amr:environment]] is created.
  !*
  !* @note
  !* As `world_context` is of type [[serial_context_t]], 
  !* any other environment configuration is not possible
  !* @endnote
  subroutine setup_context_and_environment()
    call world_context%create()
    call parameter_handler%update(environment_num_levels_key, 1)
    call parameter_handler%update(environment_num_tasks_x_level_key, [1])
    call environment%create(world_context, parameter_handler%get_values())
  end subroutine setup_context_and_environment

  !* #### Tringulation
  !* In [[tutorial_02_poisson_sharp_circular_wave_amr:setup_triangulation]] subroutine 
  !* one may readily observe that the setup triangulation helper subroutine behaves 
  !* differently depending on the value of current AMR step. 
  !* 
  !* When its value is zero the subroutine first generates a triangulation of the unit box domain which is composed of a single
  !* brick cell, i.e., a forest-of-octrees with a single adaptive octree. Then, it enters a loop
  !* in which the octree root is refined uniformly a number of steps  in order to generate \(\mathcal{T}_h\). 
  !* 
  !* The loop relies on the [[tutorial_02_poisson_sharp_circular_wave_amr:flag_all_cells_for_refinement]] helper procedure
  !*  and the [[p4est_base_triangulation_t:refine_and_coarsen]] TBP of triangulation. 
  !* 
  !* The first procedure walks through over the mesh cells, and for each cell, 
  !* sets a per-cell flag that tells the triangulation to refine the cell in order to obtain a uniformly refined mesh. 
  !* 
  !* On the other hand, the [[p4est_base_triangulation_t:refine_and_coarsen]] TBP adapts the triangulation 
  !* based on the cell flags set by the user, while transferring data that the user might have attached to the mesh objects 
  !* (e.g., cells and VEFs set identifiers) to the new mesh objects generated after mesh adaptation.
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
       end do
     else
       call triangulation%refine_and_coarsen()
     end if 
  end subroutine setup_triangulation

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:flag_all_cells_for_refinement]] 
  !* walks through over the mesh cells, and for each cell, 
  !* sets a per-cell flag that tells the triangulation to refine the cell in order to obtain a uniformly refined mesh.   
  subroutine flag_all_cells_for_refinement()
    class(cell_iterator_t), allocatable :: cell
    call triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )  
      call cell%set_for_refinement()
      call cell%next()
    end do
    call triangulation%free_cell_iterator(cell)
  end subroutine flag_all_cells_for_refinement

  !* #### Analytical expressions
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_problem_functions]] 
  !* sets up the exact solution and source term objects. 
  !* These represent the exact (analytical) solution  \(u\) and source term \(f\) of our problem. 
  !* The program does not need to implement the Dirichlet function \(g\), 
  !* as it is equivalent to \(u\) in our model problem. 
  !*
  !* While we do not cover their implementation here, 
  !* the reader is encouraged to inspect the [[tutorial_01_functions_names]] module. 
  !* Essentially, this module implements a pair of program-specific data types extending
  !* the so-called [[scalar_function_t]] **`FEMPAR`** data type. 
  !* 
  !* This latter data type represents an scalar-valued function \(h\), 
  !* and provides interfaces for the evaluation of \(h(x)\), \(\nabla h(x)\), etc., 
  !* with \(x \in \overline{\Omega}\), to be implemented by type extensions. 
  !*
  !* In particular, this tutorial requires \(u(x)\) and \(\nabla u(x)\) for the imposition of Dirichlet BCs, 
  !* and the evaluation of the energy norm, resp., and \(f(x)\) for the evaluation of the source term.
  subroutine setup_problem_functions()
      call source_term%create(triangulation%get_num_dims(),& 
                              alpha,circle_radius,circle_center)
      call exact_solution%create(triangulation%get_num_dims(),&
                                 alpha,circle_radius,circle_center)
  end subroutine setup_problem_functions

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_strong_boundary_conditions]] 
  !* sets up the strong boundary conditions object. With this object one can define
  !* the regions of the domain boundary on which to impose strong BCs, along with the function to
  !* be imposed on each of these regions. 
  !* It is required for the FE formulation, in this method, Dirichlet BCs are imposed strongly. 
  !* 
  !* The rationale behind this subroutine is as follows.
  !* Any **`FEMPAR`** triangulation handles internally (and provides on client demand) a set identifier
  !* (i.e., an integer number) per each Vertex, Edge, and Face (VEF) of the mesh. On the other
  !* hand, it is assumed that the mesh generator from which the triangulation is imported classifies
  !* the boundary of the domain into geometrical regions, and that, when the mesh is generated, it
  !* assigns the same set identifier to all VEFs which belong to the same geometrical region. 
  !* 
  !* For example, for \(d = 2\), the uniform mesh generator classifies the domain into 9 regions
  !*
  !* - the four corners of the box, which are assigned identifiers \(1, \dots, 4\), 
  !* - the four faces of the box, which are assigned identifiers \(5, \dots, 8\), 
  !* - and the interior of the box, which is assigned identifier \(9\). 
  !*
  !* For \(d = 3\), we have 27 regions, i.e., 8 corners, 12 edges, 6 faces, and the interior of the domain. (At
  !* this point, the reader should be able to grasp where the numbers 8 and 26 in
  !* [[tutorial_02_poisson_sharp_circular_wave_amr:setup_strong_boundary_conditions]] come from.) 
  !*
  !* With the aforementioned in mind, subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_strong_boundary_conditions]]
  !* sets up the strong boundary conditions
  !* instance conformally with how the VEFs of the triangulation laying on the boundary are flagged
  !* with set identifiers. In the loop, it defines a strong boundary condition to be
  !* imposed for each of the regions that span the boundary of the unit box domain, and the same
  !* function, i.e., exact solution, to be imposed on all these regions, as required by the model
  !* problem of this tutorial.
  subroutine setup_strong_boundary_conditions()
    integer(ip) :: i, boundary_ids
    if ( fe_formulation == "CG" ) then
      boundary_ids = merge(8, 26, triangulation%get_num_dims() == 2) 
      call strong_boundary_conditions%create()
      do i = 1, boundary_ids
        call strong_boundary_conditions%insert_boundary_condition(boundary_id=i, &
                                                                  field_id=1, &
                                                                  cond_type=component_1, &
                                                                  boundary_function=exact_solution)
      end do
    end if 
  end subroutine setup_strong_boundary_conditions

  !* #### System setup
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_fe_space]]
  !* follows the same pattern as [[tutorial_02_poisson_sharp_circular_wave_amr:setup_triangulation]] subroutine.
  !* When the value of current AMR step is zero, it sets up the [[tutorial_02_poisson_sharp_circular_wave_amr:fe_space]], 
  !* instance using exactly the same sequence of calls as [[tutorial_01_poisson_sharp_circular_wave:setup_fe_space]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]].
  !* 
  !* When the value of current AMR step is larger than zero, it calls the [[serial_fe_space_t:refine_and_coarsen]] TBP
  !* of [[tutorial_02_poisson_sharp_circular_wave_amr:fe_space]] to build a new global FE space \(V_h\)
  !* from the newly generated triangulation. Optionally, this TBP can by supplied with a FE function \(u_h\) 
  !* (or, more generally, an arbitrary number of them). 
  !* 
  !* In such a case, [[serial_fe_space_t:refine_and_coarsen]] injects the FE function provided into the newly generated \(V_h\)
  !* by using a suitable FE projector for the FE technology being used, the Lagrangian interpolator
  !* in the case of this tutorial. 
  !* 
  !* @note
  !* This feature is required by numerical solvers of transient and/or non-linear PDEs.
  !* @endnote
  subroutine setup_fe_space()
    type(string) :: fes_field_types(1), fes_ref_fe_types(1)
    if ( current_amr_step == 0 ) then
      call parameter_handler%update(fes_num_fields_key, 1 )
      call parameter_handler%update(fes_same_ref_fes_all_cells_key, .true.)
      call parameter_handler%update(fes_ref_fe_types_key, fes_ref_fe_types )
      fes_field_types(1) = String(field_type_scalar); 
      call parameter_handler%update(fes_field_types_key,fes_field_types)  
      fes_ref_fe_types(1) = String(fe_type_lagrangian)
      call parameter_handler%update(fes_ref_fe_types_key, fes_ref_fe_types )
      if ( fe_formulation == "CG" ) then
        call parameter_handler%update(fes_ref_fe_conformities_key, [.true.] )
        call fe_space%create( triangulation  = triangulation, &
                              conditions     = strong_boundary_conditions, &
                              parameters     = parameter_handler%get_values())
      else if ( fe_formulation == "DG" ) then
        call parameter_handler%update(fes_ref_fe_conformities_key, [.false.] )
        call fe_space%create( triangulation  = triangulation, &
                              parameters     = parameter_handler%get_values())
      end if
      call fe_space%set_up_cell_integration()
      if ( fe_formulation == "DG" ) then
        call fe_space%set_up_facet_integration()
      end if
    else
      call fe_space%refine_and_coarsen()
    end if
  end subroutine setup_fe_space

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_discrete_solution]]
  !* sets up the discrete solution object, of type [[fe_function_t]]. This FEMPAR
  !* data type represents an element of \(V_h\) , the FE solution \(u_h\) in the case of this tutorial. 
  !*
  !* This subroutine allocates room for storing the DOF values of \(u_h\) , and then it
  !* computes the values of those DOFs of \(u_h\) which lay on a region of the domain boundary which
  !* is subject to strong BCs, the whole boundary of the domain in the case of this tutorial. 
  !* 
  !* The [[serial_fe_space_t:interpolate_dirichlet_values]] TBP of [[tutorial_02_poisson_sharp_circular_wave_amr:fe_space]] 
  !* interpolates exact solution, i.e., \(u\), which is extracted from strong boundary conditions, 
  !* using a suitable FE interpolator for the FE space at hand, i.e., the Lagrangian interpolator in the case of this tutorial.
  !* 
  !* @note
  !* **`FEMPAR`** supports interpolators for _curl-_ and _div-_conforming FE spaces as well. Interpolation in such FE spaces
  !* involves the numerical evaluation of the functionals (moments) associated to the DOFs of \(V_h\)
  !* @endnote
  subroutine setup_discrete_solution()
    call discrete_solution%create(fe_space) 
    if ( fe_formulation == "CG" ) then
      call fe_space%interpolate_dirichlet_values(discrete_solution)
    end if  
  end subroutine setup_discrete_solution

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_and_assemble_fe_affine_operator]]
  !* follows the same pattern as [[tutorial_02_poisson_sharp_circular_wave_amr:setup_triangulation]] subroutine.
  !* When the value of current AMR step is zero, it sets up the [[discrete_integration_t]], 
  !* instance using exactly the same sequence of calls as [[tutorial_01_poisson_sharp_circular_wave:setup_and_assemble_fe_affine_operator]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]].
  !* 
  !* When the value of current AMR step is larger than zero, it calls the [[fe_affine_operator_t:reallocate_after_remesh]] TBP
  !* to adapt it to the new scenario after refining/coarsening the mesh.
  subroutine setup_and_assemble_fe_affine_operator()
    class(discrete_integration_t), pointer :: discrete_integration
    if ( current_amr_step == 0 ) then
      if ( fe_formulation == "CG" ) then
        call cg_discrete_integration%set_source_term(source_term)
        call cg_discrete_integration%set_boundary_function(discrete_solution)
        discrete_integration => cg_discrete_integration
      else if ( fe_formulation == "DG" ) then 
        call dg_discrete_integration%set_source_term(source_term)
        call dg_discrete_integration%set_boundary_function(exact_solution)
        discrete_integration => dg_discrete_integration
      end if
      call fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                       diagonal_blocks_symmetric_storage = [ .true. ], &
                                       diagonal_blocks_symmetric         = [ .true. ], &
                                       diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                       fe_space                          = fe_space, &
                                       discrete_integration              = discrete_integration )
    else
      call fe_affine_operator%reallocate_after_remesh()
    end if 
    call fe_affine_operator%compute()
  end subroutine setup_and_assemble_fe_affine_operator

  !* #### System solving
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:solve_system]], 
  !* follows the same pattern as [[tutorial_02_poisson_sharp_circular_wave_amr:setup_triangulation]] subroutine.
  !* When the value of current AMR step is zero, it sets up the [[discrete_integration_t]], 
  !* instance using exactly the same sequence of calls as [[tutorial_01_poisson_sharp_circular_wave:solve_system]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]].
  !* 
  !* Apart from this, the reader should note the call to the
  !* [[direct_solver_t:update_hanging_dof_values]] TBP of [[tutorial_02_poisson_sharp_circular_wave_amr:fe_space]]. 
  !* This TBP computes the values of \(u_h\) corresponding to hanging DOFs using 
  !* the algebraic constraints required to enforce the conformity of \(V_h\).  
  !* This is required each time that the values of \(u_h\) corresponding to true
  !* DOFs are updated (e.g., after linear system solution). 
  !* 
  !* @warning
  !* Not calling it results in unpredictable errors during post-processing, 
  !* or even worse, during a non-linear and/or time-stepping iterative
  !* solution loop, in which the true DOF values of \(u_h\) are updated at each iteration.
  !* @endwarning
  !* 
  !* When the value of current AMR step is larger than zero, it calls the [[direct_solver_t:reallocate_after_remesh]] TBP
  !* to adapt it to the new scenario after refining/coarsening the mesh.

  subroutine solve_system()
    class(vector_t), pointer :: dof_values
    if ( current_amr_step == 0 ) then
      call direct_solver%set_type_from_pl(parameter_handler%get_values())
      call direct_solver%set_parameters_from_pl(parameter_handler%get_values())
      call direct_solver%set_matrix(fe_affine_operator%get_matrix())
    else 
      call direct_solver%reallocate_after_remesh()
    end if
    dof_values => discrete_solution%get_free_dof_values()
    call direct_solver%solve(fe_affine_operator%get_translation(),dof_values)
    if ( fe_formulation == "CG" ) then
      call fe_space%update_hanging_dof_values(discrete_solution)
    end if  
  end subroutine solve_system

  !* #### Error computation
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:compute_error]] 
  !* computes \(e^2_k\) for each \(K \in \mathcal{T}_h\), and the global error \(e\).
  !* 
  !* This subroutine behaves exactly as [[tutorial_01_poisson_sharp_circular_wave:compute_error]] in
  !* [[tutorial_01_poisson_sharp_circular_wave]], except that, while setup of the
  !* [[tutorial_02_poisson_sharp_circular_wave_amr:error_estimator]] only occurs in AMR step zero,
  !* error calculation is performed for every AMR step.
  subroutine compute_error()
    real(rp) :: global_error_energy_norm
    if ( current_amr_step == 0 ) then
      call error_estimator%create(fe_space,parameter_handler%get_values())
      call error_estimator%set_exact_solution(exact_solution)
      call error_estimator%set_discrete_solution(discrete_solution)
      write(*,'(a54)')  repeat('=', 54)
      write(*,'(a54)') tutorial_name // ' results'
      write(*,'(a54)')  repeat('=', 54)
    end if
    call error_estimator%compute_local_true_errors()
    call error_estimator%compute_local_estimates()
    global_error_energy_norm = sqrt(sum(error_estimator%get_sq_local_true_error_entries()))
    write(*,'(a,i3)') 'AMR step:', current_amr_step
    write(*,'(a30,i24)') repeat(' ', 4) // 'NUM_CELLS:' // repeat(' ', 80), triangulation%get_num_cells()
    write(*,'(a30,i24)') repeat(' ', 4) // 'NUM_DOFS:' // repeat(' ', 80), fe_space%get_total_num_dofs()
    write(*,'(a30,e24.10)') repeat(' ', 4) // 'GLOBAL ERROR ENERGY NORM:'// repeat(' ', 80), global_error_energy_norm
    if (current_amr_step == num_amr_steps) then
      write(*,'(a54)')  repeat('=', 54)
    end if
  end subroutine compute_error

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:setup_refinement_strategy]] 
  !* forces refinement_strategy to use parameter values \(\alpha_r = 0.1\) and \(\alpha_r = 0.05\)
  subroutine setup_refinement_strategy()
    call parameter_handler%update(key = ffrs_refinement_fraction_key, value = 0.10_rp)
    call parameter_handler%update(key = ffrs_coarsening_fraction_key, value = 0.05_rp)
    call parameter_handler%update(key = ffrs_max_num_mesh_iterations_key, value = num_amr_steps )
    call refinement_strategy%create(error_estimator,parameter_handler%get_values())
  end subroutine setup_refinement_strategy

  !* #### Results plotting
  !* The handling the output writers is distributed in three different subroutines:
  !*
  !* - [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_initialize]] 
  !* - [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_write_current_amr_step]] 
  !* - [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_finalize]] 
  !* 
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_initialize]] 
  !* create the [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler]] 
  !* and add the select objects to be plotted.
  subroutine output_handler_initialize()
    if (write_postprocess_data) then
      call parameter_handler%update(output_handler_static_grid_key, value=.false.)
      call output_handler%create(parameter_handler%get_values())
      call output_handler%attach_fe_space(fe_space)
      call output_handler%add_fe_function(discrete_solution, 1, 'solution')
      call output_handler%add_cell_vector(error_estimator%get_sq_local_estimates(), 'cell_error_energy_norm_squared')
      call output_handler%open()
    end if   
  end subroutine output_handler_initialize

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_write_current_amr_step]] 
  !* is in charge of plotting every AMR step.
  subroutine output_handler_write_current_amr_step()
    if (write_postprocess_data) then
      call output_handler%append_time_step(real(current_amr_step,rp))
      call output_handler%write()
    end if  
  end subroutine output_handler_write_current_amr_step

  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:output_handler_finalize]] 
  !* finalize the writing procedure.
  subroutine output_handler_finalize()
    if (write_postprocess_data) then
      call output_handler%close()
    end if  
  end subroutine output_handler_finalize

  !* #### Environment finalization
  !* Subroutine [[tutorial_02_poisson_sharp_circular_wave_amr:free_all_objects]] 
  !* free all the objects created in the course of the program.
  subroutine free_all_objects()
    call output_handler%free()
    call error_estimator%free()
    call direct_solver%free()
    call fe_affine_operator%free()
    call discrete_solution%free()
    call fe_space%free()
    call strong_boundary_conditions%free()
    call source_term%free()
    call exact_solution%free()
    call triangulation%free()
    call environment%free()
    call world_context%free(.true.)  
  end subroutine free_all_objects
  
end program tutorial_02_poisson_sharp_circular_wave_amr
