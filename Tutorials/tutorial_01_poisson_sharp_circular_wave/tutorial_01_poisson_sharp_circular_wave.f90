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
!* The [[tutorial_01_poisson_sharp_circular_wave]],
program tutorial_01_poisson_sharp_circular_wave
!* uses the [[fempar_names]] to import all **`FEMPAR`** library symbols (i.e., derived types, parameter constants, system-wide variables, etc.),
!* and a set of tutorial-specific module units, 
!*
!* - [[tutorial_01_discrete_integration_names]]
!* - [[tutorial_01_functions_names]]
!* - [[tutorial_01_error_estimator_names]]
!*
!* These modules are not part of the **`FEMPAR`** library, but developed specifically for the problem at hand. Each of these modules defines a tutorial-specific data type and its TBPs.
  use fempar_names
  use tutorial_01_discrete_integration_names
  use tutorial_01_functions_names
  use tutorial_01_error_estimator_names

!* ### Variable definition
!* This tutorial declares:
!*
!* - a set of parameter constants, typically the tutorial name, authors, problem description, etc. (to be output on screen on demand by the user),
!* - the tutorial data type instances in charge of the FE simulation, such as the triangulation (mesh) of the computational domain or the FE space from which the approximate solution of the PDE is sought, 
!* - and a set of variables to hold the values of the Command-Line-Arguments (CLAs) which are specific to the tutorial.
  implicit none
  character(*), parameter :: tutorial_name = "tutorial_01_poisson_sharp_circular_wave"
  character(*), parameter :: tutorial_version = "v1.0.0"
  character(*), parameter :: tutorial_description = "Solves a 2D/3D Poisson problem on the unit &
                                                     square/cube, where the solution has a sharp circular & 
                                                     wave front parametrizable via command-line arguments. &
                                                     The discretization of the domain is performed using a static &
                                                     triangulation, i.e., a mesh which is generated at the beginning & 
                                                     and cannot be h-adapted during the simulation (e.g., by means of &
                                                     a-posteriori error estimators.)"
  character(*), parameter :: tutorial_authors = "Santiago Badia and Alberto F. Martín"
  

  !* The [[tutorial_01_poisson_sharp_circular_wave:world_context]] is a software abstraction 
  !* for a group of parallel tasks (a.k.a. processes) and the communication layer that orchestrates their concurrent execution.
  !* In this case is declared of type [[serial context t]] as it is designed to work in serial computing environments.
  type(serial_context_t)                  :: world_context
  !* the [[tutorial_01_poisson_sharp_circular_wave:environment]], of type [[serial_context_t]], organizes the tasks of the context
  !* from which it is set up into subgroups of tasks, referred to as levels, and builds up additional communication mechanisms 
  !* to communicate tasks belonging to different levels.
  type(environment_t)                     :: environment
  !* The [[tutorial_01_poisson_sharp_circular_wave:triangulation]], of type [[serial_tringulation_t]]. object provides the mesh. 
  !* In this case, we consider a serial static triangulation, i.e., 
  !* a triangulation that is not distributed among processors, and that cannot be h-adapted during the simulation.
  type(serial_triangulation_t)            :: triangulation
  !* The [[tutorial_01_poisson_sharp_circular_wave:fe_space]], of type [[serial_fe_space_t]] is the global finite element space to be used.
  type(serial_fe_space_t)                 :: fe_space
  !* The [[tutorial_01_poisson_sharp_circular_wave:strong_boundary_conditions]] is the object that provides the strong boundary conditions definition.
  type(strong_boundary_conditions_t)      :: strong_boundary_conditions 
  !* The [[tutorial_01_poisson_sharp_circular_wave:source_term]] is a scalar-valued function with the right hand side of the PDE problem at hand.
  type(sharp_circular_wave_source_term_t) :: source_term
  !* The [[tutorial_01_poisson_sharp_circular_wave:exact_solution]] ia a scalar-valued function with the exact (analytical) solution of the PDE problem at hand.
  type(sharp_circular_wave_solution_t)    :: exact_solution
  !* [[tutorial_01_poisson_sharp_circular_wave:fe_affine_operator]], of type [[fe_affine_operator_t]],
  !* represents the affine operator the solution of which is the one we want, i.e., \(B = Ax-f\).
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)              :: fe_affine_operator
  !* The [[tutorial_01_poisson_sharp_circular_wave:discrete_solution]], of type [[fe_function_t]], 
  !* is the object belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                     :: discrete_solution
  !* [[tutorial_01_discrete_integration_names:poisson_discrete_integration_t]] provides the definition of the bilinear and linear forms for the problem at hand.
  type(cg_discrete_integration_t), target :: cg_discrete_integration
  type(dg_discrete_integration_t), target :: dg_discrete_integration
  !* The [[tutorial_01_poisson_sharp_circular_wave:direct_solver]], of type [[direct_solver_t]], provides an interface to several external sparse direct solver packages (PARDISO, UMFPACK).
  type(direct_solver_t)                   :: direct_solver
  !* The [[tutorial_01_poisson_sharp_circular_wave:output_handler]] object, of type [[output_handler_t]], is used to generate the simulation data files for later visualization using, e.g., ParaView.
  type(output_handler_t)                  :: output_handler
  !* The [[tutorial_01_poisson_sharp_circular_wave:error_estimator]] object is used to compute the (square) of the error energy norm for all cells, i.e., \(||u-u_h||_E,K\).
  type(poisson_error_estimator_t)         :: error_estimator
  
  !* Variables to store the values of the command-line arguments which are particular to this tutorial.
  !* Execute the tutorial with the "-h" command-line flag for more details
  real(rp)                  :: alpha
  real(rp)                  :: circle_radius
  real(rp), allocatable     :: circle_center(:)
  character(:), allocatable :: fe_formulation
  logical                   :: write_postprocess_data

  !* ### Main code  
  !* Here you can see the main code of this tutorial representing the entire flow
  !* of actions/computations performed by this tutorial.
  !*
  !* As there is almost a one-to-one mapping among the data type instances and the
  !* helper procedures called in this tutorial, we will introduce them in the sequel step-by-step along
  !* with code snippets of the corresponding [helper procedures](#helper-procedures).
  call fempar_init()
  call setup_parameter_handler()
  call get_tutorial_cla_values()
  call setup_context_and_environment()
  call setup_triangulation()
  call setup_problem_functions()
  call setup_strong_boundary_conditions()
  call setup_fe_space()
  call setup_discrete_solution()
  call setup_and_assemble_fe_affine_operator()
  call solve_system()  
  call compute_error()
  call write_postprocess_data_files()
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
  !* In [[tutorial_01_poisson_sharp_circular_wave:setup_parameter_handler]] subroutine 
  !* we set up the system-wide variable referred to as [[fempar_parameter_handler_names:parameter_handler]].
  !* 
  !* **`FEMPAR`** users may define additional command-line arguments by means of a user-provided subroutine. 
  !* [[tutorial_01_poisson_sharp_circular_wave:define_tutorial_clas]] 
  !* subroutine (defined right below) is indeed the one that defines the command-line arguments for this particular 
  !* tutorial
  subroutine setup_parameter_handler()
    !> Define tutorial metadata and process CLI parameters
    call parameter_handler%process_parameters(& 
         define_user_parameters_procedure = define_tutorial_clas, & 
         progname = tutorial_name, &
         version  = tutorial_version, &
         description = tutorial_description, &
         authors     = tutorial_authors)
  end subroutine setup_parameter_handler

  !* In [[tutorial_01_poisson_sharp_circular_wave:define_tutorial_clas]] subroutine 
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
    !> Define tutorial CLAs
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
  end subroutine define_tutorial_clas
  
  !* In [[tutorial_01_poisson_sharp_circular_wave:get_tutorial_cla_values]] subroutine 
  !* we use the [[fempar_parameter_handler_names:parameter_handler]] to obtain the values of the command-line arguments 
  !* which are particular to this tutorial.
  subroutine get_tutorial_cla_values()
    !> Obtain CLA values from [[fempar_parameter_handler_names:parameter_handler]]
    call parameter_handler%getasstring("FE_FORMULATION", fe_formulation)
    call parameter_handler%get("ALPHA", alpha)
    call parameter_handler%get("CIRCLE_RADIUS", circle_radius)
    call parameter_handler%getasarray("CIRCLE_CENTER", circle_center)
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
    write(*,'(a54)')  repeat('=', 54)
  end subroutine get_tutorial_cla_values

  !* #### Environment initialization
  !* In [[tutorial_01_poisson_sharp_circular_wave:setup_context_and_environment]] subroutine 
  !* a single-task group of tasks is created.
  !* 
  !* Environment related parameters are also updated to preserve the coherence with this context
  !* (as [[tutorial_01_poisson_sharp_circular_wave:world_context]] is of type [[serial_context_t]]).
  !* 
  !* When parameters are already update, the [[tutorial_01_poisson_sharp_circular_wave:environment]] is created.
  !*
  !* @note
  !* As `world_context` is of type [[serial_context_t]], 
  !* any other environment configuration is not possible
  !* @endnote
  subroutine setup_context_and_environment()
    !> Create the context and the environment.
    call world_context%create()
    call parameter_handler%update(environment_num_levels_key, 1)
    call parameter_handler%update(environment_num_tasks_x_level_key, [1])
    call environment%create(world_context, parameter_handler%get_values())
  end subroutine setup_context_and_environment

  !* #### Tringulation
  !* In [[tutorial_01_poisson_sharp_circular_wave:setup_triangulation]] subroutine 
  !* the triangulation (of type [[triangulation_t]]) is created.  
  !* The structure hex mesh generator creates a domain that has:
  !*
  !*   - 8 boundary objects (corners+edges) in 2D
  !*   - 26 boundary objects (face+corners+edges) in 3D
  !*   - 1 interior object  (i.e., the cell interior).
  !*
  !* A graphical example for a 2D case is presented next
  !* ```text
  !* |          6
  !* |     3 - - - - 4
  !* |     -         -         1,2,3,4,5,6,7,8 are id for the boundary of the domain
  !* |   7 -    9    - 8
  !* |     -         -         9 is the id for the interior
  !* |     1 - - - - 2
  !* |          5     
  !* ```
  !* This configuration, the hexaedral uniform unit square/cube mesh (to be meshed by [[uniform_hex_mesh_t]]), 
  !* is forced by means of the [[fempar_parameter_handler_names:parameter_handler]].
  subroutine setup_triangulation()
    !> Create the triangulation
     integer(ip) :: num_dims
     call parameter_handler%update(static_triang_generate_from_key, &
                                   static_triang_generate_from_struct_hex_mesh_generator)
     
     call parameter_handler%get(struct_hex_mesh_generator_num_dims_key, num_dims)
     if ( num_dims == 2 ) then
       call parameter_handler%update(struct_hex_mesh_generator_domain_limits_key, &
                                      [0.0_rp,1.0_rp,0.0_rp,1.0_rp])
     else
       call parameter_handler%update(struct_hex_mesh_generator_domain_limits_key, &
                                     [0.0_rp,1.0_rp,0.0_rp,1.0_rp,0.0,1.0_rp])
     end if

     call triangulation%create(environment, parameter_handler%get_values())
  end subroutine setup_triangulation

  !* #### Analytical expressions
  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_problem_functions]] 
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
     !> Create the source term and the the exact solution.
     call source_term%create(triangulation%get_num_dims(),& 
                             alpha,circle_radius,circle_center)
     call exact_solution%create(triangulation%get_num_dims(),&
                             alpha,circle_radius,circle_center)
  end subroutine setup_problem_functions


  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_strong_boundary_conditions]] 
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
  !* [[tutorial_01_poisson_sharp_circular_wave:setup_strong_boundary_conditions]] come from.) 
  !*
  !* With the aforementioned in mind, subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_strong_boundary_conditions]]
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
  !* In subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_fe_space]], 
  !* this tutorial sets up the fe space instance, i.e., the computer representation of \(V_h\). 
  !* 
  !* The subroutine forces fe space to build a single-field FE space 
  !* that it is built from the same local FE for all cells \(K \in \mathcal{T}_h\). 
  !* In particular, from a scalar-valued Lagrangian-type FE, as required by our model problem. 
  !* The parameter value corresponding to the polynomial order of \(V_h|_K\) is not forced, and thus can be selected from
  !* the CLI. 
  !*
  !* Besides, [[tutorial_01_poisson_sharp_circular_wave:setup_strong_boundary_conditions]] also forces the construction 
  !* of either a conforming or non-conforming \(V_h\), depending on the FE formulation selected by the user. 
  !* In the latter case, fe space does not require strong boundary conditions, as there
  !* are not BCs to be strongly enforced in this case. The call in these lines generates a global
  !* numbering of the DOFs in \(V_h\), and (if it applies) identifies the DOFs sitting on the regions of
  !* the boundary of the domain which are subject to strong BCs (combining the data provided by
  !* the triangulation and strong boundary conditions). 
  !*
  !* Finally, [[tutorial_01_poisson_sharp_circular_wave:setup_strong_boundary_conditions]]
  !* activates the internal data structures that fe space provides for the numerical evaluation of cell
  !* and facet integrals, resp. These are required later on to evaluate the discrete bilinear and linear
  !* forms of the FE formulations implemented by this tutorial. We note that the last `call` is
  !* only required for the IP `DG` formulation, as in the `CG` formulation there are not facet integrals
  !* to be evaluated.
  subroutine setup_fe_space()
    type(string) :: fes_field_types(1), fes_ref_fe_types(1)
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
  end subroutine setup_fe_space

  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_discrete_solution]]
  !* sets up the discrete solution object, of type [[fe_function_t]]. This FEMPAR
  !* data type represents an element of \(V_h\) , the FE solution \(u_h\) in the case of this tutorial. 
  !*
  !* This subroutine allocates room for storing the DOF values of \(u_h\) , and then it
  !* computes the values of those DOFs of \(u_h\) which lay on a region of the domain boundary which
  !* is subject to strong BCs, the whole boundary of the domain in the case of this tutorial. 
  !* 
  !* The [[serial_fe_space_t:interpolate_dirichlet_values]] TBP of [[tutorial_01_poisson_sharp_circular_wave:fe_space]] 
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
  
  !* In subroutine [[tutorial_01_poisson_sharp_circular_wave:setup_and_assemble_fe_affine_operator]], 
  !* this tutorial program builds the fe affine operator instance, i.e., the computer
  !* representation of the previously operator defined. Building this instance is a two-step process. 
  !* 
  !* First, we have to call the [[fe_affine_operator_t:create]] TBP. Apart from specifying the data storage format and
  !* properties of the stiffness matrix \(A\), this call equips [[fe_affine_operator_t]] with all that it re-
  !* quires in order to evaluate the entries of the operator. 
  !* 
  !* Second, we have to call the [[fe_affine_operator_t:compute]]
  !* TBP , that triggers the actual numerical evaluation and assembly of the discrete
  !* weak form. With extensibility and flexibility in mind, this latter responsibility does not fall on
  !* [[fe_affine_operator_t]], but actually on the data type extensions of a key abstract class defined
  !* within *`FEMPAR`*, referred to as [[discrete_integration_t]]. 
  !* 
  !* This tutorial 01 implements two type extensions of discrete integration t, 
  !* namely [[tutorial_01_discrete_integration_names:cg_discrete_integration_t]],
  !* and [[tutorial_01_discrete_integration_names:dg_discrete_integration_t]]. 
  !* 
  !* First lines of this procedure passes to these instances the source term and boundary function to be
  !* imposed on the Dirichlet boundary. We note that the _IP DG_ formulation works directly with
  !* the analytical expression of the boundary function (i.e., exact solution), while the _CG_ formulation
  !* with its interpolation (i.e., discrete solution). 
  subroutine setup_and_assemble_fe_affine_operator()
    class(discrete_integration_t), pointer :: discrete_integration
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
    call fe_affine_operator%compute()
  end subroutine setup_and_assemble_fe_affine_operator

  !* #### System solving
  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:solve_system]], 
  !* finds the root of the FE operator, i.e., \(u_h \in V_h\) such that \(F_h(u_h) = 0\). For such
  !* purpose, it first sets up the [[tutorial_01_poisson_sharp_circular_wave:direct_solver]] instance 
  !* and later uses its [[direct_solver_t:solve]] TBP. This latter method is fed with the “translation” of \(F_h\) , 
  !* i.e., the right hand side \(f\) of the linear system, and the free DOF values of \(u_h\) as the unknown to be found. 
  !* The free DOFs of \(u_h\) are those whose values are not constrained, e.g., by strong BCs. direct solver t is a
  !* **`FEMPAR`**data type that offers interfaces to (non-distributed) sparse direct solvers provided by
  !* external libraries. 
  !* 
  !* This subroutine does not force any parameter value related to [[direct_solver_t]]; the default CLA values are
  !* appropriate for the Poisson problem. In any case, at this point the reader is encouraged to
  !* inspect the CLAs linked to direct solver t and play around them.
  !*
  !* @note
  !* In its first public release,**`FEMPAR`**provides interfaces to PARDISO (the
  !* version available in the Intel MKL library) and UMFPACK, although it is designed such that
  !* additional sparse direct solver implementations can be easily added.
  !* @endnote
  subroutine solve_system()
    class(vector_t), pointer :: dof_values
    call direct_solver%set_type_from_pl(parameter_handler%get_values())
    call direct_solver%set_parameters_from_pl(parameter_handler%get_values())
    call direct_solver%set_matrix(fe_affine_operator%get_matrix())
    dof_values => discrete_solution%get_free_dof_values()
    call direct_solver%solve(fe_affine_operator%get_translation(),dof_values) 
  end subroutine solve_system

  !* #### Error computation
  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:compute_error]] 
  !* computes \(e^2_k\) for each \(K \in \mathcal{T}_h\), and the global error \(e\);
  !* 
  !* This subroutine relies on the [[tutorial_01_poisson_sharp_circular_wave:error_estimator]] instance.
  !* In particular, the actual computation of \(e^2_k\) occurs at
  !* the call to the [[error_estimator_t:compute_local_true_errors]] TBP of this data type.
  !*
  !* @note
  !* At this point, the user is encouraged to inspect the implementation of this data type in order to
  !* grasp how the numerical evaluation of the integrals required for the computation of \(e^2_k\)
  !* is carried out using **`FEMPAR`**.
  !* @endnote
  subroutine compute_error()
    real(rp) :: global_error_energy_norm
    call error_estimator%create(fe_space,parameter_handler%get_values())
    call error_estimator%set_exact_solution(exact_solution)
    call error_estimator%set_discrete_solution(discrete_solution)
    call error_estimator%compute_local_true_errors()
    global_error_energy_norm = sqrt(sum(error_estimator%get_sq_local_true_error_entries()))
    write(*,'(a54)')  repeat('=', 54)
    write(*,'(a54)') tutorial_name // ' results'
    write(*,'(a54)')  repeat('=', 54)
    write(*,'(a30,i24)') 'NUM_CELLS:' // repeat(' ', 80), triangulation%get_num_cells()
    write(*,'(a30,i24)') 'NUM_DOFS:' // repeat(' ', 80), fe_space%get_total_num_dofs()
    write(*,'(a30,e24.10)') 'GLOBAL ERROR ENERGY NORM:'// repeat(' ', 80), global_error_energy_norm
    write(*,'(a54)')  repeat('=', 54)
  end subroutine compute_error

  !* #### Results plotting
  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:write_postprocess_data_files]] 
  !* write estimated error quantities along with the FE solution \(u_h\) to
  !* output data files for later visualization. The [[tutorial_01_poisson_sharp_circular_wave:output_handler]]
  !* is in charge of output data files generation.
  subroutine write_postprocess_data_files()
    if ( write_postprocess_data ) then
      call output_handler%create(parameter_handler%get_values())
      call output_handler%attach_fe_space(fe_space)
      ! call fe_space%interpolate(1, exact_solution, discrete_solution)
      call output_handler%add_fe_function(discrete_solution, 1, 'solution')
      call output_handler%add_cell_vector(error_estimator%get_sq_local_true_errors(),'cell_error_energy_norm_squared')
      call output_handler%open()
      call output_handler%write()
      call output_handler%close()
    end if  
  end subroutine write_postprocess_data_files

  !* #### Environment finalization
  !* Subroutine [[tutorial_01_poisson_sharp_circular_wave:free_all_objects]] 
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
  
end program tutorial_01_poisson_sharp_circular_wave
