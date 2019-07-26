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

program tutorial_03_poisson_sharp_circular_wave_parallel_amr
  use fempar_names
  use tutorial_03_discrete_integration_names
  use tutorial_03_functions_names
  use tutorial_03_error_estimator_names
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
  
  type(mpi_context_t)                         :: world_context
  type(environment_t)                         :: environment
  !* The triangulation_t object provides the mesh. In this case, we consider a serial p4est triangulation, i.e., 
  ! a triangulation that is distributed among processors, and that can be h-adapted and re-partitioned during the simulation.
  type(p4est_par_triangulation_t)             :: triangulation
  !* The fe_space_t is the global finite element space to be used.
  type(par_fe_space_t)                        :: fe_space
  type(strong_boundary_conditions_t)          :: strong_boundary_conditions 
  !* A scalar-valued function with the right hand side of the PDE problem at hand.
  type(sharp_circular_wave_source_term_t)     :: source_term
  !* A scalar-valued function with the exact (analytical) solution of the PDE problem at hand.
  type(sharp_circular_wave_solution_t)        :: exact_solution
  !* A fe_affine_operator_t that represents the affine operator the solution of which is the one we want, i.e., B = Ax-f.
  !* The solution is the root of this operator.
  type(fe_affine_operator_t)                  :: fe_affine_operator
  !* A fe_function_t belonging to the FE space defined above. Here, we will store the computed solution.
  type(fe_function_t)                         :: discrete_solution
  !* cg_discrete_integration_t provides the definition of the bilinear and linear forms for the problem at hand.
  type(cg_discrete_integration_t)             :: cg_discrete_integration
  !* direct_solver_t provides an interface to several external sparse direct solver packages (PARDISO, UMFPACK).
  type(iterative_linear_solver_t)             :: iterative_linear_solver
  !* The output handler type is used to generate the simulation data files for later visualization using, e.g., ParaView.
  type(output_handler_t)                      :: output_handler
  !* The following object is used to compute the (square) of the error energy norm for all cells, i.e., ||u-u_h||_E,K.
  type(poisson_error_estimator_t)             :: error_estimator
  type(fixed_fraction_refinement_strategy_t)  :: refinement_strategy
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
  
  
  !* Initialize FEMPAR library (i.e., construct system-wide variables)
  call fempar_init()
  call setup_parameter_handler()
  call setup_context_and_environment()
  call get_tutorial_cla_values()
  
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
    call redistribute_triangulation_and_fe_space()
    call setup_discrete_solution()
    call setup_and_assemble_fe_affine_operator()
    call solve_system()  
    call compute_error()
  end do
  call output_handler_finalize()
  
  call free_all_objects()
  !* Finalize FEMPAR library (i.e., destruct system-wide variables)
  call fempar_finalize()
contains

  subroutine setup_parameter_handler()
    !* In this subroutine we set up the system-wide variable referred to as "parameter_handler".
    !* This object is connected with the tutorial's command-line interface. It parses the arguments 
    !* provided to the command-line, and stores their values into a polymorphic dictionary of <key,value> pairs 
    !* (i.e., of type parameterlist_t); a pointer to such dictionary can be obtained calling "parameter_handler%get_values()", 
    !* as it will be required later on to create most of the objects provided by FEMPAR. FEMPAR users may define 
    !* additional command-line arguments by means of a user-provided subroutine. "define_tutorial_command_line_arguments" 
    !* subroutine (defined right below) is indeed the one that defines the command-line arguments for this particular 
    !* tutorial
    call parameter_handler%process_parameters(& 
         define_user_parameters_procedure = define_tutorial_clas, & 
         progname = tutorial_name, &
         version  = tutorial_version, &
         description = tutorial_description, &
         authors     = tutorial_authors)
  end subroutine setup_parameter_handler

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
  
  subroutine get_tutorial_cla_values()
    !* In this subroutine we use parameter_handler to obtain the values of the command-line arguments 
    !* which are particular to this tutorial
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

  subroutine setup_context_and_environment()
    !* Create a single-task group of tasks (as world_context is of type serial_context_t)
    call world_context%create()
    !* Force the environment to split world_context into a single group of 
    !* tasks (level) composed by a single task. As world_context is of type serial_context_t, 
    !* any other environment configuration is not possible
    !call parameter_handler%update(environment_num_levels_key, 1)
    !call parameter_handler%update(environment_num_tasks_x_level_key, [1])
    call environment%create(world_context, parameter_handler%get_values())
  end subroutine setup_context_and_environment
  
  subroutine setup_triangulation()
     integer(ip) :: i, num_dims
     if ( current_amr_step == 0 ) then
       !* Force the domain to be meshed by p4est triangulation to be the unit square/cube
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
     if ( write_postprocess_data ) then
       ! Re-adjust the size and contents of my_rank_cell_array to reflect the current status of the triangulation
       call setup_my_rank_cell_array()
     end if 
  end subroutine setup_triangulation
  
  subroutine flag_all_cells_for_refinement()
    class(cell_iterator_t), allocatable :: cell
    call triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )  
      call cell%set_for_refinement()
      call cell%next()
    end do
    call triangulation%free_cell_iterator(cell)
  end subroutine flag_all_cells_for_refinement
  
  subroutine setup_problem_functions()
      call source_term%create(triangulation%get_num_dims(),& 
                              alpha,circle_radius,circle_center)
      call exact_solution%create(triangulation%get_num_dims(),&
                                 alpha,circle_radius,circle_center)
  end subroutine setup_problem_functions
  
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
    else
       call fe_space%refine_and_coarsen()
    end if
    call fe_space%set_up_cell_integration()
 end subroutine setup_fe_space
 
 subroutine redistribute_triangulation_and_fe_space()
   call triangulation%redistribute()
   call fe_space%redistribute()
   call fe_space%set_up_cell_integration()
   if ( write_postprocess_data ) then
     ! Re-adjust the size and contents of my_rank_cell_array to reflect the current status of the triangulation
     call setup_my_rank_cell_array()
   end if
 end subroutine redistribute_triangulation_and_fe_space
  
 subroutine setup_discrete_solution()
   call discrete_solution%create(fe_space) 
   call fe_space%interpolate_dirichlet_values(discrete_solution)
 end subroutine setup_discrete_solution
  
 subroutine setup_and_assemble_fe_affine_operator()
   class(matrix_t)                  , pointer       :: matrix
   class(vector_t)                  , pointer       :: rhs
   type(par_sparse_matrix_t)        , pointer       :: sparse_matrix
 
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
   !* Now, we can compute the entries of the affine operator Ax-b by assembling the discrete weak form of the Poisson problem
   call fe_affine_operator%compute()
   rhs                => fe_affine_operator%get_translation()
   matrix             => fe_affine_operator%get_matrix()
   select type(matrix)
   class is (par_sparse_matrix_t)
     sparse_matrix => matrix
   end select
   
 end subroutine setup_and_assemble_fe_affine_operator
  
 subroutine solve_system()
   class(vector_t), pointer :: dof_values
   if ( current_amr_step == 0 ) then
      !* Iterative linear solver setup
      call iterative_linear_solver%create(environment)
      call iterative_linear_solver%set_type_from_pl(parameter_handler%get_values())
      call iterative_linear_solver%set_parameters_from_pl(parameter_handler%get_values())
      !* Next, we set the coefficient matrix and preconditioner to be used by iterative linear solver
      call iterative_linear_solver%set_operators(fe_affine_operator%get_matrix(), .identity. fe_affine_operator%get_matrix())
   else 
      call iterative_linear_solver%reallocate_after_remesh()
   end if
   !* We extract a pointer to the free nodal values of our FE function, which is the place in which we will store the result.
   dof_values => discrete_solution%get_free_dof_values()
   !* We solve the problem with the matrix already associated, the RHS obtained from the fe_affine_operator using get_translation, and 
   !* putting the result in dof_values.
   call iterative_linear_solver%solve(fe_affine_operator%get_translation(),dof_values)
   !* Once we have obtained the values corresponding to DoFs which are actually free, we have to update the values corresponding to DoFs
   !* which are hanging (constrained). This is required as, in this tutorial, we are using AMR on non-conforming meshes.
   call fe_space%update_hanging_dof_values(discrete_solution)
 end subroutine solve_system
 
 subroutine compute_error()
   real(rp) :: global_error_energy_norm
   if ( current_amr_step == 0 ) then
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
        if (current_amr_step == num_amr_steps) then
          write(*,'(a70)')  repeat('=', 70)
        end if
      end if 
   end if 
  end subroutine compute_error
  
  subroutine setup_refinement_strategy()
    real(rp) :: ref_fraction, coarse_fraction
    if ( triangulation%get_num_dims() == 2 ) then
      ref_fraction    = 0.10_rp
      coarse_fraction = 0.05_rp
    else 
      ref_fraction    = 0.10_rp
      coarse_fraction = 0.05_rp
    end if  
    call parameter_handler%update(key = ffrs_refinement_fraction_key, value = ref_fraction)
    call parameter_handler%update(key = ffrs_coarsening_fraction_key, value = coarse_fraction)
    call parameter_handler%update(key = ffrs_max_num_mesh_iterations_key, value = num_amr_steps )
    call refinement_strategy%create(error_estimator,parameter_handler%get_values())
  end subroutine setup_refinement_strategy
  
  subroutine setup_my_rank_cell_array()
    if ( environment%am_i_l1_task() ) then
      call my_rank_cell_array%resize(triangulation%get_num_local_cells())
      call my_rank_cell_array%init(real(environment%get_l1_rank(),rp))
    end if 
  end subroutine setup_my_rank_cell_array
    
  subroutine output_handler_initialize()
    if (write_postprocess_data) then
      call parameter_handler%update(output_handler_static_grid_key, value=.false.)
      call output_handler%create(parameter_handler%get_values())
      call output_handler%attach_fe_space(fe_space)
      call output_handler%add_fe_function(discrete_solution, 1, 'solution')
      call output_handler%add_cell_vector(error_estimator%get_sq_local_estimates(), 'cell_error_energy_norm_squared')
      call setup_my_rank_cell_array() 
      call output_handler%add_cell_vector(my_rank_cell_array, 'subdomain_partition')
      call output_handler%open()
    end if   
  end subroutine output_handler_initialize
  
  subroutine output_handler_write_current_amr_step()
    if (write_postprocess_data) then
      call output_handler%append_time_step(real(current_amr_step,rp))
      call setup_my_rank_cell_array() ! Re-adjust the size and contents of my_rank_cell_array to reflect the current status of the triangulation
      call output_handler%write()
    end if  
  end subroutine output_handler_write_current_amr_step
  
  subroutine output_handler_finalize()
    if (write_postprocess_data) then
      call output_handler%close()
    end if  
  end subroutine output_handler_finalize
  
  subroutine free_all_objects()
    !* Free all the objects created in the course of the program
    call output_handler%free()
    call error_estimator%free()
    call iterative_linear_solver%free()
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
  
end program tutorial_03_poisson_sharp_circular_wave_parallel_amr
