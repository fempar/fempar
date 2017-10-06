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
module par_test_hts_driver_names
  use fempar_names
  use par_test_hts_params_names
  use hts_discrete_integration_names
  use hts_conditions_names
		use maxwell_conditions_names 
  use hts_analytical_functions_names
		use maxwell_analytical_functions_names 
		use hts_theta_method_names
		use hts_nonlinear_solver_names 
# include "debug.i90"

  implicit none
  private

		integer(ip), parameter :: hts = 1
  integer(ip), parameter :: air = 2
		
  type par_test_hts_fe_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_hts_params_t)      :: test_params
     type(ParameterList_t), pointer   :: parameter_list
     
     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                          :: fe_space 
     type(p_reference_fe_t), allocatable           :: reference_fes(:) 
     class(l1_coarse_fe_handler_t), pointer        :: coarse_fe_handler 
	    type(p_l1_coarse_fe_handler_t), allocatable   :: coarse_fe_handlers(:)
	    type(hts_CG_discrete_integration_t)           :: hts_integration
					
					! Type of problem 
					type(maxwell_conditions_t)                    :: maxwell_conditions 
					type(maxwell_analytical_functions_t)          :: maxwell_analytical_functions 
     type(hts_conditions_t)                        :: hts_conditions
     type(hts_analytical_functions_t)              :: hts_analytical_functions

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)            :: fe_affine_operator
    
     ! Iterative linear solvers data type
					type(hts_nonlinear_solver_t)              :: nonlinear_solver 
     type(iterative_linear_solver_t)           :: iterative_linear_solver
					type(mlbddc_t)                            :: mlbddc
 
     ! maxwell problem solution FE function
     type(fe_function_t)                   :: H_current
					type(fe_function_t)                   :: H_previous 
     
     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                    :: environment
					
					! Time integration 
					type(theta_method_t)                   :: theta_method 
				
					! Output handler 
					real(rp), allocatable                  :: set_id_cell_vector(:)
					type(output_handler_t)                 :: oh	
				

   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
	    procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
					procedure        , private :: setup_theta_method 
					procedure        , private :: setup_nonlinear_solver 
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
					procedure        , private :: setup_fe_coarse_space 
     procedure        , private :: solve_system
					procedure        , private :: solve_nonlinear_system 
     procedure        , private :: check_solution
     procedure        , private :: initialize_output
     procedure        , private :: write_time_step_solution 
					procedure        , private :: finalize_output 
					procedure                  :: run_simulation
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
  end type par_test_hts_fe_driver_t

  ! Types
  public :: par_test_hts_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

!========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%environment%create (this%parameter_list)
  end subroutine setup_environment
   
  subroutine setup_triangulation(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(vef_iterator_t), allocatable :: vef
				real(rp)                           :: domain(6)
				integer(ip)                        :: istat 
				
				! Create a structured mesh with a custom domain 
				domain   = this%test_params%get_domain_limits()
    istat = this%parameter_list%set(key = hex_mesh_domain_limits_key , value = domain); check(istat==0)
				
    call this%triangulation%create(this%parameter_list, this%environment)
				call assign_fes_id() 
				
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       call this%triangulation%create_vef_iterator(vef)
       do while ( .not. vef%has_finished() )
          if(vef%is_at_boundary()) then 
             call vef%set_set_id(1)
          else
             call vef%set_set_id(0)
          end if
          call vef%next()
       end do
       call this%triangulation%free_vef_iterator(vef)
    end if  
	
    call this%triangulation%setup_coarse_triangulation()
				
				contains 
				subroutine assign_fes_id() 
				implicit none 
				! Locals 
    integer(ip)                 , allocatable :: cells_set(:) 
    class(cell_iterator_t)      , allocatable :: cell
    type(point_t), allocatable                :: cell_coordinates(:)
    integer(ip)                               :: inode
				integer(ip)       :: l1_rank 
    integer(ip)       :: icell, icoord 
    real(rp)          :: cx, cy, cz 
    integer(ip)       :: istat 
	   real(rp)          :: R, h, x0, y0, z0

				real(rp)          :: hts_size(0:SPACE_DIM-1)  
				real(rp)          :: domain(6)
				real(rp)          :: hts_lx, hts_ly, hts_lz
				real(rp)          :: lx, ly, lz 
				
				if ( this%environment%am_i_l1_task() ) then 
				l1_rank = this%environment%get_l1_rank() 
							
				if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
	
					hts_size = this%test_params%get_hts_domain_length() 
					hts_lx = hts_size(0) 
					hts_ly = hts_size(1) 
					hts_lz = hts_size(2) 
					
					domain = this%test_params%get_domain_limits()
					lx = domain(2)-domain(1) 
					ly = domain(4)-domain(3) 
					lz = domain(6)-domain(5) 
					
    ! Assign subset_id to different cells for the created structured mesh 
    allocate(cells_set(this%triangulation%get_num_local_cells() ), stat=istat); check(istat==0)
    call this%triangulation%create_cell_iterator(cell)
    allocate(cell_coordinates( cell%get_num_nodes() ), stat=istat); check(istat==0) 
    
    do while ( .not. cell%has_finished() )
				   if ( cell%is_local() ) then 
       call cell%get_nodes_coordinates(cell_coordinates)
       ! Compute center of the element coordinates 
       cx = 0.0_rp
       cy = 0.0_rp 
	      cz = 0.0_rp 
       do inode=1,cell%get_num_nodes()  
          cx = cx + cell_coordinates(inode)%get(1)
          cy = cy + cell_coordinates(inode)%get(2)
		        cz = cz + cell_coordinates(inode)%get(3)
       end do
       cx = cx/real(cell%get_num_nodes(),rp)
       cy = cy/real(cell%get_num_nodes(),rp)
	      cz = cz/real(cell%get_num_nodes(),rp)

       ! Select material case: HTS TAPE in the center 
       if ( ( ((lx-hts_lx)/2.0_rp<=cx) .and. (cx<=(lx+hts_lx)/2.0_rp) ) .and. &
											 ( ((ly-hts_ly)/2.0_rp<=cy) .and. (cy<=(ly+hts_ly)/2.0_rp) ) .and. & 
												( ((lz-hts_lz)/2.0_rp<=cz) .and. (cz<=(lz+hts_lz)/2.0_rp) )) then
          cells_set( cell%get_gid() ) = hts 
       else 
          cells_set( cell%get_gid() ) = air
       end if
	   
							end if 
				   call cell%next() 
    end do
				
    call this%triangulation%fill_cells_set(cells_set)  
    deallocate(cells_set, stat=istat); check(istat==0) 
    deallocate(cell_coordinates, stat=istat); check(istat==0) 
    call this%triangulation%free_cell_iterator(cell)

	     end if 
				end if 
				
				end subroutine 
  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable    :: cell
    class(reference_fe_t), pointer         :: reference_fe_geo
    
    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)
    
    if ( this%environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      reference_fe_geo => cell%get_reference_fe()
      this%reference_fes(1) =  make_reference_fe ( topology   = reference_fe_geo%get_topology(),           &
                                                   fe_type    = fe_type_nedelec,                           &
                                                   num_dims   = this%triangulation%get_num_dims(),         &
                                                   order      = this%test_params%get_reference_fe_order(), &
                                                   field_type = field_type_vector,                         &
                                                   conformity = .true. )
      call this%triangulation%free_cell_iterator(cell)
    end if  
		 
  end subroutine setup_reference_fes
  
  subroutine setup_theta_method(this) 
  implicit none 
  class(par_test_hts_fe_driver_t), intent(inout) :: this
  
    call this%theta_method%create( this%environment,                            & 
																																		 this%test_params%get_theta_value(),          &               
                                   this%test_params%get_initial_time(),         &
                                   this%test_params%get_final_time(),           & 
                                   this%test_params%get_num_time_steps(),       &
                                   this%test_params%get_max_time_step(),        & 
                                   this%test_params%get_min_time_step(),        &
                                   this%test_params%get_save_solution_n_steps() )
  
  end subroutine setup_theta_method 

  subroutine setup_nonlinear_solver (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    call this%nonlinear_solver%create( convergence_criteria = this%test_params%get_nonlinear_convergence_criteria() , &
                                       abs_tol = this%test_params%get_absolute_nonlinear_tolerance(),                 &
                                       rel_tol = this%test_params%get_relative_nonlinear_tolerance(),                 &
                                       max_iters = this%test_params%get_max_nonlinear_iterations(),                   &
                                       ideal_iters = this%test_params%get_stepping_parameter(),                       &
                                       fe_affine_operator = this%fe_affine_operator,                                  &
                                       current_dof_values = this%H_current%get_free_dof_values()                      )
    
    
  end subroutine setup_nonlinear_solver
		
  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
   	class(cell_iterator_t), allocatable        :: cell
    class(reference_fe_t), pointer             :: reference_fe_geo
	integer(ip) :: istat 
	
	if ( this%environment%am_i_l1_task() ) then
	call this%triangulation%create_cell_iterator(cell)
    reference_fe_geo => cell%get_reference_fe()
	call this%triangulation%free_cell_iterator(cell)
	
	 if (this%triangulation%get_num_dims() == 3 ) then 
	   if ( reference_fe_geo%get_topology() == topology_tet ) then
	   allocate ( tet_Hcurl_l1_coarse_fe_handler_t :: this%coarse_fe_handler )
	   elseif ( reference_fe_geo%get_topology() == topology_hex ) then 
	   allocate ( hex_Hcurl_l1_coarse_fe_handler_t :: this%coarse_fe_handler )
	   end if 
	 else 
	  allocate ( standard_l1_coarse_fe_handler_t :: this%coarse_fe_handler )
	 end if 
	 end if 
	 
    allocate(this%coarse_fe_handlers(1), stat=istat); check(istat==0)
    this%coarse_fe_handlers(1)%p => this%coarse_fe_handler
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
		
	 if (this%test_params%get_is_analytical_solution() ) then 
		call this%maxwell_conditions%set_num_dims(this%triangulation%get_num_dims())
	
		! Set-up Dirichlet boundary conditions  
 call this%maxwell_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
	call this%maxwell_conditions%set_boundary_function_Hx(this%maxwell_analytical_functions%get_boundary_function_Hx())
	call this%maxwell_conditions%set_boundary_function_Hy(this%maxwell_analytical_functions%get_boundary_function_Hy())
	if ( this%triangulation%get_num_dims() == 3) then 
	call this%maxwell_conditions%set_boundary_function_Hz(this%maxwell_analytical_functions%get_boundary_function_Hz())
	end if 
	
				! Create FE SPACE 
    call this%fe_space%create( triangulation       = this%triangulation,      &
                               reference_fes       = this%reference_fes,      &
                               coarse_fe_handlers  = this%coarse_fe_handlers, & 
							                        conditions          = this%maxwell_conditions  )			
				
	elseif ( .not. this%test_params%get_is_analytical_solution() ) then 		
	call this%hts_conditions%set_num_dims(this%triangulation%get_num_dims())
	
		! Set-up Dirichlet boundary conditions  
 call this%hts_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
	call this%hts_conditions%set_boundary_function_Hx(this%hts_analytical_functions%get_boundary_function_Hx())
	call this%hts_conditions%set_boundary_function_Hy(this%hts_analytical_functions%get_boundary_function_Hy())
	if ( this%triangulation%get_num_dims() == 3) then 
	call this%hts_conditions%set_boundary_function_Hz(this%hts_analytical_functions%get_boundary_function_Hz())
	
	    call this%hts_analytical_functions%get_parameter_values( H  = this%test_params%get_external_magnetic_field_amplitude(),  &
                                                              wH = this%test_params%get_external_magnetic_field_frequency()   ) 
	end if 
	
				! Create FE SPACE 
    call this%fe_space%create( triangulation       = this%triangulation,      &
                               reference_fes       = this%reference_fes,      &
                               coarse_fe_handlers  = this%coarse_fe_handlers, & 
							                        conditions          = this%hts_conditions  )
		end if 
    
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration()   
									
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
				! Need to initialize dof_values, no interpolation available in Nedelec
    class(vector_t) , pointer :: dof_values_current 
    class(vector_t) , pointer :: dof_values_previous
        	
				if ( this%test_params%get_is_analytical_solution() ) then 
				call this%hts_integration%set_analytical_functions(this%maxwell_analytical_functions%get_source_term())
				else
    call this%hts_integration%set_analytical_functions(this%hts_analytical_functions%get_source_term())
				end if 
				
				call this%hts_integration%set_fe_functions(this%H_previous, this%H_current)
				call this%hts_integration%set_theta_method(this%theta_method) 
				call this%hts_integration%set_parameter_values( this%test_params )
    
    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .true. ], &
                                          diagonal_blocks_symmetric         = [ .true. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                          fe_space                          = this%fe_space, &
                                          discrete_integration              = this%hts_integration )
		
	! Initialize previous time step fe_functions 
	call this%H_previous%create(this%fe_space)
 dof_values_previous => this%H_previous%get_free_dof_values()  
 call dof_values_previous%init(0.0_rp) 
	
	! Current time step fe_function 
	 call this%H_current%create(this%fe_space)
		call this%fe_space%project_dirichlet_values_curl_conforming(this%H_current, time=this%theta_method%get_current_time())
		dof_values_current => this%H_current%get_free_dof_values()
	 call dof_values_current%init(0.0_rp)
		     						
  end subroutine setup_system
  
		subroutine setup_fe_coarse_space(this) 
				implicit none 
		class(par_test_hts_fe_driver_t), intent(inout) :: this

			select type ( ch=> this%coarse_fe_handler ) 
	  class is ( Hcurl_l1_coarse_fe_handler_t ) 
	  call this%coarse_fe_handler%setup_tools( this%fe_space )
	 end select 
	
		call this%fe_space%setup_coarse_fe_space(this%parameter_list)
	end subroutine setup_fe_coarse_space
		
  subroutine setup_solver (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
	type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    integer(ip) :: FPLError
    integer(ip) :: ilev
    integer(ip) :: iparm(64)
	
    ! See https://software.intel.com/en-us/node/470298 for details
    iparm      = 0 ! Init all entries to zero
    iparm(1)   = 1 ! no solver default
    iparm(2)   = 2 ! fill-in reordering from METIS
    iparm(8)   = 2 ! numbers of iterative refinement steps
    iparm(10)  = 8 ! perturb the pivot elements with 1E-8
    iparm(11)  = 1 ! use scaling 
    iparm(13)  = 1 ! use maximum weighted matching algorithm 
    iparm(21)  = 1 ! 1x1 + 2x2 pivots

    plist => this%parameter_list 
    if ( this%environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if
    do ilev=1, this%environment%get_num_levels()-1
       ! Set current level Dirichlet solver parameters
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       
       ! Set current level Neumann solver parameters
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
     
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
       plist  => coarse 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
 
    ! Set-up MLBDDC preconditioner
    call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()

    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)

    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, this%mlbddc) 
    					
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%H_current%get_free_dof_values()
				
					call this%setup_nonlinear_solver()
				 temporal: do while ( .not. this%theta_method%finished() ) 
     call this%theta_method%print(6) 

				 ! call this%solve_nonlinear_system()
     call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), dof_values)		
					
				!	if (this%nonlinear_solver%converged() ) then  ! Theta method goes forward 
					call this%theta_method%update_solutions(this%H_current, this%H_previous)
					call this%write_time_step_solution()
					call this%theta_method%move_time_forward()
				!	elseif (.not. this%nonlinear_solver%converged()) then ! Theta method goes backwards and restarts   
    !    call this%theta_method%move_time_backwards(this%H_current, this%H_previous)
    ! end if
					
					if (.not. this%theta_method%finished() ) then 
        call this%fe_space%project_dirichlet_values_curl_conforming(this%H_current,time=this%theta_method%get_current_time(), fields_to_project=(/ 1 /) )
        call this%assemble_system() 
     end if
					
					end do temporal 
   
  end subroutine solve_system
		
		  ! -----------------------------------------------------------------------------------------------
   subroutine solve_nonlinear_system(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    
    call this%nonlinear_solver%initialize() 
								
    nonlinear: do while ( .not. this%nonlinear_solver%finished() )
    ! 0 - Update initial residual 
    call this%nonlinear_solver%start_new_iteration() 
    ! 1 - Integrate Jacobian
    call this%nonlinear_solver%compute_jacobian(this%hts_integration)
    ! 2 - Solve tangent system 
				call solve_tangent_system() 
    ! 3 - Determine step length to update solution 
    call this%nonlinear_solver%line_search%cubic_backtracking( this%nonlinear_solver )
    ! 4 - Update solution 
    call this%nonlinear_solver%update_solution() 
    ! 5 - New picard iterate with updated solution 
    call this%assemble_system()  
    ! 6 - Evaluate new residual 
    call this%nonlinear_solver%compute_residual()
    ! 7 - Print current output 
    call this%nonlinear_solver%print_current_iteration_output() 
   
    end do nonlinear 
   
    call this%nonlinear_solver%print_final_output() 
   contains 
			
			subroutine solve_tangent_system() 
			 
				call this%mlbddc%free()
			 call this%iterative_linear_solver%free()
    call this%setup_solver() 
    call this%iterative_linear_solver%solve( -this%nonlinear_solver%get_residual(),            &  
																																												  this%nonlinear_solver%get_increment_dof_values() )	 

			end subroutine solve_tangent_system
			
  end subroutine solve_nonlinear_system
   
  subroutine check_solution(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    type(error_norms_vector_t) :: error_norm 
				class(vector_function_t), pointer :: H_exact_function
				real(rp) :: l2, linfty, h1_s, hcurl 
				real(rp) :: error_tolerance 
    
				if ( this%test_params%get_is_analytical_solution() ) then 
				
				error_tolerance = 1.0e-3_rp 
    call error_norm%create(this%fe_space,1)   
				H_exact_function => this%maxwell_analytical_functions%get_solution_function()
				
    l2 = error_norm%compute(H_exact_function, this%H_current, l2_norm)     
    linfty = error_norm%compute(H_exact_function, this%H_current, linfty_norm)   
    h1_s = error_norm%compute(H_exact_function, this%H_current, h1_seminorm) 
				hcurl = error_norm%compute(H_exact_function, this%H_current, hcurl_seminorm)
    if ( this%environment%am_i_l1_root() ) then
      write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
      write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
      write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
						write(*,'(a20,e32.25)') 'hcurl_norm:', hcurl; check ( hcurl < error_tolerance )
    end if  
    call error_norm%free()
				
			 end if 
  end subroutine check_solution
  
  subroutine initialize_output(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout)    :: this
	   integer(ip)                                       :: i, istat

		 if ( .not. this%environment%am_i_l1_task() ) then
			return 
			end if 
			
	    if(this%test_params%get_write_solution()) then
				  	call build_set_id_cell_vector()
       call  this%oh%create(VTK) 
       call  this%oh%attach_fe_space(this%fe_space)
       call  this%oh%add_fe_function(this%H_current, 1, 'H')
       call  this%oh%add_fe_function(this%H_current, 1, 'J', curl_diff_operator)
							call  this%oh%add_cell_vector(this%set_id_cell_vector, 'set_id')
       call  this%oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
    endif

  contains 
      subroutine build_set_id_cell_vector()
      call memalloc(this%triangulation%get_num_local_cells(), this%set_id_cell_vector, __FILE__, __LINE__)
      do i=1, this%triangulation%get_num_local_cells()
            this%set_id_cell_vector(i) = this%environment%get_l1_rank() + 1
      enddo
    end subroutine build_set_id_cell_vector
	
	end subroutine initialize_output
	
	  ! -----------------------------------------------------------------------------------------------
  subroutine write_time_step_solution(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout)    :: this
    integer(ip)                                       :: err
    
			if ( .not. this%environment%am_i_l1_task() ) then
			return 
			end if 
			
    if( this%test_params%get_write_solution() .and. this%theta_method%print_this_step() ) then
        call this%oh%append_time_step(this%theta_method%get_current_time())
        call this%oh%write()
        call this%theta_method%update_time_to_be_printed() 
    endif
   
  end subroutine write_time_step_solution
		
	  ! -----------------------------------------------------------------------------------------------
  subroutine finalize_output(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout)    :: this
    integer(ip)                                       :: err
				
			if ( .not. this%environment%am_i_l1_task() ) then
			return 
			end if 
			
    if(this%test_params%get_write_solution()) then
    call this%oh%close()
    call this%oh%free()
				call memfree( this%set_id_cell_vector, __FILE__, __LINE__ )
    endif
  end subroutine finalize_output
		
  ! ***********************************************************************************************
  subroutine run_simulation(this) 
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    call this%setup_triangulation()
    call this%setup_reference_fes()
				call this%setup_theta_method()
	   call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
				call this%setup_fe_coarse_space() 
    call this%setup_solver()
				call this%initialize_output() 
    call this%solve_system()
				call this%finalize_output()
    call this%check_solution()
    
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%H_current%free() 
				call this%H_previous%free() 
    call this%mlbddc%free() 
    call this%iterative_linear_solver%free()
				call this%nonlinear_solver%free()
    call this%fe_affine_operator%free()
	
	if ( this%environment%am_i_l1_task() ) then 
	if (allocated(this%coarse_fe_handlers) ) then 
	select type ( ch => this%fe_space%get_coarse_fe_handler(field_id=1) )
	class is ( Hcurl_l1_coarse_fe_handler_t )
	call this%coarse_fe_handler%free() 
	deallocate( this%coarse_fe_handler ) 
	end select 
	deallocate(this%coarse_fe_handlers, stat=istat); check(istat==0) 
	end if 
	end if 
	
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
      do i=1, size(this%reference_fes)
        call this%reference_fes(i)%free()
      end do
      deallocate(this%reference_fes, stat=istat)
      check(istat==0)
    end if
    call this%triangulation%free()
		
  end subroutine free  

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    call this%environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters
  
end module par_test_hts_driver_names
