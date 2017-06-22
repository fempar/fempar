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
module fempar_sm_driver_names
  use fempar_names
  use fempar_sm_params_names
  use fempar_sm_discrete_integration_names
  use fempar_sm_conditions_names
  use fempar_sm_analytical_functions_names
  use fempar_sm_linear_solver_names
  use fempar_sm_nonlinear_solver_names
# include "debug.i90"

  implicit none
  private

  type fempar_sm_fe_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(fempar_sm_params_t)       :: test_params
     type(ParameterList_t), pointer :: parameter_list
     
     ! Cells and lower dimension objects container
     type(par_triangulation_t)                   :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                        :: fe_space 
     type(p_reference_fe_t), allocatable         :: reference_fes(:) 
     type(standard_l1_coarse_fe_handler_t)       :: coarse_fe_handler
     type(p_l1_coarse_fe_handler_t), allocatable :: coarse_fe_handlers(:)
     

     class(fempar_sm_discrete_integration_t), allocatable :: fempar_sm_integration
     type(fempar_sm_conditions_t)                :: fempar_sm_conditions
     type(fempar_sm_analytical_functions_t)      :: fempar_sm_analytical_functions
     
     ! Operators to define and solve the problem
     type(fe_affine_operator_t)                  :: fe_affine_operator
     type(mlbddc_t)                              :: mlbddc
     type(linear_solver_t)                       :: linear_solver
     type(nonlinear_solver_t)                    :: nonlinear_solver 
     
     ! Problem solution FE function
     type(fe_function_t)                         :: solution
     
     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                         :: par_environment

     ! Timers
     type(timer_t) :: timer_problem
     type(timer_t) :: timer_triangulation
     type(timer_t) :: timer_fe_space
     type(timer_t) :: timer_assemply
     type(timer_t) :: timer_solver_setup
     type(timer_t) :: timer_solver_run

   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_timers
     procedure                  :: report_timers
     procedure                  :: free_timers
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_discrete_integration     
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     !procedure        , private :: setup_system
     procedure        , private :: setup_operators
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure                  :: run_simulation
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
  end type fempar_sm_fe_driver_t

  ! Types
  public :: fempar_sm_fe_driver_t

contains

  subroutine run_simulation(this) 
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    
    ! Geometry
    call this%timer_triangulation%start()
    call this%setup_triangulation()
    call this%timer_triangulation%stop()

    ! Problem and its FE approximation
    call this%timer_fe_space%start()
    call this%setup_discrete_integration()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%timer_fe_space%stop()

    ! Algebraic operators
    !call this%timer_assemply%start()
    !call this%setup_system()
    !call this%assemble_system()
    !call this%timer_assemply%stop()

    call this%timer_solver_setup%start()
    call this%setup_operators()
    !call this%setup_solver()
    call this%timer_solver_setup%stop()
   
    call this%timer_solver_run%start()
    call this%nonlinear_solver%solve()
    !mcheck( this%nonlinear_solver%has_converged(), 'Nonlinear solver has not converged.' )
    !call this%solve_system()
    call this%timer_solver_run%stop()

    call this%write_solution()
    call this%check_solution()
    call this%free()
  end subroutine run_simulation
 
!========================================================================================

  subroutine parse_command_line_parameters(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

!========================================================================================
subroutine setup_timers(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    class(execution_context_t), pointer :: w_context
    w_context => this%par_environment%get_w_context()
    call this%timer_triangulation%create(w_context,"SETUP TRIANGULATION")
    call this%timer_fe_space%create(     w_context,"SETUP FE SPACE")
    call this%timer_assemply%create(     w_context,"FE INTEGRATION AND ASSEMBLY")
    call this%timer_solver_setup%create( w_context,"SETUP SOLVER AND PRECONDITIONER")
    call this%timer_solver_run%create(   w_context,"SOLVER RUN")
end subroutine setup_timers

!========================================================================================
subroutine report_timers(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%report(.true.)
    call this%timer_fe_space%report(.false.)
    call this%timer_assemply%report(.false.)
    call this%timer_solver_setup%report(.false.)
    call this%timer_solver_run%report(.false.)
    if (this%par_environment%get_l1_rank() == 0) then
      write(*,*)
    end if
end subroutine report_timers

!========================================================================================
subroutine free_timers(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_solver_setup%free()
    call this%timer_solver_run%free()
end subroutine free_timers

!========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%par_environment%create (this%parameter_list)
  end subroutine setup_environment

!========================================================================================
  subroutine setup_triangulation(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    type(vef_iterator_t)  :: vef
    logical :: fixed_pressure
    
    call this%triangulation%create(this%parameter_list, this%par_environment)
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       fixed_pressure = .false.
       call this%triangulation%create_vef_iterator(vef)
       do while ( .not. vef%has_finished() )
          if(vef%is_at_boundary()) then
             call vef%set_set_id(1)
             if(.not.fixed_pressure) then
                fixed_pressure = .true.
                call vef%set_set_id(2)
             end if
          else
             call vef%set_set_id(0)
          end if
          call vef%next()
       end do
       call this%triangulation%free_vef_iterator(vef)
    end if  
    call this%triangulation%setup_coarse_triangulation()
  end subroutine setup_triangulation

!========================================================================================

  subroutine setup_discrete_integration(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_discrete_integration_type() == discrete_integration_type_irreducible ) then
       allocate(irreducible_discrete_integration_t :: this%fempar_sm_integration, stat=istat); check(istat==0)
    else if(this%test_params%get_discrete_integration_type() == discrete_integration_type_mixed_u_p ) then
       allocate(mixed_u_p_discrete_integration_t :: this%fempar_sm_integration, stat=istat); check(istat==0)
   else
       write(*,*) this%test_params%get_discrete_integration_type()
       write(*,*) discrete_integration_type_irreducible
       mcheck(.false.,'Discrete integration not available')
   end if
   call this%fempar_sm_integration%create(this%triangulation%get_num_dimensions(),this%fempar_sm_analytical_functions)
   !call this%fempar_sm_integration%set_solution(this%solution)
   
  end subroutine setup_discrete_integration

!========================================================================================
  subroutine setup_reference_fes(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat, field_id
    class(cell_iterator_t), allocatable       :: cell
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    
    allocate(this%reference_fes(this%fempar_sm_integration%get_number_fields()), stat=istat)
    check(istat==0)
    
    if ( this%par_environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      reference_fe_geo => cell%get_reference_fe_geo()
      do field_id = 1, this%fempar_sm_integration%get_number_fields()
         this%reference_fes(field_id) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                             fe_type = this%fempar_sm_integration%get_fe_type(field_id), &
                                                             number_dimensions = this%triangulation%get_num_dimensions(), &
                                                             order = this%test_params%get_reference_fe_order(), &
                                                             field_type = this%fempar_sm_integration%get_field_type(field_id), &
                                                             conformity = .true. )
      end do      
      call this%triangulation%free_cell_iterator(cell)
    end if  
  end subroutine setup_reference_fes

!========================================================================================

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(fempar_sm_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat, field_id
    allocate(this%coarse_fe_handlers(this%fempar_sm_integration%get_number_fields()), stat=istat); check(istat==0)
    do field_id = 1, this%fempar_sm_integration%get_number_fields()
       this%coarse_fe_handlers(field_id)%p => this%coarse_fe_handler
    end do
  end subroutine setup_coarse_fe_handlers

!========================================================================================

  subroutine setup_fe_space(this)
    implicit none
    class(fempar_sm_fe_driver_t), target, intent(inout) :: this
    class(reference_fe_t)       , pointer       :: reference_fe
    
    call this%fempar_sm_analytical_functions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%fempar_sm_conditions%set_number_components(this%fempar_sm_integration%get_number_components())
    call this%fempar_sm_conditions%set_number_dimensions(this%triangulation%get_num_dimensions())    
    ! B.C. taken from the exact solution
    call this%fempar_sm_conditions%set_boundary_function(this%fempar_sm_analytical_functions%get_solution_function_u())

    call this%fe_space%create( triangulation       = this%triangulation, &
                               conditions          = this%fempar_sm_conditions, &
                               reference_fes       = this%reference_fes, &
                               coarse_fe_handlers  = this%coarse_fe_handlers, &
                               field_blocks        = this%fempar_sm_integration%get_field_blocks(),            &
                               field_coupling      = this%fempar_sm_integration%get_field_coupling())
    
    call this%fe_space%fill_dof_info() 
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    call this%fe_space%initialize_fe_integration()
    call this%fe_space%interpolate_dirichlet_values(this%fempar_sm_conditions)    

  end subroutine setup_fe_space

  !========================================================================================

  subroutine setup_operators (this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    type(parameterlist_t) :: linear_pl
    integer(ip) :: ilev
    integer(ip) :: FPLError
    real(rp)    :: res_norm
    class(vector_t), pointer  :: dof_values

    ! Solution and initial guess
    call this%solution%create(this%fe_space) 
    call this%fempar_sm_integration%init_solution(this%fe_space, this%solution)
    call this%solution%update_strong_dirichlet_values(this%fe_space)
    
    ! FE operator
    if(this%fempar_sm_integration%is_symmetric().and.this%fempar_sm_integration%is_coercive()) then
       call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
            &                                diagonal_blocks_symmetric_storage = [ .true. ], &
            &                                diagonal_blocks_symmetric         = [ .true. ], &
            &                                diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
            &                                fe_space                          = this%fe_space, &
            &                                discrete_integration              = this%fempar_sm_integration )
    else if(this%fempar_sm_integration%is_symmetric()) then
       call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
            &                                diagonal_blocks_symmetric_storage = [ .true. ], &
            &                                diagonal_blocks_symmetric         = [ .true. ], &
            &                                diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_INDEFINITE ], &
            &                                fe_space                          = this%fe_space, &
            &                                discrete_integration              = this%fempar_sm_integration )
    end if

    ! BDDC preconditioner
    plist => this%parameter_list 
    if ( this%par_environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    end if
    do ilev=1, this%par_environment%get_num_levels()-1
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
       plist  => coarse 
    end do
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)

    ! Set coarsest-grid solver type (currently NOT inherited from fine level matrices types)
    if(this%fempar_sm_integration%is_coercive()) then
       FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
    else
       FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
    end if

    call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    call this%mlbddc%symbolic_setup()
    !call this%mlbddc%numerical_setup()

    ! To debug
    !dof_values => this%fe_affine_operator%get_translation()
    !res_norm = dof_values%nrm2()
    !write(*,*) res_norm
    
    ! Linear solver
    call this%linear_solver%create(this%fe_space%get_environment())
    if(this%fempar_sm_integration%is_coercive()) then
       call this%linear_solver%set_type_from_string(cg_name)
    else
       call this%linear_solver%set_type_from_string(lgmres_name)
    end if
    call this%linear_solver%setup_operators(this%fe_affine_operator, this%mlbddc) 
    call linear_pl%init()
    FPLError = linear_pl%set(key = ils_rtol, value = 1.0e-12_rp); assert(FPLError == 0)
    FPLError = linear_pl%set(key = ils_max_num_iterations, value = 5000); assert(FPLError == 0)
    FPLError = linear_pl%set(key = ils_atol, value = 1.0e-9); assert(FPLError == 0)
    call this%linear_solver%set_parameters_from_pl(linear_pl) 

    ! Nonlinear solver ! fempar_sm_abs_res_norm_and_rel_inc_norm
    call this%nonlinear_solver%create(convergence_criteria = fempar_sm_abs_res_norm, & 
         &                                         abs_tol = 1.0e-6,  &
         &                                         rel_tol = 1.0e-9, &
         &                                       max_iters = 10   ,  &
         &                                   linear_solver = this%linear_solver, &
         &                                         unknown = this%solution%get_dof_values(), &
         &                                     environment = this%par_environment)    

  end subroutine setup_operators

  !========================================================================================

  subroutine setup_solver (this)
    implicit none
    class(fempar_sm_fe_driver_t), target, intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse

    integer(ip) :: ilev
    integer(ip) :: FPLError
    integer(ip) :: iparm(64)
    logical, parameter :: enable_mkl=.true. ! While the compilations system does not manage macros for drivers

    !#ifdef ENABLE_MKL  
    if(enable_mkl) then
       ! See https://software.intel.com/en-us/node/470298 for details
       iparm      = 0 ! Init all entries to zero
       !iparm(1)   = 1 ! no solver default
       !iparm(2)   = 2 ! fill-in reordering from METIS
       !iparm(8)   = 2 ! numbers of iterative refinement steps
       !iparm(10)  = 8 ! perturb the pivot elements with 1E-8
       !iparm(11)  = 1 ! use scaling 
       !iparm(13)  = 1 ! use maximum weighted matching algorithm 
       !iparm(21)  = 1 ! 1x1 + 2x2 pivots

       plist => this%parameter_list 
       if ( this%par_environment%get_l1_size() == 1 ) then
          FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
          !FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
          !FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          !FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       end if
       do ilev=1, this%par_environment%get_num_levels()-1
          ! Set current level Dirichlet solver parameters
          dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
          FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
          !FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_uns); assert(FPLError == 0)
          !FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          !FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

          ! Set current level Neumann solver parameters
          neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
          FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
          !FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_uns); assert(FPLError == 0)
          !FPLError = neumann%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
          !FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

          coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
          plist  => coarse 
       end do
       ! Set coarsest-grid solver parameters
       FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       !FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       !FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

       ! Set-up MLBDDC preconditioner
       call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
       call this%mlbddc%symbolic_setup()
       call this%mlbddc%numerical_setup()
    end if
    !#endif    

    call this%linear_solver%create(this%fe_space%get_environment())
    call this%linear_solver%set_type_from_string(lgmres_name)

    !#ifdef ENABLE_MKL
    if(enable_mkl) then
       call this%linear_solver%set_operators(this%fe_affine_operator, this%mlbddc) 
       !#else
    else
       call parameter_list%init()
       FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-12_rp)
       assert(FPLError == 0)
       FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
       assert(FPLError == 0)
       call this%linear_solver%set_parameters_from_pl(parameter_list)
       call this%linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
       call parameter_list%free()
    end if
    !#endif   

  end subroutine setup_solver
  
!========================================================================================

  subroutine assemble_system (this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
!    class(matrix_t)                  , pointer       :: matrix
!    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    !rhs                => this%fe_affine_operator%get_translation()
    !matrix             => this%fe_affine_operator%get_matrix()
    !select type(matrix)
    !class is (sparse_matrix_t)  
    !   call matrix%print_matrix_market(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
    
    !select type(rhs)
    !class is (serial_scalar_array_t)  
    !   call rhs%print(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine assemble_system
  
!========================================================================================

  subroutine solve_system(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    !class(matrix_t)                         , pointer       :: matrix
    !class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    !matrix     => this%fe_affine_operator%get_matrix()
    !rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    call this%linear_solver%solve(this%fe_affine_operator%get_translation(), &
                                            dof_values)
    
    !select type (dof_values)
    !class is (serial_scalar_array_t)  
    !   call dof_values%print(6)
    !class DEFAULT
    !   assert(.false.) 
    !end select
    
    !select type (matrix)
    !class is (sparse_matrix_t)  
    !   call this%direct_solver%update_matrix(matrix, same_nonzero_pattern=.true.)
    !   call this%direct_solver%solve(rhs , dof_values )
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine solve_system
   
  subroutine check_solution(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: scalar_error_norm 
    type(error_norms_vector_t) :: vector_error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    
    call vector_error_norm%create(this%fe_space,1)    
    mean      = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, mean_norm)   
    l1        = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, l1_norm)   
    l2        = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, l2_norm)   
    lp        = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, lp_norm)   
    linfty    = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, linfty_norm)   
    h1_s      = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, h1_seminorm) 
    h1        = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, h1_norm) 
    w1p_s     = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, w1p_seminorm)   
    w1p       = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, w1p_norm)   
    w1infty_s = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, w1infty_seminorm) 
    w1infty   = vector_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_u(), this%solution, w1infty_norm)  
    if ( this%par_environment%am_i_l1_root() ) then
      write(*,'(a20,e32.25)') 'First field (vector):'
      write(*,'(a20,e32.25)') 'mean_norm:', mean
      write(*,'(a20,e32.25)') 'l1_norm:', l1
      write(*,'(a20,e32.25)') 'l2_norm:', l2
      write(*,'(a20,e32.25)') 'lp_norm:', lp
      write(*,'(a20,e32.25)') 'linfnty_norm:', linfty
      write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s
      write(*,'(a20,e32.25)') 'h1_norm:', h1
      write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s
      write(*,'(a20,e32.25)') 'w1p_norm:', w1p
      write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s
      write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty
      !write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < 1.0e-04 )
      !write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < 1.0e-04 )
    end if  
    call vector_error_norm%free()

    !call scalar_error_norm%create(this%fe_space,2)    
    !mean      = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, mean_norm)   
    !l1        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, l1_norm)   
    !l2        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, l2_norm)   
    !lp        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, lp_norm)   
    !linfty    = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, linfty_norm)   
    !h1_s      = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, h1_seminorm) 
    !h1        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, h1_norm) 
    !w1p_s     = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1p_seminorm)   
    !w1p       = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1p_norm)   
    !w1infty_s = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1infty_seminorm) 
    !w1infty   = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1infty_norm)  
    !if ( this%par_environment%am_i_l1_root() ) then
    !  write(*,'(a20,e32.25)') 'Second field (scalar):'
    !  write(*,'(a20,e32.25)') 'mean_norm:', mean
    !  write(*,'(a20,e32.25)') 'l1_norm:', l1
    !  write(*,'(a20,e32.25)') 'l2_norm:', l2
    !  write(*,'(a20,e32.25)') 'lp_norm:', lp
    !  write(*,'(a20,e32.25)') 'linfnty_norm:', linfty
    !  write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s
    !  write(*,'(a20,e32.25)') 'h1_norm:', h1
    !  write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s
    !  write(*,'(a20,e32.25)') 'w1p_norm:', w1p
    !  write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s
    !  write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty
    !  !write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < 1.0e-04 )
    !  !write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < 1.0e-04 )
    !end if  
    !call scalar_error_norm%free()

  end subroutine check_solution
  
  !========================================================================================
   
  subroutine write_solution(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    type(output_handler_t)                   :: oh
    integer(ip) :: field_id
    character(len=:), allocatable :: name
    
    if(this%test_params%get_write_solution()) then
       if ( this%par_environment%am_i_l1_task() ) then
        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        do field_id = 1, this%fempar_sm_integration%get_number_fields()
           name =  this%fempar_sm_integration%get_field_name(field_id)
           call oh%add_fe_function(this%solution, field_id, name)
        end do
        call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
        call oh%write()
        call oh%close()
        call oh%free()
        end if
    endif
  end subroutine write_solution
  
 !========================================================================================

  subroutine free(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat

    call this%nonlinear_solver%free()
    call this%linear_solver%free()
    call this%solution%free()
    call this%mlbddc%free()
    call this%fe_affine_operator%free()
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
      do i=1, size(this%reference_fes)
        call this%reference_fes(i)%free()
      end do
      deallocate(this%reference_fes, stat=istat)
      mcheck(istat==0,'Error deallocating')
    end if
    if (allocated(this%coarse_fe_handlers)) then
      deallocate(this%coarse_fe_handlers, stat=istat)
      mcheck(istat==0,'Error deallocating')
    end if
    call this%triangulation%free()
    call this%fempar_sm_integration%free()
  end subroutine free  

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters
  
end module fempar_sm_driver_names
