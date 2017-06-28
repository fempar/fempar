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
  use fempar_sm_nonlinear_operator_names
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
     type(nonlinear_operator_t)                  :: nonlinear_operator
     type(nonlinear_solver_t)                    :: nonlinear_solver 
     
     ! Problem solution FE function
     type(fe_function_t)                         :: solution
     type(constant_vector_function_t)            :: zero_vector
     type(constant_scalar_function_t)            :: zero_scalar

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
     procedure        , private :: setup_operators
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
    class(vector_t) , pointer :: unknown

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

    ! Construct Linear and nonlinear operators
    call this%timer_solver_setup%start()
    call this%setup_operators()
    call this%timer_solver_setup%stop()
   
    ! Solve the problem
    call this%timer_solver_run%start()
    unknown => this%solution%get_dof_values()
    !write(*,*) unknown%nrm2()
    call this%nonlinear_solver%solve(this%nonlinear_operator, unknown)
    mcheck( this%nonlinear_solver%has_converged(), 'Nonlinear solver has not converged.' )
    call this%timer_solver_run%stop()

    ! Postprocess
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
    integer(ip) :: field_id
    type(vector_field_t) :: zero_vector_field
    class(vector_t), pointer  :: dof_values
    
    ! Solution and initial guess
    call this%solution%create(this%fe_space) 
    call this%zero_scalar%create(0.0_rp)
    call zero_vector_field%init(0.0_rp)
    call this%zero_vector%create(zero_vector_field)
    do field_id = 1, this%fempar_sm_integration%get_number_fields()
       if(this%fempar_sm_integration%get_field_type(field_id) == field_type_vector) then
         call this%solution%interpolate_function(this%fe_space, field_id, this%zero_vector)
       else if( this%fempar_sm_integration%get_field_type(field_id) == field_type_scalar) then
         call this%solution%interpolate_function(this%fe_space, field_id, this%zero_scalar)
       end if
    end do
    call this%solution%update_strong_dirichlet_values(this%fe_space)
    call this%fempar_sm_integration%set_solution(this%solution)
    
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

    ! allocate this%nonlinear_operator to a specific time stepping operator (e.g. BDFk)
    call this%nonlinear_operator%create(this%fe_affine_operator)

    ! Nonlinear solver ! fempar_sm_abs_res_norm_and_rel_inc_norm
    call this%nonlinear_solver%create(convergence_criteria = fempar_sm_abs_res_norm, & 
         &                                         abs_tol = 1.0e-6,  &
         &                                         rel_tol = 1.0e-9, &
         &                                       max_iters = 10   ,  &
         &                                   linear_solver = this%linear_solver, &
         &                                     environment = this%par_environment)    

  end subroutine setup_operators

   
  subroutine check_solution(this)
    implicit none
    class(fempar_sm_fe_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: scalar_error_norm 
    type(error_norms_vector_t) :: vector_error_norm 
    real(rp)    :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    integer(ip) :: field_id

    do field_id = 1, this%fempar_sm_integration%get_number_fields()

       if(this%fempar_sm_integration%get_field_type(field_id) == field_type_vector) then
          call vector_error_norm%create(this%fe_space,field_id)    
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
             write(*,'(a12,a8)')      this%fempar_sm_integration%get_field_name(field_id), ' field: '
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
          end if
          call vector_error_norm%free()
       else if( this%fempar_sm_integration%get_field_type(field_id) == field_type_scalar) then
          call scalar_error_norm%create(this%fe_space,field_id)    
          mean      = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, mean_norm)   
          l1        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, l1_norm)   
          l2        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, l2_norm)   
          lp        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, lp_norm)   
          linfty    = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, linfty_norm)   
          h1_s      = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, h1_seminorm) 
          h1        = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, h1_norm) 
          w1p_s     = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1p_seminorm)   
          w1p       = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1p_norm)   
          w1infty_s = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1infty_seminorm) 
          w1infty   = scalar_error_norm%compute(this%fempar_sm_analytical_functions%get_solution_function_p(), this%solution, w1infty_norm)  
          if ( this%par_environment%am_i_l1_root() ) then
             write(*,'(a12,a8)')      this%fempar_sm_integration%get_field_name(field_id), ' field: '
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
          end if
          call scalar_error_norm%free()

       end if
     end do
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

    call this%nonlinear_operator%free()
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
