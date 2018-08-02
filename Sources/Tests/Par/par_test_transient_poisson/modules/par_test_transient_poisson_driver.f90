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
module par_test_transient_poisson_driver_names
  use fempar_names
  use par_test_transient_poisson_params_names
  use transient_poisson_discrete_integration_names
  use transient_poisson_conditions_names
  use transient_poisson_analytical_functions_names
  use time_stepping_names

# include "debug.i90"

  implicit none
  private

  integer(ip), parameter :: PAR_TEST_POISSON_FULL = 1 ! Has to be == 1
  integer(ip), parameter :: PAR_TEST_POISSON_VOID = 2

  type par_test_transient_poisson_fe_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_transient_poisson_params_t)      :: test_params
     type(ParameterList_t), pointer                 :: parameter_list
     
     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation
     integer(ip), allocatable              :: cell_set_ids(:)
     
     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                                :: fe_space 
     type(p_reference_fe_t), allocatable                 :: reference_fes(:) 
     type(standard_l1_coarse_fe_handler_t)               :: coarse_fe_handler
     type(p_l1_coarse_fe_handler_t), allocatable         :: coarse_fe_handlers(:)
     type(transient_poisson_CG_discrete_integration_t)   :: poisson_integration
     type(transient_poisson_conditions_t)                :: poisson_conditions
     type(transient_poisson_analytical_functions_t)      :: poisson_analytical_functions

     
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                   :: fe_op
     type(nonlinear_solver_t)                     :: nl_solver
     type(time_stepping_operator_t)               :: time_operator
     type(dirk_solver_t)                          :: time_solver
     
#ifdef ENABLE_MKL     
     ! MLBDDC preconditioner
     type(mlbddc_t)                            :: mlbddc
#endif  
    
     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)           :: iterative_linear_solver
 
     ! Poisson problem solution FE function
     type(fe_function_t)                   :: solution
     
     ! Environment required for fe_op + vtk_handler
     type(environment_t)                    :: par_environment

     ! Timers
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
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure                  :: run_simulation
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
     procedure, nopass, private :: popcorn_fun => par_test_poisson_driver_popcorn_fun
  end type par_test_transient_poisson_fe_driver_t

  ! Types
  public :: par_test_transient_poisson_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

!========================================================================================
subroutine setup_timers(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
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
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
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
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_solver_setup%free()
    call this%timer_solver_run%free()
end subroutine free_timers

!========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    !istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%par_environment%create (this%parameter_list)
  end subroutine setup_environment
   
  subroutine setup_triangulation(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this

    class(cell_iterator_t), allocatable :: cell
    type(point_t), allocatable :: cell_coords(:)
    integer(ip) :: istat
    integer(ip) :: set_id
    real(rp) :: x, y
    integer(ip) :: num_void_neigs

    integer(ip)           :: ivef
    class(vef_iterator_t), allocatable  :: vef, vef_of_vef
    type(list_t), pointer :: vefs_of_vef
    type(list_t), pointer :: vertices_of_line
    type(list_iterator_t) :: vefs_of_vef_iterator
    type(list_iterator_t) :: vertices_of_line_iterator
    class(reference_fe_t), pointer :: reference_fe_geo
    integer(ip) :: ivef_pos_in_cell, vef_of_vef_pos_in_cell
    integer(ip) :: vertex_pos_in_cell, icell_arround
    integer(ip) :: inode, num


    call this%triangulation%create(this%parameter_list, this%par_environment)

    ! Set the cell ids to use void fes
    if (this%test_params%get_use_void_fes()) then
      if ( this%par_environment%am_i_l1_task() ) then
        call memalloc(this%triangulation%get_num_local_cells(),this%cell_set_ids)
        call this%triangulation%create_cell_iterator(cell)
        allocate(cell_coords(1:cell%get_num_nodes()),stat=istat); check(istat == 0)
        do while( .not. cell%has_finished() )
          if (cell%is_local()) then
            set_id = PAR_TEST_POISSON_VOID
            call cell%get_nodes_coordinates(cell_coords)
            select case (trim(this%test_params%get_use_void_fes_case()))
            case ('half')
              y = cell_coords(1)%get(2)
              if (y>=0.5) set_id = PAR_TEST_POISSON_FULL
            case ('quarter')
              x = cell_coords(1)%get(1)
              y = cell_coords(1)%get(2)
              if (x>=0.5 .and. y>=0.5) set_id = PAR_TEST_POISSON_FULL
            case ('popcorn')
              do inode = 1,cell%get_num_nodes()
                if ( this%popcorn_fun(cell_coords(inode),&
                  this%triangulation%get_num_dims()) < 0.0 ) then
                  set_id = PAR_TEST_POISSON_FULL
                  exit
                end if
              end do
            case default
              check(.false.)
            end select
            this%cell_set_ids(cell%get_gid()) = set_id
          end if
          call cell%next()
        end do
        deallocate(cell_coords, stat = istat); check(istat == 0)
        call this%triangulation%fill_cells_set(this%cell_set_ids)
      end if
        call this%triangulation%free_cell_iterator(cell)
    end if

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

    ! Set all the vefs on the interface between full/void if there are void fes
    if (this%test_params%get_use_void_fes()) then
      call this%triangulation%create_vef_iterator(vef)
      call this%triangulation%create_vef_iterator(vef_of_vef)
      call this%triangulation%create_cell_iterator(cell)
      do while ( .not. vef%has_finished() )

         ! If it is an INTERIOR face
         if( vef%get_dim() == this%triangulation%get_num_dims()-1 .and. vef%get_num_cells_around()==2 ) then

           ! Compute number of void neighbors
           num_void_neigs = 0
           do icell_arround = 1,vef%get_num_cells_around()
             call vef%get_cell_around(icell_arround,cell)
             if (cell%get_set_id() == PAR_TEST_POISSON_VOID) num_void_neigs = num_void_neigs + 1
           end do

           if(num_void_neigs==1) then ! If vef (face) is between a full and a void cell

               ! Set this face as Dirichlet boundary
               call vef%set_set_id(1)

               ! Do a loop on all edges in 3D (vertex in 2D) of the face
               ivef = vef%get_gid()
               call vef%get_cell_around(1,cell) ! There is always one cell around
               reference_fe_geo => cell%get_reference_fe()
               ivef_pos_in_cell = cell%get_vef_lid_from_gid(ivef)
               vefs_of_vef => reference_fe_geo%get_facets_n_face()
               vefs_of_vef_iterator = vefs_of_vef%create_iterator(ivef_pos_in_cell)
               do while( .not. vefs_of_vef_iterator%is_upper_bound() )

                  ! Set edge (resp. vertex) as Dirichlet
                  vef_of_vef_pos_in_cell = vefs_of_vef_iterator%get_current()
                  call cell%get_vef(vef_of_vef_pos_in_cell, vef_of_vef)
                  call vef_of_vef%set_set_id(1)

                  ! If 3D, traverse vertices of current line
                  if ( this%triangulation%get_num_dims() == 3 ) then
                    vertices_of_line          => reference_fe_geo%get_vertices_n_face()
                    vertices_of_line_iterator = vertices_of_line%create_iterator(vef_of_vef_pos_in_cell)
                    do while( .not. vertices_of_line_iterator%is_upper_bound() )

                      ! Set vertex as Dirichlet
                      vertex_pos_in_cell = vertices_of_line_iterator%get_current()
                      call cell%get_vef(vertex_pos_in_cell, vef_of_vef)
                      call vef_of_vef%set_set_id(1)

                      call vertices_of_line_iterator%next()
                    end do ! Loop in vertices in 3D only
                  end if

                  call vefs_of_vef_iterator%next()
               end do ! Loop in edges (resp. vertices)

           end if ! If face on void/full boundary
         end if ! If vef is an interior face

         call vef%next()
      end do ! Loop in vefs
      call this%triangulation%free_cell_iterator(cell)
      call this%triangulation%free_vef_iterator(vef)
      call this%triangulation%free_vef_iterator(vef_of_vef)
    end if

    call this%triangulation%setup_coarse_triangulation()

  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable       :: cell
    class(reference_fe_t), pointer :: reference_fe_geo
    
    if (this%test_params%get_use_void_fes()) then
      allocate(this%reference_fes(2), stat=istat)
    else
      allocate(this%reference_fes(1), stat=istat)
    end if
    check(istat==0)
    
    if ( this%par_environment%am_i_l1_task() ) then
      call this%triangulation%create_cell_iterator(cell)
      reference_fe_geo => cell%get_reference_fe()
      this%reference_fes(PAR_TEST_POISSON_FULL) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_lagrangian, &
                                                   num_dims = this%triangulation%get_num_dims(), &
                                                   order = this%test_params%get_reference_fe_order(), &
                                                   field_type = field_type_scalar, &
                                                   conformity = .true. )
      if (this%test_params%get_use_void_fes()) then
        this%reference_fes(PAR_TEST_POISSON_VOID) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                   fe_type = fe_type_void, &
                                                   num_dims = this%triangulation%get_num_dims(), &
                                                   order = -1, &
                                                   field_type = field_type_scalar, &
                                                   conformity = .true. )
      end if
      call this%triangulation%free_cell_iterator(cell)
    end if
  end subroutine setup_reference_fes
  
  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat
    allocate(this%coarse_fe_handlers(1), stat=istat)
    check(istat==0)
    this%coarse_fe_handlers(1)%p => this%coarse_fe_handler
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this

    integer(ip) :: set_ids_to_reference_fes(1,2)

    call this%poisson_analytical_functions%set_num_dims(this%triangulation%get_num_dims())
    call this%poisson_conditions%set_boundary_function(this%poisson_analytical_functions%get_boundary_function())

    if (this%test_params%get_use_void_fes()) then
      set_ids_to_reference_fes(1,PAR_TEST_POISSON_FULL) = PAR_TEST_POISSON_FULL
      set_ids_to_reference_fes(1,PAR_TEST_POISSON_VOID) = PAR_TEST_POISSON_VOID
      call this%fe_space%create( triangulation            = this%triangulation,       &
                                 reference_fes            = this%reference_fes,       &
                                 set_ids_to_reference_fes = set_ids_to_reference_fes, &
                                 coarse_fe_handlers       = this%coarse_fe_handlers,  &
                                 conditions               = this%poisson_conditions )
    else
      call this%fe_space%create( triangulation       = this%triangulation,      &
                                 reference_fes       = this%reference_fes,      &
                                 coarse_fe_handlers  = this%coarse_fe_handlers, &
                                 conditions          = this%poisson_conditions )
    end if
    
    call this%fe_space%set_up_cell_integration()
    !call this%fe_space%print()
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    
    call this%poisson_integration%set_analytical_functions(this%poisson_analytical_functions)
    
    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_op%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .true. ], &
                                          diagonal_blocks_symmetric         = [ .true. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                          fe_space                          = this%fe_space, &
                                          discrete_integration              = this%poisson_integration )
    
    call this%time_operator%create( fe_op                   = this%fe_op, &
                                    initial_time            = this%test_params%get_initial_time() , &
                                    final_time              = this%test_params%get_final_time() , &
                                    final_time_rel_dt_tol   = 1.0e-6_rp, &
                                    time_step               = this%test_params%get_time_step() , &
                                    time_integration_scheme = this%test_params%get_time_integration_scheme() )  
  
    call this%solution%create(this%fe_space) 
    call this%poisson_integration%set_up( fe_space = this%fe_space, fe_function = this%solution ) 
    call this%poisson_integration%set_analytical_functions(this%poisson_analytical_functions) 
    
    call this%fe_space%interpolate(field_id=1, &
                                   function = this%poisson_analytical_functions%get_solution_function(), &
                                   fe_function=this%solution, &
                                   time=this%time_operator%get_current_time())
  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), target, intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse, coarse_matrix_params

    integer(ip) :: ilev
    integer(ip) :: FPLError
    integer(ip) :: iparm(64)

#ifdef ENABLE_MKL  
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
    if ( this%par_environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if
    do ilev=1, this%par_environment%get_num_levels()-1
       ! Set current level Dirichlet solver parameters
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       if ( ilev > 1 ) then
         FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_uns); assert(FPLError == 0)
       else
         FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       end if 
       FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       
       ! Set current level Neumann solver parameters
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       if ( ilev > 1 ) then
         FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_uns); assert(FPLError == 0)
       else
         FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       end if
       FPLError = neumann%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       
       ! Set current level Coarse-grid matrix parameters
       coarse_matrix_params => plist%NewSubList(key=mlbddc_coarse_matrix_params)
       FPLError = coarse_matrix_params%set(key=mlbddc_coarse_matrix_symmetric_storage, value=.false.); assert(FPLError == 0)
       FPLError = coarse_matrix_params%set(key=mlbddc_coarse_matrix_is_symmetric     , value=.false.); assert(FPLError == 0)
       FPLError = coarse_matrix_params%set(key=mlbddc_coarse_matrix_sign             , value=SPARSE_MATRIX_SIGN_UNKNOWN); assert(FPLError == 0)
      
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
       plist  => coarse 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_uns); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

    ! Set-up MLBDDC preconditioner
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    call this%mlbddc%create(this%time_operator%get_fe_operator(), this%parameter_list)
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()
#endif    
   
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    
    call this%nl_solver%create( convergence_criteria = abs_res_norm, &
                               abs_tol = 1.0e-6_rp, &
                               rel_tol = 1.0e-6_rp, &
                               max_iters = 10_ip, &
                               linear_solver = this%iterative_linear_solver, &
                               environment = this%fe_space%get_environment(),&
                               fe_operator = this%time_operator%get_fe_operator())
   
   call this%time_solver%create( ts_op = this%time_operator, &
                                 nl_solver = this%nl_solver )

#ifdef ENABLE_MKL
    call this%iterative_linear_solver%set_operators(this%time_operator%get_matrix(), this%mlbddc) 
#else
    call parameter_list%init()
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-12_rp)
    assert(FPLError == 0)
    FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
    assert(FPLError == 0)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
    call this%iterative_linear_solver%set_operators(this%time_operator%get_matrix(), .identity. this%time_operator%get_fe_operator()) 
    call parameter_list%free()
#endif   
    
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    call this%time_operator%initial_assembly( this%solution%get_free_dof_values() )
  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    do while ( .not. this%time_operator%has_finished() )
       call this%time_operator%print(6)
       call this%time_solver%advance_fe_function(this%solution)
    end do
  end subroutine solve_system
   
  subroutine check_solution(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    real(rp)                                      :: current_time 
    type(error_norms_scalar_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    current_time = this%time_operator%get_current_time()
    call error_norm%create(this%fe_space,1)    
    mean = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, mean_norm, time=current_time)   
    l1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l1_norm, time=current_time)   
    l2 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l2_norm, time=current_time)   
    lp = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, lp_norm, time=current_time)   
    linfty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, linfty_norm, time=current_time)   
    h1_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_seminorm, time=current_time) 
    h1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_norm, time=current_time) 
    w1p_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_seminorm, time=current_time)   
    w1p = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_norm, time=current_time)   
    w1infty_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm, time=current_time) 
    w1infty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_norm, time=current_time)  
    if ( this%par_environment%am_i_l1_root() ) then
      write(*,'(a20,e32.25)') 'mean_norm:', mean; !check ( abs(mean) < 1.0e-04 )
      write(*,'(a20,e32.25)') 'l1_norm:', l1; !check ( l1 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'l2_norm:', l2; !check ( l2 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'lp_norm:', lp; !check ( lp < 1.0e-04 )
      write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; !check ( linfty < 1.0e-04 )
      write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; !check ( h1_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'h1_norm:', h1; !check ( h1 < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; !check ( w1p_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1p_norm:', w1p; !check ( w1p < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; !check ( w1infty_s < 1.0e-04 )
      write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; !check ( w1infty < 1.0e-04 )
    end if  
    call error_norm%free()
  end subroutine check_solution
  
  subroutine write_solution(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(in) :: this
    type(output_handler_t)                          :: oh
    real(rp),allocatable :: cell_vector(:)
    real(rp),allocatable :: mypart_vector(:)

    if(this%test_params%get_write_solution()) then
      if (this%par_environment%am_i_l1_task()) then

        if (this%test_params%get_use_void_fes()) then
          call memalloc(this%triangulation%get_num_local_cells(),cell_vector,__FILE__,__LINE__)
          cell_vector(:) = this%cell_set_ids(:)
        end if

        call memalloc(this%triangulation%get_num_local_cells(),mypart_vector,__FILE__,__LINE__)
        mypart_vector(:) = this%par_environment%get_l1_rank()

        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, 1, 'solution')
        call oh%add_fe_function(this%solution, 1, 'grad_solution', grad_diff_operator)
        if (this%test_params%get_use_void_fes()) then
          call oh%add_cell_vector(cell_vector,'cell_set_ids')
        end if
        call oh%add_cell_vector(mypart_vector,'l1_rank')
        call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
        call oh%write()
        call oh%close()
        call oh%free()

        if (allocated(cell_vector)) call memfree(cell_vector,__FILE__,__LINE__)
        call memfree(mypart_vector,__FILE__,__LINE__)

      end if
    endif
  end subroutine write_solution
  
  subroutine run_simulation(this) 
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this

    call this%timer_triangulation%start()
    call this%setup_triangulation()
    call this%timer_triangulation%stop()

    call this%timer_fe_space%start()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%timer_fe_space%stop()

    call this%timer_assemply%start()
    call this%setup_system()
    call this%assemble_system()
    call this%timer_assemply%stop()

    call this%timer_solver_setup%start()
    call this%setup_solver()
    call this%timer_solver_setup%stop()

    call this%timer_solver_run%start()
    call this%solve_system()
    call this%timer_solver_run%stop()

    call this%check_solution()
    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free() 
    call this%iterative_linear_solver%free()
    call this%time_operator%free()
    call this%fe_op%free()
    call this%time_solver%free()
    call this%nl_solver%free()
    call this%poisson_integration%free()
#ifdef ENABLE_MKL    
    call this%mlbddc%free()
#endif       
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
      do i=1, size(this%reference_fes)
        call this%reference_fes(i)%free()
      end do
      deallocate(this%reference_fes, stat=istat)
      check(istat==0)
    end if
    call this%triangulation%free()
    if (allocated(this%cell_set_ids)) call memfree(this%cell_set_ids,__FILE__,__LINE__)
  end subroutine free  

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_test_transient_poisson_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters

  function par_test_poisson_driver_popcorn_fun(point,num_dim) result (val)
    implicit none
    type(point_t), intent(in) :: point
    integer(ip),   intent(in) :: num_dim
    real(rp) :: val
    type(point_t) :: p
    real(rp) :: x, y, z
    real(rp) :: xk, yk, zk
    real(rp) :: r0, sg, A
    integer(ip) :: k
    p = point
    if (num_dim < 3) call p%set(3,0.62)
    x = ( 2.0*p%get(1) - 1.0 )
    y = ( 2.0*p%get(2) - 1.0 )
    z = ( 2.0*p%get(3) - 1.0 )
    r0 = 0.6
    sg = 0.2
    A  = 2.0
    val = sqrt(x**2 + y**2 + z**2) - r0
    do k = 0,11
        if (0 <= k .and. k <= 4) then
            xk = (r0/sqrt(5.0))*2.0*cos(2.0*k*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin(2.0*k*pi/5.0)
            zk = (r0/sqrt(5.0))
        else if (5 <= k .and. k <= 9) then
            xk = (r0/sqrt(5.0))*2.0*cos((2.0*(k-5)-1.0)*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin((2.0*(k-5)-1.0)*pi/5.0)
            zk =-(r0/sqrt(5.0))
        else if (k == 10) then
            xk = 0
            yk = 0
            zk = r0
        else
            xk = 0
            yk = 0
            zk = -r0
        end if
        val = val - A*exp( -( (x - xk)**2  + (y - yk)**2 + (z - zk)**2 )/(sg**2) )
    end do
  end function par_test_poisson_driver_popcorn_fun
  
end module par_test_transient_poisson_driver_names
