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
module par_pb_bddc_maxwell_driver_names
  use fempar_names
  use par_pb_bddc_maxwell_params_names
  use par_pb_bddc_maxwell_discrete_integration_names
  use par_pb_bddc_maxwell_conditions_names
  use par_pb_bddc_maxwell_analytical_functions_names
# include "debug.i90"

  implicit none
  private

  type par_pb_bddc_maxwell_fe_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(par_pb_bddc_maxwell_params_t)   :: test_params
     type(ParameterList_t), pointer       :: parameter_list

     ! Cells and lower dimension objects container
     type(par_triangulation_t)             :: triangulation
     integer(ip), allocatable              :: cells_set_id(:) 

     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                          :: fe_space 
     type(p_reference_fe_t), allocatable           :: reference_fes(:) 
     type(Hcurl_l1_coarse_fe_handler_t)            :: coarse_fe_handler 
     type(p_l1_coarse_fe_handler_t), allocatable   :: coarse_fe_handlers(:)
     type(maxwell_CG_discrete_integration_t)       :: maxwell_integration
     type(maxwell_conditions_t)                    :: maxwell_conditions
     ! Analytical functions describing parameters
     type(par_pb_bddc_maxwell_analytical_functions_t) :: maxwell_analytical_functions
     type(resistivity_holder_t), allocatable          :: resistivity_holder(:) 
     type(resistivity_function_white_t)               :: resistivity_white 
     type(resistivity_function_black_t)               :: resistivity_black 
     type(permeability_holder_t), allocatable         :: permeability_holder(:) 
     type(permeability_function_white_t)              :: permeability_white 
     type(permeability_function_black_t)              :: permeability_black 
     real(rp), allocatable                            :: average_permeability(:) 
     real(rp), allocatable                            :: average_resistivity(:) 

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)            :: fe_affine_operator
     type(fe_affine_operator_t)            :: fe_affine_prec_operator 

     !#ifdef ENABLE_MKL     
     ! MLBDDC preconditioner
     type(mlbddc_t)                            :: mlbddc
     !#endif  

     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)           :: iterative_linear_solver

     ! maxwell problem solution FE function
     type(fe_function_t)                   :: solution

     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                    :: par_environment

     ! Timers
     type(timer_t) :: timer_triangulation
     type(timer_t) :: timer_fe_space
     type(timer_t) :: timer_assemply
     type(timer_t) :: timer_solver_setup
     type(timer_t) :: timer_change_basis_setup 
     type(timer_t) :: timer_solver_run

   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_timers
     procedure                  :: report_timers
     procedure                  :: free_timers
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: set_cells_set_id 
     procedure        , private :: setup_analytical_functions 
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: compute_average_parameter_values
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: print_info 
     procedure                  :: run_simulation
     procedure        , private :: free
     procedure                  :: free_command_line_parameters
     procedure                  :: free_environment
  end type par_pb_bddc_maxwell_fe_driver_t

  ! Types
  public :: par_pb_bddc_maxwell_fe_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    call this%test_params%create()
    this%parameter_list => this%test_params%get_values()
  end subroutine parse_command_line_parameters

  !========================================================================================
  subroutine setup_timers(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(execution_context_t), pointer :: w_context
    w_context => this%par_environment%get_w_context()
    call this%timer_triangulation%create(w_context,"SETUP TRIANGULATION")
    call this%timer_fe_space%create(     w_context,"SETUP FE SPACE")
    call this%timer_assemply%create(     w_context,"FE INTEGRATION AND ASSEMBLY")
    call this%timer_change_basis_setup%create( w_context, "CHANGE BASIS SETUP")
    call this%timer_solver_setup%create( w_context,"SETUP SOLVER AND PRECONDITIONER")
    call this%timer_solver_run%create(   w_context,"SOLVER RUN")
  end subroutine setup_timers

  !========================================================================================
  subroutine report_timers(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this

    call this%timer_triangulation%report(.true.)
    call this%timer_fe_space%report(.true.)
    call this%timer_assemply%report(.true.)
    call this%timer_change_basis_setup%report(.true.)
    call this%timer_solver_setup%report(.true.)
    call this%timer_solver_run%report(.true.)

  end subroutine report_timers

  !========================================================================================
  subroutine free_timers(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_change_basis_setup%free() 
    call this%timer_solver_setup%free()
    call this%timer_solver_run%free()
  end subroutine free_timers

  !========================================================================================
  subroutine setup_environment(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    istat = this%parameter_list%set(key = execution_context_key, value = mpi_context) ; check(istat==0)
    call this%par_environment%create (this%parameter_list)
  end subroutine setup_environment

  subroutine setup_triangulation(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(vef_iterator_t), allocatable :: vef
    real(rp)                           :: domain(6)
    integer(ip)                        :: istat 

    ! Create a structured mesh 
    call this%triangulation%create(this%parameter_list, this%par_environment)

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

    call this%set_cells_set_id() 
    call this%triangulation%setup_coarse_triangulation()
  end subroutine setup_triangulation

  subroutine set_cells_set_id( this ) 
    implicit none 
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(cell_iterator_t), allocatable :: cell 
    type(point_t), allocatable    :: cell_coordinates(:)

    logical, allocatable   :: is_grav_center_within_channel(:)
    real(rp)               :: channel_size 
    real(rp)    :: grav_center(3) 
    integer(ip) :: nparts_x_dir 
    real(rp)    :: nparts  
    integer(ip) :: i, ndime
    integer(ip) :: idime, inode 
    integer(ip) :: ijk(3), aux 
    integer(ip) :: istat, dummy_val 
    real(rp)    :: contrast
    real(rp)    :: resistivity, permeability 
    real(rp)    :: resistivity_max, resistivity_min 
    real(rp)    :: permeability_max, permeability_min 
    logical     :: min_values_initialized

    if ( .not. this%par_environment%am_i_l1_task() ) return 

    ! Determine indices i,j,k for the current subdomain 
    ijk = 0
    nparts = (real(this%par_environment%get_num_tasks()-1,rp))**(1.0_rp/(real(this%triangulation%get_num_dims(),rp)))   
    nparts_x_dir = nint(nparts,ip)
    aux = this%par_environment%get_l1_rank()
    do i = 1,this%triangulation%get_num_dims()-1
       ijk(i) = mod(aux, nparts_x_dir)
       aux = aux/nparts_x_dir
    end do
    ijk(this%triangulation%get_num_dims()) = aux
    ijk = ijk+1

    call this%triangulation%create_cell_iterator(cell)
    call memalloc( this%triangulation%get_num_local_cells(), this%cells_set_id, __FILE__, __LINE__ )
    allocate(cell_coordinates( cell%get_num_nodes() ) , stat=istat); check(istat==0)
    if ( this%test_params%get_materials_distribution_case() == channels ) then 
       call memalloc( this%triangulation%get_num_dims(), is_grav_center_within_channel, __FILE__, __LINE__ )
    end if

    ! rPB-BDDC needs a first cell loop to determine the minimum value of the coefficient 
    if ( this%test_params%get_materials_distribution_case() == heterogeneous ) then  
       min_values_initialized=.false. 
       do while ( .not. cell%has_finished() ) 
          if ( cell%is_local() ) then 
             call cell%get_nodes_coordinates(cell_coordinates)

             do inode=1, cell%get_num_nodes() 
                call this%resistivity_holder(1)%p%get_value(cell_coordinates(inode), resistivity)
                call this%permeability_holder(1)%p%get_value(cell_coordinates(inode), permeability)
                if ( .not. min_values_initialized ) then 
                resistivity_min  = resistivity 
                permeability_min = permeability
                min_values_initialized = .true. 
                else 
                resistivity_min  = min(resistivity_min,  resistivity) 
                permeability_min = min(permeability_min, permeability) 
                end if 
             end do
          end if
          call cell%next()
       end do
          ! Init cell iterator 
          call cell%first() 
         
    end if

    this%cells_set_id = 0
    do while ( .not. cell%has_finished() )
       if ( cell%is_local() ) then 
          select case ( this%test_params%get_materials_distribution_case() )
          case ( homogeneous ) 
             this%cells_set_id(cell%get_gid()) = 0
          case ( checkerboard )  
             if ( mod( sum(ijk),2 ) == 0 ) then 
                this%cells_set_id(cell%get_gid()) = WHITE 
             else 
                this%cells_set_id(cell%get_gid()) = BLACK 
             end if
          case ( channels ) 
             call cell%get_nodes_coordinates(cell_coordinates)
             grav_center = 0
             do inode=1, cell%get_num_nodes() 
                do idime = 1, this%triangulation%get_num_dims()
                   grav_center(idime) = grav_center(idime) + cell_coordinates(inode)%get(idime)/cell%get_num_nodes() 
                end do
             end do

             is_grav_center_within_channel = .false. 
             do idime=1, this%triangulation%get_num_dims()
                channel_size = 1.0_rp/real(nparts_x_dir,rp)*this%test_params%get_channels_ratio()
                if ( grav_center(idime) >= (real(ijk(idime)-1,rp))/real(nparts_x_dir,rp) .and. & 
                     grav_center(idime) <= (real(ijk(idime)-1,rp))/real(nparts_x_dir,rp) + channel_size ) then  
                   is_grav_center_within_channel(idime) = .true. 
                end if
             end do

             if ( count(is_grav_center_within_channel) >= this%triangulation%get_num_dims()-1 ) then 
                this%cells_set_id(cell%get_gid()) = WHITE
             else 
                this%cells_set_id(cell%get_gid()) = BLACK 
             end if

          case ( heterogeneous ) 
             ! Extract coordinates, evaluate resistivity/permeability
             call cell%get_nodes_coordinates(cell_coordinates)

             permeability_max = 0.0_rp 
             resistivity_max  = 0.0_rp 
             do inode=1, cell%get_num_nodes() 
                call this%resistivity_holder(1)%p%get_value(cell_coordinates(inode), resistivity)
                call this%permeability_holder(1)%p%get_value(cell_coordinates(inode), permeability)
                permeability_max = max( permeability_max, permeability ) 
                resistivity_max  = max( resistivity_max, resistivity ) 
             end do
             
             contrast=(resistivity_max)/(resistivity_min)
             massert(this%test_params%get_rpb_bddc_threshold()>1.0_rp, 'Not valid Relaxed PB-BDDC threshold') 
             this%cells_set_id(cell%get_gid()) = floor( log(contrast)/log(this%test_params%get_rpb_bddc_threshold()) )
          case DEFAULT 
             massert( .false., 'Materials distribution case not valid')
          end select
       end if
       call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
    if ( this%test_params%get_materials_distribution_case() == channels ) then 
       if (allocated(cell_coordinates)) deallocate(cell_coordinates, stat=istat); check(istat==0)
       if (allocated(is_grav_center_within_channel)) call memfree(is_grav_center_within_channel, __FILE__, __LINE__ )
    end if

    call this%triangulation%fill_cells_set(this%cells_set_id) 
  end subroutine set_cells_set_id

  subroutine setup_analytical_functions(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat

    ! Initialize resistivity values for different materials 
    call this%resistivity_black%set_value(this%test_params%get_resistivity_black())
    call this%resistivity_white%set_value(this%test_params%get_resistivity_white()) 
    call this%permeability_black%set_value(this%test_params%get_permeability_black())
    call this%permeability_white%set_value(this%test_params%get_permeability_white()) 

    select case ( this%test_params%get_materials_distribution_case() ) 
    case ( checkerboard, channels ) 
    allocate(this%resistivity_holder(2), stat=istat)
    allocate(this%permeability_holder(2), stat=istat) 

    this%resistivity_holder(1)%p => this%resistivity_black 
    this%resistivity_holder(2)%p => this%resistivity_white
    this%permeability_holder(1)%p => this%permeability_black 
    this%permeability_holder(2)%p => this%permeability_white 
    case ( homogeneous, heterogeneous ) 
    
    allocate(this%resistivity_holder(1), stat=istat)
    allocate(this%permeability_holder(1), stat=istat) 
    this%resistivity_holder(1)%p => this%resistivity_white 
    this%permeability_holder(1)%p => this%permeability_white
    
    case DEFAULT 
    massert(.false., 'Not valid material distribution case') 
    end select 

    call this%maxwell_analytical_functions%set_resistivity( this%resistivity_holder )
    call this%maxwell_analytical_functions%set_permeability( this%permeability_holder ) 

    call this%maxwell_analytical_functions%set_num_dims(this%triangulation%get_num_dims())

  end subroutine setup_analytical_functions

  subroutine setup_reference_fes(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    integer(ip) :: istat
    class(cell_iterator_t), allocatable    :: cell
    class(reference_fe_t), pointer         :: reference_fe_geo

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)

    if ( this%par_environment%am_i_l1_task() ) then
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

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout), target :: this
    integer(ip) :: istat 

    allocate(this%coarse_fe_handlers(1), stat=istat); check(istat==0)
    this%coarse_fe_handlers(1)%p => this%coarse_fe_handler
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this

    call this%maxwell_conditions%set_num_dims(this%triangulation%get_num_dims())
    call this%maxwell_conditions%set_boundary_function_Hx(this%maxwell_analytical_functions%get_boundary_function_Hx())
    call this%maxwell_conditions%set_boundary_function_Hy(this%maxwell_analytical_functions%get_boundary_function_Hy())
    if ( this%triangulation%get_num_dims() == 3) then 
       call this%maxwell_conditions%set_boundary_function_Hz(this%maxwell_analytical_functions%get_boundary_function_Hz())
    end if

    call this%fe_space%create( triangulation       = this%triangulation,      &
         reference_fes       = this%reference_fes,      &
         coarse_fe_handlers  = this%coarse_fe_handlers, & 
         conditions          = this%maxwell_conditions  )

    call this%fe_space%set_up_cell_integration()
    call this%fe_space%set_up_facet_integration()   
  end subroutine setup_fe_space

  subroutine setup_system (this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(matrix_t), pointer :: matrix 

    call this%maxwell_integration%set_analytical_functions(this%maxwell_analytical_functions)
    call this%maxwell_integration%set_materials_distribution_case(this%test_params%get_materials_distribution_case())
                                                                  
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
         diagonal_blocks_symmetric_storage = [ .true. ], &
         diagonal_blocks_symmetric         = [ .true. ], &
         diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
         fe_space                          = this%fe_space, &
         discrete_integration              = this%maxwell_integration )

    if ( this%test_params%get_boundary_mass_trick() ) then 
       call this%fe_affine_prec_operator%create ( sparse_matrix_storage_format      = csr_format, &
            diagonal_blocks_symmetric_storage = [ .true. ], &
            diagonal_blocks_symmetric         = [ .true. ], &
            diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
            fe_space                          = this%fe_space, &
            discrete_integration              = this%maxwell_integration )
    end if

    ! Set-up solution with Dirichlet boundary conditions
    call this%solution%create(this%fe_space)
    call this%fe_space%interpolate_dirichlet_values(this%solution) 
    call this%maxwell_integration%set_fe_function(this%solution)

  end subroutine setup_system

  subroutine setup_solver (this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    type(parameterlist_t) :: parameter_list
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    integer(ip) :: FPLError
    integer(ip) :: ilev
    integer(ip) :: iparm(64)
    class(matrix_t), pointer :: matrix 

    call this%timer_change_basis_setup%start()
    if ( this%par_environment%am_i_l1_task() ) then       
       ! Compute average parameters to be sent to the preconditioner for weights computation 
       call this%compute_average_parameter_values() 

       matrix => this%fe_affine_operator%get_matrix() 
       select type ( matrix ) 
          class is (par_sparse_matrix_t) 
          call this%coarse_fe_handler%create( 1, this%fe_space, matrix, this%parameter_list, this%average_permeability, this%average_resistivity )
       end select
    end if
    call this%timer_change_basis_setup%stop()

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
    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    if ( this%test_params%get_boundary_mass_trick() ) then 
       call this%mlbddc%create(this%fe_affine_prec_operator, this%parameter_list)
    else 
       call this%mlbddc%create(this%fe_affine_operator, this%parameter_list)
    end if
    call this%mlbddc%symbolic_setup()
    call this%mlbddc%numerical_setup()  

    call parameter_list%init()
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-6_rp)
    FPLError = parameter_list%set(key = ils_max_num_iterations, value = 1000)
    assert(FPLError == 0)

    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)

    call this%iterative_linear_solver%set_operators(this%fe_affine_operator%get_tangent(), this%mlbddc) 
    call parameter_list%free()
  end subroutine setup_solver

  ! -----------------------------------------------------------------------------------------------
  subroutine compute_average_parameter_values(this)
    implicit none 
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this  
    class(fe_cell_iterator_t), allocatable :: fe
    ! Integration loop 
    type(quadrature_t)       , pointer     :: quad
    integer(ip)                            :: qpoin, num_quad_points
    type(point_t)            , pointer     :: quad_coords(:)
    real(rp)                               :: factor 
    real(rp), allocatable                  :: set_id_volume(:)
    real(rp), allocatable                  :: permeability(:)
    real(rp), allocatable                  :: resistivity(:) 
    integer(ip) :: max_cell_set_id, set_id, material_id 
    integer(ip) :: istat 

    max_cell_set_id = maxval(this%cells_set_id)+1 ! 1-based arrays  
    call memalloc( max_cell_set_id, this%average_permeability, __FILE__, __LINE__ )
    call memalloc( max_cell_set_id, this%average_resistivity, __FILE__, __LINE__ ) 
    call memalloc( max_cell_set_id, set_id_volume, __FILE__, __LINE__ )
    this%average_permeability = 0.0_rp 
    this%average_resistivity  = 0.0_rp 
    set_id_volume = 0.0_rp 

    ! Integrate structures needed 
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
    call fe%update_integration()
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    quad_coords      => fe%get_quadrature_points_coordinates()
    call memalloc( num_quad_points, resistivity, __FILE__, __LINE__ )
    call memalloc( num_quad_points, permeability, __FILE__, __LINE__ )

    ! Loop over elements
    do while ( .not. fe%has_finished())
       if ( fe%is_local() ) then  
          call fe%update_integration()
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Evaluate parameters on the cell 
          select case ( this%test_params%get_materials_distribution_case() )
          case ( homogeneous, heterogeneous )
          material_id = 1
          case ( checkerboard, channels ) 
          material_id = 1 + fe%get_set_id()
          end select 
          set_id = 1 + fe%get_set_id() ! 1-based arrays 
          call this%resistivity_holder(material_id)%p%get_values_set(quad_coords, resistivity)
          call this%permeability_holder(material_id)%p%get_values_set(quad_coords, permeability)  

          ! Integrate cell contribution to averages 
          do qpoin=1, num_quad_points
             factor = fe%get_det_jacobian(qpoin) * quad%get_weight(qpoin) 						         
             ! Average magnetic field 
             this%average_permeability(set_id) = this%average_permeability(set_id) + permeability(qpoin)*factor 
             this%average_resistivity(set_id)  = this%average_resistivity(set_id)  + resistivity(qpoin)*factor
             set_id_volume(set_id)             = set_id_volume(set_id) + factor
          end do
       end if
       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)

    do set_id=1, max_cell_set_id
       if ( set_id_volume(set_id) == 0.0_rp ) cycle
       this%average_permeability(set_id) = this%average_permeability(set_id)/set_id_volume(set_id)
       this%average_resistivity(set_id)  = this%average_resistivity(set_id)/set_id_volume(set_id)
    end do

    call memfree ( resistivity, __FILE__, __LINE__ ) 
    call memfree ( permeability, __FILE__, __LINE__ ) 
    call memfree(set_id_volume, __FILE__, __LINE__)
  end subroutine compute_average_parameter_values

  subroutine assemble_system (this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs

    call this%maxwell_integration%set_integration_type( problem ) 
    call this%fe_affine_operator%compute()

    if ( this%test_params%get_boundary_mass_trick() ) then 
       call this%maxwell_integration%set_integration_type( preconditioner ) 
       call this%fe_affine_prec_operator%compute() 
    end if

  end subroutine assemble_system

  subroutine solve_system(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_free_dof_values()
    call this%iterative_linear_solver%apply(this%fe_affine_operator%get_translation(), dof_values)

  end subroutine solve_system

  subroutine check_solution(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    type(error_norms_vector_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, hcurl_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: tol

    call error_norm%create(this%fe_space,1)    
    mean = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    hcurl_s = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, hcurl_seminorm)
    h1_s = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%maxwell_analytical_functions%get_solution_function(), this%solution, w1infty_norm)  
    if ( this%par_environment%am_i_l1_root() ) then
       tol=1.0e-3_rp 
       write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < tol )
       write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < tol )
       write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < tol )
       write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < tol )
       write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < tol )
       write(*,'(a20,e32.25)') 'hcurl_seminorm:', hcurl_s; check ( linfty < tol)
       write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < tol )
       write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < tol )
       write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < tol )
       write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < tol )
       write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < tol )
       write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < tol )
    end if
    call error_norm%free()
  end subroutine check_solution

  subroutine write_solution(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(in)    :: this
    type(output_handler_t)                             :: oh	
    real(rp), allocatable                              :: set_id_cell_vector(:)
    real(rp), allocatable                              :: set_id_rank(:)
    integer(ip)                                        :: i, istat

    if ( this%par_environment%am_i_l1_task() ) then
       if( this%test_params%get_write_solution() .and. (this%test_params%get_reference_fe_order()==1) ) then
          call build_set_ids()
          call oh%create()
          call oh%attach_fe_space(this%fe_space)
          call oh%add_fe_function(this%solution, 1, 'u')
          call oh%add_fe_function(this%solution, 1, 'curl(u)', curl_diff_operator)
          call oh%add_cell_vector(set_id_rank, 'rank')
          call oh%add_cell_vector(set_id_cell_vector, 'set_id')
          call oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
          call oh%write()
          call oh%close()
          call oh%free()
          call free_set_ids()
       end if
    end if

  contains 
    subroutine build_set_ids()
      implicit none 
      class(cell_iterator_t), allocatable :: cell 

      ! Allocate set ids 
      call memalloc(this%triangulation%get_num_local_cells(), set_id_cell_vector, __FILE__, __LINE__)
      call memalloc(this%triangulation%get_num_local_cells(), set_id_rank, __FILE__, __LINE__)

      ! Fill Set ids 
      call this%triangulation%create_cell_iterator(cell) 
      do while ( .not. cell%has_finished() ) 
         if ( cell%is_local() ) then 
            set_id_cell_vector(cell%get_gid()) = cell%get_set_id()
            set_id_rank(cell%get_gid())        = this%par_environment%get_l1_rank() + 1
         end if
         call cell%next() 
      enddo
      call this%triangulation%free_cell_iterator(cell) 
    end subroutine build_set_ids

    subroutine free_set_ids()
      call memfree(set_id_cell_vector, __FILE__, __LINE__)
      call memfree(set_id_rank, __FILE__, __LINE__)
    end subroutine free_set_ids

  end subroutine write_solution

  !========================================================================================
  subroutine print_info (this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this

    integer(ip) :: num_sub_domains
    real(rp)    :: num_total_cells
    real(rp)    :: num_dofs
    integer(ip) :: num_coarse_dofs

    class(environment_t), pointer :: environment
    class(coarse_fe_space_t), pointer :: coarse_fe_space

    environment => this%fe_space%get_environment()

    if (environment%am_i_l1_task()) then
       num_total_cells  = real(this%triangulation%get_num_local_cells(),kind=rp)
       num_dofs         = real(this%fe_space%get_field_num_dofs(1),kind=rp)
       call environment%l1_sum(num_total_cells )
       call environment%l1_sum(num_dofs        )
    end if

    if (environment%get_l1_rank() == 0) then
       num_sub_domains = environment%get_l1_size()
       write(*,'(a,i22)') 'num_sub_domains:          ', num_sub_domains
       write(*,'(a,i22)') 'num_total_cells:          ', nint(num_total_cells , kind=ip )
       write(*,'(a,i22)') 'num_dofs (sub-assembled): ', nint(num_dofs        , kind=ip )
    end if

    if (environment%am_i_lgt1_task()) then
       coarse_fe_space => this%fe_space%get_coarse_fe_space()
       num_coarse_dofs = coarse_fe_space%get_field_num_dofs(1)
       write(*,'(a,i22)') 'num_coarse_dofs:  ', num_coarse_dofs
    end if

  end subroutine print_info

  subroutine run_simulation(this) 
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this

    call this%timer_triangulation%start()
    call this%setup_analytical_functions()
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

    call this%print_info()
    call this%write_solution()
    ! call this%check_solution()    
    call this%free()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    integer(ip) :: i, istat

    call this%solution%free()
    call this%mlbddc%free()  
    call this%iterative_linear_solver%free()
    call this%fe_affine_operator%free()
    if ( this%test_params%get_boundary_mass_trick() ) then 
       call this%fe_affine_prec_operator%free()
    end if

    if (allocated(this%coarse_fe_handlers) ) then 
       call this%coarse_fe_handler%free() 
       deallocate(this%coarse_fe_handlers, stat=istat); check(istat==0) 
    end if

    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          call this%reference_fes(i)%free()
       end do
       deallocate(this%reference_fes, stat=istat)
       check(istat==0)
    end if

    if ( allocated(this%resistivity_holder) ) then
       deallocate(this%resistivity_holder, stat=istat)
       check(istat==0)
    end if

    if ( allocated(this%permeability_holder) ) then
       deallocate(this%permeability_holder, stat=istat)
       check(istat==0)
    end if

    call this%triangulation%free()
    if (allocated(this%cells_set_id) )        call memfree( this%cells_set_id, __FILE__ ,__LINE__ )
    if (allocated(this%average_resistivity))  call memfree( this%average_resistivity, __FILE__, __LINE__ ) 
    if (allocated(this%average_permeability)) call memfree( this%average_permeability, __FILE__, __LINE__ ) 

  end subroutine free

  !========================================================================================
  subroutine free_environment(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    call this%par_environment%free()
  end subroutine free_environment

  !========================================================================================
  subroutine free_command_line_parameters(this)
    implicit none
    class(par_pb_bddc_maxwell_fe_driver_t), intent(inout) :: this
    call this%test_params%free()
  end subroutine free_command_line_parameters

end module par_pb_bddc_maxwell_driver_names
