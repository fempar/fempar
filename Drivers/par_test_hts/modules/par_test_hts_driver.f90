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
  use hts_output_handler_field_generator_names 
# include "debug.i90"

  implicit none
  private

  type par_test_hts_fe_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(par_test_hts_params_t)      :: test_params
     type(ParameterList_t), pointer   :: parameter_list

     ! Cells and lower dimension objects container
     type(par_triangulation_t)       :: triangulation
     integer(ip), allocatable        :: cells_set_ids(:) 

     ! Discrete weak problem integration-related data type instances 
     type(par_fe_space_t)                              :: fe_space 
     type(p_reference_fe_t)             , allocatable  :: reference_fes(:) 
     type(Hcurl_l1_coarse_fe_handler_t)                :: coarse_fe_handler
     type(p_l1_coarse_fe_handler_t)     , allocatable  :: coarse_fe_handlers(:)
     type(hts_discrete_integration_t)                  :: hts_integration

     ! Type of problem 
     type(maxwell_conditions_t)                    :: maxwell_conditions 
     type(maxwell_analytical_functions_t)          :: maxwell_analytical_functions 
     type(hts_conditions_t)                        :: hts_conditions
     type(hts_analytical_functions_t)              :: hts_analytical_functions

     ! Parameter-based solver 
     real(rp), allocatable                            :: average_permeability(:) 
     real(rp), allocatable                            :: average_resistivity(:) 

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(iterative_linear_solver_t)             :: linear_solver
     type(fe_operator_t)                         :: fe_operator
     type(nonlinear_solver_t)                    :: nonlinear_solver 
     type(mlbddc_t)                              :: mlbddc

     ! maxwell problem solution FE function
     type(fe_function_t)                   :: H_current
     type(fe_function_t)                   :: H_previous 

     ! Environment required for fe_affine_operator + vtk_handler
     type(environment_t)                    :: environment

     ! Time integration 
     type(theta_method_t)                   :: theta_method 

     ! Output handler 
     type(resistivity_field_generator_t)    :: resistivity_field_generator 
     real(rp), allocatable                  :: set_id_cell_vector(:)
     real(rp), allocatable                  :: rank_cell_vector(:) 
     type(output_handler_t)                 :: oh	

     ! Timers
     type(timer_t) :: timer_triangulation
     type(timer_t) :: timer_fe_space
     type(timer_t) :: timer_assemply
     type(timer_t) :: timer_setup_system 
     type(timer_t) :: timer_solver_total
     type(timer_t) :: timer_write_sol
     type(timer_t) :: timer_run_simulation
     type(timer_t) :: timer_find_subsets 
     type(timer_t) :: timer_setup_coarse_fe_handler
     
   contains
     procedure                  :: parse_command_line_parameters
     procedure                  :: setup_timers
     procedure                  :: report_timers
     procedure                  :: free_timers
     procedure                  :: setup_environment
     procedure        , private :: setup_triangulation
     procedure        , private :: set_material_cells_set_id
     procedure        , private :: set_pb_cells_set_id 
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_coarse_fe_handlers
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_discrete_integration
     procedure        , private :: setup_operators
     procedure        , private :: setup_solver_after_pb_partition 
     procedure        , private :: setup_theta_method 
     procedure        , private :: compute_average_parameter_values
     procedure        , private :: setup_system
     procedure        , private :: solve_problem 
     procedure        , private :: compute_hysteresis_data 
     procedure        , private :: check_solution
     procedure        , private :: initialize_output
     procedure        , private :: write_time_step_solution 
     procedure        , private :: finalize_output 
     procedure                  :: run_simulation
     procedure        , private :: print_info 
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

  subroutine setup_timers(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(execution_context_t), pointer :: w_context
    w_context => this%environment%get_w_context()
    call this%timer_triangulation%create(w_context,"SETUP TRIANGULATION")
    call this%timer_fe_space%create(w_context,"SETUP FE SPACE")
    call this%timer_assemply%create(w_context,"FE INTEGRATION AND ASSEMBLY")
    call this%timer_setup_system%create(w_context,"OPERATORS SETUP")
    call this%timer_solver_total%create(w_context,"SOLVER SETUP AND RUN")
    call this%timer_write_sol%create(w_context,"WRITE SOLUTION")
    call this%timer_run_simulation%create(w_context,"RUN SIMULATION")
    call this%timer_find_subsets%create(w_context, "FIND SUBSETS")  
    call this%timer_setup_coarse_fe_handler%create(w_context, "SETUP COARSE FE HANDLER")  
  end subroutine setup_timers

  subroutine report_timers(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%report(.true.)
    call this%timer_fe_space%report(.false.)
    call this%timer_assemply%report(.false.)
    call this%timer_setup_system%report(.false.)
    call this%timer_solver_total%report(.true.)
    call this%timer_write_sol%report(.false.)
    call this%timer_run_simulation%report(.false.)
    if (this%environment%get_l1_rank() == 0) then
       write(*,*)
    end if
  end subroutine report_timers

  subroutine free_timers(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    call this%timer_triangulation%free()
    call this%timer_fe_space%free()
    call this%timer_assemply%free()
    call this%timer_setup_system%free() 
    call this%timer_solver_total%free()
    call this%timer_write_sol%free()
    call this%timer_run_simulation%free()
    call this%timer_find_subsets%free() 
    call this%timer_setup_coarse_fe_handler%free() 
  end subroutine free_timers

  subroutine setup_environment(this, world_context)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(execution_context_t)     , intent(in)    :: world_context
    integer(ip) :: istat
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       istat = this%parameter_list%set(key = environment_type_key, value = structured) ; check(istat==0)
    else
       istat = this%parameter_list%set(key = environment_type_key, value = unstructured) ; check(istat==0)
    end if
    call this%environment%create(world_context, this%parameter_list)
  end subroutine setup_environment
  
  subroutine setup_triangulation(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(vef_iterator_t), allocatable :: vef
    real(rp)                           :: domain(6)
    integer(ip)                        :: istat 
    integer(ip)                        :: i 
    class(environment_t), pointer      :: environment

    ! Create a structured mesh with a custom domain 
    domain   = this%test_params%get_domain_limits()
    istat    = this%parameter_list%set(key = hex_mesh_domain_limits_key , value = domain); check(istat==0)
    call this%triangulation%create(this%environment, this%parameter_list)

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

    call this%set_material_cells_set_id() 
    call this%triangulation%setup_coarse_triangulation()
  end subroutine setup_triangulation

  subroutine set_material_cells_set_id(this) 
    implicit none 
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    ! Locals 
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

       hts_size = this%test_params%get_hts_domain_length() 
       hts_lx = hts_size(0) 
       hts_ly = hts_size(1) 
       hts_lz = hts_size(2) 

       domain = this%test_params%get_domain_limits()
       lx = domain(2)-domain(1) 
       ly = domain(4)-domain(3) 
       lz = domain(6)-domain(5) 

       ! Assign subset_id to different cells for the created structured mesh 
       call memalloc(this%triangulation%get_num_local_cells(), this%cells_set_ids, __FILE__, __LINE__)
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
             if ( this%triangulation%get_num_dims() == 2 ) then 
                if ( ( ((lx-hts_lx)/2.0_rp<=cx) .and. (cx<=(lx+hts_lx)/2.0_rp) ) .and. &
                     ( ((ly-hts_ly)/2.0_rp<=cy) .and. (cy<=(ly+hts_ly)/2.0_rp) ) ) then
                   this%cells_set_ids( cell%get_gid() ) = hts  
                else 
                   this%cells_set_ids( cell%get_gid() ) = air
                end if
             elseif ( this%triangulation%get_num_dims() == 3) then 
                if ( ( ((lx-hts_lx)/2.0_rp<=cx) .and. (cx<=(lx+hts_lx)/2.0_rp) ) .and. &
                     ( ((ly-hts_ly)/2.0_rp<=cy) .and. (cy<=(ly+hts_ly)/2.0_rp) ) .and. & 
                     ( ((lz-hts_lz)/2.0_rp<=cz) .and. (cz<=(lz+hts_lz)/2.0_rp) )) then
                   this%cells_set_ids( cell%get_gid() ) = hts  
                else 
                   this%cells_set_ids( cell%get_gid() ) = air
                end if
             end if
          end if
          call cell%next() 
       end do

       call this%triangulation%fill_cells_set(this%cells_set_ids)   
       deallocate(cell_coordinates, stat=istat); check(istat==0) 
       call this%triangulation%free_cell_iterator(cell)
    end if

  end subroutine set_material_cells_set_id

  subroutine set_pb_cells_set_id( this ) 
    implicit none 
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(fe_cell_iterator_t), allocatable :: fe
    type(point_t)      , pointer     :: quad_coords(:)
    type(quadrature_t) , pointer     :: quad
    integer(ip)                      :: qpoin, num_quad_points
    integer(ip) :: i, ndime
    integer(ip) :: istat, dummy_val 
    real(rp)    :: contrast, cell_contrast
    ! Compute resistivity 
    type(fe_cell_function_vector_t)        :: fe_cell_function
    type(vector_field_t)                   :: H_value, H_curl
    type(tensor_field_t)                   :: resistivity_tensor 
    real(rp)                               :: resistivity
    real(rp)                               :: cell_resistivity_max, cell_resistivity_min 
    real(rp)                               :: sbd_resistivity_min, sbd_resistivity_max 
    real(rp)                               :: global_resistivity_min 
    integer(ip)                            :: num_cells_sent_to_coarse_solver  
    integer(ip)                            :: num_sets_id 
    
    if ( .not. this%environment%am_i_l1_task() ) return 

    ! Integrate structures needed 
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
    call fe%update_integration()
    call fe_cell_function%create(this%fe_space, 1)
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    quad_coords      => fe%get_quadrature_points_coordinates()
    
    ! rPB-BDDC needs a first cell loop to determine the minimum value of the coefficient 
    sbd_resistivity_min = this%test_params%get_air_resistivity() 
    sbd_resistivity_max = 0.0_rp 
    do while ( .not. fe%has_finished() ) 
       if ( fe%is_local() .and. fe%get_set_id() > AIR) then 

       call fe%update_integration()
          call fe_cell_function%update(fe, this%H_current) 
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Integrate cell contribution to averages 
          do qpoin=1, num_quad_points
            ! Evaluate parameters on the cell 
             call fe_cell_function%get_value(qpoin, H_value)
             call fe_cell_function%compute_curl(qpoin, H_curl)
             resistivity_tensor  = this%hts_integration%compute_resistivity( H_value, H_curl, fe%get_set_id() ) 
             resistivity = max( resistivity_tensor%get(1,1), resistivity_tensor%get(2,2), resistivity_tensor%get(3,3) ) 
             sbd_resistivity_min = min(sbd_resistivity_min, resistivity)
             sbd_resistivity_max = max(sbd_resistivity_max, resistivity)
          end do

       end if
       call fe%next()
    end do
    
    ! Number of subsets that will arise
    num_cells_sent_to_coarse_solver = 0 
    global_resistivity_min = 1.0_rp/sbd_resistivity_min 
    call this%environment%l1_max(global_resistivity_min)
    global_resistivity_min = 1.0_rp/global_resistivity_min 
    contrast    = sbd_resistivity_max/global_resistivity_min 
    num_sets_id = floor( log(contrast)/log(this%test_params%get_rpb_bddc_threshold()) )

    ! Init fe iterator 
    call fe%first()
    do while ( .not. fe%has_finished() )
       if ( fe%is_local() .and. fe%get_set_id() > AIR) then 
          ! Extract coordinates, evaluate resistivity/permeability
          call fe%update_integration()
          call fe_cell_function%update(fe, this%H_current) 
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Integrate cell contribution to averages 
          cell_resistivity_max  = 0.0_rp
          cell_resistivity_min  = this%test_params%get_air_resistivity()
          do qpoin=1, num_quad_points
            ! Evaluate parameters on the cell 
             call fe_cell_function%get_value(qpoin, H_value)
             call fe_cell_function%compute_curl(qpoin, H_curl)
             resistivity_tensor = this%hts_integration%compute_resistivity( H_value, H_curl, fe%get_set_id() )
             resistivity = max( resistivity_tensor%get(1,1), resistivity_tensor%get(2,2), resistivity_tensor%get(3,3) )
             cell_resistivity_max = max(cell_resistivity_max,  resistivity)
             cell_resistivity_min = min(cell_resistivity_min,  resistivity)
          end do

          contrast=cell_resistivity_max/global_resistivity_min
          cell_contrast=cell_resistivity_max/cell_resistivity_min
          massert(this%test_params%get_rpb_bddc_threshold()>1.0_rp, 'Not valid Relaxed PB-BDDC threshold')
          if ( cell_contrast < this%test_params%get_rpb_bddc_threshold() ) then 
          this%cells_set_ids(fe%get_gid()) = HTS + floor( log(contrast)/log(this%test_params%get_rpb_bddc_threshold()) )
          else 
          num_cells_sent_to_coarse_solver = num_cells_sent_to_coarse_solver +1
          this%cells_set_ids(fe%get_gid()) = HTS + num_sets_id + num_cells_sent_to_coarse_solver 
          end if 
       end if
       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)
    call this%triangulation%fill_cells_set(this%cells_set_ids) 
    
  end subroutine set_pb_cells_set_id
  
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

  subroutine setup_coarse_fe_handlers(this)
    implicit none
    class(par_test_hts_fe_driver_t), target, intent(inout) :: this
    integer(ip) :: istat 

    allocate(this%coarse_fe_handlers(1), stat=istat); check(istat==0)
    this%coarse_fe_handlers(1)%p => this%coarse_fe_handler
  end subroutine setup_coarse_fe_handlers

  subroutine setup_fe_space(this)
    implicit none
    class(par_test_hts_fe_driver_t), target, intent(inout) :: this

    if (this%test_params%get_is_analytical_solution() ) then 
       call this%maxwell_conditions%set_num_dims(this%triangulation%get_num_dims())

       ! Set-up Dirichlet boundary conditions  
       call this%maxwell_analytical_functions%set_nonlinear_exponent( this%test_params%get_nonlinear_exponent() )
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
       call this%hts_analytical_functions%get_parameter_values( H = this%test_params%get_external_magnetic_field_amplitude(),  &
            wH = this%test_params%get_external_magnetic_field_frequency()   ) 
       if ( this%triangulation%get_num_dims() == 3) then 
          call this%hts_conditions%set_boundary_function_Hz(this%hts_analytical_functions%get_boundary_function_Hz())
       end if

       ! Create FE SPACE 
       call this%fe_space%create( triangulation       = this%triangulation,      &
            reference_fes       = this%reference_fes,      &
            coarse_fe_handlers  = this%coarse_fe_handlers, & 
            conditions          = this%hts_conditions  )
    end if

    call this%fe_space%set_up_cell_integration()
  end subroutine setup_fe_space

  subroutine setup_discrete_integration(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    if ( this%test_params%get_is_analytical_solution() ) then 
       call this%hts_integration%set_analytical_functions(this%maxwell_analytical_functions%get_source_term())
    else
       call this%hts_integration%set_analytical_functions(this%hts_analytical_functions%get_source_term())
    end if

    call this%hts_integration%set_fe_functions(this%H_previous, this%H_current)
    call this%hts_integration%set_theta_method(this%theta_method) 
    call this%hts_integration%set_parameter_values( this%test_params )

  end subroutine setup_discrete_integration

  subroutine setup_system (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    ! Initialize previous time step fe_functions 
    call this%H_previous%create(this%fe_space)
    if (this%test_params%get_is_analytical_solution() ) then 
       call this%fe_space%interpolate_vector_function( 1, this%maxwell_analytical_functions%get_solution_function(), &
                                                       this%H_previous,                                              &  
                                                       time=this%theta_method%get_initial_time() ) 
    else 
       call this%fe_space%interpolate_vector_function( 1, this%hts_analytical_functions%get_initial_conditions(), &
                                                       this%H_previous,                                           &  
                                                       time=this%theta_method%get_initial_time() ) 
    end if
    call this%fe_space%interpolate_dirichlet_values( this%H_previous, time=this%theta_method%get_initial_time() ) 

    ! Current time step fe_function 
    call this%H_current%create(this%fe_space)
    if (this%test_params%get_is_analytical_solution() ) then 
       call this%fe_space%interpolate_vector_function( 1, this%maxwell_analytical_functions%get_solution_function(), &
                                                       this%H_current,                                               &  
                                                       time=0.0_rp ) 
    else 
       call this%fe_space%interpolate_vector_function( 1, this%hts_analytical_functions%get_initial_conditions(), &
                                                       this%H_current,                                            &  
                                                       time=this%theta_method%get_current_time() ) 
    end if
    call this%fe_space%interpolate_dirichlet_values( this%H_current, time=this%theta_method%get_current_time() ) 

  end subroutine setup_system

  subroutine setup_operators (this)
    implicit none
    class(par_test_hts_fe_driver_t), target, intent(inout) :: this
    class(l1_coarse_fe_handler_t), pointer :: coarse_fe_handler 
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    type(parameterlist_t) :: linear_pl
    integer(ip) :: ilev
    integer(ip) :: FPLError
    integer(ip) :: field_id
    integer(ip) :: istat
    integer(ip) :: iparm(64)
    class(matrix_t), pointer :: matrix

    ! FE operator
    call this%fe_operator%create ( sparse_matrix_storage_format      = csr_format, &
         &                                diagonal_blocks_symmetric_storage = [ .true. ], &
         &                                diagonal_blocks_symmetric         = [ .true. ], &
         &                                diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
         &                                fe_space                          = this%fe_space, &
         &                                discrete_integration              = this%hts_integration ) 


    ! BDDC preconditioner 
     ! Compute average parameters to be sent to the preconditioner for weights computation 
    call this%compute_average_parameter_values() 
    matrix => this%fe_operator%get_matrix() 
    select type ( matrix ) 
       class is (par_sparse_matrix_t) 
       call this%coarse_fe_handler%create( 1, this%fe_space, matrix, this%parameter_list, this%average_permeability, this%average_resistivity )
       class DEFAULT
       assert(.false.) 
    end select
    
    ! Prepare the internal parameter list of pardiso
    ! See https://software.intel.com/en-us/node/470298 for details
    iparm      = 0 ! Init all entries to zero
    iparm(1)   = 1  ! no solver default
    iparm(2)   = 2  ! fill-in reordering from METIS
    iparm(8)   = 2  ! numbers of iterative refinement steps
    iparm(10)  = 8  ! perturb the pivot elements with 1E-8
    iparm(21)  = 1  ! 1x1 + 2x2 pivots
    ! Customization
    iparm(11)  = 1 ! use scaling (default 0)
    iparm(13)  = 1 ! use maximum weighted matching algorithm (default 0)

    plist => this%parameter_list  
    if ( this%environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if
    do ilev=1, this%environment%get_num_levels()-1
       ! Dirichlet local problems 
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       ! Constrained Neumann problems 
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    call this%mlbddc%create(this%fe_operator, this%parameter_list)     

    ! Linear solver
    call this%linear_solver%create(this%fe_space%get_environment())
    call this%linear_solver%set_type_from_string(cg_name)
    call linear_pl%init()
    FPLError = linear_pl%set(key = ils_rtol, value = this%test_params%get_relative_linear_tolerance()); assert(FPLError == 0)
    FPLError = linear_pl%set(key = ils_max_num_iterations, value = 1000); assert(FPLError == 0)
    call this%linear_solver%set_parameters_from_pl(linear_pl)
    call this%linear_solver%set_operators( this%fe_operator%get_tangent(), this%mlbddc)
    call linear_pl%free()

    ! Nonlinear solver
    call this%nonlinear_solver%create(convergence_criteria = this%test_params%get_nonlinear_convergence_criteria(), & 
         &                                         abs_tol = this%test_params%get_absolute_nonlinear_tolerance(),   &
         &                                         rel_tol = this%test_params%get_relative_nonlinear_tolerance(),   &
         &                                       max_iters = this%test_params%get_max_nonlinear_iterations()    ,   &
         &                                   linear_solver = this%linear_solver                                 ,   &
         &                                     environment = this%environment                                   ,   &
         &                                     fe_operator = this%fe_operator                                       )   

  end subroutine setup_operators
  
   subroutine setup_solver_after_pb_partition(this)
    implicit none
    class(par_test_hts_fe_driver_t), target, intent(inout) :: this
    type(parameterlist_t)                  :: parameter_list
    class(l1_coarse_fe_handler_t), pointer :: coarse_fe_handler
    type(parameterlist_t), pointer :: plist, dirichlet, neumann, coarse
    type(parameterlist_t) :: linear_pl
    integer(ip) :: ilev
    integer(ip) :: FPLError
    integer(ip) :: field_id
    integer(ip) :: istat
    integer(ip) :: iparm(64)
    class(matrix_t), pointer :: matrix
    
    call this%timer_find_subsets%start()     
    call this%set_pb_cells_set_id()
    call this%timer_find_subsets%stop() 
    call this%timer_find_subsets%report(.true.)
    
    call this%triangulation%setup_coarse_triangulation()
    call this%fe_space%set_up_cell_integration()
   
    call this%linear_solver%free()
    call this%linear_solver%create(this%fe_space%get_environment())
    call this%linear_solver%set_type_from_string(cg_name)

    ! BDDC preconditioner 
    call this%timer_setup_coarse_fe_handler%start() 
    call this%compute_average_parameter_values() 
    matrix => this%fe_operator%get_matrix() 
    select type ( matrix ) 
       class is (par_sparse_matrix_t) 
       call this%coarse_fe_handler%free()
       call this%coarse_fe_handler%create( 1, this%fe_space, matrix, this%parameter_list, this%average_permeability, this%average_resistivity )
       class DEFAULT
       assert(.false.) 
    end select
    call this%timer_setup_coarse_fe_handler%stop() 
    call this%timer_setup_coarse_fe_handler%report(.true.) 
    
    ! Prepare the internal parameter list of pardiso
    ! See https://software.intel.com/en-us/node/470298 for details
    iparm      = 0 ! Init all entries to zero
    iparm(1)   = 1  ! no solver default
    iparm(2)   = 2  ! fill-in reordering from METIS
    iparm(8)   = 2  ! numbers of iterative refinement steps
    iparm(10)  = 8  ! perturb the pivot elements with 1E-8
    iparm(21)  = 1  ! 1x1 + 2x2 pivots
    ! Customization
    iparm(11)  = 1 ! use scaling (default 0)
    iparm(13)  = 1 ! use maximum weighted matching algorithm (default 0)

    plist => this%parameter_list  
    if ( this%environment%get_l1_size() == 1 ) then
       FPLError = plist%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = plist%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
    end if
    do ilev=1, this%environment%get_num_levels()-1
       ! Dirichlet local problems 
       dirichlet => plist%NewSubList(key=mlbddc_dirichlet_solver_params)
       FPLError = dirichlet%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_spd); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
       FPLError = dirichlet%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       ! Constrained Neumann problems 
       neumann => plist%NewSubList(key=mlbddc_neumann_solver_params)
       FPLError = neumann%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
       FPLError = neumann%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
       FPLError = neumann%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)
       coarse => plist%NewSubList(key=mlbddc_coarse_solver_params) 
    end do
    ! Set coarsest-grid solver parameters
    FPLError = coarse%set(key=direct_solver_type, value=pardiso_mkl); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_matrix_type, value=pardiso_mkl_sin); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_message_level, value=0); assert(FPLError == 0)
    FPLError = coarse%set(key=pardiso_mkl_iparm, value=iparm); assert(FPLError == 0)

    call this%fe_space%setup_coarse_fe_space(this%parameter_list)
    
    call this%mlbddc%free() 
    call this%mlbddc%create(this%fe_operator, this%parameter_list)     

    ! Linear solver
    call this%linear_solver%create(this%fe_space%get_environment())
    call this%linear_solver%set_type_from_string(cg_name)
    call linear_pl%init()
    FPLError = linear_pl%set(key = ils_rtol, value = this%test_params%get_relative_linear_tolerance()); assert(FPLError == 0)
    FPLError = linear_pl%set(key = ils_max_num_iterations, value = 1000); assert(FPLError == 0)
    call this%linear_solver%set_parameters_from_pl(linear_pl)
    call this%linear_solver%set_operators( this%fe_operator%get_tangent(), this%mlbddc)
    call linear_pl%free()

    ! Nonlinear solver
    call this%nonlinear_solver%free()
    call this%nonlinear_solver%create(convergence_criteria = this%test_params%get_nonlinear_convergence_criteria(), & 
         &                                         abs_tol = this%test_params%get_absolute_nonlinear_tolerance(),   &
         &                                         rel_tol = this%test_params%get_relative_nonlinear_tolerance(),   &
         &                                       max_iters = this%test_params%get_max_nonlinear_iterations()    ,   &
         &                                   linear_solver = this%linear_solver                                 ,   &
         &                                     environment = this%environment                                   ,   &
         &                                     fe_operator = this%fe_operator                                       )   
    
  end subroutine  setup_solver_after_pb_partition

  ! -----------------------------------------------------------------------------------------------
  subroutine compute_average_parameter_values(this)
    implicit none 
    class(par_test_hts_fe_driver_t), intent(inout) :: this  
    class(fe_cell_iterator_t), allocatable :: fe
    ! Integration loop 
    type(quadrature_t)       , pointer     :: quad
    integer(ip)                            :: qpoin, num_quad_points
    type(point_t)            , pointer     :: quad_coords(:)
    real(rp)                               :: factor 
    real(rp), allocatable                  :: set_id_volume(:)
    real(rp)                               :: permeability
    type(tensor_field_t)                   :: resistivity
    integer(ip) :: max_cell_set_id, set_id, material_id 
    integer(ip) :: istat 
    ! Compute resistivity 
    type(fe_cell_function_vector_t)        :: fe_cell_function
    type(vector_field_t)                   :: H_value, H_curl 
    logical :: reallocate 

    if ( .not. this%environment%am_i_l1_task() ) return 

    max_cell_set_id = maxval(this%cells_set_ids)+1 ! 1-based arrays  
      
    if ( allocated(this%average_permeability) ) reallocate = ( size(this%average_permeability) < max_cell_set_id ) 
    if ( .not. allocated(this%average_permeability) .or. reallocate ) then 
      if ( allocated(this%average_permeability) ) then 
          call memfree( this%average_permeability, __FILE__, __LINE__ ) 
      end if 
      call memalloc( max_cell_set_id, this%average_permeability, __FILE__, __LINE__ )
    end if 
    
    if ( allocated(this%average_resistivity) ) reallocate = ( size(this%average_resistivity) < max_cell_set_id ) 
    if ( .not. allocated(this%average_resistivity) .or. reallocate ) then 
      if ( allocated(this%average_resistivity) ) then 
          call memfree( this%average_resistivity, __FILE__, __LINE__ ) 
      end if 
      call memalloc( max_cell_set_id, this%average_resistivity, __FILE__, __LINE__ )
    end if 
    
    call memalloc( max_cell_set_id, set_id_volume, __FILE__, __LINE__ )
    
    this%average_permeability = 0.0_rp 
    this%average_resistivity  = 0.0_rp 
    set_id_volume = 0.0_rp 

    ! Integrate structures needed 
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
    call fe%update_integration()
    call fe_cell_function%create(this%fe_space, 1)
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    quad_coords      => fe%get_quadrature_points_coordinates()

    ! Loop over elements
    do while ( .not. fe%has_finished())
       if ( fe%is_local() ) then  
          set_id = 1 + fe%get_set_id() ! 1-based arrays 
          call fe%update_integration()
          call fe_cell_function%update(fe, this%H_current) 
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Integrate cell contribution to averages 
          do qpoin=1, num_quad_points
            ! Evaluate parameters on the cell 
             call fe_cell_function%get_value(qpoin, H_value)
             call fe_cell_function%compute_curl(qpoin, H_curl)
          
             factor = fe%get_det_jacobian(qpoin) * quad%get_weight(qpoin) 		
             resistivity                       = this%hts_integration%compute_resistivity( H_value, H_curl, fe%get_set_id() ) 
             this%average_permeability(set_id) = this%average_permeability(set_id) + this%test_params%get_air_permeability()*factor
             this%average_resistivity(set_id)  = this%average_resistivity(set_id)  + resistivity%get(1,1)*factor
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
    call memfree(set_id_volume, __FILE__, __LINE__)

  end subroutine compute_average_parameter_values

  subroutine solve_problem (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    class(vector_t) , allocatable  :: rhs_dof_values
    class(vector_t) , pointer      :: free_dof_values 
    integer(ip) :: istat 

    ! Solve the nonlinear, transient problem
    call this%fe_operator%create_range_vector(rhs_dof_values) 
    call rhs_dof_values%init(0.0_rp)

    temporal: do while ( .not. this%theta_method%finished() ) 
       call this%theta_method%print(6) 

       ! Solve the problem
       free_dof_values => this%H_current%get_free_dof_values() 
       call this%nonlinear_solver%apply(rhs_dof_values, free_dof_values )

       if ( this%nonlinear_solver%has_converged() ) then  ! Theta method goes forward 
          ! call this%compute_hysteresis_data() 
          call this%theta_method%update_solutions(this%H_current, this%H_previous)
          call this%print_info()
          call this%write_time_step_solution()
          if ( this%test_params%get_is_adaptive_time_stepping() ) then 
          call this%theta_method%move_time_forward(this%nonlinear_solver%get_current_iteration(), this%test_params%get_stepping_parameter())
          else 
          call this%theta_method%move_time_forward()
          end if 

       elseif (.not. this%nonlinear_solver%has_converged()) then ! Theta method goes backwards and restarts   
          call this%theta_method%move_time_backwards(this%H_current, this%H_previous)
       end if

       if (.not. this%theta_method%finished() ) then 
          call this%fe_space%interpolate_dirichlet_values( this%H_current, time=this%theta_method%get_current_time() )
          call this%setup_solver_after_pb_partition()
       end if

    end do temporal

    call rhs_dof_values%free()
    deallocate(rhs_dof_values, stat=istat); check(istat==0);

  end subroutine solve_problem

  ! -----------------------------------------------------------------------------------------------
  subroutine compute_hysteresis_data(this)
    implicit none 
    class(par_test_hts_fe_driver_t), intent(inout) :: this  
    class(fe_cell_iterator_t), allocatable :: fe
    ! Integration loop 
    type(quadrature_t)       , pointer     :: quad
    type(fe_cell_function_vector_t)        :: fe_cell_function_current
    integer(ip)                            :: qpoin, num_quad_points
    integer(ip)                            :: idime, idof 
    type(point_t)            , pointer     :: quad_coords(:)
    type(point_t)         , allocatable    :: aux_quad_coords(:)
    integer(ip)                            :: inode, num_nodes 
    real(rp)                               :: factor 
    type(vector_field_t)                   :: H_value, H_curl 
    ! Hysteresis variables for final computations 
    real(rp)                               :: Hx_average, Hz_average, Halpha_average
    real(rp)                               :: hts_volume 
    type(vector_field_t)                   :: Happ
    class(scalar_function_t) , pointer     :: boundary_function_Hx
    class(scalar_function_t) , pointer     :: boundary_function_Hy
    class(scalar_function_t) , pointer     :: boundary_function_Hz
    type(vector_field_t)                   :: e_x, e_z, e_alpha
    real(rp)                               :: Happ_x, Happ_y, Happ_z 
    real(rp)                               :: Mx, Mz, Malpha 
    real(rp)                               :: H_app, Jc, w, h_b
    real(rp)                               :: domain(6)
    type(vector_field_t)                   :: r, r0 
    real(rp)                               :: Q_JE 
    type(tensor_field_t)                   :: resistivity

    integer(ip) :: istat 

    if ( this%test_params%get_is_analytical_solution() ) return
    if ( .not. this%environment%am_i_l1_task() ) return 

    ! Integrate structures needed 
    call fe_cell_function_current%create(this%fe_space,  1)
    call this%fe_space%set_up_cell_integration()
    call this%fe_space%create_fe_cell_iterator(fe)
    call fe%update_integration()
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_num_quadrature_points()
    quad_coords      => fe%get_quadrature_points_coordinates()
    aux_quad_coords  = quad_coords

    ! Define unitary vectors to take components 
    call e_x%init(0.0_rp); call e_x%set(1, 1.0_rp) 
    call e_z%init(0.0_rp); call e_z%set(3, 1.0_rp) 
    call e_alpha%init(0.0_rp); 
    if ( this%triangulation%get_num_dims()==2) then 
       call e_alpha%set(2, 1.0_rp) 
    else 
       call e_alpha%set(1, cos(pi/6) ) 
       call e_alpha%set(3, sin(pi/6) )
    end if

    ! Define r0 vector 
    if ( this%test_params%get_triangulation_type() == triangulation_generate_structured ) then
       domain = this%test_params%get_domain_limits()
       call r0%set(1, (domain(2)-domain(1))/2.0_rp )
       call r0%set(2, (domain(4)-domain(3))/2.0_rp )
       call r0%set(3, (domain(6)-domain(5))/2.0_rp )
    else 
       call r0%init(0.0_rp) 
    end if

    ! Loop over elements
    Hx_average      = 0.0_rp; Mx     = 0.0_rp 
    Hz_average      = 0.0_rp; Mz     = 0.0_rp 
    Halpha_average  = 0.0_rp; Malpha = 0.0_rp
    hts_volume      = 0.0_rp 
    Q_JE            = 0.0_rp 
    do while ( .not. fe%has_finished())

       if ( fe%get_set_id() > air .and. fe%is_local() ) then  ! Integrate only in HTS device DOMAIN 
          ! Update FE-integration related data structures
          call fe%update_integration()
          call fe_cell_function_current%update(fe, this%H_current)

          ! Get quadrature coordinates to evaluate boundary value
          quad_coords => fe%get_quadrature_points_coordinates()

          ! Integrate cell contribution to H_y, x·J_z average 
          do qpoin=1, num_quad_points
             factor = fe%get_det_jacobian(qpoin) * quad%get_weight(qpoin) 						         
             call fe_cell_function_current%get_value(qpoin, H_value)
             call fe_cell_function_current%compute_curl(qpoin, H_curl)
             ! Average magnetic field 
             Hx_average     = Hx_average     + factor*H_value*e_x         
             Hz_average     = Hz_average     + factor*H_value*e_z 
             Halpha_average = Halpha_average + factor*H_value*e_alpha 

             ! Magnetization (Pardo definition) 
             do idime=1,3 
                call r%set(idime, quad_coords(qpoin)%get(idime) - r0%get(idime) )
             end do
             Mx     = Mx     + factor*cross_product(r, H_curl)*e_x 
             Mz     = Mz     + factor*cross_product(r, H_curl)*e_z 
             Malpha = Malpha + factor*cross_product(r, H_curl)*e_alpha

             ! AC loss 
             resistivity = this%hts_integration%compute_resistivity( H_value, H_curl, HTS ) 
             Q_JE   = Q_JE + factor * resistivity * H_curl * H_curl 
          end do

          hts_volume = hts_volume + fe%compute_volume()
       end if
       call fe%next()
    end do
    call this%fe_space%free_fe_cell_iterator(fe)

    ! Coordinates of quadrature does influence the constant value Happ(t) 
    boundary_function_Hx => this%hts_analytical_functions%get_boundary_function_Hx()
    boundary_function_Hy => this%hts_analytical_functions%get_boundary_function_Hy()
    boundary_function_Hz => this%hts_analytical_functions%get_boundary_function_Hz()
    call boundary_function_Hx%get_value_space_time( aux_quad_coords(1), this%theta_method%get_current_time() , Happ_x )
    call boundary_function_Hy%get_value_space_time( aux_quad_coords(1), this%theta_method%get_current_time() , Happ_y )
    call boundary_function_Hz%get_value_space_time( aux_quad_coords(1), this%theta_method%get_current_time() , Happ_z )
    call Happ%set(1, Happ_x); call Happ%set(2, Happ_y); call Happ%set(3, Happ_z) 

    ! Communicate and obtain average values 
    call this%environment%l1_sum(Hx_average)
    call this%environment%l1_sum(Halpha_average)
    call this%environment%l1_sum(Hz_average)

    call this%environment%l1_sum(Mx)
    call this%environment%l1_sum(Mz)
    call this%environment%l1_sum(Malpha)

    call this%environment%l1_sum(hts_volume)
    call this%environment%l1_sum(Q_JE)

    if ( this%environment%get_l1_rank() == 0 ) then 
       Jc = this%test_params%get_critical_current() 
       w  = this%test_params%get_external_magnetic_field_frequency() 
       h_b = 0.01_rp  
       
       ! User-defined H_app 
       H_app = 1.0_rp/this%test_params%get_air_permeability()*200e-3_rp*sin(2.0_rp*pi*w*this%theta_method%get_current_time())

       write(*,*) ' Hysteresis Data -----------------------------------------'
       write(*,*) ' Happ', H_app/(Jc*h_b)
       write(*,*) ' Hx', this%test_params%get_hts_permeability()*(Hx_average/hts_volume-Happ*e_x)/(Jc*h_b)
       write(*,*) ' Hz', this%test_params%get_hts_permeability()*(Hz_average/hts_volume-Happ*e_z)/(Jc*h_b)
       write(*,*) ' Halpha', this%test_params%get_hts_permeability()*(Halpha_average/hts_volume-Happ*e_alpha)/(Jc*h_b)																										                     
       write(*,*) ' Mx', (Mx/2.0_rp/hts_volume)/(Jc*h_b)
       write(*,*) ' Mz', (Mz/2.0_rp/hts_volume)/(Jc*h_b) 
       write(*,*) ' Malpha', (Malpha/2.0_rp/hts_volume)/(Jc*h_b)			
       write(*,*) ' Q_JE', Q_JE
       write(*,*) ' HTS volume', hts_volume 
       write(*,*) ' --------------------------------------------------------' 
    end if

  end subroutine compute_hysteresis_data

  subroutine check_solution(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this
    type(error_norms_vector_t) :: error_norm 
    class(vector_function_t), pointer :: H_exact_function
    real(rp) :: l2, linfty, h1_s, hcurl 
    real(rp) :: error_tolerance, current_time 

    if ( this%test_params%get_is_analytical_solution() ) then 

       error_tolerance = 1.0e-3_rp 
       call error_norm%create(this%fe_space,1)   
       H_exact_function => this%maxwell_analytical_functions%get_solution_function()

       current_time = this%theta_method%get_current_time() - this%theta_method%get_time_step() 
       l2 = error_norm%compute(H_exact_function, this%H_current, l2_norm, time=current_time )     
       linfty = error_norm%compute(H_exact_function, this%H_current, linfty_norm, time=current_time)   
       h1_s = error_norm%compute(H_exact_function, this%H_current, h1_seminorm, time=current_time) 
       hcurl = error_norm%compute(H_exact_function, this%H_current, hcurl_seminorm, time=current_time)
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
    class(cell_iterator_t) , allocatable :: cell

    if ( this%environment%am_i_l1_task() ) then
       if(this%test_params%get_write_solution()) then
          call  build_cell_vectors()
          call  initialize_field_generator() 
          call  this%oh%create(VTK) 
          call  this%oh%attach_fe_space(this%fe_space)
          call  this%oh%add_fe_function(this%H_current, 1, 'H')
          call  this%oh%add_fe_function(this%H_current, 1, 'J', curl_diff_operator)
          call  this%oh%add_field_generator('Resistivity', this%resistivity_field_generator)
          call  this%oh%add_cell_vector(this%set_id_cell_vector, 'set_id')
          call  this%oh%add_cell_vector(this%rank_cell_vector, 'rank' ) 
          call  this%oh%open(this%test_params%get_dir_path(), this%test_params%get_prefix())
       endif
    end if

  contains  
    subroutine initialize_field_generator() 
      call this%resistivity_field_generator%set_parameter_values( this%test_params%get_nonlinear_exponent(),     & 
           this%test_params%get_critical_current(),       & 
           this%test_params%get_critical_electric_field() )                                                         
      call this%resistivity_field_generator%set_magnetic_field( this%H_current ) 
    end subroutine initialize_field_generator

    subroutine build_cell_vectors()
      call memalloc(this%triangulation%get_num_local_cells(), this%set_id_cell_vector, __FILE__, __LINE__)
      call memalloc(this%triangulation%get_num_local_cells(), this%rank_cell_vector, __FILE__, __LINE__)
      this%rank_cell_vector = this%environment%get_l1_rank() 
      call this%triangulation%create_cell_iterator(cell) 
      i=0
      do while ( .not. cell%has_finished() ) 
         if (cell%is_local() ) then  
            i=i+1
            this%set_id_cell_vector(i) = cell%get_set_id()
         end if
         call cell%next() 
      enddo
      call this%triangulation%free_cell_iterator(cell) 
    end subroutine build_cell_vectors
    
  end subroutine initialize_output

  ! -----------------------------------------------------------------------------------------------
  subroutine write_time_step_solution(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout)    :: this
    integer(ip)                                       :: err
    integer(ip)                                       :: i, istat
    class(cell_iterator_t) , allocatable :: cell

    if ( this%environment%am_i_l1_task() ) then
       if( this%test_params%get_write_solution() .and. this%theta_method%print_this_step() ) then
          call update_cell_vectors_values() 
          call this%oh%append_time_step(this%theta_method%get_current_time())
          call this%oh%write()
          call this%theta_method%update_time_to_be_printed() 
       endif
    end if
    
  contains 
      subroutine update_cell_vectors_values()
      call this%triangulation%create_cell_iterator(cell) 
      i=0
      do while ( .not. cell%has_finished() ) 
         if (cell%is_local() ) then  
            i=i+1
            this%set_id_cell_vector(i) = cell%get_set_id()
         end if
         call cell%next() 
      enddo
      call this%triangulation%free_cell_iterator(cell) 
    end subroutine update_cell_vectors_values
    
  end subroutine write_time_step_solution

  ! -----------------------------------------------------------------------------------------------
  subroutine finalize_output(this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout)    :: this
    integer(ip)                                       :: err

    if ( this%environment%am_i_l1_task() ) then
       if(this%test_params%get_write_solution()) then
          call this%oh%close()
          call this%oh%free()
          call memfree( this%set_id_cell_vector, __FILE__, __LINE__ )
          call memfree( this%rank_cell_vector, __FILE__, __LINE__ ) 
       endif
    end if
  end subroutine finalize_output

  ! ***********************************************************************************************
  subroutine run_simulation(this) 
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    call this%timer_run_simulation%start()

    call this%timer_triangulation%start()
    call this%setup_triangulation()
    call this%timer_triangulation%stop()

    call this%timer_fe_space%start()
    call this%setup_reference_fes()
    call this%setup_coarse_fe_handlers()
    call this%setup_fe_space()
    call this%timer_fe_space%stop()

    call this%timer_assemply%start()
    call this%setup_discrete_integration()
    call this%timer_assemply%stop()

    call this%timer_setup_system%start() 
    call this%setup_theta_method()
    call this%setup_system()
    call this%setup_operators() 
    call this%timer_setup_system%stop() 

    call this%initialize_output() 
    call this%print_info()

    call this%timer_solver_total%start()
    call this%solve_problem ()
    call this%timer_solver_total%stop()

    call this%finalize_output()
    call this%check_solution() 
    call this%free()
    call this%timer_run_simulation%stop()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(par_test_hts_fe_driver_t), target, intent(inout) :: this							
    integer(ip) :: i, istat

    call this%H_current%free() 
    call this%H_previous%free() 
    call this%linear_solver%free()
    call this%nonlinear_solver%free()
    call this%fe_operator%free()
    call this%mlbddc%free() 
    call this%linear_solver%free()
    if (allocated(this%average_resistivity))  call memfree( this%average_resistivity, __FILE__, __LINE__ ) 
    if (allocated(this%average_permeability)) call memfree( this%average_permeability, __FILE__, __LINE__ ) 

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
    if (allocated(this%cells_set_ids)) then 
       call memfree(this%cells_set_ids, __FILE__, __LINE__)
    end if
    call this%triangulation%free()
  end subroutine free

  !========================================================================================
  subroutine print_info (this)
    implicit none
    class(par_test_hts_fe_driver_t), intent(inout) :: this

    integer(ip) :: num_sub_domains
    real(rp) :: num_total_cells
    real(rp) :: num_dofs
    real(rp) :: num_fixed_dofs
    real(rp) :: num_interface_dofs 
    integer(ip) :: num_coarse_dofs
    real(rp) :: num_hts_cells 

    class(fe_cell_iterator_t), allocatable :: fe 
    class(environment_t), pointer :: environment
    class(coarse_fe_space_t), pointer :: coarse_fe_space

    environment => this%fe_space%get_environment()

    if (environment%am_i_l1_task()) then
       call this%fe_space%create_fe_cell_iterator(fe)  
       num_hts_cells=0
       do while ( .not. fe%has_finished())
          if (fe%is_local()) then
             if ( fe%get_set_id() > 0 ) num_hts_cells = num_hts_cells+1
          end if
          call fe%next()
       enddo
       call this%fe_space%free_fe_cell_iterator(fe)

       num_total_cells      = real(this%triangulation%get_num_local_cells(),kind=rp)
       num_dofs             = real(this%fe_space%get_field_num_dofs(1),kind=rp)
       num_fixed_dofs       = real(this%fe_space%get_num_fixed_dofs(),kind=rp)
       num_interface_dofs   = real(this%fe_space%get_total_num_interface_dofs(),kind=rp)

       call environment%l1_sum(num_total_cells )
       call environment%l1_sum(num_dofs        )
       call environment%l1_sum(num_fixed_dofs  )
       call environment%l1_sum(num_interface_dofs )
       call environment%l1_sum(num_hts_cells   ) 
    end if

    if (environment%get_l1_rank() == 0) then
       num_sub_domains = environment%get_l1_size()
       write(*,'(a,i22)') 'num_sub_domains:                    ', num_sub_domains
       write(*,'(a,i22)') 'num_cells_in_hts_domain:            ', nint(num_hts_cells       , kind=ip )
       write(*,'(a,i22)') 'num_total_cells:                    ', nint(num_total_cells     , kind=ip )
       write(*,'(a,i22)') 'num_dofs (sub-assembled):           ', nint(num_dofs            , kind=ip )
       write(*,'(a,i22)') 'num_fixed_dofs (sub-assembled):     ', nint(num_fixed_dofs      , kind=ip )
       write(*,'(a,i22)') 'num_interface_dofs (sub-assembled): ', nint(num_interface_dofs  , kind=ip )
    end if

    if (environment%am_i_lgt1_task()) then
       coarse_fe_space => this%fe_space%get_coarse_fe_space()
       num_coarse_dofs = coarse_fe_space%get_field_num_dofs(1)
       write(*,'(a,i22)') 'num_coarse_dofs:  ', num_coarse_dofs
    end if

  end subroutine print_info

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
