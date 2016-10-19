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
module test_hts_nedelec_driver_names
  use fempar_names
  use hts_nedelec_params_names
  use hts_nedelec_analytical_functions_names
  use hts_nedelec_discrete_integration_names
  use hts_nedelec_conditions_names
  use hts_nonlinear_solver_names
  use hts_theta_method_names 
# include "debug.i90"

  implicit none
  private

  type test_hts_nedelec_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(hts_nedelec_params_t)               :: test_params
     type(ParameterList_t)                    :: parameter_list

     ! Cells and lower dimension objects container
     type(serial_triangulation_t)             :: triangulation

     ! Analytical functions of the problem
     type(hts_nedelec_analytical_functions_t) :: problem_functions 

     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                  :: fe_space 
     type(p_reference_fe_t) , allocatable     :: reference_fes(:) 
     type(hts_nedelec_discrete_integration_t) :: hts_nedelec_integration
     type(hts_nedelec_conditions_t)           :: hts_nedelec_conditions

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                  :: fe_affine_operator
     
     ! Temporal and Nonlinear solver data type 
     type(hts_nonlinear_solver_t)              :: nonlinear_solver 
     type(theta_method_t)                      :: theta_method 
  
     ! HTS problem solution FE function
     type(fe_function_t)            :: H_current 
     type(fe_function_t)            :: H_previous 
     
     ! Writing solution 
     type(vtk_handler_t)            :: vtk_handler

   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_nonlinear_solver
     procedure        , private :: setup_theta_method 
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: solve_nonlinear_system 
     procedure        , private :: check_solution
     procedure        , private :: initialize_output  
     procedure        , private :: write_time_step_solution
     procedure        , private :: finalize_output 
     procedure        , private :: free
  end type test_hts_nedelec_driver_t

  ! Types
  public :: test_hts_nedelec_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_hts_nedelec_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse(this%parameter_list)
  end subroutine parse_command_line_parameters

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_triangulation(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    ! Locals 
    integer(ip)                 , allocatable :: cells_set(:) 
    type(cell_iterator_t)                     :: cell_iterator
    type(cell_accessor_t)                     :: cell
    type(point_t), allocatable                :: cell_coordinates(:)
    integer(ip)                               :: inode
    integer(ip)       :: icell, icoord 
    real(rp)          :: cx, cy
    integer(ip)       :: istat 

    istat = 0
    call this%triangulation%create(this%parameter_list)

    ! Assign subset_id to different cells for the created structured mesh 
    allocate(cells_set(this%triangulation%get_num_cells() ), stat=istat); check(istat==0)
    cell_iterator = this%triangulation%create_cell_iterator()
    call cell_iterator%current(cell)
    allocate(cell_coordinates( cell%get_num_nodes() ) , stat=istat); check(istat==0) 
    
    do while ( .not. cell_iterator%has_finished() )
       call cell_iterator%current(cell)
       call cell%get_coordinates(cell_coordinates)
       ! Compute center of the element coordinates 
       cx = 0.0_rp
       cy = 0.0_rp 
       do inode=1,cell%get_num_nodes()  
          cx = cx + cell_coordinates(inode)%get(1)
          cy = cy + cell_coordinates(inode)%get(2)
       end do
       cx = cx/real(cell%get_num_nodes(),rp)
       cy = cy/real(cell%get_num_nodes(),rp)
       ! Select material case 
       !if ( ( (18e-3_rp<cx) .and. (cx<30e-3_rp) ) .and. ( (23.73e-3_rp<cy) .and. (cy<24.27e-3_rp) ) ) then 
       if ( ( (18e-3_rp<cx) .and. (cx<30e-3_rp) ) .and. ( (21e-3_rp<cy) .and. (cy<27e-3_rp) ) ) then 
          cells_set( cell%get_lid() ) = 1 ! HTS material 
       else 
          cells_set( cell%get_lid() ) = 2 ! Air material 
       end if
       call cell_iterator%next() 
    end do

    call this%triangulation%fill_cells_set(cells_set)  
    deallocate(cells_set, stat=istat); check(istat==0) 
    deallocate(cell_coordinates, stat=istat); check(istat==0) 
  end subroutine setup_triangulation

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_reference_fes(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip) :: istat
    type(vef_iterator_t)  :: vef_iterator
    type(vef_accessor_t)  :: vef

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)

    this%reference_fes(1) =  make_reference_fe ( topology = topology_hex, &
                                                 fe_type = fe_type_nedelec, &
                                                 number_dimensions = this%triangulation%get_num_dimensions(), &
                                                 order = this%test_params%get_reference_fe_order(), &
                                                 field_type = field_type_vector, &
                                                 continuity = .true. ) 
    
    if ( trim(this%test_params%get_triangulation_type()) == 'structured' ) then
       vef_iterator = this%triangulation%create_vef_iterator()
       do while ( .not. vef_iterator%has_finished() )
          call vef_iterator%current(vef)
          if(vef%is_at_boundary()) then
             call vef%set_set_id(1)
          else
             call vef%set_set_id(0)
          end if
          call vef_iterator%next()
       end do
    end if    
    
  end subroutine setup_reference_fes

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_fe_space(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this

    call this%hts_nedelec_conditions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%fe_space%create( triangulation       = this%triangulation,      &
         conditions          = this%hts_nedelec_conditions, &
         reference_fes       = this%reference_fes)
    call this%fe_space%fill_dof_info() 
    call this%fe_space%initialize_fe_integration()
    call this%fe_space%initialize_fe_face_integration() 
    call this%hts_nedelec_conditions%set_boundary_function_Hx(this%problem_functions%get_boundary_function_Hx())
    call this%hts_nedelec_conditions%set_boundary_function_Hy(this%problem_functions%get_boundary_function_Hy())
    if ( this%triangulation%get_num_dimensions() == 3) then 
       call this%hts_nedelec_conditions%set_boundary_function_Hz(this%problem_functions%get_boundary_function_Hz())
    end if
    ! Create H_previous with initial time (t0) boundary conditions 
    call this%fe_space%project_dirichlet_values_curl_conforming(this%hts_nedelec_conditions, this%theta_method%get_initial_time() )
    call this%H_previous%create(this%fe_space) 
    ! Update fe_space to the current time (t1) boundary conditions, create H_current 
    call this%fe_space%project_dirichlet_values_curl_conforming(this%hts_nedelec_conditions, this%theta_method%get_current_time() )
    call this%H_current%create(this%fe_space)

    !call this%fe_space%print()
  end subroutine setup_fe_space

  ! -----------------------------------------------------------------------------------------------
  subroutine setup_system (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this 
    ! Need to initialize dof_values, no interpolation available in Nedelec
    class(vector_t) , pointer :: dof_values_current 
    class(vector_t) , pointer :: dof_values_previous
    
    dof_values_current => this%H_current%get_dof_values() 
    dof_values_previous => this%H_previous%get_dof_values()
    call dof_values_current%init(0.0_rp) 
    call dof_values_previous%init(0.0_rp) 
    
    call this%problem_functions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%hts_nedelec_integration%create( this%theta_method, this%H_current, this%H_previous, &
                                             this%test_params, this%problem_functions%get_source_term() )
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .false.  ], &
                                          diagonal_blocks_symmetric         = [ .false. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_UNKNOWN ], &
                                          fe_space                          = this%fe_space,           &
                                          discrete_integration              = this%hts_nedelec_integration )
    
    nullify(dof_values_current) 
    nullify(dof_values_previous) 
  end subroutine setup_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_nonlinear_solver (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    
    call this%nonlinear_solver%create( abs_tol = this%test_params%get_absolute_nonlinear_tolerance(),  &
                                       rel_tol = this%test_params%get_relative_nonlinear_tolerance(),  &
                                       max_iters = this%test_params%get_max_nonlinear_iterations(),    &
                                       ideal_iters = this%test_params%get_stepping_parameter(),        &
                                       fe_affine_operator = this%fe_affine_operator,                   &
                                       current_dof_values = this%H_current%get_dof_values()            )
    
  end subroutine setup_nonlinear_solver
  
  ! -----------------------------------------------------------------------------------------------
  subroutine setup_theta_method(this) 
  implicit none 
  class(test_hts_nedelec_driver_t), intent(inout) :: this
  
    call this%theta_method%create( this%test_params%get_theta_value(),          &               
                                   this%test_params%get_initial_time(),         &
                                   this%test_params%get_final_time(),           & 
                                   this%test_params%get_number_time_steps(),    &
                                   this%test_params%get_max_time_step(),        & 
                                   this%test_params%get_min_time_step()         )
  
  
  end subroutine setup_theta_method 

  ! -----------------------------------------------------------------------------------------------
  subroutine assemble_system (this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(matrix_t), pointer       :: matrix
    class(vector_t), pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    !rhs    => this%fe_affine_operator%get_translation()
    !matrix => this%fe_affine_operator%get_matrix()
    !select type(matrix)
    !class is (sparse_matrix_t)  
    !   call matrix%print_matrix_market(6) 
    !class DEFAULT
    !   assert(.false.) 
    !end select
  end subroutine assemble_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine solve_system(this) 
  implicit none 
  class(test_hts_nedelec_driver_t), intent(inout) :: this
    
  temporal: do while ( .not. this%theta_method%finished() ) 
     call this%theta_method%print(6) 
     
     call this%solve_nonlinear_system()

     if (this%nonlinear_solver%converged() ) then  ! Theta method goes forward 
        call this%theta_method%update_solutions(this%H_current, this%H_previous)
        call this%write_time_step_solution() 
        call this%theta_method%move_time_forward( this%nonlinear_solver%get_current_iteration(), &
                                                  this%nonlinear_solver%get_ideal_num_iterations() ) 
     elseif (.not. this%nonlinear_solver%converged()) then ! Theta method goes backwards and restarts   
        call this%theta_method%move_time_backwards(this%H_current, this%H_previous)
     end if

     if (.not. this%theta_method%finished() ) then 
        call this%fe_space%project_dirichlet_values_curl_conforming(this%hts_nedelec_conditions, this%theta_method%get_current_time() )
        call this%H_current%update_strong_dirichlet_values(this%fe_space) 
        call this%assemble_system() 
     end if

  end do temporal
 
  end subroutine solve_system 
 
  ! -----------------------------------------------------------------------------------------------
   subroutine solve_nonlinear_system(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(matrix_t) , pointer :: A
  
    call this%nonlinear_solver%initialize() 
    nonlinear: do while ( .not. this%nonlinear_solver%finished() )
   
    ! 0 - Update counter
    call this%nonlinear_solver%start_new_iteration() 
    ! 1 - Evaluate initial residual 
    call this%nonlinear_solver%compute_residual(this%fe_affine_operator, this%H_current%get_dof_values() )
    ! 2 - Integrate Jacobian
    call this%nonlinear_solver%compute_jacobian(this%hts_nedelec_integration, this%fe_affine_operator)
    ! 3 - Solve tangent system 
    !A => this%fe_affine_operator%get_matrix() 
    !  select type(A)
    !   class is (sparse_matrix_t)  
    !   call A%print_matrix_market(6) 
    !   class DEFAULT
    !   assert(.false.) 
    !end select  
    call this%nonlinear_solver%solve_tangent_system(this%fe_affine_operator) 
    ! 4 - Update solution 
    call this%nonlinear_solver%update_solution(this%H_current) 
    ! 5 - New picard iterate with updated solution 
    call this%assemble_system() 
    ! 6 - Evaluate new residual 
    call this%nonlinear_solver%compute_residual(this%fe_affine_operator, this%H_current%get_dof_values() )
    ! 7 - Print current output 
    call this%nonlinear_solver%print_current_iteration_output() 
    
    end do nonlinear 
    
  end subroutine solve_nonlinear_system
  
  ! -----------------------------------------------------------------------------------------------
  subroutine check_solution(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    class(vector_function_t), pointer :: H_exact_function
    type(error_norms_vector_t) :: H_error_norm
    real(rp) :: l2, hcurl
    real(rp) :: error_tolerance
    
    H_exact_function => this%problem_functions%get_solution()
    
    call H_error_norm%create(this%fe_space,1)

    write(*,*) 'H ERROR NORMS'  
    l2 = H_error_norm%compute(H_exact_function, this%H_current, l2_norm, time=this%theta_method%get_current_time() - this%theta_method%get_time_step() )   
    hcurl = H_error_norm%compute(H_exact_function, this%H_current, hcurl_seminorm, time=this%theta_method%get_current_time() - this%theta_method%get_time_step() )    
    error_tolerance = 1.0e-04

    write(*,'(a20,f20.16)') 'l2_norm:', l2; !check ( l2 < error_tolerance )
    write(*,'(a20,f20.16)') 'hcurl_norm:', hcurl; !check ( h1 < error_tolerance )
    
    call H_error_norm%free()
  end subroutine check_solution 
  
    ! -----------------------------------------------------------------------------------------------
  subroutine initialize_output(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    if(this%test_params%get_write_solution()) then
       call  this%vtk_handler%create(this%fe_space, this%test_params%get_dir_path_out(), this%test_params%get_prefix())
    endif
  end subroutine initialize_output
  
  ! -----------------------------------------------------------------------------------------------
  subroutine write_time_step_solution(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    if(this%test_params%get_write_solution()) then
       err = this%vtk_handler%open_vtu(time_step=this%theta_method%get_current_time() ,format='ascii'); check(err==0)
       err = this%vtk_handler%write_vtu_mesh(this%H_current); check(err==0)
       err = this%vtk_handler%write_vtu_node_field(this%H_current, 1, 'H'); check(err==0)
       err = this%vtk_handler%write_vtu_node_field(this%H_current, 1, 'J'); check(err==0)
       err = this%vtk_handler%close_vtu(); check(err==0)
       err = this%vtk_handler%write_pvtu(time_step=this%theta_method%get_current_time() ); check(err==0)
       err = this%vtk_handler%write_pvd(); check(err==0)
    endif
  end subroutine write_time_step_solution
  
      ! -----------------------------------------------------------------------------------------------
  subroutine finalize_output(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip)                                      :: err
    if(this%test_params%get_write_solution()) then
    call  this%vtk_handler%free()
    endif
  end subroutine finalize_output

  ! -----------------------------------------------------------------------------------------------
  subroutine run_simulation(this) 
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_theta_method() 
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()   
    call this%setup_nonlinear_solver()
    call this%initialize_output() 
    call this%solve_system()
    call this%finalize_output()
    call this%check_solution() 
    call this%free()
  end subroutine run_simulation

  ! -----------------------------------------------------------------------------------------------
  subroutine free(this)
    implicit none
    class(test_hts_nedelec_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    call this%H_previous%free()
    call this%H_current%free() 
    call this%fe_affine_operator%free()
    call this%fe_space%free()
    if ( allocated(this%reference_fes) ) then
       do i=1, size(this%reference_fes)
          call this%reference_fes(i)%p%free()
       end do
       deallocate(this%reference_fes, stat=istat)
       check(istat==0)
    end if
    call this%triangulation%free()
    call this%test_params%free()
    call this%nonlinear_solver%free()
  end subroutine free

end module test_hts_nedelec_driver_names
