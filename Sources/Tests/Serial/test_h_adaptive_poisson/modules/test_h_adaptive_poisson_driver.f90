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
module test_h_adaptive_poisson_driver_names
  use fempar_names
  use test_poisson_params_names
  use poisson_cG_discrete_integration_names
  use poisson_conditions_names
  use poisson_analytical_functions_names
    
# include "debug.i90"

  implicit none
  private
  
  type test_h_adaptive_poisson_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(test_poisson_params_t)   :: test_params
     type(ParameterList_t)         :: parameter_list
     
     ! Cells and lower dimension objects container
     type(p4est_serial_triangulation_t)           :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(serial_hp_adaptive_fe_space_t)          :: fe_space 
     type(p_reference_fe_t), allocatable          :: reference_fes(:) 
     type(poisson_cG_discrete_integration_t)      :: poisson_cG_integration
     type(poisson_conditions_t)                   :: poisson_conditions
     type(poisson_analytical_functions_t)         :: poisson_analytical_functions
          
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                   :: fe_affine_operator
     
     ! Direct and Iterative linear solvers data type
#ifdef ENABLE_MKL     
     type(direct_solver_t)                     :: direct_solver
#else     
     type(iterative_linear_solver_t)           :: iterative_linear_solver
#endif     
 
     ! Poisson problem solution FE function
     type(fe_function_t)                       :: solution
   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: set_cells_for_refinement
     procedure        , private :: set_cells_for_coarsening
     procedure        , private :: fill_cells_set
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: write_solution
     procedure        , private :: free
  end type test_h_adaptive_poisson_driver_t

  ! Types
  public :: test_h_adaptive_poisson_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse(this%parameter_list)
  end subroutine parse_command_line_parameters
  
  subroutine setup_triangulation(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this

    class(cell_iterator_t), allocatable :: cell
    type(point_t), allocatable :: cell_coords(:)
    integer(ip) :: istat
    integer(ip) :: set_id
    real(rp) :: x, y
    integer(ip) :: num_void_neigs

    integer(ip)           :: ivef
    class(vef_iterator_t),allocatable :: vef, vef_of_vef
    type(list_t), pointer :: vefs_of_vef
    type(list_t), pointer :: vertices_of_line
    type(list_iterator_t) :: vefs_of_vef_iterator
    type(list_iterator_t) :: vertices_of_line_iterator
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    integer(ip) :: ivef_pos_in_cell, vef_of_vef_pos_in_cell
    integer(ip) :: vertex_pos_in_cell, icell_arround
    integer(ip) :: inode, num, i


    !call this%triangulation%create(this%test_params%get_dir_path(),&
    !                               this%test_params%get_prefix(),&
    !                               geometry_interpolation_order=this%test_params%get_reference_fe_geo_order())
    call this%triangulation%create(this%parameter_list)
   
    do i=1, 10
       call this%triangulation%clear_refinement_and_coarsening_flags()
       call this%set_cells_for_refinement()
       call this%fill_cells_set()
       call this%triangulation%refine_and_coarsen()

        
       if ( mod(i,3) == 0 ) then
         call this%triangulation%clear_refinement_and_coarsening_flags()
         call this%set_cells_for_coarsening()
         call this%triangulation%refine_and_coarsen()
       end if
    end do   
       
    !if ( this%test_params%get_triangulation_type() == 'structured' ) then
    !   call this%triangulation%create_vef_iterator(vef)
    !   do while ( .not. vef%has_finished() )
    !      if(vef%is_at_boundary()) then
    !         call vef%set_set_id(1)
    !      else
    !         call vef%set_set_id(0)
    !      end if
    !      call vef%next()
    !   end do
    !   call this%triangulation%free_vef_iterator(vef)
    !end if   
  end subroutine setup_triangulation
  
  subroutine set_cells_for_refinement(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    call this%triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      if ( mod(cell%get_lid()-1,2) == 0 ) then
        call cell%set_for_refinement()
      end if
      call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
  end subroutine set_cells_for_refinement
  
  subroutine set_cells_for_coarsening(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(cell_iterator_t)      , allocatable :: cell
    call this%triangulation%create_cell_iterator(cell)
    do while ( .not. cell%has_finished() )
      call cell%set_for_coarsening()
      call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
  end subroutine set_cells_for_coarsening
  
  subroutine fill_cells_set(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip), allocatable :: cell_set_ids(:)
    class(cell_iterator_t), allocatable :: cell
    
    call memalloc(this%triangulation%get_num_cells(),cell_set_ids)
    call this%triangulation%create_cell_iterator(cell)
    do while( .not. cell%has_finished() )
      if (cell%is_local()) then
         cell_set_ids(cell%get_lid()) = cell%get_lid()
      end if
      call cell%next()
    end do
    call this%triangulation%free_cell_iterator(cell)
    call this%triangulation%fill_cells_set(cell_set_ids)
    call memfree(cell_set_ids)
    
  end subroutine fill_cells_set
  
  subroutine setup_reference_fes(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    !
    !type(polytope_tree_t) :: poly, poly_old
    !integer(ip) :: topology
    
    ! Locals
    integer(ip) :: istat    
    class(cell_iterator_t)      , allocatable :: cell
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    
    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)
    
    call this%triangulation%create_cell_iterator(cell)
    reference_fe_geo => cell%get_reference_fe_geo()
    this%reference_fes(1) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                 fe_type = fe_type_lagrangian, &
                                                 number_dimensions = this%triangulation%get_num_dimensions(), &
                                                 order = this%test_params%get_reference_fe_order(), &
                                                 field_type = field_type_scalar, &
                                                 continuity = .true. )
    call this%triangulation%free_cell_iterator(cell)
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this    
    
    call this%poisson_analytical_functions%set_num_dimensions(this%triangulation%get_num_dimensions())
    call this%poisson_conditions%set_boundary_function(this%poisson_analytical_functions%get_boundary_function())
    call this%fe_space%create( triangulation       = this%triangulation, &
                               conditions          = this%poisson_conditions, &
                               reference_fes       = this%reference_fes)
        
    call this%fe_space%fill_dof_info() 
    call this%fe_space%initialize_fe_integration()    
    call this%fe_space%interpolate_dirichlet_values(this%poisson_conditions)
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this

    call this%poisson_cG_integration%set_analytical_functions(this%poisson_analytical_functions)
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                               diagonal_blocks_symmetric_storage = [ .true. ], &
                                               diagonal_blocks_symmetric         = [ .true. ], &
                                               diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                               fe_space                          = this%fe_space, &
                                               discrete_integration              = this%poisson_cG_integration )

  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer :: FPLError
    type(parameterlist_t) :: parameter_list
    integer :: iparm(64)
    class(matrix_t), pointer       :: matrix

    call parameter_list%init()
#ifdef ENABLE_MKL
    FPLError = parameter_list%set(key = direct_solver_type,        value = pardiso_mkl)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_spd)
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_message_level, value = 0)
    iparm = 0
    FPLError = FPLError + parameter_list%set(key = pardiso_mkl_iparm,         value = iparm)
    assert(FPLError == 0)
    
    call this%direct_solver%set_type_from_pl(parameter_list)
    call this%direct_solver%set_parameters_from_pl(parameter_list)
    
    matrix => this%fe_affine_operator%get_matrix()
    select type(matrix)
    class is (sparse_matrix_t)  
       call this%direct_solver%set_matrix(matrix)
    class DEFAULT
       assert(.false.) 
    end select
#else    
    FPLError = parameter_list%set(key = ils_rtol, value = 1.0e-12_rp)
    !FPLError = FPLError + parameter_list%set(key = ils_output_frequency, value = 30)
    FPLError = parameter_list%set(key = ils_max_num_iterations, value = 5000)
    assert(FPLError == 0)
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    call this%iterative_linear_solver%set_parameters_from_pl(parameter_list)
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
#endif
    call parameter_list%free()
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()
    
    select type(matrix)
    class is (sparse_matrix_t)  
       !call matrix%print_matrix_market(6) 
    class DEFAULT
       assert(.false.) 
    end select
    
    select type(rhs)
    class is (serial_scalar_array_t)  
    !call rhs%print(6) 
    class DEFAULT
       assert(.false.) 
    end select
  end subroutine assemble_system
  
  
  subroutine solve_system(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    
#ifdef ENABLE_MKL    
    call this%direct_solver%solve(this%fe_affine_operator%get_translation(), dof_values)
#else
    call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
                                            dof_values)
#endif    
    
    call this%solution%update_fixed_dof_values(this%fe_space)
    
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
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
    real(rp) :: error_tolerance
    
    call error_norm%create(this%fe_space,1)
    mean = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, mean_norm)   
    l1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l1_norm)   
    l2 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, l2_norm)   
    lp = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, lp_norm)   
    linfty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, linfty_norm)   
    h1_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_seminorm) 
    h1 = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, h1_norm) 
    w1p_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_seminorm)   
    w1p = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1p_norm)   
    w1infty_s = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_seminorm) 
    w1infty = error_norm%compute(this%poisson_analytical_functions%get_solution_function(), this%solution, w1infty_norm)

#ifdef ENABLE_MKL    
    error_tolerance = 1.0e-08
#else
    error_tolerance = 1.0e-06
#endif    
    
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < error_tolerance )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < error_tolerance )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < error_tolerance )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < error_tolerance )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < error_tolerance )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < error_tolerance )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < error_tolerance )
    call error_norm%free()
  end subroutine check_solution
  
  subroutine write_solution(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(in) :: this
    type(output_handler_t)                   :: oh
    character(len=:), allocatable            :: path
    character(len=:), allocatable            :: prefix
    if(this%test_params%get_write_solution()) then
        path = this%test_params%get_dir_path_out()
        prefix = this%test_params%get_prefix()
        call oh%create()
        call oh%attach_fe_space(this%fe_space)
        call oh%add_fe_function(this%solution, 1, 'solution')
        call oh%add_fe_function(this%solution, 1, 'grad_solution', grad_diff_operator)
        call oh%open(path, prefix)
        call oh%write()
        call oh%close()
        call oh%free()
    endif
  end subroutine write_solution
  
  subroutine run_simulation(this) 
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this    
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    call this%setup_solver()
    call this%solution%create(this%fe_space) 
    call this%solve_system()
    call this%check_solution()
    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(test_h_adaptive_poisson_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
    
#ifdef ENABLE_MKL        
    call this%direct_solver%free()
#else
    call this%iterative_linear_solver%free()
#endif
    
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
  end subroutine free  
  

  
end module test_h_adaptive_poisson_driver_names
