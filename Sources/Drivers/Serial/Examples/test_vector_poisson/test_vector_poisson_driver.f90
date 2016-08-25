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
module test_vector_poisson_driver_names
  use serial_names
!  use error_norms_names
  use test_poisson_params_names
  use vector_poisson_analytical_functions_names
  use vector_poisson_discrete_integration_names
  use vector_poisson_conditions_names
# include "debug.i90"

  implicit none
  private

  type test_vector_poisson_driver_t 
     private 

     ! Place-holder for parameter-value set provided through command-line interface
     type(test_poisson_params_t)                 :: test_params

     ! Cells and lower dimension objects container
     type(serial_triangulation_t)                :: triangulation

     ! Analytical functions of the problem
     type(vector_poisson_analytical_functions_t) :: problem_functions

     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                     :: fe_space 
     type(p_reference_fe_t) , allocatable        :: reference_fes(:) 
     type(vector_poisson_discrete_integration_t) :: vector_poisson_integration
     type(vector_poisson_conditions_t)           :: vector_poisson_conditions

     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                  :: fe_affine_operator

     ! Iterative linear solvers data type
     type(iterative_linear_solver_t)             :: iterative_linear_solver

     ! Poisson problem solution FE function
     type(fe_function_t)                         :: solution

   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     procedure        , private :: solve_system
     procedure        , private :: check_solution
     procedure        , private :: print_error_norms
     procedure        , private :: free
  end type test_vector_poisson_driver_t

  ! Types
  public :: test_vector_poisson_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_vector_poisson_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse()
  end subroutine parse_command_line_parameters

  subroutine setup_triangulation(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    call this%triangulation%create(this%test_params%get_dir_path(),&
         this%test_params%get_prefix())
  end subroutine setup_triangulation

  subroutine setup_reference_fes(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    integer(ip) :: istat

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)

    this%reference_fes(1) =  make_reference_fe ( topology = topology_tet,                                     &
                                                 fe_type = fe_type_lagrangian,                                &
                                                 number_dimensions = this%triangulation%get_num_dimensions(), &
                                                 order = 1,                                                   &
                                                 field_type = field_type_vector,                              &
                                                 continuity = .true. )
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this

    call this%fe_space%create( triangulation       = this%triangulation,      &
                               conditions          = this%vector_poisson_conditions, &
                               reference_fes       = this%reference_fes)
    call this%fe_space%fill_dof_info() 
    call this%vector_poisson_conditions%set_boundary_function(this%problem_functions%get_boundary_values())
    call this%fe_space%update_strong_dirichlet_bcs_values(this%vector_poisson_conditions)
  end subroutine setup_fe_space

  subroutine setup_system (this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    call this%vector_poisson_integration%set_source_term(this%problem_functions%get_source_term())
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = [ .true. ], &
                                          diagonal_blocks_symmetric         = [ .true. ], &
                                          diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
                                          fe_space                          = this%fe_space,           &
                                          discrete_integration              = this%vector_poisson_integration )
  end subroutine setup_system

  subroutine setup_solver (this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    integer               :: FPLError
    type(parameterlist_t) :: parameter_list
    integer               :: iparm(64)
    call this%iterative_linear_solver%create(this%fe_space%get_environment())
    call this%iterative_linear_solver%set_type_from_string(cg_name)
    call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
  end subroutine setup_solver

  subroutine assemble_system (this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    class(matrix_t), pointer       :: matrix
    class(vector_t), pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    rhs    => this%fe_affine_operator%get_translation()
    matrix => this%fe_affine_operator%get_matrix()
  end subroutine assemble_system

  subroutine solve_system(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    class(matrix_t), pointer       :: matrix
    class(vector_t), pointer       :: rhs
    class(vector_t), pointer       :: dof_values
    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
         dof_values)

    select type(dof_values)
       class is (serial_scalar_array_t)  
       call dof_values%print(6) 
       class DEFAULT
       assert(.false.) 
    end select

  end subroutine solve_system
  
    subroutine check_solution(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    class(vector_t), allocatable :: exact_solution_vector
    class(vector_t), pointer     :: computed_solution_vector

    call this%fe_affine_operator%create_range_vector(exact_solution_vector)
    call exact_solution_vector%init(1.0_rp)
    
    computed_solution_vector => this%solution%get_dof_values() 
    
    exact_solution_vector = computed_solution_vector - exact_solution_vector
   
    check ( exact_solution_vector%nrm2()/computed_solution_vector%nrm2() < 1.0e-04 )
     
    call exact_solution_vector%free()
    deallocate(exact_solution_vector)
   
  end subroutine check_solution
  
    subroutine print_error_norms(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    type(constant_vector_function_t) :: constant_function
    type(error_norms_vector_t) :: error_norm
    call constant_function%create(vector_field_t([1.0_rp,1.0_rp,0.0_rp]))
    call error_norm%create(this%fe_space,1)
    write(*,'(a20,e32.25)') 'mean_norm:', error_norm%compute(constant_function, this%solution, mean_norm)   
    write(*,'(a20,e32.25)') 'l1_norm:', error_norm%compute(constant_function, this%solution, l1_norm)   
    write(*,'(a20,e32.25)') 'l2_norm:', error_norm%compute(constant_function, this%solution, l2_norm)   
    write(*,'(a20,e32.25)') 'lp_norm:', error_norm%compute(constant_function, this%solution, lp_norm)   
    write(*,'(a20,e32.25)') 'linfnty_norm:', error_norm%compute(constant_function, this%solution, linfty_norm)   
    write(*,'(a20,e32.25)') 'h1_seminorm:', error_norm%compute(constant_function, this%solution, h1_seminorm)   
    write(*,'(a20,e32.25)') 'h1_norm:', error_norm%compute(constant_function, this%solution, h1_norm)   
    write(*,'(a20,e32.25)') 'hdiv_seminorm:', error_norm%compute(constant_function, this%solution, hdiv_seminorm)   
    write(*,'(a20,e32.25)') 'w1p_seminorm:', error_norm%compute(constant_function, this%solution, w1p_seminorm)   
    write(*,'(a20,e32.25)') 'w1p_norm:', error_norm%compute(constant_function, this%solution, w1p_norm)   
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', error_norm%compute(constant_function, this%solution, w1infty_seminorm)   
    write(*,'(a20,e32.25)') 'w1infty_norm:', error_norm%compute(constant_function, this%solution, w1infty_norm)   
    call error_norm%free()
  end subroutine print_error_norms 
  

  subroutine run_simulation(this) 
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    call this%setup_solver()
    call this%fe_space%create_fe_function(this%solution)
    call this%solve_system()
    call this%check_solution()
    call this%print_error_norms()
    call this%free()
  end subroutine run_simulation

  subroutine free(this)
    implicit none
    class(test_vector_poisson_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    call this%solution%free()
    call this%iterative_linear_solver%free()
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

end module test_vector_poisson_driver_names
