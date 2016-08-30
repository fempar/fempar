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
module test_poisson_driver_names
  use serial_names
  use test_poisson_params_names
  use poisson_cG_discrete_integration_names
  use poisson_dG_discrete_integration_names
  use poisson_conditions_names
  use poisson_analytical_functions_names
  use list_types_names
  use vtk_handler_names
# include "debug.i90"

  implicit none
  private

  type test_poisson_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(test_poisson_params_t)   :: test_params
     
     ! Cells and lower dimension objects container
     type(serial_triangulation_t)              :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(serial_fe_space_t)                   :: fe_space 
     type(p_reference_fe_t), allocatable       :: reference_fes(:) 
     type(poisson_cG_discrete_integration_t)   :: poisson_cG_integration
     type(poisson_dG_discrete_integration_t)   :: poisson_dG_integration
     type(poisson_conditions_t)                :: poisson_conditions
     type(poisson_analytical_functions_t)      :: poisson_analytical_functions
     
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(fe_affine_operator_t)                :: fe_affine_operator
     
     ! Direct and Iterative linear solvers data type
     type(direct_solver_t)                     :: direct_solver
     type(iterative_linear_solver_t)           :: iterative_linear_solver
 
     ! Poisson problem solution FE function
     type(fe_function_t)                       :: solution
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
     procedure        , private :: write_solution
     procedure        , private :: free
  end type test_poisson_driver_t

  ! Types
  public :: test_poisson_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_poisson_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse()
  end subroutine parse_command_line_parameters
  
  subroutine setup_triangulation(this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    call this%triangulation%create(this%test_params%get_dir_path(),&
                                   this%test_params%get_prefix(),&
                                   geometry_interpolation_order=this%test_params%get_reference_fe_geo_order())
    !call this%triangulation%print()
  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    ! Locals
    integer(ip) :: istat    
    logical                                   :: continuity
    type(cell_iterator_t)                     :: cell_iterator
    type(cell_accessor_t)                     :: cell
    class(lagrangian_reference_fe_t), pointer :: reference_fe_geo
    

    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)
    
    continuity = .true.
    if ( trim(this%test_params%get_fe_formulation()) == 'dG' ) then
      continuity = .false.
    end if
    
    cell_iterator = this%triangulation%create_cell_iterator()
    call cell_iterator%current(cell)
    reference_fe_geo => cell%get_reference_fe_geo()
    
    this%reference_fes(1) =  make_reference_fe ( topology = reference_fe_geo%get_topology(), &
                                                 fe_type = fe_type_lagrangian, &
                                                 number_dimensions = this%triangulation%get_num_dimensions(), &
                                                 order = this%test_params%get_reference_fe_order(), & 
                                                 field_type = field_type_scalar, &
                                                 continuity = continuity ) 
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    
    call this%fe_space%create( triangulation       = this%triangulation, &
                               conditions          = this%poisson_conditions, &
                               reference_fes       = this%reference_fes)
    call this%fe_space%fill_dof_info() 
    call this%poisson_conditions%set_boundary_function(this%poisson_analytical_functions%get_boundary_function())
    call this%fe_space%update_strong_dirichlet_bcs_values(this%poisson_conditions)
    !call this%fe_space%print()
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    
    if ( trim(this%test_params%get_fe_formulation()) == 'cG' ) then
       call this%poisson_cG_integration%set_analytical_functions(this%poisson_analytical_functions)
       call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
            diagonal_blocks_symmetric_storage = [ .true. ], &
            diagonal_blocks_symmetric         = [ .true. ], &
            diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
            fe_space                          = this%fe_space, &
            discrete_integration              = this%poisson_cG_integration )

    else
       call this%poisson_dG_integration%set_analytical_functions(this%poisson_analytical_functions)
       call this%poisson_dG_integration%set_poisson_conditions(this%poisson_conditions)
       call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
            diagonal_blocks_symmetric_storage = [ .true. ], &
            diagonal_blocks_symmetric         = [ .true. ], &
            diagonal_blocks_sign              = [ SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE ], &
            fe_space                          = this%fe_space, &
            discrete_integration              = this%poisson_dG_integration )
    end if
  end subroutine setup_system
  
  subroutine setup_solver (this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    integer :: FPLError
    type(parameterlist_t) :: parameter_list
    integer :: iparm(64)
    class(matrix_t)                         , pointer       :: matrix

    call parameter_list%init()
    FPLError = FPLError + parameter_list%set(key = direct_solver_type,        value = pardiso_mkl)
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
    
    call parameter_list%free()
    
     !call this%iterative_linear_solver%create(this%fe_space%get_environment())
     !call this%iterative_linear_solver%set_type_from_string(cg_name)
     !call this%iterative_linear_solver%set_operators(this%fe_affine_operator, .identity. this%fe_affine_operator) 
  end subroutine setup_solver
  
  
  subroutine assemble_system (this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
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
    class(test_poisson_driver_t), intent(inout) :: this
    class(matrix_t)                         , pointer       :: matrix
    class(vector_t)                         , pointer       :: rhs
    class(vector_t)                         , pointer       :: dof_values

    matrix     => this%fe_affine_operator%get_matrix()
    rhs        => this%fe_affine_operator%get_translation()
    dof_values => this%solution%get_dof_values()
    !call this%iterative_linear_solver%solve(this%fe_affine_operator%get_translation(), &
    !                                        dof_values)
    call this%direct_solver%solve(this%fe_affine_operator%get_translation(), dof_values)
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
    class(test_poisson_driver_t), intent(inout) :: this
    type(error_norms_scalar_t) :: error_norm 
    real(rp) :: mean, l1, l2, lp, linfty, h1, h1_s, w1p_s, w1p, w1infty_s, w1infty
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
    
    write(*,'(a20,e32.25)') 'mean_norm:', mean; check ( abs(mean) < 1.0e-08 )
    write(*,'(a20,e32.25)') 'l1_norm:', l1; check ( l1 < 1.0e-08 )
    write(*,'(a20,e32.25)') 'l2_norm:', l2; check ( l2 < 1.0e-08 )
    write(*,'(a20,e32.25)') 'lp_norm:', lp; check ( lp < 1.0e-08 )
    write(*,'(a20,e32.25)') 'linfnty_norm:', linfty; check ( linfty < 1.0e-08 )
    write(*,'(a20,e32.25)') 'h1_seminorm:', h1_s; check ( h1_s < 1.0e-08 )
    write(*,'(a20,e32.25)') 'h1_norm:', h1; check ( h1 < 1.0e-08 )
    write(*,'(a20,e32.25)') 'w1p_seminorm:', w1p_s; check ( w1p_s < 1.0e-08 )
    write(*,'(a20,e32.25)') 'w1p_norm:', w1p; check ( w1p < 1.0e-08 )
    write(*,'(a20,e32.25)') 'w1infty_seminorm:', w1infty_s; check ( w1infty_s < 1.0e-08 )
    write(*,'(a20,e32.25)') 'w1infty_norm:', w1infty; check ( w1infty < 1.0e-08 )
    call error_norm%free()
  end subroutine check_solution
  
    subroutine write_solution(this)
    implicit none
    class(test_poisson_driver_t), intent(in) :: this
    type(vtk_handler_t)                      :: vtk_handler
    integer(ip)                              :: err
    print*, '----------------------------------------------------------------------'
    print*, '----------------------------------------------------------------------'
    print*, '----------------------------------------------------------------------'
    print*, '----------------------------------------------------------------------'
    print*, '----------------------------------------------------------------------'
    print*, this%test_params%get_write_solution()
    if(this%test_params%get_write_solution()) then
       call  vtk_handler%create(this%fe_space, this%test_params%get_dir_path(), this%test_params%get_prefix())
       err = vtk_handler%open_vtu(); check(err==0)
       err = vtk_handler%write_vtu_mesh(); check(err==0)
       err = vtk_handler%write_vtu_node_field(this%solution, 1, 'solution'); check(err==0)
       err = vtk_handler%close_vtu(); check(err==0)
       call  vtk_handler%free()
    endif
  end subroutine write_solution
  
  subroutine run_simulation(this) 
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this    
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
    call this%write_solution()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(test_poisson_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
    
    call this%solution%free()
    call this%iterative_linear_solver%free()
    call this%direct_solver%free()
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
  
end module test_poisson_driver_names
