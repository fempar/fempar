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
module test_new_serial_fe_space_driver_names
  use serial_names
  use test_new_serial_fe_space_params_names
  use poisson_discrete_integration_names
# include "debug.i90"

  implicit none
  private

  type test_new_serial_fe_space_driver_t 
     private 
     
     ! Place-holder for parameter-value set provided through command-line interface
     type(test_new_serial_fe_space_params_t) :: test_params
     
     ! Cells and lower dimension objects container
     type(serial_triangulation_t)            :: triangulation
     
     ! Discrete weak problem integration-related data type instances 
     type(new_serial_fe_space_t)               :: fe_space 
     type(p_reference_fe_t), allocatable       :: reference_fes(:) 
     type(poisson_discrete_integration_t)      :: poisson_integration
     
     ! Place-holder for the coefficient matrix and RHS of the linear system
     type(new_fe_affine_operator_t)            :: fe_affine_operator
     
     ! Environment required for fe_affine_operator + vtk_handler
     type(serial_environment_t)                :: serial_environment
   contains
     procedure                  :: run_simulation
     procedure        , private :: parse_command_line_parameters
     procedure        , private :: setup_triangulation
     procedure        , private :: setup_reference_fes
     procedure        , private :: setup_fe_space
     procedure        , private :: setup_system
     !procedure        , private :: setup_solver
     procedure        , private :: assemble_system
     !procedure        , private :: solve_system
     procedure        , private :: free
  end type test_new_serial_fe_space_driver_t

  ! Types
  public :: test_new_serial_fe_space_driver_t

contains

  subroutine parse_command_line_parameters(this)
    implicit none
    class(test_new_serial_fe_space_driver_t ), intent(inout) :: this
    call this%test_params%create()
    call this%test_params%parse()
  end subroutine parse_command_line_parameters
  
  subroutine setup_triangulation(this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    call this%triangulation%create(this%test_params%get_dir_path(),&
                                   this%test_params%get_prefix())
    !call this%triangulation%print()
  end subroutine setup_triangulation
  
  subroutine setup_reference_fes(this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    
    integer(ip) :: istat
    
    ! if (test_single_scalar_valued_reference_fe) then
    allocate(this%reference_fes(1), stat=istat)
    check(istat==0)
    
    this%reference_fes(1) =  make_reference_fe ( topology = topology_quad, &
                                                 fe_type = fe_type_lagrangian, &
                                                 number_dimensions = this%triangulation%get_num_dimensions(), &
                                                 order = 1, &
                                                 field_type = field_type_scalar, &
                                                 continuity = .true. )
  end subroutine setup_reference_fes

  subroutine setup_fe_space(this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    call this%fe_space%create( triangulation       = this%triangulation, &
                               !boundary_conditions = this%conditions,    &
                               reference_fes       = this%reference_fes)
    call this%fe_space%fill_dof_info() 
    call this%fe_space%print()
  end subroutine setup_fe_space
  
  subroutine setup_system (this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    
    ! if (test_single_scalar_valued_reference_fe) then
    call this%fe_affine_operator%create ( sparse_matrix_storage_format      = csr_format, &
                                          diagonal_blocks_symmetric_storage = (/.false. /), &
                                          diagonal_blocks_symmetric         = (/.false. /), &
                                          diagonal_blocks_sign              = (/SPARSE_MATRIX_SIGN_INDEFINITE /), &
                                          environment                       = this%serial_environment, &
                                          fe_space                          = this%fe_space, &
                                          discrete_integration              = this%poisson_integration )
  end subroutine setup_system
  
  subroutine assemble_system (this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    class(matrix_t)                  , pointer       :: matrix
    class(vector_t)                  , pointer       :: rhs
    call this%fe_affine_operator%numerical_setup()
    rhs                => this%fe_affine_operator%get_translation()
    matrix             => this%fe_affine_operator%get_matrix()
    
    select type(matrix)
    class is (sparse_matrix_t)  
       call matrix%print_matrix_market(6) 
    class DEFAULT
       assert(.false.) 
    end select
    
    select type(rhs)
    class is (serial_scalar_array_t)  
       call rhs%print(6) 
    class DEFAULT
       assert(.false.) 
    end select
  end subroutine assemble_system
  
  
  subroutine run_simulation(this) 
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    call this%free()
    call this%parse_command_line_parameters()
    call this%setup_triangulation()
    call this%setup_reference_fes()
    call this%setup_fe_space()
    call this%setup_system()
    call this%assemble_system()
    call this%free()
  end subroutine run_simulation
  
  subroutine free(this)
    implicit none
    class(test_new_serial_fe_space_driver_t), intent(inout) :: this
    integer(ip) :: i, istat
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
  
end module test_new_serial_fe_space_driver_names
