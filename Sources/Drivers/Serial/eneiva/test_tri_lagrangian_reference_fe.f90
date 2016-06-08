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

# include "command_line_parameters.i90"
# include "poisson_discrete_integration.i90"
# include "vector_laplacian_composite_discrete_integration.i90"

!****************************************************************************************************
program test_reference_fe
  use serial_names
  use command_line_parameters_names
  implicit none
#include "debug.i90"

  ! Our data
  type(test_reference_fe_params_t) :: params
  type(mesh_t)                     :: f_mesh
  type(conditions_t)               :: f_cond
  type(triangulation_t)            :: f_trian

  call fempar_init()  
  
  call params%create()
  call params%parse()

  ! Read mesh
  call mesh_read (params%dir_path, params%prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call conditions_read (params%dir_path, params%prefix, f_mesh%npoin, f_cond)

  ! Construct triangulation
  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  if ( trim(params%laplacian_type) == 'scalar' ) then
    call test_single_scalar_valued_reference_fe()
  else  
  !  call test_single_vector_valued_reference_fe()
  !  call test_composite_reference_fe_monolithic()
  !  call test_composite_reference_fe_block()
  end if  
 
  call triangulation_free(f_trian)
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)
  call params%free()

  call fempar_finalize()

contains  
  
subroutine test_single_scalar_valued_reference_fe ()
    use poisson_discrete_integration_names
    implicit none
    type(p_reference_fe_t)               :: reference_fe_array(1)
    type(serial_fe_space_t)              :: fe_space
    type(fe_affine_operator_t)           :: fe_affine_operator
    type(poisson_discrete_integration_t) :: poisson_integration
    type(iterative_linear_solver_t)      :: iterative_linear_solver
    type(serial_environment_t)           :: senv
    class(vector_t), allocatable         :: computed_solution_vector, exact_solution_vector
    class(matrix_t)            , pointer :: matrix
    class(array_t)             , pointer :: array    
    
    ! Simple case
     reference_fe_array(1) =  make_reference_fe ( topology = topology_tet,              &
                                                  fe_type = fe_type_lagrangian,         &
                                                  number_dimensions = f_trian%num_dims, &
                                                  order = params%order,                 &
                                                  field_type = field_type_scalar,       &
                                                  continuity = .true. )
     
     call fe_space%create( triangulation = f_trian,      &
                           boundary_conditions = f_cond, &
                           reference_fe_phy = reference_fe_array )
     
     call fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
                                    bc_code = 1,                                        &
                                    fe_space_component = 1 )
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format=csr_format, &
                                     diagonal_blocks_symmetric_storage=(/.true./), &
                                     diagonal_blocks_symmetric=(/.true./), &
                                     diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                     environment=senv, &
                                     fe_space=fe_space, &
                                     discrete_integration=poisson_integration )
     
     call fe_affine_operator%symbolic_setup()
     call fe_affine_operator%numerical_setup()
  
     !!matrix => fe_affine_operator%get_matrix()
     !!select type(matrix)
     !! class is (sparse_matrix_t)
     !!   call matrix%print_matrix_market(6)
     !!end select
 
     !call fe_affine_operator%create_range_vector(computed_solution_vector)
     !call fe_affine_operator%create_range_vector(exact_solution_vector)
     !call exact_solution_vector%init(1.0_rp)

     !! Create iterative linear solver, set operators and solve linear system
     !call iterative_linear_solver%create(senv)
     !call iterative_linear_solver%set_type_from_string(cg_name)
     !call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     !call iterative_linear_solver%solve(fe_affine_operator%get_translation(), computed_solution_vector)
     !call iterative_linear_solver%free() 

     !!select type(computed_solution_vector)
     !!  class is(serial_scalar_array_t)
     !!  call computed_solution_vector%print(6)
     !!  class is(serial_block_array_t)
     !!  call computed_solution_vector%print(6)
     !!  class default
     !!  check(.false.) 
     !!end select
  
     !computed_solution_vector = computed_solution_vector - exact_solution_vector
     !check ( computed_solution_vector%nrm2()/exact_solution_vector%nrm2() < 1.0e-04 )
     
     !call computed_solution_vector%free()
     !deallocate(computed_solution_vector)

     !call exact_solution_vector%free()
     !deallocate(exact_solution_vector)
     
     call fe_affine_operator%free()
     call fe_space%free()
     
     call reference_fe_array(1)%free()
     
  end subroutine test_single_scalar_valued_reference_fe
  
  subroutine test_single_vector_valued_reference_fe ()
    use vector_laplacian_single_discrete_integration_names
    implicit none
    type(p_reference_fe_t)                        :: reference_fe_array(1)
    type(serial_fe_space_t)                       :: fe_space
    type(fe_affine_operator_t)                    :: fe_affine_operator
    type(vector_laplacian_single_discrete_integration_t) :: vector_laplacian_integration
    type(iterative_linear_solver_t)               :: iterative_linear_solver
    type(serial_environment_t)                    :: senv
    class(vector_t), allocatable                  :: computed_solution_vector, exact_solution_vector
    class(matrix_t)            , pointer          :: matrix
    class(array_t)             , pointer          :: array

    ! Simple case
     reference_fe_array(1) =  make_reference_fe ( topology = topology_quad, &
                                                  fe_type = fe_type_lagrangian, &
                                                  number_dimensions = f_trian%num_dims, &
                                                  order = params%order, &
                                                  field_type = field_type_vector, &
                                                  continuity = .true. )
     
     call fe_space%create( triangulation = f_trian, &
                           boundary_conditions = f_cond, &
                           reference_fe_phy = reference_fe_array )
     
     call fe_space%update_bc_value (vector_function=constant_vector_function_t(vector_field_t(1.0_rp)), &
                                    bc_code = 1, &
                                    fe_space_component = 1 )
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format=csr_format, &
                                     diagonal_blocks_symmetric_storage=(/.true./), &
                                     diagonal_blocks_symmetric=(/.true./), &
                                     diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                     environment=senv, &
                                     fe_space=fe_space, &
                                     discrete_integration=vector_laplacian_integration )
     
     call fe_affine_operator%symbolic_setup()
     call fe_affine_operator%numerical_setup()
  
     !matrix => fe_affine_operator%get_matrix()
     !select type(matrix)
     ! class is (sparse_matrix_t)
     !   call matrix%print_matrix_market(6)
     !end select
 
     call fe_affine_operator%create_range_vector(computed_solution_vector)
     call fe_affine_operator%create_range_vector(exact_solution_vector)
     call exact_solution_vector%init(1.0_rp)

     ! Create iterative linear solver, set operators and solve linear system
     call iterative_linear_solver%create(senv)
     call iterative_linear_solver%set_type_from_string(cg_name)
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(fe_affine_operator%get_translation(), computed_solution_vector)
     call iterative_linear_solver%free() 

     !select type(computed_solution_vector)
     !  class is(serial_scalar_array_t)
     !  call computed_solution_vector%print(6)
     !  class is(serial_block_array_t)
     !  call computed_solution_vector%print(6)
     !  class default
     !  check(.false.) 
     !end select
  
     computed_solution_vector = computed_solution_vector - exact_solution_vector
     check ( computed_solution_vector%nrm2()/exact_solution_vector%nrm2() < 1.0e-04 )
     
     call computed_solution_vector%free()
     deallocate(computed_solution_vector)

     call exact_solution_vector%free()
     deallocate(exact_solution_vector)
     
     call fe_affine_operator%free()
     call fe_space%free()
     
     call reference_fe_array(1)%free()
  end subroutine test_single_vector_valued_reference_fe  

  subroutine test_composite_reference_fe_monolithic ()
    use vector_laplacian_composite_discrete_integration_names
    implicit none
    type(p_reference_fe_t)                        :: reference_fe_array(2)
    type(serial_fe_space_t)                       :: fe_space
    type(fe_affine_operator_t)                    :: fe_affine_operator
    type(vector_laplacian_composite_discrete_integration_t) :: vector_laplacian_integration
    type(iterative_linear_solver_t)               :: iterative_linear_solver
    type(serial_environment_t)                    :: senv
    class(vector_t), allocatable                  :: computed_solution_vector, exact_solution_vector
    class(matrix_t)            , pointer          :: matrix
    class(array_t)             , pointer          :: array

    ! Simple case
    reference_fe_array(1) = make_reference_fe ( topology = topology_quad, &
                                                    fe_type = fe_type_lagrangian, &
                                                    number_dimensions = f_trian%num_dims, &
                                                    order = params%order, &
                                                    field_type = field_type_scalar, &
                                                    continuity = .true. )
     
    reference_fe_array(2) = make_reference_fe ( topology = topology_quad, &
                                                fe_type = fe_type_lagrangian, &
                                                number_dimensions = f_trian%num_dims, &
                                                order = params%order, & 
                                                field_type = field_type_scalar, &
                                                continuity = .true. )
     
     call fe_space%create( triangulation = f_trian, &
                           boundary_conditions = f_cond, &
                           reference_fe_phy = reference_fe_array, &
                           field_blocks = (/1,1/), &
                           field_coupling = reshape((/.true.,.false.,.false.,.true./),(/2,2/)) )
     
     call fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
                                    bc_code = 1, &
                                    fe_space_component = 1 )
     
     call fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
                                    bc_code = 1, &
                                    fe_space_component = 2 )
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format=csr_format, &
                                     diagonal_blocks_symmetric_storage=(/.true./), &
                                     diagonal_blocks_symmetric=(/.true./), &
                                     diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                     environment=senv, &
                                     fe_space=fe_space, &
                                     discrete_integration=vector_laplacian_integration )
     
     call fe_affine_operator%symbolic_setup()
     call fe_affine_operator%numerical_setup()
  
     !matrix => fe_affine_operator%get_matrix()
     !select type(matrix)
     ! class is (sparse_matrix_t)
     !   call matrix%print_matrix_market(6)
     !end select
 
     call fe_affine_operator%create_range_vector(computed_solution_vector)
     call fe_affine_operator%create_range_vector(exact_solution_vector)
     call exact_solution_vector%init(1.0_rp)

     ! Create iterative linear solver, set operators and solve linear system
     call iterative_linear_solver%create(senv)
     call iterative_linear_solver%set_type_from_string(cg_name)
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(fe_affine_operator%get_translation(), computed_solution_vector)
     call iterative_linear_solver%free() 

     !select type(computed_solution_vector)
     !  class is(serial_scalar_array_t)
     !  call computed_solution_vector%print(6)
     !  class is(serial_block_array_t)
     !  call computed_solution_vector%print(6)
     !  class default
     !  check(.false.) 
     !end select
  
     computed_solution_vector = computed_solution_vector - exact_solution_vector
     check ( computed_solution_vector%nrm2()/exact_solution_vector%nrm2() < 1.0e-04 )
     
     call computed_solution_vector%free()
     deallocate(computed_solution_vector)

     call exact_solution_vector%free()
     deallocate(exact_solution_vector)
     
     call fe_affine_operator%free()
     call fe_space%free()
     
     call reference_fe_array(1)%free()
     
     call reference_fe_array(1)%free()
     call reference_fe_array(2)%free()
  end subroutine test_composite_reference_fe_monolithic    

  subroutine test_composite_reference_fe_block ()
    use vector_laplacian_composite_discrete_integration_names
    implicit none
    type(p_reference_fe_t)                        :: reference_fe_array(2)
    type(serial_fe_space_t)                       :: fe_space
    type(fe_affine_operator_t)                    :: fe_affine_operator
    type(vector_laplacian_composite_discrete_integration_t) :: vector_laplacian_integration
    type(iterative_linear_solver_t)               :: iterative_linear_solver
    type(serial_environment_t)                    :: senv
    class(vector_t), allocatable                  :: computed_solution_vector, exact_solution_vector
    class(matrix_t)            , pointer          :: matrix
    class(array_t)             , pointer          :: array

    ! Simple case
    reference_fe_array(1) = make_reference_fe ( topology = topology_quad, &
                                                    fe_type = fe_type_lagrangian, &
                                                    number_dimensions = f_trian%num_dims, &
                                                    order = params%order, &
                                                    field_type = field_type_scalar, &
                                                    continuity = .true. )
     
    reference_fe_array(2) = make_reference_fe ( topology = topology_quad, &
                                                fe_type = fe_type_lagrangian, &
                                                number_dimensions = f_trian%num_dims, &
                                                order = params%order, & 
                                                field_type = field_type_scalar, &
                                                continuity = .true. )
     
     call fe_space%create( triangulation = f_trian, &
                           boundary_conditions = f_cond, &
                           reference_fe_phy = reference_fe_array, &
                           field_blocks = (/1,2/), &
                           field_coupling = reshape((/.true.,.false.,.false.,.true./),(/2,2/)) )
     
     call fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
                                    bc_code = 1, &
                                    fe_space_component = 1 )
     
     call fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
                                    bc_code = 1, &
                                    fe_space_component = 2 )
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create ( csr_format, &
                                      (/.true.,.true./), &
                                      (/.true.,.true./), &
                                      (/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE,SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/),&
                                      senv, &
                                      fe_space, &
                                      vector_laplacian_integration )
     
     call fe_affine_operator%symbolic_setup()
     call fe_affine_operator%numerical_setup()
  
     !matrix => fe_affine_operator%get_matrix()
     !select type(matrix)
     ! class is (sparse_matrix_t)
     !   call matrix%print_matrix_market(6)
     !end select
 
     call fe_affine_operator%create_range_vector(computed_solution_vector)
     call fe_affine_operator%create_range_vector(exact_solution_vector)
     call exact_solution_vector%init(1.0_rp)

     ! Create iterative linear solver, set operators and solve linear system
     call iterative_linear_solver%create(senv)
     call iterative_linear_solver%set_type_from_string(cg_name)
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(fe_affine_operator%get_translation(), computed_solution_vector)
     call iterative_linear_solver%free() 

     !select type(computed_solution_vector)
     !  class is(serial_scalar_array_t)
     !  call computed_solution_vector%print(6)
     !  class is(serial_block_array_t)
     !  call computed_solution_vector%print(6)
     !  class default
     !  check(.false.) 
     !end select
  
     computed_solution_vector = computed_solution_vector - exact_solution_vector
     check ( computed_solution_vector%nrm2()/exact_solution_vector%nrm2() < 1.0e-04 )
     
     call computed_solution_vector%free()
     deallocate(computed_solution_vector)

     call exact_solution_vector%free()
     deallocate(exact_solution_vector)
     
     call fe_affine_operator%free()
     call fe_space%free()
     
     call reference_fe_array(1)%free()
     call reference_fe_array(2)%free()
  end subroutine test_composite_reference_fe_block
  
end program test_reference_fe
