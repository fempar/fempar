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
module command_line_parameters_names
  use types_names
  use Data_Type_Command_Line_Interface
# include "debug.i90"

  implicit none
  private

  type test_reference_fe_params_t
     character(len=:), allocatable     :: default_dir_path
     character(len=1024)               :: dir_path
     character(len=:), allocatable     :: default_prefix
     character(len=2014)               :: prefix
     character(len=:), allocatable     :: default_order
     integer(ip)                       :: order
     character(len=:), allocatable     :: default_laplacian_type
     character(len=2014)               :: laplacian_type
     
     type(Type_Command_Line_Interface) :: cli
  contains
     procedure, private :: set_defaults => test_reference_fe_params_set_defaults
     procedure          :: create       => test_reference_fe_params_create
     procedure          :: parse        => test_reference_fe_params_parse
     procedure          :: free         => test_reference_fe_params_free
  end type test_reference_fe_params_t

  ! Types
  public :: test_reference_fe_params_t

contains

  subroutine test_reference_fe_params_set_defaults(this)
    implicit none
    class(test_reference_fe_params_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = 'square_4x4_scalar.gid/data/'
    this%default_prefix         = 'square_4x4_scalar'
    this%default_order          = '1' 
    this%default_laplacian_type = 'scalar'
  end subroutine test_reference_fe_params_set_defaults
  !==================================================================================================

  subroutine test_reference_fe_params_create(this)
    implicit none
    class(test_reference_fe_params_t), intent(inout) :: this
    ! Locals
    integer(ip) :: error

    call this%set_defaults()
    
    ! Initialize Command Line Interface
    call this%cli%init(progname    = 'test_reference_fe', &
         &               version     = '',                               &
         &               authors     = '',                               &
         &               license     = '',                               &
         &               description = "FEMPAR sequential test driver to solve a Scalar or Vector-Valued Laplacian problem with&
                                     & Galerkin FEM, (arbitrary order) continuous Lagrangian FEs, and a Krylov solver without preconditioning. The&
                                     & GiD input data should be set-up such that the solution is u=1 in the whole domain.", &
         &               examples    = ['test_reference_fe -h'] )
    
    
    ! Set Command Line Arguments
    ! IO parameters
    call this%cli%add(switch='--dir_path',switch_ab='-d',                              &
         &            help='Absolute or relative PATH to the data files. Must end with /',required=.false., act='store', &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Prefix for all input files (mesh, conditions, etc.). & 
                       & E.g., if these files were generated from square.gid GiD project, then --prefix square.', &
         &             required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--order',switch_ab='-o',help='Polynomial order for Lagrangian FEs',  &
         &       required=.false.,act='store',def=trim(this%default_order), error=error) 
    check(error==0)
    
    call this%cli%add(switch='--laplacian_type',switch_ab='-lt',help='Scalar or Vector Laplacian.',  &
         &       required=.false.,choices='scalar,vector', act='store',def=trim(this%default_laplacian_type), error=error) 
    check(error==0)

  end subroutine test_reference_fe_params_create
  
    subroutine test_reference_fe_params_parse(this)
    implicit none
    class(test_reference_fe_params_t), intent(inout) :: this
    ! Locals
    integer(ip) :: error
    call this%cli%parse(error=error)
    check(error==0)
    
    call this%cli%get(switch='-d',val=this%dir_path,error=error); check(error==0)
    call this%cli%get(switch='-p',val=this%prefix,error=error); check(error==0)
    call this%cli%get(switch='-o',val=this%order,error=error); check(error==0)
    call this%cli%get(switch='-lt',val=this%laplacian_type,error=error); check(error==0)
    
  end subroutine test_reference_fe_params_parse

  subroutine test_reference_fe_params_free(this)
    implicit none
    class(test_reference_fe_params_t), intent(inout) :: this
    ! IO parameters
    deallocate(this%default_dir_path)
    deallocate(this%default_prefix)
    deallocate(this%default_order)
    ! How to free the cli?
  end subroutine test_reference_fe_params_free
  
  
end module command_line_parameters_names

module poisson_discrete_integration_names
  use serial_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_discrete_integration_t
     integer(ip) :: viscosity 
   contains
     procedure :: integrate
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
  
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)          , intent(inout) :: fe_space
    class(assembler_t)                , intent(inout) :: assembler

    type(finite_element_t), pointer :: fe
    type(volume_integrator_t), pointer :: vol_int
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(quadrature_t), pointer :: quad

    integer(ip)  :: igaus,inode,jnode,ngaus
    real(rp)     :: factor

    type(vector_field_t) :: grad_test, grad_trial

    integer(ip) :: number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, iapprox, number_nodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)
    integer(ip), allocatable :: number_nodes_per_field(:)  

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )

    call fe_space%initialize_integration()
    
    quad => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map   => fe%get_fe_map()
       vol_int  => fe%get_volume_integrator(1)
       elem2dof => fe%get_elem2dof()
       bc_code  => fe%get_bc_code()
       bc_value => fe%get_bc_value()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes
             call vol_int%get_gradient(inode,igaus,grad_trial)
             do jnode = 1, number_nodes
                call vol_int%get_gradient(jnode,igaus,grad_test)
                elmat(inode,jnode) = elmat(inode,jnode) + factor * grad_test * grad_trial
             end do
          end do
       end do
       
       ! Apply boundary conditions
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
end module poisson_discrete_integration_names

module vector_laplacian_single_discrete_integration_names
use serial_names

implicit none
# include "debug.i90"
private
type, extends(discrete_integration_t) :: vector_laplacian_single_discrete_integration_t
integer(ip) :: viscosity 
contains
procedure :: integrate
end type vector_laplacian_single_discrete_integration_t

public :: vector_laplacian_single_discrete_integration_t

contains
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(vector_laplacian_single_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                   , intent(inout) :: fe_space
    class(assembler_t)                         , intent(inout) :: assembler

    type(finite_element_t), pointer :: fe
    type(volume_integrator_t), pointer :: vol_int
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(quadrature_t), pointer :: quad
    integer(ip), allocatable :: number_nodes_per_field(:)

    integer(ip)  :: igaus,inode,jnode,ioffset,joffset,ngaus
    real(rp) :: factor

    type(vector_field_t) :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t) :: grad_test_vector, grad_trial_vector
    
    integer(ip) :: i, number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, number_nodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )
    
    call fe_space%initialize_integration()
    
    quad  => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map            => fe%get_fe_map()
       vol_int           => fe%get_volume_integrator(1)
       elem2dof          => fe%get_elem2dof()
       bc_code           => fe%get_bc_code()
       bc_value          => fe%get_bc_value()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes
             call vol_int%get_gradient(inode,igaus,grad_trial_vector)
             do jnode = 1, number_nodes
                call vol_int%get_gradient(jnode,igaus,grad_test_vector)
                elmat(inode,jnode) = elmat(inode,jnode) + factor * double_contract(grad_test_vector,grad_trial_vector)
             end do
          end do
       end do
       
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )      
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate

end module vector_laplacian_single_discrete_integration_names

module vector_laplacian_composite_discrete_integration_names
use serial_names

implicit none
# include "debug.i90"
private
type, extends(discrete_integration_t) :: vector_laplacian_composite_discrete_integration_t
integer(ip) :: viscosity 
contains
procedure :: integrate
end type vector_laplacian_composite_discrete_integration_t

public :: vector_laplacian_composite_discrete_integration_t

contains
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(vector_laplacian_composite_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                   , intent(inout) :: fe_space
    class(assembler_t)                         , intent(inout) :: assembler

    type(finite_element_t), pointer :: fe
    type(volume_integrator_t), pointer :: vol_int_first_fe, vol_int_second_fe
    real(rp), allocatable :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer :: fe_map
    type(quadrature_t), pointer :: quad
    integer(ip), allocatable :: number_nodes_per_field(:)

    integer(ip)  :: igaus,inode,jnode,ioffset,joffset,ngaus
    real(rp) :: factor

    type(vector_field_t) :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t) :: grad_test_vector, grad_trial_vector
    
    integer(ip) :: i, number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical, pointer :: field_coupling(:,:)

    integer(ip) :: ielem, number_nodes
    type(i1p_t), pointer :: elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )
    
    call fe_space%initialize_integration()
    
    quad  => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    do ielem = 1, fe_space%get_number_elements()
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       
       fe_map            => fe%get_fe_map()
       vol_int_first_fe  => fe%get_volume_integrator(1)
       vol_int_second_fe => fe%get_volume_integrator(2)
       elem2dof          => fe%get_elem2dof()
       bc_code           => fe%get_bc_code()
       bc_value          => fe%get_bc_value()

       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes_per_field(1)
             call vol_int_first_fe%get_gradient(inode,igaus,grad_trial_scalar)
             do jnode = 1, number_nodes_per_field(1)
                call vol_int_first_fe%get_gradient(jnode,igaus,grad_test_scalar)
                elmat(inode,jnode) = elmat(inode,jnode) + factor * grad_test_scalar * grad_trial_scalar
             end do
          end do

          do inode = 1, number_nodes_per_field(2)
             ioffset = number_nodes_per_field(1)+inode
             call vol_int_second_fe%get_gradient(inode,igaus,grad_trial_scalar)
             ! write(*,*) inode, grad_trial_vector%value
             do jnode = 1, number_nodes_per_field(2)
                joffset = number_nodes_per_field(1)+jnode
                call vol_int_second_fe%get_gradient(jnode,igaus,grad_test_scalar)
                elmat(ioffset,joffset) = elmat(ioffset,joffset) + factor * grad_test_scalar * grad_trial_scalar
             end do
          end do
       end do
       
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )      
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate

end module vector_laplacian_composite_discrete_integration_names

!****************************************************************************************************
program test_reference_fe
  use serial_names
  use command_line_parameters_names
  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                          :: f_mesh
  type(triangulation_t)                 :: f_trian
  type(conditions_t)                    :: f_cond
  type(test_reference_fe_params_t)      :: params
  integer(ip) :: problem_id
 
  call meminit

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
    call test_single_vector_valued_reference_fe()
    call test_composite_reference_fe_monolithic()
    call test_composite_reference_fe_block()
  end if  
 
  !call fe_space%free()
  call triangulation_free(f_trian)
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)
  call params%free()

  call memstatus

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
     reference_fe_array(1) =  make_reference_fe ( topology = topology_quad, &
                                                  fe_type = fe_type_lagrangian, &
                                                  number_dimensions = f_trian%num_dims, &
                                                  order = params%order, &
                                                  field_type = field_type_scalar, &
                                                  continuity = .true. )
     
     call fe_space%create( triangulation = f_trian, &
                           boundary_conditions = f_cond, &
                           reference_fe_phy = reference_fe_array )
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format='CSR', &
                                     diagonal_blocks_symmetric_storage=(/.true./), &
                                     diagonal_blocks_symmetric=(/.true./), &
                                     diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                     environment=senv, &
                                     fe_space=fe_space, &
                                     discrete_integration=poisson_integration )
     
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
     call iterative_linear_solver%set_type_and_parameters_from_pl()
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(computed_solution_vector)
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
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format='CSR', &
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
     call iterative_linear_solver%set_type_and_parameters_from_pl()
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(computed_solution_vector)
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
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create (sparse_matrix_storage_format='CSR', &
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
     call iterative_linear_solver%set_type_and_parameters_from_pl()
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(computed_solution_vector)
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
     
     call fe_space%fill_dof_info() 
     
     call fe_affine_operator%create ( 'CSR', &
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
     call iterative_linear_solver%set_type_from_string("CG")
     call iterative_linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call iterative_linear_solver%solve(computed_solution_vector)
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
