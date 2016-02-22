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

  type test_heterogeneous_poisson_params_t
     
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
					
     ! Problem parameters
     character(len=:), allocatable :: default_number_of_values
     character(len=:), allocatable :: default_temperature_points
     character(len=:), allocatable :: default_conductivity_values
          
  end type test_heterogeneous_poisson_params_t

  ! Types
  public :: test_heterogeneous_poisson_params_t

  ! Functions
  public :: set_default_params, cli_add_params

contains

  subroutine set_default_params(params)
    implicit none
    type(test_heterogeneous_poisson_params_t), intent(inout) :: params

    ! IO parameters
    params%default_dir_path     = 'data'
    params%default_prefix       = 'UnitSquare_Heterogeneous'
    params%default_dir_path_out = 'output'

    ! Problem parameters
    params%default_number_of_values = '2'
    params%default_temperature_points = '0.0 100.0'
    params%default_conductivity_values = '1.0 2.0'
				
  end subroutine set_default_params
  !==================================================================================================

  subroutine cli_add_params(cli,params,group)
    implicit none
    type(Type_Command_Line_Interface)        , intent(inout) :: cli
    type(test_heterogeneous_poisson_params_t), intent(in)    :: params
    character(*)                             , intent(in)    :: group
    ! Locals
    integer(ip) :: error

    ! Set Command Line Arguments
    ! IO parameters
    call cli%add(group=trim(group),switch='--dir_path',switch_ab='-d',                              &
         &       help='Directory of the source files',required=.false., act='store',                &
         &       def=trim(params%default_dir_path),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--prefix',switch_ab='-pr',help='Name of the GiD files',  &
         &       required=.false.,act='store',def=trim(params%default_prefix),error=error) 
    check(error==0)
    call cli%add(group=trim(group),switch='--dir_path_out',switch_ab='-out',help='Output Directory',&
         &       required=.false.,act='store',def=trim(params%default_dir_path_out),error=error)
    check(error==0)

    ! Problem parameters
    call cli%add(group=trim(group),switch='--number_of_values',switch_ab='-nvalu',       &
         &       help='Number of values',required=.false.,act='store',                   &
         &       def=trim(params%default_number_of_values),error=error)
    call cli%add(group=trim(group),switch='--temperature_points',switch_ab='-temps',     &
         &       help='Temperature points',required=.false.,act='store',nargs='*',       &
         &       def=trim(params%default_temperature_points),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--conductivity_values',switch_ab='-conds',    &
         &       help='Conductivity values',required=.false.,act='store',nargs='*',      &
         &       def=trim(params%default_conductivity_values),error=error)
    check(error==0)
				
  end subroutine cli_add_params
  !==================================================================================================

end module command_line_parameters_names

#include "sbm_heterogeneous_poisson_discrete_integration.i90"

!****************************************************************************************************
!****************************************************************************************************

program test_heterogeneous_poisson
  use serial_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
  ! SB
  !use reference_face_names
  use reference_fe_names
  use reference_fe_factory_names
  use SB_fe_space_names
  use SB_discrete_integration_names
  use heterogeneous_poisson_discrete_integration_names
  use SB_fe_affine_operator_names
  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                         :: f_mesh
  type(conditions_t)                   :: f_cond
  type(conditions_t)                   :: f_cond_tri
  type(triangulation_t)                :: f_trian
		
  class(vector_t), allocatable, target :: dof_values, residual ! dof-stored

  type(linear_solver_t)                :: linear_solver
  type(vector_space_t) , pointer       :: fe_affine_operator_range_vector_space
  type(serial_environment_t)           :: senv

  ! Arguments
  character(len=256)                   :: dir_path, dir_path_out
  character(len=256)                   :: prefix, filename

  logical                              :: diagonal_blocks_symmetric_storage(1)
  logical                              :: diagonal_blocks_symmetric(1)
  integer(ip)                          :: diagonal_blocks_sign(1)

  integer(ip)                          :: lunio, istat

  type(Type_Command_Line_Interface)    :: cli 
  character(len=:), allocatable        :: group

  ! SB
  type(SB_serial_fe_space_t)           :: fe_space
  type(heterogeneous_poisson_discrete_integration_t) :: heterogeneous_poisson_integration
  type(SB_fe_affine_operator_t)        :: fe_affine_operator
  type(p_reference_fe_t)               :: reference_fe_array(2)
  
  integer(ip)                          :: count, max_number_iterations
  real(rp)                             :: tol, residual_nrm2
  
  count = 1
  max_number_iterations = 30
  tol = 1.0e-03_rp
  residual_nrm2 = 1.0_rp

  call meminit

  ! Read IO parameters
  call read_flap_cli_test_cdr(cli)
  call cli%parse(error=istat)
  if(cli%run_command('analytical')) then
     group = 'analytical'
  else
     group = 'analytical'
  end if
		
  call cli%get(group=trim(group),switch='-d',val=dir_path,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-pr',val=prefix,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-out',val=dir_path_out,error=istat); check(istat==0)

  ! Read mesh
  call mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read scalar conditions 
  call conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)
  
  ! Create vector conditions
  call conditions_create(3,3,f_cond%ncond,f_cond_tri)
  f_cond_tri%code(1,1:f_cond%ncond) = f_cond%code(1,1:f_cond%ncond)
  f_cond_tri%code(2,1:f_cond%ncond) = f_cond%code(1,1:f_cond%ncond)
  f_cond_tri%code(3,1:f_cond%ncond) = f_cond%code(1,1:f_cond%ncond)
  f_cond_tri%valu(1,1:f_cond%ncond) = f_cond%valu(1,1:f_cond%ncond)
  f_cond_tri%valu(2,1:f_cond%ncond) = f_cond%valu(1,1:f_cond%ncond)
  f_cond_tri%code(3,1:f_cond%ncond) = f_cond%code(1,1:f_cond%ncond)
		
  ! Read temperature - thermal conductivity mapping
  call cli%get(group=trim(group), &
               switch='-nvalu', &
               val=heterogeneous_poisson_integration%nvalu, &
               error=istat); check(istat==0)
  call cli%get_varying(group=trim(group), &
                       switch='-temps', &
                       val=heterogeneous_poisson_integration%data_points, &
                       error=istat); check(istat==0)
  call cli%get_varying(group=trim(group), &
                       switch='-conds', &
                       val=heterogeneous_poisson_integration%property_values, &
                       error=istat); check(istat==0) 

  ! Construct triangulation
  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond_tri )
  
  ! Scalar-vector problem composite space
  reference_fe_array(1) =  make_reference_fe ( topology = topology_quad, &
                                               fe_type = fe_type_lagrangian, &
                                               number_dimensions = 2, &
                                               order = 1, &
                                               field_type = field_type_scalar, &
                                               continuity = .true. )

  reference_fe_array(2) =  make_reference_fe ( topology = topology_quad, &
                                               fe_type = fe_type_lagrangian, &
                                               number_dimensions = 2, &
                                               order = 1, &
                                               field_type = field_type_vector, &
                                               continuity = .true. )
     
  call fe_space%create( triangulation = f_trian, &
                        boundary_conditions = f_cond_tri, &
                        reference_fe_phy = reference_fe_array, &
                        field_blocks = (/1,2/), &
                        field_coupling = reshape((/.true.,.false.,.false.,.true./),(/2,2/)) )

  call fe_space%fill_dof_info() 
    
  call fe_affine_operator%create ( sparse_matrix_storage_format= 'CSR', &
                                   diagonal_blocks_symmetric_storage=(/.true./), &
                                   diagonal_blocks_symmetric=(/.true./), &
                                   diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
                                   triangulation=f_trian, &
                                   fe_space=fe_space, &
                                   discrete_integration=heterogeneous_poisson_integration )
  
  fe_affine_operator_range_vector_space => fe_affine_operator%get_range_vector_space()
  call fe_affine_operator_range_vector_space%create_vector(dof_values)
  heterogeneous_poisson_integration%dof_values => dof_values
  call dof_values%init(0.0_rp)
  call fe_affine_operator_range_vector_space%create_vector(residual)   
  
  call fe_affine_operator%symbolic_setup()
  call fe_affine_operator%numerical_setup()
  

  
  do  while (residual_nrm2 > tol .and. count <= max_number_iterations)
     
     call linear_solver%create(senv)
     call linear_solver%set_type_and_parameters_from_pl()
     call linear_solver%set_operators(fe_affine_operator, .identity. fe_affine_operator)
     call linear_solver%set_initial_solution(dof_values)
     call linear_solver%solve(dof_values)
     call linear_solver%free()
     
     call fe_affine_operator%free_in_stages(free_numerical_setup)
     call fe_affine_operator%numerical_setup()
     
     call fe_affine_operator%apply(dof_values,residual)
     residual_nrm2 = residual%nrm2()
     count = count + 1
     write(*,*) 'Residual:',residual_nrm2
     
  end do
  
  write(*,*) 'Fixed-point algorithm converged in',count,'iterations'
  
  select type(dof_values)
     class is(serial_scalar_array_t)
     call dof_values%print(6)
     class is(serial_block_array_t)
     call dof_values%print(6)
     class default
     check(.false.) 
  end select
		
  call dof_values%free()
  call residual%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call triangulation_free(f_trian)
  call conditions_free(f_cond)
  call conditions_free(f_cond_tri)
  call mesh_free(f_mesh)

  call memstatus 

contains
  
  !==================================================================================================
  subroutine read_flap_cli_test_cdr(cli)
    implicit none
    type(Type_Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(test_heterogeneous_poisson_params_t) :: test_params
    logical                                   :: authors_print
    integer(ip)                               :: error

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'test_heterogeneous_poisson',                                       &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         & 'Serial FEMPAR driver for multimaterial poisson problems using Continuous-Galerkin .',   &
         &        examples    = ['test_heterogeneous_poisson            -h ',                       &
                                 'test_heterogeneous_poisson analytical -h ' ]) ! ALERT: spaces

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical dof_values')

    ! Set Command Line Arguments for each group
    call set_default_params(test_params)
    call cli_add_params(cli,test_params,'analytical')

  end subroutine read_flap_cli_test_cdr
		
end program test_heterogeneous_poisson
