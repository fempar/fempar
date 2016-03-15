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

module heterogeneous_poisson_discrete_integration_names
  use serial_names
  
  implicit none
# include "debug.i90"
  private

  type, extends(discrete_integration_t) :: heterogeneous_poisson_discrete_integration_t
     ! User-defined table
     integer(ip)                :: nvalu
     real(rp), allocatable      :: data_points(:)
     real(rp), allocatable      :: property_values(:)
     class(fe_function_t), pointer   :: dof_values => NULL()                    
   contains
     procedure                  :: integrate
     procedure, non_overridable :: interpolate_property
     procedure, non_overridable :: search_interval
  end type heterogeneous_poisson_discrete_integration_t
  
  public :: heterogeneous_poisson_discrete_integration_t
  
contains
  
  !===============================================================================================
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(heterogeneous_poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                        , intent(inout) :: fe_space
    class(assembler_t)                              , intent(inout) :: assembler
    ! Locals
    type(finite_element_t), pointer    :: fe
    type(volume_integrator_t), pointer :: vol_int_first_fe, vol_int_second_fe
    real(rp), allocatable                 :: elmat(:,:), elvec(:)
    type(fe_map_t), pointer               :: fe_map
    type(quadrature_t), pointer        :: quad

    integer(ip)                           :: igaus,inode,jnode,ioffset,joffset,ngaus
    real(rp)                              :: factor

    type(vector_field_t)                  :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t)                  :: grad_test_vector, grad_trial_vector

    integer(ip)                           :: number_fe_spaces

    integer(ip), pointer                  :: field_blocks(:)
    logical, pointer                      :: field_coupling(:,:)

    integer(ip)                           :: ielem, number_nodes
    type(i1p_t), pointer                  :: elem2dof(:)
    type(i1p_t), pointer                  :: bc_code(:)
    type(r1p_t), pointer                  :: bc_value(:)
    integer(ip), allocatable              :: number_nodes_per_field(:)
    
    type(fe_function_scalar_t)            :: fe_unknown_scalar
    real(rp)                              :: igaus_value_scalar
    
    type(fe_function_vector_t)            :: fe_unknown_vector
    type(vector_field_t)                  :: igaus_value_vector
    
    real(rp)                              :: viscosity_scalar
    type(tensor_field_t)                  :: viscosity_matrix

    class(vector_t), pointer              :: free_dof_values
				
    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks    => fe_space%get_field_blocks()
    field_coupling  => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )

    call fe_space%initialize_integration()
    
    quad => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    
    call fe_space%create_fe_function(1,fe_unknown_scalar)
    call fe_space%create_fe_function(2,fe_unknown_vector)
    
    free_dof_values => this%dof_values%free_dof_values
    
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
						
       call fe%update_values(fe_unknown_scalar, free_dof_values)
       call fe%update_values(fe_unknown_vector, free_dof_values)
       
       do igaus = 1,ngaus
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          
          call fe_unknown_scalar%get_value(igaus, igaus_value_scalar)
          viscosity_scalar = this%interpolate_property(igaus_value_scalar)
          
          do inode = 1, number_nodes_per_field(1)
             call vol_int_first_fe%get_gradient(inode,igaus,grad_trial_scalar)
             do jnode = 1, number_nodes_per_field(1)
                call vol_int_first_fe%get_gradient(jnode,igaus,grad_test_scalar)
                elmat(inode,jnode) = elmat(inode,jnode) + & 
                   & factor * viscosity_scalar * grad_test_scalar * grad_trial_scalar
             end do
          end do
          
          call fe_unknown_vector%get_value(igaus, igaus_value_vector)
          
          call viscosity_matrix%init(0.0_rp)
          call viscosity_matrix%set(1,1,this%interpolate_property(igaus_value_vector%get(1)))
          call viscosity_matrix%set(2,2,this%interpolate_property(igaus_value_vector%get(2)))
          
          do inode = 1, number_nodes_per_field(2)
             ioffset = number_nodes_per_field(1) + inode
             call vol_int_second_fe%get_gradient(inode,igaus,grad_trial_vector)
             do jnode = 1, number_nodes_per_field(2)
                joffset = number_nodes_per_field(1) + jnode
                call vol_int_second_fe%get_gradient(jnode,igaus,grad_test_vector)
                elmat(ioffset,joffset) = elmat(ioffset,joffset) + & 
                   & factor * double_contract(grad_test_vector,viscosity_matrix*grad_trial_vector)
             end do
          end do
          
       end do
       
       ! Apply boundary conditions
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )
    end do
    
    call fe_unknown_scalar%free()
    call fe_unknown_vector%free()
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
 
  !===============================================================================================
  function interpolate_property ( this, point ) result ( value )
     implicit none
     class(heterogeneous_poisson_discrete_integration_t), intent(in) :: this
     real(rp), intent(in) :: point
     real(rp)             :: value
     ! Locals
     integer(ip)          :: ileft, iright
     real(rp)             :: slope
     
     ileft = 1; iright = this%nvalu
     if (point <= this%data_points(ileft)) then
        value = this%property_values(ileft)
     else if (point >= this%property_values(iright)) then
        value = this%property_values(iright)
     else
        call this%search_interval( point, ileft, iright )
     endif
     
     slope = (this%property_values(iright)-this%property_values(ileft))/ &
             &  (this%data_points(iright)-this%data_points(ileft))
     
     value = slope*(point-this%data_points(ileft)) + this%property_values(ileft)
     
  end function interpolate_property
  
  !===============================================================================================
  recursive subroutine search_interval ( this, point, ileft, iright )
     implicit none
     class(heterogeneous_poisson_discrete_integration_t), intent(in) :: this
     real(rp)   , intent(in)    :: point
     integer(ip), intent(inout) :: ileft, iright
     ! Locals
     integer(ip)                :: imid
     
     if (iright - ileft > 1) then
        imid = (ileft + iright)/2
        if (this%data_points(imid) < point) then
           ileft = imid
           call this%search_interval( point, ileft, iright )
        else
           iright = imid
           call this%search_interval( point, ileft, iright )
        endif 
     endif
     
  end subroutine search_interval
  
end module heterogeneous_poisson_discrete_integration_names
!****************************************************************************************************

program test_heterogeneous_poisson
  use serial_names
  use command_line_parameters_names
  use heterogeneous_poisson_discrete_integration_names
  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                         :: f_mesh
  type(conditions_t)                   :: f_cond
  type(conditions_t)                   :: f_cond_tri
  type(triangulation_t)                :: f_trian
		
  class(fe_function_t), allocatable, target :: dof_values
  class(vector_t), allocatable, target      :: residual ! dof-stored

  type(iterative_linear_solver_t)                :: linear_solver
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
  type(serial_fe_space_t)           :: fe_space
  type(heterogeneous_poisson_discrete_integration_t) :: heterogeneous_poisson_integration
  type(fe_affine_operator_t)        :: fe_affine_operator
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
                                   environment=senv, &
                                   fe_space=fe_space, &
                                   discrete_integration=heterogeneous_poisson_integration )
  
  call fe_affine_operator%create_range_vector(dof_values%get_vector_dof_values())
  
  ! It must be a FE function now
  heterogeneous_poisson_integration%dof_values => dof_values
  call dof_values%init(0.0_rp)
  call fe_affine_operator%create_range_vector(residual)
  
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
  call reference_fe_array(1)%free()
  call reference_fe_array(2)%free()
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
