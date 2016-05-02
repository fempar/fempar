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

!***************************************************************************************************!
! COMMAND LINE PARAMETERS                                                                           !
! Functions to read the parameters introduced in the execution command line                         !
!***************************************************************************************************!
# include "sbm_cdr_command_line_parameters.i90"

!***************************************************************************************************!
! ANALYTICAL LINEAR-STEADY FUNCTION                                                                 !
! Definition of the linear-steady analytical function for the CDR problem.                          !
!***************************************************************************************************!
# include "sbm_cdr_analytical_linear_solution.i90"

!***************************************************************************************************!
! ANALYTICAL LINEAR-STEADY-TRANSIENT FUNCTION                                                       !
! Definition of the linear-steady analytical function for the CDR problem.                          !
!***************************************************************************************************!
# include "sbm_cdr_analytical_linear_transient_solution.i90"

!***************************************************************************************************!
! ANALYTICAL LINEAR-STEADY FUNCTION                                                                 !
! Definition of the linear-steady analytical function for the CDR problem.                          !
!***************************************************************************************************!
# include "sbm_cdr_analytical_layer_solution.i90"

!***************************************************************************************************!
! ANALYTICAL FUNCTIONS                                                                              ! 
! Definition of the analytical functions for the CDR problem.                                       !
!***************************************************************************************************!
# include "sbm_cdr_analytical.i90"

!***************************************************************************************************!
! TIME INTEGRATION                                                                                  !
! Theta method                                                                                      !
!***************************************************************************************************!
# include "sbm_theta_method.i90"

!***************************************************************************************************!
! DISCRETE INTEGRATION                                                                              !
! The method integration: cG/dG + theta-method in time for convection-diffusion problems            !
!***************************************************************************************************!
module dG_CDR_discrete_integration_names
 use serial_names
 use theta_method_names
 use cdr_analytical_functions_names 
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: dG_CDR_discrete_integration_t
     real(rp)                         :: viscosity 
     real(rp)                         :: C_IP ! Interior Penalty constant
     ! DG symmetric parameter: (0-Symmetric, 1/2-Incomplete, 1-Antisymmetric)
     real(rp)                         :: xi   
     class(fe_function_t)   , pointer :: fe_function         => NULL()
     class(theta_method_t)  , pointer :: theta_method        => NULL()
     type(CDR_analytical_functions_t) :: analytical_functions
     type(fe_function_t)    , pointer :: fe_values_previous 
   contains
     procedure, non_overridable :: set_problem
     procedure, non_overridable :: integrate
     procedure, non_overridable :: compute_analytical_force
     procedure, non_overridable :: set_initial_solution
     procedure, non_overridable :: compute_graph_laplacian_perturbed_matrix
     !procedure, non_overridable :: compute_graph_laplacian_perturbed_block_matrix
  end type DG_CDR_discrete_integration_t
  
  public :: dG_CDR_discrete_integration_t
  
contains
  
  subroutine set_problem(this, viscosity, C_IP, xi, solution_field_name)
    implicit none
    class(dG_CDR_discrete_integration_t)   , intent(inout) :: this
    real(rp)                               , intent(in)    :: viscosity 
    real(rp)                               , intent(in)    :: C_IP 
    real(rp)                               , intent(in)    :: xi
    character(len=1024)                    , intent(in)    :: solution_field_name

    this%viscosity = viscosity
    this%C_IP      = C_IP
    this%xi        = xi
    call this%analytical_functions%create(solution_field_name)
  end subroutine set_problem

  !==================================================================================================
  subroutine set_initial_solution(this,fe_space)
    implicit none
    class(dG_CDR_discrete_integration_t), intent(inout) :: this
    class(serial_fe_space_t)            , intent(in)    :: fe_space

    ! Update unknown 
    call fe_space%interpolate_fe_function(this%analytical_functions%get_solution_function(),        &
         &                                fe_space_component=1,fe_function=this%fe_values_previous, &
         &                                time=this%theta_method%initial_time)

  end subroutine set_initial_solution
  
  !==================================================================================================
  function compute_analytical_force(this,point) result(source)
    implicit none
    class(dG_cdr_discrete_integration_t), intent(in) :: this
    type(point_t)                       , intent(in) :: point
    real(rp) :: source
    real(rp) :: laplacian
    type(vector_field_t) :: gradient
    type(vector_field_t) :: convection
    real(rp) :: dt_solution
    real(rp)             :: time

    time = this%theta_method%current_time
    call this%analytical_functions%get_value_solution_laplacian(point,time,laplacian)
    call this%analytical_functions%get_value_solution_gradient(point,time,gradient)
    call this%analytical_functions%get_value_convection(point,time,convection)
    call this%analytical_functions%get_value_dt_solution(point,time,dt_solution)
    
    source = dt_solution - this%viscosity*laplacian + convection*gradient 
  end function compute_analytical_force

  !==================================================================================================
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(dG_CDR_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)            , intent(inout) :: fe_space
    class(assembler_t)                  , intent(inout) :: assembler

    type(finite_element_t)   , pointer :: fe
    type(finite_face_t)      , pointer :: face
    type(volume_integrator_t), pointer :: vol_int
    type(face_integrator_t)  , pointer :: face_int
    real(rp)             , allocatable :: elmat(:,:), elvec(:), facemat(:,:,:,:), facevec(:,:)
    type(fe_map_t)           , pointer :: fe_map
    type(face_map_t)         , pointer :: face_map
    type(quadrature_t)       , pointer :: quad

    integer(ip)                :: igaus,inode,jnode,ngaus
    real(rp)                   :: factor, h_length, bcvalue, source, outflow, time, time_factor
    real(rp)                   :: u_value_previous
    type(fe_function_scalar_t) :: fe_unknown_current, fe_unknown_previous

    real(rp)                   :: shape_test_scalar, shape_trial_scalar
    type(vector_field_t)       :: grad_test_scalar, grad_trial_scalar
    type(vector_field_t)       :: normal(2)
    type(vector_field_t)       :: convection

    integer(ip)                :: i, j, number_fe_spaces

    integer(ip)      , pointer :: field_blocks(:)
    logical          , pointer :: field_coupling(:,:)

    integer(ip)                :: ielem, iface, iapprox, number_nodes, ineigh, jneigh
    integer(ip)                :: number_neighbours
    type(i1p_t)      , pointer :: elem2dof(:),test_elem2dof(:),trial_elem2dof(:)
    type(i1p_t)      , pointer :: bc_code(:)
    type(r1p_t)      , pointer :: bc_value(:)
    integer(ip)  , allocatable :: number_nodes_per_field(:)
    type(point_t)    , pointer :: coordinates(:)

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()

    time = this%theta_method%current_time
    time_factor = this%theta_method%time_factor

    ! Initialize the integration values
    call fe_space%initialize_integration()

    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )

    ! Fill the previous array
    call fe_space%update_bc_value_scalar(this%analytical_functions%get_solution_function(),         & 
                            bc_code=1, fe_space_component=1, time=time-this%theta_method%time_step )
    call fe_space%update_global_fe_function_bcs(this%fe_values_previous)
    
    ! Initialize the finite element function for each kind of reference element
    call fe_space%create_fe_function(1,fe_unknown_previous)

    ! Update the Dirichlet boundary conditions
    call fe_space%update_bc_value_scalar(this%analytical_functions%get_solution_function(),         &
         &                               bc_code=1,fe_space_component=1, time= time )
    call fe_space%update_global_fe_function_bcs(this%fe_function)
    ! ------------------------------------ LOOP OVER THE ELEMENTS -----------------------------------
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )

    write(*,*) '*************************  Integration on elements  *************************'
    quad => fe%get_quadrature()
    ngaus = quad%get_number_quadrature_points()
    do ielem = 1, fe_space%get_number_elements()
       !write(*,*) __FILE__,__LINE__,ielem,'------------------'
       elmat = 0.0_rp
       elvec = 0.0_rp

       fe => fe_space%get_finite_element(ielem)
       call fe%update_integration()
       call fe%update_values(fe_unknown_previous, this%fe_values_previous) 

       fe_map   => fe%get_fe_map()
       vol_int  => fe%get_volume_integrator(1)
       elem2dof => fe%get_elem2dof()
       
       coordinates => fe_map%get_quadrature_coordinates() 

       do igaus = 1,ngaus
          source = this%compute_analytical_force(coordinates(igaus))
          call this%analytical_functions%get_value_convection(coordinates(igaus),time,convection)
          factor = fe_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do inode = 1, number_nodes
             call vol_int%get_gradient(inode,igaus,grad_test_scalar)
             call vol_int%get_value(inode,igaus,shape_test_scalar)
             do jnode = 1, number_nodes
                call vol_int%get_gradient(jnode,igaus,grad_trial_scalar)
                call vol_int%get_value(jnode,igaus,shape_trial_scalar)
                elmat(inode,jnode) = elmat(inode,jnode) +                                           &
                     &               factor * (this%viscosity*grad_trial_scalar*grad_test_scalar    &
                     &               - convection*grad_test_scalar*shape_trial_scalar               &
                     &               + time_factor * shape_test_scalar * shape_trial_scalar)
             end do
             call fe_unknown_previous%get_value(igaus, u_value_previous)
             elvec(inode) = elvec(inode) + factor * (source * shape_test_scalar                     &
                  &                      + time_factor * u_value_previous * shape_test_scalar)
          end do
       end do
       
       ! Apply boundary conditions
       call fe%impose_strong_dirichlet_bcs( elmat, elvec)
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,   &
            &                   field_coupling, elmat, elvec )
    end do
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    ! ---------------------------------- END LOOP OVER THE ELEMENTS ---------------------------------

    ! ------------------------------------ LOOP OVER THE FACES --------------------------------------
    call memalloc ( number_nodes, number_nodes,2,2,facemat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, 2,facevec, __FILE__, __LINE__ )

    write(*,*) '*************************   Integration on faces    *************************'
    !call fe_space%initialize_face_integration()
    do iface = 1, fe_space%get_number_interior_faces()

       facemat = 0.0_rp
       facevec = 0.0_rp

       face => fe_space%get_finite_face(iface)
       number_neighbours = face%number_neighbours()
       !write(*,*) __FILE__,__LINE__,iface,'------------------'
       call face%update_integration()
      
       quad   => face%get_quadrature()
       ngaus = quad%get_number_quadrature_points()
       face_map => face%get_map()

       j = 1
       face_int => face%get_face_integrator(j)
       coordinates => face_map%get_quadrature_coordinates()

       do igaus = 1, ngaus
          call face_map%get_normals(igaus,normal)
          h_length = face_map%compute_characteristic_length(igaus,number_neighbours)
          factor = face_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          call this%analytical_functions%get_value_convection(coordinates(igaus),0.0_rp,convection)
          do ineigh = 1, number_neighbours
             do inode = 1, number_nodes_per_field(j)
                !ioffset = number_nodes_per_field(j)*(ineigh-1) + inode
                call face_int%get_value(inode,igaus,ineigh,shape_test_scalar)
                call face_int%get_gradient(inode,igaus,ineigh,grad_test_scalar)
                do jneigh = 1, number_neighbours
                   do jnode = 1, number_nodes_per_field(j)
                      !joffset = number_nodes_per_field(j)*(jneigh-1) + jnode
                      call face_int%get_value(jnode,igaus,jneigh,shape_trial_scalar)
                      call face_int%get_gradient(jnode,igaus,jneigh,grad_trial_scalar)
                      !- mu*({{grad u}}[[v]] + xi*[[u]]{{grad v}} ) + C*mu*p^2/h * [[u]] [[v]]
                      facemat(inode,jnode,ineigh,jneigh) = facemat(inode,jnode,ineigh,jneigh) +     &
                           &  factor * (this%viscosity *                                             &
                           &  (-0.5_rp*grad_trial_scalar*normal(ineigh)*shape_test_scalar           &
                           &   -this%xi*0.5_rp*grad_test_scalar*normal(jneigh)*shape_trial_scalar   &
                           &   + this%c_IP / h_length * shape_trial_scalar*shape_test_scalar *      &
                           &   normal(ineigh)*normal(jneigh))                                       &
                           & + 0.5_rp*convection*normal(ineigh)*shape_test_scalar*shape_trial_scalar&
                           & + 0.5_rp*convection%nrm2()*shape_trial_scalar*shape_test_scalar *      &
                           &   normal(ineigh)*normal(jneigh))                         
                   end do
                end do
             end do
          end do
       end do
       do ineigh = 1, number_neighbours
          trial_elem2dof => face%get_elem2dof(ineigh)
          do jneigh = 1, number_neighbours
             test_elem2dof => face%get_elem2dof(jneigh)
             call assembler%face_assembly(number_fe_spaces,number_nodes_per_field,                  &
                  &                       number_nodes_per_field,trial_elem2dof,test_elem2dof,      &
                  &                       field_blocks,field_coupling,facemat(:,:,ineigh,jneigh),   &
                  &                       facevec(:,ineigh) )   
          end do
       end do
    end do
    
    write(*,*) '************************* Integration on boundaries *************************'
    ! Boundary Faces
    do iface = fe_space%get_number_interior_faces() + 1, fe_space%get_number_interior_faces() +     &
         &                                               fe_space%get_number_boundary_faces()

       facemat = 0.0_rp
       facevec = 0.0_rp

       face => fe_space%get_finite_face(iface)
       number_neighbours = face%number_neighbours()
       !write(*,*) __FILE__,__LINE__,iface,'------------------'
       call face%update_integration()

       face_map => face%get_map()
       quad   => face%get_quadrature()
       ngaus = quad%get_number_quadrature_points()

       j = 1
       face_int => face%get_face_integrator(j)
       coordinates => face_map%get_quadrature_coordinates()

       do igaus = 1, ngaus
          call face_map%get_normals(igaus,normal)
          h_length = face_map%compute_characteristic_length(igaus,number_neighbours)
          factor = face_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          call this%analytical_functions%get_value_solution(coordinates(igaus),time,bcvalue)
          call this%analytical_functions%get_value_convection(coordinates(igaus),time,convection)
          if (convection * normal(1) > 0 ) then
             outflow = 1.0_rp
          else
             outflow = 0.0_rp
          end if
          do ineigh = 1, number_neighbours
             do inode = 1, number_nodes_per_field(j)
                call face_int%get_value(inode,igaus,ineigh,shape_test_scalar)
                call face_int%get_gradient(inode,igaus,ineigh,grad_test_scalar)
                do jneigh = 1, number_neighbours
                   do jnode = 1, number_nodes_per_field(j)
                      call face_int%get_value(jnode,igaus,jneigh,shape_trial_scalar)
                      call face_int%get_gradient(jnode,igaus,jneigh,grad_trial_scalar)
                      facemat(inode,jnode,ineigh,jneigh) = facemat(inode,jnode,ineigh,jneigh) +     &
                           &  factor * (this%viscosity *   &
                           &  (- grad_trial_scalar*normal(ineigh)*shape_test_scalar                 &
                           &  - this%xi*grad_test_scalar*normal(jneigh)*shape_trial_scalar          &
                           &   + this%c_IP / h_length * shape_trial_scalar*shape_test_scalar *      &
                           &   normal(ineigh)*normal(jneigh))                                       &
                           &  + outflow * convection * normal(ineigh) *                             &
                           &  shape_test_scalar*shape_trial_scalar)
                   end do
                end do
                facevec(inode,ineigh) = facevec(inode,ineigh) + factor * (this%viscosity *          &
                     &                  (+ this%xi* bcvalue * grad_test_scalar*normal(jneigh)       &
                     &                  + this%c_IP/h_length * bcvalue * shape_test_scalar )        &
                     &                  - (1.0_rp-outflow) * convection *normal(ineigh) *           &
                     &                  bcvalue*shape_test_scalar)
             end do
          end do
       end do

       do ineigh = 1, number_neighbours
          trial_elem2dof => face%get_elem2dof(ineigh)
          do jneigh = 1, number_neighbours
             test_elem2dof => face%get_elem2dof(jneigh)
             call assembler%face_assembly(number_fe_spaces,number_nodes_per_field,                  &
                  &                       number_nodes_per_field,trial_elem2dof,test_elem2dof,      &
                  &                       field_blocks,field_coupling,facemat(:,:,ineigh,jneigh),   &
                  &                       facevec(:,ineigh) )   
          end do
       end do
    end do
 
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )
    ! ----------------------------------- END LOOP OVER FACES ---------------------------------------

    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
  end subroutine integrate

  subroutine compute_graph_laplacian_perturbed_matrix(this,sparse_matrix)
    use sparse_matrix_names
    use matrix_names
    implicit none
    class(dG_CDR_discrete_integration_t), intent(in)    :: this
    type(sparse_matrix_t)               , intent(inout) :: sparse_matrix

    type(sparse_matrix_iterator_t) :: matrix_entry
    integer(ip) :: i,j
    real(rp)    :: value,transposed_value
    logical     :: finished
    logical     :: found
    real(rp)    :: nu

    call sparse_matrix%get_iterator(matrix_entry)

    do while (.not. matrix_entry%has_finished())

       i = matrix_entry%get_row()
       j = matrix_entry%get_column()

       if (i < j) then
          value = matrix_entry%get_value()
          found = sparse_matrix%get_value(j,i,transposed_value); assert(found)
          nu = max(value,transposed_value)
          if (nu > 0) then
             call matrix_entry%set_value(value-nu)
             found = sparse_matrix%sum_value(j,i,-nu); assert(found)
             found = sparse_matrix%sum_value(i,i,+nu); assert(found)
             found = sparse_matrix%sum_value(j,j,+nu); assert(found)
          end if
       end if

       ! Do stuff
       call matrix_entry%next()
    end do
  end subroutine compute_graph_laplacian_perturbed_matrix
!!$  
!!$  subroutine compute_graph_laplacian_perturbed_block_matrix(this,block_sparse_matrix)
!!$    use block_sparse_matrix_names
!!$    use matrix_names
!!$    implicit none
!!$    class(dG_CDR_discrete_integration_t), intent(in)    :: this
!!$    type(block_sparse_matrix_t)         , intent(inout) :: block_sparse_matrix
!!$
!!$    class(block_sparse_matrix_iterator_t), allocatable :: iterator
!!$    integer(ip) :: i,j,iblock,jblock
!!$    real(rp)    :: value
!!$    logical     :: finished
!!$
!!$    call block_sparse_matrix%get_iterator(iterator)
!!$
!!$    do while (.not. iterator%has_finished())
!!$       iblock = iterator%get_iblock()
!!$       jblock = iterator%get_jblock()
!!$       i = iterator%get_row()
!!$       j = iterator%get_column()
!!$       value = iterator%get_value()
!!$       write(*,*)  iblock,jblock,i, j, value
!!$
!!$       if ((i==j) .and. (iblock == jblock)) then
!!$          call iterator%set_value(1.0_rp)
!!$       else
!!$          call iterator%set_value(0.0_rp)
!!$       end if
!!$          
!!$       ! Do stuff
!!$       call iterator%next()
!!$    end do
!!$
!!$    call iterator%free()
!!$  end subroutine compute_graph_laplacian_perturbed_block_matrix
end module DG_CDR_discrete_integration_names

!****************************************************************************************************
program test_cdr
  use serial_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
  use serial_names
  use theta_method_names
  use dG_CDR_discrete_integration_names
  use sparse_matrix_names
  use direct_solver_names
  use FPL
  use pardiso_mkl_direct_solver_names
  use umfpack_direct_solver_names
  use cdr_analytical_functions_names 
  use vtk_handler_names

  implicit none
#include "debug.i90"

  ! Our data
  type(mesh_t)                          :: f_mesh
  type(triangulation_t)                 :: f_trian
  type(conditions_t)                    :: f_cond
  class(matrix_t)             , pointer :: matrix
  class(array_t)              , pointer :: array
  type(serial_scalar_array_t) , pointer :: my_array
  type(serial_scalar_array_t) , target  :: feunk
  type(serial_environment_t)            :: senv

  type(Type_Command_Line_Interface)     :: cli 
  character(len=:)        , allocatable :: group

  ! Arguments
  character(len=256)       :: dir_path, dir_path_out
  character(len=256)       :: prefix, filename
  integer(ip)              :: i, j, vars_prob(1) = 1, ierror, iblock

  character(len=1024)                           :: space_solution_name
  real(rp)                                      :: viscosity
  logical                                       :: continuity

  integer(ip)                                   :: istat

  class(vector_t)        , allocatable, target  :: residual
  class(vector_t)                     , pointer :: dof_values
  class(vector_t)                     , pointer :: dof_values_previous
  type(fe_function_t)                 , target  :: fe_function
  type(fe_function_t)                 , target  :: fe_values_previous

  type(serial_fe_space_t)                       :: fe_space
  type(p_reference_fe_t)                        :: reference_fe_array_one(1)
  type(p_reference_fe_t)                        :: reference_fe_array_two(2)
  type(fe_affine_operator_t)                    :: fe_affine_operator
  type(theta_method_t)                , target  :: theta_method
  type(dG_CDR_discrete_integration_t)           :: dG_CDR_integration
  class(vector_t)        , allocatable, target  :: vector
  type(interpolation_face_restriction_t)        :: face_interpolation

  class(vector_t)                     , pointer :: rhs

  type(sparse_matrix_t), pointer   :: sparse_matrix
  type(direct_solver_t)            :: direct_solver
  type(ParameterList_t)            :: parameter_list
  type(ParameterList_t), pointer   :: direct_solver_parameters
  integer                          :: FPLError
  integer                          :: iparm(64) = 0
  logical                          :: diagonal_blocks_symmetric_storage(1)
  logical                          :: diagonal_blocks_symmetric(1)
  integer(ip)                      :: diagonal_blocks_sign(1)

  integer(ip)                      :: order, count, max_number_iterations

  real(rp)                         :: tolerance, residual_nrm2

  integer(ip)                      :: err
  type(vtk_handler_t)              :: vtk_handler
  integer(ip)                      :: vtk_error

  call meminit

  ! ParameterList: initialize
  call FPL_Init()
  call the_direct_solver_creational_methods_dictionary%init()
  call parameter_list%Init()
  direct_solver_parameters => parameter_list%NewSubList(Key=pardiso_mkl)
  
  ! ParameterList: set parameters
  FPLError = 0
  FPLError = FPLError +                                                                             &
       &     direct_solver_parameters%set(key = direct_solver_type,        value = pardiso_mkl)
  FPLError = FPLError +                                                                             &
       &      direct_solver_parameters%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_uss)
  FPLError = FPLError +                                                                             &
       &     direct_solver_parameters%set(key = pardiso_mkl_message_level, value = 0)
  FPLError = FPLError +                                                                             &
       &     direct_solver_parameters%set(key = pardiso_mkl_iparm,         value = iparm)
  check(FPLError == 0)
  ! Read IO parameters
  call read_flap_cli_test_cdr(cli)
  call cli%parse(error=istat)
  if(cli%run_command('analytical')) then
     group = 'analytical'
  elseif(cli%run_command('transient')) then
     group = 'transient'
  else
     group = 'analytical'
  end if
  call cli%get(group=trim(group),switch='-d',val=dir_path,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-pr',val=prefix,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-out',val=dir_path_out,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-cg',val= continuity,error=istat); check(istat==0)
  call cli%get(group=trim(group),switch='-p',val=order,error=istat); check(istat==0)

  ! Read mesh
  call mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)
 
  ! Construc triangulation
  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  if (.not. continuity)  call triangulation_construct_faces ( f_trian )

  ! Time integration
  call theta_method%read_and_initialize( cli, group )
  dG_CDR_integration%theta_method => theta_method

  ! Composite case
  reference_fe_array_one(1) = make_reference_fe ( topology = topology_quad,                         &
       &                                          fe_type  = fe_type_lagrangian,                    &
       &                                          number_dimensions = f_trian%num_dims,             &
       &                                          order = order,                                    &
       &                                          field_type = field_type_scalar,                   &
       &                                          continuity = continuity )

  call fe_space%create( triangulation = f_trian, boundary_conditions = f_cond,                      &
       &                reference_fe_phy = reference_fe_array_one,                                  &
       &                field_blocks = (/1/),                                                       &
       &                field_coupling = reshape((/.true./),(/1,1/)) )

  call fe_space%create_face_array()

  call fe_space%fill_dof_info() 
 
  ! Set the analytical solution and convection  
  call cli%get(group=trim(group),switch='-ssol' ,val=space_solution_name,error=istat);check(istat==0)
  call cli%get(group=trim(group),switch='-diff' ,val=viscosity,error=istat);check(istat==0)
  call dG_CDR_integration%set_problem( viscosity = viscosity, C_IP = 10.0_rp, xi = 1.0_rp,             &
       &                               solution_field_name = space_solution_name)

  ! Set the parameters for the integration in time
  dG_CDR_integration%fe_values_previous => fe_values_previous
  call fe_space%create_global_fe_function(dG_CDR_integration%fe_values_previous) 
  call dG_CDR_integration%set_initial_solution(fe_space)
  dof_values_previous => fe_values_previous%get_dof_values()
  !call dof_values_previous%init(0.0_rp)

  ! Create the operator
  diagonal_blocks_symmetric_storage = .false.
  diagonal_blocks_symmetric         = .false.
  diagonal_blocks_sign              = SPARSE_MATRIX_SIGN_INDEFINITE

  call fe_affine_operator%create ('CSR',diagonal_blocks_symmetric_storage ,                         &
       &                          diagonal_blocks_symmetric,diagonal_blocks_sign,                   &
       &                          senv, fe_space, dG_CDR_integration)

  ! Create the unknown array
  call fe_space%create_global_fe_function(fe_function)
  dG_CDR_integration%fe_function => fe_function
  dof_values => fe_function%get_dof_values()
  call dof_values%init(0.0_rp)
  !call fe_affine_operator%create_range_vector(dof_values_previous)
  !call dof_values_previous%init(0.0_rp)
  dG_CDR_integration%fe_function => fe_function

  
  ! Initialize default solver parameters
  tolerance = 1.0e-8
  call fe_affine_operator%create_range_vector(residual)
  max_number_iterations = 1
  
  temporal: do while (.not. theta_method%finished ) 

     call theta_method%print(6)

     ! Initialize nonlinear iterations parameters
     count = 1
     residual_nrm2 = 1.0_rp

     ! Nonlinear iterations
     do  while (residual_nrm2 > tolerance .and. count <= max_number_iterations)
        
        call fe_affine_operator%symbolic_setup()
        call fe_affine_operator%numerical_setup()

        ! Extract the matrix and the RHS
        matrix => fe_affine_operator%get_matrix()
        rhs    => fe_affine_operator%get_array() 
        select type (matrix)
           class is (sparse_matrix_t) 
           sparse_matrix => matrix
           class DEFAULT
           assert(.false.)
        end select

        call dG_CDR_integration%compute_graph_laplacian_perturbed_matrix(sparse_matrix)

        ! DIRECT SOLVER ===============================================
        call direct_solver%set_type_from_pl(direct_solver_parameters)
        call direct_solver%set_matrix(sparse_matrix)
        call direct_solver%set_parameters_from_pl(direct_solver_parameters)
        call direct_solver%update_matrix(sparse_matrix, same_nonzero_pattern=.true.) 

        ! Nonlinear iteration Direct solver     
        select type (rhs) 
        class is (serial_scalar_array_t)
           select type (dof_values) 
           class is (serial_scalar_array_t)
              !call rhs%print(6)
              call direct_solver%solve(rhs , dof_values)
           end select
        class DEFAULT
           assert(.false.)
        end select

        call direct_solver%log_info()
        call direct_solver%free() 

        ! Evaluate norm of the residual
        call fe_affine_operator%apply(dof_values,residual)
        residual_nrm2 = residual%nrm2()

        ! Print current norm of the residual
        write(*,*) 'Norm of the residual at iteration',count,'=',residual_nrm2                                              
        call fe_affine_operator%free_in_stages(free_numerical_setup)

        ! Update iteration counter
        count = count + 1
     end do
     
     select type (dof_values)
        class is ( serial_scalar_array_t) 
        call dof_values%print(6)
        class default
        check(.false.)
     end select
     
     ! current: u^{n+theta} <- u^{n} + update current, previous and integration
     call theta_method%update_solutions( dof_values, dof_values_previous )
     call theta_method%update_integration()

     ! Update operator with solution at current iteration
     if (.not. theta_method%finished ) then 
!!$       
!!$        select type (matrix)
!!$           class is (sparse_matrix_t)  
!!$           call direct_solver%update_matrix(matrix, same_nonzero_pattern=mod(iter,2)==0)
!!$           class DEFAULT
!!$           assert(.false.) 
!!$        end select
     end  if

  end do temporal
  call  vtk_handler%initialize(fe_space, senv, dir_path_out, prefix)
  err = vtk_handler%open_vtu() 
  err = vtk_handler%write_vtu_mesh()
  err = vtk_handler%write_vtu_node_field(fe_function, 1, 'Unkno')
  err = vtk_handler%close_vtu()

  !call fe_space%print()
  call fe_values_previous%free()
  call dG_CDR_integration%analytical_functions%free()
  call fe_affine_operator%free()
  call fe_space%free()
  call reference_fe_array_one(1)%free()
  call residual%free()
  call dof_values%free()
  call triangulation_free ( f_trian )
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)
  call memstatus
contains
 
end program test_cdr
