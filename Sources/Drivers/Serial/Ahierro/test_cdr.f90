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

  type test_cdr_params_t
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out

     character(len=:), allocatable :: default_kfl_conv 
     character(len=:), allocatable :: default_kfl_tder 
     character(len=:), allocatable :: default_kfl_react 
     character(len=:), allocatable :: default_react
     character(len=:), allocatable :: default_diffu 
     character(len=:), allocatable :: default_space_solution_flag 
     character(len=:), allocatable :: default_tempo_solution_flag 

     character(len=:), allocatable :: default_kfl_stab
     character(len=:), allocatable :: default_ftime
     character(len=:), allocatable :: default_itime
     character(len=:), allocatable :: default_theta
     character(len=:), allocatable :: default_time_step

     character(len=:), allocatable :: default_continuity
     character(len=:), allocatable :: default_enable_face_integration
     character(len=:), allocatable :: default_order

  end type test_cdr_params_t

  ! Types
  public :: test_cdr_params_t

  ! Functions
  public :: set_default_params,cli_add_params,set_default_params_analytical
  public :: set_default_params_transient

contains

  subroutine set_default_params(params)
    implicit none
    type(test_cdr_params_t), intent(inout) :: params

    ! IO parameters
    params%default_dir_path     = 'data'
    params%default_prefix       = 'square_4x4'
    params%default_dir_path_out = 'output'

    ! Problem parameters
    params%default_kfl_conv            = '0' ! Enabling advection
    params%default_kfl_tder            = '0' ! Time derivative not computed 
    params%default_kfl_react           = '0' ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '3'
    params%default_tempo_solution_flag = '0'

    ! Solver parameter
    params%default_kfl_stab = '0'   ! Stabilization of convective term (0: Off, 2: OSS)

    ! Time integration
    params%default_itime           = '0'
    params%default_ftime           = '0'
    params%default_theta           = '1.0'
    params%default_time_step       = '0.0'

    ! FE Space parameters
    params%default_continuity              = '1'
    params%default_enable_face_integration = '.false.'
    params%default_order                   = '1'
  end subroutine set_default_params
  !==================================================================================================

  subroutine cli_add_params(cli,params,group)
    implicit none
    type(Type_Command_Line_Interface)            , intent(inout) :: cli
    type(test_cdr_params_t)                      , intent(in)    :: params
    character(*)                                 , intent(in)    :: group
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
    call cli%add(group=trim(group),switch='--kfl_conv',switch_ab='-kconv',help='Convection flag',   &
         &       required=.false.,act='store',def=trim(params%default_kfl_conv),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_react',switch_ab='-kreac',help='Reaction flag',    &
         &       required=.false.,act='store',def=trim(params%default_kfl_react),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--kfl_tder',switch_ab='-ktd',                            &
         &       help='Temporal derivative computation flag',required=.false.,act='store',          &
         &       def=trim(params%default_kfl_tder),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--react',switch_ab='-reac',help='Constant Reaction Value'&
         &       ,required=.false.,act='store',def=trim(params%default_react),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--diffu',switch_ab='-diff',                              &
         &       help='Constant Diffusion Value',required=.false.,act='store',                      &
         &       def=trim(params%default_diffu),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--space_solution_flag',switch_ab='-ssol',                &
         &       help='Space analytical solution',required=.false.,act='store',nargs='*',           &
         &       def=trim(params%default_space_solution_flag),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--tempo_solution_flag',switch_ab='-tsol',                &
         &       help='Temporal analytical solution',required=.false.,act='store',nargs='*',        &
         &       def=trim(params%default_tempo_solution_flag),error=error)
    check(error==0)

    ! Solver parameters
    call cli%add(group=trim(group),switch='--kfl_stab',switch_ab='-kst',help='Stabilization flag',  &
         &       required=.false.,act='store',def=trim(params%default_kfl_stab),error=error)
    check(error==0)

    ! Time integration parameters
    call cli%add(group=trim(group),switch='--itime',switch_ab='-t0',help='Initial time',            &
         &       required=.false.,act='store',def=trim(params%default_itime),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--ftime',switch_ab='-tf',help='Final time',              &
         &       required=.false.,act='store',def=trim(params%default_ftime),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--tstep',switch_ab='-ts',help='Time step',               &
         &       required=.false.,act='store',def=trim(params%default_time_step),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--theta',switch_ab='-tht',help='Theta method',           &
         &       required=.false.,act='store',def=trim(params%default_theta),error=error)
    check(error==0)

    ! FE Space parameters 
    call cli%add(group=trim(group),switch='--continuity',switch_ab='-cg',                           &
         &       help='Flag for the continuity of the FE Space',required=.false.,act='store',       &
         &       def=trim(params%default_continuity),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--face_integ',switch_ab='-fi',                           &
         &       help='Allow face integration',required=.false.,act='store',                        &
         &       def=trim(params%default_enable_face_integration),error=error)
    check(error==0)
    call cli%add(group=trim(group),switch='--order',switch_ab='-p',                                 &
         &       help='Initial order of the approximation',required=.false.,act='store',            &
         &       def=trim(params%default_order),error=error)
    check(error==0)
  end subroutine cli_add_params
  !==================================================================================================

  subroutine set_default_params_analytical(params)
    implicit none
    type(test_cdr_params_t), intent(inout) :: params

    ! Names
    params%default_kfl_conv            = '0'    ! Enabling advection
    params%default_kfl_tder            = '0'    ! Time derivative not computed 
    params%default_kfl_react           = '0'    ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '4'
    params%default_tempo_solution_flag = '0'

  end subroutine set_default_params_analytical
  !==================================================================================================

  subroutine set_default_params_transient(params)
    implicit none
    type(test_cdr_params_t), intent(inout) :: params

    ! Problem
    params%default_kfl_conv            = '11'   ! Enabling advection
    params%default_kfl_tder            = '0'    ! Time derivative not computed 
    params%default_kfl_react           = '0'    ! Non analytical reaction
    params%default_react               = '0.0'  ! Reaction
    params%default_diffu               = '1.0'  ! Diffusion
    params%default_space_solution_flag = '3'
    params%default_tempo_solution_flag = '1'
    params%default_itime               = '0'
    params%default_ftime               = '1'

    ! Time integration method
    params%default_theta           = '1.0'
    params%default_time_step       = '0.1'


  end subroutine set_default_params_transient

end module command_line_parameters_names

!****************************************************************************************************
module CDR_discrete_integration_names
 use serial_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: CDR_discrete_integration_t
     integer(ip) :: viscosity 
   contains
     procedure :: integrate
  end type CDR_discrete_integration_t
  
  public :: CDR_discrete_integration_t
  
contains
  
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(CDR_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)      , intent(inout) :: fe_space
    class(assembler_t)            , intent(inout) :: assembler

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
       !write (*,*) 'XXXXXXXXX ELMAT XXXXXXXXX'
       !write (*,*) elmat
       
       ! Apply boundary conditions
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value, number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks,  field_coupling, elmat, elvec )
    end do
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
end module CDR_discrete_integration_names

!==================================================================================================
module vector_dG_CDR_discrete_integration_names
use serial_names

implicit none
# include "debug.i90"
private
type, extends(discrete_integration_t) :: vector_dG_CDR_discrete_integration_t
real(rp) :: viscosity 
real(rp) :: C_IP ! Interior Penalty constant
real(rp) :: xi   ! DG symmetric parameter: (0-Symmetric, 1/2-Incomplete, 1-Antisymmetric)
contains
procedure :: set_problem
procedure :: integrate
end type vector_dG_CDR_discrete_integration_t

public :: vector_dG_CDR_discrete_integration_t

contains
  subroutine set_problem(this, viscosity, C_IP, xi)
    implicit none
    class(vector_dG_CDR_discrete_integration_t), intent(inout) :: this
    real(rp)                                   , intent(in)    :: viscosity 
    real(rp)                                   , intent(in)    :: C_IP 
    real(rp)                                   , intent(in)    :: xi
    
    this%viscosity = viscosity
    this%C_IP      = C_IP
    this%xi        = xi
  end subroutine set_problem
  subroutine integrate ( this, fe_space, assembler )
    implicit none
    class(vector_dG_CDR_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)                      , intent(inout) :: assembler

    type(finite_element_t)   , pointer :: fe
    type(finite_face_t)         , pointer :: face
    type(volume_integrator_t), pointer :: vol_int_first_fe, vol_int_second_fe
    type(face_integrator_t)     , pointer :: face_int
    real(rp)                , allocatable :: elmat(:,:), elvec(:), facemat(:,:,:,:), facevec(:,:)
    type(fe_map_t)              , pointer :: fe_map
    type(face_map_t)            , pointer :: face_map
    type(quadrature_t)       , pointer :: quad
    integer(ip)             , allocatable :: number_nodes_per_field(:)

    integer(ip)  :: igaus,inode,jnode,ioffset,joffset,ngaus
    real(rp)     :: factor, h_length, bcvalue

    real(rp)             :: shape_test_scalar, shape_trial_scalar
    type(vector_field_t) :: grad_test_scalar, grad_trial_scalar
    type(tensor_field_t) :: grad_test_vector, grad_trial_vector
    type(vector_field_t) :: normal(2)
    
    integer(ip)          :: i, j, number_fe_spaces

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    integer(ip)          :: ielem, iface, number_nodes, ineigh, jneigh, number_neighbours
    type(i1p_t), pointer :: elem2dof(:),test_elem2dof(:),trial_elem2dof(:)
    type(i1p_t), pointer :: bc_code(:)
    type(r1p_t), pointer :: bc_value(:)

    number_fe_spaces = fe_space%get_number_fe_spaces()
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    fe => fe_space%get_finite_element(1)
    number_nodes = fe%get_number_nodes()

    call memalloc ( number_fe_spaces, number_nodes_per_field, __FILE__, __LINE__ )
    call fe%get_number_nodes_per_field( number_nodes_per_field )

    ! --------------------------- LOOP OVER THE ELEMENTS -------------------------
    call memalloc ( number_nodes, number_nodes, elmat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, elvec, __FILE__, __LINE__ )
    
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
                elmat(inode,jnode) = elmat(inode,jnode) + factor *                                  &
                     &               this%viscosity * grad_test_scalar*grad_trial_scalar
             end do
          end do

          do inode = 1, number_nodes_per_field(2)
             ioffset = number_nodes_per_field(1)+inode
             call vol_int_second_fe%get_gradient(inode,igaus,grad_trial_vector)
             do jnode = 1, number_nodes_per_field(2)
                joffset = number_nodes_per_field(1)+jnode
                call vol_int_second_fe%get_gradient(jnode,igaus,grad_test_vector)
                elmat(ioffset,joffset) = elmat(ioffset,joffset) + factor *                        &
                     &          this%viscosity*double_contract(grad_test_vector,grad_trial_vector)
             end do
          end do
       end do
       
       
       call this%impose_strong_dirichlet_data( elmat, elvec, bc_code, bc_value,                   &
            &                                  number_nodes_per_field, number_fe_spaces )
       call assembler%assembly( number_fe_spaces, number_nodes_per_field, elem2dof, field_blocks, &
            &                   field_coupling, elmat, elvec )      
    end do
    do inode = 1, number_nodes_per_field(1)
       write(*,*) inode, '+++++'
       write(*,*) elmat(1:4,inode)
    end do
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )


    ! --------------------------- LOOP OVER THE FACES -----------------------------
    call memalloc ( number_nodes, number_nodes,2,2,facemat, __FILE__, __LINE__ )
    call memalloc ( number_nodes, 2,facevec, __FILE__, __LINE__ )

    !call fe_space%initialize_face_integration()
    write(*,*) __FILE__,__LINE__,fe_space%get_number_interior_faces()
    do iface = 1, fe_space%get_number_interior_faces()

       facemat = 0.0_rp
       facevec = 0.0_rp

       face => fe_space%get_finite_face(iface)
       number_neighbours = face%number_neighbours()
       write(*,*) __FILE__,__LINE__,iface,'------------------'
       call face%update_integration()
      
       quad   => face%get_quadrature()
       ngaus = quad%get_number_quadrature_points()
       face_map => face%get_map()

       j = 1
       face_int => face%get_face_integrator(j)
       
       do igaus = 1, ngaus
          call face_map%get_normals(igaus,normal)
          h_length = face_map%compute_characteristic_length(igaus,number_neighbours)
          factor = face_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do ineigh = 1, number_neighbours
             do inode = 1, number_nodes_per_field(j)
                !ioffset = number_nodes_per_field(j)*(ineigh-1) + inode
                call face_int%get_value(inode,igaus,ineigh,shape_trial_scalar)
                call face_int%get_gradient(inode,igaus,ineigh,grad_trial_scalar)
                do jneigh = 1, number_neighbours
                   do jnode = 1, number_nodes_per_field(j)
                      !joffset = number_nodes_per_field(j)*(jneigh-1) + jnode
                      call face_int%get_value(jnode,igaus,jneigh,shape_test_scalar)
                      call face_int%get_gradient(jnode,igaus,jneigh,grad_test_scalar)
                      !- mu*({{grad u}}[[v]] + (1-xi)*[[u]]{{grad v}} ) + C*mu*p^2/h * [[u]] [[v]]
                      facemat(inode,jnode,ineigh,jneigh) = facemat(inode,jnode,ineigh,jneigh) +     &
                           &  factor * this%viscosity *   &
                           &  (-0.5_rp*grad_test_scalar*normal(ineigh)*shape_trial_scalar - &
                           &   0.5_rp*grad_trial_scalar*normal(jneigh)*shape_test_scalar  + &
                           &   this%c_IP / h_length * shape_test_scalar*shape_trial_scalar *        &
                           &   normal(ineigh)*normal(jneigh))
                   end do
                end do
             end do
          end do
       end do
       do ineigh = 1, number_neighbours
          test_elem2dof => face%get_elem2dof(ineigh)
          do jneigh = 1, number_neighbours
             trial_elem2dof => face%get_elem2dof(jneigh)
             call assembler%face_assembly(number_fe_spaces,number_nodes_per_field, number_nodes_per_field,&
                  &                  test_elem2dof,trial_elem2dof,field_blocks,field_coupling,       &
                  &                  facemat(:,:,ineigh,jneigh), facevec(:,ineigh) )   
          end do
       end do
    end do

    do inode = 1, number_nodes_per_field(j)
       write(*,*) inode, '-------'
       write(*,*) facemat(1:4,inode,1,1),  facemat(1:4,inode,2,1)
    end do
    do inode = 1, number_nodes_per_field(j)
       write(*,*) number_nodes_per_field(j)+inode, '-------'
       write(*,*) facemat(1:4,inode,1,2),  facemat(1:4,inode,2,2)
    end do

    do iface = fe_space%get_number_interior_faces() + 1, fe_space%get_number_interior_faces() +   &
         &                                               fe_space%get_number_boundary_faces()

       facemat = 0.0_rp
       facevec = 0.0_rp

       face => fe_space%get_finite_face(iface)
       number_neighbours = face%number_neighbours()
       write(*,*) __FILE__,__LINE__,iface,'------------------'
       call face%update_integration()

       face_map => face%get_map()
       quad   => face%get_quadrature()
       ngaus = quad%get_number_quadrature_points()

       j = 1
       face_int => face%get_face_integrator(j)
      
       do igaus = 1, ngaus
          call face_map%get_normals(igaus,normal)
          h_length = face_map%compute_characteristic_length(igaus,number_neighbours)
          factor = face_map%get_det_jacobian(igaus) * quad%get_weight(igaus)
          do ineigh = 1, number_neighbours
             do inode = 1, number_nodes_per_field(j)
                call face_int%get_value(inode,igaus,ineigh,shape_trial_scalar)
                call face_int%get_gradient(inode,igaus,ineigh,grad_trial_scalar)
                do jneigh = 1, number_neighbours
                   do jnode = 1, number_nodes_per_field(j)
                      call face_int%get_value(jnode,igaus,jneigh,shape_test_scalar)
                      call face_int%get_gradient(jnode,igaus,jneigh,grad_test_scalar)
                      facemat(inode,jnode,ineigh,jneigh) = facemat(inode,jnode,ineigh,jneigh) +     &
                           &  factor * this%viscosity *   &
                           &  (-grad_test_scalar*normal(ineigh)*shape_trial_scalar - &
                           &   grad_trial_scalar*normal(jneigh)*shape_test_scalar  + &
                           &   this%c_IP / h_length * shape_test_scalar*shape_trial_scalar)
                   end do
                end do
                bcvalue = 1.0_rp
                facevec(inode,ineigh) = facevec(inode,ineigh) + factor * this%viscosity *           &
                     &                  (bcvalue * grad_trial_scalar*normal(jneigh) +       &
                     &                  this%c_IP/h_length * bcvalue * shape_trial_scalar )
             end do
          end do
       end do
    end do
 
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )
    ! ----------------------------------------------------------------------------
    call memfree ( number_nodes_per_field, __FILE__, __LINE__ )
  end subroutine integrate

end module vector_dG_CDR_discrete_integration_names


!****************************************************************************************************
program test_cdr
  use serial_names
  use Data_Type_Command_Line_Interface
  use command_line_parameters_names
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
  integer(ip)              :: space_solution_flag(1), tempo_solution_flag(1)

  integer(ip), allocatable :: order(:,:), material(:), problem(:)

  integer(ip), allocatable :: continuity(:,:)
  logical    , allocatable :: enable_face_integration(:,:)

  logical                  :: diagonal_blocks_symmetric_storage(1)
  logical                  :: diagonal_blocks_symmetric(1)
  integer(ip)              :: diagonal_blocks_sign(1)

  integer(ip)              :: lunio, istat

  class(vector_t) , pointer :: x, y
  class(operator_t), pointer :: A

  call meminit

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

  ! Read mesh
  call mesh_read (dir_path, prefix, f_mesh, permute_c2z=.true.)

  ! Read conditions 
  call conditions_read (dir_path, prefix, f_mesh%npoin, f_cond)
  ! Construc triangulation
  call mesh_to_triangulation ( f_mesh, f_trian, gcond = f_cond )

  call triangulation_construct_faces ( f_trian )


  ! To be erased
  call test_reference_face_stuff(f_trian,f_cond)

  call triangulation_free ( f_trian )
  call conditions_free ( f_cond )
  call mesh_free (f_mesh)
  call memstatus
contains
  !==================================================================================================
  subroutine  test_reference_face_stuff(f_trian, f_cond)
    use serial_names
    use CDR_discrete_integration_names
    use vector_dG_CDR_discrete_integration_names
    use block_sparse_matrix_names

    implicit none

    type(triangulation_t), intent(inout) :: f_trian
    type(conditions_t)   , intent(in)    :: f_cond

    type(serial_fe_space_t)                       :: fe_space
    type(p_reference_fe_t)                        :: reference_fe_array_two(2)
    type(p_reference_fe_t)                        :: reference_fe_array_one(1)
    type(fe_affine_operator_t)                    :: fe_affine_operator
    type(vector_dG_CDR_discrete_integration_t)    :: vector_dG_CDR_integration
    type(CDR_discrete_integration_t)              :: CDR_integration
    class(vector_t)         , allocatable, target :: vector
    type(interpolation_face_restriction_t)                    :: face_interpolation

    class(matrix_t), pointer :: matrix
    type(sparse_matrix_t), pointer :: my_matrix
    logical                  :: diagonal_blocks_symmetric_storage(2)
    logical                  :: diagonal_blocks_symmetric(2)
    integer(ip)              :: diagonal_blocks_sign(2)
    
    ! Composite case
    reference_fe_array_two(1) = make_reference_fe ( topology = topology_quad, fe_type = fe_type_lagrangian,      &
         &                      number_dimensions = 2, order = 1, field_type = field_type_scalar,            &
         &                      continuity = .false. )

    reference_fe_array_two(2) = make_reference_fe ( topology = topology_quad, fe_type = fe_type_lagrangian,      &
         &                      number_dimensions = 2, order = 1, field_type = field_type_vector,            &
         &                      continuity = .false. )

    call fe_space%create( triangulation = f_trian, boundary_conditions = f_cond,                    &
         &                reference_fe_phy = reference_fe_array_two,                                &
         &                field_blocks = (/1,2/),                                                   &
         &                field_coupling = reshape((/.true.,.false.,.false.,.true./),(/2,2/)) )

    call fe_space%create_face_array()
    call fe_space%fill_dof_info() 

     call vector_dG_CDR_integration%set_problem( viscosity = 1.0_rp, C_IP = 10.0_rp, xi = 0.0_Rp)
    ! Create the operator
    diagonal_blocks_symmetric_storage = .false.
    diagonal_blocks_symmetric         = .false.
    diagonal_blocks_sign              = SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE
    call fe_affine_operator%create ('CSR',diagonal_blocks_symmetric_storage ,                       &
         &                          diagonal_blocks_symmetric,diagonal_blocks_sign,         &
         &                          senv, fe_space, vector_dG_CDR_integration)
    call fe_affine_operator%symbolic_setup()
    call fe_affine_operator%numerical_setup()
    matrix => fe_affine_operator%get_matrix()
    select type(matrix)
     class is(block_sparse_matrix_t)
        !my_matrix => matrix
        do i = 1,1! matrix%nblocks
           write(*,*) i,i,'+++++++++++++++++++++++++++++++'
           write(*,*) __FILE__,__LINE__
           my_matrix => matrix%blocks(i,i)%sparse_matrix
           write(*,*) __FILE__,__LINE__
           !call my_matrix%print_matrix_market(6)
           write(*,*) __FILE__,__LINE__
        end do
     class default
        check(.false.)
     end select

     !call fe_affine_operator%create_range_vector(vector)

    !call fe_space%print()
    !call reference_fe_array_two(1)%free
    !call reference_fe_array_two(2)%free
    call fe_affine_operator%free()
    call fe_space%free()
    call reference_fe_array_two(1)%free()
    call reference_fe_array_two(2)%free()
  end subroutine test_reference_face_stuff

  !==================================================================================================
  subroutine read_flap_cli_test_cdr(cli)
    implicit none
    type(Type_Command_Line_Interface), intent(out) :: cli
    ! Locals
    type(test_cdr_params_t) :: analytical_params,transient_params
    logical     :: authors_print
    integer(ip) :: error

    authors_print = .false.

    ! Initialize Command Line Interface
    call cli%init(progname    = 'test_cdr',                                                         &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 & 
         &        license     = '',                                                                 &
         &        description =                                                                     &
         &    'Serial FEMPAR driver to solve transient CDR problems using Continuous-Galerkin .',   &
         &        examples    = ['test_cdr -h            ', 'test_cdr analytical -h ' ])

    ! Set Command Line Arguments Groups, i.e. commands
    call cli%add_group(group='analytical',description='Solve a problem with an analytical solution')
    call cli%add_group(group='transient',description='Solve a problem with an transient solution')

    ! Set Command Line Arguments for each group
    call set_default_params(analytical_params)
    call set_default_params_analytical(analytical_params)
    call cli_add_params(cli,analytical_params,'analytical') 

    call set_default_params(transient_params)
    call set_default_params_transient(transient_params)
    call cli_add_params(cli,transient_params,'transient')

  end subroutine read_flap_cli_test_cdr


end program test_cdr
