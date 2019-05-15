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

!****************************************************************************************************

module stokes_discrete_integration_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: stokes_discrete_integration_t
     class(vector_function_t), pointer :: source_term          => NULL()
     type(fe_function_t), pointer :: fe_function          => NULL()
     real(rp)                     :: viscosity 
   contains
     procedure :: integrate_galerkin   => stokes_discrete_integration_integrate
     procedure :: integrate_residual   => stokes_discrete_integration_integrate_residual  
     procedure :: integrate_tangent    => stokes_discrete_integration_integrate_tangent      
     procedure :: set_source_term
     procedure :: set_viscosity
     procedure :: set_fe_function 
  end type stokes_discrete_integration_t
  
  public :: stokes_discrete_integration_t
  
contains
  subroutine set_source_term ( this, source_term )
     implicit none
     class(stokes_discrete_integration_t), intent(inout) :: this
     class(vector_function_t), target     , intent(in)    :: source_term
     this%source_term => source_term
  end subroutine set_source_term  
  
  subroutine set_viscosity ( this, viscosity )
     implicit none
     class(stokes_discrete_integration_t), intent(inout) :: this
    real(rp)                             , intent(in)    :: viscosity
     this%viscosity = viscosity
  end subroutine set_viscosity

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(stokes_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

subroutine stokes_discrete_integration_integrate ( this, fe_space, assembler )
  implicit none
  class(stokes_discrete_integration_t), intent(in)    :: this
  class(serial_fe_space_t)               , intent(inout) :: fe_space
  class(assembler_t)        , intent(inout) :: assembler
  call this%integrate_tangent(fe_space, assembler)
  call this%integrate_residual(fe_space, assembler)
end subroutine stokes_discrete_integration_integrate

subroutine stokes_discrete_integration_integrate_tangent ( this, fe_space, assembler )
  implicit none
  class(stokes_discrete_integration_t), intent(in)    :: this
  class(serial_fe_space_t)               , intent(inout) :: fe_space
  class(assembler_t)        , intent(inout) :: assembler
  ! FE space traversal-related data types
  class(fe_cell_iterator_t), allocatable :: fe
  ! FE integration-related data types
  type(quadrature_t)       , pointer :: quad
  type(point_t)            , pointer :: quad_coords(:)
  type(vector_field_t), allocatable  :: shape_p_gradients(:,:)
  real(rp)            , allocatable  :: shape_p_values(:,:)
  type(tensor_field_t), allocatable  :: shape_u_gradients(:,:)
  type(vector_field_t), allocatable  :: shape_u_values(:,:)
  ! Workspace (FE matrix and vector, assembly data), it could be allocated in the creation
  real(rp)   , allocatable :: elmat(:,:), elvec(:)
  integer(ip), allocatable :: num_dofs_per_field(:)  

  ! Locals
  integer(ip)  :: istat
  integer(ip)  :: qpoint, num_quad_points
  integer(ip)  :: idof, jdof, idof_u, jdof_u, idof_p, jdof_p , num_dofs

  type(fe_cell_function_vector_t) :: cell_solution_u
  type(fe_cell_function_scalar_t) :: cell_solution_p
  type(vector_field_t), pointer   :: solution_u(:) => null()
  type(tensor_field_t), pointer   :: solution_gradu(:) => null()
  real(rp)            , pointer   :: solution_p(:) => null()

  ! Problem variables
  type(vector_field_t) :: source_term_value
  type(tensor_field_t) :: s_u, epsd_v, epsd_u
  real(rp)     :: dV, div_v, div_u
  
  assert (associated(this%source_term)) 

  call fe_space%set_up_cell_integration()
  call fe_space%create_fe_cell_iterator(fe)
  !call cell_solution_u%create(fe_space,1)
  !call cell_solution_p%create(fe_space,2)


  num_dofs = fe%get_num_dofs()
  call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
  call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
  call memalloc ( fe_space%get_num_fields(), num_dofs_per_field, __FILE__, __LINE__ )
  num_dofs_per_field(1) = fe%get_num_dofs_field(1)
  num_dofs_per_field(2) = fe%get_num_dofs_field(2)

  quad            => fe%get_quadrature()
  num_quad_points = quad%get_num_quadrature_points()

  do while ( .not. fe%has_finished())
     if ( fe%is_local() ) then
        ! Update FE-integration related data structures
        call fe%update_integration()
        !call cell_solution_u%update(fe,this%fe_function)
        !call cell_solution_p%update(fe,this%fe_function)
        !solution_u      => cell_solution_u%get_quadrature_points_values()
        !solution_gradu => cell_solution_u%get_quadrature_points_gradients()
        !solution_p     => cell_solution_p%get_quadrature_points_values()

        ! Get quadrature coordinates to evaluate source_term
        quad_coords => fe%get_quadrature_points_coordinates()

        ! Compute element matrix and vector
        elmat = 0.0_rp

        call fe%get_gradients(shape_u_gradients,1)
        call fe%get_values(shape_u_values,1)
        call fe%get_gradients(shape_p_gradients,2)
        call fe%get_values(shape_p_values,2)

        do qpoint = 1, num_quad_points
          dV = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof_u = 1, num_dofs_per_field(1)
            idof = idof_u
            ! U-U
            epsd_v  = symmetric_part(shape_u_gradients(idof_u,qpoint)) 
            div_v = trace(epsd_v)
            do jdof_u = 1, num_dofs_per_field(1)
              jdof = jdof_u
              epsd_u  = symmetric_part(shape_u_gradients(jdof_u,qpoint))
              !s_u = 2*viscosity*epsd_u
              elmat(idof,jdof) = elmat(idof,jdof) + dV * double_contract(epsd_v,2*this%viscosity*epsd_u)
              !> Convective term $\int_\Omega u_i \cdot \partial_i u_j v_j
              !elmat(idof,jdof) = elmat(idof,jdof) + dv * (solution_u(qpoint) * shape_u_gradients(jdof_u,qpoint)) * shape_u_values(idof_u,qpoint)
            end do
            ! U-P
            do jdof_p = 1, num_dofs_per_field(2)
              jdof = num_dofs_per_field(1)+jdof_p
              elmat(idof,jdof) = elmat(idof,jdof) - dV * div_v * shape_p_values(jdof_p,qpoint)
            end do
         end do
         do idof_p = 1, num_dofs_per_field(2)
           idof = num_dofs_per_field(1)+idof_p
           ! P-U
           do jdof_u = 1, num_dofs_per_field(1)
             div_u = trace(shape_u_gradients(jdof_u,qpoint))
             jdof = jdof_u
             elmat(idof,jdof) = elmat(idof,jdof) + dV * shape_p_values(idof_p,qpoint) * div_u
           end do
         end do
       end do      

       call fe%assembly( elmat, assembler)
     end if
     call fe%next()
  end do

  !call cell_solution_u%free()
  !call cell_solution_p%free()
  call fe_space%free_fe_cell_iterator(fe)
  call memfree(shape_p_values, __FILE__, __LINE__)
  deallocate (shape_p_gradients, stat=istat); check(istat==0);
  deallocate (shape_u_values, stat=istat); check(istat==0);
  deallocate (shape_u_gradients, stat=istat); check(istat==0);
  call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
  call memfree ( elmat, __FILE__, __LINE__ )
  call memfree ( elvec, __FILE__, __LINE__ )

end subroutine stokes_discrete_integration_integrate_tangent



subroutine stokes_discrete_integration_integrate_residual ( this, fe_space, assembler )
  implicit none
  class(stokes_discrete_integration_t), intent(in)    :: this
  class(serial_fe_space_t)               , intent(inout) :: fe_space
  class(assembler_t)        , intent(inout) :: assembler

  ! FE space traversal-related data types
  class(fe_cell_iterator_t), allocatable :: fe

  ! FE integration-related data types
  type(quadrature_t)       , pointer :: quad
  type(point_t)            , pointer :: quad_coords(:)
  type(vector_field_t), allocatable  :: shape_p_gradients(:,:)
  real(rp)            , allocatable  :: shape_p_values(:,:)
  type(tensor_field_t), allocatable  :: shape_u_gradients(:,:)
  type(vector_field_t), allocatable  :: shape_u_values(:,:)

  ! Workspace (FE matrix and vector, assembly data), it could be allocated in the creation
  real(rp)   , allocatable :: elmat(:,:), elvec(:)
  integer(ip), allocatable :: num_dofs_per_field(:)  

  integer(ip)  :: istat
  integer(ip)  :: qpoint, num_quad_points
  integer(ip)  :: idof, jdof, idof_u, jdof_u, idof_p, jdof_p , num_dofs

  type(fe_cell_function_vector_t) :: cell_solution_u
  type(fe_cell_function_scalar_t) :: cell_solution_p
  type(vector_field_t), pointer   :: solution_u(:) => null()
  type(tensor_field_t), pointer   :: solution_gradu(:) => null()
  real(rp)            , pointer   :: solution_p(:) => null()

  ! Problem variables
  type(vector_field_t) :: source_term_value
  type(tensor_field_t) :: s_u, epsd_v, epsd_u
  real(rp)     :: dV, div_v, div_u
  
  assert (associated(this%source_term)) 

  call fe_space%set_up_cell_integration()
  call fe_space%create_fe_cell_iterator(fe)
  call cell_solution_u%create(fe_space,1)
  call cell_solution_p%create(fe_space,2)

  num_dofs = fe%get_num_dofs()
  call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
  call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
  call memalloc ( fe_space%get_num_fields(), num_dofs_per_field, __FILE__, __LINE__ )
  num_dofs_per_field(1) = fe%get_num_dofs_field(1)
  num_dofs_per_field(2) = fe%get_num_dofs_field(2)

  quad            => fe%get_quadrature()
  num_quad_points = quad%get_num_quadrature_points()

  do while ( .not. fe%has_finished())
     if ( fe%is_local() ) then
        ! Update FE-integration related data structures
        call fe%update_integration()
        call cell_solution_u%update(fe,this%fe_function)
        call cell_solution_p%update(fe,this%fe_function)
        solution_u      => cell_solution_u%get_quadrature_points_values()
        solution_gradu => cell_solution_u%get_quadrature_points_gradients()
        solution_p     => cell_solution_p%get_quadrature_points_values()

        ! Get quadrature coordinates to evaluate source_term
        quad_coords => fe%get_quadrature_points_coordinates()

        ! Compute element matrix and vector
        elvec = 0.0_rp

        call fe%get_gradients(shape_u_gradients,1)
        call fe%get_values(shape_u_values,1)
        call fe%get_gradients(shape_p_gradients,2)
        call fe%get_values(shape_p_values,2)

        do qpoint = 1, num_quad_points
          dV = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
         ! Residual
         call this%source_term%get_value_space(quad_coords(qpoint),source_term_value)
         epsd_u = symmetric_part(solution_gradu(qpoint))
         div_u  = trace(epsd_u)
         do idof_u = 1, num_dofs_per_field(1)
           idof = idof_u
           epsd_v  = symmetric_part(shape_u_gradients(idof_u,qpoint)) 
           div_v = trace(epsd_v)
           elvec(idof) = elvec(idof) - dV * source_term_value * shape_u_values(idof_u,qpoint)
           elvec(idof) = elvec(idof) + dV * double_contract(epsd_v,2*this%viscosity*epsd_u)
           elvec(idof) = elvec(idof) - dV * div_v * solution_p(qpoint)
           !> Convective term $\int_\Omega u_i \cdot \partial_i u_j v_j
           !elvec(idof) = elvec(idof) + dV * (solution_u(qpoint)*solution_gradu(qpoint))*shape_u_values(idof_u,qpoint)
         end do
 
         do idof_p = 1, num_dofs_per_field(2)
           idof = num_dofs_per_field(1)+idof_p
           elvec(idof) = elvec(idof) + dV * shape_p_values(idof_p,qpoint) * div_u
         end do
       end do      

       call fe%assembly( elvec, assembler)
     end if
     call fe%next()
  end do

  call cell_solution_u%free()
  call cell_solution_p%free()
  call fe_space%free_fe_cell_iterator(fe)
  call memfree(shape_p_values, __FILE__, __LINE__)
  deallocate (shape_p_gradients, stat=istat); check(istat==0);
  deallocate (shape_u_values, stat=istat); check(istat==0);
  deallocate (shape_u_gradients, stat=istat); check(istat==0);
  call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
  call memfree ( elmat, __FILE__, __LINE__ )
  call memfree ( elvec, __FILE__, __LINE__ )

end subroutine stokes_discrete_integration_integrate_residual
end module stokes_discrete_integration_names

