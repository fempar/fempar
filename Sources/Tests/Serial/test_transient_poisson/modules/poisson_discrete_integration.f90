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
module poisson_cG_discrete_integration_names
  use fempar_names
  use poisson_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  
  type, extends(linear_discrete_integration_t) :: poisson_cG_discrete_integration_t
     type(poisson_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(fe_function_t)                 , pointer :: fe_function          => NULL()
     real(rp)                                      :: current_time = 0.0_rp
   contains
     procedure :: set_analytical_functions
     procedure :: set_fe_function
     procedure :: set_current_time
     procedure :: integrate_galerkin
     procedure :: integrate_tangent
     procedure :: integrate_residual 
  end type poisson_cG_discrete_integration_t

  
  public :: poisson_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_cG_discrete_integration_t)    , intent(inout) :: this
     type(poisson_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(poisson_cG_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function
  
  subroutine set_current_time (this, current_time)
     implicit none
     class(poisson_cG_discrete_integration_t), intent(inout) :: this
     real(rp)                                , intent(in)    :: current_time
     this%current_time = current_time
  end subroutine set_current_time

  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(poisson_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(assembler_t)      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor
    real(rp)     :: source_term_value

    class(scalar_function_t), pointer :: source_term

    assert (associated(this%analytical_functions))
    assert (associated(this%fe_function)) 
    
    source_term => this%analytical_functions%get_source_term()
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
    
    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )
       
       ! Update FE-integration related data structures
       call fe%update_integration()
          
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()
       
       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe%get_quadrature_points_coordinates()

       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call fe%get_gradients(shape_gradients)
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
             end do
          end do
          
          ! Source term
          call source_term%get_value(quad_coords(qpoint),this%current_time,source_term_value)
          do idof = 1, num_dofs
             elvec(idof) = elvec(idof) + factor * source_term_value * shape_values(idof,qpoint)
          end do 
       end do
       
       call fe%assembly( this%fe_function, elmat, elvec, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)

    call memfree(shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
  subroutine integrate_tangent ( this, fe_space, assembler )
    implicit none
    class(poisson_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)                      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_gradients(:,:)

    ! FE matrix i.e., A_K 
    real(rp), allocatable              :: elmat(:,:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor

    assert (associated(this%analytical_functions))
    assert (associated(this%fe_function)) 
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    
    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )
       ! Update FE-integration related data structures
       call fe%update_integration()
          
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()
       
       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe%get_quadrature_points_coordinates()

       ! Compute element matrix and vector
       elmat = 0.0_rp
       call fe%get_gradients(shape_gradients)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
             end do
          end do
       end do
       call fe%assembly( elmat, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
  end subroutine integrate_tangent
  
  subroutine integrate_residual ( this, fe_space, assembler )
    implicit none
    class(poisson_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)                      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix i.e., A_K 
    real(rp), allocatable              :: elvec(:)
    type(fe_cell_function_scalar_t)    :: u_h
    type(vector_field_t), pointer      :: u_h_grads(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor
    class(scalar_function_t), pointer :: source_term
    real(rp)     :: source_term_value
 
    source_term => this%analytical_functions%get_source_term()
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
    
    call u_h%create(fe_space, field_id=1)
    
    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )
       ! Update FE-integration related data structures
       call fe%update_integration()
       
       call u_h%update(fe, this%fe_function)
       u_h_grads => u_h%get_quadrature_points_gradients()   
          
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()
       
       ! Get quadrature coordinates to evaluate source_term
       quad_coords => fe%get_quadrature_points_coordinates()

       ! Compute element matrix and vector
       elvec = 0.0_rp
       call fe%get_gradients(shape_gradients)
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          call source_term%get_value(quad_coords(qpoint),this%current_time,source_term_value)
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
            ! r_K(i) = (grad(phi_i),grad(u_h))-(f,phi_i)
            elvec(idof) = elvec(idof) + factor * (shape_gradients(idof,qpoint)*u_h_grads(qpoint)- &
                                                 (source_term_value*shape_values(idof,qpoint)))
          end do
       end do
       call fe%assembly( elvec, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    call memfree(shape_values, __FILE__, __LINE__)
    call memfree ( elvec, __FILE__, __LINE__ )
    call u_h%free()
  end subroutine integrate_residual
  
end module poisson_cG_discrete_integration_names
