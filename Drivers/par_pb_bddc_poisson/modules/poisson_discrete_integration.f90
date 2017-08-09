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
module pb_bddc_poisson_cG_discrete_integration_names
  use fempar_names
  use pb_bddc_poisson_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_cG_discrete_integration_t
     type(poisson_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(fe_function_t)                 , pointer :: fe_function          => NULL()
     real(rp), public :: diffusion_inclusion
   contains
     procedure :: set_analytical_functions
     procedure :: set_fe_function
     procedure :: integrate_galerkin
  end type poisson_cG_discrete_integration_t
  
  public :: poisson_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(poisson_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(poisson_cG_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)             , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine integrate_galerkin ( this, fe_space, matrix_array_assembler )
    implicit none
    class(poisson_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    class(fe_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(cell_map_t)           , pointer :: cell_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(cell_integrator_t), pointer :: cell_int
    type(vector_field_t)               :: grad_test, grad_trial
    real(rp)                           :: shape_trial


    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: factor
    real(rp)     :: source_term_value

    class(scalar_function_t), pointer :: source_term

    real(rp) :: viscosity
    
    assert (associated(this%analytical_functions))
    assert (associated(this%fe_function))

    source_term => this%analytical_functions%get_source_term()

    call fe_space%initialize_fe_integration()
    call fe_space%create_fe_iterator(fe)
    
    num_dofs = fe%get_num_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_num_quadrature_points()
    cell_map          => fe%get_cell_map()
    cell_int         => fe%get_cell_integrator(1)
    do while ( .not. fe%has_finished())
       if ( fe%is_local() ) then
          ! Update FE-integration related data structures
          call fe%update_integration()
          
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => cell_map%get_quadrature_points_coordinates()
          
          ! Get subset_id
          if ( fe%get_set_id() <= 1 ) then
             viscosity = 1.0_rp
          else 
             viscosity = this%diffusion_inclusion
          end if
          
          !if (viscosity == 0.0_rp) viscosity = 1.0_rp
          
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          do qpoint = 1, num_quad_points
             factor = cell_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
             do idof = 1, num_dofs
                call cell_int%get_gradient(idof, qpoint, grad_trial)
                do jdof = 1, num_dofs
                   call cell_int%get_gradient(jdof, qpoint, grad_test)
                   ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * grad_test * grad_trial * viscosity
                end do
             end do
             
             ! Source term
             call source_term%get_value(quad_coords(qpoint),source_term_value)
             do idof = 1, num_dofs
                call cell_int%get_value(idof, qpoint, shape_trial)
                elvec(idof) = elvec(idof) + factor * source_term_value * shape_trial !* viscosity
             end do
          end do
          
          call fe%assembly( this%fe_function, elmat, elvec, matrix_array_assembler )
       end if
       call fe%next()
    end do
    call fe_space%free_fe_iterator(fe)
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module pb_bddc_poisson_cG_discrete_integration_names
