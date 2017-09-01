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
module vector_poisson_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: vector_poisson_discrete_integration_t
     private
     class(vector_function_t), pointer :: source_term => NULL()
     type(fe_function_t)     , pointer :: fe_function => NULL()
   contains
     procedure :: set_source_term
     procedure :: set_fe_function
     procedure :: integrate_galerkin
  end type vector_poisson_discrete_integration_t
  
  public :: vector_poisson_discrete_integration_t
  
contains
   
  subroutine set_source_term (this, vector_function)
    implicit none
    class(vector_poisson_discrete_integration_t)        , intent(inout) :: this
    class(vector_function_t)                    , target, intent(in)    :: vector_function
    this%source_term => vector_function
  end subroutine set_source_term

  subroutine set_fe_function (this, fe_function)
     implicit none
     class(vector_poisson_discrete_integration_t), intent(inout) :: this
     type(fe_function_t)                          , target, intent(in)    :: fe_function
     this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(vector_poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                    , intent(inout) :: fe_space
    class(assembler_t)             , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(vector_field_t), allocatable  :: shape_values(:,:)
    type(tensor_field_t), allocatable  :: shape_gradients(:,:)

    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor, source_term_value
    
    type(vector_field_t) :: source_term

    assert (associated(this%fe_function)) 

    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )

    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished())
       
       ! Update FE-integration related data structures
       call fe%update_integration()
       
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       num_dofs = fe%get_num_dofs()

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe%get_quadrature_points_coordinates()
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call fe%get_gradients(shape_gradients)
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
       
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
          ! Diffusive term
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * double_contract(shape_gradients(jdof,qpoint),shape_gradients(idof,qpoint))
             end do
          end do
          
          ! Source term
          call this%source_term%get_value_space(quad_coords(qpoint),source_term)
          do idof = 1, num_dofs
             elvec(idof) = elvec(idof) + factor * source_term * shape_values(idof,qpoint)
          end do 
          
       end do
       
       call fe%assembly( this%fe_function, elmat, elvec, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    deallocate(shape_values, stat=istat); check(istat==0);
    deallocate(shape_gradients, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module vector_poisson_discrete_integration_names
