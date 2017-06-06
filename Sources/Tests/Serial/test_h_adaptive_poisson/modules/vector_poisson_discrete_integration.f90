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
     class(vector_function_t), pointer :: source_term
   contains
     procedure :: set_source_term
     procedure :: integrate
  end type vector_poisson_discrete_integration_t
  
  public :: vector_poisson_discrete_integration_t
  
contains
   
  subroutine set_source_term (this, vector_function)
    implicit none
    class(vector_poisson_discrete_integration_t)        , intent(inout) :: this
    class(vector_function_t)                    , target, intent(in)    :: vector_function
    this%source_term => vector_function
  end subroutine set_source_term

  subroutine integrate ( this, fe_space, matrix_array_assembler )
    implicit none
    class(vector_poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                    , intent(inout) :: fe_space
    class(matrix_array_assembler_t)             , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    class(fe_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int
    type(vector_field_t), allocatable  :: shape_values(:,:)
    type(tensor_field_t), allocatable  :: shape_gradients(:,:)

    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor, source_term_value
    
    type(vector_field_t) :: source_term

    integer(ip)  :: number_fields

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip), allocatable :: num_dofs_per_field(:)  

    
    number_fields = fe_space%get_number_fields()
    allocate( elem2dof(number_fields), stat=istat); check(istat==0);
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()
    
    max_num_dofs = fe_space%get_max_number_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )

    call fe_space%create_fe_iterator(fe)
    do while ( .not. fe%has_finished())
       
       ! Update FE-integration related data structures
       call fe%update_integration()
       
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points = quad%get_number_quadrature_points()
       fe_map          => fe%get_fe_map()
       vol_int         => fe%get_volume_integrator(1)
       num_dofs = fe%get_number_dofs()
       call fe%get_number_dofs_per_field(num_dofs_per_field)
       
       ! Get DoF numbering within current FE
       call fe%get_elem2dof(elem2dof)

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe_map%get_quadrature_coordinates()
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call vol_int%get_gradients(shape_gradients)
       call vol_int%get_values(shape_values)
       do qpoint = 1, num_quad_points
       
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
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
       
       ! Assemble and apply boundary conditions (and the rest of hanging node constraints)
       call fe%assemble( elmat, elvec, matrix_array_assembler )
       call fe%next()
    end do
    call fe_space%free_fe_iterator(fe)
    deallocate(shape_values, stat=istat); check(istat==0);
    deallocate(shape_gradients, stat=istat); check(istat==0);
    deallocate (elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
  
end module vector_poisson_discrete_integration_names
