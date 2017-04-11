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
module poisson_unfitted_discrete_integration_names
  use fempar_names
  use poisson_unfitted_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_unfitted_cG_discrete_integration_t
     type(poisson_unfitted_analytical_functions_t), pointer :: analytical_functions => NULL()
   contains
     procedure :: set_analytical_functions
     procedure :: integrate
  end type poisson_unfitted_cG_discrete_integration_t
  
  public :: poisson_unfitted_cG_discrete_integration_t
  
contains
   
  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_unfitted_cG_discrete_integration_t)    ,intent(inout)  :: this
     type(poisson_unfitted_analytical_functions_t), target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions


  subroutine integrate ( this, fe_space, matrix_array_assembler )
    implicit none
    class(poisson_unfitted_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    type(fe_iterator_t) :: fe_iterator
    type(fe_accessor_t) :: fe

    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int
    type(vector_field_t), allocatable  :: shape_gradients(:,:)
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor
    real(rp)     :: source_term_value

    integer(ip)  :: number_fields

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip), allocatable :: num_dofs_per_field(:)  
    class(scalar_function_t), pointer :: source_term
    
    assert (associated(this%analytical_functions)) 

    source_term => this%analytical_functions%get_source_term()

    number_fields = fe_space%get_number_fields()
    allocate( elem2dof(number_fields), stat=istat); check(istat==0);
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()

    max_num_dofs = fe_space%get_max_number_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )

    fe_iterator = fe_space%create_fe_iterator()
    do while ( .not. fe_iterator%has_finished() )
       ! Get current FE
       call fe_iterator%current(fe)
       if ( fe%is_local() ) then
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
          
          ! Get quadrature coordinates to evaluate source_term
          quad_coords => fe_map%get_quadrature_coordinates()
                    
          ! Compute element matrix and vector
          elmat = 0.0_rp
          elvec = 0.0_rp
          call vol_int%get_gradients(shape_gradients)
          call vol_int%get_values(shape_values)
          do qpoint = 1, num_quad_points
             factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
             do idof = 1, num_dofs
                do jdof = 1, num_dofs
                   ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                   elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients(jdof,qpoint) * shape_gradients(idof,qpoint)
                end do
             end do
             
             ! Source term
             call source_term%get_value(quad_coords(qpoint),source_term_value)
             do idof = 1, num_dofs
                elvec(idof) = elvec(idof) + factor * source_term_value * shape_values(idof,qpoint) 
             end do
          end do
          
          ! Apply boundary conditions
          call fe%impose_strong_dirichlet_bcs( elmat, elvec )
          call matrix_array_assembler%assembly( number_fields, num_dofs_per_field, elem2dof, field_blocks, field_coupling, elmat, elvec )
       end if
       call fe_iterator%next()
    end do
    call memfree(shape_values, __FILE__, __LINE__)
    deallocate (shape_gradients, stat=istat); check(istat==0);
    deallocate (elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
  
end module poisson_unfitted_discrete_integration_names
