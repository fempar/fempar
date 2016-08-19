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
  use serial_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_cG_discrete_integration_t
   contains
     procedure :: integrate
  end type poisson_cG_discrete_integration_t
  
  public :: poisson_cG_discrete_integration_t
  
contains
   
  subroutine integrate ( this, fe_space, matrix_array_assembler )
    implicit none
    class(poisson_cG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    type(fe_iterator_t) :: fe_iterator
    type(fe_accessor_t) :: fe
    
    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(volume_integrator_t), pointer :: vol_int
    type(vector_field_t)               :: grad_test, grad_trial
    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: factor

    integer(ip)  :: number_fields

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip), allocatable :: num_dofs_per_field(:)  

    
    number_fields = fe_space%get_number_fields()
    allocate( elem2dof(number_fields), stat=istat); check(istat==0);
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()
    
    call fe_space%initialize_fe_integration()
    
    fe_iterator = fe_space%create_fe_iterator()
    call fe_iterator%current(fe)
    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )
    call fe%get_number_dofs_per_field(num_dofs_per_field)
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_number_quadrature_points()
    fe_map          => fe%get_fe_map()
    vol_int         => fe%get_volume_integrator(1)
    do while ( .not. fe_iterator%has_finished() )
       ! Get current FE
       call fe_iterator%current(fe)
       
       ! Update FE-integration related data structures
       call fe%update_integration()
       
       ! Get DoF numbering within current FE
       call fe%get_elem2dof(elem2dof)

       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       do qpoint = 1, num_quad_points
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
             call vol_int%get_gradient(idof, qpoint, grad_trial)
             do jdof = 1, num_dofs
                call vol_int%get_gradient(jdof, qpoint, grad_test)
                ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * grad_test * grad_trial
             end do
          end do
       end do
       
       ! Apply boundary conditions
       call fe%impose_strong_dirichlet_bcs( elmat, elvec )
       call matrix_array_assembler%assembly( number_fields, num_dofs_per_field, elem2dof, field_blocks, field_coupling, elmat, elvec )
       call fe_iterator%next()
    end do
    deallocate (elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
  
end module poisson_cG_discrete_integration_names
