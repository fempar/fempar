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
module time_stepping_mass_discrete_integration_names
  use types_names
  use memor_names
  use reference_fe_names
  use discrete_integration_names
  use fe_space_names
  use assembler_names
  use vector_names
 
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: mass_discrete_integration_t 
  !< Mass matrix discrete integration associated with the `mass_op` in the `time_stepping_operator_t`
     type(fe_function_t) :: fe_function !< This FE function represents the values of the temporal derivative
   contains
     procedure :: integrate_galerkin
     procedure :: set_evaluation_point
     procedure :: create_fe_function
     procedure :: set_current_time
     procedure :: free
  end type mass_discrete_integration_t
  
  public :: mass_discrete_integration_t
  
contains
   
  subroutine free (this)
    class(mass_discrete_integration_t), intent(inout)    :: this
    call this%fe_function%free()
  end subroutine free
  
  !> As the mass operator is supposed to be `type(fe_affine_operator_t)` only integrate galerkin is needed to
  !> integrate the mass matrix and the restricted BC's.
  !> @note
  !> The force term is not evaluated for integrating the mass matrix, and `elvec` is kept to 0.0
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(mass_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)         , intent(inout) :: fe_space
    class(assembler_t)      , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor

    
    
    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( max_num_dofs, elvec, __FILE__, __LINE__ )
    
    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )
       if ( fe%is_local() ) then
       ! Update FE-integration related data structures
       call fe%update_integration()
          
       ! Very important: this has to be inside the loop, as different FEs can be present!
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()
       

       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) = (phi_i,phi_j)
                elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_values(jdof,qpoint) * shape_values(idof,qpoint)
             end do
          end do
       end do
       
       call fe%assembly( this%fe_function, elmat, elvec, assembler )
       end if
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)

    call memfree(shape_values, __FILE__, __LINE__)
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
  !> The `fe_function_t` associated have a composition relationship with `mass_discrete_integration_t`, 
  !> i.e., is created and freed within the scope of this module.
  subroutine create_fe_function ( this, fe_space ) 
     implicit none
     class(mass_discrete_integration_t), intent(inout)  :: this
     class(serial_fe_space_t)          ,   intent(in)   :: fe_space
     call this%fe_function%create( fe_space )
  end subroutine create_fe_function
     
  subroutine set_evaluation_point ( this, evaluation_point )
    implicit none
    class(mass_discrete_integration_t), intent(inout)  :: this
    class(vector_t)                   , intent(in)     :: evaluation_point
    call this%fe_function%set_free_dof_values(evaluation_point)
  end subroutine set_evaluation_point
  
  !> This subroutine evaluates the temporal derivative at the Dirichlet Boundary for current time
  !> @todo
  !> The name may not be well chosen
  subroutine set_current_time (this, fe_space, current_time)
   implicit none
   class(mass_discrete_integration_t), intent(inout) :: this
   class(serial_fe_space_t)          , intent(in)    :: fe_space
   real(rp)                          , intent(in)    :: current_time
   call fe_space%interpolate_dirichlet_values(this%fe_function, time=current_time , time_derivative_order = 1)
  end subroutine set_current_time
 
end module time_stepping_mass_discrete_integration_names
