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
module mass_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(linear_discrete_integration_t) :: mass_discrete_integration_t
   contains
     procedure :: integrate_galerkin
  end type mass_discrete_integration_t
  
  public :: mass_discrete_integration_t
  
contains
   
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(mass_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)          , intent(inout) :: fe_space
    class(assembler_t)                , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t), allocatable :: fe

    ! FE integration-related data types
    type(quadrature_t)       , pointer :: quad
    real(rp)            , allocatable  :: shape_values(:,:)

    ! FE matrix i.e., A_K 
    real(rp), allocatable              :: elmat(:,:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs, max_num_dofs
    real(rp)     :: factor

    max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
    call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )
    
    call fe_space%create_fe_cell_iterator(fe)
    do while ( .not. fe%has_finished() )
       
       ! Update FE-integration related data structures
       call fe%update_integration()
          
       quad            => fe%get_quadrature()
       num_quad_points =  quad%get_num_quadrature_points()
       num_dofs        =  fe%get_num_dofs()
       
       ! Compute element matrix
       elmat = 0.0_rp
       call fe%get_values(shape_values)
       do qpoint = 1, num_quad_points
          factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          do idof = 1, num_dofs
             do jdof = 1, num_dofs
                ! A_K(i,j) = ((phi_i),phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_values(jdof,qpoint) * shape_values(idof,qpoint)
             end do
          end do 
       end do
       
       call fe%assembly( elmat, assembler )
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)

    call memfree(shape_values, __FILE__, __LINE__)
    call memfree ( elmat, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module mass_discrete_integration_names
