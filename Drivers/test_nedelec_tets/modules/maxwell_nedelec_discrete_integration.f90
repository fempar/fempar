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
module maxwell_nedelec_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: maxwell_nedelec_discrete_integration_t
     private
     class(vector_function_t), pointer :: source_term        => NULL()
   contains
     procedure :: set_source_term
     procedure :: integrate
  end type maxwell_nedelec_discrete_integration_t
  
  public :: maxwell_nedelec_discrete_integration_t
  
contains
   
  subroutine set_source_term (this, vector_function)
    implicit none
    class(maxwell_nedelec_discrete_integration_t), intent(inout) :: this
    class(vector_function_t), target, intent(in)    :: vector_function
    this%source_term => vector_function
  end subroutine set_source_term

  subroutine integrate ( this, fe_space, matrix_array_assembler )
    implicit none
    class(maxwell_nedelec_discrete_integration_t), intent(in)   :: this
    class(serial_fe_space_t)                    , intent(inout) :: fe_space
    class(matrix_array_assembler_t)             , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    type(fe_iterator_t) :: fe_iterator
    type(fe_accessor_t) :: fe
    type(fe_face_iterator_t) :: fe_face_iterator
    type(fe_face_accessor_t) :: fe_face
    
    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int_H
    type(vector_field_t), allocatable  :: H_shape_values(:,:)
    type(vector_field_t), allocatable  :: H_shape_curls(:,:)
    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: factor
    type(vector_field_t), allocatable :: source_term_values(:)

    integer(ip)  :: number_fields

    integer(ip), pointer :: field_blocks(:)
    logical    , pointer :: field_coupling(:,:)

    type(i1p_t), allocatable :: elem2dof(:)
    integer(ip), allocatable :: num_dofs_per_field(:) 
    
    assert ( associated(this%source_term) )
    
    number_fields = fe_space%get_number_fields()
    allocate( elem2dof(number_fields), stat=istat); check(istat==0);
    field_blocks => fe_space%get_field_blocks()
    field_coupling => fe_space%get_field_coupling()
    
    fe_iterator = fe_space%create_fe_iterator()
    call fe_space%initialize_fe_integration()
    
    call fe_iterator%current(fe)
    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )
    call fe%get_number_dofs_per_field(num_dofs_per_field)
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_number_quadrature_points()
    fe_map           => fe%get_fe_map()
    vol_int_H => fe%get_volume_integrator(1)
    allocate (source_term_values(num_quad_points), stat=istat); check(istat==0)
    do while ( .not. fe_iterator%has_finished() )
       ! Get current FE
       call fe_iterator%current(fe)
       
       ! Update FE-integration related data structures
       call fe%update_integration()

       ! Get DoF numbering within current FE
       call fe%get_elem2dof(elem2dof)

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe_map%get_quadrature_coordinates()
       
       ! Evaluate pressure source term at quadrature points
       call this%source_term%get_values_set(quad_coords, source_term_values)
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call vol_int_H%get_values(H_shape_values)
       call vol_int_H%get_curls(H_shape_curls)
       do qpoint = 1, num_quad_points
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          ! \int_(curl(v).curl(u))
          do idof=1, num_dofs_per_field(1)
            do jdof=1, num_dofs_per_field(1)
              elmat(idof,jdof) = elmat(idof,jdof) + &
                                 ( H_shape_values(idof,qpoint)*H_shape_values(jdof,qpoint) + 0.0_rp*H_shape_curls(jdof,qpoint)*H_shape_curls(idof,qpoint) )*factor
            end do
            ! \int_(curl(v).f)
            elvec(idof) = elvec(idof) + H_shape_values(idof,qpoint) * source_term_values(qpoint) * factor
          end do
       end do
       !call fe%impose_strong_dirichlet_bcs( elmat, elvec )
       call matrix_array_assembler%assembly( number_fields, num_dofs_per_field, elem2dof, field_blocks, field_coupling, elmat, elvec )
       call fe_iterator%next()
    end do
    deallocate(H_shape_curls, stat=istat); check(istat==0);
    deallocate(H_shape_values, stat=istat); check(istat==0);
    deallocate (source_term_values, stat=istat); check(istat==0);
    deallocate (elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
  end subroutine integrate
  
end module maxwell_nedelec_discrete_integration_names
