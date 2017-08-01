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
module mixed_laplacian_rt_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: mixed_laplacian_rt_discrete_integration_t
     private
     class(scalar_function_t), pointer :: pressure_source_term        => NULL()
     class(scalar_function_t), pointer :: pressure_boundary_function  => NULL()
     type(fe_function_t)     , pointer :: fe_function                 => NULL()
   contains
     procedure :: set_pressure_source_term
     procedure :: set_pressure_boundary_function
     procedure :: set_fe_function
     procedure :: integrate_galerkin
  end type mixed_laplacian_rt_discrete_integration_t
  
  public :: mixed_laplacian_rt_discrete_integration_t
  
contains
   
  subroutine set_pressure_source_term (this, scalar_function)
    implicit none
    class(mixed_laplacian_rt_discrete_integration_t), intent(inout) :: this
    class(scalar_function_t), target, intent(in)    :: scalar_function
    this%pressure_source_term => scalar_function
  end subroutine set_pressure_source_term
  
  subroutine set_pressure_boundary_function (this, scalar_function)
    implicit none
    class(mixed_laplacian_rt_discrete_integration_t), intent(inout) :: this
    class(scalar_function_t), target, intent(in)    :: scalar_function
    this%pressure_boundary_function => scalar_function
  end subroutine set_pressure_boundary_function
  
  subroutine set_fe_function (this, fe_function)
    implicit none
    class(mixed_laplacian_rt_discrete_integration_t), intent(inout) :: this
    type(fe_function_t)                     , target, intent(in)    :: fe_function
    this%fe_function => fe_function
  end subroutine set_fe_function

  subroutine integrate_galerkin ( this, fe_space, matrix_array_assembler )
    implicit none
    class(mixed_laplacian_rt_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                    , intent(inout) :: fe_space
    class(matrix_array_assembler_t)             , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    class(fe_iterator_t)     , allocatable :: fe
    class(fe_face_iterator_t), allocatable :: fe_face
    
    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(face_map_t)         , pointer :: face_map
    type(face_integrator_t)  , pointer :: face_int_velocity
    type(vector_field_t)               :: normals(2)
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(cell_integrator_t), pointer :: cell_int_velocity, cell_int_pressure
    type(vector_field_t), allocatable  :: velocity_shape_values(:,:)
    real(rp)            , allocatable  :: velocity_shape_divs(:,:)
    real(rp)            , allocatable  :: pressure_shape_values(:,:)
    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)
    
    ! FACE matrix and vector, i.e., A_F + f_F
    real(rp), allocatable              :: facemat(:,:,:,:), facevec(:,:)
    
    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    real(rp)     :: factor
    real(rp), allocatable :: pressure_source_term_values(:)
    real(rp), allocatable :: pressure_boundary_function_values(:)

    integer(ip), pointer :: num_dofs_per_field(:) 
    
    assert ( associated(this%pressure_source_term) )
    assert ( associated(this%pressure_boundary_function) )
    
    call fe_space%initialize_fe_integration()
    call fe_space%create_fe_iterator(fe)

    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    num_dofs_per_field => fe%get_number_dofs_per_field()
    quad             => fe%get_quadrature()
    num_quad_points  = quad%get_number_quadrature_points()
    fe_map           => fe%get_fe_map()
    cell_int_velocity => fe%get_cell_integrator(1)
    cell_int_pressure => fe%get_cell_integrator(2)
    
    call memalloc ( num_quad_points, pressure_source_term_values, __FILE__, __LINE__ )
    do while ( .not. fe%has_finished())
      
       ! Update FE-integration related data structures
       call fe%update_integration()

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe_map%get_quadrature_coordinates()
       
       ! Evaluate pressure source term at quadrature points
       call this%pressure_source_term%get_values_set(quad_coords, pressure_source_term_values)
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       call cell_int_velocity%get_values(velocity_shape_values)
       call cell_int_velocity%get_divergences(velocity_shape_divs)
       call cell_int_pressure%get_values(pressure_shape_values)
       do qpoint = 1, num_quad_points
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
          ! \int_(v.u)
          do idof=1, num_dofs_per_field(1)
            do jdof=1, num_dofs_per_field(1)
              elmat(idof,jdof) = elmat(idof,jdof) + &
                                 velocity_shape_values(jdof,qpoint)*velocity_shape_values(idof,qpoint)*factor
            end do
          end do
          
          ! \int_(div(v)*p)
          do idof=1, num_dofs_per_field(1)
            do jdof=1, num_dofs_per_field(2)
              elmat(idof,jdof+num_dofs_per_field(1)) = elmat(idof,jdof+num_dofs_per_field(1)) &
                                                     - velocity_shape_divs(idof,qpoint)*pressure_shape_values(jdof,qpoint)*factor
            end do
          end do
          
          ! \int_(q*div(u))
          do idof=1, num_dofs_per_field(2)
            do jdof=1, num_dofs_per_field(1)
              elmat(idof+num_dofs_per_field(1),jdof) = elmat(idof+num_dofs_per_field(1),jdof) &
                                                     - pressure_shape_values(idof,qpoint)*velocity_shape_divs(jdof,qpoint)*factor
            end do
          end do

          do idof=1, num_dofs_per_field(2)
            elvec(idof+num_dofs_per_field(1)) = elvec(idof+num_dofs_per_field(1)) - &
                                                pressure_shape_values(idof,qpoint) * pressure_source_term_values(qpoint)*factor
          end do
       end do
       
       call fe%assemble( this%fe_function, elmat, elvec, matrix_array_assembler )
       call fe%next()
    end do
    call fe_space%free_fe_iterator(fe)
    call memfree ( pressure_source_term_values, __FILE__, __LINE__ )
    
    call fe_space%initialize_fe_face_integration()

    call memalloc ( num_dofs, num_dofs, 2, 2, facemat, __FILE__, __LINE__ )
    call memalloc ( num_dofs,              2, facevec, __FILE__, __LINE__ )
    
    ! Search for the first boundary face
    call fe_space%create_fe_face_iterator(fe_face)
    do while ( .not. fe_face%is_at_boundary() ) 
       call fe_face%next()
    end do

    quad               => fe_face%get_quadrature()
    num_quad_points    = quad%get_number_quadrature_points()
    face_map           => fe_face%get_face_map()
    face_int_velocity  => fe_face%get_face_integrator(1)
    num_dofs_per_field => fe_face%get_number_dofs_per_field(1)
    
    facemat = 0.0_rp
    call memalloc ( num_quad_points, pressure_boundary_function_values, __FILE__, __LINE__ )
    do while ( .not. fe_face%has_finished() )
       if ( fe_face%is_at_boundary() ) then
         !assert( fe_face%get_set_id() == 1 )
         facevec = 0.0_rp
         call fe_face%update_integration() 
         quad_coords => face_map%get_quadrature_coordinates()
         call this%pressure_boundary_function%get_values_set(quad_coords, pressure_boundary_function_values)
         call face_int_velocity%get_values(1,velocity_shape_values)
         do qpoint = 1, num_quad_points
            factor = face_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            call face_map%get_normals(qpoint,normals)
            do idof = 1, num_dofs_per_field(1)
              facevec(idof,1) = facevec(idof,1) - &
                                pressure_boundary_function_values(qpoint)*velocity_shape_values(idof,qpoint)*normals(1)*factor
            end do   
         end do
         call fe_face%assemble( facemat, facevec, matrix_array_assembler )
       end if
       call fe_face%next()
    end do
    call fe_space%free_fe_face_iterator(fe_face)
    call memfree ( pressure_boundary_function_values, __FILE__, __LINE__ )
    deallocate(velocity_shape_values, stat=istat); check(istat==0);
    call memfree(velocity_shape_divs, __FILE__, __LINE__)
    call memfree(pressure_shape_values, __FILE__, __LINE__)
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module mixed_laplacian_rt_discrete_integration_names
