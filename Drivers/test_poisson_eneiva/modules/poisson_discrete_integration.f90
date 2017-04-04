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
module poisson_discrete_integration_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_discrete_integration_t
     private
     class(scalar_function_t), pointer :: source_term
     class(scalar_function_t), pointer :: neumann_values
   contains
     procedure :: set_source_term
     procedure :: set_neumann_values
     procedure :: integrate
  end type poisson_discrete_integration_t
  
  public :: poisson_discrete_integration_t
  
contains
   
  subroutine set_source_term (this, scalar_function)
    implicit none
    class(poisson_discrete_integration_t)        , intent(inout) :: this
    class(scalar_function_t)             , target, intent(in)    :: scalar_function
    this%source_term => scalar_function
  end subroutine set_source_term
  
  subroutine set_neumann_values (this, scalar_function)
    implicit none
    class(poisson_discrete_integration_t)        , intent(inout) :: this
    class(scalar_function_t)             , target, intent(in)    :: scalar_function
    this%neumann_values => scalar_function
  end subroutine set_neumann_values

  subroutine integrate ( this, fe_space, matrix_array_assembler )
    implicit none
    class(poisson_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)             , intent(inout) :: fe_space
    class(matrix_array_assembler_t)      , intent(inout) :: matrix_array_assembler

    ! FE space traversal-related data types
    class(fe_accessor_t), allocatable :: fe
    type(fe_face_accessor_t) :: fe_face
    
    ! FE integration-related data types
    type(fe_map_t)           , pointer :: fe_map
    type(quadrature_t)       , pointer :: quad
    type(point_t)            , pointer :: quad_coords(:)
    type(volume_integrator_t), pointer :: vol_int
    real(rp)                           :: shape_test, shape_trial
    type(vector_field_t)               :: grad_test, grad_trial
    
    ! Face integration-related data types
    type(face_map_t)       , pointer :: face_map
    type(face_integrator_t), pointer :: face_int
    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)
    
    ! FACE matrix and vector, i.e., A_F + f_F
    real(rp), allocatable              :: facemat(:,:,:,:), facevec(:,:)

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    integer(ip)  :: ineigh, jneigh
    real(rp)     :: factor
    real(rp)     :: source_term_value, neumann_value

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
    call fe_space%create_fe_accessor(fe)
    call fe%first()    
    
    num_dofs = fe%get_number_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    call memalloc ( number_fields, num_dofs_per_field, __FILE__, __LINE__ )
    call fe%get_number_dofs_per_field(num_dofs_per_field)
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_number_quadrature_points()
    fe_map          => fe%get_fe_map()
    vol_int         => fe%get_volume_integrator(1)
    
    do while ( .not. fe%past_the_end())
       
       ! Update FE-integration related data structures
       call fe%update_integration()
       
       ! Get DoF numbering within current FE
       call fe%get_elem2dof(elem2dof)

       ! Get quadrature coordinates to evaluate boundary value
       quad_coords => fe_map%get_quadrature_coordinates()
       
       ! Compute element matrix and vector
       elmat = 0.0_rp
       elvec = 0.0_rp
       do qpoint = 1, num_quad_points
       
          factor = fe_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
          
          ! Diffusive term
          do idof = 1, num_dofs
             call vol_int%get_gradient(idof, qpoint, grad_trial)
             do jdof = 1, num_dofs
                call vol_int%get_gradient(jdof, qpoint, grad_test)
                ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                elmat(idof,jdof) = elmat(idof,jdof) + factor * grad_test * grad_trial
             end do
          end do
          
          ! Source term
          call this%source_term%get_value(quad_coords(qpoint),source_term_value)
          do idof = 1, num_dofs
             call vol_int%get_value(idof, qpoint, shape_trial)
             elvec(idof) = elvec(idof) + factor * source_term_value * shape_trial
          end do 
          
       end do
       
       ! Apply boundary conditions
       call fe%impose_strong_dirichlet_bcs( elmat, elvec )
       call matrix_array_assembler%assembly( number_fields, num_dofs_per_field, elem2dof, field_blocks, field_coupling, elmat, elvec )
       call fe%next()
    end do
    call fe%free()
    
    call fe_space%initialize_fe_face_integration()
    
    call memalloc ( num_dofs, num_dofs, 2, 2, facemat, __FILE__, __LINE__ )
    call memalloc ( num_dofs,              2, facevec, __FILE__, __LINE__ )
    
    ! Search for the first boundary face
    call fe_space%create_fe_face_accessor(fe_face)
    do while ( .not. fe_face%is_at_boundary() ) 
       call fe_face%next()
    end do
    
    quad            => fe_face%get_quadrature()
    num_quad_points = quad%get_number_quadrature_points()
    face_map        => fe_face%get_face_map()
    face_int        => fe_face%get_face_integrator(1)
    
    do while ( .not. fe_face%past_the_end() )
       facemat = 0.0_rp
       facevec = 0.0_rp
       if ( fe_face%is_at_boundary() .and. fe_face%get_set_id() == 0 ) then
         call fe_face%update_integration()    
         quad_coords => face_map%get_quadrature_coordinates()
         do qpoint = 1, num_quad_points
            call this%neumann_values%get_value(quad_coords(qpoint),neumann_value)
            factor = face_map%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            do idof = 1, num_dofs_per_field(1)
              call face_int%get_value(idof,qpoint,1,shape_trial)   
              facevec(idof,1) = facevec(idof,1) + factor * shape_trial * neumann_value
            end do   
         end do
         call fe_face%get_elem2dof(1, elem2dof)
         call matrix_array_assembler%face_assembly(number_fields, &
                                                   num_dofs_per_field, &
                                                   num_dofs_per_field, &
                                                   elem2dof, &
                                                   elem2dof, &
                                                   field_blocks, &
                                                   field_coupling, &
                                                   facemat(:,:,1,1), &
                                                   facevec(:,1) )            
       end if
       call fe_face%next()
    end do
    
    deallocate(elem2dof, stat=istat); check(istat==0);
    call memfree ( num_dofs_per_field, __FILE__, __LINE__ )
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )
  end subroutine integrate
  
end module poisson_discrete_integration_names
