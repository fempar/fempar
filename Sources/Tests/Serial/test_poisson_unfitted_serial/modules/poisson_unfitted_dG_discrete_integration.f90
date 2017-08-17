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
module poisson_unfitted_dG_discrete_integration_names
  use fempar_names
  use poisson_unfitted_analytical_functions_names
  use poisson_unfitted_conditions_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(discrete_integration_t) :: poisson_unfitted_dG_discrete_integration_t
     type(poisson_unfitted_analytical_functions_t), pointer :: analytical_functions => NULL()
     type(poisson_unfitted_conditions_t)          , pointer :: poisson_unfitted_conditions   => NULL()
   contains
     procedure :: set_analytical_functions
     procedure :: set_poisson_unfitted_conditions
     procedure :: integrate_galerkin
  end type poisson_unfitted_dG_discrete_integration_t
  
  public :: poisson_unfitted_dG_discrete_integration_t
  
contains

  subroutine set_analytical_functions ( this, analytical_functions )
     implicit none
     class(poisson_unfitted_dG_discrete_integration_t)        , intent(inout) :: this
     type(poisson_unfitted_analytical_functions_t)    , target, intent(in)    :: analytical_functions
     this%analytical_functions => analytical_functions
  end subroutine set_analytical_functions
  
  subroutine set_poisson_unfitted_conditions ( this, poisson_unfitted_conditions )
     implicit none
     class(poisson_unfitted_dG_discrete_integration_t)        , intent(inout) :: this
     type(poisson_unfitted_conditions_t)              , target, intent(in)    :: poisson_unfitted_conditions
     this%poisson_unfitted_conditions => poisson_unfitted_conditions
  end subroutine set_poisson_unfitted_conditions
  
  
  subroutine integrate_galerkin ( this, fe_space, assembler )
    implicit none
    class(poisson_unfitted_dG_discrete_integration_t), intent(in)    :: this
    class(serial_fe_space_t)                , intent(inout) :: fe_space
    class(assembler_t)         , intent(inout) :: assembler

    ! FE space traversal-related data types
    class(fe_cell_iterator_t)     , allocatable :: fe
    class(fe_facet_iterator_t), allocatable :: fe_face
    
    ! FE integration-related data types
    type(quadrature_t)       , pointer     :: quad
    type(point_t)            , pointer     :: quad_coords(:)
    type(vector_field_t)     , allocatable, target :: shape_gradients_first(:,:), shape_gradients_second(:,:)
    type(vector_field_t)     , pointer     :: shape_gradients_ineigh(:,:),shape_gradients_jneigh(:,:)
    real(rp)                 , allocatable, target :: shape_values_first(:,:), shape_values_second(:,:)
    real(rp)                 , pointer     :: shape_values_ineigh(:,:),shape_values_jneigh(:,:)
    type(i1p_t)              , allocatable :: fe_dofs(:)
    
    ! Face integration-related data types
    type(facet_integrator_t), pointer :: facet_int
    type(vector_field_t)             :: normals(2)
    real(rp)                         :: shape_test, shape_trial
    real(rp)                         :: h_length
    type(i1p_t)        , allocatable :: trial_fe_dofs(:), test_fe_dofs(:)
    
    ! FE matrix and vector i.e., A_K + f_K
    real(rp), allocatable              :: elmat(:,:), elvec(:)
    
    ! FACE matrix and vector, i.e., A_F + f_F
    real(rp), allocatable              :: facemat(:,:,:,:), facevec(:,:)
    
    ! Problem and dG discretization related parameters 
    real(rp) :: viscosity
    real(rp) :: C_IP        ! Interior Penalty constant
    
    class(scalar_function_t), pointer :: source_term, boundary_function
    type(fe_function_t)               :: boundary_fe_function
    type(fe_facet_function_scalar_t)   :: boundary_fe_facet_function
    real(rp) :: source_term_value, boundary_value, boundary_fe_function_value

    integer(ip)  :: istat
    integer(ip)  :: qpoint, num_quad_points
    integer(ip)  :: idof, jdof, num_dofs
    integer(ip)  :: ineigh, jneigh
    real(rp)     :: factor
    
    assert (associated(this%analytical_functions))
    
    source_term => this%analytical_functions%get_source_term()
    call this%poisson_unfitted_conditions%get_function(1,1,boundary_function)
    
    call boundary_fe_function%create(fe_space)
    call fe_space%interpolate(1,boundary_function,boundary_fe_function)
    call boundary_fe_facet_function%create(fe_space,1)
    
    call fe_space%set_up_cell_integration()
    call fe_space%create_fe_cell_iterator(fe)
    
    num_dofs = fe%get_num_dofs()
    call memalloc ( num_dofs, num_dofs, elmat, __FILE__, __LINE__ )
    call memalloc ( num_dofs, elvec, __FILE__, __LINE__ )
    quad            => fe%get_quadrature()
    num_quad_points = quad%get_num_quadrature_points()
    
    viscosity = 1.0_rp
    C_IP      = 10.0_rp * fe%get_order(1)**2
    
    do while ( .not. fe%has_finished())

       if ( fe%is_local() ) then
       
         ! Update FE-integration related data structures
         call fe%update_integration()
         
         ! Get quadrature coordinates to evaluate source_term
         quad_coords => fe%get_quadrature_points_coordinates()

         ! Compute element matrix and vector
         elmat = 0.0_rp
         elvec = 0.0_rp
         call fe%get_gradients(shape_gradients_first)
         call fe%get_values(shape_values_first)
         do qpoint = 1, num_quad_points
            factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            do idof = 1, num_dofs
               do jdof = 1, num_dofs
                  ! A_K(i,j) = (grad(phi_i),grad(phi_j))
                  elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_gradients_first(jdof,qpoint) * shape_gradients_first(idof,qpoint)
               end do
            end do
            
            ! Source term
            call source_term%get_value(quad_coords(qpoint),source_term_value)
            do idof = 1, num_dofs
               elvec(idof) = elvec(idof) + factor * source_term_value * shape_values_first(idof,qpoint)
            end do  
         end do        
         
         call fe%assembly( elmat, elvec, assembler )
       end if
       
       call fe%next()
    end do
    call fe_space%free_fe_cell_iterator(fe)
    
    call fe_space%set_up_facet_integration()
    
    call memalloc ( num_dofs, num_dofs, 2, 2, facemat, __FILE__, __LINE__ )
    call memalloc ( num_dofs,              2, facevec, __FILE__, __LINE__ )
    
    ! Search for the first interior face
    call fe_space%create_fe_facet_iterator(fe_face)
    do while ( fe_face%is_at_boundary() ) 
       call fe_face%next()
    end do
    
    quad            => fe_face%get_quadrature()
    num_quad_points = quad%get_num_quadrature_points()
    facet_int        => fe_face%get_facet_integrator(1)
    
    do while ( .not. fe_face%has_finished() ) 
       
       if ( .not. fe_face%is_at_boundary() ) then
         facemat = 0.0_rp
         call fe_face%update_integration()    
         
         call facet_int%get_values(1,shape_values_first)
         call facet_int%get_values(2,shape_values_second)
         call facet_int%get_gradients(1,shape_gradients_first)
         call facet_int%get_gradients(2,shape_gradients_second)
         do qpoint = 1, num_quad_points
            call fe_face%get_normals(qpoint,normals)
            h_length = fe_face%compute_characteristic_length(qpoint)
            factor = fe_face%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            do ineigh = 1, fe_face%get_num_cells_around()
               if (ineigh==1) then
                 shape_values_ineigh    => shape_values_first
                 shape_gradients_ineigh => shape_gradients_first
               else if (ineigh==2) then
                 shape_values_ineigh    => shape_values_second
                 shape_gradients_ineigh => shape_gradients_second
               end if
               
               do jneigh = 1, fe_face%get_num_cells_around()

                 if (jneigh==1) then
                   shape_values_jneigh    => shape_values_first
                   shape_gradients_jneigh => shape_gradients_first
                 else if (jneigh==2) then
                   shape_values_jneigh    => shape_values_second
                   shape_gradients_jneigh => shape_gradients_second
                 end if
 
                 do idof = 1, num_dofs
                  !call facet_int%get_value(idof,qpoint,ineigh,shape_trial)
                  !call facet_int%get_gradient(idof,qpoint,ineigh,grad_trial)
                     do jdof = 1, num_dofs
                        !call facet_int%get_value(jdof,qpoint,jneigh,shape_test)
                        !call facet_int%get_gradient(jdof,qpoint,jneigh,grad_test)
                        !- mu*({{grad u}}[[v]] + (1-xi)*[[u]]{{grad v}} ) + C*mu*p^2/h * [[u]] [[v]]
                        facemat(idof,jdof,ineigh,jneigh) = facemat(idof,jdof,ineigh,jneigh) +     &
                             &  factor * viscosity *   &
                             &  (-0.5_rp*shape_gradients_jneigh(jdof,qpoint)*normals(ineigh)*shape_values_ineigh(idof,qpoint) - &
                             &   0.5_rp*shape_gradients_ineigh(idof,qpoint)*normals(jneigh)*shape_values_jneigh(jdof,qpoint)   + &
                             &   c_IP / h_length * shape_values_jneigh(jdof,qpoint)*shape_values_ineigh(idof,qpoint) *        &
                             &   normals(ineigh)*normals(jneigh))
                     end do
                 end do
               end do
            end do
         end do
         call fe_face%assembly( facemat, assembler )
       end if
         
       call fe_face%next()
    end do
    
    ! Search for the first boundary face
    call fe_face%first()
    do while ( .not. fe_face%is_at_boundary() ) 
       call fe_face%next()
    end do

    quad            => fe_face%get_quadrature()
    num_quad_points = quad%get_num_quadrature_points()
    facet_int        => fe_face%get_facet_integrator(1)
   
    do while ( .not. fe_face%has_finished() )
       
       if ( fe_face%is_at_boundary() ) then
         facemat = 0.0_rp
         facevec = 0.0_rp
         assert( fe_face%get_set_id() == 1 )
         call fe_face%update_integration()
         call boundary_fe_facet_function%update(fe_face,boundary_fe_function)
         quad_coords => fe_face%get_quadrature_points_coordinates()
         call facet_int%get_values(1,shape_values_first)
         call facet_int%get_gradients(1,shape_gradients_first)
         do qpoint = 1, num_quad_points
            call fe_face%get_normals(qpoint,normals)
            h_length = fe_face%compute_characteristic_length(qpoint)
            factor = fe_face%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
            call boundary_function%get_value(quad_coords(qpoint),boundary_value)
            call boundary_fe_facet_function%get_value(qpoint,1,boundary_fe_function_value)
            boundary_value = 2*boundary_value - boundary_fe_function_value
            do idof = 1, num_dofs
              !call facet_int%get_value(idof,qpoint,1,shape_trial)
              !call facet_int%get_gradient(idof,qpoint,1,grad_trial)   
              do jdof = 1, num_dofs
                 !call facet_int%get_value(jdof,qpoint,1,shape_test)
                 !call facet_int%get_gradient(jdof,qpoint,1,grad_test)
                 facemat(idof,jdof,1,1) = facemat(idof,jdof,1,1) + &
                                     &  factor * viscosity *   &
                                     (-shape_gradients_first(jdof,qpoint)*normals(1)*shape_values_first(idof,qpoint) - &
                                      shape_gradients_first(idof,qpoint)*normals(1)*shape_values_first(jdof,qpoint)  + &
                                      c_IP / h_length * shape_values_first(idof,qpoint)*shape_values_first(jdof,qpoint))
              end do
              facevec(idof,1) = facevec(idof,1) + factor * viscosity * &
                                      (-boundary_value * shape_gradients_first(idof,qpoint) * normals(1) + &
                                      c_IP/h_length * boundary_value * shape_values_first(idof,qpoint) ) 
            end do   
         end do
         call fe_face%assembly( facemat, facevec, assembler )
       end if
       call fe_face%next()
    end do
    call fe_space%free_fe_facet_iterator(fe_face)
    call boundary_fe_function%free()
    call boundary_fe_facet_function%free()
    call memfree(shape_values_first, __FILE__, __LINE__) 
    call memfree(shape_values_second, __FILE__, __LINE__) 
    deallocate(shape_gradients_first, stat=istat); check(istat==0);
    deallocate(shape_gradients_second, stat=istat); check(istat==0);
    call memfree ( elmat, __FILE__, __LINE__ )
    call memfree ( elvec, __FILE__, __LINE__ )
    call memfree ( facemat, __FILE__, __LINE__ )
    call memfree ( facevec, __FILE__, __LINE__ )
  end subroutine integrate_galerkin
  
end module poisson_unfitted_dG_discrete_integration_names
