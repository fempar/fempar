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
module poisson_cG_error_estimator_names
  use fempar_names
  use analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  
  
  type, extends(error_estimator_t) :: poisson_cG_error_estimator_t
    class(base_analytical_functions_t), pointer :: analytical_functions => NULL()
    type(fe_function_t)               , pointer :: fe_function          => NULL()
   contains
    procedure :: set_analytical_functions  => pcGee_set_analytical_functions
    procedure :: set_fe_function           => pcGee_set_fe_function
    procedure :: compute_local_estimates   => pcGee_compute_local_estimates
    procedure :: compute_local_true_errors => pcGee_compute_local_true_errors
  end type poisson_cG_error_estimator_t
  
  public :: poisson_cG_error_estimator_t

contains

 subroutine pcGee_set_analytical_functions ( this, analytical_functions )
   implicit none
   class(poisson_cG_error_estimator_t)         , intent(inout) :: this
   class(base_analytical_functions_t) , target , intent(in)    :: analytical_functions
   this%analytical_functions => analytical_functions
 end subroutine pcGee_set_analytical_functions

 subroutine pcGee_set_fe_function ( this, fe_function )
   implicit none
   class(poisson_cG_error_estimator_t)         , intent(inout) :: this
   type(fe_function_t)                , target , intent(in)    :: fe_function
   this%fe_function => fe_function
 end subroutine pcGee_set_fe_function

 subroutine pcGee_compute_local_estimates(this)
   class(poisson_cG_error_estimator_t), intent(inout) :: this
   class(serial_fe_space_t)    , pointer     :: fe_space
   class(triangulation_t)      , pointer     :: triangulation
   type(std_vector_real_rp_t)  , pointer     :: sq_local_estimates
   real(rp)                    , pointer     :: sq_local_estimate_entries(:)
   class(scalar_function_t)    , pointer     :: source_term
   real(rp)                                  :: source_term_value
   class(fe_cell_iterator_t)   , allocatable :: fe
   type(fe_cell_function_scalar_t)           :: fe_cell_function
   type(fe_cell_function_duties_t)           :: fe_cell_function_duties
   type(quadrature_t)          , pointer     :: quad
   type(point_t)               , pointer     :: quad_coords(:)
   real(rp)                    , pointer     :: laplacians(:)
   real(rp)                                  :: factor, h_length
   integer(ip)                               :: qpoint, num_quad_points, ineigh
   class(fe_facet_iterator_t)  , allocatable :: fe_facet
   type(fe_facet_function_scalar_t)          :: fe_facet_function
   type(vector_field_t)                      :: normals(2)
   type(vector_field_t)                      :: fe_function_gradient
   real(rp)                                  :: sq_local_estimate_value
   real(rp)                                  :: sq_local_face_estimate_value, jump
   type(interpolation_duties_t), allocatable :: interpolation_duties(:)
   type(cell_map_duties_t)                   :: cell_map_duties
   integer(ip)                               :: istat
   
   fe_space      => this%get_fe_space()
   triangulation => fe_space%get_triangulation()
   
   sq_local_estimates         => this%get_sq_local_estimates()
   call sq_local_estimates%resize(0)
   call sq_local_estimates%resize(triangulation%get_num_local_cells(),0.0_rp)
   sq_local_estimate_entries => sq_local_estimates%get_pointer()
   
   assert (associated(this%analytical_functions))
   source_term => this%analytical_functions%get_source_term()
   
   assert (associated(this%fe_function))
   
   assert (fe_space%get_num_fields()==1)
   
   allocate(interpolation_duties(1),stat=istat); check(istat==0)
   call interpolation_duties(1)%assign_compute_first_derivatives(.true.)
   call interpolation_duties(1)%assign_compute_second_derivatives(.true.)
   
   call cell_map_duties%assign_compute_jacobian_inverse(.true.)
   call cell_map_duties%assign_compute_jacobian_derivative(.false.)
   
   call fe_space%set_up_cell_integration(interpolation_duties,cell_map_duties)
   call fe_space%create_fe_cell_iterator(fe)
   call fe_cell_function_duties%assign_evaluate_values(.false.)
   call fe_cell_function_duties%assign_evaluate_gradients(.false.)
   call fe_cell_function_duties%assign_evaluate_laplacians(.true.)
   call fe_cell_function%create(fe_space,1,fe_cell_function_duties)
   do while ( .not. fe%has_finished() )
     if ( fe%is_local() .and. .not. fe%is_void(1) ) then
       sq_local_estimate_value = 0.0_rp
       call fe%update_integration()
       call fe_cell_function%update(fe,this%fe_function)
       laplacians => fe_cell_function%get_quadrature_points_laplacians()
       quad => fe%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       quad_coords => fe%get_quadrature_points_coordinates()
       do qpoint = 1, num_quad_points
         factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
         h_length = fe%compute_characteristic_length(qpoint)
         call source_term%get_value(quad_coords(qpoint),source_term_value)
         sq_local_estimate_value = sq_local_estimate_value + & 
           ( h_length ** 2.0_rp ) * factor * ( source_term_value + laplacians(qpoint) ) ** 2.0_rp
       end do
       sq_local_estimate_entries(fe%get_gid()) = sq_local_estimate_value 
     end if
     call fe%next()
   end do
   call fe_cell_function%free()
   
   call fe_space%set_up_facet_integration()
   call fe_facet_function%create(fe_space,1)
   call fe_space%create_fe_facet_iterator(fe_facet)
   do while ( .not. fe_facet%has_finished() )
     if ( fe_facet%is_local() .and. fe_facet%is_at_field_interior(1) ) then
       sq_local_face_estimate_value = 0.0_rp
       call fe_facet%update_integration()
       quad => fe_facet%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       call fe_facet_function%update(fe_facet,this%fe_function)
       do qpoint = 1, num_quad_points
         jump = 0.0_rp
         call fe_facet%get_normals(qpoint,normals)
         do ineigh = 1, fe_facet%get_num_cells_around()
           call fe_facet_function%get_gradient(qpoint,ineigh,fe_function_gradient)
           jump = jump + normals(ineigh) * fe_function_gradient
         end do
         factor   = fe_facet%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
         h_length = fe_facet%compute_characteristic_length(qpoint)
         sq_local_face_estimate_value = sq_local_face_estimate_value + & 
                                        h_length * factor * ( jump ** 2.0_rp )
       end do
       do ineigh = 1, fe_facet%get_num_cells_around()
         call fe_facet%get_cell_around(ineigh,fe)
         sq_local_estimate_entries(fe%get_gid()) = sq_local_estimate_entries(fe%get_gid()) + &
                                                     0.5_rp * sq_local_face_estimate_value
       end do
     end if
     call fe_facet%next()
   end do
   call fe_facet_function%free()
   call fe_space%free_fe_cell_iterator(fe)
   call fe_space%free_fe_facet_iterator(fe_facet)
   
   if ( allocated(interpolation_duties) ) then
     deallocate(interpolation_duties, stat=istat)
     check(istat==0)
   end if
   
 end subroutine pcGee_compute_local_estimates

 subroutine pcGee_compute_local_true_errors(this)
   class(poisson_cG_error_estimator_t), intent(inout) :: this
   class(serial_fe_space_t)  , pointer     :: fe_space
   class(triangulation_t)    , pointer     :: triangulation
   type(std_vector_real_rp_t), pointer     :: sq_local_true_errors
   real(rp)                  , pointer     :: sq_local_true_error_entries(:)
   class(scalar_function_t)  , pointer     :: exact_solution
   real(rp)                                :: solution_value
   type(fe_cell_function_scalar_t)         :: fe_cell_function
   type(vector_field_t)      , pointer     :: fe_function_gradients(:)
   class(fe_cell_iterator_t) , allocatable :: fe
   type(quadrature_t)        , pointer     :: quad
   type(point_t)             , pointer     :: quad_coords(:)
   type(vector_field_t)      , allocatable :: work_array_gradients(:)
   real(rp)                                :: factor
   integer(ip)                             :: qpoint, num_quad_points
   real(rp)                                :: sq_local_true_error_value
   integer(ip)                             :: istat
   
   fe_space      => this%get_fe_space()
   triangulation => fe_space%get_triangulation()
   
   sq_local_true_errors        => this%get_sq_local_true_errors()
   call sq_local_true_errors%resize(0)
   call sq_local_true_errors%resize(triangulation%get_num_local_cells(),0.0_rp)
   sq_local_true_error_entries => sq_local_true_errors%get_pointer()
   
   assert (associated(this%analytical_functions))
   exact_solution => this%analytical_functions%get_solution_function()
   
   assert (associated(this%fe_function))
   allocate(work_array_gradients(fe_space%get_max_num_quadrature_points()),stat=istat)
   check(istat==0)
   
   call fe_cell_function%create(fe_space,1)
   call fe_space%create_fe_cell_iterator(fe)
   do while(.not. fe%has_finished())
     if ( fe%is_local() ) then
       sq_local_true_error_value = 0.0_rp
       call fe%update_integration()
       call fe_cell_function%update(fe,this%fe_function)
       quad => fe%get_quadrature()
       num_quad_points = quad%get_num_quadrature_points()
       quad_coords => fe%get_quadrature_points_coordinates()
       call exact_solution%get_gradients_set( quad_coords, &
                                              work_array_gradients(1:num_quad_points) )
       fe_function_gradients => fe_cell_function%get_quadrature_points_gradients()
       do qpoint = 1, num_quad_points
         factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
         work_array_gradients(qpoint) = work_array_gradients(qpoint) - fe_function_gradients(qpoint)
         sq_local_true_error_value = sq_local_true_error_value + &
           work_array_gradients(qpoint) * work_array_gradients(qpoint) * factor
       end do
       sq_local_true_error_entries(fe%get_gid()) = sq_local_true_error_value
     end if
     call fe%next()
   end do
   call fe_cell_function%free()
   call fe_space%free_fe_cell_iterator(fe)
   
   deallocate(work_array_gradients, stat=istat)
   check(istat==0)
   
 end subroutine pcGee_compute_local_true_errors

end module poisson_cG_error_estimator_names