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
module maxwell_analytical_error_estimator_names
  use fempar_names
  use maxwell_analytical_functions_names
  
  implicit none
# include "debug.i90"
  private
  
  
  type, extends(error_estimator_t) :: maxwell_analytical_error_estimator_t
    type(maxwell_analytical_functions_t), pointer     :: analytical_functions => NULL()
    type(fe_function_t)                 , pointer     :: fe_function          => NULL()
    type(interpolation_duties_t)        , allocatable :: interpolation_duties(:)
    type(cell_map_duties_t)                           :: cell_map_duties
    type(fe_cell_function_vector_t)                   :: fe_cell_function
    type(fe_cell_function_duties_t)                   :: fe_cell_function_duties
    type(vector_field_t)                , allocatable :: work_array_values(:) 
    type(vector_field_t)                , allocatable :: work_array_curls(:) 
    type(vector_field_t)                , pointer     :: fe_function_values(:)
    type(vector_field_t)                , allocatable :: fe_function_curls(:) 
   contains
    procedure :: create                    => maee_create
    procedure :: free                      => maee_free
    procedure :: set_analytical_functions  => maee_set_analytical_functions
    procedure :: set_fe_function           => maee_set_fe_function
    procedure :: compute_local_estimates   => maee_compute_local_estimates
    procedure :: compute_local_true_errors => maee_compute_local_true_errors
    procedure :: get_error_norm_exponent   => maee_get_error_norm_exponent
  end type maxwell_analytical_error_estimator_t
  
  public :: maxwell_analytical_error_estimator_t

contains

 subroutine maee_create ( this, fe_space, parameter_list )
   implicit none
   class(maxwell_analytical_error_estimator_t), intent(inout) :: this
   class(serial_fe_space_t)   , target, intent(in)    :: fe_space
   type(parameterlist_t)              , intent(in)    :: parameter_list
   class(environment_t), pointer   :: environment
   integer(ip)                   :: istat
   call ee_create(this,fe_space,parameter_list)
   assert (fe_space%get_num_fields()==1)
   environment => this%get_environment()
   if ( environment%am_i_l1_task() ) then
     allocate(this%interpolation_duties(1),stat=istat); check(istat==0)
     call this%cell_map_duties%assign_compute_jacobian_inverse(.true.)
     call this%cell_map_duties%assign_compute_jacobian_derivative(.false.)
     call this%fe_cell_function%create(fe_space,1)
     allocate(this%work_array_values(fe_space%get_max_num_quadrature_points()),stat=istat)
     allocate(this%work_array_curls(fe_space%get_max_num_quadrature_points()),stat=istat)
     allocate(this%fe_function_curls(fe_space%get_max_num_quadrature_points()),stat=istat)
     check(istat==0)
   end if
 end subroutine maee_create
 
 subroutine maee_free ( this )
   implicit none
   class(maxwell_analytical_error_estimator_t), intent(inout) :: this
   integer(ip)                   :: istat
   call this%fe_cell_function_duties%assign_evaluate_values(.false.)
   call this%fe_cell_function_duties%assign_evaluate_gradients(.false.)
   call this%fe_cell_function_duties%assign_evaluate_laplacians(.false.)
   call this%cell_map_duties%assign_compute_jacobian_inverse(.false.)
   call this%cell_map_duties%assign_compute_jacobian_derivative(.false.)
   if ( allocated(this%interpolation_duties) ) then
     deallocate(this%interpolation_duties, stat=istat)
     check(istat==0)
   end if
   call this%fe_cell_function%free()
   if ( allocated(this%work_array_values) ) then
     deallocate(this%work_array_values, stat=istat)
     check(istat==0)
   end if
   if ( allocated(this%work_array_curls) ) then
     deallocate(this%work_array_curls, stat=istat)
     check(istat==0)
   end if
   if ( allocated(this%fe_function_curls) ) then
     deallocate(this%fe_function_curls, stat=istat)
     check(istat==0)
   end if
   this%fe_function_values => NULL() 
   call ee_free(this)
 end subroutine maee_free

 subroutine maee_set_analytical_functions ( this, analytical_functions )
   implicit none
   class(maxwell_analytical_error_estimator_t)           , intent(inout) :: this
   type(maxwell_analytical_functions_t) , target , intent(in)    :: analytical_functions
   this%analytical_functions => analytical_functions
 end subroutine maee_set_analytical_functions

 subroutine maee_set_fe_function ( this, fe_function )
   implicit none
   class(maxwell_analytical_error_estimator_t)         , intent(inout) :: this
   type(fe_function_t)                , target , intent(in)    :: fe_function
   this%fe_function => fe_function
 end subroutine maee_set_fe_function

 subroutine maee_compute_local_estimates(this)
   class(maxwell_analytical_error_estimator_t), intent(inout) :: this
   class(environment_t)      , pointer     :: environment
   class(serial_fe_space_t)  , pointer     :: fe_space
   class(triangulation_t)    , pointer     :: triangulation
   type(std_vector_real_rp_t), pointer     :: sq_local_estimates
   real(rp)                  , pointer     :: sq_local_estimate_entries(:)
   class(vector_function_t)  , pointer     :: exact_solution
   real(rp)                                :: solution_value
   class(fe_cell_iterator_t) , allocatable :: fe
   type(quadrature_t)        , pointer     :: quad
   type(point_t)             , pointer     :: quad_coords(:)
   real(rp)                                :: factor
   integer(ip)                             :: qpoint, num_quad_points
   real(rp)                                :: sq_local_estimate_value
   real(rp)                                :: sq_local_curl_estimate_value 
   integer(ip) :: istat 
   
   environment => this%get_environment()
   
   if ( environment%am_i_l1_task() ) then
   
     fe_space      => this%get_fe_space()
     assert ( associated(fe_space) )
     triangulation => fe_space%get_triangulation()
     
     sq_local_estimates         => this%get_sq_local_estimates()
     call sq_local_estimates%resize(0)
     call sq_local_estimates%resize(triangulation%get_num_local_cells(), 0.0_rp)
     sq_local_estimate_entries => sq_local_estimates%get_pointer()
     
     assert (associated(this%analytical_functions))
     exact_solution => this%analytical_functions%get_solution_function()
     assert (associated(this%fe_function))
     
     !call this%interpolation_duties(1)%assign_compute_first_derivatives(.true.)
     !call this%interpolation_duties(1)%assign_compute_second_derivatives(.false.)
     !call fe_space%set_up_cell_integration(this%interpolation_duties,this%cell_map_duties)
     
     !call this%fe_cell_function_duties%assign_evaluate_values(.true.)
     !call this%fe_cell_function_duties%assign_evaluate_gradients(.false.)
     !call this%fe_cell_function_duties%assign_evaluate_laplacians(.false.)
     !call this%fe_cell_function%set_duties(this%fe_cell_function_duties)
     
     call fe_space%create_fe_cell_iterator(fe)
     do while(.not. fe%has_finished())
       if ( fe%is_local() .and. (.not. fe%is_void(1)) ) then
         sq_local_estimate_value      = 0.0_rp
         sq_local_curl_estimate_value = 0.0_rp
         call fe%update_integration()
         call this%fe_cell_function%update(fe,this%fe_function)
         quad => fe%get_quadrature()
         num_quad_points = quad%get_num_quadrature_points()
         
         if ( size(this%work_array_values) < num_quad_points ) then  
         deallocate( this%work_array_values, stat=istat); check(istat==0) 
         allocate( this%work_array_values(num_quad_points), stat=istat); check(istat==0)
         end if 
         
         if ( size(this%work_array_curls) < num_quad_points ) then  
         deallocate( this%work_array_curls, stat=istat); check(istat==0) 
         allocate( this%work_array_curls(num_quad_points), stat=istat); check(istat==0)
         end if 
         
         if ( size(this%fe_function_curls) < num_quad_points ) then  
         deallocate( this%fe_function_curls, stat=istat); check(istat==0) 
         allocate( this%fe_function_curls(num_quad_points), stat=istat); check(istat==0)
         end if
                  
         quad_coords => fe%get_quadrature_points_coordinates()
         call exact_solution%get_values_set( quad_coords, &
                                             this%work_array_values(1:num_quad_points) )         
         call exact_solution%get_curls_set( quad_coords, &
                                             this%work_array_curls(1:num_quad_points) ) 
         this%fe_function_values => this%fe_cell_function%get_quadrature_points_values()
         call this%fe_cell_function%compute_quadrature_points_curl_values(this%fe_function_curls)
         do qpoint = 1, num_quad_points
           factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
           
           this%work_array_values(qpoint) = this%work_array_values(qpoint) - this%fe_function_values(qpoint)
           this%work_array_curls(qpoint)  = this%work_array_curls(qpoint) - this%fe_function_curls(qpoint)  
                       
           sq_local_estimate_value      = sq_local_estimate_value      + factor * this%work_array_values(qpoint) * this%work_array_values(qpoint) 
           sq_local_curl_estimate_value = sq_local_curl_estimate_value + factor * this%work_array_curls(qpoint) * this%work_array_curls(qpoint) 
         end do
         sq_local_estimate_entries(fe%get_gid()) = sq_local_estimate_value + sq_local_curl_estimate_value 
       end if
       call fe%next()
     end do
     call fe_space%free_fe_cell_iterator(fe)
   end if
   
 end subroutine maee_compute_local_estimates
 
 subroutine maee_compute_local_true_errors(this)
   class(maxwell_analytical_error_estimator_t), intent(inout) :: this
   class(environment_t)      , pointer     :: environment
   class(serial_fe_space_t)  , pointer     :: fe_space
   class(triangulation_t)    , pointer     :: triangulation
   type(std_vector_real_rp_t), pointer     :: sq_local_true_errors
   real(rp)                  , pointer     :: sq_local_true_error_entries(:)
   class(vector_function_t)  , pointer     :: exact_solution
   real(rp)                                :: solution_value
   type(vector_field_t)      , pointer     :: fe_function_values(:)
   class(fe_cell_iterator_t) , allocatable :: fe
   type(quadrature_t)        , pointer     :: quad
   type(point_t)             , pointer     :: quad_coords(:)
   real(rp)                                :: factor
   integer(ip)                             :: qpoint, num_quad_points
   real(rp)                                :: sq_local_true_error_value
   
   !environment => this%get_environment()
   
   !if ( environment%am_i_l1_task() ) then
   
   !  fe_space      => this%get_fe_space()
   !  assert ( associated(fe_space) )
   !  triangulation => fe_space%get_triangulation()
   !  
   !  sq_local_true_errors        => this%get_sq_local_true_errors()
   !  call sq_local_true_errors%resize(0)
   !  call sq_local_true_errors%resize(triangulation%get_num_local_cells(), 0.0_rp)
   !  sq_local_true_error_entries => sq_local_true_errors%get_pointer()
   !  
   !  assert (associated(this%analytical_functions))
   !  exact_solution => this%analytical_functions%get_solution_function()
   !  
   !  assert (associated(this%fe_function))
   !  
   !  call this%interpolation_duties(1)%assign_compute_first_derivatives(.false.)
   !  call this%interpolation_duties(1)%assign_compute_second_derivatives(.false.)
   !  call fe_space%set_up_cell_integration(this%interpolation_duties,this%cell_map_duties)
   !  
   !  call this%fe_cell_function_duties%assign_evaluate_values(.true.)
   !  call this%fe_cell_function_duties%assign_evaluate_gradients(.false.)
   !  call this%fe_cell_function_duties%assign_evaluate_laplacians(.false.)
   !  call this%fe_cell_function%set_duties(this%fe_cell_function_duties)
   !  
   !  call fe_space%create_fe_cell_iterator(fe)
   !  do while(.not. fe%has_finished())
   !    if ( fe%is_local() ) then
   !      sq_local_true_error_value = 0.0_rp
   !      call fe%update_integration()
   !      call this%fe_cell_function%update(fe,this%fe_function)
   !      quad => fe%get_quadrature()
   !      num_quad_points = quad%get_num_quadrature_points()
   !      quad_coords => fe%get_quadrature_points_coordinates()
   !      call exact_solution%get_values_set( quad_coords, &
   !                                          this%work_array_values(1:num_quad_points) )
   !      fe_function_values => this%fe_cell_function%get_quadrature_points_values()
   !      do qpoint = 1, num_quad_points
   !        factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
   !        this%work_array_values(qpoint) = this%work_array_values(qpoint) - fe_function_values(qpoint)
   !        sq_local_true_error_value = sq_local_true_error_value + &
   !          this%work_array_values(qpoint) * this%work_array_values(qpoint) * factor
   !      end do
   !      sq_local_true_error_entries(fe%get_gid()) = sq_local_true_error_value
   !    end if
   !    call fe%next()
   !  end do
   !  call fe_space%free_fe_cell_iterator(fe)
   !end if
   
 end subroutine maee_compute_local_true_errors

 function maee_get_error_norm_exponent(this)
   class(maxwell_analytical_error_estimator_t), intent(in) :: this
   real(rp) :: maee_get_error_norm_exponent
   maee_get_error_norm_exponent = 0.5_rp
 end function maee_get_error_norm_exponent

end module maxwell_analytical_error_estimator_names
