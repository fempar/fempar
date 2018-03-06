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
module fe_cell_function_names
  use types_names
  use list_types_names
  use reference_fe_names
  use fe_space_names
  use environment_names
  use field_names
  
  implicit none
# include "debug.i90"
  private
  
#define duties fe_cell_function_duties
#define task_01 evaluate_values
#define task_02 evaluate_gradients
#define task_03 evaluate_laplacians
#include "duties_header.i90"
  
  type fe_cell_function_scalar_t
   private
   type(fe_cell_function_duties_t)   :: my_duties
   integer(ip)                       :: field_id
   integer(ip)                       :: current_num_nodes             
   integer(ip)                       :: current_num_quadrature_points
   real(rp)            , allocatable :: nodal_values(:)  
   real(rp)            , allocatable :: quadrature_points_values(:)
   type(vector_field_t), allocatable :: quadrature_points_gradients(:)
   real(rp)            , allocatable :: quadrature_points_laplacians(:)
  contains
     procedure, non_overridable :: create                                => fe_cell_function_scalar_create
     procedure, non_overridable :: update                                => fe_cell_function_scalar_update
     procedure, non_overridable :: move_alloc_values_in                  => fe_cell_function_scalar_move_alloc_values_in
     procedure, non_overridable :: move_alloc_values_out                 => fe_cell_function_scalar_move_alloc_values_out
     procedure, non_overridable :: move_alloc_gradients_in               => fe_cell_function_scalar_move_alloc_gradients_in
     procedure, non_overridable :: move_alloc_gradients_out              => fe_cell_function_scalar_move_alloc_gradients_out
     procedure, non_overridable :: get_field_id                          => fe_cell_function_scalar_get_field_id
     procedure, non_overridable :: get_nodal_values                      => fe_cell_function_scalar_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values          => fe_cell_function_scalar_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients       => fe_cell_function_scalar_get_quadrature_points_gradients
     procedure, non_overridable :: get_quadrature_points_laplacians      => fe_cell_function_scalar_get_quadrature_points_laplacians
     procedure, non_overridable :: get_value                             => fe_cell_function_scalar_get_value
     procedure, non_overridable :: get_gradient                          => fe_cell_function_scalar_get_gradient
     procedure, non_overridable :: get_laplacian                         => fe_cell_function_scalar_get_laplacian
     procedure, non_overridable :: set_current_num_nodes              => fe_cell_function_scalar_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points  => fe_cell_function_scalar_set_current_num_quadrature_points
     procedure, non_overridable :: free                                  => fe_cell_function_scalar_free
  end type fe_cell_function_scalar_t
  
  type fe_cell_function_vector_t
   private
   type(fe_cell_function_duties_t)   :: my_duties
   integer(ip)                       :: field_id
   integer(ip)                       :: current_num_nodes             
   integer(ip)                       :: current_num_quadrature_points           
   real(rp)            , allocatable :: nodal_values(:)  
   type(vector_field_t), allocatable :: quadrature_points_values(:)
   type(tensor_field_t), allocatable :: quadrature_points_gradients(:)
   type(vector_field_t), allocatable :: quadrature_points_laplacians(:)
  contains
     procedure, non_overridable :: create                                => fe_cell_function_vector_create
     procedure, non_overridable :: update                                => fe_cell_function_vector_update
     procedure, non_overridable :: move_alloc_values_in                  => fe_cell_function_vector_move_alloc_values_in
     procedure, non_overridable :: move_alloc_values_out                 => fe_cell_function_vector_move_alloc_values_out
     procedure, non_overridable :: move_alloc_gradients_in               => fe_cell_function_vector_move_alloc_gradients_in
     procedure, non_overridable :: move_alloc_gradients_out              => fe_cell_function_vector_move_alloc_gradients_out
     procedure, non_overridable :: get_field_id                          => fe_cell_function_vector_get_field_id
     procedure, non_overridable :: get_nodal_values                      => fe_cell_function_vector_get_nodal_values      
     procedure, non_overridable :: get_quadrature_points_values          => fe_cell_function_vector_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients       => fe_cell_function_vector_get_quadrature_points_gradients
     procedure, non_overridable :: get_quadrature_points_laplacians      => fe_cell_function_vector_get_quadrature_points_laplacians
     procedure, non_overridable :: compute_quadrature_points_curl_values => fe_cell_function_vector_compute_quadrature_points_curl_values
     procedure, non_overridable :: get_value                             => fe_cell_function_vector_get_value
     procedure, non_overridable :: get_gradient                          => fe_cell_function_vector_get_gradient 
     procedure, non_overridable :: get_laplacian                         => fe_cell_function_vector_get_laplacian
     procedure, non_overridable :: compute_curl                          => fe_cell_function_vector_compute_curl 
     procedure, non_overridable :: compute_divergence                    => fe_cell_function_vector_compute_divergence 
     procedure, non_overridable :: set_current_num_nodes              => fe_cell_function_vector_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points  => fe_cell_function_vector_set_current_num_quadrature_points
     procedure, non_overridable :: free                                  => fe_cell_function_vector_free
  end type fe_cell_function_vector_t
  
  type fe_cell_function_tensor_t
   private
   type(fe_cell_function_duties_t)   :: my_duties
   integer(ip)                       :: field_id
   integer(ip)                       :: current_num_nodes            
   integer(ip)                       :: current_num_quadrature_points     
   real(rp)            , allocatable :: nodal_values(:)
   type(tensor_field_t), allocatable :: quadrature_points_values(:)
  contains
     procedure, non_overridable :: create                               => fe_cell_function_tensor_create
     procedure, non_overridable :: update                               => fe_cell_function_tensor_update
     procedure, non_overridable :: move_alloc_values_in                 => fe_cell_function_tensor_move_alloc_values_in
     procedure, non_overridable :: move_alloc_values_out                => fe_cell_function_tensor_move_alloc_values_out
     procedure, non_overridable :: get_field_id                         => fe_cell_function_tensor_get_field_id
     procedure, non_overridable :: get_nodal_values                     => fe_cell_function_tensor_get_nodal_values  
     procedure, non_overridable :: get_quadrature_points_values         => fe_cell_function_tensor_get_quadrature_points_values          
     procedure, non_overridable :: get_value                            => fe_cell_function_tensor_get_value      
     procedure, non_overridable :: set_current_num_nodes             => fe_cell_function_tensor_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points => fe_cell_function_tensor_set_current_num_quadrature_points
     procedure, non_overridable :: free                                 => fe_cell_function_tensor_free
  end type fe_cell_function_tensor_t
  
 public :: fe_cell_function_scalar_t, fe_cell_function_vector_t, & 
           fe_cell_function_tensor_t, fe_cell_function_duties_t
 
contains
!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.
#include "sbm_fe_cell_function.i90"

#define duties fe_cell_function_duties
#define task_01 evaluate_values
#define task_02 evaluate_gradients
#define task_03 evaluate_laplacians
#include "duties_body.i90"

end module fe_cell_function_names
