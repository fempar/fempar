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
module face_fe_function_names
  ! Serial modules
  use types_names
  use list_types_names
  
  use reference_fe_names
  use fe_space_names
  use fe_function_names
  use cell_fe_function_names
  use field_names
  use environment_names
  
  implicit none
# include "debug.i90"
  private
 
  
  type face_fe_function_scalar_t
   private
   logical                         :: is_boundary
   type(i1p_t)                     :: quadrature_points_permutation(2)   
   type(cell_fe_function_scalar_t) :: cell_fe_function(2)
  contains
     procedure, non_overridable :: create                                => face_fe_function_scalar_create
     procedure, non_overridable :: update                                => face_fe_function_scalar_update
     procedure, non_overridable :: get_field_id                          => face_fe_function_scalar_get_field_id
     procedure, non_overridable :: get_nodal_values                      => face_fe_function_scalar_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values          => face_fe_function_scalar_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients       => face_fe_function_scalar_get_quadrature_points_gradients
     procedure, non_overridable :: get_value                             => face_fe_function_scalar_get_value
     procedure, non_overridable :: get_gradient                          => face_fe_function_scalar_get_gradient
     procedure, non_overridable :: set_current_number_nodes              => face_fe_function_scalar_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points  => face_fe_function_scalar_set_current_number_quadrature_points
     procedure, non_overridable :: free                                  => face_fe_function_scalar_free
  end type face_fe_function_scalar_t
  
  type face_fe_function_vector_t
   private
   logical                         :: is_boundary
   type(i1p_t)                     :: quadrature_points_permutation(2)  
   type(cell_fe_function_vector_t) :: cell_fe_function(2)
  contains
     procedure, non_overridable :: create                                => face_fe_function_vector_create
     procedure, non_overridable :: update                                => face_fe_function_vector_update
     procedure, non_overridable :: get_field_id                          => face_fe_function_vector_get_field_id
     procedure, non_overridable :: get_nodal_values                      => face_fe_function_vector_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values          => face_fe_function_vector_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients       => face_fe_function_vector_get_quadrature_points_gradients
     procedure, non_overridable :: get_value                             => face_fe_function_vector_get_value
     procedure, non_overridable :: get_gradient                          => face_fe_function_vector_get_gradient
     procedure, non_overridable :: set_current_number_nodes              => face_fe_function_vector_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points  => face_fe_function_vector_set_current_number_quadrature_points
     procedure, non_overridable :: free                                  => face_fe_function_vector_free
  end type face_fe_function_vector_t
  
  type face_fe_function_tensor_t
   private
   logical                         :: is_boundary
   type(i1p_t)                     :: quadrature_points_permutation(2)    
   type(cell_fe_function_tensor_t) :: cell_fe_function(2)
  contains
     procedure, non_overridable :: create                                => face_fe_function_tensor_create
     procedure, non_overridable :: update                                => face_fe_function_tensor_update
     procedure, non_overridable :: get_field_id                          => face_fe_function_tensor_get_field_id
     procedure, non_overridable :: get_nodal_values                      => face_fe_function_tensor_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values          => face_fe_function_tensor_get_quadrature_points_values
     procedure, non_overridable :: get_value                             => face_fe_function_tensor_get_value
     procedure, non_overridable :: set_current_number_nodes              => face_fe_function_tensor_set_current_number_nodes
     procedure, non_overridable :: set_current_number_quadrature_points  => face_fe_function_tensor_set_current_number_quadrature_points
     procedure, non_overridable :: free                                  => face_fe_function_tensor_free
  end type face_fe_function_tensor_t
  
 public :: fe_function_t
  
 public :: face_fe_function_scalar_t, face_fe_function_vector_t, face_fe_function_tensor_t
 
contains
!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.
#include "sbm_face_fe_function.i90"

end module face_fe_function_names
