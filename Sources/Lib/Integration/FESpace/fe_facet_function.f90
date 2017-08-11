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
module fe_facet_function_names
  ! Serial modules
  use types_names
  use list_types_names
  
  use reference_fe_names
  use fe_space_names
  use fe_cell_function_names
  use field_names
  use environment_names
  
  implicit none
# include "debug.i90"
  private
  
  
  type fe_facet_function_scalar_t
   private
   logical                                      :: is_at_boundary
   integer(ip)                                  :: active_cell_id(2)
   type(i1p_t)                                  :: quadrature_points_permutation(2)
   type(fe_cell_function_scalar_t)              :: fe_cell_function(2)
   class(fe_cell_iterator_t)      , allocatable :: fe
   class(serial_fe_space_t)       , pointer     :: fe_space => NULL()
  contains
     procedure, non_overridable :: create                             => fe_facet_function_scalar_create
     procedure, non_overridable :: update                             => fe_facet_function_scalar_update
     procedure, non_overridable :: get_field_id                       => fe_facet_function_scalar_get_field_id
     procedure, non_overridable :: get_nodal_values                   => fe_facet_function_scalar_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values       => fe_facet_function_scalar_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients    => fe_facet_function_scalar_get_quadrature_points_gradients
     procedure, non_overridable :: get_value                          => fe_facet_function_scalar_get_value
     procedure, non_overridable :: get_gradient                       => fe_facet_function_scalar_get_gradient
     procedure, non_overridable :: set_current_num_nodes              => fe_facet_function_scalar_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points  => fe_facet_function_scalar_set_current_num_quadrature_points
     procedure, non_overridable :: free                               => fe_facet_function_scalar_free
  end type fe_facet_function_scalar_t
  
  type fe_facet_function_vector_t
   private
   logical                                      :: is_at_boundary
   integer(ip)                                  :: active_cell_id(2)
   type(i1p_t)                                  :: quadrature_points_permutation(2)
   type(fe_cell_function_vector_t)              :: fe_cell_function(2)
   class(fe_cell_iterator_t)      , allocatable :: fe
   class(serial_fe_space_t)       , pointer     :: fe_space => NULL()
  contains
     procedure, non_overridable :: create                             => fe_facet_function_vector_create
     procedure, non_overridable :: update                             => fe_facet_function_vector_update
     procedure, non_overridable :: get_field_id                       => fe_facet_function_vector_get_field_id
     procedure, non_overridable :: get_nodal_values                   => fe_facet_function_vector_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values       => fe_facet_function_vector_get_quadrature_points_values
     procedure, non_overridable :: get_quadrature_points_gradients    => fe_facet_function_vector_get_quadrature_points_gradients
     procedure, non_overridable :: get_value                          => fe_facet_function_vector_get_value
     procedure, non_overridable :: get_gradient                       => fe_facet_function_vector_get_gradient
     procedure, non_overridable :: set_current_num_nodes              => fe_facet_function_vector_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points  => fe_facet_function_vector_set_current_num_quadrature_points
     procedure, non_overridable :: free                               => fe_facet_function_vector_free
  end type fe_facet_function_vector_t
  
  type fe_facet_function_tensor_t
   private
   logical                                      :: is_at_boundary
   integer(ip)                                  :: active_cell_id(2)
   type(i1p_t)                                  :: quadrature_points_permutation(2)
   type(fe_cell_function_tensor_t)              :: fe_cell_function(2)
   class(fe_cell_iterator_t)      , allocatable :: fe
   class(serial_fe_space_t)       , pointer     :: fe_space => NULL()
  contains
     procedure, non_overridable :: create                            => fe_facet_function_tensor_create
     procedure, non_overridable :: update                            => fe_facet_function_tensor_update
     procedure, non_overridable :: get_field_id                      => fe_facet_function_tensor_get_field_id
     procedure, non_overridable :: get_nodal_values                  => fe_facet_function_tensor_get_nodal_values
     procedure, non_overridable :: get_quadrature_points_values      => fe_facet_function_tensor_get_quadrature_points_values
     procedure, non_overridable :: get_value                         => fe_facet_function_tensor_get_value
     procedure, non_overridable :: set_current_num_nodes             => fe_facet_function_tensor_set_current_num_nodes
     procedure, non_overridable :: set_current_num_quadrature_points => fe_facet_function_tensor_set_current_num_quadrature_points
     procedure, non_overridable :: free                              => fe_facet_function_tensor_free
  end type fe_facet_function_tensor_t
  
 public :: fe_facet_function_scalar_t, fe_facet_function_vector_t, fe_facet_function_tensor_t
 
contains
!  ! Includes with all the TBP and supporting subroutines for the types above.
!  ! In a future, we would like to use the submodule features of FORTRAN 2008.
#include "sbm_fe_facet_function.i90"

end module fe_facet_function_names
