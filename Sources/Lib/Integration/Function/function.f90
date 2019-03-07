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
module function_names
  use types_names
  use field_names

  implicit none
# include "debug.i90"

  private
  type :: array_function_t
     private
     integer(ip)  :: num_components = 1
     integer(ip)  :: num_dims = -1
   contains
     procedure, non_overridable :: set_num_components  => array_set_num_components
     procedure, non_overridable :: get_num_components  => array_function_get_num_components
     procedure                  :: set_num_dims        => array_function_set_num_dims
     procedure                  :: get_num_dims        => array_function_get_num_dims
     procedure  :: get_component_value_space              => array_function_get_component_value_space
     procedure  :: get_component_value_space_time         => array_function_get_component_value_space_time
     generic    :: get_component_value                    => get_component_value_space, get_component_value_space_time
     procedure  :: get_component_values_set_space         => array_function_get_component_values_set_space
     procedure  :: get_component_values_set_space_time    => array_function_get_component_values_set_space_time
     generic    :: get_component_values_set               => get_component_values_set_space, get_component_values_set_space_time
     procedure  :: get_value_space                        => array_function_get_value_space
     procedure  :: get_value_space_time                   => array_function_get_value_space_time
     generic    :: get_value                              => get_value_space, get_value_space_time
     procedure  :: get_values_set_space                   => array_function_get_values_set_space
     procedure  :: get_values_set_space_time              => array_function_get_values_set_space_time
     generic    :: get_values_set                         => get_values_set_space, get_values_set_space_time     
  end type array_function_t

  type :: scalar_function_t
     private
     integer(ip)  :: num_dims = -1
   contains
     procedure                  :: set_num_dims                 => scalar_function_set_num_dims
     procedure                  :: get_num_dims                 => scalar_function_get_num_dims
     procedure                  :: get_value_space              => scalar_function_get_value_space
     procedure                  :: get_value_space_time         => scalar_function_get_value_space_time
     generic                    :: get_value                    => get_value_space, get_value_space_time
     procedure                  :: get_values_set_space         => scalar_function_get_values_set_space
     procedure                  :: get_values_set_space_time    => scalar_function_get_values_set_space_time
     generic                    :: get_values_set               => get_values_set_space, get_values_set_space_time     

     procedure                  :: get_gradient_space           => scalar_function_get_gradient_space
     procedure                  :: get_gradient_space_time      => scalar_function_get_gradient_space_time
     generic                    :: get_gradient                 => get_gradient_space, get_gradient_space_time
     procedure                  :: get_gradients_set_space      => scalar_function_get_gradients_set_space
     procedure                  :: get_gradients_set_space_time => scalar_function_get_gradients_set_space_time
     generic                    :: get_gradients_set            => get_gradients_set_space, get_gradients_set_space_time     
     
  
     procedure                  :: get_value_temporal_derivative      => scalar_function_get_value_temporal_derivative
     procedure                  :: get_values_set_temporal_derivative => scalar_function_get_values_set_temporal_derivative
  end type scalar_function_t
  
  type :: p_scalar_function_t
    class(scalar_function_t), pointer :: p => NULL()
  end type p_scalar_function_t

  type :: vector_function_t
     private
     integer(ip)  :: num_dims = -1
   contains
     procedure                  :: set_num_dims                 => vector_function_set_num_dims
     procedure                  :: get_num_dims                 => vector_function_get_num_dims
     procedure                  :: get_value_space              => vector_function_get_value_space
     procedure                  :: get_value_space_time         => vector_function_get_value_space_time
     generic                    :: get_value                    => get_value_space, get_value_space_time
     procedure                  :: get_values_set_space         => vector_function_get_values_set_space
     procedure                  :: get_values_set_space_time    => vector_function_get_values_set_space_time
     generic                    :: get_values_set               => get_values_set_space, get_values_set_space_time  
     
     procedure                  :: get_gradient_space           => vector_function_get_gradient_space
     procedure                  :: get_gradient_space_time      => vector_function_get_gradient_space_time
     generic                    :: get_gradient                 => get_gradient_space, get_gradient_space_time
     procedure                  :: get_gradients_set_space      => vector_function_get_gradients_set_space
     procedure                  :: get_gradients_set_space_time => vector_function_get_gradients_set_space_time
     generic                    :: get_gradients_set            => get_gradients_set_space, get_gradients_set_space_time   
     
     procedure                  :: get_curl_space               => vector_function_get_curl_space
     procedure                  :: get_curl_space_time          => vector_function_get_curl_space_time
     generic                    :: get_curl                     => get_curl_space, get_curl_space_time
     procedure                  :: get_curls_set_space          => vector_function_get_curls_set_space
     procedure                  :: get_curls_set_space_time     => vector_function_get_curls_set_space_time
     generic                    :: get_curls_set                => get_curls_set_space, get_curls_set_space_time 
  end type vector_function_t

  type :: p_vector_function_t
    class(vector_function_t), pointer :: p => NULL()
  end type p_vector_function_t
  
  type :: tensor_function_t
     private
     integer(ip)  :: num_dims = -1
   contains
     procedure                  :: set_num_dims              => tensor_function_set_num_dims
     procedure                  :: get_num_dims              => tensor_function_get_num_dims
     procedure                  :: get_value_space           => tensor_function_get_value_space
     procedure                  :: get_value_space_time      => tensor_function_get_value_space_time
     generic                    :: get_value                 => get_value_space, get_value_space_time
     procedure                  :: get_values_set_space      => tensor_function_get_values_set_space
     procedure                  :: get_values_set_space_time => tensor_function_get_values_set_space_time
     generic                    :: get_values_set            => get_values_set_space, get_values_set_space_time     
  end type tensor_function_t

  type, extends(scalar_function_t) :: vector_component_function_t
     private
     class(vector_function_t), pointer :: vector_function => NULL()
     integer(ip)                       :: component = 0
   contains
     procedure                  :: set                          => vector_component_function_set
     procedure                  :: get_value_space              => vector_component_function_get_value_space
     procedure                  :: get_value_space_time         => vector_component_function_get_value_space_time
     procedure                  :: get_gradient_space           => vector_component_function_get_gradient_space
     procedure                  :: get_gradient_space_time      => vector_component_function_get_gradient_space_time
   end type vector_component_function_t

  ! Functions
  public :: array_function_t, scalar_function_t, vector_function_t, tensor_function_t, vector_component_function_t
  
  ! Pointer to functions
  public :: p_scalar_function_t, p_vector_function_t

contains
! One only needs to fill array_get_component_value_space or array_get_component_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root array_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from array_function_t
#include "sbm_array_function.i90"
! One only needs to fill scalar_get__value_space or scalar_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root scalar_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from scalar_function_t
#include "sbm_scalar_function.i90"
! One only needs to fill vector_get__value_space or vector_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root vector_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from vector_function_t
#include "sbm_vector_function.i90"
! One only needs to fill tensor_get__value_space or tensor_get_value_space_time
! (depending on the type of function). The rest of TBPs are implemented based on that at the 
! root tensor_function_t. For efficiency purposes, these functions can also be filled in user
! defined functions that inherit from tensor_function_t
#include "sbm_tensor_function.i90"
! One only needs to call set to assign the vector it gets the values from
#include "sbm_vector_component_function.i90"

end module function_names
