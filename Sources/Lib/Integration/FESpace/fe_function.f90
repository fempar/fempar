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
module fe_function_names
  ! Serial modules
  use types_names
  use list_types_names
  
  use reference_fe_names
  use fe_space_names
  use function_names
  use environment_names
  use field_names
  
  ! Linear algebra
  use vector_names
  use serial_scalar_array_names
  
  implicit none
# include "debug.i90"
  private
  
  type fe_function_t
   private
   class(vector_t), allocatable  :: dof_values
   type(serial_scalar_array_t)   :: fixed_dof_values
  contains
     procedure, non_overridable          :: create                               => fe_function_create
     procedure, non_overridable          :: update_fixed_dof_values              => fe_function_update_fixed_dof_values
     procedure, non_overridable          :: gather_nodal_values_through_iterator => fe_function_gather_nodal_values_through_iterator
     procedure, non_overridable          :: gather_nodal_values_from_raw_data    => fe_function_gather_nodal_values_from_raw_data
     generic                             :: gather_nodal_values                  => gather_nodal_values_through_iterator, &
                                                                                    gather_nodal_values_from_raw_data
     procedure, non_overridable          :: scatter_nodal_values                 => fe_function_scatter_nodal_values
     procedure, private, non_overridable :: interpolate_scalar_function          => fe_function_interpolate_scalar_function
     procedure, private, non_overridable :: interpolate_vector_function          => fe_function_interpolate_vector_function
     procedure, private, non_overridable :: interpolate_tensor_function          => fe_function_interpolate_tensor_function
     generic                             :: interpolate_function                 => interpolate_scalar_function, &
                                                                                    interpolate_vector_function, &
                                                                                    interpolate_tensor_function
     procedure, non_overridable          :: copy                                 => fe_function_copy
     procedure, non_overridable          :: get_dof_values                       => fe_function_get_dof_values
     procedure, non_overridable          :: get_fixed_dof_values                 => fe_function_get_fixed_dof_values
     procedure, non_overridable          :: free                                 => fe_function_free
     generic                             :: assignment(=)                        => copy
  end type fe_function_t 
   
  public :: fe_function_t  
 
contains
! Includes with all the TBP and supporting subroutines for the types above.
! In a future, we would like to use the submodule features of FORTRAN 2008.
#include "sbm_fe_function.i90"

end module fe_function_names
