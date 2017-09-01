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

#include "debug.i90"
#include "std_vector_macros.i90"
module std_vector_real_rp_names
  use types_names
  use memor_names
  implicit none
  private
  
  STD_VECTOR_TYPE(real,rp)
    
  public :: std_vector_real_rp_t
  
contains
  
  STD_VECTOR_PUSH_BACK(real,rp)
  STD_VECTOR_ERASE(real,rp)
  STD_VECTOR_RESIZE(real,rp)
  STD_VECTOR_SHRINK_TO_FIT(real,rp)
  STD_VECTOR_COPY(real,rp)
  STD_VECTOR_FREE(real,rp)
  STD_VECTOR_GET(real,rp)
  STD_VECTOR_GET_POINTER_TO_RANGE(real,rp)
  STD_VECTOR_SET(real,rp)
  STD_VECTOR_SIZE(real,rp)
  STD_VECTOR_CAPACITY(real,rp)
  STD_VECTOR_GET_RAW_POINTER(real,rp)

end module std_vector_real_rp_names
