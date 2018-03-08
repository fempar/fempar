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
module mixed_laplacian_rt_conditions_names
  use fempar_names
  
  implicit none
# include "debug.i90"
  private
  type, extends(conditions_t) :: mixed_laplacian_rt_conditions_t
     private
     integer(ip) :: num_dims
   contains
     procedure :: set_num_dims          => mixed_laplacian_rt_conditions_set_num_dims
     procedure :: get_num_components       => mixed_laplacian_rt_conditions_get_num_components  
     procedure :: get_components_code         => mixed_laplacian_rt_conditions_get_components_code
     procedure :: get_function                => mixed_laplacian_rt_conditions_get_function
  end type mixed_laplacian_rt_conditions_t
  
  public :: mixed_laplacian_rt_conditions_t
  
contains
  
  subroutine mixed_laplacian_rt_conditions_set_num_dims (this, num_dims)
    implicit none
    class(mixed_laplacian_rt_conditions_t), intent(inout) :: this
    integer(ip)                           , intent(in)    :: num_dims
    this%num_dims = num_dims
  end subroutine mixed_laplacian_rt_conditions_set_num_dims 

  function mixed_laplacian_rt_conditions_get_num_components(this)
    implicit none
    class(mixed_laplacian_rt_conditions_t), intent(in) :: this
    integer(ip) :: mixed_laplacian_rt_conditions_get_num_components
    assert ( this%num_dims == 2 .or. this%num_dims == 3 ) 
    mixed_laplacian_rt_conditions_get_num_components = this%num_dims+1
  end function mixed_laplacian_rt_conditions_get_num_components

  subroutine mixed_laplacian_rt_conditions_get_components_code(this, boundary_id, components_code)
    implicit none
    class(mixed_laplacian_rt_conditions_t), intent(in)  :: this
    integer(ip)                       , intent(in)  :: boundary_id
    logical                           , intent(out) :: components_code(:)
    assert ( size(components_code) == 3 .or. size(components_code) == 4 )
    components_code(1:size(components_code)) = .false.
  end subroutine mixed_laplacian_rt_conditions_get_components_code
  
  subroutine mixed_laplacian_rt_conditions_get_function ( this, boundary_id, component_id, function )
    implicit none
    class(mixed_laplacian_rt_conditions_t), target     , intent(in)  :: this
    integer(ip)                                    , intent(in)  :: boundary_id
    integer(ip)                                    , intent(in)  :: component_id
    class(scalar_function_t)          , pointer    , intent(out) :: function
    assert ( component_id == 1 .or. component_id == 2 .or. component_id == 3 .or. component_id == 4 )
    nullify(function)
  end subroutine mixed_laplacian_rt_conditions_get_function 

end module mixed_laplacian_rt_conditions_names
