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

module level_set_functions_gallery_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  type, extends(scalar_function_t) :: level_set_function_t
    private
    real(rp) :: tolerance = 0.0_rp
  contains
    procedure, private :: get_level_set_value => level_set_function_get_level_set_value ! TODO: deferred better?
    procedure :: get_value_space              => level_set_function_get_value_space
    procedure :: set_tolerance                => level_set_function_set_tolerance
  end type level_set_function_t

  type, extends(level_set_function_t) :: level_set_cyllinder_t
    private
    real(rp) :: radius = 0.0_rp
  contains
    procedure :: set_radius                   => level_set_cyllinder_set_radius
    procedure, private :: get_level_set_value => level_set_cyllinder_get_level_set_value
  end type level_set_cyllinder_t

  public :: level_set_function_t
  public :: level_set_cyllinder_t

contains

  subroutine level_set_function_get_level_set_value( this, point, result )
    implicit none
    class(level_set_function_t), intent(in)    :: this
    type(point_t)              , intent(in)    :: point
    real(rp)                   , intent(inout) :: result
    check(.false.) ! Derived types MUST implement this method. TODO: deferred better?
  end subroutine level_set_function_get_level_set_value

  subroutine level_set_function_get_value_space( this, point, result )
    implicit none
    class(level_set_function_t), intent(in)    :: this
    type(point_t)              , intent(in)    :: point
    real(rp)                   , intent(inout) :: result

    call this%get_level_set_value( point, result )
    if (abs(result) < this%tolerance) then
      result = 0.0_rp
    end if

  end subroutine level_set_function_get_value_space

  subroutine level_set_function_set_tolerance ( this, tolerance_in )
    implicit none
    class(level_set_function_t), intent(inout) :: this
    real(rp)                   , intent(in)    :: tolerance_in
    this%tolerance = tolerance_in
  end subroutine level_set_function_set_tolerance

  !------------------------------------------------------------------

  subroutine level_set_cyllinder_set_radius ( this, radius_in )
    implicit none
    class(level_set_cyllinder_t), intent(inout) :: this
    real(rp)                    , intent(in)    :: radius_in
    this%radius = radius_in
  end subroutine level_set_cyllinder_set_radius


  subroutine level_set_cyllinder_get_level_set_value( this, point, result )
    implicit none
    class(level_set_cyllinder_t), intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result

    integer(ip), parameter :: x=1,y=2
    assert(this%radius > 0.0_rp)
    result = sqrt( point%get(x)**2 +  point%get(y)**2 ) - this%radius

  end subroutine level_set_cyllinder_get_level_set_value


end module level_set_functions_gallery_names
!***************************************************************************************************
