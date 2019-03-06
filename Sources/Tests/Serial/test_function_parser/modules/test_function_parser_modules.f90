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

module test_coded_function_names
  use fempar_names
  implicit none
# include "debug.i90"
  private
  type, extends(scalar_function_t) :: scalar_coded_function_add_1_2_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_add_1_2_ops_get_value_space
  end type scalar_coded_function_add_1_2_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_mul_1_2_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_mul_1_2_ops_get_value_space
  end type scalar_coded_function_mul_1_2_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_div_1_2_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_div_1_2_ops_get_value_space
  end type scalar_coded_function_div_1_2_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_pow_1_2_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_pow_1_2_ops_get_value_space
  end type scalar_coded_function_pow_1_2_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_add_3_5_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_add_3_5_ops_get_value_space
  end type scalar_coded_function_add_3_5_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_mul_3_5_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_mul_3_5_ops_get_value_space
  end type scalar_coded_function_mul_3_5_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_div_3_5_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_div_3_5_ops_get_value_space
  end type scalar_coded_function_div_3_5_ops_t

  type, extends(scalar_function_t) :: scalar_coded_function_pow_3_5_ops_t
    private 
   contains
     procedure :: get_value_space    => scalar_coded_function_pow_3_5_ops_get_value_space
  end type scalar_coded_function_pow_3_5_ops_t

  type, extends(vector_function_t) :: vector_coded_function_t
    private 
   contains
     procedure :: get_value_space    => vector_coded_function_get_value_space
  end type vector_coded_function_t

  type, extends(tensor_function_t) :: tensor_coded_function_t
    private 
   contains
     procedure :: get_value_space    => tensor_coded_function_get_value_space
  end type tensor_coded_function_t

  public :: scalar_coded_function_add_1_2_ops_t
  public :: scalar_coded_function_mul_1_2_ops_t
  public :: scalar_coded_function_div_1_2_ops_t
  public :: scalar_coded_function_pow_1_2_ops_t
  public :: scalar_coded_function_add_3_5_ops_t
  public :: scalar_coded_function_mul_3_5_ops_t
  public :: scalar_coded_function_div_3_5_ops_t
  public :: scalar_coded_function_pow_3_5_ops_t
  public :: vector_coded_function_t
  public :: tensor_coded_function_t
contains  

  !===============================================================================================
  subroutine scalar_coded_function_add_1_2_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_add_1_2_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x+y
      result = x+y 
    else if ( this%get_num_dims() == 3 ) then ! x+y+z
      z=point%get(3)
      result = x+y+z  
    end if
  end subroutine scalar_coded_function_add_1_2_ops_get_value_space

  subroutine scalar_coded_function_mul_1_2_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_mul_1_2_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x*y
      result = x*y 
    else if ( this%get_num_dims() == 3 ) then ! x*y*z
      z=point%get(3)
      result = x*y*z  
    end if
  end subroutine scalar_coded_function_mul_1_2_ops_get_value_space

  subroutine scalar_coded_function_div_1_2_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_div_1_2_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x/y
      result = x/y 
    else if ( this%get_num_dims() == 3 ) then ! x/y/z
      z=point%get(3)
      result = x/y/z  
    end if
  end subroutine scalar_coded_function_div_1_2_ops_get_value_space

  subroutine scalar_coded_function_pow_1_2_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_pow_1_2_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x^y
      result = x**y 
    else if ( this%get_num_dims() == 3 ) then ! x^y^z
      z=point%get(3)
      result = (x**y)**z  
    end if
  end subroutine scalar_coded_function_pow_1_2_ops_get_value_space

  !===============================================================================================

  subroutine scalar_coded_function_add_3_5_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_add_3_5_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x+y
      result = x+y+x+y 
    else if ( this%get_num_dims() == 3 ) then ! x+y+z
      z=point%get(3)
      result = x+y+z+x+y+z  
    end if
  end subroutine scalar_coded_function_add_3_5_ops_get_value_space

  subroutine scalar_coded_function_mul_3_5_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_mul_3_5_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x*y
      result = x*y*x*y 
    else if ( this%get_num_dims() == 3 ) then ! x*y*z
      z=point%get(3)
      result = x*y*z*x*y*z  
    end if
  end subroutine scalar_coded_function_mul_3_5_ops_get_value_space

  subroutine scalar_coded_function_div_3_5_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_div_3_5_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x/y
      result = x/y/x/y 
    else if ( this%get_num_dims() == 3 ) then ! x/y/z
      z=point%get(3)
      result = x/y/z/x/y/z  
    end if
  end subroutine scalar_coded_function_div_3_5_ops_get_value_space

  subroutine scalar_coded_function_pow_3_5_ops_get_value_space ( this, point, result )
    implicit none
    class(scalar_coded_function_pow_3_5_ops_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    real(rp),                       intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! x^y
      result = ((x**y)**x)**y 
    else if ( this%get_num_dims() == 3 ) then ! x^y^z
      z=point%get(3)
      result = ((((x**y)**z)**x)**y)**z  
    end if
  end subroutine scalar_coded_function_pow_3_5_ops_get_value_space

  !===============================================================================================
  subroutine vector_coded_function_get_value_space ( this, point, result )
    implicit none
    class(vector_coded_function_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    type(vector_field_t),           intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! (x+y, -x-y)
      call result%set(1, x+y); call result%set(2, -x-y) 
    else if ( this%get_num_dims() == 3 ) then ! (x+y+z, -x-y-z, x*y*z)
      z=point%get(3)
      call result%set(1, x+y+z); call result%set(2, -x-y-z); call result%set(3, x*y*z)  
    end if
  end subroutine vector_coded_function_get_value_space

  !===============================================================================================
  subroutine tensor_coded_function_get_value_space ( this, point, result )
    implicit none
    class(tensor_coded_function_t), intent(in)    :: this
    type(point_t),                  intent(in)    :: point
    type(tensor_field_t),           intent(inout) :: result
    real(rp)                                      :: x,y,z
    assert ( this%get_num_dims() == 2 .or. this%get_num_dims() == 3 )
    x=point%get(1); y=point%get(2)
    if ( this%get_num_dims() == 2 ) then  ! [(x+y, -x-y), (x*y, x/y)]
      call result%set(1,1, x+y); call result%set(1,2, -x-y)
      call result%set(2,1, x*y); call result%set(2,2, x/y) 
    else if ( this%get_num_dims() == 3 ) then ! [(x+y+z, -x-y-z, x*y*z), (x+y-z, -x-y+z, x/y/z), (x-y-z, -x+y-z, x^y^z)]
      z=point%get(3)
      call result%set(1,1, x+y+z); call result%set(1,2, -x-y-z); call result%set(1,3, x*y*z)
      call result%set(2,1, x+y-z); call result%set(2,2, -x-y+z); call result%set(2,3, x/y/z)
      call result%set(3,1, x-y-z); call result%set(3,2, -x+y-z); call result%set(3,3, (x**y)**z)
    end if
  end subroutine tensor_coded_function_get_value_space

end module test_coded_function_names
!***************************************************************************************************

