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

  !TODO. Add a mechanism to get the complement, the union etc.
  !Translations and complement can implemented in the base class

  character(len=*), parameter, public :: level_set_sphere_str       = 'sphere'
  character(len=*), parameter, public :: level_set_cylinder_str     = 'cylinder'
  character(len=*), parameter, public :: level_set_cheese_block_str = 'cheese_block'
  character(len=*), parameter, public :: level_set_popcorn_str      = 'popcorn'

  type :: level_set_function_factory_t
    private
  contains
    procedure, non_overridable, nopass :: create => level_set_function_factory_create
  end type level_set_function_factory_t

  type, extends(scalar_function_t) :: level_set_function_t
    private
    real(rp) :: tolerance = 1.0e-3_rp
<<<<<<< HEAD
    integer(ip) :: num_dims = SPACE_DIM
    real(rp) :: domain(6) = [0, 1, 0, 1, 0, 1]
=======
    integer(ip) :: num_dims = SPACE_DIM
>>>>>>> 893f8355af5bec5b333af960e1ff9a5e4e13b979
  contains
    procedure, private :: get_level_set_value => level_set_function_get_level_set_value
    procedure, non_overridable :: set_num_dims  => level_set_function_set_num_dims
    procedure :: get_value_space              => level_set_function_get_value_space
    procedure :: set_tolerance                => level_set_function_set_tolerance
    procedure :: set_domain                   => level_set_function_set_domain
    procedure :: get_tolerance                => level_set_function_get_tolerance
  end type level_set_function_t

  type, extends(level_set_function_t) :: level_set_sphere_t
    private
    real(rp) :: radius = 0.0_rp
  contains
    procedure :: set_radius                   => level_set_sphere_set_radius
    procedure, private :: get_level_set_value => level_set_sphere_get_level_set_value
  end type level_set_sphere_t

  type, extends(level_set_sphere_t) :: level_set_cylinder_t
    private
  contains
    procedure, private :: get_level_set_value => level_set_cylinder_get_level_set_value
  end type level_set_cylinder_t

  type, extends(level_set_function_t) :: level_set_cheese_block_t
    private
  contains
    procedure, private :: get_level_set_value => level_set_cheese_block_get_level_set_value
  end type level_set_cheese_block_t

  type, extends(level_set_function_t) :: level_set_popcorn_t
    private
  contains
    procedure, private :: get_level_set_value => level_set_popcorn_get_level_set_value
  end type level_set_popcorn_t

  public :: level_set_function_factory_t
  public :: level_set_function_t
  public :: level_set_sphere_t
  public :: level_set_cylinder_t
  public :: level_set_cheese_block_t
  public :: level_set_popcorn_t

contains

!========================================================================================
subroutine level_set_function_factory_create(level_set_str, level_set_function)

  implicit none
  character(len=*), intent(in) :: level_set_str
  class(level_set_function_t), allocatable, intent(inout) :: level_set_function

  integer(ip) :: istat

  select case (level_set_str)
    case (level_set_sphere_str)
      allocate( level_set_sphere_t::       level_set_function, stat= istat ); check(istat==0)
    case (level_set_cylinder_str)
      allocate( level_set_cylinder_t::     level_set_function, stat= istat ); check(istat==0)
    case (level_set_cheese_block_str)
      allocate( level_set_cheese_block_t:: level_set_function, stat= istat ); check(istat==0)
    case (level_set_popcorn_str)
      allocate( level_set_popcorn_t::      level_set_function, stat= istat ); check(istat==0)
    case default
      mcheck(.false., 'Unknown type of level set function `'//level_set_str//'`')
  end select

end subroutine level_set_function_factory_create

!========================================================================================
  subroutine level_set_function_get_level_set_value( this, point, result )
    implicit none
    class(level_set_function_t), intent(in)    :: this
    type(point_t)              , intent(in)    :: point
    real(rp)                   , intent(inout) :: result
    check(.false.) ! Derived types MUST implement this method.
  end subroutine level_set_function_get_level_set_value

!========================================================================================
  subroutine level_set_function_set_num_dims( this, num_dime )
    implicit none
    class(level_set_function_t), intent(inout) :: this
    integer(ip),                 intent(in)    :: num_dime
    this%num_dims = num_dime
  end subroutine level_set_function_set_num_dims

!========================================================================================
  subroutine level_set_function_get_value_space( this, point, result )
    implicit none
    class(level_set_function_t), intent(in)    :: this
    type(point_t)              , intent(in)    :: point
    real(rp)                   , intent(inout) :: result

    type(point_t) :: point2
    real(rp), parameter :: tol = 1.0e-14_rp

    ! Map to the new domain
    call point2%set(1, (this%domain(2)-this%domain(1))*point%get(1) + this%domain(1) )
    call point2%set(2, (this%domain(4)-this%domain(3))*point%get(2) + this%domain(3) )
    call point2%set(3, (this%domain(6)-this%domain(5))*point%get(3) + this%domain(5) )

    ! Restrict to the current dimensions
<<<<<<< HEAD
    if (this%num_dims < 3) call point2%set(3,0.0)
    if (this%num_dims < 2) call point2%set(2,0.0)
=======
    point2 = point
    if (this%num_dims < 3) call point2%set(3,0.0)
    if (this%num_dims < 2) call point2%set(2,0.0)
>>>>>>> 893f8355af5bec5b333af960e1ff9a5e4e13b979

    ! Call the actual level set function
    call this%get_level_set_value( point2, result )

    ! Apply tolerance
    if (abs(result) < tol) result = 0.0_rp

  end subroutine level_set_function_get_value_space

!========================================================================================
  subroutine level_set_function_set_tolerance ( this, tolerance_in )
    implicit none
    class(level_set_function_t), intent(inout) :: this
    real(rp)                   , intent(in)    :: tolerance_in
    this%tolerance = tolerance_in
  end subroutine level_set_function_set_tolerance

!========================================================================================
subroutine level_set_function_set_domain ( this, domain )
    implicit none
    class(level_set_function_t), intent(inout) :: this
    real(rp)                   , intent(in)    :: domain(6)
    this%domain(:) = domain
  end subroutine level_set_function_set_domain

!========================================================================================
  function level_set_function_get_tolerance ( this )
    implicit none
    class(level_set_function_t), intent(in) :: this
    real(rp) :: level_set_function_get_tolerance
    level_set_function_get_tolerance = this%tolerance
  end function level_set_function_get_tolerance

!========================================================================================
  subroutine level_set_sphere_set_radius ( this, radius_in )
    implicit none
    class(level_set_sphere_t), intent(inout) :: this
    real(rp)                    , intent(in)    :: radius_in
    this%radius = radius_in
  end subroutine level_set_sphere_set_radius

!========================================================================================
  subroutine level_set_sphere_get_level_set_value( this, point, result )
    implicit none
    class(level_set_sphere_t), intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    integer(ip), parameter :: x=1,y=2,z=3
    assert(this%radius > 0.0_rp)
    result = sqrt( point%get(x)**2 + point%get(y)**2 + point%get(z)**2 ) - this%radius
  end subroutine level_set_sphere_get_level_set_value

!========================================================================================
  subroutine level_set_cylinder_get_level_set_value( this, point, result )
    implicit none
    class(level_set_cylinder_t),  intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    integer(ip), parameter :: x=1,y=2
    assert(this%radius > 0.0_rp)
    result = sqrt( point%get(x)**2 + point%get(y)**2 ) - this%radius
  end subroutine level_set_cylinder_get_level_set_value

!========================================================================================
  subroutine level_set_cheese_block_get_level_set_value( this, point, result )
    implicit none
    class(level_set_cheese_block_t),  intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    real(rp) :: x, y, z
    x = 1.85*point%get(1)
    y = 1.85*point%get(2)
    z = 1.85*point%get(3)
    result = (x**2+y**2-4)**2 + (z**2-1.2)**2 + (y**2+z**2-4)**2 +&
             (x**2-1.2)**2 + (z**2+x**2-4)**2 + (y**2-1.2)**2 - 12
  end subroutine level_set_cheese_block_get_level_set_value

!========================================================================================
  subroutine level_set_popcorn_get_level_set_value( this, point, result )
    implicit none
    class(level_set_popcorn_t),  intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result

    real(rp) :: x, y, z
    real(rp) :: xk, yk, zk
    real(rp) :: r0, sg, A
    integer(ip) :: k

    x = point%get(1)
    y = point%get(2)
    z = point%get(3)

    r0 = 0.6
    sg = 0.2
    A  = 2.0
    result = sqrt(x**2 + y**2 + z**2) - r0

    do k = 0,11
        if (0 <= k .and. k <= 4) then
            xk = (r0/sqrt(5.0))*2.0*cos(2.0*k*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin(2.0*k*pi/5.0)
            zk = (r0/sqrt(5.0))
        else if (5 <= k .and. k <= 9) then
            xk = (r0/sqrt(5.0))*2.0*cos((2.0*(k-5)-1.0)*pi/5.0)
            yk = (r0/sqrt(5.0))*2.0*sin((2.0*(k-5)-1.0)*pi/5.0)
            zk =-(r0/sqrt(5.0))
        else if (k == 10) then
            xk = 0
            yk = 0
            zk = r0
        else
            xk = 0
            yk = 0
            zk = -r0
        end if
        result = result - A*exp( -( (x - xk)**2  + (y - yk)**2 + (z - zk)**2 )/(sg**2) )
    end do

  end subroutine level_set_popcorn_get_level_set_value

end module level_set_functions_gallery_names
!***************************************************************************************************
