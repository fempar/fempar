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
  use types_names 
  use field_names
  use function_names
  implicit none
# include "debug.i90"
  private

  !TODO. Add a mechanism to get the complement, the union etc.
  !Translations and complement can implemented in the base class

  character(len=*), parameter, public :: level_set_sphere_str       = 'sphere'
  character(len=*), parameter, public :: level_set_cylinder_str     = 'cylinder'
  character(len=*), parameter, public :: level_set_cheese_block_str = 'cheese_block'
  character(len=*), parameter, public :: level_set_popcorn_str      = 'popcorn'
  character(len=*), parameter, public :: level_set_bullets_str      = 'bullets'
  character(len=*), parameter, public :: level_set_spiral_str       = 'spiral'

  integer(ip), parameter, private :: NUM_SPHERES = 20

  type :: level_set_function_factory_t
    private
  contains
    procedure, non_overridable, nopass :: create => level_set_function_factory_create
  end type level_set_function_factory_t

  type, extends(scalar_function_t) :: level_set_function_t
    private
    real(rp) :: tolerance = 1.0e-3_rp
    integer(ip) :: num_dims = SPACE_DIM
    real(rp) :: domain(6) = [0, 1, 0, 1, 0, 1]
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
    type(point_t) :: center
  contains
    procedure :: set_radius                   => level_set_sphere_set_radius
    procedure :: set_center                   => level_set_sphere_set_center
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

  type, extends(level_set_function_t) :: level_set_bullets_t
    private
    type(level_set_sphere_t) ::  spheres(NUM_SPHERES)
  contains
    procedure          :: init => level_set_bullets_init
    procedure, private :: get_level_set_value => level_set_bullets_get_level_set_value
  end type level_set_bullets_t

  type, extends(level_set_function_t) :: level_set_spiral_t
    private
  contains
    procedure, private :: get_level_set_value => level_set_spiral_get_level_set_value
  end type level_set_spiral_t

  public :: level_set_function_factory_t
  public :: level_set_function_t
  public :: level_set_sphere_t
  public :: level_set_cylinder_t
  public :: level_set_cheese_block_t
  public :: level_set_popcorn_t
  public :: level_set_bullets_t
  public :: level_set_spiral_t

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
    case (level_set_bullets_str)
      allocate( level_set_bullets_t::      level_set_function, stat= istat ); check(istat==0)
      select type (level_set_function)
        type is (level_set_bullets_t)
          call level_set_function%init()
        class default
          check(.false.)
      end select
    case (level_set_spiral_str)
      allocate( level_set_spiral_t::      level_set_function, stat= istat ); check(istat==0)
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
    if (this%num_dims < 3) call point2%set(3,0.0)
    if (this%num_dims < 2) call point2%set(2,0.0)

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
  subroutine level_set_sphere_set_center ( this, center_in )
    implicit none
    class(level_set_sphere_t), intent(inout) :: this
    real(rp)                 , intent(in)    :: center_in(3)
    call this%center%set(1,center_in(1))
    call this%center%set(2,center_in(2))
    call this%center%set(3,center_in(3))
  end subroutine level_set_sphere_set_center

!========================================================================================
  subroutine level_set_sphere_get_level_set_value( this, point, result )
    implicit none
    class(level_set_sphere_t), intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    integer(ip), parameter :: x=1,y=2,z=3
    assert(this%radius > 0.0_rp)
    result = sqrt( (point%get(x)-this%center%get(x))**2 + (point%get(y)-this%center%get(y))**2 + (point%get(z)-this%center%get(z))**2 ) - this%radius
  end subroutine level_set_sphere_get_level_set_value

!========================================================================================
  subroutine level_set_cylinder_get_level_set_value( this, point, result )
    implicit none
    class(level_set_cylinder_t),  intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result
    integer(ip), parameter :: x=1,y=2
    assert(this%radius > 0.0_rp)
    result = sqrt( (point%get(x)-this%center%get(x))**2 + (point%get(y)-this%center%get(y))**2 ) - this%radius
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

!========================================================================================
  subroutine level_set_bullets_init( this )
    implicit none
    class(level_set_bullets_t)  , intent(inout)    :: this
      call this%spheres( 1)%set_radius(0.0735)
      call this%spheres( 2)%set_radius(0.0778)
      call this%spheres( 3)%set_radius(0.0866)
      call this%spheres( 4)%set_radius(0.0861)
      call this%spheres( 5)%set_radius(0.0712)
      call this%spheres( 6)%set_radius(0.0780)
      call this%spheres( 7)%set_radius(0.0805)
      call this%spheres( 8)%set_radius(0.0783)
      call this%spheres( 9)%set_radius(0.0831)
      call this%spheres(10)%set_radius(0.0826)
      call this%spheres(11)%set_radius(0.0758)
      call this%spheres(12)%set_radius(0.0786)
      call this%spheres(13)%set_radius(0.0703)
      call this%spheres(14)%set_radius(0.0897)
      call this%spheres(15)%set_radius(0.0733)
      call this%spheres(16)%set_radius(0.0721)
      call this%spheres(17)%set_radius(0.0774)
      call this%spheres(18)%set_radius(0.0740)
      call this%spheres(19)%set_radius(0.0798)
      call this%spheres(20)%set_radius(0.0768)

      call this%spheres( 1)%set_center([0.8541, 0.1336, 0.3054])
      call this%spheres( 2)%set_center([0.8295, 0.5480, 0.4257])!
      call this%spheres( 3)%set_center([0.1493, 0.7994, 0.5918])
      call this%spheres( 4)%set_center([0.6865, 0.7326, 0.3160])!
      call this%spheres( 5)%set_center([0.3190, 0.2573, 0.7543])
      call this%spheres( 6)%set_center([0.4395, 0.3972, 0.8784])
      call this%spheres( 7)%set_center([0.5375, 0.6092, 0.6805])!
      call this%spheres( 8)%set_center([0.8471, 0.8776, 0.3776])
      call this%spheres( 9)%set_center([0.4355, 0.2306, 0.5659])
      call this%spheres(10)%set_center([0.8787, 0.7787, 0.1925])
      call this%spheres(11)%set_center([0.3443, 0.6135, 0.8185])
      call this%spheres(12)%set_center([0.6577, 0.4030, 0.7976])
      call this%spheres(13)%set_center([0.6304, 0.1800, 0.7491])!
      call this%spheres(14)%set_center([0.5307, 0.4438, 0.3124])
      call this%spheres(15)%set_center([0.6553, 0.4859, 0.5740])
      call this%spheres(16)%set_center([0.6306, 0.2026, 0.1256])
      call this%spheres(17)%set_center([0.2477, 0.5702, 0.4414])
      call this%spheres(18)%set_center([0.2084, 0.2853, 0.3532])
      call this%spheres(19)%set_center([0.8913, 0.4095, 0.2346])
      call this%spheres(20)%set_center([0.2422, 0.5651, 0.2482])
  end subroutine level_set_bullets_init

!========================================================================================
  subroutine level_set_bullets_get_level_set_value( this, point, result )
    implicit none
    class(level_set_bullets_t)  , intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result

    integer(ip) :: s
    real(rp) :: sval

    result = 1.0e30
    do s=1, NUM_SPHERES
      call this%spheres(s)%get_value(point,sval)
      result = min(result, sval)
    end do

    result = -1.0*result

  end subroutine level_set_bullets_get_level_set_value

!========================================================================================
  subroutine level_set_spiral_get_level_set_value( this, point, result)

    implicit none
    class(level_set_spiral_t)   , intent(in)    :: this
    type(point_t)               , intent(in)    :: point
    real(rp)                    , intent(inout) :: result

    real(rp) :: c1, c2
    real(rp) :: x1, x2, x3

    integer(ip), parameter :: n = 2
    real(rp)   , parameter :: r = 0.15
    real(rp)   , parameter :: rd = 0.8

    x1 = point%get(1)
    x2 = point%get(2)
    x3 = point%get(3)

    c1 = rd * sin(2.0*PI*x3/real(n))
    c2 = rd * cos(2.0*PI*x3/real(n))

    result = (x1-c1)**2 + (x2-c2)**2 - r

  end subroutine level_set_spiral_get_level_set_value

end module level_set_functions_gallery_names
!***************************************************************************************************
