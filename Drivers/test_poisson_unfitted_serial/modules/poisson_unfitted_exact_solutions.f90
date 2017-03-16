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

module poisson_unfitted_exact_solutions_names
  use fempar_names
  implicit none
# include "debug.i90"
  private

  ! TODO Exact solution can be a derived type from scalar_function_t
  type :: scalar_exact_solution_t
    contains
      procedure :: u      => scalar_exact_solution_u
      procedure :: grad_u => scalar_exact_solution_grad_u
      procedure :: lapl_u => scalar_exact_solution_lapl_u
  end type scalar_exact_solution_t

  type, extends(scalar_exact_solution_t) :: sol_ex001_2d_t
    contains
      procedure :: u      => sol_ex001_2d_u
      procedure :: grad_u => sol_ex001_2d_grad_u
      procedure :: lapl_u => sol_ex001_2d_lapl_u
  end type sol_ex001_2d_t

  type, extends(scalar_exact_solution_t) :: sol_ex002_2d_t
    real(rp) :: k  = 4.0*PI
    real(rp) :: x0 = 0.3
    real(rp) :: y0 = 0.3
    contains
      procedure :: u      => sol_ex002_2d_u
      procedure :: grad_u => sol_ex002_2d_grad_u
      procedure :: lapl_u => sol_ex002_2d_lapl_u
  end type sol_ex002_2d_t

  type, extends(scalar_exact_solution_t) :: sol_ex002_3d_t
    real(rp) :: k  = 4.0*PI
    real(rp) :: x0 = 0.3
    real(rp) :: y0 = 0.3
    real(rp) :: z0 = 0.3
    contains
      procedure :: u      => sol_ex002_3d_u
      procedure :: grad_u => sol_ex002_3d_grad_u
      procedure :: lapl_u => sol_ex002_3d_lapl_u
  end type sol_ex002_3d_t

  interface scalar_exact_solution_t

    ! This is the factory constructor.
    ! It creates the desired exact solution from a string identifier
    procedure scalar_exact_solution_constructor

  end interface scalar_exact_solution_t

  public :: scalar_exact_solution_t

  contains

    !==============================================================================================
    function scalar_exact_solution_constructor(example_name) result (exact_sol)

      character(len=*), intent(in)   :: example_name
      class(scalar_exact_solution_t), pointer :: exact_sol

      integer(ip) :: istat

      select case (example_name)
      case ("ex001_2d")
        allocate(sol_ex001_2d_t:: exact_sol, stat=istat); check(istat==0)
      case ("ex002_2d")
        allocate(sol_ex002_2d_t:: exact_sol, stat=istat); check(istat==0)
      case ("ex002_3d")
        allocate(sol_ex002_3d_t:: exact_sol, stat=istat); check(istat==0)
      case ("zero")
        allocate(scalar_exact_solution_t:: exact_sol, stat=istat); check(istat==0)
      case default
        check(.false.)
      end select

    end function scalar_exact_solution_constructor

    !==============================================================================================
    !==============================================================================================
    subroutine scalar_exact_solution_u(this,point,val)
      implicit none
      class(scalar_exact_solution_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      val = 0.0
    end subroutine scalar_exact_solution_u

    !==============================================================================================
    subroutine scalar_exact_solution_grad_u(this,point,val)
      implicit none
      class(scalar_exact_solution_t), intent(in) :: this
      type(point_t), intent(in) :: point
      type(vector_field_t), intent(inout) :: val
      call val%set(1,0.0)
      call val%set(2,0.0)
      call val%set(3,0.0)
    end subroutine scalar_exact_solution_grad_u

    !==============================================================================================
    subroutine scalar_exact_solution_lapl_u(this,point,val)
      implicit none
      class(scalar_exact_solution_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      val = 0.0
    end subroutine scalar_exact_solution_lapl_u

    !==============================================================================================
    !==============================================================================================
    subroutine sol_ex001_2d_u(this,point,val)
      implicit none
      class(sol_ex001_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2
      x1 = point%get(1)
      x2 = point%get(2)
      val =  x1**2*x2**2*(x1 - 1)**2*(x2 - 1)**2
    end subroutine sol_ex001_2d_u

    !==============================================================================================
    subroutine sol_ex001_2d_grad_u(this,point,val)
      implicit none
      class(sol_ex001_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      type(vector_field_t), intent(inout) :: val
      real(rp) :: x1, x2, g1, g2
      x1 = point%get(1)
      x2 = point%get(2)
      g1 = 2*x1*x2**2*(x1 - 1)**2*(x2 - 1)**2 + x1**2*x2**2*(2*x1 - 2)*(x2 - 1)**2
      g2 = 2*x1**2*x2*(x1 - 1)**2*(x2 - 1)**2 + x1**2*x2**2*(2*x2 - 2)*(x1 - 1)**2
      call val%set(1,g1)
      call val%set(2,g2)
      call val%set(3,0.0)
    end subroutine sol_ex001_2d_grad_u

    !==============================================================================================
    subroutine sol_ex001_2d_lapl_u(this,point,val)
      implicit none
      class(sol_ex001_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2
      x1 = point%get(1)
      x2 = point%get(2)
      val = 2*x1**2*x2**2*(x1 - 1)**2 +&
            2*x1**2*x2**2*(x2 - 1)**2 +&
            2*x1**2*(x1 - 1)**2*(x2 - 1)**2 +&
            2*x2**2*(x1 - 1)**2*(x2 - 1)**2 +&
            4*x1*x2**2*(2*x1 - 2)*(x2 - 1)**2 +&
            4*x1**2*x2*(2*x2 - 2)*(x1 - 1)**2
    end subroutine sol_ex001_2d_lapl_u

    !==============================================================================================
    !==============================================================================================
    subroutine sol_ex002_2d_u(this,point,val)
      implicit none
      class(sol_ex002_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      val = sin(this%k*(x1**2 + x2**2)**(1/2.0))
    end subroutine sol_ex002_2d_u

    !==============================================================================================
    subroutine sol_ex002_2d_grad_u(this,point,val)
      implicit none
      class(sol_ex002_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      type(vector_field_t), intent(inout) :: val
      real(rp) :: x1, x2, g1, g2
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      g1 = (this%k*x1*cos(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2)**(1/2.0)
      g2 = (this%k*x2*cos(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2)**(1/2.0)
      call val%set(1,g1)
      call val%set(2,g2)
      call val%set(3,0.0)
    end subroutine sol_ex002_2d_grad_u

    !==============================================================================================
    subroutine sol_ex002_2d_lapl_u(this,point,val)
      implicit none
      class(sol_ex002_2d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      val = (2*this%k*cos(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2)**(1/2.0) -&
            (this%k**2*x1**2*sin(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2) -&
            (this%k**2*x2**2*sin(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2) -&
            (this%k*x1**2*cos(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2)**(3/2.0) -&
            (this%k*x2**2*cos(this%k*(x1**2 + x2**2)**(1/2.0)))/(x1**2 + x2**2)**(3/2.0)
    end subroutine sol_ex002_2d_lapl_u

    !==============================================================================================
    !==============================================================================================
    subroutine sol_ex002_3d_u(this,point,val)
      implicit none
      class(sol_ex002_3d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2, x3
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      x3 = point%get(3) - this%z0
      val = sin(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0))
    end subroutine sol_ex002_3d_u

    !==============================================================================================
    subroutine sol_ex002_3d_grad_u(this,point,val)
      implicit none
      class(sol_ex002_3d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      type(vector_field_t), intent(inout) :: val
      real(rp) :: x1, x2, x3, g1, g2, g3
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      x3 = point%get(3) - this%z0
      g1 = (this%k*x1*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(1/2.0)
      g2 = (this%k*x2*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(1/2.0)
      g3 = (this%k*x3*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(1/2.0)
      call val%set(1,g1)
      call val%set(2,g2)
      call val%set(3,g3)
    end subroutine sol_ex002_3d_grad_u

    !==============================================================================================
    subroutine sol_ex002_3d_lapl_u(this,point,val)
      implicit none
      class(sol_ex002_3d_t), intent(in) :: this
      type(point_t), intent(in) :: point
      real(rp), intent(inout) :: val
      real(rp) :: x1, x2, x3
      x1 = point%get(1) - this%x0
      x2 = point%get(2) - this%y0
      x3 = point%get(3) - this%z0
      val = (3*this%k*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(1/2.0) -&
            (this%k*x1**2*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(3/2.0) -&
            (this%k*x2**2*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(3/2.0) -&
            (this%k*x3**2*cos(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)**(3/2.0) -&
            (this%k**2*x1**2*sin(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2) -&
            (this%k**2*x2**2*sin(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2) -&
            (this%k**2*x3**2*sin(this%k*(x1**2 + x2**2 + x3**2)**(1/2.0)))/(x1**2 + x2**2 + x3**2)
    end subroutine sol_ex002_3d_lapl_u

end module poisson_unfitted_exact_solutions_names
!***************************************************************************************************
