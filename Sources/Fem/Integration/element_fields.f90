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
# include "debug.i90"
module element_fields_names
  use types_names
  use memor_names
  use fem_space_names
  use memory_guard_names
  use problem_names
  implicit none
  private

  type, abstract, extends(memory_guard) :: field_t
     contains
       procedure(real_product_interface), deferred :: real_product
       procedure(field_binary_interface), deferred :: field_product
       procedure(field_binary_interface), deferred :: field_sum
       procedure(field_assign_interface), deferred :: field_assign
       procedure, pass(xfield) :: field_real_right => product_field_real
       procedure, pass(yfield) :: field_real_left  => product_real_field
       procedure :: free => free_field
       generic   :: operator(*)   => field_product, field_real_right, field_real_left
       generic   :: operator(+)   => field_sum
       generic   :: assignment(=) => field_assign
  end type field_t
  type, extends(field_t) :: scalar_t
     real(rp), allocatable :: a(:) ! To circumvent IBM compiler error we need a pointer...
     !real(rp), pointer     :: a(:)               ! a(ngaus)
     contains
       procedure :: field_product => product_scalar_field
       procedure :: real_product  => product_real_scalar
       procedure :: field_assign  => scalar_assign
       procedure :: field_sum     => sum_scalar_scalar
  end type scalar_t
  type , extends(field_t) :: vector_t
     real(rp), allocatable :: a(:,:)
     !real(rp), pointer :: a(:,:)   ! a(:,ngaus)
     contains
       procedure :: field_product => product_vector_field
       procedure :: real_product  => product_real_vector
       procedure :: field_assign  => vector_assign
       procedure :: field_sum     =>  sum_vector_vector
  end type vector_t
  type , extends(field_t) :: tensor_t
     real(rp), allocatable :: a(:,:,:)
     !real(rp), pointer :: a(:,:,:) ! a(:,:,ngaus)
     contains
       procedure :: field_product => product_tensor_field
       procedure :: real_product  => product_real_tensor
       procedure :: field_assign  => tensor_assign
       procedure :: field_sum     => sum_tensor_tensor
  end type tensor_t

  abstract interface
     subroutine real_product_interface (yfield,xreal,zfield)
       import :: field_t, rp
       real(rp)    , intent(in)    :: xreal
       class(field_t), intent(in)    :: yfield
       class(field_t), intent(inout) :: zfield
     end subroutine real_product_interface
     subroutine field_assign_interface(yfield,xfield) ! y=x
       import :: field_t
       class(field_t), intent(inout) :: yfield
       class(field_t), intent(in)    :: xfield
     end subroutine field_assign_interface
     function field_binary_interface (xfield,yfield) result(zfield)
       import :: field_t
       class(field_t), intent(in)  :: xfield,yfield
       class(field_t), allocatable :: zfield
     end function field_binary_interface
  end interface

  public :: field_t, scalar_t, vector_t, tensor_t
  public :: inv, norm

contains

  subroutine free_field ( this )
    class(field_t), intent(inout) :: this
  end subroutine free_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Field functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! General

  function product_real_field(xreal,yfield) result(zfield)
    implicit none
    real(rp)    , intent(in)  :: xreal
    class(field_t), intent(in)  :: yfield
    class(field_t), allocatable :: zfield
    allocate(zfield,mold=yfield)
    call yfield%real_product(xreal,zfield)
  end function product_real_field
  function product_field_real(xfield,yreal) result(zfield)
    implicit none
    real(rp)    , intent(in)  :: yreal
    class(field_t), intent(in)  :: xfield
    class(field_t), allocatable :: zfield
    allocate(zfield,mold=xfield)
    call xfield%real_product(yreal,zfield)
  end function product_field_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scalar

  subroutine scalar_assign(yfield,xfield) ! y=x
    implicit none
    class(scalar_t), intent(inout) :: yfield
    class(field_t) , intent(in)    :: xfield
    integer(ip)  :: ng
    call xfield%GuardTemp()
    if(allocated(yfield%a)) call memfree(yfield%a,__FILE__,__LINE__)
    select type(xfield)
    class is(scalar_t)
       ng = size(xfield%a,1)
       call memalloc(ng,yfield%a,__FILE__,__LINE__)
       yfield%a = xfield%a ! Intrinsic copy
    class default
       check(1==0)
    end select
    call xfield%CleanTemp()
  end subroutine scalar_assign

  function sum_scalar_scalar(xfield,yfield) result(zfield)
    implicit none
    class(scalar_t), intent(in)  :: xfield
    class(field_t) , intent(in)  :: yfield
    class(field_t) , allocatable :: zfield
    integer(ip)  :: ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    allocate(zfield,mold=xfield)
    select type(zfield)
    class is(scalar_t)
       ng = size(xfield%a,1)
       call memalloc(ng,zfield%a,__FILE__,__LINE__)
       select type(yfield)
       class is(scalar_t)
          assert(size(yfield%a,1)==ng)
          zfield%a = xfield%a + yfield%a
       class default
          check(1==0)
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function sum_scalar_scalar

  subroutine product_real_scalar(yfield,xreal,zfield)
    implicit none
    real(rp)     , intent(in)    :: xreal
    class(scalar_t), intent(in)    :: yfield
    class(field_t) , intent(inout) :: zfield
    integer(ip)  :: ng
    call yfield%GuardTemp()
    call zfield%SetTemp()
    ng = size(yfield%a,1)
    select type(zfield)
    class is(scalar_t)
       call memalloc(ng,zfield%a,__FILE__,__LINE__)
       zfield%a = xreal * yfield%a ! Intrinsic product
    class default
       check(1==0)
    end select
    call yfield%CleanTemp()
  end subroutine product_real_scalar

  function product_scalar_field(xfield,yfield) result(zfield)
    implicit none
    class(scalar_t), intent(in)  :: xfield     ! xfield(ngaus)
    class(field_t) , intent(in)  :: yfield     ! yfield
    class(field_t) , allocatable :: zfield     ! zfield
    integer(ip)   :: i,j,k,n1,n2,ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    ng = size(xfield%a,1)
    ! allocate(zfield, mold=yfield) ! Also possible, but less explicit
    select type(yfield)
    class is(scalar_t)
       assert(size(yfield%a,1)==ng)
       allocate(scalar_t :: zfield)
       select type(zfield)
       class is(scalar_t)
          call memalloc(ng,zfield%a,__FILE__,__LINE__)
          zfield%a = xfield%a * yfield%a
       end select
    class is(vector_t)
       n1 = size(yfield%a,1)
       assert(size(yfield%a,2)==ng)
       allocate(vector_t :: zfield)
       select type(zfield)
       class is(vector_t)
          call memalloc(n1,ng,zfield%a,__FILE__,__LINE__)
          do k=1,ng
             do i=1,n1
                zfield%a(i,k) = xfield%a(k) * yfield%a(i,k)
             end do
          end do
       end select
    class is(tensor_t)
       n1 = size(yfield%a,1)
       n2 = size(yfield%a,2)
       assert(size(yfield%a,3)==ng)
       allocate(tensor_t :: zfield)
       select type(zfield)
       class is(tensor_t)
          call memalloc(n1,n2,ng,zfield%a,__FILE__,__LINE__)
          do k=1,ng
             do j=1,n2
                do i=1,n1
                   zfield%a(i,j,k) = xfield%a(k) * yfield%a(i,j,k)
                end do
             end do
          end do
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function product_scalar_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector

  subroutine vector_assign(yfield,xfield) ! y=x
    implicit none
    class(vector_t), intent(inout) :: yfield
    class(field_t) , intent(in)    :: xfield
    integer(ip)  :: n1,ng
    call xfield%GuardTemp()
    if(allocated(yfield%a)) call memfree(yfield%a,__FILE__,__LINE__)
    select type(xfield)
    class is(vector_t)
       n1 = size(xfield%a,1)
       ng = size(xfield%a,2)
       call memalloc(n1,ng,yfield%a,__FILE__,__LINE__)
       yfield%a = xfield%a ! Intrinsic copy
    class default
       check(1==0)
    end select
    call xfield%CleanTemp()
  end subroutine vector_assign

  function sum_vector_vector(xfield,yfield) result(zfield)
    implicit none
    class(vector_t), intent(in)  :: xfield
    class(field_t) , intent(in)  :: yfield
    class(field_t), allocatable :: zfield
    integer(ip)   :: n1,ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    allocate(zfield,mold=xfield)
    select type(zfield)
    class is(vector_t)
       n1 = size(xfield%a,1)
       ng = size(xfield%a,2)
       call memalloc(n1,ng,zfield%a,__FILE__,__LINE__)
       select type(yfield)
       class is(vector_t)
          assert(size(yfield%a,1)==n1)
          assert(size(yfield%a,2)==ng)
          zfield%a = xfield%a + yfield%a
       class default
          check(1==0)
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function sum_vector_vector

  subroutine product_real_vector(yfield,xreal,zfield)
    implicit none
    real(rp)     , intent(in)    :: xreal
    class(vector_t), intent(in)    :: yfield
    class(field_t), intent(inout) :: zfield
    integer(ip)  :: n1,ng
    call yfield%GuardTemp()
    call zfield%SetTemp()
    n1 = size(yfield%a,1)
    ng = size(yfield%a,2)
    select type(zfield)
    class is(vector_t)
       call memalloc(n1,ng,zfield%a,__FILE__,__LINE__)
       zfield%a = xreal * yfield%a ! Intrinsic product
    class default
       check(1==0)
    end select
    call yfield%CleanTemp()
  end subroutine product_real_vector

  function product_vector_field(xfield,yfield) result(zfield)
    implicit none
    class(vector_t), intent(in)  :: xfield  ! xfield(ngaus)
    class(field_t) , intent(in)  :: yfield  ! yfield
    class(field_t) , allocatable :: zfield  ! zfield
    integer(ip)   :: i,j,k,n1,n2,ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    n1 = size(xfield%a,1)
    ng = size(xfield%a,2)
    select type(yfield)
    class is(scalar_t)
       assert(size(yfield%a,1)==ng)
       allocate(vector_t :: zfield)
       select type(zfield)
       class is(vector_t)
          call memalloc(n1,ng,zfield%a,__FILE__,__LINE__)
          do k=1,ng
             do i=1,n1
                zfield%a(i,k) = xfield%a(i,k) * yfield%a(k)
             end do
          end do
       end select
    class is(vector_t)
       assert(size(yfield%a,1)==n1)
       assert(size(yfield%a,2)==ng)
       allocate(scalar_t :: zfield)
       select type(zfield)
       class is(scalar_t)
          call memalloc(ng,zfield%a,__FILE__,__LINE__)
          ! Change by forall, BLAS, etc...
          do k=1,ng
             zfield%a(k) = 0.0_rp
             do i=1,n1
                zfield%a(k) = zfield%a(k) + xfield%a(i,k)*yfield%a(i,k)
             end do
          end do
       end select
    class is(tensor_t)
       n2=size(yfield%a,2)
       assert(size(yfield%a,1)==n1)
       assert(size(yfield%a,3)==ng)
       allocate(vector_t :: zfield)
       select type(zfield)
       class is(vector_t)
          call memalloc(n2,ng,zfield%a,__FILE__,__LINE__)
          ! Change by forall, BLAS, etc...
          do k=1,ng
             do j=1,n2
                zfield%a(j,k) = 0.0_rp
                do i=1,n1
                   zfield%a(j,k) = zfield%a(j,k) + xfield%a(i,k)*yfield%a(i,j,k)
                end do
             end do
          end do
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function product_vector_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Tensor

  subroutine tensor_assign(yfield,xfield) ! y=x
    implicit none
    class(tensor_t), intent(inout) :: yfield
    class(field_t) , intent(in)    :: xfield
    integer(ip)  :: n1,n2,ng
    call xfield%GuardTemp()
    if(allocated(yfield%a)) call memfree(yfield%a,__FILE__,__LINE__)
    select type(xfield)
    class is(tensor_t)
       n1 = size(xfield%a,1)
       n2 = size(xfield%a,2)
       ng = size(xfield%a,3)
       call memalloc(n1,n2,ng,yfield%a,__FILE__,__LINE__)
       yfield%a = xfield%a ! Intrinsic copy
    class default
       check(1==0)
    end select
    call xfield%CleanTemp()
  end subroutine tensor_assign

  function sum_tensor_tensor(xfield,yfield) result(zfield)
    implicit none
    class(tensor_t), intent(in)  :: xfield
    class(field_t) , intent(in)  :: yfield
    class(field_t), allocatable :: zfield
    integer(ip)   :: n1,n2,ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    allocate(zfield,mold=xfield)
    select type(zfield)
    class is(tensor_t)
       n1 = size(xfield%a,1)
       n2 = size(xfield%a,2)
       ng = size(xfield%a,3)
       call memalloc(n1,n2,ng,zfield%a,__FILE__,__LINE__)
       select type(yfield)
       class is(tensor_t)
          assert(size(yfield%a,1)==n1)
          assert(size(yfield%a,2)==n2)
          assert(size(yfield%a,3)==ng)
          zfield%a = xfield%a + yfield%a
       class default
          check(1==0)
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function sum_tensor_tensor

  subroutine product_real_tensor(yfield,xreal,zfield)
    implicit none
    real(rp)     , intent(in)    :: xreal
    class(tensor_t), intent(in)    :: yfield
    class(field_t) , intent(inout) :: zfield
    integer(ip)   :: n1,n2,ng
    call yfield%GuardTemp()
    call zfield%SetTemp()
    n1 = size(yfield%a,1)
    n2 = size(yfield%a,2)
    ng = size(yfield%a,3)
    select type(zfield)
    class is(tensor_t)
       call memalloc(n1,n2,ng,zfield%a,__FILE__,__LINE__)
       zfield%a = xreal * yfield%a ! Intrinsic product
    class default
       check(1==0)
    end select
    call yfield%CleanTemp()
  end subroutine product_real_tensor

  function product_tensor_field(xfield,yfield) result(zfield)
    implicit none
    class(tensor_t), intent(in)  :: xfield    ! xfield(ngaus)
    class(field_t) , intent(in)  :: yfield     ! yfield
    class(field_t) , allocatable :: zfield     ! zfield
    integer(ip)   :: i,j,k,n1,n2,ng
    call xfield%GuardTemp()
    call yfield%GuardTemp()
    n1=size(xfield%a,1)
    n2=size(xfield%a,2)
    ng=size(xfield%a,3)

    select type(yfield)
    class is(scalar_t)
       assert(size(yfield%a,1)==ng)
       allocate(tensor_t :: zfield)
       select type(zfield)
       class is(tensor_t)
          call memalloc(n1,n2,ng,zfield%a,__FILE__,__LINE__)
          ! Change by forall, BLAS, etc...
          do k=1,ng
             do j=1,n2
                do i=1,n1
                   zfield%a(i,j,k) = xfield%a(i,j,k)*yfield%a(k)
                end do
             end do
          end do
       end select
    class is(vector_t)
       assert(size(yfield%a,1)==n2)
       assert(size(yfield%a,2)==ng)
       allocate(vector_t :: zfield)
       select type(zfield)
       class is(vector_t)
          call memalloc(n1,ng,zfield%a,__FILE__,__LINE__)
          ! Change by forall, BLAS, etc...
          do k=1,ng
             zfield%a(:,k) = 0.0_rp
             do j=1,n2
                do i=1,n1
                   zfield%a(i,k) = zfield%a(i,k) + xfield%a(i,j,k)*yfield%a(j,k)
                end do
             end do
          end do
       end select
    class is(tensor_t)
       assert(size(yfield%a,1)==n1)
       assert(size(yfield%a,2)==n2)
       assert(size(yfield%a,3)==ng)
       allocate(scalar_t :: zfield)
       select type(zfield)
       class is(scalar_t)
          call memalloc(ng,zfield%a,__FILE__,__LINE__)
          ! Change by forall, BLAS, etc...
          do k=1,ng
             zfield%a(k) = 0.0_rp
             do j=1,n2
                do i=1,n1
                   zfield%a(k) = zfield%a(k) + xfield%a(i,j,k)*yfield%a(i,j,k)
                end do
             end do
          end do
       end select
    class default
       check(1==0)
    end select
    call zfield%SetTemp()
    call xfield%CleanTemp()
    call yfield%CleanTemp()
  end function product_tensor_field

! term by term inverse
  function inv(x) result(z)
    implicit none
    !type(scalar_t) :: x    ! x(ngaus)
    class(field_t) :: x    ! x(ngaus)
    type(scalar_t) :: z    ! z(ngaus)
    integer(ip)  :: i,ng
    call x%GuardTemp()
    select type(x)
    class is(scalar_t)
       call z%SetTemp()
       ng = size(x%a,1)
       call memalloc(ng,z%a,__FILE__,__LINE__)
       forall(i=1:ng)
          z%a(i) = 1.0_rp/x%a(i)
       end forall
    class default
       check(1==0)
    end select
    call x%CleanTemp()
  end function inv

  function norm(x) result(z)
    implicit none
    type(vector_t) :: x    ! x(ndime,ngaus)
    type(scalar_t) :: z    ! z(ngaus)
    integer(ip)  :: i,j,k
    do i=1,size(x%a,2)
       z%a(i) = 0.0_rp
       do j=1,size(x%a,1)
          z%a(i) = z%a(i) + x%a(j,i)
       end do
    end do
  end function norm


end module element_fields_names
