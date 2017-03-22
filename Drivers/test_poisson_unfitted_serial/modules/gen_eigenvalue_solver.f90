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

!****************************************************************************************
module gen_eigenvalue_solver_names
  use fempar_names
  use blas77_interfaces_names
  implicit none
# include "debug.i90"
  private

  type :: gen_eigenvalue_solver_base_t
    private
    integer(ip) :: N
    real(rp), allocatable :: A(:,:)
    real(rp), allocatable :: B(:,:)
    real(rp), allocatable :: lambdas(:,:)
    real(rp), allocatable :: scaling(:)
  contains
    procedure                           :: create         => gen_eigenvalue_solver_base_create
    procedure                           :: free           => gen_eigenvalue_solver_base_free
    procedure, non_overridable          :: solve          => gen_eigenvalue_solver_base_solve
    procedure, private                  :: solve_lapack   => gen_eigenvalue_solver_base_solve_lapack
    procedure, non_overridable, private :: scale_matrices => gen_eigenvalue_solver_base_scale_matrices
  end type gen_eigenvalue_solver_base_t


  type, extends(gen_eigenvalue_solver_base_t) :: gen_eigenvalue_solver_chol_t
    private
    real(rp), allocatable :: work(:)
    integer(ip) :: lwork
  contains
    procedure          :: create        => gen_eigenvalue_solver_chol_create
    procedure          :: free          => gen_eigenvalue_solver_chol_free
    procedure, private :: solve_lapack  => gen_eigenvalue_solver_chol_solve_lapack
  end type gen_eigenvalue_solver_chol_t

  type, extends(gen_eigenvalue_solver_base_t) :: gen_eigenvalue_solver_qz_t
    private
    real(rp), allocatable :: alphar(:)
    real(rp), allocatable :: alphai(:)
    real(rp), allocatable :: beta(:)
    real(rp), allocatable :: vl(:,:)
    real(rp), allocatable :: vr(:,:)
    integer(ip) :: ldvl
    integer(ip) :: ldvr
    real(rp), allocatable :: work(:)
    integer(ip) :: lwork
  contains
    procedure          :: create        => gen_eigenvalue_solver_qz_create
    procedure          :: free          => gen_eigenvalue_solver_qz_free
    procedure, private :: solve_lapack  => gen_eigenvalue_solver_qz_solve_lapack
  end type gen_eigenvalue_solver_qz_t

  ! This is the interface class. The one to be used by the clients
  type :: gen_eigenvalue_solver_t
    private
      type(gen_eigenvalue_solver_chol_t) :: chol
      type(gen_eigenvalue_solver_qz_t)   :: qz
    contains
      procedure, non_overridable :: create => gen_eigenvalue_solver_create
      procedure, non_overridable :: free   => gen_eigenvalue_solver_free
      procedure, non_overridable :: solve  => gen_eigenvalue_solver_solve
  end type gen_eigenvalue_solver_t


  public :: gen_eigenvalue_solver_base_t
  public :: gen_eigenvalue_solver_chol_t
  public :: gen_eigenvalue_solver_qz_t
  public :: gen_eigenvalue_solver_t

contains

!========================================================================================
!========================================================================================
subroutine gen_eigenvalue_solver_base_create(this,N)
  implicit none
  class(gen_eigenvalue_solver_base_t), intent(inout) :: this
  integer(ip), intent(in) :: N
  call this%free()
  this%N = N
  call memalloc (N,N,this%A,       __FILE__, __LINE__ )
  call memalloc (N,N,this%B,       __FILE__, __LINE__ )
  call memalloc (N,2,this%lambdas, __FILE__, __LINE__ )
  call memalloc (N,  this%scaling, __FILE__, __LINE__ )
end subroutine gen_eigenvalue_solver_base_create

!========================================================================================
subroutine gen_eigenvalue_solver_base_free(this)
  implicit none
  class(gen_eigenvalue_solver_base_t), intent(inout) :: this
  this%N = 0
  if (allocated(this%A))       call memfree (this%A,       __FILE__, __LINE__)
  if (allocated(this%B))       call memfree (this%B,       __FILE__, __LINE__)
  if (allocated(this%lambdas)) call memfree (this%lambdas, __FILE__, __LINE__)
  if (allocated(this%scaling)) call memfree (this%scaling, __FILE__, __LINE__)
end subroutine gen_eigenvalue_solver_base_free

!========================================================================================
function gen_eigenvalue_solver_base_solve(this,A,B,stat,scale) result (lambdas)
  implicit none
  class(gen_eigenvalue_solver_base_t), target, intent(inout) :: this

  ! Input
  real(rp), intent(in) :: A(:,:)
  real(rp), intent(in) :: B(:,:)
  logical, optional, intent(in) :: scale

  ! Output
  integer(ip), intent(inout) :: stat
  real(rp), pointer :: lambdas(:,:)

  logical :: wescale

  ! Determine if we want to scale
  ! We do not rescale by default
  wescale = .false.
  if( present(scale)) wescale = scale

  ! Reset the lambdas
  assert(allocated(this%lambdas))
  this%lambdas(:,:) = 0.0

  ! Set the matrices
  this%A(:,:) = A(:,:)
  this%B(:,:) = B(:,:)

  ! Rescale if desired
  if (wescale) call this%scale_matrices()

  ! Call the lapack as implemented in the derived classes
  stat = this%solve_lapack()

  ! Return the pointer
  lambdas =>this%lambdas

end function gen_eigenvalue_solver_base_solve

!========================================================================================
function gen_eigenvalue_solver_base_solve_lapack(this) result (stat)
  implicit none
  class(gen_eigenvalue_solver_base_t), target, intent(inout) :: this
  integer(ip) :: stat
  ! TODO @fverdugo DRIVER PRIORITY HIGH EFFORT LOW
  ! Why I cannot call a non deferred method of an abstract class?
  ! This function has to be implemented in the derived classes
  check(.false.)
end function gen_eigenvalue_solver_base_solve_lapack

!========================================================================================
subroutine gen_eigenvalue_solver_base_scale_matrices(this)
  implicit none
  class(gen_eigenvalue_solver_base_t), target, intent(inout) :: this

  integer(ip) :: i,j

  ! Compute the scaling factor
  do i=1,this%N
    this%scaling(i) = sqrt(abs(this%B(i,i)))
    if (this%scaling(i) == 0) this%scaling(i) = 1.0
    this%scaling(i) = 1.0/this%scaling(i)
  end do

  ! Scale the matrices (this will affect the eigenvectors, but not the eigenvalues)
  do i=1,this%N
    do j=1,this%N
      this%A(i,j) = this%scaling(i)*this%A(i,j)*this%scaling(j)
      this%B(i,j) = this%scaling(i)*this%B(i,j)*this%scaling(j)
    end do
  end do

end subroutine gen_eigenvalue_solver_base_scale_matrices

!========================================================================================
!========================================================================================
subroutine gen_eigenvalue_solver_chol_create(this,N)
  implicit none
  class(gen_eigenvalue_solver_chol_t), intent(inout) :: this
  integer(ip), intent(in) :: N
  call this%free()
  call this%gen_eigenvalue_solver_base_t%create(N)
  this%lwork = max(3*N-1,0)
  call memalloc (this%lwork,this%work, __FILE__, __LINE__ )
end subroutine gen_eigenvalue_solver_chol_create

!========================================================================================
subroutine gen_eigenvalue_solver_chol_free(this)
  implicit none
  class(gen_eigenvalue_solver_chol_t), intent(inout) :: this
  call this%gen_eigenvalue_solver_base_t%free()
  this%lwork = 0
  if (allocated(this%work)) call memfree (this%work, __FILE__, __LINE__)
end subroutine gen_eigenvalue_solver_chol_free

!========================================================================================
function gen_eigenvalue_solver_chol_solve_lapack(this) result (stat)
  implicit none
  class(gen_eigenvalue_solver_chol_t), target, intent(inout) :: this

  ! Output
  integer(ip) :: stat

  ! Call the Lapack routine
  ! TODO @fverdugo FEMPAT PRIORITY HIGH EFFORT LOW
  ! How to avoid implicit interface warning?
  ! How to know if we have to call the double or single precision version?
  call DSYGV(1,'N','L',this%N,this%A,this%N,this%B,this%N,this%lambdas(:,1),this%work,this%lwork,stat)

end function gen_eigenvalue_solver_chol_solve_lapack

!========================================================================================
!========================================================================================
subroutine gen_eigenvalue_solver_qz_create(this,N)
  implicit none
  class(gen_eigenvalue_solver_qz_t), intent(inout) :: this
  integer(ip), intent(in) :: N
  call this%free()
  call this%gen_eigenvalue_solver_base_t%create(N)
  this%lwork = max(8*N,1)
  this%ldvl  = 1
  this%ldvr  = 1
  call memalloc (N,this%alphar, __FILE__, __LINE__ )
  call memalloc (N,this%alphai, __FILE__, __LINE__ )
  call memalloc (N,this%beta  , __FILE__, __LINE__ )
  call memalloc (this%lwork,this%work, __FILE__, __LINE__ )
  call memalloc (this%ldvl,N,this%vl, __FILE__, __LINE__ )
  call memalloc (this%ldvr,N,this%vr, __FILE__, __LINE__ )
end subroutine gen_eigenvalue_solver_qz_create

!========================================================================================
subroutine gen_eigenvalue_solver_qz_free(this)
  implicit none
  class(gen_eigenvalue_solver_qz_t), intent(inout) :: this
  call this%gen_eigenvalue_solver_base_t%free()
  this%lwork = 0
  this%ldvl  = 0
  this%ldvr  = 0
  if (allocated(this%alphar)) call memfree (this%alphar, __FILE__, __LINE__)
  if (allocated(this%alphai)) call memfree (this%alphai, __FILE__, __LINE__)
  if (allocated(this%beta  )) call memfree (this%beta  , __FILE__, __LINE__)
  if (allocated(this%vl    )) call memfree (this%vl,   __FILE__, __LINE__)
  if (allocated(this%vr    )) call memfree (this%vr,   __FILE__, __LINE__)
  if (allocated(this%work  )) call memfree (this%work, __FILE__, __LINE__)
end subroutine gen_eigenvalue_solver_qz_free

!========================================================================================
function gen_eigenvalue_solver_qz_solve_lapack(this) result (stat)
  implicit none
  class(gen_eigenvalue_solver_qz_t), target, intent(inout) :: this

  ! Output
  integer(ip) :: stat

  integer(ip) :: i

  ! Call the Lapack routine
  ! TODO @fverdugo FEMPAT PRIORITY HIGH EFFORT LOW
  ! How to avoid implicit interface warning?
  ! How to know if we have to call the double or single precision version?
  call DGGEV('N','N',this%N,this%A,this%N,this%B,this%N,this%alphar,this%alphai,this%beta,&
  this%vl,this%ldvl,this%vr,this%ldvr,this%work,this%lwork,stat)

  ! Compute the lambdas
  ! TODO @fverdugo FEMPAT PRIORITY HIGH EFFORT LOW
  ! How to set to Inf when dividing by zero?
  ! And a NaN when 0/0?
  do i=1,this%N
    !assert( (this%beta(i) .ne. 0.0_rp) .or. (this%alphar(i) .ne. 0.0_rp) )
    if (this%beta(i) == 0) stat = this%N + 3
    this%lambdas(i,1) = this%alphar(i)/this%beta(i)
    this%lambdas(i,2) = this%alphai(i)/this%beta(i)
  end do

end function gen_eigenvalue_solver_qz_solve_lapack

!========================================================================================
!========================================================================================
subroutine gen_eigenvalue_solver_create(this,N)
  implicit none
  class(gen_eigenvalue_solver_t), intent(inout) :: this
  integer(ip), intent(in) :: N
  call this%free()
  call this%chol%create(N)
  call this%qz%create(N)
end subroutine gen_eigenvalue_solver_create

!========================================================================================
subroutine gen_eigenvalue_solver_free(this)
  implicit none
  class(gen_eigenvalue_solver_t), intent(inout) :: this
  call this%chol%free()
  call this%qz%free()
end subroutine gen_eigenvalue_solver_free

!========================================================================================
function gen_eigenvalue_solver_solve(this,A,B,stat) result (lambdas)
  implicit none
  class(gen_eigenvalue_solver_t), target, intent(inout) :: this

  ! Input
  real(rp), intent(in) :: A(:,:)
  real(rp), intent(in) :: B(:,:)

  ! Output
  integer(ip), intent(inout) :: stat
  real(rp), pointer :: lambdas(:,:)

  ! First try to solve with chol without scaling
  lambdas => this%chol%solve(A,B,stat)
  if (stat==0) return

  ! Then try to solve with chol with scaling
  lambdas => this%chol%solve(A,B,stat,.true.)
  if (stat==0) return

  ! Then try to solve with qz without scaling
  lambdas => this%chol%solve(A,B,stat)
  if (stat==0) return

  ! Finally try to solve with qz with scaling
  lambdas => this%qz%solve(A,B,stat,.true.)

end function gen_eigenvalue_solver_solve

!****************************************************************************************
end module gen_eigenvalue_solver_names
