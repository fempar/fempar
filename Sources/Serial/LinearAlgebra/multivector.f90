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
module multivector_names
  ! Serial modules
  use types_names
  use memor_names
  
  ! Abstract modules
  use vector_names
  use environment_names
  
  implicit none
# include "debug.i90"

  private
  
  ! multivector
  type :: multivector_t
     private 
     integer(ip) :: nvectors = -1
     class(vector_t), allocatable  :: vectors(:)
     class(environment_t), pointer :: environment  
   contains     
     procedure :: create    => multivector_create
     procedure :: get       => multivector_get
     procedure :: multidot  => multivector_multidot
     procedure :: multiaxpy => multivector_multiaxpy
     procedure :: free      => multivector_free
  end type multivector_t

  ! Types
  public :: multivector_t

contains

  subroutine multivector_create(this, environment, nvectors, mold)
    implicit none
    ! Parameters
    class(multivector_t)        , intent(inout) :: this
    class(environment_t), target, intent(in)    :: environment
    integer(ip)                 , intent(in)    :: nvectors
    class(vector_t)             , intent(in)    :: mold
    integer(ip) :: i
    
    this%nvectors = nvectors
    allocate (this%vectors(this%nvectors),mold=mold)
    do i=1, this%nvectors
      call this%vectors(i)%default_initialization()
    end do
    this%environment => environment
  end subroutine multivector_create
  
  subroutine multivector_free(this)
    implicit none
    ! Parameters
    class(multivector_t), intent(inout) :: this
    integer(ip) :: i
    do i=1, this%nvectors
      call this%vectors(i)%free()
    end do
    
    this%nvectors = -1
    deallocate (this%vectors)
    nullify(this%environment)
  end subroutine multivector_free
  
  function multivector_get(this,position)
    implicit none
    class(multivector_t), target, intent(in) :: this
    integer(ip)                 , intent(in) :: position
    class(vector_t), pointer :: multivector_get
    assert ( position >= 1 .and. position <= this%nvectors )
    assert ( allocated(this%vectors) )
    multivector_get => this%vectors(position)
  end function multivector_get
  
  ! s <- this_n^T * vector, with this_n = (this(1), this(2), ..., this(n))
  subroutine multivector_multidot (this,n,vector,s)
    implicit none
    class(multivector_t), intent(in)  :: this
    integer(ip)         , intent(in)  :: n
    class(vector_t)     , intent(in)  :: vector 
    real(rp)            , intent(out) :: s(n)
    integer(ip)                       :: j
    do j=1, n
      s(j) = vector%local_dot(this%vectors(j))
    end do
    call this%environment%first_level_sum(s)
  end subroutine multivector_multidot  
  
  ! vector <- vector +  alpha*this_n^T*s, with this_n = (this(1), this(2), ..., this(n))
  subroutine multivector_multiaxpy (this,n,vector,alpha,s)
    implicit none
    class(multivector_t), intent(in)     :: this
    integer(ip)         , intent(in)     :: n 
    class(vector_t)     , intent(inout)  :: vector
    real(rp)            , intent(in)     :: alpha
    real(rp)            , intent(in)     :: s(n)
    integer(ip)                          :: j
    do j=1, n
     call vector%axpby(alpha*s(j),this%vectors(j),1.0_rp)
    end do
  end subroutine multivector_multiaxpy
  
end module multivector_names
