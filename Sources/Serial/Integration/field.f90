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
module field_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  type :: scalar_field_t
     private
     real(rp) :: value
   contains
     procedure, non_overridable :: init  => scalar_field_init
     procedure, non_overridable :: set   => scalar_field_set
  end type scalar_field_t
  
  type :: vector_field_t
     private
     real(rp) :: value(number_space_dimensions)
   contains
     procedure, non_overridable :: init  => vector_field_init
     procedure, non_overridable :: set   => vector_field_set
     procedure, non_overridable :: get   => vector_field_get	
     procedure, non_overridable :: add   => vector_field_add
     procedure, non_overridable :: nrm2  => vector_field_nrm2
  end type vector_field_t
  
  type :: tensor_field_t
     private
     real(rp)  :: value(number_space_dimensions,number_space_dimensions)
   contains
     procedure, non_overridable :: init  => tensor_field_init
     procedure, non_overridable :: set   => tensor_field_set
     procedure, non_overridable :: get   => tensor_field_get
     procedure, non_overridable :: add   => tensor_field_add
  end type tensor_field_t

  type :: symmetric_tensor_field_t
     private
     real(rp)  :: value(number_space_dimensions,number_space_dimensions)
   contains			
     procedure, non_overridable :: init  => symmetric_tensor_field_init
     procedure, non_overridable :: set   => symmetric_tensor_field_set					
  end type symmetric_tensor_field_t
  
  interface operator(*)
    module procedure single_contract_vector_vector, single_contract_tensor_vector, &
                     single_contract_vector_tensor, single_contract_tensor_tensor
    module procedure scal_left_vector, scal_right_vector, scal_left_tensor, scal_right_tensor
    module procedure scal_left_scalar, scal_right_scalar
  end interface
  
  interface double_contract
    module procedure double_contract_tensor_tensor
  end interface
  
  !public :: scalar_field_t (not actually needed, used real(rp) instead)
  public :: vector_field_t, tensor_field_t, symmetric_tensor_field_t 
  public :: operator(*)
  public :: double_contract
contains

  subroutine scalar_field_init(this,value)
    implicit none
    class(scalar_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine scalar_field_init

  subroutine scalar_field_set(this,value)
    implicit none
    class(scalar_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine scalar_field_set
		
  subroutine vector_field_init(this,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine vector_field_init

  subroutine vector_field_set(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = value
  end subroutine vector_field_set
  
  function vector_field_get(this,i) result(value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)                             :: value
    value = this%value(i)
  end function vector_field_get
  
  subroutine vector_field_add(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = this%value(i) + value
  end subroutine vector_field_add
  
  function vector_field_nrm2(this)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp) :: vector_field_nrm2
    vector_field_nrm2 = this * this
    vector_field_nrm2 = sqrt(vector_field_nrm2)
  end function vector_field_nrm2
		
  subroutine tensor_field_init(this,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine tensor_field_init

  subroutine tensor_field_set(this,i,j,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    integer(ip)          , intent(in)    :: j
    real(rp)             , intent(in)    :: value
    this%value(i,j) = value
  end subroutine tensor_field_set
  
  function tensor_field_get(this,i,j) result(value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    integer(ip)          , intent(in)    :: j
    real(rp)                             :: value
    value = this%value(i,j)
  end function tensor_field_get
  
  subroutine tensor_field_add(this,i,j,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i, j
    real(rp)             , intent(in)    :: value
    this%value(i,j) = this%value(i,j) + value
  end subroutine tensor_field_add
		
  subroutine symmetric_tensor_field_init(this,value)
    implicit none
    class(symmetric_tensor_field_t), intent(inout) :: this
    real(rp)                       , intent(in)    :: value
    this%value = value
  end subroutine symmetric_tensor_field_init

  subroutine symmetric_tensor_field_set(this,i,j,value)
    implicit none
    class(symmetric_tensor_field_t), intent(inout) :: this
    integer(ip)                    , intent(in)    :: i
    integer(ip)                    , intent(in)    :: j
    real(rp)                       , intent(in)    :: value
    assert(j>=i)
    this%value(i,j) = value
  end subroutine symmetric_tensor_field_set

  function single_contract_vector_vector(v1,v2) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v1
    type(vector_field_t), intent(in) :: v2
    real(rp)                         :: res
    integer(ip) :: k
    res=0.0_rp
    do k=1,number_space_dimensions
      res = res + v1%value(k)*v2%value(k)
    end do
   end function single_contract_vector_vector
   
   function single_contract_tensor_vector(t,v) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t
    type(vector_field_t), intent(in) :: v
    type(vector_field_t)             :: res
    integer(ip) :: i, k
    res%value=0.0_rp
    do k=1,number_space_dimensions
      do i=1,number_space_dimensions
        res%value(i) = res%value(i) + t%value(i,k) * v%value(k)
      end do
    end do
   end function single_contract_tensor_vector
   
   function single_contract_vector_tensor(v,t) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v
    type(tensor_field_t), intent(in) :: t
    type(vector_field_t)             :: res
    integer(ip) :: i, k
    res%value=0.0_rp
    do i=1,number_space_dimensions
      do k=1,number_space_dimensions
        res%value(i) = res%value(i) + v%value(k) * t%value(k,i)
      end do
    end do
   end function single_contract_vector_tensor
   
   function single_contract_tensor_tensor(t1,t2) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t1
    type(tensor_field_t), intent(in) :: t2
    type(tensor_field_t)             :: res
    integer(ip) :: i, j, k
    res%value=0.0_rp
    do i=1,number_space_dimensions
      do k=1,number_space_dimensions
        do j=1,number_space_dimensions
          res%value(i,k) = res%value(i,k) + t1%value(i,j) * t2%value(j,k)
        end do
      end do
    end do
   end function single_contract_tensor_tensor
   
   function scal_left_scalar(alpha,v) result(res)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(scalar_field_t), intent(in) :: v
    real(rp)                         :: res
    res = alpha * v%value
   end function scal_left_scalar
   
   function scal_right_scalar(v,alpha) result(res)
    implicit none
    type(scalar_field_t), intent(in) :: v
    real(rp)            , intent(in) :: alpha
    real(rp)                         :: res
    res = alpha * v%value
   end function scal_right_scalar

   function scal_left_vector(alpha,v) result(res)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(vector_field_t), intent(in) :: v
    type(vector_field_t)             :: res
    res%value = alpha * v%value
   end function scal_left_vector
   
   function scal_right_vector(v,alpha) result(res)
    implicit none
    type(vector_field_t), intent(in) :: v
    real(rp)            , intent(in) :: alpha
    type(vector_field_t)             :: res
    res%value = alpha * v%value
   end function scal_right_vector
   
   function scal_left_tensor(alpha,t) result(res)
    implicit none
    real(rp)            , intent(in) :: alpha
    type(tensor_field_t), intent(in) :: t
    type(tensor_field_t)             :: res
    res%value = alpha * t%value
   end function scal_left_tensor
   
   function scal_right_tensor(t,alpha) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t
    real(rp)            , intent(in) :: alpha
    type(tensor_field_t)             :: res
    res%value = alpha * t%value
   end function scal_right_tensor
   
   function double_contract_tensor_tensor(t1,t2) result(res)
    implicit none
    type(tensor_field_t), intent(in) :: t1
    type(tensor_field_t), intent(in) :: t2
    real(rp)                         :: res
    integer(ip) :: i, j
    res = 0.0_rp
    do j=1, number_space_dimensions
      do i=1,number_space_dimensions
        res = res + t1%value(i,j)*t2%value(i,j)
      end do
    end do
   end function double_contract_tensor_tensor
  
end module field_names
