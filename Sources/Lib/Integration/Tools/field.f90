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
#ifdef memcheck
  use iso_c_binding
#endif
  implicit none
# include "debug.i90"
  private
  
  type :: vector_field_t
     !private
     real(rp) :: value(SPACE_DIM) = 0.0_rp
   contains
     procedure, non_overridable :: vector_field_init_with_scalar
     procedure, non_overridable :: vector_field_init_with_array
     generic :: init  => vector_field_init_with_scalar, vector_field_init_with_array
     procedure, non_overridable :: set       => vector_field_set
     procedure, non_overridable :: get       => vector_field_get
     procedure, non_overridable :: vector_field_add_scalar
     procedure, non_overridable :: vector_field_add_array
     generic :: add  => vector_field_add_scalar, vector_field_add_array
     procedure, non_overridable :: nrm2      => vector_field_nrm2
     procedure, non_overridable :: get_value => vector_field_get_value
  end type vector_field_t
  
  type allocatable_array_vector_field_t
    private
    type(vector_field_t), allocatable :: a(:)
  contains
    procedure, non_overridable :: create         => allocatable_array_vector_field_create
    procedure, non_overridable :: free           => allocatable_array_vector_field_free
    procedure, non_overridable :: move_alloc_out => allocatable_array_vector_field_move_alloc_out
    procedure, non_overridable :: move_alloc_in  => allocatable_array_vector_field_move_alloc_in
    procedure, non_overridable :: get_array      => allocatable_array_vector_field_get_array
  end type allocatable_array_vector_field_t
  
  interface vector_field_t
    module procedure vector_field_constructor_with_scalar, &
                     vector_field_constructor_with_array
  end interface

  type :: tensor_field_t
     !private
     real(rp)  :: value(SPACE_DIM,SPACE_DIM) = 0.0_rp
   contains
     procedure, non_overridable :: tensor_field_init_with_scalar
     procedure, non_overridable :: tensor_field_init_with_array
     generic :: init  => tensor_field_init_with_scalar, tensor_field_init_with_array
     procedure, non_overridable :: set   => tensor_field_set
     procedure, non_overridable :: get   => tensor_field_get
     procedure, non_overridable :: add   => tensor_field_add
  end type tensor_field_t
  
  type allocatable_array_tensor_field_t
    private
    type(tensor_field_t), allocatable :: a(:)
  contains
    procedure, non_overridable :: create         => allocatable_array_tensor_field_create
    procedure, non_overridable :: free           => allocatable_array_tensor_field_free
    procedure, non_overridable :: move_alloc_out => allocatable_array_tensor_field_move_alloc_out
    procedure, non_overridable :: move_alloc_in  => allocatable_array_tensor_field_move_alloc_in
    procedure, non_overridable :: get_array      => allocatable_array_tensor_field_get_array
  end type allocatable_array_tensor_field_t
  
  type :: symmetric_tensor_field_t
     !private
     real(rp)  :: value(SPACE_DIM,SPACE_DIM)
   contains   
     procedure, non_overridable :: init  => symmetric_tensor_field_init
     procedure, non_overridable :: set   => symmetric_tensor_field_set     
  end type symmetric_tensor_field_t

  type, extends(vector_field_t) :: point_t
  end type point_t

#include "field_operators.i90"

  public :: vector_field_t, tensor_field_t, symmetric_tensor_field_t, point_t 
  public :: allocatable_array_vector_field_t, allocatable_array_tensor_field_t
  ! public :: operator(*), operator(+), operator(-), assignment(=)
  ! public :: double_contract, cross_product
  ! public :: symmetric_part, trace
  public :: init_vector_field_2D_array
  public :: fill_vector_field_2D_array_with_4D_plain_array
  public :: fill_vector_field_2D_array_with_4D_plain_array_perm
  public :: compute_point_1D_array_lin_comb_with_3D_plain_array
  
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(point_t)
# define var_size 8*SPACE_DIM
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
 
contains

# include "mem_body.i90"

  function vector_field_get_value(this)
    implicit none
    class(vector_field_t), intent(in) :: this
    real(rp) :: vector_field_get_value(SPACE_DIM)
    vector_field_get_value = this%value
  end function vector_field_get_value

  subroutine vector_field_init_with_scalar(this,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine vector_field_init_with_scalar

  subroutine vector_field_init_with_array(this,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value(SPACE_DIM)
    this%value = value
  end subroutine vector_field_init_with_array

  subroutine vector_field_set(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = value
  end subroutine vector_field_set

  function vector_field_get(this,i) result(value)
    implicit none
    class(vector_field_t), intent(in) :: this
    integer(ip)          , intent(in) :: i
    real(rp)                          :: value
    value = this%value(i)
  end function vector_field_get

  subroutine vector_field_add_scalar(this,i,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i
    real(rp)             , intent(in)    :: value
    this%value(i) = this%value(i) + value
  end subroutine vector_field_add_scalar
  
  subroutine vector_field_add_array(this,value)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value(SPACE_DIM)
    this%value = this%value + value
  end subroutine vector_field_add_array

  function vector_field_nrm2(this)
    implicit none
    class(vector_field_t), intent(inout) :: this
    real(rp) :: vector_field_nrm2
    vector_field_nrm2 = this * this
    vector_field_nrm2 = sqrt(vector_field_nrm2)
  end function vector_field_nrm2
  
  ! This constructor should be used with care. In a 2D simulation,
  ! and assuming that FEMPAR was compiled with parameter constant
  ! num_space_dims == 3, then this function will also fill
  ! with a nonzero value the third component of new_vector_field
  ! (obviously if value/= 0.0_rp). This may cause trouble if the
  ! code that consumes the resulting type(vector_field_t) also 
  ! accesses the third component, as e.g., happens with all operations
  ! among vectors and tensors (single_contration, double_contraction,etc.).
  function vector_field_constructor_with_scalar(value) result(new_vector_field)
    implicit none
    real(rp), intent(in) :: value
    type(vector_field_t) :: new_vector_field
    call new_vector_field%init(value)
  end function vector_field_constructor_with_scalar
  
  function vector_field_constructor_with_array(value) result(new_vector_field)
    implicit none
    real(rp), intent(in) :: value(SPACE_DIM)
    type(vector_field_t) :: new_vector_field
    call new_vector_field%init(value)
  end function vector_field_constructor_with_array
  
  subroutine allocatable_array_vector_field_create ( this, size )
    implicit none
    class(allocatable_array_vector_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: size
     integer(ip) :: istat
    call this%free()
    allocate(this%a(size), stat=istat); check(istat==0)
  end subroutine allocatable_array_vector_field_create
 
  subroutine allocatable_array_vector_field_free ( this )
    implicit none
    class(allocatable_array_vector_field_t), intent(inout) :: this
    integer(ip) :: istat
    if (allocated(this%a)) then
      deallocate(this%a, stat=istat); check(istat==0);
    end if
  end subroutine allocatable_array_vector_field_free
  
  subroutine allocatable_array_vector_field_move_alloc_out(this, a)
    implicit none
    class(allocatable_array_vector_field_t), intent(inout) :: this
    type(vector_field_t), allocatable      , intent(inout) :: a(:)
    assert (.not. allocated (a))
    !assert (allocated(this%a))
    call move_alloc(from=this%a, to=a) 
  end subroutine allocatable_array_vector_field_move_alloc_out
  
  subroutine allocatable_array_vector_field_move_alloc_in(this, a)
    implicit none
    class(allocatable_array_vector_field_t), intent(inout) :: this
    type(vector_field_t), allocatable      , intent(inout) :: a(:)
    !assert (allocated (a))
    assert (.not. allocated(this%a))
    call move_alloc(from=a, to=this%a) 
  end subroutine allocatable_array_vector_field_move_alloc_in
  
  function allocatable_array_vector_field_get_array(this)
    implicit none
    class(allocatable_array_vector_field_t), target , intent(in) :: this
    type(vector_field_t)                   , pointer :: allocatable_array_vector_field_get_array(:)
    allocatable_array_vector_field_get_array => this%a
  end function allocatable_array_vector_field_get_array
  
  subroutine tensor_field_init_with_scalar(this,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value
    this%value = value
  end subroutine tensor_field_init_with_scalar
  
  subroutine tensor_field_init_with_array(this,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    real(rp)             , intent(in)    :: value(SPACE_DIM,SPACE_DIM)
    this%value = value
  end subroutine tensor_field_init_with_array

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
    class(tensor_field_t), intent(in) :: this
    integer(ip)          , intent(in) :: i
    integer(ip)          , intent(in) :: j
    real(rp)                          :: value
    value = this%value(i,j)
  end function tensor_field_get

  subroutine tensor_field_add(this,i,j,value)
    implicit none
    class(tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: i, j
    real(rp)             , intent(in)    :: value
    this%value(i,j) = this%value(i,j) + value
  end subroutine tensor_field_add
  
  subroutine allocatable_array_tensor_field_create ( this, size )
    implicit none
    class(allocatable_array_tensor_field_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: size
     integer(ip) :: istat
    call this%free()
    allocate(this%a(size), stat=istat); check(istat==0)
  end subroutine allocatable_array_tensor_field_create
 
  subroutine allocatable_array_tensor_field_free ( this )
    implicit none
    class(allocatable_array_tensor_field_t), intent(inout) :: this
    integer(ip) :: istat
    if (allocated(this%a)) then
      deallocate(this%a, stat=istat); check(istat==0);
    end if
  end subroutine allocatable_array_tensor_field_free
  
  subroutine allocatable_array_tensor_field_move_alloc_out(this, a)
    implicit none
    class(allocatable_array_tensor_field_t), intent(inout) :: this
    type(tensor_field_t), allocatable      , intent(inout) :: a(:)
    assert (.not. allocated (a))
    !assert (allocated(this%a))
    call move_alloc(from=this%a, to=a) 
  end subroutine allocatable_array_tensor_field_move_alloc_out
  
  subroutine allocatable_array_tensor_field_move_alloc_in(this, a)
    implicit none
    class(allocatable_array_tensor_field_t), intent(inout) :: this
    type(tensor_field_t), allocatable      , intent(inout) :: a(:)
    !assert (allocated (a))
    assert (.not. allocated(this%a))
    call move_alloc(from=a, to=this%a) 
  end subroutine allocatable_array_tensor_field_move_alloc_in
  
  function allocatable_array_tensor_field_get_array(this)
    implicit none
    class(allocatable_array_tensor_field_t), target , intent(in) :: this
    type(tensor_field_t)                   , pointer :: allocatable_array_tensor_field_get_array(:)
    allocatable_array_tensor_field_get_array => this%a
  end function allocatable_array_tensor_field_get_array
  
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
  
  ! "Friend" subroutines of reference_fe_t implementors that are here for optimization purposes.
  ! We assume that the actual arguments passed to the dummy arguments fulfill the preconditions
  ! encompassed in the asserts(...). Thus, they cannot be considered a generally applicable 
  ! subroutines, but only in the context corresponding to the caller. 
  subroutine init_vector_field_2D_array(m,n,value,vector_field_2D_array)
    implicit none
    integer(ip)         , intent(in)    :: m
    integer(ip)         , intent(in)    :: n
    real(rp)            , intent(in)    :: value
    type(vector_field_t), intent(inout) :: vector_field_2D_array(:,:)
    integer(ip) :: i,j
    do j=1,n
      do i=1,m
        vector_field_2D_array(i,j)%value = value
      end do 
    end do
  end subroutine init_vector_field_2D_array
 
  subroutine fill_vector_field_2D_array_with_4D_plain_array ( plain_array, vector_field_2D_array )
    implicit none
    real(rp)            , intent(in)    :: plain_array(:,:,:,:)
    type(vector_field_t), intent(inout) :: vector_field_2D_array(:,:)
    integer(ip) :: i, j
    assert (size(vector_field_2D_array,1) == size(plain_array,3))
    assert (size(vector_field_2D_array,2) >= size(plain_array,4))
    assert (size(plain_array,2) == SPACE_DIM)
    assert (size(plain_array,1) >= 1)
    do j=1, size(plain_array,4)
      do i=1, size(plain_array,3)
        vector_field_2D_array(i,j)%value(:) = plain_array(1,:,i,j)
      end do
    end do
  end subroutine fill_vector_field_2D_array_with_4D_plain_array
  
  subroutine fill_vector_field_2D_array_with_4D_plain_array_perm ( plain_array, perm, vector_field_2D_array)
    implicit none
    real(rp)            , intent(in)    :: plain_array(:,:,:,:)
    integer(ip)         , intent(in)    :: perm(:)
    type(vector_field_t), intent(inout) :: vector_field_2D_array(:,:)
    integer(ip) :: i, j
    assert (size(vector_field_2D_array,1) == size(plain_array,3))
    assert (size(vector_field_2D_array,2) >= size(plain_array,4))
    assert (size(perm) >= size(plain_array,4)) 
    assert (size(plain_array,2) == SPACE_DIM)
    assert (size(plain_array,1) >= 1)
    do j=1, size(plain_array,4)
      do i=1, size(plain_array,3)
        vector_field_2D_array(i,j)%value(:) = plain_array(1,:,i,perm(j))
      end do
    end do
  end subroutine fill_vector_field_2D_array_with_4D_plain_array_perm 
  
  subroutine compute_point_1D_array_lin_comb_with_3D_plain_array ( plain_array, point_1D_array_in, point_1D_array_out )
    implicit none
    real(rp)     , intent(in)    :: plain_array(:,:,:)
    type(point_t), intent(in)    :: point_1D_array_in(:)
    type(point_t), intent(inout) :: point_1D_array_out(:)
    integer(ip) :: i, j, k
    real(rp) :: alpha
    assert ( size(plain_array,1) >= 1 )
    assert ( size(point_1D_array_out) == size(plain_array,3) )
    assert ( size(point_1D_array_in)  == size(plain_array,2) )
    do i=1, size(point_1D_array_out) 
      point_1D_array_out(i)%value(:) = 0.0_rp 
      do j=1, size(point_1D_array_in)
        alpha = plain_array(1,j,i)
        do k=1, SPACE_DIM
          point_1D_array_out(i)%value(k) = point_1D_array_out(i)%value(k) + alpha * point_1D_array_in(j)%value(k)
        end do
      end do 
    end do
  end subroutine compute_point_1D_array_lin_comb_with_3D_plain_array

#include "field_operators_overloads.i90"
  
end module field_names





