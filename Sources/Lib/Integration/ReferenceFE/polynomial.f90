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
module polynomial_names
  use types_names
  use memor_names
  use sort_names
  use allocatable_array_names

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: NUM_POLY_DERIV = 3

  type :: polynomial_t
     private
     integer(ip)           :: order
     real(rp), allocatable :: coefficients(:)
   contains
     procedure                           :: create          => polynomial_create
     procedure, non_overridable          :: copy            => polynomial_copy
     procedure, non_overridable          :: free            => polynomial_free
     procedure                           :: get_values      => polynomial_get_values
     procedure, nopass                   :: generate_basis  => polynomial_generate_basis
     ! procedure ( polynomial_assign_interface), deferred :: assign
     ! generic(=) :: assign
  end type polynomial_t
  
  type, extends(polynomial_t) :: monomial_t
   contains
     procedure                           :: get_values      => monomial_get_values
     procedure, nopass                   :: generate_basis  => monomial_generate_basis
  end type monomial_t

  type :: polynomial_allocatable_array_t
     private
     class(polynomial_t), allocatable :: polynomials(:)
   contains
     procedure, non_overridable          :: create        => polynomial_allocatable_array_create
     procedure, non_overridable          :: copy          => polynomial_allocatable_array_copy
     procedure, non_overridable          :: free          => polynomial_allocatable_array_free
  end type polynomial_allocatable_array_t

  ! type, extends(polynomial_t) :: lagrange_polynomial_t
  !    private
  !    ! For the moment, all the functionality of lagrange_polynomial_t in polynomial_t,
  !    ! by a re-interpretation of the coefficients and an if in the get_values functions
  !  contains
  ! end type lagrange_polynomial_t
  
  type :: tensor_product_polynomial_space_t
     private
     integer(ip)                          :: number_dimensions
     integer(ip)                          :: number_polynomials
     integer(ip)                          :: number_pols_dim(SPACE_DIM)
     type(polynomial_allocatable_array_t) :: polynomial_1D_basis(SPACE_DIM)
     type(allocatable_array_rp3_t)        :: work_shape_data(SPACE_DIM)
   contains
     procedure, non_overridable :: create   => tensor_product_polynomial_space_create
     procedure, non_overridable :: fill     => tensor_product_polynomial_space_fill
     procedure, non_overridable :: evaluate => tensor_product_polynomial_space_evaluate
     procedure, non_overridable :: free     => tensor_product_polynomial_space_free
     procedure, non_overridable :: get_number_polynomials => tensor_product_polynomial_space_get_number_polynomials
     
  end type tensor_product_polynomial_space_t
  ! Note: The vector space built from tensor_product_polynomial_space_t is going to be
  ! at the reference_fe implementation of RT and Nedelec.
  
  type :: tet_polynomial_prebase_t
     private
     integer(ip)                          :: number_dimensions
     integer(ip)                          :: number_polynomials
     integer(ip)                          :: number_pols_dim(SPACE_DIM)
     type(polynomial_allocatable_array_t) :: polynomial_1D_basis(SPACE_DIM)
     type(allocatable_array_rp3_t)        :: work_shape_data(SPACE_DIM)
   contains
     procedure, non_overridable :: create   => tet_polynomial_prebase_create
     procedure, non_overridable :: fill     => tet_polynomial_prebase_fill
     procedure, non_overridable :: evaluate => tet_polynomial_prebase_evaluate
     procedure, non_overridable :: free     => tet_polynomial_prebase_free
     procedure, non_overridable :: get_number_polynomials => tet_polynomial_prebase_get_number_polynomials
     
  end type tet_polynomial_prebase_t
  
  public :: polynomial_t, monomial_t, polynomial_allocatable_array_t
  public :: tensor_product_polynomial_space_t, tet_polynomial_prebase_t

contains
 
  ! polynomial_t TBPS
  !===================================================================================

! Generate the basis of 1D Lagrange polynomials for a given order of interpolation
! ===================================================================================================
subroutine polynomial_generate_basis ( order, basis )
  implicit none
  integer(ip)                         , intent(in)    :: order
  type(polynomial_allocatable_array_t), intent(inout) :: basis
  integer(ip) :: i,j,istat
  real(rp) :: node_coordinates(order+1)
  type(polynomial_t) :: mold_polynomial
  
  call basis%create (order+1, mold_polynomial )
  do i = 0,order
     node_coordinates(i+1) = i
  end do
  node_coordinates = (2.0_rp/order)*node_coordinates-1.0_rp
  do i=1,order+1
     call basis%polynomials(i)%create(order)
     basis%polynomials(i)%coefficients(1:i-1) = node_coordinates(1:i-1)
     basis%polynomials(i)%coefficients(i:order) = node_coordinates(i+1:order+1)
     basis%polynomials(i)%coefficients(order+1)  = 1.0_rp
     do j = 1,order+1
        if ( j /= i ) then
           basis%polynomials(i)%coefficients(order+1) = basis%polynomials(i)%coefficients(order+1)*(node_coordinates(i)-node_coordinates(j))
        end if
     end do
     basis%polynomials(i)%coefficients(order+1) = 1/basis%polynomials(i)%coefficients(order+1)
  end do
end subroutine polynomial_generate_basis

! Compute the 1d shape function and n-th derivatives on ALL gauss points for ALL Lagrange polynomials
! ===================================================================================================
  subroutine polynomial_get_values (this,x,p_x)
    class(polynomial_t), intent(in)    :: this
    real(rp)           , intent(in)    :: x
    real(rp)           , intent(inout) :: p_x(3)
    integer(ip) :: i,j
    p_x = 0.0_rp
    p_x(1) = 1.0_rp
    do i = 1,this%order
       do j = size(p_x),2,-1
          p_x(j) = p_x(j)*(x-this%coefficients(i))+p_x(j-1)
       end do
       p_x(1) = p_x(1)*(x-this%coefficients(i))
    end do
    p_x = p_x*this%coefficients(this%order+1)
  end subroutine polynomial_get_values
  
  subroutine polynomial_create ( this, order )
    class(polynomial_t), intent(inout) :: this
    integer(ip)        , intent(in)    :: order 
    call this%free()
    this%order = order
    call memalloc(order+1, this%coefficients, __FILE__, __LINE__)
  end subroutine polynomial_create

  subroutine polynomial_copy (lhs, rhs)
     class(polynomial_t), intent(inout) :: lhs
     type(polynomial_t),  intent(in)    :: rhs

     call lhs%free()
     if(allocated(rhs%coefficients)) then
        call lhs%create(rhs%order)
        lhs%coefficients = rhs%coefficients
     endif
  end subroutine polynomial_copy
  
  subroutine polynomial_free ( this )
    class(polynomial_t), intent(inout)    :: this
    this%order = -1
    if ( allocated(this%coefficients) ) call memfree( this%coefficients, __FILE__, __LINE__ )
  end subroutine polynomial_free
  
  ! monomial_t TBPS
  !===================================================================================
  subroutine monomial_create ( this, order )
    class(monomial_t), intent(inout) :: this
    integer(ip)      , intent(in)    :: order 
    call this%free()
    this%order = order
  end subroutine monomial_create
  
  subroutine monomial_get_values(this,x,p_x)
    implicit none
    class(monomial_t), intent(in)    :: this
    real(rp)         , intent(in)    :: x
    real(rp)         , intent(inout) :: p_x(3)

    p_x(1) = x**real(this%order,rp)
    if ( this%order >= 1 ) then
      p_x(2) = real(this%order,rp) * x**(real(this%order,rp)-1.0_rp)
    else
      p_x(2) = 0.0_rp
    end if 
    if ( this%order >= 2 ) then
      p_x(3) = (real(this%order,rp)-1.0_rp)*real(this%order,rp) * x**(real(this%order,rp)-1.0_rp)
    else
      p_x(3) = 0.0_rp
    end if  
  end subroutine monomial_get_values
  
  subroutine monomial_generate_basis ( order, basis )
    implicit none
    integer(ip)                         , intent(in)    :: order
    type(polynomial_allocatable_array_t), intent(inout) :: basis
    integer(ip) :: i
    type(monomial_t ) :: mold_monomial
  
    call basis%create(order+1, mold_monomial)
    do i = 0,order
       call basis%polynomials(i+1)%create(i)
    end do
  end subroutine monomial_generate_basis
    
  ! tensor_product_polynomial_space_t TBPS
  !===================================================================================
  subroutine tensor_product_polynomial_space_create( this, dim, polynomial_1D_basis )
    class(tensor_product_polynomial_space_t), intent(inout) :: this
    integer(ip)                             , intent(in)    :: dim
    type(polynomial_allocatable_array_t)    , intent(in)        :: polynomial_1D_basis(:)
    integer(ip) :: i
    
    this%number_dimensions = dim
    this%number_polynomials = 1
    do i = 1, dim
       call this%polynomial_1D_basis(i)%copy(polynomial_1D_basis(i))
       this%number_pols_dim(i) = size(polynomial_1D_basis(i)%polynomials)
       this%number_polynomials = this%number_polynomials * this%number_pols_dim(i)
    end do
  end subroutine tensor_product_polynomial_space_create
    
  subroutine tensor_product_polynomial_space_fill( this, points )
    implicit none
    class(tensor_product_polynomial_space_t), intent(inout) :: this
    real(rp)                                , intent(in)    :: points(:,:)
    integer(ip)                 :: n_q_points, i, j, q
    
    do i=1,this%number_dimensions
    call this%work_shape_data(i)%free()
    end do
    n_q_points = size(points,2)
    do i=1, this%number_dimensions
       call this%work_shape_data(i)%create(NUM_POLY_DERIV, &
                                          size(this%polynomial_1D_basis(i)%polynomials), &
                                          n_q_points)
    end do
    ! Can we make it more efficient having an array of points
    do i = 1,this%number_dimensions
       do j = 1,size(this%polynomial_1D_basis(i)%polynomials)
            do q = 1,n_q_points
               call this%polynomial_1D_basis(i)%polynomials(j)%get_values(points(i,q),this%work_shape_data(i)%a(:,j,q))
            end do
       end do
    end do
  end subroutine tensor_product_polynomial_space_fill

 subroutine tensor_product_polynomial_space_evaluate( this, q_point, values, gradients )
    implicit none
    class(tensor_product_polynomial_space_t), intent(in)    :: this
    integer(ip)                             , intent(in)    :: q_point
    real(rp)                                , intent(inout) :: values(:)
    real(rp)                                , intent(inout) :: gradients(:,:)
    integer(ip) :: ijk(SPACE_DIM),idime,ishape,jdime
    values = 1.0_rp
    
    !write(*,*) 'iq',q_point
    do ishape = 1, this%number_polynomials
       call index_to_ijk(ishape, this%number_dimensions, this%number_pols_dim, ijk)
          !write(*,*) 'ishape',ishape
          !write(*,*) 'ishape_ijk',ijk
       do idime = 1, this%number_dimensions
          !write(*,*) 'iq',q_point
          !write(*,*) 'dim_tensor',idime
          !write(*,*) 'value pol scalar', this%work_shape_data(idime)%a(1,ijk(idime),q_point)          
          values(ishape) = values(ishape)* &
               this%work_shape_data(idime)%a(1,ijk(idime),q_point)
          !write(*,*) 'value',values(ishape)
       end do
       !write(*,*) 'final value',values(ishape)
    end do
    do ishape = 1, this%number_polynomials
       call index_to_ijk(ishape, this%number_dimensions, this%number_pols_dim, ijk)
       do idime = 1, this%number_dimensions
          gradients(idime,ishape) = this%work_shape_data(idime)%a(2,ijk(idime),q_point)
          do jdime = 1, this%number_dimensions
             if ( jdime /= idime ) then
                gradients(idime,ishape) = & 
                   gradients(idime,ishape) * this%work_shape_data(jdime)%a(1,ijk(jdime),q_point)
             end if
          end do
       end do
    end do    
  end subroutine tensor_product_polynomial_space_evaluate   
  
  ! tet_polynomial_prebase_t TBPS
  !===================================================================================
  subroutine tet_polynomial_prebase_create ( this, dim, polynomial_1D_prebasis )
    implicit none
    class(tet_polynomial_prebase_t)     , intent(inout) :: this
    integer(ip)                         , intent(in)    :: dim
    type(polynomial_allocatable_array_t), intent(in)    :: polynomial_1D_prebasis(:)
    integer(ip) :: i, j, work_array(3), max_polynomials, ijk(SPACE_DIM)
    
    this%number_dimensions = dim
    this%number_polynomials = 1
    this%number_pols_dim = 1
    max_polynomials = 1
    do i = 1, dim
       call this%polynomial_1D_basis(i)%copy(polynomial_1D_prebasis(i))
       this%number_pols_dim(i) = size(polynomial_1D_prebasis(i)%polynomials)
       max_polynomials = max_polynomials*this%number_pols_dim(i)
    end do
    
    work_array = this%number_pols_dim-1
    call sort(SPACE_DIM, work_array)
    i = 0
    this%number_polynomials = 0
    countp : do i=1,max_polynomials
       call index_to_ijk(i,this%number_dimensions, this%number_pols_dim, ijk)
       do j=1,dim
          if (ijk(j) > this%number_pols_dim(j))  cycle countp
       end do
       if (sum(ijk)-SPACE_DIM < maxval(this%number_pols_dim)) then
          this%number_polynomials = this%number_polynomials + 1
       end if
    end do countp
  end subroutine tet_polynomial_prebase_create
  
  subroutine tet_polynomial_prebase_fill(this, points)
    implicit none
    class(tet_polynomial_prebase_t), intent(inout) :: this
    real(rp)                       , intent(in)    :: points(:,:)
    integer(ip)                 :: n_q_points, i, j, q
    
    do i=1,this%number_dimensions
       call this%work_shape_data(i)%free()
    end do
    n_q_points = size(points,2)
    do i=1, SPACE_DIM
       call this%work_shape_data(i)%create(NUM_POLY_DERIV, &
                                           this%number_pols_dim(i), &
                                           n_q_points)
    end do
    
    do i = 1,this%number_dimensions
       do j = 1,this%number_pols_dim(i)
          do q = 1,n_q_points
             call this%polynomial_1D_basis(i)%polynomials(j)%get_values(points(i,q),this%work_shape_data(i)%a(:,j,q))
          end do
       end do
    end do
    
  end subroutine tet_polynomial_prebase_fill
  
  subroutine tet_polynomial_prebase_evaluate(this, q_point, values, gradients)
    implicit none
    class(tet_polynomial_prebase_t), intent(inout) :: this
    integer(ip)                             , intent(in)    :: q_point
    real(rp)                                , intent(inout) :: values(:)
    real(rp)                                , intent(inout) :: gradients(:,:)
    integer(ip) :: ishape, i,j, idime, jdime, max_polynomials, ijk(SPACE_DIM)
    
    
    max_polynomials = 1
    do i = 1, this%number_dimensions
       max_polynomials = max_polynomials*this%number_pols_dim(i)
    end do
    values = 1.0_rp
    gradients = 0.0_rp
    gradients(1:this%number_dimensions,:) = 1.0_rp
    ishape = 1
    countp : do i=1,max_polynomials
       call index_to_ijk(i,this%number_dimensions, this%number_pols_dim, ijk)
       do j=1,this%number_dimensions
          if (ijk(j) > this%number_pols_dim(j))  cycle countp
       end do
       if (sum(ijk)-SPACE_DIM < maxval(this%number_pols_dim)) then
          do idime=1,this%number_dimensions
             values(ishape) = values(ishape)*this%work_shape_data(idime)%a(1,ijk(idime),q_point) 
             gradients(idime,ishape) = gradients(idime,ishape)*this%work_shape_data(idime)%a(2,ijk(idime),q_point)
             do jdime=1,this%number_dimensions
                if (idime /= jdime) then
                   gradients(idime,ishape) = gradients(idime,ishape)*this%work_shape_data(jdime)%a(1,ijk(jdime),q_point)
                end if
             end do
          end do
          ishape = ishape + 1
       end if
    end do countp
  end subroutine tet_polynomial_prebase_evaluate
  
  
  subroutine tet_polynomial_prebase_free( this )
    implicit none
    class(tet_polynomial_prebase_t), intent(inout)    :: this
    integer(ip) :: i
    do i=1,SPACE_DIM
    call this%work_shape_data(i)%free()
    call this%polynomial_1D_basis(i)%free()
    end do
  end subroutine tet_polynomial_prebase_free
  
  function tet_polynomial_prebase_get_number_polynomials( this ) result(num_poly)
    implicit none
    class(tet_polynomial_prebase_t), intent(inout)    :: this
    integer(ip) :: num_poly
    num_poly = this%number_polynomials    
  end function tet_polynomial_prebase_get_number_polynomials
  
  ! polynomial_allocatable_array_t TBPS
  !===================================================================================
  subroutine polynomial_allocatable_array_create ( this, number_polynomials, mold_polynomial )
   implicit none
   class(polynomial_allocatable_array_t), intent(inout)    :: this
   integer(ip), intent(in)  :: number_polynomials
   class(polynomial_t), intent(in) :: mold_polynomial
   call this%free()
   allocate ( this%polynomials(number_polynomials), mold=mold_polynomial )
  end subroutine polynomial_allocatable_array_create
  
  subroutine polynomial_allocatable_array_copy (lhs, rhs)
   implicit none
   class(polynomial_allocatable_array_t), intent(inout) :: lhs
   type(polynomial_allocatable_array_t),  intent(in)    :: rhs
   integer(ip)                                          :: idx
   type(polynomial_t)                                   :: mold_polynomial
   call lhs%free()
   if(allocated(rhs%polynomials)) then
      if (size(rhs%polynomials) > 0) then
         call lhs%create(size(rhs%polynomials), rhs%polynomials(1))
         do idx=1, size(rhs%polynomials)
           call lhs%polynomials(idx)%copy(rhs%polynomials(idx))
         enddo   
      else
         call lhs%create(0, mold_polynomial)
      end if
   endif
  end subroutine polynomial_allocatable_array_copy
  
  subroutine polynomial_allocatable_array_free( this )
    implicit none
    class(polynomial_allocatable_array_t), intent(inout)    :: this
    integer(ip) :: istat,i
    if ( allocated(this%polynomials)) then
       do i = 1, size(this%polynomials)
             call this%polynomials(i)%free()
       end do
       deallocate( this%polynomials, stat=istat )
       check( istat == 0 )
    end if
  end subroutine polynomial_allocatable_array_free
  
  subroutine tensor_product_polynomial_space_free( this )
    implicit none
    class(tensor_product_polynomial_space_t), intent(inout)    :: this
    integer(ip) :: i
    do i=1,this%number_dimensions
    call this%work_shape_data(i)%free()
    call this%polynomial_1D_basis(i)%free()
    end do
  end subroutine tensor_product_polynomial_space_free
  
  function tensor_product_polynomial_space_get_number_polynomials ( this ) result(num_poly)
    implicit none
    class(tensor_product_polynomial_space_t), intent(in)    :: this
    integer(ip) :: num_poly
    num_poly = this%number_polynomials
  end function tensor_product_polynomial_space_get_number_polynomials 

! Support subroutines
!==================================================================================================
subroutine index_to_ijk( index, ndime, n_pols_dim, ijk )
  implicit none
  integer(ip)                         , intent(in) :: index
  integer(ip)                         , intent(in) :: ndime
  integer(ip)                         , intent(in) :: n_pols_dim(SPACE_DIM)
  integer(ip)                         , intent(inout) :: ijk(SPACE_DIM)
  integer(ip) :: i,aux

  ijk = 0
  aux = (index-1)
  do i = 1,ndime-1
     ijk(i) = mod(aux, n_pols_dim(i))
     aux = aux/n_pols_dim(i)
  end do
  ijk(ndime) = aux
  ijk = ijk+1
end subroutine index_to_ijk

end module polynomial_names
