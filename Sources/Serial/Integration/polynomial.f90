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
  use allocatable_array_names

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: NUM_POLY_DERIV = 3

  type :: polynomial_t
     !private
     integer(ip)           :: order
     real(rp), allocatable :: coefficients(:)
   contains
     procedure, non_overridable :: free => polynomial_free
     procedure, non_overridable         :: get_values      => polynomial_get_values
     procedure, nopass, non_overridable :: generate_basis => polynomial_generate_basis
     ! procedure ( polynomial_assign_interface), deferred :: assign
     ! generic(=) :: assign
  end type polynomial_t

  type :: polynomial_allocatable_array_t
     !private
     class(polynomial_t), allocatable :: polynomials(:)
   contains
     procedure, non_overridable :: free => polynomial_allocatable_array_free
  end type polynomial_allocatable_array_t

  ! type, extends(polynomial_t) :: lagrange_polynomial_t
  !    private
  !    ! For the moment, all the functionality of lagrange_polynomial_t in polynomial_t,
  !    ! by a re-interpretation of the coefficients and an if in the get_values functions
  !  contains
  ! end type lagrange_polynomial_t
  
  type :: tensor_product_polynomial_space_t
     !private
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
  end type tensor_product_polynomial_space_t
  ! Note: The vector space built from tensor_product_polynomial_space_t is going to be
  ! at the reference_fe implementation of RT and Nedelec.
  
!  type :: raviart_thomas_reference_fe_t
!     private
     ! It seems that the change-of-basis matrix is auxiliary
     ! tensor_product_polynomial_space_t probably just an auxiliary variable
!   contains
!    procedure :: fill_interpolation => raviart_thomas_fill_interpolation
     ! Inside, we will probably need a pre-interpolation + change-of-basis
!     procedure :: fill_interpolation => raviart_thomas_fill_pre_interpolation
!  end type raviart_thomas_reference_fe_t

  !type :: new_lagrangian_reference_fe_t
  !   private
  !   integer(ip) :: number_shape_functions
  !   integer(ip) :: number_dimensions
  ! contains
  !   procedure :: fill_interpolation => lagrangian_fill_interpolation
  !end type new_lagrangian_reference_fe_t

  public :: polynomial_t, polynomial_allocatable_array_t
  public :: tensor_product_polynomial_space_t

contains
 
  ! polynomial_t TBPS
  !===================================================================================

! Generate the basis of Lagrange polynomials for a given order of interpolation
! ===================================================================================================
subroutine polynomial_generate_basis ( order, basis )
  implicit none
  integer(ip)                     , intent(in)    :: order
  class(polynomial_t), allocatable, intent(inout) :: basis(:)
  type(polynomial_t) , allocatable :: aux_basis(:)
  integer(ip) :: i,j,istat
  real(rp) :: node_coordinates(order+1)
  allocate( aux_basis(order+1), stat=istat)
  do i = 0,order
     node_coordinates(i+1) = i
  end do
  node_coordinates = (2.0_rp/order)*node_coordinates-1.0_rp
  do i=1,order+1
     aux_basis(i)%order = order
     call memalloc(order+1, aux_basis(i)%coefficients, __FILE__, __LINE__)
     aux_basis(i)%coefficients(1:i-1) = node_coordinates(1:i-1)
     aux_basis(i)%coefficients(i:order) = node_coordinates(i+1:order+1)
     aux_basis(i)%coefficients(order+1)  = 1.0_rp
     do j = 1,order+1
        if ( j /= i ) then
           aux_basis(i)%coefficients(order+1) = aux_basis(i)%coefficients(order+1)*(node_coordinates(i)-node_coordinates(j))
        end if
     end do
     aux_basis(i)%coefficients(order+1) = 1/aux_basis(i)%coefficients(order+1)
  end do
  call move_alloc(from=aux_basis,to=basis)
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
  
  subroutine polynomial_free ( this )
    class(polynomial_t), intent(in)    :: this
    if ( allocated(this%coefficients) ) call memfree( this%coefficients, __FILE__, __LINE__ )
  end subroutine polynomial_free
  
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
       this%polynomial_1D_basis(i) = polynomial_1D_basis(i)
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
          !associate (poly => this%polynomial_1D_basis(i)%polynomials(j))
          !  do q = 1,n_q_points
          !     call poly%get_values(points(:,q),this%work_shape_data(i)%a(:,j,q))
          !  end do
          !end associate
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
    
    write(*,*) 'iq',q_point
    do ishape = 1, this%number_polynomials
       call index_to_ijk(ishape, this%number_dimensions, this%number_pols_dim, ijk)
          write(*,*) 'ishape',ishape
          write(*,*) 'ishape_ijk',ijk
       do idime = 1, this%number_dimensions
          write(*,*) 'iq',q_point
          write(*,*) 'dim_tensor',idime
          write(*,*) 'value pol scalar', this%work_shape_data(idime)%a(1,ijk(idime),q_point)          
          values(ishape) = values(ishape)* &
               this%work_shape_data(idime)%a(1,ijk(idime),q_point)
          write(*,*) 'value',values(ishape)
       end do
       write(*,*) 'final value',values(ishape)
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
  
  !! lagrangian_reference_fe_t TBPS
  !!===================================================================================
  !subroutine lagrangian_fill_interpolation( this, dim, quadrature, interpolation )
  !  implicit none
  !  class(new_lagrangian_reference_fe_t), intent(inout) :: this
  !  type(quadrature_t)              , intent(in)    :: quadrature
  !  type(interpolation_t)           , intent(inout) :: interpolation
  !  integer(ip)                     , intent(in)    :: dim
  !  
  !  type(tensor_product_polynomial_space_t) :: tensor_product_polynomial_space
  !  type(polynomial_allocatable_array_t)    :: polynomial_1D_basis(SPACE_DIM)
  !  type(polynomial_t)           :: polynomial_1D
  !  
  !  real(rp), allocatable :: tensor_product_values(:)
  !  real(rp), allocatable :: tensor_product_gradients(:,:)
  !  real(rp) :: aux_point(SPACE_DIM)
  !  
  !  integer(ip) :: idime, q_point
  !  
  !  
  !  this%number_dimensions = dim
  !  tensor_product_polynomial_space%number_dimensions = this%number_dimensions
  !  do idime=1,this%number_dimensions
  !     call polynomial_1D%generate_basis( 5, polynomial_1D_basis(idime)%polynomials)
  !  end do
  !  
  !  call tensor_product_polynomial_space%create(this%number_dimensions,polynomial_1D_basis)
  !  call tensor_product_polynomial_space%fill( quadrature%get_coordinates() )

  !  this%number_shape_functions = tensor_product_polynomial_space%number_polynomials
  !  
  !  call memalloc(tensor_product_polynomial_space%number_polynomials,tensor_product_values,__FILE__,__LINE__)
  !  call memalloc(SPACE_DIM,tensor_product_polynomial_space%number_polynomials,tensor_product_gradients,__FILE__,__LINE__)
  !  
  !  do q_point=1, quadrature%get_number_quadrature_points()
  !     do idime = 1, this%number_dimensions
  !        call tensor_product_polynomial_space%evaluate(q_point, &
  !             tensor_product_values, &
  !             tensor_product_gradients)
  !           interpolation%shape_functions(1,:,q_point) = tensor_product_values(:)
  !           interpolation%shape_derivatives(1,:,:,q_point) = tensor_product_gradients(:,:)
  !     end do
  !  end do
  !  
  !  ! Free auxiliary memory
  !  do idime=1,this%number_dimensions
  !     call polynomial_1D_basis(idime)%free()
  !  end do
  !  call tensor_product_polynomial_space%free()
  !  
  !end subroutine lagrangian_fill_interpolation

  ! subroutine raviart_thomas_fill_interpolation( this, quadrature, interpolation )
  !   private
  !   class(raviart_thomas_reference_fe_t), intent(inout) :: this
  !   type(quadrature_t)                  , intent(in)    :: quadrature
  !   type(interpolation_t)               , intent(inout) :: interpolation
    
  !   type(tensor_product_polynomial_space_t) :: tensor_product_polynomial_space
  !   type(polynomial_allocatable_array_t)    :: polynomial_1D_basis(SPACE_DIM)
  !   type(lagrangian_polynomial_t)           :: polynomial_1D
    
  !   real(rp), allocatable :: tensor_product_polynomials_values(:)
  !   real(rp), allocatable :: tensor_product_polynomials_gradients(:)
  !   real(rp) :: aux_point(SPACE_DIM)
    
  !   call polynomial_1D%generate_basis(order=K+1, polynomial_1D_basis(1)%polynomials)
  !   do idime=2, number_dimensions
  !      call polynomial_1D%generate_basis(order=K, polynomial_1D_basis(idime)%polynomials)
  !   end do
    
  !   call tensor_product_polynomial_space%create(polynomial_1D_basis)
  !   call tensor_product_polynomial_space%fill(polynomial_1D_basis, quadrature%get_coordinates() )
    
  !   call memalloc(tensor_product_polynomial_values, this%number_polynomials)
  !   call memalloc(tensor_product_polynomial_gradients,SPACE_DIM,this%number_polynomials)
    
  !   do q_point=1, quadrature%get_number_points()
  !      c = 0
  !      do idime = 1, this%number_dimensions
  !         aux_point(1:quadrature%number_dimensions) = &
  !              quadrature%coordinates(:,q_point)
  !         ! PERMUTE AUX_POINT AS A FUNCTION OF IDEM
  !         ! aux_point = p(aux_point)
  !         ! CALL EVALUATE
  !         call tensor_product_polynomial_space%evaluate(aux_point, &
  !              tensor_product_polynomial_values, &
  !              tensor_product_polynomial_gradients)
  !         do i=1, this%number_polynomials ! scalar shape functions
  !            c = c+1
  !            interpolation%shape_values(idime,c,q_point) = tensor_product_polynomial_values(i)
  !            interpolation%shape_values(:,idime,c,q_point) = tensor_product_polynomial_gradients(:,i)
  !         end do
  !      end do
  !   end do
    
  !   ! Free auxiliary memory
  !   do idime=1, this%number_dimensions
  !      call polynomial_1D_basis(idime)%free()
  !   end do
  !   call tensor_product_polynomial_space%free()
    
  ! end subroutine raviart_thomas_fill_interpolation
  

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
