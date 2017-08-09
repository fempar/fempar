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

  type, abstract :: polynomial_t
     private
     integer(ip)           :: order
     real(rp), allocatable :: coefficients(:)
   contains
     procedure, non_overridable                  :: create => polynomial_create
     procedure, non_overridable                  :: copy   => polynomial_copy
     procedure, non_overridable                  :: free   => polynomial_free
     procedure(polynomial_get_values_interface)            , deferred :: get_values
     procedure(polynomial_generate_basis_interface), nopass, deferred :: generate_basis
     ! procedure ( polynomial_assign_interface), deferred :: assign
     ! generic(=) :: assign
  end type polynomial_t
  
  type :: polynomial_basis_t
     private
     class(polynomial_t), allocatable :: polynomials(:)
   contains
     procedure, non_overridable          :: create        => polynomial_basis_create
     procedure, non_overridable          :: copy          => polynomial_basis_copy
     procedure, non_overridable          :: free          => polynomial_basis_free
  end type polynomial_basis_t
  
  abstract interface
    subroutine polynomial_get_values_interface( this, x, p_x)
      import :: polynomial_t, rp
      implicit none
      class(polynomial_t), intent(in)    :: this
      real(rp)           , intent(in)    :: x
      real(rp)           , intent(inout) :: p_x(3)
    end subroutine polynomial_get_values_interface
    
    subroutine polynomial_generate_basis_interface( order, basis )
      import :: polynomial_basis_t, ip
      implicit none
      integer(ip)                         , intent(in)    :: order
      type(polynomial_basis_t), intent(inout) :: basis
    end subroutine polynomial_generate_basis_interface
  end interface
  
  type, extends(polynomial_t) :: lagrange_polynomial_t
   contains
     procedure          :: get_values      => lagrange_polynomial_get_values
     procedure, nopass  :: generate_basis  => lagrange_polynomial_generate_basis
  end type lagrange_polynomial_t

  type, extends(polynomial_t) :: monomial_t
   contains
     procedure          :: get_values      => monomial_get_values
     procedure, nopass  :: generate_basis  => monomial_generate_basis
  end type monomial_t
 
  type :: tensor_product_polynomial_space_t
     private
     integer(ip)                          :: num_dims
     integer(ip)                          :: num_polynomials
     integer(ip)                          :: num_pols_dim(SPACE_DIM)
     type(polynomial_basis_t) :: polynomial_1D_basis(SPACE_DIM)
     type(allocatable_array_rp3_t)        :: work_shape_data(SPACE_DIM)
   contains
     procedure                  :: create   => tensor_product_polynomial_space_create
     procedure, non_overridable :: fill     => tensor_product_polynomial_space_fill
     procedure                  :: evaluate => tensor_product_polynomial_space_evaluate
     procedure, non_overridable :: free     => tensor_product_polynomial_space_free
     procedure, non_overridable :: get_num_polynomials => tensor_product_polynomial_space_get_num_polynomials
     
  end type tensor_product_polynomial_space_t
  ! Note: The vector space built from tensor_product_polynomial_space_t is going to be
  ! at the reference_fe implementation of RT and Nedelec.
  
  type, extends(tensor_product_polynomial_space_t) :: truncated_tensor_product_polynomial_space_t
   contains
     procedure, non_overridable :: create   => truncated_tensor_product_polynomial_space_create
     procedure, non_overridable :: evaluate => truncated_tensor_product_polynomial_space_evaluate
  end type truncated_tensor_product_polynomial_space_t
  
  public :: polynomial_t, lagrange_polynomial_t, monomial_t, polynomial_basis_t
  public :: tensor_product_polynomial_space_t, truncated_tensor_product_polynomial_space_t

contains
 
#include "sbm_polynomial.i90"
#include "sbm_tensor_product_polynomial.i90"

end module polynomial_names
