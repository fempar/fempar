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
module error_norms_names
  use types_names
  use memor_names
  use fe_space_names
  use environment_names
  use reference_fe_names
  use field_names
  use function_names
# include "debug.i90"
  implicit none
  private

  character(len=*), parameter :: mean_norm         = 'mean_norm'        ! The function is integrated on each cell
  character(len=*), parameter :: l1_norm           = 'l1_norm'          ! The absolute value of the function is integrated on each cell
  character(len=*), parameter :: l2_norm           = 'l2_norm'          ! The square of the function is integrated on each cell and the square root of the sum over all cells is taken 
  character(len=*), parameter :: lp_norm           = 'lp_norm'          ! The absolute value to the pth power is integrated and the pth root of the sum over all cells is taken
  character(len=*), parameter :: linfty_norm       = 'linfty_norm'      ! The maximum absolute value of the function is computed on each cell, the maximum over all cells is taken
  character(len=*), parameter :: h1_seminorm       = 'h1_seminorm'      ! L2_norm of the gradient. 
  character(len=*), parameter :: hdiv_seminorm     = 'hdiv_seminorm'    ! L2_norm of the divergence of a vector field 
  character(len=*), parameter :: h1_norm           = 'h1_norm'          ! The square of this norm is the square of the L2_norm plus the square of the H1_seminorm
  character(len=*), parameter :: w1p_seminorm      = 'w1p_seminorm'     ! Lp_norm of the gradient
  character(len=*), parameter :: w1p_norm          = 'w1p_norm'         ! same as H1_norm for Lp
  character(len=*), parameter :: w1infty_seminorm  = 'w1infty_seminorm' ! Linfty_norm of the gradient
  character(len=*), parameter :: w1infty_norm      = 'w1infty_norm'     ! same as H1_norm for Linfty
  
  ! Data type which implements the norm of the difference of an scalar-valued (user-defined) 
  ! exact function (i.e., class(scalar_function_t)) and one of the scalar fields of a multi-field 
  ! type(fe_function_t). The norm to be evaluated should be chosen from the suite of norms declared above. 
  ! If a new norm is to be implemented,  then the module is not designed  with extensibility in mind, i.e., 
  ! the implementation of this data type's TBPs should be modified in order to accommodate for a new norm 
  ! (if it fits in the implementation pattern followed within the module).
  type error_norms_scalar_t
     private
     ! FE space over which the norm of the difference is to be computed
     class(serial_fe_space_t), pointer    :: fe_space
     
     ! Field identifier of the field for which this object is to compute
     ! the norm of the difference
     integer(ip)                           :: field_id
     
     ! Work arrays to store type(fe_function_t) values and gradients 
     ! restricted to field_id + current cell
     type(cell_fe_function_scalar_t)       :: cell_fe_function
     
     ! Work arrays to store exact function (and difference of fe - exact functions) 
     ! values and gradients. Size = (max_number_of_quadrature_points, 1). A 2-rank array 
     ! is required provided the current interface of tensor-valued function_t data types
     real(rp)                , allocatable :: work_array_values(:,:)
     type(vector_field_t)    , allocatable :: work_array_gradients(:,:)
   contains
     procedure, non_overridable          :: create                    => error_norms_scalar_create
     procedure, non_overridable          :: free                      => error_norms_scalar_free
     procedure, non_overridable          :: compute                   => error_norms_scalar_compute
     procedure, non_overridable, private :: compute_cell_contribution => error_norms_scalar_compute_cell_contribution
  end type error_norms_scalar_t
  
  ! Data type which implements the norm of the difference of an vector-valued (user-defined) 
  ! exact function (i.e., class(vector_function_t)) and one of the vector fields of a multi-field 
  ! type(fe_function_t). The norm to be evaluated should be chosen from the suite of norms declared above. 
  ! If a new norm is to be implemented,  then the module is not designed  with extensibility in mind, i.e., 
  ! the implementation of this data type's TBPs should be modified in order to accommodate for a new norm 
  ! (if it fits in the implementation pattern followed within the module).
  type error_norms_vector_t
     private
     ! FE space over which the norm of the difference is to be computed
     class(serial_fe_space_t), pointer    :: fe_space
     
     ! Field identifier of the field for which this object is to compute
     ! the norm of the difference
     integer(ip)                           :: field_id
     
     ! Work arrays to store type(fe_function_t) values and gradients 
     ! restricted to field_id + current cell
     type(cell_fe_function_vector_t)       :: cell_fe_function
     
     ! Work arrays to store exact function (and difference of fe - exact functions) 
     ! values and gradients. Size = (max_number_of_quadrature_points, 1). A 2-rank array 
     ! is required provided the current interface of tensor-valued function_t data types
     type(vector_field_t)    , allocatable :: work_array_values(:,:)
     type(tensor_field_t)    , allocatable :: work_array_gradients(:,:)
   contains
     procedure, non_overridable          :: create                    => error_norms_vector_create
     procedure, non_overridable          :: free                      => error_norms_vector_free
     procedure, non_overridable          :: compute                   => error_norms_vector_compute
     procedure, non_overridable, private :: compute_cell_contribution => error_norms_vector_compute_cell_contribution
  end type error_norms_vector_t
  

  ! Parameters
  public :: mean_norm, l1_norm, l2_norm, lp_norm, linfty_norm, h1_seminorm
  public :: hdiv_seminorm, h1_norm, w1p_seminorm, w1p_norm, w1infty_seminorm, w1infty_norm

  ! Data types
  public :: error_norms_scalar_t, error_norms_vector_t !, error_norms_tensor_t

contains

  subroutine update_norm (norm_type, cell_contribution, norm)
    implicit none
    character(*), intent(in)    :: norm_type
    real(rp)    , intent(in)    :: cell_contribution
    real(rp)    , intent(inout) :: norm
   
    select case ( trim(norm_type) )
    case (linfty_norm,W1infty_norm,w1infty_seminorm)
       norm = max(cell_contribution,norm)
    case (mean_norm, l1_norm, l2_norm, lp_norm, h1_norm, h1_seminorm, hdiv_seminorm, w1p_norm, w1p_seminorm) 
       norm = norm + cell_contribution
    end select 
  end subroutine update_norm
  
  subroutine finalize_norm (environment, norm_type, exponent, values_norm, gradients_norm, norm)
    implicit none
    class(environment_t), intent(in)  :: environment
    character(*)        , intent(in)  :: norm_type
    real(rp)            , intent(in)  :: exponent
    real(rp)            , intent(in)  :: values_norm 
    real(rp)            , intent(in)  :: gradients_norm 
    real(rp)            , intent(out) :: norm 

    real(rp) :: aux

    norm = 0.0_rp
    select case ( trim(norm_type) )
    case(mean_norm, l1_norm, l2_norm, lp_norm, h1_norm, w1p_norm) 
       norm = values_norm
       call environment%l1_sum(norm)
    case(linfty_norm, w1infty_norm)
       norm = values_norm
       call environment%l1_max(norm)
    end select

    select case ( trim(norm_type) )
    case(h1_norm, h1_seminorm, hdiv_seminorm,w1p_norm, w1p_seminorm) 
       aux = gradients_norm
       call environment%l1_sum(aux)
       norm = norm + aux
    case(linfty_norm, w1infty_norm, w1infty_seminorm)
       aux = gradients_norm
       call environment%l1_max(aux)
       norm = norm + aux
    end select

    select case( trim(norm_type) )
    case (l2_norm,h1_seminorm,hdiv_seminorm,h1_norm)
      norm = sqrt(norm)
    case (lp_norm,w1p_norm,w1p_seminorm)
      norm = norm**exponent
    end select
  end subroutine finalize_norm


  function error_norm_is_supported(norm_type) result(is_supported)
    implicit none
    character(*), intent(in) :: norm_type
    logical :: is_supported
    is_supported =  ((trim(norm_type) == mean_norm) .or. &
                     (trim(norm_type) == l1_norm) .or. &
                     (trim(norm_type) == l2_norm) .or. &
                     (trim(norm_type) == lp_norm) .or. & 
                     (trim(norm_type) == linfty_norm) .or. &
                     (trim(norm_type) == h1_norm) .or. &
                     (trim(norm_type) == w1p_norm) .or. & 
                     (trim(norm_type) == w1infty_norm) .or. &
                     (trim(norm_type) == h1_seminorm) .or. &
                     (trim(norm_type) == hdiv_seminorm) .or. &
                     (trim(norm_type) == w1p_seminorm) .or. &
                     (trim(norm_type) == w1infty_seminorm) )
  end function error_norm_is_supported
  
  ! Private helper function
  function error_norm_requires_values(norm_type) result(requires)
    implicit none
    character(*), intent(in) :: norm_type
    logical :: requires 
    assert ( error_norm_is_supported(norm_type) )
    requires = ( (trim(norm_type) == mean_norm) .or. &
                 (trim(norm_type) == l1_norm) .or. &
                 (trim(norm_type) == l2_norm) .or. &
                 (trim(norm_type) == lp_norm) .or. &
                 (trim(norm_type) == linfty_norm) .or. &
                 (trim(norm_type) == h1_norm) .or. &
                 (trim(norm_type) == w1p_norm) .or. &
                 (trim(norm_type) == w1infty_norm) )
  end function error_norm_requires_values
  
  ! Private helper function
  function error_norm_requires_gradients(norm_type) result(requires)
    implicit none
    character(*), intent(in) :: norm_type
    logical                  :: requires 
    assert ( error_norm_is_supported(norm_type) )
    requires = ( (trim(norm_type) == h1_seminorm) .or. &
                 (trim(norm_type) == hdiv_seminorm) .or. &
                 (trim(norm_type) == w1p_seminorm) .or. &
                 (trim(norm_type) == w1infty_seminorm) .or. &
                 (trim(norm_type) == h1_norm) .or. &
                 (trim(norm_type) == w1p_norm) .or. &
                 (trim(norm_type) == w1infty_norm) )
  end function error_norm_requires_gradients
  
  ! Private helper function
  function error_norm_determine_exponent(norm_type, exponent) result(exponent_)
    implicit none
    character(*)         , intent(in) :: norm_type
    integer(ip), optional, intent(in) :: exponent
    real(rp) :: exponent_

    assert ( error_norm_is_supported(norm_type) )

    if ( present(exponent) ) then
      exponent_ = exponent
    else
      exponent_ = 2.0_rp
    end if

    ! Determine "exponent_" for those norms where it 
    ! is inherently defined by the norm itself 
    select case ( trim(norm_type) )
    case (l2_norm,h1_seminorm,h1_norm,hdiv_seminorm)
      exponent_ = 2.0_rp
    case (l1_norm)   
      exponent_ = 1.0_rp
    end select 
    
  end function error_norm_determine_exponent
  
#include "sbm_error_norms_scalar.i90"
#include "sbm_error_norms_vector.i90"
  
end module error_norms_names
