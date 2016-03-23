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
  use serial_fe_space_names
  use environment_names
  use reference_fe_names
  use field_names
  use function_names
# include "debug.i90"
  implicit none
  private

  character(len=*), parameter :: mean_norm     = 'MEAN'
  character(len=*), parameter :: l1_norm       = 'L1-NORM'
  character(len=*), parameter :: l2_norm       = 'L2-NORM'
  character(len=*), parameter :: lp_norm       = 'Lp-NORM'
  character(len=*), parameter :: linfty_norm   = 'Linfty-NORM'  
  character(len=*), parameter :: h1_norm       = 'H1-NORM'
  character(len=*), parameter :: h1_seminorm   = 'H1-SEMINORM'
  character(len=*), parameter :: hdiv_seminorm = 'Hdiv-SEMINORM'

  type norm_t
   contains
     procedure :: norm_compute_scalar
     procedure :: norm_compute_vector
     procedure :: norm_compute_tensor
     generic   :: compute  => norm_compute_scalar,norm_compute_vector,norm_compute_tensor
     procedure :: norm_finalize_scalar
     procedure :: norm_finalize_vector
     procedure :: norm_finalize_tensor
     generic   :: finalize  => norm_finalize_scalar,norm_finalize_vector,norm_finalize_tensor
  end type norm_t

  type, extends(norm_t) :: l2_norm_t
   contains
     procedure :: norm_compute_scalar  => l2_norm_compute_scalar
     procedure :: norm_compute_vector  => l2_norm_compute_vector
     procedure :: norm_compute_tensor  => l2_norm_compute_tensor
     procedure :: norm_finalize_scalar => l2_norm_finalize_scalar
     procedure :: norm_finalize_vector => l2_norm_finalize_vector
     procedure :: norm_finalize_tensor => l2_norm_finalize_tensor
  end type l2_norm_t

  type error_norm_t
     private
     class(serial_fe_space_t), pointer :: fe_space => NULL()
     class(environment_t)    , pointer :: environment => NULL()
     class(norm_t)       , allocatable :: norm
     integer(ip)                       :: fe_space_id
     integer(ip)                       :: exponent
   contains
     procedure :: create  => error_norm_create
     procedure :: error_norm_compute_scalar
     procedure :: error_norm_compute_vector
     procedure :: error_norm_compute_tensor
     generic   :: compute => error_norm_compute_scalar, error_norm_compute_vector, error_norm_compute_tensor
     procedure :: free    => error_norm_free
  end type error_norm_t

  ! Parameters
  public :: mean_norm, l1_norm, l2_norm, lp_norm, linfty_norm, h1_norm, h1_seminorm, hdiv_seminorm

  ! Functions
  public :: error_norm_t

contains

  !===================================================================================================
  subroutine norm_compute_scalar(this,dvolume,value_in,value_out)
    implicit none
    class(norm_t), intent(in)    :: this
    real(rp)     , intent(in)    :: dvolume
    real(rp)     , intent(in)    :: value_in
    real(rp)     , intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_scalar

  !===================================================================================================
  subroutine norm_compute_vector(this,dvolume,value_in,value_out)
    implicit none
    class(norm_t)       , intent(in)    :: this
    real(rp)            , intent(in)    :: dvolume
    type(vector_field_t), intent(in)    :: value_in
    type(vector_field_t), intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_vector

  !===================================================================================================
  subroutine norm_compute_tensor(this,dvolume,value_in,value_out)
    implicit none
    class(norm_t)       , intent(in)    :: this
    real(rp)            , intent(in)    :: dvolume
    type(tensor_field_t), intent(in)    :: value_in
    type(tensor_field_t), intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_tensor

  !===================================================================================================
  function norm_finalize_scalar(this,value_in) result(value_out)
    implicit none
    class(norm_t), intent(in)    :: this
    real(rp)     , intent(in)    :: value_in
    real(rp) :: value_out
    check(.false.)
  end function norm_finalize_scalar

  !===================================================================================================
  function norm_finalize_vector(this,value_in) result(value_out)
    implicit none
    class(norm_t)       , intent(in)    :: this
    type(vector_field_t), intent(in)    :: value_in
    real(rp)    :: value_out
    check(.false.)
  end function norm_finalize_vector

  !===================================================================================================
  function norm_finalize_tensor(this,value_in) result(value_out)
    implicit none
    class(norm_t)       , intent(in)    :: this
    type(tensor_field_t), intent(in)    :: value_in
    real(rp)    :: value_out
    check(.false.)
  end function norm_finalize_tensor

  !===================================================================================================
  subroutine l2_norm_compute_scalar(this,dvolume,value_in,value_out)
    implicit none
    class(l2_norm_t), intent(in)    :: this
    real(rp)        , intent(in)    :: dvolume
    real(rp)        , intent(in)    :: value_in
    real(rp)        , intent(inout) :: value_out
    value_out = value_out + dvolume*(value_in**2)
  end subroutine l2_norm_compute_scalar

  !===================================================================================================
  subroutine l2_norm_compute_vector(this,dvolume,value_in,value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    real(rp)            , intent(in)    :: dvolume
    type(vector_field_t), intent(in)    :: value_in
    type(vector_field_t), intent(inout) :: value_out
    type(vector_field_t) :: aux_vector
    aux_vector = dvolume*(value_in*value_in)
    value_out = (aux_vector + value_out)
  end subroutine l2_norm_compute_vector

  !===================================================================================================
  subroutine l2_norm_compute_tensor(this,dvolume,value_in,value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    real(rp)            , intent(in)    :: dvolume
    type(tensor_field_t), intent(in)    :: value_in
    type(tensor_field_t), intent(inout) :: value_out
    value_out = (value_out + (dvolume*(value_in*value_in)))
  end subroutine l2_norm_compute_tensor

  !===================================================================================================
  function l2_norm_finalize_scalar(this,value_in) result(value_out)
    implicit none
    class(l2_norm_t), intent(in)    :: this
    real(rp)        , intent(in)    :: value_in
    real(rp) :: value_out
    value_out = sqrt(value_in)
  end function l2_norm_finalize_scalar  
  
  !===================================================================================================
  function l2_norm_finalize_vector(this,value_in) result(value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    type(vector_field_t), intent(in)    :: value_in
    real(rp)    :: value_out
    integer(ip) :: idime
    value_out = 0.0_rp
    do idime=1,number_space_dimensions
       value_out = value_out + value_in%get(idime)
    end do
    value_out = sqrt(value_out)
  end function l2_norm_finalize_vector
  
  !===================================================================================================
  function l2_norm_finalize_tensor(this,value_in) result(value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    type(tensor_field_t), intent(in)    :: value_in
    real(rp)    :: value_out
    integer(ip) :: idime,jdime
    value_out = 0.0_rp
    do idime=1,number_space_dimensions
       do jdime=1,number_space_dimensions
          value_out = value_out + value_in%get(idime,jdime)
       end do
    end do
    value_out = sqrt(value_out)
  end function l2_norm_finalize_tensor

  !====================================================================================================
  subroutine error_norm_create(this,fe_space,environment,fe_space_id,norm_type,exponent)
    implicit none
    class(error_norm_t)             , intent(inout) :: this
    class(serial_fe_space_t), target, intent(in)    :: fe_space
    class(environment_t)    , target, intent(in)    :: environment
    integer(ip)                     , intent(in)    :: fe_space_id
    character(len=*)                , intent(in)    :: norm_type
    integer(ip)           , optional, intent(in)    :: exponent

    call this%free()
    this%fe_space    => fe_space
    this%environment => environment
    this%fe_space_id = fe_space_id
    this%exponent    = 2
    if(present(exponent)) this%exponent = exponent
    select case(norm_type)
    case(l2_norm)
       allocate(l2_norm_t :: this%norm)
    case default
       check(.false.)
    end select

  end subroutine error_norm_create

  !====================================================================================================
  subroutine error_norm_free(this)
    implicit none
    class(error_norm_t), intent(inout) :: this

    if(associated(this%fe_space))    nullify(this%fe_space)
    if(associated(this%environment)) nullify(this%environment)
    if(allocated(this%norm))         deallocate(this%norm)
    this%fe_space_id = 0
    this%exponent    = 0
  end subroutine error_norm_free
  
  !====================================================================================================
  function error_norm_compute_scalar(this,fe_function,scalar_function,time) result(norm)
    implicit none
    class(error_norm_t)     , intent(in) :: this
    type(fe_function_t)     , intent(in) :: fe_function
    class(scalar_function_t), intent(in) :: scalar_function
    real(rp)   , optional   , intent(in) :: time
    ! Result
    real(rp) :: norm
    ! Locals
    integer(ip) :: ielem, qpoin
    real(rp)    :: dvolume, time_(1)
    real(rp)    :: unknown_at_qpoin, error_at_qpoin
    real(rp), allocatable              :: function_at_qpoin(:,:)
    type(fe_function_scalar_t)         :: fe_unknown
    type(finite_element_t)   , pointer :: fe
    type(quadrature_t)       , pointer :: quad
    type(fe_map_t)           , pointer :: fe_map
    type(point_t)            , pointer :: quadrature_coordinate(:)

    ! Check fine task
    if(this%environment%am_i_fine_task()) then
    
       ! Initialize
       norm = 0.0_rp
       call this%fe_space%initialize_integration()
       call this%fe_space%create_fe_function(this%fe_space_id,fe_unknown)
       if ( present(time) ) time_(1) = time

       ! Loop over elements
       do ielem = 1, this%fe_space%get_number_elements()

          ! Getters
          fe                    => this%fe_space%get_finite_element(ielem)
          quad                  => fe%get_quadrature()
          fe_map                => fe%get_fe_map()
          quadrature_coordinate => fe_map%get_quadrature_coordinates() 

          ! Update finite structures
          call fe%update_integration()
          call fe%update_values(fe_unknown, fe_function)

          ! Evaluate function on quadrature points
          call memalloc(quad%get_number_quadrature_points(),1,function_at_qpoin,__FILE__,__LINE__)
          if(present(time)) then
             call scalar_function%get_values_set_space_time(quadrature_coordinate,time_,function_at_qpoin)
          else
             call scalar_function%get_values_set_space(quadrature_coordinate,function_at_qpoin(:,1))
          end if

          ! Loop over quadrature points
          do qpoin = 1, quad%get_number_quadrature_points()

             ! |J]*wg
             dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)

             ! Unknown at Gauss point
             call fe_unknown%get_value(qpoin,unknown_at_qpoin)

             ! Compute norm
             error_at_qpoin = unknown_at_qpoin-function_at_qpoin(qpoin,1)
             call this%norm%compute(dvolume,error_at_qpoin,norm)

          end do

          ! Deallocate 
          call memfree(function_at_qpoin,__FILE__,__LINE__)

       end do

       ! Finalize norm
       norm = this%norm%finalize(norm)

       ! Free
       call fe_unknown%free()

    end if

    ! All reduce
    call this%environment%first_level_sum(norm)

  end function error_norm_compute_scalar
  
  !====================================================================================================
  function error_norm_compute_vector(this,fe_function,vector_function,time) result(norm)
    implicit none
    class(error_norm_t)     , intent(in) :: this
    type(fe_function_t)     , intent(in) :: fe_function
    class(vector_function_t), intent(in) :: vector_function
    real(rp)   , optional   , intent(in) :: time
    ! Result
    real(rp) :: norm
    ! Locals
    integer(ip) :: ielem, qpoin, istat
    real(rp)    :: dvolume, time_(1)
    type(vector_field_t)               :: unknown_at_qpoin
    type(vector_field_t)               :: error_at_qpoin
    type(vector_field_t)               :: vector_norm
    type(vector_field_t), allocatable  :: function_at_qpoin(:,:)
    type(fe_function_vector_t)         :: fe_unknown
    type(finite_element_t)   , pointer :: fe
    type(quadrature_t)       , pointer :: quad
    type(fe_map_t)           , pointer :: fe_map
    type(point_t)            , pointer :: quadrature_coordinate(:)

    ! Check fine task
    if(this%environment%am_i_fine_task()) then
    
       ! Initialize
       norm = 0.0_rp
       call vector_norm%init(0.0_rp)
       call this%fe_space%initialize_integration()
       call this%fe_space%create_fe_function(this%fe_space_id,fe_unknown)
       if ( present(time) ) time_(1) = time

       ! Loop over elements
       do ielem = 1, this%fe_space%get_number_elements()

          ! Getters
          fe                    => this%fe_space%get_finite_element(ielem)
          quad                  => fe%get_quadrature()
          fe_map                => fe%get_fe_map()
          quadrature_coordinate => fe_map%get_quadrature_coordinates() 

          ! Update finite structures
          call fe%update_integration()
          call fe%update_values(fe_unknown, fe_function)

          ! Evaluate function on quadrature points
          allocate(function_at_qpoin(quad%get_number_quadrature_points(),1),stat=istat)
          check(istat==0)
          if(present(time)) then
             call vector_function%get_values_set_space_time(quadrature_coordinate,time_,function_at_qpoin)
          else
             call vector_function%get_values_set_space(quadrature_coordinate,function_at_qpoin(:,1))
          end if

          ! Loop over quadrature points
          do qpoin = 1, quad%get_number_quadrature_points()

             ! |J]*wg
             dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)

             ! Unknown at Gauss point
             call fe_unknown%get_value(qpoin,unknown_at_qpoin)

             ! Compute norm
             error_at_qpoin = unknown_at_qpoin + ((-1.0_rp) * function_at_qpoin(qpoin,1))
             call this%norm%compute(dvolume,error_at_qpoin,vector_norm)

          end do

          ! Deallocate 
          deallocate(function_at_qpoin)

       end do

       ! Finalize norm
       norm = this%norm%finalize(vector_norm)

       ! Free
       call fe_unknown%free()

    end if

    ! All reduce
    call this%environment%first_level_sum(norm)

  end function error_norm_compute_vector
  
  !====================================================================================================
  function error_norm_compute_tensor(this,fe_function,vector_function,time) result(norm)
    implicit none
    class(error_norm_t)     , intent(in) :: this
    type(fe_function_t)     , intent(in) :: fe_function
    class(tensor_function_t), intent(in) :: vector_function
    real(rp)   , optional   , intent(in) :: time
    ! Result
    real(rp) :: norm
    check(.false.)
  end function error_norm_compute_tensor

end module error_norms_names
