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
     procedure, private :: norm_compute_scalar
     procedure, private :: norm_compute_vector
     procedure, private :: norm_compute_tensor
     generic            :: compute  => norm_compute_scalar,norm_compute_vector,norm_compute_tensor
     procedure          :: finalize => norm_finalize
  end type norm_t

  type, extends(norm_t) :: l1_norm_t
   contains
     procedure :: norm_compute_scalar  => l1_norm_compute_scalar
     procedure :: norm_compute_vector  => l1_norm_compute_vector
     procedure :: norm_compute_tensor  => l1_norm_compute_tensor
     procedure :: finalize             => l1_norm_finalize
  end type l1_norm_t

  type, extends(norm_t) :: l2_norm_t
   contains
     procedure :: norm_compute_scalar  => l2_norm_compute_scalar
     procedure :: norm_compute_vector  => l2_norm_compute_vector
     procedure :: norm_compute_tensor  => l2_norm_compute_tensor
     procedure :: finalize             => l2_norm_finalize
  end type l2_norm_t

  type error_norm_scalar_t
     private
     class(serial_fe_space_t), pointer     :: fe_space => NULL()
     class(environment_t)    , pointer     :: environment => NULL()
     class(norm_t)           , allocatable :: norm
     integer(ip)                           :: fe_space_id
     real(rp)                , allocatable :: qpoints_work_array(:,:)
   contains
     procedure          :: create   => error_norm_scalar_create
     procedure          :: compute  => error_norm_scalar_compute
     procedure          :: free     => error_norm_scalar_free
  end type error_norm_scalar_t

  type error_norm_vector_t
     private
     class(serial_fe_space_t), pointer     :: fe_space => NULL()
     class(environment_t)    , pointer     :: environment => NULL()
     class(norm_t)           , allocatable :: norm
     integer(ip)                           :: fe_space_id
     type(vector_field_t)    , allocatable :: qpoints_work_array(:,:)
   contains
     procedure          :: create   => error_norm_vector_create
     procedure          :: compute  => error_norm_vector_compute
     procedure          :: free     => error_norm_vector_free
  end type error_norm_vector_t


  ! Parameters
  public :: mean_norm, l1_norm, l2_norm, lp_norm, linfty_norm, h1_norm, h1_seminorm, hdiv_seminorm

  ! Functions
  public :: error_norm_scalar_t, error_norm_vector_t

contains

  !===================================================================================================
  subroutine norm_compute_scalar(this,quad,fe_map,value_in,value_out)
    implicit none
    class(norm_t)     , intent(in)    :: this
    type(quadrature_t), intent(in)    :: quad
    type(fe_map_t)    , intent(in)    :: fe_map
    real(rp)          , intent(in)    :: value_in(:)
    real(rp)          , intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_scalar

  !===================================================================================================
  subroutine norm_compute_vector(this,quad,fe_map,value_in,value_out)
    implicit none
    class(norm_t)     , intent(in)      :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(vector_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_vector

  !===================================================================================================
  subroutine norm_compute_tensor(this,quad,fe_map,value_in,value_out)
    implicit none
    class(norm_t)       , intent(in)    :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(tensor_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    check(.false.)
  end subroutine norm_compute_tensor

  !===================================================================================================
  function norm_finalize(this,value_in) result(value_out)
    implicit none
    class(norm_t), intent(in)    :: this
    real(rp)     , intent(in)    :: value_in
    real(rp) :: value_out
    check(.false.)
  end function norm_finalize

  !===================================================================================================
  subroutine l1_norm_compute_scalar(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l1_norm_t)  , intent(in)    :: this
    type(quadrature_t), intent(in)    :: quad
    type(fe_map_t)    , intent(in)    :: fe_map
    real(rp)          , intent(in)    :: value_in(:)
    real(rp)          , intent(inout) :: value_out
    real(rp)    :: dvolume
    integer(ip) :: qpoin
    
    ! Loop over quadrature points
    do qpoin = 1, quad%get_number_quadrature_points()
       ! |J]*wg
       dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)
       ! Compute norm ( abs(f(u)) )
       value_out = value_out + dvolume*abs(value_in(qpoin))
    end do

  end subroutine l1_norm_compute_scalar

  !===================================================================================================
  subroutine l1_norm_compute_vector(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l1_norm_t)    , intent(in)    :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(vector_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    real(rp)                            :: dvolume
    integer(ip)                         :: qpoin,idime

    ! Loop over quadrature points
    do qpoin = 1, quad%get_number_quadrature_points()
       ! |J]*wg
       dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)
       do idime = 1,number_space_dimensions

          ! Compute norm ( sum_{i=1}^{ndime} abs(f(u)_i) )
          value_out = value_out + dvolume*abs(value_in(qpoin)%get(idime))

       end do
    end do
  end subroutine l1_norm_compute_vector

  !===================================================================================================
  subroutine l1_norm_compute_tensor(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l1_norm_t)    , intent(in)    :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(tensor_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    check(.false.)
  end subroutine l1_norm_compute_tensor

  !===================================================================================================
  function l1_norm_finalize(this,value_in) result(value_out)
    implicit none
    class(l1_norm_t), intent(in)    :: this
    real(rp)        , intent(in)    :: value_in
    real(rp) :: value_out
    value_out = value_in
  end function l1_norm_finalize  

  !===================================================================================================
  subroutine l2_norm_compute_scalar(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l2_norm_t)  , intent(in)    :: this
    type(quadrature_t), intent(in)    :: quad
    type(fe_map_t)    , intent(in)    :: fe_map
    real(rp)          , intent(in)    :: value_in(:)
    real(rp)          , intent(inout) :: value_out
    real(rp)    :: dvolume
    integer(ip) :: qpoin
    
    ! Loop over quadrature points
    do qpoin = 1, quad%get_number_quadrature_points()
       ! |J]*wg
       dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)
       ! Compute norm ( f(u)**2 )
       value_out = value_out + dvolume*(value_in(qpoin)**2)
    end do

  end subroutine l2_norm_compute_scalar

  !===================================================================================================
  subroutine l2_norm_compute_vector(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(vector_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    real(rp)                            :: dvolume
    integer(ip)                         :: qpoin

    ! Loop over quadrature points
    do qpoin = 1, quad%get_number_quadrature_points()
       ! |J]*wg
       dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)

       ! Compute norm ( f(u)**2 )
       value_out = value_out + dvolume*(value_in(qpoin)*value_in(qpoin))
    end do
  end subroutine l2_norm_compute_vector

  !===================================================================================================
  subroutine l2_norm_compute_tensor(this,quad,fe_map,value_in,value_out)
    implicit none
    class(l2_norm_t)    , intent(in)    :: this
    type(quadrature_t)  , intent(in)    :: quad
    type(fe_map_t)      , intent(in)    :: fe_map
    type(tensor_field_t), intent(in)    :: value_in(:)
    real(rp)            , intent(inout) :: value_out
    check(.false.)
  end subroutine l2_norm_compute_tensor

  !===================================================================================================
  function l2_norm_finalize(this,value_in) result(value_out)
    implicit none
    class(l2_norm_t), intent(in)    :: this
    real(rp)        , intent(in)    :: value_in
    real(rp) :: value_out
    value_out = sqrt(value_in)
  end function l2_norm_finalize  
  
  !====================================================================================================
  subroutine error_norm_scalar_create(this,fe_space,environment,fe_space_id,norm_type,exponent)
    implicit none
    class(error_norm_scalar_t)      , intent(inout) :: this
    class(serial_fe_space_t), target, intent(in)    :: fe_space
    class(environment_t)    , target, intent(in)    :: environment
    integer(ip)                     , intent(in)    :: fe_space_id
    character(len=*)                , intent(in)    :: norm_type
    integer(ip)           , optional, intent(in)    :: exponent

    call this%free()
    this%fe_space    => fe_space
    this%environment => environment
    this%fe_space_id = fe_space_id
    select case(norm_type)
    case(l2_norm)
       allocate(l2_norm_t :: this%norm)
    case default
       check(.false.)
    end select
    
    call memalloc(fe_space%get_max_number_quadrature_points(),1,this%qpoints_work_array,__FILE__,__LINE__)
  end subroutine error_norm_scalar_create

  !====================================================================================================
  subroutine error_norm_scalar_free(this)
    implicit none
    class(error_norm_scalar_t), intent(inout) :: this

    if(associated(this%fe_space))    nullify(this%fe_space)
    if(associated(this%environment)) nullify(this%environment)
    if(allocated(this%norm))         deallocate(this%norm)
    this%fe_space_id = 0
    if (allocated(this%qpoints_work_array)) call memfree(this%qpoints_work_array,__FILE__,__LINE__)
  end subroutine error_norm_scalar_free
  
  !====================================================================================================
  function error_norm_scalar_compute(this,fe_function,scalar_function,time) result(norm)
    implicit none
    class(error_norm_scalar_t), intent(inout) :: this
    type(fe_function_t)       , intent(in)    :: fe_function
    class(scalar_function_t)  , intent(in)    :: scalar_function
    real(rp)   , optional     , intent(in)    :: time
    ! Result
    real(rp) :: norm

    ! Locals
    integer(ip)                            :: ielem
    real(rp)                               :: time_(1)
    type(fe_function_scalar_t)             :: fe_unknown
    type(finite_element_t)   , pointer     :: fe
    type(quadrature_t)       , pointer     :: quad
    type(fe_map_t)           , pointer     :: fe_map
    type(point_t)            , pointer     :: quadrature_coordinate(:)
    real(rp)                 , pointer     :: unknown_at_qpoints(:)
    integer(ip)                            :: number_qpoints

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
          number_qpoints        = quad%get_number_quadrature_points()

          ! Update finite structures
          call fe%update_integration()
          call fe%update_values(fe_unknown, fe_function)

          ! Evaluate function on quadrature points
          if(present(time)) then
             call scalar_function%get_values_set_space_time(quadrature_coordinate,time_,this%qpoints_work_array)
          else
             call scalar_function%get_values_set_space(quadrature_coordinate,this%qpoints_work_array(:,1))
          end if

          unknown_at_qpoints => fe_unknown%get_quadrature_points_values()
          this%qpoints_work_array(1:number_qpoints,1) = this%qpoints_work_array(1:number_qpoints,1) - unknown_at_qpoints(1:number_qpoints)

          ! Compute element contribution to the global norm
          call this%norm%compute(quad,fe_map,this%qpoints_work_array(:,1),norm)
       end do
       
       ! Free
       call fe_unknown%free()

       ! All reduce
       call this%environment%first_level_sum(norm)

       ! Finalize norm
       norm = this%norm%finalize(norm)
    end if
  end function error_norm_scalar_compute


  !====================================================================================================
  subroutine error_norm_vector_create(this,fe_space,environment,fe_space_id,norm_type,exponent)
    implicit none
    class(error_norm_vector_t)      , intent(inout) :: this
    class(serial_fe_space_t), target, intent(in)    :: fe_space
    class(environment_t)    , target, intent(in)    :: environment
    integer(ip)                     , intent(in)    :: fe_space_id
    character(len=*)                , intent(in)    :: norm_type
    integer(ip)           , optional, intent(in)    :: exponent

    call this%free()
    this%fe_space    => fe_space
    this%environment => environment
    this%fe_space_id = fe_space_id
    select case(norm_type)
    case(l2_norm)
       allocate(l2_norm_t :: this%norm)
    case default
       check(.false.)
    end select
    
    allocate(this%qpoints_work_array(fe_space%get_max_number_quadrature_points(),1))
  end subroutine error_norm_vector_create

  !====================================================================================================
  subroutine error_norm_vector_free(this)
    implicit none
    class(error_norm_vector_t), intent(inout) :: this

    if(associated(this%fe_space))    nullify(this%fe_space)
    if(associated(this%environment)) nullify(this%environment)
    if(allocated(this%norm))         deallocate(this%norm)
    this%fe_space_id = 0
    if (allocated(this%qpoints_work_array)) deallocate(this%qpoints_work_array)
  end subroutine error_norm_vector_free
  
  !====================================================================================================
  function error_norm_vector_compute(this,fe_function,vector_function,time) result(norm)
    implicit none
    class(error_norm_vector_t), intent(inout) :: this
    type(fe_function_t)       , intent(in)    :: fe_function
    class(vector_function_t)  , intent(in)    :: vector_function
    real(rp)   , optional     , intent(in)    :: time
    ! Result
    real(rp) :: norm

    ! Locals
    integer(ip)                            :: ielem
    real(rp)                               :: time_(1)
    type(fe_function_vector_t)             :: fe_unknown
    type(finite_element_t)   , pointer     :: fe
    type(quadrature_t)       , pointer     :: quad
    type(fe_map_t)           , pointer     :: fe_map
    type(point_t)            , pointer     :: quadrature_coordinate(:)
    type(vector_field_t)     , pointer     :: unknown_at_qpoints(:)
    integer(ip)                            :: qpoin, number_qpoints

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
          number_qpoints        = quad%get_number_quadrature_points()

          ! Update finite structures
          call fe%update_integration()
          call fe%update_values(fe_unknown, fe_function)

          ! Evaluate function on quadrature points
          if(present(time)) then
             call vector_function%get_values_set_space_time(quadrature_coordinate,time_,this%qpoints_work_array)
          else
             call vector_function%get_values_set_space(quadrature_coordinate,this%qpoints_work_array(:,1))
          end if

          unknown_at_qpoints => fe_unknown%get_quadrature_points_values()
          do qpoin=1, number_qpoints
             this%qpoints_work_array(qpoin,1) = this%qpoints_work_array(qpoin,1) - unknown_at_qpoints(qpoin)
          end do
          
          ! Compute element contribution to the global norm
          call this%norm%compute(quad,fe_map,this%qpoints_work_array(:,1),norm)
       end do
       
       ! Free
       call fe_unknown%free()

       ! All reduce
       call this%environment%first_level_sum(norm)

       ! Finalize norm
       norm = this%norm%finalize(norm)
    end if
  end function error_norm_vector_compute

end module error_norms_names
