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

  type, extends(norm_t) :: h1_seminorm_t
   contains
     procedure :: norm_compute_scalar  => h1_seminorm_compute_scalar
     procedure :: norm_compute_vector  => h1_seminorm_compute_vector
     procedure :: norm_compute_tensor  => h1_seminorm_compute_tensor
     procedure :: finalize             => h1_seminorm_finalize
  end type h1_seminorm_t

  type error_norm_scalar_t
     private
     class(serial_fe_space_t), pointer     :: fe_space    => NULL()
     class(environment_t)    , pointer     :: environment => NULL()
     class(norm_t)           , allocatable :: norm
     integer(ip)                           :: field_id
     real(rp)                , allocatable :: quadrature_points_values(:,:)
     type(vector_field_t)    , allocatable :: quadrature_points_gradients(:,:)
   contains
     procedure          :: create            => error_norm_scalar_create
     procedure          :: compute_values    => error_norm_scalar_compute_values
     procedure          :: compute_gradients => error_norm_scalar_compute_gradients
     generic            :: compute           => compute_values, compute_gradients
     procedure          :: free              => error_norm_scalar_free
  end type error_norm_scalar_t

  type error_norm_vector_t
     private
     class(serial_fe_space_t), pointer     :: fe_space    => NULL()
     class(environment_t)    , pointer     :: environment => NULL()
     class(norm_t)           , allocatable :: norm
     integer(ip)                           :: field_id
     type(vector_field_t)    , allocatable :: quadrature_points_values(:,:)
     type(tensor_field_t)    , allocatable :: quadrature_points_gradients(:,:)
   contains
     procedure          :: create            => error_norm_vector_create
     procedure          :: compute_values    => error_norm_vector_compute_values
     procedure          :: compute_gradients => error_norm_vector_compute_gradients
     generic            :: compute           => compute_values, compute_gradients
     procedure          :: free              => error_norm_vector_free
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
  real(rp)    :: dvolume
  integer(ip) :: qpoin,idime

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
  real(rp)    :: dvolume
  integer(ip) :: qpoin

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

!===================================================================================================
subroutine h1_seminorm_compute_scalar(this,quad,fe_map,value_in,value_out)
  implicit none
  class(h1_seminorm_t), intent(in)    :: this
  type(quadrature_t)  , intent(in)    :: quad
  type(fe_map_t)      , intent(in)    :: fe_map
  real(rp)            , intent(in)    :: value_in(:)
  real(rp)            , intent(inout) :: value_out
  check(.false.)
end subroutine h1_seminorm_compute_scalar

!===================================================================================================
subroutine h1_seminorm_compute_vector(this,quad,fe_map,value_in,value_out)
  implicit none
  class(h1_seminorm_t), intent(in)    :: this
  type(quadrature_t)  , intent(in)    :: quad
  type(fe_map_t)      , intent(in)    :: fe_map
  type(vector_field_t), intent(in)    :: value_in(:)
  real(rp)            , intent(inout) :: value_out
  real(rp)    :: dvolume
  integer(ip) :: qpoin

  ! Loop over quadrature points
  do qpoin = 1, quad%get_number_quadrature_points()
     ! |J]*wg
     dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)

     ! Compute norm ( f(u)**2 )
     value_out = value_out + dvolume*(value_in(qpoin)*value_in(qpoin))
  end do
end subroutine h1_seminorm_compute_vector

!===================================================================================================
subroutine h1_seminorm_compute_tensor(this,quad,fe_map,value_in,value_out)
  implicit none
  class(h1_seminorm_t), intent(in)    :: this
  type(quadrature_t)  , intent(in)    :: quad
  type(fe_map_t)      , intent(in)    :: fe_map
  type(tensor_field_t), intent(in)    :: value_in(:)
  real(rp)            , intent(inout) :: value_out
  real(rp)    :: dvolume
  integer(ip) :: qpoin

  ! Loop over quadrature points
  do qpoin = 1, quad%get_number_quadrature_points()
     ! |J]*wg
     dvolume = fe_map%get_det_jacobian(qpoin) * quad%get_weight(qpoin)

     ! Compute norm ( f(u)**2 )
     value_out = value_out + dvolume*double_contract(value_in(qpoin),value_in(qpoin))
  end do
end subroutine h1_seminorm_compute_tensor

!===================================================================================================
function h1_seminorm_finalize(this,value_in) result(value_out)
  implicit none
  class(h1_seminorm_t), intent(in)    :: this
  real(rp)            , intent(in)    :: value_in
  real(rp) :: value_out
  value_out = sqrt(value_in)
end function h1_seminorm_finalize

!====================================================================================================
subroutine error_norm_scalar_create(this,fe_space,environment,field_id,norm_type,exponent)
  implicit none
  class(error_norm_scalar_t)      , intent(inout) :: this
  class(serial_fe_space_t), target, intent(in)    :: fe_space
  class(environment_t)    , target, intent(in)    :: environment
  integer(ip)                     , intent(in)    :: field_id
  character(len=*)                , intent(in)    :: norm_type
  integer(ip)           , optional, intent(in)    :: exponent

  call this%free()
  this%fe_space    => fe_space
  this%environment => environment
  this%field_id = field_id
  select case(norm_type)
     case(l1_norm)
     allocate(l1_norm_t :: this%norm)
     case(l2_norm)
     allocate(l2_norm_t :: this%norm)
     case(h1_seminorm)
     allocate(h1_seminorm_t :: this%norm)
     case default
     check(.false.)
  end select

  call memalloc(fe_space%get_max_number_quadrature_points(),1,this%quadrature_points_values,__FILE__,__LINE__)
  allocate(this%quadrature_points_gradients(fe_space%get_max_number_quadrature_points(),1))

end subroutine error_norm_scalar_create

!====================================================================================================
subroutine error_norm_scalar_free(this)
  implicit none
  class(error_norm_scalar_t), intent(inout) :: this

  if(associated(this%fe_space))    nullify(this%fe_space)
  if(associated(this%environment)) nullify(this%environment)
  if(allocated(this%norm))         deallocate(this%norm)
  this%field_id = 0
  if (allocated(this%quadrature_points_values)) then
  call memfree(this%quadrature_points_values,__FILE__,__LINE__)
  end if
  if (allocated(this%quadrature_points_gradients)) deallocate(this%quadrature_points_gradients)
end subroutine error_norm_scalar_free

!====================================================================================================
function error_norm_scalar_compute_values(this,fe_function,scalar_function,time) result(norm)
  implicit none
  class(error_norm_scalar_t), intent(inout) :: this
  type(fe_function_t)       , intent(in)    :: fe_function
  class(scalar_function_t)  , intent(in)    :: scalar_function
  real(rp)   , optional     , intent(in)    :: time
  ! Result
  real(rp) :: norm

  ! Locals
  real(rp)                                     :: time_(1)
  type(cell_fe_function_scalar_t)              :: fe_unknown
  type(fe_iterator_t)                          :: fe_iterator
  type(fe_accessor_t)                          :: fe
  type(quadrature_t)             , pointer     :: quad
  type(fe_map_t)                 , pointer     :: fe_map
  type(point_t)                  , pointer     :: quadrature_coordinates(:)
  real(rp)                       , pointer     :: unknown_values_at_qpoints(:)
  integer(ip)                                  :: number_qpoints

  if(this%environment%am_i_l1_task()) then
     ! Initialize
     norm = 0.0_rp
     fe_iterator = this%fe_space%create_fe_iterator()
     call fe_iterator%current(fe)
     call this%fe_space%initialize_fe_integration()
     call this%fe_space%create_cell_fe_function(this%field_id,fe_unknown)
     if ( present(time) ) time_(1) = time

     call fe_iterator%current(fe)
     quad                   => fe%get_quadrature()
     fe_map                 => fe%get_fe_map()
     number_qpoints         = quad%get_number_quadrature_points()

     do while ( .not. fe_iterator%has_finished() )

        ! Get current FE
        call fe_iterator%current(fe)

        ! Update FE-integration related data structures
        call fe%update_integration()

        ! Get quadrature coordinates to evaluate boundary value
        quadrature_coordinates => fe_map%get_quadrature_coordinates()          

        ! Update cell FE function
        call fe%update_cell_fe_function(fe_unknown,fe_function)

        ! Evaluate function on quadrature points
        if(present(time)) then
           call scalar_function%get_values_set_space_time(quadrature_coordinates,time_, &
                                                          this%quadrature_points_values)
        else
           call scalar_function%get_values_set_space(quadrature_coordinates, &
                                                     this%quadrature_points_values(:,1))
        end if

        unknown_values_at_qpoints => fe_unknown%get_quadrature_points_values()
        this%quadrature_points_values(1:number_qpoints,1) =    &
             this%quadrature_points_values(1:number_qpoints,1) &
             - unknown_values_at_qpoints(1:number_qpoints)

        ! Compute element contribution to the global norm
        call this%norm%compute(quad,fe_map,this%quadrature_points_values(:,1),norm)

        call fe_iterator%next()
     end do

     ! Free
     call fe_unknown%free()

     ! All reduce
     call this%environment%l1_sum(norm)

     ! Finalize norm
     norm = this%norm%finalize(norm)
  end if
end function error_norm_scalar_compute_values

!====================================================================================================
function error_norm_scalar_compute_gradients(this,fe_function,vector_function,time) result(norm)
  implicit none
  class(error_norm_scalar_t), intent(inout) :: this
  type(fe_function_t)       , intent(in)    :: fe_function
  class(vector_function_t)  , intent(in)    :: vector_function
  real(rp)   , optional     , intent(in)    :: time
  ! Result
  real(rp) :: norm

  ! Locals
  real(rp)                                     :: time_(1)
  type(cell_fe_function_scalar_t)              :: fe_unknown
  type(fe_iterator_t)                          :: fe_iterator
  type(fe_accessor_t)                          :: fe
  type(quadrature_t)             , pointer     :: quad
  type(fe_map_t)                 , pointer     :: fe_map
  type(point_t)                  , pointer     :: quadrature_coordinates(:)
  type(vector_field_t)           , pointer     :: unknown_gradients_at_qpoints(:)
  integer(ip)                                  :: qpoint, number_qpoints

  if(this%environment%am_i_l1_task()) then
     ! Initialize
     norm = 0.0_rp
     fe_iterator = this%fe_space%create_fe_iterator()
     call fe_iterator%current(fe)
     call this%fe_space%initialize_fe_integration()
     call this%fe_space%create_cell_fe_function(this%field_id,fe_unknown)
     if ( present(time) ) time_(1) = time

     call fe_iterator%current(fe)
     quad                   => fe%get_quadrature()
     fe_map                 => fe%get_fe_map()
     number_qpoints         = quad%get_number_quadrature_points()

     do while ( .not. fe_iterator%has_finished() )

        ! Get current FE
        call fe_iterator%current(fe)

        ! Update FE-integration related data structures
        call fe%update_integration()

        ! Get quadrature coordinates to evaluate boundary value
        quadrature_coordinates => fe_map%get_quadrature_coordinates()          

        ! Update cell FE function
        call fe%update_cell_fe_function(fe_unknown,fe_function)

        ! Evaluate function on quadrature points
        if(present(time)) then
           call vector_function%get_values_set_space_time(quadrature_coordinates,time_, &
                                                          this%quadrature_points_gradients)
        else
           call vector_function%get_values_set_space(quadrature_coordinates, &
                                                     this%quadrature_points_gradients(:,1))
        end if

        unknown_gradients_at_qpoints => fe_unknown%get_quadrature_points_gradients()
        do qpoint = 1, number_qpoints
           this%quadrature_points_gradients(qpoint,1) =    &
                this%quadrature_points_gradients(qpoint,1) &
                - unknown_gradients_at_qpoints(qpoint)
        end do

        ! Compute element contribution to the global norm
        call this%norm%compute(quad,fe_map,this%quadrature_points_gradients(:,1),norm)

        call fe_iterator%next()
     end do

     ! Free
     call fe_unknown%free()

     ! All reduce
     call this%environment%l1_sum(norm)

     ! Finalize norm
     norm = this%norm%finalize(norm)
  end if
end function error_norm_scalar_compute_gradients

!====================================================================================================
subroutine error_norm_vector_create(this,fe_space,environment,field_id,norm_type,exponent)
  implicit none
  class(error_norm_vector_t)      , intent(inout) :: this
  class(serial_fe_space_t), target, intent(in)    :: fe_space
  class(environment_t)    , target, intent(in)    :: environment
  integer(ip)                     , intent(in)    :: field_id
  character(len=*)                , intent(in)    :: norm_type
  integer(ip)           , optional, intent(in)    :: exponent

  call this%free()
  this%fe_space    => fe_space
  this%environment => environment
  this%field_id = field_id
  select case(norm_type)
     case(l1_norm)
     allocate(l1_norm_t :: this%norm)
     case(l2_norm)
     allocate(l2_norm_t :: this%norm)
     case(h1_seminorm)
     allocate(h1_seminorm_t :: this%norm)
     case default
     check(.false.)
  end select

  allocate(this%quadrature_points_values(fe_space%get_max_number_quadrature_points(),1))
  allocate(this%quadrature_points_gradients(fe_space%get_max_number_quadrature_points(),1))

end subroutine error_norm_vector_create

!====================================================================================================
subroutine error_norm_vector_free(this)
implicit none
class(error_norm_vector_t), intent(inout) :: this

  if(associated(this%fe_space))    nullify(this%fe_space)
  if(associated(this%environment)) nullify(this%environment)
  if(allocated(this%norm))         deallocate(this%norm)
  this%field_id = 0
  if (allocated(this%quadrature_points_values)) deallocate(this%quadrature_points_values)
  if (allocated(this%quadrature_points_gradients)) deallocate(this%quadrature_points_gradients)
end subroutine error_norm_vector_free

!====================================================================================================
function error_norm_vector_compute_values(this,fe_function,vector_function,time) result(norm)
  implicit none
  class(error_norm_vector_t), intent(inout) :: this
  type(fe_function_t)       , intent(in)    :: fe_function
  class(vector_function_t)  , intent(in)    :: vector_function
  real(rp)   , optional     , intent(in)    :: time
  ! Result
  real(rp) :: norm

  ! Locals
  real(rp)                                     :: time_(1)
  type(cell_fe_function_vector_t)              :: fe_unknown
  type(fe_iterator_t)                          :: fe_iterator
  type(fe_accessor_t)                          :: fe
  type(quadrature_t)             , pointer     :: quad
  type(fe_map_t)                 , pointer     :: fe_map
  type(point_t)                  , pointer     :: quadrature_coordinates(:)
  type(vector_field_t)           , pointer     :: unknown_values_at_qpoints(:)
  integer(ip)                                  :: qpoint, number_qpoints

  if(this%environment%am_i_l1_task()) then
     ! Initialize
     norm = 0.0_rp
     fe_iterator = this%fe_space%create_fe_iterator()
     call fe_iterator%current(fe)
     call this%fe_space%initialize_fe_integration()
     call this%fe_space%create_cell_fe_function(this%field_id,fe_unknown)
     if ( present(time) ) time_(1) = time

     call fe_iterator%current(fe)
     quad                   => fe%get_quadrature()
     fe_map                 => fe%get_fe_map()
     number_qpoints         = quad%get_number_quadrature_points()

     do while ( .not. fe_iterator%has_finished() )

        ! Get current FE
        call fe_iterator%current(fe)

        ! Update FE-integration related data structures
        call fe%update_integration()

        ! Get quadrature coordinates to evaluate boundary value
        quadrature_coordinates => fe_map%get_quadrature_coordinates()          

        ! Update cell FE function
        call fe%update_cell_fe_function(fe_unknown,fe_function)

        ! Evaluate function on quadrature points
        if(present(time)) then
           call vector_function%get_values_set_space_time(quadrature_coordinates,time_, &
                                                          this%quadrature_points_values)
        else
           call vector_function%get_values_set_space(quadrature_coordinates, &
                                                     this%quadrature_points_values(:,1))
        end if

        unknown_values_at_qpoints => fe_unknown%get_quadrature_points_values()
        do qpoint = 1, number_qpoints
           this%quadrature_points_values(qpoint,1) =    &
                this%quadrature_points_values(qpoint,1) &
                - unknown_values_at_qpoints(qpoint)
        end do

        ! Compute element contribution to the global norm
        call this%norm%compute(quad,fe_map,this%quadrature_points_values(:,1),norm)

        call fe_iterator%next()
     end do

     ! Free
     call fe_unknown%free()

     ! All reduce
     call this%environment%l1_sum(norm)

     ! Finalize norm
     norm = this%norm%finalize(norm)
  end if
end function error_norm_vector_compute_values

!====================================================================================================
function error_norm_vector_compute_gradients(this,fe_function,tensor_function,time) result(norm)
  implicit none
  class(error_norm_vector_t), intent(inout) :: this
  type(fe_function_t)       , intent(in)    :: fe_function
  class(tensor_function_t)  , intent(in)    :: tensor_function
  real(rp)   , optional     , intent(in)    :: time
  ! Result
  real(rp) :: norm

  ! Locals
  real(rp)                                     :: time_(1)
  type(cell_fe_function_vector_t)              :: fe_unknown
  type(fe_iterator_t)                          :: fe_iterator
  type(fe_accessor_t)                          :: fe
  type(quadrature_t)             , pointer     :: quad
  type(fe_map_t)                 , pointer     :: fe_map
  type(point_t)                  , pointer     :: quadrature_coordinates(:)
  type(tensor_field_t)           , pointer     :: unknown_gradients_at_qpoints(:)
  integer(ip)                                  :: qpoint, number_qpoints

  if(this%environment%am_i_l1_task()) then
     ! Initialize
     norm = 0.0_rp
     fe_iterator = this%fe_space%create_fe_iterator()
     call fe_iterator%current(fe)
     call this%fe_space%initialize_fe_integration()
     call this%fe_space%create_cell_fe_function(this%field_id,fe_unknown)
     if ( present(time) ) time_(1) = time

     call fe_iterator%current(fe)
     quad                   => fe%get_quadrature()
     fe_map                 => fe%get_fe_map()
     number_qpoints         = quad%get_number_quadrature_points()

     do while ( .not. fe_iterator%has_finished() )

        ! Get current FE
        call fe_iterator%current(fe)

        ! Update FE-integration related data structures
        call fe%update_integration()

        ! Get quadrature coordinates to evaluate boundary value
        quadrature_coordinates => fe_map%get_quadrature_coordinates()          

        ! Update cell FE function
        call fe%update_cell_fe_function(fe_unknown,fe_function)

        ! Evaluate function on quadrature points
        if(present(time)) then
           call tensor_function%get_values_set_space_time(quadrature_coordinates,time_, &
                                                          this%quadrature_points_gradients)
        else
           call tensor_function%get_values_set_space(quadrature_coordinates, &
                                                     this%quadrature_points_gradients(:,1))
        end if
        
        unknown_gradients_at_qpoints => fe_unknown%get_quadrature_points_gradients()        
        do qpoint = 1, number_qpoints
           this%quadrature_points_gradients(qpoint,1) =    &
                this%quadrature_points_gradients(qpoint,1) &
                - unknown_gradients_at_qpoints(qpoint)
        end do

        ! Compute element contribution to the global norm
        call this%norm%compute(quad,fe_map,this%quadrature_points_gradients(:,1),norm)

        call fe_iterator%next()
     end do

     ! Free
     call fe_unknown%free()

     ! All reduce
     call this%environment%l1_sum(norm)

     ! Finalize norm
     norm = this%norm%finalize(norm)
  end if
end function error_norm_vector_compute_gradients

end module error_norms_names
