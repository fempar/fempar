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
  use serial_environment_names
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

  ! Parameters
  public :: mean_norm, l1_norm, l2_norm, lp_norm, linfty_norm, h1_seminorm
  public :: hdiv_seminorm, h1_norm, w1p_seminorm, w1p_norm, w1infty_seminorm, w1infty_norm

  ! Data types
  public :: error_norms_scalar_t !, error_norms_vector_t, error_norms_tensor_t

contains

  !===================================================================================================
  subroutine error_norms_scalar_create (this, fe_space, field_id)
    implicit none
    class(error_norms_scalar_t)     , intent(inout) :: this
    class(serial_fe_space_t), target, intent(in)    :: fe_space
    integer(ip)                     , intent(in)    :: field_id
    
    integer(ip) :: istat
    
    call this%free()
    this%field_id = field_id
    this%fe_space => fe_space
    call this%fe_space%create_cell_fe_function(field_id, this%cell_fe_function)
    
    call memalloc(this%fe_space%get_max_number_quadrature_points(), 1, this%work_array_values, __FILE__, __LINE__)
    allocate(this%work_array_gradients(this%fe_space%get_max_number_quadrature_points(),1), stat=istat)
    check(istat==0)
  end subroutine error_norms_scalar_create
  
  !===================================================================================================
  subroutine error_norms_scalar_free (this)
    implicit none
    class(error_norms_scalar_t), intent(inout) :: this
    integer(ip) :: istat
    this%field_id = 0
    nullify(this%fe_space)
    call this%cell_fe_function%free()
    if (allocated(this%work_array_values))    call memfree (this%work_array_values   , __FILE__, __LINE__)
    if (allocated(this%work_array_gradients)) then
      deallocate(this%work_array_gradients, stat=istat)
      check(istat==0)
    end if
  end subroutine error_norms_scalar_free
  
  !===================================================================================================
  function error_norms_scalar_compute (this, exact_solution, fe_function, norm_type, exponent, time ) result(norm)
    implicit none
    class(error_norms_scalar_t), intent(inout) :: this
    class(scalar_function_t)   , intent(in)    :: exact_solution
    type(fe_function_t)        , intent(in)    :: fe_function
    character(*)               , intent(in)    :: norm_type
    integer(ip), optional      , intent(in)    :: exponent
    real(rp)   , optional      , intent(in)    :: time
    real(rp)                                   :: norm

    type(serial_environment_t) :: environment
    
    ! Local variables
    real(rp)            :: exponent_
    real(rp)            :: time_(1)
    type(fe_iterator_t) :: fe_iterator
    type(fe_accessor_t) :: fe
    real(rp)            :: values_cell_contribution, values_norm
    real(rp)            :: gradients_cell_contribution, gradients_norm
    
    assert ( error_norm_is_supported(norm_type) )
    assert ( trim(norm_type) /= hdiv_seminorm)
    if ( environment%am_i_l1_task() ) then 
       exponent_ = error_norm_determine_exponent(norm_type, exponent)
    
       values_norm    = 0.0_rp
       gradients_norm = 0.0_rp
       fe_iterator = this%fe_space%create_fe_iterator()
       do while ( .not. fe_iterator%has_finished() ) 
         call fe_iterator%current(fe)
         if ( fe%is_local() ) then
           call fe%update_integration()
           call fe%update_cell_fe_function(fe_function, this%cell_fe_function)
           call this%compute_cell_contribution( fe, &
                                                exact_solution, & 
                                                norm_type, &
                                                exponent_, &
                                                values_cell_contribution, &
                                                gradients_cell_contribution, & 
                                                time)
           call update_norm(norm_type, values_cell_contribution   , values_norm)
           call update_norm(norm_type, gradients_cell_contribution, gradients_norm)
         end if  
         call fe_iterator%next() 
       end do
       call finalize_norm(environment, norm_type, exponent_, values_norm, gradients_norm, norm)
    end if
  end function error_norms_scalar_compute
  
  subroutine error_norms_scalar_compute_cell_contribution (this, &
                                                           fe, & 
                                                           exact_solution, & 
                                                           norm_type, & 
                                                           exponent, &
                                                           values_cell_contribution, & 
                                                           gradients_cell_contribution, &
                                                           time )
    implicit none
    class(error_norms_scalar_t), intent(inout) :: this
    type(fe_accessor_t)        , intent(in)    :: fe
    class(scalar_function_t)   , intent(in)    :: exact_solution
    character(*)               , intent(in)    :: norm_type
    real(rp)                   , intent(in)    :: exponent
    real(rp)                   , intent(out)   :: values_cell_contribution
    real(rp)                   , intent(out)   :: gradients_cell_contribution
    real(rp)   , optional      , intent(in)    :: time 

    ! Locals 
    logical                       :: norm_requires_values
    logical                       :: norm_requires_gradients
    type(fe_map_t)    , pointer   :: fe_map 
    type(quadrature_t), pointer   :: quadrature 
    type(point_t)     , pointer   :: coordinates(:)
    real(rp)                      :: time_(1) 
    integer(ip)                   :: q_point, n_q_points
    integer(ip)                   :: idime

    real(rp)            , pointer :: fe_function_values(:)
    type(vector_field_t), pointer :: fe_function_gradients(:)
    real(rp)                      :: gradient_component

    if ( present(time) ) time_ = time

    norm_requires_values    = error_norm_requires_values(norm_type)
    norm_requires_gradients = error_norm_requires_gradients(norm_type) 
    fe_map                  => fe%get_fe_map()
    quadrature              => fe%get_quadrature()
    coordinates             => fe_map%get_quadrature_coordinates()
    n_q_points              = quadrature%get_number_quadrature_points() 

    if (norm_requires_values) then
      ! First evaluate exact solution at the quadrature points.
      if(present(time)) then
         call exact_solution%get_values_set_space_time(coordinates, &
                                                       time_, &
                                                       this%work_array_values(1:n_q_points,1:1))
      else
         call exact_solution%get_values_set_space(coordinates, &
                                                  this%work_array_values(1:n_q_points,1))
      end if
      ! Then, subtract fe_function.
      fe_function_values    => this%cell_fe_function%get_quadrature_points_values()
      do q_point=1, n_q_points
        this%work_array_values(q_point,1) = this%work_array_values(q_point,1) - fe_function_values(q_point) 
      end do
    end if
    
    ! Do the same for gradients, if required by norm type 
    if (norm_requires_gradients) then  
      if(present(time)) then
         call exact_solution%get_gradients_set_space_time(coordinates, &
                                                       time_, &
                                                       this%work_array_gradients(1:n_q_points,1:1))
      else
         call exact_solution%get_gradients_set_space(coordinates, &
                                                     this%work_array_gradients(1:n_q_points,1))
      end if
      fe_function_gradients => this%cell_fe_function%get_quadrature_points_gradients()
      do q_point=1, n_q_points
        this%work_array_gradients(q_point,1) = this%work_array_gradients(q_point,1) - fe_function_gradients(q_point) 
      end do
    end if


    values_cell_contribution = 0.0_rp 
    select case ( trim(norm_type) )
    case (mean_norm)
      do q_point=1, n_q_points
        values_cell_contribution = values_cell_contribution + & 
                                    this%work_array_values(q_point,1)*fe_map%get_det_jacobian(q_point)*quadrature%get_weight(q_point)
      end do
    case (l1_norm)
      do q_point=1, n_q_points
        values_cell_contribution = values_cell_contribution + & 
                                    abs(this%work_array_values(q_point,1))*fe_map%get_det_jacobian(q_point) * quadrature%get_weight(q_point)
      end do
    case (l2_norm,h1_norm)
      do q_point=1, n_q_points
        values_cell_contribution = values_cell_contribution + & 
                                    (this%work_array_values(q_point,1)*this%work_array_values(q_point,1))*&
                                    fe_map%get_det_jacobian(q_point) * quadrature%get_weight(q_point)
      end do
    case (lp_norm,w1p_norm)
      do q_point=1, n_q_points
        values_cell_contribution = values_cell_contribution + & 
                                    (abs(this%work_array_values(q_point,1))**exponent)*&
                                   fe_map%get_det_jacobian(q_point) * quadrature%get_weight(q_point)
      end do
    case (linfty_norm,w1infty_norm)
      do q_point=1, n_q_points
        values_cell_contribution = max(values_cell_contribution, abs(this%work_array_values(q_point,1)))
      end do
    end select

    gradients_cell_contribution = 0.0_rp 
    select case ( trim(norm_type) )
    case (h1_norm, h1_seminorm)
      do q_point=1, n_q_points
        gradients_cell_contribution = gradients_cell_contribution + & 
                                      (this%work_array_gradients(q_point,1)*this%work_array_gradients(q_point,1))*&
                                       fe_map%get_det_jacobian(q_point) * quadrature%get_weight(q_point)
      end do
    case (w1p_norm, w1p_seminorm)
      do q_point=1, n_q_points
        gradients_cell_contribution = gradients_cell_contribution + & 
                                      (this%work_array_gradients(q_point,1)*this%work_array_gradients(q_point,1))**(exponent/2.0_rp)*&
                                       fe_map%get_det_jacobian(q_point) * quadrature%get_weight(q_point)
      end do
    case (w1infty_norm,w1infty_seminorm)
      do q_point=1, n_q_points
        do idime=1, SPACE_DIM
           gradient_component          = this%work_array_gradients(q_point,1)%get(idime)
           gradients_cell_contribution = max(gradients_cell_contribution, abs(gradient_component))
        end do 
      end do
    end select

  end subroutine error_norms_scalar_compute_cell_contribution 
 
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
  
  !function update_norm (norm_type, cell_contribution, norm)
  !  implicit none
  !  character(*), intent(in)    :: norm_type
  !  real(rp)    , intent(in)    :: cell_contribution
  !  real(rp)    , intent(inout) :: norm
  !  select case ( trim(norm_type) )
  !  case (linfty_norm,W1infty_norm,w1infty_seminorm)
  !     norm = max(cell_contribution,norm)
  !  case default:
  !    check(.false.)
  !  end select case
  !end function  update_norm

end module error_norms_names
