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
module error_estimator_names
  use types_names
  use fe_space_names
  use triangulation_names
  use std_vector_real_rp_names
  use environment_names
  use FPL
  
  implicit none
# include "debug.i90"
  private
  
  
  type, abstract :: error_estimator_t
    private
    class(serial_fe_space_t)  , pointer :: fe_space => NULL()
    type(std_vector_real_rp_t)          :: sq_local_estimates
    real(rp)                            :: global_estimate
    type(std_vector_real_rp_t)          :: sq_local_true_errors
    real(rp)                            :: global_true_error
    type(std_vector_real_rp_t)          :: local_effectivities
    real(rp)                            :: global_effectivity
    integer(ip)                         :: output_mode
    integer(ip)                         :: logical_unit
   contains
    procedure (compute_local_estimates_interface)  , deferred :: compute_local_estimates
    procedure (compute_local_true_errors_interface), deferred :: compute_local_true_errors
    procedure (get_error_norm_exponent_interface)  , deferred :: get_error_norm_exponent
    procedure                           :: create                           => ee_create
    procedure                           :: free                             => ee_free
    procedure, non_overridable          :: compute_global_estimate          => ee_compute_global_estimate
    procedure, non_overridable          :: compute_global_true_error        => ee_compute_global_true_error
    procedure, non_overridable          :: compute_local_effectivities      => ee_compute_local_effectivities
    procedure, non_overridable          :: compute_global_effectivity       => ee_compute_global_effectivity
    procedure, non_overridable          :: get_fe_space                     => ee_get_fe_space
    procedure, non_overridable          :: get_environment                  => ee_get_environment
    procedure, non_overridable          :: get_sq_local_estimates           => ee_get_sq_local_estimates
    procedure, non_overridable          :: get_sq_local_estimate_entries    => ee_get_sq_local_estimate_entries
    procedure, non_overridable          :: get_global_estimate              => ee_get_global_estimate
    procedure, non_overridable          :: get_sq_local_true_errors         => ee_get_sq_local_true_errors
    procedure, non_overridable          :: get_sq_local_true_error_entries  => ee_get_sq_local_true_error_entries
    procedure, non_overridable          :: get_global_true_error            => ee_get_global_true_error
    procedure, non_overridable          :: get_local_effectivities          => ee_get_local_effectivities
    procedure, non_overridable          :: get_local_effectivity_entries    => ee_get_local_effectivity_entries
    procedure, non_overridable          :: get_global_effectivity           => ee_get_global_effectivity
    !procedure, non_overridable          :: print_estimator_history_header   => ee_print_estimator_history_header
    !procedure, non_overridable          :: print_estimator_history_new_line => ee_print_estimator_history_new_line
    !procedure, non_overridable          :: print_estimator_history_footer   => ee_print_estimator_history_footer
  end type error_estimator_t
  
  abstract interface

    subroutine compute_local_estimates_interface(this)
      import :: error_estimator_t
      class(error_estimator_t), intent(inout) :: this
    end subroutine compute_local_estimates_interface

    subroutine compute_local_true_errors_interface(this)
      import :: error_estimator_t
      class(error_estimator_t), intent(inout) :: this
    end subroutine compute_local_true_errors_interface

    function get_error_norm_exponent_interface(this)
      import :: error_estimator_t, rp
      class(error_estimator_t), intent(in) :: this
      real(rp) :: get_error_norm_exponent_interface
    end function get_error_norm_exponent_interface

  end interface
  
  public :: error_estimator_t, ee_create, ee_free

contains

 subroutine ee_create ( this, fe_space, parameter_list )
   implicit none
   class(error_estimator_t)        , intent(inout) :: this
   class(serial_fe_space_t), target, intent(in)    :: fe_space
   type(parameterlist_t)           , intent(in)    :: parameter_list
   call this%free()
   this%fe_space => fe_space
 end subroutine ee_create
 
 subroutine ee_free ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   nullify(this%fe_space)
   this%output_mode        = 0_ip
   this%logical_unit       = 0_ip
   this%global_estimate    = 0.0_rp
   this%global_true_error  = 0.0_rp
   this%global_effectivity = 0.0_rp
   call this%sq_local_estimates%free()
   call this%sq_local_true_errors%free()
   call this%local_effectivities%free()
 end subroutine ee_free
 
 subroutine ee_compute_global_estimate ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   class(environment_t), pointer :: environment
   this%global_estimate = sum(this%sq_local_estimates%get_pointer())
   environment => this%get_environment()
   if ( environment%get_l1_size() > 1 ) then
     call environment%l1_sum(this%global_estimate)
   end if
   this%global_estimate = this%global_estimate ** this%get_error_norm_exponent()
 end subroutine ee_compute_global_estimate
 
 subroutine ee_compute_global_true_error ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   class(environment_t), pointer :: environment
   this%global_true_error = sum(this%sq_local_true_errors%get_pointer())
   environment => this%get_environment()
   if ( environment%get_l1_size() > 1 ) then
     call environment%l1_sum(this%global_true_error)
   end if
   this%global_true_error = this%global_true_error ** this%get_error_norm_exponent()
 end subroutine ee_compute_global_true_error
 
 subroutine ee_compute_local_effectivities ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   class(environment_t)  , pointer :: environment
   class(triangulation_t), pointer :: triangulation
   real(rp)              , pointer :: sq_local_estimate_entries(:)
   real(rp)              , pointer :: sq_local_true_error_entries(:)
   real(rp)              , pointer :: local_effectivity_entries(:)
   integer(ip)                     :: icell, num_local_cells
   real(rp)                        :: norm_exponent
   environment => this%get_environment()
   if ( environment%am_i_l1_task() ) then
     triangulation               => this%fe_space%get_triangulation()
     num_local_cells             =  triangulation%get_num_local_cells()
     assert ( this%sq_local_estimates%size()   == num_local_cells )
     assert ( this%sq_local_true_errors%size() == num_local_cells )
     sq_local_estimate_entries   => this%sq_local_estimates%get_pointer()
     sq_local_true_error_entries => this%sq_local_true_errors%get_pointer()
     call this%local_effectivities%resize(0)
     call this%local_effectivities%resize(num_local_cells,0.0_rp)
     local_effectivity_entries   => this%local_effectivities%get_pointer()
     norm_exponent = this%get_error_norm_exponent()
     do icell = 1, num_local_cells
       local_effectivity_entries(icell) =                   &
         ( sq_local_estimate_entries(icell)   ** norm_exponent ) / &
         ( sq_local_true_error_entries(icell) ** norm_exponent )
     end do
   end if
 end subroutine ee_compute_local_effectivities
 
 subroutine ee_compute_global_effectivity ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   class(environment_t), pointer :: environment
   this%global_effectivity = this%global_estimate / this%global_true_error
 end subroutine ee_compute_global_effectivity
 
 function ee_get_fe_space ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   class(serial_fe_space_t), pointer :: ee_get_fe_space
   ee_get_fe_space => this%fe_space
 end function ee_get_fe_space
 
 function ee_get_environment ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   class(environment_t), pointer :: ee_get_environment
   ee_get_environment => this%fe_space%get_environment()
 end function ee_get_environment
 
 function ee_get_sq_local_estimates ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: ee_get_sq_local_estimates
   ee_get_sq_local_estimates => this%sq_local_estimates
 end function ee_get_sq_local_estimates
 
 function ee_get_sq_local_estimate_entries ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   real(rp), pointer :: ee_get_sq_local_estimate_entries(:)
   ee_get_sq_local_estimate_entries => this%sq_local_estimates%get_pointer()
 end function ee_get_sq_local_estimate_entries
 
 function ee_get_global_estimate ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   real(rp) :: ee_get_global_estimate
   ee_get_global_estimate = this%global_estimate
 end function ee_get_global_estimate
 
 function ee_get_sq_local_true_errors ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: ee_get_sq_local_true_errors
   ee_get_sq_local_true_errors => this%sq_local_true_errors
 end function ee_get_sq_local_true_errors
 
 function ee_get_sq_local_true_error_entries ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   real(rp), pointer :: ee_get_sq_local_true_error_entries(:)
   ee_get_sq_local_true_error_entries => this%sq_local_true_errors%get_pointer()
 end function ee_get_sq_local_true_error_entries
 
 function ee_get_global_true_error ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   real(rp) :: ee_get_global_true_error
   ee_get_global_true_error = this%global_true_error
 end function ee_get_global_true_error
 
 function ee_get_local_effectivities ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: ee_get_local_effectivities
   ee_get_local_effectivities => this%local_effectivities
 end function ee_get_local_effectivities
 
 function ee_get_local_effectivity_entries ( this )
   implicit none
   class(error_estimator_t), target, intent(inout) :: this
   real(rp), pointer :: ee_get_local_effectivity_entries(:)
   ee_get_local_effectivity_entries => this%local_effectivities%get_pointer()
 end function ee_get_local_effectivity_entries
 
 function ee_get_global_effectivity ( this )
   implicit none
   class(error_estimator_t), intent(inout) :: this
   real(rp) :: ee_get_global_effectivity
   ee_get_global_effectivity = this%global_effectivity
 end function ee_get_global_effectivity

end module error_estimator_names