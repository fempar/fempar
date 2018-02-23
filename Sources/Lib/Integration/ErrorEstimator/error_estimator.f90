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
  use std_vector_real_rp_names
  use FPL
  
  implicit none
# include "debug.i90"
  private
  
  
  type, abstract :: explicit_error_estimator_t
    class(serial_fe_space_t)    , pointer     :: fe_space
    !class(refinement_strategy_t), allocatable :: refinement_strategy
    type(std_vector_real_rp_t)  , allocatable :: local_estimates
    real(rp)                                  :: global_estimate
    type(std_vector_real_rp_t)  , allocatable :: local_true_errors
    real(rp)                                  :: global_true_error
    type(std_vector_real_rp_t)  , allocatable :: local_effectivities
    real(rp)                                  :: global_effectivity
    integer(ip)                               :: output_mode
    integer(ip)                               :: logical_unit
   contains
    procedure (compute_local_estimates_interface)  , deferred :: compute_local_estimates
    procedure (compute_local_true_errors_interface), deferred :: compute_local_true_errors
    procedure                           :: create                           => eee_create
    procedure                           :: free                             => eee_free
    !procedure, non_overridable          :: setup_adaptive_strategy          => eee_setup_adaptive_strategy
    !procedure, non_overridable          :: has_finished_refinement          => eee_has_finished_adaptive_refinement
    !procedure, non_overridable          :: update_refinement_flags          => eee_update_refinement_flags
    !procedure, non_overridable          :: compute_global_estimate          => eee_compute_global_estimate
    !procedure, non_overridable          :: compute_global_true_error        => eee_compute_global_true_error
    !procedure, non_overridable          :: compute_local_effectivities      => eee_get_local_effectivities
    !procedure, non_overridable          :: compute_global_effectivity       => eee_get_global_effectivity
    procedure, non_overridable          :: get_local_estimates              => eee_get_local_estimates
    procedure, non_overridable          :: get_global_estimate              => eee_get_global_estimate
    procedure, non_overridable          :: get_local_true_errors            => eee_get_local_true_errors
    procedure, non_overridable          :: get_global_true_error            => eee_get_global_true_error
    procedure, non_overridable          :: get_local_effectivities          => eee_get_local_effectivities
    procedure, non_overridable          :: get_global_effectivity           => eee_get_global_effectivity
    !procedure, non_overridable          :: print_estimator_history_header   => eee_print_estimator_history_header
    !procedure, non_overridable          :: print_estimator_history_new_line => eee_print_estimator_history_new_line
    !procedure, non_overridable          :: print_estimator_history_footer   => eee_print_estimator_history_footer
  end type explicit_error_estimator_t
  
  abstract interface

    subroutine compute_local_estimates_interface(this)
      import :: explicit_error_estimator_t
      class(explicit_error_estimator_t), intent(inout) :: this
    end subroutine compute_local_estimates_interface

    subroutine compute_local_true_errors_interface(this)
      import :: explicit_error_estimator_t
      class(explicit_error_estimator_t), intent(inout) :: this
    end subroutine compute_local_true_errors_interface

  end interface
  
  public :: explicit_error_estimator_t

contains

 subroutine eee_create ( this, fe_space, parameter_list )
   implicit none
   class(explicit_error_estimator_t)        , intent(inout) :: this
   class(serial_fe_space_t)         , target, intent(in)    :: fe_space
   type(parameterlist_t)                    , intent(in)    :: parameter_list
   call this%free()
   this%fe_space => fe_space
 end subroutine eee_create

 subroutine eee_free ( this )
   implicit none
   class(explicit_error_estimator_t), intent(inout) :: this
   nullify(this%fe_space)
   this%output_mode        = 0_ip
   this%logical_unit       = 0_ip
   this%global_estimate    = 0.0_rp
   this%global_true_error  = 0.0_rp
   this%global_effectivity = 0.0_rp
   call this%local_estimates%free()
   call this%local_true_errors%free()
   call this%local_effectivities%free()
 end subroutine eee_free
 
 function eee_get_local_estimates ( this )
   implicit none
   class(explicit_error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: eee_get_local_estimates
   eee_get_local_estimates => this%local_estimates
 end function eee_get_local_estimates
 
 function eee_get_global_estimate ( this )
   implicit none
   class(explicit_error_estimator_t), intent(inout) :: this
   real(rp) :: eee_get_global_estimate
   eee_get_global_estimate = this%global_estimate
 end function eee_get_global_estimate
 
 function eee_get_local_true_errors ( this )
   implicit none
   class(explicit_error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: eee_get_local_true_errors
   eee_get_local_true_errors => this%local_true_errors
 end function eee_get_local_true_errors
 
 function eee_get_global_true_error ( this )
   implicit none
   class(explicit_error_estimator_t), intent(inout) :: this
   real(rp) :: eee_get_global_true_error
   eee_get_global_true_error = this%global_true_error
 end function eee_get_global_true_error
 
 function eee_get_local_effectivities ( this )
   implicit none
   class(explicit_error_estimator_t), target, intent(inout) :: this
   type(std_vector_real_rp_t), pointer :: eee_get_local_effectivities
   eee_get_local_effectivities => this%local_effectivities
 end function eee_get_local_effectivities
 
 function eee_get_global_effectivity ( this )
   implicit none
   class(explicit_error_estimator_t), intent(inout) :: this
   real(rp) :: eee_get_global_effectivity
   eee_get_global_effectivity = this%global_effectivity
 end function eee_get_global_effectivity

end module error_estimator_names