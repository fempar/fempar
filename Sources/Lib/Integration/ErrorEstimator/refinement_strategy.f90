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
module refinement_strategy_names
  use types_names
  use fe_space_names
  use p4est_triangulation_names
  use std_vector_integer_ip_names
  use error_estimator_names
  use FPL
  
  implicit none
# include "debug.i90"
  private
  
  character(len=*), parameter :: error_objective_key         = 'error_objective'
  character(len=*), parameter :: objective_tolerance_key     = 'objective_tolerance'
  character(len=*), parameter :: max_num_mesh_iterations_key = 'max_num_mesh_iterations'
  
  type :: error_objective_refinement_strategy_t
    private
    class(error_estimator_t), pointer :: error_estimator => NULL()
    real(rp)                          :: error_objective
    real(rp)                          :: objective_tolerance
    integer(ip)                       :: current_mesh_iteration
    integer(ip)                       :: max_num_mesh_iterations
   contains
    procedure :: create                  => eors_create
    procedure :: update_refinement_flags => eors_update_refinement_flags
    procedure :: has_finished_refinement => eors_has_finished_refinement
  end type error_objective_refinement_strategy_t
  
  public :: error_objective_refinement_strategy_t
  
  public :: error_objective_key, objective_tolerance_key, max_num_mesh_iterations_key
  
contains
  
  subroutine eors_create(this,error_estimator,parameter_list)
    class(error_objective_refinement_strategy_t)         , intent(inout) :: this
    class(error_estimator_t)                    , target , intent(in)    :: error_estimator
    type(parameterlist_t)                                , intent(in)    :: parameter_list
    integer(ip) :: FPLerror
    this%error_estimator => error_estimator
    assert(parameter_list%isPresent(error_objective_key))
    assert(parameter_list%isAssignable(error_objective_key,this%error_objective))
    FPLerror = parameter_list%get(key = error_objective_key, value = this%error_objective)
    assert(FPLerror==0)
    assert(parameter_list%isPresent(objective_tolerance_key))
    assert(parameter_list%isAssignable(objective_tolerance_key,this%objective_tolerance))
    FPLerror = parameter_list%get(key = objective_tolerance_key, value = this%objective_tolerance)
    assert(FPLerror==0)
    if ( parameter_list%isPresent(max_num_mesh_iterations_key) ) then
      FPLerror = parameter_list%get(key = max_num_mesh_iterations_key, value = this%max_num_mesh_iterations)
      assert(FPLerror==0)
    else
      this%max_num_mesh_iterations = 100
    end if
    this%current_mesh_iteration = 0
  end subroutine eors_create
  
  subroutine eors_update_refinement_flags(this,refinement_and_coarsening_flags)
    class(error_objective_refinement_strategy_t), intent(inout) :: this
    type(std_vector_integer_ip_t)               , intent(inout) :: refinement_and_coarsening_flags
    integer(ip), pointer :: refinement_and_coarsening_flags_entries(:)
    integer(ip)          :: i, num_of_entries
    real(rp)   , pointer :: sq_local_estimate_entries(:)
    real(rp)             :: sq_error_upper_bound, sq_error_lower_bound
    refinement_and_coarsening_flags_entries => refinement_and_coarsening_flags%get_pointer()
    sq_local_estimate_entries => this%error_estimator%get_sq_local_estimate_entries()
    num_of_entries        = size(sq_local_estimate_entries)
    sq_error_upper_bound  = ( this%error_objective * ( 1.0_rp + this%objective_tolerance ) ) ** 2.0_rp
    sq_error_lower_bound  = ( this%error_objective * ( 1.0_rp - this%objective_tolerance ) ) ** 2.0_rp
    do i = 1, num_of_entries
      if ( sq_local_estimate_entries(i) > sq_error_upper_bound ) then
        refinement_and_coarsening_flags_entries(i) = refinement
      else if ( sq_local_estimate_entries(i) < sq_error_lower_bound ) then
        !refinement_and_coarsening_flags_entries(i) = coarsening
      else
        refinement_and_coarsening_flags_entries(i) = do_nothing
      end if
    end do
    this%current_mesh_iteration = this%current_mesh_iteration + 1
  end subroutine eors_update_refinement_flags
  
  function eors_has_finished_refinement(this)
    class(error_objective_refinement_strategy_t), intent(inout) :: this
    logical :: eors_has_finished_refinement
    real(rp), pointer :: sq_local_estimate_entries(:)
    real(rp)          :: max_local_estimate
    real(rp)          :: sq_error_upper_bound
    sq_local_estimate_entries => this%error_estimator%get_sq_local_estimate_entries()
    max_local_estimate = maxval(sq_local_estimate_entries)
    sq_error_upper_bound  = ( this%error_objective * ( 1.0_rp + this%objective_tolerance ) ) ** 2.0_rp
    eors_has_finished_refinement = ( ( this%current_mesh_iteration > 1             .and. &
                                       max_local_estimate < sq_error_upper_bound )  .or. & 
                                     ( this%current_mesh_iteration > this%max_num_mesh_iterations ) )
    if ( this%current_mesh_iteration > this%max_num_mesh_iterations ) &
      write(*,*) 'Error objective mesh refinement strategy exceeded the maximum number of iterations'
  end function eors_has_finished_refinement
  
end module refinement_strategy_names