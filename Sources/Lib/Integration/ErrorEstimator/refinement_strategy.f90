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
  use environment_names
  use triangulation_names
  use fe_space_names
  use p4est_triangulation_names
  use std_vector_integer_ip_names
  use error_estimator_names
  use FPL
  
  implicit none
# include "debug.i90"
  private
  
  character(len=*), parameter :: num_uniform_refinements_key        = 'num_uniform_refinements'
  
  character(len=*), parameter :: error_objective_key                = 'error_objective'
  character(len=*), parameter :: objective_tolerance_key            = 'objective_tolerance'
  character(len=*), parameter :: max_num_mesh_iterations_key        = 'max_num_mesh_iterations'
  
  character(len=*), parameter :: refinement_fraction_key            = 'refinement_fraction'
  character(len=*), parameter :: coarsening_fraction_key            = 'coarsening_fraction'
  character(len=*), parameter :: threshold_tolerance_key            = 'threshold_tolerance_key'
  
  type, abstract :: refinement_strategy_t
    private
    class(error_estimator_t), pointer :: error_estimator => NULL()
    integer(ip)                       :: current_mesh_iteration
   contains
    procedure ( set_parameters_interface )         , deferred :: set_parameters
    procedure ( update_refinement_flags_interface ), deferred :: update_refinement_flags
    procedure ( has_finished_refinement_interface ), deferred :: has_finished_refinement
    procedure :: create                     => refinement_strategy_create
    procedure :: free                       => refinement_strategy_free
    procedure :: get_error_estimator        => refinement_strategy_get_error_estimator
    procedure :: get_current_mesh_iteration => refinement_strategy_get_current_mesh_iteration
  end type refinement_strategy_t
  
  abstract interface
  
    subroutine set_parameters_interface(this,parameter_list)
      import :: refinement_strategy_t, parameterlist_t
      class(refinement_strategy_t) , intent(inout) :: this
      type(parameterlist_t)        , intent(in)    :: parameter_list
    end subroutine set_parameters_interface
  
    subroutine update_refinement_flags_interface(this,refinement_and_coarsening_flags)
      import :: refinement_strategy_t, std_vector_integer_ip_t
      class(refinement_strategy_t) , intent(inout) :: this
      type(std_vector_integer_ip_t), intent(inout) :: refinement_and_coarsening_flags
    end subroutine update_refinement_flags_interface
  
    function has_finished_refinement_interface(this)
      import :: refinement_strategy_t
      class(refinement_strategy_t), intent(inout) :: this
      logical :: has_finished_refinement_interface
    end function has_finished_refinement_interface
  
  end interface
  
  type, extends(refinement_strategy_t) :: uniform_refinement_strategy_t
    private
    integer(ip)                       :: num_uniform_refinements
   contains
    procedure :: set_parameters             => urs_set_parameters
    procedure :: update_refinement_flags    => urs_update_refinement_flags
    procedure :: has_finished_refinement    => urs_has_finished_refinement
  end type uniform_refinement_strategy_t
  
  type, extends(refinement_strategy_t) :: error_objective_refinement_strategy_t
    private
    real(rp)                          :: error_objective
    real(rp)                          :: objective_tolerance
    integer(ip)                       :: max_num_mesh_iterations
   contains
    procedure :: set_parameters             => eors_set_parameters
    procedure :: update_refinement_flags    => eors_update_refinement_flags
    procedure :: has_finished_refinement    => eors_has_finished_refinement
  end type error_objective_refinement_strategy_t
  
  ! Implements Algorithm for calculating coarsening and refinement thresholds
  ! in Fig 5. of the following paper:
  ! Algorithms and data structures for massively parallel generic adaptive finite element codes
  ! W Bangerth, C Burstedde, T Heister, M Kronbichler
  ! ACM Transactions on Mathematical Software 38 (2)
  type, extends(refinement_strategy_t) :: fixed_fraction_refinement_strategy_t
    private
    ! Refine all cells s.t. #{e_i > \theta_r} \approx this%refinement_fraction * N_cells
    ! Coarsen all cells s.t. #{e_i < \theta_c} \approx this%coarsening_fraction * N_cells 
    real(rp)                          :: refinement_fraction
    real(rp)                          :: coarsening_fraction
    real(rp)                          :: threshold_tolerance
    real(rp)                          :: sq_refinement_threshold
    real(rp)                          :: sq_coarsening_threshold
    integer(ip)                       :: max_num_mesh_iterations
   contains
    procedure :: set_parameters             => ffrs_set_parameters
    procedure :: update_refinement_flags    => ffrs_update_refinement_flags
    procedure :: has_finished_refinement    => ffrs_has_finished_refinement
  end type fixed_fraction_refinement_strategy_t
  
  
  public :: num_uniform_refinements_key
  public :: error_objective_key, objective_tolerance_key, max_num_mesh_iterations_key
  public :: refinement_fraction_key, coarsening_fraction_key
  
  public :: refinement_strategy_t
  public :: uniform_refinement_strategy_t, error_objective_refinement_strategy_t, fixed_fraction_refinement_strategy_t
  
contains
  
  subroutine refinement_strategy_create(this,error_estimator,parameter_list)
    class(refinement_strategy_t)         , intent(inout) :: this
    class(error_estimator_t)    , target , intent(in)    :: error_estimator
    type(parameterlist_t)                , intent(in)    :: parameter_list
    call this%free()
    this%error_estimator => error_estimator
    call this%set_parameters(parameter_list)
    this%current_mesh_iteration = 0
  end subroutine refinement_strategy_create
  
  subroutine refinement_strategy_free(this)
    class(refinement_strategy_t), intent(inout) :: this
    nullify(this%error_estimator)
    this%current_mesh_iteration = 0
  end subroutine refinement_strategy_free
  
  function refinement_strategy_get_error_estimator(this)
    class(refinement_strategy_t), intent(inout) :: this
    class(error_estimator_t), pointer :: refinement_strategy_get_error_estimator
    refinement_strategy_get_error_estimator => this%error_estimator
  end function refinement_strategy_get_error_estimator
  
  function refinement_strategy_get_current_mesh_iteration(this)
    class(refinement_strategy_t), intent(inout) :: this
    integer(ip) :: refinement_strategy_get_current_mesh_iteration
    refinement_strategy_get_current_mesh_iteration = this%current_mesh_iteration
  end function refinement_strategy_get_current_mesh_iteration
  
  subroutine urs_set_parameters(this,parameter_list)
    class(uniform_refinement_strategy_t), intent(inout) :: this
    type(parameterlist_t)               , intent(in)    :: parameter_list
    integer(ip) :: FPLerror
    assert(parameter_list%isPresent(num_uniform_refinements_key))
    assert(parameter_list%isAssignable(num_uniform_refinements_key,this%num_uniform_refinements))
    FPLerror = parameter_list%get(key = num_uniform_refinements_key, value = this%num_uniform_refinements)
    assert(FPLerror==0)
  end subroutine urs_set_parameters
  
  subroutine urs_update_refinement_flags(this,refinement_and_coarsening_flags)
    class(uniform_refinement_strategy_t), intent(inout) :: this
    type(std_vector_integer_ip_t)       , intent(inout) :: refinement_and_coarsening_flags
    integer(ip), pointer :: refinement_and_coarsening_flags_entries(:)
    integer(ip)          :: i, num_of_entries
    refinement_and_coarsening_flags_entries => refinement_and_coarsening_flags%get_pointer()
    num_of_entries = size(refinement_and_coarsening_flags_entries)
    do i = 1, num_of_entries
      refinement_and_coarsening_flags_entries(i) = refinement
    end do
    this%current_mesh_iteration = this%current_mesh_iteration + 1
  end subroutine urs_update_refinement_flags
  
  function urs_has_finished_refinement(this)
    class(uniform_refinement_strategy_t), intent(inout) :: this
    logical :: urs_has_finished_refinement
    urs_has_finished_refinement = ( this%current_mesh_iteration > this%num_uniform_refinements )
  end function urs_has_finished_refinement
  
  subroutine eors_set_parameters(this,parameter_list)
    class(error_objective_refinement_strategy_t), intent(inout) :: this
    type(parameterlist_t)                       , intent(in)    :: parameter_list
    integer(ip) :: FPLerror
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
  end subroutine eors_set_parameters
  
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
        ! TO-DO: It seems that the cell is coarsened, even if
        !        one of its siblings is marked as `do_nothing`
        !        or `refinement`. This kills the algorithm.
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
  
  subroutine ffrs_set_parameters(this,parameter_list)
    implicit none
    class(fixed_fraction_refinement_strategy_t), intent(inout) :: this
    type(parameterlist_t)                       , intent(in)    :: parameter_list
    integer(ip) :: FPLerror
    ! Parse refinement/coarsening fraction
    assert(parameter_list%isPresent(refinement_fraction_key))
    assert(parameter_list%isAssignable(refinement_fraction_key,this%refinement_fraction))
    FPLerror = parameter_list%get(key = refinement_fraction_key, value = this%refinement_fraction)
    assert(FPLerror==0)
    assert ( this%refinement_fraction >= 0.0_rp .and. this%refinement_fraction <= 1.0_rp )
    
    assert(parameter_list%isPresent(coarsening_fraction_key))
    assert(parameter_list%isAssignable(coarsening_fraction_key,this%coarsening_fraction))
    FPLerror = parameter_list%get(key = coarsening_fraction_key, value = this%coarsening_fraction)
    assert(FPLerror==0)
    assert ( this%coarsening_fraction >= 0.0_rp .and. this%coarsening_fraction <= 1.0_rp )
    
    ! Parse refinement_threshold
    if ( parameter_list%isPresent(threshold_tolerance_key) ) then
      FPLerror = parameter_list%get(key = threshold_tolerance_key, value = this%threshold_tolerance)
      assert(FPLerror==0)
    else
      this%threshold_tolerance = 3.0_rp*1.0e-08_rp ! \approx 1/(2**25), with 25 being the 
                                                              ! number of recursive bisection steps
    end if
    
    if ( parameter_list%isPresent(max_num_mesh_iterations_key) ) then
      FPLerror = parameter_list%get(key = max_num_mesh_iterations_key, value = this%max_num_mesh_iterations)
      assert(FPLerror==0)
    else
      this%max_num_mesh_iterations = 10
    end if
    
  end subroutine ffrs_set_parameters
  
  subroutine ffrs_update_refinement_flags(this,refinement_and_coarsening_flags)
    implicit none
    class(fixed_fraction_refinement_strategy_t), intent(inout) :: this
    type(std_vector_integer_ip_t)               , intent(inout) :: refinement_and_coarsening_flags
    integer(ip), pointer :: refinement_and_coarsening_flags_entries(:)
    integer(ip)          :: i, num_local_cells
    integer(igp)         :: num_global_cells
    integer(igp)         :: target_num_cells_to_be_refined_coarsened(2)
    integer(igp)         :: current_num_cells_to_be_refined_coarsened(2)
    real(rp)   , pointer :: sq_local_estimate_entries(:)
    class(serial_fe_space_t), pointer :: fe_space
    class(triangulation_t), pointer :: triangulation
    class(environment_t), pointer :: environment
    real(rp) :: ref_sq_min_estimate, ref_min_estimate
    real(rp) :: ref_sq_max_estimate, ref_max_estimate
    real(rp) :: coarsening_sq_max_estimate, coarsening_max_estimate
    real(rp) :: coarsening_sq_min_estimate, coarsening_min_estimate
    real(rp) :: coarsening_split_estimate, ref_split_estimate
    integer(ip) :: num_iterations
    real(rp) :: aux(2)
    real(rp) :: initial_interval_size
    
    assert ( associated(this%error_estimator) )
    
    fe_space      => this%error_estimator%get_fe_space()
    triangulation => fe_space%get_triangulation()
    environment   => triangulation%get_environment()
    
    if ( environment%am_i_l1_task() ) then
      refinement_and_coarsening_flags_entries => refinement_and_coarsening_flags%get_pointer()
      sq_local_estimate_entries => this%error_estimator%get_sq_local_estimate_entries()
    
      num_local_cells  = triangulation%get_num_local_cells()
      num_global_cells = num_local_cells
      call environment%l1_sum(num_global_cells)
    
      target_num_cells_to_be_refined_coarsened(1) = int(real(num_global_cells,rp)*this%refinement_fraction)
      target_num_cells_to_be_refined_coarsened(2) = int(real(num_global_cells,rp)*this%coarsening_fraction)
    
      ref_sq_min_estimate = minval(sq_local_estimate_entries(1:num_local_cells))
      ref_sq_max_estimate = maxval(sq_local_estimate_entries(1:num_local_cells))
    
      aux(1) = -sqrt(ref_sq_min_estimate)
      aux(2) = sqrt(ref_sq_max_estimate)
      call environment%l1_max(aux)
      ref_min_estimate = -aux(1)
      ref_max_estimate = aux(2)
      coarsening_min_estimate = ref_min_estimate
      coarsening_max_estimate = ref_max_estimate

      initial_interval_size = abs(ref_max_estimate-ref_min_estimate)
      num_iterations=0
      do while ( abs(ref_max_estimate-ref_min_estimate)/initial_interval_size > this%threshold_tolerance) 
        ! Compute interval split point using the fact that the log of error estimators
        ! is much better uniformly scattered than the error estimators themselves. This is required
        ! in order to have faster convergence whenever the error estimators are scattered across very
        ! different orders of magnitude
        ! avg_estimate = exp(1/2*(log(min_estimate)+log(max_estimate))) = sqrt(min_estimate*max_estimate)
        if (ref_min_estimate == 0.0_rp) then
          ref_split_estimate = sqrt(1.0e-10_rp*ref_max_estimate)
        else
          ref_split_estimate = sqrt(ref_min_estimate*ref_max_estimate)
        end if
        this%sq_refinement_threshold = ref_split_estimate*ref_split_estimate
        if (coarsening_min_estimate == 0.0_rp) then
          coarsening_split_estimate = sqrt(1.0e-10_rp*coarsening_max_estimate)
        else
          coarsening_split_estimate = sqrt(coarsening_min_estimate*coarsening_max_estimate)
        end if
        this%sq_coarsening_threshold = coarsening_split_estimate*coarsening_split_estimate
        
         
        ! Count how many cells have local error estimate larger or equal to avg_estimate
        ! count = #{ i: e_i >= avg_estimate }
        current_num_cells_to_be_refined_coarsened(1) = 0 
        current_num_cells_to_be_refined_coarsened(2) = 0 
        do i=1, num_local_cells
          if ( sq_local_estimate_entries(i) >= this%sq_refinement_threshold ) then
            current_num_cells_to_be_refined_coarsened(1) = current_num_cells_to_be_refined_coarsened(1) + 1 
          end if
          if ( sq_local_estimate_entries(i) < this%sq_coarsening_threshold ) then
            current_num_cells_to_be_refined_coarsened(2) = current_num_cells_to_be_refined_coarsened(2) + 1 
          end if 
        end do 
        call environment%l1_sum(current_num_cells_to_be_refined_coarsened)
       
        if ( current_num_cells_to_be_refined_coarsened(1) > &
             target_num_cells_to_be_refined_coarsened(1) ) then
           ref_min_estimate = ref_split_estimate
        else
           ref_max_estimate = ref_split_estimate
        end if
        
        if ( current_num_cells_to_be_refined_coarsened(2) > &
             target_num_cells_to_be_refined_coarsened(2) ) then
           coarsening_max_estimate = coarsening_split_estimate
        else
           coarsening_min_estimate = coarsening_split_estimate
        end if
        
        num_iterations = num_iterations + 1 
      end do 
      if ( environment%am_i_l1_root() ) then
        write(*,*) "Converged to ", this%threshold_tolerance, " in ", num_iterations, " iterations"
        write(*,*) "Computed refinement threshold squared = ", this%sq_refinement_threshold
        write(*,*) "% cells to be refined = ", real(current_num_cells_to_be_refined_coarsened(1),rp)/real(num_global_cells,rp)
        write(*,*) "Computed coarsening threshold squared = ", this%sq_coarsening_threshold
        write(*,*) "% cells to be coarsened = ", real(current_num_cells_to_be_refined_coarsened(2),rp)/real(num_global_cells,rp)
      end if
      do i = 1, num_local_cells
        if ( sq_local_estimate_entries(i) > this%sq_refinement_threshold ) then
          refinement_and_coarsening_flags_entries(i) = refinement
        else if ( sq_local_estimate_entries(i) < this%sq_coarsening_threshold ) then
         refinement_and_coarsening_flags_entries(i) = coarsening  
        else
         refinement_and_coarsening_flags_entries(i) = do_nothing
        end if
      end do
    end if
    this%current_mesh_iteration = this%current_mesh_iteration + 1
  end subroutine ffrs_update_refinement_flags
  
  function ffrs_has_finished_refinement(this)
    class(fixed_fraction_refinement_strategy_t), intent(inout) :: this
    logical :: ffrs_has_finished_refinement
    real(rp), pointer :: sq_local_estimate_entries(:)
    real(rp)          :: max_local_estimate
    real(rp)          :: sq_error_upper_bound
    ffrs_has_finished_refinement = (this%current_mesh_iteration > this%max_num_mesh_iterations)
  end function ffrs_has_finished_refinement
  
end module refinement_strategy_names
