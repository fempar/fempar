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
module refinement_strategy_parameters_names

  use types_names

implicit none

  ! Parameter keys

  character(len=*), parameter :: urs_num_uniform_refinements_key   = 'UNIFORM_REFINEMENT_STRATEGY_NUM_UNIFORM_REFINEMENTS'
  
  character(len=*), parameter :: eors_error_objective_key          = 'ERROR_OBJECTIVE_REFINEMENT_STRATEGY_ERROR_OBJECTIVE'
  character(len=*), parameter :: eors_objective_tolerance_key      = 'ERROR_OBJECTIVE_REFINEMENT_STRATEGY_OBJETIVE_TOLERANCE'
  character(len=*), parameter :: eors_max_num_mesh_iterations_key  = 'ERROR_OBJECTIVE_REFINEMENT_STRATEGY_MAX_NUM_MESH_ITERATIONS'
  
  character(len=*), parameter :: ffrs_refinement_fraction_key      = 'FIXED_FRACTION_REFINEMENT_STRATEGY_REFINEMENT_FRACTION'
  character(len=*), parameter :: ffrs_coarsening_fraction_key      = 'FIXED_FRACTION_REFINEMENT_STRATEGY_COARSENING_FRACTION'
  character(len=*), parameter :: ffrs_max_num_mesh_iterations_key  = 'FIXED_FRACTION_REFINEMENT_STRATEGY_MAX_NUM_MESH_ITERATIONS'
  character(len=*), parameter :: ffrs_print_info_key               = 'FIXED_FRACTION_REFINEMENT_STRATEGY_PRINT_INFO'

  ! Parameter CLA names

  character(len=*), parameter :: urs_num_uniform_refinements_cla_name = '--'//urs_num_uniform_refinements_key
  
  character(len=*), parameter :: eors_error_objective_cla_name         = '--'//eors_error_objective_key
  character(len=*), parameter :: eors_objective_tolerance_cla_name     = '--'//eors_objective_tolerance_key
  character(len=*), parameter :: eors_max_num_mesh_iterations_cla_name = '--'//eors_max_num_mesh_iterations_key
  
  character(len=*), parameter :: ffrs_refinement_fraction_cla_name     = '--'//ffrs_refinement_fraction_key
  character(len=*), parameter :: ffrs_coarsening_fraction_cla_name     = '--'//ffrs_coarsening_fraction_key
  character(len=*), parameter :: ffrs_max_num_mesh_iterations_cla_name = '--'//ffrs_max_num_mesh_iterations_key
  character(len=*), parameter :: ffrs_print_info_cla_name              = '--'//ffrs_print_info_key

  ! Parameter help messages

  character(len=*), parameter :: urs_num_uniform_refinements_help  = 'Number of uniform mesh refinement steps to be' // BRK_LINE // &
                                                                     'performed using the uniform refinement strategy.'
  
  character(len=*), parameter :: eors_error_objective_help         = 'Real quantity expressing the ABSOLUTE global error value' // BRK_LINE //  &
                                                                     'to be achieved by the error objective refinement strategy.'
  character(len=*), parameter :: eors_objective_tolerance_help     = 'Real quantity that increments ERROR_OBJECTIVE_REFINEMENT_STRATEGY_ERROR_OBJECTIVE' // BRK_LINE //  &
                                                                     'to facilitate convergence of the error objective refinement strategy.'
  character(len=*), parameter :: eors_max_num_mesh_iterations_help = 'Maximum number of iterations allowed to the error objective refinement' // BRK_LINE //  &
                                                                     'strategy to find a mesh with ABSOLUTE global error below' // BRK_LINE //  &
                                                                     'ERROR_OBJECTIVE_REFINEMENT_STRATEGY_ERROR_OBJECTIVE.' 
  
  character(len=*), parameter :: ffrs_refinement_fraction_help     = 'Fraction of total number of cells to be set for REFINEMENT' // BRK_LINE //  &
                                                                     'in a single mesh adaptation step of the fixed fraction refinement strategy.'
  character(len=*), parameter :: ffrs_coarsening_fraction_help     = 'Fraction of total number of cells to be set for COARSENING' // BRK_LINE //  &
                                                                     'in a single mesh adaptation step of the fixed fraction refinement strategy.'
  character(len=*), parameter :: ffrs_max_num_mesh_iterations_help = 'Number of mesh adaptation steps to be performed' // BRK_LINE //  &
                                                                     'within the fixed fraction refinement strategy.'
  character(len=*), parameter :: ffrs_print_info_help              = 'Print convergence status and computed thresholds of' // BRK_LINE //  &
                                                                     'a single step of the fixed fraction refinement strategy.'

  ! Parameter default values

  integer(ip), parameter :: urs_num_uniform_refinements_default   = 0
  
  real(rp),    parameter :: eors_error_objective_default          = 0.0_rp
  real(rp),    parameter :: eors_objective_tolerance_default      = 0.0_rp
  integer(ip), parameter :: eors_max_num_mesh_iterations_default  = 100
  
  real(rp),    parameter :: ffrs_refinement_fraction_default      = 0.5_rp
  real(rp),    parameter :: ffrs_coarsening_fraction_default      = 0.5_rp
  integer(ip), parameter :: ffrs_max_num_mesh_iterations_default  = 10
  logical,     parameter :: ffrs_print_info_default               = .false.


end module refinement_strategy_parameters_names
