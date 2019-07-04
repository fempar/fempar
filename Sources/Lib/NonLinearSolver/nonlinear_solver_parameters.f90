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
module nonlinear_solver_parameters_names
  use types_names
  implicit none

  !-------------------------------------------------------------------
  ! FPL keys of the parameters for nonlinear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: nls_rtol_key                   = 'NLS_RTOL'
  character(len=*), parameter :: nls_atol_key                   = 'NLS_ATOL'
  character(len=*), parameter :: nls_stopping_criterium_key     = 'NLS_STOPPING_CRITERIUM'
  character(len=*), parameter :: nls_max_num_iterations_key     = 'NLS_MAX_NUM_ITERATIONS'
  character(len=*), parameter :: nls_print_iteration_output_key = 'NLS_PRINT_ITERATION_OUTPUT'

  !-------------------------------------------------------------------
  ! Parameter handler CLA names of the parameters for nonlinear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: nls_rtol_cla_name                   = '--'//nls_rtol_key
  character(len=*), parameter :: nls_atol_cla_name                   = '--'//nls_atol_key
  character(len=*), parameter :: nls_stopping_criterium_cla_name     = '--'//nls_stopping_criterium_key
  character(len=*), parameter :: nls_max_num_iterations_cla_name     = '--'//nls_max_num_iterations_key
  character(len=*), parameter :: nls_print_iteration_output_cla_name = '--'//nls_print_iteration_output_key


  character(len=*), parameter :: abs_res_norm                   = 'abs_res_norm'                  ! |r(x_i)| <= abs_tol
  character(len=*), parameter :: rel_inc_norm                   = 'rel_inc_norm'                  ! |dx_i|   <= rel_norm*|x_i|
  character(len=*), parameter :: rel_r0_res_norm                = 'rel_r0_res_norm'               ! |r(x_i)| <= rel_tol*|r(x_0)|+abs_tol
  character(len=*), parameter :: abs_res_norm_and_rel_inc_norm  = 'abs_res_norm_and_rel_inc_norm' ! |r(x_i)| <= abs_tol & |dx_i| <= rel_norm*|x_i|

  !------------------------------------------------------------------
  ! Default values of the parameters for nonlinear solvers 
  !------------------------------------------------------------------
  real    (rp), parameter :: default_nls_rtol                     = 1.0e-09_rp
  real    (rp), parameter :: default_nls_atol                     = 1.0e-14_rp
  character(*), parameter :: default_nls_stopping_criterium       = rel_r0_res_norm
  integer (ip), parameter :: default_nls_max_iter                 = 10 
  logical     , parameter :: default_nls_print_iteration_output   = .true. 

  !-------------------------------------------------------------------
  ! Parameter handler CLA help messages
  !-------------------------------------------------------------------
  character(len=*), parameter :: nls_rtol_cla_help                   = 'Relative tolerance for the nonlinear solvers stopping criteria'
  character(len=*), parameter :: nls_atol_cla_help                   = 'Absolute tolerance for the nonlinear solvers stopping criteria'
  character(len=*), parameter :: nls_stopping_criterium_cla_help     = 'Nonlinear solvers stopping criterium type'
  character(len=*), parameter :: nls_max_num_iterations_cla_help     = 'Nonlinear solvers maximum number of iterations'
  character(len=*), parameter :: nls_print_iteration_output_cla_help = 'Print output per nonlinear solver iteration'


end module nonlinear_solver_parameters_names
