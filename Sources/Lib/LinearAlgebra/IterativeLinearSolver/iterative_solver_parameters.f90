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
module iterative_linear_solver_parameters_names

  use types_names

  implicit none

  !-------------------------------------------------------------------
  ! List of convergence criteria available for iterative solvers 
  !-------------------------------------------------------------------
  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration

  !-------------------------------------------------------------------
  ! String parameters with the names of the parameters for iterative linear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: ils_type                      = 'iterative_linear_solver_type'
  character(len=*), parameter :: ils_rtol                      = 'iterative_linear_solver_rtol'
  character(len=*), parameter :: ils_atol                      = 'iterative_linear_solver_atol'
  character(len=*), parameter :: ils_stopping_criteria         = 'iterative_linear_solver_stopping_criteria'
  character(len=*), parameter :: ils_output_frequency          = 'iterative_linear_solver_output_frequency'
  character(len=*), parameter :: ils_max_num_iterations        = 'iterative_linear_solver_max_num_iterations'
  character(len=*), parameter :: ils_track_convergence_history = 'iterative_linear_solver_track_convergence_history'
  character(len=*), parameter :: ils_luout                      = 'iterative_linear_solver_luout'

  !-----------------------------------------------------------------
  ! Iterative linear solver names
  !-----------------------------------------------------------------
  character(len=*), parameter :: cg_name         = 'CG'         ! CG iterative linear solver
  character(len=*), parameter :: fgmres_name     = 'FGMRES'     ! FGMREs iterative linear solver
  character(len=*), parameter :: icg_name        = 'ICG'        ! ICG iterative linear solver
  character(len=*), parameter :: lfom_name       = 'LFOM'       ! LFOM iterative linear solver
  character(len=*), parameter :: lgmres_name     = 'LGMRES'     ! LGMRES iterative linear solver
  character(len=*), parameter :: minres_name     = 'MINRES'     ! MINRES iterative linear solver
  character(len=*), parameter :: rgmres_name     = 'RGMRES'     ! RGMRES iterative linear solver
  character(len=*), parameter :: richardson_name = 'RICHARDSON' ! RICHARDSON iterative linear solver


  !-----------------------------------------------------------------
  ! Some common parameters to FGMRES, LFOM, LGMRES and RGMRES iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_dkrymax             = 'iterative_linear_solver_dkrymax'
  character(len=*), parameter :: ils_orthonorm_strat     = 'iterative_linear_solver_orthonorm_strat'
  character(len=*), parameter :: orthonorm_strat_icgsro  = 'ICGSRO' 
  character(len=*), parameter :: orthonorm_strat_mgsro   = 'MGSRO'
  
  integer (ip), parameter :: mgsro  = 1 ! mgs : Modified Gram-Schmidt (appropriate for serial GMRES)
  integer (ip), parameter :: icgsro = 2 ! icgs: Iterative Classical Gram-Schmidt (appropriate for distributed GMRES)

  !-----------------------------------------------------------------
  ! Parameters used in RICHARDSON iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_relaxation = 'iterative_linear_solver_relaxation'
  
  
  !------------------------------------------------------------------------------------
  ! Default values for implementors of class(base_iterative_linear_solver_t) parameters
  !------------------------------------------------------------------------------------
  integer (ip), parameter :: default_luout                        = 6
  real    (rp), parameter :: default_rtol                         = 1.0e-06_rp
  real    (rp), parameter :: default_atol                         = 0.0_rp
  integer (ip), parameter :: default_output_frequency             = 1 
  integer (ip), parameter :: default_max_num_iterations           = 1000
  logical,      parameter :: default_track_convergence_history    = .false.
  integer (ip), parameter :: default_fgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  integer (ip), parameter :: default_lfom_stopping_criteria       = res_res
  integer (ip), parameter :: default_lgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  integer (ip), parameter :: default_rgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  integer (ip), parameter :: default_minres_stopping_criteria     = res_res
  integer (ip), parameter :: default_cg_stopping_criteria         = res_res
  integer (ip), parameter :: default_icg_stopping_criteria        = res_res
  integer (ip), parameter :: default_richardson_stopping_criteria = res_res
  real (rp),    parameter :: default_richardson_relaxation        = 1.0_rp
  integer (ip), parameter :: default_dkrymax                      = 1000
  integer (ip), parameter :: default_orthonorm_strat              = icgsro
  
end module iterative_linear_solver_parameters_names
