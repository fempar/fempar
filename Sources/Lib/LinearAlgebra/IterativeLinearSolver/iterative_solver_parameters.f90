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
  character(len=*), parameter :: res_nrmgiven_rhs_nrmgiven  = "RES_NRMGIVEN_RHS_NRMGIVEN"  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  character(len=*), parameter :: res_nrmgiven_res_nrmgiven  = "RES_NRMGIVEN_RES_NRMGIVEN"  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  character(len=*), parameter :: delta_rhs                  = "DELTA_RHS"                  ! || dx(i) ||  <= rtol*||  b  || + atol
  character(len=*), parameter :: delta_delta                = "DELTA_DELTA"                ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  character(len=*), parameter :: res_res                    = "RES_RES"                    ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  character(len=*), parameter :: res_rhs                    = "RES_RHS"                    ! ||  r(i) ||  <= rtol*||  b  || + atol
  character(len=*), parameter :: delta_rhs_and_res_res      = "DELTA_RHS_AND_RES_RES"      ! delta_rhs    AND res_res
  character(len=*), parameter :: delta_rhs_and_res_rhs      = "DELTA_RHS_AND_RES_RHS"      ! delta_rhs    AND res_rhs
  character(len=*), parameter :: delta_delta_and_res_res    = "DELTA_DELTA_AND_RES_RES"    ! delta_delta  AND res_res
  character(len=*), parameter :: delta_delta_and_res_rhs    = "DELTA_DELTA_AND_RES_RHS"    ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration

  character(len=*), parameter, public :: ils_stopping_criterium_cla_choices  = res_nrmgiven_rhs_nrmgiven // "," // &
                                                                               res_nrmgiven_res_nrmgiven // "," // &
                                                                               delta_rhs  // "," // &
                                                                               delta_delta // "," // &
                                                                               res_res // "," // &
                                                                               res_rhs // "," // &
                                                                               delta_rhs_and_res_res // "," // &
                                                                               delta_rhs_and_res_rhs // "," // &
                                                                               delta_delta_and_res_res // "," // &
                                                                               delta_delta_and_res_rhs
                                                                               
   character(len=*), parameter, public :: ils_stopping_criterium_cla_help    = "stopping criterium type" // BRK_LINE // & 
                                                                               BULLET_FLAP_HELP_MESSAGE // res_nrmgiven_rhs_nrmgiven // ": ||r(i)||_g <= rtol*||b||_g + atol; r(i)=b-A*x(i)" // BRK_LINE // & 
                                                                               BULLET_FLAP_HELP_MESSAGE // res_nrmgiven_res_nrmgiven // ": ||r(i)||_g <= rtol*||r(0)||_g + atol; r(i)=b-A*x(i)" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE //  delta_rhs  // ": ||dx(i)|| <= rtol*||b|| + atol; dx(i)=x(i)-x(i-1)" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // delta_delta // ": ||dx(i)|| <= rtol*||dx(1)|| + atol; dx(i)=x(i)-x(i-1)" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // res_res // ": ||r(i)|| <= rtol*||r(0)|| + atol" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // res_rhs // ": ||r(i)||  <= rtol*||b|| + atol" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // delta_rhs_and_res_res // ": delta_rhs AND res_res" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // delta_rhs_and_res_rhs // ": delta_rhs AND res_rhs" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // delta_delta_and_res_res // ": delta_delta AND res_res" // BRK_LINE // &
                                                                               BULLET_FLAP_HELP_MESSAGE // delta_delta_and_res_rhs // ": delta_delta AND res_rhs"
                                                                               
  !-------------------------------------------------------------------
  ! String parameters with the names of the parameters for iterative linear solvers
  !-------------------------------------------------------------------
  character(len=*), parameter :: ils_type_key                      = 'ILS_TYPE'
  character(len=*), parameter :: ils_rtol_key                      = 'ILS_RTOL'
  character(len=*), parameter :: ils_atol_key                      = 'ILS_ATOL'
  character(len=*), parameter :: ils_stopping_criterium_key        = 'ILS_STOPPING_CRITERIUM'
  character(len=*), parameter :: ils_output_frequency_key          = 'ILS_OUTPUT_FREQUENCY'
  character(len=*), parameter :: ils_max_num_iterations_key        = 'ILS_MAX_NUM_ITERATIONS'
  character(len=*), parameter :: ils_track_convergence_history_key = 'ILS_TRACK_CONVERGENCE_HISTORY'
  character(len=*), parameter :: ils_luout_key                     = 'ILS_LUOUT'

  character(len=*), parameter :: ils_type_cla_name                      = '--'//ils_type_key
  character(len=*), parameter :: ils_rtol_cla_name                      = '--'//ils_rtol_key
  character(len=*), parameter :: ils_atol_cla_name                      = '--'//ils_atol_key
  character(len=*), parameter :: ils_stopping_criterium_cla_name        = '--'//ils_stopping_criterium_key
  character(len=*), parameter :: ils_output_frequency_cla_name          = '--'//ils_output_frequency_key
  character(len=*), parameter :: ils_max_num_iterations_cla_name        = '--'//ils_max_num_iterations_key
  character(len=*), parameter :: ils_track_convergence_history_cla_name = '--'//ils_track_convergence_history_key
  character(len=*), parameter :: ils_luout_cla_name                     = '--'//ils_luout_key


  !-----------------------------------------------------------------
  ! Iterative linear solver names
  !-----------------------------------------------------------------
  character(len=*), parameter :: cg_name         = 'CG'         ! CG iterative linear solver
  character(len=*), parameter :: fgmres_name     = 'FGMRES'     ! FGMRES iterative linear solver
  character(len=*), parameter :: icg_name        = 'ICG'        ! ICG iterative linear solver
  character(len=*), parameter :: lfom_name       = 'LFOM'       ! LFOM iterative linear solver
  character(len=*), parameter :: lgmres_name     = 'LGMRES'     ! LGMRES iterative linear solver
  character(len=*), parameter :: minres_name     = 'MINRES'     ! MINRES iterative linear solver
  character(len=*), parameter :: rgmres_name     = 'RGMRES'     ! RGMRES iterative linear solver
  character(len=*), parameter :: richardson_name = 'RICHARDSON' ! RICHARDSON iterative linear solver
  
  character(len=*), parameter, public :: ils_type_cla_choices  = cg_name         // "," // &
                                                                 icg_name        // "," // &
                                                                 minres_name     // "," // &
                                                                 lfom_name       // "," // &
                                                                 lgmres_name     // "," // &
                                                                 rgmres_name     // "," // &
                                                                 fgmres_name     // "," // &
                                                                 richardson_name
 
  character(len=*), parameter, public :: ils_type_cla_help    = "Iterative linear solver type" // BRK_LINE // & 
                                                               BULLET_FLAP_HELP_MESSAGE // cg_name // ": Preconditioned Conjugate Gradients Solver" // BRK_LINE // & 
                                                               BULLET_FLAP_HELP_MESSAGE // icg_name // ": Preconditioned Inexact Conjugate Gradients Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // minres_name // ": Preconditioned Minimum Residual Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // lfom_name // ": Left-Preconditioned Full Orthogonalization Method Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // lgmres_name // ": Left-Preconditioned Generalized Minimum Residual Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // rgmres_name // ": Right-Preconditioned Generalized Minimum Residual Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // fgmres_name // ": Preconditioned Flexible Generalized Minimum Residual Solver" // BRK_LINE // &
                                                               BULLET_FLAP_HELP_MESSAGE // richardson_name // ": (Relaxed) Preconditioned Richardson Fixed-Point Solver"
                                                               
  !-----------------------------------------------------------------
  ! Some common parameters to FGMRES, LFOM, LGMRES and RGMRES iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_max_dim_krylov_basis_key = 'ILS_MAX_DIM_KRYLOV_BASIS'
  character(len=*), parameter :: ils_orthonorm_strategy_key   = 'ILS_ORTHONORM_STRATEGY'

  character(len=*), parameter :: ils_max_dim_krylov_basis_cla_name = '--'//ils_max_dim_krylov_basis_key
  character(len=*), parameter :: ils_orthonorm_strategy_cla_name   = '--'//ils_orthonorm_strategy_key

  character(len=*), parameter :: orthonorm_strat_icgsro  = 'ICGSRO'  ! icgs: Iterative Classical Gram-Schmidt (appropriate for distributed GMRES)
  character(len=*), parameter :: orthonorm_strat_mgsro   = 'MGSRO'   ! mgs : Modified Gram-Schmidt (appropriate for serial GMRES)
  
  character(len=*), parameter, public :: ils_orthonorm_strategy_cla_choices  = orthonorm_strat_icgsro // "," // &
                                                                               orthonorm_strat_mgsro
 
  character(len=*), parameter, public :: ils_orthonorm_strategy_cla_help    = "Orthonormalization strategy (for XGMRES/LFOM)" // BRK_LINE // & 
                                                               BULLET_FLAP_HELP_MESSAGE // orthonorm_strat_icgsro // ": Iterative Classical Gram-Schmidt (appropriate for MPI-parallel variants)" // BRK_LINE // & 
                                                               BULLET_FLAP_HELP_MESSAGE // orthonorm_strat_mgsro // ": Modified Gram-Schmidt (appropriate for serial variants)"
  
  !-----------------------------------------------------------------
  ! Parameters used in RICHARDSON iterative linear solvers
  !-----------------------------------------------------------------
  character(len=*), parameter :: ils_relaxation_key = 'ILS_RELAXATION'
  character(len=*), parameter :: ils_relaxation_cla_name = '--'//ils_relaxation_key
  
  
  !------------------------------------------------------------------------------------
  ! Default values for implementors of class(base_iterative_linear_solver_t) parameters
  !------------------------------------------------------------------------------------
  integer (ip), parameter :: default_luout                            = 6
  real    (rp), parameter :: default_rtol                             = 1.0e-06_rp
  real    (rp), parameter :: default_atol                             = 0.0_rp
  integer (ip), parameter :: default_output_frequency                 = 1 
  integer (ip), parameter :: default_max_num_iterations               = 1000
  logical,      parameter :: default_track_convergence_history        = .false.
  character(len=*), parameter :: default_fgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  character(len=*), parameter :: default_lfom_stopping_criteria       = res_res
  character(len=*), parameter :: default_lgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  character(len=*), parameter :: default_rgmres_stopping_criteria     = res_nrmgiven_res_nrmgiven
  character(len=*), parameter :: default_minres_stopping_criteria     = res_res
  character(len=*), parameter :: default_cg_stopping_criteria         = res_res
  character(len=*), parameter :: default_icg_stopping_criteria        = res_res
  character(len=*), parameter :: default_richardson_stopping_criteria = res_res
  real (rp),    parameter :: default_richardson_relaxation            = 1.0_rp
  integer (ip), parameter :: default_dkrymax                          = 1000
  character(*), parameter :: default_orthonorm_strat                  = orthonorm_strat_icgsro
  
end module iterative_linear_solver_parameters_names
