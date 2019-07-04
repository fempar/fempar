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

module direct_solver_parameters_names
  use types_names
implicit none

    !-----------------------------------------------------------------
    ! Parameters used in DIRECT SOLVERS
    !-----------------------------------------------------------------
    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: dls_type_key       = 'DLS_TYPE'
    character(len=*), parameter :: dls_type_cla_name  = '--'//dls_type_key
    character(len=*), parameter :: pardiso_mkl        = 'PARDISO_MKL'                      ! Name of the PARDISO MKL direct solver type
    character(len=*), parameter :: umfpack            = 'UMFPACK'                          ! Name of the UMFPACK direct solver type

    !-----------------------------------------------------------------
    ! Parameters used in PARDISO_MKL direct solver
    !-----------------------------------------------------------------

    ! PARDISO MKL matrix types
    integer,          parameter :: pardiso_mkl_type_guess_from_matrix_properties = 0       ! Fempar defined matrix type to define type from matrix properties
    integer,          parameter :: pardiso_mkl_spd =  2                                    ! Real Symmetric positive definite 
    integer,          parameter :: pardiso_mkl_sin = -2                                    ! Real Symmetric indefinite
    integer,          parameter :: pardiso_mkl_uss = 1                                     ! Real Unsymmetric, structurally symmetric
    integer,          parameter :: pardiso_mkl_uns = 11                                    ! Real Unsymmetric, structurally unsymmetric

    ! PARDISO MKL default values
    integer,          parameter :: pardiso_mkl_default_message_level = 0                   ! Default pardiso_mkl message level
    integer,          parameter :: pardiso_mkl_default_iparm         = 0                   ! Default pardiso_mkl_iparam(64) = 0
    integer,          parameter :: pardiso_mkl_default_matrix_type   = pardiso_mkl_type_guess_from_matrix_properties ! Default pardiso_mkl matrix type

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: pardiso_mkl_iparm         = 'PARDISO_PARAMS'           ! PARDISO MKL control parameters array
    character(len=*), parameter :: pardiso_mkl_matrix_type   = 'PARDISO_MAT_TYPE'         ! PARDISO MKL matrix type
    character(len=*), parameter :: pardiso_mkl_message_level = 'PARDISO_MSG_LEVEL'        ! PARDISO MKL verbosity level

    character(len=*), parameter :: pardiso_mkl_iparm_cla_name         = '--'//pardiso_mkl_iparm
    character(len=*), parameter :: pardiso_mkl_matrix_type_cla_name   = '--'//pardiso_mkl_matrix_type
    character(len=*), parameter :: pardiso_mkl_message_level_cla_name = '--'//pardiso_mkl_message_level


    ! Parameter choices
    character(len=*), parameter :: dls_type_cla_choices  = pardiso_mkl//','//umfpack
    character(len=*), parameter :: pardiso_mkl_message_level_cla_choices  = '0,1'
    character(len=*), parameter :: pardiso_mkl_matrix_type_cla_choices = '0,1,2,-2,11'

    ! Parameter help

    character(len=*), parameter :: dls_type_cla_help = "Direct solver type" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // pardiso_mkl // ": Intel MKL PARDISO - Parallel Direct Sparse Solver Interface" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // umfpack // ": UMFPACK - Unsymmetric MultiFrontal PACKage"

    character(len=*), parameter :: pardiso_mkl_iparm_cla_help         = "PARDISO iparm array parameter (see PARDISO users' manual for additional details)"
    character(len=*), parameter :: pardiso_mkl_message_level_cla_help = "PARDISO Message level (see PARDISO users' manual for additional details)" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // "0: No output is generated" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // "1: Prints statistical information"

    character(len=*), parameter :: pardiso_mkl_matrix_type_cla_help = "PARDISO mtype parameter (see PARDISO users' manual for additional details)" // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // " 0: Guess matrix type from matrix properties" // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 1: Real and structurally symmetric"          // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // " 2: Real and symmetric positive definite"     // BRK_LINE // &
                   BULLET_FLAP_HELP_MESSAGE // "-2: Real and symmetric indefinite"            // BRK_LINE // & 
                   BULLET_FLAP_HELP_MESSAGE // "11: Real and unsymmetric matrix"


    !-----------------------------------------------------------------
    ! Parameters used in UMFPACK direct solver
    !-----------------------------------------------------------------

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: umfpack_control_params            = 'UMFPACK_CONTROL_PARAMS' ! UMFPACK real array of 20 UMFPACK parameters
    character(len=*), parameter :: umfpack_control_params_cla_name   = '--'//umfpack_control_params
    character(len=*), parameter :: umfpack_control_params_cla_help   = "UMFPACK control array parameter (see UMFPACK users' manual for additional details)"

end module direct_solver_parameters_names
