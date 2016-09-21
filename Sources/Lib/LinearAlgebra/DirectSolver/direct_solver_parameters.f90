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
implicit none

    !-----------------------------------------------------------------
    ! Parameters used in DIRECT SOLVERS
    !-----------------------------------------------------------------
    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: direct_solver_type = 'direct_solver_type'
    character(len=*), parameter :: pardiso_mkl        = 'pardiso_mkl'                      ! Name of the PARDISO MKL direct solver type
    character(len=*), parameter :: umfpack            = 'umfpack'                          ! Name of the UMFPACK direct solver type

    !-----------------------------------------------------------------
    ! Parameters used in PARDISO_MKL direct solver
    !-----------------------------------------------------------------

    ! PARDISO MKL matrix types
    integer,          parameter :: pardiso_mkl_spd =  2                                    ! Real Symmetric positive definite 
    integer,          parameter :: pardiso_mkl_sin = -2                                    ! Real Symmetric indefinite
    integer,          parameter :: pardiso_mkl_uss = 1                                     ! Real Unsymmetric, structurally symmetric
    integer,          parameter :: pardiso_mkl_uns = 11                                    ! Real Unsymmetric, structurally unsymmetric

    ! PARDISO MKL default values
    integer,          parameter :: pardiso_mkl_default_message_level = 0                   ! Default pardiso_mkl message level
    integer,          parameter :: pardiso_mkl_default_iparm         = 0                   ! Default pardiso_mkl_iparam(64) = 0
    integer,          parameter :: pardiso_mkl_default_matrix_type   = pardiso_mkl_uns     ! Default pardiso_mkl matrix type

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: pardiso_mkl_iparm         = 'pardiso_mkl_iparm'         ! PARDISO MKL control parameters array
    character(len=*), parameter :: pardiso_mkl_matrix_type   = 'pardiso_mkl_matrix_type'   ! PARDISO MKL matrix type
    character(len=*), parameter :: pardiso_mkl_message_level = 'pardiso_mkl_message_level' ! PARDISO MKL verbosity level


    !-----------------------------------------------------------------
    ! Parameters used in UMFPACK direct solver
    !-----------------------------------------------------------------

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: umfpack_control_params   = 'umfpack_control_params'      ! UMFPACK real array of 20 UMFPACK parameters

end module direct_solver_parameters_names
