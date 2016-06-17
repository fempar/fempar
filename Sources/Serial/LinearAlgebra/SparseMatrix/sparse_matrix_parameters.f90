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
module sparse_matrix_parameters_names
   use types_names

   implicit none

    ! Sparse matrix format ID's
    character(len=3), parameter :: coo_format = 'COO'
    character(len=3), parameter :: csr_format = 'CSR'

    ! Transition diagram states
    integer(ip), parameter :: SPARSE_MATRIX_STATE_START              = 0
    integer(ip), parameter :: SPARSE_MATRIX_STATE_PROPERTIES_SET     = 1
    integer(ip), parameter :: SPARSE_MATRIX_STATE_CREATED            = 2
    integer(ip), parameter :: SPARSE_MATRIX_STATE_BUILD_SYMBOLIC     = 3
    integer(ip), parameter :: SPARSE_MATRIX_STATE_BUILD_NUMERIC      = 4
    integer(ip), parameter :: SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC = 5
    integer(ip), parameter :: SPARSE_MATRIX_STATE_ASSEMBLED          = 6
    integer(ip), parameter :: SPARSE_MATRIX_STATE_UPDATE             = 7

    ! Matrix sign
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE     = 0
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE = 1
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_INDEFINITE            = 2 ! Both positive and negative eigenvalues
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_UNKNOWN               = 3 ! No info

    ! COO matrix sort state
    integer(ip),      parameter :: COO_SPARSE_MATRIX_SORTED_NONE    = 20
    integer(ip),      parameter :: COO_SPARSE_MATRIX_SORTED_BY_ROWS = 21
    integer(ip),      parameter :: COO_SPARSE_MATRIX_SORTED_BY_COLS = 22
end module sparse_matrix_parameters_names
