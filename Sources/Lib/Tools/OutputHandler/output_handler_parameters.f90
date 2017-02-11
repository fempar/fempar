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
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Public parameters of the [[output_handler_names(module)]] module
!
! Contains the following public entities:
! [[output_handler_parameters_names(module)]]
!
!---------------------------------------------------------------------
module output_handler_parameters_names
!---------------------------------------------------------------------
!* Author: Víctor Sande Veiga
! Date: 2016-11-29
! Version: 0.0.1
! Category: IO
!
!--------------------------------------------------------------------- 
!### Public parameters of the [[output_handler_names(module)]] module
!
! Contains the following public parameters:
! [[output_handler_parameters_names(module):oh_staticgrid(variable)]], 
! [[output_handler_parameters_names(module):no_diff_operator(variable)]], 
! [[output_handler_parameters_names(module):grad_diff_operator(variable)]], 
! [[output_handler_parameters_names(module):div_diff_operator(variable)]], 
! [[output_handler_parameters_names(module):curl_diff_operator(variable)]], 
!---------------------------------------------------------------------

implicit none
private

    character(*), parameter, public :: oh_staticgrid = 'OH_StaticGrid'

    ! Diff operators
    character(len=*), parameter, public :: no_diff_operator = 'no_diff_operator'
    character(len=*), parameter, public :: grad_diff_operator = 'grad_diff_operator'
    character(len=*), parameter, public :: div_diff_operator  = 'div_diff_operator'
    character(len=*), parameter, public :: curl_diff_operator = 'curl_diff_operator'

end module output_handler_parameters_names
