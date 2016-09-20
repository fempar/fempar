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

module wsmp_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE serial_scalar_array_names
    USE base_direct_solver_names
    USE FPL

    implicit none
# include "debug.i90"
  
    private

    type, extends(base_direct_solver_t) :: wsmp_direct_solver_t
    private

    contains
    private
        procedure, public :: is_linear               => wsmp_direct_solver_is_linear
        procedure, public :: free_clean_body         => wsmp_direct_solver_free_clean_body
        procedure, public :: free_symbolic_body      => wsmp_direct_solver_free_symbolic_body
        procedure, public :: free_numerical_body     => wsmp_direct_solver_free_numerical_body
        procedure         :: initialize              => wsmp_direct_solver_initialize
        procedure         :: set_defaults            => wsmp_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl  => wsmp_direct_solver_set_parameters_from_pl
        procedure, public :: symbolic_setup_body     => wsmp_direct_solver_symbolic_setup_body
        procedure, public :: numerical_setup_body    => wsmp_direct_solver_numerical_setup_body
        procedure, public :: solve_single_rhs_body   => wsmp_direct_solver_solve_single_rhs_body
        procedure, public :: solve_several_rhs_body  => wsmp_direct_solver_solve_several_rhs_body
    end type

contains

    subroutine wsmp_direct_solver_initialize(this)
        class(wsmp_direct_solver_t),   intent(inout) :: this
    end subroutine wsmp_direct_solver_initialize

    subroutine wsmp_direct_solver_set_defaults(this)
        class(wsmp_direct_solver_t),   intent(inout) :: this
    end subroutine wsmp_direct_solver_set_defaults

    subroutine wsmp_direct_solver_set_parameters_from_pl(this, parameter_list)
        class(wsmp_direct_solver_t),  intent(inout) :: this
        type(ParameterList_t),        intent(in)    :: parameter_list
    end subroutine wsmp_direct_solver_set_parameters_from_pl

    subroutine wsmp_direct_solver_symbolic_setup_body(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_symbolic_setup_body

    subroutine wsmp_direct_solver_numerical_setup_body(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_numerical_setup_body

    function wsmp_direct_solver_is_linear(this) result(is_linear)
        class(wsmp_direct_solver_t), intent(inout) :: this
        logical                                            :: is_linear
    end function wsmp_direct_solver_is_linear

    subroutine wsmp_direct_solver_solve_single_rhs_body(op, x, y)
        class(wsmp_direct_solver_t),  intent(inout) :: op
        type(serial_scalar_array_t),  intent(in)    :: x
        type(serial_scalar_array_t),  intent(inout) :: y
    end subroutine wsmp_direct_solver_solve_single_rhs_body

    subroutine wsmp_direct_solver_solve_several_rhs_body(op, x, y)
        class(wsmp_direct_solver_t), intent(inout) :: op
        real(rp),                    intent(inout) :: x(:, :)
        real(rp),                    intent(inout) :: y(:, :)
    end subroutine wsmp_direct_solver_solve_several_rhs_body

    subroutine wsmp_direct_solver_free_clean_body(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_clean_body

    subroutine wsmp_direct_solver_free_symbolic_body(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_symbolic_body

    subroutine wsmp_direct_solver_free_numerical_body(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_numerical_body

end module wsmp_direct_solver_names
