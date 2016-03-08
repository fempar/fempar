module umfpack_direct_solver_names
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

    type, extends(base_direct_solver_t) :: umfpack_direct_solver_t
    private

    contains
    private
        procedure, public :: is_linear               => umfpack_direct_solver_is_linear
        procedure, public :: free_clean              => umfpack_direct_solver_free_clean
        procedure, public :: free_symbolic           => umfpack_direct_solver_free_symbolic
        procedure, public :: free_numerical          => umfpack_direct_solver_free_numerical
        procedure         :: initialize              => umfpack_direct_solver_initialize
        procedure, public :: set_defaults            => umfpack_direct_solver_set_defaults
        procedure, public :: set_from_parameter_list => umfpack_direct_solver_set_from_parameter_list
        procedure, public :: symbolic_setup          => umfpack_direct_solver_symbolic_setup
        procedure, public :: numerical_setup         => umfpack_direct_solver_numerical_setup
        procedure, public :: solve                   => umfpack_direct_solver_solve
        procedure, public :: log_info                => umfpack_direct_solver_log_info
    end type

contains

    subroutine umfpack_direct_solver_initialize(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_initialize

    subroutine umfpack_direct_solver_set_defaults(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_set_defaults

    subroutine umfpack_direct_solver_set_from_parameter_list(this, parameter_list)
        class(umfpack_direct_solver_t),  intent(inout) :: this
        type(ParameterList_t),           intent(in)    :: parameter_list
    end subroutine umfpack_direct_solver_set_from_parameter_list

    subroutine umfpack_direct_solver_symbolic_setup(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_symbolic_setup

    subroutine umfpack_direct_solver_numerical_setup(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_numerical_setup

    function umfpack_direct_solver_is_linear(this) result(is_linear)
        class(umfpack_direct_solver_t), intent(inout) :: this
        logical                                       :: is_linear
    end function umfpack_direct_solver_is_linear

    subroutine umfpack_direct_solver_solve(op, x, y)
        class(umfpack_direct_solver_t), intent(inout) :: op
        class(serial_scalar_array_t),   intent(in)    :: x
        class(serial_scalar_array_t),   intent(inout) :: y
    end subroutine umfpack_direct_solver_solve

    subroutine umfpack_direct_solver_log_info(this)
        class(umfpack_direct_solver_t), intent(in) :: this
    end subroutine umfpack_direct_solver_log_info

    subroutine umfpack_direct_solver_free_clean(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_free_clean

    subroutine umfpack_direct_solver_free_symbolic(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_free_symbolic

    subroutine umfpack_direct_solver_free_numerical(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
    end subroutine umfpack_direct_solver_free_numerical

end module umfpack_direct_solver_names
