module hsl_ma87_direct_solver_names
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

    type, extends(base_direct_solver_t) :: hsl_ma87_direct_solver_t
    private

    contains
    private
        procedure, public :: free_clean              => hsl_ma87_direct_solver_free_clean
        procedure, public :: free_symbolic           => hsl_ma87_direct_solver_free_symbolic
        procedure, public :: free_numerical          => hsl_ma87_direct_solver_free_numerical
        procedure         :: initialize              => hsl_ma87_direct_solver_initialize
        procedure, public :: set_defaults            => hsl_ma87_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl  => hsl_ma87_direct_solver_set_parameters_from_pl
        procedure, public :: symbolic_setup          => hsl_ma87_direct_solver_symbolic_setup
        procedure, public :: numerical_setup         => hsl_ma87_direct_solver_numerical_setup
        procedure, public :: solve                   => hsl_ma87_direct_solver_solve
    end type

contains

    subroutine hsl_ma87_direct_solver_initialize(this)
        class(hsl_ma87_direct_solver_t),      intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_initialize

    subroutine hsl_ma87_direct_solver_set_defaults(this)
        class(hsl_ma87_direct_solver_t),      intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_set_defaults

    subroutine hsl_ma87_direct_solver_set_parameters_from_pl(this, parameter_list)
        class(hsl_ma87_direct_solver_t),  intent(inout) :: this
        type(ParameterList_t),            intent(in)    :: parameter_list
    end subroutine hsl_ma87_direct_solver_set_parameters_from_pl

    subroutine hsl_ma87_direct_solver_symbolic_setup(this)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_symbolic_setup

    subroutine hsl_ma87_direct_solver_numerical_setup(this)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_numerical_setup

    function hsl_ma87_direct_solver_is_linear(this) result(is_linear)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
        logical                                 :: is_linear
    end function hsl_ma87_direct_solver_is_linear

    subroutine hsl_ma87_direct_solver_solve(op, x, y)
        class(hsl_ma87_direct_solver_t), intent(inout) :: op
        class(serial_scalar_array_t),    intent(in)    :: x
        class(serial_scalar_array_t),    intent(inout) :: y
    end subroutine hsl_ma87_direct_solver_solve

    subroutine hsl_ma87_direct_solver_free_clean(this)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_free_clean

    subroutine hsl_ma87_direct_solver_free_symbolic(this)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_free_symbolic

    subroutine hsl_ma87_direct_solver_free_numerical(this)
        class(hsl_ma87_direct_solver_t), intent(inout) :: this
    end subroutine hsl_ma87_direct_solver_free_numerical
end module hsl_ma87_direct_solver_names
