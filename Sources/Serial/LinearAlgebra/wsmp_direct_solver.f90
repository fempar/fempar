module wsmp_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE serial_scalar_array_names
    USE base_direct_solver_names

    implicit none
# include "debug.i90"
  
    private

    type, extends(base_direct_solver_t) :: wsmp_direct_solver_t
    private

    contains
    private
        procedure, public :: is_linear          => wsmp_direct_solver_is_linear
        procedure, public :: free_clean         => wsmp_direct_solver_free_clean
        procedure, public :: free_symbolic      => wsmp_direct_solver_free_symbolic
        procedure, public :: free_numerical     => wsmp_direct_solver_free_numerical
        procedure         :: initialize         => wsmp_direct_solver_initialize
        procedure, public :: set_defaults       => wsmp_direct_solver_set_defaults
        procedure, public :: symbolic_setup     => wsmp_direct_solver_symbolic_setup
        procedure, public :: numerical_setup    => wsmp_direct_solver_numerical_setup
        procedure, public :: solve              => wsmp_direct_solver_solve
    end type

contains

    subroutine wsmp_direct_solver_initialize(this)
        class(wsmp_direct_solver_t),   intent(inout) :: this
    end subroutine wsmp_direct_solver_initialize

    subroutine wsmp_direct_solver_set_defaults(this)
        class(wsmp_direct_solver_t),   intent(inout) :: this
    end subroutine wsmp_direct_solver_set_defaults

    subroutine wsmp_direct_solver_symbolic_setup(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_symbolic_setup

    subroutine wsmp_direct_solver_numerical_setup(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_numerical_setup

    function wsmp_direct_solver_is_linear(this) result(is_linear)
        class(wsmp_direct_solver_t), intent(inout) :: this
        logical                                            :: is_linear
    end function wsmp_direct_solver_is_linear

    subroutine wsmp_direct_solver_solve(op, x, y)
        class(wsmp_direct_solver_t),  intent(inout) :: op
        class(serial_scalar_array_t), intent(in)    :: x
        class(serial_scalar_array_t), intent(inout) :: y
    end subroutine wsmp_direct_solver_solve

    subroutine wsmp_direct_solver_free_clean(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_clean

    subroutine wsmp_direct_solver_free_symbolic(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_symbolic

    subroutine wsmp_direct_solver_free_numerical(this)
        class(wsmp_direct_solver_t), intent(inout) :: this
    end subroutine wsmp_direct_solver_free_numerical

end module wsmp_direct_solver_names
