module base_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE serial_scalar_array_names

    implicit none
# include "debug.i90"
  
    private

    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_START    = 0
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_INIT     = 1
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_SYMBOLIC = 2 ! Symbolic data already computed
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_NUMERIC  = 3 ! Numerical data already computed

    type, abstract :: base_direct_solver_t
    private
        character(len=:), allocatable          :: name
        integer(ip)                            :: state = BASE_DIRECT_SOLVER_STATE_START
        type(sparse_matrix_t), public, pointer :: matrix
        ! Direct solvers info
        integer(ip)                            :: mem_peak_symb
        integer(ip)                            :: mem_perm_symb
        integer(ip)                            :: nz_factors   
        integer(ip)                            :: mem_peak_num 
        real(rp)                               :: Mflops
    contains
    private
        procedure(base_direct_solver_free_clean),         public, deferred :: free_clean
        procedure(base_direct_solver_free_symbolic),      public, deferred :: free_symbolic
        procedure(base_direct_solver_free_numerical),     public, deferred :: free_numerical
        procedure(base_direct_solver_initialize),                 deferred :: initialize
        procedure(base_direct_solver_set_defaults),       public, deferred :: set_defaults
        procedure(base_direct_solver_symbolic_setup),     public, deferred :: symbolic_setup
        procedure(base_direct_solver_numerical_setup),    public, deferred :: numerical_setup
        procedure(base_direct_solver_solve),              public, deferred :: solve
        procedure, public :: set_name                  => base_direct_solver_set_name
        procedure, public :: set_matrix                => base_direct_solver_set_matrix
        procedure, public :: matrix_is_set             => base_direct_solver_matrix_is_set
        procedure, public :: set_state_start           => base_direct_solver_set_state_start
        procedure, public :: set_state_init            => base_direct_solver_set_state_init
        procedure, public :: set_state_symbolic        => base_direct_solver_set_state_symbolic
        procedure, public :: set_state_numeric         => base_direct_solver_set_state_numeric
        procedure, public :: state_is_start            => base_direct_solver_state_is_start
        procedure, public :: state_is_init             => base_direct_solver_state_is_init
        procedure, public :: state_is_symbolic         => base_direct_solver_state_is_symbolic
        procedure, public :: state_is_numeric          => base_direct_solver_state_is_numeric
        procedure, public :: set_mem_peak_symb         => base_direct_solver_set_mem_peak_symb
        procedure, public :: set_mem_perm_symb         => base_direct_solver_set_mem_perm_symb
        procedure, public :: set_nz_factors            => base_direct_solver_set_nz_factors
        procedure, public :: set_Mflops                => base_direct_solver_set_Mflops
    end type

    interface
        subroutine base_direct_solver_initialize(this)
            import base_direct_solver_t
            class(base_direct_solver_t),   intent(inout) :: this
        end subroutine base_direct_solver_initialize

        subroutine base_direct_solver_set_defaults(this)
            import base_direct_solver_t
            class(base_direct_solver_t),   intent(inout) :: this
        end subroutine base_direct_solver_set_defaults

        subroutine base_direct_solver_symbolic_setup(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_symbolic_setup

        subroutine base_direct_solver_numerical_setup(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_numerical_setup

        subroutine base_direct_solver_solve(op, x, y)
            import base_direct_solver_t
            import serial_scalar_array_t
            class(base_direct_solver_t),  intent(inout) :: op
            class(serial_scalar_array_t), intent(in)    :: x
            class(serial_scalar_array_t), intent(inout) :: y
        end subroutine base_direct_solver_solve

        subroutine  base_direct_solver_free_clean(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_clean

        subroutine  base_direct_solver_free_symbolic(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_symbolic

        subroutine  base_direct_solver_free_numerical(this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(inout) :: this
        end subroutine base_direct_solver_free_numerical

        subroutine  base_direct_solver_log_info (this)
            import base_direct_solver_t
            class(base_direct_solver_t), intent(in) :: this
        end subroutine base_direct_solver_log_info
    end interface

public :: base_direct_solver_t

contains

    subroutine base_direct_solver_set_name(this, name)
        class(base_direct_solver_t),   intent(inout) :: this
        character(len=*),              intent(in)    :: name
        assert(.not. this%state_is_symbolic() .or. .not. this%state_is_numeric())
        this%name = name
    end subroutine base_direct_solver_set_name

    subroutine base_direct_solver_set_matrix(this, matrix)
        class(base_direct_solver_t),   intent(inout) :: this
        type(sparse_matrix_t), target, intent(in)    :: matrix
        assert(.not. this%state_is_symbolic() .or. .not. this%state_is_numeric())
        this%matrix => matrix
    end subroutine base_direct_solver_set_matrix

    function base_direct_solver_matrix_is_set(this) result(matrix_is_set)
        class(base_direct_solver_t), intent(inout) :: this
        logical                                    :: matrix_is_set
        matrix_is_set = associated(this%matrix)
    end function base_direct_solver_matrix_is_set

    subroutine base_direct_solver_set_state_start(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_START
    end subroutine base_direct_solver_set_state_start

    subroutine base_direct_solver_set_state_init(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_INIT
    end subroutine base_direct_solver_set_state_init

    subroutine base_direct_solver_set_state_symbolic(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_SYMBOLIC
    end subroutine base_direct_solver_set_state_symbolic

    subroutine base_direct_solver_set_state_numeric(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_NUMERIC
    end subroutine base_direct_solver_set_state_numeric

    function base_direct_solver_state_is_start(this) result(is_start)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_start
        is_start = this%state == BASE_DIRECT_SOLVER_STATE_START
    end function base_direct_solver_state_is_start

    function base_direct_solver_state_is_init(this) result(is_init)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_init
        is_init = this%state == BASE_DIRECT_SOLVER_STATE_INIT
    end function base_direct_solver_state_is_init

    function base_direct_solver_state_is_symbolic(this) result(is_symbolic_setup)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_symbolic_setup
        is_symbolic_setup = this%state == BASE_DIRECT_SOLVER_STATE_SYMBOLIC
    end function base_direct_solver_state_is_symbolic

    function base_direct_solver_state_is_numeric(this) result(is_numerical_setup)
        class(base_direct_solver_t), intent(in) :: this
        logical                                 :: is_numerical_setup
        is_numerical_setup= this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC
    end function base_direct_solver_state_is_numeric

    subroutine base_direct_solver_set_mem_peak_symb(this, mem_peak_symb)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_peak_symb
        this%mem_peak_symb = mem_peak_symb
    end subroutine base_direct_solver_set_mem_peak_symb

    subroutine base_direct_solver_set_mem_perm_symb(this, mem_perm_symb)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_perm_symb
        this%mem_perm_symb = mem_perm_symb
    end subroutine base_direct_solver_set_mem_perm_symb

    subroutine base_direct_solver_set_nz_factors(this, nz_factors)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz_factors
        this%nz_factors = nz_factors
    end subroutine base_direct_solver_set_nz_factors

    subroutine base_direct_solver_set_Mflops(this, Mflops)
        class(base_direct_solver_t), intent(inout) :: this
        real(rp),                    intent(in)    :: Mflops
        this%Mflops = Mflops
    end subroutine base_direct_solver_set_Mflops

end module base_direct_solver_names
