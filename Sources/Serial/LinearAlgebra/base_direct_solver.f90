module base_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE serial_scalar_array_names
    USE FPL

    implicit none
# include "debug.i90"
  
    private

    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_START    = 0
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_SYMBOLIC = 1 ! Symbolic data already computed
    integer(ip), public, parameter :: BASE_DIRECT_SOLVER_STATE_NUMERIC  = 2 ! Numerical data already computed

  !-----------------------------------------------------------------
  ! State transition diagram for type(base_direct_solver_t)
  !-----------------------------------------------------------------
  ! Note: all the child classes must implement this state diagram
  !-----------------------------------------------------------------
  ! Input State         | Action                | Output State 
  !-----------------------------------------------------------------
  ! Start               | symbolic_setup        | Symbolic
  ! Start               | numerical_setup       | Numeric           ! perform symbolic_setup()
  ! Start               | solve                 | Numeric           ! perform numerical_setup()
  ! Start               | free_clean            | Start
  ! Start               | free_symbolic         | Start             ! it does nothing
  ! Start               | free_numeric          | Start             ! it does nothing

  ! Symbolic            | symbolic_setup        | Symbolic
  ! Symbolic            | numerical_setup       | Numeric           
  ! Symbolic            | solve                 | Numeric           ! perform numerical_setup()
  ! Symbolic            | free_clean            | Start
  ! Symbolic            | free_symbolic         | Start
  ! Symbolic            | free_numeric          | Start             ! it does nothing

  ! Numeric             | symbolic_setup        | Symbolic
  ! Numeric             | numeric_setup         | Numeric
  ! Numeric             | solve                 | Numeric
  ! Numeric             | free_numeric          | Symbolic
  ! Numeric             | free_symbolic         | Start
  ! Numeric             | free_clean            | Start

  

    type, abstract :: base_direct_solver_t
    private
        character(len=:), allocatable          :: name
        integer(ip)                            :: state  = BASE_DIRECT_SOLVER_STATE_START
        type(sparse_matrix_t), public, pointer :: matrix => NULL()
        ! Direct solvers info
        integer(ip)                            :: mem_peak_symb
        integer(ip)                            :: mem_perm_symb
        integer(ip)                            :: nz_factors   
        integer(ip)                            :: mem_peak_num 
        real(rp)                               :: Mflops
    contains
    private
        procedure(base_direct_solver_free_clean),              public, deferred :: free_clean
        procedure(base_direct_solver_free_symbolic),           public, deferred :: free_symbolic
        procedure(base_direct_solver_free_numerical),          public, deferred :: free_numerical
        procedure(base_direct_solver_set_defaults),            public, deferred :: set_defaults
        procedure(base_direct_solver_set_parameters_from_pl),  public, deferred :: set_parameters_from_pl
        procedure(base_direct_solver_symbolic_setup),          public, deferred :: symbolic_setup
        procedure(base_direct_solver_numerical_setup),         public, deferred :: numerical_setup
        procedure(base_direct_solver_solve),                   public, deferred :: solve
        procedure, non_overridable, public :: reset              => base_direct_solver_reset
        procedure, non_overridable, public :: set_name           => base_direct_solver_set_name
        procedure, non_overridable, public :: set_matrix         => base_direct_solver_set_matrix
        procedure, non_overridable, public :: update_matrix      => base_direct_solver_update_matrix
        procedure, non_overridable, public :: matrix_is_set      => base_direct_solver_matrix_is_set
        procedure, non_overridable, public :: set_state_start    => base_direct_solver_set_state_start
        procedure, non_overridable, public :: set_state_symbolic => base_direct_solver_set_state_symbolic
        procedure, non_overridable, public :: set_state_numeric  => base_direct_solver_set_state_numeric
        procedure, non_overridable, public :: state_is_start     => base_direct_solver_state_is_start
        procedure, non_overridable, public :: state_is_symbolic  => base_direct_solver_state_is_symbolic
        procedure, non_overridable, public :: state_is_numeric   => base_direct_solver_state_is_numeric
        procedure, non_overridable, public :: set_mem_peak_symb  => base_direct_solver_set_mem_peak_symb
        procedure, non_overridable, public :: set_mem_perm_symb  => base_direct_solver_set_mem_perm_symb
        procedure, non_overridable, public :: set_mem_peak_num   => base_direct_solver_set_mem_peak_num
        procedure, non_overridable, public :: set_nz_factors     => base_direct_solver_set_nz_factors
        procedure, non_overridable, public :: set_Mflops         => base_direct_solver_set_Mflops
        procedure, non_overridable, public :: get_mem_peak_symb  => base_direct_solver_get_mem_peak_symb
        procedure, non_overridable, public :: get_mem_perm_symb  => base_direct_solver_get_mem_perm_symb
        procedure, non_overridable, public :: get_mem_peak_num   => base_direct_solver_get_mem_peak_num
        procedure, non_overridable, public :: get_nz_factors     => base_direct_solver_get_nz_factors
        procedure, non_overridable, public :: get_Mflops         => base_direct_solver_get_Mflops
        procedure, non_overridable, public :: log_info           => base_direct_solver_log_info
    end type

    interface
        subroutine base_direct_solver_set_defaults(this)
            import base_direct_solver_t
            class(base_direct_solver_t),   intent(inout) :: this
        end subroutine base_direct_solver_set_defaults

        subroutine base_direct_solver_set_parameters_from_pl(this, parameter_list)
            import base_direct_solver_t
            import ParameterList_t
            class(base_direct_solver_t),   intent(inout) :: this
            type(ParameterList_t),         intent(in)    :: parameter_list
        end subroutine base_direct_solver_set_parameters_from_pl

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
    end interface

public :: base_direct_solver_t

contains

    subroutine base_direct_solver_reset(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state  = BASE_DIRECT_SOLVER_STATE_START
        nullify(this%matrix)
        this%mem_peak_symb = 0
        this%mem_perm_symb = 0
        this%nz_factors    = 0
        this%mem_peak_num  = 0
        this%Mflops        = 0._rp
    end subroutine base_direct_solver_reset

    subroutine base_direct_solver_set_name(this, name)
        class(base_direct_solver_t),   intent(inout) :: this
        character(len=*),              intent(in)    :: name
        assert(this%state_is_start())
        this%name = name
    end subroutine base_direct_solver_set_name

    subroutine base_direct_solver_set_matrix(this, matrix)
        class(base_direct_solver_t),   intent(inout) :: this
        type(sparse_matrix_t), target, intent(in)    :: matrix
        assert(.not. this%state_is_symbolic())
        assert(.not. this%state_is_numeric())
        this%matrix => matrix
    end subroutine base_direct_solver_set_matrix

    subroutine base_direct_solver_update_matrix(this, matrix, same_nonzero_pattern)
    !-----------------------------------------------------------------
    !< Update matrix pointer 
    !< If same_nonzero_pattern performs returns to SYMBOLIC state
    !< If not same_nonzero_pattern returns to START state
    !-----------------------------------------------------------------
        class(base_direct_solver_t),   intent(inout) :: this
        type(sparse_matrix_t), target, intent(in)    :: matrix
        logical,                       intent(in)    :: same_nonzero_pattern
    !-----------------------------------------------------------------
        this%matrix => matrix
        if(same_nonzero_pattern .and. this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC) then
            call this%free_numerical()
        elseif(.not. same_nonzero_pattern .and. &
               (this%state == BASE_DIRECT_SOLVER_STATE_SYMBOLIC .or. &
                this%state == BASE_DIRECT_SOLVER_STATE_NUMERIC)) then
                call this%free_symbolic()
        endif
    end subroutine base_direct_solver_update_matrix

    function base_direct_solver_matrix_is_set(this) result(matrix_is_set)
        class(base_direct_solver_t), intent(inout) :: this
        logical                                    :: matrix_is_set
        matrix_is_set = associated(this%matrix)
    end function base_direct_solver_matrix_is_set

    subroutine base_direct_solver_set_state_start(this)
        class(base_direct_solver_t), intent(inout) :: this
        this%state = BASE_DIRECT_SOLVER_STATE_START
    end subroutine base_direct_solver_set_state_start

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

    subroutine base_direct_solver_set_mem_peak_num(this, mem_peak_num)
        class(base_direct_solver_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: mem_peak_num
        this%mem_peak_num = mem_peak_num
    end subroutine base_direct_solver_set_mem_peak_num

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

    function base_direct_solver_get_mem_peak_symb(this) result(mem_peak_symb)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_peak_symb
        mem_peak_symb = this%mem_peak_symb
    end function base_direct_solver_get_mem_peak_symb

    function base_direct_solver_get_mem_perm_symb(this) result(mem_perm_symb)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_perm_symb
        mem_perm_symb = this%mem_perm_symb
    end function base_direct_solver_get_mem_perm_symb

    function base_direct_solver_get_mem_peak_num(this) result(mem_peak_num)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: mem_peak_num
        mem_peak_num = this%mem_peak_num
    end function base_direct_solver_get_mem_peak_num

    function base_direct_solver_get_nz_factors(this) result(nz_factors)
        class(base_direct_solver_t), intent(in) :: this
        integer(ip)                             :: nz_factors
        nz_factors = this%nz_factors
    end function base_direct_solver_get_nz_factors

    function base_direct_solver_get_Mflops(this) result(Mflops)
        class(base_direct_solver_t), intent(in) :: this
        real(rp)                                :: Mflops
        Mflops = this%Mflops
    end function base_direct_solver_get_Mflops


    subroutine base_direct_solver_log_info(this)
    !-----------------------------------------------------------------
    !< Print info
    !-----------------------------------------------------------------
        class(base_direct_solver_t), intent(in) :: this
    !-----------------------------------------------------------------
        write (*,'(a)') ''
        write (*,'(a)') '---------------------------------------'
        write (*,'(a)') ' '//this%name//' DIRECT SOLVER LOG INFO'
        write (*,'(a)') '---------------------------------------'
        write (*,'(a,i10)') 'Peak mem.      in KBytes (symb fact) = ', this%get_mem_peak_symb()
        write (*,'(a,i10)') 'Permanent mem. in KBytes (symb fact) = ', this%get_mem_perm_symb()
        write (*,'(a,i10)') 'Size of factors (thousands)          = ', this%get_nz_factors()
        write (*,'(a,i10)') 'Peak mem.      in KBytes (num fact)  = ', this%get_mem_peak_num()
        write (*,'(a,f10.2)') 'MFlops for factorization             = ', this%get_Mflops() 
        write (*,'(a)') ''
    end subroutine base_direct_solver_log_info

end module base_direct_solver_names
