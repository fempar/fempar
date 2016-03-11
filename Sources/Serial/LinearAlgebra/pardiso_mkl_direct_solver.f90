#ifdef ENABLE_MKL
  include 'mkl_pardiso.f90'
#endif
module pardiso_mkl_direct_solver_names
    ! Serial modules
    USE types_names
    USE memor_names
    USE sparse_matrix_names
    USE base_sparse_matrix_names, only: base_sparse_matrix_t
    USE csr_sparse_matrix_names, only: csr_sparse_matrix_t
    USE serial_scalar_array_names
    USE base_direct_solver_names
    USE FPL
#ifdef ENABLE_MKL
    USE mkl_pardiso
#endif

    implicit none

# include "debug.i90"

    private
    ! Parameters used in pardiso_mkl direct_solver
    integer, parameter :: dp              = 8   ! kind (1.0D0)
    integer, parameter :: pardiso_mkl_spd =  2  ! Real Symmetric positive definite 
    integer, parameter :: pardiso_mkl_sin = -2  ! Real Symmetric indefinite
    integer, parameter :: pardiso_mkl_uss = 1   ! Real Unsymmetric, structurally symmetric
    integer, parameter :: pardiso_mkl_uns = 11  ! Real Unsymmetric, structurally unsymmetric

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: pardiso_mkl_name             = 'PARDISO_MKL'
    character(len=*), parameter :: pardiso_mkl_iparm            = 'pardiso_mkl_iparm'
    character(len=*), parameter :: pardiso_mkl_matrix_type      = 'pardiso_mkl_matrix_type'
    character(len=*), parameter :: pardiso_mkl_message_level    = 'pardiso_mkl_message_level'

    integer, parameter :: default_pardiso_mkl_message_level = 0                ! Default pardiso_mkl parameters array
    integer, parameter :: default_pardiso_mkl_iparm(64)     = 0                ! Default pardiso_mkl parameters array
    integer, parameter :: default_pardiso_mkl_matrix_type   = pardiso_mkl_uns  ! Default pardiso_mkl parameters array

    type, extends(base_direct_solver_t) :: pardiso_mkl_direct_solver_t
    private
#ifdef ENABLE_MKL
        type ( mkl_pardiso_handle ) :: pardiso_mkl_pt(64)               !< Solver internal address pointer
#endif
        integer                     :: pardiso_mkl_iparm(64)
        integer                     :: phase                 = -1
        integer                     :: matrix_type           = -1500
        integer                     :: max_number_of_factors = 1
        integer                     :: actual_matrix         = 1
        integer                     :: number_of_rhs         = 1
        integer                     :: message_level         = 0
    contains
    private
        procedure, public :: free_clean              => pardiso_mkl_direct_solver_free_clean
        procedure, public :: free_symbolic           => pardiso_mkl_direct_solver_free_symbolic
        procedure, public :: free_numerical          => pardiso_mkl_direct_solver_free_numerical
        procedure         :: initialize              => pardiso_mkl_direct_solver_initialize
        procedure, public :: set_defaults            => pardiso_mkl_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl  => pardiso_mkl_direct_solver_set_parameters_from_pl
        procedure, public :: symbolic_setup          => pardiso_mkl_direct_solver_symbolic_setup
        procedure, public :: numerical_setup         => pardiso_mkl_direct_solver_numerical_setup
        procedure, public :: solve                   => pardiso_mkl_direct_solver_solve
#ifndef ENABLE_MKL
        procedure         :: not_enabled_error  => pardiso_mkl_direct_solver_not_enabled_error
#endif
    end type

public :: create_pardiso_mkl_direct_solver
public :: pardiso_mkl_name, pardiso_mkl_iparm, pardiso_mkl_matrix_type, pardiso_mkl_message_level
public :: pardiso_mkl_spd, pardiso_mkl_sin, pardiso_mkl_uss, pardiso_mkl_uns

contains

    function create_pardiso_mkl_direct_solver() result(pardiso_mkl_direct_solver)
    !-----------------------------------------------------------------
    !< Creational function for pardiso_mkl direct solver
    !-----------------------------------------------------------------
        class(base_direct_solver_t),       pointer :: pardiso_mkl_direct_solver
        type(pardiso_mkl_direct_solver_t), pointer :: pardiso_mkl
    !-----------------------------------------------------------------
        allocate(pardiso_mkl)
        call pardiso_mkl%set_name(pardiso_mkl_name)
        call pardiso_mkl%initialize()
        call pardiso_mkl%set_defaults()
        pardiso_mkl_direct_solver => pardiso_mkl
    end function create_pardiso_mkl_direct_solver


    subroutine pardiso_mkl_direct_solver_initialize(this)
    !-----------------------------------------------------------------
    !< Initiliaze the internal solver memory pointer.
    !< Set PARDISO parameters to the initial state
    !< Change STATE to direct_solver_STATE_START
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        integer(ip)                                        :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert(this%state_is_start())
        call this%reset()
        ! Initiliaze the internal solver memory pointer. This is only
        ! necessary before FIRST call of PARDISO.
        do i = 1, 64
            this%pardiso_mkl_pt(i)%dummy = 0 
            this%pardiso_mkl_iparm(i)    = 0
        end do
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_initialize


    subroutine pardiso_mkl_direct_solver_set_defaults(this)
    !-----------------------------------------------------------------
    !< Set PARDISO default parameters
    !< Choose pardiso_mkl matrix according to matrix if is set.
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        integer(ip)                                        :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        this%max_number_of_factors = 1
        this%actual_matrix         = 1
        this%number_of_rhs         = 1
        this%message_level         = default_pardiso_mkl_message_level
        this%pardiso_mkl_iparm     = default_pardiso_mkl_iparm
        this%matrix_type           = default_pardiso_mkl_matrix_type

        this%pardiso_mkl_iparm         = 0
        this%pardiso_mkl_iparm(18)     = -1 ! Output: number of nonzeros in the factor LU
        this%pardiso_mkl_iparm(19)     = -1 ! Output: Mflops for LU factorization
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_set_defaults


    subroutine pardiso_mkl_direct_solver_set_parameters_from_pl(this, parameter_list)
    !-----------------------------------------------------------------
    !< Set PARDISO parameters from a given ParameterList
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        type(ParameterList_t),               intent(in)    :: parameter_list
        integer(ip)                                        :: FPLError
        integer(ip)                                        :: matrix_type
        logical                                            :: is_present
        logical                                            :: same_data_type
        integer(ip), allocatable                           :: shape(:)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL  
        ! Matrix type
#ifdef DEBUG
        is_present     = parameter_list%isPresent(Key=pardiso_mkl_matrix_type)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=pardiso_mkl_matrix_type, mold=this%matrix_type)
            shape          = parameter_list%getshape(Key=pardiso_mkl_matrix_type)
            if(same_data_type .and. size(shape) == 0) then
#endif
                FPLError   = parameter_list%Get(Key=pardiso_mkl_matrix_type, Value=matrix_type)
                assert(FPLError == 0)
                if(this%state_is_symbolic()) then
                    ! Matrix cannot change in symbolic to numeric transition
                    assert(matrix_type == this%matrix_type)
                endif
                this%matrix_type = matrix_type
#ifdef DEBUG
            else
                write(*,'(a)') ' Warning! pardiso_mkl_matrix_type ignored. Wrong data type or shape. '
            endif
        endif

         ! iparm
        is_present     = parameter_list%isPresent(Key=pardiso_mkl_iparm)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=pardiso_mkl_iparm, mold=this%pardiso_mkl_iparm)
            shape          = parameter_list%getshape(Key=pardiso_mkl_iparm)
            if(same_data_type .and. size(shape) == 1) then
                if(shape(1) == 64) then
#endif
                    FPLError =  parameter_list%Get(Key=pardiso_mkl_iparm, Value=this%pardiso_mkl_iparm)
                    assert(FPLError == 0)
#ifdef DEBUG
                else
                    write(*,'(a)') ' Warning! pardiso_mkl_iparam ignored. Expected size (64). '
                endif
            else
                write(*,'(a)') ' Warning! pardiso_mkl_iparam ignored. Wrong data type or shape. '
            endif
        endif

         ! Message level
        is_present     = parameter_list%isPresent(Key=pardiso_mkl_message_level)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=pardiso_mkl_message_level, mold=this%message_level)
            shape          = parameter_list%getshape(Key=pardiso_mkl_message_level)
            if(same_data_type .and. size(shape) == 0) then
#endif
                FPLError   = parameter_list%Get(Key=pardiso_mkl_message_level, Value=this%message_level)
                assert(FPLError == 0)
#ifdef DEBUG
            else
                write(*,'(a)') ' Warning! pardiso_mkl_message_level ignored. Wrong data type or shape. '
            endif
        endif
#endif
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_set_parameters_from_pl


    subroutine pardiso_mkl_direct_solver_symbolic_setup(this)
    !-----------------------------------------------------------------
    !< Rewind to direct_solver_STATE_START state if necessary
    !< Perform PARDISO analysis step. Reordering and symbolic factorization, 
    !< this step also allocates all memory that is necessary for the factorization
    !< Set STATE direct_solver_STATE_SYMBOLIC
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer              :: matrix
        integer                                           :: error
        integer, target                                   :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        ! Check pre-conditions
        assert (this%state_is_start() .or. this%state_is_symbolic() .or. this%state_is_numeric())
        check(this%matrix_is_set())

        this%phase               = 11 ! only reordering and symbolic factorization
        matrix                   => this%matrix%get_pointer_to_base_matrix()

        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Reordering and symbolic factorization, this step also allocates 
                ! all memory that is necessary for the factorization
                call pardiso(pt     = this%pardiso_mkl_pt,         & !< Handle to internal data structure. The entries must be set to zero prior to the first call to pardiso
                             maxfct = this%max_number_of_factors,  & !< Maximum number of factors with identical sparsity structure that must be kept in memory at the same time
                             mnum   = this%actual_matrix,          & !< Actual matrix for the solution phase. The value must be: 1 <= mnum <= maxfct. 
                             mtype  = this%matrix_type,            & !< Defines the matrix type, which influences the pivoting method
                             phase  = this%phase,                  & !< Controls the execution of the solver (11 == Analysis)
                             n      = this%matrix%get_num_rows(),  & !< Number of equations in the sparse linear systems of equations
                             a      = ddum,                        & !< Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja
                             ia     = matrix%irp,                  & !< Pointers to columns in CSR format
                             ja     = matrix%ja,                   & !< Column indices of the CSR sparse matrix
                             perm   = idum,                        & !< Permutation vector
                             nrhs   = this%number_of_rhs,          & !< Number of right-hand sides that need to be solved for
                             iparm  = this%pardiso_mkl_iparm,      & !< This array is used to pass various parameters to Intel MKL PARDISO 
                             msglvl = this%message_level,          & !< Message level information
                             b      = ddum,                        & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                             x      = ddum,                        & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                             error  = error )

                if (error /= 0) then
                    write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                            error, 'during stage', this%phase
                    check(.false.)
                end if
            class DEFAULT
                check(.false.)
        end select

        call this%set_mem_peak_symb(this%pardiso_mkl_iparm(15))
        call this%set_mem_perm_symb(this%pardiso_mkl_iparm(16))
        call this%set_nz_factors(int(this%pardiso_mkl_iparm(18)/1e3))
        call this%set_state_symbolic()
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_symbolic_setup


    subroutine pardiso_mkl_direct_solver_numerical_setup(this)
    !-----------------------------------------------------------------
    !< Rewind to direct_solver_STATE_SYMBOLIC_SETUP state if necessary
    !< Perform numerical factorization
    !< set direct_solver_STATE_NUMERICAL_SETUP state
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer              :: matrix
        integer                                           :: error
        integer, target                                   :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert (this%state_is_start() .or. this%state_is_symbolic() .or. this%state_is_numeric())

        if(this%state_is_start()) call this%symbolic_setup()

        ! Factorization.
        this%phase = 22 ! only numerical factorization
        matrix => this%matrix%get_pointer_to_base_matrix()

        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Reordering and symbolic factorization, this step also allocates 
                ! all memory that is necessary for the factorization
                call pardiso(pt     = this%pardiso_mkl_pt,         & !< Handle to internal data structure. The entries must be set to zero prior to the first call to pardiso
                             maxfct = this%max_number_of_factors,  & !< Maximum number of factors with identical sparsity structure that must be kept in memory at the same time
                             mnum   = this%actual_matrix,          & !< Actual matrix for the solution phase. The value must be: 1 <= mnum <= maxfct. 
                             mtype  = this%matrix_type,            & !< Defines the matrix type, which influences the pivoting method
                             phase  = this%phase,                  & !< Controls the execution of the solver (22 == Numerical factorization)
                             n      = matrix%get_num_rows(),       & !< Number of equations in the sparse linear systems of equations
                             a      = matrix%val,                  & !< Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja
                             ia     = matrix%irp,                  & !< Pointers to columns in CSR format
                             ja     = matrix%ja,                   & !< Column indices of the CSR sparse matrix
                             perm   = idum,                        & !< Permutation vector
                             nrhs   = this%number_of_rhs,          & !< Number of right-hand sides that need to be solved for
                             iparm  = this%pardiso_mkl_iparm,      & !< This array is used to pass various parameters to Intel MKL PARDISO 
                             msglvl = this%message_level,          & !< Message level information
                             b      = ddum,                        & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                             x      = ddum,                        & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                             error  = error )

                if (error /= 0) then
                    write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                            error, 'during stage', this%phase, 'in row', this%pardiso_mkl_iparm(30)
                    check(.false.)
                end if

                ! Set pardiso mkl info
                call this%set_mem_peak_num(this%pardiso_mkl_iparm(16)+this%pardiso_mkl_iparm(17))
                call this%set_Mflops(real(this%pardiso_mkl_iparm(19))/1.0e3_rp)
            class DEFAULT
                check(.false.)
        end select
        call this%set_state_numeric()
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_numerical_setup


    subroutine pardiso_mkl_direct_solver_solve(op, x, y)
    !-----------------------------------------------------------------
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: op
        class(serial_scalar_array_t),       intent(in)    :: x
        class(serial_scalar_array_t),       intent(inout) :: y
        integer                                           :: error
        integer, target                                   :: idum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert (op%state_is_start() .or. op%state_is_symbolic() .or. op%state_is_numeric())
        assert (op%number_of_rhs == 1)

        if(.not. op%state_is_numeric()) call op%numerical_setup()

        ! (c) y  <- A^-1 * x
        op%phase = 33 ! only Fwd/Bck substitution

        select type (matrix => op%matrix%get_pointer_to_base_matrix())
            type is (csr_sparse_matrix_t)
                ! Solve, iterative refinement
                call pardiso(pt     = op%pardiso_mkl_pt,           & !< Handle to internal data structure. The entries must be set to zero prior to the first call to pardiso
                             maxfct = op%max_number_of_factors,    & !< Maximum number of factors with identical sparsity structure that must be kept in memory at the same time
                             mnum   = op%actual_matrix,            & !< Actual matrix for the solution phase. The value must be: 1 <= mnum <= maxfct. 
                             mtype  = op%matrix_type,              & !< Defines the matrix type, which influences the pivoting method
                             phase  = op%phase,                    & !< Controls the execution of the solver (33 == Solve, iterative refinement)
                             n      = matrix%get_num_rows(),       & !< Number of equations in the sparse linear systems of equations
                             a      = matrix%val,                  & !< Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja
                             ia     = matrix%irp,                  & !< Pointers to columns in CSR format
                             ja     = matrix%ja,                   & !< Column indices of the CSR sparse matrix
                             perm   = idum,                        & !< Permutation vector
                             nrhs   = op%number_of_rhs,            & !< Number of right-hand sides that need to be solved for
                             iparm  = op%pardiso_mkl_iparm,        & !< This array is used to pass various parameters to Intel MKL PARDISO 
                             msglvl = op%message_level,            & !< Message level information
                             b      = x%b,                         & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                             x      = y%b,                         & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                             error  = error )

                if (error /= 0) then
                    write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                            error, 'during stage', op%phase
                    check(.false.)
                end if
            class DEFAULT
                check(.false.)
        end select
#else
        call op%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_solve


    subroutine pardiso_mkl_direct_solver_free_clean(this)
    !-----------------------------------------------------------------
    !< Deallocate PARDISO internal data structure
    !< Reset state to direct_solver_STATE_SYMBOLIC_START
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert(this%state_is_start() .or. this%state_is_symbolic() .or. this%state_is_numeric())
        this%matrix_type         = -1500
        this%matrix              => NULL()
        call this%set_state_start()
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_clean


    subroutine pardiso_mkl_direct_solver_free_symbolic(this)
    !-----------------------------------------------------------------
    !< Release all internal memory for all matrices
    !< Set state to direct_solver_STATE_START
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        integer                                           :: error
        integer                                           :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert(this%state_is_symbolic() .or. this%state_is_numeric())
        this%phase = -1 ! Release all internal memory for all matrices
        call pardiso(pt     = this%pardiso_mkl_pt,             & !< Handle to internal data structure. The entries must be set to zero prior to the first call to pardiso
                     maxfct = this%max_number_of_factors,      & !< Maximum number of factors with identical sparsity structure that must be kept in memory at the same time
                     mnum   = this%actual_matrix,              & !< Actual matrix for the solution phase. The value must be: 1 <= mnum <= maxfct. 
                     mtype  = this%matrix_type,                & !< Defines the matrix type, which influences the pivoting method
                     phase  = this%phase,                      & !< Controls the execution of the solver (-1 == Release all internal memory for all matrices)
                     n      = this%matrix%get_num_rows(),      & !< Number of equations in the sparse linear systems of equations
                     a      = ddum,                            & !< Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja
                     ia     = idum,                            & !< Pointers to columns in CSR format
                     ja     = idum,                            & !< Column indices of the CSR sparse matrix
                     perm   = idum,                            & !< Permutation vector
                     nrhs   = this%number_of_rhs,              & !< Number of right-hand sides that need to be solved for
                     iparm  = this%pardiso_mkl_iparm,          & !< This array is used to pass various parameters to Intel MKL PARDISO 
                     msglvl = this%message_level,              & !< Message level information
                     b      = ddum,                            & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                     x      = ddum,                            & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                     error  = error )

        if (error /= 0) then
            write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                error, 'during stage', this%phase
            check(.false.)
        end if
        call this%set_state_start()

#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_symbolic


    subroutine pardiso_mkl_direct_solver_free_numerical(this)
    !-----------------------------------------------------------------
    !< Release internal memory only for L and U factors
    !< Set state to direct_solver_STATE_SYMBOLIC_SETUP
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        integer                                           :: error
        integer                                           :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert(this%state_is_numeric())
        ! Release internal memory only for L and U factors
        this%phase = 0 ! Release internal memory for L and U matrix number mnum
        call pardiso(pt     = this%pardiso_mkl_pt,         & !< Handle to internal data structure. The entries must be set to zero prior to the first call to pardiso
                     maxfct = this%max_number_of_factors,  & !< Maximum number of factors with identical sparsity structure that must be kept in memory at the same time
                     mnum   = this%actual_matrix,          & !< Actual matrix for the solution phase. The value must be: 1 <= mnum <= maxfct. 
                     mtype  = this%matrix_type,            & !< Defines the matrix type, which influences the pivoting method
                     phase  = this%phase,                  & !< Controls the execution of the solver (0 == Release internal memory for L and U matrix number mnum)
                     n      = this%matrix%get_num_rows(),  & !< Number of equations in the sparse linear systems of equations
                     a      = ddum,                        & !< Contains the non-zero elements of the coefficient matrix A corresponding to the indices in ja
                     ia     = idum,                        & !< Pointers to columns in CSR format
                     ja     = idum,                        & !< Column indices of the CSR sparse matrix
                     perm   = idum,                        & !< Permutation vector
                     nrhs   = this%number_of_rhs,          & !< Number of right-hand sides that need to be solved for
                     iparm  = this%pardiso_mkl_iparm,      & !< This array is used to pass various parameters to Intel MKL PARDISO 
                     msglvl = this%message_level,          & !< Message level information
                     b      = ddum,                        & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                     x      = ddum,                        & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                     error  = error )
        if (error /= 0) then
            write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                    error, 'during stage', this%phase
            check(.false.)
        end if

        call this%set_state_symbolic()
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_numerical


#ifndef ENABLE_MKL
    subroutine pardiso_mkl_direct_solver_not_enabled_error(this)
    !-----------------------------------------------------------------
    !< Show NOT_ENABLED error and stops execution
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(in) :: this
    !-----------------------------------------------------------------
        write (0,*) 'Error: Fempar was not compiled with -DENABLE_MKL.'
        write (0,*) "Error: You must activate this cpp macro in order to use Intel MKL's interface to PARDISO"
        check(.false.)
    end subroutine
#endif

end module pardiso_mkl_direct_solver_names
