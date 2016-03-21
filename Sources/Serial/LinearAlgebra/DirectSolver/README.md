# Direct Solvers: Basic usage


```fortran
...
    type(sparse_matrix_t)          :: sparse_matrix
    type(ParameterList_t)          :: parameter_list
    type(ParameterList_t), pointer :: direct_solver_params
    type(direct_solver_t)          :: direct_solver
    type(serial_scalar_array_t)    :: x
    type(serial_scalar_array_t)    :: y
    integer                        :: FPLError
...
    call FPL_Init()                                                                        !
    call TheDirectSolverCreationalMethodsDictionary%Init()                                 ! Initialization. In a FEMPAR_INIT subroutine in the future
    call parameter_list%Init()                                                             !
...
    direct_solver_params => parameter_list%NewSubList(Key=PARDISO_MKL)                     ! Create a new parameter sublist to store Direct solver parameters
    direct_solver_params%set(key = DIRECT_SOLVER_TYPE, Value = PARDISO_MKL)                ! Set the DIRECT_SOLVER_TYPE parameter into the list
...
    call direct_solver%set_type_from_pl(direct_solver_params)                              ! Create direct solver type from parameter list
    call direct_solver%set_matrix(sparse_matrix)                                           ! Associate matrix and direct solver entities
    call direct_solver%set_parameters_from_pl(direct_solver_params)                        ! Set parameters from parameter list
    call direct_solver%update_matrix(sparse_matrix, same_nonzero_pattern=.true.)           ! Notify matrix changes to the direct solver
    call direct_solver%solve(x,y)                                                          ! Solve
    call direct_solver%log_info()
...

```

# Direct Solvers: Public Parameters

```fortran

    !-----------------------------------------------------------------
    ! Parameters used in DIRECT SOLVERS
    !-----------------------------------------------------------------
    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: DIRECT_SOLVER_TYPE = 'DIRECT_SOLVER_TYPE'
    character(len=*), parameter :: PARDISO_MKL        = 'PARDISO_MKL'                      ! Name of the PARDISO MKL direct solver type
    character(len=*), parameter :: UMFPACK            = 'UMFPACK'                          ! Name of the UMFPACK direct solver type

    !-----------------------------------------------------------------
    ! Parameters used in PARDISO_MKL direct solver
    !-----------------------------------------------------------------

    ! PARDISO MKL matrix types
    integer,          parameter :: PARDISO_MKL_SPD =  2                                    ! Real Symmetric positive definite 
    integer,          parameter :: PARDISO_MKL_SIN = -2                                    ! Real Symmetric indefinite
    integer,          parameter :: PARDISO_MKL_USS = 1                                     ! Real Unsymmetric, structurally symmetric
    integer,          parameter :: PARDISO_MKL_UNS = 11                                    ! Real Unsymmetric, structurally unsymmetric

    ! PARDISO MKL default values
    integer,          parameter :: PARDISO_MKL_DEFAULT_MESSAGE_LEVEL = 0                   ! Default pardiso_mkl message level
    integer,          parameter :: PARDISO_MKL_DEFAULT_IPARM         = 0                   ! Default pardiso_mkl_iparam(64) = 0
    integer,          parameter :: PARDISO_MKL_DEFAULT_MATRIX_TYPE   = pardiso_mkl_uns     ! Default pardiso_mkl matrix type

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: PARDISO_MKL_IPARM         = 'PARDISO_MKL_IPARM'         ! PARDISO MKL control parameters array
    character(len=*), parameter :: PARDISO_MKL_MATRIX_TYPE   = 'PARDISO_MKL_MATRIX_TYPE'   ! PARDISO MKL matrix type
    character(len=*), parameter :: PARDISO_MKL_MESSAGE_LEVEL = 'PARDISO_MKL_MESSAGE_LEVEL' ! PARDISO MKL verbosity level


    !-----------------------------------------------------------------
    ! Parameters used in UMFPACK direct solver
    !-----------------------------------------------------------------

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: UMFPACK_CONTROL_PARAMS   = 'UMFPACK_CONTROL_PARAMS'      ! UMFPACK real array of 20 UMFPACK parameters
```
