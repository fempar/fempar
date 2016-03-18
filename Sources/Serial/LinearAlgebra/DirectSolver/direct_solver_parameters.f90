module direct_solver_parameters
implicit none

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

end module direct_solver_parameters
