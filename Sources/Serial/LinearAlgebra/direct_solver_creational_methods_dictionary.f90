module direct_solver_creational_methods_dictionary_names

USE direct_solver_parameters,        only: PARDISO_MKL, UMFPACK
USE ProcedureDictionary,             only: ProcedureDictionary_t
USE base_direct_solver_names,        only: create_direct_solver_interface
USE pardiso_mkl_direct_solver_names, only: create_pardiso_mkl_direct_solver
USE umfpack_direct_solver_names,     only: create_umfpack_direct_solver

implicit none
private

    type :: DirectSolverCreationalMethodsDictionary_t
    private
        type(ProcedureDictionary_t) :: CreationalMethods
    contains
    private
        procedure, public :: Init          => DirectSolverCreationalMethodsDictionary_Init
        procedure, public :: isInitialized => DirectSolverCreationalMethodsDictionary_isInitialized
        procedure, public :: Set           => DirectSolverCreationalMethodsDictionary_Set
        procedure, public :: Get           => DirectSolverCreationalMethodsDictionary_Get
        procedure, public :: Free          => DirectSolverCreationalMethodsDictionary_Free
    end type

    type(DirectSolverCreationalMethodsDictionary_t) :: TheDirectSolverCreationalMethodsDictionary

public :: TheDirectSolverCreationalMethodsDictionary

contains
    subroutine DirectSolverCreationalMethodsDictionary_Init(this,Size)
    !-----------------------------------------------------------------
        class(DirectSolverCreationalMethodsDictionary_t), intent(INOUT) :: this
        integer, optional,                                intent(IN)    :: Size
    !-----------------------------------------------------------------
        call this%CreationalMethods%Free()
        call this%CreationalMethods%Init(Size = Size)
        call this%Set(Key=PARDISO_MKL, Proc=create_pardiso_mkl_direct_solver)
        call this%Set(Key=UMFPACK,     Proc=create_umfpack_direct_solver)
    end subroutine

    function DirectSolverCreationalMethodsDictionary_isInitialized(this) result(isInitialized)
    !-----------------------------------------------------------------
        class(DirectSolverCreationalMethodsDictionary_t), intent(INOUT) :: this
        logical                                                         :: isInitialized
    !-----------------------------------------------------------------
        isInitialized = this%CreationalMethods%isInitialized()
    end function

    subroutine DirectSolverCreationalMethodsDictionary_Set(this,Key,Proc)
    !-----------------------------------------------------------------
        class(DirectSolverCreationalMethodsDictionary_t), intent(INOUT) :: this
        character(len=*),                                 intent(IN)    :: Key
        procedure(create_direct_solver_interface)                       :: Proc
        procedure(create_direct_solver_interface), pointer              :: ProcPointer
    !-----------------------------------------------------------------
        ProcPointer => Proc
        call this%CreationalMethods%Set(Key=Key,Proc=ProcPointer)
    end subroutine

    subroutine DirectSolverCreationalMethodsDictionary_Get(this,Key,Proc)
    !-----------------------------------------------------------------
        class(DirectSolverCreationalMethodsDictionary_t),  intent(IN) :: this
        character(len=*),                                  intent(IN) :: Key
        procedure(create_direct_solver_interface), pointer            :: Proc
    !-----------------------------------------------------------------
        call this%CreationalMethods%Get(Key=Key,Proc=Proc)
    end subroutine

    subroutine DirectSolverCreationalMethodsDictionary_Free(this)
    !-----------------------------------------------------------------
        class(DirectSolverCreationalMethodsDictionary_t), intent(INOUT) :: this
    !-----------------------------------------------------------------
        call this%CreationalMethods%Free()
    end subroutine
end module direct_solver_creational_methods_dictionary_names

