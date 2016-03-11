module umfpack_direct_solver_names

    USE types_names
    USE memor_names
!    USE iso_c_binding
    USE sparse_matrix_names
    USE base_sparse_matrix_names
    USE csr_sparse_matrix_names
    USE serial_scalar_array_names
    USE base_direct_solver_names
    USE umfpack_interface_names
    USE FPL

implicit none
# include "debug.i90"
  
private

    ! Parameter strings to be used in the Parameter List
    character(len=*), parameter :: umfpack_name             = 'UMFPACK'
    character(len=*), parameter :: umfpack_control_params   = 'umfpack_control_params'

    !< Numbering format of matrix column and row indices
    integer(ip),      parameter :: C_NUMBERING       = 0
    integer(ip),      parameter :: FORTRAN_NUMBERING = 1

    type, extends(base_direct_solver_t) :: umfpack_direct_solver_t
    private
        type(c_ptr)    :: Symbolic
        type(c_ptr)    :: Numeric
#ifdef ENABLE_UMFPACK
        real(c_double) :: Control(0:UMFPACK_CONTROL-1)            
        real(c_double) :: Info(0:UMFPACK_INFO-1)            
#endif
        integer(ip)    :: Matrix_numbering = FORTRAN_NUMBERING
    contains
    private
        procedure, public :: free_clean              => umfpack_direct_solver_free_clean
        procedure, public :: free_symbolic           => umfpack_direct_solver_free_symbolic
        procedure, public :: free_numerical          => umfpack_direct_solver_free_numerical
        procedure         :: initialize              => umfpack_direct_solver_initialize
        procedure, public :: set_defaults            => umfpack_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl  => umfpack_direct_solver_set_parameters_from_pl
        procedure         :: Fortran_to_C_numbering  => umfpack_direct_solver_Fortran_to_C_numbering
        procedure         :: C_to_Fortran_numbering  => umfpack_direct_solver_C_to_Fortran_numbering
        procedure, public :: symbolic_setup          => umfpack_direct_solver_symbolic_setup
        procedure, public :: numerical_setup         => umfpack_direct_solver_numerical_setup
        procedure, public :: solve                   => umfpack_direct_solver_solve
#ifndef ENABLE_UMFPACK
        procedure         :: not_enabled_error       => umfpack_direct_solver_not_enabled_error
#endif
    end type

public :: create_umfpack_direct_solver, umfpack_name, umfpack_control_params

contains

    function create_umfpack_direct_solver() result(umfpack_direct_solver)
    !-----------------------------------------------------------------
    !< Creational function for umfpack direct solver
    !-----------------------------------------------------------------
        class(base_direct_solver_t),   pointer :: umfpack_direct_solver
        type(umfpack_direct_solver_t), pointer :: umfpack
    !-----------------------------------------------------------------
        allocate(umfpack)
        call umfpack%set_name(umfpack_name)
        call umfpack%initialize()
        call umfpack%set_defaults()
        umfpack_direct_solver => umfpack
    end function create_umfpack_direct_solver

    subroutine umfpack_direct_solver_initialize(this)
    !-----------------------------------------------------------------
    !< Change STATE to direct_solver_STATE_START
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert(this%state_is_start())
        call this%reset()
        this%Matrix_numbering = FORTRAN_NUMBERING
        call this%set_state_start()
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_initialize


    subroutine umfpack_direct_solver_set_defaults(this)
    !-----------------------------------------------------------------
    !< Set UMFPACK default parameters
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        !  Set the default control parameters.
        call umfpack_di_defaults (this%Control)
#else
        call this%not_enabled_error
#endif
    end subroutine umfpack_direct_solver_set_defaults


    subroutine umfpack_direct_solver_set_parameters_from_pl(this, parameter_list)
    !-----------------------------------------------------------------
    !< Set UMPACK parameters from a given ParameterList
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t),  intent(inout) :: this
        type(ParameterList_t),           intent(in)    :: parameter_list
        integer(ip)                                    :: FPLError
        logical                                        :: is_present
        logical                                        :: same_data_type
        integer(ip), allocatable                       :: shape(:)
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
#ifdef DEBUG
        is_present     = parameter_list%isPresent(Key=umfpack_control_params)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=umfpack_control_params, mold=this%Control)
            shape          = parameter_list%getshape(Key=umfpack_control_params)
            if(same_data_type .and. size(shape) == 1) then
                if(shape(1) == UMFPACK_CONTROL) then
#endif
                    ! UMFPACK control parameters
                    FPLError = parameter_list%Get(Key=umfpack_control_params, Value=this%Control)
                    assert(FPLError == 0)
#ifdef DEBUG
                else
                    write(*,'(a)') ' Warning! pardiso_mkl_iparam ignored. Expected size (20). '
                endif
            else
                write(*,'(a)') ' Warning! pardiso_mkl_iparam ignored. Wrong data type or shape. '
            endif
        endif
#endif
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_set_parameters_from_pl


    subroutine umfpack_direct_solver_symbolic_setup(this)
    !-----------------------------------------------------------------
    !< Pre-orders the columns of the matrix to reduce fill-in
    !< adn performs a symbolic analysis
    !< Set STATE direct_solver_STATE_SYMBOLIC
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
        integer(ip)                                   :: status
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert ( .not. this%matrix%get_symmetric_storage() )
        assert (this%state_is_START() .or. this%state_is_symbolic() .or. this%state_is_numeric())
        check(this%matrix_is_set())

        matrix => this%matrix%get_pointer_to_base_matrix()

        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Fortran to C numbering 
                call this%Fortran_to_C_numbering()
                !  From the matrix data, create the symbolic factorization information.
                status = umfpack_di_symbolic (matrix%get_num_rows(), &  !< Number of rows of the transposed CSC (CSR) matrix
                                              matrix%get_num_cols(), &  !< Number of columns of the transposed CSC (CSR) matrix
                                              matrix%irp,            &  !< Transposed CSC column (CSR rows) pointers
                                              matrix%ja,             &  !< Transposed CSC row (CSR columns) indices
                                              C_NULL_PTR,            &  !< Numerical values
                                              this%Symbolic,         &  !< Opaque symbolic object
                                              this%Control,          &  !< UMPACK control parameters array
                                              this%Info)                !< UMPACK info array
                if ( status < 0 ) then
                    write ( *, '(a)' ) ''
                    write ( *, '(a)' ) 'UMFPACK - Fatal error!'
                    write ( *, '(a,i10)' ) '  umfpack_di_symbolic returns STATUS = ', status
                    check ( status == UMFPACK_OK )  
                end if
                ! C to Fortran numbering 
                call this%C_to_Fortran_numbering()
                ! Memmory usage statistics
                call this%set_mem_peak_symb( int((this%Info(UMFPACK_SYMBOLIC_PEAK_MEMORY)*this%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp) )
                call this%set_mem_perm_symb( int((this%Info(UMFPACK_SYMBOLIC_SIZE)*this%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp) )
            class DEFAULT
                check(.false.)
        end select
        call this%set_state_symbolic()    
#else
        call this%not_enabled_error
#endif
    end subroutine umfpack_direct_solver_symbolic_setup


    subroutine umfpack_direct_solver_numerical_setup(this)
    !-----------------------------------------------------------------
    !< Numerically scales and then factorizes a sparse matrix 
    !< into the produc LU
    !< set direct_solver_STATE_NUMERICAL_SETUP state
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
        real(rp),                    pointer          :: val(:)
        integer(ip)                                   :: status
#ifdef ENABLE_UMFPACK
        assert ( .not. this%matrix%get_symmetric_storage() )
        assert (this%state_is_symbolic() .or. this%state_is_numeric())

        matrix => this%matrix%get_pointer_to_base_matrix()
        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Fortran to C numbering 
                call this%Fortran_to_C_numbering()
                val => matrix%val(:)
                !  From the symbolic factorization information, carry out the numeric factorization.
                status = umfpack_di_numeric ( matrix%irp,    &  !< Transposed CSC column (CSR rows) pointers
                                              matrix%ja,     &  !< Transposed CSC row (CSR columns) indices
                                              val,           &  !< Numerical values
                                              this%Symbolic, &  !< Opaque symbolic object
                                              this%Numeric,  &  !< Opaque numeric object
                                              this%Control,  &  !< UMPACK control parameters array
                                              this%Info)        !< UMPACK info array
          
                if ( status < 0 ) then
                    write ( *, '(a)' ) ''
                    write ( *, '(a)' ) 'UMFPACK - Fatal error!'
                    write ( *, '(a,i10)' ) '  umfpack_di_numeric returns status = ', status
                    check ( status == UMFPACK_OK ) 
                end if
                ! C to Fortran numbering 
                call this%C_to_Fortran_numbering()
                ! Memory usage statistics
                call this%set_mem_peak_num( int((this%Info(UMFPACK_PEAK_MEMORY)*this%Info(UMFPACK_SIZE_OF_UNIT))/1024.0_rp) )
                call this%set_Mflops( this%Info(UMFPACK_FLOPS)/(1.0e+06_rp * this%Info(UMFPACK_NUMERIC_TIME)) )
                call this%set_nz_factors( int((this%Info(UMFPACK_UNZ)+this%Info(UMFPACK_LNZ))/1.0e+03_rp) )
            class DEFAULT
                check(.false.)
        end select
        call this%set_state_numeric()    
#else
        call this%not_enabled_error
#endif

    end subroutine umfpack_direct_solver_numerical_setup


    subroutine umfpack_direct_solver_solve(op, x, y)
    !-----------------------------------------------------------------
    ! Computes y <- AT^-1 * x, using previously computed LU factorization
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: op
        class(serial_scalar_array_t),   intent(in)    :: x
        class(serial_scalar_array_t),   intent(inout) :: y
        class(base_sparse_matrix_t), pointer          :: matrix
        real(rp),                    pointer          :: val(:)
        real(rp),                    pointer          :: x_b(:)
        real(rp),                    pointer          :: y_b(:)
        integer(ip)                                   :: status
#ifdef ENABLE_UMFPACK
        assert ( .not. op%matrix%get_symmetric_storage() )
        assert (op%state_is_numeric())

        matrix => op%matrix%get_pointer_to_base_matrix()
        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Fortran to C numbering 
                call op%Fortran_to_C_numbering()
                val => matrix%val(:)
                x_b => x%b
                y_b => y%b
                ! Solve the linear system.
                status = umfpack_di_solve ( UMFPACK_At, &  !< A'x=b
                                            matrix%irp, &  !< Transposed CSC column (CSR rows) pointers
                                            matrix%ja,  &  !< Transposed CSC row (CSR columns) indices
                                            val,        &  !< Numerical values
                                            y_b,        &  !< Output X array
                                            x_b,        &  !< Input B array
                                            op%Numeric, &  !< Opaque numeric object
                                            op%Control, &  !< UMPACK control parameters array
                                            op%Info )      !< UMPACK info array

                if ( status < 0 ) then
                    write ( *, '(a)' ) ''
                    write ( *, '(a)' ) 'UMFPACK - Fatal error!'
                    write ( *, '(a,i10)' ) '  umfpack_di_solve returns status = ', status
                    check ( status == UMFPACK_OK ) 
                end if

                ! C to Fortran numbering 
                call op%C_to_Fortran_numbering()
        end select
#else
        call op%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_solve


    subroutine umfpack_direct_solver_free_clean(this)
    !-----------------------------------------------------------------
    !< Deallocate UMFPACK internal data structure
    !< Reset state to direct_solver_STATE_SYMBOLIC_START
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert(this%state_is_start() .or. this%state_is_symbolic() .or. this%state_is_numeric())
        call this%set_state_start()
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_free_clean


    subroutine umfpack_direct_solver_free_symbolic(this)
    !-----------------------------------------------------------------
    !< Release all internal memory for all matrices
    !< Set state to direct_solver_STATE_START
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert(this%state_is_symbolic() .or. this%state_is_numeric())
        call umfpack_di_free_symbolic ( this%Symbolic )
        call this%set_state_start()
#else
        call this%not_enabled_error()
#endif
    
    end subroutine umfpack_direct_solver_free_symbolic


    subroutine umfpack_direct_solver_free_numerical(this)
    !-----------------------------------------------------------------
    !< Release internal memory only for L and U factors
    !< Set state to direct_solver_STATE_SYMBOLIC_SETUP
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert(this%state_is_numeric())
        call umfpack_di_free_numeric ( this%Numeric )
        call this%set_state_symbolic()
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_free_numerical


    subroutine umfpack_direct_solver_Fortran_to_C_numbering(this)
    !-----------------------------------------------------------------
    !< Change matrix column and row indices to C style (starting in 0)
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        if(this%Matrix_Numbering == FORTRAN_NUMBERING) then
            matrix => this%matrix%get_pointer_to_base_matrix()
            select type (matrix)
                type is (csr_sparse_matrix_t)
                    matrix%irp = matrix%irp - 1
                    matrix%ja  = matrix%ja  - 1
                    this%Matrix_Numbering = C_NUMBERING
            end select
        endif
#else
        call this%not_enabled_error()
#endif
    end subroutine


    subroutine umfpack_direct_solver_C_to_Fortran_numbering(this)
    !-----------------------------------------------------------------
    !< Change matrix column and row indices to Fortran style (starting in 1)
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        if(this%Matrix_Numbering == C_NUMBERING) then
            matrix => this%matrix%get_pointer_to_base_matrix()
            select type (matrix)
                type is (csr_sparse_matrix_t)
                    matrix%irp = matrix%irp + 1
                    matrix%ja  = matrix%ja  + 1
                    this%Matrix_Numbering = FORTRAN_NUMBERING
            end select
        endif
#else
        call this%not_enabled_error()
#endif
    end subroutine


#ifndef ENABLE_UMFPACK
    subroutine umfpack_direct_solver_not_enabled_error(this)
        class(umfpack_direct_solver_t), intent(inout) :: this
        write (0,*) 'Error: Fem was not compiled with -DENABLE_UMFPACK.'
        write (0,*) "Error: You must activate this cpp macro in order to use UMFPACK"
        check(.false.)
  end subroutine
#endif

end module umfpack_direct_solver_names
