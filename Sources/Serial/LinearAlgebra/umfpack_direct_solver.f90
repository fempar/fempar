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
        procedure, public :: free_clean_body         => umfpack_direct_solver_free_clean_body
        procedure, public :: free_symbolic_body      => umfpack_direct_solver_free_symbolic_body
        procedure, public :: free_numerical_body     => umfpack_direct_solver_free_numerical_body
        procedure         :: initialize              => umfpack_direct_solver_initialize
        procedure         :: set_defaults            => umfpack_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl  => umfpack_direct_solver_set_parameters_from_pl
        procedure         :: Fortran_to_C_numbering  => umfpack_direct_solver_Fortran_to_C_numbering
        procedure         :: C_to_Fortran_numbering  => umfpack_direct_solver_C_to_Fortran_numbering
        procedure, public :: symbolic_setup_body     => umfpack_direct_solver_symbolic_setup_body
        procedure, public :: numerical_setup_body    => umfpack_direct_solver_numerical_setup_body
        procedure, public :: solve_body              => umfpack_direct_solver_solve_body
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
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        call this%reset()
        this%Matrix_numbering = FORTRAN_NUMBERING
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
            FPLError       = parameter_list%getshape(Key=umfpack_control_params, shape=shape)
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


    subroutine umfpack_direct_solver_symbolic_setup_body(this)
    !-----------------------------------------------------------------
    !< Pre-orders the columns of the matrix to reduce fill-in
    !< adn performs a symbolic analysis
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
        integer(ip)                                   :: status
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        ! Check pre-conditions
        assert ( .not. this%matrix%get_symmetric_storage() )

!        print*, '(1) --> symbolic_setup'
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
#else
        call this%not_enabled_error
#endif
    end subroutine umfpack_direct_solver_symbolic_setup_body


    subroutine umfpack_direct_solver_numerical_setup_body(this)
    !-----------------------------------------------------------------
    !< Numerically scales and then factorizes a sparse matrix 
    !< into the produc LU
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer          :: matrix
        real(rp),                    pointer          :: val(:)
        integer(ip)                                   :: status
#ifdef ENABLE_UMFPACK
        assert ( .not. this%matrix%get_symmetric_storage() )
!        print*, '(2) --> numerical_setup'
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
#else
        call this%not_enabled_error
#endif

    end subroutine umfpack_direct_solver_numerical_setup_body


    subroutine umfpack_direct_solver_solve_body(op, x, y)
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
!        print*, '(3) --> solve'
        matrix => op%matrix%get_pointer_to_base_matrix()
        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Fortran to C numbering 
                call op%Fortran_to_C_numbering()
                val => matrix%val(:)
                x_b => x%get_entries()
                y_b => y%get_entries()
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
    end subroutine umfpack_direct_solver_solve_body


    subroutine umfpack_direct_solver_free_clean_body(this)
    !-----------------------------------------------------------------
    !< Deallocate UMFPACK internal data structure
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
!        print*, '(4) --> free_clean'
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_free_clean_body


    subroutine umfpack_direct_solver_free_symbolic_body(this)
    !-----------------------------------------------------------------
    !< Release all internal memory for all matrices
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
!        print*, '(5) --> free_symbolic'
        call umfpack_di_free_symbolic ( this%Symbolic )
#else
        call this%not_enabled_error()
#endif
    
    end subroutine umfpack_direct_solver_free_symbolic_body


    subroutine umfpack_direct_solver_free_numerical_body(this)
    !-----------------------------------------------------------------
    !< Release internal memory only for L and U factors
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
!        print*, '(6) --> free_numeric'
        call umfpack_di_free_numeric ( this%Numeric )
#else
        call this%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_free_numerical_body


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
