! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module umfpack_direct_solver_names

    USE types_names
    USE memor_names
!    USE iso_c_binding
    USE sparse_matrix_names
    USE base_sparse_matrix_names
    USE csr_sparse_matrix_names
    USE serial_scalar_array_names
    USE base_direct_solver_names
    USE direct_solver_parameters_names
    USE umfpack_interface_names
    USE FPL

implicit none
# include "debug.i90"
  
private
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
        procedure, public :: free_clean_body           => umfpack_direct_solver_free_clean_body
        procedure, public :: free_symbolic_body        => umfpack_direct_solver_free_symbolic_body
        procedure, public :: free_numerical_body       => umfpack_direct_solver_free_numerical_body
        procedure         :: initialize                => umfpack_direct_solver_initialize
        procedure         :: set_defaults              => umfpack_direct_solver_set_defaults
        procedure, public :: set_parameters_from_pl    => umfpack_direct_solver_set_parameters_from_pl
        procedure         :: Fortran_to_C_numbering    => umfpack_direct_solver_Fortran_to_C_numbering
        procedure         :: C_to_Fortran_numbering    => umfpack_direct_solver_C_to_Fortran_numbering
        procedure, public :: symbolic_setup_body       => umfpack_direct_solver_symbolic_setup_body
        procedure, public :: numerical_setup_body      => umfpack_direct_solver_numerical_setup_body
        procedure, public :: solve_single_rhs_body     => umfpack_direct_solver_solve_single_rhs_body
        procedure, public :: solve_several_rhs_body    => umfpack_direct_solver_solve_several_rhs_body
#ifndef ENABLE_UMFPACK
        procedure         :: not_enabled_error       => umfpack_direct_solver_not_enabled_error
#endif
    end type

public :: create_umfpack_direct_solver

contains

    subroutine create_umfpack_direct_solver(umfpack_direct_solver)
    !-----------------------------------------------------------------
    !< Creational function for umfpack direct solver
    !-----------------------------------------------------------------
        class(base_direct_solver_t),   pointer, intent(inout) :: umfpack_direct_solver
        type(umfpack_direct_solver_t), pointer                :: umfpack_instance
    !-----------------------------------------------------------------
        allocate(umfpack_instance)
        call umfpack_instance%set_name(umfpack)
        call umfpack_instance%initialize()
        call umfpack_instance%set_defaults()
        umfpack_direct_solver => umfpack_instance
    end subroutine create_umfpack_direct_solver

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
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        if( parameter_list%isPresent(umfpack_control_params)) then
            ! UMFPACK control parameters
            assert(parameter_list%isAssignable(umfpack_control_params, this%control))
            FPLError = parameter_list%Get(Key=umfpack_control_params, Value=this%Control)
            assert(FPLError == 0)
        endif
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
                                              matrix%get_irp(),      &  !< Transposed CSC column (CSR rows) pointers
                                              matrix%get_ja(),       &  !< Transposed CSC row (CSR columns) indices
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
                val => matrix%get_val()
                !  From the symbolic factorization information, carry out the numeric factorization.
                status = umfpack_di_numeric ( matrix%get_irp(), &  !< Transposed CSC column (CSR rows) pointers
                                              matrix%get_ja(),  &  !< Transposed CSC row (CSR columns) indices
                                              val,              &  !< Numerical values
                                              this%Symbolic,    &  !< Opaque symbolic object
                                              this%Numeric,     &  !< Opaque numeric object
                                              this%Control,     &  !< UMPACK control parameters array
                                              this%Info)           !< UMPACK info array
          
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


    subroutine umfpack_direct_solver_solve_single_rhs_body(op, x, y)
    !-----------------------------------------------------------------
    ! Computes y <- AT^-1 * x, using previously computed LU factorization
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: op
        type(serial_scalar_array_t),    intent(in)    :: x
        type(serial_scalar_array_t),    intent(inout) :: y
        class(base_sparse_matrix_t), pointer          :: matrix
        real(rp),                    pointer          :: val(:)
        real(rp),                    pointer          :: x_b(:)
        real(rp),                    pointer          :: y_b(:)
        integer(ip)                                   :: status
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert ( .not. op%matrix%get_symmetric_storage() )
!        print*, '(3) --> solve'
        matrix => op%matrix%get_pointer_to_base_matrix()
        select type (matrix)
            type is (csr_sparse_matrix_t)
                ! Fortran to C numbering 
                call op%Fortran_to_C_numbering()
                val => matrix%get_val()
                x_b => x%get_entries()
                y_b => y%get_entries()
                ! Solve the linear system.
                status = umfpack_di_solve ( UMFPACK_At,       &  !< A'x=b
                                            matrix%get_irp(), &  !< Transposed CSC column (CSR rows) pointers
                                            matrix%get_ja(),  &  !< Transposed CSC row (CSR columns) indices
                                            val,              &  !< Numerical values
                                            y_b,              &  !< Output X array
                                            x_b,              &  !< Input B array
                                            op%Numeric,       &  !< Opaque numeric object
                                            op%Control,       &  !< UMPACK control parameters array
                                            op%Info )            !< UMPACK info array

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
    end subroutine umfpack_direct_solver_solve_single_rhs_body


    subroutine umfpack_direct_solver_solve_several_rhs_body(op, x, y)
    !-----------------------------------------------------------------
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    !-----------------------------------------------------------------
        class(umfpack_direct_solver_t), intent(inout) :: op
        real(rp),                       intent(inout) :: x(:, :)
        real(rp),                       intent(inout) :: y(:, :)
        class(base_sparse_matrix_t), pointer          :: matrix
        real(rp),                    pointer          :: val(:)
        integer(ip)                                   :: number_rows
        integer(ip)                                   :: number_rhs
        integer(ip)                                   :: status
        integer(ip)                                   :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        assert ( .not. op%matrix%get_symmetric_storage() )
!        print*, '(3) --> solve'
        matrix      => op%matrix%get_pointer_to_base_matrix()
        number_rows = size(x,1)
        number_rhs  = size(x,2)
        select type (matrix)
            type is (csr_sparse_matrix_t)
                assert(matrix%get_num_rows()==number_rows .and. size(y,1) == number_rows)
                assert(size(y,2) == number_rhs)
                ! Fortran to C numbering 
                call op%Fortran_to_C_numbering()
                val => matrix%get_val()

                do i=1, number_rhs
                    ! Solve the linear system.
                    status = umfpack_di_solve ( UMFPACK_At, &        !< A'x=b
                                                matrix%get_irp(), &  !< Transposed CSC column (CSR rows) pointers
                                                matrix%get_ja(),  &  !< Transposed CSC row (CSR columns) indices
                                                val,        &        !< Numerical values
                                                y(:,i),     &        !< Output X array
                                                x(:,i),     &        !< Input B array
                                                op%Numeric, &        !< Opaque numeric object
                                                op%Control, &        !< UMPACK control parameters array
                                                op%Info )            !< UMPACK info array

                    if ( status < 0 ) then
                        write ( *, '(a)' ) ''
                        write ( *, '(a)' ) 'UMFPACK - Fatal error!'
                        write ( *, '(a,i10)' ) '  umfpack_di_solve returns status = ', status
                        check ( status == UMFPACK_OK ) 
                    end if
                enddo

                ! C to Fortran numbering 
                call op%C_to_Fortran_numbering()
        end select
#else
        call op%not_enabled_error()
#endif
    end subroutine umfpack_direct_solver_solve_several_rhs_body


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
        integer(ip),                 pointer          :: irp(:)
        integer(ip),                 pointer          :: ja(:)
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        if(this%Matrix_Numbering == FORTRAN_NUMBERING) then
            matrix => this%matrix%get_pointer_to_base_matrix()
            select type (matrix)
                type is (csr_sparse_matrix_t)
                    irp => matrix%get_irp()
                    ja => matrix%get_ja()
                    irp = irp - 1
                    ja  = ja  - 1
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
        integer(ip),                 pointer          :: irp(:)
        integer(ip),                 pointer          :: ja(:)
    !-----------------------------------------------------------------
#ifdef ENABLE_UMFPACK
        if(this%Matrix_Numbering == C_NUMBERING) then
            matrix => this%matrix%get_pointer_to_base_matrix()
            select type (matrix)
                type is (csr_sparse_matrix_t)
                    irp => matrix%get_irp()
                    ja => matrix%get_ja()
                    irp = irp + 1
                    ja  = ja  + 1
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
