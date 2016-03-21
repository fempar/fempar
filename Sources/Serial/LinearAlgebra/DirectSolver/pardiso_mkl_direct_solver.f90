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
    USE direct_solver_parameters_names
    USE FPL
#ifdef ENABLE_MKL
    USE mkl_pardiso
#endif

    implicit none

# include "debug.i90"

    private
    ! Parameters used in pardiso_mkl direct_solver
    integer, parameter :: dp              = 8   ! kind (1.0D0)

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
        logical                     :: forced_matrix_type    = .false.
    contains
    private
        procedure, public :: free_clean_body             => pardiso_mkl_direct_solver_free_clean_body
        procedure, public :: free_symbolic_body          => pardiso_mkl_direct_solver_free_symbolic_body
        procedure, public :: free_numerical_body         => pardiso_mkl_direct_solver_free_numerical_body
        procedure         :: initialize                  => pardiso_mkl_direct_solver_initialize
        procedure         :: set_defaults                => pardiso_mkl_direct_solver_set_defaults
        procedure         :: set_matrix_type_from_matrix => pardiso_mkl_direct_solver_set_matrix_type_from_matrix
        procedure, public :: set_parameters_from_pl      => pardiso_mkl_direct_solver_set_parameters_from_pl
        procedure, public :: symbolic_setup_body         => pardiso_mkl_direct_solver_symbolic_setup_body
        procedure, public :: numerical_setup_body        => pardiso_mkl_direct_solver_numerical_setup_body
        procedure, public :: solve_body                  => pardiso_mkl_direct_solver_solve_body
#ifndef ENABLE_MKL
        procedure         :: not_enabled_error           => pardiso_mkl_direct_solver_not_enabled_error
#endif
    end type

public :: create_pardiso_mkl_direct_solver

contains

    subroutine create_pardiso_mkl_direct_solver(base_direct_solver)
    !-----------------------------------------------------------------
    !< Creational function for pardiso_mkl direct solver
    !-----------------------------------------------------------------
        class(base_direct_solver_t),       pointer, intent(inout) :: base_direct_solver
        type(pardiso_mkl_direct_solver_t), pointer                :: pardiso_mkl_instance
    !-----------------------------------------------------------------
        assert(.not. associated(base_direct_solver))
        allocate(pardiso_mkl_instance)
        call pardiso_mkl_instance%set_name(pardiso_mkl)
        call pardiso_mkl_instance%initialize()
        call pardiso_mkl_instance%set_defaults()
        base_direct_solver => pardiso_mkl_instance
    end subroutine create_pardiso_mkl_direct_solver


    subroutine pardiso_mkl_direct_solver_initialize(this)
    !-----------------------------------------------------------------
    !< Initiliaze the internal solver memory pointer.
    !< Set PARDISO parameters to the initial state
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        integer(ip)                                        :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        call this%reset()
        this%forced_matrix_type = .false.
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
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        integer(ip)                                        :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        this%max_number_of_factors = 1
        this%actual_matrix         = 1
        this%number_of_rhs         = 1
        this%message_level         = PARDISO_MKL_DEFAULT_MESSAGE_LEVEL
        this%matrix_type           = PARDISO_MKL_DEFAULT_MATRIX_TYPE
        this%pardiso_mkl_iparm     = PARDISO_MKL_DEFAULT_IPARM

        this%pardiso_mkl_iparm(18)     = -1 ! Output: number of nonzeros in the factor LU
        this%pardiso_mkl_iparm(19)     = -1 ! Output: Mflops for LU factorization
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_set_defaults


    subroutine pardiso_mkl_direct_solver_set_matrix_type_from_matrix(this)
    !-----------------------------------------------------------------
    !< Choose pardiso_mkl matrix according to matrix if is set.
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t),  intent(inout) :: this
        integer(ip)                                        :: i
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        if(this%matrix_is_set()) then
            ! Choose pardiso_mkl matrix according to matrix
            if(this%matrix%get_symmetric_storage() .and. this%matrix%get_sign() == SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE) then
               this%matrix_type = pardiso_mkl_spd
            elseif(this%matrix%get_symmetric_storage() .and. this%matrix%get_sign() /= SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE) then
               this%matrix_type = pardiso_mkl_sin
            elseif(.not. this%matrix%get_symmetric_storage() .and. this%matrix%is_symmetric() ) then
               this%matrix_type = pardiso_mkl_uss
            else
               this%matrix_type = pardiso_mkl_uns
            end if
        endif
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_set_matrix_type_from_matrix


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
        is_present     = parameter_list%isPresent(Key=PARDISO_MKL_MATRIX_TYPE)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=PARDISO_MKL_MATRIX_TYPE, mold=this%matrix_type)
            FPLError       = parameter_list%getshape(Key=PARDISO_MKL_MATRIX_TYPE, shape=shape)
            if(same_data_type .and. size(shape) == 0) then
#endif
                FPLError   = parameter_list%Get(Key=PARDISO_MKL_MATRIX_TYPE, Value=matrix_type)
                assert(FPLError == 0)
                if(this%state_is_start()) then
                    ! Matrix cannot change in symbolic to numeric transition
                    this%matrix_type = matrix_type
                    this%forced_matrix_type = .true.
                else
                    write(*,'(a)') ' Warning! pardiso_mkl_matrix_type ignored. It cannot be changed after analysis phase'
                endif

#ifdef DEBUG
            else
                write(*,'(a)') ' Warning! pardiso_mkl_matrix_type ignored. Wrong data type or shape. '
            endif
        endif

         ! iparm
        is_present     = parameter_list%isPresent(Key=PARDISO_MKL_IPARM)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=PARDISO_MKL_IPARM, mold=this%pardiso_mkl_iparm)
            FPLError       = parameter_list%getshape(Key=PARDISO_MKL_IPARM, shape=shape)
            if(same_data_type .and. size(shape) == 1) then
                if(shape(1) == 64) then
#endif
                    FPLError =  parameter_list%Get(Key=PARDISO_MKL_IPARM, Value=this%pardiso_mkl_iparm)
                    assert(FPLError == 0)
#ifdef DEBUG
                else
                    write(*,'(a)') ' Warning! pardiso_mkl_iparm ignored. Expected size (64). '
                endif
            else
                write(*,'(a)') ' Warning! pardiso_mkl_iparm ignored. Wrong data type or shape. '
            endif
        endif

         ! Message level
        is_present     = parameter_list%isPresent(Key=PARDISO_MKL_MESSAGE_LEVEL)
        if(is_present) then
            same_data_type = parameter_list%isOfDataType(Key=PARDISO_MKL_MESSAGE_LEVEL, mold=this%message_level)
            FPLError       = parameter_list%getshape(Key=PARDISO_MKL_MESSAGE_LEVEL, shape=shape)
            if(same_data_type .and. size(shape) == 0) then
#endif
                FPLError   = parameter_list%Get(Key=PARDISO_MKL_MESSAGE_LEVEL, Value=this%message_level)
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


    subroutine pardiso_mkl_direct_solver_symbolic_setup_body(this)
    !-----------------------------------------------------------------
    !< Perform PARDISO analysis step. Reordering and symbolic factorization, 
    !< this step also allocates all memory that is necessary for the factorization
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer              :: matrix
        integer                                           :: error
        integer, target                                   :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
!        print*, '(1) --> symbolic_setup'
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
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_symbolic_setup_body


    subroutine pardiso_mkl_direct_solver_numerical_setup_body(this)
    !-----------------------------------------------------------------
    !< Perform numerical factorization
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        class(base_sparse_matrix_t), pointer              :: matrix
        integer                                           :: error
        integer, target                                   :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
!        print*, '(2) --> numerical_setup'
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
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_numerical_setup_body


    subroutine pardiso_mkl_direct_solver_solve_body(op, x, y)
    !-----------------------------------------------------------------
    ! Computes y <- A^-1 * x, using previously computed LU factorization
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: op
        type(serial_scalar_array_t),        intent(in)    :: x
        type(serial_scalar_array_t),        intent(inout) :: y
        real(rp), pointer                                 :: x_b(:)
        real(rp), pointer                                 :: y_b(:)
        integer                                           :: error
        integer,  target                                  :: idum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        assert (op%number_of_rhs == 1)
        call x%GuardTemp()
!        print*, '(3) --> solve'
        ! (c) y  <- A^-1 * x
        op%phase = 33 ! only Fwd/Bck substitution
        x_b => x%get_entries()
        y_b => y%get_entries()

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
                             b      = x_b,                         & !< Array, size (n, nrhs). On entry, contains the right-hand side vector/matrix
                             x      = y_b,                         & !< Array, size (n, nrhs). If iparm(6)=0 it contains solution vector/matrix X
                             error  = error )

                if (error /= 0) then
                    write (0,*) 'Error, PARDISO_MKL: the following ERROR was detected: ', & 
                            error, 'during stage', op%phase
                    check(.false.)
                end if
            class DEFAULT
                check(.false.)
        end select
        call x%CleanTemp()
#else
        call op%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_solve_body


    subroutine pardiso_mkl_direct_solver_free_clean_body(this)
    !-----------------------------------------------------------------
    !< Deallocate PARDISO internal data structure
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
!        print*, '(4) --> free_clean' 
        this%matrix_type         = -1500
        this%matrix              => NULL()
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_clean_body


    subroutine pardiso_mkl_direct_solver_free_symbolic_body(this)
    !-----------------------------------------------------------------
    !< Release all internal memory for all matrices
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        integer                                           :: error
        integer                                           :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
        if(.not. this%forced_matrix_type) call this%set_matrix_type_from_matrix()
!        print*, '(5) --> free_symbolic'
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
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_symbolic_body


    subroutine pardiso_mkl_direct_solver_free_numerical_body(this)
    !-----------------------------------------------------------------
    !< Release internal memory only for L and U factors
    !-----------------------------------------------------------------
        class(pardiso_mkl_direct_solver_t), intent(inout) :: this
        integer                                           :: error
        integer                                           :: idum(1)
        real(dp)                                          :: ddum(1)
    !-----------------------------------------------------------------
#ifdef ENABLE_MKL
!        print*, '(6) --> free_numerical'
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
#else
        call this%not_enabled_error()
#endif
    end subroutine pardiso_mkl_direct_solver_free_numerical_body


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
