module base_sparse_matrix_names

USE types_names
USE memor_names
USE vector_names

implicit none

# include "debug.i90"

private

!---------------------------------------------------------------------
!< BASE SPARSE MATRIX DERIVED  TYPE
!---------------------------------------------------------------------

    ! States
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_START           = 0 
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_CREATED         = 1 
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_BUILD_SYMBOLIC  = 2 
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_BUILD_NUMERIC   = 3  
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_ASSEMBLED       = 4 
    integer(ip), public, parameter :: SPARSE_MATRIX_STATE_UPDATE          = 5

    !-----------------------------------------------------------------
    ! State transition diagram for type(base_sparse_matrix_t)
    !-----------------------------------------------------------------
    ! Note: it is desireable that the state management occurs only
    !       inside this class to get a cleaner implementation
    !       of the son classes
    !-----------------------------------------------------------------
    ! Input State         | Action                | Output State 
    !-----------------------------------------------------------------
    ! Start               | Create                | Created
    ! Start               | Free                  | Start

    ! Created             | Insert (2-values)     | Build_symbolic
    ! Created             | Insert (3-values)     | Build_numeric
    ! Created             | Free                  | Start

    ! Build_symbolic      | Free                  | Start
    ! Build_symbolic      | Insert (2-values)     | Build_symbolic
    ! Build_symbolic      | Insert (2-values)     | * Error
    ! Build_symbolic      | convert               | Assembled

    ! Build_numeric       | Free                  | Start
    ! Build_numeric       | Insert (2-values)     | * Error
    ! Build_numeric       | Insert (3-values)     | Build_numeric
    ! Build_numeric       | convert               | Assembled

    ! Assembled           | Free                  | Start
    ! Assembled           | Set                   | Update

    ! Update              | Free                  | Start

  
    ! Matrix sign
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE     = 10
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE = 11
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_INDEFINITE            = 12 ! Both positive and negative eigenvalues
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_UNKNOWN               = 13 ! No info
    
    type, abstract :: base_sparse_matrix_t
    private 
        integer(ip) :: num_rows                          !< Number of rows
        integer(ip) :: num_cols                          !< Number of colums
        integer(ip) :: state = SPARSE_MATRIX_STATE_START !< Matrix state (one of SPARSE_MATRIX_STATE_XXX parameters)
        integer(ip) :: sign                              !< Matrix sign (one of SPARSE_MATRIX_SIGN_XXX parameters)
        logical     :: symmetric                         !< Matrix is symmetric (.true.) or not (.false.)
    contains
    private
        procedure(base_sparse_matrix_is_by_rows),            public, deferred :: is_by_rows
        procedure(base_sparse_matrix_is_by_cols),            public, deferred :: is_by_cols
        procedure(base_sparse_matrix_set_symmetric_storage), public, deferred :: set_symmetric_storage
        procedure(base_sparse_matrix_get_symmetric_storage), public, deferred :: get_symmetric_storage
        procedure(base_sparse_matrix_copy_to_coo),           public, deferred :: copy_to_coo
        procedure(base_sparse_matrix_copy_from_coo),         public, deferred :: copy_from_coo
        procedure(base_sparse_matrix_move_to_coo),           public, deferred :: move_to_coo
        procedure(base_sparse_matrix_move_from_coo),         public, deferred :: move_from_coo
        procedure(base_sparse_matrix_move_to_fmt),           public, deferred :: move_to_fmt
        procedure(base_sparse_matrix_move_from_fmt),         public, deferred :: move_from_fmt
        procedure(base_sparse_matrix_free_arrays),           public, deferred :: free_arrays
        procedure(base_sparse_matrix_get_nnz),               public, deferred :: get_nnz
        procedure(base_sparse_matrix_print),                 public, deferred :: print
        procedure         ::                             base_sparse_matrix_create_square
        procedure         ::                             base_sparse_matrix_create_rectangular
        procedure         ::                             base_sparse_matrix_append_entries
        procedure         ::                             base_sparse_matrix_append_values
        procedure         :: append_entries_body      => base_sparse_matrix_append_entries_body
        procedure         :: append_values_body       => base_sparse_matrix_append_values_body
        procedure         :: is_valid_sign            => base_sparse_matrix_is_valid_sign
        procedure         :: apply_body               => base_sparse_matrix_apply_body
        procedure, public :: copy_to_fmt              => base_sparse_matrix_copy_to_fmt
        procedure, public :: copy_from_fmt            => base_sparse_matrix_copy_from_fmt
        procedure, public :: set_sign                 => base_sparse_matrix_set_sign
        procedure, public :: get_sign                 => base_sparse_matrix_get_sign
        procedure, public :: set_num_rows             => base_sparse_matrix_set_num_rows
        procedure, public :: get_num_rows             => base_sparse_matrix_get_num_rows
        procedure, public :: set_num_cols             => base_sparse_matrix_set_num_cols
        procedure, public :: get_num_cols             => base_sparse_matrix_get_num_cols
        procedure, public :: set_symmetry             => base_sparse_matrix_set_symmetry
        procedure, public :: is_symmetric             => base_sparse_matrix_is_symmetric
        procedure, public :: set_state                => base_sparse_matrix_set_state
        procedure, public :: set_state_start          => base_sparse_matrix_set_state_start
        procedure, public :: set_state_created        => base_sparse_matrix_set_state_created
        procedure, public :: set_state_build_symbolic => base_sparse_matrix_set_state_build_symbolic
        procedure, public :: set_state_build_numeric  => base_sparse_matrix_set_state_build_numeric
        procedure, public :: set_state_assembled      => base_sparse_matrix_set_state_assembled
        procedure, public :: set_state_update         => base_sparse_matrix_set_state_update
        procedure, public :: get_state                => base_sparse_matrix_get_state
        procedure, public :: allocate_arrays          => base_sparse_matrix_allocate_arrays
        procedure, public :: apply                    => base_sparse_matrix_apply
        procedure, public :: free                     => base_sparse_matrix_free
        generic,   public :: create                   => base_sparse_matrix_create_square, &
                                                         base_sparse_matrix_create_rectangular
        generic,   public :: set_values               => base_sparse_matrix_append_entries, &
                                                         base_sparse_matrix_append_values
        generic           :: append_body              => append_entries_body, &
                                                         append_values_body
    end type


!---------------------------------------------------------------------
!< COO SPARSE MATRIX DERIVED TYPE
!---------------------------------------------------------------------

    integer(ip), parameter :: COO_SPARSE_MATRIX_SORTED_NONE    = 20
    integer(ip), parameter :: COO_SPARSE_MATRIX_SORTED_BY_ROWS = 21
    integer(ip), parameter :: COO_SPARSE_MATRIX_SORTED_BY_COLS = 22

    type, extends(base_sparse_matrix_t) :: coo_sparse_matrix_t
        integer(ip), private       :: sort_status = COO_SPARSE_MATRIX_SORTED_NONE ! Not sorted
        integer(ip), private       :: nnz = 0                     !< Number of non zeros
        logical,     private       :: symmetric_storage = .false. !< .True.   Implicitly assumes that G=(V,E) is such that 
                                                                  !<          (i,j) \belongs E <=> (j,i) \belongs E, forall i,j \belongs V.
                                                                  !<          Only edges (i,j) with j>=i are stored.
                                                                  !< .False.  All (i,j) \belongs E are stored.  
        integer(ip), allocatable   :: ia(:)                       !< Row indices
        integer(ip), allocatable   :: ja(:)                       !< Column indices        
        real(rp),    allocatable   :: val(:)                      !< Values
    contains
    private
        procedure         :: append_values_body      => coo_sparse_matrix_append_values_body
        procedure         :: append_entries_body     => coo_sparse_matrix_append_entries_body
        procedure, public :: is_by_rows              => coo_sparse_matrix_is_by_rows
        procedure, public :: is_by_cols              => coo_sparse_matrix_is_by_cols
        procedure, public :: set_nnz                 => coo_sparse_matrix_set_nnz
        procedure, public :: get_nnz                 => coo_sparse_matrix_get_nnz
        procedure, public :: sort_and_compress       => coo_sparse_matrix_sort_and_compress
        procedure, public :: set_sort_status_none    => coo_sparse_matrix_set_sort_status_none
        procedure, public :: set_sort_status_by_rows => coo_sparse_matrix_set_sort_status_by_rows
        procedure, public :: set_sort_status_by_cols => coo_sparse_matrix_set_sort_status_by_cols
        procedure, public :: get_sort_status         => coo_sparse_matrix_get_sort_status
        procedure, public :: set_symmetric_storage   => coo_sparse_matrix_set_symmetric_storage
        procedure, public :: get_symmetric_storage   => coo_sparse_matrix_get_symmetric_storage
        procedure, public :: allocate_arrays         => coo_sparse_matrix_allocate_arrays
        procedure, public :: copy_to_coo             => coo_sparse_matrix_copy_to_coo
        procedure, public :: copy_from_coo           => coo_sparse_matrix_copy_from_coo
        procedure, public :: copy_to_fmt             => coo_sparse_matrix_copy_to_fmt
        procedure, public :: copy_from_fmt           => coo_sparse_matrix_copy_from_fmt
        procedure, public :: move_to_coo             => coo_sparse_matrix_move_to_coo
        procedure, public :: move_from_coo           => coo_sparse_matrix_move_from_coo
        procedure, public :: move_to_fmt             => coo_sparse_matrix_move_to_fmt
        procedure, public :: move_from_fmt           => coo_sparse_matrix_move_from_fmt
        procedure, public :: free_arrays             => coo_sparse_matrix_free_arrays
        procedure, public :: print                   => coo_sparse_matrix_print
    end type coo_sparse_matrix_t

!---------------------------------------------------------------------
!< BASE SPARSE MATRIX INTERFACES
!---------------------------------------------------------------------

    interface
        function base_sparse_matrix_is_by_rows(this) result(by_rows)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(in) :: this
            logical                                 :: by_rows
        end function base_sparse_matrix_is_by_rows

        function base_sparse_matrix_is_by_cols(this) result(by_cols)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(in) :: this
            logical                                 :: by_cols
        end function base_sparse_matrix_is_by_cols

        function base_sparse_matrix_get_nnz(this) result(nnz)
            import base_sparse_matrix_t
            import ip    
            class(base_sparse_matrix_t), intent(in) :: this
            integer(ip)                             :: nnz
        end function base_sparse_matrix_get_nnz

        subroutine base_sparse_matrix_set_symmetric_storage(this, symmetric_storage)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            logical,                     intent(in)    :: symmetric_storage
        end subroutine base_sparse_matrix_set_symmetric_storage


        function base_sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(in) :: this
            logical                                 :: symmetric_storage
        end function base_sparse_matrix_get_symmetric_storage

        subroutine base_sparse_matrix_copy_to_coo(this, to)
            import base_sparse_matrix_t
            import coo_sparse_matrix_t
            class(base_sparse_matrix_t), intent(in)    :: this
            class(coo_sparse_matrix_t),  intent(inout) :: to
        end subroutine base_sparse_matrix_copy_to_coo

        subroutine base_sparse_matrix_copy_from_coo(this, from)
            import base_sparse_matrix_t
            import coo_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            class(coo_sparse_matrix_t), intent(in)     :: from
        end subroutine base_sparse_matrix_copy_from_coo


        subroutine base_sparse_matrix_move_to_coo(this, to)
            import base_sparse_matrix_t
            import coo_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            class(coo_sparse_matrix_t),  intent(inout) :: to
        end subroutine base_sparse_matrix_move_to_coo

        subroutine base_sparse_matrix_move_from_coo(this, from)
            import base_sparse_matrix_t
            import coo_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            class(coo_sparse_matrix_t),  intent(inout) :: from
        end subroutine base_sparse_matrix_move_from_coo

        subroutine base_sparse_matrix_move_to_fmt(this, to)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            class(base_sparse_matrix_t), intent(inout) :: to
        end subroutine base_sparse_matrix_move_to_fmt

        subroutine base_sparse_matrix_move_from_fmt(this, from)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t), intent(inout) :: this
            class(base_sparse_matrix_t), intent(inout) :: from
        end subroutine base_sparse_matrix_move_from_fmt

        subroutine base_sparse_matrix_free_arrays(this)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t),  intent(inout) :: this
        end subroutine base_sparse_matrix_free_arrays

        subroutine base_sparse_matrix_print(this,lunou, only_graph)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),  intent(in) :: this
            integer(ip),                  intent(in) :: lunou
            logical, optional,            intent(in) :: only_graph
        end subroutine base_sparse_matrix_print
    end interface

!---------------------------------------------------------------------
!< COO SPARSE MATRIX INTERFACES
!---------------------------------------------------------------------

    interface 
        subroutine duplicates_operation(input, output)
            import rp
            real(rp), intent(in)    :: input
            real(rp), intent(inout) :: output
        end subroutine duplicates_operation
    end interface

!---------------------------------------------------------------------
!< PUBLIC TYPES
!---------------------------------------------------------------------

public :: base_sparse_matrix_t
public :: coo_sparse_matrix_t

contains

!---------------------------------------------------------------------
!< BASE SPARSE MATRIX PROCEDURES
!---------------------------------------------------------------------

    subroutine base_sparse_matrix_set_state(this, state)
    !-----------------------------------------------------------------
    !< Set the matrix state. See state diagram
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: state
    !-----------------------------------------------------------------
        this%state = state
    end subroutine base_sparse_matrix_set_state


    subroutine base_sparse_matrix_set_state_start(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_START
    !< The matrix can be in this state in the initial stage
    !< of after a free
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_START
    end subroutine base_sparse_matrix_set_state_start


    subroutine base_sparse_matrix_set_state_created(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after SPARSE_MATRIX_STATE_START
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        check(this%state == SPARSE_MATRIX_STATE_START)
        this%state = SPARSE_MATRIX_STATE_CREATED
    end subroutine base_sparse_matrix_set_state_created


    subroutine base_sparse_matrix_set_state_build_symbolic(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after SPARSE_MATRIX_STATE_CREATE
    !< or SPARSE_MATRIX_STATE_BUILD_SYMBOLIC
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        check(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        this%state = SPARSE_MATRIX_STATE_BUILD_SYMBOLIC
    end subroutine base_sparse_matrix_set_state_build_symbolic


    subroutine base_sparse_matrix_set_state_build_numeric(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after SPARSE_MATRIX_STATE_CREATED
    !< or SPARSE_MATRIX_STATE_BUILD_NUMERIC
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        check(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC)
        this%state = SPARSE_MATRIX_STATE_BUILD_NUMERIC
    end subroutine base_sparse_matrix_set_state_build_numeric


    subroutine base_sparse_matrix_set_state_assembled(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_ASSEMBLED
    !< The matrix can be in this state after SPARSE_MATRIX_STATE_UPDATE
    !< or SPARSE_MATRIX_STATE_BUILD_SYMBOLIC or SPARSE_MATRIX_STATE_BUILD_NUMERIC
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        check(this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        this%state = SPARSE_MATRIX_STATE_ASSEMBLED
    end subroutine base_sparse_matrix_set_state_assembled


    subroutine base_sparse_matrix_set_state_update(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_ASSEMBLED
    !< The matrix can be in this state after SPARSE_MATRIX_STATE_ASSEMBLED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        check(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        this%state = SPARSE_MATRIX_STATE_update
    end subroutine base_sparse_matrix_set_state_update


    function base_sparse_matrix_get_state(this) result(state)
    !-----------------------------------------------------------------
    !< Get the matrix state. See state diagram
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip)                             :: state
    !-----------------------------------------------------------------
        state = this%state
    end function base_sparse_matrix_get_state


    subroutine base_sparse_matrix_set_sign(this, sign)
    !-----------------------------------------------------------------
    !< Set the sign of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: sign
    !-----------------------------------------------------------------
        this%sign = sign
    end subroutine base_sparse_matrix_set_sign


    function base_sparse_matrix_get_sign(this) result(sign)
    !-----------------------------------------------------------------
    !< Get the matrix sign
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip)                             :: sign
    !-----------------------------------------------------------------
        sign = this%sign
    end function base_sparse_matrix_get_sign


    subroutine base_sparse_matrix_set_symmetry(this, symmetric)
    !-----------------------------------------------------------------
    !< Set symmetry of the matrix. .true. for a symmetric matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        logical,                     intent(in)    :: symmetric
    !-----------------------------------------------------------------
        this%symmetric = symmetric
    end subroutine base_sparse_matrix_set_symmetry


    function base_sparse_matrix_is_symmetric(this) result(is_symmetric)
    !-----------------------------------------------------------------
    !< Return .true. for a symmetric matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: is_symmetric
    !-----------------------------------------------------------------
        is_symmetric = this%symmetric
    end function base_sparse_matrix_is_symmetric


    function base_sparse_matrix_is_valid_sign(this, sign) result(is_valid_sign)
    !-----------------------------------------------------------------
    !< Return .true. if sign is one of the allowed
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip),                 intent(in) :: sign
        logical                                 :: is_valid_sign
    !-----------------------------------------------------------------
        is_valid_sign = (sign == SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE .or. &
                         sign == SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE  .or. &
                         sign == SPARSE_MATRIX_SIGN_INDEFINITE .or. &
                         sign == SPARSE_MATRIX_SIGN_UNKNOWN )
    end function base_sparse_matrix_is_valid_sign


    subroutine base_sparse_matrix_set_num_rows(this, num_rows)
    !-----------------------------------------------------------------
    !< Set the number of rows of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
    !-----------------------------------------------------------------
        this%num_rows = num_rows
    end subroutine base_sparse_matrix_set_num_rows


    subroutine base_sparse_matrix_set_num_cols(this, num_cols)
    !-----------------------------------------------------------------
    !< Set the number of colums of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_cols
    !-----------------------------------------------------------------
        this%num_cols = num_cols
    end subroutine base_sparse_matrix_set_num_cols


    function base_sparse_matrix_get_num_rows(this) result(num_rows)
    !-----------------------------------------------------------------
    !< Get the number of rows of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip)                             :: num_rows
    !-----------------------------------------------------------------
        num_rows = this%num_rows
    end function base_sparse_matrix_get_num_rows


    function base_sparse_matrix_get_num_cols(this) result(num_cols)
    !-----------------------------------------------------------------
    !< Get the number of colums of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip)                             :: num_cols
    !-----------------------------------------------------------------
        num_cols = this%num_cols
    end function base_sparse_matrix_get_num_cols


    subroutine base_sparse_matrix_create_square(this,num_rows_and_cols,symmetric_storage,is_symmetric,sign, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a square matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows_and_cols
        logical,                     intent(in)    :: symmetric_storage
        logical,                     intent(in)    :: is_symmetric
        integer(ip),                 intent(in)    :: sign
        integer(ip), optional,       intent(in)    :: nz
    !-----------------------------------------------------------------
        call this%set_state_created()
        assert(this%is_valid_sign(sign))
        if(symmetric_storage) then
            check(is_symmetric)
        endif
        call this%set_symmetric_storage(symmetric_storage)
        this%symmetric = is_symmetric
        this%sign = sign    
        this%num_rows = num_rows_and_cols
        this%num_cols = num_rows_and_cols
        if(present(nz)) then
            call this%allocate_arrays(nz)
        else
            call this%allocate_arrays()
        endif
    end subroutine base_sparse_matrix_create_square


    subroutine base_sparse_matrix_create_rectangular(this,num_rows,num_cols, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: num_cols
        integer(ip), optional,       intent(in)    :: nz
    !-----------------------------------------------------------------
        call this%set_state_created()
        call this%set_symmetric_storage(.false.)
        this%symmetric = .false.
        this%sign = SPARSE_MATRIX_SIGN_UNKNOWN
        this%num_rows = num_rows
        this%num_cols = num_cols
        if(present(nz)) then
            call this%allocate_arrays(nz)
        else
            call this%allocate_arrays()
        endif
    end subroutine base_sparse_matrix_create_rectangular


    subroutine base_sparse_matrix_append_values(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< That must be overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, val, imin, imax, jmin, jmax)
        call this%set_state_build_numeric()
    end subroutine base_sparse_matrix_append_values


    subroutine base_sparse_matrix_append_entries(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< That must be overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, imin, imax, jmin, jmax)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_append_entries


    subroutine base_sparse_matrix_append_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_values_body


    subroutine base_sparse_matrix_append_entries_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_entries_body


    subroutine base_sparse_matrix_allocate_arrays(this, nz)
    !-----------------------------------------------------------------
    !< Allocate arrays
    !< Must be overloaded 
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip), optional,       intent(in)    :: nz
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_allocate_arrays


    subroutine base_sparse_matrix_apply(op, x, y)
    !-----------------------------------------------------------------
    !< Matrix vector product
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in)    :: op
        class(vector_t),             intent(in)    :: x
        class(vector_t),             intent(inout) :: y
    !-----------------------------------------------------------------
        check(op%state == SPARSE_MATRIX_STATE_ASSEMBLED)
        call op%apply_body(x, y)
    end subroutine base_sparse_matrix_apply


    subroutine base_sparse_matrix_apply_body(op, x, y)
    !-----------------------------------------------------------------
    !< Allocate arrays
    !< Must be overloaded 
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in)    :: op
        class(vector_t),             intent(in)    :: x
        class(vector_t),             intent(inout) :: y
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_apply_body


    subroutine base_sparse_matrix_copy_to_fmt(this, to)
    !-----------------------------------------------------------------
    !< Copy this (FTM) -> to (FMT)
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in)    :: this
        class(base_sparse_matrix_t), intent(inout) :: to
        type(coo_sparse_matrix_t)                  :: tmp
    !-----------------------------------------------------------------
        select type(to)
            type is (coo_sparse_matrix_t)
                call this%copy_to_coo(to)
            class default
                call this%copy_to_coo(tmp)
                call to%move_from_coo(tmp)
        end select
    end subroutine base_sparse_matrix_copy_to_fmt


    subroutine base_sparse_matrix_copy_from_fmt(this, from)
    !-----------------------------------------------------------------
    !< Copy from (FMT) -> this (FMT)
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        class(base_sparse_matrix_t), intent(in)    :: from
        type(coo_sparse_matrix_t)                  :: tmp
    !-----------------------------------------------------------------
        select type(from)
            type is (coo_sparse_matrix_t)
                call this%copy_from_coo(from)
            class default
                call from%copy_to_coo(tmp)
                call this%move_from_coo(tmp)
        end select
    end subroutine base_sparse_matrix_copy_from_fmt


    subroutine base_sparse_matrix_free (this)
    !-----------------------------------------------------------------
    !< Clean the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%sign = SPARSE_MATRIX_SIGN_UNKNOWN
        this%num_rows = 0
        this%num_cols = 0
        this%symmetric = .false.
        call this%set_symmetric_storage(.false.)
        call this%free_arrays()
        call this%set_state_start()
    end subroutine base_sparse_matrix_free


!---------------------------------------------------------------------
!< COO SPARSE MATRIX PROCEDURES
!---------------------------------------------------------------------

    subroutine coo_sparse_matrix_set_nnz(this, nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nnz
    !-----------------------------------------------------------------
        this%nnz = nnz
    end subroutine coo_sparse_matrix_set_nnz


    function coo_sparse_matrix_get_nnz(this) result(nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        integer(ip)                            :: nnz
    !-----------------------------------------------------------------
        nnz = this%nnz
    end function coo_sparse_matrix_get_nnz


    subroutine coo_sparse_matrix_set_symmetric_storage(this, symmetric_storage)
    !-----------------------------------------------------------------
    !< Set symmetry storage property of the matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        logical,                    intent(in)    :: symmetric_storage
    !-----------------------------------------------------------------
        this%symmetric_storage = symmetric_storage
    end subroutine coo_sparse_matrix_set_symmetric_storage


    function coo_sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
    !-----------------------------------------------------------------
    !< Get symmetric storage property of the matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        logical                                :: symmetric_storage
    !-----------------------------------------------------------------
        symmetric_storage = this%symmetric_storage
    end function coo_sparse_matrix_get_symmetric_storage


    subroutine coo_sparse_matrix_set_sort_status_none(this)
    !-----------------------------------------------------------------
    !< Set sort status to COO_SPARSE_MATRIX_SORT_STATUS_NONE (unsorted)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_set_sort_status_none


    subroutine coo_sparse_matrix_set_sort_status_by_rows(this)
    !-----------------------------------------------------------------
    !< Set sort status to COO_SPARSE_MATRIX_SORT_STATUS_BY_ROWS (sorted by rows)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_ROWS
    end subroutine coo_sparse_matrix_set_sort_status_by_rows


    subroutine coo_sparse_matrix_set_sort_status_by_cols(this)
    !-----------------------------------------------------------------
    !< Set sort status to COO_SPARSE_MATRIX_SORT_STATUS_BY_COLS (sorted by columns)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_COLS
    end subroutine coo_sparse_matrix_set_sort_status_by_cols


    function coo_sparse_matrix_get_sort_status(this) result(sort_status)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        integer(ip)                            :: sort_status
    !-----------------------------------------------------------------
        sort_status = this%sort_status
    end function coo_sparse_matrix_get_sort_status


    function coo_sparse_matrix_is_by_rows(this) result(is_by_rows)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        logical                                :: is_by_rows 
    !-----------------------------------------------------------------
        is_by_rows = (this%sort_status == COO_SPARSE_MATRIX_SORTED_BY_ROWS)
    end function coo_sparse_matrix_is_by_rows


    function coo_sparse_matrix_is_by_cols(this) result(is_by_cols)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by columns
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        logical                                :: is_by_cols
    !-----------------------------------------------------------------
        is_by_cols = (this%sort_status == COO_SPARSE_MATRIX_SORTED_BY_COLS)
    end function coo_sparse_matrix_is_by_cols


    subroutine coo_sparse_matrix_allocate_arrays(this, nz)
    !-----------------------------------------------------------------
    !< Allocate COO arrays
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        integer(ip), optional,      intent(in)     :: nz
    !-----------------------------------------------------------------
        if(present(nz)) then
            call memalloc(nz, this%ia,  __FILE__, __LINE__)
            call memalloc(nz, this%ja,  __FILE__, __LINE__)
            call memalloc(nz, this%val, __FILE__, __LINE__)
        else
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%ia,  __FILE__, __LINE__)
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%ja,  __FILE__, __LINE__)
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%val, __FILE__, __LINE__)
        endif
    end subroutine coo_sparse_matrix_allocate_arrays


    subroutine coo_sparse_matrix_append_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0) return

        nnz = this%nnz
        newnnz = nnz + nz 

        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Not yet used positions initialized to 0
        this%ia(nnz+1:) = 0
        this%ja(nnz+1:) = 0
        this%val(nnz+1:) = 0.0_rp
        ! Append the new entries and values
        do i=1, nz
            ir = ia(i); ic = ja(i)
            if(ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or. &
               (this%symmetric_storage .and. ir>ic) ) cycle
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ic
            this%val(nnz) = val(i)
            !write(*,*) 'Inserted:', ir, ic, val(i)
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_values_body


    subroutine coo_sparse_matrix_append_entries_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the COO sparse matrix
    !< It allows duplicates entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0) return

        nnz = this%nnz
        newnnz = nnz+nz

        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Not yet used positions initialized to 0
        this%ia(nnz+1:) = 0
        this%ja(nnz+1:) = 0
        ! Append the new entries
        do i=1, nz
            ir = ia(i); ic = ja(i)
            if(ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or. &
               (this%symmetric_storage .and. ir>ic) ) cycle
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ic
            !write(*,*) 'Inserted:', ir, ic
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_entries_body


    subroutine coo_sparse_matrix_sort_and_compress(this, by_cols, sum_duplicates)
    !-----------------------------------------------------------------
    !< Sort ia, ja and val by rows as default or by_cols if forced
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        logical,     optional,      intent(in)     :: by_cols
        logical,     optional,      intent(in)     :: sum_duplicates
        logical                                    :: sum_dupl
        integer(ip), allocatable                   :: ias(:)
        integer(ip), allocatable                   :: jas(:)
        real(rp),    allocatable                   :: vs(:)
        integer(ip), allocatable                   :: iaux(:)
        integer(ip), allocatable                   :: ix2(:)
        logical                                    :: by_rows 
        logical                                    :: sorted
        logical                                    :: use_buffers
        integer(ip)                                :: i, j, k, nnz, nzl, imx, iret, ip, is, irw, icl
        integer                                    :: info
        procedure(duplicates_operation), pointer   :: apply_duplicates => null ()
    !-----------------------------------------------------------------
        by_rows     = .true.
        sum_dupl    = .true.
        sorted      = .true.
        use_buffers = .true.
        if(present(by_cols)) by_rows = .not. by_cols

        if(present(sum_duplicates)) sum_dupl = sum_duplicates
        if(sum_dupl) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        nnz = this%nnz

        call memalloc(max(this%num_rows, this%num_cols, nnz) + 2, iaux, __FILE__, __LINE__)
        allocate(ias(nnz),  &
                 jas(nnz),  &
                 vs(nnz),   &
                 ix2(max(this%num_rows, this%num_cols, nnz) + 2), stat=info)

        use_buffers = (info == 0)
        iaux = 0

        if(by_rows) then
        ! By-rows sorting
            if(use_buffers) then
            ! If there is enough memory space
                if(this%ia(1) >= 1 .and. this%ia(1) <= this%num_rows) then
                    iaux(this%ia(1)) = iaux(this%ia(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        if(this%ia(i) < 1 .or. this%ia(i) > this%num_rows) then
                            use_buffers = .false.
                            exit
                        endif
                        iaux(this%ia(i)) = iaux(this%ia(i)) + 1
                        sorted = sorted .and. (this%ia(i-1) <= this%ia(i))
                    enddo
                else
                    use_buffers = .false.
                endif
            endif

            if(use_buffers) then
            ! If there is enough memory space
                if(sorted) then
                ! If rows are already sorted
                    k = 0
                    i = 1
                    do j = 1, this%num_rows
                        nzl = iaux(j)
                        imx = i + nzl - 1
                        if(nzl > 0)  then
                            ! Sort the colums of a particular row
                            call mergesort(nzl,this%ja(i:imx),ix2, iret)
                            if(iret == 0) call reorder(nzl, this%val(i:imx), this%ia(i:imx), this%ja(i:imx), ix2)
                            k = k + 1
                            this%ia(k)  = this%ia(i)
                            this%ja(k)  = this%ja(i)
                            this%val(k) = this%val(i)
                            irw = this%ia(k)
                            icl = this%ja(k)
                            do 
                                i=i+1
                                if(i > imx) exit
                                ! If symmetric_storage, only upper triangle is stored
                                if(this%symmetric_storage .and. this%ia(i)>this%ja(i)) cycle
                                if(this%ia(i) == irw .and. this%ja(i) == icl) then
                                    ! Duplicated values action: assign the last value
                                    call apply_duplicates(input=this%val(i), output=this%val(k))
                                else
                                    k = k + 1
                                    this%ia(k)  = this%ia(i)
                                    this%ja(k)  = this%ja(i)
                                    this%val(k) = this%val(i)
                                    irw = this%ia(k)
                                    icl = this%ja(k)
                                endif
                            enddo
                        endif
                    enddo
                else
                ! If rows are NOT sorted
                    ip = iaux(1)
                    iaux(1) = 0
                    do i=2, this%num_rows
                        is = iaux(i)
                        iaux(i) = ip
                        ip = ip + is
                    enddo
                    iaux(this%num_rows + 1) = ip

                    do i=1, nnz
                        irw = this%ia(i)
                        ip = iaux(irw) + 1
                        ias(ip) = this%ia(i)
                        jas(ip) = this%ja(i)
                        vs(ip)  = this%val(i)
                        iaux(irw) = ip
                    enddo

                    k = 0
                    i = 1
                    do j=1, this%num_rows
                        nzl = iaux(j) - i + 1
                        imx = i + nzl - 1

                        if(nzl > 0) then
                            ! Sort the colums of a particular row
                            call mergesort(nzl,jas(i:imx),ix2, iret)
                            if(iret == 0) call reorder(nzl, vs(i:imx), ias(i:imx), jas(i:imx), ix2)
                            k = k + 1
                            this%ia(k)  = ias(i)
                            this%ja(k)  = jas(i)
                            this%val(k) = vs(i)
                            irw = this%ia(k)
                            icl = this%ja(k)
                            do 
                                i=i+1
                                if(i > imx) exit
                                ! If symmetric_storage, only upper triangle is stored
                                if(this%symmetric_storage .and. ias(i)>jas(i)) cycle
                                if(ias(i) == irw .and. jas(i) == icl) then
                                    ! Duplicated values: assign the last value
                                    call apply_duplicates(input=vs(i), output=this%val(k))
                                else
!                                    this%val(k) = vs(i)
                                    k = k + 1
                                    this%ia(k)  = ias(i)
                                    this%ja(k)  = jas(i)
                                    this%val(k) = vs(i)
                                    irw = this%ia(k)
                                    icl = this%ja(k)
                                endif
                            enddo
                        endif
                    enddo
                endif
                i = k
            else
            ! If there is'n enough memory space
                ! Sort the rows
                call mergesort(nnz,this%ia, iaux, iret)
                if(iret == 0) call reorder(nnz, this%val, this%ia, this%ja, iaux)
                i = 1
                j = i
                do while (i <= nnz)
                    do while (this%ia(j) == this%ia(i))
                        j = j + 1
                        if(j > nnz) exit
                    enddo
                    nzl = j - i
                    ! Sort the colums of a particular row
                    call mergesort(nzl, this%ja(i:), iaux, iret)
                    if(iret == 0) call reorder(nzl, this%val(i:i+nzl-1), this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                    i = j
                enddo

                i = 1
                j = 1
                irw = this%ia(i)
                icl = this%ja(i)

                do 
                    j = j + 1
                    if(j > nnz) exit
                    ! If symmetric_storage, only upper triangle is stored
                    if(this%symmetric_storage .and. this%ia(j)>this%ja(j)) cycle
                    if(this%ia(j) == irw .and. this%ja(j) == icl) then
                        ! Duplicated values: assign the last value
                        call apply_duplicates(input=this%val(j), output=this%val(i))
!                        this%val(i) = this%val(j)
                    else
                        i = i + 1
                        this%ia(i)  = this%ia(j)
                        this%ja(i)  = this%ja(j)
                        this%val(i) = this%val(j)
                        irw = this%ia(i)
                        icl = this%ja(i)
                    endif
                enddo
            endif
            this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_ROWS
        else
        ! By-columns sorting
            if(use_buffers) then
            ! If there is enough memory space
                if(this%ja(1) >= 1 .and. this%ja(1) <= this%num_cols) then
                    iaux(this%ja(1)) = iaux(this%ja(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        if(this%ja(i) < 1 .or. this%ja(i) > this%num_cols) then
                            use_buffers = .false.
                            exit
                        endif
                        iaux(this%ja(i)) = iaux(this%ja(i)) + 1
                        sorted = sorted .and. (this%ja(i-1) <= this%ja(i))
                    enddo
                else
                    use_buffers = .false.
                endif
            endif

            if(use_buffers) then
            ! If there is enough memory space
                if(sorted) then
                ! If columns are already sorted
                    k = 0
                    i = 1
                    do j = 1, this%num_cols
                        nzl = iaux(j)
                        imx = i + nzl - 1
                        if(nzl > 0)  then
                            ! Sort the rows of a particular columns
                            call mergesort(nzl,this%ia(i:imx),ix2, iret)
                            if(iret == 0) call reorder(nzl, this%val(i:imx), this%ia(i:imx), this%ja(i:imx), ix2)
                            k = k + 1
                            this%ia(k)  = this%ia(i)
                            this%ja(k)  = this%ja(i)
                            this%val(k) = this%val(i)
                            irw = this%ia(k)
                            icl = this%ja(k)
                            do 
                                i=i+1
                                if(i > imx) exit
                                ! If symmetric_storage, only lower triangle is stored
                                if(this%symmetric_storage .and. this%ja(i)>this%ia(i)) cycle
                                if(this%ia(i) == irw .and. this%ja(i) == icl) then
                                    ! Duplicated values action: assign the last value or sum the values
                                    call apply_duplicates(input=this%val(i), output=this%val(k))
                                else
                                    k = k + 1
                                    this%ia(k)  = this%ia(i)
                                    this%ja(k)  = this%ja(i)
                                    this%val(k) = this%val(i)
                                    irw = this%ia(k)
                                    icl = this%ja(k)
                                endif
                            enddo
                        endif
                    enddo
                else
                ! If columns are NOT sorted
                    ip = iaux(1)
                    iaux(1) = 0
                    do i=2, this%num_cols
                        is = iaux(i)
                        iaux(i) = ip
                        ip = ip + is
                    enddo
                    iaux(this%num_cols + 1) = ip

                    do i=1, nnz
                        icl = this%ja(i)
                        ip = iaux(icl) + 1
                        ias(ip) = this%ia(i)
                        jas(ip) = this%ja(i)
                        vs(ip)  = this%val(i)
                        iaux(icl) = ip
                    enddo

                    k = 0
                    i = 1
                    do j=1, this%num_cols
                        nzl = iaux(j) - i + 1
                        imx = i + nzl - 1

                        if(nzl > 0) then
                            ! Sort the rows of a particular column
                            call mergesort(nzl,ias(i:imx),ix2, iret)
                            if(iret == 0) call reorder(nzl, vs(i:imx), ias(i:imx), jas(i:imx), ix2)
                            k = k + 1
                            this%ia(k)  = ias(i)
                            this%ja(k)  = jas(i)
                            this%val(k) = vs(i)
                            irw = this%ia(k)
                            icl = this%ja(k)
                            do 
                                i=i+1
                                if(i > imx) exit
                                ! If symmetric_storage, only lower triangle is stored 
                                if(this%symmetric_storage .and. jas(i)>ias(i)) cycle
                                if(ias(i) == irw .and. jas(i) == icl) then
                                    ! Duplicated values: assign the last value or sum the values
                                    call apply_duplicates(input=vs(i), output=this%val(k))
                                else
                                    k = k + 1
                                    this%ia(k)  = ias(i)
                                    this%ja(k)  = jas(i)
                                    this%val(k) = vs(i)
                                    irw = this%ia(k)
                                    icl = this%ja(k)
                                endif
                            enddo
                        endif
                    enddo
                endif
                i = k
            else
            ! If there is'n enough memory space
                ! Sort the columns
                call mergesort(nnz,this%ja, iaux, iret)
                if(iret == 0) call reorder(nnz, this%val, this%ia, this%ja, iaux)
                i = 1
                j = i
                do while (i <= nnz)
                    do while (this%ja(j) == this%ja(i))
                        j = j + 1
                        if(j > nnz) exit
                    enddo
                    nzl = j - i
                    ! Sort the rows of a particular column
                    call mergesort(nzl, this%ia(i:), iaux, iret)
                    if(iret == 0) call reorder(nzl, this%val(i:i+nzl-1), this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                    i = j
                enddo

                i = 1
                j = 1
                irw = this%ia(i)
                icl = this%ja(i)

                do 
                    j = j + 1
                    if(j > nnz) exit
                    ! If symmetric_storage, only lower triangle is stored
                    if(this%symmetric_storage .and. this%ja(j)>this%ia(j)) cycle
                    if(this%ia(j) == irw .and. this%ja(j) == icl) then
                        ! Duplicated values: assign the last value or sum the values
                        call apply_duplicates(input=this%val(j), output=this%val(i))
                    else
                        i = i + 1
                        this%ia(i)  = this%ia(j)
                        this%ja(i)  = this%ja(j)
                        this%val(i) = this%val(j)
                        irw = this%ia(i)
                        icl = this%ja(i)
                    endif
                enddo
            endif
            this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_COLS
        endif

        if(allocated(iaux)) call memfree(iaux, __FILE__, __LINE__)
        if(allocated(ias)) then; deallocate(ias, stat=info); assert(info == 0); endif
        if(allocated(jas)) then; deallocate(jas, stat=info); assert(info == 0); endif
        if(allocated(vs))  then; deallocate(vs,  stat=info); assert(info == 0); endif
        this%nnz = i
        call memrealloc(nnz, this%ia,  __FILE__, __LINE__)
        call memrealloc(nnz, this%ja,  __FILE__, __LINE__)
        call memrealloc(nnz, this%val, __FILE__, __LINE__)
                    
    contains

        subroutine assign_value(input, output)
        !-------------------------------------------------------------
        ! Assign an input value to the autput
        !-------------------------------------------------------------
            use types_names
            real(rp), intent(in)    :: input
            real(rp), intent(inout) :: output
        !-------------------------------------------------------------
            output = input
        end subroutine assign_value


        subroutine sum_value(input, output)
        !-------------------------------------------------------------
        ! Sum an input value to the autput
        !-------------------------------------------------------------
            use types_names
            real(rp), intent(in)    :: input
            real(rp), intent(inout) :: output
        !-------------------------------------------------------------
            output = output + input
        end subroutine sum_value


        subroutine mergesort(n,k,l,iret)
        !-------------------------------------------------------------
        !   This subroutine sorts an integer array into ascending order.
        !
        ! Arguments:
        !   n    -  integer           Input: size of the array 
        !   k    -  integer(*)        input: array of keys to be sorted
        !   l    -  integer(0:n+1)   output: link list 
        !   iret -  integer          output: 0 Normal termination
        !                                    1 the array was already sorted 
        !
        ! REFERENCES  = (1) D. E. Knuth
        !                   The Art of Computer Programming,
        !                     vol.3: Sorting and Searching
        !                   Addison-Wesley, 1973
        !-------------------------------------------------------------
            use types_names
            integer(ip) :: n, iret
            integer(ip) :: k(n),l(0:n+1)
            integer(ip) :: p,q,s,t
        !-------------------------------------------------------------
            iret = 0
            !  first step: we are preparing ordered sublists, exploiting
            !  what order was already in the input data; negative links
            !  mark the end of the sublists
            l(0) = 1
            t = n + 1
            do  p = 1,n - 1
                if (k(p) <= k(p+1)) then
                    l(p) = p + 1
                else
                    l(t) = - (p+1)
                    t = p
                end if
            end do
            l(t) = 0
            l(n) = 0
            ! see if the input was already sorted
            if (l(n+1) == 0) then
                iret = 1
                return 
            else
                l(n+1) = abs(l(n+1))
            end if

            mergepass: do 
                ! otherwise, begin a pass through the list.
                ! throughout all the subroutine we have:
                !  p, q: pointing to the sublists being merged
                !  s: pointing to the most recently processed record
                !  t: pointing to the end of previously completed sublist
                s = 0
                t = n + 1
                p = l(s)
                q = l(t)
                if (q == 0) exit mergepass

                outer: do 

                    if (k(p) > k(q)) then 
                        l(s) = sign(q,l(s))
                        s = q
                        q = l(q)
                        if (q > 0) then 
                            do 
                                if (k(p) <= k(q)) cycle outer
                                s = q
                                q = l(q)
                                if (q <= 0) exit
                            end do
                        end if
                        l(s) = p
                        s = t
                        do 
                            t = p
                            p = l(p)
                            if (p <= 0) exit
                        end do

                    else 
                        l(s) = sign(p,l(s))
                        s = p
                        p = l(p)
                        if (p>0) then 
                            do 
                                if (k(p) > k(q)) cycle outer 
                                s = p
                                p = l(p)
                                if (p <= 0) exit
                            end do
                        end if
                        !  otherwise, one sublist ended, and we append to it the rest
                        !  of the other one.
                        l(s) = q
                        s = t
                        do 
                            t = q
                            q = l(q)
                            if (q <= 0) exit
                        end do
                    end if

                    p = -p
                    q = -q
                    if (q == 0) then
                        l(s) = sign(p,l(s))
                        l(t) = 0
                        exit outer 
                    end if
                end do outer
            end do mergepass

        end subroutine mergesort


        subroutine reorder(n,x,i1,i2,iaux)
        !-------------------------------------------------------------
        !  Reorder (an) input vector(s) based on a list sort output.
        !  Based on: D. E. Knuth: The Art of Computer Programming
        !            vol. 3: Sorting and Searching, Addison Wesley, 1973
        !            ex. 5.2.12
        !-------------------------------------------------------------
            use types_names
            integer(ip), intent(in) :: n
            integer(ip) :: iaux(0:*) 
            real(rp)   :: x(*)
            integer(ip) :: i1(*), i2(*)     
            integer(ip) :: lswap, lp, k, isw1, isw2
            complex(rp)  :: swap
        !-------------------------------------------------------------
            lp = iaux(0)
            k  = 1
            do 
                if ((lp == 0).or.(k>n)) exit
                do 
                    if (lp >= k) exit
                    lp = iaux(lp)
                end do
                swap     = x(lp)
                x(lp)    = x(k)
                x(k)     = swap
                isw1     = i1(lp)
                i1(lp)   = i1(k)
                i1(k)    = isw1
                isw2     = i2(lp)
                i2(lp)   = i2(k)
                i2(k)    = isw2
                lswap    = iaux(lp)
                iaux(lp) = iaux(k)
                iaux(k)  = lp
                lp = lswap 
                k  = k + 1
            enddo
            return
        end subroutine reorder

    end subroutine coo_sparse_matrix_sort_and_compress


    subroutine coo_sparse_matrix_copy_to_coo(this, to)
    !-----------------------------------------------------------------
    !< Copy this (COO) -> to (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in)    :: this
        class(coo_sparse_matrix_t), intent(inout) :: to
        integer(ip)                               :: nnz
    !-----------------------------------------------------------------
        nnz = this%nnz
        call to%free()
        if(this%num_rows == this%num_cols) then
            call to%create(num_rows_and_cols = this%num_rows,          &
                           symmetric_storage = this%symmetric_storage, &
                           is_symmetric      = this%symmetric,         &
                           sign              = this%sign,              &
                           nz                = nnz)
        else
            call to%create(num_rows = this%num_rows, &
                           num_cols = this%num_cols )
        endif
        to%nnz = nnz
        to%ia(1:nnz)  = this%ia(1:nnz)
        to%ja(1:nnz)  = this%ja(1:nnz)
        to%val(1:nnz) = this%val(1:nnz)
        to%state = this%state
    end subroutine coo_sparse_matrix_copy_to_coo


    subroutine coo_sparse_matrix_copy_from_coo(this, from)
    !-----------------------------------------------------------------
    !< Copy from (COO) -> this (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(in)    :: from
        integer(ip)                               :: nnz
    !-----------------------------------------------------------------
        nnz = from%get_nnz()
        call this%free()
        if(from%get_num_rows() == from%get_num_cols()) then
            call this%create(num_rows_and_cols = from%get_num_rows(),          &
                             symmetric_storage = from%get_symmetric_storage(), &
                             is_symmetric      = from%is_symmetric(),          &
                             sign              = from%get_sign(),              &
                             nz                = nnz)
        else
            call this%create(num_rows = from%get_num_rows(), &
                             num_cols = from%get_num_cols() )
        endif
        this%nnz = nnz
        this%ia(1:nnz)  = from%ia(1:nnz)
        this%ja(1:nnz)  = from%ja(1:nnz)
        this%val(1:nnz) = from%val(1:nnz)
        this%state = from%state
    end subroutine coo_sparse_matrix_copy_from_coo


    subroutine coo_sparse_matrix_copy_to_fmt(this, to)
    !-----------------------------------------------------------------
    !< Copy this (FTM) -> to (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),  intent(in)    :: this
        class(base_sparse_matrix_t), intent(inout) :: to
    !-----------------------------------------------------------------
        call to%copy_from_coo(from=this)
    end subroutine coo_sparse_matrix_copy_to_fmt


    subroutine coo_sparse_matrix_copy_from_fmt(this, from)
    !-----------------------------------------------------------------
    !< Copy from (FMT) -> this (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),  intent(inout) :: this
        class(base_sparse_matrix_t), intent(in)    :: from
    !-----------------------------------------------------------------
        call from%copy_to_coo(to=this)
    end subroutine coo_sparse_matrix_copy_from_fmt


    subroutine coo_sparse_matrix_move_to_coo(this, to)
    !-----------------------------------------------------------------
    !< Move this (COO) -> to (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(inout) :: to
        integer(ip)                               :: nnz
    !-----------------------------------------------------------------
        nnz = this%nnz
        call to%free()
        to%num_cols  = this%num_cols
        to%num_rows  = this%num_rows
        to%symmetric = this%symmetric
        to%symmetric_storage = this%symmetric_storage
        to%sign      = this%sign
        to%nnz       = this%nnz
        to%state = this%state
        call move_alloc(from=this%ia, to=to%ia)
        call move_alloc(from=this%ja, to=to%ja)
        call move_alloc(from=this%val, to=to%val)
        call this%free()
    end subroutine coo_sparse_matrix_move_to_coo


    subroutine coo_sparse_matrix_move_from_coo(this, from)
    !-----------------------------------------------------------------
    !< Move from (COO) -> this (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        class(coo_sparse_matrix_t), intent(inout) :: from
        integer(ip)                               :: nnz
    !-----------------------------------------------------------------
        nnz = this%nnz
        call this%free()
        this%num_cols = from%num_cols
        this%num_rows = from%num_rows
        this%symmetric = from%symmetric
        this%symmetric_storage = from%symmetric_storage
        this%sign = from%sign
        this%nnz = from%nnz
        this%state = from%state
        call move_alloc(from=from%ia, to=this%ia)
        call move_alloc(from=from%ja, to=this%ja)
        call move_alloc(from=from%val, to=this%val)
        call from%free()
    end subroutine coo_sparse_matrix_move_from_coo


    subroutine coo_sparse_matrix_move_to_fmt(this, to)
    !-----------------------------------------------------------------
    !< Move this (CSR) -> to (FMT)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),  intent(inout) :: this
        class(base_sparse_matrix_t), intent(inout) :: to
    !-----------------------------------------------------------------
        call to%move_from_coo(from=this)
    end subroutine coo_sparse_matrix_move_to_fmt


    subroutine coo_sparse_matrix_move_from_fmt(this, from)
    !-----------------------------------------------------------------
    !< Move from (FMT) -> this (COO)
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),  intent(inout) :: this
        class(base_sparse_matrix_t), intent(inout) :: from
    !-----------------------------------------------------------------
        call from%move_to_coo(to=this)
    end subroutine coo_sparse_matrix_move_from_fmt


    subroutine coo_sparse_matrix_free_arrays(this)
    !-----------------------------------------------------------------
    !< Clean COO sparse matrix format derived type
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
    !-----------------------------------------------------------------
        if(allocated(this%ia))  call memfree (this%ia, __FILE__, __LINE__)
        if(allocated(this%ja))  call memfree (this%ja, __FILE__, __LINE__)
        if(allocated(this%val)) call memfree (this%val, __FILE__, __LINE__)
    end subroutine coo_sparse_matrix_free_arrays


    subroutine coo_sparse_matrix_print(this,lunou, only_graph)
    !-----------------------------------------------------------------
    !< Print a COO matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),  intent(in) :: this
        integer(ip),                 intent(in) :: lunou
        logical, optional,           intent(in) :: only_graph
        logical                                 :: print_vals
        integer(ip)                             :: i,j
    !-----------------------------------------------------------------
        print_vals = .true.; if(present(only_graph)) print_vals = .not. only_graph
        write (lunou, '(a)')     '********************************************'
        write (lunou, '(a)')     '************* COO data structure ***********'
        write (lunou, '(a)')     '********************************************'
        write (lunou, '(a,i10)') 'Number of rows:', this%num_rows
        write (lunou, '(a,i10)') 'Number of cols:', this%num_cols
        write (lunou, '(a,i10)') 'Number of non zeros (nnz):', this%nnz
    
        write (lunou, '(a)')     'Rows list (ia):'
        if(allocated(this%ia)) then
            write (lunou, *)    this%ia(1:this%nnz)
        else
            write (lunou,'(A)') 'Not allocated'
        endif
    
        write (lunou, '(a)')      'Columns list (ja):'
        if(allocated(this%ja)) then
            write (lunou, *)    this%ja(1:this%nnz)
        else
            write (lunou,'(A)') 'Not allocated'
        endif

        if(print_vals) then
            write (lunou, '(a)')      'Values list (val):'
            if(allocated(this%val)) then
                write (lunou, *)    this%val(1:this%nnz)
            else
                write (lunou,'(A)') 'Not allocated'
            endif
        endif
    end subroutine coo_sparse_matrix_print

end module base_sparse_matrix_names
