module sparse_matrix_names

USE types_names
USE memor_names
USE vector_names
USE base_sparse_matrix_names, only: base_sparse_matrix_t
USE coo_sparse_matrix_names
USE csr_sparse_matrix_names

implicit none

# include "debug.i90"

private

    ! Matrix sign
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE     = 0
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE = 1
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_INDEFINITE            = 2 ! Both positive and negative eigenvalues
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_UNKNOWN               = 3 ! No info

    type :: sparse_matrix_t
    private
        class(base_sparse_matrix_t), allocatable :: State
    contains
        procedure         ::                              sparse_matrix_create_square
        procedure         ::                              sparse_matrix_create_rectangular
        procedure         ::                              sparse_matrix_append_entries
        procedure         ::                              sparse_matrix_append_values
        procedure         ::                              sparse_matrix_convert
        procedure         ::                              sparse_matrix_convert_string
        procedure         ::                              sparse_matrix_convert_sparse_matrix_mold
        procedure         ::                              sparse_matrix_convert_base_sparse_matrix_mold
        procedure, public :: get_nnz                   => sparse_matrix_get_nnz
        procedure, public :: get_sign                  => sparse_matrix_get_sign
        procedure, public :: get_num_rows              => sparse_matrix_get_num_rows
        procedure, public :: get_num_cols              => sparse_matrix_get_num_cols
        procedure, public :: get_symmetric_storage     => sparse_matrix_get_symmetric_storage
        procedure, public :: is_by_rows                => sparse_matrix_is_by_rows
        procedure, public :: is_by_cols                => sparse_matrix_is_by_cols
        procedure, public :: is_symmetric              => sparse_matrix_is_symmetric
        procedure, public :: get_default_sparse_matrix => sparse_matrix_get_default_sparse_matrix
        generic,   public :: create                    => sparse_matrix_create_square, &
                                                          sparse_matrix_create_rectangular
        generic,   public :: set_values                => sparse_matrix_append_entries, &
                                                          sparse_matrix_append_values
        generic,   public :: convert                   => sparse_matrix_convert,                         &
                                                          sparse_matrix_convert_string,                  &
                                                          sparse_matrix_convert_sparse_matrix_mold,      &
                                                          sparse_matrix_convert_base_sparse_matrix_mold
        procedure, public :: free                      => sparse_matrix_free
        procedure, public :: apply                     => sparse_matrix_apply
        procedure, public :: print                     => sparse_matrix_print
    end type sparse_matrix_t

    class(base_sparse_matrix_t), allocatable, target, save :: default_sparse_matrix

public :: sparse_matrix_t

contains

    function sparse_matrix_is_symmetric(this) result(is_symmetric)
    !-----------------------------------------------------------------
    !< Get the symmetry property of the matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        logical                               :: is_symmetric
    !-----------------------------------------------------------------
        is_symmetric = this%State%is_symmetric()
    end function sparse_matrix_is_symmetric


    function sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
    !-----------------------------------------------------------------
    !< Get the symmetry storage property of the concrete matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        logical                               :: symmetric_storage
    !-----------------------------------------------------------------
        symmetric_storage = this%State%get_symmetric_storage()
    end function sparse_matrix_get_symmetric_storage


    function sparse_matrix_is_by_rows(this) result(is_by_rows)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        logical                               :: is_by_rows
    !-----------------------------------------------------------------
        is_by_rows = this%State%is_by_rows()
    end function sparse_matrix_is_by_rows


    function sparse_matrix_is_by_cols(this) result(is_by_cols)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by cols
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        logical                               :: is_by_cols
    !-----------------------------------------------------------------
        is_by_cols = this%State%is_by_cols()
    end function sparse_matrix_is_by_cols


    function sparse_matrix_get_num_rows(this) result( num_rows)
    !-----------------------------------------------------------------
    !< Get the number of rows
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip)                           :: num_rows
    !-----------------------------------------------------------------
        num_rows = this%State%get_num_rows()
    end function sparse_matrix_get_num_rows


    function sparse_matrix_get_num_cols(this) result( num_cols)
    !-----------------------------------------------------------------
    !< Get the number of columns
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip)                           :: num_cols
    !-----------------------------------------------------------------
        num_cols = this%State%get_num_cols()
    end function sparse_matrix_get_num_cols


    function sparse_matrix_get_nnz(this) result(nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip)                           :: nnz
    !-----------------------------------------------------------------
        nnz = this%State%get_nnz()
    end function sparse_matrix_get_nnz


    function sparse_matrix_get_sign(this) result( sign)
    !-----------------------------------------------------------------
    !< Get the sign of the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip)                           :: sign
    !-----------------------------------------------------------------
        sign = this%State%get_sign()
    end function sparse_matrix_get_sign


    subroutine sparse_matrix_create_square(this, num_rows_and_cols, symmetric_storage, is_symmetric, sign)
    !-----------------------------------------------------------------
    !< Set the properties and size of a square matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows_and_cols
        logical,                intent(in)    :: symmetric_storage
        logical,                intent(in)    :: is_symmetric
        integer(ip),            intent(in)    :: sign
    !-----------------------------------------------------------------
        if(.not. allocated(this%State)) allocate(coo_sparse_matrix_t :: this%State)
        call this%State%create(num_rows_and_cols, symmetric_storage, is_symmetric, sign)
    end subroutine sparse_matrix_create_square
  

    subroutine sparse_matrix_create_rectangular(this, num_rows, num_cols)
    !-----------------------------------------------------------------
    !< Set the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: num_cols
    !-----------------------------------------------------------------
        if(.not. allocated(this%State)) allocate(coo_sparse_matrix_t :: this%State)
        call this%State%create(num_rows, num_cols)
    end subroutine sparse_matrix_create_rectangular


    subroutine sparse_matrix_append_values(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        call this%State%set_values(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_append_values


    subroutine sparse_matrix_append_entries(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        call this%State%set_values(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_append_entries


    subroutine sparse_matrix_convert(this)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to the default concrete implementation
    !-----------------------------------------------------------------
        class(sparse_matrix_t),    intent(inout) :: this
        class(base_sparse_matrix_t), allocatable :: tmp
        integer                                  :: error
    !-----------------------------------------------------------------
        call this%State%set_state_assembled()
        allocate(tmp, mold=this%get_default_sparse_matrix(), stat=error)
        check(error==0)
        call tmp%move_from_fmt(from=this%State)
        if(allocated(this%State)) deallocate(this%State)
        call move_alloc(from=tmp, to=this%State)
    end subroutine sparse_matrix_convert


    subroutine sparse_matrix_convert_string(this, string)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< The new sparse matrix format is specified using a character array
    !< Valid format strings are 'CSR', 'csr', 'COO' and 'coo'
    !-----------------------------------------------------------------
        class(sparse_matrix_t),    intent(inout) :: this
        character(len=*),          intent(in)    :: string
        class(base_sparse_matrix_t), allocatable :: tmp
        integer                                  :: error
    !-----------------------------------------------------------------
        call this%State%set_state_assembled()
        error = 0
        select case (string)
            case ('CSR', 'csr')
                allocate(csr_sparse_matrix_t :: tmp, stat=error)
            case ('COO', 'coo')
                allocate(coo_sparse_matrix_t :: tmp, stat=error) 
            case default
                check(.false.)
        end select
        check(error==0)
        call tmp%move_from_fmt(from=this%State)
        if(allocated(this%State)) deallocate(this%State)
        call move_alloc(from=tmp, to=this%State)
    end subroutine sparse_matrix_convert_string


    subroutine sparse_matrix_convert_sparse_matrix_mold(this, mold)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< given by a mold
    !-----------------------------------------------------------------
        class(sparse_matrix_t),    intent(inout) :: this
        class(sparse_matrix_t),    intent(in)    :: mold
        class(base_sparse_matrix_t), allocatable :: tmp
        integer                                  :: error
    !-----------------------------------------------------------------
        call this%State%set_state_assembled()
        allocate(tmp, mold=mold%State, stat=error)
        check(error==0)
        call tmp%move_from_fmt(from=this%State)
        if(allocated(this%State)) deallocate(this%State)
        call move_alloc(from=tmp, to=this%State)
    end subroutine sparse_matrix_convert_sparse_matrix_mold


    subroutine sparse_matrix_convert_base_sparse_matrix_mold(this, mold)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< given by a mold
    !-----------------------------------------------------------------
        class(sparse_matrix_t),      intent(inout) :: this
        class(base_sparse_matrix_t), intent(in)    :: mold
        class(base_sparse_matrix_t), allocatable   :: tmp
        integer                                    :: error
    !-----------------------------------------------------------------
        call this%State%set_state_assembled()
        allocate(tmp, mold=mold, stat=error)
        check(error==0)
        call tmp%move_from_fmt(from=this%State)
        if(allocated(this%State)) deallocate(this%State)
        call move_alloc(from=tmp, to=this%State)
    end subroutine sparse_matrix_convert_base_sparse_matrix_mold


    subroutine sparse_matrix_set_default_sparse_matrix(this, mold) 
    !-----------------------------------------------------------------
    !< Allocate the default sparse matrix to a given mold
    !-----------------------------------------------------------------
        class(sparse_matrix_t),      intent(in) :: this
        class(base_sparse_matrix_t), intent(in) :: mold
    !-----------------------------------------------------------------
        if (allocated(default_sparse_matrix)) deallocate(default_sparse_matrix)
        allocate(default_sparse_matrix, mold=mold)
    end subroutine sparse_matrix_set_default_sparse_matrix


    function sparse_matrix_get_default_sparse_matrix(this) result(default_sparse_matrix_pointer)
    !-----------------------------------------------------------------
    !< Get a pointer to the default sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t),   intent(in) :: this
        class(base_sparse_matrix_t), pointer :: default_sparse_matrix_pointer
    !-----------------------------------------------------------------
        if (.not.allocated(default_sparse_matrix)) then 
            allocate(csr_sparse_matrix_t :: default_sparse_matrix)
        end if
        default_sparse_matrix_pointer => default_sparse_matrix
    end function sparse_matrix_get_default_sparse_matrix


    subroutine sparse_matrix_apply(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in)    :: op
        class(vector_t),        intent(in)    :: x
        class(vector_t),        intent(inout) :: y 
    !-----------------------------------------------------------------
        call op%State%apply(x,y)
    end subroutine sparse_matrix_apply


    subroutine sparse_matrix_free(this)
    !-----------------------------------------------------------------
    !< Clean the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(allocated(default_sparse_matrix)) then
            call default_sparse_matrix%free()
            deallocate(default_sparse_matrix)
        endif
        if(allocated(this%State)) then
            call this%State%free()
            deallocate(this%State)
        endif
    end subroutine sparse_matrix_free


    subroutine sparse_matrix_print(this,lunou, only_graph)
    !-----------------------------------------------------------------
    !< Print a Sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t),  intent(in) :: this
        integer(ip),             intent(in) :: lunou
        logical,     optional,   intent(in) :: only_graph
    !-----------------------------------------------------------------
        if(present(only_graph)) then
            call this%State%print(lunou, only_graph=only_graph)
        else
            call this%State%print(lunou)
        endif
    end subroutine sparse_matrix_print

end module sparse_matrix_names
