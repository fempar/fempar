module sparse_matrix_names

USE types_names
USE memor_names
USE vector_names
USE matrix_names
USE vector_space_names
USE serial_scalar_array_names
USE base_sparse_matrix_names, only: base_sparse_matrix_t, coo_sparse_matrix_t, base_sparse_matrix_iterator_t
USE csr_sparse_matrix_names

implicit none

# include "debug.i90"

private

    ! Matrix sign
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE     = 0
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE = 1
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_INDEFINITE            = 2 ! Both positive and negative eigenvalues
    integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_UNKNOWN               = 3 ! No info

    type, extends(matrix_t) :: sparse_matrix_t
    private
        class(base_sparse_matrix_t), allocatable :: State
    contains
        procedure, non_overridable ::                                           sparse_matrix_create_square
        procedure, non_overridable ::                                           sparse_matrix_create_rectangular
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_coords
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_values
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_coords_by_row
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_coords_by_col
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_values_by_row
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_values_by_col
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_single_coord
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_single_value
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_dense_values
        procedure, non_overridable ::                                           sparse_matrix_insert_bounded_square_dense_values
        procedure, non_overridable ::                                           sparse_matrix_insert_dense_values
        procedure, non_overridable ::                                           sparse_matrix_insert_square_dense_values
        procedure, non_overridable ::                                           sparse_matrix_insert_coords
        procedure, non_overridable ::                                           sparse_matrix_insert_values
        procedure, non_overridable ::                                           sparse_matrix_insert_coords_by_row
        procedure, non_overridable ::                                           sparse_matrix_insert_coords_by_col
        procedure, non_overridable ::                                           sparse_matrix_insert_values_by_row
        procedure, non_overridable ::                                           sparse_matrix_insert_values_by_col
        procedure, non_overridable ::                                           sparse_matrix_insert_single_coord
        procedure, non_overridable ::                                           sparse_matrix_insert_single_value
        procedure, non_overridable ::                                           sparse_matrix_convert
        procedure, non_overridable ::                                           sparse_matrix_convert_string
        procedure, non_overridable ::                                           sparse_matrix_convert_sparse_matrix_mold
        procedure, non_overridable ::                                           sparse_matrix_convert_base_sparse_matrix_mold
        procedure, non_overridable ::         create_vector_spaces           => sparse_matrix_create_vector_spaces
        procedure, non_overridable, public :: get_pointer_to_base_matrix     => sparse_matrix_get_pointer_to_base_matrix
        procedure, non_overridable, public :: get_nnz                        => sparse_matrix_get_nnz
        procedure, non_overridable, public :: get_sign                       => sparse_matrix_get_sign
        procedure, non_overridable, public :: get_num_rows                   => sparse_matrix_get_num_rows
        procedure, non_overridable, public :: get_num_cols                   => sparse_matrix_get_num_cols
        procedure, non_overridable, public :: get_symmetric_storage          => sparse_matrix_get_symmetric_storage 
        procedure, non_overridable, public :: is_by_rows                     => sparse_matrix_is_by_rows
        procedure, non_overridable, public :: is_by_cols                     => sparse_matrix_is_by_cols
        procedure, non_overridable, public :: is_symmetric                   => sparse_matrix_is_symmetric
        procedure, non_overridable, public :: get_default_sparse_matrix      => sparse_matrix_get_default_sparse_matrix
        procedure,                  public :: allocate                       => sparse_matrix_allocate
        procedure,                  public :: free_in_stages                 => sparse_matrix_free_in_stages  
        generic,                    public :: create                         => sparse_matrix_create_square, &
                                                                                sparse_matrix_create_rectangular
        generic,                    public :: insert                         => sparse_matrix_insert_bounded_coords,              &
                                                                                sparse_matrix_insert_bounded_values,              &
                                                                                sparse_matrix_insert_bounded_coords_by_row,       &
                                                                                sparse_matrix_insert_bounded_coords_by_col,       &
                                                                                sparse_matrix_insert_bounded_values_by_row,       &
                                                                                sparse_matrix_insert_bounded_values_by_col,       &
                                                                                sparse_matrix_insert_bounded_single_coord,        &
                                                                                sparse_matrix_insert_bounded_single_value,        &
                                                                                sparse_matrix_insert_bounded_dense_values,        &
                                                                                sparse_matrix_insert_bounded_square_dense_values, &
                                                                                sparse_matrix_insert_coords,                      &
                                                                                sparse_matrix_insert_values,                      &
                                                                                sparse_matrix_insert_dense_values,                &
                                                                                sparse_matrix_insert_square_dense_values,         &
                                                                                sparse_matrix_insert_coords_by_row,               & 
                                                                                sparse_matrix_insert_coords_by_col,               &
                                                                                sparse_matrix_insert_values_by_row,               &
                                                                                sparse_matrix_insert_values_by_col,               &
                                                                                sparse_matrix_insert_single_coord,                &
                                                                                sparse_matrix_insert_single_value
        generic,                     public :: convert                       => sparse_matrix_convert,                         &
                                                                                sparse_matrix_convert_string,                  &
                                                                                sparse_matrix_convert_sparse_matrix_mold,      &
                                                                                sparse_matrix_convert_base_sparse_matrix_mold
        procedure,                  public :: split_2x2_numeric              => sparse_matrix_split_2x2_numeric
        procedure,                  public :: split_2x2_symbolic             => sparse_matrix_split_2x2_symbolic
        procedure,                  public :: permute_and_split_2x2_numeric  => sparse_matrix_permute_and_split_2x2_numeric
        procedure,                  public :: permute_and_split_2x2_symbolic => sparse_matrix_permute_and_split_2x2_symbolic
        procedure,                  public :: expand_matrix_numeric          => sparse_matrix_expand_matrix_numeric
        procedure,                  public :: expand_matrix_symbolic         => sparse_matrix_expand_matrix_symbolic
        procedure,                  public :: extract_diagonal               => sparse_matrix_extract_diagonal
        procedure,                  public :: free                           => sparse_matrix_free
        procedure,                  public :: apply                          => sparse_matrix_apply
        procedure, non_overridable, public :: print                          => sparse_matrix_print
        procedure, non_overridable, public :: print_matrix_market            => sparse_matrix_print_matrix_market
        procedure                 , public :: create_iterator                => sparse_matrix_create_iterator
        procedure, non_overridable, public :: get_entry                      => sparse_matrix_get_entry
        procedure, non_overridable, public :: set_value                      => sparse_matrix_set_value
        procedure, non_overridable, public :: add_to_entry                   => sparse_matrix_add_to_entry
    end type sparse_matrix_t

    class(base_sparse_matrix_t), allocatable, target, save :: default_sparse_matrix

    type, extends(matrix_iterator_t) :: sparse_matrix_iterator_t
       private
       class(base_sparse_matrix_iterator_t), allocatable :: base_iterator
     contains
       procedure, non_overridable :: free         => sparse_matrix_iterator_free
       procedure, non_overridable :: next         => sparse_matrix_iterator_next
       procedure, non_overridable :: has_finished => sparse_matrix_iterator_has_finished
       procedure, non_overridable :: get_row      => sparse_matrix_iterator_get_row
       procedure, non_overridable :: get_column   => sparse_matrix_iterator_get_column
       procedure, non_overridable :: get_entry    => sparse_matrix_iterator_get_entry
       procedure, non_overridable :: set_value    => sparse_matrix_iterator_set_value
    end type sparse_matrix_iterator_t

public :: sparse_matrix_t, sparse_matrix_iterator_t

contains

    function sparse_matrix_get_pointer_to_base_matrix(this) result(pointer_to_base_matrix)
    !-----------------------------------------------------------------
    !< Get the symmetry property of the matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t),      target, intent(in) :: this
        class(base_sparse_matrix_t), pointer            :: pointer_to_base_matrix
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        pointer_to_base_matrix => this%State
    end function sparse_matrix_get_pointer_to_base_matrix


    function sparse_matrix_is_symmetric(this) result(is_symmetric)
    !-----------------------------------------------------------------
    !< Get the symmetry property of the matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        logical                            :: is_symmetric
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        is_symmetric = this%State%is_symmetric()
    end function sparse_matrix_is_symmetric


    function sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
    !-----------------------------------------------------------------
    !< Get the symmetry storage property of the concrete matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        logical                            :: symmetric_storage
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        symmetric_storage = this%State%get_symmetric_storage()
    end function sparse_matrix_get_symmetric_storage


    function sparse_matrix_is_by_rows(this) result(is_by_rows)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        logical                            :: is_by_rows
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        is_by_rows = this%State%is_by_rows()
    end function sparse_matrix_is_by_rows


    function sparse_matrix_is_by_cols(this) result(is_by_cols)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by cols
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        logical                            :: is_by_cols
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        is_by_cols = this%State%is_by_cols()
    end function sparse_matrix_is_by_cols


    function sparse_matrix_get_num_rows(this) result( num_rows)
    !-----------------------------------------------------------------
    !< Get the number of rows
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: num_rows
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        num_rows = this%State%get_num_rows()
    end function sparse_matrix_get_num_rows


    function sparse_matrix_get_num_cols(this) result( num_cols)
    !-----------------------------------------------------------------
    !< Get the number of columns
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: num_cols
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        num_cols = this%State%get_num_cols()
    end function sparse_matrix_get_num_cols


    function sparse_matrix_get_nnz(this) result(nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: nnz
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        nnz = this%State%get_nnz()
    end function sparse_matrix_get_nnz


    function sparse_matrix_get_sign(this) result( sign)
    !-----------------------------------------------------------------
    !< Get the sign of the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: sign
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        sign = this%State%get_sign()
    end function sparse_matrix_get_sign


    subroutine sparse_matrix_allocate(this)
    !-----------------------------------------------------------------
    !< Allocate matrix values only if is in a assembled symbolic stage
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%allocate_values()
    end subroutine sparse_matrix_allocate


    subroutine sparse_matrix_create_vector_spaces(this)
    !-----------------------------------------------------------------
    !< Create vector spaces
    !-----------------------------------------------------------------
        class(sparse_matrix_t),      intent(inout) :: this
        type(serial_scalar_array_t)                :: range_vector
        type(serial_scalar_array_t)                :: domain_vector
        type(vector_space_t), pointer              :: range_vector_space
        type(vector_space_t), pointer              :: domain_vector_space
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call range_vector%create(this%get_num_rows())
        call domain_vector%create(this%get_num_cols())
        range_vector_space => this%get_range_vector_space()
        call range_vector_space%create(range_vector)
        domain_vector_space => this%get_domain_vector_space()
        call domain_vector_space%create(domain_vector)
        call range_vector%free()
        call domain_vector%free()
    end subroutine sparse_matrix_create_vector_spaces


    subroutine sparse_matrix_create_square(this, num_rows_and_cols, symmetric_storage, is_symmetric, sign, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a square matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows_and_cols
        logical,                intent(in)    :: symmetric_storage
        logical,                intent(in)    :: is_symmetric
        integer(ip),            intent(in)    :: sign
        integer(ip), optional,  intent(in)    :: nz
    !-----------------------------------------------------------------
        if(.not. allocated(this%State)) allocate(coo_sparse_matrix_t :: this%State)
        call this%State%create(num_rows_and_cols, symmetric_storage, is_symmetric, sign, nz)
        call this%create_vector_spaces()
    end subroutine sparse_matrix_create_square
  

    subroutine sparse_matrix_create_rectangular(this, num_rows, num_cols, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: num_cols
        integer(ip), optional,  intent(in)    :: nz
    !-----------------------------------------------------------------
        if(.not. allocated(this%State)) allocate(coo_sparse_matrix_t :: this%State)
        call this%State%create(num_rows, num_cols, nz)
        call this%create_vector_spaces()
    end subroutine sparse_matrix_create_rectangular


    subroutine sparse_matrix_insert_bounded_values(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
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
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_values


    subroutine sparse_matrix_insert_bounded_coords(this, nz, ia, ja, imin, imax, jmin, jmax) 
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
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_coords


    subroutine sparse_matrix_insert_bounded_values_by_row(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_values_by_row


    subroutine sparse_matrix_insert_bounded_values_by_col(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_values_by_col


    subroutine sparse_matrix_insert_bounded_coords_by_row(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_coords_by_row


    subroutine sparse_matrix_insert_bounded_coords_by_col(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_coords_by_col


    subroutine sparse_matrix_insert_bounded_single_value(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entry and value to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(ia, ja, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_single_value


    subroutine sparse_matrix_insert_bounded_single_coord(this, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entry to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(ia, ja, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_single_coord


    subroutine sparse_matrix_insert_bounded_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: num_cols
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_cols)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_dense_values


    subroutine sparse_matrix_insert_bounded_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_rows)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
    end subroutine sparse_matrix_insert_bounded_square_dense_values


    subroutine sparse_matrix_insert_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: num_cols
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_cols)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(num_rows, num_cols, ia, ja, ioffset, joffset, val)
    end subroutine sparse_matrix_insert_dense_values


    subroutine sparse_matrix_insert_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_rows)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(num_rows, ia, ja, ioffset, joffset, val)
    end subroutine sparse_matrix_insert_square_dense_values


    subroutine sparse_matrix_insert_values(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val)
    end subroutine sparse_matrix_insert_values


    subroutine sparse_matrix_insert_coords(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja)
    end subroutine sparse_matrix_insert_coords


    subroutine sparse_matrix_insert_values_by_row(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val)
    end subroutine sparse_matrix_insert_values_by_row


    subroutine sparse_matrix_insert_values_by_col(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja, val)
    end subroutine sparse_matrix_insert_values_by_col


    subroutine sparse_matrix_insert_coords_by_row(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja)
    end subroutine sparse_matrix_insert_coords_by_row


    subroutine sparse_matrix_insert_coords_by_col(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(nz, ia, ja)
    end subroutine sparse_matrix_insert_coords_by_col


    subroutine sparse_matrix_insert_single_value(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entry and value to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(ia, ja, val)
    end subroutine sparse_matrix_insert_single_value


    subroutine sparse_matrix_insert_single_coord(this, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entry to the sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%insert(ia, ja)
    end subroutine sparse_matrix_insert_single_coord


    subroutine sparse_matrix_convert(this)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to the default concrete implementation
    !-----------------------------------------------------------------
        class(sparse_matrix_t),    intent(inout) :: this
        class(base_sparse_matrix_t), allocatable :: tmp
        integer                                  :: error
    !-----------------------------------------------------------------
        assert(allocated(this%State))
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
        assert(allocated(this%State))
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
        type(sparse_matrix_t) ,    intent(in)    :: mold
        class(base_sparse_matrix_t), allocatable :: tmp
        integer                                  :: error
    !-----------------------------------------------------------------
        assert(allocated(this%State))
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
        assert(allocated(this%State))
        allocate(tmp, mold=mold, stat=error)
        check(error==0)
        call this%State%convert_body()
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


    subroutine sparse_matrix_split_2x2_numeric(this, num_row, num_col, A_II, A_IG, A_GI, A_GG, symmetric_storage, symmetric, sign)
    !-----------------------------------------------------------------
    !< Split matrix in a 2x2 submatrix
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !-----------------------------------------------------------------
        class(sparse_matrix_t),          intent(in)    :: this
        integer,                         intent(in)    :: num_row
        integer,                         intent(in)    :: num_col
        type(sparse_matrix_t),           intent(inout) :: A_II
        type(sparse_matrix_t),           intent(inout) :: A_IG
        type(sparse_matrix_t), optional, intent(inout) :: A_GI
        type(sparse_matrix_t),           intent(inout) :: A_GG
        logical,               optional, intent(in)    :: symmetric_storage
        logical,               optional, intent(in)    :: symmetric
        integer(ip),           optional, intent(in)    :: sign
        logical                                        :: symmetric_storage_aux
        logical                                        :: symmetric_aux
        integer(ip)                                    :: sign_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        assert (.not. present(A_GI) .or. (.not. this%get_symmetric_storage()) )

        if(present(A_GI)) then
            if(.not. allocated(A_GI%State)) allocate(A_GI%State, mold=this%State)
        endif

        if(.not. allocated(A_II%State)) allocate(A_II%State, mold=this%State)
        if(.not. allocated(A_IG%State)) allocate(A_IG%State, mold=this%State)
        if(.not. allocated(A_GG%State)) allocate(A_GG%State, mold=this%State)

        sign_aux = this%State%get_sign()
        symmetric_aux = this%State%is_symmetric()
        symmetric_storage_aux = this%State%get_symmetric_storage()

        if(present(sign)) sign_aux = sign
        if(present(symmetric)) symmetric_aux = symmetric
        if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

        if(A_II%State%state_is_start()) call A_II%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)
        if(A_IG%State%state_is_start()) call A_IG%State%set_properties(.false., .false., SPARSE_MATRIX_SIGN_UNKNOWN)
        if(A_GG%State%state_is_start()) call A_GG%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)

        if(present(A_GI)) then
            if(A_GI%State%state_is_start()) call A_GI%State%set_properties(.false., .false., SPARSE_MATRIX_SIGN_UNKNOWN)
            call this%State%split_2x2_numeric(num_row, num_col, A_II%State, A_IG%State, A_GI%State, A_GG%State)
            call A_GI%create_vector_spaces()
        else
            call this%State%split_2x2_numeric(num_row, num_col, A_II=A_II%State, A_IG=A_IG%State, A_GG=A_GG%State)
        endif
        call A_II%create_vector_spaces()
        call A_IG%create_vector_spaces()
        call A_GG%create_vector_spaces()
    end subroutine sparse_matrix_split_2x2_numeric


    subroutine sparse_matrix_split_2x2_symbolic(this, num_row, num_col, A_II, A_IG, A_GI, A_GG, symmetric_storage, symmetric, sign)
    !-----------------------------------------------------------------
    !< Split matrix in a 2x2 submatrix
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !-----------------------------------------------------------------
        class(sparse_matrix_t),          intent(in)    :: this
        integer,                         intent(in)    :: num_row
        integer,                         intent(in)    :: num_col
        type(sparse_matrix_t),           intent(inout) :: A_II
        type(sparse_matrix_t),           intent(inout) :: A_IG
        type(sparse_matrix_t), optional, intent(inout) :: A_GI
        type(sparse_matrix_t),           intent(inout) :: A_GG
        integer,               optional, intent(in)    :: sign
        logical,               optional, intent(in)    :: symmetric
        logical,               optional, intent(in)    :: symmetric_storage
        integer                                        :: sign_aux
        logical                                        :: symmetric_aux
        logical                                        :: symmetric_storage_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        assert (.not. present(A_GI) .or. (.not. this%get_symmetric_storage()) )
        assert(.not. allocated(A_II%State))
        assert(.not. allocated(A_IG%State))
        assert(.not. allocated(A_GG%State))

        if(present(A_GI)) then
            assert(.not. allocated(A_GI%State))
            allocate(A_GI%State, mold=this%State)
        endif

        allocate(A_II%State, mold=this%State)
        allocate(A_IG%State, mold=this%State)
        allocate(A_GG%State, mold=this%State)

        sign_aux = this%State%get_sign()
        symmetric_aux = this%State%is_symmetric()
        symmetric_storage_aux = this%State%get_symmetric_storage()

        if(present(sign)) sign_aux = sign
        if(present(symmetric)) symmetric_aux = symmetric
        if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

        call A_II%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)
        call A_IG%State%set_properties(.false., .false., SPARSE_MATRIX_SIGN_UNKNOWN)
        call A_GG%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)

        if(present(A_GI)) then
            call A_GI%State%set_properties(.false., .false., SPARSE_MATRIX_SIGN_UNKNOWN)
            call this%State%split_2x2_symbolic(num_row, num_col, A_II%State, A_IG%State, A_GI%State, A_GG%State)
            call A_GI%create_vector_spaces()
        else
            call this%State%split_2x2_symbolic(num_row, num_col, A_II=A_II%State, A_IG=A_IG%State, A_GG=A_GG%State)
        endif
        call A_II%create_vector_spaces()
        call A_IG%create_vector_spaces()
        call A_GG%create_vector_spaces()
    end subroutine sparse_matrix_split_2x2_symbolic


    subroutine sparse_matrix_permute_and_split_2x2_numeric(this, num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR, symmetric_storage, symmetric, sign) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows 
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_CC, A_RC, A_CR and A_RR given the global 
    !< matrix A. Note that A_CC, A_RC, A_CR are dense and A_RR is sparse
    !-----------------------------------------------------------------
        class(sparse_matrix_t),           intent(in)    :: this
        integer(ip),                      intent(in)    :: num_row
        integer(ip),                      intent(in)    :: num_col
        integer(ip),                      intent(in)    :: perm(:)
        integer(ip),                      intent(in)    :: iperm(:)
        real(rp),    allocatable,         intent(inout) :: A_CC(:,:)
        real(rp),    allocatable,         intent(inout) :: A_CR(:,:)
        real(rp),    allocatable,         intent(inout) :: A_RC(:,:)
        class(sparse_matrix_t),           intent(inout) :: A_RR
        logical,                optional, intent(in)    :: symmetric_storage
        logical,                optional, intent(in)    :: symmetric
        integer(ip),            optional, intent(in)    :: sign
        logical                                         :: symmetric_storage_aux
        logical                                         :: symmetric_aux
        integer(ip)                                     :: sign_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))

        if(allocated(A_CC)) call memfree(A_CC, __FILE__, __LINE__)
        if(allocated(A_CR)) call memfree(A_CR, __FILE__, __LINE__)
        if(allocated(A_RC)) call memfree(A_RC, __FILE__, __LINE__)

        if(.not. allocated(A_RR%State)) then
            allocate(A_RR%State, mold=this%State)

            sign_aux = this%State%get_sign()
            symmetric_aux = this%State%is_symmetric()
            symmetric_storage_aux = this%State%get_symmetric_storage()

            if(present(sign)) sign_aux = sign
            if(present(symmetric)) symmetric_aux = symmetric
            if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

            call A_RR%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)
        endif

        call this%State%permute_and_split_2x2_numeric(num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR%State)
        call A_RR%create_vector_spaces()
    end subroutine sparse_matrix_permute_and_split_2x2_numeric


    subroutine sparse_matrix_permute_and_split_2x2_symbolic(this, num_row, num_col, perm, iperm, A_RR, symmetric_storage, symmetric, sign) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !<
    !< this routine computes A_RR from the global matrix A
    !< A_CC, ACR and A_RC sparsity pattern calculation is not
    !< performed because they are dense matrices
    !-----------------------------------------------------------------
        class(sparse_matrix_t),           intent(in)    :: this
        integer(ip),                      intent(in)    :: num_row
        integer(ip),                      intent(in)    :: num_col
        integer(ip),                      intent(in)    :: perm(:)
        integer(ip),                      intent(in)    :: iperm(:)
        class(sparse_matrix_t),           intent(inout) :: A_RR
        logical,                optional, intent(in)    :: symmetric_storage
        logical,                optional, intent(in)    :: symmetric
        integer(ip),            optional, intent(in)    :: sign
        logical                                         :: symmetric_storage_aux
        logical                                         :: symmetric_aux
        integer(ip)                                     :: sign_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        assert(.not. allocated(A_RR%State))

        allocate(A_RR%State, mold=this%State)

        sign_aux = this%State%get_sign()
        symmetric_aux = this%State%is_symmetric()
        symmetric_storage_aux = this%State%get_symmetric_storage()

        if(present(sign)) sign_aux = sign
        if(present(symmetric)) symmetric_aux = symmetric
        if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

        call A_RR%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)

        call this%State%permute_and_split_2x2_symbolic(num_row, num_col, perm, iperm, A_RR%State)
        call A_RR%create_vector_spaces()
    end subroutine sparse_matrix_permute_and_split_2x2_symbolic


    subroutine sparse_matrix_expand_matrix_numeric(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to, symmetric_storage, symmetric, sign)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(sparse_matrix_t),          intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        real(rp),                        intent(in)    :: C_T_val(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        real(rp),                        intent(in)    :: I_val(C_T_nz)
        class(sparse_matrix_t),          intent(inout) :: to
        logical,               optional, intent(in)    :: symmetric_storage
        logical,               optional, intent(in)    :: symmetric
        integer(ip),           optional, intent(in)    :: sign
        logical                                        :: symmetric_storage_aux
        logical                                        :: symmetric_aux
        integer(ip)                                    :: sign_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        if(.not. allocated(to%State)) allocate(to%State, mold=this%State)

        if(to%State%state_is_start()) then
            sign_aux = this%State%get_sign()
            symmetric_aux = this%State%is_symmetric()
            symmetric_storage_aux = this%State%get_symmetric_storage()

            if(present(sign)) sign_aux = sign
            if(present(symmetric)) symmetric_aux = symmetric
            if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

            call to%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)
        endif

        call this%State%expand_matrix_numeric(C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to%State)
        call to%create_vector_spaces()
    end subroutine sparse_matrix_expand_matrix_numeric


    subroutine sparse_matrix_expand_matrix_symbolic(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to, symmetric_storage, symmetric, sign)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !-----------------------------------------------------------------
        class(sparse_matrix_t),          intent(inout) :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        class(sparse_matrix_t),          intent(inout) :: to
        logical,               optional, intent(in)    :: symmetric_storage
        logical,               optional, intent(in)    :: symmetric
        integer(ip),           optional, intent(in)    :: sign
        logical                                        :: symmetric_storage_aux
        logical                                        :: symmetric_aux
        integer(ip)                                    :: sign_aux
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        if(.not. allocated(to%State)) allocate(to%State, mold=this%State)
        assert(to%State%state_is_start())

        sign_aux = this%State%get_sign()
        symmetric_aux = this%State%is_symmetric()
        symmetric_storage_aux = this%State%get_symmetric_storage()

        if(present(sign)) sign_aux = sign
        if(present(symmetric)) symmetric_aux = symmetric
        if(present(symmetric_storage)) symmetric_storage_aux = symmetric_storage

        call to%State%set_properties(symmetric_storage_aux, symmetric_aux, sign_aux)

        call this%State%expand_matrix_symbolic(C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to%State)
        call to%create_vector_spaces()
    end subroutine sparse_matrix_expand_matrix_symbolic


    subroutine sparse_matrix_extract_diagonal(this, diagonal) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in)    :: this
        real(rp), allocatable,  intent(inout) :: diagonal(:)
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        call this%State%extract_diagonal(diagonal)
    end subroutine sparse_matrix_extract_diagonal


    subroutine sparse_matrix_apply(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in)    :: op
        class(vector_t),        intent(in)    :: x
        class(vector_t),        intent(inout) :: y 
    !-----------------------------------------------------------------
        assert(allocated(op%State))
        call op%abort_if_not_in_domain(x)
        call op%abort_if_not_in_range(y)
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
        call this%free_vector_spaces()
    end subroutine sparse_matrix_free


    subroutine sparse_matrix_free_in_stages(this, action)
    !-----------------------------------------------------------------
    !< free_in_stages procedure.
    !< As it extends from matrix_t, it must be implemented
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: action
        integer(ip) :: nnz
        integer(ip) :: num_rows                        
        integer(ip) :: num_cols
        integer(ip) :: state
        integer(ip) :: sign
        logical     :: symmetric
        logical     :: symmetric_storage
    !-----------------------------------------------------------------
        if(allocated(this%State)) then
            if(action == free_numerical_setup) then
                call this%State%free_numeric()
            elseif(action == free_symbolic_setup) then
                call this%State%free_symbolic()
            elseif(action == free_clean) then
                call this%State%free_clean()
                call this%free_vector_spaces()
                deallocate(this%State)
                if(allocated(default_sparse_matrix)) then
                    call default_sparse_matrix%free()
                    deallocate(default_sparse_matrix)
                endif
            else
                call this%free()
            endif
        endif
    end subroutine sparse_matrix_free_in_stages


    subroutine sparse_matrix_print(this,lunou, only_graph)
    !-----------------------------------------------------------------
    !< Print a Sparse matrix
    !-----------------------------------------------------------------
        class(sparse_matrix_t),  intent(in) :: this
        integer(ip),             intent(in) :: lunou
        logical,     optional,   intent(in) :: only_graph
    !-----------------------------------------------------------------
        assert(allocated(this%State))
        if(present(only_graph)) then
            call this%State%print(lunou, only_graph=only_graph)
        else
            call this%State%print(lunou)
        endif
    end subroutine sparse_matrix_print

    subroutine sparse_matrix_print_matrix_market (this, lunou, ng, l2g)
    !-----------------------------------------------------------------
    !< Print a Sparse matrix in matrix market format
    !-----------------------------------------------------------------
        class(sparse_matrix_t), intent(in) :: this
        integer(ip),            intent(in) :: lunou
        integer(ip), optional,  intent(in) :: ng
        integer(ip), optional,  intent(in) :: l2g (*)
    !-----------------------------------------------------------------
        call this%State%print_matrix_market(lunou, ng, l2g)
    end subroutine sparse_matrix_print_matrix_market
    
    subroutine sparse_matrix_create_iterator(this, iblock, jblock, iterator)
      !-----------------------------------------------------------------
      !< Get a pointer to an iterator over the matrix entries
      !-----------------------------------------------------------------
      class(sparse_matrix_t)               , intent(in)    :: this
      integer(ip)                          , intent(in)    :: iblock 
      integer(ip)                          , intent(in)    :: jblock 
      class(matrix_iterator_t), allocatable, intent(inout) :: iterator
      !-----------------------------------------------------------------
      assert(iblock == 1)
      assert(jblock == 1)
      assert(allocated(this%State))
      assert(this%State%state_is_assembled())
      if (allocated(iterator)) deallocate(iterator)
      allocate(sparse_matrix_iterator_t :: iterator)
      select type ( iterator)
      class is (sparse_matrix_iterator_t) 
         call this%State%create_iterator(iterator%base_iterator)
      class DEFAULT
         assert(.false.)
      end select
    end subroutine sparse_matrix_create_iterator

    function sparse_matrix_get_entry(this, ia, ja, val)
      !-----------------------------------------------------------------
      !< Get the value in the (ia,ja) entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_t), intent(in)  :: this
      integer(ip)           , intent(in)  :: ia
      integer(ip)           , intent(in)  :: ja
      real(rp)              , intent(out) :: val
      logical                             :: sparse_matrix_get_entry
      !-----------------------------------------------------------------
      assert(allocated(this%State))
      assert(this%State%state_is_assembled())
      sparse_matrix_get_entry = this%State%get_entry(ia, ja, val)
    end function sparse_matrix_get_entry

    function sparse_matrix_set_value(this, ia, ja, val)
      !-----------------------------------------------------------------
      !< Set the value in the (ia,ja) entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_t), intent(inout) :: this
      integer(ip)           , intent(in)    :: ia
      integer(ip)           , intent(in)    :: ja
      real(rp)              , intent(in)    :: val
      logical                               :: sparse_matrix_set_value
      !-----------------------------------------------------------------
      assert(allocated(this%State))
      assert(this%State%state_is_assembled())
      sparse_matrix_set_value = this%State%set_value(ia, ja, val)
    end function sparse_matrix_set_value

    function sparse_matrix_add_to_entry(this, ia, ja, val)
      !-----------------------------------------------------------------
      !< Sum the value in the (ia,ja) entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_t), intent(inout) :: this
      integer(ip)           , intent(in)    :: ia
      integer(ip)           , intent(in)    :: ja
      real(rp)              , intent(in)    :: val
      logical                               :: sparse_matrix_add_to_entry
      !-----------------------------------------------------------------
      assert(allocated(this%State))
      assert(this%State%state_is_assembled())
      sparse_matrix_add_to_entry = this%State%add_to_entry(ia, ja, val)
    end function sparse_matrix_add_to_entry
    !-----------------------------------------------------------------
    !< SPARSE_MATRIX_ITERATOR SUBROUTINES
    !-----------------------------------------------------------------
     subroutine sparse_matrix_iterator_free(this)
      !-----------------------------------------------------------------
      !< Free the information in the iterator
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(inout) :: this
      call this%base_iterator%free()
    end subroutine sparse_matrix_iterator_free

    subroutine sparse_matrix_iterator_next(this)
      !-----------------------------------------------------------------
      !< Set the pointer to the following entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(inout) :: this
      call this%base_iterator%next()
    end subroutine sparse_matrix_iterator_next

    function sparse_matrix_iterator_has_finished(this)
      !-----------------------------------------------------------------
      !< Check if the pointer of the matrix has reached the end
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(in) :: this
      logical :: sparse_matrix_iterator_has_finished
      sparse_matrix_iterator_has_finished = this%base_iterator%has_finished()
    end function sparse_matrix_iterator_has_finished

    function sparse_matrix_iterator_get_row(this)
      !-----------------------------------------------------------------
      !< Get the row index of the entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(in) :: this
      integer(ip) :: sparse_matrix_iterator_get_row
      sparse_matrix_iterator_get_row = this%base_iterator%get_row()
    end function sparse_matrix_iterator_get_row
    
    function sparse_matrix_iterator_get_column(this)
      !-----------------------------------------------------------------
      !<  Get the column index of the entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(in) :: this
      integer(ip) :: sparse_matrix_iterator_get_column
      sparse_matrix_iterator_get_column = this%base_iterator%get_column()
    end function sparse_matrix_iterator_get_column

    function sparse_matrix_iterator_get_entry(this)
      !-----------------------------------------------------------------
      !< Get the value of the entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(in) :: this
      real(rp) :: sparse_matrix_iterator_get_entry
      sparse_matrix_iterator_get_entry = this%base_iterator%get_entry()
    end function sparse_matrix_iterator_get_entry

    subroutine sparse_matrix_iterator_set_value(this,new_value)
      !-----------------------------------------------------------------
      !< Set the value of the entry of the matrix
      !-----------------------------------------------------------------
      class(sparse_matrix_iterator_t), intent(inout) :: this
      real(rp)                       , intent(in)    :: new_value
      call this%base_iterator%set_value(new_value)
    end subroutine sparse_matrix_iterator_set_value

  end module sparse_matrix_names
