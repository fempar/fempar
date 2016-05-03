module par_sparse_matrix_names
USE types_names
USE memor_names
USE vector_names
USE matrix_names
use base_sparse_matrix_names
use sparse_matrix_names
USE vector_space_names
use dof_import_names
use par_environment_names
use serial_scalar_array_names
USE par_scalar_array_names

implicit none

# include "debug.i90"

private

    type, extends(matrix_t) :: par_sparse_matrix_t
      private
      type(sparse_matrix_t)            :: sparse_matrix
      type(par_environment_t), pointer :: p_env             => NULL()
      type(dof_import_t)     , pointer :: dof_import_domain => NULL()
      type(dof_import_t)     , pointer :: dof_import_range  => NULL()
    contains
        procedure, non_overridable ::                                           par_sparse_matrix_create_square
        procedure, non_overridable ::                                           par_sparse_matrix_create_rectangular
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_coords
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_coords_by_row
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_coords_by_col
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_values_by_row
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_values_by_col
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_single_coord
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_single_value
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_dense_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_bounded_square_dense_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_dense_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_square_dense_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_coords
        procedure, non_overridable ::                                           par_sparse_matrix_insert_values
        procedure, non_overridable ::                                           par_sparse_matrix_insert_coords_by_row
        procedure, non_overridable ::                                           par_sparse_matrix_insert_coords_by_col
        procedure, non_overridable ::                                           par_sparse_matrix_insert_values_by_row
        procedure, non_overridable ::                                           par_sparse_matrix_insert_values_by_col
        procedure, non_overridable ::                                           par_sparse_matrix_insert_single_coord
        procedure, non_overridable ::                                           par_sparse_matrix_insert_single_value
        procedure, non_overridable ::                                           par_sparse_matrix_convert
        procedure, non_overridable ::                                           par_sparse_matrix_convert_string
        procedure, non_overridable ::                                           par_sparse_matrix_convert_par_sparse_matrix_mold
        procedure, non_overridable ::                                           par_sparse_matrix_convert_base_par_sparse_matrix_mold
        procedure, non_overridable ::         create_vector_spaces           => par_sparse_matrix_create_vector_spaces
        procedure, non_overridable, public :: get_nnz                        => par_sparse_matrix_get_nnz
        procedure, non_overridable, public :: get_sign                       => par_sparse_matrix_get_sign
        procedure, non_overridable, public :: get_num_rows                   => par_sparse_matrix_get_num_rows
        procedure, non_overridable, public :: get_num_cols                   => par_sparse_matrix_get_num_cols
        procedure, non_overridable, public :: get_symmetric_storage          => par_sparse_matrix_get_symmetric_storage 
        procedure, non_overridable, public :: is_by_rows                     => par_sparse_matrix_is_by_rows
        procedure, non_overridable, public :: is_by_cols                     => par_sparse_matrix_is_by_cols
        procedure, non_overridable, public :: is_symmetric                   => par_sparse_matrix_is_symmetric
        procedure,                  public :: allocate                       => par_sparse_matrix_allocate
        procedure,                  public :: free_in_stages                 => par_sparse_matrix_free_in_stages  
        generic,                    public :: create                         => par_sparse_matrix_create_square, &
                                                                                par_sparse_matrix_create_rectangular
        generic,                    public :: insert                         => par_sparse_matrix_insert_bounded_coords,              &
                                                                                par_sparse_matrix_insert_bounded_values,              &
                                                                                par_sparse_matrix_insert_bounded_coords_by_row,       &
                                                                                par_sparse_matrix_insert_bounded_coords_by_col,       &
                                                                                par_sparse_matrix_insert_bounded_values_by_row,       &
                                                                                par_sparse_matrix_insert_bounded_values_by_col,       &
                                                                                par_sparse_matrix_insert_bounded_single_coord,        &
                                                                                par_sparse_matrix_insert_bounded_single_value,        &
                                                                                par_sparse_matrix_insert_bounded_dense_values,        &
                                                                                par_sparse_matrix_insert_bounded_square_dense_values, &
                                                                                par_sparse_matrix_insert_coords,                      &
                                                                                par_sparse_matrix_insert_values,                      &
                                                                                par_sparse_matrix_insert_dense_values,                &
                                                                                par_sparse_matrix_insert_square_dense_values,         &
                                                                                par_sparse_matrix_insert_coords_by_row,               & 
                                                                                par_sparse_matrix_insert_coords_by_col,               &
                                                                                par_sparse_matrix_insert_values_by_row,               &
                                                                                par_sparse_matrix_insert_values_by_col,               &
                                                                                par_sparse_matrix_insert_single_coord,                &
                                                                                par_sparse_matrix_insert_single_value
        generic,                     public :: convert                       => par_sparse_matrix_convert,                         &
                                                                                par_sparse_matrix_convert_string,                  &
                                                                                par_sparse_matrix_convert_par_sparse_matrix_mold,      &
                                                                                par_sparse_matrix_convert_base_par_sparse_matrix_mold
        procedure,                  public :: extract_diagonal               => par_sparse_matrix_extract_diagonal
        procedure,                  public :: apply                          => par_sparse_matrix_apply
        procedure, non_overridable, public :: print                          => par_sparse_matrix_print
        procedure, non_overridable, public :: print_matrix_market            => par_sparse_matrix_print_matrix_market
        procedure, non_overridable, public :: create_iterator                => par_sparse_matrix_create_iterator
    end type par_sparse_matrix_t

public :: par_sparse_matrix_t

contains

    function par_sparse_matrix_is_symmetric(this) result(is_symmetric)
    !-----------------------------------------------------------------
    !< Get the symmetry property of the matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        logical                            :: is_symmetric
    !-----------------------------------------------------------------
        is_symmetric = this%sparse_matrix%is_symmetric()
    end function par_sparse_matrix_is_symmetric


    function par_sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
    !-----------------------------------------------------------------
    !< Get the symmetry storage property of the concrete matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        logical                            :: symmetric_storage
    !-----------------------------------------------------------------
        symmetric_storage = this%sparse_matrix%get_symmetric_storage()
    end function par_sparse_matrix_get_symmetric_storage


    function par_sparse_matrix_is_by_rows(this) result(is_by_rows)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by rows
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        logical                            :: is_by_rows
    !-----------------------------------------------------------------
        is_by_rows = this%sparse_matrix%is_by_rows()
    end function par_sparse_matrix_is_by_rows


    function par_sparse_matrix_is_by_cols(this) result(is_by_cols)
    !-----------------------------------------------------------------
    !< Check if the matrix is sorted by cols
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        logical                            :: is_by_cols
    !-----------------------------------------------------------------
        is_by_cols = this%sparse_matrix%is_by_cols()
    end function par_sparse_matrix_is_by_cols


    function par_sparse_matrix_get_num_rows(this) result( num_rows)
    !-----------------------------------------------------------------
    !< Get the number of rows
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: num_rows
    !-----------------------------------------------------------------
        num_rows = this%sparse_matrix%get_num_rows()
    end function par_sparse_matrix_get_num_rows


    function par_sparse_matrix_get_num_cols(this) result( num_cols)
    !-----------------------------------------------------------------
    !< Get the number of columns
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: num_cols
    !-----------------------------------------------------------------
        num_cols = this%sparse_matrix%get_num_cols()
    end function par_sparse_matrix_get_num_cols


    function par_sparse_matrix_get_nnz(this) result(nnz)
    !-----------------------------------------------------------------
    !< Get the number of non zeros of the matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: nnz
    !-----------------------------------------------------------------
        nnz = this%sparse_matrix%get_nnz()
    end function par_sparse_matrix_get_nnz


    function par_sparse_matrix_get_sign(this) result( sign)
    !-----------------------------------------------------------------
    !< Get the sign of the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        integer(ip)                        :: sign
    !-----------------------------------------------------------------
        sign = this%sparse_matrix%get_sign()
    end function par_sparse_matrix_get_sign


    subroutine par_sparse_matrix_allocate(this)
    !-----------------------------------------------------------------
    !< Allocate matrix values only if is in a assembled symbolic stage
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%allocate()
    end subroutine par_sparse_matrix_allocate


    subroutine par_sparse_matrix_create_vector_spaces(this)
    !-----------------------------------------------------------------
    !< Create vector spaces
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t),      intent(inout) :: this
        type(par_scalar_array_t)                :: range_vector
        type(par_scalar_array_t)                :: domain_vector
        type(vector_space_t), pointer           :: range_vector_space
        type(vector_space_t), pointer           :: domain_vector_space
    !-----------------------------------------------------------------
        call range_vector%create(this%p_env, this%dof_import_range)
        call domain_vector%create(this%p_env, this%dof_import_domain)
        range_vector_space => this%get_range_vector_space()
        call range_vector_space%create(range_vector)
        domain_vector_space => this%get_domain_vector_space()
        call domain_vector_space%create(domain_vector)
        call range_vector%free()
        call domain_vector%free()
    end subroutine par_sparse_matrix_create_vector_spaces


    subroutine par_sparse_matrix_create_square(this, p_env, dof_import, symmetric_storage, is_symmetric, sign, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a square matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t)     , intent(inout) :: this
        type(par_environment_t), target, intent(in)    :: p_env
        type(dof_import_t)     , target, intent(in)    :: dof_import
        logical                        , intent(in)    :: symmetric_storage
        logical                        , intent(in)    :: is_symmetric
        integer(ip)                    , intent(in)    :: sign
        integer(ip)          , optional,  intent(in)   :: nz
    !-----------------------------------------------------------------
        this%p_env             => p_env
        this%dof_import_domain => dof_import
        this%dof_import_range  => dof_import
        if (this%p_env%p_context%iam>=0) then
          call this%sparse_matrix%create(dof_import%get_number_dofs(), symmetric_storage, is_symmetric, sign, nz)
        end if
        call this%create_vector_spaces()        
    end subroutine par_sparse_matrix_create_square
  

    subroutine par_sparse_matrix_create_rectangular(this, p_env, dof_import_range, dof_import_domain, nz)
    !-----------------------------------------------------------------
    !< Set the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t)     , intent(inout) :: this
        type(par_environment_t), target, intent(in)    :: p_env
        type(dof_import_t)     , target, intent(in)    :: dof_import_range
        type(dof_import_t)     , target, intent(in)    :: dof_import_domain
        integer(ip)          , optional,  intent(in)   :: nz
    !-----------------------------------------------------------------
        this%p_env             => p_env
        this%dof_import_domain => dof_import_domain
        this%dof_import_range  => dof_import_range
        if (this%p_env%p_context%iam>=0) then
          call this%sparse_matrix%create(dof_import_range%get_number_dofs(), dof_import_domain%get_number_dofs(),nz)
        end if
        call this%create_vector_spaces()
    end subroutine par_sparse_matrix_create_rectangular


    subroutine par_sparse_matrix_insert_bounded_values(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_values


    subroutine par_sparse_matrix_insert_bounded_coords(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_coords


    subroutine par_sparse_matrix_insert_bounded_values_by_row(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_values_by_row


    subroutine par_sparse_matrix_insert_bounded_values_by_col(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_values_by_col


    subroutine par_sparse_matrix_insert_bounded_coords_by_row(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_coords_by_row


    subroutine par_sparse_matrix_insert_bounded_coords_by_col(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_coords_by_col


    subroutine par_sparse_matrix_insert_bounded_single_value(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entry and value to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(ia, ja, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_single_value


    subroutine par_sparse_matrix_insert_bounded_single_coord(this, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entry to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        integer(ip),            intent(in)    :: imin
        integer(ip),            intent(in)    :: imax
        integer(ip),            intent(in)    :: jmin
        integer(ip),            intent(in)    :: jmax
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(ia, ja, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_single_coord


    subroutine par_sparse_matrix_insert_bounded_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
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
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_dense_values


    subroutine par_sparse_matrix_insert_bounded_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
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
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
    end subroutine par_sparse_matrix_insert_bounded_square_dense_values


    subroutine par_sparse_matrix_insert_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: num_cols
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_cols)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(num_rows, num_cols, ia, ja, ioffset, joffset, val)
    end subroutine par_sparse_matrix_insert_dense_values


    subroutine par_sparse_matrix_insert_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: num_rows
        integer(ip),            intent(in)    :: ia(num_rows)
        integer(ip),            intent(in)    :: ja(num_rows)
        integer(ip),            intent(in)    :: ioffset
        integer(ip),            intent(in)    :: joffset
        real(rp),               intent(in)    :: val(:,:)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(num_rows, ia, ja, ioffset, joffset, val)
    end subroutine par_sparse_matrix_insert_square_dense_values


    subroutine par_sparse_matrix_insert_values(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val)
    end subroutine par_sparse_matrix_insert_values


    subroutine par_sparse_matrix_insert_coords(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja)
    end subroutine par_sparse_matrix_insert_coords


    subroutine par_sparse_matrix_insert_values_by_row(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val)
    end subroutine par_sparse_matrix_insert_values_by_row


    subroutine par_sparse_matrix_insert_values_by_col(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja, val)
    end subroutine par_sparse_matrix_insert_values_by_col


    subroutine par_sparse_matrix_insert_coords_by_row(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja)
    end subroutine par_sparse_matrix_insert_coords_by_row


    subroutine par_sparse_matrix_insert_coords_by_col(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: nz
        integer(ip),            intent(in)    :: ia(nz)
        integer(ip),            intent(in)    :: ja
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(nz, ia, ja)
    end subroutine par_sparse_matrix_insert_coords_by_col


    subroutine par_sparse_matrix_insert_single_value(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entry and value to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
        real(rp),               intent(in)    :: val
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(ia, ja, val)
    end subroutine par_sparse_matrix_insert_single_value


    subroutine par_sparse_matrix_insert_single_coord(this, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entry to the sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(inout) :: this
        integer(ip),            intent(in)    :: ia
        integer(ip),            intent(in)    :: ja
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%insert(ia, ja)
    end subroutine par_sparse_matrix_insert_single_coord


    subroutine par_sparse_matrix_convert(this)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to the default concrete implementation
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t),    intent(inout) :: this
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%convert()
    end subroutine par_sparse_matrix_convert


    subroutine par_sparse_matrix_convert_string(this, string)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< The new sparse matrix format is specified using a character array
    !< Valid format strings are 'CSR', 'csr', 'COO' and 'coo'
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t),    intent(inout) :: this
        character(len=*),          intent(in)    :: string
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%convert(string)
    end subroutine par_sparse_matrix_convert_string


    subroutine par_sparse_matrix_convert_par_sparse_matrix_mold(this, mold)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< given by a mold
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t),    intent(inout) :: this
        type(par_sparse_matrix_t) ,    intent(in)    :: mold
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%convert(mold%sparse_matrix)
    end subroutine par_sparse_matrix_convert_par_sparse_matrix_mold


    subroutine par_sparse_matrix_convert_base_par_sparse_matrix_mold(this, mold)
    !-----------------------------------------------------------------
    !< Change the state of the matrix to different concrete implementation
    !< given by a mold
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t) , intent(inout) :: this
        class(base_sparse_matrix_t), intent(in)    :: mold
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%convert(mold)
    end subroutine par_sparse_matrix_convert_base_par_sparse_matrix_mold

    subroutine par_sparse_matrix_extract_diagonal(this, diagonal) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in)    :: this
        real(rp), allocatable,  intent(inout) :: diagonal(:)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%extract_diagonal(diagonal)
    end subroutine par_sparse_matrix_extract_diagonal


    subroutine par_sparse_matrix_apply(op,x,y) 
    !-----------------------------------------------------------------
    !< Apply matrix vector product y=op*x
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in)    :: op
        class(vector_t),        intent(in)    :: x
        class(vector_t),        intent(inout) :: y 
    !-----------------------------------------------------------------
        type(serial_scalar_array_t), pointer :: y_serial
        call op%abort_if_not_in_domain(x)
        call op%abort_if_not_in_range(y)
        call x%GuardTemp()
        select type(x)
          class is (par_scalar_array_t)
          select type(y)
            class is(par_scalar_array_t)
               if(op%p_env%p_context%iam<0) return
               y_serial=> y%get_serial_scalar_array()
               call op%sparse_matrix%apply(x%get_serial_scalar_array(), y_serial)
               call y%comm()
           end select
        end select
        call x%CleanTemp()
    end subroutine par_sparse_matrix_apply

    subroutine par_sparse_matrix_free_in_stages(this, action)
    !-----------------------------------------------------------------
    !< free_in_stages procedure.
    !< As it extends from matrix_t, it must be implemented
    !-----------------------------------------------------------------
       class(par_sparse_matrix_t), intent(inout) :: this
       integer(ip)               ,   intent(in)  :: action
       if ( action == free_clean ) then
          !if ( this%state == created ) then
            if(this%p_env%p_context%iam>=0) then
              call this%sparse_matrix%free_in_stages(action)
            end if  
            nullify ( this%dof_import_domain )
            nullify ( this%dof_import_range )
            nullify ( this%p_env )
            call this%free_vector_spaces()
          !this%state = start
         !end if   
       else if ( action == free_symbolic_setup ) then
          if(this%p_env%p_context%iam>=0) then
            call this%sparse_matrix%free_in_stages(action)
          end if 
       else if ( action == free_numerical_setup ) then
          if(this%p_env%p_context%iam>=0) then
            call this%sparse_matrix%free_in_stages(action)
          end if     
       end if
    end subroutine par_sparse_matrix_free_in_stages

    subroutine par_sparse_matrix_print(this,lunou, only_graph)
    !-----------------------------------------------------------------
    !< Print a Sparse matrix
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t),  intent(in) :: this
        integer(ip),             intent(in) :: lunou
        logical,     optional,   intent(in) :: only_graph
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        if(present(only_graph)) then
            call this%sparse_matrix%print(lunou, only_graph=only_graph)
        else
            call this%sparse_matrix%print(lunou)
        endif
    end subroutine par_sparse_matrix_print

    subroutine par_sparse_matrix_print_matrix_market (this, lunou, ng, l2g)
    !-----------------------------------------------------------------
    !< Print a Sparse matrix in matrix market format
    !-----------------------------------------------------------------
        class(par_sparse_matrix_t), intent(in) :: this
        integer(ip),            intent(in) :: lunou
        integer(ip), optional,  intent(in) :: ng
        integer(ip), optional,  intent(in) :: l2g (*)
    !-----------------------------------------------------------------
        if(this%p_env%p_context%iam<0) return
        call this%sparse_matrix%print_matrix_market(lunou, ng, l2g)
    end subroutine par_sparse_matrix_print_matrix_market
    
    subroutine par_sparse_matrix_create_iterator(this, iblock, jblock, iterator)
      !-----------------------------------------------------------------
      !< Get a pointer to an iterator over the matrix entries
      !-----------------------------------------------------------------
      class(par_sparse_matrix_t)            , intent(in)    :: this
      integer(ip)                           , intent(in)    :: iblock 
      integer(ip)                           , intent(in)    :: jblock 
      class(matrix_iterator_t), allocatable , intent(inout) :: iterator
      !-----------------------------------------------------------------
      ! NOT IMPLEMENTED YET
      assert(.not. allocated(iterator))
      assert(.false.)
    end subroutine par_sparse_matrix_create_iterator

end module par_sparse_matrix_names
