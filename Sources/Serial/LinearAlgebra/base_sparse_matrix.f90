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
  integer(ip), parameter :: SPARSE_MATRIX_STATE_START              = 0
  integer(ip), parameter :: SPARSE_MATRIX_STATE_PROPERTIES_SET     = 1
  integer(ip), parameter :: SPARSE_MATRIX_STATE_CREATED            = 2
  integer(ip), parameter :: SPARSE_MATRIX_STATE_BUILD_SYMBOLIC     = 3
  integer(ip), parameter :: SPARSE_MATRIX_STATE_BUILD_NUMERIC      = 4
  integer(ip), parameter :: SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC = 5
  integer(ip), parameter :: SPARSE_MATRIX_STATE_ASSEMBLED          = 6
  integer(ip), parameter :: SPARSE_MATRIX_STATE_UPDATE             = 7

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
  ! Start               | Free_clean            | Start
  ! Start               | Free_symbolic         | Start
  ! Start               | Free_numeric          | Start
  
  ! Created             | Insert (2-values)     | Build_symbolic
  ! Created             | Insert (3-values)     | Build_numeric
  ! Created             | Free_clean            | Start
  ! Created             | Free_symbolic         | Created
  ! Created             | Free_numeric          | Assembled_symbolic
  
  ! Build_symbolic      | Insert (2-values)     | Build_symbolic
  ! Build_symbolic      | Insert (2-values)     | * Error
  ! Build_symbolic      | convert               | Assembled_symbolic
  ! Build_symbolic      | Free_clean            | Start
  ! Build_symbolic      | Free_symbolic         | Created
  ! Build_symbolic      | Free_numeric          | Created
  
  ! Build_numeric       | Insert (2-values)     | * Error
  ! Build_numeric       | Insert (3-values)     | Build_numeric
  ! Build_numeric       | convert               | Assembled
  ! Build_numeric       | Free_clean            | Start
  ! Build_numeric       | Free_symbolic         | Created
  ! Build_numeric       | Free_numeric          | Created
  !-----------------------------------------------------------------
  ! Under development
  !-----------------------------------------------------------------
  ! Assembled           | Free_clean            | Start
  ! Assembled           | Free_symbolic         | Created
  ! Assembled           | Free_numeric          | Assembled_numeric
  ! Assembled           | Set                   | Update
  
  ! Assembled_symbolic  | Free_clean            | Start
  ! Assembled_symbolic  | Set                   | Update
  
  ! Update              | Free_clean            | Start
  
  
  ! Matrix sign
  integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE     = 0
  integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_POSITIVE_SEMIDEFINITE = 1
  integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_INDEFINITE            = 2 ! Both positive and negative eigenvalues
  integer(ip), public, parameter :: SPARSE_MATRIX_SIGN_UNKNOWN               = 3 ! No info
  
  type, abstract :: base_sparse_matrix_t
     private 
     integer(ip) :: num_rows                          !< Number of rows
     integer(ip) :: num_cols                          !< Number of colums
     integer(ip) :: state = SPARSE_MATRIX_STATE_START !< Matrix state (one of SPARSE_MATRIX_STATE_XXX parameters)
     integer(ip) :: sign                              !< Matrix sign (one of SPARSE_MATRIX_SIGN_XXX pa
     logical     :: sum_duplicates = .true.           !< If .false. overwrites duplicated values, else perform sum
     logical     :: symmetric                         !< Matrix is symmetric (.true.) or not (.false.)
     logical     :: symmetric_storage = .false.       !< .True.   Implicitly assumes that G=(V,E) is such that 
     !<          (i,j) \belongs E <=> (j,i) \belongs E, forall i,j \belongs V.
     !<          Only edges (i,j) with j>=i are stored.
     !< .False.  All (i,j) \belongs E are stored.  
   contains
     private
        procedure(base_sparse_matrix_is_by_rows),                public, deferred :: is_by_rows
        procedure(base_sparse_matrix_is_by_cols),                public, deferred :: is_by_cols
        procedure(base_sparse_matrix_copy_to_coo),               public, deferred :: copy_to_coo
        procedure(base_sparse_matrix_copy_from_coo),             public, deferred :: copy_from_coo
        procedure(base_sparse_matrix_move_to_coo),               public, deferred :: move_to_coo
        procedure(base_sparse_matrix_move_from_coo),             public, deferred :: move_from_coo
        procedure(base_sparse_matrix_move_to_fmt),               public, deferred :: move_to_fmt
        procedure(base_sparse_matrix_move_from_fmt),             public, deferred :: move_from_fmt
        procedure(base_sparse_matrix_initialize_values),         public, deferred :: initialize_values
        procedure(base_sparse_matrix_allocate_values_body),      public, deferred :: allocate_values_body
        procedure(base_sparse_matrix_update_bounded_values_body),              &
                                                                 public, deferred :: update_bounded_values_body
        procedure(base_sparse_matrix_update_bounded_value_body),               &
                                                                 public, deferred :: update_bounded_value_body
        procedure(base_sparse_matrix_update_bounded_values_by_row_body),       &
                                                                 public, deferred :: update_bounded_values_by_row_body
        procedure(base_sparse_matrix_update_bounded_values_by_col_body),       &
                                                                 public, deferred :: update_bounded_values_by_col_body
        procedure(base_sparse_matrix_update_bounded_dense_values_body),        &
                                                                 public, deferred :: update_bounded_dense_values_body
        procedure(base_sparse_matrix_update_bounded_square_dense_values_body), &
                                                                 public, deferred :: update_bounded_square_dense_values_body
        procedure(base_sparse_matrix_update_dense_values_body),                &
                                                                 public, deferred :: update_dense_values_body
        procedure(base_sparse_matrix_update_square_dense_values_body),         &
                                                                 public, deferred :: update_square_dense_values_body
        procedure(base_sparse_matrix_update_values_body),        public, deferred :: update_values_body
        procedure(base_sparse_matrix_update_values_by_row_body), public, deferred :: update_values_by_row_body
        procedure(base_sparse_matrix_update_values_by_col_body), public, deferred :: update_values_by_col_body
        procedure(base_sparse_matrix_update_value_body),         public, deferred :: update_value_body
        procedure(base_sparse_matrix_split_2x2_symbolic),        public, deferred :: split_2x2_symbolic
        procedure(base_sparse_matrix_split_2x2_numeric),         public, deferred :: split_2x2_numeric
        procedure(base_sparse_matrix_permute_and_split_2x2_numeric),           &
                                                                 public, deferred :: permute_and_split_2x2_numeric
        procedure(base_sparse_matrix_permute_and_split_2x2_symbolic),          &
                                                                 public, deferred :: permute_and_split_2x2_symbolic
        procedure(base_sparse_matrix_expand_matrix_numeric),     public, deferred :: expand_matrix_numeric
        procedure(base_sparse_matrix_expand_matrix_symbolic),    public, deferred :: expand_matrix_symbolic
        procedure(base_sparse_matrix_extract_diagonal),          public, deferred :: extract_diagonal
        procedure(base_sparse_matrix_print_matrix_market_body),  public, deferred :: print_matrix_market_body
        procedure(base_sparse_matrix_free_coords),               public, deferred :: free_coords
        procedure(base_sparse_matrix_free_val),                  public, deferred :: free_val
        procedure(base_sparse_matrix_set_nnz),                   public, deferred :: set_nnz
        procedure(base_sparse_matrix_get_nnz),                   public, deferred :: get_nnz
        procedure(base_sparse_matrix_print),                     public, deferred :: print
        procedure         ::                                     base_sparse_matrix_create_square
        procedure         ::                                     base_sparse_matrix_create_rectangular
        procedure         :: insert_bounded_coords            => base_sparse_matrix_insert_bounded_coords
        procedure         :: insert_bounded_values            => base_sparse_matrix_insert_bounded_values
        procedure         :: insert_bounded_coords_by_row     => base_sparse_matrix_insert_bounded_coords_by_row
        procedure         :: insert_bounded_coords_by_col     => base_sparse_matrix_insert_bounded_coords_by_col
        procedure         :: insert_bounded_values_by_row     => base_sparse_matrix_insert_bounded_values_by_row
        procedure         :: insert_bounded_values_by_col     => base_sparse_matrix_insert_bounded_values_by_col
        procedure         :: insert_bounded_single_coord      => base_sparse_matrix_insert_bounded_single_coord
        procedure         :: insert_bounded_single_value      => base_sparse_matrix_insert_bounded_single_value
        procedure         :: insert_bounded_dense_values      => base_sparse_matrix_insert_bounded_dense_values
        procedure         :: insert_bounded_square_dense_values=> base_sparse_matrix_insert_bounded_square_dense_values
        procedure         :: insert_coords                    => base_sparse_matrix_insert_coords
        procedure         :: insert_values                    => base_sparse_matrix_insert_values
        procedure         :: insert_dense_values              => base_sparse_matrix_insert_dense_values
        procedure         :: insert_square_dense_values       => base_sparse_matrix_insert_square_dense_values
        procedure         :: insert_coords_by_row             => base_sparse_matrix_insert_coords_by_row
        procedure         :: insert_coords_by_col             => base_sparse_matrix_insert_coords_by_col
        procedure         :: insert_values_by_row             => base_sparse_matrix_insert_values_by_row
        procedure         :: insert_values_by_col             => base_sparse_matrix_insert_values_by_col
        procedure         :: insert_single_coord              => base_sparse_matrix_insert_single_coord
        procedure         :: insert_single_value              => base_sparse_matrix_insert_single_value
        procedure         :: append_bounded_coords_body       => base_sparse_matrix_append_bounded_coords_body
        procedure         :: append_bounded_values_body       => base_sparse_matrix_append_bounded_values_body
        procedure         :: append_bounded_coords_by_row_body=> base_sparse_matrix_append_bounded_coords_by_row_body
        procedure         :: append_bounded_coords_by_col_body=> base_sparse_matrix_append_bounded_coords_by_col_body
        procedure         :: append_bounded_values_by_row_body=> base_sparse_matrix_append_bounded_values_by_row_body
        procedure         :: append_bounded_values_by_col_body=> base_sparse_matrix_append_bounded_values_by_col_body
        procedure         :: append_bounded_single_coord_body => base_sparse_matrix_append_bounded_single_coord_body
        procedure         :: append_bounded_single_value_body => base_sparse_matrix_append_bounded_single_value_body
        procedure         :: append_bounded_dense_values_body => base_sparse_matrix_append_bounded_dense_values_body
        procedure         :: append_bounded_square_dense_values_body => base_sparse_matrix_append_bounded_square_dense_values_body
        procedure         :: append_coords_body               => base_sparse_matrix_append_coords_body
        procedure         :: append_values_body               => base_sparse_matrix_append_values_body
        procedure         :: append_dense_values_body         => base_sparse_matrix_append_dense_values_body
        procedure         :: append_square_dense_values_body  => base_sparse_matrix_append_square_dense_values_body
        procedure         :: append_coords_by_row_body        => base_sparse_matrix_append_coords_by_row_body
        procedure         :: append_coords_by_col_body        => base_sparse_matrix_append_coords_by_col_body
        procedure         :: append_values_by_row_body        => base_sparse_matrix_append_values_by_row_body
        procedure         :: append_values_by_col_body        => base_sparse_matrix_append_values_by_col_body
        procedure         :: append_single_coord_body         => base_sparse_matrix_append_single_coord_body
        procedure         :: append_single_value_body         => base_sparse_matrix_append_single_value_body
        procedure         :: is_valid_sign                    => base_sparse_matrix_is_valid_sign
        procedure         :: apply_body                       => base_sparse_matrix_apply_body
        procedure, public :: is_symbolic                      => base_sparse_matrix_is_symbolic
        procedure, public :: copy_to_fmt                      => base_sparse_matrix_copy_to_fmt
        procedure, public :: copy_from_fmt                    => base_sparse_matrix_copy_from_fmt
        procedure, public :: set_sign                         => base_sparse_matrix_set_sign
        procedure, public :: get_sign                         => base_sparse_matrix_get_sign
        procedure, public :: set_num_rows                     => base_sparse_matrix_set_num_rows
        procedure, public :: get_num_rows                     => base_sparse_matrix_get_num_rows
        procedure, public :: set_num_cols                     => base_sparse_matrix_set_num_cols
        procedure, public :: get_num_cols                     => base_sparse_matrix_get_num_cols
        procedure, public :: set_sum_duplicates               => base_sparse_matrix_set_sum_duplicates
        procedure, public :: get_sum_duplicates               => base_sparse_matrix_get_sum_duplicates
        procedure, public :: set_symmetry                     => base_sparse_matrix_set_symmetry
        procedure, public :: set_symmetric_storage            => base_sparse_matrix_set_symmetric_storage
        procedure, public :: get_symmetric_storage            => base_sparse_matrix_get_symmetric_storage
        procedure, public :: is_symmetric                     => base_sparse_matrix_is_symmetric
        procedure, public :: set_state                        => base_sparse_matrix_set_state
        procedure, public :: set_state_start                  => base_sparse_matrix_set_state_start
        procedure, public :: set_state_properties_set         => base_sparse_matrix_set_state_properties_set
        procedure, public :: set_state_created                => base_sparse_matrix_set_state_created
        procedure, public :: set_state_build_symbolic         => base_sparse_matrix_set_state_build_symbolic
        procedure, public :: set_state_build_numeric          => base_sparse_matrix_set_state_build_numeric
        procedure, public :: set_state_assembled              => base_sparse_matrix_set_state_assembled
        procedure, public :: set_state_assembled_symbolic     => base_sparse_matrix_set_state_assembled_symbolic
        procedure, public :: set_state_update                 => base_sparse_matrix_set_state_update
        procedure, public :: state_is_start                   => base_sparse_matrix_state_is_start
        procedure, public :: state_is_properties_setted       => base_sparse_matrix_state_is_properties_setted
        procedure, public :: state_is_created                 => base_sparse_matrix_state_is_created
        procedure, public :: state_is_build_symbolic          => base_sparse_matrix_state_is_build_symbolic
        procedure, public :: state_is_build_numeric           => base_sparse_matrix_state_is_build_numeric
        procedure, public :: state_is_assembled               => base_sparse_matrix_state_is_assembled
        procedure, public :: state_is_assembled_symbolic      => base_sparse_matrix_state_is_assembled_symbolic
        procedure, public :: state_is_update                  => base_sparse_matrix_state_is_update
        procedure, public :: get_state                        => base_sparse_matrix_get_state
        procedure, public :: allocate_coords                  => base_sparse_matrix_allocate_coords
        procedure, public :: allocate_values                  => base_sparse_matrix_allocate_values
        procedure, public :: convert_body                     => base_sparse_matrix_convert_body
        procedure, public :: apply                            => base_sparse_matrix_apply
        procedure, public :: print_matrix_market              => base_sparse_matrix_print_matrix_market
        procedure, public :: free                             => base_sparse_matrix_free
        procedure, public :: free_clean                       => base_sparse_matrix_free_clean
        procedure, public :: free_symbolic                    => base_sparse_matrix_free_symbolic
        procedure, public :: free_numeric                     => base_sparse_matrix_free_numeric
        procedure, public :: set_properties                   => base_sparse_matrix_set_properties
        generic,   public :: create                           => base_sparse_matrix_create_square, &
                                                                 base_sparse_matrix_create_rectangular
        generic,   public :: insert                           => insert_bounded_coords,              &
                                                                 insert_bounded_values,              &
                                                                 insert_bounded_coords_by_row,       &
                                                                 insert_bounded_coords_by_col,       &
                                                                 insert_bounded_values_by_row,       &
                                                                 insert_bounded_values_by_col,       &
                                                                 insert_bounded_single_value,        &
                                                                 insert_bounded_single_coord,        &
                                                                 insert_bounded_dense_values,        &
                                                                 insert_bounded_square_dense_values, &
                                                                 insert_coords,                      &
                                                                 insert_values,                      &
                                                                 insert_dense_values,                &
                                                                 insert_square_dense_values,         &
                                                                 insert_coords_by_row,               &
                                                                 insert_coords_by_col,               &
                                                                 insert_values_by_row,               &
                                                                 insert_values_by_col,               &
                                                                 insert_single_value,                &
                                                                 insert_single_coord
        generic           :: append_body                      => append_bounded_coords_body,             &  
                                                                 append_bounded_values_body,             &
                                                                 append_bounded_coords_by_row_body,      &
                                                                 append_bounded_coords_by_col_body,      &
                                                                 append_bounded_values_by_row_body,      &
                                                                 append_bounded_values_by_col_body,      &
                                                                 append_bounded_single_value_body,       &
                                                                 append_bounded_single_coord_body,       &
                                                                 append_bounded_dense_values_body,         &
                                                                 append_bounded_square_dense_values_body,  &
                                                                 append_coords_body,                     &
                                                                 append_values_body,                     &
                                                                 append_dense_values_body,               &
                                                                 append_square_dense_values_body,        &
                                                                 append_coords_by_row_body,              &
                                                                 append_coords_by_col_body,              &
                                                                 append_values_by_row_body,              &
                                                                 append_values_by_col_body,              &
                                                                 append_single_value_body,               &
                                                                 append_single_coord_body
        generic,   public :: update_body                      => update_bounded_values_body ,              &
                                                                 update_bounded_values_by_row_body,        &
                                                                 update_bounded_values_by_col_body,        &
                                                                 update_bounded_value_body,                &
                                                                 update_bounded_dense_values_body ,        &
                                                                 update_bounded_square_dense_values_body , &
                                                                 update_dense_values_body ,                &
                                                                 update_square_dense_values_body ,         &
                                                                 update_values_body ,                      &
                                                                 update_values_by_row_body ,               &
                                                                 update_values_by_col_body ,               &
                                                                 update_value_body
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
        integer(ip), allocatable   :: ia(:)                       !< Row indices
        integer(ip), allocatable   :: ja(:)                       !< Column indices        
        real(rp),    allocatable   :: val(:)                      !< Values
    contains
    private
        procedure         :: append_bounded_values_body              => coo_sparse_matrix_append_bounded_values_body
        procedure         :: append_bounded_coords_body              => coo_sparse_matrix_append_bounded_coords_body
        procedure         :: append_bounded_values_by_row_body       => coo_sparse_matrix_append_bounded_values_by_row_body
        procedure         :: append_bounded_values_by_col_body       => coo_sparse_matrix_append_bounded_values_by_col_body
        procedure         :: append_bounded_coords_by_row_body       => coo_sparse_matrix_append_bounded_coords_by_row_body
        procedure         :: append_bounded_coords_by_col_body       => coo_sparse_matrix_append_bounded_coords_by_col_body
        procedure         :: append_bounded_single_value_body        => coo_sparse_matrix_append_bounded_single_value_body
        procedure         :: append_bounded_single_coord_body        => coo_sparse_matrix_append_bounded_single_coord_body
        procedure         :: append_bounded_dense_values_body        => coo_sparse_matrix_append_bounded_dense_values_body
        procedure         :: append_bounded_square_dense_values_body => coo_sparse_matrix_append_bounded_square_dense_values_body
        procedure         :: append_values_body                      => coo_sparse_matrix_append_values_body
        procedure         :: append_coords_body                      => coo_sparse_matrix_append_coords_body
        procedure         :: append_dense_values_body                => coo_sparse_matrix_append_dense_values_body
        procedure         :: append_square_dense_values_body         => coo_sparse_matrix_append_square_dense_values_body
        procedure         :: append_values_by_row_body               => coo_sparse_matrix_append_values_by_row_body
        procedure         :: append_values_by_col_body               => coo_sparse_matrix_append_values_by_col_body
        procedure         :: append_coords_by_row_body               => coo_sparse_matrix_append_coords_by_row_body
        procedure         :: append_coords_by_col_body               => coo_sparse_matrix_append_coords_by_col_body
        procedure         :: append_single_value_body                => coo_sparse_matrix_append_single_value_body
        procedure         :: append_single_coord_body                => coo_sparse_matrix_append_single_coord_body
        procedure, public :: update_bounded_values_body              => coo_sparse_matrix_update_bounded_values_body
        procedure, public :: update_bounded_values_by_row_body       => coo_sparse_matrix_update_bounded_values_by_row_body
        procedure, public :: update_bounded_values_by_col_body       => coo_sparse_matrix_update_bounded_values_by_col_body
        procedure, public :: update_bounded_value_body               => coo_sparse_matrix_update_bounded_value_body
        procedure, public :: update_bounded_dense_values_body        => coo_sparse_matrix_update_bounded_dense_values_body
        procedure, public :: update_bounded_square_dense_values_body => coo_sparse_matrix_update_bounded_square_dense_values_body
        procedure, public :: update_dense_values_body                => coo_sparse_matrix_update_dense_values_body
        procedure, public :: update_square_dense_values_body         => coo_sparse_matrix_update_square_dense_values_body
        procedure, public :: update_values_body                      => coo_sparse_matrix_update_values_body
        procedure, public :: update_values_by_row_body               => coo_sparse_matrix_update_values_by_row_body
        procedure, public :: update_values_by_col_body               => coo_sparse_matrix_update_values_by_col_body
        procedure, public :: update_value_body                       => coo_sparse_matrix_update_value_body
        procedure, public :: split_2x2_symbolic                      => coo_sparse_matrix_split_2x2_symbolic
        procedure, public :: split_2x2_numeric                       => coo_sparse_matrix_split_2x2_numeric
        procedure, public :: permute_and_split_2x2_numeric           => coo_sparse_matrix_permute_and_split_2x2_numeric
        procedure, public :: permute_and_split_2x2_symbolic          => coo_sparse_matrix_permute_and_split_2x2_symbolic
        procedure, public :: expand_matrix_numeric                   => coo_sparse_matrix_expand_matrix_numeric
        procedure, public :: expand_matrix_symbolic                  => coo_sparse_matrix_expand_matrix_symbolic
        procedure, public :: extract_diagonal                        => coo_sparse_matrix_extract_diagonal
        procedure, public :: is_by_rows                              => coo_sparse_matrix_is_by_rows
        procedure, public :: is_by_cols                              => coo_sparse_matrix_is_by_cols
        procedure, public :: set_nnz                                 => coo_sparse_matrix_set_nnz
        procedure, public :: get_nnz                                 => coo_sparse_matrix_get_nnz
        procedure, public :: sort_and_compress                       => coo_sparse_matrix_sort_and_compress
        procedure, public :: set_sort_status_none                    => coo_sparse_matrix_set_sort_status_none
        procedure, public :: set_sort_status_by_rows                 => coo_sparse_matrix_set_sort_status_by_rows
        procedure, public :: set_sort_status_by_cols                 => coo_sparse_matrix_set_sort_status_by_cols
        procedure, public :: get_sort_status                         => coo_sparse_matrix_get_sort_status
        procedure, public :: allocate_coords                         => coo_sparse_matrix_allocate_coords
        procedure, public :: allocate_values_body                    => coo_sparse_matrix_allocate_values_body
        procedure, public :: initialize_values                       => coo_sparse_matrix_initialize_values
        procedure, public :: copy_to_coo                             => coo_sparse_matrix_copy_to_coo
        procedure, public :: copy_from_coo                           => coo_sparse_matrix_copy_from_coo
        procedure, public :: copy_to_fmt                             => coo_sparse_matrix_copy_to_fmt
        procedure, public :: copy_from_fmt                           => coo_sparse_matrix_copy_from_fmt
        procedure, public :: move_to_coo                             => coo_sparse_matrix_move_to_coo
        procedure, public :: move_from_coo                           => coo_sparse_matrix_move_from_coo
        procedure, public :: move_to_fmt                             => coo_sparse_matrix_move_to_fmt
        procedure, public :: move_from_fmt                           => coo_sparse_matrix_move_from_fmt
        procedure, public :: free_coords                             => coo_sparse_matrix_free_coords
        procedure, public :: free_val                                => coo_sparse_matrix_free_val
        procedure, public :: print_matrix_market_body                => coo_sparse_matrix_print_matrix_market_body
        procedure, public :: print                                   => coo_sparse_matrix_print
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

        subroutine base_sparse_matrix_set_nnz(this, nnz)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nnz
        end subroutine base_sparse_matrix_set_nnz

        function base_sparse_matrix_get_nnz(this) result(nnz)
            import base_sparse_matrix_t
            import ip    
            class(base_sparse_matrix_t), intent(in) :: this
            integer(ip)                             :: nnz
        end function base_sparse_matrix_get_nnz

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

        subroutine base_sparse_matrix_allocate_values_body(this, nz)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),  intent(inout) :: this
            integer(ip), optional,        intent(in)    :: nz
        end subroutine base_sparse_matrix_allocate_values_body

        subroutine base_sparse_matrix_initialize_values(this, val)
            import base_sparse_matrix_t
            import rp
            class(base_sparse_matrix_t),  intent(inout) :: this
            real(rp),                     intent(in)    :: val
        end subroutine base_sparse_matrix_initialize_values

        subroutine base_sparse_matrix_update_bounded_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia(nz)
            integer(ip),                 intent(in)    :: ja(nz)
            real(rp),                    intent(in)    :: val(nz)
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_values_body

        subroutine base_sparse_matrix_update_bounded_value_body(this, ia, ja, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: ia
            integer(ip),                 intent(in)    :: ja
            real(rp),                    intent(in)    :: val
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_value_body

        subroutine base_sparse_matrix_update_bounded_values_by_row_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia
            integer(ip),                 intent(in)    :: ja(nz)
            real(rp),                    intent(in)    :: val(nz)
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_values_by_row_body

        subroutine base_sparse_matrix_update_bounded_values_by_col_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia(nz)
            integer(ip),                 intent(in)    :: ja
            real(rp),                    intent(in)    :: val(nz)
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_values_by_col_body

        subroutine base_sparse_matrix_update_bounded_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: num_rows
            integer(ip),                 intent(in)    :: num_cols
            integer(ip),                 intent(in)    :: ia(num_rows)
            integer(ip),                 intent(in)    :: ja(num_cols)
            integer(ip),                 intent(in)    :: ioffset
            integer(ip),                 intent(in)    :: joffset
            real(rp),                    intent(in)    :: val(:, :)
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_dense_values_body

        subroutine base_sparse_matrix_update_bounded_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: num_rows
            integer(ip),                 intent(in)    :: ia(num_rows)
            integer(ip),                 intent(in)    :: ja(num_rows)
            integer(ip),                 intent(in)    :: ioffset
            integer(ip),                 intent(in)    :: joffset
            real(rp),                    intent(in)    :: val(:, :)
            integer(ip),                 intent(in)    :: imin
            integer(ip),                 intent(in)    :: imax
            integer(ip),                 intent(in)    :: jmin
            integer(ip),                 intent(in)    :: jmax
        end subroutine base_sparse_matrix_update_bounded_square_dense_values_body

        subroutine base_sparse_matrix_update_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: num_rows
            integer(ip),                 intent(in)    :: num_cols
            integer(ip),                 intent(in)    :: ia(num_rows)
            integer(ip),                 intent(in)    :: ja(num_cols)
            integer(ip),                 intent(in)    :: ioffset
            integer(ip),                 intent(in)    :: joffset
            real(rp),                    intent(in)    :: val(:, :)
        end subroutine base_sparse_matrix_update_dense_values_body

        subroutine base_sparse_matrix_update_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: num_rows
            integer(ip),                 intent(in)    :: ia(num_rows)
            integer(ip),                 intent(in)    :: ja(num_rows)
            integer(ip),                 intent(in)    :: ioffset
            integer(ip),                 intent(in)    :: joffset
            real(rp),                    intent(in)    :: val(:, :)
        end subroutine base_sparse_matrix_update_square_dense_values_body


        subroutine base_sparse_matrix_update_values_body(this, nz, ia, ja, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia(nz)
            integer(ip),                 intent(in)    :: ja(nz)
            real(rp),                    intent(in)    :: val(nz)
        end subroutine base_sparse_matrix_update_values_body


        subroutine base_sparse_matrix_update_value_body(this, ia, ja, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: ia
            integer(ip),                 intent(in)    :: ja
            real(rp),                    intent(in)    :: val
        end subroutine base_sparse_matrix_update_value_body

        subroutine base_sparse_matrix_update_values_by_row_body(this, nz, ia, ja, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia
            integer(ip),                 intent(in)    :: ja(nz)
            real(rp),                    intent(in)    :: val(nz)
        end subroutine base_sparse_matrix_update_values_by_row_body

        subroutine base_sparse_matrix_update_values_by_col_body(this, nz, ia, ja, val) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t), intent(inout) :: this
            integer(ip),                 intent(in)    :: nz
            integer(ip),                 intent(in)    :: ia(nz)
            integer(ip),                 intent(in)    :: ja
            real(rp),                    intent(in)    :: val(nz)
        end subroutine base_sparse_matrix_update_values_by_col_body

        subroutine base_sparse_matrix_split_2x2_symbolic(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),           intent(in)    :: this
            integer(ip),                           intent(in)    :: num_row
            integer(ip),                           intent(in)    :: num_col
            class(base_sparse_matrix_t),           intent(inout) :: A_II
            class(base_sparse_matrix_t),           intent(inout) :: A_IG
            class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
            class(base_sparse_matrix_t),           intent(inout) :: A_GG
        end subroutine base_sparse_matrix_split_2x2_symbolic

        subroutine base_sparse_matrix_split_2x2_numeric(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),           intent(in)    :: this
            integer(ip),                           intent(in)    :: num_row
            integer(ip),                           intent(in)    :: num_col
            class(base_sparse_matrix_t),           intent(inout) :: A_II
            class(base_sparse_matrix_t),           intent(inout) :: A_IG
            class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
            class(base_sparse_matrix_t),           intent(inout) :: A_GG
        end subroutine base_sparse_matrix_split_2x2_numeric

        subroutine base_sparse_matrix_permute_and_split_2x2_numeric(this, num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR) 
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t),           intent(in)    :: this
            integer(ip),                           intent(in)    :: num_row
            integer(ip),                           intent(in)    :: num_col
            integer(ip),                           intent(in)    :: perm(:)
            integer(ip),                           intent(in)    :: iperm(:)
            real(rp),    allocatable,              intent(out)   :: A_CC(:,:)
            real(rp),    allocatable,              intent(out)   :: A_CR(:,:)
            real(rp),    allocatable,              intent(out)   :: A_RC(:,:)
            class(base_sparse_matrix_t),           intent(inout) :: A_RR
        end subroutine base_sparse_matrix_permute_and_split_2x2_numeric

        subroutine base_sparse_matrix_permute_and_split_2x2_symbolic(this, num_row, num_col, perm, iperm, A_RR) 
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),           intent(in)    :: this
            integer(ip),                           intent(in)    :: num_row
            integer(ip),                           intent(in)    :: num_col
            integer(ip),                           intent(in)    :: perm(:)
            integer(ip),                           intent(in)    :: iperm(:)
            class(base_sparse_matrix_t),           intent(inout) :: A_RR
        end subroutine base_sparse_matrix_permute_and_split_2x2_symbolic

        subroutine base_sparse_matrix_expand_matrix_numeric(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to)
            import base_sparse_matrix_t
            import ip
            import rp
            class(base_sparse_matrix_t),     intent(in)    :: this
            integer,                         intent(in)    :: C_T_num_cols
            integer,                         intent(in)    :: C_T_nz
            integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
            integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
            real(rp),                        intent(in)    :: C_T_val(C_T_nz)
            integer(ip),                     intent(in)    :: I_nz
            integer(ip),                     intent(in)    :: I_ia(I_nz)
            integer(ip),                     intent(in)    :: I_ja(I_nz)
            real(rp),                        intent(in)    :: I_val(C_T_nz)
            class(base_sparse_matrix_t),     intent(inout) :: to
        end subroutine base_sparse_matrix_expand_matrix_numeric

        subroutine base_sparse_matrix_expand_matrix_symbolic(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),     intent(in)    :: this
            integer,                         intent(in)    :: C_T_num_cols
            integer,                         intent(in)    :: C_T_nz
            integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
            integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
            integer(ip),                     intent(in)    :: I_nz
            integer(ip),                     intent(in)    :: I_ia(I_nz)
            integer(ip),                     intent(in)    :: I_ja(I_nz)
            class(base_sparse_matrix_t),     intent(inout) :: to
        end subroutine base_sparse_matrix_expand_matrix_symbolic

        subroutine base_sparse_matrix_extract_diagonal(this, diagonal)
            import base_sparse_matrix_t
            import  rp
            class(base_sparse_matrix_t),  intent(in)    :: this
            real(rp), allocatable,        intent(inout) :: diagonal(:)
        end subroutine base_sparse_matrix_extract_diagonal

        subroutine base_sparse_matrix_free_coords(this)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t),  intent(inout) :: this
        end subroutine base_sparse_matrix_free_coords

        subroutine base_sparse_matrix_free_val(this)
            import base_sparse_matrix_t
            class(base_sparse_matrix_t),  intent(inout) :: this
        end subroutine base_sparse_matrix_free_val

        subroutine base_sparse_matrix_print_matrix_market_body (this, lunou, ng, l2g)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t), intent(in) :: this
            integer(ip),                 intent(in) :: lunou
            integer(ip), optional,       intent(in) :: ng
            integer(ip), optional,       intent(in) :: l2g (*)
        end subroutine base_sparse_matrix_print_matrix_market_body

        subroutine base_sparse_matrix_print(this,lunou, only_graph)
            import base_sparse_matrix_t
            import ip
            class(base_sparse_matrix_t),  intent(in) :: this
            integer(ip),                  intent(in) :: lunou
            logical, optional,            intent(in) :: only_graph
        end subroutine base_sparse_matrix_print
    end interface

!---------------------------------------------------------------------
!< AUX PROCEDURES INTERFACES
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
public :: duplicates_operation
public :: assign_value
public :: sum_value
public :: binary_search
public :: mergesort_link_list
public :: reorder_ip_rp_from_link_list
public :: reorder_ip_from_link_list


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
    !< This state can be reached after free() calls
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_START
    end subroutine base_sparse_matrix_set_state_start


    subroutine base_sparse_matrix_set_state_properties_set(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_PROPERTIES_SET
    !< The matrix can jump to this state after set_properties() calls
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_PROPERTIES_SET
    end subroutine base_sparse_matrix_set_state_properties_set


    subroutine base_sparse_matrix_set_state_created(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after create() calls
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_CREATED
    end subroutine base_sparse_matrix_set_state_created


    subroutine base_sparse_matrix_set_state_build_symbolic(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after insert(ia,ja) calls
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_BUILD_SYMBOLIC
    end subroutine base_sparse_matrix_set_state_build_symbolic


    subroutine base_sparse_matrix_set_state_build_numeric(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_CREATED
    !< The matrix can jump to this state after inser(ia,ja,val) calls
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_BUILD_NUMERIC
    end subroutine base_sparse_matrix_set_state_build_numeric


    subroutine base_sparse_matrix_set_state_assembled_symbolic(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_ASSEMBLED
    !< The matrix can be in this state after convert() calls from a
    !< SPARSE_MATRIX_STATE_BUILD_SYMBOLIC state or callfree_numeric()
    !< from SPARSE_MATRIX_STATE_ASSEMBLED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC
    end subroutine base_sparse_matrix_set_state_assembled_symbolic


    subroutine base_sparse_matrix_set_state_assembled(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_ASSEMBLED
    !< The matrix can be in this state after convert() calls from a
    !< SPARSE_MATRIX_STATE_BUILD_NUMERIC state or 
    !< from SPARSE_MATRIX_STATE_UPDATE
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_ASSEMBLED
    end subroutine base_sparse_matrix_set_state_assembled


    subroutine base_sparse_matrix_set_state_update(this)
    !-----------------------------------------------------------------
    !< Set the matrix state to SPARSE_MATRIX_STATE_UPDATE
    !< The matrix can be in this state after SPARSE_MATRIX_STATE_ASSEMBLED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%state = SPARSE_MATRIX_STATE_UPDATE
    end subroutine base_sparse_matrix_set_state_update


    function base_sparse_matrix_state_is_start(this) result(state_start)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_START
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_start
    !-----------------------------------------------------------------
        state_start = (this%state == SPARSE_MATRIX_STATE_START)
    end function base_sparse_matrix_state_is_start


    function base_sparse_matrix_state_is_properties_setted(this) result(state_properties_setted)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_PROPERTIES_SET
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_properties_setted
    !-----------------------------------------------------------------
        state_properties_setted = (this%state == SPARSE_MATRIX_STATE_PROPERTIES_SET)
    end function base_sparse_matrix_state_is_properties_setted


    function base_sparse_matrix_state_is_created(this) result(state_created)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_CREATED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_created
    !-----------------------------------------------------------------
        state_created = (this%state == SPARSE_MATRIX_STATE_CREATED)
    end function base_sparse_matrix_state_is_created


    function base_sparse_matrix_state_is_build_symbolic(this) result(state_build_symbolic)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_CREATED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_build_symbolic
    !-----------------------------------------------------------------
        state_build_symbolic = (this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
    end function base_sparse_matrix_state_is_build_symbolic


    function base_sparse_matrix_state_is_build_numeric(this) result(state_build_numeric)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_CREATED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_build_numeric
    !-----------------------------------------------------------------
        state_build_numeric = (this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC)
    end function base_sparse_matrix_state_is_build_numeric


    function base_sparse_matrix_state_is_assembled_symbolic(this) result(state_assembled_symbolic)
    !-----------------------------------------------------------------
    !< check if the matrix state is SPARSE_MATRIX_STATE_ASSEMBLED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_assembled_symbolic
    !-----------------------------------------------------------------
        state_assembled_symbolic = (this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC)
    end function base_sparse_matrix_state_is_assembled_symbolic


    function base_sparse_matrix_state_is_assembled(this) result(state_assembled)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_ASSEMBLED
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_assembled
    !-----------------------------------------------------------------
        state_assembled = (this%state == SPARSE_MATRIX_STATE_ASSEMBLED)
    end function base_sparse_matrix_state_is_assembled


    function base_sparse_matrix_state_is_update(this) result(state_update)
    !-----------------------------------------------------------------
    !< Check if the matrix state is SPARSE_MATRIX_STATE_UPDATE
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: state_update
    !-----------------------------------------------------------------
        state_update = (this%state == SPARSE_MATRIX_STATE_UPDATE)
    end function base_sparse_matrix_state_is_update


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


    subroutine base_sparse_matrix_set_sum_duplicates(this, sum_duplicates)
    !-----------------------------------------------------------------
    !< Set sum_duplicates value (.true. => sum, .false. => overwrite)
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        logical,                     intent(in)    :: sum_duplicates
    !-----------------------------------------------------------------
        this%sum_duplicates = sum_duplicates
    end subroutine base_sparse_matrix_set_sum_duplicates


    function base_sparse_matrix_get_sum_duplicates(this) result(sum_duplicates)
    !-----------------------------------------------------------------
    !< Get sum_duplicates value (.true. => sum, .false. => overwrite)
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: sum_duplicates
    !-----------------------------------------------------------------
        sum_duplicates = this%sum_duplicates
    end function base_sparse_matrix_get_sum_duplicates


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


    subroutine base_sparse_matrix_set_symmetric_storage(this, symmetric_storage)
    !-----------------------------------------------------------------
    !< Set symmetric storage property of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        logical,                     intent(in)    :: symmetric_storage
    !-----------------------------------------------------------------
        this%symmetric_storage = symmetric_storage
    end subroutine base_sparse_matrix_set_symmetric_storage


    function base_sparse_matrix_get_symmetric_storage(this) result(symmetric_storage)
    !-----------------------------------------------------------------
    !< Get symmetric storage property of the matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: symmetric_storage
    !-----------------------------------------------------------------
        symmetric_storage = this%symmetric_storage
    end function base_sparse_matrix_get_symmetric_storage


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


    function base_sparse_matrix_is_symbolic(this) result(symbolic)
    !-----------------------------------------------------------------
    !< Check if values are allocated/managed
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        logical                                 :: symbolic
    !-----------------------------------------------------------------
        symbolic = (this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC)
    end function base_sparse_matrix_is_symbolic


    subroutine base_sparse_matrix_set_properties(this,symmetric_storage,is_symmetric,sign)
    !-----------------------------------------------------------------
    !< Set the sign and symmetry properties of a sparse matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        logical,                     intent(in)    :: symmetric_storage
        logical,                     intent(in)    :: is_symmetric
        integer(ip),                 intent(in)    :: sign
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_START)
        assert(this%is_valid_sign(sign))
        if(symmetric_storage) then
            assert(is_symmetric)
        endif
        call this%set_symmetric_storage(symmetric_storage)
        this%symmetric = is_symmetric
        this%sign = sign
        call this%set_state_properties_set()
    end subroutine base_sparse_matrix_set_properties


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
        assert(this%state == SPARSE_MATRIX_STATE_START .or. this%state == SPARSE_MATRIX_STATE_PROPERTIES_SET)
        assert(this%is_valid_sign(sign))
        if(symmetric_storage) then
            assert(is_symmetric)
        endif
        call this%set_symmetric_storage(symmetric_storage)
        this%symmetric = is_symmetric
        this%sign = sign
        this%num_rows = num_rows_and_cols
        this%num_cols = num_rows_and_cols
        call this%allocate_coords(nz)
        call this%set_state_created()
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
        assert(this%state == SPARSE_MATRIX_STATE_START .or. this%state == SPARSE_MATRIX_STATE_PROPERTIES_SET)
        if(.not. this%state_is_properties_setted()) then
            call this%set_symmetric_storage(.false.)
            this%symmetric = .false.
            this%sign = SPARSE_MATRIX_SIGN_UNKNOWN
        endif
        this%num_rows = num_rows
        this%num_cols = num_cols
        call this%allocate_coords(nz)
        call this%set_state_created()
    end subroutine base_sparse_matrix_create_rectangular


    subroutine base_sparse_matrix_insert_bounded_values(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
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
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_values


    subroutine base_sparse_matrix_insert_bounded_coords(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
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
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja, imin, imax, jmin, jmax)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_bounded_coords


    subroutine base_sparse_matrix_insert_bounded_values_by_row(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_values_by_row


    subroutine base_sparse_matrix_insert_bounded_values_by_col(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_values_by_col


    subroutine base_sparse_matrix_insert_bounded_coords_by_row(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja, imin, imax, jmin, jmax)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_bounded_coords_by_row


    subroutine base_sparse_matrix_insert_bounded_coords_by_col(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja, imin, imax, jmin, jmax)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_bounded_coords_by_col

    
    subroutine base_sparse_matrix_insert_bounded_single_value(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(ia, ja, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(ia, ja, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_single_value


    subroutine base_sparse_matrix_insert_bounded_single_coord(this, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(ia, ja, imin, imax, jmin, jmax)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_bounded_single_coord


    subroutine base_sparse_matrix_insert_values(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_values


    subroutine base_sparse_matrix_insert_coords(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_coords


    subroutine base_sparse_matrix_insert_bounded_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: num_cols
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_cols)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_dense_values


    subroutine base_sparse_matrix_insert_bounded_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_rows)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_bounded_square_dense_values


    subroutine base_sparse_matrix_insert_dense_values(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: num_cols
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_cols)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(num_rows, num_cols, ia, ja, ioffset, joffset, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(num_rows, num_cols, ia, ja, ioffset, joffset, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_dense_values


    subroutine base_sparse_matrix_insert_square_dense_values(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_rows)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(num_rows, ia, ja, ioffset, joffset, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(num_rows, ia, ja, ioffset, joffset, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_square_dense_values


    subroutine base_sparse_matrix_insert_values_by_row(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_values_by_row


    subroutine base_sparse_matrix_insert_values_by_col(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(nz, ia, ja, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(nz, ia, ja, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_values_by_col


    subroutine base_sparse_matrix_insert_coords_by_row(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_coords_by_row


    subroutine base_sparse_matrix_insert_coords_by_col(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(nz, ia, ja)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_coords_by_col
    
    
    subroutine base_sparse_matrix_insert_single_value(this, ia, ja, val)
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) then
            call this%append_body(ia, ja, val)
            call this%set_state_build_numeric()
        else
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC) then
                call this%allocate_values()
                call this%initialize_values(val=0.0_rp)
            endif
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) then
                call this%update_body(ia, ja, val)
                call this%set_state_update()
            endif
        endif
    end subroutine base_sparse_matrix_insert_single_value


    subroutine base_sparse_matrix_insert_single_coord(this, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries to the sparse matrix
    !< This is a common interface to control the state diagram
    !< It delegates the insert of new entries on the append procedures
    !< Append procedures must be only overloaded on the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_CREATED .or. this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)
        call this%append_body(ia, ja)
        call this%set_state_build_symbolic()
    end subroutine base_sparse_matrix_insert_single_coord

    
    subroutine base_sparse_matrix_append_bounded_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
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
    end subroutine base_sparse_matrix_append_bounded_values_body


    subroutine base_sparse_matrix_append_bounded_coords_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
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
    end subroutine base_sparse_matrix_append_bounded_coords_body


    subroutine base_sparse_matrix_append_bounded_values_by_row_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_values_by_row_body


    subroutine base_sparse_matrix_append_bounded_values_by_col_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_values_by_col_body


    subroutine base_sparse_matrix_append_bounded_coords_by_row_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_coords_by_row_body


    subroutine base_sparse_matrix_append_bounded_coords_by_col_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_coords_by_col_body
    
    
    subroutine base_sparse_matrix_append_bounded_single_value_body(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append a new entry and value to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_single_value_body


    subroutine base_sparse_matrix_append_bounded_single_coord_body(this, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append a new entrie to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_single_coord_body


    subroutine base_sparse_matrix_append_values_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_values_body


   subroutine base_sparse_matrix_append_coords_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_coords_body


    subroutine base_sparse_matrix_append_values_by_row_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_values_by_row_body


    subroutine base_sparse_matrix_append_values_by_col_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val(nz)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_values_by_col_body


   subroutine base_sparse_matrix_append_coords_by_row_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja(nz)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_coords_by_row_body


   subroutine base_sparse_matrix_append_coords_by_col_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: nz
        integer(ip),                 intent(in)    :: ia(nz)
        integer(ip),                 intent(in)    :: ja
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_coords_by_col_body


    subroutine base_sparse_matrix_append_bounded_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: num_cols
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_cols)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_dense_values_body


    subroutine base_sparse_matrix_append_bounded_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin ,jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_rows)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
        integer(ip),                 intent(in)    :: imin
        integer(ip),                 intent(in)    :: imax
        integer(ip),                 intent(in)    :: jmin
        integer(ip),                 intent(in)    :: jmax
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_bounded_square_dense_values_body


    subroutine base_sparse_matrix_append_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: num_cols
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_cols)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_dense_values_body


    subroutine base_sparse_matrix_append_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: num_rows
        integer(ip),                 intent(in)    :: ia(num_rows)
        integer(ip),                 intent(in)    :: ja(num_rows)
        integer(ip),                 intent(in)    :: ioffset
        integer(ip),                 intent(in)    :: joffset
        real(rp),                    intent(in)    :: val(:, :)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_square_dense_values_body
    
    
    subroutine base_sparse_matrix_append_single_value_body(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append a new entry and value to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
        real(rp),                    intent(in)    :: val
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_single_value_body


    subroutine base_sparse_matrix_append_single_coord_body(this, ia, ja) 
    !-----------------------------------------------------------------
    !< Append a new entrie to the sparse matrix
    !< Must be overloaded only in the COO format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip),                 intent(in)    :: ia
        integer(ip),                 intent(in)    :: ja
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_append_single_coord_body


    subroutine base_sparse_matrix_allocate_coords(this, nz)
    !-----------------------------------------------------------------
    !< Allocate coords
    !< Must be overloaded 
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
        integer(ip), optional,       intent(in)    :: nz
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine base_sparse_matrix_allocate_coords


    subroutine base_sparse_matrix_allocate_values(this)
    !-----------------------------------------------------------------
    !< Allocate val
    !< Must be overloaded 
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED)
        if ( this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC ) then
          call this%allocate_values_body(this%get_nnz())
          call this%set_state_assembled()
        end if  
    end subroutine base_sparse_matrix_allocate_values

    subroutine base_sparse_matrix_apply(op, x, y)
    !-----------------------------------------------------------------
    !< Matrix vector product
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in)    :: op
        class(vector_t),             intent(in)    :: x
        class(vector_t),             intent(inout) :: y
    !-----------------------------------------------------------------
        assert(op%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. op%state == SPARSE_MATRIX_STATE_UPDATE)
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


    subroutine base_sparse_matrix_convert_body(this)
    !-----------------------------------------------------------------
    !< Convert procedure from the base_sparse_matrix point of view
    !< State diagram management. Sort and compress?
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED_SYMBOLIC .or. this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        if(this%is_symbolic()) then
            call this%set_state_assembled_symbolic()
        else
            call this%set_state_assembled()
        endif
    end subroutine base_sparse_matrix_convert_body


    subroutine base_sparse_matrix_free (this)
    !-----------------------------------------------------------------
    !< Clean the properties and size of a rectangular matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%free_clean()
    end subroutine base_sparse_matrix_free


    subroutine base_sparse_matrix_free_clean (this)
    !-----------------------------------------------------------------
    !< Clean the properties, size and arrays of a matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        this%sign = SPARSE_MATRIX_SIGN_UNKNOWN
        this%num_rows = 0
        this%num_cols = 0
        call this%set_nnz(nnz=0)
        this%symmetric = .false.
        this%symmetric_storage = .false.
        call this%free_coords()
        call this%free_val()
        call this%set_state_start()
    end subroutine base_sparse_matrix_free_clean


    subroutine base_sparse_matrix_free_symbolic (this)
    !-----------------------------------------------------------------
    !< Clean the coords and values of a matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        call this%free_coords()
        call this%free_val()
        call this%set_nnz(nnz=0)
        if(this%state /= SPARSE_MATRIX_STATE_START) call this%set_state_created()
    end subroutine base_sparse_matrix_free_symbolic


    subroutine base_sparse_matrix_free_numeric (this)
    !-----------------------------------------------------------------
    !< Clean the values of a matrix
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(inout) :: this
    !-----------------------------------------------------------------
        if(this%state > SPARSE_MATRIX_STATE_CREATED) then
            call this%free_val()
            if(this%state == SPARSE_MATRIX_STATE_BUILD_NUMERIC) call this%set_state_build_symbolic()
            if(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE) call this%set_state_assembled_symbolic()
        endif
    end subroutine base_sparse_matrix_free_numeric


    subroutine base_sparse_matrix_print_matrix_market (this, lunou, ng, l2g)
    !-----------------------------------------------------------------
    !< Print the Sparse matrix in matrix market format
    !-----------------------------------------------------------------
        class(base_sparse_matrix_t), intent(in) :: this
        integer(ip),                 intent(in) :: lunou
        integer(ip), optional,       intent(in) :: ng
        integer(ip), optional,       intent(in) :: l2g (*)
    !-----------------------------------------------------------------
        assert(this%state == SPARSE_MATRIX_STATE_ASSEMBLED .or. this%state == SPARSE_MATRIX_STATE_UPDATE)
        call this%print_matrix_market_body(lunou, ng, l2g)
    end subroutine base_sparse_matrix_print_matrix_market

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


    subroutine coo_sparse_matrix_allocate_coords(this, nz)
    !-----------------------------------------------------------------
    !< Allocate COO arrays
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        integer(ip), optional,      intent(in)     :: nz
    !-----------------------------------------------------------------
        if(present(nz)) then
            call memalloc(nz, this%ia,  __FILE__, __LINE__)
            call memalloc(nz, this%ja,  __FILE__, __LINE__)
        else
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%ia,  __FILE__, __LINE__)
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%ja,  __FILE__, __LINE__)
        endif
    end subroutine coo_sparse_matrix_allocate_coords


    subroutine coo_sparse_matrix_allocate_values_body(this, nz)
    !-----------------------------------------------------------------
    !< Allocate COO arrays
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        integer(ip), optional,      intent(in)     :: nz
    !-----------------------------------------------------------------
        check(.not. allocated(this%val))
        if(present(nz)) then
            call memalloc(nz, this%val, __FILE__, __LINE__)
        else
            call memalloc(max(7*this%num_rows, 7*this%num_cols, 1), this%val,  __FILE__, __LINE__)
        endif
    end subroutine coo_sparse_matrix_allocate_values_body


    subroutine coo_sparse_matrix_initialize_values(this, val)
    !-----------------------------------------------------------------
    !< Initialize COO values
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        real(rp),                   intent(in)     :: val
    !-----------------------------------------------------------------
        if(allocated(this%val)) this%val(1:this%get_nnz()) = val
    end subroutine coo_sparse_matrix_initialize_values


    subroutine coo_sparse_matrix_append_bounded_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
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
        if(.not. allocated(this%val)) call this%allocate_values_body(size(this%ia))
        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Append the new entries and values
        do i=1, nz
            ir = ia(i); ic = ja(i)
            if(ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or.             & ! Check imposed bounds
               ir<1 .or. ir>this%num_rows .or. ic<1 .or. ic>this%num_cols .or. & ! Check matrix bounds
               (this%symmetric_storage .and. ir>ic) ) cycle                      ! Check if has symmetric_storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ic
            this%val(nnz) = val(i)
            !write(*,*) 'Inserted:', ir, ic, val(i)
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_values_body


    subroutine coo_sparse_matrix_append_bounded_coords_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new coords to the COO sparse matrix
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
        ! Append the new entries
        do i=1, nz
            ir = ia(i); ic = ja(i)
            if(ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or.             & ! Check imposed bounds
               ir<1 .or. ir>this%num_rows .or. ic<1 .or. ic>this%num_cols .or. & ! Check matrix bounds
               (this%symmetric_storage .and. ir>ic) ) cycle                      ! Check if has symmetric_storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ic
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_coords_body


    subroutine coo_sparse_matrix_append_bounded_values_by_row_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0 .or. ia<imin .or. ia>imax .or. ia<1 .or. ia>this%num_rows) return
        nnz = this%nnz
        newnnz = nnz + nz 
        if(.not. allocated(this%val)) call this%allocate_values_body(size(this%ia))
        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Append the new entries and values
        do i=1, nz
            ic = ja(i)
            if(ic<jmin .or. ic>jmax .or.             &      ! Check imposed column bounds
               ic<1 .or. ic>this%num_cols .or.       &      ! Check matrix column bounds
               (this%symmetric_storage .and. ia>ic) ) cycle ! Check if has symmetric_storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ia
            this%ja(nnz) = ic
            this%val(nnz) = val(i)
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_values_by_row_body


    subroutine coo_sparse_matrix_append_bounded_values_by_col_body(this, nz, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ir, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0 .or. ja<jmin .or. ja>jmax .or. ja<1 .or. ja>this%num_cols) return
        nnz = this%nnz
        newnnz = nnz + nz 
        if(.not. allocated(this%val)) call this%allocate_values_body(size(this%ia))
        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Append the new entries and values
        do i=1, nz
            ir = ia(i)
            if(ir<imin .or. ir>imax .or.             &      ! Check imposed row bounds
               ir<1 .or. ir>this%num_rows .or.       &      ! Check matrix row bounds
               (this%symmetric_storage .and. ir>ja) ) cycle ! Check if has symmetric storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ja
            this%val(nnz) = val(i)
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_values_by_col_body


    subroutine coo_sparse_matrix_append_bounded_coords_by_row_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0 .or. ia<imin .or. ia>imax .or. ia<1 .or. ia>this%num_rows) return
        nnz = this%nnz
        newnnz = nnz + nz 
        ! Realloc this%ia, this%ja to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
        endif
        ! Append the new entries
        do i=1, nz
            ic = ja(i)
            if(ic<jmin .or. ic>jmax .or.            &       ! Check imposed column bounds
               ic<1 .or. ic>this%num_cols .or.      &       ! Check matrix column bounds
               (this%symmetric_storage .and. ia>ic) ) cycle ! Check if has symmetric storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ia
            this%ja(nnz) = ic
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_coords_by_row_body


    subroutine coo_sparse_matrix_append_bounded_coords_by_col_body(this, nz, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, ir, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        if(nz == 0 .or. ja<jmin .or. ja>jmax .or. ja<1 .or. ja>this%num_cols) return
        nnz = this%nnz
        newnnz = nnz + nz 
        ! Realloc this%ia, this%ja to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
        endif
        ! Append the new entries
        do i=1, nz
            ir = ia(i)
            if(ir<imin .or. ir>imax .or.            &       ! Check imposed row bounds
               ir<1 .or. ir>this%num_rows .or.      &       ! Check matrix row bounds
               (this%symmetric_storage .and. ir>ja) ) cycle ! Check if has symmetric storage
                !If symmetric_storage is .true. only the upper triangle is stored
            nnz = nnz + 1
            this%ia(nnz) = ir
            this%ja(nnz) = ja
        enddo
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_coords_by_col_body

    
    subroutine coo_sparse_matrix_append_bounded_single_value_body(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entrie and value to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: nnz, newnnz, newsize
    !-----------------------------------------------------------------
        nnz = this%nnz
        newnnz = nnz + 1
        if(.not. allocated(this%val)) call this%allocate_values_body(size(this%ia))
        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Append the new entries and values
        if(ia<imin .or. ia>imax .or. ja<jmin .or. ja>jmax .or.             & ! Check imposed bounds
           ia<1 .or. ia>this%num_rows .or. ja<1 .or. ja>this%num_cols .or. & ! Check matrix bounds
               (this%symmetric_storage .and. ia>ja) ) return                 ! Check if has symmetric storage
        !If symmetric_storage is .true. only the upper triangle is stored
        nnz = nnz+1
        this%ia(nnz) = ia
        this%ja(nnz) = ja
        this%val(nnz) = val
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_single_value_body


    subroutine coo_sparse_matrix_append_bounded_single_coord_body(this, ia, ja, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new coord to the COO sparse matrix
    !< It allows duplicates entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: nnz, newnnz, newsize
    !-----------------------------------------------------------------
        nnz = this%nnz
        newnnz = nnz+1

        ! Realloc this%ia, this%ja, this%val to the right size if is needed
        if( size(this%ia) < newnnz) then
            newsize = max(newnnz, int(1.5*size(this%ia)))
            call memrealloc(newsize, this%ia,  __FILE__, __LINE__)
            call memrealloc(newsize, this%ja,  __FILE__, __LINE__)
            call memrealloc(newsize, this%val, __FILE__, __LINE__)
        endif
        ! Append the new entries
        if(ia<imin .or. ia>imax .or. ja<jmin .or. ja>jmax .or.             & ! Check imposed bounds
           ia<1 .or. ia>this%num_rows .or. ja<1 .or. ja>this%num_cols .or. & ! Check matrix bounds
               (this%symmetric_storage .and. ia>ja) ) return                 ! Check if has symmetric storage
        !If symmetric_storage is .true. only the upper triangle is stored
        nnz = nnz + 1
        this%ia(nnz) = ia
        this%ja(nnz) = ja
        this%nnz = nnz
        this%sort_status = COO_SPARSE_MATRIX_SORTED_NONE
    end subroutine coo_sparse_matrix_append_bounded_single_coord_body    


    subroutine coo_sparse_matrix_append_values_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip)                               :: nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_values_body


    subroutine coo_sparse_matrix_append_coords_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new coords to the COO sparse matrix
    !< It allows duplicates entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        integer(ip)                               :: nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_coords_body


    subroutine coo_sparse_matrix_append_bounded_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols
            do i=1, num_rows
                call this%append_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
            enddo
        enddo
    end subroutine coo_sparse_matrix_append_bounded_dense_values_body


    subroutine coo_sparse_matrix_append_bounded_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
          do i=1, num_rows
            call this%append_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
          enddo
        enddo
    end subroutine coo_sparse_matrix_append_bounded_square_dense_values_body


    subroutine coo_sparse_matrix_append_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols
            do i=1, num_rows
                call this%append_body(ia(i), ja(j), val(i+ioffset,j+joffset),1 , this%num_rows, 1, this%num_cols)
            enddo
        enddo
    end subroutine coo_sparse_matrix_append_dense_values_body


    subroutine coo_sparse_matrix_append_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
            do i=1, num_rows
                call this%append_body(ia(i), ja(j), val(i+ioffset,j+joffset),1 , this%num_rows, 1, this%num_cols)
            enddo
        enddo
    end subroutine coo_sparse_matrix_append_square_dense_values_body


    subroutine coo_sparse_matrix_append_values_by_row_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_values_by_row_body


    subroutine coo_sparse_matrix_append_values_by_col_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_values_by_col_body


    subroutine coo_sparse_matrix_append_coords_by_row_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_coords_by_row_body


    subroutine coo_sparse_matrix_append_coords_by_col_body(this, nz, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new entries and values to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(nz, ia, ja, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_coords_by_col_body

    
    subroutine coo_sparse_matrix_append_single_value_body(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Append new entrie and value to the COO sparse matrix
    !< It allows duplicated entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        integer(ip)                               :: nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_single_value_body


    subroutine coo_sparse_matrix_append_single_coord_body(this, ia, ja) 
    !-----------------------------------------------------------------
    !< Append new coord to the COO sparse matrix
    !< It allows duplicates entries
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        integer(ip)                               :: i, ir, ic, nnz, newnnz, newsize
    !-----------------------------------------------------------------
        call this%append_body(ia, ja, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_append_single_coord_body    


    subroutine coo_sparse_matrix_update_bounded_values_body(this, nz, ia, ja, val, imin, imax, jmin, jmax)
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
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
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ic, ilr, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        if(nz==0) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        nnz = this%nnz

        if(this%is_by_rows()) then
            ilr = -1

            do i=1, nz
                ir = ia(i)
                ic = ja(i) 
                ! Ignore out of bounds entries
                if (ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or. &
                    ir<1 .or. ir>this%num_rows .or. ic<1 .or. ic>this%num_cols .or. &
                    (this%symmetric_storage .and. ic>ir)) cycle
                if (ir /= ilr) then 
                    i1 = binary_search(ir,nnz,this%ia)
                    i2 = i1
                    do 
                        if (i2+1 > nnz) exit
                        if (this%ia(i2+1) /= this%ia(i2)) exit
                        i2 = i2 + 1
                    end do
                    do 
                        if (i1-1 < 1) exit
                        if (this%ia(i1-1) /= this%ia(i1)) exit
                        i1 = i1 - 1
                    end do
                    ilr = ir
                    nc = i2-i1+1
                    ipaux = binary_search(ic,nc,this%ja(i1:i2))
                    assert(ipaux>0) ! Entry not found
                    if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
                end if
            end do
        elseif(this%is_by_cols()) then
            !Not tested yet!
            ilc = -1

            do i=1, nz
                ir = ia(i)
                ic = ja(i) 
                ! Ignore out of bounds entries
                if (ir<imin .or. ir>imax .or. ic<jmin .or. ic>jmax .or. &
                    ir<1 .or. ir>this%num_rows .or. ic<1 .or. ic>this%num_cols .or. &
                    (this%symmetric_storage .and. ic>ir)) cycle
                if (ic /= ilc) then 
                    i1 = binary_search(ic,nnz,this%ja)
                    i2 = i1
                    do 
                        if (i2+1 > nnz) exit
                        if (this%ja(i2+1) /= this%ja(i2)) exit
                        i2 = i2 + 1
                    end do
                    do 
                        if (i1-1 < 1) exit
                        if (this%ja(i1-1) /= this%ja(i1)) exit
                        i1 = i1 - 1
                    end do
                    ilc = ic
                    nr = i2-i1+1
                    ipaux = binary_search(ir,nc,this%ia(i1:i2))
                    assert(ipaux>0) ! Entry not found
                    if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
                end if
            end do
        endif
    end subroutine coo_sparse_matrix_update_bounded_values_body


    subroutine coo_sparse_matrix_update_bounded_value_body(this, ia, ja, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        ! Ignore out of bounds entriy
        if (ia<imin .or. ia>imax .or. ja<jmin .or. ja>jmax .or. &
            ia<1 .or. ia>this%num_rows .or. ja<1 .or. ja>this%num_cols .or. &
            (this%symmetric_storage .and. ja>ia)) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        nnz = this%nnz
        if(this%is_by_rows()) then
            i1 = binary_search(ia,nnz,this%ia)
            i2 = i1
            do 
                if (i2+1 > nnz) exit
                if (this%ia(i2+1) /= this%ia(i2)) exit
                i2 = i2 + 1
            end do
            do 
                if (i1-1 < 1) exit
                if (this%ia(i1-1) /= this%ia(i1)) exit
                i1 = i1 - 1
            end do
            nc = i2-i1+1
            ipaux = binary_search(ja,nc,this%ja(i1:i2))
            assert(ipaux>0) ! Entry not found
            if (ipaux>0) call apply_duplicates(input=val, output=this%val(i1+ipaux-1))

        elseif(this%is_by_rows()) then
            ! Not tested yet! 
            i1 = binary_search(ja,nnz,this%ja)
            i2 = i1
            do 
                if (i2+1 > nnz) exit
                if (this%ja(i2+1) /= this%ja(i2)) exit
                i2 = i2 + 1
            end do
            do 
                if (i1-1 < 1) exit
                if (this%ja(i1-1) /= this%ja(i1)) exit
                i1 = i1 - 1
            end do
            nr = i2-i1+1
            ipaux = binary_search(ia,nc,this%ia(i1:i2))
            assert(ipaux>0) ! Entry not found
            if (ipaux>0) call apply_duplicates(input=val, output=this%val(i1+ipaux-1))
        endif
    end subroutine coo_sparse_matrix_update_bounded_value_body


    subroutine coo_sparse_matrix_update_bounded_values_by_row_body(this, nz, ia, ja, val, imin, imax, jmin, jmax)
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i, ic, ilr, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        if(nz==0 .or. ia<imin .or. ia<1 .or. ia>imax .or. ia>this%num_rows ) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        nnz = this%nnz

        if(this%is_by_rows()) then
            ! Search the range for this row number
            i1 = binary_search(ia,nnz,this%ia)
            if(i1<1) return ! Row not found
            i2 = i1
            do 
                if (i2+1 > nnz) exit
                if (this%ia(i2+1) /= this%ia(i2)) exit
                i2 = i2 + 1
            end do
            do 
                if (i1-1 < 1) exit
                if (this%ia(i1-1) /= this%ia(i1)) exit
                i1 = i1 - 1
            end do
            ! Search each column
            do i=1, nz
                ic = ja(i) 
                ! Ignore out of bounds entries
                if (ic<jmin .or. ic>jmax .or. ic<1 .or. ic>this%num_cols .or. &
                    (this%symmetric_storage .and. ic>ia)) cycle
                nc = i2-i1+1
                ipaux = binary_search(ic,nc,this%ja(i1:i2))
                assert(ipaux>0) ! Entry not found
                if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
            end do
        elseif(this%is_by_cols()) then
            !Not tested yet!
            ilc = -1

            do i=1, nz
                ic = ja(i) 
                ! Ignore out of bounds entries
                if (ic<jmin .or. ic>jmax .or. ic<1 .or. ic>this%num_cols .or. &
                    (this%symmetric_storage .and. ic>ia)) cycle
                if (ic /= ilc) then 
                    i1 = binary_search(ic,nnz,this%ja)
                    i2 = i1
                    do 
                        if (i2+1 > nnz) exit
                        if (this%ja(i2+1) /= this%ja(i2)) exit
                        i2 = i2 + 1
                    end do
                    do 
                        if (i1-1 < 1) exit
                        if (this%ja(i1-1) /= this%ja(i1)) exit
                        i1 = i1 - 1
                    end do
                    ilc = ic
                end if
                nr = i2-i1+1
                ipaux = binary_search(ia,nc,this%ia(i1:i2))
                assert(ipaux>0) ! Entry not found
                if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
            end do
        endif
    end subroutine coo_sparse_matrix_update_bounded_values_by_row_body


    subroutine coo_sparse_matrix_update_bounded_values_by_col_body(this, nz, ia, ja, val, imin, imax, jmin, jmax)
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir, ilr, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        if(nz==0 .or. ja<jmin .or. ja<1 .or. ja>jmax .or. ja>this%num_cols ) return

        if(this%get_sum_duplicates()) then
            apply_duplicates => sum_value
        else
            apply_duplicates => assign_value
        endif

        nnz = this%nnz

        if(this%is_by_rows()) then
            ilr = -1

            do i=1, nz
                ir = ia(i) 
                ! Ignore out of bounds entries
                if (ir<imin .or. ir>imax .or. ir<1 .or. ir>this%num_rows .or. &
                    (this%symmetric_storage .and. ja>ir)) cycle
                if (ir /= ilr) then 
                    i1 = binary_search(ir,nnz,this%ia)
                    i2 = i1
                    do 
                        if (i2+1 > nnz) exit
                        if (this%ia(i2+1) /= this%ia(i2)) exit
                        i2 = i2 + 1
                    end do
                    do 
                        if (i1-1 < 1) exit
                        if (this%ia(i1-1) /= this%ia(i1)) exit
                        i1 = i1 - 1
                    end do
                    ilr = ir
                end if
                nc = i2-i1+1
                ipaux = binary_search(ja,nc,this%ja(i1:i2))
                assert(ipaux>0) ! Entry not found
                if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
            end do
        elseif(this%is_by_cols()) then
            ! Not tested yet!

            ! Search the range for this column number
            i1 = binary_search(ja,nnz,this%ja)
            if(i1<1) return ! Column not found
            i2 = i1
            do 
                if (i2+1 > nnz) exit
                if (this%ja(i2+1) /= this%ja(i2)) exit
                i2 = i2 + 1
            end do
            do 
                if (i1-1 < 1) exit
                if (this%ja(i1-1) /= this%ja(i1)) exit
                i1 = i1 - 1
            end do
            ! Search each row
            do i=1, nz
                ir = ia(i) 
                ! Ignore out of bounds entries
                if (ir<imin .or. ir>imax .or. ir<1 .or. ir>this%num_rows .or. &
                    (this%symmetric_storage .and. ja>ir)) cycle
                nr = i2-i1+1
                ipaux = binary_search(ir,nc,this%ia(i1:i2))
                assert(ipaux>0) ! Entry not found
                if (ipaux>0) call apply_duplicates(input=val(i), output=this%val(i1+ipaux-1))
            end do
        endif
    end subroutine coo_sparse_matrix_update_bounded_values_by_col_body



    subroutine coo_sparse_matrix_update_values_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir,ic, ilr, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_update_values_body


    subroutine coo_sparse_matrix_update_bounded_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
            enddo
        enddo
    end subroutine coo_sparse_matrix_update_bounded_dense_values_body


    subroutine coo_sparse_matrix_update_bounded_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val, imin, imax, jmin, jmax) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip),                intent(in)    :: imin
        integer(ip),                intent(in)    :: imax
        integer(ip),                intent(in)    :: jmin
        integer(ip),                intent(in)    :: jmax
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), imin, imax, jmin, jmax)
            enddo
        enddo
    end subroutine coo_sparse_matrix_update_bounded_square_dense_values_body


    subroutine coo_sparse_matrix_update_dense_values_body(this, num_rows, num_cols, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: num_cols
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_cols)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1 .or. num_cols<1) return
        do j=1, num_cols
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%num_rows, 1, this%num_cols)
            enddo
        enddo
    end subroutine coo_sparse_matrix_update_dense_values_body


    subroutine coo_sparse_matrix_update_square_dense_values_body(this, num_rows, ia, ja, ioffset, joffset, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: num_rows
        integer(ip),                intent(in)    :: ia(num_rows)
        integer(ip),                intent(in)    :: ja(num_rows)
        integer(ip),                intent(in)    :: ioffset
        integer(ip),                intent(in)    :: joffset
        real(rp),                   intent(in)    :: val(:, :)
        integer(ip)                               :: i, j
    !-----------------------------------------------------------------
        if(num_rows<1) return
        do j=1, num_rows
            do i=1, num_rows
                call this%update_body(ia(i), ja(j), val(i+ioffset,j+joffset), 1, this%num_rows, 1, this%num_cols)
            enddo
        enddo
    end subroutine coo_sparse_matrix_update_square_dense_values_body


    subroutine coo_sparse_matrix_update_values_by_row_body(this, nz, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja(nz)
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i, ic, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_update_values_by_row_body


    subroutine coo_sparse_matrix_update_values_by_col_body(this, nz, ia, ja, val)
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: nz
        integer(ip),                intent(in)    :: ia(nz)
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val(nz)
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i,ir, ilr, ilc, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        call this%update_body(nz, ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_update_values_by_col_body


    subroutine coo_sparse_matrix_update_value_body(this, ia, ja, val) 
    !-----------------------------------------------------------------
    !< Update the values and entries in the sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout) :: this
        integer(ip),                intent(in)    :: ia
        integer(ip),                intent(in)    :: ja
        real(rp),                   intent(in)    :: val
        procedure(duplicates_operation), pointer  :: apply_duplicates => null ()
        integer(ip)                               :: i, ipaux,i1,i2,nr,nc,nnz
    !-----------------------------------------------------------------
        call this%update_body(ia, ja, val, 1, this%num_rows, 1, this%num_cols)
    end subroutine coo_sparse_matrix_update_value_body


    subroutine coo_sparse_matrix_sort_and_compress(this, by_cols)
    !-----------------------------------------------------------------
    !< Sort ia, ja and val by rows as default or by_cols if forced
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
        logical,     optional,      intent(in)     :: by_cols
        logical                                    :: by_rows 
    !-----------------------------------------------------------------
        assert(this%state >= SPARSE_MATRIX_STATE_CREATED)
               
        by_rows     = .true.

        if(present(by_cols)) by_rows = .not. by_cols

        if(this%nnz == 0) then
            if(by_rows) then
                this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_ROWS
            else
                this%sort_status = COO_SPARSE_MATRIX_SORTED_BY_COLS
            endif
            if(.not. allocated(this%val)) allocate(this%val(this%nnz))
            call this%set_state_assembled()
            return
        endif

        if(this%is_symbolic()) then
            call sort_and_compress_symbolic(this, .not. by_rows)
            call this%set_state_assembled_symbolic()
        else
            call sort_and_compress_numeric(this, .not. by_rows)
            call this%set_state_assembled()
        endif
    contains


        subroutine sort_and_compress_symbolic(this, by_cols)
        !-------------------------------------------------------------
        !< Sort ia, ja and val by rows as default or by_cols if forced
        !-------------------------------------------------------------
            use types_names
            class(coo_sparse_matrix_t), intent(inout)  :: this
            logical,     optional,      intent(in)     :: by_cols
            integer(ip), allocatable                   :: ias(:)
            integer(ip), allocatable                   :: jas(:)
            integer(ip), allocatable                   :: iaux(:)
            integer(ip), allocatable                   :: ix2(:)
            logical                                    :: by_rows 
            logical                                    :: sorted
            logical                                    :: use_buffers
            integer(ip)                                :: i, j, k, nnz, nzl, imx, iret, ipaux, is, irw, icl
            integer                                    :: info
            procedure(duplicates_operation), pointer   :: apply_duplicates => null ()
        !-------------------------------------------------------------
            by_rows     = .true.
            sorted      = .true.
            use_buffers = .true.
            if(present(by_cols)) by_rows = .not. by_cols

            if(this%get_sum_duplicates()) then
                apply_duplicates => sum_value
            else
                apply_duplicates => assign_value
            endif

            nnz = this%nnz

            call memalloc(max(this%num_rows, this%num_cols, nnz) + 2, iaux, __FILE__, __LINE__)
            allocate(ias(nnz),  &
                     jas(nnz),  &
                     ix2(max(this%num_rows, this%num_cols, nnz) + 2), stat=info)

            use_buffers = (info == 0)
            iaux = 0

            if(by_rows) then
            ! By-rows sorting
                if(use_buffers) then
                ! If there is enough memory space
                    iaux(this%ia(1)) = iaux(this%ia(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        iaux(this%ia(i)) = iaux(this%ia(i)) + 1
                        sorted = sorted .and. (this%ia(i-1) <= this%ia(i))
                    enddo
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
                                call mergesort_link_list(nzl,this%ja(i:imx),ix2, iret)
                                if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, this%ia(i:imx), this%ja(i:imx), ix2)
                                k = k + 1
                                this%ia(k)  = this%ia(i)
                                this%ja(k)  = this%ja(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
                                    if(this%ia(i) == irw .and. this%ja(i) == icl) then
                                        ! Duplicated values action: assign the last value
                                    else
                                        k = k + 1
                                        this%ia(k)  = this%ia(i)
                                        this%ja(k)  = this%ja(i)
                                        irw = this%ia(k)
                                        icl = this%ja(k)
                                    endif
                                enddo
                            endif
                        enddo
                    else
                    ! If rows are NOT sorted
                        ipaux = iaux(1)
                        iaux(1) = 0
                        do i=2, this%num_rows
                            is = iaux(i)
                            iaux(i) = ipaux
                            ipaux = ipaux + is
                        enddo
                        iaux(this%num_rows + 1) = ipaux

                        do i=1, nnz
                            irw = this%ia(i)
                            ipaux = iaux(irw) + 1
                            ias(ipaux) = this%ia(i)
                            jas(ipaux) = this%ja(i)
                            iaux(irw) = ipaux
                        enddo

                        k = 0
                        i = 1
                        do j=1, this%num_rows
                            nzl = iaux(j) - i + 1
                            imx = i + nzl - 1

                            if(nzl > 0) then
                                ! Sort the colums of a particular row
                                call mergesort_link_list(nzl,jas(i:imx),ix2, iret)
                                if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, ias(i:imx), jas(i:imx), ix2)
                                k = k + 1
                                this%ia(k)  = ias(i)
                                this%ja(k)  = jas(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
                                    if(ias(i) == irw .and. jas(i) == icl) then
                                    else
                                        k = k + 1
                                        this%ia(k)  = ias(i)
                                        this%ja(k)  = jas(i)
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
                    call mergesort_link_list(nnz,this%ia, iaux, iret)
                    if(iret == 0) call reorder_symbolic_coo_from_link_list(nnz, this%ia, this%ja, iaux)
                    i = 1
                    j = i
                    do while (i <= nnz)
                        do while (this%ia(j) == this%ia(i))
                            j = j + 1
                            if(j > nnz) exit
                        enddo
                        nzl = j - i
                        ! Sort the colums of a particular row
                        call mergesort_link_list(nzl, this%ja(i:), iaux, iret)
                        if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                        i = j
                    enddo

                    i = 1
                    j = 1
                    irw = this%ia(i)
                    icl = this%ja(i)

                    do 
                        j = j + 1
                        if(j > nnz) exit
                        if(this%ia(j) == irw .and. this%ja(j) == icl) then
                            ! Duplicated values: assign the last value
                        else
                            i = i + 1
                            this%ia(i)  = this%ia(j)
                            this%ja(i)  = this%ja(j)
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
                    iaux(this%ja(1)) = iaux(this%ja(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        iaux(this%ja(i)) = iaux(this%ja(i)) + 1
                        sorted = sorted .and. (this%ja(i-1) <= this%ja(i))
                    enddo
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
                                call mergesort_link_list(nzl,this%ia(i:imx),ix2, iret)
                                if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, this%ia(i:imx), this%ja(i:imx), ix2)
                                k = k + 1
                                this%ia(k)  = this%ia(i)
                                this%ja(k)  = this%ja(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
                                    if(this%ia(i) == irw .and. this%ja(i) == icl) then
                                        ! Duplicated values action: assign the last value or sum the values
                                    else
                                        k = k + 1
                                        this%ia(k)  = this%ia(i)
                                        this%ja(k)  = this%ja(i)
                                        irw = this%ia(k)
                                        icl = this%ja(k)
                                    endif
                                enddo
                            endif
                        enddo
                    else
                    ! If columns are NOT sorted
                        ipaux = iaux(1)
                        iaux(1) = 0
                        do i=2, this%num_cols
                            is = iaux(i)
                            iaux(i) = ipaux
                            ipaux = ipaux + is
                        enddo
                        iaux(this%num_cols + 1) = ipaux

                        do i=1, nnz
                            icl = this%ja(i)
                            ipaux = iaux(icl) + 1
                            ias(ipaux) = this%ia(i)
                            jas(ipaux) = this%ja(i)
                            iaux(icl) = ipaux
                        enddo

                        k = 0
                        i = 1
                        do j=1, this%num_cols
                            nzl = iaux(j) - i + 1
                            imx = i + nzl - 1

                            if(nzl > 0) then
                                ! Sort the rows of a particular column
                                call mergesort_link_list(nzl,ias(i:imx),ix2, iret)
                                if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, ias(i:imx), jas(i:imx), ix2)
                                k = k + 1
                                this%ia(k)  = ias(i)
                                this%ja(k)  = jas(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
                                    if(ias(i) == irw .and. jas(i) == icl) then
                                        ! Duplicated values: assign the last value or sum the values
                                    else
                                        k = k + 1
                                        this%ia(k)  = ias(i)
                                        this%ja(k)  = jas(i)
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
                    call mergesort_link_list(nnz,this%ja, iaux, iret)
                    if(iret == 0) call reorder_symbolic_coo_from_link_list(nnz, this%ia, this%ja, iaux)
                    i = 1
                    j = i
                    do while (i <= nnz)
                        do while (this%ja(j) == this%ja(i))
                            j = j + 1
                            if(j > nnz) exit
                        enddo
                        nzl = j - i
                        ! Sort the rows of a particular column
                        call mergesort_link_list(nzl, this%ia(i:), iaux, iret)
                        if(iret == 0) call reorder_symbolic_coo_from_link_list(nzl, this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                        i = j
                    enddo

                    i = 1
                    j = 1
                    irw = this%ia(i)
                    icl = this%ja(i)

                    do 
                        j = j + 1
                        if(j > nnz) exit
                        if(this%ia(j) == irw .and. this%ja(j) == icl) then
                            ! Duplicated values: assign the last value or sum the values
                        else
                            i = i + 1
                            this%ia(i)  = this%ia(j)
                            this%ja(i)  = this%ja(j)
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
            this%nnz = i
            call memrealloc(nnz, this%ia,  __FILE__, __LINE__)
            call memrealloc(nnz, this%ja,  __FILE__, __LINE__)

        end subroutine sort_and_compress_symbolic


        subroutine sort_and_compress_numeric(this, by_cols)
        !-------------------------------------------------------------
        !< Sort ia, ja and val by rows as default or by_cols if forced
        !-------------------------------------------------------------
            use types_names
            class(coo_sparse_matrix_t), intent(inout)  :: this
            logical,     optional,      intent(in)     :: by_cols
            integer(ip), allocatable                   :: ias(:)
            integer(ip), allocatable                   :: jas(:)
            real(rp),    allocatable                   :: vs(:)
            integer(ip), allocatable                   :: iaux(:)
            integer(ip), allocatable                   :: ix2(:)
            logical                                    :: by_rows 
            logical                                    :: sorted
            logical                                    :: use_buffers
            integer(ip)                                :: i, j, k, nnz, nzl, imx, iret, ipaux, is, irw, icl
            integer                                    :: info
            procedure(duplicates_operation), pointer   :: apply_duplicates => null ()
        !-------------------------------------------------------------
            by_rows     = .true.
            sorted      = .true.
            use_buffers = .true.
            if(present(by_cols)) by_rows = .not. by_cols

            if(this%get_sum_duplicates()) then
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
                    iaux(this%ia(1)) = iaux(this%ia(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        iaux(this%ia(i)) = iaux(this%ia(i)) + 1
                        sorted = sorted .and. (this%ia(i-1) <= this%ia(i))
                    enddo
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
                                call mergesort_link_list(nzl,this%ja(i:imx),ix2, iret)
                                if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, this%val(i:imx), this%ia(i:imx), this%ja(i:imx), ix2)
                                ! If row/column index out of range ignore it
                                if(i > imx) exit
                                k = k + 1
                                this%ia(k)  = this%ia(i)
                                this%ja(k)  = this%ja(i)
                                this%val(k) = this%val(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
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
                        ipaux = iaux(1)
                        iaux(1) = 0
                        do i=2, this%num_rows
                            is = iaux(i)
                            iaux(i) = ipaux
                            ipaux = ipaux + is
                        enddo
                        iaux(this%num_rows + 1) = ipaux

                        do i=1, nnz
                            irw = this%ia(i)
                            ipaux = iaux(irw) + 1
                            ias(ipaux) = this%ia(i)
                            jas(ipaux) = this%ja(i)
                            vs(ipaux)  = this%val(i)
                            iaux(irw) = ipaux
                        enddo

                        k = 0
                        i = 1
                        do j=1, this%num_rows
                            nzl = iaux(j) - i + 1
                            imx = i + nzl - 1

                            if(nzl > 0) then
                                ! Sort the colums of a particular row
                                call mergesort_link_list(nzl,jas(i:imx),ix2, iret)
                                if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, vs(i:imx), ias(i:imx), jas(i:imx), ix2)
                                if(i > imx) exit
                                k = k + 1
                                this%ia(k)  = ias(i)
                                this%ja(k)  = jas(i)
                                this%val(k) = vs(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
                                    if(ias(i) == irw .and. jas(i) == icl) then
                                        ! Duplicated values: assign the last value
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
                    ! Sort the rows
                    call mergesort_link_list(nnz,this%ia, iaux, iret)
                    if(iret == 0) call reorder_numeric_coo_from_link_list(nnz, this%val, this%ia, this%ja, iaux)
                    i = 1
                    j = i
                    do while (i <= nnz)
                        do while (this%ia(j) == this%ia(i))
                            j = j + 1
                            if(j > nnz) exit
                        enddo
                        nzl = j - i
                        ! Sort the colums of a particular row
                        call mergesort_link_list(nzl, this%ja(i:), iaux, iret)
                        if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, this%val(i:i+nzl-1), this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                        i = j
                    enddo

                    i = 1
                    j = 1
                    irw = this%ia(i)
                    icl = this%ja(i)
                    do 
                        j = j + 1
                        if(j > nnz) exit
                        if(this%ia(j) == irw .and. this%ja(j) == icl) then
                            ! Duplicated values: assign the last value
                            call apply_duplicates(input=this%val(j), output=this%val(i))
                        else
                            i = i + 1
                            this%ia(i)  = this%ia(j)
                            this%ja(i)  = this%ja(j)
                            this%val(i)  = this%val(j)
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
                    iaux(this%ja(1)) = iaux(this%ja(1))+1
                    sorted = .true.
                    do i = 2, nnz
                        iaux(this%ja(i)) = iaux(this%ja(i)) + 1
                        sorted = sorted .and. (this%ja(i-1) <= this%ja(i))
                    enddo
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
                                call mergesort_link_list(nzl,this%ia(i:imx),ix2, iret)
                                if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, this%val(i:imx), this%ia(i:imx), this%ja(i:imx), ix2)
                                ! If row/column index out of range ignore it
                                if(i > imx) exit
                                k = k + 1
                                this%ia(k)  = this%ia(i)
                                this%ja(k)  = this%ja(i)
                                this%val(k) = this%val(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
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
                        ipaux = iaux(1)
                        iaux(1) = 0
                        do i=2, this%num_cols
                            is = iaux(i)
                            iaux(i) = ipaux
                            ipaux = ipaux + is
                        enddo
                        iaux(this%num_cols + 1) = ipaux

                        do i=1, nnz
                            icl = this%ja(i)
                            ipaux = iaux(icl) + 1
                            ias(ipaux) = this%ia(i)
                            jas(ipaux) = this%ja(i)
                            vs(ipaux)  = this%val(i)
                            iaux(icl) = ipaux
                        enddo

                        k = 0
                        i = 1
                        do j=1, this%num_cols
                            nzl = iaux(j) - i + 1
                            imx = i + nzl - 1

                            if(nzl > 0) then
                                ! Sort the rows of a particular column
                                call mergesort_link_list(nzl,ias(i:imx),ix2, iret)
                                if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, vs(i:imx), ias(i:imx), jas(i:imx), ix2)
                                if(i > imx) exit
                                k = k + 1
                                this%ia(k)  = ias(i)
                                this%ja(k)  = jas(i)
                                this%val(k) = vs(i)
                                irw = this%ia(k)
                                icl = this%ja(k)
                                do 
                                    i=i+1
                                    if(i > imx) exit
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
                    call mergesort_link_list(nnz,this%ja, iaux, iret)
                    if(iret == 0) call reorder_numeric_coo_from_link_list(nnz, this%val, this%ia, this%ja, iaux)
                    i = 1
                    j = i
                    do while (i <= nnz)
                        do while (this%ja(j) == this%ja(i))
                            j = j + 1
                            if(j > nnz) exit
                        enddo
                        nzl = j - i
                        ! Sort the rows of a particular column
                        call mergesort_link_list(nzl, this%ia(i:), iaux, iret)
                        if(iret == 0) call reorder_numeric_coo_from_link_list(nzl, this%val(i:i+nzl-1), this%ia(i:i+nzl-1), this%ja(i:i+nzl-1), iaux)
                        i = j
                    enddo

                    i = 1
                    j = 1
                    irw = this%ia(i)
                    icl = this%ja(i)

                    do 
                        j = j + 1
                        if(j > nnz) exit
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

        end subroutine sort_and_compress_numeric

    end subroutine coo_sparse_matrix_sort_and_compress


    subroutine coo_sparse_matrix_split_2x2_numeric(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !<
    !< this routine computes A_II, A_IG, A_GI and A_GG given the global 
    !< matrix A. Note that A_II, A_IG, A_GI and A_GG are all optional.
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        class(base_sparse_matrix_t),           intent(inout) :: A_II
        class(base_sparse_matrix_t),           intent(inout) :: A_IG
        class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
        class(base_sparse_matrix_t),           intent(inout) :: A_GG
        integer(ip)                                          :: i
        integer(ip)                                          :: nz
        integer(ip)                                          :: total_cols
        integer(ip)                                          :: total_rows
    !-----------------------------------------------------------------
        assert(this%get_symmetric_storage() .eqv. A_II%get_symmetric_storage()) 
        assert(this%get_symmetric_storage() .eqv. A_GG%get_symmetric_storage())
        total_rows = this%get_num_rows()
        total_cols = this%get_num_cols()
        assert((A_II%get_num_rows()==num_row) .and. (A_II%get_num_cols()==num_col))
        assert((A_IG%get_num_rows()==num_row) .and. (A_IG%get_num_cols()==total_cols-num_col))
        assert((A_GG%get_num_rows()==total_rows-num_row) .and. (A_GG%get_num_cols()==total_cols-num_col))
        if(present(A_GI)) then
            assert((A_GI%get_num_rows()==this%get_num_rows()-num_row) .and. (A_GI%get_num_cols()==num_col))
        endif

        ! algorithm for split_2x2 a non sorted COO matrix. Not tested yet!
        ! TODO: algorithm for split_2x2 a sorted by rows COO matrix
        ! TODO: algorithm for split_2x2 a sorted by cols COO matrix
        if(this%is_symbolic()) then
            call this%split_2x2_symbolic(num_row, num_col, A_II, A_IG, A_GI, A_GG)
        else
            do i=1, this%get_nnz()
                if(this%ia(i)<=num_row) then
                    if(this%ja(i)<=num_col) then
                        call A_II%insert(ia  = this%ia(i),   &
                                         ja  = this%ja(i),   &
                                         val = this%val(i),  &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_IG%insert(ia  = this%ia(i),         &
                                         ja  = this%ja(i)-num_col, &
                                         val = this%val(i),        &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                else
                    if(this%ja(i)<=num_col) then
                        if(present(A_GI)) call A_GI%insert(ia  = this%ia(i)-num_row, &
                                                           ja  = this%ja(i),         &
                                                           imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    else
                        call A_GG%insert(ia  = this%ia(i)-num_row, &
                                         ja  = this%ja(i)-num_col, &
                                         val = this%val(i),        &
                                         imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                    endif
                endif
            enddo
        endif
    end subroutine coo_sparse_matrix_split_2x2_numeric



    subroutine coo_sparse_matrix_split_2x2_symbolic(this, num_row, num_col, A_II, A_IG, A_GI, A_GG) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2
    !< A = [A_II A_IG]
    !<     [A_GI A_GG]
    !<
    !< this routine computes A_II, A_IG, A_GI and A_GG given the global 
    !< matrix A. Note that A_II, A_IG, A_GI and A_GG are all optional.
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        class(base_sparse_matrix_t),           intent(inout) :: A_II
        class(base_sparse_matrix_t),           intent(inout) :: A_IG
        class(base_sparse_matrix_t), optional, intent(inout) :: A_GI
        class(base_sparse_matrix_t),           intent(inout) :: A_GG
        integer(ip)                                          :: i
        integer(ip)                                          :: nz
        integer(ip)                                          :: total_cols
        integer(ip)                                          :: total_rows
    !-----------------------------------------------------------------
        assert(this%get_symmetric_storage() .eqv. A_II%get_symmetric_storage()) 
        assert(this%get_symmetric_storage() .eqv. A_GG%get_symmetric_storage())
        total_rows = this%get_num_rows()
        total_cols = this%get_num_cols()
        assert((A_II%get_num_rows()==num_row) .and. (A_II%get_num_cols()==num_col))
        assert((A_IG%get_num_rows()==num_row) .and. (A_IG%get_num_cols()==total_cols-num_col))
        assert((A_GG%get_num_rows()==total_rows-num_row) .and. (A_GG%get_num_cols()==total_cols-num_col))
        if(present(A_GI)) then
            assert((A_GI%get_num_rows()==this%get_num_rows()-num_row) .and. (A_GI%get_num_cols()==num_col))
        endif

        ! algorithm for split_2x2 a non sorted COO matrix. Not tested yet!
        ! TODO: algorithm for split_2x2 a sorted by rows COO matrix
        ! TODO: algorithm for split_2x2 a sorted by cols COO matrix
        do i=1, this%get_nnz()
            if(this%ia(i)<=num_row) then
                if(this%ja(i)<=num_col) then
                    call A_II%insert(ia  = this%ia(i),   &
                                     ja  = this%ja(i),   &
                                     imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                else
                    call A_IG%insert(ia  = this%ia(i),         &
                                     ja  = this%ja(i)-num_col, &
                                     imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                endif
            else
                if(this%ja(i)<=num_col) then
                    if(present(A_GI)) call A_GI%insert(ia  = this%ia(i)-num_row, &
                                                       ja  = this%ja(i),         &
                                                       imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                else
                    call A_GG%insert(ia  = this%ia(i)-num_row, &
                                     ja  = this%ja(i)-num_col, &
                                     imin = 1, imax = num_row, jmin = 1, jmax = num_col)
                endif
            endif
        enddo
    end subroutine coo_sparse_matrix_split_2x2_symbolic


    subroutine coo_sparse_matrix_permute_and_split_2x2_numeric(this, num_row, num_col, perm, iperm, A_CC, A_CR, A_RC, A_RR) 
    !-----------------------------------------------------------------
    !< Split matrix in 2x2 and permute some columns and rows
    !< given 2 permutation arrays (perm and iperm)
    !< 
    !< A = [A_CC A_RC]
    !<     [A_CR A_RR]
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        real(rp),    allocatable,              intent(out)   :: A_CC(:,:)
        real(rp),    allocatable,              intent(out)   :: A_CR(:,:)
        real(rp),    allocatable,              intent(out)   :: A_RC(:,:)
        class(base_sparse_matrix_t),           intent(inout) :: A_RR
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine coo_sparse_matrix_permute_and_split_2x2_numeric


    subroutine coo_sparse_matrix_permute_and_split_2x2_symbolic(this, num_row, num_col, perm, iperm, A_RR) 
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
        class(coo_sparse_matrix_t),            intent(in)    :: this
        integer(ip),                           intent(in)    :: num_row
        integer(ip),                           intent(in)    :: num_col
        integer(ip),                           intent(in)    :: perm(:)
        integer(ip),                           intent(in)    :: iperm(:)
        class(base_sparse_matrix_t),           intent(inout) :: A_RR
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine coo_sparse_matrix_permute_and_split_2x2_symbolic


    subroutine coo_sparse_matrix_expand_matrix_numeric(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, C_T_val, I_nz, I_ia, I_ja, I_val, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input) sparse matrix must be in ASSEMBLED state
    !<  - TO (output) sparse matrix must be in START state
    !<  - C_T coordinate arrays (C_T_ia, C_T_ja and C_T_val) must 
    !<    have the same size (C_T_nz)
    !<  - I coordinate arrays (I_ia, I_ja and I_val) must 
    !<    have the same size (I_nz)
    !<  - Row index arrays (X_ia) must be in ascendent order
    !<  - Column index arrays (X_ja) must be in ascendent order for 
    !<    each row
    !<  - For each C_T row index (C_T_ia): 1<=C_T_ia(i)<=this%get_num_rows()
    !<  - For each C_T column index (C_T_ja): 1<=C_T_ja(i)<=C_T_num_cols
    !<  - For each I row and column index (I_ia and I_ja): 
    !<    1<=I_ia(i) and I_ia(i)<=C_T_num_cols
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),      intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        real(rp),                        intent(in)    :: C_T_val(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        real(rp),                        intent(in)    :: I_val(I_nz)
        class(base_sparse_matrix_t),     intent(inout) :: to
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine coo_sparse_matrix_expand_matrix_numeric


    subroutine coo_sparse_matrix_expand_matrix_symbolic(this, C_T_num_cols, C_T_nz, C_T_ia, C_T_ja, I_nz, I_ia, I_ja, to)
    !-----------------------------------------------------------------
    !< Expand matrix A given a (by_row) sorted C_T and I in COO
    !< A = [A C_T]
    !<     [C  I ]
    !< Some considerations:
    !<  - C = transpose(C_T)
    !<  - I is a square matrix
    !<  - THIS (input) sparse matrix must be in ASSEMBLED or 
    !<    ASSEMBLED_SYMBOLIC state
    !<  - TO (output) sparse matrix must be in START state
    !<  - C_T coordinate arrays (C_T_ia, C_T_ja and C_T_val) must 
    !<    have the same size (C_T_nz)
    !<  - I coordinate arrays (I_ia, I_ja and I_val) must 
    !<    have the same size (I_nz)
    !<  - Row index arrays (X_ia) must be in ascendent order
    !<  - Column index arrays (X_ja) must be in ascendent order for 
    !<    each row
    !<  - For each C_T row index (C_T_ia): 1<=C_T_ia(i)<=this%get_num_rows()
    !<  - For each C_T column index (C_T_ja): 1<=C_T_ja(i)<=C_T_num_cols
    !<  - For each I row and column index (I_ia and I_ja): 
    !<    1<=I_ia(i) and I_ia(i)<=C_T_num_cols
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t),      intent(in)    :: this
        integer,                         intent(in)    :: C_T_num_cols
        integer,                         intent(in)    :: C_T_nz
        integer(ip),                     intent(in)    :: C_T_ia(C_T_nz)
        integer(ip),                     intent(in)    :: C_T_ja(C_T_nz)
        integer,                         intent(in)    :: I_nz
        integer(ip),                     intent(in)    :: I_ia(I_nz)
        integer(ip),                     intent(in)    :: I_ja(I_nz)
        class(base_sparse_matrix_t),     intent(inout) :: to
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine coo_sparse_matrix_expand_matrix_symbolic


    subroutine coo_sparse_matrix_extract_diagonal(this, diagonal)
    !-----------------------------------------------------------------
    !< Return the diagonal of a coo sparse matrix
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in)    :: this
        real(rp), allocatable,      intent(inout) :: diagonal(:)
    !-----------------------------------------------------------------
        check(.false.)
    end subroutine coo_sparse_matrix_extract_diagonal


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
        if(.not. this%is_symbolic() .and. .not. this%state_is_created()) then
            call to%allocate_values_body(nnz)
            to%val(1:nnz) = this%val(1:nnz)
        endif
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
        if(.not. from%is_symbolic() .and. .not. from%state_is_created()) then
            call this%allocate_values_body(nnz)
            this%val(1:nnz) = from%val(1:nnz)
        endif
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
        if(.not. this%is_symbolic()) call move_alloc(from=this%val, to=to%val)
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
        if(.not. from%is_symbolic()) call move_alloc(from=from%val, to=this%val)
        call from%free()
    end subroutine coo_sparse_matrix_move_from_coo


    subroutine coo_sparse_matrix_move_to_fmt(this, to)
    !-----------------------------------------------------------------
    !< Move this (COO) -> to (FMT)
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


    subroutine coo_sparse_matrix_free_coords(this)
    !-----------------------------------------------------------------
    !< Clean coords of COO sparse matrix format derived type
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
    !-----------------------------------------------------------------
        if(allocated(this%ia))  call memfree (this%ia, __FILE__, __LINE__)
        if(allocated(this%ja))  call memfree (this%ja, __FILE__, __LINE__)
    end subroutine coo_sparse_matrix_free_coords


    subroutine coo_sparse_matrix_free_val(this)
    !-----------------------------------------------------------------
    !< Clean values of COO sparse matrix format derived type
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(inout)  :: this
    !-----------------------------------------------------------------
        if(allocated(this%val))  call memfree (this%val, __FILE__, __LINE__)
    end subroutine coo_sparse_matrix_free_val


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


    subroutine coo_sparse_matrix_print_matrix_market_body (this, lunou, ng, l2g)
    !-----------------------------------------------------------------
    !< Print a COO matrix to matrix market format
    !-----------------------------------------------------------------
        class(coo_sparse_matrix_t), intent(in) :: this
        integer(ip),                intent(in) :: lunou
        integer(ip), optional,      intent(in) :: ng
        integer(ip), optional,      intent(in) :: l2g (*)
        integer(ip) :: i, j
        integer(ip) :: nr, nc, nnz
    !-----------------------------------------------------------------
        if ( present(ng) ) then 
            nr = ng
            nc = ng
        else
            nr = this%num_rows
            nc = this%num_cols
        end if
        nnz = this%nnz

        write (lunou,'(a)') '%%MatrixMarket matrix coordinate real general'
        if (.not. this%get_symmetric_storage()) then
            write (lunou,*) nr,nc,nnz
            do i=1,this%get_nnz()
                if (present(l2g)) then
                    write(lunou,'(i12, i12, e32.25)') l2g(this%ia(i)), l2g(this%ja(i)), this%val(i)
                else
                    write(lunou,'(i12, i12, e32.25)') this%ja(i), this%ja(i), this%val(i)
                end if
            end do
        else 
            if(nnz>0) nnz = 2*(nnz) - this%num_rows
            write (lunou,*) nr,nc,nnz

            do i=1,this%get_nnz()
                if (present(l2g)) then
                    write(lunou,'(i12, i12, e32.25)') l2g(this%ia(i)), l2g(this%ja(i)), this%val(i)
                else
                    write(lunou,'(i12, i12, e32.25)') this%ia(i), this%ja(i), this%val(i)
                end if
                if (this%ia(i) /= this%ja(i)) then
                    if (present(l2g)) then
                        write(lunou,'(i12, i12, e32.25)') l2g(this%ja(j)), l2g(this%ia(i)), this%val(i)
                    else
                        write(lunou,'(i12, i12, e32.25)') this%ja(i), this%ia(i), this%val(i)
                    end if
                end if
            end do
        end if
    end subroutine coo_sparse_matrix_print_matrix_market_body

!---------------------------------------------------------------------
!< AUX PROCEDURES
!---------------------------------------------------------------------

    subroutine assign_value(input, output)
    !-----------------------------------------------------------------
    ! Assign an input value to the autput
    !-----------------------------------------------------------------
        real(rp), intent(in)    :: input
        real(rp), intent(inout) :: output
    !-------------------------------------------------------------
        output = input
    end subroutine assign_value


    subroutine sum_value(input, output)
    !-----------------------------------------------------------------
    ! Sum an input value to the autput
    !-----------------------------------------------------------------
        real(rp), intent(in)    :: input
        real(rp), intent(inout) :: output
    !-----------------------------------------------------------------
        output = output + input
    end subroutine sum_value


    function  binary_search(key,size,vector) result(ipos)
    !-----------------------------------------------------------------
    ! Perform binary search in a vector
    !-----------------------------------------------------------------
        integer(ip), intent(in) :: key
        integer(ip), intent(in) :: size
        integer(ip), intent(in) :: vector(size)
        integer(ip)             :: ipos
        integer(ip)             :: lowerbound
        integer(ip)             :: upperbound
        integer(ip)             :: midpoint
    !-----------------------------------------------------------------
        lowerbound = 1 
        upperbound = size
        ipos = -1 

        do while (lowerbound.le.upperbound) 
            midpoint = (lowerbound+upperbound)/2
            if (key.eq.vector(midpoint))  then
                ipos = midpoint
                lowerbound  = upperbound + 1
            else if (key < vector(midpoint))  then
                upperbound = midpoint-1
            else 
                lowerbound = midpoint + 1
            end if
        enddo
        return
    end function binary_search

    subroutine mergesort_link_list(n,k,l,iret)
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
        integer(ip), intent(in)  :: n
        integer(ip), intent(in)  :: k(n)
        integer(ip), intent(out) :: l(0:n+1)
        integer(ip), intent(out) :: iret
        integer(ip)              :: p,q,s,t
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

    end subroutine mergesort_link_list


    subroutine reorder_numeric_coo_from_link_list(n,x,i1,i2,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        real(rp),    intent(inout) :: x(*)
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: i2(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1, isw2
        real(rp)    :: swap
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
    end subroutine reorder_numeric_coo_from_link_list


    subroutine reorder_symbolic_coo_from_link_list(n,i1,i2,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: i2(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1, isw2
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
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
    end subroutine reorder_symbolic_coo_from_link_list
    

    subroutine reorder_ip_from_link_list(n,i1,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1
    !-------------------------------------------------------------
        lp = iaux(0)
        k  = 1
        do 
            if ((lp == 0).or.(k>n)) exit
            do 
                if (lp >= k) exit
                lp = iaux(lp)
            end do
            isw1     = i1(lp)
            i1(lp)   = i1(k)
            i1(k)    = isw1
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_ip_from_link_list

    
    subroutine reorder_ip_rp_from_link_list(n,x,i1,iaux)
    !-------------------------------------------------------------
    !  Reorder (an) input vector(s) based on a list sort output.
    !  Based on: D. E. Knuth: The Art of Computer Programming
    !            vol. 3: Sorting and Searching, Addison Wesley, 1973
    !            ex. 5.2.12
    !-------------------------------------------------------------
        use types_names
        integer(ip), intent(in)    :: n
        real(rp),    intent(inout) :: x(*)
        integer(ip), intent(inout) :: i1(*)
        integer(ip), intent(inout) :: iaux(0:*) 
        integer(ip) :: lswap, lp, k, isw1
        real(rp)    :: swap
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
            lswap    = iaux(lp)
            iaux(lp) = iaux(k)
            iaux(k)  = lp
            lp = lswap 
            k  = k + 1
        enddo
        return
    end subroutine reorder_ip_rp_from_link_list


end module base_sparse_matrix_names
