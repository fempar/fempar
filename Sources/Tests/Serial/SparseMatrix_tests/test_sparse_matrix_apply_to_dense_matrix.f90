program test_sparse_matrix

USE types_names
USE memor_names
USE serial_names
USE sparse_matrix_names
USE csr_sparse_matrix_names


implicit none

# include "debug.i90"

    type(sparse_matrix_t)       :: sparse_matrix
    integer                     :: n        = 5
    real                        :: alpha    = 1.0
    integer                     :: LDB      = 5
    real                        :: b(5,5)   = 1.0
    real                        :: beta     = 0.0
    integer                     :: LDC      = 5
    real                        :: c(5,5)   = 0.0
    real                        :: val(4,4) = 9

    call meminit()


!------------------------------------------------------------------
! GENERAL - EMPTY MATRIX (NNZ==0)
!------------------------------------------------------------------
    ! Create: START=>CREATE
    call sparse_matrix%create(num_rows=5,num_cols=5)

    ! Convert: CREATE=>ASSEMBLED
    call sparse_matrix%convert(csr_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
#ifdef ENABLE_MKL
    ! Convert: ASSEMBLED=>ASSEMBLED
    call sparse_matrix%convert(coo_format)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
#endif

    ! Free: CREATED=>START
    call sparse_matrix%free_in_stages(free_clean)

!------------------------------------------------------------------
! GENERAL - NUMERIC
!------------------------------------------------------------------
    ! Create: START=>CREATE
    call sparse_matrix%create(num_rows=5,num_cols=5)

    ! Append: CREATE=>BUILD_NUMERIC
    call sparse_matrix%insert(nz=10,                                 &
                              ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                              ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                              val=(/1.,2.,3.,2.,1.,5.,4.,3.,4.,5./))

    ! Convert: BUILD_NUMERIC=>ASSEMBLED
    call sparse_matrix%convert(csr_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
    c = 0.0
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
    c = 0.0
#ifdef ENABLE_MKL
    ! Convert: ASSEMBLED=>ASSEMBLED
    call sparse_matrix%convert(coo_format)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
#endif

    ! Free: CREATED=>START
    call sparse_matrix%free_in_stages(free_clean)


!------------------------------------------------------------------
! SYMMETRIC STORGE - EMPTY MATRIX (NNZ==0)
!------------------------------------------------------------------
    ! Create: START=>CREATE
    call sparse_matrix%create(5, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    ! Convert: CREATE=>ASSEMBLED
    call sparse_matrix%convert(csr_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))

#ifdef ENABLE_MKL
    ! Convert: ASSEMBLED=>ASSEMBLED
    call sparse_matrix%convert(coo_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 0.0))
#endif

    ! Free: CREATED=>START
    call sparse_matrix%free_in_stages(free_clean)

!------------------------------------------------------------------
! SYMMETRIC STORAGE - NUMERIC
!------------------------------------------------------------------
    ! Create: START=>CREATE
    call sparse_matrix%create(5, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    ! Append: CREATE=>BUILD_NUMERIC
    call sparse_matrix%insert(nz=7,                                 &
                              ia=(/1,2,3,4,5,1,2/),            &
                              ja=(/1,2,3,4,5,5,4/),            &
                              val=(/1.,2.,6.,2.,1.,5.,4./))

    ! Convert: BUILD_NUMERIC=>ASSEMBLED
    call sparse_matrix%convert(csr_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
    c = 0.0
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) s
    check(all(c == 6.0))
    c = 0.0

#ifdef ENABLE_MKL
    ! Convert: BUILD_NUMERIC=>ASSEMBLED
    call sparse_matrix%convert(coo_format)
    call sparse_matrix%print( 6)

    ! Apply to dense matrix
    call sparse_matrix%apply_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
    ! Apply transpose to dense matrix
    call sparse_matrix%apply_transpose_to_dense_matrix(n, alpha, LDB, b, beta, LDC, c) 
    check(all(c == 6.0))
#endif

    ! Free: CREATED=>START
    call sparse_matrix%free_in_stages(free_clean)

    call memstatus()

end program
