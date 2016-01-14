program test_coo_sort_and_compress

USE types_names
USE memor_names
USE base_sparse_matrix_names

implicit none

# include "debug.i90"

    type(coo_sparse_matrix_t) :: coo_matrix

    call meminit()

!------------------------------------------------------------------
! NUMERIC
!------------------------------------------------------------------

    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_START)

    call coo_matrix%create(num_rows_and_cols=5,             &
                           symmetric_storage=.true.,        &
                           is_symmetric=.true.,             &
                           sign=SPARSE_MATRIX_SIGN_UNKNOWN, &
                           nz=12)
    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_CREATED)

    call coo_matrix%insert(nz=12,                                          &
                           ia=(/1,2,3,4,5,1,6,1,2,3,4,5/),                 &
                           ja=(/1,2,3,4,5,6,1,5,4,3,2,1/),                 &
                           val=(/1.,2.,3.,4.,5.,99.,99.,5.,4.,3.,2.,1./))
    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_BUILD_NUMERIC)

    call coo_matrix%print( 6 )

    call coo_matrix%sort_and_compress(by_cols=.false.)

    call coo_matrix%print( 6 )

    call coo_matrix%sort_and_compress(by_cols=.true.)

    call coo_matrix%print( 6 )

    call coo_matrix%free()

    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_START)

!------------------------------------------------------------------
! SYMBOLIC
!------------------------------------------------------------------

    call coo_matrix%create(num_rows_and_cols=5,             &
                           symmetric_storage=.true.,        &
                           is_symmetric=.true.,             &
                           sign=SPARSE_MATRIX_SIGN_UNKNOWN, &
                           nz=12)
    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_CREATED)

    call coo_matrix%insert(nz=12,                                 &
                           ia=(/1,2,3,4,5,1,6,1,2,3,4,5/),        &
                           ja=(/1,2,3,4,5,6,1,5,4,3,2,1/))
    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_BUILD_SYMBOLIC)

    call coo_matrix%print( 6 )

    call coo_matrix%sort_and_compress(by_cols=.false.)

    call coo_matrix%print( 6 )

    call coo_matrix%sort_and_compress(by_cols=.true.)

    call coo_matrix%print( 6 )

    call coo_matrix%free()

    check(coo_matrix%get_state() == SPARSE_MATRIX_STATE_START)


    call memstatus()

end program
