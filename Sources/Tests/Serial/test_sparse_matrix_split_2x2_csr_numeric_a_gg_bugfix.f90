program test_coo_append

USE types_names
USE memor_names
USE fempar_names
USE sparse_matrix_names
USE csr_sparse_matrix_names


implicit none

# include "debug.i90"

    type(sparse_matrix_t)       :: sparse_matrix
    type(sparse_matrix_t)       :: A_II, A_IG, A_GI, A_GG

    call meminit()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FULL MATRIX (NNZ>0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------
! SYMMETRIC STORAGE
!------------------------------------------------------------------

    call sparse_matrix%create(9, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%insert(nz=6,          &
                              ia=1,          &
                              ja=(/1,2,3,4,6,9/),  &
                              val=(/1.,2.,3.,4.,6.,9./))

    call sparse_matrix%insert(nz=8,          &
                              ia=2,          &
                              ja=(/2,3,4,5,6,7,8,9/),  &
                              val=(/2.,3.,4.,5.,6.,7.,8.,9./))

    call sparse_matrix%insert(nz=2,          &
                              ia=3,          &
                              ja=(/3,4/),  &
                              val=(/3.,4./))

    call sparse_matrix%insert(nz=4,          &
                              ia=4,          &
                              ja=(/4,5,8/),  &
                              val=(/4.,5.,8./))

    call sparse_matrix%insert(nz=2,          &
                              ia=5,          &
                              ja=(/5,8/),  &
                              val=(/5.,8./))

    call sparse_matrix%insert(nz=2,          &
                              ia=6,          &
                              ja=(/6,9/),  &
                              val=(/6.,9./))


    call sparse_matrix%insert(nz=3,          &
                              ia=7,          &
                              ja=(/7,8,9/),  &
                              val=(/7.,8.,9./))

    call sparse_matrix%insert(nz=2,          &
                              ia=8,          &
                              ja=(/8,9/),  &
                              val=(/8.,9./))

    call sparse_matrix%insert(nz=1,          &
                              ia=9,          &
                              ja=(/9/),  &
                              val=(/9./))

    call sparse_matrix%convert('CSR')

    call sparse_matrix%split_2x2_symbolic(num_row=4, num_col=4, A_II=A_II, A_IG=A_IG, A_GG=A_GG)
    call sparse_matrix%split_2x2_numeric(num_row=4, num_col=4, A_II=A_II, A_IG=A_IG, A_GG=A_GG)

    print*, '--------------------------------------------'
    print*, ' Original matrix (SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'
    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_II'
    print*, '--------------------------------------------'
    call A_II%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_IG'
    print*, '--------------------------------------------'
    call A_IG%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_GG'
    print*, '--------------------------------------------'
    call A_GG%print_matrix_market(6)

    call A_II%Free()
    call A_IG%Free()
    call A_GG%Free()

!------------------------------------------------------------------
! SYMMETRIC STORAGE - NON SYMMETRIC STORAGE FOR SUBMATRICES
!------------------------------------------------------------------

    call sparse_matrix%split_2x2_symbolic(num_row=4, num_col=4, A_II=A_II, A_IG=A_IG, A_GG=A_GG, symmetric_storage=.false.)
    call sparse_matrix%split_2x2_numeric(num_row=4, num_col=4, A_II=A_II, A_IG=A_IG, A_GG=A_GG, symmetric_storage=.false.)

    print*, '--------------------------------------------'
    print*, ' Original matrix (SYMMETRIC STORAGE)'
    print*, ' Submatrices (NON SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'

    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_II'
    print*, '--------------------------------------------'
    call A_II%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_IG'
    print*, '--------------------------------------------'
    call A_IG%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_GG'
    print*, '--------------------------------------------'
    call A_GG%print_matrix_market(6)

    call A_II%Free()
    call A_IG%Free()
    call A_GG%Free()

    call sparse_matrix%free_in_stages(free_clean)

    call memstatus()

end program
