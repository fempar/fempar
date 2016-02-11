program test_coo_append

USE types_names
USE memor_names
USE serial_names
USE sparse_matrix_names
USE csr_sparse_matrix_names


implicit none

# include "debug.i90"

    type(sparse_matrix_t)       :: sparse_matrix
    integer(ip), allocatable    :: perm(:)
    integer(ip), allocatable    :: iperm(:)
    real(rp),    allocatable    :: A_CC(:,:)
    real(rp),    allocatable    :: A_CR(:,:)
    real(rp),    allocatable    :: A_RC(:,:)
    type(sparse_matrix_t)       :: A_RR

    call meminit()

!------------------------------------------------------------------
! NON SYMMETRIC STORAGE
!------------------------------------------------------------------

    call memalloc(4, perm, __FILE__, __LINE__)
    call memalloc(4, iperm, __FILE__, __LINE__)

    perm  = (/4,2,1,3/)
    iperm = (/3,2,4,1/)

    call sparse_matrix%create(4, .false., .false., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%insert(nz=8,                    &
                              ia=(/1,1,2,3,3,3,4,4/),  &
                              ja=(/1,3,2,1,3,4,3,4/),  &
                              val=(/1.,2.,3.,4.,5.,6.,7.,8./))

    call sparse_matrix%convert('CSR')

    call sparse_matrix%permute_and_split_2x2_symbolic(num_row=2, num_col=2, perm=perm, iperm=iperm, A_RR=A_RR)
    call sparse_matrix%permute_and_split_2x2_numeric(num_row=2, num_col=2, perm=perm, iperm=iperm, A_CC=A_CC, A_CR=A_CR, A_RC=A_RC, A_RR=A_RR)

    print*, '--------------------------------------------'
    print*, ' Original matrix (NON SYMMETRIC STORAGE)'
    print*, ' Submatrices (NON SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'
    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_CC'
    print*, '--------------------------------------------'
    print*, A_CC
    print*, '--------------------------------------------'
    print*, ' A_CR'
    print*, '--------------------------------------------'
    print*, A_CR
    print*, '--------------------------------------------'
    print*, ' A_RC'
    print*, '--------------------------------------------'
    print*, A_RC
    print*, '--------------------------------------------'
    print*, ' A_RR'
    print*, '--------------------------------------------'
    call A_RR%print_matrix_market(6)

    call A_RR%Free()
    call memfree(A_CC, __FILE__, __LINE__)
    call memfree(A_CR, __FILE__, __LINE__)
    call memfree(A_RC, __FILE__, __LINE__)

!------------------------------------------------------------------
! NON SYMMETRIC STORAGE
!------------------------------------------------------------------

    call sparse_matrix%permute_and_split_2x2_symbolic(num_row=2, num_col=2, perm=perm, iperm=iperm, A_RR=A_RR, symmetric_storage=.true., symmetric=.true.)
    call sparse_matrix%permute_and_split_2x2_numeric(num_row=2, num_col=2, perm=perm, iperm=iperm, A_CC=A_CC, A_CR=A_CR, A_RC=A_RC, A_RR=A_RR, symmetric_storage=.true., symmetric=.true.)

    print*, '--------------------------------------------'
    print*, ' Original matrix (NON SYMMETRIC STORAGE)'
    print*, ' Submatrices (SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'
    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_CC'
    print*, '--------------------------------------------'
    print*, A_CC
    print*, '--------------------------------------------'
    print*, ' A_CR'
    print*, '--------------------------------------------'
    print*, A_CR
    print*, '--------------------------------------------'
    print*, ' A_RC'
    print*, '--------------------------------------------'
    print*, A_RC
    print*, '--------------------------------------------'
    print*, ' A_RR'
    print*, '--------------------------------------------'
    call A_RR%print_matrix_market(6)

    call A_RR%Free()
    call sparse_matrix%free()
    call memfree(A_CC, __FILE__, __LINE__)
    call memfree(A_CR, __FILE__, __LINE__)
    call memfree(A_RC, __FILE__, __LINE__)

!------------------------------------------------------------------
! SYMMETRIC STORAGE
!------------------------------------------------------------------

    call sparse_matrix%create(4, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%insert(nz=8,                    &
                              ia=(/1,1,2,3,3,3,4,4/),  &
                              ja=(/1,4,2,1,3,4,3,4/),  &
                              val=(/1.,2.,3.,4.,5.,6.,7.,8./))

    call sparse_matrix%convert('CSR')

    call sparse_matrix%permute_and_split_2x2_symbolic(num_row=2, num_col=2, perm=perm, iperm=iperm, A_RR=A_RR)
!    call sparse_matrix%permute_and_split_2x2_numeric(num_row=2, num_col=2, perm=perm, iperm=iperm, A_CC=A_CC, A_CR=A_CR, A_RC=A_RC, A_RR=A_RR)

    print*, '--------------------------------------------'
    print*, ' Original matrix (SYMMETRIC STORAGE)'
    print*, ' Submatrices (SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'
    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_CC'
    print*, '--------------------------------------------'
    print*, A_CC
    print*, '--------------------------------------------'
    print*, ' A_CR'
    print*, '--------------------------------------------'
    print*, A_CR
    print*, '--------------------------------------------'
    print*, ' A_RC'
    print*, '--------------------------------------------'
    print*, A_RC
    print*, '--------------------------------------------'
    print*, ' A_RR'
    print*, '--------------------------------------------'
    call A_RR%print(6)

    call A_RR%Free()
!    call memfree(A_CC, __FILE__, __LINE__)
!    call memfree(A_CR, __FILE__, __LINE__)
!    call memfree(A_RC, __FILE__, __LINE__)

!------------------------------------------------------------------
! SYMMETRIC STORAGE
!------------------------------------------------------------------

    call sparse_matrix%permute_and_split_2x2_symbolic(num_row=2, num_col=2, perm=perm, iperm=iperm, A_RR=A_RR, symmetric_storage=.false.)
    call sparse_matrix%permute_and_split_2x2_numeric(num_row=2, num_col=2, perm=perm, iperm=iperm, A_CC=A_CC, A_CR=A_CR, A_RC=A_RC, A_RR=A_RR, symmetric_storage=.false.)

    print*, '--------------------------------------------'
    print*, ' Original matrix (SYMMETRIC STORAGE)'
    print*, ' Submatrices (NON SYMMETRIC STORAGE)'
    print*, '--------------------------------------------'
    call sparse_matrix%print_matrix_market(6)
    print*, '--------------------------------------------'
    print*, ' A_CC'
    print*, '--------------------------------------------'
    print*, A_CC
    print*, '--------------------------------------------'
    print*, ' A_CR'
    print*, '--------------------------------------------'
    print*, A_CR
    print*, '--------------------------------------------'
    print*, ' A_RC'
    print*, '--------------------------------------------'
    print*, A_RC
    print*, '--------------------------------------------'
    print*, ' A_RR'
    print*, '--------------------------------------------'
    call A_RR%print_matrix_market(6)

    call A_RR%Free()

    call sparse_matrix%free()

    if(allocated(A_CC)) call memfree(A_CC, __FILE__, __LINE__)
    if(allocated(A_CR)) call memfree(A_CR, __FILE__, __LINE__)
    if(allocated(A_RC)) call memfree(A_RC, __FILE__, __LINE__)
    call memfree(perm, __FILE__, __LINE__)
    call memfree(iperm,__FILE__, __LINE__)
    call memstatus()

end program
