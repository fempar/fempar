program test_coo_csr_copy_move

USE types_names
USE memor_names
USE base_sparse_matrix_names
USE csr_sparse_matrix_names

implicit none

# include "debug.i90"

    type(coo_sparse_matrix_t) :: coo_matrix
    type(csr_sparse_matrix_t) :: csr_matrix

    call meminit()


!------------------------------------------------------------------
! EMPTY MATRIX (NNZ==0)
!------------------------------------------------------------------

    check(coo_matrix%state_is_start())

    call coo_matrix%create(num_rows_and_cols = 5,                          &
                           symmetric_storage = .true.,                     &
                           is_symmetric      = .true.,                     &
                           sign              = SPARSE_MATRIX_SIGN_UNKNOWN, &
                           nz                = 0)

    check(coo_matrix%state_is_created())

    call coo_matrix%sort_and_compress()

    ! Copy coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%copy_to_fmt(csr_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_fmt(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) <- coo_matrix (COO)
    call csr_matrix%copy_from_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)


    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%move_to_fmt(csr_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- csr_matrix_copy (CSR)
    call coo_matrix%move_from_fmt(csr_matrix)
    check(csr_matrix%state_is_start())
    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call csr_matrix%move_from_coo(coo_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (CSR)
    call csr_matrix%move_to_coo(coo_matrix)
    check(csr_matrix%state_is_start())

    call coo_matrix%free()

    check(coo_matrix%state_is_start())

!------------------------------------------------------------------
! NUMERIC
!------------------------------------------------------------------

    check(coo_matrix%state_is_start())

    call coo_matrix%create(num_rows_and_cols = 5,                          &
                           symmetric_storage = .true.,                     &
                           is_symmetric      = .true.,                     &
                           sign              = SPARSE_MATRIX_SIGN_UNKNOWN, &
                           nz                = 10)

    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                           ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                           ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                           val=(/1.,2.,3.,4.,5.,5.,4.,3.,2.,1./), &
                           imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_numeric())

    call coo_matrix%sort_and_compress()

    ! Copy coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%copy_to_fmt(csr_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_fmt(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) <- coo_matrix (COO)
    call csr_matrix%copy_from_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)


    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%move_to_fmt(csr_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- csr_matrix_copy (CSR)
    call coo_matrix%move_from_fmt(csr_matrix)
    check(csr_matrix%state_is_start())
    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call csr_matrix%move_from_coo(coo_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (CSR)
    call csr_matrix%move_to_coo(coo_matrix)
    check(csr_matrix%state_is_start())

    call coo_matrix%free()

    check(coo_matrix%state_is_start())

!------------------------------------------------------------------
! SYMBOLIC
!------------------------------------------------------------------

    call coo_matrix%create(num_rows_and_cols = 5,                          &
                           symmetric_storage = .true.,                     &
                           is_symmetric      = .true.,                     &
                           sign              = SPARSE_MATRIX_SIGN_UNKNOWN, &
                           nz                = 10)

    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                           ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                           ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                           imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_symbolic())

    call coo_matrix%sort_and_compress()

    ! Copy coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%copy_to_fmt(csr_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_fmt(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) <- coo_matrix (COO)
    call csr_matrix%copy_from_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)
    ! Copy csr_matrix (CSR) -> coo_matrix (COO)
    call csr_matrix%copy_to_coo(coo_matrix)
    call compare_coo_csr_matrix(coo_matrix, csr_matrix)


    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call coo_matrix%move_to_fmt(csr_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- csr_matrix_copy (CSR)
    call coo_matrix%move_from_fmt(csr_matrix)
    check(csr_matrix%state_is_start())
    ! Move coo_matrix (COO) -> csr_matrix (CSR)
    call csr_matrix%move_from_coo(coo_matrix)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (CSR)
    call csr_matrix%move_to_coo(coo_matrix)
    check(csr_matrix%state_is_start())

    call coo_matrix%free()

    check(coo_matrix%state_is_start())



    call memstatus()

contains

    subroutine compare_coo_csr_matrix(a, b)
        class(coo_sparse_matrix_t), intent(in) :: a
        class(csr_sparse_matrix_t), intent(in) :: b
    
        check(a%get_num_rows()            == b%get_num_rows())
        check(a%get_num_cols()            == b%get_num_cols())
        check(a%get_nnz()                 == b%get_nnz())
        check(a%is_symmetric()           .eqv. b%is_symmetric())
        check(a%get_symmetric_storage()  .eqv. b%get_symmetric_storage())
        check(a%get_sign()                == b%get_sign())
        check(a%is_by_rows()             .eqv. b%is_by_rows())
        check(a%get_state()               == b%get_state())
        check(a%get_nnz()+1               == b%irp(a%get_num_rows()+1))
    end subroutine compare_coo_csr_matrix

end program
