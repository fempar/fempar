program test_coo_copy_move

USE types_names
USE memor_names
USE base_sparse_matrix_names

implicit none

# include "debug.i90"

    type(coo_sparse_matrix_t) :: coo_matrix
    type(coo_sparse_matrix_t) :: coo_matrix_copy

    call meminit()

    check(coo_matrix%state_is_start())

!------------------------------------------------------------------
! EMPTY MATRIX (NNZ==0)
!------------------------------------------------------------------

    call coo_matrix%create(num_rows=0,num_cols=0, nz=0)
    check(coo_matrix%state_is_created())

    ! Copy coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%copy_to_coo(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix_copy (COO) <- coo_matrix (COO)
    call coo_matrix_copy%copy_from_coo(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) -> coo_matrix_fmt (FMT)
    call coo_matrix%copy_to_fmt(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix_copy%copy_from_fmt(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)

    ! Move coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%move_to_coo(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (COO)
    call coo_matrix%move_from_coo(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())
    ! Move coo_matrix (COO) -> coo_matrix_copy (FMT)
    call coo_matrix%move_to_fmt(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix%move_from_fmt(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())


    call coo_matrix%free()

    check(coo_matrix%state_is_start())

    call coo_matrix_copy%free()

    check(coo_matrix_copy%state_is_start())

!------------------------------------------------------------------
! NUMERIC
!------------------------------------------------------------------

    call coo_matrix%create(num_rows=5,num_cols=5, nz=10)
    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                           ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                           ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                           val=(/1.,2.,3.,4.,5.,5.,4.,3.,2.,1./), &
                           imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_numeric())

    ! Copy coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%copy_to_coo(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix_copy (COO) <- coo_matrix (COO)
    call coo_matrix_copy%copy_from_coo(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) -> coo_matrix_fmt (FMT)
    call coo_matrix%copy_to_fmt(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix_copy%copy_from_fmt(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)

    ! Move coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%move_to_coo(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (COO)
    call coo_matrix%move_from_coo(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())
    ! Move coo_matrix (COO) -> coo_matrix_copy (FMT)
    call coo_matrix%move_to_fmt(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix%move_from_fmt(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())


    call coo_matrix%free()

    check(coo_matrix%state_is_start())

    call coo_matrix_copy%free()

    check(coo_matrix_copy%state_is_start())

!------------------------------------------------------------------
! SYMBOLIC
!------------------------------------------------------------------

    call coo_matrix%create(num_rows=5,num_cols=5, nz=10)
    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                           ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                           ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                           imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_symbolic())

    ! Copy coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%copy_to_coo(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix_copy (COO) <- coo_matrix (COO)
    call coo_matrix_copy%copy_from_coo(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) -> coo_matrix_fmt (FMT)
    call coo_matrix%copy_to_fmt(coo_matrix_copy)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)
    ! Copy coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix_copy%copy_from_fmt(coo_matrix)
    call compare_coo_matrix(coo_matrix, coo_matrix_copy)

    ! Move coo_matrix (COO) -> coo_matrix_copy (COO)
    call coo_matrix%move_to_coo(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (COO)
    call coo_matrix%move_from_coo(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())
    ! Move coo_matrix (COO) -> coo_matrix_copy (FMT)
    call coo_matrix%move_to_fmt(coo_matrix_copy)
    check(coo_matrix%state_is_start())
    ! Move coo_matrix (COO) <- coo_matrix_copy (FMT)
    call coo_matrix%move_from_fmt(coo_matrix_copy)
    check(coo_matrix_copy%state_is_start())


    call coo_matrix%free()

    check(coo_matrix%state_is_start())

    call coo_matrix_copy%free()

    check(coo_matrix_copy%state_is_start())


    call memstatus()

contains

    subroutine compare_coo_matrix(a, b)
        class(coo_sparse_matrix_t), intent(in) :: a
        class(coo_sparse_matrix_t), intent(in) :: b
    
        check(a%get_num_rows()            == b%get_num_rows())
        check(a%get_num_cols()            == b%get_num_cols())
        check(a%get_nnz()                 == b%get_nnz())
        check(a%is_symmetric()           .eqv. b%is_symmetric())
        check(a%get_symmetric_storage()  .eqv. b%get_symmetric_storage())
        check(a%get_sign()                == b%get_sign())
        check(a%get_sort_status()         == b%get_sort_status())
!        check(a%get_state()               == b%get_state()) State must change after conversion
        check(all(a%ia                    == b%ia))
        check(all(a%ja                    == b%ja))
        check(allocated(a%val)           .eqv. allocated(b%val))

    end subroutine compare_coo_matrix

end program
