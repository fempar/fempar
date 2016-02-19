program test_coo_append

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
    check(coo_matrix%state_is_start())

    call coo_matrix%create(num_rows=5,num_cols=5)
    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                               ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                               ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                               val=(/1.,2.,3.,4.,5.,5.,4.,3.,2.,1./), &
                               imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_numeric())

    call coo_matrix%print( 6 )

    call coo_matrix%free()

    check(coo_matrix%state_is_start())

!------------------------------------------------------------------
! SYMBOLIC
!------------------------------------------------------------------

    call coo_matrix%create(num_rows=5,num_cols=5)
    check(coo_matrix%state_is_created())

    call coo_matrix%insert(nz=10,                                 &
                               ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                               ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                               imin=1, imax=5, jmin=1, jmax=5 )
    check(coo_matrix%state_is_build_symbolic())

    call coo_matrix%print( 6 )

    call coo_matrix%free()

    check(coo_matrix%state_is_start())

    call memstatus()

end program
