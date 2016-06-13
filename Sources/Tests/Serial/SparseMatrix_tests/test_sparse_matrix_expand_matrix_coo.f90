program test_expand_matrix

USE types_names
USE memor_names
USE serial_names
USE sparse_matrix_names
USE base_sparse_matrix_names, only: coo_sparse_matrix_t


implicit none

# include "debug.i90"

    type(sparse_matrix_t)       :: sparse_matrix
    type(coo_sparse_matrix_t)   :: C_T
    type(coo_sparse_matrix_t)   :: I
    type(sparse_matrix_t)       :: expanded_sparse_matrix
    integer(ip)                 :: num_rows_and_cols_to_expand = 5
    integer(ip)                 :: nz_to_expand = 5
    integer(ip), allocatable    :: C_T_ia(:)
    integer(ip), allocatable    :: C_T_ja(:)
    real(rp),    allocatable    :: C_T_val(:)
    integer(ip), allocatable    :: I_ia(:)
    integer(ip), allocatable    :: I_ja(:)
    real(rp),    allocatable    :: i_val(:)


    call meminit()

    call C_t%create(5, 5)
    call C_T%insert(nz=5,                      &
                    ia=(/1,1,3,5,5/),          &
                    ja=(/3,4,1,2,5/),          &
                    val=(/3.0,4.0,1.0,2.0,5.0/))
    call C_T%sort_and_compress()

print*, '!------------------------------------------------------------------'
print*, '! C_T SPARSE MATRIX'
print*, '!------------------------------------------------------------------'
    call C_T%print_matrix_market(6)

    call I%create(5, 5)
    call I%insert(nz=5,                      &
                  ia=(/1,1,2,3,4/),          &
                  ja=(/2,3,1,1,4/),          &
                  val=(/2.0,3.0,1.0,1.0,4.0/))
    call I%sort_and_compress()

print*, '!------------------------------------------------------------------'
print*, '! I SPARSE MATRIX'
print*, '!------------------------------------------------------------------'
    call I%print_matrix_market(6)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FULL MATRIX (NNZ>0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, '!------------------------------------------------------------------'
print*, '! ORIGINAL EMPTY (NNZ==0) SPARSE MATRIX'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%create(5, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%convert('CSR')

    call sparse_matrix%print_matrix_market(6)

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - NUMERIC - OPTIONAL I'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - NUMERIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - NUMERIC - NON SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I,                      &
                                             symmetric_storage = .false.)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - SYMBOLIC - OPTIONAL I'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - SYMBOLIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - SYMMETRIC STORAGE - SYMBOLIC - NON SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I,                      &
                                              symmetric_storage = .false.)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()


print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - NON SYMMETRIC STORAGE - NUMERIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%free()
    call sparse_matrix%create(5, .false., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%convert('CSR')

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - NON SYMMETRIC STORAGE - NUMERIC - SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I,                      &
                                             symmetric_storage = .true.)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()


print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - NON SYMMETRIC STORAGE - SYMBOLIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! EMPTY (NNZ==0) - NON SYMMETRIC STORAGE - SYMBOLIC - SYMMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I,                      &
                                              symmetric_storage = .true.)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()
    call sparse_matrix%free()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FULL MATRIX (NNZ>0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print*, '!------------------------------------------------------------------'
print*, '! ORIGINAL FULL (NNZ>0) SPARSE MATRIX'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%create(5, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%insert(nz=3,          &
                              ia=1,          &
                              ja=(/1,2,3/),  &
                              val=(/1.,2.,3./))
    call sparse_matrix%insert(nz=3,          &
                              ia=2,          &
                              ja=(/2,3,4/),  &
                              val=(/2.,3.,4./))
    call sparse_matrix%insert(nz=3,          &
                              ia=3,          &
                              ja=(/3,4,5/),  &
                              val=(/3.,4.,5./))
    call sparse_matrix%insert(nz=3,          &
                              ia=4,          &
                              ja=(/1,2,4/),  &
                              val=(/1.,2.,4./))
    call sparse_matrix%insert(nz=3,          &
                              ia=5,          &
                              ja=(/3,4,5/),  &
                              val=(/3.,4.,5./))

    call sparse_matrix%convert('CSR')

    call sparse_matrix%print_matrix_market(6)

print*, '!------------------------------------------------------------------'
print*, '! SYMMETRIC STORAGE - NUMERIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! SYMMETRIC STORAGE - NUMERIC - NON SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I,                      &
                                             symmetric_storage = .false.)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! SYMMETRIC STORAGE - SYMBOLIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! SYMMETRIC STORAGE - SYMBOLIC - NON SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I,                      &
                                              symmetric_storage = .false.)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()


print*, '!------------------------------------------------------------------'
print*, '! NON SYMMETRIC STORAGE - NUMERIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%free()
    call sparse_matrix%create(5, .false., .true., SPARSE_MATRIX_SIGN_UNKNOWN )

    call sparse_matrix%insert(nz=3,          &
                              ia=1,          &
                              ja=(/1,2,3/),  &
                              val=(/1.,2.,3./))
    call sparse_matrix%insert(nz=3,          &
                              ia=2,          &
                              ja=(/1,2,4/),  &
                              val=(/1.,2.,4./))
    call sparse_matrix%insert(nz=3,          &
                              ia=3,          &
                              ja=(/1,3,5/),  &
                              val=(/1.,3.,5./))
    call sparse_matrix%insert(nz=2,          &
                              ia=4,          &
                              ja=(/2,4/),    &
                              val=(/2.,4./))
    call sparse_matrix%insert(nz=2,          &
                              ia=5,          &
                              ja=(/3,5/),    &
                              val=(/3.,5./))

    call sparse_matrix%convert('CSR')


    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I) 

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! NON SYMMETRIC STORAGE - NUMERIC - SYMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_numeric(C_T = C_T,                    &
                                             to  = expanded_sparse_matrix, &
                                             I   = I,                      &
                                             symmetric_storage = .true.)

    call expanded_sparse_matrix%print_matrix_market(6)
    call expanded_sparse_matrix%free()


print*, '!------------------------------------------------------------------'
print*, '! NON SYMMETRIC STORAGE - SYMBOLIC'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I) 

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()

print*, '!------------------------------------------------------------------'
print*, '! NON SYMMETRIC STORAGE - SYMBOLIC - SYMMETRIC STORAGE FOR EXPANDED'
print*, '!------------------------------------------------------------------'

    call sparse_matrix%expand_matrix_symbolic(C_T = C_T,                    &
                                              to  = expanded_sparse_matrix, &
                                              I   = I,                      &
                                              symmetric_storage = .true.)

    call expanded_sparse_matrix%print(6)
    call expanded_sparse_matrix%free()
    call sparse_matrix%free()
    call C_T%free()
    call I%free()

    call memstatus()

end program
