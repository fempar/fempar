# Sparse Matrix: State transition diagram

```fortran
  !-----------------------------------------------------------------
  ! State transition diagram for type(base_sparse_matrix_t)
  !-----------------------------------------------------------------
  ! Input State         | Action                | Output State 
  !-----------------------------------------------------------------
  ! Start               | Create                | Created
  ! Start               | Free_clean            | Start
  ! Start               | Free_symbolic         | Start
  ! Start               | Free_numeric          | Start
  !-----------------------------------------------------------------
  ! Created             | Insert (2-values)     | Build_symbolic
  ! Created             | Insert (3-values)     | Build_numeric
  ! Created             | Free_clean            | Start
  ! Created             | Free_symbolic         | Created
  ! Created             | Free_numeric          | Assembled_symbolic
  !-----------------------------------------------------------------
  ! Build_symbolic      | Insert (2-values)     | Build_symbolic
  ! Build_symbolic      | Insert (3-values)     | * Error
  ! Build_symbolic      | convert               | Assembled_symbolic
  ! Build_symbolic      | Free_clean            | Start
  ! Build_symbolic      | Free_symbolic         | Created
  ! Build_symbolic      | Free_numeric          | Created
  !-----------------------------------------------------------------
  ! Build_numeric       | Insert (2-values)     | * Error
  ! Build_numeric       | Insert (3-values)     | Build_numeric
  ! Build_numeric       | convert               | Assembled
  ! Build_numeric       | Free_clean            | Start
  ! Build_numeric       | Free_symbolic         | Created
  ! Build_numeric       | Free_numeric          | Created
  !-----------------------------------------------------------------
  ! Assembled           | Free_clean            | Start
  ! Assembled           | Free_symbolic         | Created
  ! Assembled           | Free_numeric          | Assembled_symbolic
  ! Assembled           | Insert (2-values)     | * Error
  ! Assembled           | Insert (3-values)     | Update
  !-----------------------------------------------------------------
  ! Assembled_symbolic  | Free_clean            | Start
  ! Assembled_symbolic  | Free_symbolic         | Created
  ! Assembled_symbolic  | Free_numeric          | Assembed_symbolic
  ! Assembled_symbolic  | Insert (2-values)     | * Error
  ! Assembled_symbolic  | Insert (3-values)     | Update
  !-----------------------------------------------------------------
  ! Assembled           | Free_clean            | Start
  ! Assembled           | Free_symbolic         | Created
  ! Assembled           | Free_numeric          | Assembled_symbolic
  ! Assembled           | Insert (2-values)     | * Error
  ! Assembled           | Insert (3-values)     | Update
  !-----------------------------------------------------------------```
```

# Sparse Matrix: Basic usage

```fortran
...
    type(sparse_matrix_t)       :: sparse_matrix
...
    ! Create: START=>CREATE
    call sparse_matrix%create(num_rows_and_cols=5,             &
                              symmetric_storage=.false.,       &
                              is_symmetric=.false.,            &
                              sign=SPARSE_MATRIX_SIGN_UNKNOWN, &
                              nz=10)

    ! Append: CREATE=>BUILD_SYMBOLIC
    call sparse_matrix%insert(ia=1, ja=1, imin=1, imax=5, jmin=1, jmax=5 )
    call sparse_matrix%insert(nz=10,                                 &
                              ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                              ja=(/1,2,3,4,5,5,4,3,2,1/))
    call sparse_matrix%insert(ia=1, ja=1)

    ! Convert: BUILD_SYMBOLIC=>ASSEMBLED_SYMBOLIC
    call sparse_matrix%convert('CSR')

    ! Update: ASSEMBLED_SYMBOLIC=>UPDATE
    call sparse_matrix%insert(ia=1, ja=1, val=1., imin=1, imax=5, jmin=1, jmax=5 )
    call sparse_matrix%insert(nz=10,                                 &
                              ia=(/1,2,3,4,5,1,2,3,4,5/),            &
                              ja=(/1,2,3,4,5,5,4,3,2,1/),            &
                              val=(/1.,2.,3.,4.,5.,5.,4.,3.,2.,1./), &
                              imin=1, imax=5, jmin=1, jmax=5 )
    call sparse_matrix%insert(ia=1, ja=1, val=1., imin=1, imax=5, jmin=1, jmax=5 )

    ! Update: UPDATE=>ASSEMBLED
    call sparse_matrix%convert('CSR')

    ! Free: ASSEMBLED=>ASSEMBLED_SYMBOLIC
    call sparse_matrix%free_in_stages(free_numerical_setup)
    ! Free: ASSEMBLED_SYMBOLIC=>CREATED
    call sparse_matrix%free_in_stages(free_symbolic_setup)
    ! Free: CREATED=>START
    call sparse_matrix%free_in_stages(free_clean)
...
```


