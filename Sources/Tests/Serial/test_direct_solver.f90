program test_direct_solver

USE types_names
USE memor_names
USE serial_names
USE sparse_matrix_names
USE direct_solver_names
USE FPL
USE IR_Precision



implicit none

# include "debug.i90"

    type(sparse_matrix_t)          :: sparse_matrix
    type(direct_solver_t)          :: direct_solver
    type(ParameterList_t)          :: parameter_list
    type(ParameterList_t), pointer :: direct_solver_parameters
    type(serial_scalar_array_t)    :: x
    type(serial_scalar_array_t)    :: y
    integer                        :: FPLError
    integer                        :: iparm(64) = 0

    call meminit()
    ! ParameterList: initialize
    call FPL_Init()
    call parameter_list%Init()
    direct_solver_parameters => parameter_list%NewSubList(Key=pardiso_mkl_name)

    ! ParameterList: set parameters
    FPLError = 0
    FPLError = FPLError + direct_solver_parameters%set(key = direct_solver_type,        value = pardiso_mkl_name)
    FPLError = FPLError + direct_solver_parameters%set(key = pardiso_mkl_matrix_type,   value = pardiso_mkl_sin)
    FPLError = FPLError + direct_solver_parameters%set(key = pardiso_mkl_message_level, value = 0)
    FPLError = FPLError + direct_solver_parameters%set(key = pardiso_mkl_iparm,         value = iparm)
    assert(FPLError == 0)

    call parameter_list%print()
    call direct_solver_parameters%print()

    ! Sparse matrix: create
    call sparse_matrix%create(4, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )
    call sparse_matrix%insert(nz=8,                    &
                              ia=(/1,1,2,3,3,3,4,4/),  &
                              ja=(/1,4,2,1,3,4,3,4/),  &
                              val=(/1.,2.,3.,4.,5.,6.,7.,8./))
    call sparse_matrix%convert('CSR')
    call sparse_matrix%print(6)

    ! Serial scalar array: create and initialize
    call x%create_and_allocate(sparse_matrix%get_num_rows())
    call x%init(1.0_rp)
    call y%create_and_allocate(sparse_matrix%get_num_cols())
    call x%print(6)

    ! Direct solver: create and set properties
    call direct_solver%set_type_from_parameter_list(direct_solver_parameters)
    call direct_solver%set_matrix(sparse_matrix)
    call direct_solver%set_defaults()
    call direct_solver%set_from_parameter_list(direct_solver_parameters)

    ! Direct solver: analisys, factorization and solve
    call direct_solver%symbolic_setup()
    call direct_solver%numerical_setup()
    call direct_solver%solve(x,y)
    call direct_solver%log_info()
    call y%print(6)

    ! Free
    call parameter_list%free()
    call direct_solver%free()
    call sparse_matrix%free()
    call x%free()
    call y%free()

    call memstatus()

end program test_direct_solver
