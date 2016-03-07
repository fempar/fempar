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

    type(sparse_matrix_t) :: sparse_matrix
    type(direct_solver_t) :: direct_solver
    type(ParameterList_t) :: parameter_list
    integer               :: FPLError

    call meminit()

    call FPL_Init()

    call parameter_list%Init()
    FPLError = parameter_list%set(key=direct_solver_type, value='PARDISO_MKL')
    assert(FPLError == 0)

    call sparse_matrix%create(4, .true., .true., SPARSE_MATRIX_SIGN_UNKNOWN )
    call sparse_matrix%insert(nz=8,                    &
                              ia=(/1,1,2,3,3,3,4,4/),  &
                              ja=(/1,4,2,1,3,4,3,4/),  &
                              val=(/1.,2.,3.,4.,5.,6.,7.,8./))
    call sparse_matrix%convert('CSR')

    call direct_solver%set_type_from_parameter_list(parameter_list)
    call direct_solver%set_matrix(sparse_matrix)
    call direct_solver%set_defaults()
    call direct_solver%set_from_parameter_list(parameter_list)

    call parameter_list%print()

    call parameter_list%free()
    call direct_solver%free()
    call sparse_matrix%free()

    call memstatus()

end program test_direct_solver
