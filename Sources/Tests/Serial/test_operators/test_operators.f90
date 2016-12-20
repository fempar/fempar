! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_operators
    use fempar_names
implicit none
#include "debug.i90"
    type(sparse_matrix_t)            :: Mat
    type(serial_scalar_array_t)      :: Vec1
    type(serial_scalar_array_t)      :: Vec2
    type(lvalue_operator_t)          :: Op
    type(serial_scalar_array_t)      :: op_result
    real(rp)                         :: vector_values(3)
    integer                          :: tam = 3
    integer                          :: iters = 999
    integer                          :: i

    call FEMPAR_INIT()

    call Mat%create(num_rows_and_cols=tam, symmetric_storage=.false., is_symmetric=.false., sign=SPARSE_MATRIX_SIGN_UNKNOWN, nz=6)
    call Mat%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
    call Mat%convert('CSR')

    call Vec1%create_and_allocate(tam)
    call Vec2%create_and_allocate(tam)
    call Vec1%init(1.0_rp)

    vector_values = [4._rp, 3._rp,4._rp] 

    call op_result%create(3)
    call op_result%allocate()

    do i=1, iters
        Op = Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], vector_values)
        call check_vector_value(Vec2, op_result)

        Op = Mat + Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 2.0*vector_values)
        call check_vector_value(Vec2, op_result)

        Op = Mat - Mat
        Vec2 = Op * Vec1
        call op_result%insert_subvector(1, tam, [1,2,3], 0.0*vector_values)
        call check_vector_value(Vec2, op_result)

        Op = Mat * Mat 
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], [16._rp, 10._rp, 16._rp])

        Op = .minus. Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], -vector_values)
        call check_vector_value(Vec2, op_result)

        Op = .identity. Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 0.0*vector_values-1.0)
        call check_vector_value(Vec2, op_result)

        Op = 3.0*Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 3.0*vector_values)
        call check_vector_value(Vec2, op_result)

        Op = Mat*3.0 
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 3.0*vector_values)
        call check_vector_value(Vec2, op_result)

        Op =  .minus. Mat + Mat * .identity. Mat * 2.0 - Mat
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 0.0*vector_values)
        call check_vector_value(Vec2, op_result)

        Op = Mat * .identity. Mat * 2.0 - Mat
        Op = .minus. Mat + Op
        call Op%apply(Vec1,Vec2)
        call op_result%insert_subvector(1, tam, [1,2,3], 0.0*vector_values)
        call check_vector_value(Vec2, op_result)

    enddo

    call Mat%free()
    call Vec1%free()
    call Vec2%free()
    call Op%free()
    call Op_result%free()

    call FEMPAR_FINALIZE()

contains

    subroutine check_vector_value(vector, result_vector)
        type(serial_scalar_array_t), intent(in) :: vector
        type(serial_scalar_array_t), intent(in) :: result_vector
        check(abs(vector%nrm2()-result_vector%nrm2())<EPSILON(1.0_rp)*1.e3)
    end subroutine
  
end program test_operators
