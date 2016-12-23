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
    type(sparse_matrix_t), pointer   :: Block
    type(block_sparse_matrix_t)      :: BlockMat
    type(serial_scalar_array_t)      :: Vec1
    type(serial_scalar_array_t)      :: Vec2
    type(serial_block_array_t)       :: BlockVec1
    type(serial_block_array_t)       :: BlockVec2
    type(lvalue_operator_t)          :: Op
    type(serial_scalar_array_t)      :: op_result
    real(rp)                         :: vector_values(3)
    integer                          :: tam = 3
    integer                          :: iters = 999
    integer                          :: i, j

    call FEMPAR_INIT()

    call Mat%create(num_rows_and_cols=tam, symmetric_storage=.false., is_symmetric=.false., sign=SPARSE_MATRIX_SIGN_UNKNOWN, nz=6)
    call Mat%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
    call Mat%convert('CSR')

    call Vec1%create_and_allocate(tam)
    call Vec2%create_and_allocate(tam)
    call Vec1%init(1.0_rp)

    vector_values = [4._rp, 3._rp,4._rp] 

    call test_operators_apply(Mat, Vec1, Vec2, 999)
    call test_operators_apply_Add(Mat, Vec1, Vec2, 999)

    call Mat%free()
    call Vec1%free()
    call Vec2%free()

    call BlockMat%create(2, [tam,tam], [tam, tam], [.false.,.false.], [.false.,.false.], [SPARSE_MATRIX_SIGN_UNKNOWN,SPARSE_MATRIX_SIGN_UNKNOWN])
    do i=1, BlockMat%get_nblocks()
        do j=1, BlockMat%get_nblocks()
            Block => BlockMat%get_block(i,j)
            call Block%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
        enddo
    enddo
    call BlockMat%compress_storage('CSR')

    call BlockVec1%create_and_allocate(2, [tam, tam])
    call BlockVec2%create_and_allocate(2, [tam, tam])
    call BlockVec1%init(1.0_rp)

    call test_operators_apply(BlockMat, BlockVec1, BlockVec2, 999)
    call test_operators_apply_add(BlockMat, BlockVec1, BlockVec2, 999)

    call BlockMat%free()
    call BlockVec1%free()
    call BlockVec2%free()
    call Op%free()
    call Op_result%free()

    call FEMPAR_FINALIZE()

contains

    subroutine check_vector_value(vector, result_vector)
        class(vector_t), intent(in) :: vector
        class(vector_t), intent(in) :: result_vector
        class(vector_t), allocatable :: error
        allocate(error,mold=vector)
        error = vector-result_vector
        check(error%nrm2()<=EPSILON(1.0_rp)*1.e3*result_vector%nrm2())
        call error%free()
        deallocate(error) 
    end subroutine

    subroutine test_operators_apply(Mat, x, y, iters)
        class(matrix_t), intent(in)      :: Mat
        class(vector_t), intent(in)      :: x
        class(vector_t), intent(inout)   :: y
        integer, intent(in)              :: iters
        type(lvalue_operator_t)          :: Op
        class(vector_t), allocatable     :: op_result
        real(rp)                         :: vector_values(3)
        integer                          :: tam = 3
        integer                          :: i, j

        vector_values = [4._rp, 3._rp,4._rp] 

        allocate(op_result, mold=y)
        select type(op_result)
            type is (serial_scalar_array_t)
                call op_result%create_and_allocate(tam)
            type is (serial_block_array_t)
                call op_result%create_and_allocate(y%get_number_blocks(), [tam, tam])
            class default
                check(.false.)
        end select

        do i=1, iters
            call Op%assign(Mat)
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat + Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 2.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat - Mat
            y = Op * x
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat * Mat 
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp])
            enddo
            call check_vector_value(y, op_result)

            Op = .minus. Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], -y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = .identity. Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values+1.0)
            enddo
            call check_vector_value(y, op_result)

            Op = 3.0*Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 3.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = 0.0*Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat*3.0 
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 3.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op =  .minus. Mat + Mat * .identity. Mat * 2.0 - Mat
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat * .identity. Mat * 2.0 - Mat
            Op = .minus. Mat + Op
            call Op%apply(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values)
            enddo
            call check_vector_value(y, op_result)
        enddo

        call Op%free()
        call Op_result%free()
        deallocate(Op_result)

    end subroutine


    subroutine test_operators_apply_add(Mat, x, y, iters)
        class(matrix_t), intent(in)      :: Mat
        class(vector_t), intent(in)      :: x
        class(vector_t), intent(inout)   :: y
        integer, intent(in)              :: iters
        type(lvalue_operator_t)          :: Op
        class(vector_t), allocatable     :: op_result
        real(rp)                         :: vector_values(3)
        integer                          :: tam = 3
        integer                          :: i, j

        vector_values = [4._rp, 3._rp,4._rp] 

        allocate(op_result, mold=y)
        select type(op_result)
            type is (serial_scalar_array_t)
                call op_result%create_and_allocate(tam)
            type is (serial_block_array_t)
                call op_result%create_and_allocate(y%get_number_blocks(), [tam, tam])
            class default
                check(.false.)
        end select

        do i=1, iters
            call y%init(0.0_rp)

            Op = 0.0*Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 0.0*vector_values)
            enddo
            call check_vector_value(y, op_result)

            call Op%assign(Mat)
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 1.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat + Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 3.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat - Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 3.0*y%get_number_blocks()*vector_values)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat * Mat 
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 3.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp])
            enddo
            call check_vector_value(y, op_result)

            Op = .minus. Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 2.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp])
            enddo
            call check_vector_value(y, op_result)

            Op = .identity. Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 2.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp]+1.0)
            enddo
            call check_vector_value(y, op_result)

            Op = 3.0*Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 5.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp]+1.0)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat*3.0 
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 8.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp]+1.0)
            enddo
            call check_vector_value(y, op_result)

            Op =  .minus. Mat + Mat * .identity. Mat * 2.0 - Mat
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 8.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp]+1.0)
            enddo
            call check_vector_value(y, op_result)

            Op = Mat * .identity. Mat * 2.0 - Mat
            Op = .minus. Mat + Op
            call Op%apply_add(x,y)
            do j=1, y%get_number_blocks()
                call op_result%insert_subvector(j, tam, [1,2,3], 8.0*y%get_number_blocks()*vector_values+y%get_number_blocks()*y%get_number_blocks()*[16._rp, 10._rp, 16._rp]+1.0)
            enddo
            call check_vector_value(y, op_result)
        enddo

        call Op%free()
        call Op_result%free()
        deallocate(Op_result)

    end subroutine

  
end program test_operators
