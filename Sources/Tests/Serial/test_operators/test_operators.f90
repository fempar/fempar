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
    real(rp)                         :: vector_values(3)
    integer                          :: i, iters = 999

    call FEMPAR_INIT()

    call Mat%create(num_rows_and_cols=3, symmetric_storage=.false., is_symmetric=.false., sign=SPARSE_MATRIX_SIGN_UNKNOWN, nz=6)
    call Mat%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
    call Mat%convert('CSR')

    call Vec1%create_and_allocate(3)
    call Vec2%create_and_allocate(3)
    call Vec1%init(1.0_rp)

    vector_values = [4._rp, 3._rp,4._rp] 

    do i=1, iters
        Op = Mat
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, vector_values)

        Op = Mat + Mat
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, 2*vector_values)

        Op = Mat - Mat
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, 0*vector_values)

        Op = Mat * Mat 
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, [16._rp, 10._rp, 16._rp])

        Op = .minus. Mat
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, -vector_values)

        Op = 3.0*Mat
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, 3*vector_values)

        Op = Mat*3.0 
        call Op%apply(Vec1,Vec2)
        call check_vector_value(Vec2, 3*vector_values)
    enddo

    call Mat%free()
    call Vec1%free()
    call Vec2%free()
    call Op%free()

    call FEMPAR_FINALIZE()

contains

    subroutine check_vector_value(vector, values)
        class(serial_scalar_array_t), intent(in) :: vector
        real(rp),                     intent(in) :: values(:)
        real(rp), pointer                        :: vector_entries(:)
        nullify(vector_entries)
        vector_entries => vector%get_entries()
        check(associated(vector_entries))
        check(size(vector_entries) == size(values))
        check(all(vector_entries == values))
    end subroutine
  
end program test_operators
