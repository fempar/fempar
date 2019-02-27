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
  integer                          :: itest
  integer                          :: tam
  integer                          :: tam_init = 3
  integer                          :: tam_adapted = 4
  integer                          :: iters = 999
  integer                          :: iter_reallocate = 500
  integer                          :: i, j

  call FEMPAR_INIT()

  assert (iters > iter_reallocate)

  ! Working with sparse matrices   
  do itest = 1,2
     tam = tam_init
     call create_matrix(Mat,tam)
     call create_and_allocate_arrays(Vec1,Vec2,tam)
     if (itest==1) then
        call test_operators_apply(Mat, Vec1, Vec2, 999)
     else
        call test_operators_apply_Add(Mat, Vec1, Vec2, 999)
     end if
  enddo
  call Mat%free()
  call Vec1%free()
  call Vec2%free()

  ! Working with Block sparse matrices   
  do itest = 1,2  
     tam = tam_init 
     call create_matrix(BlockMat,tam)
     call create_and_allocate_arrays(BlockVec1,BlockVec2,tam)
     if (itest==1) then
        call test_operators_apply(BlockMat, BlockVec1, BlockVec2, 999)
     else      
        call test_operators_apply_add(BlockMat, BlockVec1, BlockVec2, 999)
     end if
  enddo
  call BlockMat%free()
  call BlockVec1%free()
  call BlockVec2%free()

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
  end subroutine check_vector_value

  subroutine test_operators_apply(Mat, x, y, iters)
    class(matrix_t), intent(inout)   :: Mat
    class(vector_t), intent(inout)   :: x
    class(vector_t), intent(inout)   :: y
    integer,  intent(in)             :: iters
    type(lvalue_operator_t)          :: Op
    class(vector_t), allocatable     :: Op_result
    integer, allocatable             :: indices(:)
    integer                          :: i, j
    real(rp), allocatable            :: vector_values(:)        
    real(rp), allocatable            :: Mat_x_Mat_values(:)        

    call memalloc(tam,indices,__FILE__,__LINE__)
    indices = [1,2,3]

    call memalloc(tam,vector_values, __FILE__, __LINE__ ) 
    vector_values  = [4._rp, 3._rp, 4._rp]         

    call memalloc(tam,Mat_x_Mat_values, __FILE__, __LINE__ )
    Mat_x_Mat_values = [16._rp, 10._rp, 16._rp]

    call create_and_allocate_result_operator(Op_result,y,tam)

    do i=1, iters
       call Op%assign(Mat)
       call Op%apply(x,y)           
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat + Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 2.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat - Mat
       y = Op * x
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat * Mat 
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = .minus. Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, -y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = .identity. Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       Op = 3.0*Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 3.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = 0.0*Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat*3.0 
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 3.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op =  .minus. Mat + Mat * .identity. Mat * 2.0 - Mat
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat * .identity. Mat * 2.0 - Mat
       Op = .minus. Mat + Op
       call Op%apply(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       if (i == iter_reallocate) then
          ! Test update_after_remesh TBP of lvalue_operator_t
          ! 1. Destroy sparse matrix and associated arrays and
          !    generate new ones of different size
          tam = tam_adapted
          call create_adapted_sparse_matrix(Mat,tam)
          call create_and_allocate_arrays(x,y,tam)
          call create_and_allocate_result_operator(Op_result,y,tam)
          call adapt_vector_values_array(indices, vector_values, Mat_x_Mat_values, tam)              

          ! 2. Update direct solver to reflect that the sparse matrix has changed 
          !    (typically after AMR mesh adaptation, although not necessarily)
          call Op%reallocate_after_remesh()
       end if

    enddo
    call memfree(indices, __FILE__, __LINE__)
    call memfree(vector_values, __FILE__, __LINE__)
    call memfree(Mat_x_Mat_values, __FILE__, __LINE__)
    call Op%free()
    call Op_result%free()
    deallocate(Op_result)
  end subroutine test_operators_apply

  subroutine test_operators_apply_add(Mat, x, y, iters)
    class(matrix_t), intent(inout)   :: Mat
    class(vector_t), intent(inout)   :: x
    class(vector_t), intent(inout)   :: y
    integer, intent(in)              :: iters
    type(lvalue_operator_t)          :: Op
    class(vector_t), allocatable     :: Op_result
    integer, allocatable             :: indices(:)
    integer                          :: i, j
    real(rp), allocatable            :: vector_values(:)        
    real(rp), allocatable            :: Mat_x_Mat_values(:)  

    call memalloc(tam,indices,__FILE__,__LINE__)
    indices = [1,2,3]

    call memalloc(tam,vector_values, __FILE__, __LINE__ ) 
    vector_values  = [4._rp, 3._rp, 4._rp]         

    call memalloc(tam,Mat_x_Mat_values, __FILE__, __LINE__ )
    Mat_x_Mat_values = [16._rp, 10._rp, 16._rp]

    call create_and_allocate_result_operator(Op_result,y,tam)

    do i=1, iters
       call y%init(0.0_rp)

       Op = 0.0*Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 0.0*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       call Op%assign(Mat)
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 1.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat + Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 3.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat - Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 3.0*y%get_num_blocks()*vector_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat * Mat 
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 3.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = .minus. Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 2.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values)
       enddo
       call check_vector_value(y, Op_result)

       Op = .identity. Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 2.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       Op = 3.0*Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 5.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat*3.0 
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 8.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       Op =  .minus. Mat + Mat * .identity. Mat * 2.0 - Mat
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 8.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       Op = Mat * .identity. Mat * 2.0 - Mat
       Op = .minus. Mat + Op
       call Op%apply_add(x,y)
       do j=1, y%get_num_blocks()
          call Op_result%insert_subvector(j, tam, indices, 8.0*y%get_num_blocks()*vector_values+y%get_num_blocks()*y%get_num_blocks()*Mat_x_Mat_values+1.0)
       enddo
       call check_vector_value(y, Op_result)

       if (i == iter_reallocate) then
          ! Test update_after_remesh TBP of lvalue_operator_t
          ! 1. Destroy sparse matrix and associated arrays and
          !    generate new ones of different size
          tam = tam_adapted
          call create_adapted_sparse_matrix(Mat,tam)
          call create_and_allocate_arrays(x,y,tam)
          call create_and_allocate_result_operator(Op_result,y,tam)
          call adapt_vector_values_array(indices, vector_values, Mat_x_Mat_values, tam)              

          ! 2. Update direct solver to reflect that the sparse matrix has changed 
          !    (typically after AMR mesh adaptation, although not necessarily)
          call Op%reallocate_after_remesh()
       end if

    enddo
    call memfree(indices, __FILE__, __LINE__)
    call memfree(vector_values, __FILE__, __LINE__)
    call memfree(Mat_x_Mat_values, __FILE__, __LINE__)
    call Op%free()
    call Op_result%free()
    deallocate(Op_result)
  end subroutine test_operators_apply_add

  subroutine create_and_allocate_arrays(array1,array2,tam)
    implicit none
    class(vector_t),  intent(inout) :: array1
    class(vector_t),  intent(inout) :: array2
    integer,          intent(in)    :: tam     
    select type(array1)
    type is (serial_scalar_array_t)     
       call Vec1%free()
       call Vec1%create_and_allocate(Mat%get_num_rows())
       call Vec1%init(1.0_rp)
    type is (serial_block_array_t)
       call BlockVec1%free()
       call BlockVec1%create_and_allocate(2, [tam, tam])
       call BlockVec1%init(1.0_rp)
       class default
       check(.false.)
    end select
    select type(array2)
    type is (serial_scalar_array_t)
       call Vec2%free()
       call Vec2%create_and_allocate(Mat%get_num_cols())
    type is (serial_block_array_t)
       call BlockVec2%free()
       call BlockVec2%create_and_allocate(2, [tam, tam])
       class default
       check(.false.)
    end select
  end subroutine create_and_allocate_arrays

  subroutine create_and_allocate_result_operator(Op_result,y,tam)
    implicit none
    class(vector_t), allocatable,  intent(inout) :: Op_result
    class(vector_t),               intent(in)    :: y      
    integer,                       intent(in)    :: tam
    if (allocated(Op_result)) then
       call Op_result%free()
       deallocate(Op_result)
    end if
    allocate(Op_result, mold=y)
    select type(Op_result)
    type is (serial_scalar_array_t)
       call Op_result%create_and_allocate(tam)
    type is (serial_block_array_t)
       call Op_result%create_and_allocate(y%get_num_blocks(), [tam, tam])
       class default
       check(.false.)
    end select
  end subroutine create_and_allocate_result_operator

  subroutine adapt_vector_values_array(indices,vector_values,Mat_x_Mat_values,tam)
    implicit none
    integer,  allocatable,  intent(inout) :: indices(:)
    real,     allocatable,  intent(inout) :: vector_values(:)
    real,     allocatable,  intent(inout) :: Mat_x_Mat_values(:)
    integer,                intent(in)    :: tam     
    call memrealloc(tam,indices,__FILE__,__LINE__)
    indices = [1,2,3,4]
    call memrealloc(tam,vector_values,__FILE__,__LINE__)
    vector_values = [8._rp, 4._rp, 6._rp, 5._rp]    
    call memrealloc(tam,Mat_x_Mat_values,__FILE__,__LINE__)
    Mat_x_Mat_values = [42._rp, 17._rp, 40._rp, 34._rp]        
  end subroutine  adapt_vector_values_array

  subroutine create_matrix(Mat,tam)
    implicit none
    class(matrix_t),   intent(inout)  :: Mat
    integer,           intent(in)     :: tam
    call Mat%free()  
    select type(Mat)
    type is (sparse_matrix_t)
       call Mat%create(num_rows_and_cols=tam, symmetric_storage=.false., is_symmetric=.false., sign=SPARSE_MATRIX_SIGN_UNKNOWN, nz=6)
       call Mat%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
       call Mat%convert('CSR')
    type is (block_sparse_matrix_t)
       call BlockMat%create(2, [tam,tam], [tam, tam], [.false.,.false.], [.false.,.false.], [SPARSE_MATRIX_SIGN_UNKNOWN,SPARSE_MATRIX_SIGN_UNKNOWN])
       do i=1, BlockMat%get_nblocks()
          do j=1, BlockMat%get_nblocks()
             Block => BlockMat%get_block(i,j)
               call Block%insert(nz=6, ia=[1,1,2,2,3,3], ja=[1,3,1,2,1,3], val=[1.0,3.0,1.0,2.0,1.0,3.0])
            enddo
         enddo
         call BlockMat%compress_storage('CSR')
         class default
         check(.false.)
      end select
    end subroutine create_matrix

    subroutine create_adapted_sparse_matrix(Mat,tam)
      implicit none
      class(matrix_t),       intent(inout)  :: Mat
      integer,               intent(in)     :: tam     
      call Mat%free()  
      select type(Mat)
      type is (sparse_matrix_t)
         call Mat%create(num_rows_and_cols=tam, &
              symmetric_storage=.false., &
              is_symmetric=.false., &
              sign=SPARSE_MATRIX_SIGN_UNKNOWN, &
              nz=9)
         call Mat%insert(nz=9, &
              ia=[1,1,1,2,2,3,3,4,4], &
              ja=[1,2,4,2,4,1,3,1,4], &
              val=[1.0,1.0,6.0,3.0,1.0,2.0,4.0,3.0,2.0])
         call Mat%convert('CSR')
      type is (block_sparse_matrix_t)
         call Mat%create(2, [tam,tam], [tam, tam], [.false.,.false.], [.false.,.false.], [SPARSE_MATRIX_SIGN_UNKNOWN,SPARSE_MATRIX_SIGN_UNKNOWN])
         do i=1, BlockMat%get_nblocks()
            do j=1, BlockMat%get_nblocks()
               Block => BlockMat%get_block(i,j)
                 call Block%insert(nz=9, &
                      ia=[1,1,1,2,2,3,3,4,4], &
                      ja=[1,2,4,2,4,1,3,1,4], &
                      val=[1.0,1.0,6.0,3.0,1.0,2.0,4.0,3.0,2.0])
              enddo
           enddo
           call BlockMat%compress_storage('CSR')
           class default
           check(.false.)
        end select
      end subroutine create_adapted_sparse_matrix    

    end program test_operators
