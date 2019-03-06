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

!* The `[[test_function_parser]]`,
program test_function_parser

  use fempar_names
  use test_coded_function_names
  implicit none
# include "debug.i90"

  type(fempar_parameter_handler_t)            :: parameter_handler
  type(ParameterList_t), pointer              :: parameter_list

  type(scalar_coded_function_add_1_2_ops_t)   :: scalar_coded_function_add_1_2_ops
  type(scalar_coded_function_mul_1_2_ops_t)   :: scalar_coded_function_mul_1_2_ops
  type(scalar_coded_function_div_1_2_ops_t)   :: scalar_coded_function_div_1_2_ops
  type(scalar_coded_function_pow_1_2_ops_t)   :: scalar_coded_function_pow_1_2_ops

  type(scalar_coded_function_add_3_5_ops_t)   :: scalar_coded_function_add_3_5_ops
  type(scalar_coded_function_mul_3_5_ops_t)   :: scalar_coded_function_mul_3_5_ops
  type(scalar_coded_function_div_3_5_ops_t)   :: scalar_coded_function_div_3_5_ops
  type(scalar_coded_function_pow_3_5_ops_t)   :: scalar_coded_function_pow_3_5_ops

  type(scalar_function_parser_t)              :: scalar_function_parser

  type(vector_coded_function_t)               :: vector_coded_function
  type(tensor_coded_function_t)               :: tensor_coded_function

  type(scalar_function_parser_t)              :: scalar_function_parser_1
  type(scalar_function_parser_t)              :: scalar_function_parser_2
  type(scalar_function_parser_t)              :: scalar_function_parser_3
  type(scalar_function_parser_t)              :: scalar_function_parser_4
  type(scalar_function_parser_t)              :: scalar_function_parser_5
  type(scalar_function_parser_t)              :: scalar_function_parser_6
  type(scalar_function_parser_t)              :: scalar_function_parser_7
  type(scalar_function_parser_t)              :: scalar_function_parser_8
  type(scalar_function_parser_t)              :: scalar_function_parser_9
  type(vector_function_parser_t)              :: vector_function_parser
  type(tensor_function_parser_t)              :: tensor_function_parser

  type(scalar_function_and_gradient_parser_t) :: scalar_function_and_gradient_parser
  type(vector_function_and_gradient_parser_t) :: vector_function_and_gradient_parser

  type(point_t), allocatable                  :: points(:)
  real(rp), allocatable                       :: parser_result(:)
  real(rp), allocatable                       :: result(:)
  type(vector_field_t), allocatable           :: vector_parser_result(:)
  type(vector_field_t), allocatable           :: vector_result(:)
  type(tensor_field_t), allocatable           :: tensor_parser_result(:)
  type(tensor_field_t), allocatable           :: tensor_result(:)
  character(len=:), allocatable               :: ftype, tname, op, snum, sops
  integer                                     :: num, ops, i, j, k
  real(rp)                                    :: error


  call fempar_init()

  call parameter_handler%process_parameters(define_test_function_parameters)
  parameter_list => parameter_handler%get_values()
  error = parameter_list%get(key='num_points',    value=num)
  error = parameter_list%get(key='num_ops',       value=ops)
  error = parameter_list%getAsString(key='num_points',    string=snum)
  error = parameter_list%getAsString(key='num_ops',       string=sops)
  error = parameter_list%getAsString(key='function_type', string=ftype)
  error = parameter_list%getAsString(key='operator',      string=op)
  error = parameter_list%getAsString(key='test_name',     string=tname)
  assert(error == 0)
  assert(ops >= 2 .and. ops <= 5)
  assert(ftype=='scalar' .or. ftype=='vector' .or. ftype=='tensor' .or. ftype=='scalar_function_and_gradient' .or. ftype=='vector_function_and_gradient')
  if(ftype=='scalar') then
      assert(op == '+' .or. op == '/' .or. op == '*' .or. op == '^')
      write(*,*) '[TEST] ', tname, ', type: ', ftype, ', operator: ', op, ', points: ', snum, ', num operators: ', sops
  else
      write(*,*) '[TEST] ', tname, ', type: ', ftype, ', points: ', snum
  endif

    error = 1.e-14
    call create_random_points(num, points)

    if(ftype == 'scalar') then
      allocate(parser_result(num))
      allocate(result(num))
      if(ops == 2) then
        select case(op)
          case('+')
            call scalar_function_parser%create("x+y", num_dims=2)
            call scalar_coded_function_add_1_2_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_add_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('*')
            call scalar_function_parser%create("x*y", num_dims=2)
            call scalar_coded_function_mul_1_2_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_mul_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('/')
            call scalar_function_parser%create("x/y", num_dims=2)
            call scalar_coded_function_div_1_2_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_div_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('^')
            call scalar_function_parser%create("x^y", num_dims=2)
            call scalar_coded_function_pow_1_2_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_pow_1_2_ops%get_value_space(points(i), result(i))
            enddo
        end select
      elseif(ops == 3) then
        select case(op)
          case('+')
            call scalar_function_parser%create("x+y+z", num_dims=3)
            call scalar_coded_function_add_1_2_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_add_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('*')
            call scalar_function_parser%create("x*y*z", num_dims=3)
            call scalar_coded_function_mul_1_2_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_mul_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('/')
            call scalar_function_parser%create("x/y/z", num_dims=3)
            call scalar_coded_function_div_1_2_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_div_1_2_ops%get_value_space(points(i), result(i))
            enddo
          case('^')
            call scalar_function_parser%create("x^y^z", num_dims=3)
            call scalar_coded_function_pow_1_2_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_pow_1_2_ops%get_value_space(points(i), result(i))
            enddo
        end select
      elseif(ops == 4) then
        select case(op)
          case('+')
            call scalar_function_parser%create("x+y+x+y", num_dims=2)
            call scalar_coded_function_add_3_5_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_add_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('*')
            call scalar_function_parser%create("x*y*x*y", num_dims=2)
            call scalar_coded_function_mul_3_5_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_mul_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('/')
            call scalar_function_parser%create("x/y/x/y", num_dims=2)
            call scalar_coded_function_div_3_5_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_div_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('^')
            call scalar_function_parser%create("x^y^x^y", num_dims=2)
            call scalar_coded_function_pow_3_5_ops%set_num_dims(2)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_pow_3_5_ops%get_value_space(points(i), result(i))
            enddo
        end select
      elseif(ops == 5) then
        select case(op)
          case('+')
            call scalar_function_parser%create("x+y+z+x+y+z", num_dims=3)
            call scalar_coded_function_add_3_5_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_add_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('*')
            call scalar_function_parser%create("x*y*z*x*y*z", num_dims=3)
            call scalar_coded_function_mul_3_5_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_mul_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('/')
            call scalar_function_parser%create("x/y/z/x/y/z", num_dims=3)
            call scalar_coded_function_div_3_5_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_div_3_5_ops%get_value_space(points(i), result(i))
            enddo
          case('^')
            call scalar_function_parser%create("x^y^z^x^y^z", num_dims=3)
            call scalar_coded_function_pow_3_5_ops%set_num_dims(3)
            do i=1, num
              call scalar_function_parser%get_value_space(points(i), parser_result(i))
            enddo
            do i=1, num
              call scalar_coded_function_pow_3_5_ops%get_value_space(points(i), result(i))
            enddo
        end select
      endif
     assert(all(abs(result - parser_result) < error))
    elseif(ftype == 'vector') then
      allocate(vector_parser_result(num))
      allocate(vector_result(num))
      if(ops == 2) then
        call scalar_function_parser_1%create("x+y", num_dims=2)
        call scalar_function_parser_2%create("-x-y", num_dims=2)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2)
        call vector_coded_function%set_num_dims(2)
        do i=1, num
          call vector_function_parser%get_value_space(points(i), vector_parser_result(i))
        enddo
        do i=1, num
          call vector_coded_function%get_value_space(points(i), vector_result(i))
        enddo
      elseif(ops == 3) then
        call scalar_function_parser_1%create("x+y+z", num_dims=3)
        call scalar_function_parser_2%create("-x-y-z", num_dims=3)
        call scalar_function_parser_3%create("x*y*z", num_dims=3)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2, scalar_function_parser_3)
        call vector_coded_function%set_num_dims(3)
        do i=1, num
          call vector_function_parser%get_value_space(points(i), vector_parser_result(i))
        enddo
        do i=1, num
          call vector_coded_function%get_value_space(points(i), vector_result(i))
        enddo
      endif
      assert(all((/ (abs(vector_result(i)%get_value() - vector_parser_result(i)%get_value()) < error , i=1, num) /)))
    elseif(ftype == 'tensor') then
      allocate(tensor_parser_result(num))
      allocate(tensor_result(num))
      if(ops == 2) then
         ! [(x+y, -x-y), (x*y, x/y)]
        call scalar_function_parser_1%create("x+y", num_dims=2)
        call scalar_function_parser_2%create("-x-y", num_dims=2)
        call scalar_function_parser_3%create("x*y", num_dims=2)
        call scalar_function_parser_4%create("x/y", num_dims=2)
        call tensor_function_parser%create(scalar_function_parser_1, scalar_function_parser_2, &
                                           scalar_function_parser_3, scalar_function_parser_4)
        call tensor_coded_function%set_num_dims(2)
        do i=1, num
          call tensor_function_parser%get_value_space(points(i), tensor_parser_result(i))
        enddo
        do i=1, num
          call tensor_coded_function%get_value_space(points(i), tensor_result(i))
        enddo
      elseif(ops == 3) then
        ! [(x+y+z, -x-y-z, x*y*z), (x+y-z, -x-y+z, x/y/z), (x-y-z, -x+y-z, x^y^z)]
        call scalar_function_parser_1%create("x+y+z", num_dims=3)
        call scalar_function_parser_2%create("-x-y-z", num_dims=3)
        call scalar_function_parser_3%create("x*y*z", num_dims=3)
        call scalar_function_parser_4%create("x+y-z", num_dims=3)
        call scalar_function_parser_5%create("-x-y+z", num_dims=3)
        call scalar_function_parser_6%create("x/y/z", num_dims=3)
        call scalar_function_parser_7%create("x-y-z", num_dims=3)
        call scalar_function_parser_8%create("-x+y-z", num_dims=3)
        call scalar_function_parser_9%create("x^y^z", num_dims=3)
        call tensor_function_parser%create(scalar_function_parser_1, scalar_function_parser_2, scalar_function_parser_3, &
                                           scalar_function_parser_4, scalar_function_parser_5, scalar_function_parser_6, &
                                           scalar_function_parser_7, scalar_function_parser_8, scalar_function_parser_9)
        call tensor_coded_function%set_num_dims(3)
        do i=1, num
          call tensor_function_parser%get_value_space(points(i), tensor_parser_result(i))
        enddo
        do i=1, num
          call tensor_coded_function%get_value_space(points(i), tensor_result(i))
        enddo
      endif
      assert(all((/(((abs(tensor_result(i)%get(j,k)-tensor_parser_result(i)%get(j,k))<error,k=1,ops),j=1,ops),i=1,num)/)))
    elseif(ftype == 'scalar_function_and_gradient') then
      allocate(parser_result(num))
      allocate(result(num))
      allocate(vector_parser_result(num))
      allocate(vector_result(num))
      if(ops == 2) then
        call scalar_function_parser%create("x^2", num_dims=2)
        call scalar_function_parser_1%create("2*x", num_dims=2)
        call scalar_function_parser_2%create("0", num_dims=2)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2)
        call scalar_function_and_gradient_parser%create(scalar_function_parser, vector_function_parser, scalar_function_parser_2)
        do i=1, num
          call scalar_function_and_gradient_parser%get_value_temporal_derivative(points(i), 0._rp, 1, parser_result(i))
          call scalar_function_and_gradient_parser%get_value_space(points(i), parser_result(i))
          call scalar_function_and_gradient_parser%get_gradient_space(points(i), vector_parser_result(i))
          result(i) = points(i)%get(1)**2
          call vector_result(i)%set(1, 2*points(i)%get(1))
        enddo
      elseif(ops == 3) then
        call scalar_function_parser%create("x^2", num_dims=3)
        call scalar_function_parser_1%create("2*x", num_dims=3)
        call scalar_function_parser_2%create("0", num_dims=3)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2, scalar_function_parser_2)
        call scalar_function_and_gradient_parser%create(scalar_function_parser, vector_function_parser)
        do i=1, num
          call scalar_function_and_gradient_parser%get_value_space(points(i), parser_result(i))
          call scalar_function_and_gradient_parser%get_gradient_space(points(i), vector_parser_result(i))
          result(i) = points(i)%get(1)**2
          call vector_result(i)%set(1, 2*points(i)%get(1))
        enddo
      endif
      assert(all(abs(result - parser_result) < error))
      assert(all((/ (abs(vector_result(i)%get_value() - vector_parser_result(i)%get_value()) < error , i=1, num) /)))
    elseif(ftype == 'vector_function_and_gradient') then
      allocate(tensor_parser_result(num))
      allocate(tensor_result(num))
      allocate(vector_parser_result(num))
      allocate(vector_result(num))
      if(ops == 2) then
        call scalar_function_parser_1%create("x^2", num_dims=2)
        call scalar_function_parser_2%create("0",   num_dims=2)
        call scalar_function_parser_3%create("2*x", num_dims=2)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2)
        call tensor_function_parser%create(scalar_function_parser_3, scalar_function_parser_2, scalar_function_parser_3, scalar_function_parser_2)
        call vector_function_and_gradient_parser%create(vector_function_parser, tensor_function_parser)
        do i=1, num
          call vector_function_and_gradient_parser%get_value_space(points(i), vector_parser_result(i))
          call vector_function_and_gradient_parser%get_gradient_space(points(i), tensor_parser_result(i))
          call vector_result(i)%set(1, points(i)%get(1)**2)
          call tensor_result(i)%set(1, 1, 2*points(i)%get(1))
          call tensor_result(i)%set(2, 1, 2*points(i)%get(1))
        enddo
      elseif(ops == 3) then
        call scalar_function_parser_1%create("x^2", num_dims=3)
        call scalar_function_parser_2%create("0",   num_dims=3)
        call scalar_function_parser_3%create("2*x", num_dims=3)
        call vector_function_parser%create(scalar_function_parser_1, scalar_function_parser_2, scalar_function_parser_2)
        call tensor_function_parser%create(scalar_function_parser_3, scalar_function_parser_2, scalar_function_parser_2, &
                                           scalar_function_parser_3, scalar_function_parser_2, scalar_function_parser_2, &
                                           scalar_function_parser_3, scalar_function_parser_2, scalar_function_parser_2)
        call vector_function_and_gradient_parser%create(vector_function_parser, tensor_function_parser)
        do i=1, num
          call vector_function_and_gradient_parser%get_value_space(points(i), vector_parser_result(i))
          call vector_function_and_gradient_parser%get_gradient_space(points(i), tensor_parser_result(i))
          call vector_result(i)%set(1, points(i)%get(1)**2)
          call tensor_result(i)%set(1, 1, 2*points(i)%get(1))
          call tensor_result(i)%set(2, 1, 2*points(i)%get(1))
          call tensor_result(i)%set(3, 1, 2*points(i)%get(1))
        enddo
      endif
      assert(all((/ (abs(vector_result(i)%get_value() - vector_parser_result(i)%get_value()) < error , i=1, num) /)))
      assert(all((/(((abs(tensor_result(i)%get(j,k)-tensor_parser_result(i)%get(j,k))<error,k=1,ops),j=1,ops),i=1,num)/)))
    endif

  call scalar_function_and_gradient_parser%free()
  call vector_function_and_gradient_parser%free()
  call scalar_function_parser%free()
  call scalar_function_parser_1%free()
  call scalar_function_parser_2%free()
  call scalar_function_parser_3%free()
  call scalar_function_parser_4%free()
  call scalar_function_parser_5%free()
  call scalar_function_parser_6%free()
  call scalar_function_parser_7%free()
  call scalar_function_parser_8%free()
  call scalar_function_parser_9%free()
  call vector_function_parser%free()
  call tensor_function_parser%free()

  if(allocated(parser_result))        deallocate(parser_result)
  if(allocated(result))               deallocate(result)
  if(allocated(vector_parser_result)) deallocate(vector_parser_result)
  if(allocated(vector_result))        deallocate(vector_result)
  if(allocated(tensor_parser_result)) deallocate(tensor_parser_result)
  if(allocated(tensor_result))        deallocate(tensor_result)

  call parameter_handler%free()
  call fempar_finalize()

contains

  subroutine define_test_function_parameters(this)
    class(fempar_parameter_handler_t), intent(inout) :: this
    type(parameterlist_t), pointer :: values, switches, switches_ab, helpers, required 
    integer(ip) :: error

    values      => this%get_values()
    switches    => this%get_switches()
    switches_ab => this%get_switches_ab()
    helpers     => this%get_helpers()

    error = 0
    error = error + helpers%set(key = 'num_ops',        Value= 'Dimension of the function to evaluate (2/3)')
    error = error + switches%set(key = 'num_ops',       Value= '--TEST_NUM_OPERATORS')
    error = error + values%set(key = 'num_ops',         Value= 2)

    error = error + helpers%set(key = 'function_type',  Value= 'Type of the function to evaluate (scalar)')
    error = error + switches%set(key = 'function_type', Value= '--TEST_FUNCTION_TYPE')
    error = error + values%set(key = 'function_type',   Value= 'scalar')

    error = error + helpers%set(key = 'num_points',     Value= 'Number of points to evaluate')
    error = error + switches%set(key = 'num_points',    Value= '--TEST_NUM_POINTS')
    error = error + values%set(key = 'num_points',      Value= 9999)

    error = error + helpers%set(key = 'operator',       Value= 'Operator of the function to evaluate')
    error = error + switches%set(key = 'operator',      Value= '--TEST_OPERATOR')
    error = error + values%set(key = 'operator',        Value= '+')

    error = error + helpers%set(key = 'test_name',      Value= 'Type of the function to evaluate (none/pointers)')
    error = error + switches%set(key = 'test_name',     Value= '--TEST_FUNCTION_NAME')
    error = error + values%set(key = 'test_name',       Value= 'none')

    assert(error == 0)
  end subroutine  define_test_function_parameters

  subroutine create_random_points(num, points)
    integer,                    intent(in)    :: num
    type(point_t), allocatable, intent(inout) :: points(:)
    real(rp),      allocatable                :: values(:,:)
    integer                                   :: i

    if(allocated(points)) deallocate(points)
    allocate(points(num))
    allocate(values(SPACE_DIM, num))
    CALL RANDOM_NUMBER(values)
    do i=1,num
        call points(i)%init(values(:,i))
    enddo
  end subroutine create_random_points

end program test_function_parser
