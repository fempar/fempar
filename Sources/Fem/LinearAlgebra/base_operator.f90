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
module base_operator_names
  use types
  use memory_guard_names
  use base_operand_names
  implicit none

  private

  ! Abstract operator (and its pure virtual function apply)
  type, abstract :: base_operator
     private
     logical :: allocated = .false.
   contains
     procedure :: is_allocated
     procedure :: set_allocated
     procedure :: free_expression
     procedure (apply_interface)       , deferred :: apply
     procedure (apply_fun_interface)   , deferred :: apply_fun
     procedure (info_interface)        , deferred :: info
     procedure (am_i_fine_task_interface), deferred :: am_i_fine_task
     procedure (bcast_interface)       , deferred :: bcast
     procedure  :: sum       => sum_operator_constructor
     procedure  :: sub       => sub_operator_constructor
     procedure  :: mult      => mult_operator_constructor
     procedure  :: minus     => minus_operator_constructor
     procedure, pass(op_left)  :: scal_left => scal_left_operator_constructor
     procedure, pass(op_right) :: scal_right => scal_right_operator_constructor
     generic    :: operator(+) => sum
     generic    :: operator(*) => mult, scal_right, scal_left, apply_fun
     generic    :: operator(-) => minus, sub
  end type base_operator

  ! Son class abstract operator
  type, extends(base_operator) :: abs_operator
     class(base_operator), pointer :: op => null()
   contains
     procedure  :: apply     => abs_operator_apply
     procedure  :: apply_fun => abs_operator_apply_fun
     procedure  :: info => abs_operator_info
     procedure  :: am_i_fine_task => abs_operator_am_i_fine_task
     procedure  :: bcast => abs_operator_bcast
     procedure  :: free_expression => abs_operator_free_expression
     procedure  :: abs  => abs_operator_constructor
     generic    :: assignment(=) => abs
  end type abs_operator

  ! Son class sum
  type, extends(base_operator) :: sum_operator
     class(base_operator), pointer :: op1 => null(), op2 => null()
   contains
     procedure  :: apply => sum_operator_apply
     procedure  :: apply_fun => sum_operator_apply_fun 
     procedure  :: info => sum_operator_info
     procedure  :: am_i_fine_task => sum_operator_am_i_fine_task
     procedure  :: bcast => sum_operator_bcast
     procedure  :: free_expression => sum_operator_free_expression
  end type sum_operator

  ! Son class sub
  type, extends(base_operator) :: sub_operator
     class(base_operator), pointer :: op1 => null(), op2 => null()
   contains
     procedure  :: apply => sub_operator_apply
     procedure  :: apply_fun => sub_operator_apply_fun 
     procedure  :: info => sub_operator_info
     procedure  :: am_i_fine_task => sub_operator_am_i_fine_task
     procedure  :: bcast => sub_operator_bcast
     procedure  :: free_expression => sub_operator_free_expression
  end type sub_operator

  ! Son class mult
  type, extends(base_operator) :: mult_operator
     class(base_operator), pointer :: op1 => null() , op2 => null()
   contains
     procedure  :: apply => mult_operator_apply
     procedure  :: apply_fun => mult_operator_apply_fun 
     procedure  :: info => mult_operator_info
     procedure  :: am_i_fine_task => mult_operator_am_i_fine_task
     procedure  :: bcast => mult_operator_bcast
     procedure  :: free_expression => mult_operator_free_expression
  end type mult_operator

  ! Son class scal
  type, extends(base_operator) :: scal_operator
     class(base_operator), pointer :: op => null()
     real(rp)                     :: alpha
   contains
     procedure  :: apply => scal_operator_apply
     procedure  :: apply_fun => scal_operator_apply_fun
     procedure  :: info => scal_operator_info
     procedure  :: am_i_fine_task => scal_operator_am_i_fine_task
     procedure  :: bcast => scal_operator_bcast
     procedure  :: free_expression => scal_operator_free_expression
  end type scal_operator

  ! Son class minus
  type, extends(base_operator) :: minus_operator
     class(base_operator), pointer :: op => null()
     real(rp)                     :: alpha
   contains
     procedure  :: apply => minus_operator_apply
     procedure  :: apply_fun => minus_operator_apply_fun 
     procedure  :: info => minus_operator_info
     procedure  :: am_i_fine_task => minus_operator_am_i_fine_task
     procedure  :: bcast => minus_operator_bcast
     procedure  :: free_expression => minus_operator_free_expression
  end type minus_operator

  ! Abstract interfaces
  abstract interface
     ! op%apply(x,y) <=> y <- op*x
     ! Implicitly assumes that y is already allocated
     subroutine apply_interface(op,x,y) 
       import :: base_operator, base_operand
       implicit none
       class(base_operator), intent(in)    :: op
       class(base_operand) , intent(in)    :: x
       class(base_operand) , intent(inout) :: y 
     end subroutine apply_interface

     ! op%apply(x)
     ! Allocates room for (temporary) y
     function apply_fun_interface(op,x) result(y)
       import :: base_operator, base_operand
       implicit none
       class(base_operator), intent(in)  :: op
       class(base_operand) , intent(in)  :: x
       class(base_operand) , allocatable :: y 
     end function apply_fun_interface

     subroutine info_interface(op,me,np) 
       import :: base_operator, ip
       implicit none
       class(base_operator), intent(in)  :: op
       integer(ip)         , intent(out) :: me
       integer(ip)         , intent(out) :: np
     end subroutine info_interface

     function am_i_fine_task_interface(op) 
       import :: base_operator, ip
       implicit none
       class(base_operator), intent(in)  :: op
       logical                           :: am_i_fine_task_interface 
     end function am_i_fine_task_interface

     subroutine bcast_interface (op, condition)
       import :: base_operator
       implicit none
       class(base_operator), intent(in) :: op
       logical, intent(inout) :: condition
     end subroutine bcast_interface

  end interface

  public :: abs_operator, base_operator

contains
  
  subroutine free_expression(this)
    implicit none
    class(base_operator), intent(inout) :: this
  end subroutine free_expression

  function is_allocated(op)
    implicit none
    class(base_operator), intent(in) :: op
    logical :: is_allocated
    is_allocated = op%allocated
  end function is_allocated

  subroutine set_allocated(op)
    implicit none
    class(base_operator), intent(inout) :: op
    op%allocated = .true.
  end subroutine set_allocated

  !--------------------------------------------------------------------!
  ! Construction and deallocation functions/subroutines of the nodes of! 
  ! the tree that represents an expression among matrix operators      !
  ! -------------------------------------------------------------------!
  function sum_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator)    , intent(in), target :: op1, op2
    type(sum_operator), pointer             :: res
    allocate(res)
    res%op1 => op1
    res%op2 => op2
    call res%set_allocated()
  end function sum_operator_constructor

  subroutine sum_operator_free_expression(this)
    implicit none
    class(sum_operator), intent(inout) :: this 
    call this%op1%free_expression()
    if ( this%op1%is_allocated() ) deallocate(this%op1)   
    call this%op2%free_expression()
    if ( this%op2%is_allocated() ) deallocate(this%op2)   
  end subroutine sum_operator_free_expression

  function sub_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator)    , intent(in), target :: op1, op2
    type(sub_operator), pointer             :: res
    allocate(res)
    res%op1 => op1
    res%op2 => op2
    call res%set_allocated()
  end function sub_operator_constructor

  subroutine sub_operator_free_expression(this)
    implicit none
    class(sub_operator), intent(inout) :: this 
    call this%op1%free_expression()
    if ( this%op1%is_allocated() ) deallocate(this%op1)   
    call this%op2%free_expression()
    if ( this%op2%is_allocated() ) deallocate(this%op2)
  end subroutine sub_operator_free_expression

  function mult_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator)    , intent(in), target :: op1, op2
    type(mult_operator), pointer            :: res
    allocate(res)
    res%op1 => op1
    res%op2 => op2
    call res%set_allocated()
  end function mult_operator_constructor

  subroutine mult_operator_free_expression(this) 
    implicit none
    class(mult_operator), intent(inout) :: this
    call this%op1%free_expression()
    if ( this%op1%is_allocated() ) deallocate(this%op1)   
    call this%op2%free_expression()
    if ( this%op2%is_allocated() ) deallocate(this%op2)
  end subroutine mult_operator_free_expression

  function scal_left_operator_constructor(alpha, op_left) result (res)
    implicit none
    class(base_operator)     , intent(in), target :: op_left
    real(rp)                , intent(in)          :: alpha
    type(scal_operator), pointer             :: res
    allocate(res)
    res%op => op_left
    res%alpha = alpha
    call res%set_allocated()
  end function scal_left_operator_constructor

  subroutine scal_operator_free_expression(this)
    implicit none
    class(scal_operator), intent(inout) :: this
    call this%op%free_expression()
    if ( this%op%is_allocated() ) deallocate(this%op)   
  end subroutine scal_operator_free_expression

  function scal_right_operator_constructor(op_right,alpha) result (res)
    implicit none
    class(base_operator)     , intent(in), target  :: op_right
    real(rp)                , intent(in)           :: alpha
    class(scal_operator), pointer              :: res
    allocate(res)
    res%op => op_right
    res%alpha = alpha
    call res%set_allocated()
  end function scal_right_operator_constructor

  function minus_operator_constructor(op) result (res)
    implicit none
    class(base_operator), intent(in), target :: op
    type(minus_operator), pointer            :: res
    allocate(res)
    res%op => op
    call res%set_allocated()
  end function minus_operator_constructor

  subroutine minus_operator_free_expression(this)
    implicit none
    class(minus_operator), intent(inout) :: this
    call this%op%free_expression()
    if ( this%op%is_allocated() ) deallocate(this%op)  
  end subroutine minus_operator_free_expression

  subroutine abs_operator_constructor(res,op)
    implicit none
    class(base_operator), intent(in), target :: op
    class(abs_operator), intent(inout)      :: res
    res%op => op
  end subroutine abs_operator_constructor

  subroutine abs_operator_free_expression(this)
    implicit none
    class(abs_operator), intent(inout) :: this
    call this%op%free_expression()
    if ( this%op%is_allocated() ) deallocate(this%op)  
  end subroutine abs_operator_free_expression

  !-------------------------------------!
  ! apply_fun and apply implementations !
  !-------------------------------------!
  function sum_operator_apply_fun(op,x) result(y)
    implicit none
    class(sum_operator), intent(in)       :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y = op%op2 * x
    call y%axpby( 1.0, op%op1*x, 1.0 )
    call x%CleanTemp()
    call y%SetTemp()
  end function sum_operator_apply_fun

  subroutine sum_operator_apply(op,x,y)
    implicit none
    class(sum_operator), intent(in)    :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 

    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    ! y <- 1.0 * op1*x + 1.0*y
    call y%axpby( 1.0, op%op1*x, 1.0 )
    call x%CleanTemp()
  end subroutine sum_operator_apply
  
  function sub_operator_apply_fun(op,x) result(y)
    implicit none
    class(sub_operator), intent(in)       :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y = op%op2 * x
    call y%axpby( 1.0, op%op1*x, -1.0 )
    call x%CleanTemp()
    call y%SetTemp()
  end function sub_operator_apply_fun

  subroutine sub_operator_apply(op,x,y)
    implicit none
    class(sub_operator), intent(in)    :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 
    
    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    ! y <- 1.0 * op1*x - 1.0*y
    call y%axpby( 1.0, op%op1*x, -1.0 )
    call x%CleanTemp()
  end subroutine sub_operator_apply
  
  function mult_operator_apply_fun(op,x) result(y)
    implicit none
    class(mult_operator), intent(in)       :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y = op%op1 * ( op%op2 * x )
    call x%CleanTemp()
    call y%SetTemp()
  end function mult_operator_apply_fun

  subroutine mult_operator_apply(op,x,y)
    implicit none
    class(mult_operator), intent(in)    :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 

    call x%GuardTemp()
    call op%op2%apply ( op%op1*x, y )
    call x%CleanTemp()
  end subroutine mult_operator_apply

  function  scal_operator_apply_fun(op,x) result(y)
    implicit none
    class(scal_operator), intent(in)      :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y =  op%op * x
    call y%scal ( op%alpha, y)
    call x%CleanTemp()
    call y%SetTemp()
  end function scal_operator_apply_fun

  subroutine scal_operator_apply(op,x,y)
    implicit none
    class(scal_operator), intent(in)   :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 

    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( op%alpha, y )
    call x%CleanTemp()
  end subroutine scal_operator_apply

  function  minus_operator_apply_fun(op,x) result(y)
    implicit none
    class(minus_operator), intent(in)       :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y = op%op*x
    call y%scal( -1.0, y )
    call x%CleanTemp()
    call y%SetTemp()
  end function minus_operator_apply_fun

  subroutine minus_operator_apply(op,x,y)
    implicit none
    class(minus_operator), intent(in)  :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 

    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( -1.0, y )
    call x%CleanTemp()
  end subroutine minus_operator_apply

  function  abs_operator_apply_fun(op,x) result(y)
    implicit none
    class(abs_operator), intent(in)       :: op
    class(base_operand)     , intent(in)  :: x
    class(base_operand)     , allocatable :: y 

    call x%GuardTemp()
    allocate(y, mold=x)
    y = op%op*x
    call x%CleanTemp()
    call y%SetTemp()
  end function abs_operator_apply_fun

  subroutine abs_operator_apply(op,x,y)
    implicit none
    class(abs_operator), intent(in)    :: op
    class(base_operand), intent(in)    :: x
    class(base_operand), intent(inout) :: y 

    call x%GuardTemp()
    call op%op%apply(x,y)
    call x%CleanTemp()
  end subroutine abs_operator_apply

  !-------------------------------------!
  ! info implementations                !
  !-------------------------------------!
  subroutine sum_operator_info(op,me,np)
    implicit none
    class(sum_operator), intent(in)    :: op
    integer(ip)        , intent(out)   :: me
    integer(ip)        , intent(out)   :: np

    call op%op1%info(me,np)
  end subroutine sum_operator_info

  subroutine sub_operator_info(op,me,np)
    implicit none
    class(sub_operator), intent(in)    :: op
    integer(ip)        , intent(out)   :: me
    integer(ip)        , intent(out)   :: np

    call op%op1%info(me,np)
  end subroutine sub_operator_info
  
  subroutine mult_operator_info(op,me,np)
    implicit none
    class(mult_operator), intent(in)    :: op
    integer(ip)         , intent(out)   :: me
    integer(ip)         , intent(out)   :: np

    call op%op1%info(me,np)
  end subroutine mult_operator_info

  subroutine minus_operator_info(op,me,np)
    implicit none
    class(minus_operator), intent(in)   :: op
    integer(ip)         , intent(out)   :: me
    integer(ip)         , intent(out)   :: np

    call op%op%info(me,np)
  end subroutine minus_operator_info
  
  subroutine scal_operator_info(op,me,np)
    implicit none
    class(scal_operator), intent(in)    :: op
    integer(ip)         , intent(out)   :: me
    integer(ip)         , intent(out)   :: np

    call op%op%info(me,np)
  end subroutine scal_operator_info

  subroutine abs_operator_info(op,me,np)
    implicit none
    class(abs_operator), intent(in)    :: op
    integer(ip)        , intent(out)   :: me
    integer(ip)        , intent(out)   :: np
    
    call op%op%info(me,np)
  end subroutine abs_operator_info

  !-------------------------------------!
  ! am_i_fine_task implementations      !
  !-------------------------------------!
  function sum_operator_am_i_fine_task(op)
    implicit none
    class(sum_operator), intent(in)    :: op
    logical :: sum_operator_am_i_fine_task
    
    sum_operator_am_i_fine_task = op%op1%am_i_fine_task()
  end function sum_operator_am_i_fine_task
  
  function sub_operator_am_i_fine_task(op)
    implicit none
    class(sub_operator), intent(in)    :: op
    logical :: sub_operator_am_i_fine_task
    
    sub_operator_am_i_fine_task = op%op1%am_i_fine_task()
  end function sub_operator_am_i_fine_task

  function mult_operator_am_i_fine_task(op)
    implicit none
    class(mult_operator), intent(in)    :: op
    logical :: mult_operator_am_i_fine_task
    
    mult_operator_am_i_fine_task = op%op1%am_i_fine_task()
  end function mult_operator_am_i_fine_task
  
  function minus_operator_am_i_fine_task(op)
    implicit none
    class(minus_operator), intent(in)    :: op
    logical :: minus_operator_am_i_fine_task
    
    minus_operator_am_i_fine_task = op%op%am_i_fine_task()
  end function minus_operator_am_i_fine_task

  function scal_operator_am_i_fine_task(op)
    implicit none
    class(scal_operator), intent(in)    :: op
    logical :: scal_operator_am_i_fine_task
    
    scal_operator_am_i_fine_task = op%op%am_i_fine_task()
  end function scal_operator_am_i_fine_task

  function abs_operator_am_i_fine_task(op)
    implicit none
    class(abs_operator), intent(in)    :: op
    logical :: abs_operator_am_i_fine_task
    
    abs_operator_am_i_fine_task = op%op%am_i_fine_task()
  end function abs_operator_am_i_fine_task

  !-------------------------------------!
  ! bcast implementations               !
  !-------------------------------------!
  subroutine sum_operator_bcast(op,condition)
    implicit none
    class(sum_operator), intent(in)    :: op
    logical            , intent(inout)   :: condition

    call op%op1%bcast(condition)
  end subroutine sum_operator_bcast

  subroutine sub_operator_bcast(op,condition)
    implicit none
    class(sub_operator), intent(in)    :: op
    logical            , intent(inout)   :: condition

    call op%op1%bcast(condition)
  end subroutine sub_operator_bcast

  subroutine mult_operator_bcast(op,condition)
    implicit none
    class(mult_operator), intent(in)   :: op
    logical            , intent(inout)   :: condition

    call op%op1%bcast(condition)
  end subroutine mult_operator_bcast

  subroutine minus_operator_bcast(op,condition)
    implicit none
    class(minus_operator), intent(in)   :: op
    logical            , intent(inout)   :: condition

    call op%op%bcast(condition)
  end subroutine minus_operator_bcast

  subroutine scal_operator_bcast(op,condition)
    implicit none
    class(scal_operator), intent(in)   :: op
    logical            , intent(inout)   :: condition

    call op%op%bcast(condition)
  end subroutine scal_operator_bcast

  subroutine abs_operator_bcast(op,condition)
    implicit none
    class(abs_operator), intent(in)     :: op
    logical            , intent(inout)   :: condition

    call op%op%bcast(condition)
  end subroutine abs_operator_bcast

end module base_operator_names
