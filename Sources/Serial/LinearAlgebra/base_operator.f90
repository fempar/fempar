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
  use types_names
  use memory_guard_names
  use integrable_names
  use base_operand_names
  implicit none
# include "debug.i90"

  private

  ! Abstract operator (and its pure virtual function apply)
  type, abstract, extends(integrable_t) :: base_operator_t
     integer(ip) :: fill_values_stage = update_nonlinear
     integer(ip) :: free_values_stage = update_nonlinear
   contains
     procedure (apply_interface)         , deferred :: apply
     procedure (apply_fun_interface)     , deferred :: apply_fun
     procedure (fill_values_interface)   , deferred :: fill_values
     procedure (free_values_interface)   , deferred :: free_values
     procedure  :: sum       => sum_operator_constructor
     procedure  :: sub       => sub_operator_constructor
     procedure  :: mult      => mult_operator_constructor
     procedure  :: minus     => minus_operator_constructor
     procedure, pass(op_left)  :: scal_left => scal_left_operator_constructor
     procedure, pass(op_right) :: scal_right => scal_right_operator_constructor
     generic    :: operator(+) => sum
     generic    :: operator(*) => mult, scal_right, scal_left, apply_fun
     generic    :: operator(-) => minus, sub
  end type base_operator_t

  ! Son class expression_operator_t. These operators are always temporary
  ! and therefore an assignment is needed to make copies. The gfortran
  ! compiler only supports A=B when A and B are polymorphic if the assignment 
  ! is overwritten.
  type, abstract, extends(base_operator_t) :: expression_operator_t 
   contains
     procedure (expression_operator_assign_interface), deferred :: assign
     generic  :: assignment(=) => assign
  end type expression_operator_t

  ! Derived class binary
  type, abstract, extends(expression_operator_t) :: binary_operator_t
     class(base_operator_t), pointer :: op1 => null(), op2 => null()
   contains
     procedure :: default_initialization => binary_operator_default_init
     procedure :: fill_values => binary_operator_fill_values
     procedure :: free_values => binary_operator_free_values
     procedure :: free    => binary_operator_destructor
     procedure :: assign  => binary_operator_copy
  end type binary_operator_t


  ! Son class abstract operator
  type, extends(base_operator_t) :: abs_operator_t
     class(base_operator_t), pointer :: op_stored => null()
     class(base_operator_t), pointer :: op        => null()
   contains
     procedure  :: default_initialization => abs_operator_default_init
     procedure  :: apply     => abs_operator_apply
     procedure  :: apply_fun => abs_operator_apply_fun
     procedure  :: fill_values => abs_operator_fill_values
     procedure  :: free_values => abs_operator_free_values
     procedure  :: free  => abs_operator_destructor
     procedure  :: assign => abs_operator_constructor
     generic    :: assignment(=) => assign
  end type abs_operator_t

  ! Son class sum
  type, extends(binary_operator_t) :: sum_operator_t
  contains
     procedure  :: apply => sum_operator_apply
     procedure  :: apply_fun => sum_operator_apply_fun 
  end type sum_operator_t


  ! Son class sub
  type, extends(binary_operator_t) :: sub_operator_t
   contains
     procedure  :: apply => sub_operator_apply
     procedure  :: apply_fun => sub_operator_apply_fun 
  end type sub_operator_t


  ! Son class mult
  type, extends(binary_operator_t) :: mult_operator_t
   contains
     procedure  :: apply => mult_operator_apply
     procedure  :: apply_fun => mult_operator_apply_fun 
  end type mult_operator_t


  ! Son class scal
  type, extends(expression_operator_t) :: scal_operator_t
     class(base_operator_t), pointer :: op => null()
     real(rp)                     :: alpha
   contains
     procedure  :: default_initialization => scal_operator_default_init
     procedure  :: apply => scal_operator_apply
     procedure  :: apply_fun => scal_operator_apply_fun
     procedure  :: fill_values => scal_operator_fill_values
     procedure  :: free_values => scal_operator_free_values
     procedure  :: free => scal_operator_destructor
     procedure  :: assign => scal_operator_copy
  end type scal_operator_t


  ! Son class minus
  type, extends(base_operator_t) :: minus_operator_t
     class(base_operator_t), pointer :: op => null()
   contains
     procedure  :: default_initialization => minus_operator_default_init
     procedure  :: apply => minus_operator_apply
     procedure  :: apply_fun => minus_operator_apply_fun 
     procedure  :: fill_values => minus_operator_fill_values
     procedure  :: free_values => minus_operator_free_values
     procedure  :: free => minus_operator_destructor
     procedure  :: assign => minus_operator_copy
  end type minus_operator_t

  ! Abstract interfaces
  abstract interface
     ! op%apply(x,y) <=> y <- op*x
     ! Implicitly assumes that y is already allocated
     subroutine apply_interface(op,x,y) 
       import :: base_operator_t, base_operand_t
       implicit none
       class(base_operator_t), intent(in)    :: op
       class(base_operand_t) , intent(in)    :: x
       class(base_operand_t) , intent(inout) :: y 
     end subroutine apply_interface

     ! op%apply(x)
     ! Allocates room for (temporary) y
     function apply_fun_interface(op,x) result(y)
       import :: base_operator_t, base_operand_t
       implicit none
       class(base_operator_t), intent(in)  :: op
       class(base_operand_t) , intent(in)  :: x
       class(base_operand_t) , allocatable :: y 
     end function apply_fun_interface

     ! op1%fill_values()
     ! Fill preconditioner values
     subroutine fill_values_interface(op,stage)
       import :: base_operator_t,ip
       implicit none
       class(base_operator_t), intent(inout) :: op
       integer(ip), optional , intent(in)    :: stage
     end subroutine fill_values_interface

     ! op1%free_values()
     ! Free preconditioner values
     subroutine free_values_interface(op,stage)
       import :: base_operator_t,ip
       implicit none
       class(base_operator_t), intent(inout) :: op
       integer(ip), optional , intent(in)    :: stage
     end subroutine free_values_interface

     subroutine expression_operator_assign_interface(op1,op2)
       import :: base_operator_t, expression_operator_t
       implicit none
       class(base_operator_t)      , intent(in)    :: op2
       class(expression_operator_t), intent(inout) :: op1
     end subroutine expression_operator_assign_interface

  end interface

  public :: abs_operator_t, base_operator_t

contains

  subroutine binary_operator_default_init(this)
    implicit none
    class(binary_operator_t), intent(inout) :: this
    nullify(this%op1)
    nullify(this%op2)
    call this%NullifyTemporary()
  end subroutine binary_operator_default_init

  subroutine binary_operator_destructor(this)
    implicit none
    class(binary_operator_t), intent(inout) :: this
    ! out_verbosity10_2( 'Destructing binary',this%id )

    select type(that => this%op1)
    class is(expression_operator_t)
       call that%CleanTemp()
    class is(abs_operator_t)
       call that%CleanTemp()
    class default
       check(1==0)
       ! out_verbosity10_1( 'Why I am here?' )
    end select
    ! out_verbosity10_2( 'Deallocating',this%op1%id )
    deallocate(this%op1)

    select type(that => this%op2)
    class is(expression_operator_t)
       call that%CleanTemp()
    class is(abs_operator_t)
       call that%CleanTemp()
    class default
       check(1==0)
       ! out_verbosity10_1( 'Why I am here?')
    end select
    ! out_verbosity10_2( 'Deallocating',this%op2%id )
    deallocate(this%op2)
  end subroutine binary_operator_destructor

  subroutine binary_operator_copy(op1,op2)
    implicit none
    class(base_operator_t)  , intent(in)    :: op2
    class(binary_operator_t), intent(inout) :: op1

    ! global_id=global_id+1
    ! op1%id=global_id
    ! out_verbosity10_4( 'Copy binary_operator', op1%id, ' from ',op2%id)
    select type(op2)
    class is(binary_operator_t)
       call binary_operator_constructor(op2%op1,op2%op2,op1)
    class default
       check(1==0)
    end select
  end subroutine binary_operator_copy

  subroutine binary_operator_constructor(op1,op2,res) 
    implicit none
    class(base_operator_t)  , intent(in)    :: op1, op2
    class(binary_operator_t), intent(inout) :: res

    call op1%GuardTemp()
    call op2%GuardTemp()

    ! out_verbosity10_1( 'Creating binary left operator')

    ! Allocate op1
    select type(op1)
    class is(expression_operator_t)
       allocate(res%op1,mold=op1); call res%op1%default_initialization()
    class default
       allocate(abs_operator_t::res%op1)
    end select
    ! Assign op1
    select type(this => res%op1)
    class is(expression_operator_t)
       ! out_verbosity10_1( 'binary left is an expression' )
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
    class is(abs_operator_t)
       ! out_verbosity10_1( 'binary left is not an expression' )
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select

    ! out_verbosity10_1( 'Creating binary right operator')

    ! Allocate op2
    select type(op2)
    class is(expression_operator_t)
       allocate(res%op2,mold=op2); call res%op2%default_initialization()
    class default
       allocate(abs_operator_t::res%op2)
    end select
    ! Assign op2
    select type(that => res%op2)
    class is(expression_operator_t)
       ! out_verbosity10_1( 'binary right is an expression' )
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%GuardTemp()
    class is(abs_operator_t)
       ! out_verbosity10_1( 'binary right is not an expression' )
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%SetTemp()
       call that%GuardTemp()
    end select
    call res%SetTemp()
    call op1%CleanTemp()
    call op2%CleanTemp()
  end subroutine binary_operator_constructor

  subroutine abs_operator_default_init(this)
    implicit none
    class(abs_operator_t), intent(inout) :: this
    nullify(this%op)
    nullify(this%op_stored)
    call this%NullifyTemporary()
  end subroutine abs_operator_default_init

  recursive subroutine abs_operator_constructor(op1,op2)
    implicit none
    class(abs_operator_t) , intent(inout) :: op1
    class(base_operator_t), intent(in), target  :: op2

    call op1%free()
    ! global_id=global_id+1
    ! op1%id=global_id

    call op2%GuardTemp()
    select type(op2)
    class is(abs_operator_t) ! Can be temporary (or not)
       if(associated(op2%op_stored)) then
          assert(.not.associated(op2%op))
          allocate(op1%op_stored, mold = op2%op_stored); call op1%op_stored%default_initialization()
          select type(this => op1%op_stored)
          class is(expression_operator_t)
             ! out_verbosity10_4( 'Creating abs ', op1%id, ' from abs (copy content, which is expression)', op2%id)
             this = op2%op_stored
          class is(abs_operator_t)
             ! out_verbosity10_4( 'Creating abs ', op1%id, ' from abs (copy content, which is abs)', op2%id)
             this = op2%op_stored
          class default
             ! out_verbosity10_1( 'How this is possible????')
             check(1==0)
          end select
          call op1%op_stored%GuardTemp()
       else if(associated(op2%op)) then
          assert(.not.associated(op2%op_stored))
          ! out_verbosity10_4( 'Creating abs ', op1%id, ' from abs (reassign pointer)', op2%id)
          op1%op => op2%op
          call op1%op%GuardTemp()
       else
          ! out_verbosity10_1( 'How this is possible????')
          check(1==0)
       end if
    class is(expression_operator_t) ! Temporary
       ! out_verbosity10_2( 'Creating abs from expression (copy expression)', op1%id)
       allocate(op1%op_stored,mold=op2); call op1%op_stored%default_initialization()
       select type(this => op1%op_stored)
       class is(expression_operator_t)
          this = op2              ! Here = overloaded
       end select
       call op1%op_stored%GuardTemp()
   class default                 ! Cannot be temporary (I don't know how to copy it!)
       !assert(.not.op2%IsTemp())
       !out_verbosity10_2( 'Creating abs from base_operator (point to a permanent base_operator)', op1%id)
       op1%op => op2
       call op1%op%GuardTemp()
    end select
    ! out_verbosity10_2( 'Creating abs rhs argument has temporary counter', op2%GetTemp())
    call op2%CleanTemp()
  end subroutine abs_operator_constructor

  subroutine abs_operator_destructor(this)
    implicit none
    class(abs_operator_t), intent(inout) :: this

    if(associated(this%op)) then
       assert(.not.associated(this%op_stored))
       ! out_verbosity10_2( 'Destructing abs association', this%id)
       ! Nothing to free, the pointer points to permanent data
       this%op => null()
    else if(associated(this%op_stored)) then
       assert(.not.associated(this%op))
       ! out_verbosity10_2( 'Destructing abs allocation', this%id)
       ! Any of the following two lines should give the same result
       call this%op_stored%CleanTemp()
       ! out_verbosity10_2( 'Deallocating',this%op_stored%id )
       deallocate(this%op_stored)
    end if
    ! else
       ! out_verbosity10_2( 'Called when initialized',this%id)
    ! end if
  end subroutine abs_operator_destructor

  subroutine scal_operator_default_init(this)
    implicit none
    class(scal_operator_t), intent(inout) :: this
    nullify(this%op)
    call this%NullifyTemporary()
  end subroutine scal_operator_default_init

  subroutine scal_operator_constructor(alpha,op,res)
    implicit none
    class(base_operator_t), intent(in)    :: op
    real(rp)            , intent(in)    :: alpha
    type(scal_operator_t) , intent(inout) :: res

    call op%GuardTemp()
    res%alpha = alpha
    ! Allocate op
    select type(op)
    class is(expression_operator_t)
       allocate(res%op,mold=op); call res%op%default_initialization()
    class default
       allocate(abs_operator_t::res%op)
    end select
    ! Assign op
    select type(this => res%op)
    class is(expression_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
    class is(abs_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select
    call res%SetTemp()
    call op%CleanTemp()
  end subroutine scal_operator_constructor

  subroutine scal_operator_destructor(this)
    implicit none
    class(scal_operator_t), intent(inout) :: this 
    ! out_verbosity10_2( 'Destructing scal ', this%id)
    select type(that => this%op)
    class is(expression_operator_t)
       call that%CleanTemp()
    class is(abs_operator_t)
       call that%CleanTemp()
    class default
       ! out_verbosity10_1( 'Why I am here?' )
       check(1==0)
    end select
    ! out_verbosity10_2( 'Deallocating',this%op%id )
    deallocate(this%op)
  end subroutine scal_operator_destructor

  subroutine scal_operator_copy(op1,op2)
    implicit none
    class(base_operator_t), intent(in)    :: op2
    class(scal_operator_t), intent(inout) :: op1
    !global_id=global_id+1
    !op1%id=global_id
    !out_verbosity10_4( 'Copy scal_operator ',op1%id,' from ', op2%id)
    select type(op2)
    class is(scal_operator_t)
       call scal_operator_constructor(op2%alpha,op2%op,op1) ! Not the default constructor
    class default
       check(1==0)
       !out_verbosity10_1( 'Error assigning scal operators')
       !stop
    end select
  end subroutine scal_operator_copy

  subroutine minus_operator_default_init(this)
    implicit none
    class(minus_operator_t), intent(inout) :: this
    nullify(this%op)
    call this%NullifyTemporary()
  end subroutine minus_operator_default_init

  function minus_operator_constructor(op) result (res)
    implicit none
    class(base_operator_t)    , intent(in)  :: op
    type(minus_operator_t) :: res
    !type(scal_operator_t), allocatable :: res
    !allocate(res)
    !global_id=global_id+1
    !res%id=global_id
    !out_verbosity10_2('Creating scal left', res%id )
    call minus_operator_constructor_sub(op,res)
  end function minus_operator_constructor

  subroutine minus_operator_constructor_sub(op,res)
    implicit none
    class(base_operator_t) , intent(in)    :: op
    type(minus_operator_t) , intent(inout) :: res

    call op%GuardTemp()
    ! Allocate op
    select type(op)
    class is(expression_operator_t)
       allocate(res%op,mold=op); call res%op%default_initialization()
    class default
       allocate(abs_operator_t::res%op)
    end select
    ! Assign op
    select type(this => res%op)
    class is(expression_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
    class is(abs_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select
    call res%SetTemp()
    call op%CleanTemp()
  end subroutine minus_operator_constructor_sub

  subroutine minus_operator_destructor(this)
    implicit none
    class(minus_operator_t), intent(inout) :: this 
    ! out_verbosity10_2( 'Destructing minus ', this%id)
    select type(that => this%op)
    class is(expression_operator_t)
       call that%CleanTemp()
    class is(abs_operator_t)
       call that%CleanTemp()
    class default
       ! out_verbosity10_1( 'Why I am here?' )
       check(1==0)
    end select
    ! out_verbosity10_2( 'Deallocating',this%op%id )
    deallocate(this%op)
  end subroutine minus_operator_destructor

  subroutine minus_operator_copy(op1,op2)
    implicit none
    class(base_operator_t), intent(in)    :: op2
    class(minus_operator_t), intent(inout) :: op1
    !global_id=global_id+1
    !op1%id=global_id
    !out_verbosity10_4( 'Copy minus_operator ',op1%id,' from ', op2%id)
    select type(op2)
    class is(minus_operator_t)
       call minus_operator_constructor_sub(op2%op,op1) ! Not the default constructor
    class default
       check(1==0)
       !out_verbosity10_1( 'Error assigning minus operators')
       !stop
    end select
  end subroutine minus_operator_copy

!!$  !--------------------------------------------------------------------!
!!$  ! Construction and deallocation functions/subroutines of the nodes of! 
!!$  ! the tree that represents an expression among matrix operators      !
!!$  ! -------------------------------------------------------------------!
  function sum_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator_t), intent(in)  :: op1, op2
    type(sum_operator_t)  :: res
    !type(sum_operator_t)  , allocatable :: res
    !allocate(res)
    ! global_id=global_id+1
    ! res%id=global_id
    ! out_verbosity10_2( 'Creating sum', res%id)
    call binary_operator_constructor(op1,op2,res) 
  end function sum_operator_constructor

  function sub_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator_t), intent(in)  :: op1, op2
    type(sub_operator_t)  :: res
    !type(sub_operator_t)  , allocatable :: res
    !allocate(res)
    ! global_id=global_id+1
    ! res%id=global_id
    ! out_verbosity10_2( 'Creating sub', res%id)
    call binary_operator_constructor(op1,op2,res) 
  end function sub_operator_constructor

  function mult_operator_constructor(op1,op2) result (res)
    implicit none
    class(base_operator_t), intent(in)  :: op1, op2
    type(mult_operator_t) :: res
    !type(mult_operator_t) , allocatable :: res
    !allocate(res)
    !global_id=global_id+1
    !res%id=global_id
    !out_verbosity10_2( 'Creating mult', res%id)
    call binary_operator_constructor(op1,op2,res)
  end function mult_operator_constructor

  function scal_left_operator_constructor(alpha, op_left) result (res)
    implicit none
    class(base_operator_t)    , intent(in)  :: op_left
    real(rp)                , intent(in)  :: alpha
    type(scal_operator_t) :: res
    !type(scal_operator_t), allocatable :: res
    !allocate(res)
    !global_id=global_id+1
    !res%id=global_id
    !out_verbosity10_2('Creating scal left', res%id )
    call scal_operator_constructor(alpha,op_left,res)
  end function scal_left_operator_constructor
  
  function scal_right_operator_constructor(op_right, alpha) result (res)
    implicit none
    class(base_operator_t)    , intent(in)  :: op_right
    real(rp)                , intent(in)  :: alpha
    type(scal_operator_t) :: res
    !type(scal_operator_t), allocatable :: res
    !allocate(res)
    !global_id=global_id+1
    !res%id=global_id
    !out_verbosity10_2( 'Creating scal right', res%id )
    call scal_operator_constructor(alpha,op_right,res)
  end function scal_right_operator_constructor

  !-------------------------------------!
  ! apply_fun and apply implementations !
  !-------------------------------------!
  function sum_operator_apply_fun(op,x) result(y)
    implicit none
    class(sum_operator_t), intent(in)       :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    y = op%op2 * x
    call y%axpby( 1.0, op%op1*x, 1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
    call y%SetTemp()
  end function sum_operator_apply_fun

  subroutine sum_operator_apply(op,x,y)
    implicit none
    class(sum_operator_t), intent(in)    :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    ! y <- 1.0 * op1*x + 1.0*y
    call y%axpby( 1.0, op%op1*x, 1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine sum_operator_apply

  
  function sub_operator_apply_fun(op,x) result(y)
    implicit none
    class(sub_operator_t), intent(in)       :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    y = op%op2 * x
    call y%axpby( 1.0, op%op1*x, -1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
    call y%SetTemp()
  end function sub_operator_apply_fun

  subroutine sub_operator_apply(op,x,y)
    implicit none
    class(sub_operator_t), intent(in)    :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    ! y <- 1.0 * op1*x - 1.0*y
    call y%axpby( 1.0, op%op1*x, -1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine sub_operator_apply
  
  function mult_operator_apply_fun(op,x) result(y)
    implicit none
    class(mult_operator_t), intent(in)       :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    y = op%op1 * ( op%op2 * x )
    call x%CleanTemp()
    call op%CleanTemp()
    call y%SetTemp()
  end function mult_operator_apply_fun

  subroutine mult_operator_apply(op,x,y)
    implicit none
    class(mult_operator_t), intent(in)    :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op2%apply ( op%op1*x, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine mult_operator_apply

  function  scal_operator_apply_fun(op,x) result(y)
    implicit none
    class(scal_operator_t), intent(in)      :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    y =  op%op * x
    call y%scal ( op%alpha, y)
    call x%CleanTemp()
    call op%CleanTemp()
    call y%SetTemp()
  end function scal_operator_apply_fun

  subroutine scal_operator_apply(op,x,y)
    implicit none
    class(scal_operator_t), intent(in)   :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( op%alpha, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine scal_operator_apply

  function  minus_operator_apply_fun(op,x) result(y)
    implicit none
    class(minus_operator_t), intent(in)       :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    y = op%op*x
    call y%scal( -1.0, y )
    call x%CleanTemp()
    call y%SetTemp()
    call op%CleanTemp()
  end function minus_operator_apply_fun

  subroutine minus_operator_apply(op,x,y)
    implicit none
    class(minus_operator_t), intent(in)  :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( -1.0, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine minus_operator_apply

  function  abs_operator_apply_fun(op,x) result(y)
    implicit none
    class(abs_operator_t), intent(in)       :: op
    class(base_operand_t)     , intent(in)  :: x
    class(base_operand_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()

    if(associated(op%op_stored)) then
       assert(.not.associated(op%op))
       y = op%op_stored*x
    else if(associated(op%op)) then
       assert(.not.associated(op%op_stored))
       y = op%op*x
    else
       check(1==0)
    end if

    call x%CleanTemp()
    call op%CleanTemp()
    call y%SetTemp()
  end function abs_operator_apply_fun

  subroutine abs_operator_apply(op,x,y)
    implicit none
    class(abs_operator_t), intent(in)    :: op
    class(base_operand_t), intent(in)    :: x
    class(base_operand_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()

    if(associated(op%op_stored)) then
       assert(.not.associated(op%op))
       call op%op_stored%apply(x,y)
    else if(associated(op%op)) then
       assert(.not.associated(op%op_stored))
       call op%op%apply(x,y)
    else
       check(1==0)
    end if

    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine abs_operator_apply

  !-----------------------------!
  ! fill_values implementations !
  !-----------------------------!

  subroutine binary_operator_fill_values(op,stage)
    implicit none
    class(binary_operator_t), intent(inout) :: op
    integer(ip), optional   , intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage
    
    call op%op1%fill_values(stage_)
    call op%op2%fill_values(stage_)

  end subroutine binary_operator_fill_values

  subroutine abs_operator_fill_values(op,stage)
    implicit none
    class(abs_operator_t), intent(inout) :: op
    integer(ip), optional, intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    if(associated(op%op_stored)) then
       assert(.not.associated(op%op))
       call op%op_stored%fill_values(stage_)
    else if(associated(op%op)) then
       assert(.not.associated(op%op_stored))
       call op%op%fill_values(stage_)
    else
       check(1==0)
    end if

  end subroutine abs_operator_fill_values

  subroutine minus_operator_fill_values(op,stage)
    implicit none
    class(minus_operator_t), intent(inout) :: op
    integer(ip), optional  , intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    call op%op%fill_values(stage_)

  end subroutine minus_operator_fill_values

  subroutine scal_operator_fill_values(op,stage)
    implicit none
    class(scal_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    call op%op%fill_values(stage_)

  end subroutine scal_operator_fill_values

  !-----------------------------!
  ! free_values implementations !
  !-----------------------------!

  subroutine binary_operator_free_values(op,stage)
    implicit none
    class(binary_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)      :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    call op%op1%free_values(stage_)
    call op%op2%free_values(stage_)

  end subroutine binary_operator_free_values

  subroutine abs_operator_free_values(op,stage)
    implicit none
    class(abs_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)   :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    if(associated(op%op_stored)) then
       assert(.not.associated(op%op))
       call op%op_stored%free_values(stage_)
    else if(associated(op%op)) then
       assert(.not.associated(op%op_stored))
       call op%op%free_values(stage_)
    else
       check(1==0)
    end if

  end subroutine abs_operator_free_values

  subroutine minus_operator_free_values(op,stage)
    implicit none
    class(minus_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)     :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    call op%op%free_values(stage_)

  end subroutine minus_operator_free_values

  subroutine scal_operator_free_values(op,stage)
    implicit none
    class(scal_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)      :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    call op%op%free_values(stage_)

  end subroutine scal_operator_free_values

end module base_operator_names
