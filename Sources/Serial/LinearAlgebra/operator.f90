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
module operator_names
  use types_names
  use memory_guard_names
  use vector_names
  use vector_space_names
  implicit none
# include "debug.i90"

  private

  ! Abstract operator (and its deferred TBPs)
  type, abstract, extends(memory_guard_t) :: operator_t
     private
     ! Which is the relationship among operator_t and vector_space_t?
     ! Option 1: always association relationship
     ! Option 2: always composition relationship
     ! Option 3: either association or composition depending on the subclass of operator_t
     ! I finally took Option 2.
     type(vector_space_t) :: domain_vector_space
     type(vector_space_t) :: range_vector_space
   contains
     procedure (apply_interface)         , deferred :: apply
     procedure (apply_fun_interface)     , deferred :: apply_fun

     procedure :: free_vector_spaces
     procedure :: get_domain_vector_space
     procedure :: get_range_vector_space

     procedure  :: sum       => sum_operator_constructor
     procedure  :: sub       => sub_operator_constructor
     procedure  :: mult      => mult_operator_constructor
     procedure  :: minus     => minus_operator_constructor
     procedure, pass(op_left)  :: scal_left => scal_left_operator_constructor
     procedure, pass(op_right) :: scal_right => scal_right_operator_constructor
     generic    :: operator(+) => sum
     generic    :: operator(*) => mult, scal_right, scal_left, apply_fun
     generic    :: operator(-) => minus, sub
  end type operator_t

  ! Son class expression_operator_t. These operators are always temporary
  ! and therefore an assignment is needed to make copies. The gfortran
  ! compiler only supports A=B when A and B are polymorphic if the assignment 
  ! is overwritten.
  type, abstract, extends(operator_t) :: expression_operator_t 
   contains
     procedure (expression_operator_assign_interface), deferred :: assign
     generic  :: assignment(=) => assign
  end type expression_operator_t

  type, abstract, extends(expression_operator_t) :: binary_operator_t
     class(operator_t), pointer :: op1 => null(), op2 => null()
   contains
     procedure :: default_initialization => binary_operator_default_init
     procedure :: free    => binary_operator_destructor
     procedure :: assign  => binary_operator_copy
  end type binary_operator_t

  type, abstract, extends(expression_operator_t) :: unary_operator_t
     class(operator_t), pointer :: op => null()
   contains
     procedure :: default_initialization => unary_operator_default_init
     procedure :: free    => unary_operator_destructor
     procedure :: assign  => unary_operator_copy
  end type unary_operator_t


  type, extends(operator_t) :: dynamic_state_operator_t
     class(operator_t), pointer :: op_stored => null()
     class(operator_t), pointer :: op        => null()
   contains
     procedure  :: default_initialization => dynamic_state_operator_default_init
     procedure  :: apply     => dynamic_state_operator_apply
     procedure  :: apply_fun => dynamic_state_operator_apply_fun
     procedure  :: free  => dynamic_state_operator_destructor
     procedure  :: assign => dynamic_state_operator_constructor
     generic    :: assignment(=) => assign
  end type dynamic_state_operator_t

  type, extends(binary_operator_t) :: sum_operator_t
   contains
     procedure  :: apply => sum_operator_apply
     procedure  :: apply_fun => sum_operator_apply_fun 
  end type sum_operator_t

  type, extends(binary_operator_t) :: sub_operator_t
   contains
     procedure  :: apply => sub_operator_apply
     procedure  :: apply_fun => sub_operator_apply_fun 
  end type sub_operator_t

  type, extends(binary_operator_t) :: mult_operator_t
   contains
     procedure  :: apply => mult_operator_apply
     procedure  :: apply_fun => mult_operator_apply_fun 
  end type mult_operator_t

  type, extends(unary_operator_t) :: scal_operator_t
     real(rp) :: alpha
   contains
     procedure  :: apply => scal_operator_apply
     procedure  :: apply_fun => scal_operator_apply_fun
  end type scal_operator_t

  type, extends(unary_operator_t) :: minus_operator_t
   contains
     procedure  :: apply => minus_operator_apply
     procedure  :: apply_fun => minus_operator_apply_fun
  end type minus_operator_t

  abstract interface
     ! op%apply(x,y) <=> y <- op*x
     ! Implicitly assumes that y is already allocated
     subroutine apply_interface(op,x,y) 
       import :: operator_t, vector_t
       implicit none
       class(operator_t), intent(in)    :: op
       class(vector_t) , intent(in)    :: x
       class(vector_t) , intent(inout) :: y 
     end subroutine apply_interface

     ! op%apply(x)
     ! Allocates room for (temporary) y
     function apply_fun_interface(op,x) result(y)
       import :: operator_t, vector_t
       implicit none
       class(operator_t), intent(in)  :: op
       class(vector_t) , intent(in)  :: x
       class(vector_t) , allocatable :: y 
     end function apply_fun_interface

     subroutine expression_operator_assign_interface(op1,op2)
       import :: operator_t, expression_operator_t
       implicit none
       class(operator_t)      , intent(in)    :: op2
       class(expression_operator_t), intent(inout) :: op1
     end subroutine expression_operator_assign_interface
  end interface

  public :: dynamic_state_operator_t, operator_t, sum_operator_t, scal_operator_t

contains

  function get_domain_vector_space ( this )
    implicit none
    class(operator_t), target, intent(in) :: this
    type(vector_space_t), pointer :: get_domain_vector_space
    get_domain_vector_space => this%domain_vector_space
  end function get_domain_vector_space

  function get_range_vector_space ( this )
    implicit none
    class(operator_t), target, intent(in) :: this
    type(vector_space_t), pointer :: get_range_vector_space
    get_range_vector_space => this%range_vector_space
  end function get_range_vector_space

  subroutine free_vector_spaces ( this )
    implicit none
    class(operator_t), intent(inout) :: this
    call this%domain_vector_space%free()
    call this%range_vector_space%free()
  end subroutine free_vector_spaces

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

    select type(that => this%op1)
       class is(expression_operator_t)
       call that%CleanTemp()
       class is(dynamic_state_operator_t)
       call that%CleanTemp()
       class default
       check(1==0)
    end select
    deallocate(this%op1)

    select type(that => this%op2)
       class is(expression_operator_t)
       call that%CleanTemp()
       class is(dynamic_state_operator_t)
       call that%CleanTemp()
       class default
       check(1==0)
    end select
    deallocate(this%op2)
  end subroutine binary_operator_destructor

  subroutine binary_operator_copy(op1,op2)
    implicit none
    class(operator_t)  , intent(in)    :: op2
    class(binary_operator_t), intent(inout) :: op1

    select type(op2)
       class is(binary_operator_t)
       call binary_operator_constructor(op2%op1,op2%op2,op1)
       class default
       check(1==0)
    end select
  end subroutine binary_operator_copy

  subroutine binary_operator_constructor(op1,op2,res) 
    implicit none
    class(operator_t)  , intent(in)    :: op1, op2
    class(binary_operator_t), intent(inout) :: res

    call op1%GuardTemp()
    call op2%GuardTemp()

    ! Allocate op1
    select type(op1)
       class is(expression_operator_t)
       allocate(res%op1,mold=op1); call res%op1%default_initialization()
       class default
       allocate(dynamic_state_operator_t::res%op1)
    end select
    ! Assign op1
    select type(this => res%op1)
       class is(expression_operator_t)
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
       class is(dynamic_state_operator_t)
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select

    ! Allocate op2
    select type(op2)
       class is(expression_operator_t)
       allocate(res%op2,mold=op2); call res%op2%default_initialization()
       class default
       allocate(dynamic_state_operator_t::res%op2)
    end select
    ! Assign op2
    select type(that => res%op2)
       class is(expression_operator_t)
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%GuardTemp()
       class is(dynamic_state_operator_t)
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%SetTemp()
       call that%GuardTemp()
    end select
    call res%SetTemp()
    call op1%CleanTemp()
    call op2%CleanTemp()
  end subroutine binary_operator_constructor

  subroutine unary_operator_default_init(this)
    implicit none
    class(unary_operator_t), intent(inout) :: this
    nullify(this%op)
    call this%NullifyTemporary()
  end subroutine unary_operator_default_init

  subroutine unary_operator_destructor(this)
    implicit none
    class(unary_operator_t), intent(inout) :: this

    select type(that => this%op)
       class is(expression_operator_t)
       call that%CleanTemp()
       class is(dynamic_state_operator_t)
       call that%CleanTemp()
       class default
       check(1==0)
    end select
    deallocate(this%op)
  end subroutine unary_operator_destructor

  subroutine unary_operator_copy(op1,op2)
    implicit none
    class(unary_operator_t), intent(inout) :: op1
    class(operator_t)  , intent(in)    :: op2

    select type(op2)
       class is(unary_operator_t)
       call unary_operator_constructor(op2%op,op1)
       class default
       check(1==0)
    end select
  end subroutine unary_operator_copy

  subroutine unary_operator_constructor(op,res) 
    implicit none
    class(operator_t)  , intent(in)    :: op
    class(unary_operator_t), intent(inout) :: res

    call op%GuardTemp()

    ! Allocate op1
    select type(op)
       class is(expression_operator_t)
       allocate(res%op,mold=op); call res%op%default_initialization()
       class default
       allocate(dynamic_state_operator_t::res%op)
    end select

    ! Assign op1
    select type(this => res%op)
       class is(expression_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
       class is(dynamic_state_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select

    call res%SetTemp()
    call op%CleanTemp()
  end subroutine unary_operator_constructor

  subroutine dynamic_state_operator_default_init(this)
    implicit none
    class(dynamic_state_operator_t), intent(inout) :: this
    nullify(this%op)
    nullify(this%op_stored)
    call this%NullifyTemporary()
  end subroutine dynamic_state_operator_default_init

  recursive subroutine dynamic_state_operator_constructor(op1,op2)
    implicit none
    class(dynamic_state_operator_t) , intent(inout) :: op1
    class(operator_t), intent(in), target  :: op2

    call op1%free()

    call op2%GuardTemp()
    select type(op2)
       class is(dynamic_state_operator_t) ! Can be temporary (or not)
       if(associated(op2%op_stored)) then
          assert(.not.associated(op2%op))
          allocate(op1%op_stored, mold = op2%op_stored); call op1%op_stored%default_initialization()
          select type(this => op1%op_stored)
             class is(expression_operator_t)
             this = op2%op_stored
             class is(dynamic_state_operator_t)
             this = op2%op_stored
             class default
             check(1==0)
          end select
          call op1%op_stored%GuardTemp()
       else if(associated(op2%op)) then
          assert(.not.associated(op2%op_stored))
          op1%op => op2%op
          call op1%op%GuardTemp()
       else
          check(1==0)
       end if
       class is(expression_operator_t) ! Temporary
       allocate(op1%op_stored,mold=op2); call op1%op_stored%default_initialization()
       select type(this => op1%op_stored)
          class is(expression_operator_t)
          this = op2              ! Here = overloaded
       end select
       call op1%op_stored%GuardTemp()
       class default                 ! Cannot be temporary (I don't know how to copy it!)
       op1%op => op2
       call op1%op%GuardTemp()
    end select
    call op2%CleanTemp()
  end subroutine dynamic_state_operator_constructor

  subroutine dynamic_state_operator_destructor(this)
    implicit none
    class(dynamic_state_operator_t), intent(inout) :: this

    if(associated(this%op)) then
       assert(.not.associated(this%op_stored))
       ! Nothing to free, the pointer points to permanent data
       this%op => null()
    else if(associated(this%op_stored)) then
       assert(.not.associated(this%op))
       call this%op_stored%CleanTemp()
       deallocate(this%op_stored)
    end if
  end subroutine dynamic_state_operator_destructor

  function minus_operator_constructor(op) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op
    type(minus_operator_t) :: res
    call unary_operator_constructor(op,res)
  end function minus_operator_constructor

!!$  !--------------------------------------------------------------------!
!!$  ! Construction and deallocation functions/subroutines of the nodes of! 
!!$  ! the tree that represents an expression among matrix operators      !
!!$  ! -------------------------------------------------------------------!
  function sum_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sum_operator_t)  :: res
    call binary_operator_constructor(op1,op2,res) 
  end function sum_operator_constructor

  function sub_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sub_operator_t)  :: res
    call binary_operator_constructor(op1,op2,res) 
  end function sub_operator_constructor

  function mult_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(mult_operator_t) :: res
    call binary_operator_constructor(op1,op2,res)
  end function mult_operator_constructor

  function scal_left_operator_constructor(alpha, op_left) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op_left
    real(rp)             , intent(in)  :: alpha
    type(scal_operator_t)              :: res
    res%alpha=alpha
    call unary_operator_constructor(op_left,res)
  end function scal_left_operator_constructor

  function scal_right_operator_constructor(op_right, alpha) result (res)
    implicit none
    class(operator_t)       , intent(in)  :: op_right
    real(rp)                , intent(in)  :: alpha
    type(scal_operator_t)                 :: res
    res%alpha=alpha
    call unary_operator_constructor(op_right,res)
  end function scal_right_operator_constructor

  !-------------------------------------!
  ! apply_fun and apply implementations !
  !-------------------------------------!
  function sum_operator_apply_fun(op,x) result(y)
    implicit none
    class(sum_operator_t), intent(in)       :: op
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
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
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
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
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op2%apply ( op%op1*x, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine mult_operator_apply

  function  scal_operator_apply_fun(op,x) result(y)
    implicit none
    class(scal_operator_t), intent(in)      :: op
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
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
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( -1.0, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine minus_operator_apply

  function  dynamic_state_operator_apply_fun(op,x) result(y)
    implicit none
    class(dynamic_state_operator_t), intent(in)       :: op
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
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
  end function dynamic_state_operator_apply_fun

  subroutine dynamic_state_operator_apply(op,x,y)
    implicit none
    class(dynamic_state_operator_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
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
  end subroutine dynamic_state_operator_apply

end module operator_names
