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
     ! Deferred methods
     procedure (apply_interface)    , deferred :: apply
     procedure (is_linear_interface), deferred :: is_linear 
     
     procedure :: get_tangent                             => operator_get_tangent
     procedure :: get_translation                         => operator_get_translation
     procedure :: free_vector_spaces                      => operator_free_vector_spaces
     procedure :: get_domain_vector_space                 => operator_get_domain_vector_space
     procedure :: get_range_vector_space                  => operator_get_range_vector_space
     procedure :: abort_if_not_in_range                   => operator_abort_if_not_in_range
     procedure :: abort_if_not_in_domain                  => operator_abort_if_not_in_domain
     procedure :: apply_fun                               => operator_apply_fun

     procedure  :: sum       => sum_operator_constructor
     procedure  :: sub       => sub_operator_constructor
     procedure  :: mult      => mult_operator_constructor
     procedure  :: minus     => minus_operator_constructor
     procedure, pass(op_left)  :: scal_left => scal_left_operator_constructor
     procedure, pass(op_right) :: scal_right => scal_right_operator_constructor
     generic    :: operator(+) => sum
     generic    :: operator(*) => mult, scal_right, scal_left, apply_fun
     generic    :: operator(-) => sub
     generic    :: operator(.minus.) => minus
  end type operator_t

  ! Son class expression_operator_t. These operators are always temporary
  ! and therefore an assignment is needed to make copies. The gfortran
  ! compiler only supports A=B when A and B are polymorphic if the assignment 
  ! is overwritten.
  type, abstract, extends(operator_t) :: expression_operator_t 
   contains
     procedure (expression_operator_assign_interface), deferred :: assign
     generic                                                    :: assignment(=) => assign
  end type expression_operator_t

  type, abstract, extends(expression_operator_t) :: binary_operator_t
     class(operator_t), pointer :: op1 => null(), op2 => null()
   contains
     procedure :: default_initialization => binary_operator_default_init
     procedure :: free    => binary_operator_destructor
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
     procedure  :: is_linear => dynamic_state_operator_is_linear
     procedure  :: free  => dynamic_state_operator_destructor
     procedure  :: assign => dynamic_state_operator_constructor
     generic    :: assignment(=) => assign
  end type dynamic_state_operator_t

  type, extends(binary_operator_t) :: sum_operator_t
   contains
     procedure  :: apply     => sum_operator_apply
     procedure  :: is_linear => sum_operator_is_linear
     procedure  :: assign => sum_operator_copy
  end type sum_operator_t

  type, extends(binary_operator_t) :: sub_operator_t
   contains
     procedure  :: apply => sub_operator_apply
     procedure  :: is_linear => sub_operator_is_linear
     procedure  :: assign => sub_operator_copy
  end type sub_operator_t

  type, extends(binary_operator_t) :: mult_operator_t
   contains
     procedure  :: apply => mult_operator_apply
     procedure  :: is_linear => mult_operator_is_linear
     procedure  :: assign => mult_operator_copy
  end type mult_operator_t

  type, extends(unary_operator_t) :: scal_operator_t
     real(rp) :: alpha
   contains
     procedure  :: apply => scal_operator_apply
     procedure  :: is_linear => scal_operator_is_linear
     ! scal_operator must overwrite assign for unary_operator_t
     ! as it adds a new member variable "alpha" to unary_operator_t
     procedure  :: assign => scal_operator_copy
  end type scal_operator_t

  type, extends(unary_operator_t) :: minus_operator_t
   contains
     procedure  :: apply => minus_operator_apply
     procedure  :: is_linear => minus_operator_is_linear
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

     subroutine expression_operator_assign_interface(op1,op2)
       import :: operator_t, expression_operator_t
       implicit none
       class(operator_t)      , intent(in)    :: op2
       class(expression_operator_t), intent(inout) :: op1
     end subroutine expression_operator_assign_interface
     
     function is_linear_interface(op) 
       import :: operator_t
       implicit none
       class(operator_t), intent(in)    :: op
       logical                          :: is_linear_interface
     end function is_linear_interface
  end interface

  public :: dynamic_state_operator_t, operator_t !, sum_operator_t, scal_operator_t
  public :: operator_get_domain_vector_space, operator_get_range_vector_space, &
            operator_abort_if_not_in_domain, operator_abort_if_not_in_range

contains

  function operator_get_domain_vector_space ( this )
    implicit none
    class(operator_t), target, intent(in) :: this
    type(vector_space_t), pointer :: operator_get_domain_vector_space
    operator_get_domain_vector_space => this%domain_vector_space
  end function operator_get_domain_vector_space

  function operator_get_range_vector_space ( this )
    implicit none
    class(operator_t), target, intent(in) :: this
    type(vector_space_t), pointer :: operator_get_range_vector_space
    operator_get_range_vector_space => this%range_vector_space
  end function operator_get_range_vector_space

  subroutine operator_free_vector_spaces ( this )
    implicit none
    class(operator_t), intent(inout) :: this
    call this%domain_vector_space%free()
    call this%range_vector_space%free()
  end subroutine operator_free_vector_spaces
  
  function  operator_apply_fun(op,x) result(y)
    implicit none
    class(operator_t), intent(in)       :: op
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
    call op%GuardTemp()
    call x%GuardTemp()
    call op%range_vector_space%create_vector(y)
    call op%apply(x,y)
    call x%CleanTemp()
    call y%SetTemp()
    call op%CleanTemp()
  end function operator_apply_fun
  
  subroutine operator_abort_if_not_in_domain ( this, vector )
    implicit none
    class(operator_t), intent(in)  :: this
    class(vector_t)  , intent(in)  :: vector
    if ( .not. this%domain_vector_space%belongs_to(vector) ) then
       write(0,'(a)') 'operator_t%abort_if_not_in_domain: input vector does not belong to domain vector space' 
       check(.false.)
    end if
  end subroutine operator_abort_if_not_in_domain
  
  subroutine operator_abort_if_not_in_range ( this, vector )
    implicit none
    class(operator_t), intent(in) :: this
    class(vector_t)  , intent(in) :: vector
    if ( .not. this%range_vector_space%belongs_to(vector) ) then
       write(0,'(a)') 'operator_t%abort_if_not_in_range: input vector does not belong to range vector space'
       check(.false.)
    end if
  end subroutine operator_abort_if_not_in_range
  
  function operator_get_tangent(op,x) result(tangent)
    implicit none
    class(operator_t)          , intent(in) :: op
    class(vector_t)  , optional, intent(in) :: x
    type(dynamic_state_operator_t)          :: tangent 
    
    if (op%is_linear()) then
      tangent = op
      call tangent%SetTemp()
    else
      write(0,'(a)') 'Error: operator_t%get_tangent(x) :: tangent unknown, your MUST override operator_t%get_tangent(x)'
      check(.false.)
    end if
  end function operator_get_tangent
  
  function operator_get_translation(op) result(translation)
    implicit none
    class(operator_t) , intent(in) :: op
    class(vector_t)   , pointer    :: translation
    
    if (op%is_linear()) then
      ! Linear operators do not have any translation
      ! This situation is signaled by returning a nullified pointer to the caller
      ! The caller should be aware that this TBP may return a nullified pointer
      nullify(translation)
    else
      write(0,'(a)') 'Error: operator_t%get_translation() :: translation unknown, your MUST override operator_t%get_translation()'
      check(.false.)
    end if
  end function operator_get_translation
  
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
    call this%free_vector_spaces()
  end subroutine binary_operator_destructor

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
    call this%free_vector_spaces()
  end subroutine unary_operator_destructor

  subroutine unary_operator_copy(op1,op2)
    implicit none
    class(unary_operator_t), intent(inout) :: op1
    class(operator_t)  , intent(in)    :: op2

    select type(op2)
       class is(unary_operator_t)
       call op2%op%domain_vector_space%clone(op1%domain_vector_space)
       call op2%op%range_vector_space%clone(op1%range_vector_space)
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
    
    if ( associated(op1%op) ) then
      call op1%op%domain_vector_space%clone(op1%domain_vector_space)
      call op1%op%range_vector_space%clone(op1%range_vector_space)
    else if ( associated(op1%op_stored) ) then
      call op1%op_stored%domain_vector_space%clone(op1%domain_vector_space)
      call op1%op_stored%range_vector_space%clone(op1%range_vector_space)
    else
      check(.false.)
    end if
    
  end subroutine dynamic_state_operator_constructor

  subroutine dynamic_state_operator_destructor(this)
    implicit none
    class(dynamic_state_operator_t), intent(inout) :: this

    if(associated(this%op)) then
       assert(.not.associated(this%op_stored))
       ! Nothing to free, the pointer points to permanent data
       this%op => null()
       call this%free_vector_spaces()
    else if(associated(this%op_stored)) then
       assert(.not.associated(this%op))
       call this%op_stored%CleanTemp()
       deallocate(this%op_stored)
       call this%free_vector_spaces()
    end if
  end subroutine dynamic_state_operator_destructor

  function minus_operator_constructor(op) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op
    type(minus_operator_t) :: res
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op
    type(vector_space_t), pointer :: domain_res
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op
    type(vector_space_t), pointer :: range_res
    
    domain_op => op%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()
    range_op => op%get_range_vector_space()
    range_res => res%get_range_vector_space()
    call unary_operator_constructor(op,res)
    call range_op%clone(range_res)
    call domain_op%clone(domain_res)
  end function minus_operator_constructor

!!$  !--------------------------------------------------------------------!
!!$  ! Construction and deallocation functions/subroutines of the nodes of! 
!!$  ! the tree that represents an expression among matrix operators      !
!!$  ! -------------------------------------------------------------------!
  function sum_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sum_operator_t)  :: res
    
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op1
    type(vector_space_t), pointer :: domain_op2
    type(vector_space_t), pointer :: domain_res
    
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op1
    type(vector_space_t), pointer :: range_op2
    type(vector_space_t), pointer :: range_res
        
    domain_op1 => op1%get_domain_vector_space()
    domain_op2 => op2%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()

    range_op1 => op1%get_range_vector_space()
    range_op2 => op2%get_range_vector_space()
    range_res => res%get_range_vector_space()
    
    if ( .not. domain_op1%equal_to(domain_op2) ) then
       write(0,'(a)') 'sum_operator_t%constructor: domain(op1)/=domain(op2)'
       check(.false.)
    end if
        
    if ( .not. range_op1%equal_to(range_op2) ) then
       write(0,'(a)') 'sum_operator_t%constructor: range(op1)/=range(op2)'
       check(.false.)
    end if
     
    call binary_operator_constructor(op1,op2,res) 
    
    call range_op1%clone(range_res)
    call domain_op1%clone(domain_res)
  end function sum_operator_constructor
  
  subroutine sum_operator_copy(op1,op2)
    implicit none
    class(sum_operator_t), intent(inout) :: op1
    class(operator_t)  , intent(in)    :: op2 
    select type(op2)
       class is(sum_operator_t)
       call op2%op1%domain_vector_space%clone(op1%domain_vector_space)
       call op2%op1%range_vector_space%clone(op1%range_vector_space)
       call binary_operator_constructor(op2%op1,op2%op2,op1)
       class default
       check(1==0)
    end select
  end subroutine sum_operator_copy

  function sub_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sub_operator_t)  :: res
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op1
    type(vector_space_t), pointer :: domain_op2
    type(vector_space_t), pointer :: domain_res
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op1
    type(vector_space_t), pointer :: range_op2
    type(vector_space_t), pointer :: range_res
    
    domain_op1 => op1%get_domain_vector_space()
    domain_op2 => op2%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()
    range_op1 => op1%get_range_vector_space()
    range_op2 => op2%get_range_vector_space()
    range_res => res%get_range_vector_space()
    if ( .not. domain_op1%equal_to(domain_op2) ) then
       write(0,'(a)') 'sub_operator_t%constructor: domain(op1)/=domain(op2)'
       check(.false.)
    end if
        
    if ( .not. range_op1%equal_to(range_op2) ) then
       write(0,'(a)') 'sub_operator_t%constructor: range(op1)/=range(op2)'
       check(.false.)
    end if
    call binary_operator_constructor(op1,op2,res)
    call range_op1%clone(range_res)
    call domain_op1%clone(domain_res)
  end function sub_operator_constructor

  subroutine sub_operator_copy(op1,op2)
    implicit none
    class(sub_operator_t), intent(inout) :: op1
    class(operator_t)  , intent(in)    :: op2 
    select type(op2)
       class is(sub_operator_t)
       call op2%op1%domain_vector_space%clone(op1%domain_vector_space)
       call op2%op1%range_vector_space%clone(op1%range_vector_space)
       call binary_operator_constructor(op2%op1,op2%op2,op1)
       class default
       check(1==0)
    end select
  end subroutine sub_operator_copy
  
  function mult_operator_constructor(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(mult_operator_t) :: res
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op1
    type(vector_space_t), pointer :: domain_op2
    type(vector_space_t), pointer :: domain_res
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op1
    type(vector_space_t), pointer :: range_op2
    type(vector_space_t), pointer :: range_res
    
    domain_op1 => op1%get_domain_vector_space()
    domain_op2 => op2%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()
    range_op1 => op1%get_range_vector_space()
    range_op2 => op2%get_range_vector_space()
    range_res => res%get_range_vector_space()
    if ( .not. domain_op1%equal_to(range_op2) ) then
       write(0,'(a)') 'mult_operator_t%constructor: domain(op1)/=range(op2)'
       check(.false.)
    end if
    call binary_operator_constructor(op1,op2,res)
    call range_op1%clone(range_res)
    call domain_op2%clone(domain_res)
  end function mult_operator_constructor
  
   subroutine mult_operator_copy(op1,op2)
    implicit none
    class(mult_operator_t), intent(inout) :: op1
    class(operator_t)     , intent(in)    :: op2 
    select type(op2)
       class is(mult_operator_t)
       call op2%op1%range_vector_space%clone(op1%range_vector_space)
       call op2%op2%domain_vector_space%clone(op1%domain_vector_space)
       call binary_operator_constructor(op2%op1,op2%op2,op1)
       class default
       check(1==0)
    end select
  end subroutine mult_operator_copy

  function scal_left_operator_constructor(alpha, op_left) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op_left
    real(rp)             , intent(in)  :: alpha
    type(scal_operator_t)              :: res
    
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op_left
    type(vector_space_t), pointer :: domain_res
    
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op_left
    type(vector_space_t), pointer :: range_res
    
    domain_op_left => op_left%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()
    range_op_left => op_left%get_range_vector_space()
    range_res => res%get_range_vector_space()
    res%alpha=alpha
    call unary_operator_constructor(op_left,res)
    call range_op_left%clone(range_res)
    call domain_op_left%clone(domain_res)
  end function scal_left_operator_constructor

  function scal_right_operator_constructor(op_right, alpha) result (res)
    implicit none
    class(operator_t)       , intent(in)  :: op_right
    real(rp)                , intent(in)  :: alpha
    type(scal_operator_t)                 :: res
    ! Pointers to domain vector spaces
    type(vector_space_t), pointer :: domain_op_right
    type(vector_space_t), pointer :: domain_res
    ! Pointers to range vector spaces
    type(vector_space_t), pointer :: range_op_right
    type(vector_space_t), pointer :: range_res
    
    domain_op_right => op_right%get_domain_vector_space()
    domain_res => res%get_domain_vector_space()
    range_op_right => op_right%get_range_vector_space()
    range_res => res%get_range_vector_space()    
    res%alpha=alpha
    call unary_operator_constructor(op_right,res)
    call range_op_right%clone(range_res)
    call domain_op_right%clone(domain_res)
  end function scal_right_operator_constructor
  
    subroutine scal_operator_copy(op1,op2)
    implicit none
    class(scal_operator_t), intent(inout) :: op1
    class(operator_t)       , intent(in)    :: op2

    select type(op2)
       class is(scal_operator_t)
       op1%alpha = op2%alpha
       call op2%op%domain_vector_space%clone(op1%domain_vector_space)
       call op2%op%range_vector_space%clone(op1%range_vector_space)
       call unary_operator_constructor(op2%op,op1)
       class default
       check(1==0)
    end select
  end subroutine scal_operator_copy


  !-------------------------------------!
  ! apply implementations               !
  !-------------------------------------!
  subroutine sum_operator_apply(op,x,y)
    implicit none
    class(sum_operator_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    class(vector_t), allocatable :: w

    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call op%GuardTemp()
    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    call op%op1%range_vector_space%create_vector(w)
    call op%op1%apply(x,w)
    ! y <- 1.0 * op1*x + 1.0*y
    call y%axpby( 1.0, w, 1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
    call w%free()
    deallocate(w)
  end subroutine sum_operator_apply

  subroutine sub_operator_apply(op,x,y)
    implicit none
    class(sub_operator_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    class(vector_t), allocatable :: w
    
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call op%GuardTemp()
    call x%GuardTemp()
    ! y <- op2*x
    call op%op2%apply(x,y)
    call op%op1%range_vector_space%create_vector(w)
    call op%op1%apply(x,w)
    ! y <- 1.0 * op1*x - 1.0*y
    call y%axpby( 1.0, w, -1.0 )
    call x%CleanTemp()
    call op%CleanTemp()
    call w%free()
    deallocate(w)
  end subroutine sub_operator_apply
  
  subroutine mult_operator_apply(op,x,y)
    implicit none
    class(mult_operator_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    class(vector_t), allocatable :: w
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op1%range_vector_space%create_vector(w)
    call op%op1%apply(x,w)
    call op%op2%apply(w, y )
    call x%CleanTemp()
    call op%CleanTemp()
    call w%free()
    deallocate(w)
  end subroutine mult_operator_apply

  subroutine scal_operator_apply(op,x,y)
    implicit none
    class(scal_operator_t), intent(in)   :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    call op%GuardTemp()
    call x%GuardTemp()
    call op%op%apply(x,y)
    call y%scal( op%alpha, y )
    call x%CleanTemp()
    call op%CleanTemp()
  end subroutine scal_operator_apply

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


  subroutine dynamic_state_operator_apply(op,x,y)
    implicit none
    class(dynamic_state_operator_t), intent(in)    :: op
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
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
  
  !-------------------------------------!
  ! is_linear implementations           !
  !-------------------------------------!
  function sum_operator_is_linear(op)
    implicit none
    class(sum_operator_t), intent(in)    :: op
    logical :: sum_operator_is_linear
    sum_operator_is_linear = .false.
  end function sum_operator_is_linear

  function sub_operator_is_linear(op)
    implicit none
    class(sub_operator_t), intent(in)    :: op
    logical :: sub_operator_is_linear
    sub_operator_is_linear = .false.
  end function sub_operator_is_linear
  
  function mult_operator_is_linear(op)
    implicit none
    class(mult_operator_t), intent(in)    :: op
    logical :: mult_operator_is_linear
    mult_operator_is_linear = .false.
  end function mult_operator_is_linear

  function scal_operator_is_linear(op)
    implicit none
    class(scal_operator_t), intent(in)   :: op
    logical :: scal_operator_is_linear
    scal_operator_is_linear = .false.
  end function scal_operator_is_linear

  function minus_operator_is_linear(op)
    implicit none
    class(minus_operator_t), intent(in)  :: op
    logical :: minus_operator_is_linear
    minus_operator_is_linear = .false.
  end function minus_operator_is_linear

  function dynamic_state_operator_is_linear(op)
    implicit none
    class(dynamic_state_operator_t), intent(in) :: op
    logical :: dynamic_state_operator_is_linear
    dynamic_state_operator_is_linear = .false.
  end function dynamic_state_operator_is_linear
  
end module operator_names
