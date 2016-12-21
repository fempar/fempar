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
     procedure (apply_add_interface), deferred :: apply_add
     procedure (is_linear_interface), deferred :: is_linear 
     
     procedure :: get_tangent                             => operator_get_tangent
     procedure :: get_translation                         => operator_get_translation
     procedure :: free_vector_spaces                      => operator_free_vector_spaces
     procedure :: create_domain_vector                    => operator_create_domain_vector
     procedure :: create_range_vector                     => operator_create_range_vector
     procedure :: vector_spaces_are_created               => operator_vector_spaces_are_created
     procedure :: get_domain_vector_space                 => operator_get_domain_vector_space
     procedure :: get_range_vector_space                  => operator_get_range_vector_space
     procedure :: abort_if_not_in_range                   => operator_abort_if_not_in_range
     procedure :: abort_if_not_in_domain                  => operator_abort_if_not_in_domain
     procedure :: apply_fun                               => operator_apply_fun

     procedure  :: sum       => sum_operator_create
     procedure  :: sub       => sub_operator_create
     procedure  :: mult      => mult_operator_create
     procedure  :: minus     => minus_operator_create
     procedure, pass(op_left)  :: scal_left => scal_left_operator_create
     procedure, pass(op_right) :: scal_right => scal_right_operator_create
     procedure  :: identity => identity_operator_create
     generic    :: operator(+)          => sum
     generic    :: operator(*)          => mult, scal_right, scal_left, apply_fun
     generic    :: operator(-)          => sub
     generic    :: operator(.minus.)    => minus
     generic    :: operator(.identity.) => identity
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
     procedure :: free    => binary_operator_free
     procedure :: assign  => binary_operator_assign
  end type binary_operator_t

  type, abstract, extends(expression_operator_t) :: unary_operator_t
     class(operator_t), pointer :: op => null()
   contains
     procedure :: default_initialization => unary_operator_default_init
     procedure :: free    => unary_operator_free
     procedure :: assign  => unary_operator_assign
  end type unary_operator_t

  type, extends(operator_t) :: lvalue_operator_t
     class(operator_t), pointer :: op_stored     => null()
     class(operator_t), pointer :: op_referenced => null()
   contains
     procedure  :: default_initialization => lvalue_operator_default_init
     procedure  :: apply     => lvalue_operator_apply
     procedure  :: apply_add => lvalue_operator_apply_add
     procedure  :: is_linear => lvalue_operator_is_linear
     procedure  :: free      => lvalue_operator_free
     procedure  :: assign    => lvalue_operator_create
     generic    :: assignment(=) => assign
  end type lvalue_operator_t
  
  type, extends(binary_operator_t) :: sum_operator_t
   contains
     procedure  :: apply     => sum_operator_apply
     procedure  :: apply_add => sum_operator_apply_add
     procedure  :: is_linear => sum_operator_is_linear
  end type sum_operator_t

  type, extends(binary_operator_t) :: sub_operator_t
   contains
     procedure  :: apply     => sub_operator_apply
     procedure  :: apply_add => sub_operator_apply_add
     procedure  :: is_linear => sub_operator_is_linear
  end type sub_operator_t

  type, extends(binary_operator_t) :: mult_operator_t
   contains
     procedure  :: apply     => mult_operator_apply
     procedure  :: apply_add => mult_operator_apply_add
     procedure  :: is_linear => mult_operator_is_linear
  end type mult_operator_t

  type, extends(unary_operator_t) :: scal_operator_t
     real(rp) :: alpha
   contains
     procedure  :: apply     => scal_operator_apply
     procedure  :: apply_Add => scal_operator_apply_add
     procedure  :: is_linear => scal_operator_is_linear
     ! scal_operator must overwrite assign for unary_operator_t
     ! as it adds a new member variable "alpha" to unary_operator_t
     procedure  :: assign => scal_operator_assign
  end type scal_operator_t
  
  type, extends(unary_operator_t) :: minus_operator_t
   contains
     procedure  :: apply     => minus_operator_apply
     procedure  :: apply_add => minus_operator_apply_add
     procedure  :: is_linear => minus_operator_is_linear
  end type minus_operator_t
    
  type, extends(unary_operator_t) :: identity_operator_t
  contains
     procedure :: apply     => identity_operator_apply
     procedure :: apply_add => identity_operator_apply_add
     procedure :: is_linear => identity_operator_is_linear
  end type

  abstract interface
     ! op%apply(x,y) <=> y <- op*x
     ! Implicitly assumes that y is already allocated
     subroutine apply_interface(this,x,y) 
       import :: operator_t, vector_t
       implicit none
       class(operator_t), intent(in)    :: this
       class(vector_t) , intent(in)    :: x
       class(vector_t) , intent(inout) :: y 
     end subroutine apply_interface 
     
     ! op%apply_add(x,y) <=> y <- op*x+y
     ! Implicitly assumes that y is already allocated
     subroutine apply_add_interface(this,x,y) 
       import :: operator_t, vector_t
       implicit none
       class(operator_t), intent(in)    :: this
       class(vector_t) , intent(in)    :: x
       class(vector_t) , intent(inout) :: y 
     end subroutine apply_add_interface 

     subroutine expression_operator_assign_interface(this,rvalue)
       import :: operator_t, expression_operator_t
       implicit none
       class(operator_t)      , intent(in)    :: rvalue
       class(expression_operator_t), intent(inout) :: this
     end subroutine expression_operator_assign_interface
     
     function is_linear_interface(this) 
       import :: operator_t
       implicit none
       class(operator_t), intent(in)    :: this
       logical                          :: is_linear_interface
     end function is_linear_interface
  end interface

  public :: lvalue_operator_t, operator_t !, sum_operator_t, scal_operator_t
  public :: operator_get_domain_vector_space, operator_get_range_vector_space, &
            operator_abort_if_not_in_domain, operator_abort_if_not_in_range

contains

  subroutine operator_create_domain_vector ( this, vector )
    implicit none
    class(operator_t), target, intent(in) :: this
    class(vector_t), allocatable, intent(inout) :: vector
    assert(this%domain_vector_space%is_created())
    call this%domain_vector_space%create_vector(vector)
  end subroutine operator_create_domain_vector

  subroutine operator_create_range_vector ( this, vector )
    implicit none
    class(operator_t), target, intent(in) :: this
    class(vector_t), allocatable, intent(inout) :: vector
    assert(this%range_vector_space%is_created())
    call this%range_vector_space%create_vector(vector)
  end subroutine operator_create_range_vector

  function operator_vector_spaces_are_created( this ) result(vector_spaces_are_created)
    implicit none
    class(operator_t), target, intent(in) :: this
    logical                               :: vector_spaces_are_created
    vector_spaces_are_created = &
        (this%domain_vector_space%is_created() .and. this%range_vector_space%is_created())
  end function operator_vector_spaces_are_created

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
  
  function  operator_apply_fun(this,x) result(y)
    implicit none
    class(operator_t), intent(in)       :: this
    class(vector_t)     , intent(in)  :: x
    class(vector_t)     , allocatable :: y 
    call this%GuardTemp()
    call x%GuardTemp()
    call this%range_vector_space%create_vector(y)
    call this%apply(x,y)
    call x%CleanTemp()
    call y%SetTemp()
    call this%CleanTemp()
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
  
  function operator_get_tangent(this,x) result(tangent)
    implicit none
    class(operator_t)          , intent(in) :: this
    class(vector_t)  , optional, intent(in) :: x
    type(lvalue_operator_t)          :: tangent 
    
    if (this%is_linear()) then
      tangent = this
      call tangent%SetTemp()
    else
      write(0,'(a)') 'Error: operator_t%get_tangent(x) :: tangent unknown, your MUST override operator_t%get_tangent(x)'
      check(.false.)
    end if
  end function operator_get_tangent
  
  function operator_get_translation(this) result(translation)
    implicit none
    class(operator_t) , intent(in) :: this
    class(vector_t)   , pointer    :: translation
    
    if (this%is_linear()) then
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

  subroutine binary_operator_free(this)
    implicit none
    class(binary_operator_t), intent(inout) :: this
    integer(ip)                             :: istat
    if(associated(this%op1)) then
      call this%op1%CleanTemp()
      deallocate(this%op1, stat=istat); check(istat==0)
    endif
    if(associated(this%op2)) then
      call this%op2%CleanTemp()
      deallocate(this%op2, stat=istat); check(istat==0)
    endif
    call this%free_vector_spaces()
  end subroutine binary_operator_free

  recursive subroutine binary_operator_create(op1,op2,res) 
    implicit none
    class(operator_t)  , intent(in)    :: op1, op2
    class(binary_operator_t), intent(inout) :: res
    integer(ip)                             :: istat

    call op1%GuardTemp()
    call op2%GuardTemp()

    ! Allocate op1
    select type(op1)
       class is(expression_operator_t)
       allocate(res%op1,mold=op1, stat=istat); check(istat==0)
       call res%op1%default_initialization()
       class default
       allocate(lvalue_operator_t::res%op1, stat=istat); check(istat==0)
    end select
    ! Assign op1
    select type(this => res%op1)
       class is(expression_operator_t)
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
       class is(lvalue_operator_t)
       this = op1 ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select

    ! Allocate op2
    select type(op2)
       class is(expression_operator_t)
       allocate(res%op2,mold=op2, stat=istat); check(istat==0)
       call res%op2%default_initialization()
       class default
       allocate(lvalue_operator_t::res%op2, stat=istat); check(istat==0)
    end select
    ! Assign op2
    select type(that => res%op2)
       class is(expression_operator_t)
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%GuardTemp()
       class is(lvalue_operator_t)
       that = op2 ! Here = is overloaded (and potentially recursive)
       call that%SetTemp()
       call that%GuardTemp()
    end select
    call res%SetTemp()
    call op1%CleanTemp()
    call op2%CleanTemp()
  end subroutine binary_operator_create
  
    recursive subroutine binary_operator_assign(this,rvalue)
    implicit none
    class(binary_operator_t), intent(inout) :: this
    class(operator_t)  , intent(in)    :: rvalue 
    assert(same_type_as(this, rvalue))
    select type(rvalue)
       class is(binary_operator_t)
       call rvalue%op1%domain_vector_space%clone(this%domain_vector_space)
       call rvalue%op1%range_vector_space%clone(this%range_vector_space)
       call binary_operator_create(rvalue%op1,rvalue%op2,this)
       class default
       check(1==0)
    end select
  end subroutine binary_operator_assign

  subroutine unary_operator_default_init(this)
    implicit none
    class(unary_operator_t), intent(inout) :: this
    nullify(this%op)
    call this%NullifyTemporary()
  end subroutine unary_operator_default_init

  subroutine unary_operator_free(this)
    implicit none
    class(unary_operator_t), intent(inout) :: this
    integer(ip)                            :: istat

    if(associated(this%op)) then
       call this%op%CleanTemp()
       deallocate(this%op, stat=istat); check(istat==0)
    end if
    call this%free_vector_spaces()
  end subroutine unary_operator_free

  recursive subroutine unary_operator_assign(this,rvalue)
    implicit none
    class(unary_operator_t), intent(inout) :: this
    class(operator_t)  , intent(in)    :: rvalue
    assert(same_type_as(this, rvalue))
    select type(rvalue)
       class is(unary_operator_t)
       call rvalue%op%domain_vector_space%clone(this%domain_vector_space)
       call rvalue%op%range_vector_space%clone(this%range_vector_space)
       call unary_operator_create(rvalue%op,this)
       class default
       check(1==0)
    end select
  end subroutine unary_operator_assign

  recursive subroutine unary_operator_create(op,res) 
    implicit none
    class(operator_t)  , intent(in)    :: op
    class(unary_operator_t), intent(inout) :: res
    integer(ip)                            :: istat

    call op%GuardTemp()

    ! Allocate op1
    select type(op)
       class is(expression_operator_t)
       allocate(res%op,mold=op, stat=istat); check(istat==0)
       call res%op%default_initialization()
       class default
       allocate(lvalue_operator_t::res%op, stat=istat); check(istat==0)
    end select

    ! Assign op1
    select type(this => res%op)
       class is(expression_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%GuardTemp()
       class is(lvalue_operator_t)
       this = op ! Here = is overloaded (and potentially recursive)
       call this%SetTemp()
       call this%GuardTemp()
    end select

    call res%SetTemp()
    call op%CleanTemp()
  end subroutine unary_operator_create

  subroutine lvalue_operator_default_init(this)
    implicit none
    class(lvalue_operator_t), intent(inout) :: this
    nullify(this%op_referenced)
    nullify(this%op_stored)
    call this%NullifyTemporary()
  end subroutine lvalue_operator_default_init

  recursive subroutine lvalue_operator_create(this,op)
    implicit none
    class(lvalue_operator_t) , intent(inout) :: this
    class(operator_t), intent(in), target    :: op
    integer(ip)                              :: istat
    
    call this%free()
    call op%GuardTemp()
    select type(op)
       class is(lvalue_operator_t) ! Can be temporary (or not)
          if(associated(op%op_stored)) then
             assert(.not. associated(op%op_referenced))
             allocate(this%op_stored, mold = op%op_stored, stat=istat); check(istat==0)
             call this%op_stored%default_initialization()
             select type(this => this%op_stored)
                class is(expression_operator_t)
                   this = op%op_stored
                class is(lvalue_operator_t)
                   this = op%op_stored
                class default
                   check(1==0)
             end select
             call this%op_stored%GuardTemp()
          else if(associated(op%op_referenced)) then
             assert(.not. associated(op%op_stored))
             this%op_referenced => op%op_referenced
          else
             check(1==0)
          end if
       class is(expression_operator_t) ! Temporary
          allocate(this%op_stored,mold=op, stat=istat); check(istat==0)
          call this%op_stored%default_initialization()
          select type(this => this%op_stored)
             class is(expression_operator_t)
                this = op              ! Here = overloaded
          end select
          call this%op_stored%GuardTemp()
       class default                 ! Cannot be temporary (I don't know how to copy it!)
          this%op_referenced => op
    end select
    call op%CleanTemp()
    
    if ( associated(this%op_referenced) ) then
      call this%op_referenced%domain_vector_space%clone(this%domain_vector_space)
      call this%op_referenced%range_vector_space%clone(this%range_vector_space)
    else if ( associated(this%op_stored) ) then
      call this%op_stored%domain_vector_space%clone(this%domain_vector_space)
      call this%op_stored%range_vector_space%clone(this%range_vector_space)
    else
      check(.false.)
    end if
  end subroutine lvalue_operator_create

  subroutine lvalue_operator_free(this)
    implicit none
    class(lvalue_operator_t), intent(inout) :: this
    integer(ip)                             :: istat

    if(associated(this%op_referenced)) then
       assert(.not.associated(this%op_stored))
       ! Nothing to free, the pointer points to permanent data
       this%op_referenced => null()
       call this%free_vector_spaces()
    else if(associated(this%op_stored)) then
       assert(.not.associated(this%op_referenced))
       call this%op_stored%CleanTemp()
       deallocate(this%op_stored, stat=istat); check(istat==0)
       call this%free_vector_spaces()
    end if
  end subroutine lvalue_operator_free

  recursive function minus_operator_create(op) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op
    type(minus_operator_t) :: res
    
    call op%GuardTemp()
    
    call unary_operator_create(op,res)
    call op%range_vector_space%clone(res%range_vector_space)
    call op%domain_vector_space%clone(res%domain_vector_space)
    
    call op%CleanTemp()

  end function minus_operator_create
  
  recursive function identity_operator_create(op) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op
    type(identity_operator_t) :: res
    
    call op%GuardTemp()
    
    call unary_operator_create(op,res)
    call op%range_vector_space%clone(res%range_vector_space)
    call op%domain_vector_space%clone(res%domain_vector_space)
    
    call op%CleanTemp()
    
  end function identity_operator_create

!!$  !--------------------------------------------------------------------!
!!$  ! Construction and deallocation functions/subroutines of the nodes of! 
!!$  ! the tree that represents an expression among matrix operators      !
!!$  ! -------------------------------------------------------------------!
  recursive function sum_operator_create(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sum_operator_t)  :: res
    
    call op1%GuardTemp()
    call op2%GuardTemp()
    
    if ( .not. op1%domain_vector_space%equal_to(op2%domain_vector_space) ) then
       write(0,'(a)') 'sum_operator_t%create: domain(op1)/=domain(op2)'
       check(.false.)
    end if
        
    if ( .not. op1%range_vector_space%equal_to(op2%range_vector_space) ) then
       write(0,'(a)') 'sum_operator_t%create: range(op1)/=range(op2)'
       check(.false.)
    end if
     
    call binary_operator_create(op1,op2,res) 
    
    call op1%range_vector_space%clone(res%range_vector_space)
    call op1%domain_vector_space%clone(res%domain_vector_space)
  end function sum_operator_create

  recursive function sub_operator_create(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(sub_operator_t)  :: res
    
    call op1%GuardTemp()
    call op2%GuardTemp()

    if ( .not. op1%domain_vector_space%equal_to(op2%domain_vector_space) ) then
       write(0,'(a)') 'sub_operator_t%create: domain(op1)/=domain(op2)'
       check(.false.)
    end if
        
    if ( .not. op1%range_vector_space%equal_to(op2%range_vector_space) ) then
       write(0,'(a)') 'sub_operator_t%create: range(op1)/=range(op2)'
       check(.false.)
    end if
    call binary_operator_create(op1,op2,res)
    call op1%range_vector_space%clone(res%range_vector_space)
    call op1%domain_vector_space%clone(res%domain_vector_space)
    
    call op1%CleanTemp()
    call op2%CleanTemp()
    
  end function sub_operator_create
  
  recursive function mult_operator_create(op1,op2) result (res)
    implicit none
    class(operator_t), intent(in)  :: op1, op2
    type(mult_operator_t) :: res
    
    call op1%GuardTemp()
    call op2%GuardTemp()
    
    if ( .not. op1%domain_vector_space%equal_to(op2%range_vector_space) ) then
       write(0,'(a)') 'mult_operator_t%create: domain(op1)/=range(op2)'
       check(.false.)
    end if
    call binary_operator_create(op1,op2,res)
    call op1%range_vector_space%clone(res%range_vector_space)
    call op2%domain_vector_space%clone(res%domain_vector_space)
    
    call op1%CleanTemp()
    call op2%CleanTemp()
    
  end function mult_operator_create

  recursive function scal_left_operator_create(alpha, op_left) result (res)
    implicit none
    class(operator_t)    , intent(in)  :: op_left
    real(rp)             , intent(in)  :: alpha
    type(scal_operator_t)              :: res

    call op_left%GuardTemp()

    res%alpha=alpha
    call unary_operator_create(op_left,res)
    call op_left%range_vector_space%clone(res%range_vector_space)
    call op_left%domain_vector_space%clone(res%domain_vector_space)
    
    call op_left%CleanTemp()

  end function scal_left_operator_create

  recursive function scal_right_operator_create(op_right, alpha) result (res)
    implicit none
    class(operator_t)       , intent(in)  :: op_right
    real(rp)                , intent(in)  :: alpha
    type(scal_operator_t)                 :: res

    call op_right%GuardTemp()
  
    res%alpha=alpha
    call unary_operator_create(op_right,res)
    call op_right%range_vector_space%clone(res%range_vector_space)
    call op_right%domain_vector_space%clone(res%domain_vector_space)
    
    call op_right%CleanTemp()

  end function scal_right_operator_create
  
  recursive subroutine scal_operator_assign(this,rvalue)
    implicit none
    class(scal_operator_t), intent(inout) :: this
    class(operator_t)       , intent(in)    :: rvalue
    assert(same_type_as(this, rvalue))
    select type(rvalue)
       class is(scal_operator_t)
       this%alpha = rvalue%alpha
       call rvalue%op%domain_vector_space%clone(this%domain_vector_space)
       call rvalue%op%range_vector_space%clone(this%range_vector_space)
       call unary_operator_create(rvalue%op,this)
       class default
       check(1==0)
    end select
  end subroutine scal_operator_assign


  !-------------------------------------!
  ! apply implementations               !
  !-------------------------------------!
  recursive subroutine identity_operator_apply(this,x,y)
    implicit none
    class(identity_operator_t), intent(in) :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y

    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call y%copy(x)
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine identity_operator_apply
  
  recursive subroutine sum_operator_apply(this,x,y)
    implicit none
    class(sum_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    integer(ip)                    :: istat

    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op2%apply(x,y)
    call this%op1%apply_add(x,y)
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine sum_operator_apply

  recursive subroutine sub_operator_apply(this,x,y)
    implicit none
    class(sub_operator_t), intent(in)    :: this
    class(vector_t),       intent(in)    :: x
    class(vector_t),       intent(inout) :: y 
    type(lvalue_operator_t)              :: op
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op1%apply(x,y)
    op = .minus. this%op2
    call op%apply_add(x,y)
    call op%free()
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine sub_operator_apply
  
  recursive subroutine mult_operator_apply(this,x,y)
    implicit none
    class(mult_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    class(vector_t), allocatable   :: w
    integer(ip)                    :: istat
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op1%range_vector_space%create_vector(w)
    call this%op1%apply(x,w)
    call this%op2%apply(w,y)
    call x%CleanTemp()
    call this%CleanTemp()
    call w%free()
    deallocate(w, stat=istat); check(istat==0)
  end subroutine mult_operator_apply

  recursive subroutine scal_operator_apply(this,x,y)
    implicit none
    class(scal_operator_t), intent(in)   :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op%apply(x,y)
    call y%scal( this%alpha, y )
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine scal_operator_apply

  recursive subroutine minus_operator_apply(this,x,y)
    implicit none
    class(minus_operator_t), intent(in)  :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op%apply(x,y)
    call y%scal( -1.0, y )
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine minus_operator_apply

  recursive subroutine lvalue_operator_apply(this,x,y)
    implicit none
    class(lvalue_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    
    call this%GuardTemp()
    call x%GuardTemp()
    if(associated(this%op_stored)) then
       assert(.not. associated(this%op_referenced))
       call this%op_stored%apply(x,y)
    else if(associated(this%op_referenced)) then
       assert(.not. associated(this%op_stored))
       call this%op_referenced%apply(x,y)
    else
       check(1==0)
    end if
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine lvalue_operator_apply
  
  !-------------------------------------!
  ! apply_add implementations           !
  !-------------------------------------!
  recursive subroutine identity_operator_apply_add(this,x,y)
    implicit none
    class(identity_operator_t), intent(in) :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y

    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call y%axpby( 1.0, x, 1.0 )
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine identity_operator_apply_add
  
  recursive subroutine sum_operator_apply_add(this,x,y)
    implicit none
    class(sum_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op2%apply_add(x,y)
    call this%op1%apply_add(x,y)
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine sum_operator_apply_add

  recursive subroutine sub_operator_apply_add(this,x,y)
    implicit none
    class(sub_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y 
    type(lvalue_operator_t)        :: op
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op1%apply_add(x,y)
    op = .minus. this%op2
    call op%apply_add(x,y)
    call op%free()
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine sub_operator_apply_add
  
  recursive subroutine mult_operator_apply_add(this,x,y)
    implicit none
    class(mult_operator_t), intent(in)    :: this
    class(vector_t), intent(in)    :: x
    class(vector_t), intent(inout) :: y
    class(vector_t), allocatable   :: w
    integer(ip)                    :: istat
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    call this%op1%range_vector_space%create_vector(w)
    call this%op1%apply(x,w)
    call this%op2%apply_add(w,y)
    call x%CleanTemp()
    call this%CleanTemp()
    call w%free()
    deallocate(w, stat=istat); check(istat==0)
  end subroutine mult_operator_apply_add

  recursive subroutine scal_operator_apply_add(this,x,y)
    implicit none
    class(scal_operator_t), intent(in)    :: this
    class(vector_t),        intent(in)    :: x
    class(vector_t),        intent(inout) :: y 
    
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    call this%GuardTemp()
    call x%GuardTemp()
    if(this%alpha /= 0.0_rp) then
       call y%scal( 1.0/this%alpha, y )
       call this%op%apply_add(x,y)
       call y%scal( this%alpha, y )
    endif
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine scal_operator_apply_add

  recursive subroutine minus_operator_apply_add(this,x,y)
    implicit none
    class(minus_operator_t), intent(in)    :: this
    class(vector_t),         intent(in)    :: x
    class(vector_t),         intent(inout) :: y 
    call this%GuardTemp()
    call x%GuardTemp()
    call y%scal(-1.0, y)
    call this%op%apply_add(x,y)
    call y%scal(-1.0, y)
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine minus_operator_apply_add

  recursive subroutine lvalue_operator_apply_add(this,x,y)
    implicit none
    class(lvalue_operator_t), intent(in)    :: this
    class(vector_t),          intent(in)    :: x
    class(vector_t),          intent(inout) :: y 
    
    call this%abort_if_not_in_domain(x)
    call this%abort_if_not_in_range(y)
    
    call this%GuardTemp()
    call x%GuardTemp()
    if(associated(this%op_stored)) then
       assert(.not. associated(this%op_referenced))
       call this%op_stored%apply_add(x,y)
    else if(associated(this%op_referenced)) then
       assert(.not. associated(this%op_stored))
       call this%op_referenced%apply_add(x,y)
    else
       check(1==0)
    end if
    call x%CleanTemp()
    call this%CleanTemp()
  end subroutine lvalue_operator_apply_add  
  
  !-------------------------------------!
  ! is_linear implementations           !
  !-------------------------------------!
  function identity_operator_is_linear(this)
    implicit none
    class(identity_operator_t), intent(in)    :: this
    logical :: identity_operator_is_linear
    identity_operator_is_linear = .false.
  end function identity_operator_is_linear
  
  function sum_operator_is_linear(this)
    implicit none
    class(sum_operator_t), intent(in)    :: this
    logical :: sum_operator_is_linear
    sum_operator_is_linear = .false.
  end function sum_operator_is_linear

  function sub_operator_is_linear(this)
    implicit none
    class(sub_operator_t), intent(in)    :: this
    logical :: sub_operator_is_linear
    sub_operator_is_linear = .false.
  end function sub_operator_is_linear
  
  function mult_operator_is_linear(this)
    implicit none
    class(mult_operator_t), intent(in)    :: this
    logical :: mult_operator_is_linear
    mult_operator_is_linear = .false.
  end function mult_operator_is_linear

  function scal_operator_is_linear(this)
    implicit none
    class(scal_operator_t), intent(in)   :: this
    logical :: scal_operator_is_linear
    scal_operator_is_linear = .false.
  end function scal_operator_is_linear

  function minus_operator_is_linear(this)
    implicit none
    class(minus_operator_t), intent(in)  :: this
    logical :: minus_operator_is_linear
    minus_operator_is_linear = .false.
  end function minus_operator_is_linear

  function lvalue_operator_is_linear(this)
    implicit none
    class(lvalue_operator_t), intent(in) :: this
    logical :: lvalue_operator_is_linear
    lvalue_operator_is_linear = .false.
  end function lvalue_operator_is_linear
  
end module operator_names
