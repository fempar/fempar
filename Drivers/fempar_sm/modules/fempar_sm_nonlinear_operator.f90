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
module fempar_sm_nonlinear_operator_names
  use fempar_names
  use fempar_sm_discrete_integration_names
  implicit none
# include "debug.i90"
  private

  type, extends(operator_t) :: nonlinear_operator_t
     private
     type(fe_affine_operator_t)             , pointer :: residual  => null()
     type(fe_function_t)                              :: argument
     class(fempar_sm_discrete_integration_t), pointer :: discrete_integration
   contains
     procedure :: create           => nonlinear_operator_create
     procedure :: apply            => nonlinear_operator_apply
     procedure :: apply_add        => nonlinear_operator_apply_add
     procedure :: is_linear        => nonlinear_operator_is_linear
     procedure :: get_tangent      => nonlinear_operator_get_tangent
     procedure :: get_translation  => nonlinear_operator_get_translation
     procedure :: free             => nonlinear_operator_free
  end type nonlinear_operator_t

  public :: nonlinear_operator_t

contains

  subroutine nonlinear_operator_create(this, residual)
    implicit none
    class(nonlinear_operator_t)          , intent(inout) :: this
    type(fe_affine_operator_t)   , target, intent(in)    :: residual
    class(discrete_integration_t), pointer :: discrete_integration
    class(serial_fe_space_t), pointer :: fe_space

    ! Create a nonlinear fe operator
    this%residual => residual
    fe_space => this%residual%get_fe_space()
    call this%argument%create(fe_space)     ! tmp workspace

    discrete_integration => this%residual%get_discrete_integration()
    select type(discrete_integration)
    class is(fempar_sm_discrete_integration_t)
       this%discrete_integration => discrete_integration
       call this%discrete_integration%set_solution(this%argument)
    class default
       mcheck(.false.,'nonlinear operator only works with fempar_sm_discrete_integration')
    end select

  end subroutine nonlinear_operator_create

  subroutine nonlinear_operator_free(this)
    implicit none
    class(nonlinear_operator_t), intent(inout) :: this
    call this%argument%free()
    this%residual => null()
    this%discrete_integration => null()
  end subroutine nonlinear_operator_free
  
  !=================================================================================================
  
  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine nonlinear_operator_apply(this,x,y) 
    implicit none
    class(nonlinear_operator_t), intent(in)   :: this
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 
    class(vector_t) , pointer :: tmp
    tmp => this%argument%get_dof_values()
    tmp = x
    call this%discrete_integration%set_terms_to_integrate(translation_terms)
    call this%residual%numerical_setup()
    tmp => this%residual%get_translation()
    y = tmp 
  end subroutine nonlinear_operator_apply

  ! op%apply_add(x,y) <=> y <- op*x+y
  ! Implicitly assumes that y is already allocated
  subroutine nonlinear_operator_apply_add(this,x,y) 
    implicit none
    class(nonlinear_operator_t), intent(in)    :: this
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 
    class(vector_t) , pointer :: tmp
    tmp => this%argument%get_dof_values()
    tmp = x
    call this%discrete_integration%set_terms_to_integrate(translation_terms)
    call this%residual%numerical_setup()
    tmp => this%residual%get_translation()
    y = y + tmp 
  end subroutine nonlinear_operator_apply_add

  function nonlinear_operator_is_linear(this) 
    implicit none
    class(nonlinear_operator_t), intent(in)    :: this
    logical                          :: nonlinear_operator_is_linear
    nonlinear_operator_is_linear = .false.
  end function nonlinear_operator_is_linear

  function nonlinear_operator_get_tangent(this,x) result(tangent)
    implicit none
    class(nonlinear_operator_t)          , intent(in) :: this
    class(vector_t)  , optional, intent(in) :: x
    type(lvalue_operator_t)          :: tangent 
    class(vector_t) , pointer :: tmp
    if(present(x)) then
       tmp => this%argument%get_dof_values()
       tmp = x
       call this%discrete_integration%set_terms_to_integrate(tangent_terms)
       call this%residual%numerical_setup()
    end if
    tangent = this%residual%get_tangent()
  end function nonlinear_operator_get_tangent

  function nonlinear_operator_get_translation(this) result(translation)
    implicit none
    class(nonlinear_operator_t) , intent(in) :: this
    class(vector_t)   , pointer    :: translation
    translation => this%residual%get_translation()
  end function nonlinear_operator_get_translation
     
 end module fempar_sm_nonlinear_operator_names
