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
module inverse_operator_names
  use types_names
  use base_operator_names
  use base_operand_names
  use abstract_solver_names
  use abstract_environment_names
  implicit none
# include "debug.i90"
  private
  
  type, extends(base_operator_t) :: inverse_operator_t
     class(base_operator_t)       , pointer :: A => NULL()     ! System operator
     class(base_operator_t)       , pointer :: M => NULL()     ! Preconditioner operator
     type(solver_control_t)       , pointer :: sctrl => NULL() ! Solver parameters
     class(abstract_environment_t), pointer :: env => NULL()   ! Serial/parallel environment 
   contains
     procedure :: create      => inverse_operator_create     
     procedure :: apply       => inverse_operator_apply      
     procedure :: apply_fun   => inverse_operator_apply_fun  
     procedure :: fill_values => inverse_operator_fill_values
     procedure :: free_values => inverse_operator_free_values
     procedure :: free        => inverse_operator_free       
  end type inverse_operator_t

  ! Types
  public :: inverse_operator_t

contains

  !==================================================================================================
  subroutine inverse_operator_create(op,A,M,sctrl,env)
    implicit none
    class(inverse_operator_t)            , intent(inout) :: op
    class(base_operator_t)       , target, intent(in)    :: A
    class(base_operator_t)       , target, intent(in)    :: M
    type(solver_control_t)       , target, intent(in)    :: sctrl
    class(abstract_environment_t), target, intent(in)    :: env

    ! Assign pointers
    op%A     => A
    op%M     => M
    op%sctrl => sctrl
    op%env   => env

  end subroutine inverse_operator_create

  !==================================================================================================
  subroutine inverse_operator_free(this)
    implicit none
    class(inverse_operator_t), intent(inout) :: this

    ! Unassign pointers
    this%A     => null()
    this%M     => null()
    this%sctrl => null()
    this%env   => null()

  end subroutine inverse_operator_free
    
  !==================================================================================================
  subroutine inverse_operator_apply(op,x,y)
    implicit none
    class(inverse_operator_t), intent(in)    :: op
    class(base_operand_t)    , intent(in)    :: x
    class(base_operand_t)    , intent(inout) :: y

    ! Checks
    check(associated(op%A))
    check(associated(op%M))
    check(associated(op%sctrl))
    check(associated(op%env))

    call x%GuardTemp()
    call abstract_solve(op%A,op%M,x,y,op%sctrl,op%env)
    call x%CleanTemp()    

  end subroutine inverse_operator_apply
  
  !==================================================================================================
  function inverse_operator_apply_fun(op,x) result(y)
    implicit none
    class(inverse_operator_t), intent(in) :: op
    class(base_operand_t)    , intent(in) :: x
    ! Locals
    class(base_operand_t), allocatable :: y

    call x%GuardTemp()
    allocate(y, mold=x); call y%default_initialization()
    call op%apply(x,y)
    call x%CleanTemp()    
    call y%SetTemp()
    
  end function inverse_operator_apply_fun

  !==================================================================================================
  subroutine inverse_operator_fill_values(op,stage)
    implicit none
    class(inverse_operator_t), intent(inout) :: op
    integer(ip), optional    , intent(in)    :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage

    if(op%fill_values_stage==stage_) call op%M%fill_values(stage_)

  end subroutine inverse_operator_fill_values

  !==================================================================================================
  subroutine inverse_operator_free_values(op,stage)
    implicit none
    class(inverse_operator_t), intent(inout) :: op
    integer(ip), optional , intent(in)       :: stage
    ! Locals
    integer(ip) :: stage_
    
    stage_ = update_nonlinear
    if(present(stage)) stage_ = stage
    
    if(op%free_values_stage==stage_) call op%M%free_values(stage_)

  end subroutine inverse_operator_free_values

end module inverse_operator_names

  
