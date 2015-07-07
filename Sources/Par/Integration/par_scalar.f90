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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module par_scalar_names
  use types_names
  use integrable_names
  use scalar_names
  use par_environment_names
  use psb_penv_mod_names
# include "debug.i90"
  implicit none
  private

  type, extends(integrable_t) :: par_scalar_t
     private
     type(scalar_t) :: f_scalar
     type(par_environment_t), pointer :: p_env => NULL()
  contains
    procedure :: create => par_scalar_create
    procedure :: free   => par_scalar_free
    procedure :: init   => par_scalar_init
    procedure :: sum    => par_scalar_sum
    procedure :: get    => par_scalar_get
    procedure :: reduce => par_scalar_reduce
  end type par_scalar_t

  ! Types
  public :: par_scalar_t

contains

  !==================================================================================================
  subroutine par_scalar_create (this,p_env)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine initialize the par_scalar.                                                  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t)            , intent(inout) :: this
    type(par_environment_t), target, intent(in)    :: p_env
    
    this%p_env => p_env

  end subroutine par_scalar_create

  !==================================================================================================
  subroutine par_scalar_free (this)
    !-----------------------------------------------------------------------------------------------!
    !   Dummy subroutine to specialize deferred procedure.                                          !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t), intent(inout) :: this
  end subroutine par_scalar_free

  !==================================================================================================
  subroutine par_scalar_init (this,b)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine initialize the par_scalar.                                                  !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t), intent(inout) :: this
    real(rp), optional , intent(in)    :: b
    
    assert(associated(this%p_env))

    if(this%p_env%am_i_fine_task()) then
       call this%f_scalar%init(b)
    end if

  end subroutine par_scalar_init

  !==================================================================================================
  subroutine par_scalar_sum (this,b)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine sum a real to the par_scalar.                                               !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t), intent(inout) :: this
    real(rp)       , intent(in)    :: b
    
    assert(associated(this%p_env))

    if(this%p_env%am_i_fine_task()) then
       call this%f_scalar%sum(b)
    end if

  end subroutine par_scalar_sum

  !==================================================================================================
  function par_scalar_get (this)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine returns the par_scalar value.                                               !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t), intent(inout) :: this
    real(rp)                       :: par_scalar_get
    
    assert(associated(this%p_env))

    if(this%p_env%am_i_fine_task()) then
       par_scalar_get = this%f_scalar%get()
    else
       par_scalar_get = 0.0_rp
    end if

  end function par_scalar_get

  !==================================================================================================
  subroutine par_scalar_reduce (this)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine returns the par_scalar value reduced at all fine tasks.                     !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(par_scalar_t), intent(inout) :: this
    ! Locals
    real(rp) :: a
    
    assert(associated(this%p_env))

    if(this%p_env%am_i_fine_task()) then
       a = this%f_scalar%get()
       call psb_sum ( this%p_env%p_context%icontxt, a )
       call this%f_scalar%init(a)
    end if

  end subroutine par_scalar_reduce

end module par_scalar_names
