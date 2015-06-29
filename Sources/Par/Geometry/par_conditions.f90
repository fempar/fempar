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
module par_conditions_names
  ! Serial modules
  use types_names
  use memor_names
  use stdio_names
  use renumbering_names
  use conditions_names
  use conditions_io_names

  ! Parallel modules
  use par_environment_names
  use par_io_names

# include "debug.i90"
  implicit none
  private

  type par_conditions_t
     ! Data structure which stores the local part
     ! of the BC's mapped to the current processor
     type(conditions_t) :: f_conditions

     ! Parallel environment control
     type(par_environment_t), pointer :: p_env => NULL()
  end type par_conditions_t

  ! Types
  public :: par_conditions_t

  ! Methods
  public :: par_conditions_create, par_conditions_free     , & 
            par_conditions_copy, par_conditions_apply_renumbering, &
            par_conditions_read

contains

  !===============================================================================================
  subroutine par_conditions_create(ncode,nvalu,ncond,p_env,cnd)
    implicit none
    integer(ip)                  , intent(in)    :: ncode, nvalu, ncond
    type(par_environment_t), target, intent(in)    :: p_env
    type(par_conditions_t)         , intent(inout) :: cnd

    ! Parallel environment MUST BE already created
    assert ( p_env%created )

    cnd%p_env => p_env

    if( p_env%p_context%iam >= 0 ) then
       call conditions_create ( ncode,nvalu,ncond,cnd%f_conditions)
    end if

  end subroutine par_conditions_create

  !===============================================================================================
  subroutine par_conditions_copy(cnd_old,cnd_new)
    implicit none
    type(par_conditions_t), target, intent(in)    :: cnd_old
    type(par_conditions_t)        , intent(inout) :: cnd_new

    ! Parallel environment MUST BE already created
    assert ( cnd_old%p_env%created )
    
    cnd_new%p_env => cnd_old%p_env

    if( cnd_new%p_env%p_context%iam >= 0 ) then
       call conditions_copy ( cnd_old%f_conditions, cnd_new%f_conditions )
    end if

  end subroutine par_conditions_copy

  !===============================================================================================
  subroutine par_conditions_apply_renumbering(renumbering, cnd)
    implicit none
    type(renumbering_t)         ,    intent(in)  :: renumbering
    type(par_conditions_t), intent(inout)  :: cnd
    
    ! Parallel environment MUST BE already created
    assert ( cnd%p_env%created )
    
    if( cnd%p_env%p_context%iam >= 0 ) then
       call conditions_apply_renumbering ( renumbering, cnd%f_conditions )
    end if

  end subroutine par_conditions_apply_renumbering

  !===============================================================================================
  subroutine par_conditions_free(cnd)
    implicit none
    type(par_conditions_t), intent(inout) :: cnd

    ! Parallel environment MUST BE already created
    assert ( cnd%p_env%created )
    
    if( cnd%p_env%p_context%iam >= 0 ) then
       call conditions_free ( cnd%f_conditions )
    end if
    
  end subroutine par_conditions_free

  !=============================================================================
  subroutine par_conditions_read ( dir_path, prefix, npoin, p_env, p_conditions )
    implicit none 
    ! Parameters
    character (*)                , intent(in)  :: dir_path
    character (*)                , intent(in)  :: prefix
    integer(ip)                  , intent(in)  :: npoin
    type(par_environment_t), target, intent(in)  :: p_env
    type(par_conditions_t)         , intent(out) :: p_conditions
    
    ! Locals
    character(len=:), allocatable  :: name
    integer(ip) :: lunio
    
    ! Parallel environment MUST BE already created
    assert ( p_env%created )
    
    p_conditions%p_env => p_env
    if(p_env%p_context%iam>=0) then
       call conditions_compose_name ( prefix, name )
       call par_filename( p_conditions%p_env%p_context, name )
       
       ! Read conditions
       lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' )
       call conditions_read_file ( lunio, npoin, p_conditions%f_conditions )
       call io_close(lunio)
    end if
    
  end subroutine par_conditions_read
  
end module par_conditions_names
