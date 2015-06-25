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
module par_mesh_names
  ! Serial modules
use types_names
  use fem_mesh_names
  use fem_mesh_distribution_names
use fem_mesh_io_names
use stdio_names

  ! Parallel modules
use par_io_names
  use par_environment_names

# include "debug.i90"
  implicit none
  private

  ! Distributed mesh
  type par_mesh
     type(fem_mesh)                 :: f_mesh
     type(fem_mesh_distribution)    :: f_mesh_dist
     type(par_environment), pointer :: p_env
  end type par_mesh

  interface par_mesh_free
     module procedure par_mesh_free_progressively, par_mesh_free_one_shot
  end interface par_mesh_free

  ! Types
  public :: par_mesh

  ! Functions
  public :: par_mesh_free, par_mesh_create, par_mesh_read

contains

  !=============================================================================
  subroutine par_mesh_free_one_shot (p_mesh)
    !-----------------------------------------------------------------------
    ! This routine deallocates a partition object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(par_mesh), intent(inout)  :: p_mesh
    call par_mesh_free_progressively(p_mesh, free_only_struct)
    call par_mesh_free_progressively(p_mesh, free_clean)
  end subroutine par_mesh_free_one_shot


  !============================================================================!
  ! This subroutine is provided for those contexts where the parallel mesh was ! 
  ! created incrementally (as, e.g., coarse-grid mesh in MLBDDC)               !
  !============================================================================!
  subroutine par_mesh_free_progressively (p_mesh, mode)
    implicit none

    ! Parameters
    type(par_mesh), intent(inout)  :: p_mesh
    integer(ip)   , intent(in)     :: mode

    assert ( mode == free_clean .or. mode == free_only_struct )

    ! Parallel environment MUST BE already created
    assert ( associated(p_mesh%p_env) )
    assert ( p_mesh%p_env%created )

    if(p_mesh%p_env%p_context%iam<0) return

    if ( mode == free_clean ) then
       nullify (p_mesh%p_env)
    else if ( mode == free_only_struct ) then
       call fem_mesh_free ( p_mesh%f_mesh )
       call fem_mesh_distribution_free ( p_mesh%f_mesh_dist )
    end if
  end subroutine par_mesh_free_progressively
  
  !============================================================================!
  ! This subroutine is provided for those contexts where the parallel mesh is  ! 
  ! created incrementally (as, e.g., coarse-grid mesh in MLBDDC)               !
  !============================================================================!
  subroutine par_mesh_create ( p_env, p_mesh )
    implicit none 
    ! Parameters
    type(par_environment), target, intent(in)  :: p_env
    type(par_mesh)               , intent(out) :: p_mesh
    
    ! Parallel environment MUST BE already created
    assert ( p_env%created )
    
    p_mesh%p_env => p_env
  end subroutine par_mesh_create

  !=============================================================================
  subroutine par_mesh_read ( dir_path, prefix, p_env, p_mesh )
    implicit none 
    ! Parameters
    character (*)                , intent(in)  :: dir_path
    character (*)                , intent(in)  :: prefix
    type(par_environment), target, intent(in)  :: p_env
    type(par_mesh)               , intent(out) :: p_mesh

    ! Locals
    integer                        :: iam, num_procs
    integer(ip)                    :: j, ndigs_iam, ndigs_num_procs, lunio
    character(len=:), allocatable  :: name 

    ! Parallel environment MUST BE already created
    assert ( p_env%created )

    p_mesh%p_env => p_env
    if(p_env%p_context%iam>=0) then
       call fem_mesh_compose_name ( prefix, name )
       call par_filename( p_mesh%p_env%p_context, name )
       ! Read mesh
       lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' )
       call fem_mesh_read_file ( lunio, p_mesh%f_mesh, permute_c2z = .false. )
       call io_close(lunio)

       call fem_mesh_distribution_compose_name ( prefix, name )
       call par_filename( p_mesh%p_env%p_context, name )
       ! Read mesh distribution control data
       lunio = io_open (trim(dir_path) // '/' // trim(name))
       call fem_mesh_distribution_read ( lunio, p_mesh%f_mesh_dist )
       call io_close(lunio)
    end if

    ! Transfer ndime from fine-grid tasks to coarse-grid task
    ! *** MISSING ISSUE HERE: IMPORTANT !!!! call par_mesh_bcast (p_mesh, p_mesh%f_mesh%ndime)
  end subroutine par_mesh_read

end module par_mesh_names
