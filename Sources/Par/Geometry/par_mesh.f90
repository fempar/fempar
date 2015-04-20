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
  use types
  use fem_mesh_names
  use fem_mesh_io
  use stdio
  use psb_penv_mod

  use par_io
  use par_partition_names
  use par_context_names

# include "debug.i90"
  implicit none
  private

  ! Distributed mesh
  type par_mesh
     type(fem_mesh)                 :: f_mesh
     type(par_partition), pointer   :: p_part
  end type par_mesh

  interface par_mesh_create
     module procedure par_mesh_create_new, par_mesh_create_old
  end interface par_mesh_create

  interface par_mesh_free
     module procedure par_mesh_free_progressively, par_mesh_free_one_shot
  end interface par_mesh_free

  ! Types
  public :: par_mesh

  ! Functions
  public :: par_mesh_free, par_mesh_create, par_mesh_bcast

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

  
  subroutine par_mesh_free_progressively (p_mesh, mode)
    !-----------------------------------------------------------------------
    ! This routine deallocates a partition object
    !-----------------------------------------------------------------------
    implicit none

    ! Parameters
    type(par_mesh), intent(inout)  :: p_mesh
    integer(ip)   , intent(in)     :: mode

    assert ( associated(p_mesh%p_part%p_context) )
    assert ( p_mesh%p_part%p_context%created .eqv. .true.)
    assert ( mode == free_clean .or. mode == free_only_struct )

    if(p_mesh%p_part%p_context%iam<0) return

    if ( mode == free_clean ) then
       nullify (p_mesh%p_part)
    else if ( mode == free_only_struct ) then
       call fem_mesh_free ( p_mesh%f_mesh )
    end if
  end subroutine par_mesh_free_progressively
  


  !=============================================================================
  subroutine par_mesh_create_new ( p_part, p_mesh )
    implicit none 
    ! Parameters
    type(par_partition), target, intent(in)  :: p_part
    type(par_mesh)             , intent(out) :: p_mesh
    p_mesh%p_part => p_part
  end subroutine par_mesh_create_new

  subroutine par_mesh_read ( dir_path, prefix, p_mesh )
    implicit none 
    ! Parameters
    character(*), intent(in)      :: dir_path
    character(*), intent(in)      :: prefix
    type(par_mesh), intent(inout) :: p_mesh

    ! Locals
    integer(ip)                   :: lunio
    character(len=:), allocatable :: name 

    assert ( associated(p_mesh%p_part%p_context) )
    assert ( p_mesh%p_part%p_context%created .eqv. .true.)

    if(p_mesh%p_part%p_context%iam>=0) then
       call fem_mesh_compose_name ( prefix, name )
       call par_filename( p_mesh%p_part%p_context, name)
       ! Read mesh
       lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' )
       ! Read fem_partition data from path_file file
       call fem_mesh_read ( lunio, p_mesh%f_mesh )
       call io_close(lunio)
    end if
    ! Transfer ndime from fine-grid tasks to coarse-grid task
    call par_mesh_bcast (p_mesh, p_mesh%f_mesh%ndime)

  end subroutine par_mesh_read

  !=============================================================================
  subroutine par_mesh_create_old ( dir_path, prefix, p_part, p_mesh )
    implicit none 
    ! Parameters
    character (*)              , intent(in)  :: dir_path
    character (*)              , intent(in)  :: prefix
    type(par_partition), target, intent(in)  :: p_part
    type(par_mesh)             , intent(out) :: p_mesh

    ! Locals
    integer         :: iam, num_procs
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, lunio
    character(len=:), allocatable  :: name 
    logical                        :: flag

    assert ( associated(p_part%p_context) )
    assert ( p_part%p_context%created .eqv. .true.)

    p_mesh%p_part => p_part
    if(p_part%p_context%iam>=0) then
       call fem_mesh_compose_name ( prefix, name )
       call par_filename( p_part%p_context, name )
       
       ! Read mesh
       lunio = io_open( trim(dir_path) // '/' // trim(name), 'read' )
       
       ! Read fem_partition data from path_file file
       flag = .true. ! Assuming that the partitioned meshes are already permuted
       call fem_mesh_read ( lunio, p_mesh%f_mesh, flag )
       
       call io_close(lunio)
    end if

    ! Transfer ndime from fine-grid tasks to coarse-grid task
    call par_mesh_bcast (p_mesh, p_mesh%f_mesh%ndime)

  end subroutine par_mesh_create_old

  !=============================================================================
  subroutine par_mesh_bcast (p_mesh, value)
    implicit none
    ! Parameters 
    type(par_mesh), intent(in)    :: p_mesh    
    integer(ip)   , intent(inout) :: value
    assert ( associated(p_mesh%p_part))
    call par_partition_bcast(p_mesh%p_part,value)
  end subroutine par_mesh_bcast

end module par_mesh_names
