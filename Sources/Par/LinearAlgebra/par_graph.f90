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
module par_graph_names
  ! Serial modules
  use types
  use memor
  use fem_graph_names
  use fem_partition_names
#ifdef memcheck       
  use iso_c_binding
#endif

  ! Module associated with the F90 interface to Trilinos.
  ! Remember: the F90 interface to Trilinos requires C
  ! interoperability (i.e., iso_c_binding module)
  !use for_trilinos_shadow_interfaces

  ! Parallel modules
  use par_partition_names
  use par_context_names

# include "debug.i90"

! *** IMPORTANT NOTE: This cpp macro should go to a 
! common include file or it should be a program 
! subroutine otherwise
#define blk2scal(iv,idof,ndof) (((iv)-1)*(ndof)+(idof))
  
  implicit none
  private

  ! Distributed Graph (eb or vb)
  type par_graph
       ! GENERAL COMMENT: Pointers and/or instances ? 
       ! I will follow the following rules:
       ! 
       !    x I will declare as a pointer, data which
       !      are potentially common to several objects,
       !      such as the partition or the communicator (to
       !      avoid data replication).

       !    x I will declare as instances objects which are 
       !      particular to the containing instance, such as graph or 
       !      epetra_crsgraph below 
       !
       ! Agreed ? Pros, cons ... ?
       !type( epetra_crsgraph ) :: epg

       ! Data structure which stores the local part 
       ! of the graph mapped to the current processor.
       ! This is required for both eb and vb data 
       ! distributions
       type( fem_graph )       :: f_graph

       ! Parallel partition control info.
       type ( par_partition ), pointer  :: p_part      => NULL()
       type ( par_partition ), pointer  :: p_part_cols => NULL()

  end type par_graph    

  interface par_graph_create
     module procedure par_graph_create_square, par_graph_create_rectangular
  end interface par_graph_create

  interface par_graph_free
     module procedure par_graph_free_one_shot, par_graph_free_progressively
  end interface par_graph_free

  ! Types
  public :: par_graph

  ! Functions
  public :: par_graph_create, par_graph_free, par_graph_print !, par_graph_epetra_crsgraph_create

!***********************************************************************
! Allocatable arrays of type(par_graph)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_graph)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !=============================================================================
  subroutine par_graph_create_square ( p_part, p_graph )
    implicit none 
    ! Parameters
    type(par_partition), target, intent(in)  :: p_part
    type(par_graph)            , intent(out) :: p_graph
    p_graph%p_part => p_part
    p_graph%p_part_cols => p_part
  end subroutine par_graph_create_square

    !=============================================================================
  subroutine par_graph_create_rectangular ( p_part, p_part_cols, p_graph )
    implicit none 
    ! Parameters
    type(par_partition), target, intent(in)  :: p_part
    type(par_partition), target, intent(in)  :: p_part_cols
    type(par_graph)            , intent(out) :: p_graph
    p_graph%p_part => p_part
    p_graph%p_part_cols => p_part_cols
  end subroutine par_graph_create_rectangular

  subroutine par_graph_free_one_shot(p_graph)
    !-----------------------------------------------------------------------
    ! This routine
    !-----------------------------------------------------------------------
    implicit none
    type(par_graph), intent(inout)  :: p_graph
    
    call par_graph_free_progressively (p_graph, free_only_struct)
    call par_graph_free_progressively (p_graph, free_clean )
  end subroutine par_graph_free_one_shot

  subroutine par_graph_free_progressively(p_graph, mode)
    !-----------------------------------------------------------------------
    ! This routine
    !-----------------------------------------------------------------------
    implicit none
    type(par_graph), intent(inout)  :: p_graph
    integer(ip)    , intent(in)     :: mode
    
    ! p_graph%p_part is required within this subroutine
    assert ( associated(p_graph%p_part) )
    
    ! p_graph%p_part%p_context is required within this subroutine
    assert ( associated(p_graph%p_part%p_context) )
    
    assert ( mode == free_only_struct .or. mode == free_clean  ) 
    
    if(p_graph%p_part%p_context%iam<0) return
    
    if ( mode == free_only_struct ) then
       call fem_graph_free ( p_graph%f_graph )
    else if ( mode == free_clean ) then
       ! Nullify parallel partition
       nullify( p_graph%p_part ) 
       nullify( p_graph%p_part_cols ) 
    end if
  end subroutine par_graph_free_progressively

  !=============================================================================
  subroutine par_graph_print(lunou, p_graph)
    implicit none
    type(par_graph)  ,  intent(in) :: p_graph
    integer(ip)      ,  intent(in) :: lunou
    
    ! p_graph%p_part is required within this subroutine
    assert ( associated(p_graph%p_part) )
    
    ! p_graph%p_part%p_context is required within this subroutine
    assert ( associated(p_graph%p_part%p_context) )
    
    call fem_graph_print (lunou, p_graph%f_graph)
  end subroutine par_graph_print

end module par_graph_names
