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
       if ( p_graph%p_part%p_context%handler == trilinos ) then 
          !call epetra_crsgraph_destruct  ( p_graph%epg )
       end if
       
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
    
    if ( p_graph%p_part%p_context%handler == trilinos ) then 
        !call epetra_crsgraph_print ( p_graph%epg )
    end if

  end subroutine par_graph_print

  ! !=============================================================================
  ! ! Creates and fills epetra_crsgraph p_graph%epg 
  ! ! as a view of the local graph f_graph
  ! !============================================================================= 
  ! subroutine par_graph_epetra_crsgraph_create (p_graph, storage, ndof1, ndof2)
  !   implicit none

  !   ! Parameters 
  !   type(par_graph), intent(inout) :: p_graph
  !   integer(ip)    , intent(in)    :: storage, ndof1, ndof2

  !   ! Local variables
  !   integer (c_int)              :: ierrc
  !   integer (c_int), allocatable :: nnzs(:)
  !   integer (c_int), allocatable :: rc_map_values(:)
  !   integer (ip)                 :: id, idof


  !   ! Pointer to part/context object is required
  !   assert ( associated(p_graph%p_part) )
  !   assert ( associated(p_graph%p_part%p_context) )
  !   ! Vertex-based partitioning/trilinos handler required for Epetra
  !   assert ( p_graph%p_part%f_part%ptype == vertex_based )
  !   assert ( p_graph%p_part%p_context%handler == trilinos )

  !   ! Zero-based indexing for the column identifiers required by Epetra
  !   call fem_graph_ja_one_to_zero_indexing ( p_graph%f_graph )

  !   if ( storage == blk ) then
  !      if ( p_graph%p_part%maps_state (1) == map_non_created ) then

  !         call memalloc ( p_graph%p_part%f_part%nmap%nl, rc_map_values,          __FILE__,__LINE__)

  !         do id=1, p_graph%p_part%f_part%nmap%nl
  !            rc_map_values(id) = p_graph%p_part%f_part%nmap%l2g(id)-1
  !         end do

  !         ! Create row-map
  !         call epetra_map_construct ( p_graph%p_part%row_map(1), -1, & 
  !              &     p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb, & 
  !              &     rc_map_values, 0, p_graph%p_part%p_context%epcomm)
  !         ! Create col-map
  !         call epetra_map_construct ( p_graph%p_part%col_map(1), -1, p_graph%p_part%f_part%nmap%nl, &
  !              &       rc_map_values, 0, p_graph%p_part%p_context%epcomm) 
  !         ! Create importer                                                 TARGET                 SOURCE
  !         call epetra_import_construct ( p_graph%p_part%importer(1), p_graph%p_part%col_map(1), p_graph%p_part%row_map(1) )

  !         call memfree ( rc_map_values,__FILE__,__LINE__)

  !         p_graph%p_part%maps_state(1) = map_created
  !      end if

  !      ! Allocate space for the number of neighbours of each vertex
  !      call memalloc ( p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb, nnzs,          __FILE__,__LINE__)

  !      ! Compute the number of neighbours of each vertex
  !      call compute_nnzs ( p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb, & 
  !           &   p_graph%f_graph%ia, nnzs )

  !      ! Create epetra_crsgraph object with static profile=true, mode='View', with
  !      ! provided row and column maps, and variable nnzs per row
  !      call epetra_crsgraph_construct  ( p_graph%epg, p_graph%p_part%row_map(1), p_graph%p_part%col_map(1), nnzs )

  !      ! Signal to inform that we are going to provide local indices
  !      call epetra_crsgraph_fill_complete ( p_graph%epg, ierrc ) 
  !      assert (ierrc == 0)

  !      ! Insert row/column entries information via calls to InsertMyValues
  !      call insert_my_indices ( p_graph%epg, & 
  !           &   p_graph%p_part%f_part%nmap%ni + & 
  !           &   p_graph%p_part%f_part%nmap%nb,  &
  !           &   p_graph%f_graph%ia,             &
  !           &   p_graph%f_graph%ja )  

  !      ! Signal to inform that we are done 
  !      call epetra_crsgraph_fill_complete ( p_graph%epg, ierrc ); assert (ierrc == 0)

  !      call memfree (nnzs,__FILE__,__LINE__)

  !   else if ( storage == scal ) then
  !      assert ( ndof1 <= max_ndofs )
  !      if ( p_graph%p_part%maps_state (ndof1) == map_non_created ) then
  !         call memalloc ( p_graph%p_part%f_part%nmap%nl*ndof1, rc_map_values,         __FILE__,__LINE__)

  !         do id=1, p_graph%p_part%f_part%nmap%nl 
  !            do idof=1, ndof1
  !               rc_map_values(blk2scal(id,idof,ndof1)) = blk2scal(p_graph%p_part%f_part%nmap%l2g(id),idof,ndof1)-1
  !            end do
  !         end do

  !         ! Create row-map
  !         call epetra_map_construct ( p_graph%p_part%row_map(ndof1), -1, & 
  !              &     (p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb)*ndof1, & 
  !              &     rc_map_values, 0, p_graph%p_part%p_context%epcomm)
  !         ! Create col-map
  !         call epetra_map_construct ( p_graph%p_part%col_map(ndof1), -1, p_graph%p_part%f_part%nmap%nl*ndof1, &
  !              &       rc_map_values, 0, p_graph%p_part%p_context%epcomm) 
  !         ! Create importer                                                      TARGET                            SOURCE
  !         call epetra_import_construct ( p_graph%p_part%importer(ndof1), p_graph%p_part%col_map(ndof1), p_graph%p_part%row_map(ndof1) )

  !         call memfree ( rc_map_values,__FILE__,__LINE__)

  !         p_graph%p_part%maps_state(ndof1) = map_created
  !      end if
  !      assert ( ndof2 <= max_ndofs )
  !      if ( p_graph%p_part%maps_state (ndof2) == map_non_created ) then
  !         call memalloc ( p_graph%p_part%f_part%nmap%nl*ndof2, rc_map_values,     __FILE__,__LINE__)

  !         do id=1, p_graph%p_part%f_part%nmap%nl 
  !            do idof=1, ndof2
  !               rc_map_values(blk2scal(id,idof,ndof2)) = blk2scal(p_graph%p_part%f_part%nmap%l2g(id),idof,ndof2)-1
  !            end do
  !         end do

  !         ! Create row-map
  !         call epetra_map_construct ( p_graph%p_part%row_map(ndof2), -1, & 
  !              &     (p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb)*ndof2, & 
  !              &     rc_map_values, 0, p_graph%p_part%p_context%epcomm)
  !         ! Create col-map
  !         call epetra_map_construct ( p_graph%p_part%col_map(ndof2), -1, p_graph%p_part%f_part%nmap%nl*ndof2, &
  !              &       rc_map_values, 0, p_graph%p_part%p_context%epcomm) 
  !         ! Create importer                                                        
  !         call epetra_import_construct ( p_graph%p_part%importer(ndof2), &
  !         !                                           TARGET                          SOURCE
  !                                      & p_graph%p_part%col_map(ndof2), p_graph%p_part%row_map(ndof2) )

  !         call memfree ( rc_map_values,__FILE__,__LINE__)

  !         p_graph%p_part%maps_state(ndof2) = map_created
  !      end if

  !      ! Allocate space for the number of neighbours of each vertex
  !      call memalloc ( (p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb)*ndof1, nnzs, __FILE__,__LINE__)

  !      ! Compute the number of neighbours of each vertex
  !      call compute_nnzs ( (p_graph%p_part%f_part%nmap%ni + p_graph%p_part%f_part%nmap%nb)*ndof1, & 
  !           &   p_graph%f_graph%ia, nnzs )

  !      ! Create epetra_crsgraph object with static profile=true, mode='View', with
  !      ! provided row and column maps, and variable nnzs per row
  !      call epetra_crsgraph_construct  ( p_graph%epg, p_graph%p_part%row_map(ndof1), & 
  !                                      & p_graph%p_part%col_map(ndof2),  nnzs )

  !      ! Signal to inform that we are going to provide local indices
  !      call epetra_crsgraph_fill_complete ( p_graph%epg, p_graph%p_part%row_map(ndof2), & 
  !                                         & p_graph%p_part%row_map(ndof1), ierrc ) 
  !      assert (ierrc == 0)

  !      ! Insert row/column entries information via calls to InsertMyValues
  !      call insert_my_indices ( p_graph%epg,           & 
  !           &   (p_graph%p_part%f_part%nmap%ni +       &  
  !           &   p_graph%p_part%f_part%nmap%nb)*ndof1,  &
  !           &   p_graph%f_graph%ia,                    &
  !           &   p_graph%f_graph%ja )  

  !      ! Signal to inform that we are done                       
  !      call epetra_crsgraph_fill_complete ( p_graph%epg, p_graph%p_part%row_map(ndof2), &
  !                                         & p_graph%p_part%row_map(ndof1), ierrc ); assert (ierrc == 0)

  !      call memfree (nnzs,__FILE__,__LINE__)

  !   end if

  ! end subroutine par_graph_epetra_crsgraph_create

  !  subroutine compute_nnzs (n, ia, nnzs)
  !    ! Parameters
  !    integer (c_int) , intent(in)  :: n
  !    integer (c_int) , intent(in)  :: ia(n+1)
  !    integer (c_int) , intent(out) :: nnzs(n)
  !
  !    ! Local variables
  !    integer (c_int) :: i
  !
  !    ! Compute the number of neighbours of each vertex
  !    do i=1,n
  !       nnzs(i) = ia(i+1) - ia(i)
  !    end do
  !    
  !  end subroutine compute_nnzs

  ! ! Private routine which unpacks allocatable derived data type 
  ! ! members as explicit size arrays
  ! subroutine insert_my_indices (epg, n, ia, ja)
  !   ! Parameters
  !   !type(epetra_crsgraph) , intent(inout) :: epg
  !   integer (c_int)       , intent(in)    :: n
  !   integer (c_int)       , intent(in)    :: ia(n+1)
  !   integer (c_int)       , intent(in)    :: ja(ia(n+1)-1)

  !   ! Local variables
  !   integer (c_int) :: i, ierrc

  !   ! Insert row and column entries information via calls 
  !   ! to InsertMyIndices
  !   do i=1, n
  !       call epetra_crsgraph_insert_my_indices ( epg, i-1,              & 
  !            &                                   ia(i+1)-ia(i),         &  
  !            &                                   ja(ia(i):(ia(i+1)-1)), & 
  !            &                                   ierrc )
  !       assert (ierrc == 0)
  !   end do
  ! end subroutine insert_my_indices

end module par_graph_names
