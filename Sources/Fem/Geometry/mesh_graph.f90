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
module mesh_graph
  use types
  use memor
  use fem_mesh_names
  use fem_graph_names
  use renum_names
  implicit none
# include "debug.i90"
  private

  interface mesh_to_graph
     module procedure mesh_to_graph_part, mesh_to_graph_part_new
  end interface mesh_to_graph

  ! Functions
  public :: mesh_to_graph, mesh_graph_compute_connected_components

contains


  !================================================================================================
  subroutine mesh_to_graph_part (primal_mesh, dual_mesh, min_freq_neig, primal_graph)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! The routine exploits the duality among primal/dual meshes. In order to 
    ! work properly pnods has to be an array of pointers to address lnods
    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
    !
    ! min_freq_neig is an input integer which determines which are
    ! the neighbours of a particular vertex of the graph. In particular, 
    ! neighbours are those vertices which are listed at least min_freq_neig times 
    ! when traversing the list of primal elements around a primal point
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, mesh_to_graph_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
    integer(ip), intent(in)      :: min_freq_neig
    type(fem_graph), intent(out) :: primal_graph

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(4)            ! Pointers to work space

    ! Allocate space for ia on the primal graph
    primal_graph%type = part
    primal_graph%nv   = primal_mesh%npoin
    call memalloc (primal_graph%nv+1, primal_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_primal_graph and list_primal_graph routines
    ! (TOTAL WS SIZE = primal mesh npoin + 2*maximum number of neighbours of any primal
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + primal_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    pwork(4) = pwork(3) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(4), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_primal_graph' ! DBG:
    ! Calculate ia array
    call count_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
         &                           iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)), &
         &                           iwork(pwork(3):pwork(4)) ) 
    ! write(*,*) 'end count_primal_graph'   ! DBG:

    ! Allocate space for ja on the primal graph 
    call memalloc (primal_graph%ia(primal_graph%nv+1)-1, primal_graph%ja,           __FILE__,__LINE__)

    ! write(*,*) 'begin list_primal_graph' ! DBG:
    ! Calculate ja array
    call list_primal_graph_part  (primal_mesh, dual_mesh, min_freq_neig, primal_graph, & 
         &                        iwork(pwork(1):pwork(2)), iwork(pwork(2):pwork(3)),  &
         &                        iwork(pwork(3):pwork(4)) )
    ! write(*,*) 'end list_primal_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine mesh_to_graph_part

  !============================================================================================
  subroutine  count_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                  ws_position, ws_freq, ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    integer(ip), intent(in)         :: min_freq_neig
    type(fem_graph), intent(inout)  :: primal_graph
    integer(ip), intent(out)        :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)        :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                 :: ineigh  
    integer(ip)                 :: ipoindm ! Dual   mesh point identifier
    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
    integer(ip)                 :: ipoinpg ! Primal graph point identifier 
    integer(ip)                 :: p_ipoindm   
    integer(ip)                 :: p_ipoinpm   
    integer(ip)                 :: inods1 , inods2   
    integer(ip)                 :: inods1d, inods2d

    integer(ip)                 :: first_free_pos
    integer(ip)                 :: num_neighbors

    
    ! I do not agree with the following
    ! assertion. dual_mesh can actually be
    ! a primal finite element mesh when this
    ! routine has to compute dual graphs (i.e., elements 
    ! around elements). This last  situation was not taken 
    ! into account in the code below.
    ! I have modified the code below so that it
    ! is also able to compute dual graphs. 
    
    ! !!!! IMPORTANT NOTE: INCORRECT ASSERTION !!!
    ! assert(dual_mesh%nelty/=1)

    ! CORRECT ASSERTION: Either primal_mesh or dual_mesh has
    !                    to be a dual_mesh
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_freq)                                              ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    primal_graph%ia(1) = 1

    do ipoinpg=1, primal_mesh%npoin
       primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg)
       first_free_pos = 1
       num_neighbors  = 0
       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       if(dual_mesh%nelty==1) then
          inods1d=(ipoinpg-1)*dual_mesh%nnode+1
          inods2d=ipoinpg*dual_mesh%nnode
       else
          inods1d=dual_mesh%pnods(ipoinpg)
          inods2d=dual_mesh%pnods(ipoinpg+1)-1
       end if
       do p_ipoindm = inods1d, inods2d
       !do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm) 

             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq (ws_position(ipoinpm)) = ws_freq( ws_position(ipoinpm) ) + 1
             end if

          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 
       ! write (*, '(10i10)')    ws_freq     (1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if (.not.ws_neighbors(ineigh)==ipoinpg) then ! Exclude self-edges
                primal_graph%ia(ipoinpg+1) = primal_graph%ia(ipoinpg+1) + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

  end subroutine count_primal_graph_part

  !============================================================================================
  subroutine  list_primal_graph_part ( primal_mesh, dual_mesh, min_freq_neig, primal_graph,  &
       &                                 ws_position, ws_freq, ws_neighbors )
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: primal_graph
    integer(ip), intent(in)        :: min_freq_neig
    integer(ip), intent(out)       :: ws_position  (primal_mesh%npoin)
    integer(ip), intent(out)       :: ws_freq      (primal_mesh%nnode*dual_mesh%nnode)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoindm ! Dual   mesh point identifier
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoinpg ! Primal graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    ! I do not agree with the following
    ! assertion. dual_mesh can actually be
    ! a primal finite element mesh when this
    ! routine has to compute dual graphs (i.e., elements 
    ! around elements). This last  situation was not taken 
    ! into account in the code below.
    ! I have modified the code below so that it
    ! is also able to compute dual graphs. 
    
    ! !!!! IMPORTANT NOTE: INCORRECT ASSERTION !!!
    ! assert(dual_mesh%nelty/=1)

    ! CORRECT ASSERTION: Either primal_mesh or dual_mesh has
    !                    to be a dual_mesh
    assert ( primal_mesh%nelty /= 1 .or. dual_mesh%nelty /= 1 ) 


    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize work space of freqs. 
    ! associated to neighbours
    ws_freq       = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoinpg=1, primal_mesh%npoin       
       first_free_pos = 1
       num_neighbors  = 0

       ! Traverse the dual nodes of the dual element number ipoinpg
       ! (i.e., primal elements around primal node number ipoinpg)
       if(dual_mesh%nelty==1) then
          inods1d=(ipoinpg-1)*dual_mesh%nnode+1
          inods2d=ipoinpg*dual_mesh%nnode
       else
          inods1d=dual_mesh%pnods(ipoinpg)
          inods2d=dual_mesh%pnods(ipoinpg+1)-1
       end if

       do p_ipoindm = inods1d, inods2d
       ! do p_ipoindm = dual_mesh%pnods(ipoinpg), dual_mesh%pnods(ipoinpg+1)-1
          ipoindm = dual_mesh%lnods(p_ipoindm)

          ! Traverse the primal nodes of the primal element number ipoindm
          ! (i.e., dual elements around dual node number ipoindm)
          if(primal_mesh%nelty==1) then
             inods1=(ipoindm-1)*primal_mesh%nnode+1
             inods2=ipoindm*primal_mesh%nnode
          else
             inods1=primal_mesh%pnods(ipoindm)
             inods2=primal_mesh%pnods(ipoindm+1)-1
          end if
          do p_ipoinpm = inods1,inods2
             ipoinpm = primal_mesh%lnods(p_ipoinpm)
             ! If ipoinpm not visited yet
             if ( ws_position(ipoinpm) == 0 ) then
                ws_position(ipoinpm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoinpm
                ws_freq(first_free_pos) = 1
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             else
                ws_freq ( ws_position(ipoinpm) ) = ws_freq( ws_position(ipoinpm) ) + 1
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (ws_freq(ineigh).ge.min_freq_neig) then
             if ( .not.ws_neighbors(ineigh)==ipoinpg ) then ! Exclude self-edges
                primal_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
                first_free_ja_pos = first_free_ja_pos + 1
             end if
          end if
          ws_freq(ineigh) = 0
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_primal_graph_part

  !============================================================================================
  subroutine mesh_to_graph_part_new (primal_mesh, dual_mesh, dual_graph, only_faces)
    !---------------------------------------------------------------------------
    ! This routine generates a graph of a given primal mesh (this graph 
    ! is a list of primal points around primal points)
    !
    ! The routine exploits the duality among primal/dual meshes. In order to 
    ! work properly pnods has to be an array of pointers to address lnods
    ! (!!! this does not currently happen with FE meshes with eltyp==1 !!!)
    !
    ! min_freq_neig is an input integer which determines which are
    ! the neighbours of a particular vertex of the graph. In particular, 
    ! neighbours are those vertices which are listed at least min_freq_neig times 
    ! when traversing the list of primal elements around a primal point
    !
    ! IMPORTANT NOTE:  This routine DOES NOT generate self-edges for both primal 
    !                  and dual graphs. However, mesh_to_graph_dual
    !                  generates self-edges. METIS graphs do NOT have self-edges.
    !---------------------------------------------------------------------------
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)  :: primal_mesh, dual_mesh
    type(fem_graph), intent(out) :: dual_graph
    integer(ip), intent(in)      :: only_faces 

    ! Local variables
    integer(ip), allocatable :: iwork(:)            ! Integer ip working array
    integer(ip)              :: pwork(3)            ! Pointers to work space

    assert(only_faces==1)

    ! Allocate space for ia on the dual graph
    dual_graph%nv = dual_mesh%npoin
    call memalloc (dual_graph%nv+1, dual_graph%ia, __FILE__,__LINE__)

    ! Allocate working space for count_dual_graph and list_dual_graph routines
    ! (TOTAL WS SIZE = dual mesh npoin + 2*maximum number of neighbours of any dual
    ! graph node)
    pwork(1) = 1
    pwork(2) = pwork(1) + dual_mesh%npoin
    pwork(3) = pwork(2) + dual_mesh%nnode*primal_mesh%nnode
    call memalloc (pwork(3), iwork, __FILE__,__LINE__)

    ! write(*,*) 'begin count_dual_graph' ! DBG:
    ! Calculate ia array
    call count_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph,iwork(pwork(1):pwork(2)), &
         &                           iwork(pwork(2):pwork(3)))
    ! write(*,*) 'end count_dual_graph'   ! DBG:

    ! Allocate space for ja on the dual graph 
    call memalloc (dual_graph%ia(dual_graph%nv+1)-1, dual_graph%ja,__FILE__,__LINE__)

    ! write(*,*) 'begin list_dual_graph' ! DBG:
    ! Calculate ja array
    call list_dual_graph_part_new (primal_mesh, dual_mesh, dual_graph, iwork(pwork(1):pwork(2)), & 
         &                         iwork(pwork(2):pwork(3)))
    ! write(*,*) 'end list_dual_graph'   ! DBG: 

    ! Deallocate working space
    call memfree (iwork,__FILE__,__LINE__)
    return
  end subroutine mesh_to_graph_part_new

  !============================================================================================
  subroutine  count_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph, ws_position, &
       &                                  ws_neighbors)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)     :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout)  :: dual_graph
    integer(ip), intent(out)        :: ws_position  (dual_mesh%npoin)
    integer(ip), intent(out)        :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                 :: ineigh  
    integer(ip)                 :: ipoinpm ! Primal mesh point identifier
    integer(ip)                 :: ipoindm ! Dual mesh point identifier
    integer(ip)                 :: ipoindg ! Dual graph point identifier 
    integer(ip)                 :: p_ipoindm   
    integer(ip)                 :: p_ipoinpm   
    integer(ip)                 :: inods1 , inods2   
    integer(ip)                 :: inods1d, inods2d

    integer(ip)                 :: first_free_pos
    integer(ip)                 :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions
    ws_position   = 0

    ! Initialize work space of neighbour 
    ! identifiers 
    ws_neighbors  = 0

    !primal_graph%ia    = -1                                                ! DBG:
    !write (*,*) size(primal_graph%ia)                                      ! DBG: 
    !write (*,*) size(primal_mesh%pnods)                                    ! DBG: 
    !write (*,*) size(primal_mesh%lnods)                                    ! DBG: 
    !write (*,*) size(dual_mesh%pnods)                                      ! DBG: 
    !write (*,*) size(dual_mesh%lnods)                                      ! DBG: 
    !write (*,*) size(ws_position)                                          ! DBG: 
    !write (*,*) size(ws_neighbors)                                         ! DBG: 
    !write (*, '(10i10)')    primal_mesh%pnods  (1:(primal_mesh%nelem+1))   ! DBG: 
    !write (*, '(10i10)')    primal_graph%ia    (1:(primal_graph%nv+1))     ! DBG:  
    !primal_graph%ia    = -1 
    dual_graph%ia(1) = 1

    do ipoindg=1, dual_mesh%npoin
       dual_graph%ia(ipoindg+1) = dual_graph%ia(ipoindg)
       first_free_pos = 1
       num_neighbors  = 0
       inods1=primal_mesh%p_face(ipoindg)
       inods2=primal_mesh%pnods(ipoindg+1)-1
       do p_ipoinpm = inods1,inods2
          ipoinpm = primal_mesh%lnods(p_ipoinpm) 

          inods1d=dual_mesh%pnods(ipoinpm)
          inods2d=dual_mesh%pnods(ipoinpm+1)-1
          do p_ipoindm = inods1d, inods2d
             ipoindm = dual_mesh%lnods(p_ipoindm)

             ! write (*,*) ipoinpg, ipoindm, ipoinpm ! DBG: 
             ! If ipoinpm not visited yet
             if ( ws_position(ipoindm) == 0 ) then
                ws_position(ipoindm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoindm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
                ! If ipoinpm already visited
             end if

          end do
       end do

       ! write (*, '(10i10)')    ws_position (1:primal_graph%nv)  ! DBG:
       ! write (*, '(10i10)')    ws_neighbors(1:num_neighbors)    ! DBG: 

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if (.not.ws_neighbors(ineigh)==ipoindg) then ! Exclude self-edges
              dual_graph%ia(ipoindg+1) = dual_graph%ia(ipoindg+1) + 1
          end if
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors(ineigh) = 0
       end do
    end do

  end subroutine count_dual_graph_part_new

  !============================================================================================
  subroutine  list_dual_graph_part_new ( primal_mesh, dual_mesh, dual_graph, ws_position,  &
       &                                 ws_neighbors )
    implicit none
    ! Parameters
    type(fem_mesh) , intent(in)    :: primal_mesh, dual_mesh
    type(fem_graph), intent(inout) :: dual_graph
    integer(ip), intent(out)       :: ws_position  (dual_mesh%npoin)
    integer(ip), intent(out)       :: ws_neighbors (primal_mesh%nnode*dual_mesh%nnode)

    ! Local variables
    integer(ip)                     :: ineigh  
    integer(ip)                     :: ipoinpm ! Primal mesh point identifier
    integer(ip)                     :: ipoindm ! Dual mesh point identifier
    integer(ip)                     :: ipoindg ! Dual graph point identifier 
    integer(ip)                     :: p_ipoindm   
    integer(ip)                     :: p_ipoinpm   
    integer(ip)                     :: inods1,inods2
    integer(ip)                     :: inods1d, inods2d  
  

    integer(ip)                     :: first_free_ja_pos
    integer(ip)                     :: first_free_pos
    integer(ip)                     :: num_neighbors

    ! Initialize work space of filled-up 
    ! positions to zeros
    ws_position   = 0

    ! Initialize working space of neighbour 
    ! identifiers 
    ws_neighbors = 0

    first_free_ja_pos = 1

    do ipoindg=1, dual_mesh%npoin       
       first_free_pos = 1
       num_neighbors  = 0

       inods1=primal_mesh%p_face(ipoindg)
       inods2=primal_mesh%pnods(ipoindg+1)-1
       do p_ipoinpm = inods1,inods2
          ipoinpm = primal_mesh%lnods(p_ipoinpm)

          inods1d=dual_mesh%pnods(ipoinpm)
          inods2d=dual_mesh%pnods(ipoinpm+1)-1
          do p_ipoindm = inods1d, inods2d
             ipoindm = dual_mesh%lnods(p_ipoindm)

             ! If ipoinpm not visited yet
             if ( ws_position(ipoindm) == 0 ) then
                ws_position(ipoindm) = first_free_pos
                ws_neighbors(first_free_pos) = ipoindm
                first_free_pos = first_free_pos + 1
                num_neighbors = num_neighbors + 1 
             end if
          end do
       end do

       ! Extract those neighbours which are listed at least min_freq_neig times
       ! while restoring working space to initial state
       do ineigh=1, num_neighbors
          if ( .not.ws_neighbors(ineigh)==ipoindg ) then ! Exclude self-edges
             dual_graph%ja(first_free_ja_pos) = ws_neighbors(ineigh)
             first_free_ja_pos = first_free_ja_pos + 1
          end if
          ws_position ( ws_neighbors(ineigh) ) = 0
          ws_neighbors( ineigh               ) = 0
       end do
    end do

  end subroutine list_dual_graph_part_new


  ! Inspired on http://en.wikipedia.org/wiki/Breadth-first_search.
  ! Given a mesh (m) and its dual graph (g), it computes the list 
  ! of nodes (lconn) of each connected component in m
  subroutine mesh_graph_compute_connected_components (m, g, lconn)
    implicit none

    ! Parameters
    type(fem_mesh) , intent(in)   :: m   
    type(fem_graph), intent(in)   :: g
    type(list)     , intent(out)  :: lconn
 
    ! Locals
    integer(ip), allocatable :: auxv(:), auxe(:), e(:)
    integer(ip), allocatable :: emarked(:), vmarked(:)
    integer(ip), allocatable :: q(:)
    integer(ip)              :: head, tail, i, esize, vsize, current, & 
                                j, l, k, inods1d, inods2d, p_ipoin, ipoin

    call memalloc ( g%nv   , auxe     , __FILE__,__LINE__)
    call memalloc ( g%nv   , auxv     , __FILE__,__LINE__)
    call memalloc ( g%nv   , q        , __FILE__,__LINE__)
    call memalloc ( g%nv   , emarked  , __FILE__,__LINE__)
    call memalloc ( m%npoin, vmarked  , __FILE__,__LINE__)
    call memalloc ( g%nv  ,  e        , __FILE__,__LINE__)

    lconn%n  = 0
    emarked  = 0
    current  = 1 

    do i=1, g%nv
       if (emarked(i) == 0) then
          ! New connected component
          lconn%n = lconn%n +1
          esize   = 0
          vsize   = 0
          vmarked = 0 
!!$1  procedure BFS(G,v):
!!$2      create a queue Q
          head=1
          tail=1
!!$3      enqueue v onto Q
          q(tail)=i
          tail=tail+1
!!$4      mark v
          emarked(i)=1
          e(current)=i
          esize  = esize + 1
          current = current + 1  

!!$5      while Q is not empty:
          do while (head/=tail)
!!$6         t ← Q.dequeue()
             j=q(head)
             head = head + 1

             ! Traverse the nodes of the element number j
             if( m%nelty == 1 ) then
                inods1d = (j-1) * m%nnode+1
                inods2d = j     * m%nnode
             else
                inods1d = m%pnods(j)
                inods2d = m%pnods(j+1)-1
             end if

             do p_ipoin = inods1d, inods2d
                ipoin = m%lnods(p_ipoin)
                if (vmarked(ipoin)==0) then
                   vmarked(ipoin)=1
                   vsize = vsize+1
                end if
             end do

!!$9         for all edges e in G.adjacentEdges(t) do
             do k=g%ia(j), g%ia(j+1)-1
!!$12           u ← G.adjacentVertex(t,e)
                l=g%ja(k)
!!$13           if u is not emarked:
                if (emarked(l)==0) then
!!$14              mark u
                   emarked(l)=1
                   e(current)=l
                   esize  = esize + 1
                   current = current + 1  

!!$15              enqueue u onto Q
                   q(tail)=l
                   tail=tail+1
                end if
             end do
          end do
          auxe(lconn%n) = esize
          auxv(lconn%n) = vsize
       end if
    end do
    
    call memalloc( lconn%n+1, lconn%p, __FILE__,__LINE__)

    lconn%p(1) = 1
    do i=1, lconn%n
       lconn%p(i+1) = lconn%p(i) + auxv(i)
    end do

    call memfree( auxv   ,__FILE__,__LINE__)
    call memfree( q      ,__FILE__,__LINE__)
    call memfree( emarked,__FILE__,__LINE__)

    call memalloc( lconn%p(lconn%n+1)-1, lconn%l, __FILE__,__LINE__)

    current=1
    l=1
    do i=1, lconn%n
       vmarked = 0
       ! Traverse elements of current connected component  
       do current=current,current+auxe(i)-1
          j=e(current)

          ! Traverse the nodes of the element number j
          if( m%nelty == 1 ) then
             inods1d = (j-1) * m%nnode+1
             inods2d = j     * m%nnode
          else
             inods1d = m%pnods(j)
             inods2d = m%pnods(j+1)-1
          end if

          do p_ipoin = inods1d, inods2d
             ipoin = m%lnods(p_ipoin)
             if (vmarked(ipoin)==0) then
                vmarked(ipoin)=1
                lconn%l(l)=ipoin
                l=l+1
             end if
          end do
          
       end do
    end do

    ! write(*,*) 'ZZ', g%nv, m%npoin, l
    ! write(*,*) 'XX', lconn%n
    ! write(*,*) 'YY', lconn%p
    ! write(*,*) 'PP', lconn%l

    call memfree( auxe,__FILE__,__LINE__)
    call memfree( e,__FILE__,__LINE__)
    call memfree( vmarked,__FILE__,__LINE__)

  end subroutine mesh_graph_compute_connected_components

end module mesh_graph
