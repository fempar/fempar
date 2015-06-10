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
module dof_distribution_names
  ! Serial modules
  use types
  use memor
  use maps_names
  use dof_import_names
    
  implicit none
# include "debug.i90"
  private
  
  type dof_distribution
     integer(ip)              :: ipart
     integer(ip)              :: nparts
     
     integer(ip)              :: nl ! Number of local DoFs (formely in f_part%nmap%nl)
     integer(ip)              :: ni ! Number of interior DoFs (formely in f_part%nmap%ni)
     integer(ip)              :: nb ! Number of interface DoFs (formely in f_part%nmap%nb)

     integer(ip)              :: npadj       ! Number of adjacent parts
     integer(ip), allocatable :: lpadj(:)    ! List of adjacent parts
     
     ! To add to the documentation, the contents of lobjs are as follows:
     ! lobjs(1,*) : pysical variable associated to the DoF
     ! lobjs(2,*) : first DoF of the object in local ordering
     ! lobjs(3,*) : last DoF of the object in local ordering
     ! lobjs(4,*) : number of parts around this object
     ! lobjs(5:,*) : list of parts around this object

     integer(ip)              :: max_nparts  ! Maximum number of parts around any DoF communication object
     integer(ip)              :: nobjs       ! Number of local DoF communication objects
     integer(ip), allocatable :: lobjs(:,:)  ! List of local DoF communication objects
     
     type(list)               :: int_objs    ! List of objects on each edge to an adjacent part / Interface_objects
     type(map)                :: omap        ! Objects local to global map

     type(dof_import)         :: dof_import  ! Object which contains the control data to drive DoF nearest neigbour exchanges
  end type dof_distribution

  ! Types
  public :: dof_distribution
  
  ! Functions
  public :: dof_distribution_compute_import, dof_distribution_free, dof_distribution_print
  
contains

  subroutine dof_distribution_free(dof_dist)
    implicit none
    type(dof_distribution), intent(inout) :: dof_dist
    
    call memfree ( dof_dist%lpadj, __FILE__,__LINE__  )  
    call memfree ( dof_dist%lobjs, __FILE__,__LINE__  )
    call memfree ( dof_dist%int_objs%p, __FILE__, __LINE__ )
    call memfree ( dof_dist%int_objs%l, __FILE__, __LINE__ )
    call map_free ( dof_dist%omap )
    call dof_import_free ( dof_dist%dof_import )

  end subroutine dof_distribution_free

  ! This module uses a naive algorithm to map
  ! shared DoFs on the interface to owner parts.
  ! It is based on the following two rules:
  ! 1. Balance shared DoFs among two parts.
  ! 2. Map arbitrarily shared DoFs among several
  !    parts to the part with minimum identifer.
  subroutine dof_distribution_compute_import (dof_dist)
    implicit none

    ! Parameters
    type(dof_distribution), intent(inout)  :: dof_dist
    
    ! Locals
    integer(ip), allocatable :: ws_parts_visited_rcv_pos(:), ws_parts_visited_snd_pos(:), &
         &  ws_parts_visited_rcv_iedge(:), ws_parts_visited_snd_iedge(:)

    
    dof_dist%dof_import%ipart  = dof_dist%ipart
    dof_dist%dof_import%nparts = dof_dist%nparts
    
    call memalloc ( dof_dist%nb, dof_dist%dof_import%owner, __FILE__,__LINE__)
    
    call compute_owner ( dof_dist%lobjs, dof_dist%max_nparts, dof_dist%nobjs, & 
                         dof_dist%dof_import%owner, dof_dist%nb )
      
    ! Allocate workspace 
    call memalloc ( dof_dist%npadj, ws_parts_visited_rcv_pos, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj, ws_parts_visited_snd_pos, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj, ws_parts_visited_rcv_iedge, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj, ws_parts_visited_snd_iedge, __FILE__, __LINE__)
    

    call memalloc ( dof_dist%npadj, dof_dist%dof_import%list_rcv, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj+1, dof_dist%dof_import%rcv_ptrs, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj, dof_dist%dof_import%list_snd, __FILE__, __LINE__)
    call memalloc ( dof_dist%npadj+1, dof_dist%dof_import%snd_ptrs, __FILE__, __LINE__)

    call compute_snd_rcv_control_data_count ( dof_dist%ipart, &  
         dof_dist%lpadj, dof_dist%npadj, &
         dof_dist%lobjs, dof_dist%max_nparts, dof_dist%nobjs, &
         dof_dist%int_objs%l, dof_dist%int_objs%p, dof_dist%int_objs%n, &
         dof_dist%dof_import%owner, dof_dist%nb, &
         dof_dist%dof_import%num_rcv, dof_dist%dof_import%list_rcv, dof_dist%dof_import%rcv_ptrs, &
         dof_dist%dof_import%num_snd, dof_dist%dof_import%list_snd, dof_dist%dof_import%snd_ptrs, &
         ws_parts_visited_rcv_pos, &
         ws_parts_visited_snd_pos, &
         ws_parts_visited_rcv_iedge, & 
         ws_parts_visited_snd_iedge )

    ! Realloc list/snd info to its actual size. This is
    ! required, e.g., by psb_snd and psb_rcv comm routines
    call memrealloc ( dof_dist%dof_import%num_rcv  , dof_dist%dof_import%list_rcv, __FILE__, __LINE__)
    call memrealloc ( dof_dist%dof_import%num_rcv+1, dof_dist%dof_import%rcv_ptrs, __FILE__, __LINE__)
    call memrealloc ( dof_dist%dof_import%num_snd  , dof_dist%dof_import%list_snd, __FILE__, __LINE__)
    call memrealloc ( dof_dist%dof_import%num_snd+1, dof_dist%dof_import%snd_ptrs, __FILE__, __LINE__)
    
    call memalloc ( dof_dist%dof_import%rcv_ptrs(dof_dist%dof_import%num_rcv+1)-dof_dist%dof_import%rcv_ptrs(1), & 
         dof_dist%dof_import%unpack_idx, __FILE__,__LINE__)
    
    call memalloc ( dof_dist%dof_import%snd_ptrs(dof_dist%dof_import%num_snd+1)-dof_dist%dof_import%snd_ptrs(1), & 
         dof_dist%dof_import%pack_idx, __FILE__,__LINE__)
    
    ! Free workspace 
    call memfree ( ws_parts_visited_rcv_pos,__FILE__,__LINE__)
    call memfree ( ws_parts_visited_snd_pos,__FILE__,__LINE__)
    
    call compute_snd_rcv_control_data_list  ( dof_dist%ipart, & 
         dof_dist%lpadj, dof_dist%npadj, &
         dof_dist%lobjs, dof_dist%max_nparts, dof_dist%nobjs, &
         dof_dist%int_objs%l, dof_dist%int_objs%p, dof_dist%int_objs%n, &
         dof_dist%dof_import%owner, dof_dist%nb, &
         dof_dist%dof_import%num_rcv, dof_dist%dof_import%list_rcv, dof_dist%dof_import%rcv_ptrs, dof_dist%dof_import%unpack_idx, &
         dof_dist%dof_import%num_snd, dof_dist%dof_import%list_snd, dof_dist%dof_import%snd_ptrs, dof_dist%dof_import%pack_idx, &
         ws_parts_visited_rcv_iedge, &
         ws_parts_visited_snd_iedge )

    ! Free workspace 
    call memfree ( ws_parts_visited_rcv_iedge,__FILE__,__LINE__)
    call memfree ( ws_parts_visited_snd_iedge,__FILE__,__LINE__)
  end subroutine dof_distribution_compute_import

  ! Auxiliar subroutine
  subroutine compute_owner ( lobjs, max_nparts, nobjs, & 
                             owner, nowner )
    implicit none 

    ! Parameters
    integer(ip), intent(in)  :: max_nparts, nobjs
    integer(ip), intent(in)  :: lobjs(max_nparts+4, nobjs)
    integer(ip), intent(in)  :: nowner
    integer(ip), intent(out) :: owner(nowner)

    ! Local variables
    integer(ip) :: iobj, ni, i, ipart_min, ipart_max, ipoin
    
    ! Determine number of interior nodes
    ni = ( lobjs(3,1)-lobjs(2,1) + 1 )

    ! Determine owner of each shared node
    ! Loop on part (boundary) objects
    do iobj=2,nobjs
       ! If object shared among two parts
       if ( lobjs(4, iobj) == 2 ) then
          ipart_min = min( lobjs(5, iobj), lobjs(6, iobj) )
          do ipoin=lobjs(2,iobj), lobjs(2,iobj) + ( lobjs(3,iobj)-lobjs(2,iobj) )/2
             owner( ipoin-ni ) = ipart_min
          end do
          
          ipart_max = max( lobjs(5, iobj), lobjs(6, iobj) )
          
!!$          write(*,*) 'ZZZ', ni 
!!$          write(*,*) 'XXX', lobjs(2,iobj), lobjs(2,iobj) + ( lobjs(3,iobj)-lobjs(2,iobj) )/2
!!$          write(*,*) 'YYY', lobjs(2,iobj) + ( lobjs(3,iobj)-lobjs(2,iobj) )/2 + 1, lobjs(3,iobj)
!!$          write(*,*) 'ipart_min', ipart_min, 'ipart_max', ipart_max, lobjs(:,iobj) ! DBG:
!!$          write(*,*) iobj
          
          do ipoin=lobjs(2,iobj) + ( lobjs(3,iobj)-lobjs(2,iobj) )/2 + 1, lobjs(3,iobj)
             owner( ipoin-ni ) = ipart_max
          end do
          
          ! If objects shared among more than two parts
       else if ( lobjs(4, iobj) > 2 ) then
          ipart_min = lobjs(5, iobj)
!!$          ipart_max = lobjs(5, iobj) ! DBG: 
          do i=2, lobjs(4, iobj)
             ipart_min = min( ipart_min, lobjs(4+i, iobj) )
!!$             ipart_max = max( ipart_max, lobjs(4+i, iobj) ) ! DBG: 
          end do
!!$          write (*,*) 'ipart_min', ipart_min, 'ipart_max', ipart_max
!!$          write(*,*) 'ipart_min', ipart_min, lobjs(:,iobj) ! DBG:
          
          do ipoin=lobjs(2,iobj), lobjs(3,iobj)
             owner( ipoin-ni ) = ipart_min
          end do
       end if
    end do

  end subroutine compute_owner

  ! Auxiliar subroutine
  subroutine compute_snd_rcv_control_data_count ( ipart,                              &
                                                  lpadj, npadj,                       & 
                                                  lobjs, max_nparts, nobjs,           &
                                                  int_objs_l, int_objs_p, n_int_objs, &
                                                  owner, nowner,                      &
                                                  num_rcv, list_rcv, rcv_ptrs,        &
                                                  num_snd, list_snd, snd_ptrs,        &
                                                  ws_parts_visited_rcv_pos,           &
                                                  ws_parts_visited_snd_pos,           &
                                                  ws_parts_visited_rcv_iedge,         &
                                                  ws_parts_visited_snd_iedge     )
    implicit none

    ! Parameters
    integer(ip), intent(in)  :: ipart
    integer(ip), intent(in)  :: npadj
    integer(ip), intent(in)  :: lpadj(npadj)
    integer(ip), intent(in)  :: max_nparts, nobjs
    integer(ip), intent(in)  :: lobjs(max_nparts+4, nobjs)
    integer(ip), intent(in)  :: n_int_objs
    integer(ip), intent(in)  :: int_objs_p(n_int_objs+1)
    integer(ip), intent(in)  :: int_objs_l(int_objs_p(n_int_objs+1)-1)
    integer(ip), intent(in)  :: owner( nowner )
    integer(ip), intent(in)  :: nowner

    ! Control info to receive
    integer, intent(out) :: num_rcv, list_rcv(npadj), rcv_ptrs(npadj+1)

    ! Control info to send
    integer, intent(out) :: num_snd, list_snd(npadj), snd_ptrs(npadj+1)

    ! Work space
    integer(ip), intent(inout) :: ws_parts_visited_rcv_pos (npadj)
    integer(ip), intent(inout) :: ws_parts_visited_snd_pos (npadj)
    integer(ip), intent(inout) :: ws_parts_visited_rcv_iedge (npadj)
    integer(ip), intent(inout) :: ws_parts_visited_snd_iedge (npadj)

   
    ! Locals
    integer(ip) :: iedge, neighbour_part_id, iposobj, ipoin, ni
    integer(ip) :: num_neighbours_snd, i, idobj

    ! Determine number of interior nodes
    ni = ( lobjs(3,1)-lobjs(2,1) + 1 )
    
    num_rcv = 0
    num_snd = 0 
    ws_parts_visited_snd_pos    =  0 
    ws_parts_visited_rcv_pos    =  0
    ws_parts_visited_rcv_iedge  = -1
    ws_parts_visited_snd_iedge  = -1
    

    ! Traverse ipart's neighbours in the graph of parts
    do iedge = 1, npadj
       neighbour_part_id = lpadj (iedge)

       ! write(*,*) 'XXX', ipart, neighbour_part_id, int_objs_p(iedge), int_objs_p(iedge+1) ! DBG:

       ! Traverse list of objects on the current edge
       do iposobj=int_objs_p(iedge), int_objs_p(iedge+1)-1
          idobj = int_objs_l (iposobj)
          ! write(*,*) 'XXX', ipart, neighbour_part_id, idobj ! DBG:

          ! if ( ipart == 13 .and. neighbour_part_id == 3 ) then
          !  write(*,*) 'XXX', lobjs(3,idobj)-lobjs(2,idobj)
          ! else if ( ipart == 3 .and. neighbour_part_id == 13) then
          !  write(*,*) 'YYY', lobjs(3,idobj)-lobjs(2,idobj)
          ! end if
          ! Traverse mesh points of the current object
          do ipoin=lobjs(2,idobj), lobjs(3,idobj)
             !if(ipart==1) then
             !   write(*,*) ipart, neighbour_part_id, iedge, iposobj, idobj, (ipoin-ni)
             !end if
             if ( owner(ipoin-ni) == ipart ) then
                ! Did I already visit neighbour_part_id ? 
                if ( ws_parts_visited_rcv_pos( iedge ) == 0 ) then
                   ! No 
                   num_rcv = num_rcv + 1 
                   list_rcv ( num_rcv ) = neighbour_part_id
                   ws_parts_visited_rcv_pos  (iedge)    = num_rcv
                   ws_parts_visited_rcv_iedge (num_rcv) = iedge
                   rcv_ptrs (num_rcv+1) = 1
                else
                   ! Yes
                   rcv_ptrs (ws_parts_visited_rcv_pos(iedge)+1) = rcv_ptrs(ws_parts_visited_rcv_pos(iedge)+1) + 1
                end if
             else

                if ( owner(ipoin-ni) == neighbour_part_id ) then
                    ! Did I already visit neighbour_part_id ? 
                    if ( ws_parts_visited_snd_pos( iedge ) == 0 ) then
                       ! No 
                       num_snd = num_snd + 1 
                       list_snd ( num_snd ) = neighbour_part_id
                       ws_parts_visited_snd_pos  (iedge)    = num_snd
                       ws_parts_visited_snd_iedge (num_snd) = iedge
                       snd_ptrs (num_snd+1) = 1
                    else
                       ! Yes
                       snd_ptrs (ws_parts_visited_snd_pos(iedge)+1) = snd_ptrs(ws_parts_visited_snd_pos(iedge)+1) + 1
                    end if
                end if
             end if
          end do ! ipoin
       end do   ! iposobj
    end do     ! iedge

    ! write(*,*) 'KKK', num_rcv, num_snd ! DBG:


    ! Transform rcv_ptrs from size to pointers
    rcv_ptrs(1) = 1
    do i = 1, num_rcv
       rcv_ptrs (i+1) = rcv_ptrs (i+1) + rcv_ptrs (i)  
    end do

    ! Transform snd_ptrs from size to pointers
    snd_ptrs(1) = 1
    do i = 1, num_snd
       snd_ptrs (i+1) = snd_ptrs(i+1) + snd_ptrs (i)  
    end do

  end subroutine compute_snd_rcv_control_data_count

  ! Auxiliar routine
  subroutine compute_snd_rcv_control_data_list ( ipart,                                   & 
                                                 lpadj, npadj,                            & 
                                                 lobjs, max_nparts, nobjs,                &
                                                 int_objs_l, int_objs_p, n_int_objs,      &
                                                 owner, nowner,                           &
                                                 num_rcv, list_rcv, rcv_ptrs, unpack_idx, &
                                                 num_snd, list_snd, snd_ptrs, pack_idx  , &
                                                 ws_parts_visited_rcv_edge_id,            &
                                                 ws_parts_visited_snd_edge_id             )
    use types
    implicit none

    ! Parameters
    integer(ip), intent(in)  :: ipart
    integer(ip), intent(in)  :: npadj  
    integer(ip), intent(in)  :: lpadj(npadj)
    integer(ip), intent(in)  :: max_nparts, nobjs
    integer(ip), intent(in)  :: lobjs(max_nparts+4, nobjs)
    integer(ip), intent(in)  :: n_int_objs
    integer(ip), intent(in)  :: int_objs_p(n_int_objs+1)
    integer(ip), intent(in)  :: int_objs_l(int_objs_p(n_int_objs+1)-1)
    integer(ip), intent(in)  :: nowner 
    integer(ip), intent(in)  :: owner( nowner )

    ! Control info to receive
    integer    , intent(in)  :: num_rcv, list_rcv(num_rcv), rcv_ptrs(num_rcv+1)
    integer(ip), intent(out) :: unpack_idx ( rcv_ptrs(num_rcv+1)-rcv_ptrs(1)  )

    ! Control info to send
    integer    , intent(in)  :: num_snd, list_snd(num_snd), snd_ptrs(num_snd+1)
    integer(ip), intent(out) :: pack_idx (snd_ptrs(num_snd+1)-snd_ptrs(1))

    ! Work-space
    integer(ip), intent(in) :: ws_parts_visited_rcv_edge_id (num_rcv)
    integer(ip), intent(in) :: ws_parts_visited_snd_edge_id (num_snd)

    ! Locals
    integer(ip) :: ipos_rcv, ipos_snd, iposobj, ipoin, ni
    integer(ip) :: i, idobj, neighbour_part_id, cur_pos, iedge

    ! Determine number of interior nodes
    ni = ( lobjs(3,1)-lobjs(2,1) + 1 )    


!!$    write(*,'(a)') 'List of parts I have to receive from:' 
!!$    write(*,'(10i10)') list_rcv(1:num_rcv)
!!$
!!$    write(*,'(a)') 'Rcv_ptrs:'
!!$    write(*,'(10i10)') rcv_ptrs(1:num_rcv+1)
!!$
!!$    write(*,'(a,i10)') 'Number of parts I have send to:', &
!!$          &  num_snd
!!$
!!$    write(*,'(a)') 'List of parts I have to send to:'
!!$    write(*,'(10i10)') list_snd(1:num_snd)
!!$
!!$    write(*,'(a)') 'Snd_ptrs::'
!!$    write(*,'(10i10)') snd_ptrs(1:num_snd+1)



    do ipos_rcv = 1, num_rcv
       iedge   = ws_parts_visited_rcv_edge_id ( ipos_rcv )
       cur_pos = rcv_ptrs (ipos_rcv)
       ! write(*,*) 'KKKR', iedge, cur_pos, lpadj(iedge) ! DBG:
       ! Traverse list of objects on the current edge
       do iposobj=int_objs_p(iedge), int_objs_p(iedge+1)-1
          idobj = int_objs_l (iposobj)
          ! write(*,*) 'XXXRA', iedge, cur_pos, lpadj(iedge), idobj ! DBG:
          ! Traverse mesh points of the current object
          do ipoin=lobjs(2,idobj), lobjs(3,idobj)
             if ( owner(ipoin-ni) == ipart ) then
                unpack_idx (cur_pos) = ipoin - ni 
                cur_pos = cur_pos + 1          
             end if
          end do ! ipoin
          ! write(*,*) 'XXXRD', iedge, cur_pos, lpadj(iedge), idobj ! DBG:
       end do   ! iposobj
    end do     ! ipos_rcv

    do ipos_snd = 1, num_snd
       iedge   = ws_parts_visited_snd_edge_id ( ipos_snd )
       neighbour_part_id = lpadj (iedge) 
       cur_pos = snd_ptrs (ipos_snd)
       ! write(*,*) 'KKKS', iedge, cur_pos, lpadj(iedge) ! DBG:  
       ! Traverse list of objects on the current edge
       do iposobj=int_objs_p(iedge), int_objs_p(iedge+1)-1
          idobj = int_objs_l (iposobj)
          ! write(*,*) 'XXXSA', iedge, cur_pos, lpadj(iedge), idobj ! DBG:

          ! Traverse mesh points of the current object
          do ipoin=lobjs(2,idobj), lobjs(3,idobj)
             if ( owner(ipoin-ni) == neighbour_part_id  ) then
                pack_idx (cur_pos) = ipoin - ni
                cur_pos = cur_pos + 1          
             end if
          end do ! ipoin
          ! write(*,*) 'XXXSD', iedge, cur_pos, lpadj(iedge), idobj ! DBG:
       end do   ! iposobj
    end do     ! ipos_snd
    ! write(*,*) 'END'

  end subroutine compute_snd_rcv_control_data_list


  !=============================================================================
  subroutine dof_distribution_print (lu_out, dof_dist)
    !-----------------------------------------------------------------------
    ! This routine prints a dof_distribution object
    !-----------------------------------------------------------------------
    use types
    implicit none


    ! Parameters
    type(dof_distribution), intent(in)  :: dof_dist
    integer(ip)        , intent(in)  :: lu_out

    ! Local variables
    integer (ip) :: i, j

    if(lu_out>0) then

       write(lu_out,'(a)') '*** begin dof_distribution data structure ***'

       write(lu_out,'(a,i10)') 'Number of parts:', &
          &  dof_dist%nparts

       write(lu_out,'(a,i10)') 'Number of adjacent parts:', &
          &  dof_dist%npadj
       
       write(lu_out,'(10i10)') dof_dist%lpadj(1:dof_dist%npadj)

       write(lu_out,'(a,i10)') 'Number of local objects:', &
          &  dof_dist%nobjs

       write(lu_out,'(a)') 'List of local objects:'
       do i=1,dof_dist%omap%nl
          write(lu_out,'(10i10)') dof_dist%omap%l2g(i),dof_dist%lobjs(:,i)
       end do

       write(lu_out,'(a)') 'List of interface objects:'
       do i=1,dof_dist%npadj
          write(lu_out,'(10i10)') i, &
             & (dof_dist%int_objs%l(j),j=dof_dist%int_objs%p(i),dof_dist%int_objs%p(i+1)-1)
       end do

       write(lu_out,'(a)') '*** end dof_distribution data structure ***'

    end if
 
  end subroutine dof_distribution_print

  
end module dof_distribution_names
