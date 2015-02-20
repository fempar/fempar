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
module graph_renum
  use types
  use memor
  use fem_graph_class
  use renum_class
  use metis_interface
  use fem_mesh_partition_base
  use rcm_renum
#include "debug.i90"
  implicit none
  private

  ! Functions
  public :: graph_nd_renumbering, graph_pt_renumbering

contains
  !=================================================================================================
  ! Methods defined for fem_graph are:
  !
  ! - graph_nd_renumbering
  ! - graph_pt_renumbering
  !
  ! These methods call auxiliar (private) routines
  !
  ! - parts_renum
  ! - metis_nodendextractseparatortree
  ! - metis_partgraphrecursive
  !
  ! Routines without interface
  !
  ! - subroutine permut_graph
  ! - subroutine rcm_renum_masked
  ! - subroutine rcm_renum
  !
  ! TO DO:
  !    * interfaces
  !
  !=================================================================================================
  subroutine graph_nd_renumbering(prt_parts,gp,nleve,ren,ldomn)
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    type(part_params), intent(in) :: prt_parts
    type(fem_graph), target, intent(in) :: gp
    integer(ip), target,  intent(in)     :: nleve
    type(renum), target,  intent(inout)  :: ren
    integer(ip), target,  intent(out)    :: ldomn(gp%nv)

    assert(ren%n==gp%nv)
    assert(nleve>0)
    
    if ( gp%nv == 1 ) then
       ren%lperm(1) = 1
       ren%iperm(1) = 1
       ldomn(1) = 1
    else
#ifdef ENABLE_METIS5
       ierr = fp_metis_setdefaultoptions(c_loc(options))
       assert(ierr == METIS_OK) 
       
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ierr = fp_metis_nodendextractseparatortree( c_loc(gp%nv),c_loc(gp%ia),c_loc(gp%ja),C_NULL_PTR,c_loc(options), &
            &                                      c_loc(ren%iperm),c_loc(ren%lperm),c_loc(nleve),c_loc(ldomn))
       
       assert(ierr == METIS_OK)
#else
#ifdef ENABLE_METIS4
       optmt_nd(5) = prt_parts%metis_option_debug
       call fp_metis_nodendextractseparatortree(gp%nv,gp%ia,gp%ja,1,optmt_nd, &
            &                                   ren%iperm,ren%lperm,nleve,ldomn)
#else
       call enable_metis_error_message
#endif
#endif
    end if

  end subroutine graph_nd_renumbering

  !=================================================================================================
  subroutine graph_pt_renumbering(prt_parts,gp,ren,ldomn)
    !-----------------------------------------------------------------------
    ! This routine computes a nparts-way-partitioning of graph
    !-----------------------------------------------------------------------
    implicit none
    type(part_params), target, intent(in)    :: prt_parts
    type(fem_graph)  , target, intent(inout) :: gp
    type(renum)      , target, intent(inout)   :: ren
    integer(ip)      , target, intent(out)     :: ldomn(gp%nv)

    ! Local variables 
    integer(ip), target      :: kedge
    integer(ip)              :: idumm,iv
    integer(ip), allocatable :: lwork(:)
    integer(ip)              :: i, j, m, k, ipart
    integer(ip), allocatable :: iperm(:)

   
#ifdef ENABLE_METIS5
    ierr = fp_metis_setdefaultoptions(c_loc(options))
    assert(ierr == METIS_OK) 

!!$      From METIS 5.0 manual:
!!$
!!$      The following options are valid for METIS PartGraphRecursive:
!!$      
!!$      METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE, METIS_OPTION_RTYPE,
!!$      METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS, METIS_OPTION_NITER,
!!$      METIS_OPTION_SEED, METIS_OPTION_UFACTOR, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL
!!$     
!!$      The following options are valid for METIS PartGraphKway:
!!$ 
!!$      METIS_OPTION_OBJTYPE, METIS_OPTION_CTYPE, METIS_OPTION_IPTYPE,
!!$      METIS_OPTION_RTYPE, METIS_OPTION_NO2HOP, METIS_OPTION_NCUTS,
!!$      METIS_OPTION_NITER, METIS_OPTION_UFACTOR, METIS_OPTION_MINCONN,
!!$      METIS_OPTION_CONTIG, METIS_OPTION_SEED, METIS_OPTION_NUMBERING,
!!$      METIS_OPTION_DBGLVL

    if ( prt_parts%strat == part_kway ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       
       ! Enforce contiguous partititions
       options(METIS_OPTION_CONTIG)    = prt_parts%metis_option_contig
       
       ! Explicitly minimize the maximum degree of the subdomain graph
       options(METIS_OPTION_MINCONN)   = prt_parts%metis_option_minconn

       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor
       
       ncon = 1 
       ierr = fp_metis_partgraphkway( c_loc(gp%nv), c_loc(ncon), c_loc(gp%ia)  , c_loc(gp%ja) , & 
            C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
            C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
       
       assert(ierr == METIS_OK) 
       
    else if ( prt_parts%strat == part_recursive ) then
       options(METIS_OPTION_NUMBERING) = 1
       options(METIS_OPTION_DBGLVL)    = prt_parts%metis_option_debug
       options(METIS_OPTION_UFACTOR)   = prt_parts%metis_option_ufactor

       ncon = 1 
       ierr = fp_metis_partgraphrecursive( c_loc(gp%nv), c_loc(ncon), c_loc(gp%ia)  , c_loc(gp%ja) , & 
            C_NULL_PTR  , C_NULL_PTR , C_NULL_PTR    , c_loc(prt_parts%nparts), &
            C_NULL_PTR  , C_NULL_PTR , c_loc(options), c_loc(kedge), c_loc(ldomn) )
    end if
    
    ! write(*,*) 'Cut edges:', kedge
#else
#ifdef ENABLE_METIS4
    optmt_pr(5) = prt_parts%metis_option_debug
    if ( prt_parts%strat == part_kway ) then
       call fp_metis_partgraphkway(gp%nv,gp%ia,gp%ja,idumm,idumm,0,1,prt_parts%nparts,optmt_pr, &
            &                           kedge,ldomn)
    else if ( prt_parts%strat == part_recursive ) then
       call fp_metis_partgraphrecursive(gp%nv,gp%ia,gp%ja,idumm,idumm,0,1,prt_parts%nparts,optmt_pr, &
            &                          kedge,ldomn)
    end if
#else
    call enable_metis_error_message
#endif
#endif

    if ( prt_parts%strat == part_strip ) then
       j = gp%nv
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(m+i) = ipart
          end do
          m = m + k
          j = j - k
       end do
    else if ( prt_parts%strat == part_rcm_strip ) then
       call memalloc ( gp%nv, iperm, __FILE__,__LINE__ )
       call genrcm ( gp%nv, gp%ia(gp%nv+1)-1, gp%ia, gp%ja, iperm )
       j = gp%nv
       m = 0
       do ipart=1,prt_parts%nparts
          k = j / (prt_parts%nparts-ipart+1)
          do i = 1, k
             ldomn(iperm(m+i)) = ipart
          end do
          m = m + k
          j = j - k
       end do
       call memfree ( iperm,__FILE__,__LINE__)
    end if

    if(ren%n>0) call renum_by_sets(prt_parts%nparts,ldomn,ren)

  end subroutine graph_pt_renumbering

  !!================================================================================================
  !subroutine permut_graph(np,np2,nz,ia,ja,is,lrord,lrord2)
  !  !-----------------------------------------------------------------------
  !  ! This routine 
  !  !-----------------------------------------------------------------------
  !  implicit none
  !  integer(ip), intent(in)    :: np,np2,nz
  !  integer(ip), intent(in)    :: lrord(np),lrord2(np2)
  !  integer(ip), intent(inout) :: ia(np+1),is(np+1),ja(nz)
  !  integer(ip), allocatable   :: iaw(:),jaw(:)
  !  integer(ip)                :: i,j,k,iz,istat

  !  ! Reorder ia
  !  allocate(iaw(np+1),stat=istat)
  !  do i=1,np
  !     iaw(lrord(i)+1)=ia(i+1)-ia(i)
  !  end do
  !  iaw(1)=1
  !  do i=1,np
  !     iaw(i+1)=iaw(i+1)+iaw(i)
  !  end do

  !  ! Save ja
  !  allocate(jaw(nz),stat=istat)
  !  do i=1,np
  !     jaw(iaw(lrord(i)):iaw(lrord(i)+1)-1)=lrord2(ja(ia(i):ia(i+1)-1))
  !  end do

  !  ! Copy iaw and deallocate
  !  ia=iaw
  !  deallocate(iaw)

  !  ! Reorder ja and rebuild is
  !  is=0
  !  do i=1,np
  !     k=0
  !     do iz=ia(i),ia(i+1)-1
  !        if(jaw(iz)<=i) then
  !           ja(ia(i)+k)=jaw(iz)
  !           k=k+1
  !        else
  !           is(i+1)=is(i+1)+1
  !           ja(ia(i+1)-is(i+1))=jaw(iz)
  !        end if
  !     end do
  !  end do
  !  deallocate(jaw)

  !  ! Compress is
  !  is(1)=1
  !  do i=1,np
  !     is(i+1)=is(i+1)+is(i)
  !  end do

  !end subroutine permut_graph

!!$  !=================================================================================================
!!$  subroutine rcm_renum_masked (np,nz,ia,ja,lmask,imask,lr,nlv,lv)
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    ! This routine renumbers a graph using RCM algorithm.
!!$    !
!!$    ! np = number of vertices
!!$    ! nz = size of the graph
!!$    ! ia = pointers to adjacencies
!!$    ! ja = list of adjacencies
!!$    ! lr = On input first level nodes (pseudo-peripheral nodes) 
!!$    !      On output permutation vector
!!$    !
!!$    ! In this routine lr(new)=old
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)    :: np,nz,imask
!!$    integer(ip), intent(in)    :: lmask(np)
!!$    integer(ip), intent(inout) :: nlv
!!$    integer(ip), intent(inout) :: ia(np+1),ja(nz),lr(np),lv(np)
!!$    integer(ip)                :: i,j,ingb,iz,izbeg,izend,lvbeg,lvend,lsize,npc
!!$
!!$    ! First level already defined by the user
!!$    nlv=1
!!$    lvbeg=lv(nlv)
!!$    lvend=lv(nlv+1)-1
!!$    npc=lvend              ! Number of nodes already considered
!!$    lsize=lvend-lvbeg+1    ! Level size
!!$
!!$    do while(lsize>0)
!!$
!!$       nlv=nlv+1
!!$
!!$       ! Loop on nodes or current level
!!$       do j = lvbeg,lvend
!!$          i=lr(j)
!!$
!!$          ! Loop on the neighbors
!!$          izbeg = -ia(i)
!!$          izend =  abs(ia(i+1))-1
!!$
!!$          do iz=izbeg,izend
!!$             ingb = ja(iz)
!!$
!!$             ! Only on untouched and masked nodes
!!$             if( lmask(ingb)==imask .and. ia(ingb) > 0) then
!!$                npc=npc+1
!!$                lr(npc) = ingb
!!$                ia(ingb) = - ia(ingb) ! Touch it.
!!$             end if
!!$
!!$          end do
!!$
!!$       end do
!!$
!!$       ! Define next level
!!$       lvbeg=lvend+1
!!$       lvend=npc
!!$       lsize=lvend-lvbeg+1
!!$
!!$       lv(nlv+1) = npc+1
!!$
!!$    end do
!!$
!!$  end subroutine rcm_renum_masked
!!$
!!$  !=================================================================================================
!!$  subroutine rcm_renum (np,nz,ia,ja,lr,nlv,lv)
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    ! This routine renumber mesh graph using RCM algorithm.
!!$    !
!!$    ! np = number of vertices
!!$    ! nz = size of the graph
!!$    ! ia = pointers to adjacencies
!!$    ! ja = list of adjacencies
!!$    ! lr = On input first level nodes (pseudo-peripheral nodes) 
!!$    !      On output permutation vector
!!$    !
!!$    ! In this routine lr(new)=old
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    implicit none
!!$    integer(ip), intent(in)    :: np,nz
!!$    integer(ip), intent(inout) :: nlv
!!$    integer(ip), intent(inout) :: ia(np+1),ja(nz),lr(np),lv(np)
!!$    integer(ip)                :: i,j,ingb,iz,izbeg,izend,lvbeg,lvend,lsize,npc
!!$
!!$    ! First level defined by permutation array lr
!!$    lvbeg=1              ! Begining of the level
!!$    i=lvbeg
!!$    do while(lr(i)/=0)
!!$       ia(lr(i)) = - ia(lr(i))
!!$       lv(i)=1
!!$       i=i+1
!!$    end do
!!$    lvend=i-1            ! End of the level
!!$    npc=lvend            ! Number of nodes already considered
!!$    lsize=lvend-lvbeg    ! Level size
!!$    nlv=1                ! Number of levels found
!!$
!!$    lv(1) = 1            ! Pointers to the begining of the level
!!$    lv(nlv+1) = npc+1
!!$
!!$    do while(lsize>0)
!!$
!!$       nlv=nlv+1
!!$
!!$       ! Loop on nodes or current level
!!$       do j = lvbeg,lvend
!!$          i=lr(j)
!!$
!!$          ! Loop on the neighbors
!!$          izbeg = -ia(i)
!!$          izend =  abs(ia(i+1))-1
!!$
!!$          do iz=izbeg,izend
!!$             ingb = ja(iz)
!!$
!!$             ! Only on untouched nodes.
!!$             if( ia(ingb) > 0) then
!!$                npc=npc+1
!!$                lr(npc) = ingb
!!$                ia(ingb) = - ia(ingb) ! Touch it.
!!$             end if
!!$
!!$          end do
!!$
!!$       end do
!!$
!!$       ! Define next level
!!$       lvbeg=lvend+1
!!$       lvend=npc
!!$       lsize=lvend-lvbeg
!!$
!!$       lv(nlv+1) = npc+1
!!$
!!$    end do
!!$
!!$    if(npc/=np) then
!!$       write(*,*) 'More than one conected component!!'
!!$       stop
!!$    endif
!!$
!!$    do i=1,np
!!$       ia(i)=-ia(i)
!!$    end do
!!$
!!$    ! Reverse numbering
!!$    !do i = 1, np/2
!!$    !   j=lr(i)
!!$    !   lr(i)=lr(np+1-i)
!!$    !   lr(np+1-i)=j
!!$    !end do
!!$
!!$    return
!!$
!!$  end subroutine rcm_renum


end module graph_renum
