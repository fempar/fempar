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
module map_apply_names
  use types_names
  use memor_names
  use map_names
  use renumbering_names
  use mesh_names
  use hash_table_names
# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: l2g_add=0
  integer(ip), parameter :: l2g_copy=1

  interface map_apply_l2g
     module procedure map_apply_l2g_r1, map_apply_l2g_r2, map_apply_l2g_i1
  end interface map_apply_l2g

  interface map_apply_g2l
       module procedure map_apply_g2l_r1, map_apply_g2l_r2, map_apply_g2l_i1, map_apply_g2l_i2, &
            map_igp_apply_g2l_r1, map_igp_apply_g2l_r2, map_igp_apply_g2l_i1, map_igp_apply_g2l_i2
  end interface map_apply_g2l

  interface mesh_g2l
     module procedure mesh_g2l_emap_ip, mesh_g2l_emap_igp, mesh_g2l_nmap_igp_emap_igp
  end interface

  ! Constants
  public :: l2g_add, l2g_copy

  ! Functions
  public :: map_apply_l2g, map_apply_g2l, mesh_g2l, mesh_l2l

contains

  !=============================================================================
  subroutine map_apply_l2g_r1(ope,mp,xl,xg,renumbering)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map_t)  , intent(in)  :: mp
    real(rp)   , intent(in)  :: xl(mp%nl)
    real(rp)   , intent(out) :: xg(mp%ng)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(renumbering%iperm(mp%l2g(i)))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(renumbering%iperm(mp%l2g(i)))=xg(renumbering%iperm(mp%l2g(i)))+xl(i)
          end do
       end if
    else
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(mp%l2g(i))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(mp%l2g(i))=xg(mp%l2g(i))+xl(i)
          end do
       end if
    end if

  end subroutine map_apply_l2g_r1

  !================================================================================================
  subroutine map_apply_l2g_i1(ope,mp,xl,xg,renumbering)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: xl(mp%nl)
    integer(ip), intent(out) :: xg(mp%ng)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(renumbering%iperm(mp%l2g(i)))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(renumbering%iperm(mp%l2g(i)))=xg(renumbering%iperm(mp%l2g(i)))+xl(i)
          end do
       end if
    else
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(mp%l2g(i))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(mp%l2g(i))=xg(mp%l2g(i))+xl(i)
          end do
       end if
    end if

  end subroutine map_apply_l2g_i1

  !================================================================================================
  subroutine map_apply_l2g_r2(ope,mp,ld,xl,xg,renumbering)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    real(rp)   , intent(in)  :: xl(ld,mp%nl)
    real(rp)   , intent(out) :: xg(ld,mp%ng)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip)                   :: i

    if(present(renumbering)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(:,renumbering%iperm(mp%l2g(i)))=xl(:,i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(:,renumbering%iperm(mp%l2g(i)))=xg(:,renumbering%iperm(mp%l2g(i)))+xl(:,i)
          end do
       end if
    else
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(:,mp%l2g(i))=xl(:,i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(:,mp%l2g(i))=xg(:,mp%l2g(i))+xl(:,i)
          end do
       end if
    end if

  end subroutine map_apply_l2g_r2

  !================================================================================================
  subroutine map_apply_g2l_r1(mp,xg,xl,renumbering)
    implicit none
    type(map_t)  , intent(in)  :: mp
    real(rp)   , intent(in)  :: xg(mp%ng)
    real(rp)   , intent(out) :: xl(mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(i)=xg(renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_r1

  !================================================================================================
  subroutine map_apply_g2l_r2(mp,ld,xg,xl,renumbering)
    implicit none
    type(map_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    real(rp)   , intent(in)  :: xg(ld,mp%ng)
    real(rp)   , intent(out) :: xl(ld,mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if
  end subroutine map_apply_g2l_r2

  !================================================================================================
  subroutine map_apply_g2l_i1(mp,xg,xl,renumbering)
    implicit none
    type(map_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: xg(mp%ng)
    integer(ip), intent(out) :: xl(mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(i)=xg(renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_i1

  !================================================================================================
  subroutine map_apply_g2l_i2(mp,ld,xg,xl,renumbering)
    implicit none
    type(map_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    integer(ip), intent(in)  :: xg(ld,mp%ng)
    integer(ip), intent(out) :: xl(ld,mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_i2

   !================================================================================================
  subroutine map_igp_apply_g2l_r1(mp,xg,xl,renumbering)
    implicit none
    type(map_igp_t)  , intent(in)  :: mp
    real(rp)   , intent(in)  :: xg(mp%ng)
    real(rp)   , intent(out) :: xl(mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(i)=xg(renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_igp_apply_g2l_r1

  !================================================================================================
  subroutine map_igp_apply_g2l_r2(mp,ld,xg,xl,renumbering)
    implicit none
    type(map_igp_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    real(rp)   , intent(in)  :: xg(ld,mp%ng)
    real(rp)   , intent(out) :: xl(ld,mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if
  end subroutine map_igp_apply_g2l_r2

  !================================================================================================
  subroutine map_igp_apply_g2l_i1(mp,xg,xl,renumbering)
    implicit none
    type(map_igp_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: xg(mp%ng)
    integer(ip), intent(out) :: xl(mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(i)=xg(renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_igp_apply_g2l_i1

  !================================================================================================
  subroutine map_igp_apply_g2l_i2(mp,ld,xg,xl,renumbering)
    implicit none
    type(map_igp_t)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    integer(ip), intent(in)  :: xg(ld,mp%ng)
    integer(ip), intent(out) :: xl(ld,mp%nl)
    type(renumbering_t), intent(in), optional :: renumbering
    integer(ip) :: i

    if(present(renumbering)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,renumbering%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if

  end subroutine map_igp_apply_g2l_i2


  !================================================================================================
  subroutine mesh_g2l_nmap_igp_emap_igp(nmap,emap,bmap,gmesh,lmesh,nrenumbering,erenumbering)
    implicit none
    type(map_igp_t)        , intent(in)  :: nmap, emap
    type(map_t)            , intent(in)  :: bmap
    type(mesh_t)       , intent(in)  :: gmesh
    type(mesh_t)       , intent(out) :: lmesh
    type(renumbering_t), optional, intent(in)  :: nrenumbering,erenumbering
    type(hash_table_igp_ip_t)      :: ws_inmap
    type(hash_table_igp_ip_t)      :: el_inmap
    integer(ip)                  :: aux, ipoin,inode,knode,ielem_lmesh,ielem_gmesh,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_gmesh,p_ipoin_gmesh,p_ipoin_lmesh, istat

    assert(nmap%ng == gmesh%npoin)
    assert(emap%ng == gmesh%nelem)

    ! both or none
    assert( (present(nrenumbering).and.present(erenumbering)) .or. ( (.not.present(nrenumbering)).and.(.not.present(erenumbering)) ) )


    lmesh%ndime=gmesh%ndime
    lmesh%nnode=gmesh%nnode
    lmesh%npoin=nmap%nl
    lmesh%nelem=emap%nl

    call memalloc (lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    call memalloc (lmesh%nnode*lmesh%nelem, lmesh%lnods, __FILE__,__LINE__)

    ! call memalloc (nmap%ng , ws_inmap, __FILE__,__LINE__)
    ! ws_inmap=-1
    call ws_inmap%init(max(int(nmap%nl*0.25,ip),10))
    do ipoin=1,nmap%nl
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call ws_inmap%put(key=nmap%l2g(ipoin),val=aux,stat=istat) 
       ! ws_inmap(nmap%l2g(ipoin))=ipoin
    end do

    ! call memalloc (emap%ng , el_inmap, __FILE__,__LINE__)
    ! el_inmap=-1
    call el_inmap%init(max(int(emap%nl*0.25,ip),10))
    do ipoin=1,emap%nl
       ! aux is used to avoid compiler warning related to val being an intent(inout) argument
       aux = ipoin
       call el_inmap%put(key=emap%l2g(ipoin),val=aux,stat=istat) 
       ! el_inmap(emap%l2g(ipoin))=ipoin
    end do

    if(present(nrenumbering).and.present(erenumbering)) then
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = erenumbering%iperm(emap%l2g(ielem_lmesh))
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(nrenumbering%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
             call ws_inmap%get(key=int(gmesh%lnods(p_ipoin_gmesh+inode),igp),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    else
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = emap%l2g(ielem_lmesh)
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(gmesh%lnods(p_ipoin_gmesh+inode))
             call ws_inmap%get(key=int(gmesh%lnods(p_ipoin_gmesh+inode),igp),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    end if
    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord, nrenumbering)

  end subroutine mesh_g2l_nmap_igp_emap_igp


  !================================================================================================
  subroutine mesh_g2l_emap_igp(nmap,emap,bmap,gmesh,lmesh,nrenumbering,erenumbering)
    implicit none
    type(map_t)            , intent(in)  :: nmap,bmap
    type(map_igp_t)        , intent(in)  :: emap
    type(mesh_t)       , intent(in)  :: gmesh
    type(mesh_t)       , intent(out) :: lmesh
    type(renumbering_t), optional, intent(in)  :: nrenumbering,erenumbering
    type(hash_table_ip_ip_t)       :: ws_inmap
    type(hash_table_igp_ip_t)      :: el_inmap
    integer(ip)                  :: ipoin,inode,knode,ielem_lmesh,ielem_gmesh,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_gmesh,p_ipoin_gmesh,p_ipoin_lmesh, istat

    assert(nmap%ng == gmesh%npoin)
    assert(emap%ng == gmesh%nelem)

    ! both or none
    assert( (present(nrenumbering).and.present(erenumbering)) .or. ( (.not.present(nrenumbering)).and.(.not.present(erenumbering)) ) )

    lmesh%ndime=gmesh%ndime
    lmesh%nnode=gmesh%nnode
    lmesh%npoin=nmap%nl
    lmesh%nelem=emap%nl

    call memalloc (lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    call memalloc (lmesh%nnode*lmesh%nelem, lmesh%lnods, __FILE__,__LINE__)

    ! call memalloc (nmap%ng , ws_inmap, __FILE__,__LINE__)
    ! ws_inmap=-1
    call ws_inmap%init(max(int(nmap%nl*0.25,ip),10))
    do ipoin=1,nmap%nl
       call ws_inmap%put(key=nmap%l2g(ipoin),val=ipoin,stat=istat) 
       ! ws_inmap(nmap%l2g(ipoin))=ipoin
    end do

    ! call memalloc (emap%ng , el_inmap, __FILE__,__LINE__)
    ! el_inmap=-1
    call el_inmap%init(max(int(emap%nl*0.25,ip),10))
    do ipoin=1,emap%nl
       call el_inmap%put(key=emap%l2g(ipoin),val=ipoin,stat=istat) 
       ! el_inmap(emap%l2g(ipoin))=ipoin
    end do

    if(present(nrenumbering).and.present(erenumbering)) then
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = erenumbering%iperm(emap%l2g(ielem_lmesh))
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(nrenumbering%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
             call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    else
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = emap%l2g(ielem_lmesh)
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(gmesh%lnods(p_ipoin_gmesh+inode))
             call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    end if


    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord,nrenumbering)

  end subroutine mesh_g2l_emap_igp

  !================================================================================================
  subroutine mesh_g2l_emap_ip(nmap,emap,bmap,gmesh,lmesh,nrenumbering,erenumbering)
    implicit none
    type(map_t)            , intent(in)  :: nmap,bmap
    type(map_t)            , intent(in)  :: emap
    type(mesh_t)       , intent(in)  :: gmesh
    type(mesh_t)       , intent(out) :: lmesh
    type(renumbering_t), optional, intent(in)  :: nrenumbering,erenumbering
    type(hash_table_ip_ip_t)       :: ws_inmap
    type(hash_table_ip_ip_t)       :: el_inmap
    integer(ip)                  :: ipoin,inode,knode,ielem_lmesh,ielem_gmesh,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_gmesh,p_ipoin_gmesh,p_ipoin_lmesh, istat

    assert(nmap%ng == gmesh%npoin)
    assert(emap%ng == gmesh%nelem)

    ! both or none
    assert( (present(nrenumbering).and.present(erenumbering)) .or. ( (.not.present(nrenumbering)).and.(.not.present(erenumbering)) ) )


    lmesh%ndime=gmesh%ndime
    lmesh%nnode=gmesh%nnode
    lmesh%npoin=nmap%nl
    lmesh%nelem=emap%nl

    call memalloc (lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    call memalloc (lmesh%nnode*lmesh%nelem, lmesh%lnods, __FILE__,__LINE__)

    ! call memalloc (nmap%ng , ws_inmap, __FILE__,__LINE__)
    ! ws_inmap=-1
    call ws_inmap%init(max(int(nmap%nl*0.25,ip),10))
    do ipoin=1,nmap%nl
       call ws_inmap%put(key=nmap%l2g(ipoin),val=ipoin,stat=istat) 
       ! ws_inmap(nmap%l2g(ipoin))=ipoin
    end do

    ! call memalloc (emap%ng , el_inmap, __FILE__,__LINE__)
    ! el_inmap=-1
    call el_inmap%init(max(int(emap%nl*0.25,ip),10))
    do ipoin=1,emap%nl
       call el_inmap%put(key=emap%l2g(ipoin),val=ipoin,stat=istat) 
       ! el_inmap(emap%l2g(ipoin))=ipoin
    end do

    if(present(nrenumbering).and.present(erenumbering)) then
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = erenumbering%iperm(emap%l2g(ielem_lmesh))
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(nrenumbering%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
             call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    else
       lmesh%pnods=0
       lmesh%pnods(1)=1
       do ielem_lmesh=1,lmesh%nelem
          ielem_gmesh = emap%l2g(ielem_lmesh)
          p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
          p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
          knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
          lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
          do inode=1,knode
             ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(gmesh%lnods(p_ipoin_gmesh+inode))
             call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
          end do
       end do
    end if
    
    call ws_inmap%free
    call el_inmap%free

    call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
    call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord,nrenumbering)

  end subroutine mesh_g2l_emap_ip



  !================================================================================================
  subroutine mesh_l2l(nrenumbering,erenumbering,lmeshin,lmeshout)
    implicit none
    type(renumbering_t)    , intent(in)  :: nrenumbering,erenumbering
    type(mesh_t) , intent(in)  :: lmeshin
    type(mesh_t) , intent(out) :: lmeshout

    integer(ip)                  :: ipoin,inode,knode,ielem_lmeshout,ielem_lmeshin,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_lmeshin,p_ipoin_lmeshin,p_ipoin_lmeshout

    lmeshout%ndime=lmeshin%ndime
    lmeshout%nnode=lmeshin%nnode
    lmeshout%npoin=lmeshin%npoin
    lmeshout%nelem=lmeshin%nelem
  
    call memalloc (lmeshout%nelem+1, lmeshout%pnods, __FILE__,__LINE__)
    call memalloc (lmeshin%pnods(lmeshin%nelem+1)-1, lmeshout%lnods, __FILE__,__LINE__)

    lmeshout%pnods=0
    lmeshout%pnods(1)=1
    do ielem_lmeshout=1,lmeshout%nelem
       ielem_lmeshin = erenumbering%iperm(ielem_lmeshout)
       p_ipoin_lmeshin = lmeshin%pnods(ielem_lmeshin)-1
       p_ipoin_lmeshout = lmeshout%pnods(ielem_lmeshout)-1
       knode = lmeshin%pnods(ielem_lmeshin+1)-lmeshin%pnods(ielem_lmeshin)
       lmeshout%pnods(ielem_lmeshout+1)=lmeshout%pnods(ielem_lmeshout)+knode
       do inode=1,knode
          lmeshout%lnods(p_ipoin_lmeshout+inode) =  nrenumbering%lperm(lmeshin%lnods(p_ipoin_lmeshin+inode))
       end do
    end do

    if ( allocated(lmeshin%coord) ) then 
      call memalloc(lmeshout%ndime, lmeshout%npoin, lmeshout%coord, __FILE__,__LINE__)
      call renumbering_apply (lmeshin%ndime, nrenumbering, lmeshin%coord,  lmeshout%coord)
    end if
  end subroutine mesh_l2l

end module map_apply_names
