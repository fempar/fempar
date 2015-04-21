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
module map_apply
  use types
  use memor
  use maps_names
  use renum_names
  use fem_mesh_names
  use hash_table_names
# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: l2g_add=0
  integer(ip), parameter :: l2g_copy=1

  interface map_apply_l2g
     module procedure map_apply_l2g_r1, map_apply_l2g_r2, map_apply_l2g_i1
  end interface

  interface map_apply_g2l
     module procedure map_apply_g2l_r1, map_apply_g2l_r2, map_apply_g2l_i1, map_apply_g2l_i2
  end interface

  interface fem_mesh_g2l
     module procedure fem_mesh_g2l_emap_ip, fem_mesh_g2l_emap_igp
  end interface

  ! Constants
  public :: l2g_add, l2g_copy

  ! Functions
  public :: map_apply_l2g, map_apply_g2l, fem_mesh_g2l, fem_mesh_l2l

contains

  !=============================================================================
  subroutine map_apply_l2g_r1(ope,mp,xl,xg,ren)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map)  , intent(in)  :: mp
    real(rp)   , intent(in)  :: xl(mp%nl)
    real(rp)   , intent(out) :: xg(mp%ng)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(ren%iperm(mp%l2g(i)))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(ren%iperm(mp%l2g(i)))=xg(ren%iperm(mp%l2g(i)))+xl(i)
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
  subroutine map_apply_l2g_i1(ope,mp,xl,xg,ren)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map)  , intent(in)  :: mp
    integer(ip), intent(in)  :: xl(mp%nl)
    integer(ip), intent(out) :: xg(mp%ng)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(ren%iperm(mp%l2g(i)))=xl(i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(ren%iperm(mp%l2g(i)))=xg(ren%iperm(mp%l2g(i)))+xl(i)
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
  subroutine map_apply_l2g_r2(ope,mp,ld,xl,xg,ren)
    implicit none
    integer(ip), intent(in)  :: ope
    type(map)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    real(rp)   , intent(in)  :: xl(ld,mp%nl)
    real(rp)   , intent(out) :: xg(ld,mp%ng)
    type(renum), intent(in), optional :: ren
    integer(ip)                   :: i

    if(present(ren)) then
       if(ope==l2g_copy) then
          do i=1,mp%nl
             xg(:,ren%iperm(mp%l2g(i)))=xl(:,i)
          end do
       else if(ope==l2g_add) then
          do i=1,mp%nl
             xg(:,ren%iperm(mp%l2g(i)))=xg(:,ren%iperm(mp%l2g(i)))+xl(:,i)
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
  subroutine map_apply_g2l_r1(mp,xg,xl,ren)
    implicit none
    type(map)  , intent(in)  :: mp
    real(rp)   , intent(in)  :: xg(mp%ng)
    real(rp)   , intent(out) :: xl(mp%nl)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       do i=1,mp%nl
          xl(i)=xg(ren%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_r1

  !================================================================================================
  subroutine map_apply_g2l_r2(mp,ld,xg,xl,ren)
    implicit none
    type(map)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    real(rp)   , intent(in)  :: xg(ld,mp%ng)
    real(rp)   , intent(out) :: xl(ld,mp%nl)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,ren%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if
  end subroutine map_apply_g2l_r2

  !================================================================================================
  subroutine map_apply_g2l_i1(mp,xg,xl,ren)
    implicit none
    type(map)  , intent(in)  :: mp
    integer(ip), intent(in)  :: xg(mp%ng)
    integer(ip), intent(out) :: xl(mp%nl)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       do i=1,mp%nl
          xl(i)=xg(ren%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(i)=xg(mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_i1

  !================================================================================================
  subroutine map_apply_g2l_i2(mp,ld,xg,xl,ren)
    implicit none
    type(map)  , intent(in)  :: mp
    integer(ip), intent(in)  :: ld
    integer(ip), intent(in)  :: xg(ld,mp%ng)
    integer(ip), intent(out) :: xl(ld,mp%nl)
    type(renum), intent(in), optional :: ren
    integer(ip) :: i

    if(present(ren)) then
       do i=1,mp%nl
          xl(:,i)=xg(:,ren%iperm(mp%l2g(i)))
       end do
    else
       do i=1,mp%nl
          xl(:,i)=xg(:,mp%l2g(i))
       end do
    end if

  end subroutine map_apply_g2l_i2

  !================================================================================================
  subroutine fem_mesh_g2l_emap_igp(nmap,emap,bmap,gmesh,lmesh,nren,eren)
    implicit none
    type(map)            , intent(in)  :: nmap,bmap
    type(map_igp)        , intent(in)  :: emap
    type(fem_mesh)       , intent(in)  :: gmesh
    type(fem_mesh)       , intent(out) :: lmesh
    type(renum), optional, intent(in)  :: nren,eren
    type(hash_table_ip_ip)       :: ws_inmap
    type(hash_table_igp_ip)      :: el_inmap
    integer(ip)                  :: ipoin,inode,knode,ielem_lmesh,ielem_gmesh,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_gmesh,p_ipoin_gmesh,p_ipoin_lmesh, istat

    assert(nmap%ng == gmesh%npoin)
    assert(emap%ng == gmesh%nelem)

    ! both or none
    assert( (present(nren).and.present(eren)) .or. ( (.not.present(nren)).and.(.not.present(eren)) ) )

    ! We could be more precise and check how many
    ! element types we have in the local mesh...
    lmesh%nelty=gmesh%nelty

    lmesh%ndime=gmesh%ndime
    lmesh%nnode=gmesh%nnode
    lmesh%npoin=nmap%nl
    lmesh%nelem=emap%nl

    if(lmesh%nelty==1) then
       call memalloc (            1, lmesh%pnods, __FILE__,__LINE__)
    else
       call memalloc (lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    endif
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

    if(present(nren).and.present(eren)) then
       if(lmesh%nelty==1) then
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = eren%iperm(emap%l2g(ielem_lmesh))
             p_ipoin_gmesh = (ielem_gmesh-1)*gmesh%nnode
             knode = gmesh%nnode
             do inode=1,knode
                ! lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode) =  ws_inmap(nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
                call ws_inmap%get(key=nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)),val=lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode),stat=istat) 
             end do
          end do
       else 
          lmesh%pnods=0
          lmesh%pnods(1)=1
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = eren%iperm(emap%l2g(ielem_lmesh))
             p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
             p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
             knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
             lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
             do inode=1,knode
                ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
             end do
          end do
       end if
    else
       if(lmesh%nelty==1) then
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = emap%l2g(ielem_lmesh)
             p_ipoin_gmesh = (ielem_gmesh-1)*gmesh%nnode
             knode = gmesh%nnode
             do inode=1,knode
                ! lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode) =  ws_inmap(gmesh%lnods(p_ipoin_gmesh+inode))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode),stat=istat) 
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
    end if

    if ( gmesh%nboun > 0 ) then 
       ! Boundary elements
       lmesh%nboun = bmap%nl
       lmesh%nnodb = gmesh%nnodb
       call memalloc(1,lmesh%pboun,__FILE__,__LINE__)
       call memalloc(lmesh%nnodb*lmesh%nboun,lmesh%lboun,__FILE__,__LINE__)
       call memalloc(lmesh%nnodb+1,lmesh%nboun,lmesh%lboel,__FILE__,__LINE__)

       ! g2l
       lmesh%pboun = 1
       if(present(nren).and.present(eren)) then
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             do inode=1,lmesh%nnodb
                gnode = gmesh%lboun((gelem-1)*lmesh%nnodb+inode)
                ! lmesh%lboun((iboun-1)*lmesh%nnodb+inode) = ws_inmap(nren%lperm(gnode))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
             end do
          end do
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             velem = gmesh%lboel(lmesh%nnodb+1,gelem)
             do inode=1,lmesh%nnodb
                lmesh%lboel(inode,iboun) = gmesh%lboel(inode,gelem)
             end do
             ! lmesh%lboel(lmesh%nnodb+1,iboun) = el_inmap(eren%lperm(velem))
             call ws_inmap%get(key=eren%lperm(velem),val=lmesh%lboel(lmesh%nnodb+1,iboun),stat=istat) 
          end do
       else
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             do inode=1,lmesh%nnodb
                gnode = gmesh%lboun((gelem-1)*lmesh%nnodb+inode)
                ! lmesh%lboun((iboun-1)*lmesh%nnodb+inode) = ws_inmap(gnode)
                call ws_inmap%get(key=gnode,val=lmesh%lboun((iboun-1)*lmesh%nnodb+inode),stat=istat) 
             end do
          end do
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             velem = gmesh%lboel(lmesh%nnodb+1,gelem)
             do inode=1,lmesh%nnodb
                lmesh%lboel(inode,iboun) = gmesh%lboel(inode,gelem)
             end do
             ! lmesh%lboel(lmesh%nnodb+1,iboun) = el_inmap(velem)
             call el_inmap%get(key=int(velem,igp),val=lmesh%lboel(lmesh%nnodb+1,iboun),stat=istat) 
          end do
       end if
    else
       lmesh%nboun = 0
       lmesh%nnodb = 0
    end if
    lmesh%nelpo = 0
    
    call ws_inmap%free
    call el_inmap%free

    if ( allocated(gmesh%coord) ) then
       call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
       call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord,nren)
    end if

  end subroutine fem_mesh_g2l_emap_igp

  !================================================================================================
  subroutine fem_mesh_g2l_emap_ip(nmap,emap,bmap,gmesh,lmesh,nren,eren)
    implicit none
    type(map)            , intent(in)  :: nmap,bmap
    type(map)            , intent(in)  :: emap
    type(fem_mesh)       , intent(in)  :: gmesh
    type(fem_mesh)       , intent(out) :: lmesh
    type(renum), optional, intent(in)  :: nren,eren
    type(hash_table_ip_ip)       :: ws_inmap
    type(hash_table_ip_ip)       :: el_inmap
    integer(ip)                  :: ipoin,inode,knode,ielem_lmesh,ielem_gmesh,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_gmesh,p_ipoin_gmesh,p_ipoin_lmesh, istat

    assert(nmap%ng == gmesh%npoin)
    assert(emap%ng == gmesh%nelem)

    ! both or none
    assert( (present(nren).and.present(eren)) .or. ( (.not.present(nren)).and.(.not.present(eren)) ) )

    ! We could be more precise and check how many
    ! element types we have in the local mesh...
    lmesh%nelty=gmesh%nelty

    lmesh%ndime=gmesh%ndime
    lmesh%nnode=gmesh%nnode
    lmesh%npoin=nmap%nl
    lmesh%nelem=emap%nl

    if(lmesh%nelty==1) then
       call memalloc (            1, lmesh%pnods, __FILE__,__LINE__)
    else
       call memalloc (lmesh%nelem+1, lmesh%pnods, __FILE__,__LINE__)
    endif
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

    if(present(nren).and.present(eren)) then
       if(lmesh%nelty==1) then
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = eren%iperm(emap%l2g(ielem_lmesh))
             p_ipoin_gmesh = (ielem_gmesh-1)*gmesh%nnode
             knode = gmesh%nnode
             do inode=1,knode
                ! lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode) =  ws_inmap(nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
                call ws_inmap%get(key=nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)),val=lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode),stat=istat) 
             end do
          end do
       else 
          lmesh%pnods=0
          lmesh%pnods(1)=1
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = eren%iperm(emap%l2g(ielem_lmesh))
             p_ipoin_gmesh = gmesh%pnods(ielem_gmesh)-1
             p_ipoin_lmesh = lmesh%pnods(ielem_lmesh)-1
             knode = gmesh%pnods(ielem_gmesh+1)-gmesh%pnods(ielem_gmesh)
             lmesh%pnods(ielem_lmesh+1)=lmesh%pnods(ielem_lmesh)+knode
             do inode=1,knode
                ! lmesh%lnods(p_ipoin_lmesh+inode) =  ws_inmap(nren%lperm(gmesh%lnods(p_ipoin_gmesh+inode)))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
             end do
          end do
       end if
    else
       if(lmesh%nelty==1) then
          do ielem_lmesh=1,lmesh%nelem
             ielem_gmesh = emap%l2g(ielem_lmesh)
             p_ipoin_gmesh = (ielem_gmesh-1)*gmesh%nnode
             knode = gmesh%nnode
             do inode=1,knode
                ! lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode) =  ws_inmap(gmesh%lnods(p_ipoin_gmesh+inode))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods((ielem_lmesh-1)*lmesh%nnode +inode),stat=istat) 
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
    end if

    if ( gmesh%nboun > 0 ) then 
       ! Boundary elements
       lmesh%nboun = bmap%nl
       lmesh%nnodb = gmesh%nnodb
       call memalloc(1,lmesh%pboun,__FILE__,__LINE__)
       call memalloc(lmesh%nnodb*lmesh%nboun,lmesh%lboun,__FILE__,__LINE__)
       call memalloc(lmesh%nnodb+1,lmesh%nboun,lmesh%lboel,__FILE__,__LINE__)

       ! g2l
       lmesh%pboun = 1
       if(present(nren).and.present(eren)) then
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             do inode=1,lmesh%nnodb
                gnode = gmesh%lboun((gelem-1)*lmesh%nnodb+inode)
                ! lmesh%lboun((iboun-1)*lmesh%nnodb+inode) = ws_inmap(nren%lperm(gnode))
                call ws_inmap%get(key=gmesh%lnods(p_ipoin_gmesh+inode),val=lmesh%lnods(p_ipoin_lmesh+inode),stat=istat) 
             end do
          end do
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             velem = gmesh%lboel(lmesh%nnodb+1,gelem)
             do inode=1,lmesh%nnodb
                lmesh%lboel(inode,iboun) = gmesh%lboel(inode,gelem)
             end do
             ! lmesh%lboel(lmesh%nnodb+1,iboun) = el_inmap(eren%lperm(velem))
             call ws_inmap%get(key=eren%lperm(velem),val=lmesh%lboel(lmesh%nnodb+1,iboun),stat=istat) 
          end do
       else
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             do inode=1,lmesh%nnodb
                gnode = gmesh%lboun((gelem-1)*lmesh%nnodb+inode)
                ! lmesh%lboun((iboun-1)*lmesh%nnodb+inode) = ws_inmap(gnode)
                call ws_inmap%get(key=gnode,val=lmesh%lboun((iboun-1)*lmesh%nnodb+inode),stat=istat) 
             end do
          end do
          do iboun=1,lmesh%nboun
             gelem = bmap%l2g(iboun)
             velem = gmesh%lboel(lmesh%nnodb+1,gelem)
             do inode=1,lmesh%nnodb
                lmesh%lboel(inode,iboun) = gmesh%lboel(inode,gelem)
             end do
             ! lmesh%lboel(lmesh%nnodb+1,iboun) = el_inmap(velem)
             call el_inmap%get(key=int(velem,igp),val=lmesh%lboel(lmesh%nnodb+1,iboun),stat=istat) 
          end do
       end if
    else
       lmesh%nboun = 0
       lmesh%nnodb = 0
    end if
    lmesh%nelpo = 0
    
    call ws_inmap%free
    call el_inmap%free

    if ( allocated(gmesh%coord) ) then
       call memalloc(lmesh%ndime, lmesh%npoin, lmesh%coord, __FILE__,__LINE__)
       call map_apply_g2l(nmap, gmesh%ndime, gmesh%coord, lmesh%coord,nren)
    end if

  end subroutine fem_mesh_g2l_emap_ip

  !================================================================================================
  subroutine fem_mesh_l2l(nren,eren,lmeshin,lmeshout)
    implicit none
    type(renum)    , intent(in)  :: nren,eren
    type(fem_mesh) , intent(in)  :: lmeshin
    type(fem_mesh) , intent(out) :: lmeshout

    integer(ip)                  :: ipoin,inode,knode,ielem_lmeshout,ielem_lmeshin,iboun,gelem,velem,gnode
    integer(ip)                  :: p_ielem_lmeshin,p_ipoin_lmeshin,p_ipoin_lmeshout

    ! We could be more precise and check how many
    ! element types we have in the local mesh...
    lmeshout%nelty=lmeshin%nelty

    lmeshout%ndime=lmeshin%ndime
    lmeshout%nnode=lmeshin%nnode
    lmeshout%npoin=lmeshin%npoin
    lmeshout%nelem=lmeshin%nelem

    if(lmeshout%nelty==1) then
       call memalloc (            1, lmeshout%pnods, __FILE__,__LINE__)
       call memalloc (lmeshout%nnode*lmeshout%nelem, lmeshout%lnods, __FILE__,__LINE__)
    else
       call memalloc (lmeshout%nelem+1, lmeshout%pnods, __FILE__,__LINE__)
       call memalloc (lmeshin%pnods(lmeshin%nelem+1)-1, lmeshout%lnods, __FILE__,__LINE__)
    endif

    if(lmeshout%nelty==1) then
       ! write(*,*) nren%lperm
       do ielem_lmeshout=1,lmeshout%nelem
          ielem_lmeshin = eren%iperm(ielem_lmeshout)
          p_ipoin_lmeshin = (ielem_lmeshin-1)*lmeshin%nnode
          knode = lmeshin%nnode
          do inode=1,knode
             lmeshout%lnods((ielem_lmeshout-1)*lmeshout%nnode +inode) =  nren%lperm(lmeshin%lnods(p_ipoin_lmeshin+inode))
             ! write (*,*) 'PPP', nren%lperm(lmeshin%lnods(p_ipoin_lmeshin+inode))
          end do
       end do
    else 
       lmeshout%pnods=0
       lmeshout%pnods(1)=1
       do ielem_lmeshout=1,lmeshout%nelem
          ielem_lmeshin = eren%iperm(ielem_lmeshout)
          p_ipoin_lmeshin = lmeshin%pnods(ielem_lmeshin)-1
          p_ipoin_lmeshout = lmeshout%pnods(ielem_lmeshout)-1
          knode = lmeshin%pnods(ielem_lmeshin+1)-lmeshin%pnods(ielem_lmeshin)
          lmeshout%pnods(ielem_lmeshout+1)=lmeshout%pnods(ielem_lmeshout)+knode
          do inode=1,knode
             lmeshout%lnods(p_ipoin_lmeshout+inode) =  nren%lperm(lmeshin%lnods(p_ipoin_lmeshin+inode))
          end do
       end do
    end if

!!$    ****** PENDING ******
!!$    ! Boundary elements
!!$    lmeshout%nboun = bmap%nl
!!$    lmeshout%nnodb = lmeshin%nnodb
!!$    call memalloc(1,lmeshout%pboun,__FILE__,__LINE__)
!!$    call memalloc(lmeshout%nnodb*lmeshout%nboun,lmeshout%lboun,__FILE__,__LINE__)
!!$    call memalloc(lmeshout%nnodb+1,lmeshout%nboun,lmeshout%lboel,__FILE__,__LINE__)
!!$
!!$    ! l2l
!!$    lmeshout%pboun = 1
!!$    do iboun=1,lmeshout%nboun
!!$       gelem = bmap%l2g(iboun)
!!$       do inode=1,lmeshout%nnodb
!!$          gnode = lmeshin%lboun((gelem-1)*lmeshout%nnodb+inode)
!!$          lmeshout%lboun((iboun-1)*lmeshout%nnodb+inode) = ws_inmap(nren%lperm(gnode))
!!$       end do
!!$    end do
!!$    do iboun=1,lmeshout%nboun
!!$       gelem = bmap%l2g(iboun)
!!$       velem = lmeshin%lboel(lmeshout%nnodb+1,gelem)
!!$       do inode=1,lmeshout%nnodb
!!$          lmeshout%lboel(inode,iboun) = lmeshin%lboel(inode,gelem)
!!$       end do
!!$       lmeshout%lboel(lmeshout%nnodb+1,iboun) = el_inmap(eren%lperm(velem))
!!$    end do

    if ( allocated ( lmeshin%coord ) ) then
       call memalloc(lmeshout%ndime, lmeshout%npoin, lmeshout%coord, __FILE__,__LINE__)
       call renum_apply (lmeshin%ndime, nren, lmeshin%coord,  lmeshout%coord)
    end if

  end subroutine fem_mesh_l2l

end module map_apply
