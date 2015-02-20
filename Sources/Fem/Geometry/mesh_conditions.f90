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
module mesh_conditions
  use types
  use stdio
  use memor
  use fem_mesh_class
  use fem_conditions_class
  use fem_per_trans_class
  implicit none
  private

  ! Functions
  public :: fem_mesh_conditions_per_trans

contains
  !=============================================================================
  subroutine fem_mesh_conditions_per_trans (isper,gmesh,gnodes,umesh,unodes,f_p_tr)
    !-----------------------------------------------------------------------
    !  This routine computes the periodic mesh and boundary conditions
    !-----------------------------------------------------------------------
    implicit none
    type(fem_mesh)      , intent(in)    :: gmesh
    integer(ip)         , intent(in)    :: isper(3)
    type(fem_conditions), intent(in)    :: gnodes
    type(fem_mesh)      , intent(out)   :: umesh
    type(fem_conditions), intent(out)   :: unodes
    type(fem_per_trans) , intent(out)   :: f_p_tr

    integer(ip)   :: ind,jnd,knd,indb,indc,fic_coor(2),coord1,coord2,nfict,upoin,kfl_perio
    integer(ip), allocatable :: touched(:)
    real(rp)      :: amin,amax,apos,bmax,bmin,bpos,cpos

    ! Periodic cases
    if(gmesh%ndime==2) then
       if(isper(1)==1.and.isper(2)==0) then
          kfl_perio=1
       else if(isper(1)==0.and.isper(2)==1) then
          kfl_perio=2
       else if(isper(1)==1.and.isper(2)==1) then
          kfl_perio=4
       end if
    else if(gmesh%ndime==3) then
       ! FAENA PARA ORIOL !
    end if

    ! Set the dimension with periodic bc's
    fic_coor(2)=0
    if(kfl_perio==1)then              !a=x,b=y,c=z
       fic_coor(1)=1          !X-coord
       coord1=2               !Y-coord
       if(gmesh%ndime==3) then
          coord2=3            !Z-coord
       end if
    else if(kfl_perio==2) then        !a=y,b=x,c=z
       fic_coor(1)=2          !Y-coord
       coord1=1               !X-coord
       if(gmesh%ndime==3) then
          coord2=3            !Z-coord
       end if
    else if(kfl_perio==3) then        !a=z,b=x,c=y
       fic_coor(1)=3          !Z-coord
       coord1=1               !X-coord
       coord2=2               !Y-coord
    else if(kfl_perio==4) then        !a=x,b=y,c=z
       fic_coor(1)=1          !X-coord
       fic_coor(2)=2
       coord1=2               !Y-coord
       if(gmesh%ndime==3) then
          coord2=3            !Z-coord
       end if
    end if

    call fem_per_trans_alloc (gmesh%npoin, f_p_tr)

    ! Find amin, amax (the faces to be linked) & initialize f_p_tr%fictio
    amin = 9999.9_rp
    amax = -9999.9_rp
    do ind=1,gmesh%npoin
       amin=min(amin,gmesh%coord(fic_coor(1),ind))
       amax=max(amax,gmesh%coord(fic_coor(1),ind))
       f_p_tr%fictio(ind)=ind
    end do
    if(kfl_perio==4) then
       bmin = 9999.9_rp
       bmax = -9999.9_rp
       do ind=1,gmesh%npoin
          bmin=min(bmin,gmesh%coord(fic_coor(2),ind))
          bmax=max(bmax,gmesh%coord(fic_coor(2),ind))   
       end do
    end if

    ! Find the correspondence between nodes in both faces
    nfict=0
    if (gmesh%ndime==3) then
       do ind=1,gmesh%npoin
          if(kfl_perio/=4) then
             if (abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14) then
                bpos=gmesh%coord(coord1,ind)
                cpos=gmesh%coord(coord2,ind)
                do indb=1,gmesh%npoin
                   if ((abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14).and.&
                        (abs(bpos-gmesh%coord(coord1,indb)).le.1e-6).and.&
                        (abs(cpos-gmesh%coord(coord2,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          else if(kfl_perio==4) then
             if ((abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14).and.&
                  (abs(gmesh%coord(coord1,ind)-bmax).ge.1e-14)) then
                bpos=gmesh%coord(coord1,ind)
                cpos=gmesh%coord(coord2,ind)
                do indb=1,gmesh%npoin
                   if ((abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14).and.&
                        (abs(bpos-gmesh%coord(coord1,indb)).le.1e-6).and.&
                        (abs(cpos-gmesh%coord(coord2,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          end if
       end do
    elseif (gmesh%ndime==2) then
       do ind=1,gmesh%npoin
          if(kfl_perio/=4) then
             if (abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14) then
                bpos=gmesh%coord(coord1,ind)
                do indb=1,gmesh%npoin
                   if ((abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14).and.&
                        (abs(bpos-gmesh%coord(coord1,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          else if(kfl_perio==4) then
             if ((abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14).and.&
                  (abs(gmesh%coord(coord1,ind)-bmax).ge.1e-14)) then
                bpos=gmesh%coord(coord1,ind)
                do indb=1,gmesh%npoin
                   if ((abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14).and.&
                        (abs(bpos-gmesh%coord(coord1,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          end if
       end do
    end if
    if(kfl_perio==4) then
       if (gmesh%ndime==3) then
          do ind=1,gmesh%npoin
             if ((abs(gmesh%coord(fic_coor(2),ind)-bmax).le.1e-14).and.&
                  (abs(gmesh%coord(fic_coor(1),ind)-amax).ge.1e-14))  then
                apos=gmesh%coord(fic_coor(1),ind)
                cpos=gmesh%coord(coord2,ind)
                do indb=1,gmesh%npoin
                   if ((abs(bmin-gmesh%coord(fic_coor(2),indb)).le.1e-14).and.&
                        (abs(apos-gmesh%coord(fic_coor(1),indb)).le.1e-6).and.&
                        (abs(cpos-gmesh%coord(coord2,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
             if((abs(gmesh%coord(fic_coor(2),ind)-bmax).le.1e-14).and.&
                  (abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14)) then
                cpos=gmesh%coord(coord2,ind)
                do indb=1,gmesh%npoin
                   if ((abs(bmin-gmesh%coord(fic_coor(2),indb)).le.1e-14).and.&
                        (abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14).and.&
                        (abs(cpos-gmesh%coord(coord2,indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          end do
       elseif (gmesh%ndime==2) then
          do ind=1,gmesh%npoin
             if ((abs(gmesh%coord(fic_coor(2),ind)-bmax).le.1e-14).and.&
                  (abs(gmesh%coord(fic_coor(1),ind)-amax).ge.1e-14)) then
                apos=gmesh%coord(fic_coor(1),ind)
                do indb=1,gmesh%npoin
                   if ((abs(bmin-gmesh%coord(fic_coor(2),indb)).le.1e-14).and.&
                        (abs(apos-gmesh%coord(fic_coor(1),indb)).le.1e-6)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
             if((abs(gmesh%coord(fic_coor(2),ind)-bmax).le.1e-14).and.&
                  (abs(gmesh%coord(fic_coor(1),ind)-amax).le.1e-14)) then
                do indb=1,gmesh%npoin
                   if ((abs(bmin-gmesh%coord(fic_coor(2),indb)).le.1e-14).and.&
                        (abs(amin-gmesh%coord(fic_coor(1),indb)).le.1e-14)) then
                      f_p_tr%fictio(ind)=indb
                      nfict=nfict+1
                   end if
                end do
             end if
          end do
       end if
    end if

    ! Allocate umesh
    upoin=gmesh%npoin-nfict
    call fem_mesh_alloc(gmesh%ndime,gmesh%nelty,gmesh%nnode,upoin,gmesh%nelem,umesh)

    ! Create umesh and unodes
    umesh%ndime = gmesh%ndime
    umesh%nelty = gmesh%nelty
    umesh%nnode = gmesh%nnode
    umesh%npoin = upoin
    umesh%nelem = gmesh%nelem
    
    call fem_conditions_create(gnodes%ncode,gnodes%nvalu,umesh%npoin,unodes)

    ! here, f_p_tr%fictio is the correspondence 
    ! among original mesh numbering and 
    ! numbering corresponding to the corresponding
    ! periodic mesh BEFORE making ids consecutive
    knd = 1

    call memalloc(gmesh%npoin,touched,__FILE__,__LINE__) 
    touched = 0
    do ind=1,gmesh%npoin
       if( ind == f_p_tr%fictio(ind) ) then
          f_p_tr%fictio(ind) = knd
          umesh%coord(:,knd) = gmesh%coord(:,ind)
          unodes%code(:,knd) = gnodes%code(:,ind)
          unodes%valu(:,knd) = gnodes%valu(:,ind)
          knd = knd+1
          touched(ind) = 1
       end if
    end do
  

    do ind=1,gmesh%npoin
       if (touched(ind) == 0) then 
          f_p_tr%fictio(ind) = f_p_tr%fictio( f_p_tr%fictio(ind) )
          unodes%code(:,f_p_tr%fictio(ind)) = 0
       end if
    end do
    call memfree(touched,__FILE__,__LINE__) 


    do ind=1,gmesh%nelem*gmesh%nnode
       umesh%lnods(ind) = f_p_tr%fictio(gmesh%lnods(ind))
    end do         

  end subroutine fem_mesh_conditions_per_trans

end module mesh_conditions
