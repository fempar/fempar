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
module postpro_names
use types_names
use stdio_names
#ifdef ENABLE_GIDPOST
use gidpost_names
#endif
  use renum_names
  implicit none
  private

  type post_file_t
     integer(ip)    :: form = -1
#ifdef ENABLE_GIDPOST
     type(gid_file) :: gid_luou
#endif
     integer(ip)    :: luou = -1
  end type post_file_t

  interface postpro
     module procedure possca,posvec,poisca
  end interface postpro

  interface postpro_gp
     module procedure pgpint,pgpr1p,pgpr3p
  end interface

  ! Types
  public :: post_file_t

  ! Methods
  public :: postpro_open_file, postpro_close_file, postpro, postpro_gp, &
       &    postpro_compose_result_name, postpro_compose_mesh_name, &
       &    postpro_gp_init

contains


  !=============================================================================
#ifndef ENABLE_GIDPOST
  subroutine enable_gidpost_error_message
      implicit none
      write (0,*) 'Error: Fem was not compiled with -DENABLE_GIDPOST.'
      write (0,*) "Error: You must activate this cpp macro in order to use GIDPOST"
      stop
  end subroutine
#endif

  !==================================================================================
  subroutine postpro_compose_result_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.pos.res'
  end subroutine postpro_compose_result_name

  subroutine postpro_compose_mesh_name ( prefix, name ) 
    implicit none
    character *(*), intent(in)        :: prefix 
    character *(*), intent(out)       :: name
    name = trim(prefix) // '.pos.msh'
  end subroutine postpro_compose_mesh_name

  !==================================================================================
  subroutine postpro_open_file(form,name,pos)
    implicit none
    character(*)   , intent(in)  :: name
    integer(ip)    , intent(in)  :: form
    type(post_file_t), intent(out) :: pos

    pos%form=form

    if(pos%form==1) then
       pos%luou =  io_open(adjustl(name))
       write(pos%luou,'(a)') 'GiD Post Results File 1.0'
       write(pos%luou,'(a)') ' '
    else if(pos%form==2) then
       pos%luou =  io_open(adjustl(name))
    else if(pos%form==3) then
#ifdef ENABLE_GIDPOST
       CALL GID_POSTINIT()
       pos%gid_luou = GID_fOPENPOSTRESULTFILE(trim(adjustl(name)),GID_POSTASCII)
#else
       call enable_gidpost_error_message
#endif
    else if(pos%form==4) then
#ifdef ENABLE_GIDPOST
       CALL GID_POSTINIT()
       pos%gid_luou = GID_fOPENPOSTRESULTFILE(trim(adjustl(name)),GID_POSTBINARY)
#else
       call enable_gidpost_error_message
#endif
    end if

  end subroutine postpro_open_file

  logical function postpro_opened_file(pos)
    implicit none
    type(post_file_t), intent(in) :: pos
    if(pos%form==-1) then
       postpro_opened_file = .false.
    else
       postpro_opened_file = .true.
    end if
  end function postpro_opened_file
  !==================================================================================
  subroutine postpro_close_file(pos)
    implicit none
    type(post_file_t), intent(in) :: pos

    if(pos%form==1.or.pos%form==2) then
       call io_close(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
#ifdef ENABLE_GIDPOST
       call GiD_fClosePostResultFile(pos%gid_luou)
#else
       call enable_gidpost_error_message
#endif
    end if
  end subroutine postpro_close_file

  !==================================================================================
  subroutine possca(pos,bridge,wopos,istep,ttime,nren)
    !-----------------------------------------------------------------------
    !
    ! Write a scalar in postprocess file
    !
    !-----------------------------------------------------------------------
    implicit none
    type(post_file_t), intent(in)  :: pos
    character(*)   , intent(in)  :: wopos
    real(rp)       , intent(in)  :: bridge(:)
    integer(ip)    , intent(in)  :: istep
    real(rp)       , intent(in)  :: ttime
    type(renum_t),optional,intent(in) :: nren
    integer(ip)               :: npoin,ipoin
    character(8)              :: state

    npoin=size(bridge)
    state='ANALYSIS'

    if(pos%form==1) then
       write(pos%luou,2) wopos,state,ttime,'Scalar'
       write(pos%luou,3) wopos
       write(pos%luou,1) 'Values'
!!$       do ipoin=1,npoin
!!$          write(pos%luou,4)ipoin,bridge(ipoin)
!!$       end do
!!$       write(pos%luou,1) 'End Values'
       if ( present(nren) ) then
          do ipoin=1,npoin
             write(pos%luou,4) nren%lperm(ipoin),bridge(ipoin)
          end do
       else
          do ipoin=1,npoin
             write(pos%luou,4)ipoin,bridge(ipoin)
          end do
       end if
       write(pos%luou,1) 'End Values' 
       flush(pos%luou)
    else if(pos%form==2) then
       write(pos%luou,10) 100,'C','Step  ',ttime,1,istep
       write(pos%luou,20) -4,wopos,1,1,0
       write(pos%luou,30) -5,wopos,1,1
       do ipoin=1,npoin
          write(pos%luou,40) -1,ipoin,bridge(ipoin)
       end do
       write(pos%luou,50) -3
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
#ifdef ENABLE_GIDPOST
       CALL GID_fBEGINSCALARRESULT(pos%gid_luou,wopos,state,ttime,GID_ONNODES, &
            &                        GID_NULL,GID_NULL,GID_NULL)
       do ipoin=1,npoin
          CALL GID_fWRITESCALAR(pos%gid_luou,ipoin,bridge(ipoin))
       end do
       CALL GID_fENDRESULT(pos%gid_luou)
       CALL GID_fFLUSHPOSTFILE(pos%gid_luou)
#else
       call enable_gidpost_error_message
#endif
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e14.6,' ',a,' OnNodes')
3   format('ComponentNames ',a)
4   format(i7, 3(1x,e16.8E3))

    ! Femview formats
10  format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20  format(1x,i2,2x,a5,3x,3i5)
30  format(1x,i2,2x,a5,3x,2i5)
40  format(1x,i2,i5,e12.5)
50  format(1x,i2)

  end subroutine possca

  !==================================================================================
  subroutine posvec(pos,bridge,wopos,istep,ttime,kpoin,lpoty)
    !-----------------------------------------------------------------------
    !
    ! Write a vector in postprocess file
    !
    !-----------------------------------------------------------------------
    implicit none
    type(post_file_t), intent(in)  :: pos
    character(*)   , intent(in)  :: wopos
    real(rp)       , intent(in)  :: bridge(:,:)
    integer(ip)    , intent(in)  :: istep
    real(rp)       , intent(in)  :: ttime
    integer(ip)    , optional    :: kpoin,lpoty(:)
    integer(ip)                  :: npoin,ndime,ipoin,idime,kbopo,ibopo
    real(rp)                     :: dummr
    character(8)                 :: state

    ! Dimensions
    if(present(kpoin)) then
       npoin=kpoin
    else
       npoin=size(bridge,2)       
    end if
    ndime=size(bridge,1)

    ! Check if there is arguments
    if(present(lpoty)) then
       kbopo=1               ! bridge only defined on boundary nodes
       dummr=0.0_rp
    else
       kbopo=0               ! bridge is defined on all nodes
    end if

    state='ANALYSIS'
    if(pos%form==1) then
       write(pos%luou,2) wopos,state,ttime,'Vector'
       write(pos%luou,3) wopos//'_X',wopos//'_Y',wopos//'_Z'
       write(pos%luou,1) 'Values'
       if(kbopo==0) then
          do ipoin=1,npoin
             write(pos%luou,4)ipoin,(bridge(idime,ipoin),idime=1,ndime)
          end do
       else
          do ipoin=1,npoin
             ibopo=lpoty(ipoin)
             if(ibopo/=0) then
                write(pos%luou,4)ipoin,(bridge(idime,ibopo),idime=1,ndime)
             else
                write(pos%luou,4)ipoin,(dummr,idime=1,ndime)
             end if
          end do
       end if
       write(pos%luou,1) 'End Values'
       flush(pos%luou)
    else if(pos%form==2) then
       write(pos%luou,10) 100,'C','Step  ',ttime,1,istep
       write(pos%luou,20) -4,wopos,4,1,0
       write(pos%luou,30) -5,'X-COMP  ',1,2,1,0,0
       write(pos%luou,30) -5,'Y-COMP  ',1,2,2,0,0
       write(pos%luou,30) -5,'Z-COMP  ',1,2,3,0,0
       write(pos%luou,80) -5,'ALL     ',1,2,0,0,1,'all     '
       if(kbopo==0) then
          do ipoin = 1,npoin
             write(pos%luou,50) -1,ipoin,(bridge(idime,ipoin),idime=1,ndime)
          end do
       else
          do ipoin = 1,npoin
             ibopo=lpoty(ipoin)
             if(ibopo/=0) then
                write(pos%luou,50) -1,ipoin,(bridge(idime,ibopo),idime=1,ndime)
             else
                write(pos%luou,50) -1,ipoin,(dummr,idime=1,ndime)
             end if
          end do
       end if
       write(pos%luou,70) -3
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
#ifdef ENABLE_GIDPOST
       call GID_fBEGINVECTORRESULT(pos%gid_luou,wopos,state,ttime,GID_ONNODES,&
            &          GID_NULL,GID_NULL,wopos//'_X',wopos//'_Y',wopos//'_Z',GID_NULL)
       if(ndime==2) then
          if(kbopo==0) then
             do ipoin=1,npoin
                CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,bridge(1,ipoin), &
                     &                  bridge(2,ipoin),0.0_rp)
             end do
          else
             do ipoin = 1,npoin
                ibopo=lpoty(ipoin)
                if(ibopo/=0) then
                   CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,bridge(1,ibopo), &
                        &                  bridge(2,ibopo),0.0_rp)
                else
                   CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,dummr,dummr,dummr)
                end if
             end do
          end if
       else
          if(kbopo==0) then
             do ipoin=1,npoin
                CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,bridge(1,ipoin), &
                     &                  bridge(2,ipoin),bridge(3,ipoin))
             end do
          else
             do ipoin=1,npoin
                ibopo=lpoty(ipoin)
                if(ibopo/=0) then
                   CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,bridge(1,ibopo), &
                        &                  bridge(2,ibopo),bridge(3,ibopo))
                else
                   CALL GID_fWRITEVECTOR(pos%gid_luou,ipoin,dummr,dummr,dummr)
                end if
             end do
          end if
       end if
       CALL GID_fENDRESULT(pos%gid_luou)
       CALL GID_fFLUSHPOSTFILE(pos%gid_luou)
#else
       call enable_gidpost_error_message       
#endif
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e14.6,' ',a,' OnNodes')
3   format('ComponentNames ',a,',',a,',',a)
4   format(i7, 3(1x,e16.8E3))

    ! Femview formats
10  format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20  format(1x,i2,2x,a5,3x,3i5)
30  format(1x,i2,2x,a8,5i5)
50  format(1x,i2,i5,3e12.5)
70  format(1x,i2)
80  format(1x,i2,2x,a8,5i5,a8)

  end subroutine posvec

  !==================================================================================
  subroutine poisca(pos,bridge,wopos,istep,ttime,nren)
    !-----------------------------------------------------------------------
    !
    ! Write a scalar in postprocess file
    !
    !-----------------------------------------------------------------------
    implicit none
    type(post_file_t), intent(in)  :: pos
    character(*)   , intent(in)  :: wopos
    integer(ip)    , intent(in)  :: bridge(:)
    integer(ip)    , intent(in)  :: istep
    real(rp)       , intent(in)  :: ttime
    type(renum_t)    , intent(in), optional   :: nren
    ! Locals
    integer(ip)  :: npoin,ipoin
    character(8) :: state

    npoin=size(bridge)
    state='ANALYSIS'

    if(pos%form==1) then
       write(pos%luou,2) wopos,state,ttime,'Scalar'
       write(pos%luou,3) wopos
       write(pos%luou,1) 'Values'
       if(present(nren))then
          do ipoin=1,npoin
             write(pos%luou,4) nren%lperm(ipoin),bridge(ipoin)
          end do
       else
          do ipoin=1,npoin
             write(pos%luou,4) ipoin,bridge(ipoin)
          end do
       end if
       write(pos%luou,1) 'End Values'
       flush(pos%luou)
    else if(pos%form==2) then
       write(pos%luou,10) 100,'C','Step  ',ttime,1,istep
       write(pos%luou,20) -4,wopos,1,1,0
       write(pos%luou,30) -5,wopos,1,1
       do ipoin=1,npoin
          write(pos%luou,40) -1,ipoin,bridge(ipoin)
       end do
       write(pos%luou,50) -3
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
#ifdef ENABLE_GIDPOST
       CALL GID_fBEGINSCALARRESULT(pos%gid_luou,wopos,state,ttime,GID_ONNODES, &
            & GID_NULL,GID_NULL,GID_NULL)
       do ipoin=1,npoin
          CALL GID_fWRITESCALAR(pos%gid_luou,ipoin,real(bridge(ipoin)))
       end do
       CALL GID_fENDRESULT(pos%gid_luou)
       CALL GID_fFLUSHPOSTFILE(pos%gid_luou)
#else
       call enable_gidpost_error_message
#endif
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e14.6,' ',a,' OnNodes')
3   format('ComponentNames ',a)
4   format(i7, 3(1x,i7))

    ! Femview formats
10  format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20  format(1x,i2,2x,a5,3x,3i5)
30  format(1x,i2,2x,a5,3x,2i5)
40  format(1x,i2,i5,i7)
50  format(1x,i2)

  end subroutine poisca

  !==================================================================================
  subroutine postpro_gp_init(pos,ngaus,nnode,ndime)
    implicit none
    type(post_file_t), intent(in) :: pos
    integer(ip)    , intent(in) :: ngaus, nnode, ndime

    character(13)             :: elemt

    if(ndime==2) then
       if(nnode==3.or.nnode==6.or.nnode==7) then
          elemt='Triangle' 
       else
          elemt='Quadrilateral'
       end if
    else
       if(nnode==4.or.nnode==10) then 
          elemt='Tetrahedra'
       else if(nnode==8.or.nnode==20.or.nnode==27) then 
          elemt='Hexahedra'
       else if(nnode==6.or.nnode==15) then 
          elemt='Prism'
       end if
    end if
    write(pos%luou,1) 'GaussPoints GP Elemtype '//trim(elemt)
    write(pos%luou,5) 'Number of Gauss Points: ',ngaus
    write(pos%luou,1) 'Natural Coordinates: Internal'
    write(pos%luou,1) 'End GaussPoints'

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e13.6,' ',a,' OnGaussPoints ',a)
3   format('ComponentNames ',a)
4   format(i7, 3(1x,e16.8E3))
5   format(a,1x,i2)
6   format(6(1x,e16.8E3))
7   format(i7,$)

  end subroutine postpro_gp_init

  !=================================================================================================
  subroutine pgpint(pos,ndime,nnode,bridge,wopos,istep,ttime,eren)
    !-----------------------------------------------------------------------
    !
    ! Write a constant by element integer function in postprocess file.
    ! Only for meshes with nelty=1
    !
    !-----------------------------------------------------------------------
    implicit none
    type(post_file_t), intent(in)  :: pos
    character(*)   , intent(in)  :: wopos
    integer(ip)    , intent(in)  :: ndime, nnode
    integer(ip)    , intent(in)  :: bridge(:)
    integer(ip)    , intent(in)  :: istep
    real(rp)       , intent(in)  :: ttime
    integer(ip)    , intent(in), optional   :: eren(:)
    !type(renum_t)    , intent(in), optional   :: eren

    character(8)   :: state
    character(13)  :: elemt
    integer(ip)    :: ielem


    state='ANALYSIS'

    if(pos%form==1) then
       write(pos%luou,2) wopos,state,ttime,'Scalar','GP'
       write(pos%luou,3) wopos
       write(pos%luou,1) 'Values'
       if(present(eren)) then
          do ielem=1,size(bridge)
             write(pos%luou,7) eren(ielem)
             write(pos%luou,6) real(bridge(ielem),rp)
          end do
       else
          do ielem=1,size(bridge)
             write(pos%luou,7) ielem
             write(pos%luou,6) real(bridge(ielem),rp)
          end do
       end if
       write(pos%luou,1) 'End Values'
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
       !call GIDPOST LIB 
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e14.6,' ',a,' OnGaussPoints ',a)
3   format('ComponentNames ',a)
4   format(i7, 3(1x,e16.8E3))
5   format(a,1x,i2)
6   format(1x,e16.8E3)
7   format(i7,$)

  end subroutine pgpint

  !=================================================================================================
  subroutine pgpr3p(bridge,pos,ndime,wopos,istep,ttime,eren)
    !-----------------------------------------------------------------------
    !
    ! Write a constant by element integer function in postprocess file.
    ! Only for meshes with nelty=1
    !
    !-----------------------------------------------------------------------
    implicit none
    type(post_file_t), intent(in)  :: pos
    character(*)   , intent(in)  :: wopos
    real(rp)       , intent(in)  :: bridge(:,:,:)
    integer(ip)    , intent(in)  :: istep,ndime
    real(rp)       , intent(in)  :: ttime
    integer(ip)    , intent(in), optional   :: eren(:)

    character(8)              :: state
    integer(ip)               :: icomp,ncomp
    integer(ip)               :: igaus,ielem,ngaus,nelem
    character(13)             :: elemt

    state='ANALYSIS'

    if(pos%form==1) then
       ncomp=size(bridge,1)
       ngaus=size(bridge,2)
       nelem=size(bridge,3)
       if(ncomp==ndime) then
          write(pos%luou,2) wopos,state,ttime,'Vector','GP'
          write(pos%luou,3) wopos//'_X,'//wopos//'_Y,'//wopos//'_Z'
          write(pos%luou,1) 'Values'
          do ielem=1,nelem
             write(pos%luou,7) ielem
             do igaus=1,ngaus
                write(pos%luou,6) (bridge(icomp,igaus,ielem),icomp=1,ncomp)
             end do
          end do
       else
          write(pos%luou,2) wopos,state,ttime,'Matrix','GP'
          if(ndime==2) then
             write(pos%luou,3) wopos//'_XX,'//wopos//'_YY,'//wopos//'_ZZ,'//wopos//'_XY'
          else
             write(pos%luou,3) wopos//'_XX,'//wopos//'_YY,'//wopos//'_ZZ,'//wopos//'_XY,'&
                  //wopos//'_YZ,'//wopos//'_XZ'
          end if
          write(pos%luou,1) 'Values'
          if(ndime==2) then
             do ielem=1,nelem
                write(pos%luou,7) ielem                    
                do igaus=1,ngaus
                   write(pos%luou,6) &
                        bridge(1,igaus,ielem),&
                        bridge(2,igaus,ielem),&
                        bridge(4,igaus,ielem),&
                        bridge(3,igaus,ielem)
                end do
             end do
          else
             do ielem=1,nelem
                write(pos%luou,7) ielem                    
                do igaus=1,ngaus
                   write(pos%luou,6) &
                        bridge(1,igaus,ielem),&
                        bridge(2,igaus,ielem),&
                        bridge(4,igaus,ielem),&
                        bridge(3,igaus,ielem),&
                        bridge(6,igaus,ielem),&
                        bridge(5,igaus,ielem)
                end do
             end do
          end if
       end if
       write(pos%luou,1) 'End Values'
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
       !call GIDPOST LIB 
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e13.6,' ',a,' OnGaussPoints ',a)
3   format('ComponentNames ',a)
4   format(i7, 3(1x,e16.8E3))
5   format(a,1x,i2)
6   format(6(1x,e16.8E3))
7   format(i7,$)

  end subroutine pgpr3p

  !=================================================================================================
  subroutine pgpr1p(bridge,pos,wopos,istep,ttime,eren)
    !-----------------------------------------------------------------------
    !
    ! Write a constant by element integer function in postprocess file.
    ! Only for meshes with nelty=1
    !
    !-----------------------------------------------------------------------
    implicit none
    character(*)   , intent(in)  :: wopos
    type(post_file_t), intent(in)  :: pos
    real(rp)       , intent(in)  :: bridge(:,:)
    integer(ip)    , intent(in)  :: istep
    real(rp)       , intent(in)  :: ttime
    integer(ip)    , intent(in), optional   :: eren(:)

    character(8)              :: state
    integer(ip)               :: igaus,ielem,ngaus,nelem

    state='ANALYSIS'

    ngaus = size(bridge,2)
    if(pos%form==1) then

       write(pos%luou,2) wopos,state,ttime,'Scalar','GP'
       write(pos%luou,3) wopos
       write(pos%luou,1) 'Values'
       ngaus=size(bridge,1)
       nelem=size(bridge,2)
       do ielem=1,nelem
          write(pos%luou,7) ielem
          do igaus=1,ngaus
             write(pos%luou,6) bridge(igaus,ielem)
          end do
       end do
       write(pos%luou,1) 'End Values'
       flush(pos%luou)
    else if(pos%form==3.or.pos%form==4) then
       !call GIDPOST LIB 
    end if

    ! GiD formats
1   format(a)
2   format('Result ',a,' ',a,' ',e13.6,' ',a,' OnGaussPoints ',a)
3   format('ComponentNames ',a)
4   format(i7, 3(1x,e16.8E3))
5   format(a,1x,i2)
6   format(6(1x,e16.8E3))
7   format(i7,$)

  end subroutine pgpr1p

  !   !=================================================================================================
  !   subroutine pgpr1p(bridge,wopos,itste,ttime)
  !     !-----------------------------------------------------------------------
  !     !
  !     ! Write a matrix at Gauss points in postprocess file
  !     !
  !     !-----------------------------------------------------------------------
  !     implicit none
  !     character(*), intent(in)  :: wopos
  !     type(r1p_t),    intent(in)  :: bridge(:) 
  !     integer(ip) , intent(in)  :: itste
  !     real(rp)    , intent(in)  :: ttime
  !     character(8)              :: state
  !     character(4)              :: CNUL=CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
  !     integer(ip),  save        :: ipass=0
  !     integer(ip)               :: iesta,iesto,ielty,icomp,ncomp
  !     character(13)             :: elemt

  !     state='ANALYSIS'

  !     if(pos%form==1) then

  !        if(ipagp==0) then
  !           ipagp=1
  !           do ielty=1,nelty
  !              if(ndime==2) then
  !                 if(nnode==3.or.nnode==6.or.nnode==7) then
  !                    elemt='Triangle' 
  !                 else
  !                    elemt='Quadrilateral'
  !                 end if
  !              else
  !                 if(nnode==4.or.nnode==10) then 
  !                    elemt='Tetrahedra'
  !                 else if(nnode==8.or.nnode==20.or.nnode==27) then 
  !                    elemt='Hexahedra'
  !                 else if(nnode==6.or.nnode==15) then 
  !                    elemt='Prism'
  !                 end if
  !              end if
  !              write(pos%luou,1) 'GaussPoints '//'GP_'//trim(intost)&
  !                 &             //' Elemtype '//trim(elemt)
  !              write(pos%luou,5) 'Number of Gauss Points: ',ngaus
  !              write(pos%luou,1) 'Natural Coordinates: Internal'
  !              write(pos%luou,1) 'End GaussPoints'
  !           end do
  !        end if
  !        do ielty=1,nelty
  !           write(pos%luou,2) wopos,state,ttime,'Scalar','GP_'//trim(intost)
  !           write(pos%luou,3) wopos
  !           write(pos%luou,1) 'Values'
  !           do ielem=1,nelem
  !              if(ltype(ielem)==ielty) then
  !                 write(pos%luou,7) ielem
  !                 do igaus=1,ngaus
  !                    write(pos%luou,6) bridge(ielem)%a(igaus)
  !                 end do
  !              end if
  !           end do
  !           write(pos%luou,1) 'End Values'
  !        end do
  !        call flush(pos%luou)
  !     else if(pos%form==3.or.pos%form==4) then
  !        !call GIDPOST LIB 
  !     end if

  !     ! GiD formats
  ! 1   format(a)
  ! 2   format('Result ',a,' ',a,' ',e13.6,' ',a,' OnGaussPoints ',a)
  ! 3   format('ComponentNames ',a)
  ! 4   format(i7, 3(1x,e16.8E3))
  ! 5   format(a,1x,i2)
  ! 6   format(6(1x,e16.8E3))
  ! 7   format(i7,$)

  !     ! Femview formats
  ! 10  format(1x,i4,a1,a6,e12.5,32x,i2,i5)
  ! 20  format(1x,i2,2x,a5,3x,3i5)
  ! 30  format(1x,i2,2x,a5,3x,2i5)
  ! 40  format(1x,i2,i5,e12.5)
  ! 50  format(1x,i2)

  !   end subroutine pgpr1p

  !=================================================================================================
  subroutine poserr(msg)
    implicit none
    character(*),intent(in) :: msg
    write(6,'(a)') 
  end subroutine poserr

end module postpro_names

