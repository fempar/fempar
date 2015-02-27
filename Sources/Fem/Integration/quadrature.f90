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
module quadrature_names
  use types
  use memor
  implicit none
  private

  type quadrature
     integer(ip)           :: &
        ndime,                &    ! Number of space dimensions
        ngaus                      ! Number of integration points
     real(rp), allocatable :: &
        pos(:,:),             &    ! Quadrature points position
        weight(:)                  ! Quadrature points weight
  end type quadrature

  ! Types
  public :: quadrature

  ! Functions
  public :: quadrature_create, quadrature_free

contains

  subroutine quadrature_create(nrule,ndime,ngaus,rule)
    !-----------------------------------------------------------------------
    ! 
    ! This routine creates the integration rule
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)     , intent(in)  :: nrule,ndime,ngaus
    type(quadrature), intent(out) :: rule

    rule%ndime=ndime
    rule%ngaus=ngaus
    call memalloc(ndime,ngaus,rule%pos,__FILE__,__LINE__)
    call memalloc(ngaus,rule%weight,__FILE__,__LINE__)
    rule%pos=0.0_rp
    rule%weight=0.0_rp

    select case (nrule)

    case(ruqope_id) 
       call ruqope(ndime,ngaus,rule%pos,rule%weight)
    case(ruqclo_id) 
       call ruqclo(ndime,ngaus,rule%pos,rule%weight)
    case(rutope_id) 
       call rutope(ndime,ngaus,rule%pos,rule%weight)
    case(rutclo_id) 
       call rutclo(ndime,ngaus,rule%pos,rule%weight)
    case(rupope_id) 
       call rupope(ndime,ngaus,rule%pos,rule%weight)
    case(rupclo_id) 
       call rupclo(ndime,ngaus,rule%pos,rule%weight)
    case (ruqope_tp_id)
       call ruqope(    1,ngaus,rule%pos,rule%weight)

    end select

  end subroutine quadrature_create

  subroutine quadrature_free(rule)
    !-----------------------------------------------------------------------
    ! 
    ! This routine frees the integration rule
    !
    !-----------------------------------------------------------------------
    implicit none
    type(quadrature), intent(inout) :: rule

    call memfree(rule%pos,__FILE__,__LINE__)
    call memfree(rule%weight,__FILE__,__LINE__)

  end subroutine quadrature_free

  !-----------------------------------------------------------------------
  subroutine ruqope(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    !
    !     This routine sets up the integration constants of open
    !     integration rules for brick elements:
    ! 
    !          NDIME = 1             NDIME = 2             NDIME = 3
    ! 
    !      NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
    !      -----  ----------     -----  ----------     -----  ----------
    !        1      q1           1 x 1     q1          1x1x1     q1	
    !        2      q3           2 x 2     q3          2x2x2     q3   
    !        3      q5           3 x 3     q5          3x3x3     q5
    !        4      q7           4 x 4     q7          4x4x4     q7
    !        5      q9           5 x 5     q9          5x5x5     q9
    !        6      q11          6 x 6     q11         6x6x6     q11
    !        7      q13          7 x 7     q13         7x7x7     q13
    !        8      q15          8 x 8     q15         8x8x8     q15
    !       16      q31         16 x 16    q31        16x16x16   q31
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: posgl(16),weigl(16)
    integer(ip)              :: nlocs,igaus,ilocs,jlocs,klocs

    if(ndime==1) then
       nlocs=ngaus
    else if(ndime==2) then
       nlocs=nint(sqrt(real(ngaus,rp)))
    else
       nlocs=nint(real(ngaus,rp)**(1.0_rp/3.0_rp))
    end if

    if(nlocs==1) then
       posgl(1)=0.0_rp
       weigl(1)=2.0_rp
    else if(nlocs==2) then
       posgl(1)=-0.577350269189626_rp
       posgl(2)= 0.577350269189626_rp
       weigl(1)= 1.0_rp
       weigl(2)= 1.0_rp
    else if(nlocs==3) then
       posgl(1)=-0.774596669241483_rp
       posgl(2)= 0.0_rp
       posgl(3)= 0.774596669241483_rp
       weigl(1)= 0.555555555555556_rp
       weigl(2)= 0.888888888888889_rp
       weigl(3)= 0.555555555555556_rp
    else if(nlocs==4)  then
       posgl(1)=-0.861136311594053_rp
       posgl(2)=-0.339981043584856_rp
       posgl(3)= 0.339981043584856_rp
       posgl(4)= 0.861136311594053_rp
       weigl(1)= 0.347854845137454_rp
       weigl(2)= 0.652145154862546_rp
       weigl(3)= 0.652145154862546_rp
       weigl(4)= 0.347854845137454_rp
    else if(nlocs==5)  then
       posgl(1) = -0.906179845938664_rp
       posgl(2) = -0.538469310105683_rp
       posgl(3) =  0.0_rp
       posgl(4) =  0.538469310105683_rp
       posgl(5) =  0.906179845938664_rp
       weigl(1) =  0.236926885056189_rp
       weigl(2) =  0.478628670499366_rp
       weigl(3) =  0.568888888888889_rp
       weigl(4) =  0.478628670499366_rp
       weigl(5) =  0.236926885056189_rp
    else if(nlocs==6)  then
       posgl(1) = -0.932469514203152_rp
       posgl(2) = -0.661209386466265_rp
       posgl(3) = -0.238619186083197_rp
       posgl(4) =  0.238619186083197_rp
       posgl(5) =  0.661209386466265_rp
       posgl(6) =  0.932469514203152_rp
       weigl(1) =  0.171324492379170_rp
       weigl(2) =  0.360761573048139_rp
       weigl(3) =  0.467913934572691_rp
       weigl(4) =  0.467913934572691_rp
       weigl(5) =  0.360761573048139_rp
       weigl(6) =  0.171324492379170_rp
    else if(nlocs==7)  then
       posgl(1) = -0.949107912342759_rp
       posgl(2) = -0.741531185599394_rp
       posgl(3) = -0.405845151377397_rp
       posgl(4) =  0.0_rp
       posgl(5) =  0.405845151377397_rp
       posgl(6) =  0.741531185599394_rp
       posgl(7) =  0.949107912342759_rp
       weigl(1) =  0.129484966168870_rp
       weigl(2) =  0.279705391489277_rp
       weigl(3) =  0.381830050505119_rp
       weigl(4) =  0.417959183673469_rp
       weigl(5) =  0.381830050505119_rp
       weigl(6) =  0.279705391489277_rp
       weigl(7) =  0.129484966168870_rp
    else if(nlocs==8)  then
       posgl(1) = -0.960289856497536_rp
       posgl(2) = -0.796666477413627_rp
       posgl(3) = -0.525532409916329_rp
       posgl(4) = -0.183434642495650_rp
       posgl(5) =  0.183434642495650_rp
       posgl(6) =  0.525532409916329_rp
       posgl(7) =  0.796666477413627_rp
       posgl(8) =  0.960289856497536_rp

       weigl(1) =  0.101228536290376_rp
       weigl(2) =  0.222381034453374_rp
       weigl(3) =  0.313706645877887_rp
       weigl(4) =  0.362683783378362_rp
       weigl(5) =  0.362683783378362_rp
       weigl(6) =  0.313706645877887_rp
       weigl(7) =  0.222381034453374_rp
       weigl(8) =  0.101228536290376_rp
    else if(nlocs==16)  then
       posgl( 1) =-0.98940093499165_rp
       posgl( 2) =-0.94457502307323_rp
       posgl( 3) =-0.86563120238783_rp
       posgl( 4) =-0.75540440835500_rp
       posgl( 5) =-0.61787624440264_rp
       posgl( 6) =-0.45801677765723_rp
       posgl( 7) =-0.28160355077926_rp
       posgl( 8) =-0.09501250983764_rp
       posgl( 9) = 0.09501250983764_rp
       posgl(10) = 0.28160355077926_rp
       posgl(11) = 0.45801677765723_rp
       posgl(12) = 0.61787624440264_rp
       posgl(13) = 0.75540440835500_rp
       posgl(14) = 0.86563120238783_rp
       posgl(15) = 0.94457502307323_rp
       posgl(16) = 0.98940093499165_rp

       weigl( 1) =  0.02715245941175_rp
       weigl( 2) =  0.06225352393865_rp
       weigl( 3) =  0.09515851168249_rp
       weigl( 4) =  0.12462897125553_rp
       weigl( 5) =  0.14959598881658_rp
       weigl( 6) =  0.16915651939500_rp
       weigl( 7) =  0.18260341504492_rp
       weigl( 8) =  0.18945061045507_rp
       weigl( 9) =  0.18945061045507_rp
       weigl(10) =  0.18260341504492_rp
       weigl(11) =  0.16915651939500_rp
       weigl(12) =  0.14959598881658_rp
       weigl(13) =  0.12462897125553_rp
       weigl(14) =  0.09515851168249_rp
       weigl(15) =  0.06225352393865_rp
       weigl(16) =  0.02715245941175_rp
    else
       write(*,*) __FILE__,__LINE__,'ERROR:: Quadrature not defined'
       stop
    end if

    if(ndime==1) then
       igaus=0
       do ilocs=1,nlocs
          igaus=igaus+1
          weigp(  igaus)=weigl(ilocs)
          posgp(1,igaus)=posgl(ilocs)
       end do
    else if(ndime==2) then
       igaus=0
       do jlocs=1,nlocs
          do ilocs=1,nlocs
             igaus=igaus+1
             weigp(  igaus)=weigl(ilocs)*weigl(jlocs)
             posgp(1,igaus)=posgl(ilocs)
             posgp(2,igaus)=posgl(jlocs)
          end do
       end do
    else if(ndime==3) then
       igaus=0
       do klocs=1,nlocs
          do jlocs=1,nlocs
             do ilocs=1,nlocs
                igaus=igaus+1
                weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
                posgp(1,igaus)=posgl(ilocs)
                posgp(2,igaus)=posgl(jlocs)
                posgp(3,igaus)=posgl(klocs)
             end do
          end do
       end do
    end if

  end subroutine ruqope

  !-----------------------------------------------------------------------
  subroutine ruqclo(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    !
    !    This routine sets up the integration constants of closed 
    !    integration rules for brick elements:
    !
    !         NDIME = 1             NDIME = 2             NDIME = 3
    !
    !     NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
    !     -----  ----------     -----  ----------     -----  ----------  
    !       2      q1           2 x 2     q1          2x2x2     q1   
    !       3      q3           3 x 3     q3          3x3x3     q3
    !       4      q3           4 x 4     q3          4x4x4     q3
    !                             5       p2?           9       p2?
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: posgl(4), weigl(4)
    integer(ip)              :: inoga(27),igaus,nlocs,ilocs,jlocs,klocs

    ! 2D 8-point integration
    if(ndime==2.and.ngaus==8) then
       do igaus=5,8
          weigp(igaus-4)=-1.0_rp/3.0_rp
          weigp(igaus  )= 4.0_rp/3.0_rp
       end do
       posgp(1,5)= 0.0_rp 
       posgp(1,6)= 1.0_rp
       posgp(1,7)= 0.0_rp
       posgp(1,8)=-1.0_rp
       posgp(2,5)=-1.0_rp
       posgp(2,6)= 0.0_rp
       posgp(2,7)= 1.0_rp
       posgp(2,8)= 0.0_rp
       return
    end if

    ! 3D 20-point integration
    if(ndime==3.and.ngaus==20) then
       nlocs=nint(dfloat(ngaus)**(1.0_rp/3.0_rp))
       do igaus=1,8
          weigp(igaus)=-1.0_rp
       end do
       do igaus=9,20
          weigp(igaus)= 4.0_rp/3.0_rp
       end do
       posgl(1)=-1.0_rp
       posgl(2)= 0.0_rp
       posgl(3)= 1.0_rp
       call chaord(inoga,27)
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             do klocs=1,nlocs
                igaus=igaus+1
                if(inoga(igaus).le.20) then
                   posgp(1,inoga(igaus))=posgl(ilocs)
                   posgp(2,inoga(igaus))=posgl(jlocs)
                   posgp(3,inoga(igaus))=posgl(klocs)
                end if
             end do
          end do
       end do
       return
    end if

    ! Rules obtained from one-dimensional integration
    if(ndime==1) then
       nlocs=ngaus
    else if(ndime==2) then
       nlocs=nint(sqrt(dfloat(ngaus)))
    else
       nlocs=nint(dfloat(ngaus)**(1.0_rp/3.0_rp))
    end if

    call chaord(inoga,nlocs**ndime)

    if(nlocs==2) then
       posgl(1)=-1.0_rp
       posgl(2)= 1.0_rp
       weigl(1)= 1.0_rp
       weigl(2)= 1.0_rp
    else if(nlocs==3) then
       posgl(1)=-1.0_rp
       posgl(2)= 0.0_rp
       posgl(3)= 1.0_rp
       weigl(1)= 1.0_rp/3.0_rp
       weigl(2)= 4.0_rp/3.0_rp
       weigl(3)= 1.0_rp/3.0_rp
    else if(nlocs==4) then
       posgl(1)=-1.0_rp
       posgl(2)=-1.0_rp/3.0_rp
       posgl(3)= 1.0_rp/3.0_rp
       posgl(4)= 1.0_rp
       weigl(1)= 1.0_rp/4.0_rp
       weigl(2)= 3.0_rp/4.0_rp
       weigl(3)= 3.0_rp/4.0_rp
       weigl(4)= 1.0_rp/4.0_rp
    end if

    if(ndime==1) then
       igaus=0
       do ilocs=1,nlocs
          igaus=igaus+1
          weigp(  igaus)=weigl(ilocs)
          posgp(1,igaus)=posgl(ilocs)
       end do
    else if(ndime==2) then
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             igaus=igaus+1
             weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)
             posgp(1,inoga(igaus))=posgl(ilocs)
             posgp(2,inoga(igaus))=posgl(jlocs)
          end do
       end do
       if(ngaus==5) then                                    !  4       3
          do igaus=1,4                                      !
             weigp(igaus)=1.0_rp/3.0_rp                     !      5
          end do                                            !
          weigp(  5)=8.0_rp/3.0_rp                          !  1       2
          posgp(1,5)=0.0_rp
          posgp(2,5)=0.0_rp
       end if
    else if(ndime==3) then
       igaus=0
       do ilocs=1,nlocs
          do jlocs=1,nlocs
             do klocs=1,nlocs
                igaus=igaus+1
                weigp(  inoga(igaus))=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
                posgp(1,inoga(igaus))=posgl(ilocs)
                posgp(2,inoga(igaus))=posgl(jlocs)
                posgp(3,inoga(igaus))=posgl(klocs)
             end do
          end do
       end do
       if(ngaus==9) then                                 !    8------7
          do igaus=1,8                                   !   /|     /|
             weigp(igaus)=weigp(igaus)/3.0_rp            !  5------6 |
          end do                                         !  | | 9  | |
          weigp(  9)=16.0_rp/3.0_rp                      !  | 4----|-3
          posgp(1,9)=0.0_rp                              !  |/     |/
          posgp(2,9)=0.0_rp                              !  1------2
          posgp(3,9)=0.0_rp
       end if
    end if

  end subroutine ruqclo

  !-----------------------------------------------------------------------
  subroutine rutope(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    ! 
    !     This routine sets up the integration constants of open rules for
    !     triangles and tetrahedra
    ! 
    !             NDIME = 2             NDIME = 3
    ! 
    !          NGAUS  EXACT POL.     NGAUS  EXACT POL. 
    !          -----  ----------     -----  ----------
    !            1       p1            1       p1
    !            3       p2            4       p2
    !            4       p3            5       p3
    !            6       p4           11       p4
    !            7       p5           14       p5
    !           13       p7
    !           19       p9
    !           28       p11
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: a,b,c,d,e,f,g,h,p,q,r,s,t,u,v,w,x,y,z
    real(rp)                 :: w1,w2,w3,w4,w5,w6,w7,w8
    real(rp)                 :: ex1,et1,ez1,ex2,et2,ez2
    integer(ip)              :: istop

    istop=0

    ! Line integral (the same as for brick elements)
    if(ndime==1) then
       if(ngaus==1) then
          posgp(1,1)=0.0_rp
          weigp(  1)=2.0_rp
       else if(ngaus==2) then
          posgp(1,1)=-0.577350269189626_rp
          posgp(1,2)= 0.577350269189626_rp
          weigp(  1)= 1.0_rp
          weigp(  2)= 1.0_rp
       else if(ngaus==3) then
          posgp(1,1)=-0.774596669241483_rp
          posgp(1,2)= 0.0_rp
          posgp(1,3)= 0.774596669241483_rp
          weigp(  1)= 0.555555555555556_rp
          weigp(  2)= 0.888888888888889_rp
          weigp(  3)= 0.555555555555556_rp
       else
          istop=1
       end if

       ! Area integral (triangles)
    else if(ndime==2) then
       if(ngaus==1) then
          posgp(1,1)= 1.0_rp/3.0_rp
          posgp(2,1)= 1.0_rp/3.0_rp
          weigp(  1)= 1.0_rp/2.0_rp
       else if(ngaus==3) then
          posgp(1,2)= 2.0_rp/3.0_rp
          posgp(2,2)= 1.0_rp/6.0_rp
          posgp(1,3)= 1.0_rp/6.0_rp
          posgp(2,3)= 2.0_rp/3.0_rp
          posgp(1,1)= 1.0_rp/6.0_rp
          posgp(2,1)= 1.0_rp/6.0_rp
          weigp(  2)= 1.0_rp/6.0_rp
          weigp(  3)= 1.0_rp/6.0_rp
          weigp(  1)= 1.0_rp/6.0_rp
       else if(ngaus==4) then
          posgp(1,1)= 1.0_rp/3.0_rp
          posgp(2,1)= 1.0_rp/3.0_rp
          posgp(1,2)= 1.0_rp/5.0_rp
          posgp(2,2)= 1.0_rp/5.0_rp
          posgp(1,3)= 3.0_rp/5.0_rp
          posgp(2,3)= 1.0_rp/5.0_rp
          posgp(1,4)= 1.0_rp/5.0_rp
          posgp(2,4)= 3.0_rp/5.0_rp
          weigp(  1)=-27.0_rp/96.0_rp
          weigp(  2)= 25.0_rp/96.0_rp
          weigp(  3)= 25.0_rp/96.0_rp
          weigp(  4)= 25.0_rp/96.0_rp 
       else if(ngaus==6) then
          ex1 = 0.816847572980459_rp
          et1 = 0.091576213509771_rp
          ez1 = 0.091576213509771_rp
          ex2 = 0.108103018168070_rp
          et2 = 0.445948490915965_rp
          ez2 = 0.445948490915965_rp
          posgp(1,3)= ex1
          posgp(2,3)= et1
          posgp(1,1)= et1
          posgp(2,1)= ez1
          posgp(1,6)= ez1
          posgp(2,6)= ex1
          posgp(1,4)= ex2
          posgp(2,4)= et2
          posgp(1,5)= et2
          posgp(2,5)= ez2
          posgp(1,2)= ez2
          posgp(2,2)= ex2
          a = 0.054975870996713638_rp
          b = 0.1116907969117165_rp    
          weigp(3)  = a
          weigp(1)  = a
          weigp(6)  = a
          weigp(4)  = b
          weigp(5)  = b
          weigp(2)  = b
       else if(ngaus==7) then
          a = 1.0_rp / 3.0_rp
          b = ( 9.0_rp + 2.0_rp * sqrt ( 15.0_rp ) ) / 21.0_rp
          c = ( 6.0_rp -          sqrt ( 15.0_rp ) ) / 21.0_rp
          d = ( 9.0_rp - 2.0_rp * sqrt ( 15.0_rp ) ) / 21.0_rp
          e = ( 6.0_rp +          sqrt ( 15.0_rp ) ) / 21.0_rp
          w1 = 0.1125_rp
          w2 = ( 155.0_rp - sqrt ( 15.0_rp ) ) / 2400.0_rp
          w3 = ( 155.0_rp + sqrt ( 15.0_rp ) ) / 2400.0_rp
          posgp(1,1)= a
          posgp(2,1)= a
          posgp(1,2)= b
          posgp(2,2)= c
          posgp(1,3)= c
          posgp(2,3)= b
          posgp(1,4)= c
          posgp(2,4)= c
          posgp(1,5)= d
          posgp(2,5)= e
          posgp(1,6)= e
          posgp(2,6)= d
          posgp(1,7)= e
          posgp(2,7)= e
          weigp(  1)= w1
          weigp(  2)= w2
          weigp(  3)= w2
          weigp(  4)= w2
          weigp(  5)= w3
          weigp(  6)= w3
          weigp(  7)= w3
       else if(ngaus==13) then
          a = 0.333333333333333_rp
          b = 0.479308067841923_rp
          c = 0.869739794195568_rp
          d = 0.638444188569809_rp
          e = 0.260345966079038_rp
          f = 0.065130102902216_rp
          g = 0.312865496004875_rp
          h = 0.048690315425316_rp
          w1=-0.149570044467670_rp/2.0_rp
          w2= 0.175615257433204_rp/2.0_rp
          w3= 0.053347235608839_rp/2.0_rp
          w4= 0.077113760890257_rp/2.0_rp
          posgp(1, 1)= a
          posgp(2, 1)= a         
          posgp(1, 2)= e
          posgp(2, 2)= e
          posgp(1, 3)= b
          posgp(2, 3)= e        
          posgp(1, 4)= e
          posgp(2, 4)= b        
          posgp(1, 5)= f
          posgp(2, 5)= f        
          posgp(1, 6)= c
          posgp(2, 6)= f        
          posgp(1, 7)= f
          posgp(2, 7)= c        
          posgp(1, 8)= d
          posgp(2, 8)= g        
          posgp(1, 9)= d
          posgp(2, 9)= h        
          posgp(1,10)= g
          posgp(2,10)= d        
          posgp(1,11)= g
          posgp(2,11)= h        
          posgp(1,12)= h
          posgp(2,12)= d        
          posgp(1,13)= h
          posgp(2,13)= g
          weigp( 1) = w1
          weigp( 2) = w2
          weigp( 3) = w2
          weigp( 4) = w2
          weigp( 5) = w3
          weigp( 6) = w3
          weigp( 7) = w3
          weigp( 8) = w4
          weigp( 9) = w4
          weigp(10) = w4
          weigp(11) = w4
          weigp(12) = w4
          weigp(13) = w4
       else if(ngaus.eq.19) then
          a = 1.0_rp / 3.0_rp
          b = 0.02063496160252593_rp
          c = 0.4896825191987370_rp
          d = 0.1258208170141290_rp
          e = 0.4370895914929355_rp
          f = 0.6235929287619356_rp
          g = 0.1882035356190322_rp
          r = 0.9105409732110941_rp
          s = 0.04472951339445297_rp
          t = 0.7411985987844980_rp
          u = 0.03683841205473626_rp
          v = 0.22196288916076574_rp

          w1 = 0.09713579628279610_rp/2.0_rp
          w2 = 0.03133470022713983_rp/2.0_rp
          w3 = 0.07782754100477543_rp/2.0_rp
          w4 = 0.07964773892720910_rp/2.0_rp
          w5 = 0.02557767565869810_rp/2.0_rp
          w6 = 0.04328353937728940_rp/2.0_rp

          posgp(1, 1) = a 
          posgp(1, 2) = b
          posgp(1, 3) = c 
          posgp(1, 4) = c 
          posgp(1, 5) = d 
          posgp(1, 6) = e
          posgp(1, 7) = e
          posgp(1, 8) = f 
          posgp(1, 9) = g 
          posgp(1,10) = g
          posgp(1,11) = r
          posgp(1,12) = s
          posgp(1,13) = s
          posgp(1,14) = t
          posgp(1,15) = t
          posgp(1,16) = u
          posgp(1,17) = u
          posgp(1,18) = v
          posgp(1,19) = v

          posgp(2, 1) = a
          posgp(2, 2) = c
          posgp(2, 3) = b
          posgp(2, 4) = c
          posgp(2, 5) = e
          posgp(2, 6) = d
          posgp(2, 7) = e
          posgp(2, 8) = g
          posgp(2, 9) = f
          posgp(2,10) = g
          posgp(2,11) = s
          posgp(2,12) = r
          posgp(2,13) = s
          posgp(2,14) = u
          posgp(2,15) = v
          posgp(2,16) = t
          posgp(2,17) = v
          posgp(2,18) = t
          posgp(2,19) = u

          weigp( 1) = w1
          weigp( 2) = w2
          weigp( 3) = w2
          weigp( 4) = w2
          weigp( 5) = w3
          weigp( 6) = w3
          weigp( 7) = w3
          weigp( 8) = w4
          weigp( 9) = w4
          weigp(10) = w4
          weigp(11) = w5
          weigp(12) = w5
          weigp(13) = w5
          weigp(14) = w6
          weigp(15) = w6
          weigp(16) = w6
          weigp(17) = w6
          weigp(18) = w6
          weigp(19) = w6
       else if(ngaus.eq.28) then
          a = 1.0_rp / 3.0_rp
          b = 0.9480217181434233_rp
          c = 0.02598914092828833_rp
          d = 0.8114249947041546_rp
          e = 0.09428750264792270_rp
          f = 0.01072644996557060_rp
          g = 0.4946367750172147_rp
          p = 0.5853132347709715_rp
          q = 0.2073433826145142_rp
          r = 0.1221843885990187_rp
          s = 0.4389078057004907_rp
          t = 0.6779376548825902_rp
          u = 0.04484167758913055_rp
          v = 0.27722066752827925_rp
          w = 0.8588702812826364_rp
          x = 0.0_rp
          y = 0.1411297187173636_rp

          w1 = 0.08797730116222190_rp/2.0_rp
          w2 = 0.008744311553736190_rp/2.0_rp
          w3 = 0.03808157199393533_rp/2.0_rp
          w4 = 0.01885544805613125_rp/2.0_rp
          w5 = 0.07215969754474100_rp/2.0_rp
          w6 = 0.06932913870553720_rp/2.0_rp
          w7 = 0.04105631542928860_rp/2.0_rp
          w8 = 0.007362383783300573_rp/2.0_rp

          posgp(1, 1) = a  
          posgp(1, 2) = b  
          posgp(1, 3) = c  
          posgp(1, 4) = c  
          posgp(1, 5) = d  
          posgp(1, 6) = e  
          posgp(1, 7) = e  
          posgp(1, 8) = f  
          posgp(1, 9) = g
          posgp(1,10) = g
          posgp(1,11) = p 
          posgp(1,12) = q 
          posgp(1,13) = q
          posgp(1,14) = r 
          posgp(1,15) = s 
          posgp(1,16) = s 
          posgp(1,17) = t 
          posgp(1,18) = t 
          posgp(1,19) = u 
          posgp(1,20) = u 
          posgp(1,21) = v 
          posgp(1,22) = v 
          posgp(1,23) = w 
          posgp(1,24) = w 
          posgp(1,25) = x 
          posgp(1,26) = x 
          posgp(1,27) = y 
          posgp(1,28) = y

          posgp(2, 1) = a  
          posgp(2, 2) = c  
          posgp(2, 3) = b  
          posgp(2, 4) = c  
          posgp(2, 5) = e  
          posgp(2, 6) = d  
          posgp(2, 7) = e  
          posgp(2, 8) = g  
          posgp(2, 9) = f  
          posgp(2,10) = g  
          posgp(2,11) = q  
          posgp(2,12) = p  
          posgp(2,13) = q
          posgp(2,14) = s  
          posgp(2,15) = r  
          posgp(2,16) = s  
          posgp(2,17) = u  
          posgp(2,18) = v  
          posgp(2,19) = t  
          posgp(2,20) = v  
          posgp(2,21) = t  
          posgp(2,22) = u  
          posgp(2,23) = x  
          posgp(2,24) = y  
          posgp(2,25) = w  
          posgp(2,26) = y  
          posgp(2,27) = w  
          posgp(2,28) = x

          weigp( 1) = w1 
          weigp( 2) = w2
          weigp( 3) = w2 
          weigp( 4) = w2 
          weigp( 5) = w3 
          weigp( 6) = w3 
          weigp( 7) = w3 
          weigp( 8) = w4 
          weigp( 9) = w4 
          weigp(10) = w4 
          weigp(11) = w5 
          weigp(12) = w5 
          weigp(13) = w5
          weigp(14) = w6 
          weigp(15) = w6 
          weigp(16) = w6 
          weigp(17) = w7 
          weigp(18) = w7 
          weigp(19) = w7 
          weigp(20) = w7 
          weigp(21) = w7 
          weigp(22) = w7 
          weigp(23) = w8 
          weigp(24) = w8 
          weigp(25) = w8 
          weigp(26) = w8 
          weigp(27) = w8 
          weigp(28) = w8
       else
          istop=1
       end if

       ! Volume integral ( tetrahedra )
    else if(ndime==3) then
       if(ngaus==1) then
          posgp(1,1)= 1.0_rp/4.0_rp
          posgp(2,1)= 1.0_rp/4.0_rp
          posgp(3,1)= 1.0_rp/4.0_rp
          weigp(1)  = 1.0_rp/6.0_rp
       else if(ngaus==4) then
          a=0.5854101966249685_rp
          b=0.1381966011250105_rp
          posgp(1,1)= b
          posgp(2,1)= b
          posgp(3,1)= b
          posgp(1,2)= a
          posgp(2,2)= b
          posgp(3,2)= b
          posgp(1,3)= b
          posgp(2,3)= a
          posgp(3,3)= b
          posgp(1,4)= b
          posgp(2,4)= b
          posgp(3,4)= a
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 1.0_rp/24.0_rp
       else if(ngaus==5) then
          posgp(1,1)= 1.0_rp/4.0_rp
          posgp(2,1)= 1.0_rp/4.0_rp
          posgp(3,1)= 1.0_rp/4.0_rp
          posgp(1,2)= 1.0_rp/6.0_rp
          posgp(2,2)= 1.0_rp/6.0_rp
          posgp(3,2)= 1.0_rp/6.0_rp
          posgp(1,3)= 1.0_rp/2.0_rp
          posgp(2,3)= 1.0_rp/6.0_rp
          posgp(3,3)= 1.0_rp/6.0_rp
          posgp(1,4)= 1.0_rp/6.0_rp
          posgp(2,4)= 1.0_rp/2.0_rp
          posgp(3,4)= 1.0_rp/6.0_rp
          posgp(1,5)= 1.0_rp/6.0_rp
          posgp(2,5)= 1.0_rp/6.0_rp
          posgp(3,5)= 1.0_rp/2.0_rp
          weigp(  1)=-2.0_rp/15.0_rp
          weigp(  2)= 1.5_rp/20.0_rp
          weigp(  3)= 1.5_rp/20.0_rp
          weigp(  4)= 1.5_rp/20.0_rp
          weigp(  5)= 1.5_rp/20.0_rp
       else if(ngaus==11) then
          a=0.3994035761667992_rp
          b=0.1005964238332008_rp
          c=343.0_rp/7500.0_rp/6.0_rp
          d=56.0_rp/375.0_rp/6.0_rp
          posgp(1,1) = 1.0_rp/4.0_rp
          posgp(2,1) = 1.0_rp/4.0_rp
          posgp(3,1) = 1.0_rp/4.0_rp
          posgp(1,2) = 11.0_rp/14.0_rp
          posgp(2,2) = 1.0_rp/14.0_rp
          posgp(3,2) = 1.0_rp/14.0_rp
          posgp(1,3) = 1.0_rp/14.0_rp
          posgp(2,3) = 11.0_rp/14.0_rp
          posgp(3,3) = 1.0_rp/14.0_rp
          posgp(1,4) = 1.0_rp/14.0_rp
          posgp(2,4) = 1.0_rp/14.0_rp
          posgp(3,4) = 11.0_rp/14.0_rp
          posgp(1,5) = 1.0_rp/14.0_rp
          posgp(2,5) = 1.0_rp/14.0_rp
          posgp(3,5) = 1.0_rp/14.0_rp
          posgp(1,6) = a
          posgp(2,6) = a
          posgp(3,6) = b
          posgp(1,7) = a
          posgp(2,7) = b
          posgp(3,7) = a
          posgp(1,8) = a
          posgp(2,8) = b
          posgp(3,8) = b
          posgp(1,9) = b
          posgp(2,9) = a
          posgp(3,9) = a
          posgp(1,10)= b
          posgp(2,10)= a
          posgp(3,10)= b
          posgp(1,11)= b
          posgp(2,11)= b
          posgp(3,11)= a
          weigp(1)   =-148.0_rp/1875.0_rp/6.0_rp
          weigp(2)   = c
          weigp(3)   = c
          weigp(4)   = c
          weigp(5)   = c
          weigp(6)   = d
          weigp(7)   = d
          weigp(8)   = d
          weigp(9)   = d
          weigp(10)  = d
          weigp(11)  = d
       else if(ngaus==14) then
          a=0.0673422422100983_rp
          b=0.3108859192633005_rp
          c=0.7217942490673264_rp
          d=0.0927352503108912_rp
          e=0.4544962958743506_rp
          f=0.0455037041256494_rp
          p=0.1126879257180162_rp/6.0_rp
          q=0.0734930431163619_rp/6.0_rp
          r=0.0425460207770812_rp/6.0_rp
          posgp(1,1) = a
          posgp(2,1) = b
          posgp(3,1) = b
          posgp(1,2) = b
          posgp(2,2) = a
          posgp(3,2) = b
          posgp(1,3) = b
          posgp(2,3) = b
          posgp(3,3) = a
          posgp(1,4) = b
          posgp(2,4) = b
          posgp(3,4) = b
          posgp(1,5) = c
          posgp(2,5) = d
          posgp(3,5) = d
          posgp(1,6) = d
          posgp(2,6) = c
          posgp(3,6) = d
          posgp(1,7) = d
          posgp(2,7) = d
          posgp(3,7) = c
          posgp(1,8) = d
          posgp(2,8) = d
          posgp(3,8) = d
          posgp(1,9) = e
          posgp(2,9) = e
          posgp(3,9) = f
          posgp(1,10)= e
          posgp(2,10)= f
          posgp(3,10)= e
          posgp(1,11)= e
          posgp(2,11)= f
          posgp(3,11)= f
          posgp(1,12)= f
          posgp(2,12)= e
          posgp(3,12)= e
          posgp(1,13)= f
          posgp(2,13)= e
          posgp(3,13)= f
          posgp(1,14)= f
          posgp(2,14)= f
          posgp(3,14)= e
          weigp(1)   = p
          weigp(2)   = p
          weigp(3)   = p
          weigp(4)   = p
          weigp(5)   = q
          weigp(6)   = q
          weigp(7)   = q
          weigp(8)   = q
          weigp(9)   = r
          weigp(10)  = r
          weigp(11)  = r
          weigp(12)  = r
          weigp(13)  = r
          weigp(14)  = r
       else
          istop=1
       end if
    end if

  end subroutine rutope

  !-----------------------------------------------------------------------
  subroutine rutclo(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    ! 
    !     This routine sets up the integration constants of closed rules for
    !     triangles and tetrahedra 
    ! 
    !             NDIME = 2             NDIME = 3
    ! 
    !          NGAUS  EXACT POL.     NGAUS  EXACT POL. 
    !          -----  ----------     -----  ----------
    !            3       p1            4       p1
    !            4       p2            5       p2
    !            6       p1           10       p2
    !            7       p3           11       p3
    !           10       p1           15       p3
    !                                 20       p2 & x^4,...      
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    real(rp)                 :: w1,w2,w3
    integer(ip)              :: istop

    ! Line integral
    istop=0
    if(ndime==1) then
       if(ngaus==2) then
          posgp(1,1)=-1.0_rp
          posgp(1,2)= 1.0_rp
          weigp(1)= 1.0_rp
          weigp(2)= 1.0_rp
       else if(ngaus==3) then
          posgp(1,1)=-1.0_rp
          posgp(1,2)= 0.0_rp
          posgp(1,3)= 1.0_rp
          weigp(1)= 1.0_rp/3.0_rp
          weigp(2)= 4.0_rp/3.0_rp
          weigp(3)= 1.0_rp/3.0_rp
       else
          istop=1
       end if

       ! Area integral (triangles)
    else if(ndime==2) then
       if(ngaus==3) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          weigp(  1)= 1.0_rp/6.0_rp
          weigp(  2)= 1.0_rp/6.0_rp
          weigp(  3)= 1.0_rp/6.0_rp
       else if(ngaus==4) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 1.0_rp/3.0_rp
          posgp(2,4)= 1.0_rp/3.0_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 3.0_rp/8.0_rp
       else if(ngaus==6) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 0.5_rp
          posgp(2,4)= 0.0_rp
          posgp(1,5)= 0.5_rp
          posgp(2,5)= 0.5_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 0.5_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 1.0_rp/8.0_rp
          weigp(  5)= 1.0_rp/8.0_rp
          weigp(  6)= 1.0_rp/8.0_rp
       else if(ngaus==7) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(1,4)= 0.5_rp
          posgp(2,4)= 0.0_rp
          posgp(1,5)= 0.5_rp
          posgp(2,5)= 0.5_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 0.5_rp
          posgp(1,7)= 1.0_rp/3.0_rp
          posgp(2,7)= 1.0_rp/3.0_rp
          weigp(  1)= 1.0_rp/40.0_rp
          weigp(  2)= 1.0_rp/40.0_rp
          weigp(  3)= 1.0_rp/40.0_rp
          weigp(  4)= 1.0_rp/15.0_rp
          weigp(  5)= 1.0_rp/15.0_rp
          weigp(  6)= 1.0_rp/15.0_rp
          weigp(  7)= 9.0_rp/40.0_rp
       else if(ngaus==10) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(1, 4)= 1.0_rp/3.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(1, 5)= 2.0_rp/3.0_rp
          posgp(2, 5)= 0.0_rp
          posgp(1, 6)= 2.0_rp/3.0_rp 
          posgp(2, 6)= 1.0_rp/3.0_rp
          posgp(1, 7)= 1.0_rp/3.0_rp
          posgp(2, 7)= 2.0_rp/3.0_rp
          posgp(1, 8)= 0.0_rp
          posgp(2, 8)= 2.0_rp/3.0_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 1.0_rp/3.0_rp
          posgp(1,10)= 1.0_rp/3.0_rp
          posgp(2,10)= 1.0_rp/3.0_rp
          weigp(   1)= 1.0_rp/54.0_rp
          weigp(   2)= 1.0_rp/54.0_rp
          weigp(   3)= 1.0_rp/54.0_rp
          weigp(   4)= 1.0_rp/18.0_rp
          weigp(   5)= 1.0_rp/18.0_rp
          weigp(   6)= 1.0_rp/18.0_rp
          weigp(   7)= 1.0_rp/18.0_rp
          weigp(   8)= 1.0_rp/18.0_rp 
          weigp(   9)= 1.0_rp/18.0_rp 
          weigp(  10)= 1.0_rp/9.0_rp
       else
          istop=1
       end if

       ! Volume integral ( tetrahedra )
    else if(ndime==3) then
       if(ngaus==4) then
          posgp(1,1)= 0.0_rp 
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          weigp(  1)= 1.0_rp/24.0_rp
          weigp(  2)= 1.0_rp/24.0_rp
          weigp(  3)= 1.0_rp/24.0_rp
          weigp(  4)= 1.0_rp/24.0_rp
       else if(ngaus==5) then
          posgp(1,1)= 0.0_rp 
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          posgp(1,5)= 1.0_rp/4.0_rp
          posgp(2,5)= 1.0_rp/4.0_rp
          posgp(3,5)= 1.0_rp/4.0_rp
          weigp(  1)= 1.0_rp/120.0_rp
          weigp(  2)= 1.0_rp/120.0_rp
          weigp(  3)= 1.0_rp/120.0_rp
          weigp(  4)= 1.0_rp/120.0_rp
          weigp(  5)= 2.0_rp/15.0_rp
       else if(ngaus==10) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.5_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 0.5_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.0_rp
          posgp(3,10)= 0.5_rp
          weigp(   1)=-1.0_rp/120.0_rp
          weigp(   2)=-1.0_rp/120.0_rp 
          weigp(   3)=-1.0_rp/120.0_rp 
          weigp(   4)=-1.0_rp/120.0_rp 
          weigp(   5)= 1.0_rp/30.0_rp 
          weigp(   6)= 1.0_rp/30.0_rp 
          weigp(   7)= 1.0_rp/30.0_rp 
          weigp(   8)= 1.0_rp/30.0_rp 
          weigp(   9)= 1.0_rp/30.0_rp 
          weigp(  10)= 1.0_rp/30.0_rp 
       else if(ngaus==11) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.5_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 0.5_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.0_rp
          posgp(3,10)= 0.5_rp
          posgp(1,11)= 1.0_rp/4.0_rp
          posgp(2,11)= 1.0_rp/4.0_rp
          posgp(3,11)= 1.0_rp/4.0_rp
          weigp(   1)= 1.0_rp/360.0_rp
          weigp(   2)= 1.0_rp/360.0_rp
          weigp(   3)= 1.0_rp/360.0_rp
          weigp(   4)= 1.0_rp/360.0_rp
          weigp(   5)= 1.0_rp/90.0_rp
          weigp(   6)= 1.0_rp/90.0_rp
          weigp(   7)= 1.0_rp/90.0_rp
          weigp(   8)= 1.0_rp/90.0_rp
          weigp(   9)= 1.0_rp/90.0_rp
          weigp(  10)= 1.0_rp/90.0_rp
          weigp(  11)= 4.0_rp/45.0_rp
       else if(ngaus==15) then
          posgp(1, 1)= 0.0_rp 
          posgp(2, 1)= 0.0_rp
          posgp(3, 1)= 0.0_rp
          posgp(1, 2)= 1.0_rp
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 0.5_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 0.5_rp
          posgp(2, 6)= 0.5_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 0.0_rp
          posgp(2, 7)= 0.5_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 0.0_rp
          posgp(2, 8)= 0.0_rp
          posgp(3, 8)= 0.5_rp
          posgp(1, 9)= 0.5_rp
          posgp(2, 9)= 0.0_rp
          posgp(3, 9)= 0.5_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 0.5_rp
          posgp(3,10)= 0.5_rp
          posgp(1,11)= 1.0_rp/3.0_rp
          posgp(2,11)= 1.0_rp/3.0_rp
          posgp(3,11)= 0.0_rp
          posgp(1,12)= 1.0_rp/3.0_rp
          posgp(2,12)= 0.0_rp
          posgp(3,12)= 1.0_rp/3.0_rp
          posgp(1,13)= 1.0_rp/3.0_rp
          posgp(2,13)= 1.0_rp/3.0_rp
          posgp(3,13)= 1.0_rp/3.0_rp
          posgp(1,14)= 0.0_rp
          posgp(2,14)= 1.0_rp/3.0_rp
          posgp(3,14)= 1.0_rp/3.0_rp
          posgp(1,15)= 1.0_rp/4.0_rp
          posgp(2,15)= 1.0_rp/4.0_rp
          posgp(3,15)= 1.0_rp/4.0_rp
          weigp(   1)= 17.0_rp/5040.0_rp
          weigp(   2)= 17.0_rp/5040.0_rp
          weigp(   3)= 17.0_rp/5040.0_rp
          weigp(   4)= 17.0_rp/5040.0_rp
          weigp(   5)=  2.0_rp/315.0_rp
          weigp(   6)=  2.0_rp/315.0_rp
          weigp(   7)=  2.0_rp/315.0_rp
          weigp(   8)=  2.0_rp/315.0_rp
          weigp(   9)=  2.0_rp/315.0_rp
          weigp(  10)=  2.0_rp/315.0_rp
          weigp(  11)=  9.0_rp/840.0_rp
          weigp(  12)=  9.0_rp/840.0_rp
          weigp(  13)=  9.0_rp/840.0_rp
          weigp(  14)=  9.0_rp/840.0_rp
          weigp(  15)= 16.0_rp/315.0_rp
       else if(ngaus==20) then                           ! Integrates P2
          posgp(1, 1)= 0.0_rp                                ! and quartic mono-
          posgp(2, 1)= 0.0_rp                                ! mials to avoid zero 
          posgp(3, 1)= 0.0_rp                                ! weights at the 
          posgp(1, 2)= 1.0_rp                                ! edges
          posgp(2, 2)= 0.0_rp
          posgp(3, 2)= 0.0_rp
          posgp(1, 3)= 0.0_rp
          posgp(2, 3)= 1.0_rp
          posgp(3, 3)= 0.0_rp
          posgp(1, 4)= 0.0_rp
          posgp(2, 4)= 0.0_rp
          posgp(3, 4)= 1.0_rp
          posgp(1, 5)= 1.0_rp/3.0_rp
          posgp(2, 5)= 0.0_rp
          posgp(3, 5)= 0.0_rp
          posgp(1, 6)= 2.0_rp/3.0_rp
          posgp(2, 6)= 0.0_rp
          posgp(3, 6)= 0.0_rp
          posgp(1, 7)= 2.0_rp/3.0_rp
          posgp(2, 7)= 1.0_rp/3.0_rp
          posgp(3, 7)= 0.0_rp
          posgp(1, 8)= 1.0_rp/3.0_rp
          posgp(2, 8)= 2.0_rp/3.0_rp
          posgp(3, 8)= 0.0_rp
          posgp(1, 9)= 0.0_rp
          posgp(2, 9)= 2.0_rp/3.0_rp
          posgp(3, 9)= 0.0_rp
          posgp(1,10)= 0.0_rp
          posgp(2,10)= 1.0_rp/3.0_rp
          posgp(3,10)= 0.0_rp
          posgp(1,11)= 0.0_rp
          posgp(2,11)= 0.0_rp
          posgp(3,11)= 1.0_rp/3.0_rp
          posgp(1,12)= 2.0_rp/3.0_rp
          posgp(2,12)= 0.0_rp
          posgp(3,12)= 1.0_rp/3.0_rp
          posgp(1,13)= 0.0_rp
          posgp(2,13)= 2.0_rp/3.0_rp
          posgp(3,13)= 1.0_rp/3.0_rp
          posgp(1,14)= 0.0_rp
          posgp(2,14)= 0.0_rp
          posgp(3,14)= 2.0_rp/3.0_rp
          posgp(1,15)= 1.0_rp/3.0_rp
          posgp(2,15)= 0.0_rp
          posgp(3,15)= 2.0_rp/3.0_rp
          posgp(1,16)= 0.0_rp
          posgp(2,16)= 1.0_rp/3.0_rp
          posgp(3,16)= 2.0_rp/3.0_rp
          posgp(1,17)= 1.0_rp/3.0_rp
          posgp(2,17)= 1.0_rp/3.0_rp
          posgp(3,17)= 0.0_rp
          posgp(1,18)= 1.0_rp/3.0_rp
          posgp(2,18)= 0.0_rp
          posgp(3,18)= 1.0_rp/3.0_rp
          posgp(1,19)= 1.0_rp/3.0_rp
          posgp(2,19)= 1.0_rp/3.0_rp
          posgp(3,19)= 1.0_rp/3.0_rp
          posgp(1,20)= 0.0_rp
          posgp(2,20)= 1.0_rp/3.0_rp
          posgp(3,20)= 1.0_rp/3.0_rp
          w1 = 383.0_rp/2400.0_rp
          w2 = 1.0_rp/240.0_rp - w1
          w3 = 3.0_rp/80.0_rp  - w2
          weigp(   1)= w1
          weigp(   2)= w1
          weigp(   3)= w1
          weigp(   4)= w1
          weigp(   5)= w2
          weigp(   6)= w2
          weigp(   7)= w2
          weigp(   8)= w2
          weigp(   9)= w2
          weigp(  10)= w2
          weigp(  11)= w2
          weigp(  12)= w2
          weigp(  13)= w2
          weigp(  14)= w2
          weigp(  15)= w2
          weigp(  16)= w2
          weigp(  17)= w3
          weigp(  18)= w3
          weigp(  19)= w3
          weigp(  20)= w3
       else
          istop=1
       end if
    end if

  end subroutine rutclo

  !-----------------------------------------------------------------------
  subroutine rupope(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    ! 
    !     This routine sets up the integration constants of open rules for
    !     PRYSMS
    ! 
    !             NDIME = 3    
    ! 
    !          NGAUS  EXACT POL.
    !          -----  ----------  
    !            1       p1       
    !            6       p2       
    !       
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    integer(ip)              :: istop

    istop=0

    ! Volume integral
    if(ndime==3) then
       if(ngaus==1) then
          posgp(1,1)= 1.0_rp/3.0_rp
          posgp(2,1)= 1.0_rp/3.0_rp
          posgp(3,1)= 1.0_rp/2.0_rp
          weigp(  1)= 1.0_rp/2.0_rp
       else if(ngaus==6) then
          posgp(1,1)= 2.0_rp/3.0_rp
          posgp(2,1)= 1.0_rp/6.0_rp
          posgp(3,1)= 0.21132486540518711774542560974902_rp
          posgp(1,2)= 1.0_rp/6.0_rp
          posgp(2,2)= 2.0_rp/3.0_rp
          posgp(3,2)= 0.21132486540518711774542560974902_rp
          posgp(1,3)= 1.0_rp/6.0_rp
          posgp(2,3)= 1.0_rp/6.0_rp
          posgp(3,3)= 0.21132486540518711774542560974902_rp
          posgp(1,4)= 2.0_rp/3.0_rp
          posgp(2,4)= 1.0_rp/6.0_rp
          posgp(3,4)= 0.78867513459481288225457439025098_rp
          posgp(1,5)= 1.0_rp/6.0_rp
          posgp(2,5)= 2.0_rp/3.0_rp
          posgp(3,5)= 0.78867513459481288225457439025098_rp
          posgp(1,6)= 1.0_rp/6.0_rp
          posgp(2,6)= 1.0_rp/6.0_rp
          posgp(3,6)= 0.78867513459481288225457439025098_rp
          weigp(  1)= 1.0_rp/12.0_rp
          weigp(  2)= 1.0_rp/12.0_rp
          weigp(  3)= 1.0_rp/12.0_rp
          weigp(  4)= 1.0_rp/12.0_rp
          weigp(  5)= 1.0_rp/12.0_rp
          weigp(  6)= 1.0_rp/12.0_rp
       else
          istop=1
       end if
    end if

  end subroutine rupope

  !-----------------------------------------------------------------------
  subroutine rupclo(ndime,ngaus,posgp,weigp)
    !-----------------------------------------------------------------------
    ! 
    !     This routine sets up the integration constants of closed rules for
    !     PRISMS
    ! 
    !             NDIME = 3    
    ! 
    !          NGAUS  EXACT POL.
    !          -----  ----------  
    !            3       p1       
    ! 
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ngaus
    real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
    integer(ip)              :: istop

    ! Volume integral
    if(ndime==3) then
       if(ngaus==6) then
          posgp(1,1)= 0.0_rp
          posgp(2,1)= 0.0_rp
          posgp(3,1)= 0.0_rp
          posgp(1,2)= 1.0_rp
          posgp(2,2)= 0.0_rp
          posgp(3,2)= 0.0_rp
          posgp(1,3)= 0.0_rp
          posgp(2,3)= 1.0_rp
          posgp(3,3)= 0.0_rp
          posgp(1,4)= 0.0_rp
          posgp(2,4)= 0.0_rp
          posgp(3,4)= 1.0_rp
          posgp(1,5)= 1.0_rp
          posgp(2,5)= 0.0_rp
          posgp(3,5)= 1.0_rp
          posgp(1,6)= 0.0_rp
          posgp(2,6)= 1.0_rp
          posgp(3,6)= 1.0_rp
          weigp(  1)= 1.0_rp/12.0_rp
          weigp(  2)= 1.0_rp/12.0_rp
          weigp(  3)= 1.0_rp/12.0_rp
          weigp(  4)= 1.0_rp/12.0_rp
          weigp(  5)= 1.0_rp/12.0_rp
          weigp(  6)= 1.0_rp/12.0_rp
       end if
    end if

  end subroutine rupclo

  !-----------------------------------------------------------------------
  subroutine chaord(inoga,ninte)
    !-----------------------------------------------------------------------
    !
    !     Change the ordering of integration points for closed integration
    !     rules, according to : INOGA(IGAUS)=INODE
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ninte
    integer(ip), intent(out) :: inoga(ninte)

    if(ninte==4) then                                   ! 2 x 2 
       inoga(1)= 1
       inoga(2)= 4
       inoga(3)= 2
       inoga(4)= 3
    else if(ninte==9) then                              ! 3 x 3
       inoga(1)= 1
       inoga(2)= 8
       inoga(3)= 4
       inoga(4)= 5
       inoga(5)= 9
       inoga(6)= 7
       inoga(7)= 2
       inoga(8)= 6
       inoga(9)= 3
    else if(ninte==16) then                             ! 4 x 4
       inoga( 1)= 1
       inoga( 2)=12
       inoga( 3)=11
       inoga( 4)= 4
       inoga( 5)= 5                                        ! 4  10   9   3
       inoga( 6)=13                                        !
       inoga( 7)=16                                        ! 11 16   15  8
       inoga( 8)=10                                        ! 
       inoga( 9)= 6                                        ! 12 13   14  7
       inoga(10)=14                                        ! 
       inoga(11)=15                                        ! 1   5   6   2
       inoga(12)= 9                                        
       inoga(13)= 2
       inoga(14)= 7
       inoga(15)= 8
       inoga(16)= 3
    else if(ninte==8) then                              ! 2 x 2 x 2 
       inoga(1)= 1
       inoga(2)= 5
       inoga(3)= 4
       inoga(4)= 8
       inoga(5)= 2
       inoga(6)= 6
       inoga(7)= 3
       inoga(8)= 7
    else if(ninte==27) then                             ! 3 x 3 x 3 
       inoga( 1)= 1
       inoga( 2)=13
       inoga( 3)= 5
       inoga( 4)=12
       inoga( 5)=25
       inoga( 6)=20
       inoga( 7)= 4
       inoga( 8)=16
       inoga( 9)= 8
       inoga(10)= 9
       inoga(11)=22
       inoga(12)=17
       inoga(13)=21
       inoga(14)=27
       inoga(15)=26
       inoga(16)=11
       inoga(17)=24
       inoga(18)=19
       inoga(19)= 2
       inoga(20)=14
       inoga(21)= 6
       inoga(22)=10
       inoga(23)=23
       inoga(24)=18
       inoga(25)= 3
       inoga(26)=15
       inoga(27)= 7
    else if(ninte==64) then
       inoga( 1)= 1
       inoga( 2)=12
       inoga( 3)=21
       inoga( 4)= 5
       inoga( 5)=16
       inoga( 6)=44
       inoga( 7)=52
       inoga( 8)=32
       inoga( 9)=15
       inoga(10)=43
       inoga(11)=51
       inoga(12)=31
       inoga(13)= 4
       inoga(14)=20
       inoga(15)=24
       inoga(16)= 8
       inoga(17)= 9
       inoga(18)=37
       inoga(19)=45
       inoga(20)=25
       inoga(21)=33
       inoga(22)=57
       inoga(23)=61
       inoga(24)=53
       inoga(25)=36
       inoga(26)=60
       inoga(27)=64
       inoga(28)=56
       inoga(29)=14
       inoga(30)=42
       inoga(31)=50
       inoga(32)=30
       inoga(33)=10
       inoga(34)=38
       inoga(35)=46
       inoga(36)=26
       inoga(37)=34
       inoga(38)=58
       inoga(39)=62
       inoga(40)=54
       inoga(41)=35
       inoga(42)=59
       inoga(43)=63
       inoga(44)=55
       inoga(45)=13
       inoga(46)=41
       inoga(47)=49
       inoga(48)=29
       inoga(49)= 2
       inoga(50)=18
       inoga(51)=22
       inoga(52)= 6
       inoga(53)=11
       inoga(54)=39
       inoga(55)=47
       inoga(56)=27
       inoga(57)=12
       inoga(58)=40
       inoga(59)=48
       inoga(60)=28
       inoga(61)= 3
       inoga(62)=19
       inoga(63)=23
       inoga(64)= 7
    end if

  end subroutine chaord

end module quadrature_names
