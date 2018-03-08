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
# include "debug.i90"
program test_quadratures
  use fempar_names
  implicit none
  type(quadrature_t) :: quad
  real(rp)           :: posgl(20),weigl(20), errx, errw
  integer(ip) :: i,nlocs

  call fempar_init()

  write(*,*) ' Num points   Error in x                Error in w'
  
  do nlocs=1,20
     call quad%create(1,nlocs)
     call quad%fill_hex_gauss_legendre()

     posgl = 0.0_rp
     weigl = 0.0_rp

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
     else if(nlocs== 9 )  then 
        posgl( 1 ) = 0.968160239507626_rp 
        posgl( 2 ) = 0.836031107326636_rp 
        posgl( 3 ) = 0.613371432700590_rp 
        posgl( 4 ) = 0.324253423403809_rp 
        posgl( 5 ) = 0.000000000000000_rp 
        posgl( 6 ) = -0.324253423403809_rp 
        posgl( 7 ) = -0.613371432700590_rp 
        posgl( 8 ) = -0.836031107326636_rp 
        posgl( 9 ) = -0.968160239507626_rp 

        weigl( 1 ) = 0.081274388361575_rp 
        weigl( 2 ) = 0.180648160694857_rp 
        weigl( 3 ) = 0.260610696402936_rp 
        weigl( 4 ) = 0.312347077040003_rp 
        weigl( 5 ) = 0.330239355001260_rp 
        weigl( 6 ) = 0.312347077040003_rp 
        weigl( 7 ) = 0.260610696402936_rp 
        weigl( 8 ) = 0.180648160694857_rp 
        weigl( 9 ) = 0.081274388361575_rp 
     else if(nlocs== 10 )  then 
        posgl( 1 ) = 0.973906528517172_rp 
        posgl( 2 ) = 0.865063366688985_rp 
        posgl( 3 ) = 0.679409568299024_rp 
        posgl( 4 ) = 0.433395394129247_rp 
        posgl( 5 ) = 0.148874338981631_rp 
        posgl( 6 ) = -0.148874338981631_rp 
        posgl( 7 ) = -0.433395394129247_rp 
        posgl( 8 ) = -0.679409568299024_rp 
        posgl( 9 ) = -0.865063366688985_rp 
        posgl( 10 ) = -0.973906528517172_rp 

        weigl( 1 ) = 0.066671344308688_rp 
        weigl( 2 ) = 0.149451349150581_rp 
        weigl( 3 ) = 0.219086362515982_rp 
        weigl( 4 ) = 0.269266719309996_rp 
        weigl( 5 ) = 0.295524224714753_rp 
        weigl( 6 ) = 0.295524224714753_rp 
        weigl( 7 ) = 0.269266719309996_rp 
        weigl( 8 ) = 0.219086362515982_rp 
        weigl( 9 ) = 0.149451349150581_rp 
        weigl( 10 ) = 0.066671344308688_rp 
     else if(nlocs== 11 )  then 
        posgl( 1 ) = 0.978228658146057_rp 
        posgl( 2 ) = 0.887062599768095_rp 
        posgl( 3 ) = 0.730152005574049_rp 
        posgl( 4 ) = 0.519096129206812_rp 
        posgl( 5 ) = 0.269543155952345_rp 
        posgl( 6 ) = 0.000000000000000_rp 
        posgl( 7 ) = -0.269543155952345_rp 
        posgl( 8 ) = -0.519096129206812_rp 
        posgl( 9 ) = -0.730152005574049_rp 
        posgl( 10 ) = -0.887062599768095_rp 
        posgl( 11 ) = -0.978228658146057_rp 

        weigl( 1 ) = 0.055668567116174_rp 
        weigl( 2 ) = 0.125580369464904_rp 
        weigl( 3 ) = 0.186290210927734_rp 
        weigl( 4 ) = 0.233193764591990_rp 
        weigl( 5 ) = 0.262804544510247_rp 
        weigl( 6 ) = 0.272925086777901_rp 
        weigl( 7 ) = 0.262804544510247_rp 
        weigl( 8 ) = 0.233193764591990_rp 
        weigl( 9 ) = 0.186290210927734_rp 
        weigl( 10 ) = 0.125580369464904_rp 
        weigl( 11 ) = 0.055668567116174_rp 
     else if(nlocs== 12 )  then 
        posgl( 1 ) = 0.981560634246719_rp 
        posgl( 2 ) = 0.904117256370475_rp 
        posgl( 3 ) = 0.769902674194305_rp 
        posgl( 4 ) = 0.587317954286617_rp 
        posgl( 5 ) = 0.367831498998180_rp 
        posgl( 6 ) = 0.125233408511469_rp 
        posgl( 7 ) = -0.125233408511469_rp 
        posgl( 8 ) = -0.367831498998180_rp 
        posgl( 9 ) = -0.587317954286617_rp 
        posgl( 10 ) = -0.769902674194305_rp 
        posgl( 11 ) = -0.904117256370475_rp 
        posgl( 12 ) = -0.981560634246719_rp 

        weigl( 1 ) = 0.047175336386512_rp 
        weigl( 2 ) = 0.106939325995318_rp 
        weigl( 3 ) = 0.160078328543346_rp 
        weigl( 4 ) = 0.203167426723066_rp 
        weigl( 5 ) = 0.233492536538355_rp 
        weigl( 6 ) = 0.249147045813403_rp 
        weigl( 7 ) = 0.249147045813403_rp 
        weigl( 8 ) = 0.233492536538355_rp 
        weigl( 9 ) = 0.203167426723066_rp 
        weigl( 10 ) = 0.160078328543346_rp 
        weigl( 11 ) = 0.106939325995318_rp 
        weigl( 12 ) = 0.047175336386512_rp 
     else if(nlocs== 13 )  then 
        posgl( 1 ) = 0.984183054718588_rp 
        posgl( 2 ) = 0.917598399222978_rp 
        posgl( 3 ) = 0.801578090733310_rp 
        posgl( 4 ) = 0.642349339440340_rp 
        posgl( 5 ) = 0.448492751036447_rp 
        posgl( 6 ) = 0.230458315955135_rp 
        posgl( 7 ) = 0.000000000000000_rp 
        posgl( 8 ) = -0.230458315955135_rp 
        posgl( 9 ) = -0.448492751036447_rp 
        posgl( 10 ) = -0.642349339440340_rp 
        posgl( 11 ) = -0.801578090733310_rp 
        posgl( 12 ) = -0.917598399222978_rp 
        posgl( 13 ) = -0.984183054718588_rp 

        weigl( 1 ) = 0.040484004765316_rp 
        weigl( 2 ) = 0.092121499837728_rp 
        weigl( 3 ) = 0.138873510219787_rp 
        weigl( 4 ) = 0.178145980761946_rp 
        weigl( 5 ) = 0.207816047536888_rp 
        weigl( 6 ) = 0.226283180262897_rp 
        weigl( 7 ) = 0.232551553230874_rp 
        weigl( 8 ) = 0.226283180262897_rp 
        weigl( 9 ) = 0.207816047536888_rp 
        weigl( 10 ) = 0.178145980761946_rp 
        weigl( 11 ) = 0.138873510219787_rp 
        weigl( 12 ) = 0.092121499837728_rp 
        weigl( 13 ) = 0.040484004765316_rp 
     else if(nlocs== 14 )  then 
        posgl( 1 ) = 0.986283808696812_rp 
        posgl( 2 ) = 0.928434883663574_rp 
        posgl( 3 ) = 0.827201315069765_rp 
        posgl( 4 ) = 0.687292904811685_rp 
        posgl( 5 ) = 0.515248636358154_rp 
        posgl( 6 ) = 0.319112368927890_rp 
        posgl( 7 ) = 0.108054948707344_rp 
        posgl( 8 ) = -0.108054948707344_rp 
        posgl( 9 ) = -0.319112368927890_rp 
        posgl( 10 ) = -0.515248636358154_rp 
        posgl( 11 ) = -0.687292904811685_rp 
        posgl( 12 ) = -0.827201315069765_rp 
        posgl( 13 ) = -0.928434883663574_rp 
        posgl( 14 ) = -0.986283808696812_rp 

        weigl( 1 ) = 0.035119460331752_rp 
        weigl( 2 ) = 0.080158087159760_rp 
        weigl( 3 ) = 0.121518570687903_rp 
        weigl( 4 ) = 0.157203167158194_rp 
        weigl( 5 ) = 0.185538397477938_rp 
        weigl( 6 ) = 0.205198463721296_rp 
        weigl( 7 ) = 0.215263853463158_rp 
        weigl( 8 ) = 0.215263853463158_rp 
        weigl( 9 ) = 0.205198463721296_rp 
        weigl( 10 ) = 0.185538397477938_rp 
        weigl( 11 ) = 0.157203167158194_rp 
        weigl( 12 ) = 0.121518570687903_rp 
        weigl( 13 ) = 0.080158087159760_rp 
        weigl( 14 ) = 0.035119460331752_rp 
     else if(nlocs== 15 )  then 
        posgl( 1 ) = 0.987992518020485_rp 
        posgl( 2 ) = 0.937273392400706_rp 
        posgl( 3 ) = 0.848206583410427_rp 
        posgl( 4 ) = 0.724417731360170_rp 
        posgl( 5 ) = 0.570972172608539_rp 
        posgl( 6 ) = 0.394151347077563_rp 
        posgl( 7 ) = 0.201194093997435_rp 
        posgl( 8 ) = 0.000000000000000_rp 
        posgl( 9 ) = -0.201194093997435_rp 
        posgl( 10 ) = -0.394151347077563_rp 
        posgl( 11 ) = -0.570972172608539_rp 
        posgl( 12 ) = -0.724417731360170_rp 
        posgl( 13 ) = -0.848206583410427_rp 
        posgl( 14 ) = -0.937273392400706_rp 
        posgl( 15 ) = -0.987992518020485_rp 

        weigl( 1 ) = 0.030753241996117_rp 
        weigl( 2 ) = 0.070366047488108_rp 
        weigl( 3 ) = 0.107159220467172_rp 
        weigl( 4 ) = 0.139570677926154_rp 
        weigl( 5 ) = 0.166269205816994_rp 
        weigl( 6 ) = 0.186161000015562_rp 
        weigl( 7 ) = 0.198431485327112_rp 
        weigl( 8 ) = 0.202578241925561_rp 
        weigl( 9 ) = 0.198431485327112_rp 
        weigl( 10 ) = 0.186161000015562_rp 
        weigl( 11 ) = 0.166269205816994_rp 
        weigl( 12 ) = 0.139570677926154_rp 
        weigl( 13 ) = 0.107159220467172_rp 
        weigl( 14 ) = 0.070366047488108_rp 
        weigl( 15 ) = 0.030753241996117_rp 
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
     else if(nlocs== 17 )  then 
        posgl( 1 ) = 0.990575475314417_rp 
        posgl( 2 ) = 0.950675521768768_rp 
        posgl( 3 ) = 0.880239153726986_rp 
        posgl( 4 ) = 0.781514003896801_rp 
        posgl( 5 ) = 0.657671159216691_rp 
        posgl( 6 ) = 0.512690537086477_rp 
        posgl( 7 ) = 0.351231763453876_rp 
        posgl( 8 ) = 0.178484181495848_rp 
        posgl( 9 ) = 0.000000000000000_rp 
        posgl( 10 ) = -0.178484181495848_rp 
        posgl( 11 ) = -0.351231763453876_rp 
        posgl( 12 ) = -0.512690537086477_rp 
        posgl( 13 ) = -0.657671159216691_rp 
        posgl( 14 ) = -0.781514003896801_rp 
        posgl( 15 ) = -0.880239153726986_rp 
        posgl( 16 ) = -0.950675521768768_rp 
        posgl( 17 ) = -0.990575475314417_rp 

        weigl( 1 ) = 0.024148302868548_rp 
        weigl( 2 ) = 0.055459529373987_rp 
        weigl( 3 ) = 0.085036148317179_rp 
        weigl( 4 ) = 0.111883847193404_rp 
        weigl( 5 ) = 0.135136368468525_rp 
        weigl( 6 ) = 0.154045761076810_rp 
        weigl( 7 ) = 0.168004102156450_rp 
        weigl( 8 ) = 0.176562705366993_rp 
        weigl( 9 ) = 0.179446470356207_rp 
        weigl( 10 ) = 0.176562705366993_rp 
        weigl( 11 ) = 0.168004102156450_rp 
        weigl( 12 ) = 0.154045761076810_rp 
        weigl( 13 ) = 0.135136368468525_rp 
        weigl( 14 ) = 0.111883847193404_rp 
        weigl( 15 ) = 0.085036148317179_rp 
        weigl( 16 ) = 0.055459529373987_rp 
        weigl( 17 ) = 0.024148302868548_rp 
     else if(nlocs== 18 )  then 
        posgl( 1 ) = 0.991565168420931_rp 
        posgl( 2 ) = 0.955823949571398_rp 
        posgl( 3 ) = 0.892602466497556_rp 
        posgl( 4 ) = 0.803704958972523_rp 
        posgl( 5 ) = 0.691687043060353_rp 
        posgl( 6 ) = 0.559770831073948_rp 
        posgl( 7 ) = 0.411751161462843_rp 
        posgl( 8 ) = 0.251886225691506_rp 
        posgl( 9 ) = 0.084775013041735_rp 
        posgl( 10 ) = -0.084775013041735_rp 
        posgl( 11 ) = -0.251886225691506_rp 
        posgl( 12 ) = -0.411751161462843_rp 
        posgl( 13 ) = -0.559770831073948_rp 
        posgl( 14 ) = -0.691687043060353_rp 
        posgl( 15 ) = -0.803704958972523_rp 
        posgl( 16 ) = -0.892602466497556_rp 
        posgl( 17 ) = -0.955823949571398_rp 
        posgl( 18 ) = -0.991565168420931_rp 

        weigl( 1 ) = 0.021616013526483_rp 
        weigl( 2 ) = 0.049714548894969_rp 
        weigl( 3 ) = 0.076425730254889_rp 
        weigl( 4 ) = 0.100942044106287_rp 
        weigl( 5 ) = 0.122555206711478_rp 
        weigl( 6 ) = 0.140642914670651_rp 
        weigl( 7 ) = 0.154684675126265_rp 
        weigl( 8 ) = 0.164276483745833_rp 
        weigl( 9 ) = 0.169142382963144_rp 
        weigl( 10 ) = 0.169142382963144_rp 
        weigl( 11 ) = 0.164276483745833_rp 
        weigl( 12 ) = 0.154684675126265_rp 
        weigl( 13 ) = 0.140642914670651_rp 
        weigl( 14 ) = 0.122555206711478_rp 
        weigl( 15 ) = 0.100942044106287_rp 
        weigl( 16 ) = 0.076425730254889_rp 
        weigl( 17 ) = 0.049714548894969_rp 
        weigl( 18 ) = 0.021616013526483_rp 
     else if(nlocs== 19 )  then 
        posgl( 1 ) = 0.992406843843584_rp 
        posgl( 2 ) = 0.960208152134830_rp 
        posgl( 3 ) = 0.903155903614818_rp 
        posgl( 4 ) = 0.822714656537143_rp 
        posgl( 5 ) = 0.720966177335229_rp 
        posgl( 6 ) = 0.600545304661681_rp 
        posgl( 7 ) = 0.464570741375961_rp 
        posgl( 8 ) = 0.316564099963630_rp 
        posgl( 9 ) = 0.160358645640225_rp 
        posgl( 10 ) = 0.000000000000000_rp 
        posgl( 11 ) = -0.160358645640225_rp 
        posgl( 12 ) = -0.316564099963630_rp 
        posgl( 13 ) = -0.464570741375961_rp 
        posgl( 14 ) = -0.600545304661681_rp 
        posgl( 15 ) = -0.720966177335229_rp 
        posgl( 16 ) = -0.822714656537143_rp 
        posgl( 17 ) = -0.903155903614818_rp 
        posgl( 18 ) = -0.960208152134830_rp 
        posgl( 19 ) = -0.992406843843584_rp 

        weigl( 1 ) = 0.019461788229726_rp 
        weigl( 2 ) = 0.044814226765699_rp 
        weigl( 3 ) = 0.069044542737641_rp 
        weigl( 4 ) = 0.091490021622450_rp 
        weigl( 5 ) = 0.111566645547334_rp 
        weigl( 6 ) = 0.128753962539336_rp 
        weigl( 7 ) = 0.142606702173607_rp 
        weigl( 8 ) = 0.152766042065860_rp 
        weigl( 9 ) = 0.158968843393954_rp 
        weigl( 10 ) = 0.161054449848784_rp 
        weigl( 11 ) = 0.158968843393954_rp 
        weigl( 12 ) = 0.152766042065860_rp 
        weigl( 13 ) = 0.142606702173607_rp 
        weigl( 14 ) = 0.128753962539336_rp 
        weigl( 15 ) = 0.111566645547334_rp 
        weigl( 16 ) = 0.091490021622450_rp 
        weigl( 17 ) = 0.069044542737641_rp 
        weigl( 18 ) = 0.044814226765699_rp 
        weigl( 19 ) = 0.019461788229726_rp 
     else if(nlocs== 20 )  then 
        posgl( 1 ) = 0.993128599185095_rp 
        posgl( 2 ) = 0.963971927277914_rp 
        posgl( 3 ) = 0.912234428251326_rp 
        posgl( 4 ) = 0.839116971822219_rp 
        posgl( 5 ) = 0.746331906460151_rp 
        posgl( 6 ) = 0.636053680726515_rp 
        posgl( 7 ) = 0.510867001950827_rp 
        posgl( 8 ) = 0.373706088715420_rp 
        posgl( 9 ) = 0.227785851141645_rp 
        posgl( 10 ) = 0.076526521133497_rp 
        posgl( 11 ) = -0.076526521133497_rp 
        posgl( 12 ) = -0.227785851141645_rp 
        posgl( 13 ) = -0.373706088715420_rp 
        posgl( 14 ) = -0.510867001950827_rp 
        posgl( 15 ) = -0.636053680726515_rp 
        posgl( 16 ) = -0.746331906460151_rp 
        posgl( 17 ) = -0.839116971822219_rp 
        posgl( 18 ) = -0.912234428251326_rp 
        posgl( 19 ) = -0.963971927277914_rp 
        posgl( 20 ) = -0.993128599185095_rp 

        weigl( 1 ) = 0.017614007139152_rp 
        weigl( 2 ) = 0.040601429800387_rp 
        weigl( 3 ) = 0.062672048334109_rp 
        weigl( 4 ) = 0.083276741576705_rp 
        weigl( 5 ) = 0.101930119817240_rp 
        weigl( 6 ) = 0.118194531961518_rp 
        weigl( 7 ) = 0.131688638449177_rp 
        weigl( 8 ) = 0.142096109318382_rp 
        weigl( 9 ) = 0.149172986472604_rp 
        weigl( 10 ) = 0.152753387130726_rp 
        weigl( 11 ) = 0.152753387130726_rp 
        weigl( 12 ) = 0.149172986472604_rp 
        weigl( 13 ) = 0.142096109318382_rp 
        weigl( 14 ) = 0.131688638449177_rp 
        weigl( 15 ) = 0.118194531961518_rp 
        weigl( 16 ) = 0.101930119817240_rp 
        weigl( 17 ) = 0.083276741576705_rp 
        weigl( 18 ) = 0.062672048334109_rp 
        weigl( 19 ) = 0.040601429800387_rp 
        weigl( 20 ) = 0.017614007139152_rp 
     else
        !mcheck(.false.,'ERROR:: Quadrature not defined')
     end if

     errx = 0.0_rp
     errw = 0.0_rp
     do i=1,nlocs
        errx=max(min(posgl(i)-quad%coordinates(1,i),posgl(nlocs+1-i)-quad%coordinates(1,i)),errx)
        errw=max(weigl(i)-quad%weight(i),errw)
     end do

     write(*,*) nlocs, errx, errw
     if(errx>1.0e-10_rp) then
        mcheck(.false.,'Error in quadrature coordintes')
     end if
     if(errw>1.0e-10_rp) then
        mcheck(.false.,'Error in quadrature weights')
     end if
     call quad%free()

  end do

  !nlocs = 60
  !call quad%create(1,nlocs)
  !call quad%fill_hex_gauss_legendre()
  !do i=1,nlocs
  !   write(*,*) quad%coordinates(1,i),quad%weight(i)
  !end do
  
end program test_quadratures
