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
module interpolation_names
  use types
  use memor
  use fem_space_types

  implicit none
# include "debug.i90"
  private

  type interpolation
     integer(ip)                :: &
          kfun = 1,                &    ! Shape functions calculation (1=yes 0=no)
          kder = 1,                &    ! Derivatives calculation (1=yes 0=no)
          khes = 0,                &    ! Hessian calculation (1=yes 0=no)
          nlocs = 1                     ! Number of interpolation locs 
                                        ! (usually integration points)

     integer(ip)                :: &
          nnode = 3,               &    ! Number of interpolation coefs (nodes)
          ndime = 2,               &    ! Number of space dimensions
          ntens = 3                     ! Number of tensor components 

     logical(lg)                :: &
          khie = .false.                ! Hierarchical shape functions (1=yes, 0=no)

     real(rp), allocatable      :: &
          shape(:,:),              &    ! Shape functions
          deriv(:,:,:),            &    ! Derivatives
          hessi(:,:,:)                  ! Hessian
  end type interpolation

  interface interpolation_local
     module procedure interpolation_local_tp, interpolation_local_old
  end interface interpolation_local

  interface shape1
     module procedure shape1_old, shape1_hessi
  end interface shape1

  ! Types
  public :: interpolation

  ! Functions
  public :: interpolation_create, interpolation_free, interpolation_local
  public :: shape1, shapen, shafun

contains
  !==============================================================================
  ! TODO: separate computations of shape, deriv and hessi allowing partial
  !       definition of the interpolation (shape only, for example or derivatives
  !       only, which is actually the case in the physical domain). See
  !       interpolation_local below.
  !==============================================================================  
  subroutine interpolation_create(kfun,kder,khes,ndime,nnode,nlocs,int,khie)
    implicit none
    integer(ip)          , intent(in)    :: kfun,kder,khes,nlocs
    integer(ip)          , intent(in)    :: nnode,ndime
    type(interpolation)  , intent(out)   :: int
    logical(lg), optional, intent(in)    :: khie
    integer(ip) :: iloc,ntens

    if(ndime==1) then
       ntens=1
    elseif(ndime==2) then
       ntens=3
    else if(ndime==3) then
       ntens=6
    end if
    int%kfun=kfun
    int%kder=kder
    int%khes=khes
    int%nlocs=nlocs
    int%nnode=nnode
    int%ndime=ndime
    int%ntens=ntens
    if(present(khie)) then
       int%khie=khie
    else
       int%khie=.false.
    end if

    call memalloc(nnode,nlocs,int%shape,__FILE__,__LINE__)
    if(kder==1) call memalloc(ndime,nnode,nlocs,int%deriv,   __FILE__,__LINE__)
    if(khes==1) call memalloc(ntens,nnode,nlocs,int%hessi,   __FILE__,__LINE__)

  end subroutine interpolation_create

  !==============================================================================
  subroutine interpolation_free(int)
    implicit none
    type(interpolation), intent(inout) :: int

    call memfree(int%shape,__FILE__,__LINE__)
    if(int%kder==1) call memfree(int%deriv,__FILE__,__LINE__)
    if(int%khes==1) call memfree(int%hessi,__FILE__,__LINE__)

  end subroutine interpolation_free

  !=============================================================================
  subroutine interpolation_local_old(clocs,int)
    implicit none
    type(interpolation), intent(inout) :: int
    real(rp)           , intent(in)    :: clocs(:,:)
    integer(ip) :: iloc,ndime,nnode,ntens,nlocs

    ! TODO:
    ! The next two constrains can be eliminated after
    ! appropriate modifications below:
    assert(int%kfun==1)
    assert(int%kder==1)
    assert(int%khes==1)

    ! TODO:
    ! Depending on kfun, kder and khes different arrays will be
    ! allocated and therefore some dimensions may not be needed.
    ! If subroutine shafun below is appropriately split, they need not to be stored.
    nnode=size(int%shape,dim=1)
    nlocs=size(int%shape,dim=2)
    ndime=size(int%deriv,dim=1)
    !nnode=size(int%deriv,dim=2)
    !nlocs=size(int%deriv,dim=3)
    ntens=size(int%hessi,dim=1)
    !nnode=size(int%hessi,dim=1)
    !nlocs=size(int%hessi,dim=1)

    ! Check input data (interpolation locs)
    assert(ndime==size(clocs,dim=1))
    assert(nlocs==size(clocs,dim=2))

    do iloc=1,int%nlocs
       call shafun(clocs(:,iloc),ndime,nnode,ntens,int%khie,int%kder,int%khes,int%shape(:,iloc),  &
            &        int%deriv(:,:,iloc),int%hessi(:,:,iloc))
    end do

  end subroutine interpolation_local_old

  ! ==================================================================================================
  subroutine interpolation_local_tp(clocs,int,nd,order,ngaus)
    implicit none
    ! Parameters
    type(interpolation), intent(inout) :: int
    integer(ip)        , intent(in)    :: nd, order, ngaus
    real(rp)           , intent(in)    :: clocs(ngaus)

    real(rp)    :: coord(order+1),shpe1(order+1,ngaus),shpd1(order+1,ngaus),shph1(order+1,ngaus)
    integer(ip) :: iloc,aux(ngaus)

    ! Set the coordenades of the nodal points
    call Q_coord_1d(coord,order+1)

    ! Compute the 1d shape function on the gauss points
    call shape1(clocs,order,ngaus,coord,shpe1,shpd1,shph1,int%khes)
    
    ! Compute the tensorial product
    if (int%khes == 1) then
       call shapen(int%shape,int%deriv,shpe1,shpd1,shph1,nd,order,ngaus,int%ntens,int%khes,int%hessi)
    else
       call shapen(int%shape,int%deriv,shpe1,shpd1,shph1,nd,order,ngaus,int%ntens,int%khes)
    end if
  end subroutine interpolation_local_tp

  !=========================================================================
  subroutine shafun(posgp,ndime,nnode,ntens,khie,kder,khes,shape,deriv,heslo)
    !-----------------------------------------------------------------------
    !
    ! This routine evaluates shape functions and their derivatives
    ! for linear and quadratic isoparametric elements
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)       , intent(in)  :: ndime,nnode,ntens,kder,khes
    logical(lg)       , intent(in)  :: khie
    real(rp)          , intent(in)  :: posgp(ndime)
    real(rp)          , intent(out) :: shape(nnode)
    real(rp), optional, intent(out) :: deriv(ndime,nnode)
    real(rp), optional, intent(out) :: heslo(ntens,nnode)

    ! Initializations
    shape=0.0_rp
    deriv=0.0_rp
    heslo=0.0_rp
    ! Evaluation of the shape functions
    if(ndime==1) then
       call shape1_old(posgp(1),nnode,kder,shape,deriv)
    else if(ndime==2) then
       call shape2(posgp(1),posgp(2),nnode,ntens,khie,kder,khes,shape,deriv,heslo)
    else if(ndime==3) then
       call shape3(posgp(1),posgp(2),posgp(3),nnode,ntens,khie,kder,khes,shape,deriv,heslo)
    end if

  end subroutine shafun

  ! ===================================================================================================
  ! Compute the shape function and its derivative
  subroutine shape1_hessi(xg,p,ng,xn,shpe1,shpd1,shph1,khes)
    implicit none
    integer(ip), intent(in)  :: p,ng,khes
    real(rp),    intent(in)  :: xn(p+1),xg(ng)
    real(rp),    intent(out) :: shpe1(p+1,ng),shpd1(p+1,ng),shph1(p+1,ng)
    integer(ip)              :: i,j,k,m,ig
    real(rp)                 :: aux, aux2, aux3, auxv(ng),auxv2(ng),auxv3(ng)

    shpe1 = 1.0_rp
    shpd1 = 0.0_rp
    shph1 = 0.0_rp
    do i = 1,p+1
       do j = 1,p+1
          if (j /= i) then
             aux = 1/(xn(i)-xn(j)) 
             auxv = 1/(xn(i)-xn(j))
             if (khes == 1) auxv3 = 0
             do k = 1,p+1
                if (k /= j .and. k /= i) then
                   aux2 = 1/(xn(i)-xn(k))
                   if (khes == 1) auxv2 = 1/(xn(i)-xn(k))
                   do ig = 1,ng
                      auxv(ig) = auxv(ig)*(xg(ig)-xn(k))*aux2
                   end do
                   if (khes == 1) then
                      do m = 1, p+1
                         if (m/=k .and. m/= j .and. m /= i) then
                            aux3 = 1/(xn(i)-xn(m))
                            do ig = 1,ng
                               auxv2(ig) = auxv2(ig)*(xg(ig)-xn(m))*aux3
                            end do
                         end if
                      end do
                   end if
                   auxv3 = auxv3+auxv2
                end if
             end do
             do ig = 1,ng
                shpe1(i,ig) = shpe1(i,ig)*(xg(ig)-xn(j))*aux
                shpd1(i,ig) = shpd1(i,ig) + auxv(ig)
                shph1(i,ig) = shph1(i,ig) + auxv3(ig)*aux
             end do
          end if
       end do
    end do
  end subroutine shape1_hessi

  !=========================================================================
  subroutine shape1_old(s,nnode,kder,shape,deriv)
    !-----------------------------------------------------------------------
    !
    ! This routine evaluates shape functions and their first derivates 
    ! for 1-d continuos with 2 & 3 nodes
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nnode,kder
    real(rp),    intent(in)  :: s
    real(rp),    intent(out) :: deriv(nnode), shape(nnode)

    real(rp)                 :: s1,s2,s3,s4,q,r


    if(nnode==2) then     
       shape(1)= 0.5_rp*(1.0_rp-s)
       shape(2)= 0.5_rp*(1.0_rp+s)
       if(kder == 1) then
          deriv(1)=-0.5_rp
          deriv(2)= 0.5_rp
       end if
    else if(nnode==3) then
       shape(1)= 0.5_rp*s*(s-1.0_rp)
       shape(3)= 0.5_rp*s*(s+1.0_rp)
       shape(2)=-(s+1.0_rp)*(s-1.0_rp)
       if(kder == 1) then
          deriv(1)= s-0.5_rp
          deriv(3)= s+0.5_rp
          deriv(2)=-2.0_rp*s
       end if
    else if(nnode==4) then
       s1 = (s+1.0_rp)
       s2 = (s-1.0_rp)
       s3 = (s + 1.0_rp/3.0_rp)
       s4 = (s - 1.0_rp/3.0_rp)
       q =  9.0_rp/16.0_rp
       r =  27.0_rp/16.0_rp
       shape(1)= -q*s2*s3*s4
       shape(2)=  q*s1*s3*s4
       shape(3)=  r*s1*s2*s4
       shape(4)=  -r*s1*s2*s3

       if(kder == 1) then
          deriv(1)= -q*(s2*s3 + s2*s4 + s3*s4)
          deriv(2)=  q*(s1*s3 + s1*s4 + s3*s4)
          deriv(3)=  r*(s2*s1 + s2*s4 + s1*s4)
          deriv(4)= -r*(s2*s3 + s2*s1 + s3*s1)
       end if
    end if

  end subroutine shape1_old
  
  !==================================================================================================
  subroutine shapen (shape,deriv,s1,sd1,sdd1,nd,p,ng,nt,khes,hessi)
    implicit none
    ! Parameters
    integer(ip)       , intent(in)  :: nd,p,ng,nt,khes
    real(rp)          , intent(in)  :: s1(p+1,ng),sd1(p+1,ng),sdd1(p+1,ng)
    real(rp)          , intent(out) :: shape((p+1)**nd,ng**nd)
    real(rp)          , intent(out) :: deriv(nd,(p+1)**nd,ng**nd)
    real(rp), optional, intent(out) :: hessi(nt,(p+1)**nd,ng**nd)


    ! Local variables
    integer(ip)              :: i,ig,d,d2,d3,it
    integer(ip)              :: ijk(nd),ijkg(nd),permu(nt)

    ! The presumed order of the varibles in the hessian is not the one obtained by generation
    if (nd == 2) then
       permu = (/ 1, 3, 2 /)
    elseif (nd == 3) then
       permu = (/ 1, 4, 5, 2, 6, 3/)
    end if

    ! Initialize values
    shape = 0.0_rp
    deriv = 0.0_rp
    if (khes == 1) hessi = 0.0_rp

    ! Initialize nodal coordinates vector
    ijk = 0; ijk(1) = -1
    do i = 1,(p+1)**nd

       ! Set coordinates of node i 
       call Q_ijkg(ijk,i,nd,p)

       ! Initialize Gauss point coordinates vector
       ijkg = 0; ijkg(1) = -1
       do ig = 1,ng**nd

          ! Set coordinates of Gauss point ig 
          call Q_ijkg(ijkg,ig,nd,ng-1)

          ! Initialize shape
          shape(i,ig) = 1.0_rp
          it = 0
          do d = 1,nd
             ! Shape is the tensor product 1d shape: s_ijk(x,y,z) = s_i(x)*s_j(y)*s_k(z)
             shape(i,ig) = shape(i,ig)*s1(ijk(d)+1,ijkg(d)+1)

             ! Initialize deriv and hessi
             deriv(d,i,ig) = 1.0_rp
             it = it+1
             if (khes == 1) hessi(permu(it),i,ig)= 1.0_rp

             ! Deriv: d(s_ijk)/dx (x,y,z) = s'_i(x)*s_j(y)*s_k(z)
             ! Hessi: d2(s_ijk)/dx2 (x,y,z) = s''_i(x)*s_j(y)*s_k(z)
             do d2 = 1,nd
                if (d2 /= d) then
                   deriv( d,i,ig) = deriv( d,i,ig)*s1(ijk(d2)+1,ijkg(d2)+1)
                    if (khes==1) hessi(permu(it),i,ig)=hessi(permu(it),i,ig)*s1(ijk(d2)+1,ijkg(d2)+1)
                else
                   deriv( d,i,ig) = deriv( d,i,ig)* sd1(ijk(d)+1,ijkg(d)+1)  
                   if (khes==1) hessi(permu(it),i,ig)=hessi(permu(it),i,ig)*sdd1(ijk(d)+1,ijkg(d)+1)             
                end if
             end do

             if (khes == 1) then
                ! Hessi: d2(s_ijk)/dxdy (x,y,z) = s'_i(x)*s'_j(y)*s_k(z)
                do d2 = d+1,nd
                   it = it+1
                   hessi(permu(it),i,ig) = 1.0_rp
                   do d3 = 1,nd
                      if (d3 /= d .and. d3 /= d2) then
                         hessi(permu(it),i,ig) = hessi(permu(it),i,ig)*s1(ijk(d3)+1,ijkg(d3)+1)
                      else
                         hessi(permu(it),i,ig) = hessi(permu(it),i,ig)*sd1(ijk(d3)+1,ijkg(d3)+1)             
                      end if
                   end do
                end do
             end if

          end do
       end do
    end do
  end subroutine shapen

  !=========================================================================
  subroutine shape2(s,t,nnode,ntens,khie,kder,khes,shape,deriv,heslo)
    !-----------------------------------------------------------------------
    !
    ! This routine evaluates shape functions and their first and second
    ! derivatives for 2-d continuos standar interpolation elements.
    !
    !    TRIANGLES       3   6  &  10  nodes
    !    QUADRILATERALS  4   9  &  16  nodes
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)       , intent(in)  :: nnode,ntens,kder,khes
    logical(lg)       , intent(in)  :: khie
    real(rp)          , intent(in)  :: s,t
    real(rp)          , intent(out) :: shape(nnode)
    real(rp), optional, intent(out) :: deriv(2,nnode),heslo(ntens,nnode)
    real(rp)                 :: st,a1,a2,a3,ss,tt,s1,s2,s3,s4
    real(rp)                 :: t1,t2,t3,t4,s9,t9,c,a
    if(nnode==3) then     
       shape(1)=1.0_rp-s-t                                
       shape(2)=s                                  
       shape(3)=t                                           
       if (kder == 1) then
          deriv(1,1)=-1.0_rp                                   !  3 
          deriv(1,2)= 1.0_rp                                   !
          deriv(1,3)= 0.0_rp                                   !
          deriv(2,1)=-1.0_rp                                   !         
          deriv(2,2)= 0.0_rp                                   !  1       2
          deriv(2,3)= 1.0_rp
       end if
    else if(nnode==4) then
       st=s*t                                           
       shape(1)=(1.0_rp-t-s+st)*0.25_rp                     !  3         4
       shape(2)=(1.0_rp-t+s-st)*0.25_rp                     !
       shape(4)=(1.0_rp+t+s+st)*0.25_rp                     !      
       shape(3)=(1.0_rp+t-s-st)*0.25_rp                     !
       if (kder == 1) then                                  !  1         2
          deriv(1,1)=(-1.0_rp+t)*0.25_rp                              
          deriv(1,2)=(+1.0_rp-t)*0.25_rp
          deriv(1,4)=(+1.0_rp+t)*0.25_rp
          deriv(1,3)=(-1.0_rp-t)*0.25_rp
          deriv(2,1)=(-1.0_rp+s)*0.25_rp
          deriv(2,2)=(-1.0_rp-s)*0.25_rp
          deriv(2,4)=(+1.0_rp+s)*0.25_rp
          deriv(2,3)=(+1.0_rp-s)*0.25_rp
       end if
       if (khes == 1) then
          heslo(3,1)= 0.25_rp
          heslo(3,2)=-0.25_rp
          heslo(3,4)= 0.25_rp
          heslo(3,3)=-0.25_rp
       end if
    else if(nnode==6) then
       a1=1.0_rp-s-t
       a2=s                                             
       a3=t
       shape( 1)=(2.0_rp*a1-1.0_rp)*a1                      !  6
       shape( 3)=(2.0_rp*a2-1.0_rp)*a2                      !   
       shape( 6)=(2.0_rp*a3-1.0_rp)*a3                      !   
       shape( 2)=4.0_rp*a1*a2                               !  4      5
       shape( 5)=4.0_rp*a2*a3                               !     
       shape( 4)=4.0_rp*a1*a3                               !
       if (kder == 1) then                                  !  1     2     3
          deriv(1,1)= 1.0_rp-4.0_rp*a1                        
          deriv(1,3)= 4.0_rp*a2-1.0_rp    
          deriv(1,6)= 0.0_rp           
          deriv(1,2)= 4.0_rp*(a1-a2)   
          deriv(1,5)= 4.0_rp*a3        
          deriv(1,4)=-4.0_rp*a3       
          deriv(2,1)= 1.0_rp-4.0_rp*a1    
          deriv(2,3)= 0.0_rp           
          deriv(2,6)= 4.0_rp*a3-1.0_rp    
          deriv(2,2)=-4.0_rp*a2       
          deriv(2,5)= 4.0_rp*a2        
          deriv(2,4)= 4.0_rp*(a1-a3)
       end if
       if (khes == 1) then
          heslo(1,1)= 4.0_rp
          heslo(1,3)= 4.0_rp
          heslo(1,2)=-8.0_rp
          heslo(2,1)= 4.0_rp
          heslo(2,6)= 4.0_rp
          heslo(2,4)=-8.0_rp
          heslo(3,1)= 4.0_rp
          heslo(3,2)=-4.0_rp
          heslo(3,5)= 4.0_rp
          heslo(3,4)=-4.0_rp
       end if
    else if(nnode==9) then
       ss=s*s
       st=s*t
       tt=t*t
       s1=s+1.0_rp
       t1=t+1.0_rp
       s2=s*2.0_rp
       t2=t*2.0_rp
       s9=s-1.0_rp                               
       t9=t-1.0_rp
       if(khie) then
          ! Quadratic (1 to 4)
          shape( 1)=0.25_rp*s9*st*t9                           !  4      7      3
          shape( 3)=0.25_rp*s1*st*t9                           !        
          shape( 9)=0.25_rp*s1*st*t1                           !      
          shape( 7)=0.25_rp*s9*st*t1                           !                    
          if (kder == 1) then
             deriv(1,1)= 0.25_rp*t*t9*(-1.0_rp+s2)                !  8      9      6
             deriv(1,3)= 0.25_rp*(1.0_rp+s2)*t*t9                 !
             deriv(1,9)= 0.25_rp*(1.0_rp+s2)*t*t1                 !
             deriv(1,7)= 0.25_rp*(-1.0_rp+s2)*t*t1                !
             deriv(2,1)= 0.25_rp*(-1.0_rp+t2)*s*s9                !  1      5      2
             deriv(2,3)= 0.25_rp*s*s1*(-1.0_rp+t2)
             deriv(2,9)= 0.25_rp*s*s1*(1.0_rp+t2)
             deriv(2,7)= 0.25_rp*s*s9*(1.0_rp+t2)
          end if
          if (khes == 1) then
             heslo(1,1)= 0.5_rp*t*t9
             heslo(1,3)= 0.5_rp*t*t9
             heslo(1,9)= 0.5_rp*t*t1
             heslo(1,7)= 0.5_rp*t*t1
             heslo(2,1)= 0.5_rp*s*s9
             heslo(2,3)= 0.5_rp*s*s1
             heslo(2,9)= 0.5_rp*s*s1
             heslo(2,7)= 0.5_rp*s*s9
             heslo(3,1)= 0.25_rp*(-1.0_rp+t2)*(s9+s)
             heslo(3,3)= 0.25_rp*(-1.0_rp+t2)*(s1+s)
             heslo(3,9)= 0.25_rp*(1.0_rp+t2)*(s1+s)
             heslo(3,7)= 0.25_rp*(1.0_rp+t2)*(s9+s)
          end if
       else
          ! Linear (1 to 4)
          shape(1)=(1.0_rp-t-s+st)*0.25_rp                     !  4         3
          shape(2)=(1.0_rp-t+s-st)*0.25_rp                     !
          shape(3)=(1.0_rp+t+s+st)*0.25_rp                     !      
          shape(4)=(1.0_rp+t-s-st)*0.25_rp                     !
          deriv(1,1)=(-1.0_rp+t)*0.25_rp                       !  1         2
          deriv(1,2)=(+1.0_rp-t)*0.25_rp
          deriv(1,3)=(+1.0_rp+t)*0.25_rp
          deriv(1,4)=(-1.0_rp-t)*0.25_rp
          deriv(2,1)=(-1.0_rp+s)*0.25_rp
          deriv(2,2)=(-1.0_rp-s)*0.25_rp
          deriv(2,3)=(+1.0_rp+s)*0.25_rp
          deriv(2,4)=(+1.0_rp-s)*0.25_rp
          heslo(3,1)= 0.25_rp
          heslo(3,2)=-0.25_rp
          heslo(3,3)= 0.25_rp
          heslo(3,4)=-0.25_rp                  
       end if
       ! Quadratic (5 to 9)
       shape( 2)=0.5_rp*(1.0_rp-ss)*t*t9                    !       7  
       shape( 6)=0.5_rp*s*s1*(1.0_rp-tt)                    !   
       shape( 8)=0.5_rp*(1.0_rp-ss)*t*t1                    !  8    9    6
       shape( 4)=0.5_rp*s*s9*(1.0_rp-tt)                    !   
       shape( 5)=(1.0_rp-ss)*(1.0_rp-tt)                    !       5
       if ( kder == 1 ) then
          deriv(1,2)=-st*t9
          deriv(1,6)= 0.5_rp*(1.0_rp+s2)*(1.0_rp-tt)
          deriv(1,8)=-st*t1
          deriv(1,4)= 0.5_rp*(-1.0_rp+s2)*(1.0_rp-tt)
          deriv(1,5)=-s2*(1.0_rp-tt)
          deriv(2,2)= 0.5_rp*(1.0_rp-ss)*(-1.0_rp+t2)
          deriv(2,6)=-st*s1
          deriv(2,8)= 0.5_rp*(1.0_rp-ss)*(1.0_rp+t2)
          deriv(2,4)=-st*s9
          deriv(2,5)=-t2*(1.0_rp-ss)
       end if
       if ( khes == 1 ) then
          heslo(1,2)=-t*t9
          heslo(1,6)= 1.0_rp-tt
          heslo(1,8)=-t*t1
          heslo(1,4)= 1.0_rp-tt
          heslo(1,5)=-2.0_rp*(1.0_rp-tt)
          heslo(2,2)= 1.0_rp-ss
          heslo(2,6)=-s*s1
          heslo(2,8)= 1.0_rp-ss
          heslo(2,4)=-s*s9
          heslo(2,5)=-2.0_rp*(1.0_rp-ss)
          heslo(3,2)=-s*(-1.0_rp+t2)
          heslo(3,6)=-t*s1-st
          heslo(3,8)=-s*(1.0_rp+t2)
          heslo(3,4)=-t*s9-st
          heslo(3,5)= s2*t2
       end if
    else if(nnode==10) then
       c=9.0_rp/2.0_rp
       a1=1.0_rp-s-t
       a2=2.0_rp/3.0_rp-s-t
       a3=1.0_rp/3.0_rp-s-t
       shape( 1)=c*a1*a2*a3                                 !  10
       shape( 4)=c*(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s)*s    !
       shape(10)=c*(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t)*t    !
       shape( 2)= 3.0_rp*c*a1*a2*s                          !  8    9
       shape( 3)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-s)*s           !
       shape( 7)=-3.0_rp*c*(1.0_rp/3.0_rp-s)*s*t            !
       shape( 9)=-3.0_rp*c*s*(1.0_rp/3.0_rp-t)*t            !  5    6     7
       shape( 8)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-t)*t           !
       shape( 5)= 3.0_rp*c*a1*a2*t                          !
       shape( 6)= 6.0_rp*c*a1*s*t                           !  1    2    3    4
       if ( kder == 1 ) then
          deriv(1, 1)=-c*(a1*a2+a1*a3+a2*a3)       
          deriv(1, 4)=-c*((2.0_rp/3.0_rp-s)*s&
               + (1.0_rp/3.0_rp-s)*s-(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s))
          deriv(1, 10)=0.0_rp
          deriv(1, 2)= 3.0_rp*c*(a1*a2-a1*s-a2*s)
          deriv(1, 3)=-3.0_rp*c*(a1*(1.0_rp/3.0_rp-s)&
               - a1*s-(1.0_rp/3.0_rp-s)*s)
          deriv(1, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*t-s*t)
          deriv(1, 9)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
          deriv(1, 8)= 3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
          deriv(1, 5)= 3.0_rp*c*(-a1*t-a2*t)
          deriv(1,6)= 6.0_rp*c*(a1*t-s*t)
          deriv(2, 1)=-c*(a1*a2+a1*a3+a2*a3)
          deriv(2, 4)= 0.0_rp
          deriv(2, 10)=-c*((2.0_rp/3.0_rp-t)*t&
               + (1.0_rp/3.0_rp-t)*t-(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t))
          deriv(2, 2)= 3.0_rp*c*(-a1*s-a2*s)
          deriv(2, 3)=-3.0_rp*c*(-(1.0_rp/3.0_rp-s)*s)
          deriv(2, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*s)
          deriv(2, 9)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*s-s*t)
          deriv(2, 8)=-3.0_rp*c*(-(1.0_rp/3.0_rp-t)*t&
               - a1*t+a1*(1.0_rp/3.0_rp-t))
          deriv(2, 5)= 3.0_rp*c*(-a1*t-a2*t+a1*a2)
          deriv(2,6)= 6.0_rp*c*(a1*s-s*t)
       end if
       if ( khes == 1 ) then
          heslo(1, 1)= c*2.0_rp*(a1+a2+a3) 
          heslo(1, 4)= c*(6.0_rp*s-2.0_rp) 
          heslo(1, 10)= 0.0_rp 
          heslo(1, 2)= c*( 18.0_rp*s+12.0_rp*t-10.0_rp)
          heslo(1, 3)= c*(-18.0_rp*s- 6.0_rp*t+ 8.0_rp)
          heslo(1, 7)= c*6.0_rp*t 
          heslo(1, 9)= 0.0_rp 
          heslo(1, 8)= 0.0_rp  
          heslo(1, 5)= c*6.0_rp*t 
          heslo(1,6)=-c*12.0_rp*t 
          heslo(2, 1)= c*2.0_rp*(a1+a2+a3) 
          heslo(2, 4)= 0.0_rp 
          heslo(2, 10)= c*(6.0_rp*t-2.0_rp) 
          heslo(2, 2)= c*6.0_rp*s
          heslo(2, 3)= 0.0_rp
          heslo(2, 7)= 0.0_rp
          heslo(2, 9)= c*6.0_rp*s
          heslo(2, 8)= c*( -6.0_rp*s-18.0*t+ 8.0_rp)
          heslo(2, 5)= c*( 12.0_rp*s+18.0*t-10.0_rp)
          heslo(2,6)=-c*12.0_rp*s
          heslo(3, 1)= 2.0_rp*c*(a1+a2+a3) 
          heslo(3, 4)= 0.0_rp  
          heslo(3, 10)= 0.0_rp 
          heslo(3, 2)= c*( 12.0_rp*s+6.0_rp*t-5.0_rp) 
          heslo(3, 3)= c*(- 6.0_rp*s+1.0_rp) 
          heslo(3, 7)= c*(  6.0_rp*s-1.0_rp) 
          heslo(3, 9)= c*(  6.0_rp*t-1.0_rp) 
          heslo(3, 8)= c*(- 6.0_rp*t+1.0_rp)  
          heslo(3, 5)= c*(  6.0_rp*s+12.0_rp*t-5.0_rp) 
          heslo(3,6)= c*(-12.0_rp*s-12.0_rp*t+6.0_rp) 
       end if
    else if(nnode==16) then
       a =81.0_rp/256.0_rp
       c =1.0_rp/3.0_rp
       s1=1.0_rp+s
       s2=c+s
       s3=c-s
       s4=1.0_rp-s
       t1=1.0_rp+t
       t2=c+t
       t3=c-t
       t4=1.0_rp-t
       st=s*t
       if(khie) then
          ! Cubic (1 to 4)
          shape( 1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
          shape( 2) =   a*s1*s2*s3*t2*t3*t4                   ! 
          shape( 3) =   a*s1*s2*s3*t1*t2*t3                   ! 
          shape( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
          if (kder ==1 ) then
             deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)   !
             deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)   ! 12   13   14    7
             deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)   !
             deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)   !  
             deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4) ! 1     5    6    2
             deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)    
             deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
             deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
          end if
          if ( khes == 1 ) then
             heslo(1, 1) =&
                  a *t2*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
             heslo(1, 2) =&
                  a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
             heslo(1, 3) =&
                  a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
             heslo(1, 4) =&
                  a *t1*t2*t3*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
             heslo(2, 1) =&
                  a *s2*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
             heslo(2, 2) =&
                  a *s1*s2*s3*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
             heslo(2, 3) =&
                  a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
             heslo(2, 4) =&
                  a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
             heslo(3, 1) =&
                  a*(-s2*s3-s2*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
             heslo(3, 2) =&
                  a*(-s1*s2+s1*s3+s2*s3)*(-t2*t3-t2*t4+t3*t4)       
             heslo(3, 3) =&
                  a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t3+t2*t3)
             heslo(3, 4) =&
                  a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
          end if
       else
          ! Linear (1 to 4)
          shape(1)=(1.0_rp-t-s+st)*0.25_rp                     !  4         3
          shape(2)=(1.0_rp-t+s-st)*0.25_rp                     !
          shape(3)=(1.0_rp+t+s+st)*0.25_rp                     !      
          shape(4)=(1.0_rp+t-s-st)*0.25_rp                     !
          deriv(1,1)=(-1.0_rp+t)*0.25_rp                       !  1         2
          deriv(1,2)=(+1.0_rp-t)*0.25_rp
          deriv(1,3)=(+1.0_rp+t)*0.25_rp
          deriv(1,4)=(-1.0_rp-t)*0.25_rp
          deriv(2,1)=(-1.0_rp+s)*0.25_rp
          deriv(2,2)=(-1.0_rp-s)*0.25_rp
          deriv(2,3)=(+1.0_rp+s)*0.25_rp
          deriv(2,4)=(+1.0_rp-s)*0.25_rp
          heslo(3,1)= 0.25_rp
          heslo(3,2)=-0.25_rp
          heslo(3,3)= 0.25_rp
          heslo(3,4)=-0.25_rp
       end if
       ! Cubic (5 to 16)
       shape( 5) =-3.0_rp*a*s1*s3*s4*t2*t3*t4              !      10    9    
       shape( 6) =-3.0_rp*a*s1*s2*s4*t2*t3*t4              !
       shape( 7) =-3.0_rp*a*s1*s2*s3*t1*t3*t4              !
       shape( 8) =-3.0_rp*a*s1*s2*s3*t1*t2*t4              ! 11   16   15    8
       shape( 9) =-3.0_rp*a*s1*s2*s4*t1*t2*t3              !
       shape(10) =-3.0_rp*a*s1*s3*s4*t1*t2*t3              !
       shape(11) =-3.0_rp*a*s2*s3*s4*t1*t2*t4              ! 12   13   14    7   
       shape(12) =-3.0_rp*a*s2*s3*s4*t1*t3*t4              !
       shape(13) = 9.0_rp*a*s1*s3*s4*t1*t3*t4              !
       shape(14) = 9.0_rp*a*s1*s2*s4*t1*t3*t4              !       5    6    
       shape(15) = 9.0_rp*a*s1*s2*s4*t1*t2*t4
       shape(16) = 9.0_rp*a*s1*s3*s4*t1*t2*t4
       if ( kder == 1 ) then
          deriv(1, 5)=-3.0_rp*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(1, 6)=-3.0_rp*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1, 7)=-3.0_rp*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 8)=-3.0_rp*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
          deriv(1, 9)=-3.0_rp*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
          deriv(1,10)=-3.0_rp*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
          deriv(1,11)=-3.0_rp*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
          deriv(1,12)=-3.0_rp*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
          deriv(1,13)= 9.0_rp*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(1,14)= 9.0_rp*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1,15)= 9.0_rp*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
          deriv(1,16)= 9.0_rp*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
          deriv(2, 5)= -3.0_rp*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 6)= -3.0_rp*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
          deriv(2, 7)= -3.0_rp*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
          deriv(2, 8)= -3.0_rp*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
          deriv(2, 9)= -3.0_rp*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
          deriv(2,10)= -3.0_rp*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
          deriv(2,11)= -3.0_rp*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
          deriv(2,12)= -3.0_rp*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,13)=  9.0_rp*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,14)=  9.0_rp*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
          deriv(2,15)=  9.0_rp*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
          deriv(2,16)=  9.0_rp*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)
       end if
       if ( khes == 1 ) then
          heslo(1, 5) =&
               -3.0_rp*a *t2*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1, 6) =&
               -3.0_rp*a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1, 7) =&
               -3.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 8) =&
               -3.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
          heslo(1, 9) =&
               -3.0_rp*a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,10) =&
               -3.0_rp*a *t1*t2*t3*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1,11) =&
               -3.0_rp*a *t1*t2*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1,12) =&
               -3.0_rp*a *t1*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
          heslo(1,13) =&
               9.0_rp*a *t1*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(1,14) =&
               9.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,15) =&
               9.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
          heslo(1,16) =&
               9.0_rp*a *t1*t2*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
          heslo(2, 5) =&
               -3.0_rp*a *s1*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 6) =&
               -3.0_rp*a *s1*s2*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 7) =&
               -3.0_rp*a *s1*s2*s3*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2, 8) =&
               -3.0_rp*a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2, 9) =&
               -3.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2,10) =&
               -3.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
          heslo(2,11) =&
               -3.0_rp*a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2,12) =&
               -3.0_rp*a *s2*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,13) =&
               9.0_rp*a *s1*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,14) =&
               9.0_rp*a *s1*s2*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
          heslo(2,15) =&
               9.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(2,16) =&
               9.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
          heslo(3, 5) =&
               -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 6) =&
               -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t2*t3-t2*t4+t3*t4)
          heslo(3, 7) =&
               -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t3-t1*t4+t3*t4)
          heslo(3, 8) =&
               -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t4+t2*t4)
          heslo(3, 9) =&
               -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t3+t2*t3)
          heslo(3,10) =&
               -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
          heslo(3,11) =&
               -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
          heslo(3,12) =&
               -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,13) =&
               9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,14) =&
               9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t3-t1*t4+t3*t4)
          heslo(3,15) =&
               9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t4+t2*t4)
          heslo(3,16) =&
               9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
       end if
    end if
  end subroutine shape2

  !=========================================================================
  subroutine shape3(s,t,z,nnode,ntens,khie,kder,khes,shape,deriv,heslo)
    !-----------------------------------------------------------------------
    !
    ! This routine evaluates shape functions and their first and second
    ! derivatives 3-d standar continuous interpolation elements.
    ! 
    !    TETRAHEDRA:  4  10  &  20  nodes
    !    HEXAHEDRA:   8  27  &  64  nodes
    !    PRISM:       6             nodes
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)       , intent(in)  :: nnode,ntens,kder,khes
    logical(lg)       , intent(in)  :: khie
    real(rp)          , intent(in)  :: s,t,z
    real(rp)          , intent(out) :: shape(nnode)
    real(rp), optional, intent(out) :: deriv(3,nnode),heslo(ntens,nnode)
    integer(ip)              :: istop,i
    real(rp)                 :: a1,a2,a3,a4,a,p1,p2,p3,z1,z2,z3,z4,s1,s2,s3,s4
    real(rp)                 :: t1,t2,t3,t4,sm,tm,zm,sq,tp,zp,s11,s21,s31,s41
    real(rp)                 :: t11,t21,t31,t41,z11,z21,z31,s12,s22,s32,s42
    real(rp)                 :: t12,t22,t32,t42,z41,z12,z22,z32,z42,sl,tl,zl

    istop=0

    ! Linear tetrahedron 
    if(nnode==4) then
       shape(   1) = 1.0_rp-s-t-z
       shape(   2) = s
       shape(   3) = t
       shape(   4) = z
       if ( kder == 1 ) then
          deriv(1, 1) =-1.0_rp
          deriv(2, 1) =-1.0_rp
          deriv(3, 1) =-1.0_rp
          deriv(3, 4) = 1.0_rp
          deriv(1, 2) = 1.0_rp
          deriv(2, 3) = 1.0_rp
       end if

       ! Quadratic tetrahedron 
    else if(nnode==10) then
       a1= 1.0_rp-s-t-z
       a2=s
       a3=t
       a4=z
       shape(   1) = (2.0_rp*a1-1.0_rp)*a1
       deriv(1, 1) = 1.0_rp-4.0_rp*a1
       deriv(2, 1) = 1.0_rp-4.0_rp*a1
       deriv(3, 1) = 1.0_rp-4.0_rp*a1
       do i = 1,6
          heslo(i, 1) = 4.0_rp
       enddo
       shape(   3) = (2.0_rp*a2-1.0_rp)*a2
       deriv(1, 3) = 4.0_rp*a2-1.0_rp
       heslo(1, 3) = 4.0_rp
       shape(   6) = (2.0_rp*a3-1.0_rp)*a3
       deriv(2, 6) = 4.0_rp*a3-1.0_rp
       heslo(2, 6) = 4.0_rp
       shape(  10) = (2.0_rp*a4-1.0_rp)*a4
       deriv(3,10) = 4.0_rp*a4-1.0_rp
       heslo(3,10) = 4.0_rp
       shape(   2) = 4.0_rp*a1*a2
       deriv(1, 2) = 4.0_rp*(a1-a2)
       deriv(2, 2) =-4.0_rp*a2
       deriv(3, 2) =-4.0_rp*a2
       heslo(1, 2) =-8.0_rp
       heslo(4, 2) =-4.0_rp
       heslo(5, 2) =-4.0_rp
       shape(   5) = 4.0_rp*a2*a3
       deriv(1, 5) = 4.0_rp*a3
       deriv(2, 5) = 4.0_rp*a2
       heslo(4, 5) = 4.0_rp
       shape(   4) = 4.0_rp*a1*a3
       deriv(1, 4) =-4.0_rp*a3
       deriv(2, 4) = 4.0_rp*(a1-a3)
       deriv(3, 4) =-4.0_rp*a3
       heslo(2, 4) =-8.0_rp
       heslo(4, 4) =-4.0_rp
       heslo(6, 4) =-4.0_rp
       shape(   7) = 4.0_rp*a1*a4
       deriv(1, 7) =-4.0_rp*a4
       deriv(2, 7) =-4.0_rp*a4
       deriv(3, 7) = 4.0_rp*(a1-a4)
       heslo(3, 7) =-8.0_rp
       heslo(5, 7) =-4.0_rp
       heslo(6, 7) =-4.0_rp
       shape(   8) = 4.0_rp*a2*a4
       deriv(1, 8) = 4.0_rp*a4
       deriv(3, 8) = 4.0_rp*a2
       heslo(5, 8) = 4.0_rp
       shape(   9) = 4.0_rp*a3*a4
       deriv(2, 9) = 4.0_rp*a4
       deriv(3, 9) = 4.0_rp*a3
       heslo(6, 9) = 4.0_rp

       ! Cubic tetrahedron
    else if(nnode==20) then
       write(*,*) 'interpolation:: ERROR! Cubic tethraedron not implemented!'
       check(1 == 0)
       a=4.5_rp
       p1= 1-s-t-z
       p2= 2.0_rp/3.0_rp-s-t-z
       p3= 1.0_rp/3.0_rp-s-t-z
       z1= 2.0_rp/3.0_rp-z
       z2= 1.0_rp/3.0_rp-z
       s1= 2.0_rp/3.0_rp-s
       s2= 1.0_rp/3.0_rp-s
       t1= 2.0_rp/3.0_rp-t
       t2= 1.0_rp/3.0_rp-t

       shape(   1) = a*p1*p2*p3
       deriv(1, 1) =-a*(p1*p2+p1*p3+p2*p3)
       deriv(2, 1) =-a*(p1*p2+p1*p3+p2*p3)
       deriv(3, 1) =-a*(p1*p2+p1*p3+p2*p3)
       shape(   2) = a*z*z1*z2
       deriv(1, 2) = 0.0_rp
       deriv(2, 2) = 0.0_rp
       deriv(3, 2) =-a*(z*z1+z*z2-z1*z2)
       shape(   3) = a*s*s1*s2
       deriv(1, 3) =-a*(s*s1+s*s2-s1*s2)
       deriv(2, 3) = 0.0_rp
       deriv(3, 3) = 0.0_rp
       shape(   4) = a*t*t1*t2
       deriv(1, 4) = 0.0_rp
       deriv(2, 4) =-a*(t*t1+t*t2-t1*t2)
       deriv(3, 4) = 0.0_rp
       shape(   5) = 3.0_rp*a*p1*p2*t
       deriv(1, 5) = 3.0_rp*a*(-p2*t-p1*t)
       deriv(2, 5) = 3.0_rp*a*(-p2*t-p1*t+p1*p2)
       deriv(3, 5) = 3.0_rp*a*(-p2*t-p1*t)      
       shape(   6) = 3.0_rp*a*p1*p2*z
       deriv(1, 6) = 3.0_rp*a*(-p2*z-p1*z)      
       deriv(2, 6) = 3.0_rp*a*(-p2*z-p1*z)
       deriv(3, 6) = 3.0_rp*a*(-p2*z-p1*z+p1*p2)      
       shape(   7) = 3.0_rp*a*p1*p2*s
       deriv(1, 7) = 3.0_rp*a*(-p2*s-p1*s+p1*p2)
       deriv(2, 7) = 3.0_rp*a*(-p2*s-p1*s)
       deriv(3, 7) = 3.0_rp*a*(-p2*s-p1*s)
       shape(   8) =-3.0_rp*a*p1*t2*t
       deriv(1, 8) = 3.0_rp*a*(t2*t)
       deriv(2, 8) = 3.0_rp*a*(t2*t+p1*t-p1*t2)
       deriv(3, 8) = 3.0_rp*a*(t2*t)
       shape(   9) =-3.0_rp*a*p1*s2*s
       deriv(1, 9) = 3.0_rp*a*(s2*s+p1*s-p1*s2)
       deriv(2, 9) = 3.0_rp*a*(s2*s)
       deriv(3, 9) = 3.0_rp*a*(s2*s)
       shape(  10) =-3.0_rp*a*p1*z2*z
       deriv(1,10) = 3.0_rp*a*(z2*z)
       deriv(2,10) = 3.0_rp*a*(z2*z)
       deriv(3,10) = 3.0_rp*a*(z2*z+p1*z-p1*z2)
       shape(  11) =-3.0_rp*a*t2*t*z
       deriv(1,11) = 0.0_rp
       deriv(2,11) =-3.0_rp*a*(t2*z-t*z)
       deriv(3,11) =-3.0_rp*a*t2*t
       shape(  12) =-3.0_rp*a*z2*t*z
       deriv(1,12) = 0.0_rp
       deriv(2,12) =-3.0_rp*a*z2*z
       deriv(3,12) =-3.0_rp*a*(t*z2-t*z)
       shape(  13) =-3.0_rp*a*z2*s*z
       deriv(1,13) =-3.0_rp*a*z2*z
       deriv(2,13) = 0.0_rp
       deriv(3,13) =-3.0_rp*a*(s*z2-s*z)
       shape(  14) =-3.0_rp*a*s2*s*z
       deriv(1,14) =-3.0_rp*a*(s2*z-s*z)
       deriv(2,14) = 0.0_rp
       deriv(3,14) =-3.0_rp*a*s2*s
       shape(  15) =-3.0_rp*a*s2*s*t
       deriv(1,15) =-3.0_rp*a*(s2*t-s*t)
       deriv(2,15) =-3.0_rp*a*s2*s
       deriv(3,15) = 0.0_rp
       shape(  16) =-3.0_rp*a*t2*t*s
       deriv(1,16) =-3.0_rp*a*t2*t
       deriv(2,16) =-3.0_rp*a*(t2*s-s*t)
       deriv(3,16) = 0.0_rp
       shape(  17) = 27.0_rp*p1*t*z
       deriv(1,17) =-27.0_rp*t*z
       deriv(2,17) = 27.0_rp*(p1*z-t*z)
       deriv(3,17) = 27.0_rp*(p1*t-t*z)
       shape(  18) = 27.0_rp*p1*s*z
       deriv(1,18) = 27.0_rp*(p1*z-s*z)
       deriv(2,18) =-27.0_rp*s*z
       deriv(3,18) = 27.0_rp*(p1*s-s*z)
       shape(  19) = 27.0_rp*p1*s*t
       deriv(1,19) = 27.0_rp*(p1*t-s*t)
       deriv(2,19) = 27.0_rp*(p1*s-s*t)
       deriv(3,19) =-27.0_rp*s*t
       shape(  20) = 27.0_rp*s*t*z
       deriv(1,20) = 27.0_rp*z*t
       deriv(2,20) = 27.0_rp*s*z
       deriv(3,20) = 27.0_rp*s*t

       heslo(1, 1) = 2.0_rp*a*(p1+p2+p3)
       heslo(2, 1) = 2.0_rp*a*(p1+p2+p3)
       heslo(3, 1) = 2.0_rp*a*(p1+p2+p3)
       heslo(4, 1) = 2.0_rp*a*(p1+p2+p3)
       heslo(5, 1) = 2.0_rp*a*(p1+p2+p3)
       heslo(6, 1) = 2.0_rp*a*(p1+p2+p3)

       heslo(3, 2) = 2.0_rp*a*(z-z1-z2)

       heslo(1, 3) = 2.0_rp*a*(s-s1-s2)

       heslo(2, 4) = 2.0_rp*a*(t-t1-t2)

       heslo(1, 5) = 6.0_rp*a*t
       heslo(2, 5) = 6.0_rp*a*(-p1-p2+t)
       heslo(3, 5) = 6.0_rp*a*t
       heslo(4, 5) = 3.0_rp*a*(-p2-p1+2*t)
       heslo(5, 5) = 6.0_rp*a*t
       heslo(6, 5) = 3.0_rp*a*(-p2-p1+2*t)

       heslo(1, 6) = 6.0_rp*a*z
       heslo(2, 6) = 6.0_rp*a*z
       heslo(3, 6) = 6.0_rp*a*(-p1-p2+z)
       heslo(4, 6) = 6.0_rp*a*z
       heslo(5, 6) = 3.0_rp*a*(-p2-p1+2*z)
       heslo(6, 6) = 3.0_rp*a*(-p2-p1+2*z)

       heslo(1, 7) = 6.0_rp*a*(-p1-p2+s)
       heslo(2, 7) = 6.0_rp*a*s
       heslo(3, 7) = 6.0_rp*a*s
       heslo(4, 7) = 3.0_rp*a*(-p2-p1+2*s)
       heslo(5, 7) = 3.0_rp*a*(-p2-p1+2*s)
       heslo(6, 7) = 6.0_rp*a*s

       heslo(1, 8) = 0.0_rp
       heslo(2, 8) = 6.0_rp*a*(t2-t+p1)
       heslo(3, 8) = 0.0_rp
       heslo(4, 8) = 3.0_rp*a*(t2-t)
       heslo(5, 8) = 0.0_rp
       heslo(6, 8) = 3.0_rp*a*(t2-t)

       heslo(1, 9) = 6.0_rp*a*(s2-s+p1)
       heslo(2, 9) = 0.0_rp
       heslo(3, 9) = 0.0_rp
       heslo(4, 9) = 3.0_rp*a*(s2-s)
       heslo(5, 9) = 3.0_rp*a*(s2-s)
       heslo(6, 9) = 0.0_rp

       heslo(1,10) = 0.0_rp
       heslo(2,10) = 0.0_rp
       heslo(3,10) = 6.0_rp*a*(z2-z+p1)
       heslo(4,10) = 0.0_rp
       heslo(5,10) = 3.0_rp*a*(z2-z)           
       heslo(6,10) = 3.0_rp*a*(z2-z)

       heslo(1,11) = 0.0_rp
       heslo(2,11) = 6.0_rp*a*z
       heslo(3,11) = 0.0_rp
       heslo(4,11) = 0.0_rp
       heslo(5,11) = 0.0_rp
       heslo(6,11) = 3.0_rp*a*(t-t2)

       heslo(1,12) = 0.0_rp
       heslo(2,12) = 0.0_rp
       heslo(3,12) = 6.0_rp*a*t
       heslo(4,12) = 0.0_rp
       heslo(5,12) = 0.0_rp
       heslo(6,12) = 3.0_rp*a*(z-z2)

       heslo(1,13) = 0.0_rp
       heslo(2,13) = 0.0_rp
       heslo(3,13) = 6.0_rp*a*s
       heslo(4,13) = 0.0_rp
       heslo(5,13) = 3.0_rp*a*(z-z2)
       heslo(6,13) = 0.0_rp

       heslo(1,14) = 6.0_rp*a*z
       heslo(2,14) = 0.0_rp
       heslo(3,14) = 0.0_rp
       heslo(4,14) = 0.0_rp
       heslo(5,14) = 3.0_rp*a*(s-s2)
       heslo(6,14) = 0.0_rp

       heslo(1,15) = 6.0_rp*a*t
       heslo(2,15) = 0.0_rp
       heslo(3,15) = 0.0_rp
       heslo(4,15) = 3.0_rp*a*(s-s2)
       heslo(5,15) = 0.0_rp
       heslo(6,15) = 0.0_rp

       heslo(1,16) = 0.0_rp
       heslo(2,16) = 6.0_rp*a*s
       heslo(3,16) = 0.0_rp
       heslo(4,16) = 3.0_rp*a*(t-t2)
       heslo(5,16) = 0.0_rp
       heslo(6,16) = 0.0_rp

       heslo(1,17) = 0.0_rp
       heslo(2,17) =-54.0_rp*z
       heslo(3,17) =-54.0_rp*t
       heslo(4,17) =-27.0_rp*z
       heslo(5,17) =-27.0_rp*t
       heslo(6,17) = 27.0_rp*(p1-z-t)

       heslo(1,18) =-54.0_rp*z
       heslo(2,18) = 0.0_rp
       heslo(3,18) =-54.0_rp*s
       heslo(4,18) =-27.0_rp*z
       heslo(5,18) =-27.0_rp*s
       heslo(6,18) = 27.0_rp*(p1-z-s)

       heslo(1,19) =-54.0_rp*t
       heslo(2,19) =-54.0_rp*s
       heslo(3,19) = 0.0_rp
       heslo(4,19) = 27.0_rp*(p1-t-s)
       heslo(5,19) =-27.0_rp*t
       heslo(6,19) =-27.0_rp*s

       heslo(1,20) = 0.0_rp
       heslo(2,20) = 0.0_rp
       heslo(3,20) = 0.0_rp
       heslo(4,20) = 27.0_rp*z
       heslo(5,20) = 27.0_rp*t
       heslo(6,20) = 27.0_rp*s

       ! Trilinear brick 
!!$    else if(nnode==8) then
!!$       sm = 0.5_rp*(1.0_rp-s)
!!$       tm = 0.5_rp*(1.0_rp-t)
!!$       zm = 0.5_rp*(1.0_rp-z)
!!$       sq = 0.5_rp*(1.0_rp+s)
!!$       tp = 0.5_rp*(1.0_rp+t)
!!$       zp = 0.5_rp*(1.0_rp+z)
!!$       shape(   1) = sm*tm*zm
!!$       deriv(1, 1) =-0.5_rp*tm*zm
!!$       deriv(2, 1) =-0.5_rp*sm*zm
!!$       deriv(3, 1) =-0.5_rp*sm*tm
!!$       heslo(4, 1) = 0.25_rp*zm
!!$       heslo(5, 1) = 0.25_rp*tm
!!$       heslo(6, 1) = 0.25_rp*sm
!!$       shape(   2) = sq*tm*zm
!!$       deriv(1, 2) = 0.5_rp*tm*zm
!!$       deriv(2, 2) =-0.5_rp*sq*zm
!!$       deriv(3, 2) =-0.5_rp*sq*tm
!!$       heslo(4, 2) =-0.25_rp*zm
!!$       heslo(5, 2) =-0.25_rp*tm
!!$       heslo(6, 2) = 0.25_rp*sq
!!$       shape(   3) = sq*tp*zm
!!$       deriv(1, 3) = 0.5_rp*tp*zm
!!$       deriv(2, 3) = 0.5_rp*sq*zm
!!$       deriv(3, 3) =-0.5_rp*sq*tp
!!$       heslo(4, 3) = 0.25_rp*zm
!!$       heslo(5, 3) =-0.25_rp*tp
!!$       heslo(6, 3) =-0.25_rp*sq
!!$       shape(   4) = sm*tp*zm
!!$       deriv(1, 4) =-0.5_rp*tp*zm
!!$       deriv(2, 4) = 0.5_rp*sm*zm
!!$       deriv(3, 4) =-0.5_rp*sm*tp
!!$       heslo(4, 4) =-0.25_rp*zm
!!$       heslo(5, 4) = 0.25_rp*tp
!!$       heslo(6, 4) =-0.25_rp*sm
!!$       shape(   5) = sm*tm*zp
!!$       deriv(1, 5) =-0.5_rp*tm*zp
!!$       deriv(2, 5) =-0.5_rp*sm*zp
!!$       deriv(3, 5) = 0.5_rp*sm*tm
!!$       heslo(4, 5) = 0.25_rp*zp
!!$       heslo(5, 5) =-0.25_rp*tm
!!$       heslo(6, 5) =-0.25_rp*sm
!!$       shape(   6) = sq*tm*zp 
!!$       deriv(1, 6) = 0.5_rp*tm*zp
!!$       deriv(2, 6) =-0.5_rp*sq*zp
!!$       deriv(3, 6) = 0.5_rp*sq*tm
!!$       heslo(4, 6) =-0.25_rp*zp
!!$       heslo(5, 6) = 0.25_rp*tm
!!$       heslo(6, 6) =-0.25_rp*sq
!!$       shape(   7) = sq*tp*zp
!!$       deriv(1, 7) = 0.5_rp*tp*zp
!!$       deriv(2, 7) = 0.5_rp*sq*zp
!!$       deriv(3, 7) = 0.5_rp*sq*tp
!!$       heslo(4, 7) = 0.25_rp*zp
!!$       heslo(5, 7) = 0.25_rp*tp
!!$       heslo(6, 7) = 0.25_rp*sq
!!$       shape(   8) = sm*tp*zp
!!$       deriv(1, 8) =-0.5_rp*tp*zp
!!$       deriv(2, 8) = 0.5_rp*sm*zp
!!$       deriv(3, 8) = 0.5_rp*sm*tp
!!$       heslo(4, 8) =-0.25_rp*zp
!!$       heslo(5, 8) =-0.25_rp*tp
!!$       heslo(6, 8) = 0.25_rp*sm

 else if(nnode==8) then
       sm = 0.5_rp*(1.0_rp-s)
       tm = 0.5_rp*(1.0_rp-t)
       zm = 0.5_rp*(1.0_rp-z)
       sq = 0.5_rp*(1.0_rp+s)
       tp = 0.5_rp*(1.0_rp+t)
       zp = 0.5_rp*(1.0_rp+z)
       shape(   1) = sm*tm*zm
       shape(   2) = sq*tm*zm
       shape(   4) = sq*tp*zm
       shape(   3) = sm*tp*zm
       shape(   5) = sm*tm*zp
       shape(   6) = sq*tm*zp 
       shape(   8) = sq*tp*zp
       shape(   7) = sm*tp*zp

       if ( kder == 1 ) then 
          deriv(1, 1) =-0.5_rp*tm*zm
          deriv(2, 1) =-0.5_rp*sm*zm
          deriv(3, 1) =-0.5_rp*sm*tm
          deriv(1, 2) = 0.5_rp*tm*zm
          deriv(2, 2) =-0.5_rp*sq*zm
          deriv(3, 2) =-0.5_rp*sq*tm
          deriv(1, 4) = 0.5_rp*tp*zm
          deriv(2, 4) = 0.5_rp*sq*zm
          deriv(3, 4) =-0.5_rp*sq*tp
          deriv(1, 3) =-0.5_rp*tp*zm
          deriv(2, 3) = 0.5_rp*sm*zm
          deriv(3, 3) =-0.5_rp*sm*tp
          deriv(1, 5) =-0.5_rp*tm*zp
          deriv(2, 5) =-0.5_rp*sm*zp
          deriv(3, 5) = 0.5_rp*sm*tm
          deriv(1, 6) = 0.5_rp*tm*zp
          deriv(2, 6) =-0.5_rp*sq*zp
          deriv(3, 6) = 0.5_rp*sq*tm
          deriv(1, 8) = 0.5_rp*tp*zp
          deriv(2, 8) = 0.5_rp*sq*zp
          deriv(3, 8) = 0.5_rp*sq*tp
          deriv(1, 7) =-0.5_rp*tp*zp
          deriv(2, 7) = 0.5_rp*sm*zp
          deriv(3, 7) = 0.5_rp*sm*tp
       end if

       if (khes == 1 ) then
          heslo(4, 1) = 0.25_rp*zm
          heslo(5, 1) = 0.25_rp*tm
          heslo(6, 1) = 0.25_rp*sm
          heslo(4, 2) =-0.25_rp*zm
          heslo(5, 2) =-0.25_rp*tm
          heslo(6, 2) = 0.25_rp*sq
          heslo(4, 4) = 0.25_rp*zm
          heslo(5, 4) =-0.25_rp*tp
          heslo(6, 4) =-0.25_rp*sq
          heslo(4, 3) =-0.25_rp*zm
          heslo(5, 3) = 0.25_rp*tp
          heslo(6, 3) =-0.25_rp*sm
          heslo(4, 5) = 0.25_rp*zp
          heslo(5, 5) =-0.25_rp*tm
          heslo(6, 5) =-0.25_rp*sm
          heslo(4, 6) =-0.25_rp*zp
          heslo(5, 6) = 0.25_rp*tm
          heslo(6, 6) =-0.25_rp*sq
          heslo(4, 8) = 0.25_rp*zp
          heslo(5, 8) = 0.25_rp*tp
          heslo(6, 8) = 0.25_rp*sq
          heslo(4, 7) =-0.25_rp*zp
          heslo(5, 7) =-0.25_rp*tp
          heslo(6, 7) = 0.25_rp*sm
       end if
       ! Triquadratic brick
    else if(nnode==27) then
       sl=s*(s-1.0_rp)
       tl=t*(t-1.0_rp)
       zl=z*(z-1.0_rp)
       s1= 2.0_rp*s-1.0_rp
       t1= 2.0_rp*t-1.0_rp
       z1= 2.0_rp*z-1.0_rp
       s2= 1.0_rp-s*s
       t2= 1.0_rp-t*t
       z2= 1.0_rp-z*z
       s3= 1.0_rp+2.0_rp*s
       t3= 1.0_rp+2.0_rp*t
       z3= 1.0_rp+2.0_rp*z
       s4=-2.0_rp*s
       t4=-2.0_rp*t
       z4=-2.0_rp*z
       if(khie) then
          ! Quadratic nodes (1 to 8)
          sq=s*(s+1.0_rp)
          tp=t*(t+1.0_rp)
          zp=z*(z+1.0_rp)         
          shape(   1) = 0.125_rp*sl*tl*zl
          deriv(1, 1) = 0.125_rp*s1*tl*zl
          deriv(2, 1) = 0.125_rp*sl*t1*zl
          deriv(3, 1) = 0.125_rp*sl*tl*z1
          heslo(1, 1) = 0.25_rp*tl*zl
          heslo(2, 1) = 0.25_rp*sl*zl
          heslo(3, 1) = 0.25_rp*sl*tl
          heslo(4, 1) = 0.125_rp*s1*t1*zl
          heslo(5, 1) = 0.125_rp*s1*tl*z1
          heslo(6, 1) = 0.125_rp*sl*t1*z1
          shape(   2) = 0.125_rp*sq*tl*zl
          deriv(1, 2) = 0.125_rp*s3*tl*zl
          deriv(2, 2) = 0.125_rp*sq*t1*zl
          deriv(3, 2) = 0.125_rp*sq*tl*z1
          heslo(1, 2) = 0.25_rp*tl*zl
          heslo(2, 2) = 0.25_rp*sq*zl
          heslo(3, 2) = 0.25_rp*sq*tl
          heslo(4, 2) = 0.125_rp*s3*t1*zl
          heslo(5, 2) = 0.125_rp*s3*tl*z1
          heslo(6, 2) = 0.125_rp*sq*t1*z1
          shape(   3) = 0.125_rp*sq*tp*zl
          deriv(1, 3) = 0.125_rp*s3*tp*zl
          deriv(2, 3) = 0.125_rp*sq*t3*zl
          deriv(3, 3) = 0.125_rp*sq*tp*z1
          heslo(1, 3) = 0.25_rp*tp*zl
          heslo(2, 3) = 0.25_rp*sq*zl
          heslo(3, 3) = 0.25_rp*sq*tp
          heslo(4, 3) = 0.125_rp*s3*t3*zl
          heslo(5, 3) = 0.125_rp*s3*tp*z1
          heslo(6, 3) = 0.125_rp*sq*t3*z1
          shape(   4) = 0.125_rp*sl*tp*zl
          deriv(1, 4) = 0.125_rp*s1*tp*zl
          deriv(2, 4) = 0.125_rp*sl*t3*zl
          deriv(3, 4) = 0.125_rp*sl*tp*z1
          heslo(1, 4) = 0.25_rp*tp*zl
          heslo(2, 4) = 0.25_rp*sl*zl
          heslo(3, 4) = 0.25_rp*sl*tp
          heslo(4, 4) = 0.125_rp*s1*t3*zl
          heslo(5, 4) = 0.125_rp*s1*tp*z1
          heslo(6, 4) = 0.125_rp*sl*t3*z1
          shape(   5) = 0.125_rp*sl*tl*zp
          deriv(1, 5) = 0.125_rp*s1*tl*zp
          deriv(2, 5) = 0.125_rp*sl*t1*zp
          deriv(3, 5) = 0.125_rp*sl*tl*z3
          heslo(1, 5) = 0.25_rp*tl*zp
          heslo(2, 5) = 0.25_rp*sl*zp
          heslo(3, 5) = 0.25_rp*sl*tl
          heslo(4, 5) = 0.125_rp*s1*t1*zp
          heslo(5, 5) = 0.125_rp*s1*tl*z3
          heslo(6, 5) = 0.125_rp*sl*t1*z3
          shape(   6) = 0.125_rp*sq*tl*zp
          deriv(1, 6) = 0.125_rp*s3*tl*zp
          deriv(2, 6) = 0.125_rp*sq*t1*zp
          deriv(3, 6) = 0.125_rp*sq*tl*z3
          heslo(1, 6) = 0.25_rp*tl*zp
          heslo(2, 6) = 0.25_rp*sq*zp
          heslo(3, 6) = 0.25_rp*sq*tl
          heslo(4, 6) = 0.125_rp*s3*t1*zp
          heslo(5, 6) = 0.125_rp*s3*tl*z3
          heslo(6, 6) = 0.125_rp*sq*t1*z3
          shape(   7) = 0.125_rp*sq*tp*zp
          deriv(1, 7) = 0.125_rp*s3*tp*zp
          deriv(2, 7) = 0.125_rp*sq*t3*zp
          deriv(3, 7) = 0.125_rp*sq*tp*z3
          heslo(1, 7) = 0.25_rp*tp*zp
          heslo(2, 7) = 0.25_rp*sq*zp
          heslo(3, 7) = 0.25_rp*sq*tp
          heslo(4, 7) = 0.125_rp*s3*t3*zp
          heslo(5, 7) = 0.125_rp*s3*tp*z3
          heslo(6, 7) = 0.125_rp*sq*t3*z3
          shape(   8) = 0.125_rp*sl*tp*zp
          deriv(1, 8) = 0.125_rp*s1*tp*zp
          deriv(2, 8) = 0.125_rp*sl*t3*zp
          deriv(3, 8) = 0.125_rp*sl*tp*z3
          heslo(1, 8) = 0.25_rp*tp*zp
          heslo(2, 8) = 0.25_rp*sl*zp
          heslo(3, 8) = 0.25_rp*sl*tp
          heslo(4, 8) = 0.125_rp*s1*t3*zp
          heslo(5, 8) = 0.125_rp*s1*tp*z3
          heslo(6, 8) = 0.125_rp*sl*t3*z3
       else
          ! Linear nodes (1 to 8)
          sm = 0.5_rp*(1.0_rp-s)
          tm = 0.5_rp*(1.0_rp-t)
          zm = 0.5_rp*(1.0_rp-z)
          sq = 0.5_rp*(1.0_rp+s)
          tp = 0.5_rp*(1.0_rp+t)
          zp = 0.5_rp*(1.0_rp+z)
          shape(   1) = sm*tm*zm
          deriv(1, 1) =-0.5_rp*tm*zm
          deriv(2, 1) =-0.5_rp*sm*zm
          deriv(3, 1) =-0.5_rp*sm*tm
          heslo(4, 1) = 0.25_rp*zm
          heslo(5, 1) = 0.25_rp*tm
          heslo(6, 1) = 0.25_rp*sm
          shape(   2) = sq*tm*zm
          deriv(1, 2) = 0.5_rp*tm*zm
          deriv(2, 2) =-0.5_rp*sq*zm
          deriv(3, 2) =-0.5_rp*sq*tm
          heslo(4, 2) =-0.25_rp*zm
          heslo(5, 2) =-0.25_rp*tm
          heslo(6, 2) = 0.25_rp*sq
          shape(   3) = sq*tp*zm
          deriv(1, 3) = 0.5_rp*tp*zm
          deriv(2, 3) = 0.5_rp*sq*zm
          deriv(3, 3) =-0.5_rp*sq*tp
          heslo(4, 3) = 0.25_rp*zm
          heslo(5, 3) =-0.25_rp*tp
          heslo(6, 3) =-0.25_rp*sq
          shape(   4) = sm*tp*zm
          deriv(1, 4) =-0.5_rp*tp*zm
          deriv(2, 4) = 0.5_rp*sm*zm
          deriv(3, 4) =-0.5_rp*sm*tp
          heslo(4, 4) =-0.25_rp*zm
          heslo(5, 4) = 0.25_rp*tp
          heslo(6, 4) =-0.25_rp*sm
          shape(   5) = sm*tm*zp
          deriv(1, 5) =-0.5_rp*tm*zp
          deriv(2, 5) =-0.5_rp*sm*zp
          deriv(3, 5) = 0.5_rp*sm*tm
          heslo(4, 5) = 0.25_rp*zp
          heslo(5, 5) =-0.25_rp*tm
          heslo(6, 5) =-0.25_rp*sm
          shape(   6) = sq*tm*zp 
          deriv(1, 6) = 0.5_rp*tm*zp
          deriv(2, 6) =-0.5_rp*sq*zp
          deriv(3, 6) = 0.5_rp*sq*tm
          heslo(4, 6) =-0.25_rp*zp
          heslo(5, 6) = 0.25_rp*tm
          heslo(6, 6) =-0.25_rp*sq
          shape(   7) = sq*tp*zp
          deriv(1, 7) = 0.5_rp*tp*zp
          deriv(2, 7) = 0.5_rp*sq*zp
          deriv(3, 7) = 0.5_rp*sq*tp
          heslo(4, 7) = 0.25_rp*zp
          heslo(5, 7) = 0.25_rp*tp
          heslo(6, 7) = 0.25_rp*sq
          shape(   8) = sm*tp*zp
          deriv(1, 8) =-0.5_rp*tp*zp
          deriv(2, 8) = 0.5_rp*sm*zp
          deriv(3, 8) = 0.5_rp*sm*tp
          heslo(4, 8) =-0.25_rp*zp
          heslo(5, 8) =-0.25_rp*tp
          heslo(6, 8) = 0.25_rp*sm
       end if
       ! Quadratic nodes (9 to 27)
       sq=s*(s+1.0_rp)
       tp=t*(t+1.0_rp)
       zp=z*(z+1.0_rp)
       shape(   9) = 0.25_rp*s2*tl*zl
       deriv(1, 9) = 0.25_rp*s4*tl*zl
       deriv(2, 9) = 0.25_rp*s2*t1*zl
       deriv(3, 9) = 0.25_rp*s2*tl*z1
       heslo(1, 9) =-0.5_rp*tl*zl
       heslo(2, 9) = 0.5_rp*s2*zl
       heslo(3, 9) = 0.5_rp*s2*tl
       heslo(4, 9) = 0.25_rp*s4*t1*zl
       heslo(5, 9) = 0.25_rp*s4*tl*z1
       heslo(6, 9) = 0.25_rp*s2*t1*z1
       shape(  10) = 0.25_rp*sq*t2*zl
       deriv(1,10) = 0.25_rp*s3*t2*zl
       deriv(2,10) = 0.25_rp*sq*t4*zl
       deriv(3,10) = 0.25_rp*sq*t2*z1
       heslo(1,10) = 0.5_rp*t2*zl
       heslo(2,10) =-0.5_rp*sq*zl
       heslo(3,10) = 0.5_rp*sq*t2
       heslo(4,10) = 0.25_rp*s3*t4*zl
       heslo(5,10) = 0.25_rp*s3*t2*z1
       heslo(6,10) = 0.25_rp*sq*t4*z1
       shape(  11) = 0.25_rp*s2*tp*zl
       deriv(1,11) = 0.25_rp*s4*tp*zl
       deriv(2,11) = 0.25_rp*s2*t3*zl
       deriv(3,11) = 0.25_rp*s2*tp*z1
       heslo(1,11) =-0.5_rp*tp*zl
       heslo(2,11) = 0.5_rp*s2*zl
       heslo(3,11) = 0.5_rp*s2*tp
       heslo(4,11) = 0.25_rp*s4*t3*zl
       heslo(5,11) = 0.25_rp*s4*tp*z1
       heslo(6,11) = 0.25_rp*s2*t3*z1
       shape(  12) = 0.25_rp*sl*t2*zl
       deriv(1,12) = 0.25_rp*s1*t2*zl
       deriv(2,12) = 0.25_rp*sl*t4*zl
       deriv(3,12) = 0.25_rp*sl*t2*z1
       heslo(1,12) = 0.5_rp*t2*zl
       heslo(2,12) =-0.5_rp*sl*zl
       heslo(3,12) = 0.5_rp*sl*t2
       heslo(4,12) = 0.25_rp*s1*t4*zl
       heslo(5,12) = 0.25_rp*s1*t2*z1
       heslo(6,12) = 0.25_rp*sl*t4*z1
       shape(  13) = 0.25_rp*sl*tl*z2
       deriv(1,13) = 0.25_rp*s1*tl*z2
       deriv(2,13) = 0.25_rp*sl*t1*z2
       deriv(3,13) = 0.25_rp*sl*tl*z4
       heslo(1,13) = 0.5_rp*tl*z2
       heslo(2,13) = 0.5_rp*sl*z2
       heslo(3,13) =-0.5_rp*sl*tl
       heslo(4,13) = 0.25_rp*s1*t1*z2
       heslo(5,13) = 0.25_rp*s1*tl*z4
       heslo(6,13) = 0.25_rp*sl*t1*z4
       shape(  14) = 0.25_rp*sq*tl*z2
       deriv(1,14) = 0.25_rp*s3*tl*z2
       deriv(2,14) = 0.25_rp*sq*t1*z2
       deriv(3,14) = 0.25_rp*sq*tl*z4
       heslo(1,14) = 0.5_rp*tl*z2
       heslo(2,14) = 0.5_rp*sq*z2
       heslo(3,14) =-0.5_rp*sq*tl
       heslo(4,14) = 0.25_rp*s3*t1*z2
       heslo(5,14) = 0.25_rp*s3*tl*z4
       heslo(6,14) = 0.25_rp*sq*t1*z4
       shape(  15) = 0.25_rp*sq*tp*z2
       deriv(1,15) = 0.25_rp*s3*tp*z2
       deriv(2,15) = 0.25_rp*sq*t3*z2
       deriv(3,15) = 0.25_rp*sq*tp*z4
       heslo(1,15) = 0.5_rp*tp*z2
       heslo(2,15) = 0.5_rp*sq*z2
       heslo(3,15) =-0.5_rp*sq*tp
       heslo(4,15) = 0.25_rp*s3*t3*z2
       heslo(5,15) = 0.25_rp*s3*tp*z4
       heslo(6,15) = 0.25_rp*sq*t3*z4
       shape(  16) = 0.25_rp*sl*tp*z2
       deriv(1,16) = 0.25_rp*s1*tp*z2
       deriv(2,16) = 0.25_rp*sl*t3*z2
       deriv(3,16) = 0.25_rp*sl*tp*z4
       heslo(1,16) = 0.5_rp*tp*z2
       heslo(2,16) = 0.5_rp*sl*z2
       heslo(3,16) =-0.5_rp*sl*tp
       heslo(4,16) = 0.25_rp*s1*t3*z2
       heslo(5,16) = 0.25_rp*s1*tp*z4
       heslo(6,16) = 0.25_rp*sl*t3*z4
       shape(  17) = 0.25_rp*s2*tl*zp
       deriv(1,17) = 0.25_rp*s4*tl*zp
       deriv(2,17) = 0.25_rp*s2*t1*zp
       deriv(3,17) = 0.25_rp*s2*tl*z3
       heslo(1,17) =-0.5_rp*tl*zp
       heslo(2,17) = 0.5_rp*s2*zp
       heslo(3,17) = 0.5_rp*s2*tl
       heslo(4,17) = 0.25_rp*s4*t1*zp
       heslo(5,17) = 0.25_rp*s4*tl*z3
       heslo(6,17) = 0.25_rp*s2*t1*z3
       shape(  18) = 0.25_rp*sq*t2*zp
       deriv(1,18) = 0.25_rp*s3*t2*zp
       deriv(2,18) = 0.25_rp*sq*t4*zp
       deriv(3,18) = 0.25_rp*sq*t2*z3
       heslo(1,18) = 0.5_rp*t2*zp
       heslo(2,18) =-0.5_rp*sq*zp
       heslo(3,18) = 0.5_rp*sq*t2
       heslo(4,18) = 0.25_rp*s3*t4*zp
       heslo(5,18) = 0.25_rp*s3*t2*z3
       heslo(6,18) = 0.25_rp*sq*t4*z3
       shape(  19) = 0.25_rp*s2*tp*zp
       deriv(1,19) = 0.25_rp*s4*tp*zp
       deriv(2,19) = 0.25_rp*s2*t3*zp
       deriv(3,19) = 0.25_rp*s2*tp*z3
       heslo(1,19) =-0.5_rp*tp*zp
       heslo(2,19) = 0.5_rp*s2*zp
       heslo(3,19) = 0.5_rp*s2*tp
       heslo(4,19) = 0.25_rp*s4*t3*zp
       heslo(5,19) = 0.25_rp*s4*tp*z3
       heslo(6,19) = 0.25_rp*s2*t3*z3
       shape(  20) = 0.25_rp*sl*t2*zp
       deriv(1,20) = 0.25_rp*s1*t2*zp
       deriv(2,20) = 0.25_rp*sl*t4*zp
       deriv(3,20) = 0.25_rp*sl*t2*z3
       heslo(1,20) = 0.5_rp*t2*zp
       heslo(2,20) =-0.5_rp*sl*zp
       heslo(3,20) = 0.5_rp*sl*t2
       heslo(4,20) = 0.25_rp*s1*t4*zp
       heslo(5,20) = 0.25_rp*s1*t2*z3
       heslo(6,20) = 0.25_rp*sl*t4*z3
       shape(  21) = 0.5_rp*s2*t2*zl
       deriv(1,21) = 0.5_rp*s4*t2*zl
       deriv(2,21) = 0.5_rp*s2*t4*zl
       deriv(3,21) = 0.5_rp*s2*t2*z1
       heslo(1,21) =-t2*zl
       heslo(2,21) =-s2*zl
       heslo(3,21) = s2*t2
       heslo(4,21) = 0.5_rp*s4*t4*zl
       heslo(5,21) = 0.5_rp*s4*t2*z1
       heslo(6,21) = 0.5_rp*s2*t4*z1
       shape(  22) = 0.5_rp*s2*tl*z2
       deriv(1,22) = 0.5_rp*s4*tl*z2
       deriv(2,22) = 0.5_rp*s2*t1*z2
       deriv(3,22) = 0.5_rp*s2*tl*z4
       heslo(1,22) =-tl*z2
       heslo(2,22) = s2*z2
       heslo(3,22) =-s2*tl
       heslo(4,22) = 0.5_rp*s4*t1*z2
       heslo(5,22) = 0.5_rp*s4*tl*z4
       heslo(6,22) = 0.5_rp*s2*t1*z4
       shape(  23) = 0.5_rp*sq*t2*z2
       deriv(1,23) = 0.5_rp*s3*t2*z2
       deriv(2,23) = 0.5_rp*sq*t4*z2
       deriv(3,23) = 0.5_rp*sq*t2*z4
       heslo(1,23) = t2*z2
       heslo(2,23) =-sq*z2
       heslo(3,23) =-sq*t2
       heslo(4,23) = 0.5_rp*s3*t4*z2
       heslo(5,23) = 0.5_rp*s3*t2*z4
       heslo(6,23) = 0.5_rp*sq*t4*z4
       shape(  24) = 0.5_rp*s2*tp*z2
       deriv(1,24) = 0.5_rp*s4*tp*z2
       deriv(2,24) = 0.5_rp*s2*t3*z2
       deriv(3,24) = 0.5_rp*s2*tp*z4
       heslo(1,24) =-tp*z2
       heslo(2,24) = s2*z2
       heslo(3,24) =-s2*tp
       heslo(4,24) = 0.5_rp*s4*t3*z2
       heslo(5,24) = 0.5_rp*s4*tp*z4
       heslo(6,24) = 0.5_rp*s2*t3*z4
       shape(  25) = 0.5_rp*sl*t2*z2
       deriv(1,25) = 0.5_rp*s1*t2*z2
       deriv(2,25) = 0.5_rp*sl*t4*z2
       deriv(3,25) = 0.5_rp*sl*t2*z4
       heslo(1,25) = t2*z2
       heslo(2,25) =-sl*z2
       heslo(3,25) =-sl*t2
       heslo(4,25) = 0.5_rp*s1*t4*z2
       heslo(5,25) = 0.5_rp*s1*t2*z4
       heslo(6,25) = 0.5_rp*sl*t4*z4
       shape(  26) = 0.5_rp*s2*t2*zp
       deriv(1,26) = 0.5_rp*s4*t2*zp
       deriv(2,26) = 0.5_rp*s2*t4*zp
       deriv(3,26) = 0.5_rp*s2*t2*z3
       heslo(1,26) =-t2*zp
       heslo(2,26) =-s2*zp
       heslo(3,26) = s2*t2
       heslo(4,26) = 0.5_rp*s4*t4*zp
       heslo(5,26) = 0.5_rp*s4*t2*z3
       heslo(6,26) = 0.5_rp*s2*t4*z3
       shape(  27) = s2*t2*z2
       deriv(1,27) = s4*t2*z2
       deriv(2,27) = s2*t4*z2
       deriv(3,27) = s2*t2*z4
       heslo(1,27) =-2.0_rp*t2*z2
       heslo(2,27) =-2.0_rp*s2*z2
       heslo(3,27) =-2.0_rp*s2*t2
       heslo(4,27) = s4*t4*z2
       heslo(5,27) = s4*t2*z4
       heslo(6,27) = s2*t4*z4

       ! Tricubic brick
    else if (nnode==64) then
       a=729.0_rp/4096.0_rp                             
       s1=-(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)*(1.0_rp-s)         
       s2=(1.0_rp+s)*(1.0_rp/3.0_rp-s)*(1.0_rp-s)
       s3=(1.0_rp+s)*(1.0_rp/3.0_rp+s)*(1.0_rp-s)                         
       s4=-(1.0_rp+s)*(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)
       t1=-(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)*(1.0_rp-t)         
       t2=(1.0_rp+t)*(1.0_rp/3.0_rp-t)*(1.0_rp-t)
       t3=(1.0_rp+t)*(1.0_rp/3.0_rp+t)*(1.0_rp-t)                         
       t4=-(1.0_rp+t)*(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)
       z1=-(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)*(1.0_rp-z)         
       z2=(1.0_rp+z)*(1.0_rp/3.0_rp-z)*(1.0_rp-z)
       z3=(1.0_rp+z)*(1.0_rp/3.0_rp+z)*(1.0_rp-z)                         
       z4=-(1.0_rp+z)*(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)
       s11=-((1.0_rp/3.0_rp-s)*(1.0_rp-s)-(1.0_rp/3.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s))
       s21=(1.0_rp/3.0_rp-s)*(1.0_rp-s)-(1.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp-s)    
       s31=(1.0_rp/3.0_rp+s)*(1.0_rp-s)+(1.0_rp+s)*(1.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp+s)    
       s41=-((1.0_rp/3.0_rp+s)*(1.0_rp/3.0_rp-s)+(1.0_rp+s)*(1.0_rp/3.0_rp-s)&
            -(1.0_rp+s)*(1.0_rp/3.0_rp+s))
       t11=-((1.0_rp/3.0_rp-t)*(1.0_rp-t)-(1.0_rp/3.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t))
       t21=(1.0_rp/3.0_rp-t)*(1.0_rp-t)-(1.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp-t)    
       t31=(1.0_rp/3.0_rp+t)*(1.0_rp-t)+(1.0_rp+t)*(1.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp+t)    
       t41=-((1.0_rp/3.0_rp+t)*(1.0_rp/3.0_rp-t)+(1.0_rp+t)*(1.0_rp/3.0_rp-t)&
            -(1.0_rp+t)*(1.0_rp/3.0_rp+t))
       z11=-((1.0_rp/3.0_rp-z)*(1.0_rp-z)-(1.0_rp/3.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z))
       z21=(1.0_rp/3.0_rp-z)*(1.0_rp-z)-(1.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp-z)    
       z31=(1.0_rp/3.0_rp+z)*(1.0_rp-z)+(1.0_rp+z)*(1.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp+z)    
       z41=-((1.0_rp/3.0_rp+z)*(1.0_rp/3.0_rp-z)+(1.0_rp+z)*(1.0_rp/3.0_rp-z)&
            -(1.0_rp+z)*(1.0_rp/3.0_rp+z))
       s12= 2.0_rp*((1.0_rp-s)+(1.0_rp/3.0_rp-s)-(1.0_rp/3.0_rp+s))
       s22=-2.0_rp*((1.0_rp-s)+(1.0_rp/3.0_rp-s)-(1.0_rp+s))
       s32=-2.0_rp*((1.0_rp-s)-(1.0_rp/3.0_rp+s)-(1.0_rp+s))
       s42= 2.0_rp*((1.0_rp/3.0_rp-s)-(1.0_rp/3.0_rp+s)-(1.0_rp+s))
       t12= 2.0_rp*((1.0_rp-t)+(1.0_rp/3.0_rp-t)-(1.0_rp/3.0_rp+t))
       t22=-2.0_rp*((1.0_rp-t)+(1.0_rp/3.0_rp-t)-(1.0_rp+t)) 
       t32=-2.0_rp*((1.0_rp-t)-(1.0_rp/3.0_rp+t)-(1.0_rp+t))  
       t42= 2.0_rp*((1.0_rp/3.0_rp-t)-(1.0_rp/3.0_rp+t)-(1.0_rp+t))
       z12= 2.0_rp*((1.0_rp-z)+(1.0_rp/3.0_rp-z)-(1.0_rp/3.0_rp+z))
       z22=-2.0_rp*((1.0_rp-z)+(1.0_rp/3.0_rp-z)-(1.0_rp+z)) 
       z32=-2.0_rp*((1.0_rp-z)-(1.0_rp/3.0_rp+z)-(1.0_rp+z))  
       z42= 2.0_rp*((1.0_rp/3.0_rp-z)-(1.0_rp/3.0_rp+z)-(1.0_rp+z))
       if(khie) then
          ! Cubic nodes (1 to 8)
          shape(   1) =   a*s1*t1*z1
          deriv(1, 1) =   a*s11*t1*z1
          deriv(2, 1) =   a*s1*t11*z1
          deriv(3, 1) =   a*s1*t1*z11
          shape(   2) =   a*s4*t1*z1
          deriv(1, 2) =   a*s41*t1*z1
          deriv(2, 2) =   a*s4*t11*z1
          deriv(3, 2) =   a*s4*t1*z11
          shape(   3) =   a*s4*t4*z1
          deriv(1, 3) =   a*s41*t4*z1
          deriv(2, 3) =   a*s4*t41*z1
          deriv(3, 3) =   a*s4*t4*z11
          shape(   4) =   a*s1*t4*z1
          deriv(1, 4) =   a*s11*t4*z1
          deriv(2, 4) =   a*s1*t41*z1
          deriv(3, 4) =   a*s1*t4*z11
          shape(   5) =   a*s1*t1*z4
          deriv(1, 5) =   a*s11*t1*z4
          deriv(2, 5) =   a*s1*t11*z4
          deriv(3, 5) =   a*s1*t1*z41
          shape(   6) =   a*s4*t1*z4
          deriv(1, 6) =   a*s41*t1*z4
          deriv(2, 6) =   a*s4*t11*z4
          deriv(3, 6) =   a*s4*t1*z41
          shape(   7) =   a*s4*t4*z4
          deriv(1, 7) =   a*s41*t4*z4
          deriv(2, 7) =   a*s4*t41*z4
          deriv(3, 7) =   a*s4*t4*z41
          shape(   8) =   a*s1*t4*z4
          deriv(1, 8) =   a*s11*t4*z4
          deriv(2, 8) =   a*s1*t41*z4
          deriv(3, 8) =   a*s1*t4*z41
          heslo(1, 1) =   a*s12*t1*z1
          heslo(2, 1) =   a*s1*t12*z1
          heslo(3, 1) =   a*s1*t1*z12
          heslo(4, 1) =   a*s11*t11*z1
          heslo(5, 1) =   a*s11*t1*z11
          heslo(6, 1) =   a*s1*t11*z11
          heslo(1, 2) =   a*s42*t1*z1  
          heslo(2, 2) =   a*s4*t12*z1 
          heslo(3, 2) =   a*s4*t1*z12 
          heslo(4, 2) =   a*s41*t11*z1
          heslo(5, 2) =   a*s41*t1*z11
          heslo(6, 2) =   a*s4*t11*z11          
          heslo(1, 3) =   a*s42*t4*z1           
          heslo(2, 3) =   a*s4*t42*z1  
          heslo(3, 3) =   a*s4*t4*z12      
          heslo(4, 3) =   a*s41*t41*z1
          heslo(5, 3) =   a*s41*t4*z11
          heslo(6, 3) =   a*s4*t41*z11      
          heslo(1, 4) =   a*s12*t4*z1 
          heslo(2, 4) =   a*s1*t42*z1 
          heslo(3, 4) =   a*s1*t4*z12 
          heslo(4, 4) =   a*s11*t41*z1
          heslo(5, 4) =   a*s11*t4*z11
          heslo(6, 4) =   a*s1*t41*z11
          heslo(1, 5) =   a*s12*t1*z4 
          heslo(2, 5) =   a*s1*t12*z4 
          heslo(3, 5) =   a*s1*t1*z42 
          heslo(4, 5) =   a*s11*t11*z4
          heslo(5, 5) =   a*s11*t1*z41
          heslo(6, 5) =   a*s1*t11*z41
          heslo(1, 6) =   a*s42*t4*z4 
          heslo(2, 6) =   a*s4*t42*z4 
          heslo(3, 6) =   a*s4*t4*z42 
          heslo(4, 6) =   a*s41*t41*z4
          heslo(5, 6) =   a*s41*t4*z41
          heslo(6, 6) =   a*s4*t41*z41
          heslo(1, 7) =   a*s42*t4*z4 
          heslo(2, 7) =   a*s4*t42*z4 
          heslo(3, 7) =   a*s4*t4*z42 
          heslo(4, 7) =   a*s41*t41*z4
          heslo(5, 7) =   a*s41*t4*z41
          heslo(6, 7) =   a*s4*t41*z41
          heslo(1, 8) =   a*s12*t4*z4 
          heslo(2, 8) =   a*s1*t42*z4 
          heslo(3, 8) =   a*s1*t4*z42 
          heslo(4, 8) =   a*s11*t41*z4
          heslo(5, 8) =   a*s11*t4*z41
          heslo(6, 8) =   a*s1*t41*z41
       else
          ! Linear nodes (1 to 8)
          sm = 0.5_rp*(1.0_rp-s)
          tm = 0.5_rp*(1.0_rp-t)
          zm = 0.5_rp*(1.0_rp-z)
          sq = 0.5_rp*(1.0_rp+s)
          tp = 0.5_rp*(1.0_rp+t)
          zp = 0.5_rp*(1.0_rp+z)
          shape(   1) = sm*tm*zm
          deriv(1, 1) =-0.5_rp*tm*zm
          deriv(2, 1) =-0.5_rp*sm*zm
          deriv(3, 1) =-0.5_rp*sm*tm
          heslo(4, 1) = 0.25_rp*zm
          heslo(5, 1) = 0.25_rp*tm
          heslo(6, 1) = 0.25_rp*sm
          shape(   2) = sq*tm*zm
          deriv(1, 2) = 0.5_rp*tm*zm
          deriv(2, 2) =-0.5_rp*sq*zm
          deriv(3, 2) =-0.5_rp*sq*tm
          heslo(4, 2) =-0.25_rp*zm
          heslo(5, 2) =-0.25_rp*tm
          heslo(6, 2) = 0.25_rp*sq
          shape(   3) = sq*tp*zm
          deriv(1, 3) = 0.5_rp*tp*zm
          deriv(2, 3) = 0.5_rp*sq*zm
          deriv(3, 3) =-0.5_rp*sq*tp
          heslo(4, 3) = 0.25_rp*zm
          heslo(5, 3) =-0.25_rp*tp
          heslo(6, 3) =-0.25_rp*sq
          shape(   4) = sm*tp*zm
          deriv(1, 4) =-0.5_rp*tp*zm
          deriv(2, 4) = 0.5_rp*sm*zm
          deriv(3, 4) =-0.5_rp*sm*tp
          heslo(4, 4) =-0.25_rp*zm
          heslo(5, 4) = 0.25_rp*tp
          heslo(6, 4) =-0.25_rp*sm
          shape(   5) = sm*tm*zp
          deriv(1, 5) =-0.5_rp*tm*zp
          deriv(2, 5) =-0.5_rp*sm*zp
          deriv(3, 5) = 0.5_rp*sm*tm
          heslo(4, 5) = 0.25_rp*zp
          heslo(5, 5) =-0.25_rp*tm
          heslo(6, 5) =-0.25_rp*sm
          shape(   6) = sq*tm*zp 
          deriv(1, 6) = 0.5_rp*tm*zp
          deriv(2, 6) =-0.5_rp*sq*zp
          deriv(3, 6) = 0.5_rp*sq*tm
          heslo(4, 6) =-0.25_rp*zp
          heslo(5, 6) = 0.25_rp*tm
          heslo(6, 6) =-0.25_rp*sq
          shape(   7) = sq*tp*zp
          deriv(1, 7) = 0.5_rp*tp*zp
          deriv(2, 7) = 0.5_rp*sq*zp
          deriv(3, 7) = 0.5_rp*sq*tp
          heslo(4, 7) = 0.25_rp*zp
          heslo(5, 7) = 0.25_rp*tp
          heslo(6, 7) = 0.25_rp*sq
          shape(   8) = sm*tp*zp
          deriv(1, 8) =-0.5_rp*tp*zp
          deriv(2, 8) = 0.5_rp*sm*zp
          deriv(3, 8) = 0.5_rp*sm*tp
          heslo(4, 8) =-0.25_rp*zp
          heslo(5, 8) =-0.25_rp*tp
          heslo(6, 8) = 0.25_rp*sm
       end if
       ! Cubic nodes (9 to 64)
       shape(   9) = 3.0_rp*a*s2*t1*z1
       deriv(1, 9) = 3.0_rp*a*s21*t1*z1
       deriv(2, 9) = 3.0_rp*a*s2*t11*z1
       deriv(3, 9) = 3.0_rp*a*s2*t1*z11
       shape(  10) = 3.0_rp*a*s3*t1*z1
       deriv(1,10) = 3.0_rp*a*s31*t1*z1
       deriv(2,10) = 3.0_rp*a*s3*t11*z1
       deriv(3,10) = 3.0_rp*a*s3*t1*z11        
       shape(  11) = 3.0_rp*a*s4*t2*z1
       deriv(1,11) = 3.0_rp*a*s41*t2*z1
       deriv(2,11) = 3.0_rp*a*s4*t21*z1
       deriv(3,11) = 3.0_rp*a*s4*t2*z11        
       shape(  12) = 3.0_rp*a*s4*t3*z1
       deriv(1,12) = 3.0_rp*a*s41*t3*z1        
       deriv(2,12) = 3.0_rp*a*s4*t31*z1        
       deriv(3,12) = 3.0_rp*a*s4*t3*z11                
       shape(  13) = 3.0_rp*a*s3*t4*z1
       deriv(1,13) = 3.0_rp*a*s31*t4*z1
       deriv(2,13) = 3.0_rp*a*s3*t41*z1
       deriv(3,13) = 3.0_rp*a*s3*t4*z11        
       shape(  14) = 3.0_rp*a*s2*t4*z1
       deriv(1,14) = 3.0_rp*a*s21*t4*z1
       deriv(2,14) = 3.0_rp*a*s2*t41*z1
       deriv(3,14) = 3.0_rp*a*s2*t4*z11        
       shape(  15) = 3.0_rp*a*s1*t3*z1
       deriv(1,15) = 3.0_rp*a*s11*t3*z1
       deriv(2,15) = 3.0_rp*a*s1*t31*z1
       deriv(3,15) = 3.0_rp*a*s1*t3*z11        
       shape(  16) = 3.0_rp*a*s1*t2*z1
       deriv(1,16) = 3.0_rp*a*s11*t2*z1
       deriv(2,16) = 3.0_rp*a*s1*t21*z1
       deriv(3,16) = 3.0_rp*a*s1*t2*z11        
       shape(  17) = 3.0_rp*a*s1*t1*z2
       deriv(1,17) = 3.0_rp*a*s11*t1*z2
       deriv(2,17) = 3.0_rp*a*s1*t11*z2
       deriv(3,17) = 3.0_rp*a*s1*t1*z21        
       shape(  18) = 3.0_rp*a*s4*t1*z2
       deriv(1,18) = 3.0_rp*a*s41*t1*z2
       deriv(2,18) = 3.0_rp*a*s4*t11*z2
       deriv(3,18) = 3.0_rp*a*s4*t1*z21        
       shape(  19) = 3.0_rp*a*s4*t4*z2
       deriv(1,19) = 3.0_rp*a*s41*t4*z2
       deriv(2,19) = 3.0_rp*a*s4*t41*z2
       deriv(3,19) = 3.0_rp*a*s4*t4*z21        
       shape(  20) = 3.0_rp*a*s1*t4*z2
       deriv(1,20) = 3.0_rp*a*s11*t4*z2
       deriv(2,20) = 3.0_rp*a*s1*t41*z2
       deriv(3,20) = 3.0_rp*a*s1*t4*z21        
       shape(  21) = 3.0_rp*a*s1*t1*z3
       deriv(1,21) = 3.0_rp*a*s11*t1*z3
       deriv(2,21) = 3.0_rp*a*s1*t11*z3
       deriv(3,21) = 3.0_rp*a*s1*t1*z31        
       shape(  22) = 3.0_rp*a*s4*t1*z3
       deriv(1,22) = 3.0_rp*a*s41*t1*z3
       deriv(2,22) = 3.0_rp*a*s4*t11*z3
       deriv(3,22) = 3.0_rp*a*s4*t1*z31
       shape(  23) = 3.0_rp*a*s4*t4*z3
       deriv(1,23) = 3.0_rp*a*s41*t4*z3
       deriv(2,23) = 3.0_rp*a*s4*t41*z3
       deriv(3,23) = 3.0_rp*a*s4*t4*z31        
       shape(  24) = 3.0_rp*a*s1*t4*z3
       deriv(1,24) = 3.0_rp*a*s11*t4*z3
       deriv(2,24) = 3.0_rp*a*s1*t41*z3
       deriv(3,24) = 3.0_rp*a*s1*t4*z31        
       shape(  25) = 3.0_rp*a*s2*t1*z4
       deriv(1,25) = 3.0_rp*a*s21*t1*z4
       deriv(2,25) = 3.0_rp*a*s2*t11*z4
       deriv(3,25) = 3.0_rp*a*s2*t1*z41        
       shape(  26) = 3.0_rp*a*s3*t1*z4
       deriv(1,26) = 3.0_rp*a*s31*t1*z4
       deriv(2,26) = 3.0_rp*a*s3*t11*z4
       deriv(3,26) = 3.0_rp*a*s3*t1*z41        
       shape(  27) = 3.0_rp*a*s4*t2*z4
       deriv(1,27) = 3.0_rp*a*s41*t2*z4
       deriv(2,27) = 3.0_rp*a*s4*t21*z4
       deriv(3,27) = 3.0_rp*a*s4*t2*z41        
       shape(  28) = 3.0_rp*a*s4*t3*z4
       deriv(1,28) = 3.0_rp*a*s41*t3*z4
       deriv(2,28) = 3.0_rp*a*s4*t31*z4
       deriv(3,28) = 3.0_rp*a*s4*t3*z41        
       shape(  29) = 3.0_rp*a*s3*t4*z4
       deriv(1,29) = 3.0_rp*a*s31*t4*z4
       deriv(2,29) = 3.0_rp*a*s3*t41*z4
       deriv(3,29) = 3.0_rp*a*s3*t4*z41        
       shape(  30) = 3.0_rp*a*s2*t4*z4
       deriv(1,30) = 3.0_rp*a*s21*t4*z4  
       deriv(2,30) = 3.0_rp*a*s2*t41*z4  
       deriv(3,30) = 3.0_rp*a*s2*t4*z41          
       shape(  31) = 3.0_rp*a*s1*t3*z4
       deriv(1,31) = 3.0_rp*a*s11*t3*z4
       deriv(2,31) = 3.0_rp*a*s1*t31*z4
       deriv(3,31) = 3.0_rp*a*s1*t3*z41        
       shape(  32) = 3.0_rp*a*s1*t2*z4
       deriv(1,32) = 3.0_rp*a*s11*t2*z4
       deriv(2,32) = 3.0_rp*a*s1*t21*z4
       deriv(3,32) = 3.0_rp*a*s1*t2*z41        
       shape(  33) = 9.0_rp*a*s2*t2*z1
       deriv(1,33) = 9.0_rp*a*s21*t2*z1
       deriv(2,33) = 9.0_rp*a*s2*t21*z1
       deriv(3,33) = 9.0_rp*a*s2*t2*z11        
       shape(  34) = 9.0_rp*a*s3*t2*z1
       deriv(1,34) = 9.0_rp*a*s31*t2*z1
       deriv(2,34) = 9.0_rp*a*s3*t21*z1
       deriv(3,34) = 9.0_rp*a*s3*t2*z11        
       shape(  35) = 9.0_rp*a*s3*t3*z1
       deriv(1,35) = 9.0_rp*a*s31*t3*z1
       deriv(2,35) = 9.0_rp*a*s3*t31*z1
       deriv(3,35) = 9.0_rp*a*s3*t3*z11        
       shape(  36) = 9.0_rp*a*s2*t3*z1
       deriv(1,36) = 9.0_rp*a*s21*t3*z1
       deriv(2,36) = 9.0_rp*a*s2*t31*z1
       deriv(3,36) = 9.0_rp*a*s2*t3*z11        
       shape(  37) = 9.0_rp*a*s2*t1*z2
       deriv(1,37) = 9.0_rp*a*s21*t1*z2
       deriv(2,37) = 9.0_rp*a*s2*t11*z2
       deriv(3,37) = 9.0_rp*a*s2*t1*z21        
       shape(  38) = 9.0_rp*a*s3*t1*z2
       deriv(1,38) = 9.0_rp*a*s31*t1*z2        
       deriv(2,38) = 9.0_rp*a*s3*t11*z2        
       deriv(3,38) = 9.0_rp*a*s3*t1*z21                
       shape(  39) = 9.0_rp*a*s4*t2*z2
       deriv(1,39) = 9.0_rp*a*s41*t2*z2        
       deriv(2,39) = 9.0_rp*a*s4*t21*z2        
       deriv(3,39) = 9.0_rp*a*s4*t2*z21                
       shape(  40) = 9.0_rp*a*s4*t3*z2
       deriv(1,40) = 9.0_rp*a*s41*t3*z2        
       deriv(2,40) = 9.0_rp*a*s4*t31*z2        
       deriv(3,40) = 9.0_rp*a*s4*t3*z21                
       shape(  41) = 9.0_rp*a*s3*t4*z2
       deriv(1,41) = 9.0_rp*a*s31*t4*z2        
       deriv(2,41) = 9.0_rp*a*s3*t41*z2        
       deriv(3,41) = 9.0_rp*a*s3*t4*z21                
       shape(  42) = 9.0_rp*a*s2*t4*z2
       deriv(1,42) = 9.0_rp*a*s21*t4*z2             
       deriv(2,42) = 9.0_rp*a*s2*t41*z2             
       deriv(3,42) = 9.0_rp*a*s2*t4*z21                     
       shape(  43) = 9.0_rp*a*s1*t3*z2
       deriv(1,43) = 9.0_rp*a*s11*t3*z2             
       deriv(2,43) = 9.0_rp*a*s1*t31*z2             
       deriv(3,43) = 9.0_rp*a*s1*t3*z21                     
       shape(  44) = 9.0_rp*a*s1*t2*z2
       deriv(1,44) = 9.0_rp*a*s11*t2*z2             
       deriv(2,44) = 9.0_rp*a*s1*t21*z2             
       deriv(3,44) = 9.0_rp*a*s1*t2*z21                     
       shape(  45) = 9.0_rp*a*s2*t1*z3
       deriv(1,45) = 9.0_rp*a*s21*t1*z3             
       deriv(2,45) = 9.0_rp*a*s2*t11*z3             
       deriv(3,45) = 9.0_rp*a*s2*t1*z31                     
       shape(  46) = 9.0_rp*a*s3*t1*z3
       deriv(1,46) = 9.0_rp*a*s31*t1*z3
       deriv(2,46) = 9.0_rp*a*s3*t11*z3
       deriv(3,46) = 9.0_rp*a*s3*t1*z31        
       shape(  47) = 9.0_rp*a*s4*t2*z3
       deriv(1,47) = 9.0_rp*a*s41*t2*z3
       deriv(2,47) = 9.0_rp*a*s4*t21*z3
       deriv(3,47) = 9.0_rp*a*s4*t2*z31
       shape(  48) = 9.0_rp*a*s4*t3*z3
       deriv(1,48) = 9.0_rp*a*s41*t3*z3
       deriv(2,48) = 9.0_rp*a*s4*t31*z3
       deriv(3,48) = 9.0_rp*a*s4*t3*z31        
       shape(  49) = 9.0_rp*a*s3*t4*z3
       deriv(1,49) = 9.0_rp*a*s31*t4*z3
       deriv(2,49) = 9.0_rp*a*s3*t41*z3
       deriv(3,49) = 9.0_rp*a*s3*t4*z31        
       shape(  50) = 9.0_rp*a*s2*t4*z3
       deriv(1,50) = 9.0_rp*a*s21*t4*z3
       deriv(2,50) = 9.0_rp*a*s2*t41*z3
       deriv(3,50) = 9.0_rp*a*s2*t4*z31        
       shape(  51) = 9.0_rp*a*s1*t3*z3
       deriv(1,51) = 9.0_rp*a*s11*t3*z3
       deriv(2,51) = 9.0_rp*a*s1*t31*z3
       deriv(3,51) = 9.0_rp*a*s1*t3*z31        
       shape(  52) = 9.0_rp*a*s1*t2*z3
       deriv(1,52) = 9.0_rp*a*s11*t2*z3
       deriv(2,52) = 9.0_rp*a*s1*t21*z3
       deriv(3,52) = 9.0_rp*a*s1*t2*z31        
       shape(  53) = 9.0_rp*a*s2*t2*z4
       deriv(1,53) = 9.0_rp*a*s21*t2*z4
       deriv(2,53) = 9.0_rp*a*s2*t21*z4
       deriv(3,53) = 9.0_rp*a*s2*t2*z41        
       shape(  54) = 9.0_rp*a*s3*t2*z4
       deriv(1,54) = 9.0_rp*a*s31*t2*z4
       deriv(2,54) = 9.0_rp*a*s3*t21*z4
       deriv(3,54) = 9.0_rp*a*s3*t2*z41        
       shape(  55) = 9.0_rp*a*s3*t3*z4
       deriv(1,55) = 9.0_rp*a*s31*t3*z4
       deriv(2,55) = 9.0_rp*a*s3*t31*z4
       deriv(3,55) = 9.0_rp*a*s3*t3*z41        
       shape(  56) = 9.0_rp*a*s2*t3*z4
       deriv(1,56) = 9.0_rp*a*s21*t3*z4
       deriv(2,56) = 9.0_rp*a*s2*t31*z4
       deriv(3,56) = 9.0_rp*a*s2*t3*z41        
       shape(  57) = 27.0_rp*a*s2*t2*z2
       deriv(1,57) = 27.0_rp*a*s21*t2*z2
       deriv(2,57) = 27.0_rp*a*s2*t21*z2
       deriv(3,57) = 27.0_rp*a*s2*t2*z21        
       shape(  58) = 27.0_rp*a*s3*t2*z2
       deriv(1,58) = 27.0_rp*a*s31*t2*z2
       deriv(2,58) = 27.0_rp*a*s3*t21*z2
       deriv(3,58) = 27.0_rp*a*s3*t2*z21        
       shape(  59) = 27.0_rp*a*s3*t3*z2
       deriv(1,59) = 27.0_rp*a*s31*t3*z2
       deriv(2,59) = 27.0_rp*a*s3*t31*z2
       deriv(3,59) = 27.0_rp*a*s3*t3*z21        
       shape(  60) = 27.0_rp*a*s2*t3*z2
       deriv(1,60) = 27.0_rp*a*s21*t3*z2
       deriv(2,60) = 27.0_rp*a*s2*t31*z2
       deriv(3,60) = 27.0_rp*a*s2*t3*z21        
       shape(  61) = 27.0_rp*a*s2*t2*z3
       deriv(1,61) = 27.0_rp*a*s21*t2*z3
       deriv(2,61) = 27.0_rp*a*s2*t21*z3
       deriv(3,61) = 27.0_rp*a*s2*t2*z31        
       shape(  62) = 27.0_rp*a*s3*t2*z3
       deriv(1,62) = 27.0_rp*a*s31*t2*z3
       deriv(2,62) = 27.0_rp*a*s3*t21*z3
       deriv(3,62) = 27.0_rp*a*s3*t2*z31        
       shape(  63) = 27.0_rp*a*s3*t3*z3
       deriv(1,63) = 27.0_rp*a*s31*t3*z3
       deriv(2,63) = 27.0_rp*a*s3*t31*z3
       deriv(3,63) = 27.0_rp*a*s3*t3*z31        
       shape(  64) = 27.0_rp*a*s2*t3*z3
       deriv(1,64) = 27.0_rp*a*s21*t3*z3
       deriv(2,64) = 27.0_rp*a*s2*t31*z3
       deriv(3,64) = 27.0_rp*a*s2*t3*z31

       heslo(1, 9) = 3.0_rp*a*s22*t1*z1 
       heslo(2, 9) = 3.0_rp*a*s2*t12*z1 
       heslo(3, 9) = 3.0_rp*a*s2*t1*z12 
       heslo(4, 9) = 3.0_rp*a*s21*t11*z1
       heslo(5, 9) = 3.0_rp*a*s21*t1*z11
       heslo(6, 9) = 3.0_rp*a*s2*t11*z11

       heslo(1,10) = 3.0_rp*a*s32*t1*z1 
       heslo(2,10) = 3.0_rp*a*s3*t12*z1 
       heslo(3,10) = 3.0_rp*a*s3*t1*z12 
       heslo(4,10) = 3.0_rp*a*s31*t11*z1
       heslo(5,10) = 3.0_rp*a*s31*t1*z11
       heslo(6,10) = 3.0_rp*a*s3*t11*z11

       heslo(1,11) = 3.0_rp*a*s42*t2*z1 
       heslo(2,11) = 3.0_rp*a*s4*t22*z1 
       heslo(3,11) = 3.0_rp*a*s4*t2*z12 
       heslo(4,11) = 3.0_rp*a*s41*t21*z1
       heslo(5,11) = 3.0_rp*a*s41*t2*z11
       heslo(6,11) = 3.0_rp*a*s4*t21*z11

       heslo(1,12) = 3.0_rp*a*s42*t3*z1 
       heslo(2,12) = 3.0_rp*a*s4*t32*z1 
       heslo(3,12) = 3.0_rp*a*s4*t3*z12 
       heslo(4,12) = 3.0_rp*a*s41*t31*z1
       heslo(5,12) = 3.0_rp*a*s41*t3*z11
       heslo(6,12) = 3.0_rp*a*s4*t31*z11

       heslo(1,13) = 3.0_rp*a*s32*t4*z1 
       heslo(2,13) = 3.0_rp*a*s3*t42*z1 
       heslo(3,13) = 3.0_rp*a*s3*t4*z12 
       heslo(4,13) = 3.0_rp*a*s31*t41*z1
       heslo(5,13) = 3.0_rp*a*s31*t4*z11
       heslo(6,13) = 3.0_rp*a*s3*t41*z11

       heslo(1,14) = 3.0_rp*a*s22*t4*z1 
       heslo(2,14) = 3.0_rp*a*s2*t42*z1 
       heslo(3,14) = 3.0_rp*a*s2*t4*z12 
       heslo(4,14) = 3.0_rp*a*s21*t41*z1
       heslo(5,14) = 3.0_rp*a*s21*t4*z11
       heslo(6,14) = 3.0_rp*a*s2*t41*z11

       heslo(1,15) = 3.0_rp*a*s12*t3*z1 
       heslo(2,15) = 3.0_rp*a*s1*t32*z1 
       heslo(3,15) = 3.0_rp*a*s1*t3*z12 
       heslo(4,15) = 3.0_rp*a*s11*t31*z1
       heslo(5,15) = 3.0_rp*a*s11*t3*z11
       heslo(6,15) = 3.0_rp*a*s1*t31*z11

       heslo(1,16) = 3.0_rp*a*s12*t2*z1 
       heslo(2,16) = 3.0_rp*a*s1*t22*z1 
       heslo(3,16) = 3.0_rp*a*s1*t2*z12 
       heslo(4,16) = 3.0_rp*a*s11*t21*z1
       heslo(5,16) = 3.0_rp*a*s11*t2*z11
       heslo(6,16) = 3.0_rp*a*s1*t21*z11

       heslo(1,17) = 3.0_rp*a*s12*t1*z2 
       heslo(2,17) = 3.0_rp*a*s1*t12*z2 
       heslo(3,17) = 3.0_rp*a*s1*t1*z22 
       heslo(4,17) = 3.0_rp*a*s11*t11*z2
       heslo(5,17) = 3.0_rp*a*s11*t1*z21
       heslo(6,17) = 3.0_rp*a*s1*t11*z21

       heslo(1,18) = 3.0_rp*a*s42*t1*z2 
       heslo(2,18) = 3.0_rp*a*s4*t12*z2 
       heslo(3,18) = 3.0_rp*a*s4*t1*z22 
       heslo(4,18) = 3.0_rp*a*s41*t11*z2
       heslo(5,18) = 3.0_rp*a*s41*t1*z21
       heslo(6,18) = 3.0_rp*a*s4*t11*z21

       heslo(1,19) = 3.0_rp*a*s42*t4*z2 
       heslo(2,19) = 3.0_rp*a*s4*t42*z2 
       heslo(3,19) = 3.0_rp*a*s4*t4*z22 
       heslo(4,19) = 3.0_rp*a*s41*t41*z2
       heslo(5,19) = 3.0_rp*a*s41*t4*z21
       heslo(6,19) = 3.0_rp*a*s4*t41*z21

       heslo(1,20) = 3.0_rp*a*s12*t4*z2 
       heslo(2,20) = 3.0_rp*a*s1*t42*z2 
       heslo(3,20) = 3.0_rp*a*s1*t4*z22 
       heslo(4,20) = 3.0_rp*a*s11*t41*z2
       heslo(5,20) = 3.0_rp*a*s11*t4*z21
       heslo(6,20) = 3.0_rp*a*s1*t41*z21

       heslo(1,21) = 3.0_rp*a*s12*t1*z3 
       heslo(2,21) = 3.0_rp*a*s1*t12*z3 
       heslo(3,21) = 3.0_rp*a*s1*t1*z32 
       heslo(4,21) = 3.0_rp*a*s11*t11*z3
       heslo(5,21) = 3.0_rp*a*s11*t1*z31
       heslo(6,21) = 3.0_rp*a*s1*t11*z31

       heslo(1,22) = 3.0_rp*a*s42*t1*z3 
       heslo(2,22) = 3.0_rp*a*s4*t12*z3 
       heslo(3,22) = 3.0_rp*a*s4*t1*z32 
       heslo(4,22) = 3.0_rp*a*s41*t11*z3
       heslo(5,22) = 3.0_rp*a*s41*t1*z31
       heslo(6,22) = 3.0_rp*a*s4*t11*z31

       heslo(1,23) = 3.0_rp*a*s42*t4*z3 
       heslo(2,23) = 3.0_rp*a*s4*t42*z3 
       heslo(3,23) = 3.0_rp*a*s4*t4*z32 
       heslo(4,23) = 3.0_rp*a*s41*t41*z3
       heslo(5,23) = 3.0_rp*a*s41*t4*z31
       heslo(6,23) = 3.0_rp*a*s4*t41*z31

       heslo(1,24) = 3.0_rp*a*s12*t4*z3 
       heslo(2,24) = 3.0_rp*a*s1*t42*z3 
       heslo(3,24) = 3.0_rp*a*s1*t4*z32 
       heslo(4,24) = 3.0_rp*a*s11*t41*z3
       heslo(5,24) = 3.0_rp*a*s11*t4*z31
       heslo(6,24) = 3.0_rp*a*s1*t41*z31

       heslo(1,25) = 3.0_rp*a*s22*t4*z4 
       heslo(2,25) = 3.0_rp*a*s2*t42*z4 
       heslo(3,25) = 3.0_rp*a*s2*t4*z42 
       heslo(4,25) = 3.0_rp*a*s21*t41*z4
       heslo(5,25) = 3.0_rp*a*s21*t4*z41
       heslo(6,25) = 3.0_rp*a*s2*t41*z41

       heslo(1,26) = 3.0_rp*a*s32*t1*z4 
       heslo(2,26) = 3.0_rp*a*s3*t12*z4 
       heslo(3,26) = 3.0_rp*a*s3*t1*z42 
       heslo(4,26) = 3.0_rp*a*s31*t11*z4
       heslo(5,26) = 3.0_rp*a*s31*t1*z41
       heslo(6,26) = 3.0_rp*a*s3*t11*z41

       heslo(1,27) = 3.0_rp*a*s42*t2*z4 
       heslo(2,27) = 3.0_rp*a*s4*t22*z4 
       heslo(3,27) = 3.0_rp*a*s4*t2*z42 
       heslo(4,27) = 3.0_rp*a*s41*t21*z4
       heslo(5,27) = 3.0_rp*a*s41*t2*z41
       heslo(6,27) = 3.0_rp*a*s4*t21*z41

       heslo(1,28) = 3.0_rp*a*s42*t3*z4 
       heslo(2,28) = 3.0_rp*a*s4*t32*z4 
       heslo(3,28) = 3.0_rp*a*s4*t3*z42 
       heslo(4,28) = 3.0_rp*a*s41*t31*z4
       heslo(5,28) = 3.0_rp*a*s41*t3*z41
       heslo(6,28) = 3.0_rp*a*s4*t31*z41

       heslo(1,29) = 3.0_rp*a*s32*t4*z4 
       heslo(2,29) = 3.0_rp*a*s3*t42*z4 
       heslo(3,29) = 3.0_rp*a*s3*t4*z42 
       heslo(4,29) = 3.0_rp*a*s31*t41*z4
       heslo(5,29) = 3.0_rp*a*s31*t4*z41
       heslo(6,29) = 3.0_rp*a*s3*t41*z41

       heslo(1,30) = 3.0_rp*a*s22*t4*z4 
       heslo(2,30) = 3.0_rp*a*s2*t42*z4 
       heslo(3,30) = 3.0_rp*a*s2*t4*z42 
       heslo(4,30) = 3.0_rp*a*s21*t41*z4
       heslo(5,30) = 3.0_rp*a*s21*t4*z41
       heslo(6,30) = 3.0_rp*a*s2*t41*z31

       heslo(1,31) = 3.0_rp*a*s12*t3*z4 
       heslo(2,31) = 3.0_rp*a*s1*t32*z4 
       heslo(3,31) = 3.0_rp*a*s1*t3*z42 
       heslo(4,31) = 3.0_rp*a*s11*t31*z4
       heslo(5,31) = 3.0_rp*a*s11*t3*z41
       heslo(6,31) = 3.0_rp*a*s1*t31*z41

       heslo(1,32) = 3.0_rp*a*s12*t2*z4 
       heslo(2,32) = 3.0_rp*a*s1*t22*z4 
       heslo(3,32) = 3.0_rp*a*s1*t2*z42 
       heslo(4,32) = 3.0_rp*a*s11*t21*z4
       heslo(5,32) = 3.0_rp*a*s11*t2*z41
       heslo(6,32) = 3.0_rp*a*s1*t21*z41

       heslo(1,33) = 9.0_rp*a*s22*t2*z1 
       heslo(2,33) = 9.0_rp*a*s2*t22*z1 
       heslo(3,33) = 9.0_rp*a*s2*t2*z12 
       heslo(4,33) = 9.0_rp*a*s21*t21*z1
       heslo(5,33) = 9.0_rp*a*s21*t2*z11
       heslo(6,33) = 9.0_rp*a*s2*t21*z11

       heslo(1,34) = 9.0_rp*a*s32*t2*z1 
       heslo(2,34) = 9.0_rp*a*s3*t22*z1 
       heslo(3,34) = 9.0_rp*a*s3*t2*z12 
       heslo(4,34) = 9.0_rp*a*s31*t21*z1
       heslo(5,34) = 9.0_rp*a*s31*t2*z11
       heslo(6,34) = 9.0_rp*a*s3*t21*z11

       heslo(1,35) = 9.0_rp*a*s32*t3*z1 
       heslo(2,35) = 9.0_rp*a*s3*t32*z1 
       heslo(3,35) = 9.0_rp*a*s3*t3*z12 
       heslo(4,35) = 9.0_rp*a*s31*t31*z1
       heslo(5,35) = 9.0_rp*a*s31*t3*z11
       heslo(6,35) = 9.0_rp*a*s3*t31*z11

       heslo(1,36) = 9.0_rp*a*s22*t3*z1 
       heslo(2,36) = 9.0_rp*a*s2*t32*z1 
       heslo(3,36) = 9.0_rp*a*s2*t3*z12 
       heslo(4,36) = 9.0_rp*a*s21*t31*z1
       heslo(5,36) = 9.0_rp*a*s21*t3*z11
       heslo(6,36) = 9.0_rp*a*s2*t31*z11

       heslo(1,37) = 9.0_rp*a*s22*t1*z2 
       heslo(2,37) = 9.0_rp*a*s2*t12*z2 
       heslo(3,37) = 9.0_rp*a*s2*t1*z22 
       heslo(4,37) = 9.0_rp*a*s21*t11*z2
       heslo(5,37) = 9.0_rp*a*s21*t1*z21
       heslo(6,37) = 9.0_rp*a*s2*t11*z21

       heslo(1,38) = 9.0_rp*a*s32*t1*z2 
       heslo(2,38) = 9.0_rp*a*s3*t12*z2 
       heslo(3,38) = 9.0_rp*a*s3*t1*z22 
       heslo(4,38) = 9.0_rp*a*s31*t11*z2
       heslo(5,38) = 9.0_rp*a*s31*t1*z21
       heslo(6,38) = 9.0_rp*a*s3*t11*z21

       heslo(1,39) = 9.0_rp*a*s42*t2*z2 
       heslo(2,39) = 9.0_rp*a*s4*t22*z2 
       heslo(3,39) = 9.0_rp*a*s4*t2*z22 
       heslo(4,39) = 9.0_rp*a*s41*t21*z2
       heslo(5,39) = 9.0_rp*a*s41*t2*z21
       heslo(6,39) = 9.0_rp*a*s4*t21*z21

       heslo(1,40) = 9.0_rp*a*s42*t3*z2 
       heslo(2,40) = 9.0_rp*a*s4*t32*z2 
       heslo(3,40) = 9.0_rp*a*s4*t3*z22 
       heslo(4,40) = 9.0_rp*a*s41*t31*z2
       heslo(5,40) = 9.0_rp*a*s41*t3*z21
       heslo(6,40) = 9.0_rp*a*s4*t31*z21

       heslo(1,41) = 9.0_rp*a*s32*t4*z2 
       heslo(2,41) = 9.0_rp*a*s3*t42*z2 
       heslo(3,41) = 9.0_rp*a*s3*t4*z22 
       heslo(4,41) = 9.0_rp*a*s31*t41*z2
       heslo(5,41) = 9.0_rp*a*s31*t4*z21
       heslo(6,41) = 9.0_rp*a*s3*t41*z21

       heslo(1,42) = 9.0_rp*a*s22*t4*z2 
       heslo(2,42) = 9.0_rp*a*s2*t42*z2 
       heslo(3,42) = 9.0_rp*a*s2*t4*z22 
       heslo(4,42) = 9.0_rp*a*s21*t41*z2
       heslo(5,42) = 9.0_rp*a*s21*t4*z21
       heslo(6,42) = 9.0_rp*a*s2*t41*z21

       heslo(1,43) = 9.0_rp*a*s12*t3*z2 
       heslo(2,43) = 9.0_rp*a*s1*t32*z2 
       heslo(3,43) = 9.0_rp*a*s1*t3*z22 
       heslo(4,43) = 9.0_rp*a*s11*t31*z2
       heslo(5,43) = 9.0_rp*a*s11*t3*z21
       heslo(6,43) = 9.0_rp*a*s1*t31*z21

       heslo(1,44) = 9.0_rp*a*s12*t2*z2 
       heslo(2,44) = 9.0_rp*a*s1*t22*z2 
       heslo(3,44) = 9.0_rp*a*s1*t2*z22 
       heslo(4,44) = 9.0_rp*a*s11*t21*z2
       heslo(5,44) = 9.0_rp*a*s11*t2*z21
       heslo(6,44) = 9.0_rp*a*s1*t21*z21

       heslo(1,45) = 9.0_rp*a*s22*t1*z3 
       heslo(2,45) = 9.0_rp*a*s2*t12*z3 
       heslo(3,45) = 9.0_rp*a*s2*t1*z32 
       heslo(4,45) = 9.0_rp*a*s21*t11*z3
       heslo(5,45) = 9.0_rp*a*s21*t1*z31
       heslo(6,45) = 9.0_rp*a*s2*t11*z31

       heslo(1,46) = 9.0_rp*a*s32*t1*z3 
       heslo(2,46) = 9.0_rp*a*s3*t12*z3 
       heslo(3,46) = 9.0_rp*a*s3*t1*z32 
       heslo(4,46) = 9.0_rp*a*s31*t11*z3
       heslo(5,46) = 9.0_rp*a*s31*t1*z31
       heslo(6,46) = 9.0_rp*a*s3*t11*z31

       heslo(1,47) = 9.0_rp*a*s42*t2*z3 
       heslo(2,47) = 9.0_rp*a*s4*t22*z3 
       heslo(3,47) = 9.0_rp*a*s4*t2*z32 
       heslo(4,47) = 9.0_rp*a*s41*t21*z3
       heslo(5,47) = 9.0_rp*a*s41*t2*z31
       heslo(6,47) = 9.0_rp*a*s4*t21*z31

       heslo(1,48) = 9.0_rp*a*s42*t3*z3 
       heslo(2,48) = 9.0_rp*a*s4*t32*z3 
       heslo(3,48) = 9.0_rp*a*s4*t3*z32 
       heslo(4,48) = 9.0_rp*a*s41*t31*z3
       heslo(5,48) = 9.0_rp*a*s41*t3*z31
       heslo(6,48) = 9.0_rp*a*s4*t31*z31

       heslo(1,49) = 9.0_rp*a*s32*t4*z3 
       heslo(2,49) = 9.0_rp*a*s3*t42*z3 
       heslo(3,49) = 9.0_rp*a*s3*t4*z32 
       heslo(4,49) = 9.0_rp*a*s31*t41*z3
       heslo(5,49) = 9.0_rp*a*s31*t4*z31
       heslo(6,49) = 9.0_rp*a*s3*t41*z31

       heslo(1,50) = 9.0_rp*a*s22*t4*z3 
       heslo(2,50) = 9.0_rp*a*s2*t42*z3 
       heslo(3,50) = 9.0_rp*a*s2*t4*z32 
       heslo(4,50) = 9.0_rp*a*s21*t41*z3
       heslo(5,50) = 9.0_rp*a*s21*t4*z31
       heslo(6,50) = 9.0_rp*a*s2*t41*z31

       heslo(1,51) = 9.0_rp*a*s12*t3*z3 
       heslo(2,51) = 9.0_rp*a*s1*t32*z3 
       heslo(3,51) = 9.0_rp*a*s1*t3*z32 
       heslo(4,51) = 9.0_rp*a*s11*t31*z3
       heslo(5,51) = 9.0_rp*a*s11*t3*z31
       heslo(6,51) = 9.0_rp*a*s1*t31*z31

       heslo(1,52) = 9.0_rp*a*s12*t2*z3 
       heslo(2,52) = 9.0_rp*a*s1*t22*z3 
       heslo(3,52) = 9.0_rp*a*s1*t2*z32 
       heslo(4,52) = 9.0_rp*a*s11*t21*z3
       heslo(5,52) = 9.0_rp*a*s11*t2*z31
       heslo(6,52) = 9.0_rp*a*s1*t21*z31

       heslo(1,53) = 9.0_rp*a*s22*t2*z4 
       heslo(2,53) = 9.0_rp*a*s2*t22*z4 
       heslo(3,53) = 9.0_rp*a*s2*t2*z42 
       heslo(4,53) = 9.0_rp*a*s21*t21*z4
       heslo(5,53) = 9.0_rp*a*s21*t2*z41
       heslo(6,53) = 9.0_rp*a*s2*t21*z41

       heslo(1,54) = 9.0_rp*a*s32*t2*z4 
       heslo(2,54) = 9.0_rp*a*s3*t22*z4 
       heslo(3,54) = 9.0_rp*a*s3*t2*z42 
       heslo(4,54) = 9.0_rp*a*s31*t21*z4
       heslo(5,54) = 9.0_rp*a*s31*t2*z41
       heslo(6,54) = 9.0_rp*a*s3*t21*z41

       heslo(1,55) = 9.0_rp*a*s32*t3*z4 
       heslo(2,55) = 9.0_rp*a*s3*t32*z4 
       heslo(3,55) = 9.0_rp*a*s3*t3*z42 
       heslo(4,55) = 9.0_rp*a*s31*t31*z4
       heslo(5,55) = 9.0_rp*a*s31*t3*z41
       heslo(6,55) = 9.0_rp*a*s3*t31*z41

       heslo(1,56) = 9.0_rp*a*s22*t3*z4 
       heslo(2,56) = 9.0_rp *a*s2*t32*z4 
       heslo(3,56) = 9.0_rp*a*s2*t3*z42 
       heslo(4,56) = 9.0_rp*a*s21*t31*z4
       heslo(5,56) = 9.0_rp*a*s21*t3*z41
       heslo(6,56) = 9.0_rp*a*s2*t31*z41



       heslo(1,57) = 27.0_rp*a*s22*t2*z2 
       heslo(2,57) = 27.0_rp*a*s2*t22*z2 
       heslo(3,57) = 27.0_rp*a*s2*t2*z22 
       heslo(4,57) = 27.0_rp*a*s21*t21*z2
       heslo(5,57) = 27.0_rp*a*s21*t2*z21
       heslo(6,57) = 27.0_rp*a*s2*t21*z21

       heslo(1,58) = 27.0_rp*a*s32*t2*z2 
       heslo(2,58) = 27.0_rp*a*s3*t22*z2 
       heslo(3,58) = 27.0_rp*a*s3*t2*z22 
       heslo(4,58) = 27.0_rp*a*s31*t21*z2
       heslo(5,58) = 27.0_rp*a*s31*t2*z21
       heslo(6,58) = 27.0_rp*a*s3*t21*z21

       heslo(1,59) = 27.0_rp*a*s32*t3*z2 
       heslo(2,59) = 27.0_rp*a*s3*t32*z2 
       heslo(3,59) = 27.0_rp*a*s3*t3*z22 
       heslo(4,59) = 27.0_rp*a*s31*t31*z2
       heslo(5,59) = 27.0_rp*a*s31*t3*z21
       heslo(6,59) = 27.0_rp*a*s3*t31*z21

       heslo(1,60) = 27.0_rp*a*s22*t3*z2 
       heslo(2,60) = 27.0_rp*a*s2*t32*z2 
       heslo(3,60) = 27.0_rp*a*s2*t3*z22 
       heslo(4,60) = 27.0_rp*a*s21*t31*z2
       heslo(5,60) = 27.0_rp*a*s21*t3*z21
       heslo(6,60) = 27.0_rp*a*s2*t31*z21

       heslo(1,61) = 27.0_rp*a*s22*t2*z3 
       heslo(2,61) = 27.0_rp*a*s2*t22*z3 
       heslo(3,61) = 27.0_rp*a*s2*t2*z32 
       heslo(4,61) = 27.0_rp*a*s21*t21*z3
       heslo(5,61) = 27.0_rp*a*s21*t2*z31
       heslo(6,61) = 27.0_rp*a*s2*t21*z31

       heslo(1,62) = 27.0_rp*a*s32*t2*z3 
       heslo(2,62) = 27.0_rp*a*s3*t22*z3 
       heslo(3,62) = 27.0_rp*a*s3*t2*z32 
       heslo(4,62) = 27.0_rp*a*s31*t21*z3
       heslo(5,62) = 27.0_rp*a*s31*t2*z31
       heslo(6,62) = 27.0_rp*a*s3*t21*z31

       heslo(1,63) = 27.0_rp*a*s32*t3*z3 
       heslo(2,63) = 27.0_rp*a*s3*t32*z3 
       heslo(3,63) = 27.0_rp*a*s3*t3*z32 
       heslo(4,63) = 27.0_rp*a*s31*t31*z3
       heslo(5,63) = 27.0_rp*a*s31*t3*z31
       heslo(6,63) = 27.0_rp*a*s3*t31*z31

       heslo(1,64) = 27.0_rp*a*s22*t3*z3 
       heslo(2,64) = 27.0_rp*a*s2*t32*z3 
       heslo(3,64) = 27.0_rp*a*s2*t3*z32 
       heslo(4,64) = 27.0_rp*a*s21*t31*z3
       heslo(5,64) = 27.0_rp*a*s21*t3*z31
       heslo(6,64) = 27.0_rp*a*s2*t31*z31

       ! Linear Prism
    else if(nnode==6) then
       shape(   1) = (1.0_rp-s-t)*(1.0_rp-z)
       deriv(1, 1) = z-1.0_rp
       deriv(2, 1) = z-1.0_rp
       deriv(3, 1) = s+t-1.0_rp
       shape(   2) = s*(1.0_rp-z)
       deriv(1, 2) = 1.0_rp-z
       deriv(3, 2) = -s
       shape(   3) = t*(1.0_rp-z)
       deriv(2, 3) = 1.0_rp-z
       deriv(3, 3) = -t
       shape(   4) = (1.0_rp-s-t)*z
       deriv(1, 4) = -z
       deriv(2, 4) = -z
       deriv(3, 4) = 1.0_rp-s-t
       shape(   5) = s*z
       deriv(1, 5) = z
       deriv(3, 5) = s
       shape(   6) = t*z
       deriv(2, 6) = z
       deriv(3, 6) = t

    end if

  end subroutine shape3
  
end module interpolation_names
