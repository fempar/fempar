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
module sort_class
  use types
  implicit none
# include "debug.i90"
  private

  interface icomp
     module procedure icomp_ip, icomp_igp
  end interface icomp

  interface sort
     module procedure sort_ip, sort_igp
  end interface sort

  interface sort_array_cols_by_row_section
     module procedure sort_array_cols_by_row_section_ip, sort_array_cols_by_row_section_igp
  end interface sort_array_cols_by_row_section

  interface sort_array_cols_by_row_element
     module procedure sort_array_cols_by_row_element_ip, sort_array_cols_by_row_element_igp
  end interface sort_array_cols_by_row_element

  public :: sort_array_cols_by_row_element
  public :: sort_array_cols_by_row_section         ! would replace intsort if everything is ok.
  public :: sort
  public :: intsort
  public :: icomp

!-----------------------------------------------------------------------
! 
! sortix is used in mesh_partition.f90 line 2312
! 
! intsort is used in 
! * matrix.f90 line 431
! * mesh_partitipon.f90 lines 850 y 2240
! * par_migration_mesh_partition.f90 lines 558 and 845
! * par_dof_handler_mesh_partition.f90 lines 517 and 824
!
! icomp is used in
! * mesh_partitipon.f90 line 2262
! * par_migration_mesh_partition.f90 line 599
! * par_dof_handler_mesh_partition.f90 lines 579
!
!-----------------------------------------------------------------------

contains

! Specialization to a sort by the first k elements of an array 
! (old use of sortix in mesh_partition.f90)
# define name sort_array_cols_by_row_element_ip
# define interface (k,ld,n,data,index)
# define data_size_def  integer(ip) :: k,ld
# define data_input_def integer(ip) :: data(ld,n)
# define data_temp_def  integer(ip) :: datap(ld), datat(ld)
# define data_acces(n)  data(:,n)
! Any function could be used here
# define greater(a,b) greater_by_row_element_ip(k,ld,a,b)
# include "sort.i90"
  logical(lg) function greater_by_row_element_ip(k,n,ia,ib)
    implicit none
    integer(ip) :: k,n,ia(n),ib(n)
    greater_by_row_element_ip = ia(k) > ib(k)
  end function greater_by_row_element_ip

# define name sort_array_cols_by_row_element_igp
# define interface (k,ld,n,data,index)
# define data_size_def  integer(ip)  :: k,ld
# define data_input_def integer(igp) :: data(ld,n)
# define data_temp_def  integer(igp) :: datap(ld), datat(ld)
# define data_acces(n)  data(:,n)
! Any function could be used here
# define greater(a,b) greater_by_row_element_igp(k,ld,a,b)
# include "sort.i90"
  logical(lg) function greater_by_row_element_igp(k,n,ia,ib)
    implicit none
    integer(ip)  :: k,n
    integer(igp) :: ia(n), ib(n)
    greater_by_row_element_igp = ia(k) > ib(k)
  end function greater_by_row_element_igp

# define in_place

! Specialization to a simple ip array
# define name sort_ip
# define interface (n,data,index)
# define data_size_def
# define data_input_def integer(ip) :: data(n)
# define data_temp_def  integer(ip) :: datap, datat
# define data_acces(n)  data(n)
# define greater(a,b) (a)>(b)
# include "sort.i90"

! Specialization to a simple igp array
# define name sort_igp
# define interface (n,data,index)
# define data_size_def
# define data_input_def integer(igp) :: data(n)
# define data_temp_def  integer(igp) :: datap, datat
# define data_acces(n)  data(n)
# define greater(a,b) (a)>(b)
# include "sort.i90"

! Specialization to a sort by the first k elements of an array
# define name sort_array_cols_by_row_section_ip
# define interface (k,ld,n,data,index,datap,datat)
# define data_size_def  integer(ip) :: k,ld
# define data_input_def integer(ip) :: data(ld,n)
# define data_temp_def  integer(ip) :: datap(ld), datat(ld)
# define data_acces(n)  data(:,n)
! Any function could be used here
# define greater(a,b) greater_by_section_ip(k,ld,a,b)
# include "sort.i90"
  logical(lg) function greater_by_section_ip(k,n,ia,ib)
    implicit none
    integer(ip) :: k,n,ia(n),ib(n)
    integer(ip) :: i 
    greater_by_section_ip = .false.

    do i=1,k
       if (ia(i)==ib(i)) then
          cycle
       else
         greater_by_section_ip = ia(i) > ib(i)
         exit
       end if
    end do
  end function greater_by_section_ip

! Specialization to a sort by the first k elements of an array
# define name sort_array_cols_by_row_section_igp
# define interface (k,ld,n,data,index,datap,datat)
# define data_size_def  integer(ip)  :: k,ld
# define data_input_def integer(igp) :: data(ld,n)
# define data_temp_def  integer(igp) :: datap(ld), datat(ld)
# define data_acces(n)  data(:,n)
! Any function could be used here
# define greater(a,b) greater_by_section_igp(k,ld,a,b)
# include "sort.i90"
  logical(lg) function greater_by_section_igp(k,n,ia,ib)
    implicit none
    integer(ip) :: k, n
    integer(igp):: ia(n), ib(n)
    integer(ip) :: i 
    greater_by_section_igp = .false.
    
    do i=1,k
       if (ia(i)==ib(i)) then
          cycle
       else
         greater_by_section_igp = ia(i) > ib(i)
         exit
       end if
    end do
  end function greater_by_section_igp

!-----------------------------------------------------------------------
! The next routine is taken from internet (where?). It is a standard
! quicksort algorithm modified introducing a function icomp which 
! implements a comparison of the whole array IX(:,I) and IX(:,J)
!
! Here:
! N1, intent(in) : number of elementes of IX(:,I) to be compared (see icomp below) 
! N0, intent(in) : leading dimension of IX
! N , intent(in) : last dimension of IX
! IX, intent(in) : data to sort
! IY, intent(out): permutation in the form IY(new)=old
! T , intent(in) : working space of size N0
! TT, intent(in) : working space of size N0
!
!-----------------------------------------------------------------------
         SUBROUTINE INTSORT (N1,N0,N,IX,IY,T,TT)
           implicit none
           INTEGER N0,N1,N
           INTEGER IX(N0,N),IY(N),T(N0),TT(N0)
           REAL R
           INTEGER I, IJ, J, K, KK, L, M, NN, TY, TTY
           INTEGER IL(21), IU(21)
           !INTEGER icomp
           !external icomp
           
           IF (N .EQ. 0) RETURN

           M = 1
           I = 1
           J = N
           R = 0.375E0

20         IF (I .EQ. J) GO TO 60
           IF (R .LE. 0.5898437E0) THEN
              R = R+3.90625E-2
           ELSE
              R = R-0.21875E0
           ENDIF

30         K = I
           !
           !     Select a central element of the array and save it in location T
           !
           IJ = I + INT((J-I)*R)
           T = IX(:,IJ)
           TY = IY(IJ)
           !
           !     If first element of array is greater than T, interchange with T
           !
           !      IF (IX(I1,I) .GT. T(I1)) THEN
           IF (icomp(N1,IX(:,I),T).gt.0) THEN
              IX(:,IJ) = IX(:,I)
              IX(:,I) = T
              T = IX(:,IJ)
              IY(IJ) = IY(I)
              IY(I) = TY
              TY = IY(IJ)
           ENDIF
           L = J
           !
           !     If last element of array is less than than T, interchange with T
           !
           !      IF (IX(I1,J) .LT. T(I1)) THEN
           IF (icomp(N1,IX(:,J),T).lt.0) THEN
              IX(:,IJ) = IX(:,J)
              IX(:,J) = T
              T = IX(:,IJ)
              IY(IJ) = IY(J)
              IY(J) = TY
              TY = IY(IJ)
              !
              !        If first element of array is greater than T, interchange with T
              !
              !         IF (IX(I1,I) .GT. T(I1)) THEN
              IF (icomp(N1,IX(:,I),T).gt.0) THEN
                 IX(:,IJ) = IX(:,I)
                 IX(:,I) = T
                 T = IX(:,IJ)
                 IY(IJ) = IY(I)
                 IY(I) = TY
                 TY = IY(IJ)
              ENDIF
           ENDIF
           !
           !     Find an element in the second half of the array which is smaller
           !     than T
           !
40         L = L-1
           !      IF (IX(I1,L) .GT. T(I1)) GO TO 40
           IF (icomp(N1,IX(:,L),T).gt.0) GO TO 40
           !
           !     Find an element in the first half of the array which is greater
           !     than T
           !
50         K = K+1
           !      IF (IX(I1,K) .LT. T(I1)) GO TO 50
           IF (icomp(N1,IX(:,K),T).lt.0) GO TO 50
           !
           !     Interchange these elements
           !
           IF (K .LE. L) THEN
              TT = IX(:,L)
              IX(:,L) = IX(:,K)
              IX(:,K) = TT
              TTY = IY(L)
              IY(L) = IY(K)
              IY(K) = TTY
              GO TO 40
           ENDIF
           !
           !     Save upper and lower subscripts of the array yet to be sorted
           !
           IF (L-I .GT. J-K) THEN
              IL(M) = I
              IU(M) = L
              I = K
              M = M+1
           ELSE
              IL(M) = K
              IU(M) = J
              J = L
              M = M+1
           ENDIF
           GO TO 70
           !
           !     Begin again on another portion of the unsorted array
           !
60         M = M-1
           IF (M .EQ. 0) GO TO 190
           I = IL(M)
           J = IU(M)

70         IF (J-I .GE. 1) GO TO 30
           IF (I .EQ. 1) GO TO 20
           I = I-1

80         I = I+1
           IF (I .EQ. J) GO TO 60
           T = IX(:,I+1)
           TY = IY(I+1)
           !      IF (IX(I1,I) .LE. T(I1)) GO TO 80
           IF (icomp(N1,IX(:,I),T).le.0) GO TO 80
           K = I

90         IX(:,K+1) = IX(:,K)
           IY(K+1) = IY(K)
           K = K-1
           !      IF (T(I1) .LT. IX(I1,K)) GO TO 90
           IF (icomp(N1,T,IX(:,K)).lt.0) GO TO 90
           IX(:,K+1) = T
           IY(K+1) = TY
           GO TO 80
190        RETURN
         END SUBROUTINE INTSORT
         !-----------------------------------------------------------------------
         function icomp_ip(n,ia,ib)
           implicit none
           integer(ip) :: n,ia(n),ib(n)
           integer(ip) :: icomp_ip
           integer(ip) :: i

           do i=1,n
              if(ia(i).gt.ib(i)) then
                 icomp_ip=1
              else if(ia(i).lt.ib(i)) then
                 icomp_ip=-1
              else
                 icomp_ip=0
              end if
              if(icomp_ip/=0) exit
           end do

         end function icomp_ip

         !-----------------------------------------------------------------------
         function icomp_igp(n,ia,ib)
           implicit none
           integer(ip)  :: n
           integer(igp) :: ia(n),ib(n)
           integer(ip)  :: icomp_igp
           integer(ip)  :: i

           do i=1,n
              if(ia(i).gt.ib(i)) then
                 icomp_igp=1
              else if(ia(i).lt.ib(i)) then
                 icomp_igp=-1
              else
                 icomp_igp=0
              end if
              if(icomp_igp/=0) exit
           end do

         end function icomp_igp



end module sort_class
