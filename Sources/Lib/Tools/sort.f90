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
module sort_names
use types_names
  implicit none
# include "debug.i90"
  private

  interface sort
     module procedure sort_ip, sort_igp, sort_rp
  end interface sort

  interface sort_array_cols_by_row_section
     module procedure sort_array_cols_by_row_section_ip, sort_array_cols_by_row_section_igp
  end interface sort_array_cols_by_row_section

  interface sort_array_cols_by_row_element
     module procedure sort_array_cols_by_row_element_ip, sort_array_cols_by_row_element_igp
  end interface sort_array_cols_by_row_element

  public :: sort_array_cols_by_row_element
  public :: sort_array_cols_by_row_section
  public :: sort

contains

# define in_place

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
  logical function greater_by_row_element_ip(k,n,ia,ib)
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
  logical function greater_by_row_element_igp(k,n,ia,ib)
    implicit none
    integer(ip)  :: k,n
    integer(igp) :: ia(n), ib(n)
    greater_by_row_element_igp = ia(k) > ib(k)
  end function greater_by_row_element_igp

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

! Specialization to a simple rp array
# define name sort_rp
# define interface (n,data,index)
# define data_size_def
# define data_input_def real(rp) :: data(n)
# define data_temp_def  real(rp) :: datap, datat
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
  logical function greater_by_section_ip(k,n,ia,ib)
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
  logical function greater_by_section_igp(k,n,ia,ib)
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
end module sort_names
