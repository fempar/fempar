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
#include "debug.i90"
program test_std_vector
  use fempar_names
  implicit none
  type(std_vector_integer_ip_t) :: v
  type(std_vector_integer_ip_t) :: v2
  integer(ip) :: i
  
  call fempar_init()
 
  call v%resize(0,1)
  call v%resize(1,1)
  
  !Set some initial content:
  do i=1, 100
    call v%push_back(i+1)
    write(*,'(a,i4,a,i4)') 'Size=', v%size(), ' Capacity=', v%capacity() 
    check ( v%get(i) == i )
  end do
  
  call v%resize(1000,34)
  do i=1, 101
    check ( v%get(i) == i )
  end do
  do i=102,1000
    check ( v%get(i) == 34 )
  end do
  
  call v%resize(38,1)
  do i=1, v%size()
    check ( v%get(i) == i )
  end do
  
  call v2%copy(v)
  check (v2%size() == v%size())
  do i=1, v2%size()
    check ( v2%get(i) == v%get(i) )
  end do
  
  call v%erase(10)
  do i=1, 9
    check ( v%get(i) == i )
  end do
  
  do i=10, v%size()
    check ( v%get(i) == i+1)
  end do
  
  call v%shrink_to_fit()
  check ( v%size() == v%capacity() )
  
  do i=1, 9
    check ( v%get(i) == i )
  end do
  
  do i=10, v%size()
    check ( v%get(i) == i+1)
  end do
  
  
   
  call v%free()
  call v2%free()
  call fempar_finalize()
end program test_std_vector
