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
program test_hash_table
use types_names
  use hash_table_names
  implicit none
  type(hash_table_igp_ip_t) :: table
  integer(ip), parameter :: tbl_length = 100
  integer(ip), parameter :: tbl_min    = 18
  integer(ip), parameter :: tbl_max    = 1300
  integer(ip)            :: out,istat
  integer(igp)           :: input, igp_tmp
!
!  The initial status of the table is, e.g. 
!    18     (-1) => NULL
!    25     (-1) => NULL
!    80     (-1) => NULL
!
!  After elements are added:
!    18     (18) <=> (218) <=> (318) <=> (418) <=> (-1) => NULL
!    25     (25) <=> (125) <=> (225) <=> (325) <=> (-1) => NULL
!    80     (180) <=> (380) <=> (480) <=> (-1) => NULL
!
  call table%init(tbl_length)

  ! Set global to local values
  !input=18
  !call table%put(input , val=1 , stat=istat) ; if(istat/=now_stored) call table%status(6,istat)
      input = 2320000000_igp
  out =  1; call table%put(key=input+18 , val=out , stat=istat) ; if(istat/=now_stored) call table%status(6,istat)
  out =  2; call table%put(key=input+380, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out =  3; call table%put(key=input+25 , val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out =  4; call table%put(key=input+318, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out =  5; call table%put(key=input+418, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  ! Here is a non increasing one
  out =  6; call table%put(key=input+218, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  !call table%print
  out =  7; call table%put(key=input+180, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out =  8; call table%put(key=input+480, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out =  9; call table%put(key=input+125, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out = 10; call table%put(key=input+325, val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  igp_tmp = 225
  out = 11; call table%put(key=igp_tmp, val=out, stat=istat); if(istat/=now_stored) call table%status(6,istat)

  write(*,*) 'Two repeated keys :'
  out = 12; call table%put(key=input+18 , val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)
  out = 12; call table%put(key=input+25 , val=out , stat=istat); if(istat/=now_stored) call table%status(6,istat)

  ! Delete a key
  call table%print
  call table%del(key=input+318 , stat=istat); if(istat/=deleted) call table%status(6,istat)
  write(*,*) 'Deleted 318'
  call table%del(key=input+180 , stat=istat); if(istat/=deleted) call table%status(6,istat)
  write(*,*) 'Deleted 180'
  call table%del(key=input+380 , stat=istat); if(istat/=deleted) call table%status(6,istat)
  write(*,*) 'Deleted 380'
  call table%del(key=input+480 , stat=istat); if(istat/=deleted) call table%status(6,istat)
  write(*,*) 'Deleted 480'
  call table%del(key=input+225 , stat=istat); if(istat/=deleted) call table%status(6,istat)
  write(*,*) 'Deleted 225'

  ! Get some values
  call table%get(key=input+181,val=out, stat=istat); if(istat/=key_found) call table%status(6,istat)
  write(*,*) 'Get result 1: ', out
  call table%get(key=input+25 ,val=out, stat=istat); if(istat/=key_found) call table%status(6,istat)
  write(*,*) 'Get result 2: ', out
  call table%get(key=input+218,val=out, stat=istat); if(istat/=key_found) call table%status(6,istat)
  write(*,*) 'Get result 3: ', out


  ! Print table
  print*, 'indices of the hash table with content:'
  call table%print

  ! Free table
  call table%free

end program test_hash_table

