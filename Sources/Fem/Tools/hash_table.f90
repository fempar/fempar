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
module hash_table_class

  ! A simple hash table for INTEGER key val pairs. 
  ! - if a minimum and a maximum values are given tbl_size is computed
  !   as (tbl_max-tbl_min)/dll_size in order to have (a mean of) sll_size elements 
  !   in each linked list.
  ! - if not given min=0 is assumed.
  ! - a default size of 1000 is assumed for the table (corresponding to
  !   tbl_max=10^5)

!$ use omp_lib
  use types
  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: tbl_size = 1000
  integer(ip), parameter :: nod_size = 20

  ! Public constants
  integer(ip), parameter :: was_stored     = 1
  integer(ip), parameter :: now_stored     = 2
  integer(ip), parameter :: bad_keyval     = 3 ! A given key already stored with a different value
  integer(ip), parameter :: key_found      = 4
  integer(ip), parameter :: key_not_found  = 5
  integer(ip), parameter :: deleted        = 6
  integer(ip), parameter :: error          = 10

  ! Public status
  character(13), parameter :: stat(7)=(/'was_stored   ',&  
                                        'now_stored   ',&
                                        'bad_keyval   ',&
                                        'key_found    ',&
                                        'key_not_found',&
                                        'deleted      ',&
                                        'error        '/)

#define hash_table hash_table_ip_ip
#define hash_node  hash_node_ip_ip
#define key_type   integer(ip)
#define val_type   integer(ip)
#define put_hash_node     put_hash_node_ip_ip   
#define get_hash_node     get_hash_node_ip_ip   
#define del_hash_node     del_hash_node_ip_ip   
#define free_hash_node    free_hash_node_ip_ip  
#define print_hash_node   print_hash_node_ip_ip 
#define init_hash_table   init_hash_table_ip_ip 
#define put_hash_table    put_hash_table_ip_ip  
#define get_hash_table    get_hash_table_ip_ip  
#define del_hash_table    del_hash_table_ip_ip  
#define free_hash_table   free_hash_table_ip_ip 
#define print_hash_table  print_hash_table_ip_ip
#define status_hash_table  status_hash_table_ip_ip
#include "hash_table_header.i90"

#define hash_table hash_table_igp_ip
#define hash_node  hash_node_igp_ip
#define key_type   integer(igp)
#define val_type   integer(ip)
#define put_hash_node     put_hash_node_igp_ip   
#define get_hash_node     get_hash_node_igp_ip   
#define del_hash_node     del_hash_node_igp_ip   
#define free_hash_node    free_hash_node_igp_ip  
#define print_hash_node   print_hash_node_igp_ip 
#define init_hash_table   init_hash_table_igp_ip 
#define put_hash_table    put_hash_table_igp_ip  
#define get_hash_table    get_hash_table_igp_ip  
#define del_hash_table    del_hash_table_igp_ip  
#define free_hash_table   free_hash_table_igp_ip 
#define print_hash_table  print_hash_table_igp_ip
#define status_hash_table  status_hash_table_igp_ip
#include "hash_table_header.i90"

  public :: hash_table_ip_ip, hash_table_igp_ip

  public :: was_stored, now_stored, bad_keyval, key_found, key_not_found, deleted, error

  public :: tbl_size, nod_size, stat

contains

#define hash_table hash_table_ip_ip
#define hash_node  hash_node_ip_ip
#define key_type   integer(ip)
#define key_size   ip 
#define val_type   integer(ip)
#define print_key_val     write(*,'(i10,2x,i10)') list%key, list%val
#define put_hash_node     put_hash_node_ip_ip   
#define get_hash_node     get_hash_node_ip_ip   
#define del_hash_node     del_hash_node_ip_ip   
#define free_hash_node    free_hash_node_ip_ip  
#define print_hash_node   print_hash_node_ip_ip 
#define init_hash_table   init_hash_table_ip_ip 
#define put_hash_table    put_hash_table_ip_ip  
#define del_hash_table    del_hash_table_ip_ip  
#define get_hash_table    get_hash_table_ip_ip  
#define free_hash_table   free_hash_table_ip_ip 
#define print_hash_table  print_hash_table_ip_ip
#define status_hash_table  status_hash_table_ip_ip
#include "hash_table_body.i90"

#define hash_table hash_table_igp_ip
#define hash_node  hash_node_igp_ip
#define key_type   integer(igp)
#define key_size   igp
#define val_type   integer(ip)
#define print_key_val     write(*,'(i10,2x,i10)') list%key, list%val
#define put_hash_node     put_hash_node_igp_ip   
#define get_hash_node     get_hash_node_igp_ip   
#define del_hash_node     del_hash_node_igp_ip   
#define free_hash_node    free_hash_node_igp_ip  
#define print_hash_node   print_hash_node_igp_ip 
#define init_hash_table   init_hash_table_igp_ip 
#define put_hash_table    put_hash_table_igp_ip  
#define del_hash_table    del_hash_table_igp_ip  
#define get_hash_table    get_hash_table_igp_ip  
#define free_hash_table   free_hash_table_igp_ip 
#define print_hash_table  print_hash_table_igp_ip
#define status_hash_table  status_hash_table_igp_ip
#include "hash_table_body.i90"

end module hash_table_class
