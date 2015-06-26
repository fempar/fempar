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
module fem_mesh_partition_base_names
use types_names

  ! Basic parameters and type definitions for serial partitioning routines
  ! All PUBLIC 

  integer(ip), parameter :: part_kway      = 0
  integer(ip), parameter :: part_recursive = 1
  integer(ip), parameter :: part_strip     = 2
  integer(ip), parameter :: part_rcm_strip = 3
  
  integer(ip), parameter :: dual=0
  integer(ip), parameter :: primal=1

  type part_params_t
     integer(ip) :: nparts      = 2             ! nparts
     integer(ip) :: debug       = 1             ! Print info partition

     integer(ip) :: strat = part_kway  ! Partitioning algorithm 
                                       ! (for vertex_based+primal or element_based+dual)
                                       ! part_kway, part_recursive, part_strip, part_rcm_strip

     ! Only applicable to metis 5.0 for both part_kway and part_recursive
     ! Use METIS defaults (i.e., == -1) 30 for part_kway, and 1 for part_recursive
     integer(ip) :: metis_option_ufactor = -1 ! Imbalance tol of x/1000 + 1

     ! Only applicable to metis 5.0 and part_kway
     integer(ip) :: metis_option_minconn = 1 ! (Try to) Minimize maximum degree of subdomain graph
     integer(ip) :: metis_option_contig  = 1 ! (Try to) Produce partitions that are contiguous
     
     ! Applicable to both metis 4.0 and metis 5.0
     integer(ip) :: metis_option_debug  =  0 
  end type part_params_t


end module fem_mesh_partition_base_names
