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
module dof_distribution_names
  ! Serial modules
  use types
  use memor
  use maps_names
  use dof_import_names
    
  implicit none
# include "debug.i90"
  private
  
  type dof_distribution
     integer(ip)              :: ipart
     integer(ip)              :: nparts
     
     integer(ip)              :: nl ! Number of local DoFs (formely in f_part%nmap%nl)
     integer(ip)              :: ni ! Number of interior DoFs (formely in f_part%nmap%ni)
     integer(ip)              :: nb ! Number of interface DoFs (formely in f_part%nmap%nb)

     integer(ip)              :: npadj       ! Number of adjacent parts
     integer(ip), allocatable :: lpadj(:)    ! List of adjacent parts
     
     ! To add to the documentation, the contents of lobjs are as follows:
     ! lobjs(1,*) : pysical variable associated to the DoF
     ! lobjs(2,*) : first DoF of the object in local ordering
     ! lobjs(3,*) : last DoF of the object in local ordering
     ! lobjs(4,*) : number of parts around this object
     ! lobjs(5:,*) : list of parts around this object

     integer(ip)              :: max_nparts  ! Maximum number of parts around any DoF communication object
     integer(ip)              :: nobjs       ! Number of local DoF communication objects
     integer(ip), allocatable :: lobjs(:,:)  ! List of local DoF communication objects
     
     type(list)               :: int_objs    ! List of objects on each edge to an adjacent part / Interface_objects
     type(map)                :: omap        ! Objects local to global map

     type(dof_import) :: dof_import ! Object which contains the control data to drive DoF nearest neigbour exchanges
  end type dof_distribution

  ! Types
  public :: dof_distribution
  
  ! Functions
  public :: dof_distribution_free
  
contains

  subroutine dof_distribution_free(dof_dist)
    implicit none
    type(dof_distribution), intent(inout) :: dof_dist
    
    call memfree ( dof_dist%lpadj, __FILE__,__LINE__  )  
    call memfree ( dof_dist%lobjs, __FILE__,__LINE__  )
    call memfree ( dof_dist%int_objs%p, __FILE__, __LINE__ )
    call memfree ( dof_dist%int_objs%l, __FILE__, __LINE__ )
    call map_free ( dof_dist%omap )
    call dof_import_free ( dof_dist%dof_import )

  end subroutine dof_distribution_free
  
end module dof_distribution_names
