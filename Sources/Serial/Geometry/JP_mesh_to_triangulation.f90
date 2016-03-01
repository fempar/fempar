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
module JP_mesh_to_triangulation_names
  use types_names
  use memor_names
  use mesh_names
  use JP_triangulation_names
  use serial_triangulation_names
  use JP_element_topology_names
  !use fe_space_types_names
  use generate_vefs_mesh_conditions_names
  use conditions_names

  implicit none
# include "debug.i90"
  private

  public :: JP_mesh_to_triangulation, serial_mesh_to_triangulation

contains

  !*********************************************************************************
  ! This subroutine takes as input a 'plain' mesh and creates a triangulation,
  ! filling both element and vef (dual mesh) info. The triangulation includes
  ! vertices, edges, and vertices (in 3D). 
  !*********************************************************************************
  subroutine serial_mesh_to_triangulation (gmesh,trian)
    implicit none
    type(mesh_t), intent(in)                    :: gmesh ! Geometry mesh
    type(serial_triangulation_t), intent(inout) :: trian 
    call trian%create(gmesh%nelem)
    call JP_mesh_to_triangulation(gmesh,trian)
    call create_reference_elements(trian)
    call trian%to_dual()
  end subroutine serial_mesh_to_triangulation

  !*********************************************************************************
  ! This subroutine takes as input a 'plain' mesh and creates a triangulation,
  ! filling elements only and assuming the triangulation already allocated. 
  !*********************************************************************************
  subroutine JP_mesh_to_triangulation (gmesh, trian)
    implicit none
    ! Parameters
    type(mesh_t)             , intent(in)    :: gmesh ! Geometry mesh
    class(JP_triangulation_t), intent(inout) :: trian 

    ! Locals
    integer(ip)              :: ielem, count, g_node, inode, num_vefs
    class(JP_element_topology_t), pointer :: elem
    real(rp)   , allocatable :: coordinates(:,:)

    assert ( allocated(gmesh%coord) )
    call memalloc(gmesh%ndime, gmesh%nnode, coordinates, __FILE__, __LINE__)

    trian%num_elems = gmesh%nelem
    trian%num_dims  = gmesh%ndime
    trian%num_vefs  = gmesh%npoin

    ! trian and mesh are traversed at the same time, the first one using
    ! the iterator and the second one directly.
    call trian%element_iterator%begin()
    ielem=1
    do while( (.not.trian%element_iterator%finished()) .and. (ielem<=gmesh%nelem) )
       elem => downcast_to_element_topology( trian%element_iterator%current() )
       num_vefs = gmesh%pnods(ielem+1)-gmesh%pnods(ielem)
       ! Gather coordinates
       count = 0
       do inode = gmesh%pnods(ielem),gmesh%pnods(ielem+1)-1
          count = count+1
          g_node = gmesh%lnods(inode)
          coordinates(1:trian%num_dims, count) = gmesh%coord(1:trian%num_dims, g_node)
       end do
       ! Fill element
       call elem%fill(trian%num_dims, num_vefs, trian%num_vefs, &
            &         gmesh%lnods(gmesh%pnods(ielem):gmesh%pnods(ielem+1)-1), coordinates)
       ! Keep going
       call trian%element_iterator%next()
       ielem = ielem+1
    end do

  end subroutine JP_mesh_to_triangulation

end module JP_mesh_to_triangulation_names



















