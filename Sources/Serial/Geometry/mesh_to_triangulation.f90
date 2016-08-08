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
module mesh_to_triangulation_names
  use types_names
  use memor_names
  use mesh_names
  use triangulation_names
  !use generate_vefs_mesh_conditions_names
  use conditions_names
  use reference_fe_names
  use reference_fe_factory_names
  use elem_to_subset_id_names

  implicit none
# include "debug.i90"
  private

  public :: mesh_to_triangulation, mesh_to_triangulation_fill_elements

contains

  !*********************************************************************************
  ! This subroutine takes as input a 'plain' mesh and creates a triangulation,
  ! filling both element and vef (dual mesh) info. The triangulation includes
  ! vertices, edges, and vertices (in 3D). The length_trian optional arguments 
  ! is used when used from par_mesh_triangulation.f90, because in a parallel 
  ! environment ghost elements are also included in the list of vefs of the
  ! triangulation, even though the triangulation is not aware of that, i.e., the
  ! number of elements in triangulation does not include ghost elements, only
  ! local elements
  !*********************************************************************************
  subroutine mesh_to_triangulation (gmesh,trian,gcond,elem_to_subset_id)
    implicit none
    ! Parameters
    type(mesh_t), intent(in)                           :: gmesh ! Geometry mesh
    type(triangulation_t), intent(inout)               :: trian 
    type(conditions_t), optional, intent(inout)        :: gcond
    type(elem_to_subset_id_t), optional, intent(inout) :: elem_to_subset_id

    call mesh_to_triangulation_fill_elements( gmesh, trian, gcond = gcond, &
                                              elem_to_subset_id = elem_to_subset_id)
    call triangulation_to_dual ( trian )

  end subroutine mesh_to_triangulation

subroutine mesh_to_triangulation_fill_elements (gmesh, trian, length_trian, gcond, &
     &                                           elem_to_subset_id)
    implicit none
    ! Parameters
    type(mesh_t)                       , intent(in)    :: gmesh ! Geometry mesh
    type(triangulation_t)              , intent(inout) :: trian 
    integer(ip)              , optional, intent(in)    :: length_trian
    type(conditions_t)       , optional, intent(inout) :: gcond
    type(elem_to_subset_id_t), optional, intent(inout) :: elem_to_subset_id

    ! Locals
    type(mesh_t)            :: tmesh ! Topological mesh
    integer(ip)             :: istat, ielem, iobj
    integer(ip)             :: count, g_node, inode, p, length_trian_
    type(conditions_t)      :: tcond

    assert(trian%state == triangulation_not_created .or. trian%state == triangulation_filled)

    if ( trian%state == triangulation_filled ) call triangulation_free(trian)

    if (present(length_trian)) then
       length_trian_ = length_trian
    else
       length_trian_ = gmesh%nelem 
    endif

    call triangulation_create ( length_trian_, trian )
    
    trian%num_elems = gmesh%nelem
    trian%num_dims  = gmesh%ndime

    !! Variable values depending of the element ndime
    ! nvef = trian%elems(ielem)%num_vefs 
    ! ndime = trian%num_dims
    ! etype = 0
    ! if(ndime == 2) then       ! 2D
    !   if(nvef == 6) then     ! Linear triangles (P1)
    !      etype = P_type_id
    !   elseif(nvef == 8) then ! Linear quads (Q1)
    !      etype = Q_type_id
    !   end if
    ! elseif(ndime == 3) then    ! 3D
    !   if(nvef == 14) then     ! Linear tetrahedra (P1)
    !      etype = P_type_id
    !   elseif(nvef == 26) then ! Linear hexahedra (Q1)
    !      etype = Q_type_id
    !   end if
    ! end if
    ! assert( etype /= 0 )
    
    ! This will not be here in the future
    trian%reference_fe_geo_list(1) = make_reference_fe ( topology = topology_hex, fe_type = fe_type_lagrangian, &
         &                                               number_dimensions = trian%num_dims, order = 1,         &
         &                                               field_type = field_type_vector, continuity = .true. )

    if ( trian%num_elems > trian%elem_array_len ) then
       if (allocated(trian%elems)) then
          deallocate(trian%elems, stat=istat)
          check(istat==0)
       end if
       allocate(trian%elems(trian%num_elems), stat=istat)
       check(istat==0)
       trian%elem_array_len = trian%num_elems
    end if

    if (present(gcond)) then
       !call generate_vefs_mesh_conditions(gmesh, tmesh, trian%reference_fe_geo_list(1)%p, gcond, tcond)
       call conditions_free( gcond )
       call conditions_copy( tcond, gcond )
       call conditions_free( tcond )
    else
       !call generate_vefs_mesh_conditions(gmesh, tmesh, trian%reference_fe_geo_list(1)%p)
    end if
       
    do ielem=1, trian%num_elems
       trian%elems(ielem)%num_vefs = tmesh%pnods(ielem+1)-tmesh%pnods(ielem)
       call memalloc(trian%elems(ielem)%num_vefs, trian%elems(ielem)%vefs, __FILE__, __LINE__)
       trian%elems(ielem)%vefs(1:trian%elems(ielem)%num_vefs) = tmesh%lnods(tmesh%pnods(ielem):tmesh%pnods(ielem+1)-1)
       call put_topology_element_triangulation ( ielem, trian )
    end do

    call tmesh%free()

    assert ( allocated(gmesh%coord) )
    do ielem = 1, trian%num_elems
       call memalloc( trian%num_dims, gmesh%pnods(ielem+1)-gmesh%pnods(ielem), &
                   &  trian%elems(ielem)%coordinates, __FILE__, __LINE__ )
       count = 0
       do inode = gmesh%pnods(ielem),gmesh%pnods(ielem+1)-1
          count = count+1
          g_node = gmesh%lnods(inode)
          trian%elems(ielem)%coordinates(1:trian%num_dims, count) = gmesh%coord(1:trian%num_dims, g_node)
       end do
    end do

    if (present(elem_to_subset_id)) then
       do ielem = 1, trian%num_elems
          trian%elems(ielem)%subset_id = elem_to_subset_id%get_subset_id(ielem)
       end do
    end if

  end subroutine mesh_to_triangulation_fill_elements


end module mesh_to_triangulation_names



















