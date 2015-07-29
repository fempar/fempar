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
module fe_space_faces_names
  use types_names
  use memor_names
  use hash_table_names
  use mesh_names
  !use integration_names
  use integration_tools_names
  !  use mesh_faces
  use fe_space_names
  use fe_space_types_names
  use triangulation_names
  use array_names

# include "debug.i90"
  implicit none

  private
  integer(ip), parameter :: max_subfaces = 4
  ! Functions
  public :: fe_space_faces_list_create

contains

  !===================================================================================================
  subroutine fe_space_faces_list_create( trian, fe_space )
    implicit none
    ! Parameters
    type(fe_space_t)      , intent(inout),target  :: fe_space
    type(triangulation_t)   , intent(in), target   :: trian   

    ! Local variables
    integer(ip) :: iface, iobje, max_order, ielem, nvars, ivars, ndofs, pos_elmatvec
    integer(ip) :: i, l_faci, pos_faint, iprob, istat
    integer(ip) :: gtype, utype, g_ord, u_ord, v_key, iface_l, max_elmat, ndime
    type(reference_element_pointer_t) :: geo_reference_element, unk_reference_element
    integer(ip) :: aux_val

    !allocate(fe_space%interior_faces( fe_space%num_interior_faces ), stat=istat)
    !allocate(fe_space%boundary_faces( fe_space%num_boundary_faces ), stat=istat)

    ndime = fe_space%g_trian%num_dims

    if (fe_space%num_interior_faces /= 0) then

       ! Fill the list of faces with the neighbouring elements
       do iface_l = 1,fe_space%num_interior_faces
          iobje = fe_space%interior_faces(iface_l)%face_vef
          max_order = 0
          ndofs = 0
          do i = 1,2
             ielem = trian%vefs(iobje)%elems_around(i)
             iprob = fe_space%finite_elements(ielem)%problem
             nvars = fe_space%dof_descriptor%problems(iprob)%p%nvars
             do ivars = 1, nvars
                max_order = max(max_order,fe_space%finite_elements(ielem)%order(ivars))
             end do
             ndofs = ndofs + fe_space%finite_elements(ielem)%reference_element_vars(ivars)%p%nnode
          end do

          ! Create elemental matrix and vectors
          call fe_space%pos_elmatvec%get(key=ndofs,val=pos_elmatvec,stat=istat)
          if ( istat == new_index ) then 
             call array_create ( ndofs, ndofs, fe_space%lelmat(pos_elmatvec) )
             call array_create ( ndofs, fe_space%lelvec(pos_elmatvec) )
          end if
          fe_space%fe_faces(iface)%p_mat => fe_space%lelmat(pos_elmatvec)
          fe_space%fe_faces(iface)%p_vec => fe_space%lelvec(pos_elmatvec)

          do i = 1,2
             ielem = trian%vefs(iobje)%elems_around(i)
             iprob = fe_space%finite_elements(ielem)%problem
             fe_space%interior_faces(iface_l)%neighbor_element(i) = ielem 
             l_faci = local_position(iobje, trian%elems(ielem)%vefs, &
                  & trian%elems(ielem)%num_vefs )
             fe_space%interior_faces(iface_l)%local_face(i) = l_faci - &
                  fe_space%g_trian%elems(ielem)%geo_reference_element%nvef_dim(ndime) + 1 
             nvars = fe_space%dof_descriptor%problems(iprob)%p%nvars

             call memalloc( nvars, fe_space%fe_faces(iface)%integ(i)%p, __FILE__, __LINE__ )

             do ivars = 1, nvars
                gtype = fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%p_geo_reference_element%ftype
                utype = fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%reference_element_vars(ivars)%p%ftype
                !assert( utype == gtype )
                g_ord = fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%p_geo_reference_element%order
                u_ord = fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%reference_element_vars(ivars)%p%order
                ! SB.alert : The last part to include gauss points being used
                v_key =  utype + (max_FE_types+1)*u_ord + (max_FE_types+1)*(max_order+1)*max_order

                call fe_space%pos_face_integrator%get(key=v_key, val=pos_faint, stat = istat)
                if ( istat == new_index ) then 
                   geo_reference_element%p => fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%p_geo_reference_element
                   unk_reference_element%p => fe_space%finite_elements(fe_space%fe_faces(iface)%neighbor_element(i))%reference_element_vars(ivars)%p
                   call face_integrator_create(geo_reference_element,unk_reference_element,ndime,fe_space%lfaci(pos_faint))
                end if
                fe_space%fe_faces(iface)%integ(i)%p(ivars)%p => fe_space%lfaci(pos_faint)

             end do
          end do
       end do

    end if

    ! SB.alert : o2n pending... to be reconsidered. What does the key depend on? Hash table, etc.

  end subroutine fe_space_faces_list_create

  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)

    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( 0 == 1 )

  end function local_position

end module fe_space_faces_names
