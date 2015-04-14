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
module fem_space_faces
  use types
  use memor
  use hash_table_names
  use fem_mesh_names
  use integration_names
  !  use fem_mesh_faces
  use fem_space_names
  use fem_space_types
  use fem_triangulation_names
  use array_names

# include "debug.i90"
  implicit none

  private
  integer(ip), parameter :: max_subfaces = 4
  ! Functions
  public :: fem_space_faces_list_create

contains

  !===================================================================================================
  subroutine fem_space_faces_list_create( trian, femsp )
    implicit none
    ! Parameters
    type(fem_space)      , intent(inout),target  :: femsp
    type(fem_triangulation)   , intent(in), target   :: trian   

    ! Local variables
    integer(ip) :: iface, iobje, max_order, ielem, nvars, ivars, ndofs, pos_elmat, pos_elvec
    integer(ip) :: i, l_faci, pos_faint, iprob, istat
    integer(ip) :: gtype, utype, g_ord, u_ord, v_key, iface_l, max_elmat, ndime
    type(fem_fixed_info_pointer) :: gfinf, ufinf

    !allocate(femsp%interior_faces( femsp%num_interior_faces ), stat=istat)
    !allocate(femsp%boundary_faces( femsp%num_boundary_faces ), stat=istat)

    ndime = femsp%g_trian%num_dims

    if (femsp%num_interior_faces /= 0) then

       ! Fill the list of faces with the neighbouring elements
       do iface_l = 1,femsp%num_interior_faces
          iobje = femsp%interior_faces(iface_l)%face_object
          max_order = 0
          ndofs = 0
          do i = 1,2
             ielem = trian%objects(iobje)%elems_around(i)
             iprob = femsp%lelem(ielem)%problem
             nvars = femsp%dof_handler%problems(iprob)%nvars
             do ivars = 1, nvars
                max_order = max(max_order,femsp%lelem(ielem)%order(ivars))
             end do
             ndofs = ndofs + femsp%lelem(ielem)%f_inf(ivars)%p%nnode
          end do

          ! Create elemental matrix and vectors
          call femsp%ht_pos_elmat%put(key=ndofs,val=femsp%cur_elmat,stat=istat)!
          call femsp%ht_pos_elvec%put(key=ndofs,val=femsp%cur_elvec,stat=istat)!
          if ( istat == now_stored ) then
             call array_create ( ndofs, ndofs, femsp%lelmat(femsp%cur_elmat) )
             pos_elmat = femsp%cur_elmat!
             femsp%cur_elmat = femsp%cur_elmat + 1!
             call array_create ( ndofs, femsp%lelvec(femsp%cur_elmat) )
             pos_elvec = femsp%cur_elvec!
             femsp%cur_elvec = femsp%cur_elvec + 1!
          else if ( istat == was_stored ) then
             call femsp%ht_pos_elmat%get(key=ndofs,val=pos_elmat,stat=istat)!
             assert ( istat == key_found )
             call femsp%ht_pos_elvec%get(key=ndofs,val=pos_elvec,stat=istat)!
             assert ( istat == key_found )
          end if
          femsp%lface(iface)%p_mat        => femsp%lelmat(pos_elmat)!
          femsp%lface(iface)%p_vec        => femsp%lelvec(pos_elvec)!

          do i = 1,2
             ielem = trian%objects(iobje)%elems_around(i)
             iprob = femsp%lelem(ielem)%problem
             femsp%interior_faces(iface_l)%neighbor_element(i) = ielem 
             l_faci = local_position(iobje, trian%elems(ielem)%objects, &
                  & trian%elems(ielem)%num_objects )
             femsp%interior_faces(iface_l)%local_face(i) = l_faci - &
                  femsp%g_trian%elems(ielem)%topology%nobje_dim(ndime) + 1 
             nvars = femsp%dof_handler%problems(iprob)%nvars

             call memalloc( nvars, femsp%lface(iface)%integ(i)%p, __FILE__, __LINE__ )

             do ivars = 1, nvars
                gtype = femsp%lelem(femsp%lface(iface)%neighbor_element(i))%p_geo_info%ftype
                utype = femsp%lelem(femsp%lface(iface)%neighbor_element(i))%f_inf(ivars)%p%ftype
                !assert( utype == gtype )
                g_ord = femsp%lelem(femsp%lface(iface)%neighbor_element(i))%p_geo_info%order
                u_ord = femsp%lelem(femsp%lface(iface)%neighbor_element(i))%f_inf(ivars)%p%order
                ! SB.alert : The last part to include gauss points being used
                v_key =  utype + (max_FE_types+1)*u_ord + (max_FE_types+1)*(max_order+1)*max_order
                ! Put in hash table
                call femsp%ht_pos_face_integ%put(key=v_key, val=femsp%cur_lfaci, stat = istat)
                if ( istat == now_stored ) then 
                   gfinf%p => femsp%lelem(femsp%lface(iface)%neighbor_element(i))%p_geo_info
                   ufinf%p => femsp%lelem(femsp%lface(iface)%neighbor_element(i))%f_inf(ivars)%p
                   call integ_create(gfinf,ufinf,ndime,femsp%lfaci(femsp%cur_lfaci))
                   pos_faint       = femsp%cur_lfaci
                   femsp%cur_lfaci = femsp%cur_lfaci + 1
                else if ( istat == was_stored ) then
                   call femsp%ht_pos_face_integ%get(key=v_key,val=pos_faint,stat=istat)
                   assert ( istat == key_found )
                end if
                femsp%lface(iface)%integ(i)%p(ivars)%p => femsp%lfaci(pos_faint)
             end do

          end do
       end do

    end if

    ! SB.alert : o2n pending... to be reconsidered. What does the key depend on? Hash table, etc.

  end subroutine fem_space_faces_list_create

  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)

    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( 0 == 1 )

  end function local_position

end module fem_space_faces
