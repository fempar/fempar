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
module integration_names
  use types
  use integrable_names
  use problem_names
  use integration_tools_names
  use fem_space_names
  use assembly_names
  implicit none
  private

contains

  subroutine volume_integral(femsp,int1,int2)
    implicit none
    ! Parameters
    type(fem_space)  , intent(inout) :: femsp
    class(integrable), intent(inout) :: int1
    class(integrable), optional, intent(inout) :: int2

    ! Locals
    integer(ip) :: ielem,ivar 
    type(physical_problem), pointer :: problem

    ! Main element loop
    do ielem=1,femsp%g_trian%num_elems

       problem => femsp%dof_handler%problems(femsp%lelem(ielem)%problem)

       ! Compute integration tools on ielem
       do ivar=1,femsp%lelem(ielem)%num_vars
          !call volume_integrator_fill_val(ielem,femsp%g_mesh,femsp%lelem(ielem)%integ(ivar)%p)
       end do

       ! Compute element matrix and rhs
       !call problem%matvec(femsp%lelem(ielem))

       ! Apply boundary conditions
       ! The origional one is in common_tools.f90, e.g.
       ! scalar_conditions_fem_element_new
       ! and takes blocks into account
       !call generic_common_conditions_apply_w_dof_handler(femsp%lelem(ielem),femsp%dofh)  

       ! Assembly contribution
       ! It must be a type bound procedure implementing the abstract interface
       !call assembly(femsp%lelem(ielem), dhand)

       !if(present(int2)) call int2%assembly(femsp%lelem(ielem), dhand)

    end do

  end subroutine volume_integral


  ! SB.alert : integ_element and integ_faces would depend on the fem_space geom, and it would
  ! create a cycle. To be put out ot here. For the moment I have commented it.

  ! ! =================================================================================================
  ! subroutine volume_integrator_fill_val(ielem,geom,integ)
  !   implicit none
  !   ! Parameters
  !   integer(ip)           , intent(in)    :: ielem      ! Element ID
  !   type(fem_space)        , intent(in)    :: geom      ! Geometry interpolation space
  !   type(volume_integrator)       , intent(inout) :: integ      ! Volume integrator of ielem

  !   integer(ip)              :: nnode
  !   real(rp), allocatable    :: elcod(:,:)

  !   ! SB.alert :  I have modified this subroutine assuming that we have a fem_space-like
  !   ! geometry descriptor. I have declared and used geom as a fem_space, but I am not sure
  !   ! we want a full fem_space for the geometry. To be decided.

  !   nnode = femsp%lelem(ielem)%f_inf(ivars)%p%nnode

  !   call memalloc(femsp%g_trian%num_dims,nnode,elcod, __FILE__,__LINE__ )

  !   call gather (femsp%g_trian%num_dims,nnode, geom%lelem(ielem)%elem2dof(1,:), &
  !        & gmesh%lelem(ielem)%unkno,elcod)

  !   ! Define fe map  by interpolation
  !   call femap_from_interp(integ%gint_ref,elcod,integ%femap)

  !   ! Obtain physical interpolation
  !   call femap_apply_to_interp(integ%femap,integ%uint_ref,integ%uint_phy)

  !   call memfree(elcod,__FILE__,__LINE__)

  ! end subroutine integ_element

  ! !==================================================================================================
  ! subroutine integ_faces(elem,face,gmesh,integ,ntxob,nobje)
  !   implicit none
  !   ! Parameters
  !   integer(ip)              , intent(in)    :: elem(2),face(2)
  !   type(fem_space)           , intent(in)    :: geom     ! Geometry interpolation space
  !   type(face_integrator)         , intent(inout) :: integ
  !   type(list)               , intent(in)    :: ntxob
  !   integer(ip)              , intent(in)    :: nobje

  !   integer(ip)                      :: i,gfnod,poins(integ%gint_ref%nnode),nelem
  !   real(rp)                         :: TOL=1e-8
  !   real(rp)                         :: facod(gmesh%ndime,integ%gint_ref%nnode)
  !   real(rp)   , allocatable         :: elcod(:,:)
  !   integer(ip), allocatable         :: elpos(:)
  !   integer(ip)                      :: nnode,inode,iobje,jobje

  !   ! Define face map (from reference to physical space) by interpolation
  !   ! TODO: Assert that the first element is not the one which the subfaces are in

  !   nnode = 0
  !   do inode = ntxob%p(nobje+face),ntxob%p(nobje+face+1)-1
  !      nnode = nnode +1
  !      poins(nnode) = geom%lelem(ielem)%elem2dof(1,ntxob%l(inode))
  !   end do

  !   call gather (femsp%g_trian%num_dims,integ%gint_ref%nnode,poins,geom%lelem(ielem)%unkno,facod)
  !   call bomap_from_interp(integ%gint_ref,facod,integ%bomap)

  !   ! Define elements map  by interpolation
  !   call memalloc(femsp%g_trian%num_dims,integ%gfint_ref%nnode,elcod,__FILE__,__LINE__)
  !   call memalloc(integ%gfint_ref%nnode,elpos,__FILE__,__LINE__)
  !   elpos = geom%lelem(ielem)%elem2dof(1,:)
  !   call gather (femsp%g_trian%num_dims,integ%gfint_ref%nnode,elpos,geom%lelem(ielem)%unkno,elcod)
  !   call femap_from_face_interp(integ%gfint_ref,elcod,face,integ%femap)
  !   call femap_face_to_interp(integ%femap,integ%ufint_ref,face,integ%ufint_phy)
  !   call femap_face_to_interp(integ%femap,integ%gfint_ref,face,integ%gfint_phy)
  !   call memfree(elcod,__FILE__,__LINE__)
  !   call memfree(elpos,__FILE__,__LINE__)

  ! end subroutine integ_faces

end module integration_names
