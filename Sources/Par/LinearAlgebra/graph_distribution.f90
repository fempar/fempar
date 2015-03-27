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
module graph_distribution_names
  ! Serial modules
  use types
  use memor
  use sort_names
  use maps_names
  use dof_handler_names
  use fem_space_names
  use hash_table_names
  
  ! Parallel modules
  use par_partition_names
  use par_triangulation_names
  
  implicit none
# include "debug.i90"
  private

  type graph_distribution
     
     integer(ip) :: part
     integer(ip) :: num_parts
     
     integer(ip)              :: npadj       ! Number of adjacent parts
     integer(ip), allocatable :: lpadj(:)    ! List of adjacent parts
     
     integer(ip)              :: max_nparts  ! Maximum number of parts around any vef communication object
     integer(ip)              :: nobjs       ! Number of local vef communication objects
     integer(ip), allocatable :: lobjs(:,:)  ! List of local vef communication objects
     
     type(list)  :: int_objs  ! List of objects on each edge to an adjacent part / Interface_objects
     type(map)   :: omap      ! Objects local to global map

  end type graph_distribution

  ! Types
  public :: graph_distribution
  
  ! Functions
  public :: dof_lobjs_create

  contains

  !=============================================================================
    subroutine dof_lobjs_create(p_trian, femsp, dhand )
      implicit none
      ! Parameters
      type(par_triangulation), intent(inout) :: p_trian

      ! Locals
      integer(ip)               :: i, j, k, iobj, ielem, jelem, iprob, nvapb, l_var, g_var, mater, obje_l, idof, g_dof, g_mat, l_pos, iblock, ivars
      integer(ip)               :: est_max_nparts, est_max_itf_dofs
      integer(ip)               :: ipart, istat, count 
      integer(igp), allocatable :: lst_parts_per_dof_obj (:,:)
      integer(igp), allocatable :: ws_lobjs_temp (:,:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l1 (:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l2 (:)
      type(hash_table_ip_ip)    :: ws_parts_visited
      integer(ip), parameter    :: tbl_length = 100
      type(fem_space)           :: femsp
      type(dof_handler)         :: dhand
      integer(ip), allocatable  :: touch(:,:,:)

      ipart = p_trian%p_context%iam + 1

      ! Compute an estimation (upper bound) of the maximum number of parts around any local interface vef.
      ! This estimation assumes that all elements around all local interface vefs are associated to different parts.
      est_max_nparts = 0
      est_max_itf_dofs = 0
      do i=1, p_trian%num_itfc_objs
         iobj = p_trian%lst_itfc_objs(i)
         est_max_nparts = max(p_trian%f_trian%objects(iobj)%num_elems_around, est_max_nparts)
         est_max_itf_dofs = est_max_itf_dofs + femsp%object2dof%p(iobj+1) - femsp%object2dof%p(iobj)         
      end do

      call memalloc ( est_max_nparts+5, est_max_itf_dofs, lst_parts_per_dof_obj, __FILE__, __LINE__ )

      call memalloc (  femsp%num_materials, dhand%nvars_global, est_max_nparts+2 , touch, __FILE__, __LINE__ )

      do iblock = 1, dhand%nblocks  

         count = 0

         do i = 1, p_trian%num_itfc_objs
            iobj = p_trian%lst_itfc_objs(i)

            touch = 0

            call ws_parts_visited%init(tbl_length)

            do ielem = 1, p_trian%f_trian%objects(iobj)%num_elems_around
               jelem = p_trian%f_trian%objects(iobj)%elems_around(ielem)

               iprob = femsp%lelem(jelem)%problem
               nvapb = dhand%prob_block(iblock,iprob)%nd1
               do ivars = 1, nvapb
                  !l_var = g2l(ivars,iprob)
                  l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                  g_var = dhand%problems(iprob)%l2g_var(l_var)
                  mater = femsp%lelem(jelem)%material ! SB.alert : material can be used as p 
                  do obje_l = 1, p_trian%f_trian%elems(jelem)%num_objects
                     if ( p_trian%f_trian%elems(jelem)%objects(obje_l) == iobj ) exit
                  end do

                  call ws_parts_visited%put(key=p_trian%elems(jelem)%mypart,val=1,stat=istat)
                  if ( istat == now_stored ) then
                     touch(mater,g_var,1) = touch(mater,g_var,1) + 1 ! New part in the counter             
                     touch(mater,g_var,touch(mater,g_var,1)+2) = p_trian%elems(jelem)%mypart           ! Put new part
                  end if
                  touch(mater,g_var,2) = max(touch(mater,g_var,2),p_trian%elems(jelem)%globalID) ! Max elem GID   
               end do

            end do

            ! Sort list of parts in increasing order by part identifiers
            ! This is required by the call to icomp subroutine below 
            do mater = 1, femsp%num_materials
               do g_var = 1, dhand%nvars_global
                  call sort ( touch(mater,g_var,1), touch(mater,g_var,3:(touch(mater,g_var,1)+2)) )
               end do
            end do

            !p_trian%max_nparts = max(p_trian%max_nparts, count)
            call ws_parts_visited%free

            do idof = femsp%object2dof%p(iobj), femsp%object2dof%p(iobj+1)-1

               g_dof = femsp%object2dof%l(idof,1)
               g_var = femsp%object2dof%l(idof,2)  
               g_mat = femsp%object2dof%l(idof,3)

               ! parts
               if ( touch(g_mat,g_var,1) > 1 ) then ! Interface dof
                  
                  count = count + 1
                  lst_parts_per_dof_obj (1,count) = g_dof ! Dof local ID
                  lst_parts_per_dof_obj (2,count) = g_var ! Variable
                  lst_parts_per_dof_obj (3,count) = g_mat ! Material
                  
                  ! Use the local pos of dof in elem w/ max GID to sort
                  l_var = dhand%g2l_vars(g_var,femsp%lelem(touch(g_mat,g_var,2))%problem)
                  l_pos =  local_node( g_dof, iobj, femsp%lelem(touch(g_mat,g_var,2)), l_var, &
                       & p_trian%f_trian%elems(touch(g_mat,g_var,2))%num_objects, &
                       & p_trian%f_trian%elems(touch(g_mat,g_var,2))%objects )
                  lst_parts_per_dof_obj (4,count) = l_pos ! Local pos in Max elem GID

                  lst_parts_per_dof_obj (5,count) = touch(g_mat,g_var,1) ! Number parts 
                  lst_parts_per_dof_obj (6:(touch(g_mat,g_var,1)+5),count) = touch(g_mat,g_var,3:(touch(mater,g_var,1)+2)) ! List parts

               end if
            end do
         end do
      end do

      ! Here all verbatim !!!

    end subroutine dof_lobjs_create


    integer(ip) function local_node( g_dof, iobj, elem, l_var, nobje, objects )
      
      implicit none
      integer(ip) :: g_dof, iobj, l_var, nobje, objects(:)
      type(fem_element) :: elem

      integer(ip) :: inode, obje_l

      do obje_l = 1, nobje
         if ( objects(obje_l) == iobj ) exit
      end do

      do inode = elem%nodes_object(l_var)%p%p(obje_l), &
           &     elem%nodes_object(l_var)%p%p(obje_l+1)-1  
         local_node = elem%nodes_object(l_var)%p%l(inode)
         if ( elem%elem2dof(local_node,l_var) == g_dof ) exit
      end do

    end function local_node

  end module graph_distribution_names
