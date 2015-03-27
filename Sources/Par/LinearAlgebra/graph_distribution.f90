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
     integer(ip)              :: part
     integer(ip)              :: num_parts    
     integer(ip)              :: npadj       ! Number of adjacent parts
     integer(ip), allocatable :: lpadj(:)    ! List of adjacent parts
     
     integer(ip)              :: max_nparts  ! Maximum number of parts around any dof communication object
     integer(ip)              :: nobjs       ! Number of local dof communication objects
     integer(ip), allocatable :: lobjs(:,:)  ! List of local dof communication objects
     
     type(list)               :: int_objs    ! List of objects on each edge to an adjacent part / Interface_objects
     type(map)                :: omap        ! Objects local to global map
  end type graph_distribution

  ! Types
  public :: graph_distribution
  
  ! Functions
  public :: graph_distribution_create

  contains

  !=============================================================================
    subroutine graph_distribution_create(p_trian, femsp, dhand, gdist)
      implicit none
      ! Parameters
      type(par_triangulation)                , intent(in)  :: p_trian
      type(fem_space)                        , intent(in)  :: femsp
      type(dof_handler)                      , intent(in)  :: dhand
      type(graph_distribution) , allocatable , intent(out) :: gdist(:)

      ! Locals
      integer(ip)               :: i, j, k, iobj, ielem, jelem, iprob, nvapb, l_var, g_var, mater, obje_l, idof, g_dof, g_mat, l_pos, iblock, ivars
      integer(ip)               :: est_max_nparts, est_max_itf_dofs
      integer(ip)               :: ipart, nparts, istat, count, nparts_around
      integer(igp), allocatable :: lst_parts_per_dof_obj (:,:)
      integer(igp), allocatable :: ws_lobjs_temp (:,:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l1 (:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l2 (:)
      type(hash_table_ip_ip)    :: ws_parts_visited, ws_parts_visited_all
      integer(ip), parameter    :: tbl_length = 100

      integer(ip), allocatable :: touch(:,:,:)
      integer(ip), allocatable :: ws_parts_visited_list_all (:)

      integer(ip) :: max_nparts, nobjs, npadj

      ! ** IMPORTANT NOTE
      ! An allocatable array storing the local to local new to old
      ! permutation of DoFs. It might be required that we somehow
      ! provide this vector to the caller subroutine, as the latter might
      ! require to permute other data structures as well (e.g., boundary
      ! conditions). Temporarily declared it here as local allocatable array to
      ! let the code compile.
      integer(ip), allocatable :: l2ln2o(:), l2ln2o_ext(:) 

      ! ** IMPORTANT NOTE
      ! nint and nboun will have to store the number of actual interior and boundary
      ! DoFs for each block. They have not been yet computed within this subroutine.
      ! Declared to let the code compile
      integer(ip) :: nint, nboun

      ipart  = p_trian%p_context%iam + 1
      nparts = p_trian%p_context%np

      ! Compute an estimation (upper bound) of the maximum number of parts around any local interface vef.
      ! This estimation assumes that all elements around all local interface vefs are associated to different parts.

      est_max_nparts = 0
      est_max_itf_dofs = 0
      do i=1, p_trian%num_itfc_objs
         iobj = p_trian%lst_itfc_objs(i)
         est_max_nparts = max(p_trian%f_trian%objects(iobj)%num_elems_around, est_max_nparts)       
      end do

      call memalloc ( est_max_nparts+3, est_max_itf_dofs, lst_parts_per_dof_obj, __FILE__, __LINE__ )
      call memalloc ( femsp%num_materials, dhand%nvars_global, est_max_nparts+2 , touch, __FILE__, __LINE__ )
      call memalloc ( est_max_nparts+3, sort_parts_per_itfc_obj_l1, __FILE__,__LINE__  )  
      call memalloc ( est_max_nparts+3, sort_parts_per_itfc_obj_l2, __FILE__,__LINE__  )

      ! The size of the following array does not scale 
      ! properly with nparts. We should think in the meantime
      ! a strategy to keep bounded (below a local quantity) its size 
      call memalloc ( nparts, ws_parts_visited_list_all, __FILE__, __LINE__ )

      allocate(gdist(dhand%nblocks), stat=istat)
      check(istat==0)

      do iblock = 1, dhand%nblocks  

         est_max_itf_dofs = 0
         do i=1, p_trian%num_itfc_objs
            iobj = p_trian%lst_itfc_objs(i)
            est_max_itf_dofs = est_max_itf_dofs + femsp%object2dof(iblock)%p(iobj+1) &
                 & - femsp%object2dof(iblock)%p(iobj)
         end do
         ! l2ln2o interior vefs
         nint = 0
         nboun = 0
         do iobj = 1, p_trian%f_trian%num_objects
            if ( p_trian%objects(iobj)%interface /= -1 ) then
               do idof = femsp%object2dof(iblock)%p(iobj), femsp%object2dof(iblock)%p(iobj+1)-1
                  nint = nint + 1
                  l2ln2o(nint) = femsp%object2dof(iblock)%l(idof,1)                  
               end do
            end if
         end do

         call ws_parts_visited_all%init(tbl_length)

         ! Initialize trivial components of gdist(iblock)
         gdist(iblock)%part      = ipart
         gdist(iblock)%num_parts = nparts
         npadj     = 0

         count = 0
         max_nparts = 0
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
                     touch(mater,g_var,touch(mater,g_var,1)+2) = p_trian%elems(jelem)%mypart ! Put new part
                  end if
                  touch(mater,g_var,2) = max(touch(mater,g_var,2),p_trian%elems(jelem)%globalID) ! Max elem GID 

                  call ws_parts_visited_all%put(key=p_trian%elems(jelem)%mypart,val=1,stat=istat)
                  if ( istat == now_stored ) then
                     npadj = npadj + 1
                     ws_parts_visited_list_all(npadj) = p_trian%elems(jelem)%mypart
                  end if
               end do
            end do
            max_nparts = max(max_nparts, touch(mater,g_var,1))

            ! Sort list of parts in increasing order by part identifiers
            ! This is required by the call to icomp subroutine below 
            do mater = 1, femsp%num_materials
               do g_var = 1, dhand%nvars_global
                  call sort ( touch(mater,g_var,1), touch(mater,g_var,3:(touch(mater,g_var,1)+2)) )
               end do
            end do

            !p_trian%max_nparts = max(p_trian%max_nparts, count)
            call ws_parts_visited%free

            do idof = femsp%object2dof(iblock)%p(iobj), femsp%object2dof(iblock)%p(iobj+1)-1
               ! parts
               if ( touch(g_mat,g_var,1) > 1 ) then ! Interface dof
                  g_dof = femsp%object2dof(iblock)%l(idof,1)
                  g_var = femsp%object2dof(iblock)%l(idof,2)  
                  g_mat = femsp%object2dof(iblock)%l(idof,3)
                  nparts_around = touch(g_mat,g_var,1)

                  count = count + 1
                  lst_parts_per_dof_obj (1,count) = g_var ! Variable

                  ! Use the local pos of dof in elem w/ max GID to sort
                  l_var = dhand%g2l_vars(g_var,femsp%lelem(touch(g_mat,g_var,2))%problem)
                  l_pos =  local_node( g_dof, iobj, femsp%lelem(touch(g_mat,g_var,2)), l_var, &
                       & p_trian%f_trian%elems(touch(g_mat,g_var,2))%num_objects, &
                       & p_trian%f_trian%elems(touch(g_mat,g_var,2))%objects )

                  lst_parts_per_dof_obj (2,count) = nparts_around ! Number parts 
                  lst_parts_per_dof_obj (3:(nparts_around+2),count) = touch(g_mat,g_var,3:(touch(mater,g_var,1)+2)) ! List parts
                  lst_parts_per_dof_obj (nparts_around+3:est_max_nparts+2,count) = 0 ! Zero the rest of entries in current col except last
                  lst_parts_per_dof_obj (est_max_nparts+3,count) = l_pos ! Local pos in Max elem GID

                  ! l2ln2o interface vefs to interface dofs
                  nboun = nboun + 1
                  l2ln2o_ext(nboun) = femsp%object2dof(iblock)%l(idof,1)                  
               else
                  ! l2ln2o interface vefs to interior dofs
                  nint = nint + 1
                  l2ln2o(nint) = femsp%object2dof(iblock)%l(idof,1)  
               end if
            end do
         end do

         ! l2ln2o interface vefs to interface dofs from l2ln20_ext
         l2ln2o(nint+1:nint+nboun-1) = l2ln2o_ext

         nint = nint
         nboun = nboun
         call memfree( l2ln2o_ext, __FILE__, __LINE__ )

         call ws_parts_visited_all%free()

         ! Initialize max_nparts for gdist(iblock)%npadj
         gdist(iblock)%max_nparts = max_nparts

         ! Initialize npadj/lpadj for gdist(iblock)%npadj
         gdist(iblock)%npadj = npadj
         call memalloc ( gdist(iblock)%npadj, gdist(iblock)%lpadj, __FILE__, __LINE__ )
         gdist(iblock)%lpadj = ws_parts_visited_list_all(1:gdist(iblock)%npadj)

         ! Re-number boundary DoFs in increasing order by physical unknown identifier, the 
         ! number of parts they belong and, for DoFs sharing the same number of parts,
         ! in increasing order by the list of parts shared by each DoF.
         ! ** IMPORTANT NOTE
         !    The computation of nint/nboun is pending. nboun is actually equal to count, right?
         !    l2ln2o_dofs has to be allocated and properly initialized as well
         call sort_array_cols_by_row_section( est_max_nparts+3,            & ! #Rows 
              &                                 est_max_nparts+3,            & ! Leading dimension
              &                                 nboun,                       & ! #Cols
              &                                 lst_parts_per_dof_obj,       &
              &                                 l2ln2o(nint+1:),             &
              &                                 sort_parts_per_itfc_obj_l1,  &
              &                                 sort_parts_per_itfc_obj_l2)

         ! Identify interface communication objects 
         call memalloc ( max_nparts+4, nboun, ws_lobjs_temp, __FILE__,__LINE__ )
         ws_lobjs_temp = 0
         if (nboun == 0) then
            nobjs=0
         else
            nobjs=1
            ws_lobjs_temp(2,nobjs) = nint+1 ! begin obj

            k = nint+1
            ! Loop over vertices of the current separator (i.e., inode)
            do i=1,nboun-1
               if(icomp( max_nparts+2,lst_parts_per_dof_obj(:,i),lst_parts_per_dof_obj(:,i+1)) /= 0) then
                  ! Complete current object
                  ws_lobjs_temp(1,nobjs)=lst_parts_per_dof_obj(1,i)
                  ws_lobjs_temp(3,nobjs)=k ! end obj
                  ws_lobjs_temp(4,nobjs)=lst_parts_per_dof_obj(2,i)
                  ws_lobjs_temp(5:4+lst_parts_per_dof_obj(2,i),nobjs) = &
                       & lst_parts_per_dof_obj(3:2+lst_parts_per_dof_obj(2,i),i)

                  ! Prepare next object
                  nobjs=nobjs+1
                  ws_lobjs_temp(2,nobjs)=k+1 ! begin obj
               end if
               k=k+1
            end do

            ! Complete last object
            ws_lobjs_temp(1,nobjs)=lst_parts_per_dof_obj(1,i)
            ws_lobjs_temp(3,nobjs)=k ! end obj
            ws_lobjs_temp(4,nobjs)=lst_parts_per_dof_obj(2,i)
            ws_lobjs_temp(5:4+lst_parts_per_dof_obj(2,i),nobjs) = &
                 & lst_parts_per_dof_obj(3:2+lst_parts_per_dof_obj(2,i),i)
         end if

         ! Reallocate lobjs and add internal object first
         nobjs = nobjs + 1
         call memalloc (max_nparts+4, nobjs, gdist(iblock)%lobjs,__FILE__,__LINE__)
         gdist(iblock)%lobjs = 0


         ! Internal object
         gdist(iblock)%lobjs (1,1) = -1
         gdist(iblock)%lobjs (2,1) = 1
         gdist(iblock)%lobjs (3,1) = nint
         gdist(iblock)%lobjs (4,1) = 1
         gdist(iblock)%lobjs (5,1) = ipart
         gdist(iblock)%lobjs (6:max_nparts+4,1) = 0

         ! Copy ws_lobjs_temp to lobjs ...
         do i=1,nobjs-1
            gdist(iblock)%lobjs(:,i+1)=ws_lobjs_temp(1:(max_nparts+4), i)
         end do

         call memfree ( ws_lobjs_temp,__FILE__,__LINE__)

!!$         ** IMPORTANT NOTE
!!$         The update of object2dof(iblock) and elem2dof conformally with l2ln2o and aux is pending (see code commented below)
!!$         ! Auxiliary inverse of l2ln2o
!!$         call memalloc (ndofs, aux, __FILE__,__LINE__)
!!$                  
!!$         do i=1,ndofs
!!$            aux(l2ln2o(i)) = i
!!$         end do   
!!$         ! Update object2dof(iblock)
!!$         do i = 1,f_mesh%npoin
!!$            do j=object2dof(iblock)_p(i),object2dof(iblock)_p(i+1)-1
!!$               object2dof(iblock)_l(j,1) = aux(object2dof(iblock)_l(j,1))
!!$               dof2object_l(object2dof(iblock)_l(j,1)) = i
!!$            end do
!!$         end do
!!$         
!!$         ! Update elem2dof
!!$         do ielem = 1,f_mesh%nelem
!!$            do k = i_varsxprob(femsp%lelem(ielem)%prob), i_varsxprob(femsp%lelem(ielem)%prob+1)-1
!!$               do j = 1,femsp%lelem(ielem)%f_inf(femsp%lelem(ielem)%iv(j_varsxprob(k)))%p%nnode 
!!$                  if ( femsp%lelem(ielem)%elem2dof(j,j_varsxprob(k)) /= 0 ) then
!!$                     femsp%lelem(ielem)%elem2dof(j,j_varsxprob(k)) = aux(femsp%lelem(ielem)%elem2dof(j,j_varsxprob(k)))
!!$                  else
!!$                     femsp%lelem(ielem)%elem2dof(j,j_varsxprob(k)) = 0
!!$                  end if
!!$               end do
!!$            end do
!!$         end do    
!!$         call memfree ( aux,__FILE__,__LINE__)

!!$    call create_int_objs_new ( me+1, &
!!$                               new_f_part%npadj, &
!!$                               new_f_part%lpadj, &
!!$                               new_f_part%max_nparts , &
!!$                               new_f_part%nobjs      , &
!!$                               new_f_part%lobjs      , &
!!$                               new_f_part%int_objs%n , &
!!$                               new_f_part%int_objs%p , &
!!$                               new_f_part%int_objs%l ) 
!!$
!!$        call create_omap ( icontxt    , & ! Communication context
!!$                       me         , &
!!$                       np         , &
!!$                       new_f_part%npadj, &
!!$                       new_f_part%lpadj, & 
!!$                       new_f_part%int_objs%p, &  
!!$                       new_f_part%int_objs%l, &
!!$                       new_f_part%max_nparts , & 
!!$                       new_f_part%nobjs      , & 
!!$                       new_f_part%lobjs      , &
!!$                       new_f_part%omap%nl,     &
!!$                       new_f_part%omap%ng,     &
!!$                       new_f_part%omap%ni,     &
!!$                       new_f_part%omap%nb,     &
!!$                       new_f_part%omap%ne,     &
!!$                       new_f_part%omap%l2g )
      end do


      call memfree ( sort_parts_per_itfc_obj_l1, __FILE__,__LINE__  )  
      call memfree ( sort_parts_per_itfc_obj_l2, __FILE__,__LINE__  )
      call memfree ( lst_parts_per_dof_obj, __FILE__, __LINE__ )
      call memfree ( touch, __FILE__, __LINE__ )
      call memfree ( ws_parts_visited_list_all, __FILE__, __LINE__ )

    end subroutine graph_distribution_create


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
