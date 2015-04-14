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
  use fem_space_types
  use fem_space_names
  use hash_table_names
  
  ! Parallel modules
  use par_partition_names
  use par_triangulation_names
  use psb_penv_mod

  
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
      type(fem_space)                        , intent(inout)  :: femsp
      type(dof_handler)                      , intent(in)  :: dhand
      type(graph_distribution) , allocatable , intent(out) :: gdist(:)

      ! Locals
      integer(ip)               :: i, j, k, iobj, ielem, jelem, iprob, nvapb, l_var, g_var, mater, obje_l, idof, g_dof, g_mat, l_pos, iblock, ivars
      integer(ip)               :: est_max_nparts, est_max_itf_dofs
      integer(ip)               :: ipart, nparts, istat, count, nparts_around, touching
      integer(igp), allocatable :: lst_parts_per_dof_obj (:,:)
      integer(igp), allocatable :: ws_lobjs_temp (:,:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l1 (:)
      integer(igp), allocatable :: sort_parts_per_itfc_obj_l2 (:)
      type(hash_table_ip_ip)    :: ws_parts_visited, ws_parts_visited_all
      integer(ip), parameter    :: tbl_length = 100

      integer(ip), allocatable :: touch(:,:,:)
      integer(ip), allocatable :: ws_parts_visited_list_all (:)

      integer(ip) :: max_nparts, nobjs, npadj

      integer(ip), allocatable :: l2ln2o(:), l2ln2o_ext(:), l2lo2n(:) 

      integer(ip) :: nint, nboun, l_node, inode, iobje

      integer(ip) :: elem_ghost, elem_local, l_var_ghost, l_var_local, nnode, obje_ghost, obje_local, order 
      integer(ip) :: o2n(max_nnode)

      ipart  = p_trian%p_context%iam + 1
      nparts = p_trian%p_context%np

      ! Compute an estimation (upper bound) of the maximum number of parts around any local interface vef.
      ! This estimation assumes that all elements around all local interface vefs are associated to different parts.


      ! if ( p_trian%p_context%iam > 0 ) then
      !    !pause
      !    do while ( 1 > 0)
      !       i = i + 1
      !    end do
      ! else
      !    write (*,*) 'Processor 0 not stopped'
      !    !i = 1
      !    !do while ( 1 > 0)
      !    !   i = i + 1
      !    !end do
      ! end if



      est_max_nparts = 0
      est_max_itf_dofs = 0
      do i=1, p_trian%num_itfc_objs
         iobj = p_trian%lst_itfc_objs(i)
         est_max_nparts = max(p_trian%f_trian%objects(iobj)%num_elems_around, est_max_nparts)       
      end do

      ! The size of the following array does not scale 
      ! properly with nparts. We should think in the meantime
      ! a strategy to keep bounded (below a local quantity) its size 

      allocate(gdist(dhand%nblocks), stat=istat)
      check(istat==0)

      !write (*,*) 'GRAPH DISTRIBUTION****'

      do iblock = 1, dhand%nblocks  

         est_max_itf_dofs = 0
         do i=1, p_trian%num_itfc_objs
            iobj = p_trian%lst_itfc_objs(i)
            est_max_itf_dofs = est_max_itf_dofs + femsp%object2dof(iblock)%p(iobj+1) &
                 & - femsp%object2dof(iblock)%p(iobj)
         end do

         !write (*,*) 'est_max_nparts****',est_max_nparts
         !write (*,*) 'est_max_itf_dofs****',est_max_itf_dofs
         call memalloc ( nparts, ws_parts_visited_list_all, __FILE__, __LINE__ )

         call memalloc ( est_max_nparts+4, est_max_itf_dofs, lst_parts_per_dof_obj, __FILE__, __LINE__ )
         call memalloc ( femsp%num_continuity, dhand%nvars_global, est_max_nparts+5 , touch, __FILE__, __LINE__ )
         call memalloc ( est_max_nparts+4, sort_parts_per_itfc_obj_l1, __FILE__,__LINE__  )  
         call memalloc ( est_max_nparts+4, sort_parts_per_itfc_obj_l2, __FILE__,__LINE__  )

         call memalloc ( femsp%ndofs(iblock), l2ln2o, __FILE__, __LINE__ )
         call memalloc ( est_max_itf_dofs, l2ln2o_ext, __FILE__, __LINE__ )

         ! l2ln2o interior vefs
         nint = 0
         nboun = 0
         do iobj = 1, p_trian%f_trian%num_objects
            if ( p_trian%objects(iobj)%interface == -1 ) then
               do idof = femsp%object2dof(iblock)%p(iobj), femsp%object2dof(iblock)%p(iobj+1)-1
                  nint = nint + 1
                  l2ln2o(nint) = femsp%object2dof(iblock)%l(idof,1) 
                  !write(*,*) 'l2ln2o',l2ln2o
               end do
            end if
         end do

         !write (*,*) 'nint',nint

         ! interior
         if ( .not. femsp%static_condensation ) then 
            do ielem = 1, p_trian%f_trian%num_elems
               iprob = femsp%lelem(ielem)%problem
               nvapb = dhand%prob_block(iblock,iprob)%nd1
               do ivars = 1, nvapb
                  l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                  g_var = dhand%problems(iprob)%l2g_var(l_var)  
                  iobje = p_trian%f_trian%elems(ielem)%num_objects+1
                  do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                       &     femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1 
                     l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                     nint = nint + 1
                     l2ln2o(nint) = femsp%lelem(ielem)%elem2dof(l_node,l_var) 
                  end do
               end do
            end do
         end if

         !write (*,*) 'nint',nint

         call ws_parts_visited_all%init(tbl_length)

         ! Initialize trivial components of gdist(iblock)
         gdist(iblock)%part      = ipart
         gdist(iblock)%num_parts = nparts
         npadj     = 0

         touching = 1
         count = 0
         max_nparts = 0
         do i = 1, p_trian%num_itfc_objs
            iobj = p_trian%lst_itfc_objs(i)
            !write (*,*) '***** OBJECT *****',iobj

            touch = 0
            call ws_parts_visited%init(tbl_length)

            do ielem = 1, p_trian%f_trian%objects(iobj)%num_elems_around
               jelem = p_trian%f_trian%objects(iobj)%elems_around(ielem)

               iprob = femsp%lelem(jelem)%problem
               nvapb = dhand%prob_block(iblock,iprob)%nd1
               do ivars = 1, nvapb
                  !l_var = g2l(ivars,iprob)
                  l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                  if ( femsp%lelem(jelem)%continuity(l_var) /= 0 ) then
                     g_var = dhand%problems(iprob)%l2g_var(l_var)
                     mater = femsp%lelem(jelem)%material ! SB.alert : material can be used as p 
                     do obje_l = 1, p_trian%f_trian%elems(jelem)%num_objects
                        if ( p_trian%f_trian%elems(jelem)%objects(obje_l) == iobj ) exit
                     end do

                     !call ws_parts_visited%put(key=p_trian%elems(jelem)%mypart,val=1,stat=istat)
                     call ws_parts_visited%put(key=p_trian%elems(jelem)%mypart,val=touching,stat=istat)
                     if ( istat == now_stored ) then
                        touch(mater,g_var,1) = touch(mater,g_var,1) + 1 ! New part in the counter  
                        !write(*,*) 'touch',touch
                        touch(mater,g_var,touch(mater,g_var,1)+5) = p_trian%elems(jelem)%mypart ! Put new part
                     end if
                     if( p_trian%elems(jelem)%globalID > touch(mater,g_var,2) ) then
                        touch(mater,g_var,2) = jelem 
                        touch(mater,g_var,3) = obje_l
                     end if
                     if ( jelem <= p_trian%f_trian%num_elems ) then
                        touch(mater,g_var,4) = jelem 
                        touch(mater,g_var,5) = obje_l
                     end if
                     if ( p_trian%elems(jelem)%mypart /= ipart ) then
                        !call ws_parts_visited_all%put(key=p_trian%elems(jelem)%mypart,val=1,stat=istat)
                        call ws_parts_visited_all%put(key=p_trian%elems(jelem)%mypart,val=touching,stat=istat)
                        if ( istat == now_stored ) then
                           npadj = npadj + 1
                           ws_parts_visited_list_all(npadj) = p_trian%elems(jelem)%mypart
                        end if
                     end if
                  end if
               end do
            end do

            max_nparts = max(max_nparts, touch(mater,g_var,1))

            !write(*,*) '**** TOUCH ****',touch


            ! Sort list of parts in increasing order by part identifiers
            ! This is required by the call to icomp subroutine below 
            do mater = 1, femsp%num_continuity
               do g_var = 1, dhand%nvars_global
                  call sort ( touch(mater,g_var,1), touch(mater,g_var,6:(touch(mater,g_var,1)+5)) )
               end do
            end do

            !p_trian%max_nparts = max(p_trian%max_nparts, count)
            call ws_parts_visited%free

            do idof = femsp%object2dof(iblock)%p(iobj), femsp%object2dof(iblock)%p(iobj+1)-1
               ! parts
               g_var = femsp%object2dof(iblock)%l(idof,2)  
               g_mat = femsp%object2dof(iblock)%l(idof,3)
               !write(*,*) 'g_mat',g_mat
               if ( touch(g_mat,g_var,1) > 1 ) then ! Interface dof
                  g_dof = femsp%object2dof(iblock)%l(idof,1)
                  nparts_around = touch(g_mat,g_var,1)

                  count = count + 1
                  lst_parts_per_dof_obj (1,count) = g_var ! Variable

                  ! Use the local pos of dof in elem w/ max GID to sort
                  l_var = dhand%g2l_vars(g_var,femsp%lelem(touch(g_mat,g_var,2))%problem)

                  if ( touch(g_mat,g_var,2) <= p_trian%f_trian%num_elems ) then
                     l_pos =  local_node( g_dof, iobj, femsp%lelem(touch(g_mat,g_var,2)), l_var, &
                          & p_trian%f_trian%elems(touch(g_mat,g_var,2))%num_objects, &
                          & p_trian%f_trian%elems(touch(g_mat,g_var,2))%objects )
                  else
                     l_pos =  local_node( g_dof, iobj, femsp%lelem(touch(g_mat,g_var,4)), l_var, &
                          & p_trian%f_trian%elems(touch(g_mat,g_var,4))%num_objects, &
                          & p_trian%f_trian%elems(touch(g_mat,g_var,4))%objects )
                     ! SB.alert :  This part is difficult... can I simplify it ?
                     elem_ghost = touch(g_mat,g_var,2)
                     elem_local = touch(g_mat,g_var,4)
                     l_var_ghost = dhand%g2l_vars(g_var,femsp%lelem(elem_ghost)%problem)
                     l_var_local = dhand%g2l_vars(g_var,femsp%lelem(elem_local)%problem)
                     obje_ghost = touch(g_mat,g_var,3)
                     obje_local = touch(g_mat,g_var,5)
                     order = femsp%lelem(elem_local)%order(l_var_local)
                     ! SB.alert : We must think about it and the relation between material and interpolation
                     nnode = femsp%lelem(elem_local)%nodes_object(l_var_local)%p%p(obje_local+1) &
                          &  -femsp%lelem(elem_local)%nodes_object(l_var_local)%p%p(obje_local)

                     if ( p_trian%f_trian%objects(iobj)%dimension == p_trian%f_trian%num_dims .and. &
                          & nnode ==  (order+1)**p_trian%f_trian%num_dims ) then
                        order = order    ! hdG case
                     elseif ( nnode ==  (order-1)**p_trian%f_trian%objects(iobj)%dimension ) then
                        order = order -2 ! cG case
                     else
                        assert ( 0 == 1) ! SB.alert : Other situations possible when dG_material, cdG, hp-adaptivity ?
                     end if
                     assert( femsp%lelem(elem_local)%order(l_var_local) == femsp%lelem(elem_ghost)%order(l_var_ghost)  )

                     call permute_nodes_object(                                                                 &
                          & femsp%lelem(elem_local)%f_inf(l_var_local)%p,                                                  &
                          & femsp%lelem(elem_ghost)%f_inf(l_var_ghost)%p,                                           &
                          & o2n,obje_local,obje_ghost,                                                                &
                          & p_trian%f_trian%elems(elem_local)%objects,                                                         &
                          & p_trian%f_trian%elems(elem_ghost)%objects,                                                      &
                          & p_trian%f_trian%objects(iobj)%dimension,                                                     &
                          & order )
                     l_pos = o2n(l_pos)
                  end if
                  lst_parts_per_dof_obj (2,count) = nparts_around ! Number parts 
                  lst_parts_per_dof_obj (3:(nparts_around+2),count) = touch(g_mat,g_var,6:(touch(g_mat,g_var,1)+5)) ! List parts
                  lst_parts_per_dof_obj (nparts_around+3:est_max_nparts+2,count) = 0 ! Zero the rest of entries in current col except last
                  lst_parts_per_dof_obj (est_max_nparts+3,count) = p_trian%objects(iobj)%globalID
                  lst_parts_per_dof_obj (est_max_nparts+4,count) = l_pos ! Local pos in Max elem GID

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
         l2ln2o(nint+1:nint+nboun) = l2ln2o_ext(1:nboun)
         write (*,*) 'nint,nboun,ndofs',nint,nboun,femsp%ndofs(iblock) 
         assert( nint + nboun == femsp%ndofs(iblock) )  ! check
         call memfree( l2ln2o_ext, __FILE__, __LINE__ )

         write (*,*) 'l2ln2o:',l2ln2o

         call ws_parts_visited_all%free()

         ! Initialize max_nparts for gdist(iblock)%npadj
         gdist(iblock)%max_nparts = max_nparts

         ! Initialize npadj/lpadj for gdist(iblock)%npadj
         gdist(iblock)%npadj = npadj
         call memalloc ( gdist(iblock)%npadj, gdist(iblock)%lpadj, __FILE__, __LINE__ )
         gdist(iblock)%lpadj = ws_parts_visited_list_all(1:gdist(iblock)%npadj)

         ! write (*,*) 'npadj=', gdist(iblock)%npadj, 'lpadj=', gdist(iblock)%lpadj

         ! Re-number boundary DoFs in increasing order by physical unknown identifier, the 
         ! number of parts they belong and, for DoFs sharing the same number of parts,
         ! in increasing order by the list of parts shared by each DoF.
         call sort_array_cols_by_row_section( est_max_nparts+4,            & ! #Rows 
              &                                 est_max_nparts+4,            & ! Leading dimension
              &                                 nboun,                       & ! #Cols
              &                                 lst_parts_per_dof_obj,       &
              &                                 l2ln2o(nint+1:),             &
              &                                 sort_parts_per_itfc_obj_l1,  &
              &                                 sort_parts_per_itfc_obj_l2)

         ! write (*,*) 'l2ln2o:',l2ln2o

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
         gdist(iblock)%nobjs = nobjs
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

         ! Auxiliary inverse of l2ln2o
         call memalloc ( femsp%ndofs(iblock), l2lo2n, __FILE__,__LINE__)

         do i=1, femsp%ndofs(iblock)
            l2lo2n(l2ln2o(i)) = i
         end do

         ! Update object2dof(iblock)
         do i = 1,femsp%object2dof(iblock)%p(p_trian%f_trian%num_objects+1)-1
            femsp%object2dof(iblock)%l(i,1) = l2lo2n(femsp%object2dof(iblock)%l(i,1))
         end do

         do ielem = 1, p_trian%f_trian%num_elems
            iprob = femsp%lelem(ielem)%problem
            nvapb = dhand%prob_block(iblock,iprob)%nd1
            do ivars = 1, nvapb
               l_var = dhand%prob_block(iblock,iprob)%a(ivars)
               do inode = 1,femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                  if ( femsp%lelem(ielem)%elem2dof(inode,l_var) > 0 ) then 
                     femsp%lelem(ielem)%elem2dof(inode,l_var) = l2lo2n(femsp%lelem(ielem)%elem2dof(inode,l_var))
                  end if
               end do
            end do
         end do

         call memfree ( l2lo2n,__FILE__,__LINE__)


    call create_int_objs ( ipart, &
                           gdist(iblock)%npadj, &
                           gdist(iblock)%lpadj, &
                           gdist(iblock)%max_nparts , &
                           gdist(iblock)%nobjs      , &
                           gdist(iblock)%lobjs      , &
                           gdist(iblock)%int_objs%n , &
                           gdist(iblock)%int_objs%p , &
                           gdist(iblock)%int_objs%l ) 

    call create_omap ( p_trian%p_context%icontxt    , & ! Communication context
                       p_trian%p_context%iam         , &
                       p_trian%p_context%np         , &
                       gdist(iblock)%npadj, &
                       gdist(iblock)%lpadj, & 
                       gdist(iblock)%int_objs%p, &  
                       gdist(iblock)%int_objs%l, &
                       gdist(iblock)%max_nparts , & 
                       gdist(iblock)%nobjs      , & 
                       gdist(iblock)%lobjs      , &
                       gdist(iblock)%omap%nl,     &
                       gdist(iblock)%omap%ng,     &
                       gdist(iblock)%omap%ni,     &
                       gdist(iblock)%omap%nb,     &
                       gdist(iblock)%omap%ne,     &
                       gdist(iblock)%omap%l2g )

         call memfree ( sort_parts_per_itfc_obj_l1, __FILE__,__LINE__  )  
         call memfree ( sort_parts_per_itfc_obj_l2, __FILE__,__LINE__  )
         call memfree ( lst_parts_per_dof_obj, __FILE__, __LINE__ )
         call memfree ( touch, __FILE__, __LINE__ )
         call memfree ( ws_parts_visited_list_all, __FILE__, __LINE__ )

      end do

 end subroutine graph_distribution_create



    

    integer(ip) function local_node( g_dof, iobj, elem, l_var, nobje, objects )
      
      implicit none
      integer(ip) :: g_dof, iobj, l_node, l_var, nobje, objects(:)
      type(fem_element) :: elem

      integer(ip) :: inode, obje_l

      do obje_l = 1, nobje
         if ( objects(obje_l) == iobj ) exit
      end do

      do inode = elem%nodes_object(l_var)%p%p(obje_l), &
           &     elem%nodes_object(l_var)%p%p(obje_l+1)-1  
         l_node = elem%nodes_object(l_var)%p%l(inode)
         if ( elem%elem2dof(l_node,l_var) == g_dof ) then
            local_node = inode - elem%nodes_object(l_var)%p%p(obje_l) + 1
            exit
         end if
      end do

    end function local_node

    subroutine create_int_objs ( ipart      , &
                                 npadj      , &
                                 lpadj      , & 
                                 max_nparts , &
                                 nobjs      , & 
                                 lobjs      , &
                                 n          , &
                                 p          , &
                                 l )
      implicit none
      integer(ip)   , intent(in)  :: ipart 
      integer(ip)   , intent(in)  :: npadj 
      integer(ip)   , intent(in)  :: lpadj(npadj)
      integer(ip)   , intent(in)  :: max_nparts
      integer(ip)   , intent(in)  :: nobjs
      integer(ip)   , intent(in)  :: lobjs(max_nparts+4,nobjs)
      integer(ip)   , intent(out) :: n
      integer(ip)   , intent(out), allocatable :: p(:)
      integer(ip)   , intent(out), allocatable :: l(:) 


!!$
!!$  ! ==========================================
!!$  ! BEGIN Compute int_objs from npadj/lpadj
!!$  ! ==========================================

      ! Locals
      integer(ip) :: i, iedge, iobj, j, jpart, iposobj
      integer(ip), allocatable :: ws_elems_list (:,:)
      integer(ip), allocatable :: ws_sort_l1 (:), ws_sort_l2 (:)

      n = npadj   
      call memalloc ( n+1, p,         __FILE__,__LINE__)

      ! Count objects on each edge of the graph of parts
      p = 0

      do iobj=1,nobjs
         do j=1,lobjs(4,iobj)
            jpart = lobjs(4+j,iobj)   ! SB.alert : Error index '9' max 8
            if (ipart /= jpart) then  ! Exclude self-edges 
               ! Locate edge identifier of ipart => jpart on the list of edges of ipart 
               iedge = 1
               do while ((lpadj(iedge) /= jpart).and. &
                    &      (iedge < npadj) )
                  iedge = iedge + 1
               end do

               ! Increment by one the number of objects on the edge ipart => jpart
               p(iedge+1) = p(iedge+1) + 1 
            end if
         end do
      end do

      ! Generate pointers to the list of objects on each edge
      p(1) = 1
      do iedge=1, n
         p(iedge+1)=p(iedge+1)+p(iedge)
      end do

      ! Generate list of objects on each edge
      call memalloc (p(n+1)-1, l,            __FILE__,__LINE__)

      do iobj=1,nobjs
         do j=1,lobjs(4,iobj)
            jpart = lobjs(4+j,iobj)
            if (ipart /= jpart) then  ! Exclude self-edges 
               ! Locate edge identifier of ipart => jpart on the list of edges of ipart 
               iedge = 1
               do while ((lpadj(iedge) /= jpart).and. &
                    &      (iedge < npadj) )
                  iedge = iedge + 1
               end do

               ! Insert current object on the edge ipart => jpart
               l(p(iedge)) = iobj
               p(iedge) = p(iedge) + 1 
            end if
         end do
      end do

      ! Recover p
      do iedge=n+1, 2, -1
         p(iedge) = p(iedge-1)
      end do
      p(iedge) = 1



      call memalloc ( 2, p(n+1)-1, ws_elems_list,         __FILE__,__LINE__ )

      ws_elems_list = 0 
      do iedge=1, n
         do iposobj=p(iedge),p(iedge+1)-1
            ws_elems_list(1,iposobj) = iedge
            ws_elems_list(2,iposobj) = l(iposobj)
         end do
      end do

      call memalloc ( 2, ws_sort_l1,         __FILE__,__LINE__ )

      call memalloc ( 2, ws_sort_l2,         __FILE__,__LINE__ )

      call intsort( 2,   & ! Rows of ws_parts_list_sep
           &          2,   & ! LD of ws_parts_list_sep
           &          p(n+1)-1,      & ! Cols of ws_parts_list_sep
           &          ws_elems_list,  &
           &          l, &
           &          ws_sort_l1, ws_sort_l2)



      call memfree ( ws_sort_l1,__FILE__,__LINE__)

      call memfree ( ws_sort_l2,__FILE__,__LINE__)

      call memfree ( ws_elems_list,__FILE__,__LINE__)

      write(*,'(a)') 'List of interface objects:'
      do i=1,npadj
         write(*,'(10i10)') i, &
              & (l(j),j=p(i),p(i+1)-1)
      end do

      ! ========================================
      ! END Compute int_objs from npadj/lpadj
      ! ========================================

    end subroutine create_int_objs

    subroutine create_omap ( icontxt    , & ! Communication context
                             me         , &
                             np         , &
                             npadj      , &
                             lpadj      , &
                             int_objs_p , &
                             int_objs_l , &
                             max_nparts , & 
                             nobjs      , & 
                             lobjs      , &
                             onl        , &
                             ong        , &
                             oni        , &
                             onb        , &
                             one        , &
                             ol2g       )

#ifdef MPI_MOD
      use mpi
#endif     
      implicit none
#ifdef MPI_H
      include 'mpif.h'
#endif
      
    integer(ip), intent(in)               :: icontxt, me, np 
    integer(ip), intent(in)               :: npadj 
    integer(ip), intent(in)               :: lpadj(npadj)
    integer(ip), intent(in)               :: int_objs_p (npadj+1)
    integer(ip), intent(in)               :: int_objs_l (int_objs_p(npadj+1)-1)
    integer(ip), intent(in)               :: max_nparts, nobjs
    integer(ip), intent(in)               :: lobjs(max_nparts+4,nobjs)
    integer(ip), intent(out)              :: onl, ong, oni, onb, one
    integer(ip), intent(out), allocatable :: ol2g(:) 

    
    integer(ip) :: num_rcv             ! From how many neighbours does the part receive data ?
    integer(ip), allocatable :: &           
         list_rcv(:),           &      ! From which neighbours does the part receive data ?
         rcv_ptrs(:),           &      ! How much data does the part receive from each neighbour ?
         lids_rcv(:)                 

    integer(ip) :: num_snd             ! To how many neighbours does the part send data ? 
    integer(ip), allocatable :: &           
         list_snd(:),           &      ! To which neighbours does the part send data ?
         snd_ptrs(:),           &      ! How much data does the part send to each neighbour?
         lids_snd(:)


    integer(ip) :: ipart ! Which part am I ?
    integer(ip) :: iedge, idobj, neighbour_part_id
    integer(ip) :: iposobj, visited_snd, visited_rcv, i, j
    integer(ip) :: cur_snd, cur_rcv, sizmsg

    integer(ip)              :: nobjs_with_gid, iobj
    integer(ip), allocatable :: lobjs_with_gid(:)
    integer(ip), allocatable :: ptr_gids(:)
    integer(ip)              :: start_gid

    integer(ip), allocatable :: sndbuf(:)
    integer(ip), allocatable :: rcvbuf(:)

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: rcvhd

    ! Request handlers for non-blocking receives
    integer, allocatable, dimension(:) :: sndhd

    integer :: mpi_comm, iret, proc_to_comm
    integer :: p2pstat(mpi_status_size)
    integer, parameter :: root_pid = 0

    onl = nobjs
    oni = 1 
    onb = onl - oni
    one = 0 

    ipart = me + 1 

    num_rcv = 0
    num_snd = 0 

    cur_snd = 1
    cur_rcv = 1


    call memalloc ( npadj  , list_rcv, __FILE__,__LINE__ )
    call memalloc ( npadj+1, rcv_ptrs, __FILE__,__LINE__ )
    
    call memalloc ( npadj  , list_snd, __FILE__,__LINE__ )
    call memalloc ( npadj+1, snd_ptrs, __FILE__,__LINE__ )

    call memalloc ( int_objs_p(npadj+1)-1, lids_rcv, __FILE__,__LINE__ )
    call memalloc ( int_objs_p(npadj+1)-1, lids_snd, __FILE__,__LINE__ )
  
    ! Traverse ipart's neighbours in the graph of parts
    do iedge = 1, npadj
       neighbour_part_id = lpadj (iedge)
       
       visited_snd = 0 
       visited_rcv = 0 
     
       ! write(*,*) 'PPP', ipart, neighbour_part_id, int_objs_p(iedge), int_objs_p(iedge+1)

       ! Traverse list of objects on the current edge
       do iposobj=int_objs_p(iedge), int_objs_p(iedge+1)-1
          idobj = int_objs_l (iposobj)

          ! write (*,*) 'ZZZ', idobj, ipart, lobjs(5,idobj)
          ! write(*,*) 'XXX', ipart, neighbour_part_id, idobj ! DBG:
          if ( ipart == lobjs(5,idobj) ) then                       ! ipart has to send global object id to neighbour_part_id
             if ( visited_snd == 0 ) then
                ! No 
                num_snd = num_snd + 1 
                list_snd (num_snd) = neighbour_part_id              
                snd_ptrs (num_snd+1) = 1
                visited_snd = 1 
             else
                ! Yes
                snd_ptrs (num_snd+1) = snd_ptrs(num_snd+1) + 1
             end if
             lids_snd(cur_snd)=idobj
             cur_snd = cur_snd + 1
          else if ( neighbour_part_id == lobjs(5,idobj) ) then      ! ipart has to receive global object id from neighbour_part_id
             if ( visited_rcv == 0 ) then
                ! No 
                num_rcv = num_rcv + 1 
                list_rcv (num_rcv) = neighbour_part_id
                rcv_ptrs (num_rcv+1) = 1
                visited_rcv = 1 
             else
                ! Yes
                rcv_ptrs (num_rcv+1) = rcv_ptrs(num_rcv+1) + 1
             end if
             lids_rcv(cur_rcv)=idobj
             cur_rcv = cur_rcv + 1
             ! else ! Nothing to be sent from ipart/received from neighbour_part_id 
          end if
       end do   ! iposobj
    end do     ! iedge

    ! write(*,*) 'XXX', list_snd(1:num_snd)
    ! write(*,*) 'YYY', list_rcv(1:num_rcv)


    ! Transform rcv_ptrs from size to pointers
    rcv_ptrs(1) = 1
    do i = 1, num_rcv
       rcv_ptrs (i+1) = rcv_ptrs (i+1) + rcv_ptrs (i)  
    end do
  
    ! Transform snd_ptrs from size to pointers
    snd_ptrs(1) = 1
    do i = 1, num_snd
       snd_ptrs (i+1) = snd_ptrs(i+1) + snd_ptrs (i)  
    end do

!!$    write (*,*) 'num_rcv' , num_rcv
!!$    write (*,*) 'rcv_ptrs', rcv_ptrs(1:num_rcv+1)
!!$    write (*,*) 'lst_rcv' , list_rcv (1:num_rcv)
!!$    write (*,*) 'lids_rcv', lids_rcv(1:rcv_ptrs(num_rcv+1)-1)
!!$
!!$    write (*,*) 'num_snd' , num_snd
!!$    write (*,*) 'snd_ptrs', snd_ptrs(1:num_snd+1)
!!$    write (*,*) 'lst_snd' , list_snd (1:num_snd)
!!$    write (*,*) 'lids_snd', lids_snd(1:snd_ptrs(num_snd+1)-1)


    ! 1. Count/List how many local objects i have to identify a global id
    nobjs_with_gid = 0
    call memalloc (nobjs, lobjs_with_gid, __FILE__,__LINE__)

    do iobj = 1, nobjs
       if ( ipart == lobjs(5, iobj) ) then
          nobjs_with_gid = nobjs_with_gid + 1
          lobjs_with_gid(nobjs_with_gid) = iobj
       end if
    end do

    !! write (*,*) 'XX', lobjs_with_gid(1:nobjs_with_gid)

    ! 2. Gather + Scatter
    ! Get MPI communicator associated to icontxt (in
    ! the current implementation of our wrappers
    ! to the MPI library icontxt and mpi_comm are actually 
    ! the same)
    call psb_get_mpicomm (icontxt, mpi_comm)

    if ( me == root_pid ) then
       call memalloc( np+1, ptr_gids, __FILE__,__LINE__ )
       call mpi_gather( nobjs_with_gid  , 1, psb_mpi_integer,  &
                        ptr_gids(2:np+1), 1, psb_mpi_integer, root_pid, mpi_comm, iret)
    else
       call memalloc( 0, ptr_gids, __FILE__,__LINE__ )
       call mpi_gather( nobjs_with_gid  , 1, psb_mpi_integer,  &
                        ptr_gids        , 1, psb_mpi_integer, root_pid, mpi_comm, iret)
    end if

    if ( iret /= mpi_success ) then
       write (0,*) 'Error: mpi_gather returned != mpi_success'
       call psb_abort (icontxt)    
    end if
    
    if ( me == root_pid ) then
       ! Transform length to header
       ptr_gids(1)=1 
       do i=1, np
          ptr_gids(i+1) = ptr_gids(i) + ptr_gids(i+1) 
       end do
       ong = ptr_gids(np+1)-1 
     ! write (*,*) 'XXX', ptr_coarse_dofs(1:np+1) ! DBG:
    end if

    call psb_bcast ( icontxt, ong, root=root_pid)

    call mpi_scatter  ( ptr_gids , 1, psb_mpi_integer, &
                        start_gid, 1, psb_mpi_integer, &
                        root_pid , mpi_comm, iret )

    
    call memalloc( onl, ol2g, __FILE__,__LINE__ )


    ol2g = -1 

    ! 3. Assign global id to my "local" objects
    do i=1, nobjs_with_gid
       ol2g(lobjs_with_gid(i))=start_gid
       start_gid = start_gid + 1
    end do


  ! 4. Nearest neighbour comm
  call memalloc (num_rcv, rcvhd, __FILE__,__LINE__)
  call memalloc (num_snd, sndhd, __FILE__,__LINE__)

  call memalloc (snd_ptrs(num_snd+1)-snd_ptrs(1), sndbuf, __FILE__,__LINE__)
  call memalloc (rcv_ptrs(num_rcv+1)-rcv_ptrs(1), rcvbuf, __FILE__,__LINE__)

  ! Pack sndbuf
  do i=1, snd_ptrs(num_snd+1)-1
     sndbuf(i) = ol2g(lids_snd(i))
  end do

  ! First post all the non blocking receives   
  do i=1, num_rcv
     proc_to_comm = list_rcv(i)
     
     ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
     call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
     
     ! Message size to be received
     sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)
     
     call mpi_irecv(  rcvbuf(rcv_ptrs(i)), sizmsg,        &
             &        psb_mpi_integer, proc_to_comm, &
             &        psb_int_swap_tag, mpi_comm, rcvhd(i), iret)
             
     if ( iret /= mpi_success ) then
        write (0,*) 'Error: mpi_irecv returned != mpi_success'
        call psb_abort (icontxt)    
     end if
  end do

  ! Secondly post all non-blocking sends
  do i=1, num_snd
     proc_to_comm = list_snd(i)

     ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
     call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

     ! Message size to be sent
     sizmsg = snd_ptrs(i+1)-snd_ptrs(i)

     call mpi_isend(sndbuf(snd_ptrs(i)), sizmsg, &
          & psb_mpi_integer, proc_to_comm,    &
          & psb_int_swap_tag, mpi_comm, sndhd(i), iret)

     if ( iret /= mpi_success ) then
        write (0,*) 'Error: mpi_isend returned != mpi_success'
        call psb_abort (icontxt)    
     end if
  end do

  ! Wait on all non-blocking receives
  do i=1, num_rcv
     proc_to_comm = list_rcv(i)

     ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
     call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)

     ! Message size to be received
     sizmsg = rcv_ptrs(i+1)-rcv_ptrs(i)

     call mpi_wait(rcvhd(i), p2pstat, iret)
     
     if ( iret /= mpi_success ) then
        write (0,*) 'Error: mpi_wait returned != mpi_success'
        call psb_abort (icontxt)    
     end if
     
  end do

  ! Finally wait on all non-blocking sends
  do i=1, num_snd
     proc_to_comm = list_snd(i)
     
     ! Get MPI rank id associated to proc_to_comm - 1 in proc_to_comm
     call psb_get_rank (proc_to_comm, icontxt, proc_to_comm-1)
     
     ! Message size to be received
     sizmsg = snd_ptrs(i+1)-snd_ptrs(i)
     
     call mpi_wait(sndhd(i), p2pstat, iret)
     if ( iret /= mpi_success ) then
        write (0,*) 'Error: mpi_wait returned != mpi_success'
        call psb_abort (icontxt)    
     end if
  end do


  ! 5. Assign global id to "remote"  objects
  ! Pack sndbuf
  do i=1, rcv_ptrs(num_rcv+1)-1
     ol2g(lids_rcv(i)) = rcvbuf(i) 
  end do

  !write (*,*) 'AQUI JODEr AQUI'

  ! do i=1, npadj
  !   write (*,*) 'neig', lpadj(i)
  !   do j=int_objs_p(i),int_objs_p(i+1)-1
  !      write (*,*) 'PPP', int_objs_l(j), lobjs(2,int_objs_l(j)), ol2g(int_objs_l(j))
  !   end do
  ! end do
  
  call memfree( ptr_gids, __FILE__,__LINE__ )

  call memfree (rcvhd,__FILE__,__LINE__) 
  call memfree (sndhd,__FILE__,__LINE__)

  call memfree (sndbuf,__FILE__,__LINE__)
  call memfree (rcvbuf,__FILE__,__LINE__)

  call memfree ( list_rcv,__FILE__,__LINE__)
  call memfree ( rcv_ptrs,__FILE__,__LINE__)
    
  call memfree ( list_snd,__FILE__,__LINE__)
  call memfree ( snd_ptrs,__FILE__,__LINE__)

  call memfree ( lids_rcv,__FILE__,__LINE__)
  call memfree ( lids_snd,__FILE__,__LINE__)

  call memfree (lobjs_with_gid,__FILE__,__LINE__)

!!$  do i=1, npadj
!!$     write (*,*) 'neig', lpadj(i)
!!$     do j=int_objs_p(i),int_objs_p(i+1)-1
!!$        write (*,*) ol2g(int_objs_l(j))
!!$     end do
!!$  end do

  end subroutine create_omap


  end module graph_distribution_names
