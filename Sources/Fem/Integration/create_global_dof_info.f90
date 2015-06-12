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
module create_global_dof_info_names
  use types
  use array_names
  use memor
  use fem_triangulation_names
  use fem_space_names
  use dof_handler_names
  use fem_space_types
  use hash_table_names
  use fem_graph_names
  use sort_names
  implicit none
# include "debug.i90"
  private

  public :: create_dof_info, create_element_to_dof_and_ndofs, create_object2dof, &
          & create_dof_graph_block


contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_dof_info ( dhand, trian, femsp, dof_graph, gtype ) ! graph
    implicit none
    ! Dummy arguments
    type(dof_handler), intent(in)             :: dhand
    type(fem_triangulation), intent(in)       :: trian 
    type(fem_space), intent(inout)            :: femsp 
    type(fem_graph), allocatable, intent(out) :: dof_graph(:,:) 
    integer(ip), optional, intent(in)         :: gtype(dhand%nblocks) 

    ! Locals
    integer(ip) :: iblock, jblock, istat


    call create_element_to_dof_and_ndofs( dhand, trian, femsp )

    call create_object2dof( dhand, trian, femsp )

    ! Create graph
    allocate( dof_graph(dhand%nblocks,dhand%nblocks), stat = istat )
    check( istat == 0)

    ! To be called after the reordering of dofs
    do iblock = 1, dhand%nblocks
       do jblock = 1, dhand%nblocks
          if ( iblock == jblock .and. present(gtype) ) then
             call create_dof_graph_block ( iblock, jblock, dhand, trian, femsp, dof_graph(iblock,jblock), gtype(iblock) )
          else
             call create_dof_graph_block ( iblock, jblock, dhand, trian, femsp, dof_graph(iblock,jblock) )
          end if
       end do
    end do

  end subroutine create_dof_info

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the triangulation 
  ! and fills the element2dof structure at every finite element, i.e., it labels all
  ! dofs related to local elements (not ghost), after a count-list procedure, and
  ! puts the number of dofs in ndofs structure (per block).
  ! Note 1: The numbering is per every block independently, where the blocks are 
  ! defined at the dof_handler. A global dof numbering is not needed in the code, 
  ! when blocks are being used.
  !*********************************************************************************
  subroutine create_element_to_dof_and_ndofs( dhand, trian, femsp ) 
    implicit none
    ! Parameters
    type(dof_handler), intent(in)             :: dhand
    type(fem_triangulation), intent(in)       :: trian 
    type(fem_space), intent(inout)            :: femsp 

    ! Local variables
    integer(ip) :: iprob, l_var, iblock, count, iobje, ielem, jelem, nvapb, ivars, g_var
    integer(ip) :: obje_l, inode, l_node, elem_ext, obje_ext, prob_ext, l_var_ext, inode_ext, inode_l
    integer(ip) :: mater, order, nnode
    integer(ip) :: touch(femsp%num_continuity,dhand%nvars_global,2)

    integer(ip)     :: o2n(max_nnode)

    call memalloc ( dhand%nblocks, femsp%ndofs, __FILE__, __LINE__ )

    do iblock = 1, dhand%nblocks  

       count = 0
       ! interface
       do iobje = 1, trian%num_objects          
          touch = 0
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then 
                iprob = femsp%lelem(jelem)%problem
                nvapb = dhand%prob_block(iblock,iprob)%nd1
                do ivars = 1, nvapb
                   !l_var = g2l(ivars,iprob)
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)

                   if ( femsp%lelem(jelem)%continuity(g_var) /= 0 ) then
                      mater = femsp%lelem(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
                      do obje_l = 1, trian%elems(jelem)%num_objects
                         if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                      end do
                      if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then 
                         if ( touch(mater,g_var,1) == 0) then
                            touch(mater,g_var,1) = jelem
                            touch(mater,g_var,2) = obje_l
                            !do inode = 1, femsp%lelem(jelem)%nodes_object(inter,obje_l)%nd1
                            do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l), &
                                 &     femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1)-1 
                               l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                               !femsp%lelem(jelem)%nodes_object(inter,obje_l)%a(inode)
                               count = count + 1

                               ! write (*,*) '****PUT DOF**** (elem,obj_l,obj_g,node,idof) ',jelem,obje_l,iobje,l_node,count

                               femsp%lelem(jelem)%elem2dof(l_node,l_var) = count
                            end do
                         else
                            elem_ext = touch(mater,g_var,1)
                            obje_ext = touch(mater,g_var,2)
                            prob_ext = femsp%lelem(elem_ext)%problem
                            l_var_ext = dhand%g2l_vars(g_var,prob_ext)

                            !write (*,*) '****EXTRACT DOF**** (object)', iobje, ' FROM: (elem,obj_l) ',elem_ext,obje_ext, ' TO  : (elem,obj_l)', jelem,obje_l

                            assert ( l_var_ext > 0 )

                            nnode = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext+1) &
                                 &  -femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext) 

                            !write (*,*) 'nnode',nnode
                            !write (*,*) 'trian%objects(iobje)%dimension ',trian%objects(iobje)%dimension
                            !write (*,*) '(order-1)**trian%objects(iobje)%dimension ',(order-1)**trian%objects(iobje)%dimension 


                            !if ( nnode ==  (order-1)**trian%objects(iobje)%dimension ) then
                            !write (*,*) 'nnode XXX',nnode ! cG case
                            !end if
                            if ( nnode > 0) then  
                               order = femsp%lelem(elem_ext)%f_inf(l_var_ext)%p%order
                               !write (*,*) 'order',order
                               if ( trian%objects(iobje)%dimension == trian%num_dims .and. &
                                    & nnode ==  (order+1)**trian%num_dims ) then
                                  order = order    ! hdG case
                               elseif ( nnode ==  (order-1)**trian%objects(iobje)%dimension ) then
                                  order = order -2 ! cG case
                               else
                                  assert ( 0 == 1) ! SB.alert : Other situations possible when dG_continuity, cdG, hp-adaptivity ?
                               end if

                               !write (*,*) 'order',order

                               call permute_nodes_object(                                                                 &
                                    & femsp%lelem(elem_ext)%f_inf(l_var_ext)%p,                                           &
                                    & femsp%lelem(jelem)%f_inf(l_var)%p,                                                  &
                                    & o2n,obje_ext,obje_l,                                                                &
                                    & trian%elems(elem_ext)%objects,                                                      &
                                    & trian%elems(jelem)%objects,                                                         &
                                    & trian%objects(iobje)%dimension,                                                     &
                                    & order )

                               !write (*,*) 'PERMUTATION ARRAY:',o2n, '*'
                               ! do inode_ext = femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext), &
                               !      &         femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext+1)-1
                               !    inode_l = femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext) &
                               !         &    + femsp%lelem(jelem)%nodes_object(interl_var)%p%l(o2n(inode_ext))-1
                               do inode = 1, femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext+1) - &
                                    femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext)
                                  l_node = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext) + inode - 1
                                  inode_ext = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%l(l_node )
                                  l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l) + o2n(inode) - 1
                                  inode_l = femsp%lelem(jelem)%nodes_object(l_var)%p%l(l_node)
                                  !inode_l = femsp%lelem(jelem)%nodes_object(inter,obje_l)%a(o2n(inode_ext))
                                  !write (*,*) '****COPY FROM NODE: ',inode_ext,' TO NODE: ',inode_l, 'DOF', &
                                  !     & femsp%lelem(elem_ext)%elem2dof(inode_ext,l_var_ext)
                                  femsp%lelem(jelem)%elem2dof(inode_l,l_var) = femsp%lelem(elem_ext)%elem2dof(inode_ext,l_var_ext)
                               end do ! SB.alert : 1) face object for cG and hdG, where the face must have all their nodes
                               !                   2) corner / edge only for cG
                               !            * Never here for dG, continuity interface, hanging objects, etc.
                            end if
                         end if
                      end if
                   ! else ! dG case
                   !    !write (*,*) 'XXXXXXXXXXX DG XXXXXXXXXX'
                   !    !write(*,*) 'l_var',l_var
                   !    !write(*,*) 'iobje',iobje
                   !    !write(*,*) 'jelem',jelem                      
                   !    do obje_l = 1, trian%elems(jelem)%num_objects
                   !       if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                   !    end do
                   !    !write(*,*) 'obje_l',obje_l
                   !    !write(*,*) 'femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0',femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0
                   !    if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then    
                   !       !write(*,*) 'femsp%lelem(jelem)%nodes_object(l_var)%p%p',femsp%lelem(jelem)%nodes_object(l_var)%p%p   
                   !       !write(*,*) 'femsp%lelem(jelem)%f_inf(l_var)%p%ndxob%p',femsp%lelem(jelem)%f_inf(l_var)%p%ndxob%p
                   !       do inode = femsp%lelem(jelem)%f_inf(l_var)%p%ndxob%p(obje_l), &
                   !            &     femsp%lelem(jelem)%f_inf(l_var)%p%ndxob%p(obje_l+1)-1
                   !       !do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l), &
                   !       !     &     femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1)-1 
                   !          !write(*,*) 'inode',inode
                   !          !l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                   !          l_node = femsp%lelem(jelem)%f_inf(l_var)%p%ndxob%l(inode)
                   !          !write(*,*) 'l_node',l_node
                   !          !femsp%lelem(jelem)%nodes_object(inter,obje_l)%a(inode)
                   !          count = count + 1
                   !          write (*,*) '****PUT DOF**** (elem,obj_l,obj_g,node,idof) ',jelem,obje_l,iobje,l_node,count
                   !          femsp%lelem(jelem)%elem2dof(l_node,l_var) = count
                   !       end do
                   !    end if
                   end if
                end do
             end if
          end do
       end do

       ! interior
       if ( .not. femsp%static_condensation ) then 
          do ielem = 1, trian%num_elems
             iprob = femsp%lelem(ielem)%problem
             nvapb = dhand%prob_block(iblock,iprob)%nd1
             do ivars = 1, nvapb
                l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                g_var = dhand%problems(iprob)%p%l2g_var(l_var)  
                iobje = trian%elems(ielem)%num_objects+1
                do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                     &     femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1 
                   l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                   !l_node = femsp%lelem(ielem)%nodes_object(inter,iobje)%a(inode)
                   count = count +1
                   femsp%lelem(ielem)%elem2dof(l_node,l_var) = count
                end do
             end do
          end do
       end if

       femsp%ndofs(iblock) = count

       !write (*,*) 'NDOFS:',femsp%ndofs(iblock)
    end do

  end subroutine create_element_to_dof_and_ndofs

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the triangulation 
  ! and fills the object2dof structure. The object2dof structure puts on top of vefs
  ! the dofs that are meant to be continuous between elements with the same continuity
  ! label. As an example, when using dG only, object2dof is void. It is more an 
  ! acceleration array than a really needed structure, but it is convenient when
  ! creating the dof graph. 
  !*********************************************************************************
  subroutine create_object2dof ( dhand, trian, femsp ) 
    implicit none
    ! Parameters
    type(dof_handler), intent(in)             :: dhand
    type(fem_triangulation), intent(in)       :: trian 
    type(fem_space), intent(inout)            :: femsp 

    ! Local variables
    integer(ip) :: iprob, l_var, iblock, count, iobje, ielem, jelem, nvapb, ivars, g_var
    integer(ip) :: obje_l, inode, l_node, mater, istat
    integer(ip) :: touch(dhand%nvars_global,femsp%num_continuity)

    allocate( femsp%object2dof(dhand%nblocks), stat = istat )
    check( istat == 0)

    do iblock = 1, dhand%nblocks  
       ! Create object to dof
       !call memalloc ( dhand%nvars_global, touch, __FILE__, __LINE__ ) 
       femsp%object2dof(iblock)%n1 = trian%num_objects
       femsp%object2dof(iblock)%n2 = 3
       
       call memalloc ( trian%num_objects+1, femsp%object2dof(iblock)%p, __FILE__, __LINE__, 0 )
       do iobje = 1, trian%num_objects
          touch = 0
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then 
                iprob = femsp%lelem(jelem)%problem
                nvapb = dhand%prob_block(iblock,iprob)%nd1
                !write(*,*) 'nvapb:',nvapb
                do ivars = 1, nvapb
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                   mater = femsp%lelem(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p
                   if ( mater /= 0 ) then
                      if ( touch(g_var,mater) == 0 ) then
                         touch(g_var,mater) = 1
                         do obje_l = 1, trian%elems(jelem)%num_objects
                            if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                         end do
                         if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then
                            !write (*,*) 'ADD TO OBJECT',iobje,' #DOFS',femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1) &
                            ! & - femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l)
                            femsp%object2dof(iblock)%p(iobje+1) = femsp%object2dof(iblock)%p(iobje+1) + femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1) &
                                 & - femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l)
                            !do inode = 1, femsp%lelem(jelem)%nodes_object(l_var)%p%(iobje+1) - &
                            !     & femsp%lelem(jelem)%nodes_object(l_var)%p%(iobje)
                            !femsp%object2dof(iblock)%p(iobje+1) = femsp%object2dof(iblock)%p(iobje+1) + femsp%lelem(jelem)%nodes_object(l_var)%p%(iobje)%nd1
                            !end do
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do
       !call memfree( touch, __FILE__, __LINE__ )
       !
       !write (*,*) 'femsp%object2dof(iblock)%p', femsp%object2dof(iblock)%p
       femsp%object2dof(iblock)%p(1) = 1
       do iobje = 2, trian%num_objects+1
          femsp%object2dof(iblock)%p(iobje) = femsp%object2dof(iblock)%p(iobje) + femsp%object2dof(iblock)%p(iobje-1)
       end do
       !write (*,*) 'femsp%object2dof(iblock)%p', femsp%object2dof(iblock)%p
       !write (*,*) 'femsp%object2dof(iblock)%p(trian%num_objects+1)-1',femsp%object2dof(iblock)%p(trian%num_objects+1)-1
       call memalloc ( femsp%object2dof(iblock)%p(trian%num_objects+1)-1, 3, femsp%object2dof(iblock)%l, __FILE__, __LINE__ )
       ! 
       count = 0
       do iobje = 1, trian%num_objects
          touch = 0
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then 
                iprob = femsp%lelem(jelem)%problem
                nvapb = dhand%prob_block(iblock,iprob)%nd1
                do ivars = 1, nvapb
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                   mater = femsp%lelem(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p
                   if ( mater /= 0) then
                      if ( touch(g_var,mater) == 0 ) then
                         touch(g_var,mater) = 1
                         do obje_l = 1, trian%elems(jelem)%num_objects
                            if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                         end do
                         if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then
                            do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l), &
                                 &     femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1)-1 

                               !1, femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1) - &
                               !                        & femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l)
                               l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                               count = count + 1
                               femsp%object2dof(iblock)%l(count,1) = femsp%lelem(jelem)%elem2dof(l_node,l_var)
                               femsp%object2dof(iblock)%l(count,2) = dhand%problems(iprob)%p%l2g_var(l_var)
                               femsp%object2dof(iblock)%l(count,3) = mater
                            end do
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do
       !write (*,*) 'femsp%object2dof(iblock)%p', femsp%object2dof(iblock)%p
       !write (*,*) 'femsp%object2dof(iblock)%l', femsp%object2dof(iblock)%l
       !call print_list_2d(6,femsp%object2dof(iblock))

    end do
  end subroutine create_object2dof

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the triangulation 
  ! and creates a dof_graph. The dof_graph includes both the coupling by continuity
  ! like in continuous Galerkin methods, and the coupling by face terms (of 
  ! discontinuous Galerkin type). The algorithm considers both the case with static 
  ! condensation and without it. In order to call this subroutine, we need to compute 
  ! first element2dof and object2dof arrays.
  !*********************************************************************************
  subroutine create_dof_graph_block( iblock, jblock, dhand, trian, femsp, dof_graph, gtype ) 
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(out)                :: dof_graph
    integer(ip), optional, intent(in)           :: gtype


    ! Local variables
    integer(ip) :: iprob, l_var, count, iobje, ielem, jelem, nvapb, ivars, g_var, inter, inode, l_node
    integer(ip) :: ltype, idof, jdof, int_i, int_j, istat, jnode, job_g, jobje, jvars, k_var , touch
    integer(ip) :: l_dof, m_dof, m_node, m_var, posi, posf, l_mat, m_mat, knode

    integer(ip) :: nvapbi, nvapbj, nnode, i, iface, jprob, l_faci, l_facj, ic


    integer(ip), allocatable :: aux_ia(:)
    type(hash_table_ip_ip) :: visited

    if ( iblock == jblock ) then
       if (present(gtype) ) then 
          dof_graph%type = gtype
       else
          dof_graph%type = csr
       end if
    else ! iblock /= jblock
       dof_graph%type = csr
    end if


    touch = 1
    ltype = dof_graph%type
    assert ( ltype == csr_symm .or. ltype == csr )

    ! Initialize
    dof_graph%type = ltype
    dof_graph%nv  = femsp%ndofs(iblock) ! SB.alert : not stored there anymore
    dof_graph%nv2 = femsp%ndofs(jblock)
    call memalloc( dof_graph%nv+1, dof_graph%ia, __FILE__,__LINE__ )
    dof_graph%ia = 0

    ! COUNT

    do iobje = 1, trian%num_objects             
       if ( femsp%object2dof(iblock)%p(iobje+1)-femsp%object2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then
                do jobje = 1, trian%elems(jelem)%num_objects
                   job_g = trian%elems(jelem)%objects(jobje)
                   !write (*,*) 'job_g',job_g
                   !call visited%put(key=job_g, val=1, stat=istat)
                   call visited%put(key=job_g, val=touch, stat=istat)
                   !write (*,*) 'istat',istat
                   if ( istat == now_stored ) then   ! interface-interface
                      do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                         l_dof = femsp%object2dof(iblock)%l(idof,1)
                         l_var = femsp%object2dof(iblock)%l(idof,2)
                         l_mat = femsp%object2dof(iblock)%l(idof,3)
                         do jdof = femsp%object2dof(jblock)%p(job_g), femsp%object2dof(jblock)%p(job_g+1)-1
                            m_dof = femsp%object2dof(jblock)%l(jdof,1)
                            m_var = femsp%object2dof(jblock)%l(jdof,2)
                            m_mat = femsp%object2dof(jblock)%l(jdof,3)
                            if ( dhand%dof_coupl(l_var,m_var) == 1 .and. l_mat == m_mat ) then
                               if ( ltype == csr ) then
                                  dof_graph%ia(l_dof+1) = &
                                       & dof_graph%ia(l_dof+1) + 1
                               else ! ltype == csr_symm 
                                  if ( m_dof >= l_dof ) then
                                     dof_graph%ia(l_dof+1) = &
                                          & dof_graph%ia(l_dof+1) + 1
                                  end if
                               end if
                            end if
                         end do
                      end do
                   end if
                end do
                !end do
                if (.not.femsp%static_condensation) then  ! interface-interior
                   ! jobje = jobje 
                   iprob = femsp%lelem(jelem)%problem
                   do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                      l_dof = femsp%object2dof(iblock)%l(idof,1)
                      l_var = femsp%object2dof(iblock)%l(idof,2)
                      l_mat = femsp%object2dof(iblock)%l(idof,3)
                      nvapb = dhand%prob_block(iblock,iprob)%nd1
                      do ivars = 1, nvapb
                         k_var = dhand%prob_block(iblock,iprob)%a(ivars)
                         m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                         ! do ivars = 1, dhand%problems(femsp%lelem(jelem)%problem)%nvars
                         if ( dhand%dof_coupl(l_var, m_var) == 1 .and. l_mat == m_mat ) then                
                            if ( ltype == csr ) then
                               dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                                    & + femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje+1) &
                                    & - femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje)
                            else ! ltype == csr_symm
                               do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,m_var)
                                  if ( m_dof >= l_dof ) then
                                     dof_graph%ia(l_dof+1) = &
                                          & dof_graph%ia(l_dof+1) + 1
                                  end if
                               end do
                            end if
                         end if
                      end do
                   end do
                end if
             end if
          end do
          call visited%free
       end if
    end do
    ! Interior nodes
    if (.not.femsp%static_condensation) then
       do ielem  = 1, trian%num_elems
          iobje = trian%elems(ielem)%num_objects+1
          iprob = femsp%lelem(ielem)%problem
          nvapb = dhand%prob_block(iblock,iprob)%nd1 
          !write(*,*) 'ielem,iobje,iprob,nvapb',ielem,iobje,iprob,nvapb 
          do ivars = 1, nvapb
             !l_var = g2l(ivars,iprob)
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             int_i = l_var
             !write(*,*) 'ivars,l_var,g_var',ivars,l_var,g_var
             ! Interior - interior (inside element)
             do jvars = 1, nvapb
                m_var = dhand%prob_block(iblock,iprob)%a(jvars)
                k_var = dhand%problems(iprob)%p%l2g_var(m_var)
                !write(*,*) 'jvars,m_var,k_var',jvars,m_var,k_var
                !write(*,*) ' dhand%dof_coupl(g_var,k_var) ', dhand%dof_coupl(g_var,k_var)
                if ( dhand%dof_coupl(g_var,k_var) == 1 ) then
                   !write (*,*) 'femsp%lelem(ielem)%nodes_object(int_i)%p%p',femsp%lelem(ielem)%nodes_object(int_i)%p%p
                   do inode = femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje), &
                        & femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje+1)-1
                      l_node = femsp%lelem(ielem)%nodes_object(int_i)%p%l(inode)
                      l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                      ! write (*,*) 'l_node,l_dof',l_node,l_dof
                      if ( ltype == csr ) then
                         dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                              &  + femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje+1) &
                              & - femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje)
                      else ! ltype == csr_symm 
                         do jnode = femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(int_i)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,m_var)
                            !write (*,*) 'm_node,m_dof',m_node,m_dof
                            if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             l_mat = femsp%lelem(ielem)%continuity(g_var)
             if ( l_mat /= 0 ) then
                ! Interior - interface (inside element)
                do jobje = 1, trian%elems(ielem)%num_objects
                   job_g = trian%elems(ielem)%objects(jobje)
                   do jdof = femsp%object2dof(jblock)%p(job_g), femsp%object2dof(jblock)%p(job_g+1)-1
                      m_dof = femsp%object2dof(jblock)%l(jdof,1)
                      m_var = femsp%object2dof(jblock)%l(jdof,2)   
                      m_mat = femsp%object2dof(jblock)%l(jdof,3)                      
                      if ( dhand%dof_coupl(g_var,m_var) == 1 .and. l_mat == m_mat ) then
                         do inode = femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(int_i)%p%p(iobje+1)-1
                            l_node = femsp%lelem(ielem)%nodes_object(int_i)%p%l(inode)
                            l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                            if ( ltype == csr ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1 
                            else if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end do
       end do
    end if

    !write(*,*) 'AFTER INTERIOR PART ONLY COUPLING INTRA ELEMENT'
    !write(*,*) '****START****'
    !write(*,*) 'dof_graph%ia',dof_graph%ia
    !write(*,*) '****END****'



    ! coupling by face integration (discontinuous Galerkin terms)
    ! count

    !write(*,*) 'femsp%num_interior_faces',femsp%num_interior_faces
    do iface = 1,femsp%num_interior_faces
       !write(*,*) 'iface',iface
       iobje = femsp%interior_faces(iface)%face_object
       !write (*,*) 'iobje',iobje
       assert ( trian%objects(iobje)%num_elems_around == 2 ) 
       !write (*,*) 'trian%objects(iobje)%elems_around',trian%objects(iobje)%elems_around
       do i=1,2
          ielem = trian%objects(iobje)%elems_around(i)
          jelem = trian%objects(iobje)%elems_around(3-i)
          l_faci = local_position(femsp%interior_faces(iface)%face_object,trian%elems(ielem)%objects, &
               & trian%elems(ielem)%num_objects )
          l_facj = local_position(femsp%interior_faces(iface)%face_object,trian%elems(jelem)%objects, &
               & trian%elems(jelem)%num_objects )
          !write(*,*) 'ielem,jelem:',ielem,jelem
          !write(*,*) 'ielem femsp%lelem(ielem)%elem2dof :',femsp%lelem(ielem)%elem2dof
          !write(*,*) 'ielem femsp%lelem(jelem)%elem2dof :',femsp%lelem(jelem)%elem2dof
          !write(*,*) 'lfaces:',l_faci,l_facj
          iprob = femsp%lelem(ielem)%problem
          jprob = femsp%lelem(jelem)%problem
          !write(*,*) 'probs:',iprob,jprob
          nvapbi = dhand%prob_block(iblock,iprob)%nd1 
          nvapbj = dhand%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                m_var = dhand%prob_block(iblock,jprob)%a(jvars)
                k_var = dhand%problems(jprob)%p%l2g_var(m_var)
                if ( dhand%dof_coupl(g_var,k_var) == 1 ) then
                   if ( ltype == csr ) then
                      nnode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1) &
                           &  -femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj)
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         !femsp%lelem(ielem)%nodes_object(l_var)%p%p(l_faci), &
                         !  & femsp%lelem(ielem)%nodes_object(l_var)%p%p(l_faci+1)-1
                         !l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                              & + nnode
                      end do
                      nnode = femsp%lelem(jelem)%f_inf(m_var)%p%nnode - nnode
                      assert ( nnode > 0)
                      do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                           &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                         l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                         l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                         dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                              & + nnode
                      end do
                   else ! ltype == csr_symm 
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)
                            if ( m_dof >= l_dof ) then
                               !write(*,*) 'INSERTION 1 DG'
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end do

                      !write(*,*) 'XXXXX'
                      write(*,*) 'femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p',femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p
                      write(*,*) 'femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l',femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l
                      count = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                      end if

                      do jnode = 1, femsp%lelem(jelem)%f_inf(m_var)%p%nnode
                         if ( jnode == knode) then
                            !write (*,*) 'jump node jnode',jnode
                            count = count+1
                            !write (*,*) 'count',count
                            if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                            end if

                            !write (*,*) 'knode',knode

                         else
                            m_node = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)

                            do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                                 &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                               l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)

                               !write(*,*) 'l_node,l_dof',l_node,l_dof

                               if ( m_dof >= l_dof ) then

                                  write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  !write (*,*) 'ielem,jelem,iobje,l_faci,l_facj',ielem,jelem,iobje,l_faci,l_facj
                                  !write (*,*) 'ielem,jelem'
                                  write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  !write (*,*) 'NOW',l_dof,'COUPLED TO',dof_graph%ia(l_dof+1) + 1
                                  !write (*,*) 'CHECK COLISION:count,knode,jnode',count,knode,jnode
                                  dof_graph%ia(l_dof+1) = &
                                       & dof_graph%ia(l_dof+1) + 1
                               end if
                            end do
                         end if
                      end do

                   end if
                end if
             end do
          end do
       end do
    end do!


    !
    dof_graph%ia(1) = 1
    do idof = 2, femsp%ndofs(iblock)+1
       dof_graph%ia(idof) = dof_graph%ia(idof) + dof_graph%ia(idof-1)
    end do

    ! write(*,*) 'DOF_GRAPH%IA'
    ! write(*,*) '****START****'
    ! write(*,*) 'dof_graph%ia',dof_graph%ia
    ! write(*,*) '****END****'

    call memalloc ( dof_graph%ia(femsp%ndofs(iblock)+1)-1, dof_graph%ja, __FILE__, __LINE__ )

    ! LIST
    call memalloc( dof_graph%nv+1, aux_ia, __FILE__,__LINE__ )
    aux_ia = dof_graph%ia

    !write(*,*) 'graph%ia : ',dof_graph%ia

    !write(*,*) '******* LIST *******'
    count = 0
    do iobje = 1, trian%num_objects 
       !write(*,*) 'LOOP OBJECTS **** IOBJE:',iobje    
       if ( femsp%object2dof(iblock)%p(iobje+1)-femsp%object2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             !write(*,*) '**** LOOP ELEMENTS CONTAIN OBJECT **** ELEM:',jelem
             if ( jelem <= trian%num_elems ) then 
                do jobje = 1, trian%elems(jelem)%num_objects
                   job_g = trian%elems(jelem)%objects(jobje)
                   !write (*,*) 'job_g',job_g
                   !call visited%put(key=job_g, val=1, stat=istat)
                   call visited%put(key=job_g, val=touch, stat=istat)
                   !write (*,*) 'istat',istat
                   if ( istat == now_stored ) then  ! interface-interface
                      !write(*,*) '******** LOOP  OBJECTS IN ELEM **** JOBJE:',job_g
                      do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                         l_dof = femsp%object2dof(iblock)%l(idof,1)
                         l_var = femsp%object2dof(iblock)%l(idof,2)
                         do jdof = femsp%object2dof(jblock)%p(job_g), femsp%object2dof(jblock)%p(job_g+1)-1
                            m_dof = femsp%object2dof(jblock)%l(jdof,1)
                            m_var = femsp%object2dof(jblock)%l(jdof,2)
                            if ( dhand%dof_coupl(l_var,m_var) == 1 ) then
                               if ( ltype == csr ) then
                                  !write(*,*) '************INSERT IN IDOF: ',l_dof,' JDOF: ',m_dof
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof)+1
                                  !& dof_graph%ia(l_dof+1) + 1
                               else ! ltype == csr_symm 
                                  if ( m_dof >= l_dof ) then
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                                     !& dof_graph%ia(l_dof+1) + 1
                                  end if
                               end if
                            end if
                         end do
                      end do
                   end if
                end do
                !end do
                if (.not.femsp%static_condensation) then  ! interface-interior
                   ! jobje = jobje 
                   iprob = femsp%lelem(jelem)%problem
                   do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                      l_dof = femsp%object2dof(iblock)%l(idof,1)
                      l_var = femsp%object2dof(iblock)%l(idof,2)
                      nvapb = dhand%prob_block(iblock,iprob)%nd1
                      do ivars = 1, nvapb
                         k_var = dhand%prob_block(iblock,iprob)%a(ivars)
                         m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                         ! do ivars = 1, dhand%problems(femsp%lelem(jelem)%problem)%nvars
                         if ( dhand%dof_coupl(l_var, m_var) == 1 ) then                
                            if ( ltype == csr ) then
                               do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,m_var)
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof)+1
                               end do
                            else ! ltype == csr_symm
                               do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(l_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,m_var)
                                  if ( m_dof >= l_dof ) then
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                                  end if
                               end do
                            end if
                         end if
                      end do
                   end do
                end if
             end if
          end do
          call visited%free
       end if
    end do


    !write(*,*) 'aux_ia',aux_ia

    ! Interior nodes
    if (.not.femsp%static_condensation) then   
       do ielem  = 1, trian%num_elems
          iobje = trian%elems(ielem)%num_objects+1
          iprob = femsp%lelem(ielem)%problem
          nvapb = dhand%prob_block(iblock,iprob)%nd1  
          do ivars = 1, nvapb
             !l_var = g2l(ivars,iprob)
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             ! Interior - interior (inside element)
             do jvars = 1, nvapb
                m_var = dhand%prob_block(iblock,iprob)%a(jvars)
                k_var = dhand%problems(iprob)%p%l2g_var(m_var)
                if ( dhand%dof_coupl(g_var,k_var) == 1 ) then
                   do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                        & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                      l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                      l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                      if ( ltype == csr ) then
                         do jnode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,m_var)
                            i= aux_ia(l_dof)
                            dof_graph%ja(i) = m_dof
                            aux_ia(l_dof) = aux_ia(l_dof)+1
                         end do
                      else ! ltype == csr_symm 
                         do jnode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,m_var)
                            if ( m_dof >= l_dof ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            end if
                         end do
                      end if
                   end do
                end if
             end do
             if ( femsp%lelem(ielem)%continuity(g_var) /= 0 ) then
                ! Interior - border (inside element)
                do jobje = 1, trian%elems(ielem)%num_objects
                   job_g = trian%elems(ielem)%objects(jobje)
                   do jdof = femsp%object2dof(jblock)%p(job_g), femsp%object2dof(jblock)%p(job_g+1)-1
                      m_dof = femsp%object2dof(jblock)%l(jdof,1)
                      m_var = femsp%object2dof(jblock)%l(jdof,2)                         
                      if ( dhand%dof_coupl(g_var,m_var) == 1 ) then
                         do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                            l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                            l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                            if ( ltype == csr ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            else if ( m_dof >= l_dof ) then
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof)+1
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end do
       end do
    end if

    !write(*,*) 'aux_ia',aux_ia




    ! coupling by face integration (discontinuous Galerkin terms)
    ! list
    !write(*,*) 'femsp%num_interior_faces',femsp%num_interior_faces
    do iface = 1, femsp%num_interior_faces
       !write(*,*) 'iface',iface
       iobje = femsp%interior_faces(iface)%face_object
       !write (*,*) 'iobje',iobje
       assert ( trian%objects(iobje)%num_elems_around == 2 ) 
       !write (*,*) 'trian%objects(iobje)%elems_around',trian%objects(iobje)%elems_around
       do i=1,2
          ielem = trian%objects(iobje)%elems_around(i)
          jelem = trian%objects(iobje)%elems_around(3-i)
          l_faci = local_position(femsp%interior_faces(iface)%face_object,trian%elems(ielem)%objects, &
               & trian%elems(ielem)%num_objects )
          l_facj = local_position(femsp%interior_faces(iface)%face_object,trian%elems(jelem)%objects, &
               & trian%elems(jelem)%num_objects )
          !write(*,*) 'ielem,jelem:',ielem,jelem
          !write(*,*) 'ielem femsp%lelem(ielem)%elem2dof :',femsp%lelem(ielem)%elem2dof
          !write(*,*) 'ielem femsp%lelem(jelem)%elem2dof :',femsp%lelem(jelem)%elem2dof
          !write(*,*) 'lfaces:',l_faci,l_facj
          iprob = femsp%lelem(ielem)%problem
          jprob = femsp%lelem(jelem)%problem
          !write(*,*) 'probs:',iprob,jprob
          nvapbi = dhand%prob_block(iblock,iprob)%nd1 
          nvapbj = dhand%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                m_var = dhand%prob_block(iblock,jprob)%a(jvars)
                k_var = dhand%problems(jprob)%p%l2g_var(m_var)
                if ( dhand%dof_coupl(g_var,k_var) == 1 ) then
                   if ( ltype == csr ) then

                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%nodes_object(m_var)%p%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)
                            ic = aux_ia(l_dof)
                            dof_graph%ja(ic) = m_dof
                            aux_ia(l_dof) = aux_ia(l_dof) + 1
                         end do
                      end do


                      count = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                      end if


                      do jnode = 1, femsp%lelem(jelem)%f_inf(m_var)%p%nnode
                         if ( jnode == knode) then
                            !write (*,*) 'jump node jnode',jnode
                            count = count+1
                            !write (*,*) 'count',count
                            if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                            end if


                         else
                            m_node = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)

                            do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                                 &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                               l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)

                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof) + 1

                            end do
                         end if
                      end do

                   else ! ltype == csr_symm 

                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)
                            if ( m_dof >= l_dof ) then



                               write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'
                               write (*,*) 'INSERTION DUE TO COUPLING BY (ielem)-FACE(jelem)'
                               write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                               write (*,*) 'IN POSITION',aux_ia(l_dof)
                               write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'

                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof) + 1
                            end if
                         end do
                      end do



                      count = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                      end if


                      do jnode = 1, femsp%lelem(jelem)%f_inf(m_var)%p%nnode
                         if ( jnode == knode) then
                            !write (*,*) 'jump node jnode',jnode
                            count = count+1
                            !write (*,*) 'count',count
                            if (count <= femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(count)
                            end if


                         else
                            m_node = femsp%lelem(jelem)%f_inf(m_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,m_var)

                            do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                                 &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                               l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)

                               if ( m_dof >= l_dof ) then

                                  write (*,*) 'KKKKKKKKKKKKKK'
                                  write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  write (*,*) 'IN POSITION',aux_ia(l_dof)
                                  write (*,*) 'KKKKKKKKKKKKKK'
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof) + 1
                               end if

                            end do
                         end if
                      end do




                   end if
                end if
             end do
          end do
       end do
    end do

    ! write(*,*) 'DOF_GRAPH%JA'
    ! write(*,*) '****START****'
    ! write(*,*) 'dof_graph%ja',dof_graph%ja
    ! write(*,*) '****END****'

    do idof = 1, femsp%ndofs(iblock)
       ! Order increasingly column identifiers of current row 
       ! using heap sort algorithm
       posi = dof_graph%ia(idof)
       posf = dof_graph%ia(idof+1)-1
       call sort(posf-posi+1,dof_graph%ja(posi:posf))
    end do
    do idof = 1, femsp%ndofs(iblock)
       !write (*,*) 'DOFS COUPLED TO IDOF:',idof
       !write (*,*) '****** START:'
       do l_dof = dof_graph%ia(idof),dof_graph%ia(idof+1)-1
          !write(*,'(I5,$)') dof_graph%ja(l_dof)
       end do
       !write (*,*) '****** END'
    end do
    !call fem_graph_print( 6, dof_graph )
    call memfree (aux_ia,__FILE__,__LINE__)


  end subroutine create_dof_graph_block
  
  integer(ip) function local_position(key,list,size)
    implicit none
    integer(ip) :: key, size, list(size)
    
    do local_position = 1,size
       if ( list(local_position) == key) exit
    end do
    assert ( local_position < size + 1 )
    
  end function local_position
  




   end module create_global_dof_info_names
