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
  use fem_block_graph_names
  use sort_names
  implicit none
# include "debug.i90"
  private

  public :: create_dof_info, create_element_to_dof_and_ndofs, create_object2dof, &
       & create_dof_graph_block


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_dof_info ( dhand, trian, femsp, f_blk_graph, gtype ) ! graph
    implicit none
    ! Dummy arguments
    type(dof_handler)              , intent(in)    :: dhand
    type(fem_triangulation)        , intent(in)    :: trian 
    type(fem_space)                , intent(inout) :: femsp 
    type(fem_block_graph)          , intent(inout) :: f_blk_graph 
    integer(ip)          , optional, intent(in)    :: gtype(dhand%nblocks) 

    ! Locals
    integer(ip) :: iblock, jblock
    type(fem_graph), pointer :: f_graph


    call create_element_to_dof_and_ndofs( dhand, trian, femsp )

    call create_object2dof( dhand, trian, femsp )

    ! Create block graph
    call f_blk_graph%alloc(dhand%nblocks)

    ! To be called after the reordering of dofs
    do iblock = 1, dhand%nblocks
       do jblock = 1, dhand%nblocks
          f_graph => f_blk_graph%get_block(iblock,jblock)
          if ( iblock == jblock .and. present(gtype) ) then
             call create_dof_graph_block ( iblock, jblock, dhand, trian, femsp, f_graph, gtype(iblock) )
          else
             call create_dof_graph_block ( iblock, jblock, dhand, trian, femsp, f_graph )
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

       ! Part 1: Put DOFs on VEFs, taking into account that DOFs only belong to VEFs when we do not
       ! enforce continuity (continuity(ielem) /= 0). We go through all objects, elements around the
       ! object, variables of the element, and if for the value of continuity of this element no 
       ! dofs have already been added, we add them and touch this object for this continuity value.
       ! In FEMPAR, continuity is an elemental value. If it is different from 0, the nodes/DOFs 
       ! geometrically on the interface belong to the interface objects (VEFs). Next, we only
       ! enforce continuity for elements with same continuity value (mater below), in order to 
       ! allow for situations in which we want to have continuity in patches and discontinuity among 
       ! patches based on physical arguments (every patch would involve its own value of continuity).
       ! For hp-adaptivity, we could consider the value in continuity to be p (order) and as a result
       ! not to enforce continuity among elements with different order SINCE it would lead to ERROR
       ! to enforce continuity among elements of different order.
       do iobje = 1, trian%num_objects          
          touch = 0
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             iprob = femsp%lelem(jelem)%problem
             nvapb = dhand%prob_block(iblock,iprob)%nd1
             if ( jelem <= trian%num_elems ) then 
                do ivars = 1, nvapb
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                   if ( femsp%lelem(jelem)%continuity(g_var) /= 0 ) then
                      mater = femsp%lelem(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
                      


                      do obje_l = 1, trian%elems(jelem)%num_objects
                         if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                      end do

                      !obje_l = local_position( iobje, trian%elems(jelem)%objects, trian%elems(jelem)%num_objects)
                      if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then
                         if ( touch(mater,g_var,1) == 0 ) then                            
                            touch(mater,g_var,1) = jelem
                            touch(mater,g_var,2) = obje_l
                            call put_new_vefs_dofs_in_vef_of_element ( dhand, trian, femsp, g_var, jelem, l_var, &
                                 count, obje_l )
                         else
                            call put_existing_vefs_dofs_in_vef_of_element ( dhand, trian, femsp, touch, mater, g_var, iobje, &
                                 &                                          jelem, l_var, o2n, obje_l )
                         end if
                      end if
                   end if
                end do
             else
                do ivars = 1, nvapb
                   l_var = dhand%prob_block(iblock,iprob)%a(ivars)
                   g_var = dhand%problems(iprob)%p%l2g_var(l_var)
                   if ( femsp%lelem(jelem)%continuity(g_var) /= 0 ) then
                      mater = femsp%lelem(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
                      do obje_l = 1, trian%elems(jelem)%num_objects
                         if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                      end do
                      !if ( femsp%lelem(jelem)%bc_code(l_var,obje_l) == 0 ) then

                      !write (*,*) 'INSERT DOF IN GHOST ELEMENT ELEM2DOF'
                      !write (*,*) 'OF iobje FROM ielem VALUE ',iobje,touch(mater,g_var,1)
                      !if ( touch(mater,g_var,1) /= 0) 
                      !write(*,*) 'WITH ELEM2DOF(',touch(mater,g_var,1),')=', femsp%lelem(touch(mater,g_var,1))%elem2dof
                      !write (*,*) 'INTO ghost elem',jelem

                      if ( touch(mater,g_var,1) /= 0) then
                         call put_existing_vefs_dofs_in_vef_of_element ( dhand, trian, femsp, touch, mater, g_var, iobje, &
                              &                                          jelem, l_var, o2n, obje_l )
                      end if

                      !write(*,*) 'RESULTING ELEM2DOF(',jelem,')=', femsp%lelem(jelem)%elem2dof
                      !end if
                   end if
                end do
             end if
          end do
       end do

       ! Part 2: Put DOFs on nodes belonging to the volume object (element). For cG we only do that when 
       ! static condensation is not active. Static condensation is for all variables, elements, etc. BUT
       ! it cannot be used with dG. The following algorithm is ASSUMING that this is the case, and we are
       ! not using dG + static condensations. In any case, when creating the fem_space there is an 
       ! automatic check for satisfying that.
       ! No check about strong Dirichlet boundary conditions, because they are imposed weakly in dG, and
       ! never appear in interior nodes in cG.
       if ( ( .not. femsp%static_condensation )  ) then
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
                   count = count +1
                   femsp%lelem(ielem)%elem2dof(l_node,l_var) = count
                end do
             end do
          end do
       end if

       ! Part 3: Assign total number of dofs created to fem space object
       femsp%ndofs(iblock) = count
    end do

  end subroutine create_element_to_dof_and_ndofs

  !*********************************************************************************
  ! This subroutine takes the triangulation, the dof handler, and the triangulation 
  ! and fills the object2dof structure. The object2dof structure puts on top of VEFs
  ! the DOFs that are meant to be continuous between elements with the same continuity
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

    ! Part 1: Count DOFs on VEFs, using the notion of continuity described above (in elem2dof)
    do iblock = 1, dhand%nblocks  
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
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do

       femsp%object2dof(iblock)%p(1) = 1
       do iobje = 2, trian%num_objects+1
          femsp%object2dof(iblock)%p(iobje) = femsp%object2dof(iblock)%p(iobje) + femsp%object2dof(iblock)%p(iobje-1)
       end do

       call memalloc ( femsp%object2dof(iblock)%p(trian%num_objects+1)-1, 3, femsp%object2dof(iblock)%l, __FILE__, __LINE__ )

       ! Part 2: List DOFs on VEFs, using the notion of continuity described above (in elem2dof)
       ! We note that the object2dof(iblock)%l(:,X) is defined for X = 1,2,3
       ! object2dof(iblock)%l(:,1) : DOF LID
       ! object2dof(iblock)%l(:,2) : Variable GID associated to that DOF
       ! object2dof(iblock)%l(:,3) : Continuity value associated to that DOF (to enforce continuity)
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
       ! call print_list_2d(6,femsp%object2dof(iblock))

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

    ! COUNT PART
    call count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph )  
    call count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph )  
    call count_nnz_all_dofs_vs_all_dofs_by_face_integration( iblock, jblock, dhand, trian, femsp, dof_graph )  

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

    ! LIST PART

    call memalloc( dof_graph%nv+1, aux_ia, __FILE__,__LINE__ )
    aux_ia = dof_graph%ia

    call list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia ) 
    call list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia )
    call list_nnz_all_dofs_vs_all_dofs_by_face_integration( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia ) 

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

    ! call fem_graph_print( 6, dof_graph )

    call memfree (aux_ia,__FILE__,__LINE__)

  end subroutine create_dof_graph_block



  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)                :: dof_graph

    ! Local variables
    type(hash_table_ip_ip) :: visited
    integer(ip) :: idof, ielem, inode, iobje, iprob, istat, ivars 
    integer(ip) :: jdof, jelem, job_g, jobje, k_var, l_dof, l_mat
    integer(ip) :: l_node, l_var, m_dof, m_mat, m_var, nvapb, touch, ltype

    ltype = dof_graph%type

    do iobje = 1, trian%num_objects             
       if ( femsp%object2dof(iblock)%p(iobje+1)-femsp%object2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then
                do jobje = 1, trian%elems(jelem)%num_objects
                   job_g = trian%elems(jelem)%objects(jobje)
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
                   iprob = femsp%lelem(jelem)%problem
                   nvapb = dhand%prob_block(iblock,iprob)%nd1
                   do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                      l_dof = femsp%object2dof(iblock)%l(idof,1)
                      l_var = femsp%object2dof(iblock)%l(idof,2)
                      l_mat = femsp%object2dof(iblock)%l(idof,3)
                      do ivars = 1, nvapb
                         k_var = dhand%prob_block(iblock,iprob)%a(ivars)
                         m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                         m_mat = femsp%lelem(jelem)%continuity(m_var)
                         if ( dhand%dof_coupl(l_var, m_var) == 1 .and. l_mat == m_mat ) then                
                            if ( ltype == csr ) then
                               dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                                    & + femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje+1) &
                                    & - femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje)
                            else ! ltype == csr_symm
                               do inode = femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(k_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,k_var)
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
  end subroutine count_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)                :: dof_graph
    integer(ip), intent(inout)                :: aux_ia(:)

    ! Local variables
    type(hash_table_ip_ip) :: visited
    integer(ip) :: idof, ielem, inode, iobje, iprob, istat, ivars 
    integer(ip) :: jdof, jelem, job_g, jobje, k_var, l_dof, l_mat
    integer(ip) :: l_node, l_var, m_dof, m_mat, m_var, nvapb, touch, ltype, count, ic

    ltype = dof_graph%type

    count = 0
    do iobje = 1, trian%num_objects 
       if ( femsp%object2dof(iblock)%p(iobje+1)-femsp%object2dof(iblock)%p(iobje) > 0) then
          call visited%init(100) 
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             if ( jelem <= trian%num_elems ) then 
                do jobje = 1, trian%elems(jelem)%num_objects
                   job_g = trian%elems(jelem)%objects(jobje)
                   call visited%put(key=job_g, val=touch, stat=istat)
                   if ( istat == now_stored ) then  ! interface-interface
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
                               else ! ltype == csr_symm 
                                  if ( m_dof >= l_dof ) then
                                     ic = aux_ia(l_dof)
                                     dof_graph%ja(ic) = m_dof
                                     aux_ia(l_dof) = aux_ia(l_dof)+1
                                  end if
                               end if
                            end if
                         end do
                      end do
                   end if
                end do
                !end do
                if (.not.femsp%static_condensation) then  ! interface-interior
                   iprob = femsp%lelem(jelem)%problem
                   nvapb = dhand%prob_block(iblock,iprob)%nd1
                   do idof = femsp%object2dof(iblock)%p(iobje), femsp%object2dof(iblock)%p(iobje+1)-1
                      l_dof = femsp%object2dof(iblock)%l(idof,1)
                      l_var = femsp%object2dof(iblock)%l(idof,2)
                      l_mat = femsp%object2dof(iblock)%l(idof,3)
                      do ivars = 1, nvapb
                         k_var = dhand%prob_block(iblock,iprob)%a(ivars)
                         m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                         m_mat = femsp%lelem(jelem)%continuity(m_var)
                         if ( dhand%dof_coupl(l_var, m_var) == 1 .and. l_mat == m_mat ) then                
                            if ( ltype == csr ) then
                               do inode = femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(k_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,k_var)
                                  ic = aux_ia(l_dof)
                                  dof_graph%ja(ic) = m_dof
                                  aux_ia(l_dof) = aux_ia(l_dof)+1
                               end do
                            else ! ltype == csr_symm
                               do inode = femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje), &
                                    & femsp%lelem(jelem)%nodes_object(k_var)%p%p(jobje+1)-1
                                  l_node = femsp%lelem(jelem)%nodes_object(k_var)%p%l(inode)
                                  m_dof = femsp%lelem(jelem)%elem2dof(l_node,k_var)
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

  end subroutine list_nnz_dofs_vefs_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)                :: dof_graph

    ! Local variables
    integer(ip) :: g_var, ielem, inode, int_i, iobje, iprob, ivars, jdof, jnode, job_g
    integer(ip) :: jobje, jvars, k_var, l_dof, l_mat, l_node, l_var, ltype, m_dof, m_mat
    integer(ip) :: m_node, m_var, nvapb

    ltype = dof_graph%type

    ! As commented for elem2dof, static condensation is false for dG, by construction of the 
    ! fem space.
    if (.not.femsp%static_condensation) then
       do ielem  = 1, trian%num_elems
          iobje = trian%elems(ielem)%num_objects+1
          iprob = femsp%lelem(ielem)%problem
          nvapb = dhand%prob_block(iblock,iprob)%nd1 
          do ivars = 1, nvapb
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             ! Interior - interior 
             do jvars = 1, nvapb
                k_var = dhand%prob_block(iblock,iprob)%a(jvars)
                m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                if ( dhand%dof_coupl(g_var,m_var) == 1 ) then
                   do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                        & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                      l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                      l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                      if ( ltype == csr ) then
                         dof_graph%ia(l_dof+1) =  dof_graph%ia(l_dof+1) &
                              &  + femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje+1) &
                              & - femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje)
                      else ! ltype == csr_symm 
                         do jnode = femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(k_var)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,k_var)
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
                ! Interior - interface 
                do jobje = 1, trian%elems(ielem)%num_objects
                   job_g = trian%elems(ielem)%objects(jobje)
                   do jdof = femsp%object2dof(jblock)%p(job_g), femsp%object2dof(jblock)%p(job_g+1)-1
                      m_dof = femsp%object2dof(jblock)%l(jdof,1)
                      m_var = femsp%object2dof(jblock)%l(jdof,2)   
                      m_mat = femsp%object2dof(jblock)%l(jdof,3)                      
                      if ( dhand%dof_coupl(g_var,m_var) == 1 .and. l_mat == m_mat ) then
                         do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                            l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
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

  end subroutine count_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! List NNZ (number of nonzero entries) for DOFs on the interface (VEFs) of elements against
  ! both interior and interface nodes.
  !*********************************************************************************
  subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity ( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)              :: dof_graph
    integer(ip), intent(inout)                  :: aux_ia(:) 

    ! Local variables
    integer(ip) :: g_var, ielem, inode, iobje, iprob, ivars, jdof, jnode, job_g
    integer(ip) :: jobje, jvars, k_var, l_dof, l_mat, l_node, l_var, ltype, m_dof, m_mat
    integer(ip) :: m_node, m_var, nvapb, i, ic

    ltype = dof_graph%type

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
                k_var = dhand%prob_block(iblock,iprob)%a(jvars)
                m_var = dhand%problems(iprob)%p%l2g_var(k_var)
                if ( dhand%dof_coupl(g_var,m_var) == 1 ) then
                   do inode = femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje), &
                        & femsp%lelem(ielem)%nodes_object(l_var)%p%p(iobje+1)-1
                      l_node = femsp%lelem(ielem)%nodes_object(l_var)%p%l(inode)
                      l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                      if ( ltype == csr ) then
                         do jnode = femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(k_var)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,k_var)
                            i= aux_ia(l_dof)
                            dof_graph%ja(i) = m_dof
                            aux_ia(l_dof) = aux_ia(l_dof)+1
                         end do
                      else ! ltype == csr_symm 
                         do jnode = femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje), &
                              & femsp%lelem(ielem)%nodes_object(k_var)%p%p(iobje+1)-1
                            m_node = femsp%lelem(ielem)%nodes_object(k_var)%p%l(jnode)
                            m_dof = femsp%lelem(ielem)%elem2dof(m_node,k_var)
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


  end subroutine list_nnz_dofs_vol_vs_dofs_vefs_vol_by_continuity

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the element interior (for dG) being 
  ! coupled due to integration on faces. *** This part requires more ellaboration for cdG
  ! generalization***
  !*********************************************************************************
  ! Note: Here we take a face in which we want to integrate dG terms for pairs of unknowns
  ! and couple all DOFs a la DG. Note that we are not using AT ALL the continuity value,
  ! since the integration in faces is driven by the faces selected for integration, which
  ! are being built accordingly to what one wants to do (cG, dG, dG for jump of cont value, etc.).
  ! One could think that it could happen that two elements K1 and K2 with two different 
  ! values of continuity that share a face where we want to integrate could put more than
  ! once the coupling among two nodes. It can never happen AS SOON AS one never creates
  ! an integration face between two elements with same continuity value, which is the
  ! expected usage. 
  ! *** We could put an assert about it when creating the integration list.
  !*********************************************************************************
  subroutine count_nnz_all_dofs_vs_all_dofs_by_face_integration ( iblock, jblock, dhand, trian, femsp, dof_graph )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)              :: dof_graph

    ! Local variables
    integer(ip) :: count, g_var, i, ielem, iface, inode, iobje, iprob, ivars, jelem
    integer(ip) :: jnode, jprob, jvars, k_var, knode, l_dof, l_faci, l_facj, l_node
    integer(ip) :: l_var, ltype, m_dof, m_var, m_node, nnode, nvapbi, nvapbj

    ltype = dof_graph%type

    ! Loop over all interior faces (boundary faces do not include additional coupling)
    do iface = 1,femsp%num_interior_faces
       iobje = femsp%interior_faces(iface)%face_object
       assert ( trian%objects(iobje)%num_elems_around == 2 ) 
       do i=1,2
          ielem = trian%objects(iobje)%elems_around(i)
          jelem = trian%objects(iobje)%elems_around(3-i)
          l_faci = local_position(femsp%interior_faces(iface)%face_object,trian%elems(ielem)%objects, &
               & trian%elems(ielem)%num_objects )
          l_facj = local_position(femsp%interior_faces(iface)%face_object,trian%elems(jelem)%objects, &
               & trian%elems(jelem)%num_objects )
          iprob = femsp%lelem(ielem)%problem
          jprob = femsp%lelem(jelem)%problem
          nvapbi = dhand%prob_block(iblock,iprob)%nd1 
          nvapbj = dhand%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                k_var = dhand%prob_block(iblock,jprob)%a(jvars)
                m_var = dhand%problems(jprob)%p%l2g_var(k_var)
                if ( dhand%dof_coupl(g_var,m_var) == 1 ) then
                   if ( ltype == csr ) then
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      nnode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1) &
                           &  -femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj)
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                              & + nnode
                      end do
                      nnode = femsp%lelem(jelem)%f_inf(k_var)%p%nnode - nnode
                      assert ( nnode > 0)
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                           &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                         l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                         l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                         dof_graph%ia(l_dof+1) = dof_graph%ia(l_dof+1) &
                              & + nnode
                      end do
                   else ! ltype == csr_symm 
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)
                            if ( m_dof >= l_dof ) then
                               dof_graph%ia(l_dof+1) = &
                                    & dof_graph%ia(l_dof+1) + 1
                            end if
                         end do
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, femsp%lelem(jelem)%f_inf(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)

                            do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                                 &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                               l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                               if ( m_dof >= l_dof ) then
                                  !write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  !write (*,*) 'ielem,jelem,iobje,l_faci,l_facj',ielem,jelem,iobje,l_faci,l_facj
                                  !write (*,*) 'ielem,jelem'
                                  !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
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
    end do

  end subroutine count_nnz_all_dofs_vs_all_dofs_by_face_integration

  !*********************************************************************************
  ! Count NNZ (number of nonzero entries) for DOFs on the element interior (for dG) being 
  ! coupled due to integration on faces. *** This part requires more ellaboration for cdG
  ! generalization***
  !*********************************************************************************
  subroutine list_nnz_all_dofs_vs_all_dofs_by_face_integration ( iblock, jblock, dhand, trian, femsp, dof_graph, aux_ia )  
    implicit none
    ! Parameters
    integer(ip), intent(in)                     :: iblock, jblock
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(in)                 :: femsp 
    type(fem_graph), intent(inout)              :: dof_graph
    integer(ip), intent(inout)                  :: aux_ia(:) 

    ! Local variables
    integer(ip) :: count, g_var, i, ielem, iface, inode, iobje, iprob, ivars, jelem
    integer(ip) :: jnode, jprob, jvars, k_var, knode, l_dof, l_faci, l_facj, l_node
    integer(ip) :: l_var, ltype, m_dof, m_var, m_node, nnode, nvapbi, nvapbj, ic

    ltype = dof_graph%type

    ! Loop over all interior faces (boundary faces do not include additional coupling)
    do iface = 1, femsp%num_interior_faces
       iobje = femsp%interior_faces(iface)%face_object
       assert ( trian%objects(iobje)%num_elems_around == 2 ) 
       do i=1,2
          ielem = trian%objects(iobje)%elems_around(i)
          jelem = trian%objects(iobje)%elems_around(3-i)
          l_faci = local_position(femsp%interior_faces(iface)%face_object,trian%elems(ielem)%objects, &
               & trian%elems(ielem)%num_objects )
          l_facj = local_position(femsp%interior_faces(iface)%face_object,trian%elems(jelem)%objects, &
               & trian%elems(jelem)%num_objects )
          iprob = femsp%lelem(ielem)%problem
          jprob = femsp%lelem(jelem)%problem
          nvapbi = dhand%prob_block(iblock,iprob)%nd1 
          nvapbj = dhand%prob_block(iblock,jprob)%nd1 
          do ivars = 1, nvapbi
             l_var = dhand%prob_block(iblock,iprob)%a(ivars)
             g_var = dhand%problems(iprob)%p%l2g_var(l_var)
             do jvars = 1, nvapbj
                k_var = dhand%prob_block(iblock,jprob)%a(jvars)
                m_var = dhand%problems(jprob)%p%l2g_var(k_var)
                if ( dhand%dof_coupl(g_var,m_var) == 1 ) then
                   if ( ltype == csr ) then
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%nodes_object(k_var)%p%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)
                            ic = aux_ia(l_dof)
                            dof_graph%ja(ic) = m_dof
                            aux_ia(l_dof) = aux_ia(l_dof) + 1
                         end do
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, femsp%lelem(jelem)%f_inf(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)
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
                      ! Couple all DOFs in ielem with face DOFs in jelem and viceversa (i=1,2)
                      do inode = 1, femsp%lelem(ielem)%f_inf(l_var)%p%nnode
                         l_dof = femsp%lelem(ielem)%elem2dof(inode,l_var)
                         do jnode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj), &
                              & femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1
                            m_node = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)
                            if ( m_dof >= l_dof ) then
                               !write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'
                               !write (*,*) 'INSERTION DUE TO COUPLING BY (ielem)-FACE(jelem)'
                               !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                               !write (*,*) 'IN POSITION',aux_ia(l_dof)
                               !write (*,*) 'VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV'
                               ic = aux_ia(l_dof)
                               dof_graph%ja(ic) = m_dof
                               aux_ia(l_dof) = aux_ia(l_dof) + 1
                            end if
                         end do
                      end do
                      ! Couple all face DOFs in ielem with DOFs in jelem NOT in the face and viceversa (i=1,2)
                      count = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj)
                      knode = -1
                      if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                         knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                      end if
                      do jnode = 1, femsp%lelem(jelem)%f_inf(k_var)%p%nnode
                         if ( jnode == knode) then
                            count = count+1
                            if (count <= femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%p(l_facj+1)-1) then
                               knode = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(count)
                            end if
                         else
                            m_node = femsp%lelem(jelem)%f_inf(k_var)%p%ntxob%l(jnode)
                            m_dof = femsp%lelem(jelem)%elem2dof(m_node,k_var)
                            do inode = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci), &
                                 &     femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%p(l_faci+1)-1
                               l_node = femsp%lelem(ielem)%f_inf(l_var)%p%ntxob%l(inode)
                               l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
                               if ( m_dof >= l_dof ) then
                                  !write (*,*) 'KKKKKKKKKKKKKK'
                                  !write (*,*) 'INSERTION DUE TO COUPLING BY FACE(ielem)-INTERIOR(jelem)'
                                  !write (*,*) 'IN DOF',l_dof,'BEING COUPLED TO DOF',m_dof
                                  !write (*,*) 'IN POSITION',aux_ia(l_dof)
                                  !write (*,*) 'KKKKKKKKKKKKKK'
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

  end subroutine list_nnz_all_dofs_vs_all_dofs_by_face_integration

  subroutine put_new_vefs_dofs_in_vef_of_element ( dhand, trian, femsp, g_var, jelem, l_var, &
                                                   count, obje_l )
    implicit none
    ! Parameters
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(inout)              :: femsp 
    integer(ip), intent(inout)                  :: count
    integer(ip), intent(in)                     :: g_var, jelem, l_var, obje_l

    ! Local variables
    integer(ip) :: inode, l_node

    do inode = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l), &
         &     femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l+1)-1 
       l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%l(inode)
       count = count + 1
       !write (*,*) '****PUT DOF**** (elem,obj_l,obj_g,node,idof) ',jelem,obje_l,iobje,l_node,count
       femsp%lelem(jelem)%elem2dof(l_node,l_var) = count
    end do

  end subroutine put_new_vefs_dofs_in_vef_of_element

  subroutine put_existing_vefs_dofs_in_vef_of_element ( dhand, trian, femsp, touch, mater, g_var, iobje, jelem, l_var, &
                                                        o2n, obje_l )
    implicit none
    ! Parameters
    type(dof_handler), intent(in)               :: dhand
    type(fem_triangulation), intent(in)         :: trian 
    type(fem_space), intent(inout)              :: femsp
    integer(ip), intent(in)                     :: touch(:,:,:), mater, g_var, iobje, jelem, l_var, obje_l
    integer(ip), intent(out)                    :: o2n(:)

    ! Local variables
    integer(ip) :: elem_ext, obje_ext, prob_ext, l_var_ext
    integer(ip) :: nnode, order, inode, l_node, inode_ext, inode_l

    elem_ext = touch(mater,g_var,1)
    obje_ext = touch(mater,g_var,2)
    prob_ext = femsp%lelem(elem_ext)%problem
    l_var_ext = dhand%g2l_vars(g_var,prob_ext)
    !write (*,*) '****EXTRACT DOF**** (object)', iobje, ' FROM: (elem,obj_l) ',elem_ext,obje_ext, ' TO  : (elem,obj_l)', jelem,obje_l
    assert ( l_var_ext > 0 )
    nnode = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext+1) &
         &  -femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext) 
    if ( nnode > 0) then  
       order = femsp%lelem(elem_ext)%f_inf(l_var_ext)%p%order
       if ( trian%objects(iobje)%dimension == trian%num_dims .and. &
            & nnode ==  (order+1)**trian%num_dims ) then
          order = order    ! hdG case
       elseif ( nnode ==  (order-1)**trian%objects(iobje)%dimension ) then
          order = order -2 ! cG case
       else
          assert ( 0 == 1) ! SB.alert : Other situations possible when dG_continuity, cdG, hp-adaptivity ?
       end if
       call permute_nodes_object(                                                                 &
            & femsp%lelem(elem_ext)%f_inf(l_var_ext)%p,                                           &
            & femsp%lelem(jelem)%f_inf(l_var)%p,                                                  &
            & o2n,obje_ext,obje_l,                                                                &
            & trian%elems(elem_ext)%objects,                                                      &
            & trian%elems(jelem)%objects,                                                         &
            & trian%objects(iobje)%dimension,                                                     &
            & order )
       do inode = 1, femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext+1) - &
            femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext)
          l_node = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%p(obje_ext) + inode - 1
          inode_ext = femsp%lelem(elem_ext)%nodes_object(l_var_ext)%p%l(l_node )
          l_node = femsp%lelem(jelem)%nodes_object(l_var)%p%p(obje_l) + o2n(inode) - 1
          inode_l = femsp%lelem(jelem)%nodes_object(l_var)%p%l(l_node)
          femsp%lelem(jelem)%elem2dof(inode_l,l_var) = femsp%lelem(elem_ext)%elem2dof(inode_ext,l_var_ext)
       end do ! SB.alert : 1) face object for cG and hdG, where the face must have all their nodes
       !                   2) corner / edge only for cG
       !            * Never here for dG, continuity interface, hanging objects, etc.


    end if
  end subroutine put_existing_vefs_dofs_in_vef_of_element



    !*********************************************************************************
    ! Auxiliary function that returns the position of key in list
    !*********************************************************************************
    integer(ip) function local_position(key,list,size)
      implicit none
      integer(ip) :: key, size, list(size)

      do local_position = 1,size
         if ( list(local_position) == key) exit
      end do
      assert ( local_position < size + 1 )

    end function local_position



  end module create_global_dof_info_names
