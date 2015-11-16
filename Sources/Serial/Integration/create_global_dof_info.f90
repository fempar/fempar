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
  use types_names
  use memor_names
  use serial_triangulation_names
  use serial_fe_space_names
  use dof_descriptor_names
  use fe_space_types_names
  use hash_table_names
  use graph_names
  use sort_names
  implicit none
# include "debug.i90"
  private

  public :: create_dof_info, create_element_to_dof_and_ndofs, create_vef2dof
  
contains
  !*********************************************************************************
  ! This subroutine generates the following DOF-related data for serial runs.
  ! 1) It generates DOF numbering and fills the elem2dof arrays (see explanation of 
  !    the subroutine below)
  ! 2) It generates the vef2dof array ( see explanation of the subroutine below)
  ! NOTE: In order to understand the following subroutines, it is useful to know that 
  ! currently, the *par_fe_space* (parallel) has a pointer to a *fe_space*, which 
  ! has not info about the distributed environment (since it is a serial object),
  ! but it includes the *num_elems* local elements + *num_ghosts* ghost elements.
  ! The ghost part is filled in *par_fe_space_create*. Thus, we know if we have a
  ! local or ghost element below checking whether *ielem* <=  *num_elems* (local)
  ! or not. We are using this trick below in some cases, to be able to use the same
  ! subroutines for parallel runs, and fill properly ghost info. For serial runs, 
  ! the only thing is that *ielem* <=  *num_elems* always, so it never goes to the
  ! ghost element part.
  !*********************************************************************************
  subroutine create_dof_info ( fe_space )
    implicit none
    type(serial_fe_space_t), intent(inout) :: fe_space 
    call create_element_to_dof_and_ndofs( fe_space )
    call create_vef2dof( fe_space )
  end subroutine create_dof_info

  !*********************************************************************************
  ! This subroutine takes the triangulation and fills the element2dof structure at every 
  ! finite element, i.e., it labels all dofs related to local elements (not ghost), after a 
  ! count-list procedure, and puts the number of dofs in ndofs structure (per block).
  ! Note 1: The numbering is per every block independently, where the blocks are 
  ! defined at the dof_descriptor. A global dof numbering is not needed in the code, 
  ! when blocks are being used.
  !*********************************************************************************
  subroutine create_element_to_dof_and_ndofs( fe_space ) 
    implicit none
    ! Parameters
    type(serial_fe_space_t)     , intent(inout) :: fe_space 

    ! Local variables
    integer(ip) :: iprob, l_var, iblock, count, iobje, ielem, jelem, nvapb, ivars, g_var
    integer(ip) :: obje_l, inode, l_node, elem_ext, obje_ext, prob_ext, l_var_ext, inode_ext, inode_l
    integer(ip) :: mater, order, nnode
    integer(ip) :: touch(fe_space%num_continuity,fe_space%dof_descriptor%nvars_global,2)

    integer(ip)     :: o2n(max_nnode)

    call memalloc ( fe_space%dof_descriptor%nblocks, fe_space%ndofs, __FILE__, __LINE__ )

    do iblock = 1, fe_space%dof_descriptor%nblocks  
       count = 0

       ! Part 1: Put DOFs on VEFs, taking into account that DOFs only belong to VEFs when we do not
       ! enforce continuity (continuity(ielem) /= 0). We go through all objects, elements around the
       ! object, variables of the element, and if for the value of continuity of this element no 
       ! DOFs have already been added, we add them and touch this object for this continuity value.
       ! In FEMPAR, continuity is an elemental value. If it is different from 0, the nodes/DOFs 
       ! geometrically on the interface belong to the interface objects (VEFs). Next, we only
       ! enforce continuity for elements with same continuity value (mater below), in order to 
       ! allow for situations in which we want to have continuity in patches and discontinuity among 
       ! patches based on physical arguments (every patch would involve its own value of continuity).
       ! For hp-adaptivity, we could consider the value in continuity to be p (order) and as a result
       ! not to enforce continuity among elements with different order SINCE it would lead to ERROR
       ! to enforce continuity among elements of different order.
       do iobje = 1, fe_space%g_trian%num_vefs          
          touch = 0
          do ielem = 1, fe_space%g_trian%vefs(iobje)%num_elems_around
             jelem = fe_space%g_trian%vefs(iobje)%elems_around(ielem)
             iprob = fe_space%finite_elements(jelem)%problem
             nvapb = fe_space%dof_descriptor%prob_block(iblock,iprob)%nd1
             if ( jelem <= fe_space%g_trian%num_elems ) then ! Local elements
                do ivars = 1, nvapb
                   l_var = fe_space%dof_descriptor%prob_block(iblock,iprob)%a(ivars)
                   g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var)
                   if ( fe_space%finite_elements(jelem)%continuity(g_var) /= 0 ) then
                      mater = fe_space%finite_elements(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
                      do obje_l = 1, fe_space%g_trian%elems(jelem)%num_vefs
                         if ( fe_space%g_trian%elems(jelem)%vefs(obje_l) == iobje ) exit
                      end do
                      if ( fe_space%finite_elements(jelem)%bc_code(l_var,obje_l) == 0 ) then
                         if ( touch(mater,g_var,1) == 0 ) then                            
                            touch(mater,g_var,1) = jelem
                            touch(mater,g_var,2) = obje_l
                            call put_new_vefs_dofs_in_vef_of_element ( fe_space%dof_descriptor, fe_space%g_trian, fe_space, g_var, jelem, l_var, &
                                 count, obje_l )
                         else
                            call put_existing_vefs_dofs_in_vef_of_element ( fe_space%dof_descriptor, fe_space%g_trian, fe_space, touch, mater, g_var, iobje, &
                                 &                                          jelem, l_var, o2n, obje_l )
                         end if
                      end if
                   end if
                end do
             else ! Ghost elements
                do ivars = 1, nvapb
                   l_var = fe_space%dof_descriptor%prob_block(iblock,iprob)%a(ivars)
                   g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var)
                   if ( fe_space%finite_elements(jelem)%continuity(g_var) /= 0 ) then
                      mater = fe_space%finite_elements(jelem)%continuity(g_var) ! SB.alert : continuity can be used as p 
                      do obje_l = 1, fe_space%g_trian%elems(jelem)%num_vefs
                         if ( fe_space%g_trian%elems(jelem)%vefs(obje_l) == iobje ) exit
                      end do
                      if ( touch(mater,g_var,1) /= 0) then
                         call put_existing_vefs_dofs_in_vef_of_element ( fe_space%dof_descriptor, fe_space%g_trian, fe_space, touch, mater, g_var, iobje, &
                              &                                          jelem, l_var, o2n, obje_l )
                      end if
                   end if
                end do
             end if
          end do
       end do

       ! Part 2: Put DOFs on nodes belonging to the volume object (element). For cG we only do that when 
       ! static condensation is not active. Static condensation is for all variables, elements, etc. BUT
       ! it cannot be used with dG. The following algorithm is ASSUMING that this is the case, and we are
       ! not using dG + static condensations. In any case, when creating the fe_space there is an 
       ! automatic check for satisfying that.
       ! No check about strong Dirichlet boundary conditions, because they are imposed weakly in dG, and
       ! never appear in interior nodes in cG.
       if ( ( .not. fe_space%static_condensation )  ) then
          do ielem = 1, fe_space%g_trian%num_elems
             iprob = fe_space%finite_elements(ielem)%problem
             nvapb = fe_space%dof_descriptor%prob_block(iblock,iprob)%nd1
             do ivars = 1, nvapb
                l_var = fe_space%dof_descriptor%prob_block(iblock,iprob)%a(ivars)
                g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var) 
                iobje = fe_space%g_trian%elems(ielem)%num_vefs+1
                do inode = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje), &
                     &     fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%p(iobje+1)-1 
                   l_node = fe_space%finite_elements(ielem)%nodes_per_vef(l_var)%p%l(inode)
                   count = count +1
                   fe_space%finite_elements(ielem)%elem2dof(l_node,l_var) = count
                end do
             end do
          end do
       end if

       ! Part 3: Assign total number of dofs created to fem space object
       fe_space%ndofs(iblock) = count
    end do

  end subroutine create_element_to_dof_and_ndofs

  !*********************************************************************************
  ! This subroutine takes the finite element space and fills the vef2dof structure. 
  ! The vef2dof structure puts on top of VEFs the DOFs that are meant to be continuous 
  ! between elements with the same continuity label. As an example, when using dG only, 
  ! vef2dof is void. It is more an acceleration array than a really needed structure, 
  ! but it is convenient when creating the dof graph. 
  !*********************************************************************************
  subroutine create_vef2dof ( fe_space ) 
    implicit none
    ! Parameters
    type(serial_fe_space_t), intent(inout) :: fe_space 

    ! Local variables
    integer(ip) :: iprob, l_var, iblock, count, iobje, ielem, jelem, nvapb, ivars, g_var
    integer(ip) :: obje_l, inode, l_node, mater, istat
    integer(ip) :: touch(fe_space%dof_descriptor%nvars_global,fe_space%num_continuity)

    allocate( fe_space%vef2dof(fe_space%dof_descriptor%nblocks), stat = istat )
    check( istat == 0)

    ! Part 1: Count DOFs on VEFs, using the notion of continuity described above (in elem2dof)
    do iblock = 1, fe_space%dof_descriptor%nblocks  
       fe_space%vef2dof(iblock)%n1 = fe_space%g_trian%num_vefs
       fe_space%vef2dof(iblock)%n2 = 3
       call memalloc ( fe_space%g_trian%num_vefs+1, fe_space%vef2dof(iblock)%p, __FILE__, __LINE__, 0 )
       do iobje = 1, fe_space%g_trian%num_vefs
          touch = 0
          do ielem = 1, fe_space%g_trian%vefs(iobje)%num_elems_around
             jelem = fe_space%g_trian%vefs(iobje)%elems_around(ielem)
             if ( jelem <= fe_space%g_trian%num_elems ) then 
                iprob = fe_space%finite_elements(jelem)%problem
                nvapb = fe_space%dof_descriptor%prob_block(iblock,iprob)%nd1
                do ivars = 1, nvapb
                   l_var = fe_space%dof_descriptor%prob_block(iblock,iprob)%a(ivars)
                   g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var)
                   mater = fe_space%finite_elements(jelem)%continuity(g_var) 
                   if ( mater /= 0 ) then
                      if ( touch(g_var,mater) == 0 ) then
                         touch(g_var,mater) = 1
                         do obje_l = 1, fe_space%g_trian%elems(jelem)%num_vefs
                            if ( fe_space%g_trian%elems(jelem)%vefs(obje_l) == iobje ) exit
                         end do
                         if ( fe_space%finite_elements(jelem)%bc_code(l_var,obje_l) == 0 ) then
                            fe_space%vef2dof(iblock)%p(iobje+1) = fe_space%vef2dof(iblock)%p(iobje+1) + fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l+1) &
                                 & - fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l)
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do

       fe_space%vef2dof(iblock)%p(1) = 1
       do iobje = 2, fe_space%g_trian%num_vefs+1
          fe_space%vef2dof(iblock)%p(iobje) = fe_space%vef2dof(iblock)%p(iobje) + fe_space%vef2dof(iblock)%p(iobje-1)
       end do

       call memalloc ( fe_space%vef2dof(iblock)%p(fe_space%g_trian%num_vefs+1)-1, 3, fe_space%vef2dof(iblock)%l, __FILE__, __LINE__ )

       ! Part 2: List DOFs on VEFs, using the notion of continuity described above (in elem2dof)
       ! We note that the vef2dof(iblock)%l(:,X) is defined for X = 1,2,3
       ! vef2dof(iblock)%l(:,1) : DOF LID
       ! vef2dof(iblock)%l(:,2) : Variable GID associated to that DOF
       ! vef2dof(iblock)%l(:,3) : Continuity value associated to that DOF (to enforce continuity)
       count = 0
       do iobje = 1, fe_space%g_trian%num_vefs
          touch = 0
          do ielem = 1, fe_space%g_trian%vefs(iobje)%num_elems_around
             jelem = fe_space%g_trian%vefs(iobje)%elems_around(ielem)
             if ( jelem <= fe_space%g_trian%num_elems ) then 
                iprob = fe_space%finite_elements(jelem)%problem
                nvapb = fe_space%dof_descriptor%prob_block(iblock,iprob)%nd1
                do ivars = 1, nvapb
                   l_var = fe_space%dof_descriptor%prob_block(iblock,iprob)%a(ivars)
                   g_var = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var)
                   mater = fe_space%finite_elements(jelem)%continuity(g_var)
                   if ( mater /= 0) then
                      if ( touch(g_var,mater) == 0 ) then
                         touch(g_var,mater) = 1
                         do obje_l = 1, fe_space%g_trian%elems(jelem)%num_vefs
                            if ( fe_space%g_trian%elems(jelem)%vefs(obje_l) == iobje ) exit
                         end do
                         if ( fe_space%finite_elements(jelem)%bc_code(l_var,obje_l) == 0 ) then
                            do inode = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l), &
                                 &     fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l+1)-1 
                               l_node = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%l(inode)
                               count = count + 1
                               fe_space%vef2dof(iblock)%l(count,1) = fe_space%finite_elements(jelem)%elem2dof(l_node,l_var)
                               fe_space%vef2dof(iblock)%l(count,2) = fe_space%dof_descriptor%problems(iprob)%p%l2g_var(l_var)
                               fe_space%vef2dof(iblock)%l(count,3) = mater
                            end do
                         end if
                      end if
                   end if
                end do
             end if
          end do
       end do
    end do
  end subroutine create_vef2dof
  
    !*********************************************************************************
  ! Auxiliary function that generates new DOFs and put them in a particular VEF of a given element
  !*********************************************************************************
  subroutine put_new_vefs_dofs_in_vef_of_element ( dof_descriptor, trian, fe_space, g_var, jelem, l_var, &
       count, obje_l )
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(serial_triangulation_t), intent(in)         :: trian 
    type(serial_fe_space_t), intent(inout)              :: fe_space 
    integer(ip), intent(inout)                  :: count
    integer(ip), intent(in)                     :: g_var, jelem, l_var, obje_l

    ! Local variables
    integer(ip) :: inode, l_node

    do inode = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l), &
         &     fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l+1)-1 
       l_node = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%l(inode)
       count = count + 1
       !write (*,*) '****PUT DOF**** (elem,obj_l,obj_g,node,idof) ',jelem,obje_l,iobje,l_node,count
       fe_space%finite_elements(jelem)%elem2dof(l_node,l_var) = count
    end do

  end subroutine put_new_vefs_dofs_in_vef_of_element

  !*********************************************************************************
  ! Auxiliary function that puts existing DOFs in a particular VEF of a given element
  !*********************************************************************************
  subroutine put_existing_vefs_dofs_in_vef_of_element ( dof_descriptor, trian, fe_space, touch, mater, g_var, iobje, jelem, l_var, &
       o2n, obje_l )
    implicit none
    ! Parameters
    type(dof_descriptor_t), intent(in)               :: dof_descriptor
    type(serial_triangulation_t), intent(in)         :: trian 
    type(serial_fe_space_t), intent(inout)              :: fe_space
    integer(ip), intent(in)                     :: touch(:,:,:), mater, g_var, iobje, jelem, l_var, obje_l
    integer(ip), intent(out)                    :: o2n(:)

    ! Local variables
    integer(ip) :: elem_ext, obje_ext, prob_ext, l_var_ext
    integer(ip) :: nnode, order, inode, l_node, inode_ext, inode_l

    elem_ext = touch(mater,g_var,1)
    obje_ext = touch(mater,g_var,2)
    prob_ext = fe_space%finite_elements(elem_ext)%problem
    l_var_ext = dof_descriptor%g2l_vars(g_var,prob_ext)
    !write (*,*) '****EXTRACT DOF**** (object)', iobje, ' FROM: (elem,obj_l) ',elem_ext,obje_ext, ' TO  : (elem,obj_l)', jelem,obje_l
    assert ( l_var_ext > 0 )
    nnode = fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%p(obje_ext+1) &
         &  -fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%p(obje_ext) 
    if ( nnode > 0) then  
       order = fe_space%finite_elements(elem_ext)%reference_element_vars(l_var_ext)%p%order
       if ( trian%vefs(iobje)%dimension == trian%num_dims .and. &
            & nnode ==  (order+1)**trian%num_dims ) then
          order = order    ! hdG case
       elseif ( nnode ==  (order-1)**trian%vefs(iobje)%dimension ) then
          order = order -2 ! cG case
       else
          assert ( 0 == 1) ! SB.alert : Other situations possible when dG_continuity, cdG, hp-adaptivity ?
       end if
       call permute_nodes_per_vef(                                                                 &
            & fe_space%finite_elements(elem_ext)%reference_element_vars(l_var_ext)%p,                                           &
            & fe_space%finite_elements(jelem)%reference_element_vars(l_var)%p,                                                  &
            & o2n,obje_ext,obje_l,                                                                &
            & trian%elems(elem_ext)%vefs,                                                      &
            & trian%elems(jelem)%vefs,                                                         &
            & trian%vefs(iobje)%dimension,                                                     &
            & order )
       do inode = 1, fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%p(obje_ext+1) - &
            fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%p(obje_ext)
          l_node = fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%p(obje_ext) + inode - 1
          inode_ext = fe_space%finite_elements(elem_ext)%nodes_per_vef(l_var_ext)%p%l(l_node )
          l_node = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%p(obje_l) + o2n(inode) - 1
          inode_l = fe_space%finite_elements(jelem)%nodes_per_vef(l_var)%p%l(l_node)
          fe_space%finite_elements(jelem)%elem2dof(inode_l,l_var) = fe_space%finite_elements(elem_ext)%elem2dof(inode_ext,l_var_ext)
       end do ! SB.alert : 1) face object for cG and hdG, where the face must have all their nodes
       !                   2) corner / edge only for cG
       !            * Never here for dG, continuity interface, hanging objects, etc.


    end if
  end subroutine put_existing_vefs_dofs_in_vef_of_element

end module create_global_dof_info_names
