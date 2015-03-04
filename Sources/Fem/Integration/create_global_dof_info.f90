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
module create_global_dof_info_module
  use types
  use array_class
  use memor
  use fem_triangulation_class
  use fem_space_class
  use dof_handler_class
  use fem_space_types

  
  implicit none
# include "debug.i90"
  private

  public :: create_global_dof_info

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine create_global_dof_info ( dhand, trian, femsp ) ! graph
    implicit none
    ! Parameters
    type(dof_handler), intent(in)          :: dhand
    type(fem_triangulation), intent(in)    :: trian 
    type(fem_space), intent(inout)         :: femsp 

    ! Local variables
    integer(ip) :: iprob, l_var, iblock, count, iobje, ielem, jelem, nvapb, ivars, g_var, inter
    integer(ip) :: obje_l, inode, l_node, elem_ext, obje_ext, prob_ext, l_var_ext, inter_ext, inode_ext, inode_l
    integer(ip) :: mater, order, nnode
    integer(ip) :: g2l_vars(dhand%nvars_global,dhand%nprobs), touch(femsp%num_materials,dhand%nvars_global,2)

    type(array_ip1) :: prob_block(dhand%nprobs)
    integer(ip)     :: o2n(max_nnode)

    ! Global to local of variables per problem

    do iprob = 1, dhand%nprobs
       do l_var = 1, dhand%problems(iprob)%nvars
          g2l_vars( dhand%problems(iprob)%l2g_var(l_var), iprob) = l_var
       end do
    end do

    ! Create elem2dof and obj2dof

    do iblock = 1, dhand%nblocks  

       ! variables of a given problem for current block (prob_block)
       do iprob = 1, dhand%nprobs
          count = 0
          do l_var = 1, dhand%problems(iprob)%nvars
             if ( dhand%vars_block(dhand%problems(iprob)%l2g_var(l_var)) == iblock ) then
                count = count + 1 
             end if
          end do
          call array_create( count, prob_block(iprob))
          count = 0 
          do l_var = 1, dhand%problems(iprob)%nvars
             if ( dhand%vars_block(dhand%problems(iprob)%l2g_var(l_var)) == iblock ) then
                count = count + 1 
                prob_block(iprob)%a(count) = l_var!dhand%problems(iprob)%l2g_var(l_var)
             end if
          end do
       end do

       count = 0
       ! interface
       do iobje = 1, trian%num_objects
          do ielem = 1, trian%objects(iobje)%num_elems_around
             jelem = trian%objects(iobje)%elems_around(ielem)
             iprob = femsp%lelem(jelem)%prob
             nvapb = prob_block(iprob)%nd1
             do ivars = 1, nvapb
                !l_var = g2l(ivars,iprob)
                l_var = prob_block(iprob)%a(ivars)
                g_var = dhand%problems(iprob)%l2g_var(l_var)
                inter = femsp%lelem(jelem)%iv(l_var)
                mater = femsp%lelem(jelem)%material(inter) ! SB.alert : material can be used as p 
                do obje_l = 1, trian%elems(jelem)%num_objects
                   if ( trian%elems(jelem)%objects(obje_l) == iobje ) exit
                end do
                if ( touch(mater,g_var,1) == 0) then
                   touch(mater,g_var,1) = jelem
                   touch(mater,g_var,1) = obje_l
                   !do inode = 1, femsp%lelem(jelem)%nodes_object(inter,obje_l)%nd1
                   do inode = femsp%lelem(jelem)%nodes_object(inter)%p%p(obje_l), &
                        &     femsp%lelem(jelem)%nodes_object(inter)%p%p(obje_l+1)-1 
                      l_node = femsp%lelem(jelem)%nodes_object(inter)%p%l(inode)
                      !femsp%lelem(jelem)%nodes_object(inter,obje_l)%a(inode)
                      count = count + 1
                      femsp%lelem(jelem)%elem2dof(l_node,l_var) = count
                   end do
                else
                   elem_ext = touch(mater,g_var,1)
                   obje_ext = touch(mater,g_var,2)
                   prob_ext = femsp%lelem(elem_ext)%prob
                   l_var_ext = g2l_vars(g_var,prob_ext)
                   assert ( l_var_ext == 0 )
                   inter_ext = femsp%lelem(elem_ext)%iv(l_var_ext)

                   nnode = femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext+1) &
                        &  -femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext) 
                   order = femsp%lelem(elem_ext)%f_inf(inter_ext)%p%order
                   if ( nnode ==  (order-2)**2 ) then
                      order = order -2 ! cG case
                   else if ( nnode ==  order**2 ) then
                      !order = order    ! hdG case
                   else
                      assert ( 0 == 1) ! SB.alert : Other situations possible when dG_material, cdG, hp-adaptivity ?
                   end if

                   call permute_nodes_object(                                                                 &
                        & femsp%lelem(elem_ext)%f_inf(inter_ext)%p,                                           &
                        & femsp%lelem(jelem)%f_inf(inter)%p,                                                  &
                        & o2n,obje_ext,obje_l,                                                                &
                        & trian%elems(elem_ext)%objects,                                                      &
                        & trian%elems(jelem)%objects,                                                         &
                        & trian%objects(iobje)%dimension,                                                     &
                        & order ) 
                   do inode_ext = femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext), &
                        &         femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext)-1
                      inode_l = femsp%lelem(elem_ext)%nodes_object(inter_ext)%p%p(obje_ext) &
                           &    + femsp%lelem(jelem)%nodes_object(inter)%p%l(o2n(inode_ext))-1
                      !inode_l = femsp%lelem(jelem)%nodes_object(inter,obje_l)%a(o2n(inode_ext))
                      femsp%lelem(jelem)%elem2dof(inode_l,l_var) = femsp%lelem(elem_ext)%elem2dof(inode_ext,l_var_ext)
                   end do ! SB.alert : 1) face object for cG and hdG, where the face must have all their nodes
                   !                   2) corner / edge only for cG
                   !            * Never here for dG, material interface, hanging objects, etc.
                end if
             end do
          end do
       end do
       ! interior
       if (femsp%kfl_scond == 0) then 
          do ielem = 1, trian%num_elems
             iprob = femsp%lelem(ielem)%prob
             nvapb = prob_block(iprob)%nd1
             do ivars = 1, nvapb
                l_var = prob_block(iprob)%a(ivars)
                g_var = dhand%problems(iprob)%l2g_var(l_var)  
                iobje = trian%elems(ielem)%num_objects+1
                do inode = femsp%lelem(ielem)%nodes_object(inter)%p%p(iobje), &
                     &     femsp%lelem(ielem)%nodes_object(inter)%p%p(iobje+1)-1 
                   l_node = femsp%lelem(ielem)%nodes_object(inter)%p%l(inode)
                   !l_node = femsp%lelem(ielem)%nodes_object(inter,iobje)%a(inode)
                   femsp%lelem(ielem)%elem2dof(l_node,l_var) = count
                   count = count +1
                end do
             end do
          end do
       end if

       ! Create object to dof

       ! call memalloc ( dofh%nvars_global, touch, __FILE, __LINE__ ) 
       ! call memalloc ( trian%num_objects+1, object2dof%p, __FILE, __LINE__ )
       ! do iobje = 1, trian%num_objects
       !    do ielem = 1, trian%objects(iobje)%num_elems_around
       !       jelem = trian%objects(iobje)%elems_around(ielem)
       !       nvapb = prob_block(iprob)%nd1
       !       do ivars = 1, nvapb
       !          l_var = prob_block(iprob)%a(ivars)
       !          g_var = dhand%problems(iprob)%l2g_var(l_var)
       !          inter = femsp%lelem(jelem)%iv(l_var)
       !          if ( touch(g_var) == 0 ) then
       !             touch(g_var) = 1
       !             inter = femsp%lelem(jelem)%iv(l_var)
       !             do inode = 1, femsp%lelem(jelem)%nodes_object(inter,iobje)%nd1
       !                object2dof%p(iobje+1) = object2dof%p(iobje+1) + femsp%lelem(jelem)%nodes_object(inter,iobje)%nd1
       !             end do
       !          end if
       !       end do
       !    end do
       ! end do
       ! call memfree( touch, __FILE__, __LINE__ )
       ! !
       ! object2dof%p(1) = 1
       ! do iobje = 2, trian%num_objects+1
       !    object2dof%p(iobje) = object2dof%p(iobje) + object2dof%p(iobje-1)
       ! end do
       ! ! 
       ! touch = 0
       ! count = 0
       ! do iobje = 1, trian%num_objects
       !    do ielem = 1, trian%objects(iobje)%num_elems_around
       !       jelem = trian%objects(iobje)%elems_around(ielem)
       !       nvapb = prob_block(iprob)%nd1
       !       do ivars = 1, nvapb
       !          l_var = prob_block(iprob)%a(ivars)
       !          g_var = dhand%problems(iprob)%l2g_var(l_var)
       !          inter = femsp%lelem(jelem)%iv(l_var)
       !          if ( touch(g_var) == 0 ) then
       !             touch(g_var) = 1
       !             inter = femsp%lelem(jelem)%iv(l_var)
       !             do inode = 1, femsp%lelem(jelem)%nodes_object(inter,iobje)%nd1
       !                l_node = nodes_object(inter,iobje)%a(inode)
       !                count = count + 1
       !                object2dof%l(count,1) = femsp%lelem(jelem)%elem2dof(l_node,l_var)
       !                object2dof%l(count,2) = l2g_var(l_var)
       !             end do
       !          end if
       !       end do
       !    end do
       ! end do

    end do

    ! Create graph

    ! forall ( iblock = 1, dhand%nblocks, jblock = 1, dhand%nblocks )

    !    assert ( gtype == csr_symm .or. gtype == csr )

    !    ! Initialize
    !    dof_graph%type = gtype
    !    dof_graph%nv  = dhand%ndofs(ibloc) ! SB.alert : not stored there anymore
    !    dof_graph%nv2 = dhand%ndofs(jbloc)
    !    call memalloc( dof_graph%nv+1, dof_graph%ia, __FILE__,__LINE__ )
    !    dof_graph%ia = 0

    !    do iobje = 1, trian%num_objects
    !       if ( object2dof%p(iobje+1)-object2dof%p(iobje) > 0) then
    !          call visited%init(100) 
    !          do ielem = 1, trian%objects(iobje)%num_elems_around
    !             jelem = trian%objects(iobje)%elems_around
    !             do jobje = 1, trian%elems(ielem)%num_objects
    !                job_g = trian%elems(jelem)%objects(jobje)
    !                call visited%put(key=job_g, val=1, stat=istat)
    !                if ( istat == now_stored ) then
    !                   do idof = object2dof%p(iobje), object2dof%p(iobje+1)-1
    !                      l_dof = object2dof%l(idof,1)
    !                      l_var = object2dof%l(idof,2)
    !                      do jdof = object2dof%p(job_g), object2dof%p(job_g+1)-1
    !                         m_dof = object2dof%l(jdof,2)
    !                         m_var = object2dof%l(jdof,2)
    !                         if ( dofh%dof_coupl(l_var,m_var) == 1 ) then
    !                            if ( gtype == csr ) then
    !                               dof_graph(iblock,jblock)%ia(l_dof+1) = &
    !                                    & dof_graph(iblock,jblock)%ia(l_dof+1) +
    !                            else ! gtype == csr_symm 
    !                               if ( m_dof >= l_dof ) then
    !                                  dof_graph(iblock,jblock)%ia(dof_g+1) = &
    !                                       & dof_graph(iblock,jblock)%ia(dof_g+1) + 1
    !                               end if
    !                            end if
    !                         end if
    !                      end do
    !                   end do
    !                end if
    !             end do
    !             !end do
    !             call visited%free
    !             if (femsp%kfl_scond == 0) then
    !                jobje = jobje + 1
    !                iprob = femsp%lelem(jelem)%prob
    !                do idof = object2dof%p(iobje), object2dof%p(iobje+1)-1
    !                   l_dof = object2dof%l(idof,1)
    !                   l_var = object2dof%l(idof,2)
    !                   nvapb = prob_block(iprob)%nd1
    !                   do ivars = 1, nvapb
    !                      k_var = prob_block(iprob)%a(ivars)
    !                      m_var = dhand%problems(iprob)%l2g_var(k_var)
    !                      inter = femsp%lelem(jelem)%iv(l_var)
    !                      ! do ivars = 1, dhand%problems(femsp%lelem(jelem)%prob)%nvars
    !                      if ( dofh%dof_coupl(l_var, m_var) == 1 ) then                
    !                         if ( gtype == csr ) then
    !                            dof_graph(iblock,jblock)%ia(l_dof+1) =  dof_graph(iblock,jblock)%ia(l_dof+1) &
    !                                 &  + femsp%lelem(jelem)%nodes_object(inter,jobje)%nd1
    !                         else ! gtype == csr_symm 
    !                            do inode = 1, femsp%lelem(jelem)%nodes_object(inter,jobje)%nd1
    !                               ind_l = femsp%lelem(jelem)%nodes_object(inter,jobje)%a(inode)
    !                               m_dof = femsp%lelem(jelem)%elem2dof(ind_l,m_var)
    !                               if ( m_dof >= l_dof ) then
    !                                  dof_graph(iblock,jblock)%ia(l_dof+1) = &
    !                                       & dof_graph(iblock,jblock)%ia(l_dof+1) + 1
    !                               end if
    !                            end if
    !                         end do
    !                      end if
    !                   end do
    !                end do
    !             end if
    !          end do
    !       end if
    !    end do
    !    ! Interior nodes
    !    if (femsp%kfl_scond == 0) then
    !       do ielem  = 1, trian%num_elems
    !          jobje = trian%elems(ielem)%num_objects+1
    !          iprob = femsp%lelem(ielem)%prob
    !          nvapb = prob_block(iprob)%nd1
    !          do ivars = 1, nvapb
    !             !l_var = g2l(ivars,iprob)
    !             l_var = prob_block(iprob)%a(ivars)
    !             g_var = dhand%problems(iprob)%l2g_var(l_var)
    !             int_i = femsp%lelem(ielem)%iv(l_var)
    !             !do ivars = 1, dhand%problems(femsp%lelem(ielem)%prob)%nvars
    !             !int_i = iv(g2l_var(ivars))    
    !             do jvars = 1, nvapb
    !                m_var = prob_block(iprob)%a(jvars)
    !                g_var = dhand%problems(iprob)%l2g_var(m_var)
    !                if ( dofh%dof_coupl(l_var,m_var) == 1 ) then
    !                   int_j = femsp%lelem(ielem)%iv(m_var)
    !                   do inode = 1, nodes_object(int_i,jobje)%nd1
    !                      l_node = nodes_object(int_i,jobje)%a(inode)
    !                      l_dof = femsp%lelem(ielem)%elem2dof(l_node,l_var)
    !                      if ( gtype == csr ) then
    !                         dof_graph(iblock,jblock)%ia(l_dof+1) =  dof_graph(iblock,jblock)%ia(l_dof+1) &
    !                              &  + femsp%lelem(Ielem)%nodes_object(int_i,jobje)%nd1
    !                      else ! gtype == csr_symm 
    !                         do inode = 1, femsp%lelem(jelem)%nodes_object(int_j,jobje)%nd1
    !                            ind_l = femsp%lelem(jelem)%nodes_object(int_j,jobje)%a(inode)
    !                            m_dof = femsp%lelem(jelem)%elem2dof(ind_l,m_var)
    !                            if ( m_dof >= l_dof ) then
    !                               dof_graph(iblock,jblock)%ia(l_dof+1) = &
    !                                    & dof_graph(iblock,jblock)%ia(l_dof+1) + 1
    !                            end if
    !                         end if
    !                      end do
    !                   end if
    !                end do
    !             end do
    !          end do
    !       end do
    !    end if
    !    !
    !    dof_graph(iblock,jblock)%ia(1) = 1
    !    do idof = 2, dof_graph(iblock,jblock)%nd+1
    !       dof_graph(iblock,jblock)%ia(idof) = dof_graph(iblock,jblock)%ia(idof) + dof_graph(iblock,jblock)%ia(idof-1)
    !    end do
    !    ! 
    !    ! FALTA LIST... DESPUES DE DEPURAR dof_graph(iblock,jblock)%ia(l_dof+1)+1 ->
    ! end forall

    ! TO BE DONE IN THE NEAR FUTURE FOR DG THINGS (SB.alert)
    !    ! Interface nodes coupling via face integration
    !    if (dg) then
    !       do iface = 1,nface
    !          do other element in face
    !             couple with nodes on that face () BOTH DIRECTIONS
    !          end do
    !       end do
    !    end if
    ! end do

  end subroutine create_global_dof_info

end module create_global_dof_info_module
