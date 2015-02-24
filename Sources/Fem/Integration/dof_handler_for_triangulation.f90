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
module dof_handler_for_triangulation_class
  use types
  use memor
  use psb_sort_mod
  use fem_mesh_class
  use fem_partition_class
  use fem_graph_class
  use fem_conditions_class
  use fem_materials_class
  use fem_space_class
  use fem_space_types
  use maps_class
  use mesh_graph
  use fem_blocks_class
  use array_class
  use adaptivity

  implicit none
# include "debug.i90"
  private
  ! Data needed to go from geometrical info (import, partition, mesh) to 
  ! a dof info (dof_import, dof_partition, dof_mesh)

  type dof_handler_per_block

     integer(ip) :: ndofs                    ! Total number of dofs (nprob)

     integer(ip), allocatable ::     &
          ndofs,                     &       ! Total number of dofs (nprob)
          nvarsxprob(:),             &       ! Number of vars per problem vars_prob(nprob)
          i_varsxprob(:),            &       ! Physical vars associated to a problem (ia)
          j_varsxprob(:)                     ! (ja)

  end type dof_handler_per_block

  type dof_handler

     integer(ip) ::     &
          nblocks,                   &       ! Number of blocks
          nprob,                     &       ! Number of different problems
          nvars                              ! Total number of different physical variables

     integer(ip), allocatable ::     &
          phys_prob(:),              &       ! List of physical problems (nprob) 
          dof_coupl(:,:)                     ! Dof_coupling(nvar,nvar) for avoiding allocation & assembly of zero blocks

     type(dof_handler), allocatable :: dof_handler_block(:)

     type(fem_blocks)  ::           &
          blocks                             ! blocks ordering

  end type dof_handler

  interface triangulation_to_graph
     module procedure triangulation_to_graph
     module procedure triangulation_to_grah_adaptivity
  end interface triangulation_to_graph

  ! Types
  public :: dof_handler

  ! Functions
  public :: dof_handler_print, dof_handler_create, dof_handler_fill, dof_handler_free

  public :: triangulation_to_graph, element2dof_create, object2dof_create ! two different needed?

  !integer(ip), parameter       :: cG_cont   = 1
  !integer(ip), parameter       :: dG_cont   = 2
  !integer(ip), parameter       :: cdG_cont  = 3

  ! Parameters
  !public :: cG_cont, dG_cont, cdG_cont

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine element2dof_create
    ! XXX

    ! We assume already allocated
    ! Allocate ndofs
    !if(.not.allocated(dhand%ndofs)) then
    !   call memalloc(dhand%blocks%nb,dhand%ndofs,__FILE__,__LINE__)
    !end if

    do ielem = 1, mesh%nelem
       femsp%lelem(ielem)%elem2dof = 0
    end do

    do ibloc=1,dhand%nblocks
       idofs = 1
       do iobje = 1,trian%num_objects
          touch = 0
          nd = 0
          do lelem = 1,trian%objects(iobje)%num_elems_around
             gelem = trian%objects(iobje)%elems_around(lelem)
             iprob = femsp%lelem(gelem)prob
             kvars = 0
             do lvars = dhand%dof_handler_block(ibloc)%i_varsxprob(iprob), &
                  & dhand%dof_handler_block(ibloc)%i_varsxprob(iprob+1)-1
                kvars = kvars+1
                gv = dhand%dof_handler_block(ibloc)%j_varsxprob(lvars)
                ig = femsp%lelem(gelem)%iv(gv)
                lv = lvars - dhand%dof_handler_block(ibloc)%i_varsxprob(iprob) + 1
                do obje_loc = 1,femsp%lelem(gelem)%f_inf(ig)%p%nobje
                   if (trian%elems(gelem)%objects(obje_loc) == iobje ) then
                      exit
                   end if
                end do
                ! Determine the dimension of the object                
                if ( nd == 0 ) then
                   nd = trian%objects(iobje)%dimension
                end if
                obje_sta = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc)
                obje_end = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc+1)-1
                !
                if ( fixities%code(gv,iobje) /= 1 ) then
                   ! TODO: For 3D cdG_cont on the edges may stablish continuity between elems
                   ! that should not be paired and to make discontinuous some edges that should 
                   ! be continuous. In any case it should work and solve the problem properly.
                   ! SB: I do not understand this comment

                   if ( touch(gv,1) == 0 .or. dhand%continuity == dG_cont ) then
                      add_new_DOFs = .true.
                   elseif (nd == 1 .or. femsp%lelem(gelem)%f_inf(ig)%p%order==touch(gv,4)) then
                      add_new_DOFs = .false.
                   else
                      assert(dhand%continuity == cdG_cont )
                      add_new_DOFs = .true.
                   end if
                   !
                   if (add_new_DOFs) then
                      ! BE CAREFUL: Think about boundary conditions (...)
                      ! Create element_dofs(ielem,ivar,inode)
                      do jnode = obje_sta, obje_end
                         inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                         femsp%lelem(gelem)%p_nod(inode+1) = femsp%lelem(gelem)%p_nod(inode+1) + 1
                         ! put a dof into the node
                         femsp%lelem(gelem)%elem2dof(inode,gv) = idofs
                         idofs = idofs + 1
                      end do
                      if (touch(gv,1) == 0) then
                         touch(gv,1) = gelem
                         touch(gv,2) = obje_loc
                         touch(gv,3) = obje_sta
                         touch(gv,4) = femsp%lelem(gelem)%f_inf(ig)%p%order
                         touch(gv,5) = lv
                      end if
                      !
                   else
                      call permute_nodes_object(                                                    &
                           & femsp%lelem(touch(gv,1))%f_inf(touch(gv,5))%p,                         &
                           & femsp%lelem(gelem)%f_inf(ig)%p,o2n,touch(gv,2),obje_loc,               &
                           & trian%elems(touch(gv,1))%objects,                                      &
                           & trian%elems(gelem)%objects,nd-1,                                       &
                           & femsp%lelem(gelem)%f_inf(ig)%p%order-2)
                      i = 1
                      do jnode = obje_sta, obje_end
                         inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                         knode = femsp%lelem(touch(gv,1))%f_inf(touch(gv,5))%   &
                              & p%ndxob_j(touch(gv,3)+o2n(i)-1)
                         femsp%lelem(gelem)%elem2dof(inode,gv) = &
                              femsp%lelem(touch(gv,1))%elem2dof(knode,gv)
                         femsp%lelem(gelem)%p_nod(inode+1) = femsp%lelem(gelem)%p_nod(inode+1) + 1
                         i = i+1
                      end do
                   end if
                else
                   ! Initialize the boundary conditions array
                   if (.not. allocated(femsp%lelem(gelem)%bc_value)) then
                      nnode = size(femsp%lelem(gelem)%unkno,dim=1)
                      call memalloc(nnode,femsp%lelem(gelem)%bc_value,__FILE__,__LINE__)
                   end if
                   do jnode = obje_sta, obje_end
                      inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                      if (.not. allocated( femsp%lelem(gelem)%bc_value(inode)%a)) then
                         nvars = size(femsp%lelem(gelem)%unkno,dim=2)
                         call array_create(nvars, femsp%lelem(gelem)%bc_value(inode))
                      end if
                      femsp%lelem(gelem)%p_nod(inode+1) = femsp%lelem(gelem)%p_nod(inode+1) + 1
                   end do
                end if
             end do
          end do
       end do
       ! Interior dofs
       if (femsp%kfl_scond == scond_off) then
          do gelem = 1,trian%num_elems
             iprob = femsp%lelem(gelem)%prob
             kvars = 0
             do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                kvars = kvars+1
                gv = dhand%j_varsxprob(ibloc)%a(lvars)
                ig = femsp%lelem(gelem)%iv(gv)
                obje_loc = femsp%lelem(gelem)%f_inf(ig)%p%nobje_dim(trian%num_dims+1)
                obje_sta = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc)
                obje_end = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc+1)-1
                do jnode = obje_sta, obje_end
                   inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                   femsp%lelem(gelem)%elem2dof(inode,gv) = idofs 
                   femsp%lelem(gelem)%p_nod(inode+1) = femsp%lelem(gelem)%p_nod(inode+1) + 1
                   idofs = idofs + 1
                end do
             end do
          end do
       end if
       dhand%ndofs(ibloc) = idofs - 1       
    end do

    do ielem = 1,trian%num_elems
       femsp%lelem(ielem)%p_nod(1) = 1
       do i=2,size(femsp%lelem(ielem)%p_nod,1)
          femsp%lelem(ielem)%p_nod(i) = femsp%lelem(ielem)%p_nod(i-1) + femsp%lelem(ielem)%p_nod(i)
       end do
    end do

  end subroutine element2dof_create






       ! Interior dofs
       do ielem = 1,nelem
          do inode OWNED by ielem
             PUT new dofs on interior nodes
          end do
       end do

       ! Interface dofs
       do iobje = 1,nobje
          do ielem neighboring iobje
             do iobje in ielem
                if (iobje.not.touched) then
                   label NEW dofs on iobje (elem2dof,object2dof)
                   if (method.not.dG) then
                      TOUCH iobje
                   end if
                else
                   TAKE dofs and PUT on iobje (elem2dof,object2dof)
                end if
             end do
          end do
       end do
       ! In order for this algorithm to work for dG we must consider dofs on vertices, edges, and faces
       ! I would eliminate ob2dof and dof2ob. The first one is being used in par_dof_handler_partition 
       ! and the second one in par_partition_corner_detection

     end subroutine element2dof_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine object2dof_create
    ! I would merge element2dof and object2dof computations
  end subroutine object2dof_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine triangulation_to_graph( trian, dhand, femsp, dof_graph, gtype, ibloc, jbloc  )
    implicit none
    type(fem_triangulation)  , intent(in)  :: trian
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space), optional, intent(in)  :: femsp

    ! Not graph, not mesh
    ! ibloc and jbloc ?

    ! CG case
    ! DG case
    ! Adaptivity case

    ! Interior nodes
    do ielem = 1,nelem
       do inode IN interior (and interface for dg and hdg)
          couple with all nodes in the interior/interface of the elements
       end do
    end do
    
    ! Interface nodes via coupling continuity
    do iobje = 1,nobje
       do ielem contains iobje
          do ivar in iobje/ielem
             if (continuity) then
                do jobje in ielem
                   do jvar in iobje/ielem
                      if (jvar > ivar) then
                         couple ivar with jvar and viceversa 
                         touch jvar
                      end if
                   end do
                end do
             end if
          end do
       end do
    end do

    ! Interface nodes coupling via face integration
       if (dg) then
          do iface = 1,nface
             do other element in face
                couple with nodes on that face () BOTH DIRECTIONS
             end do
          end do
       end if
    end do

    ! 1. cG graph continuity in first and second loop
    ! 2. dG graph coupling in first and third loop
    ! 3. hdG graph coupling in first and third loop, taking into account that for hdG only nodes on faces,
    !    not defined for vertices/edges.
    ! 4. cdG for hanging nodes solved due to the fact that hanging faces are different objects
    ! 5. cdG coupling in hanging faces via third loop
    ! 6. cdG w/ material in continuity w/ material (touch(ndofs,nmats)
    
  end subroutine triangulation_to_graph

end module dof_handler_for_triangulation_class
