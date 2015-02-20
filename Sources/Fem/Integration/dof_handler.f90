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
module dof_handler_class
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

  type dof_handler
     integer(ip)        ::           &
          nprob,                     &       ! Number of different problems
          nvars,                     &       ! Total number of different physical variables 
          continuity                         ! Type of continuity (cG,dG,dG-material...)
     ! In a future, continuity can be an object-based quantity
     integer(ip), allocatable ::     &
          ndofs(:),                  &       ! Total number of dofs
          phys_prob(:),              &       ! List of physical problems (nprob)
          nvarsxprob(:),             &       ! Number of vars per problem vars_prob(nprob)
          dof_coupl(:,:)!,             &     ! Dof_coupling(nvar,nvar) for avoiding allocation & assembly of zero blocks
          !i_varsxprob(:,:),          &       ! Physical vars associated to a problem (ia)
          !j_varsxprob(:,:)!,         &       ! (ja)
          !ob2dof_l(:,:),             &       ! (nobje,2) : provides the dofs x object and associated variable
          !ob2dof_p(:)!,              &       ! pointer to ob2dof_l
          !dof2node(:)                        ! (ndofs)

     ! if object is constrained, a list of constraining dofs is stored
     ! if some of the constraining nodes is on essential bc, instead of the dof (which does not exist)
     ! a negative value of constraining node is stored. It has to be allways checked for negativity...
     type(array_ip2), allocatable :: &
          ob2dof_l(:)                         ! (nobje,2) : provides the dofs x object and associated variable

     type(array_ip1), allocatable :: &
          ob2dof_p(:),               &       ! pointer to ob2dof_l
          dof2ob_l(:),               &       ! provides the object associated at each dof
          blknvarsxprob(:),          &       ! Number of vars per problem vars_prob(nprob)
          i_varsxprob(:),            &       ! Physical vars associated to a problem (ia)
          j_varsxprob(:)!,           &       ! (ja)

     type(fem_blocks)  ::           &
          blocks                             ! blocks ordering

     !type(dofh_elem), allocatable :: &
     !     elem2dof(:)                       ! Map from elem to dof

     !logical(lg), allocatable     ::     &
     !     continuity(:)                      ! Continuity of dofs on a given node continuity (npoin)
     
     ! this pointer is at various places: here, in fem_space and in mesh
     ! it is mainly for convenience, since it is needed in multiple places in the code
     ! todo: maybe it should be rewised
     type(constraint_list), pointer     :: constraint_list => NULL()

  end type dof_handler

  interface graph_to_dof_graph
     module procedure  graph_to_dof_graph_block !graph_to_dof_graph_mono,
     module procedure graph_to_dof_graph_serial_mono, graph_to_dof_graph_serial_block
  end interface graph_to_dof_graph

  interface dG_dof_graph_create
     module procedure dG_dof_graph_mono_create, dG_dof_graph_block_create
  end interface dG_dof_graph_create

  ! Types
  public :: dof_handler

  ! Functions
  public ::  graph_to_dof_graph, dof_handler_print, dof_handler_create_one_physics, &
       &     dof_handler_free, elem2dof_create, dof_handler_fill, object2dof_create,               &
       &     dg_dof_graph_create, dof_handler_create_multiphysics

  integer(ip), parameter       :: cG_cont   = 1
  integer(ip), parameter       :: dG_cont   = 2
  integer(ip), parameter       :: cdG_cont  = 3

  ! Parameters
  public :: cG_cont, dG_cont, cdG_cont

contains

  subroutine dof_handler_fill( f_dof_handler, femsp, mesh, fixities )
    implicit none
    type(dof_handler),   intent(inout)            :: f_dof_handler
    type(fem_space),     intent(inout)            :: femsp
    type(fem_mesh),      intent(in)               :: mesh
    type(fem_conditions),intent(in)               :: fixities

    call elem2dof_create( f_dof_handler, femsp, mesh, fixities )

    call jvars_create( f_dof_handler, femsp, mesh ) 

  end subroutine dof_handler_fill   

  !=============================================================================
  subroutine elem2dof_create( dhand, femsp, mesh, fixities)
    implicit none
    type(dof_handler),   intent(inout)                    :: dhand
    type(fem_space),     intent(inout)                    :: femsp
    type(fem_mesh),      intent(in)                       :: mesh
    type(fem_conditions),intent(in)                       :: fixities

    ! Local variables
    type(fem_mesh)     :: dual_mesh
    integer(ip)        :: idofs,ielem,inode,jnode,knode,ipoin,inod0,ivars,ivar0,mvars,iobje,gv
    integer(ip)        :: lelem,gelem,obje_loc,ibloc
    integer(ip)        :: nnode, nvars, touch(dhand%nvars,4), kinte
    integer(ip)        :: obje_sta,obje_end, iprob, kvars, lvars, i, nd, ig
    integer(ip)        :: o2n(max_nnode)
    integer(ip)        :: topo_mesh_node_idx, clist_id
    logical            :: node_is_constrained, add_new_DOFs

    ! Create dual mesh
    call mesh_to_dual( mesh, dual_mesh ) ! IT COULD BE STORED IN MESH

    ! Allocate ndofs
    if(.not.allocated(dhand%ndofs)) then
       call memalloc(dhand%blocks%nb,dhand%ndofs,__FILE__,__LINE__)
    end if
    
    do ielem = 1, mesh%nelem
       femsp%lelem(ielem)%elem2dof = 0
    end do

    do ibloc=1,dhand%blocks%nb
       idofs = 1
       
       ! indices to constrains array
       ! todo: this should be done only once for all blocks, since it is more geometrical property
       ! todo: it may happen, however, that close to the boundary of 1 field inside another field, 
       ! it might be more complex
       ! todo: investigate further!
       do iobje = 1,mesh%npoin ! THIS WAY, FIRST INTERIOR, NEXT INTERFACE
          touch = 0
          nd = 0           
          do lelem = dual_mesh%pnods(iobje), dual_mesh%pnods(iobje+1)-1
             gelem = dual_mesh%lnods(lelem)
             iprob = femsp%lelem(gelem)%prob
             kvars = 0
             do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                kvars = kvars+1
                gv = dhand%j_varsxprob(ibloc)%a(lvars)
                ig = femsp%lelem(gelem)%iv(gv)
                do obje_loc = 1,femsp%lelem(gelem)%f_inf(ig)%p%nobje
                   if (mesh%lnods(mesh%pnods(gelem)+obje_loc-1) == iobje ) then
                      exit
                   end if
                end do
                ! Determine the dimension of the object
                if ( nd == 0 ) then
                   do nd = 1,mesh%ndime
                      if ( obje_loc < femsp%lelem(gelem)%f_inf(ig)%p%nobje_dim(nd+1) ) exit
                   end do
                end if
                obje_sta = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc)
                obje_end = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_i(obje_loc+1)-1
                if ( fixities%code(gv,iobje) /= 1 ) then
                   ! TODO: For 3D cdG_cont on the edges may stablish continuity between elems
                   ! that should not be paired and to make discontinuous some edges that should 
                   ! be continuous. In any case it should work and solve the problem properly.
                   if ( touch(gv,1) == 0 .or. dhand%continuity == dG_cont ) then
                      add_new_DOFs = .true.
                   elseif (nd == 1 .or. femsp%lelem(gelem)%f_inf(ig)%p%order==touch(gv,4)) then
                      add_new_DOFs = .false.
                   else
                      assert(dhand%continuity == cdG_cont )
                      add_new_DOFs = .true.
                   end if

                   if (add_new_DOFs) then
                      ! BE CAREFUL: Think about boundary conditions (...)
                      ! Create element_dofs(ielem,ivar,inode)
                      do jnode = obje_sta, obje_end
                         inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                         femsp%lelem(gelem)%p_nod(inode+1) = femsp%lelem(gelem)%p_nod(inode+1) + 1
                         
                         node_is_constrained = .false.
                         if ( fem_space_has_constraints(femsp) .and. (dhand%continuity == cG_cont)) then
                            topo_mesh_node_idx = mesh%lnods(mesh%pnods(gelem) + inode - 1) 
                            node_is_constrained = is_constrained_nodal(femsp%constraint_list,       &
                                 topo_mesh_node_idx )
                         end if
                         if(node_is_constrained) then
                            ! there is no dof in this node
                            ! put there negative value of index to the constrains table
                            clist_id = femsp%constraint_list%map_node2list(topo_mesh_node_idx)
                            femsp%lelem(gelem)%elem2dof(inode, gv) = - clist_id
                         else
                            ! put a dof into the node
                            femsp%lelem(gelem)%elem2dof(inode,gv) = idofs
                            idofs = idofs + 1
                         end if
                      end do
                      if (touch(gv,1) == 0) then
                         touch(gv,1) = gelem
                         touch(gv,2) = obje_loc
                         touch(gv,3) = obje_sta
                         touch(gv,4) = femsp%lelem(gelem)%f_inf(ig)%p%order
                      end if
                      !
                   else
                      call permute_nodes_object(                                                    &
                           & femsp%lelem(touch(gv,1))%f_inf(femsp%lelem(touch(gv,1))%iv(gv))%p,     &
                           & femsp%lelem(gelem)%f_inf(ig)%p,o2n,touch(gv,2),obje_loc,               &
                           & mesh%lnods(mesh%pnods(touch(gv,1)):mesh%pnods(touch(gv,1)+1)-1),       &
                           & mesh%lnods(mesh%pnods(gelem):mesh%pnods(gelem+1)-1),nd-1,              &
                           &  femsp%lelem(gelem)%f_inf(ig)%p%order-2)
                      i = 1
                      do jnode = obje_sta, obje_end
                         inode = femsp%lelem(gelem)%f_inf(ig)%p%ndxob_j(jnode)
                         knode = femsp%lelem(touch(gv,1))%f_inf(femsp%lelem(touch(gv,1))%iv(gv))%   &
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
       if (femsp%kfl_scond == scond_off) then
          do gelem = 1,femsp%o_mesh%nelem
             iprob = femsp%lelem(gelem)%prob
             kvars = 0
             do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                kvars = kvars+1
                gv = dhand%j_varsxprob(ibloc)%a(lvars)
                ig = femsp%lelem(gelem)%iv(gv)
                obje_loc = femsp%lelem(gelem)%f_inf(ig)%p%nobje_dim(mesh%ndime+1)
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

    do ielem = 1, femsp%o_mesh%nelem
!!$       kinte = 0
!!$       do i=1,femsp%lelem(ielem)%nint
!!$          if (femsp%lelem(ielem)%f_inf(femsp%lelem(ielem)%iv(i))%p%order > 2 ) kinte = kinte + 1
!!$       end do
!!$       if(kinte > 0) then
!!$          call permute_elem_dof_to_global(femsp%lelem(ielem),femsp%o_mesh%lnods                     &
!!$               &                         (femsp%o_mesh%pnods(ielem):femsp%o_mesh%pnods(ielem+1)-1), &
!!$               &                          dhand%nvars)
!!$       end if
       femsp%lelem(ielem)%p_nod(1) = 1
       do i=2,size(femsp%lelem(ielem)%p_nod,1)
          femsp%lelem(ielem)%p_nod(i) = femsp%lelem(ielem)%p_nod(i-1) + femsp%lelem(ielem)%p_nod(i)
       end do
    end do

    ! ????????????????????????????
    ! if (present(f_part) .and. present(new_f_part)) then
    !    if( iobje == f_part%nmap%ni ) then
    !       new_f_part%nmap%ni = idofs - 1
    !    end if
    ! end if
    ! idofs = idofs - 1

    ! Here, I would put dof2obje(i,j), where i is the object
    ! and j is the node in the object
    ! TO DO
    ! if (present(f_part) .and. present(new_f_part)) then
    !        new_f_part%nmap%nl = idofs 
    !        new_f_part%nmap%nb = idofs - new_f_part%nmap%ni
    !        new_f_part%nmap%ne = 0
    !     end if
    !     dhand%ndofs = idofs

    !     call memalloc (idofs, dhand%dof2node, __FILE__,__LINE__ )
    !     ! Store the amount of DOF
    !     dhand%ndofs = idofs


    call fem_mesh_free ( dual_mesh )

  end subroutine elem2dof_create

 !=============================================================================
  subroutine jvars_create( dhand, femsp, mesh )
    implicit none
    type(dof_handler),   intent(in)    :: dhand
    type(fem_space),     intent(inout) :: femsp
    type(fem_mesh),      intent(in)    :: mesh

    ! Locals
    integer(ip) :: ibloc,ielem,inode,iint,maxnnode,ivars,iprob,gvars,jvars

    do ielem=1,mesh%nelem
       iprob = femsp%lelem(ielem)%prob
       maxnnode=0
       do iint=1,femsp%lelem(ielem)%nint
          maxnnode=max(maxnnode,femsp%lelem(ielem)%f_inf(iint)%p%nnode)
       end do
       do inode=1,maxnnode
          jvars = 0
          do ibloc=1,dhand%blocks%nb
             do ivars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                gvars = dhand%j_varsxprob(ibloc)%a(ivars)
                iint  = femsp%lelem(ielem)%iv(gvars)
                if(femsp%lelem(ielem)%f_inf(iint)%p%nnode>=inode) then
                   jvars = jvars + 1
                   femsp%lelem(ielem)%jvars(inode,gvars) = jvars
                end if
             end do
          end do
       end do
    end do
    
  end subroutine jvars_create

  !=============================================================================
  subroutine graph_to_dof_graph_serial_mono( gtype, dhand, mesh, dof_graph,femsp)
    implicit none
    integer(ip)              , intent(in)  :: gtype
    type(fem_mesh)           , intent(in)  :: mesh
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space), optional, intent(in)  :: femsp
    ! Locals
    type(fem_graph) :: graph

    assert ( gtype == csr_symm .or. gtype == csr ) 

    call mesh_to_graph( scal, csr, 1, 1, mesh, graph)
    call graph_to_dof_graph_block(gtype,1,1,dhand,mesh,graph,dof_graph,femsp)
    call fem_graph_free( graph )

  end subroutine graph_to_dof_graph_serial_mono

  !=============================================================================
  subroutine graph_to_dof_graph_serial_block( gtype, ibloc, jbloc, dhand, mesh, dof_graph,femsp)
    implicit none
    integer(ip)              , intent(in)  :: gtype
    integer(ip)              , intent(in)  :: ibloc,jbloc
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_mesh)           , intent(in)  :: mesh
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space), optional, intent(in)  :: femsp
    ! Locals
    type(fem_graph) :: graph

    assert ( gtype == csr_symm .or. gtype == csr ) 

    call mesh_to_graph( scal, csr, 1, 1, mesh, graph)
    call graph_to_dof_graph_block(gtype,ibloc,jbloc,dhand,mesh,graph,dof_graph,femsp)
    call fem_graph_free( graph )

  end subroutine graph_to_dof_graph_serial_block

  !=============================================================================
  subroutine graph_to_dof_graph_block(gtype,ibloc,jbloc,dhand,mesh,graph,dof_graph,femsp)
    implicit none
    integer(ip)              , intent(in)  :: gtype
    integer(ip)              , intent(in)  :: ibloc,jbloc
    type(fem_mesh)           , intent(in)  :: mesh
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_graph)          , intent(in)  :: graph
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space), optional, intent(in)  :: femsp

    ! Local variables
    integer(ip)                    :: dof_nelty, dof_nnode, dof_gtype, iobje, idx
    integer(ip)                    :: ivars, jvars, idofs, jdofs, count !ipoin, lpoin, jpoin,
    integer(ip)                    :: ndofi, dof_g, lobje, jobje, tobje, gvars, lvars, mtype,ntype
    integer(ip)                    :: ielem, jelem, gelem, iprob, jprob, inode, jnode, gnode
    type(fem_mesh)                 :: dmesh
    integer(ip)                    :: var_idof, var_jdof

    assert ( gtype == csr_symm .or. gtype == csr )

    ! Initialize
    dof_graph%type = gtype
    dof_graph%nv  = dhand%ndofs(ibloc)
    dof_graph%nv2 = dhand%ndofs(jbloc)
    call memalloc( dof_graph%nv+1, dof_graph%ia, __FILE__,__LINE__ )
    dof_graph%ia = 0

    if ( present(femsp) ) then
       call mesh_to_dual(mesh,dmesh)
    end if

    !write(*,*) "dual mesh with constraints "
    !do idx = 1, dmesh%nelem
    !   write(*,*) "node ", idx, " -> elems : ", dmesh%lnods(dmesh%pnods(idx):dmesh%pnods(idx+1) - 1)
    !end do

    ! Count iDOFs on objects
    do iobje = 1,mesh%npoin
       do idofs = dhand%ob2dof_p(ibloc)%a(iobje),dhand%ob2dof_p(ibloc)%a(iobje+1)-1
          dof_g = dhand%ob2dof_l(ibloc)%a(idofs,1)
          ! Check the dof's are coupled in the system
          var_idof = dhand%ob2dof_l(ibloc)%a(idofs,2)
          ! jDOFs on objects
          do lobje = graph%ia(iobje),graph%ia(iobje+1)-1
             jobje = graph%ja(lobje)
             ! If jobje < femsp%o_msh%npoin
             do jdofs = dhand%ob2dof_p(jbloc)%a(jobje),dhand%ob2dof_p(jbloc)%a(jobje+1)-1
                var_jdof = dhand%ob2dof_l(jbloc)%a(jdofs,2)
                if(dhand%dof_coupl(var_idof,var_jdof)==0) cycle
                if ( gtype == csr ) then
                   dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                else ! gtype == csr_symm 
                   if ( dhand%ob2dof_l(jbloc)%a(jdofs,1) >= dof_g ) then
                      dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                   end if
                end if
             end do
             ! Else
             ! jelem = iobje-umesh%npoin
             ! check jelem is not in the loop below and add th einterior nodes
             ! this can wait a
          end do
          ! jDOFs in elements
          if ( present(femsp) ) then
             if (femsp%kfl_scond == scond_off) then
                do jelem = dmesh%pnods(iobje), dmesh%pnods(iobje+1)-1
                   gelem = dmesh%lnods(jelem)
                   iprob = femsp%lelem(gelem)%prob
                   do ivars = dhand%i_varsxprob(jbloc)%a(iprob),dhand%i_varsxprob(jbloc)%a(iprob+1)-1
                      gvars = dhand%j_varsxprob(jbloc)%a(ivars)
                      mtype = femsp%lelem(gelem)%iv(gvars)
                      if(dhand%dof_coupl(var_idof,gvars)==0) cycle
                      if ( gtype == csr ) then
                         dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + femsp%lelem(gelem)%        &
                              & f_inf(mtype)%p%nodes_obj(mesh%ndime+1)
                      else ! gtype == csr_symm
                         lobje = femsp%lelem(gelem)%f_inf(mtype)%p%nobje_dim(mesh%ndime+1)
                         do jnode = femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_i(lobje),               &
                              &     femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_i(lobje+1)-1
                            gnode = femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_j(jnode)
                            if ( dof_g <= femsp%lelem(gelem)%elem2dof(gnode,gvars)) then
                               dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                            end if
                         end do
                      end if
                   end do
                end do
             end if
          end if
       end do
    end do

    ! Count iDOFs in elements
    if ( present(femsp) ) then
       if (femsp%kfl_scond == scond_off) then
          do ielem = 1,mesh%nelem
             iprob = femsp%lelem(ielem)%prob
             do ivars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                gvars = dhand%j_varsxprob(ibloc)%a(ivars)
                mtype = femsp%lelem(ielem)%iv(gvars) 
                tobje = femsp%lelem(ielem)%f_inf(mtype)%p%nobje_dim(mesh%ndime+1)
                do inode = femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_i(tobje),               &
                     &     femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_i(tobje+1)-1
                   gnode = femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_j(inode)
                   dof_g = femsp%lelem(ielem)%elem2dof(gnode,gvars)
                   ! jDOFs on objects
                   do lobje = femsp%o_mesh%pnods(ielem), femsp%o_mesh%pnods(ielem+1)-1
                      jobje = femsp%o_mesh%lnods(lobje)
                      do jdofs = dhand%ob2dof_p(jbloc)%a(jobje),dhand%ob2dof_p(jbloc)%a(jobje+1)-1
                         var_jdof = dhand%ob2dof_l(jbloc)%a(jdofs,2)
                         if(dhand%dof_coupl(gvars,var_jdof)==0) cycle
                         if ( gtype == csr ) then
                            dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                         else ! gtype == csr_symm
                            if (dhand%ob2dof_l(jbloc)%a(jdofs,1) >= dof_g) then
                               dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                            end if
                         end if
                      end do
                   end do
                   ! jDOFs in ielem
                   do jvars = dhand%i_varsxprob(jbloc)%a(iprob),dhand%i_varsxprob(jbloc)%a(iprob+1)-1
                      lvars = dhand%j_varsxprob(jbloc)%a(jvars)
                      ntype = femsp%lelem(ielem)%iv(lvars)
                      if(dhand%dof_coupl(gvars,lvars)==0) cycle
                      if ( gtype == csr ) then
                         dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + &
                              femsp%lelem(ielem)%f_inf(ntype)%p%nodes_obj(mesh%ndime+1)
                      else ! gtype == csr_symm
                         jobje = femsp%lelem(ielem)%f_inf(ntype)%p%nobje_dim(mesh%ndime+1)
                         do jnode = femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_i(jobje),               &
                              &     femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_i(jobje+1)-1
                            gnode = femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_j(jnode)
                            if ( dof_g <= femsp%lelem(ielem)%elem2dof(gnode,lvars) ) then
                               dof_graph%ia(dof_g+1) = dof_graph%ia(dof_g+1) + 1
                            end if
                         end do
                      end if
                   end do
                end do
             end do
          end do
       end if
    end if

    dof_graph%ia(1) = 1
    do idofs = 2,dof_graph%nv+1
       dof_graph%ia(idofs) = dof_graph%ia(idofs-1) + dof_graph%ia(idofs)
    end do

    call memalloc( dof_graph%ia(dof_graph%nv+1)-1, dof_graph%ja, __FILE__,__LINE__ )

    ! List iDOFs on objects
    do iobje = 1,mesh%npoin
       do idofs = dhand%ob2dof_p(ibloc)%a(iobje),dhand%ob2dof_p(ibloc)%a(iobje+1)-1
          dof_g = dhand%ob2dof_l(ibloc)%a(idofs,1)
          count = dof_graph%ia(dof_g)
          ! Check the dof's are coupled in the system
          var_idof = dhand%ob2dof_l(ibloc)%a(idofs,2)
          ! jDOFs on objects
          do lobje = graph%ia(iobje),graph%ia(iobje+1)-1
             jobje = graph%ja(lobje)
             do jdofs = dhand%ob2dof_p(jbloc)%a(jobje),dhand%ob2dof_p(jbloc)%a(jobje+1)-1
                var_jdof = dhand%ob2dof_l(jbloc)%a(jdofs,2)
                if(dhand%dof_coupl(var_idof,var_jdof)==0) cycle
                if ( gtype == csr_symm ) then 
                   if (dhand%ob2dof_l(jbloc)%a(jdofs,1) >= dof_g) then
                      dof_graph%ja(count) = dhand%ob2dof_l(jbloc)%a(jdofs,1)
                      count = count+1
                   end if
                else
                   dof_graph%ja(count) = dhand%ob2dof_l(jbloc)%a(jdofs,1)
                   count = count+1
                end if
             end do
          end do
          ! jDOFs in elements
          if ( present(femsp) ) then
             if (femsp%kfl_scond == scond_off) then
                do jelem = dmesh%pnods(iobje), dmesh%pnods(iobje+1)-1
                   gelem = dmesh%lnods(jelem)
                   jprob = femsp%lelem(gelem)%prob
                   do jvars = dhand%i_varsxprob(jbloc)%a(jprob),dhand%i_varsxprob(jbloc)%a(jprob+1)-1
                      lvars = dhand%j_varsxprob(jbloc)%a(jvars)
                      mtype = femsp%lelem(gelem)%iv(lvars)
                      if(dhand%dof_coupl(var_idof,lvars)==0) cycle
                      lobje = femsp%lelem(gelem)%f_inf(mtype)%p%nobje_dim(mesh%ndime+1)
                      do jnode = femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_i(lobje),               &
                           &     femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_i(lobje+1)-1
                         gnode = femsp%lelem(gelem)%f_inf(mtype)%p%ndxob_j(jnode)
                         if ( gtype == csr_symm ) then
                            if ( dof_g <= femsp%lelem(gelem)%elem2dof(gnode,lvars)) then    
                               dof_graph%ja(count) = femsp%lelem(gelem)%elem2dof(gnode,lvars)
                               count = count+1
                            end if
                         else
                            dof_graph%ja(count) = femsp%lelem(gelem)%elem2dof(gnode,lvars)
                            count = count+1
                         end if
                      end do
                   end do
                end do
             end if
          end if
          ! Order increasingly column identifiers of current row 
          ! using heap sort algorithm
          call psb_hsort(dof_graph%ja( dof_graph%ia(dof_g):dof_graph%ia(dof_g+1)-1 ))
       end do
    end do

    ! List iDOFs in elements
    if ( present(femsp) ) then
       if (femsp%kfl_scond == scond_off) then
          do ielem = 1,mesh%nelem
             iprob = femsp%lelem(ielem)%prob
             do ivars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                gvars = dhand%j_varsxprob(ibloc)%a(ivars)
                mtype = femsp%lelem(ielem)%iv(gvars)
                tobje = femsp%lelem(ielem)%f_inf(mtype)%p%nobje_dim(mesh%ndime+1)
                do inode = femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_i(tobje),                        &
                     &     femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_i(tobje+1)-1
                   gnode = femsp%lelem(ielem)%f_inf(mtype)%p%ndxob_j(inode)
                   dof_g = femsp%lelem(ielem)%elem2dof(gnode,gvars)
                   count = dof_graph%ia(dof_g)
                   ! jDOFs on objects
                   do lobje = femsp%o_mesh%pnods(ielem), femsp%o_mesh%pnods(ielem+1)-1
                      jobje = femsp%o_mesh%lnods(lobje)
                      do jdofs = dhand%ob2dof_p(jbloc)%a(jobje),dhand%ob2dof_p(jbloc)%a(jobje+1)-1
                         var_jdof = dhand%ob2dof_l(jbloc)%a(jdofs,2)
                         if(dhand%dof_coupl(gvars,var_jdof)==0) cycle
                         if ( gtype == csr_symm ) then
                            if (dhand%ob2dof_l(jbloc)%a(jdofs,1) >= dof_g) then
                               dof_graph%ja(count) = dhand%ob2dof_l(jbloc)%a(jdofs,1)
                               count = count+1
                            end if
                         else
                            dof_graph%ja(count) = dhand%ob2dof_l(jbloc)%a(jdofs,1)
                            count = count+1
                         end if
                      end do
                   end do
                   ! jDOFs in ielem
                   jprob = femsp%lelem(ielem)%prob
                   do jvars = dhand%i_varsxprob(jbloc)%a(jprob),dhand%i_varsxprob(jbloc)%a(jprob+1)-1
                      lvars = dhand%j_varsxprob(jbloc)%a(jvars)
                      ntype = femsp%lelem(ielem)%iv(lvars)
                      if(dhand%dof_coupl(gvars,lvars)==0) cycle
                      tobje = femsp%lelem(ielem)%f_inf(ntype)%p%nobje_dim(mesh%ndime+1)
                      do jnode = femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_i(tobje),               &
                           &     femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_i(tobje+1)-1
                         gnode = femsp%lelem(ielem)%f_inf(ntype)%p%ndxob_j(jnode)
                         if ( gtype == csr_symm ) then
                            if ( dof_g <= femsp%lelem(ielem)%elem2dof(gnode,lvars) ) then
                               dof_graph%ja(count) = femsp%lelem(ielem)%elem2dof(gnode,lvars)
                               count = count+1
                            end if
                         else
                            dof_graph%ja(count) = femsp%lelem(ielem)%elem2dof(gnode,lvars)
                            count = count+1
                         end if
                      end do
                   end do
                   ! Order increasingly column identifiers of current row 
                   ! using heap sort algorithm
                   call psb_hsort(dof_graph%ja( dof_graph%ia(dof_g):dof_graph%ia(dof_g+1)-1 ))
                end do
             end do
          end do
       end if
    end if

    if ( present(femsp) ) then
       call fem_mesh_free(dmesh)
    end if

  end subroutine graph_to_dof_graph_block

  !=============================================================================
  subroutine dG_dof_graph_mono_create( dhand, dof_graph,femsp)
    implicit none
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space)          , intent(in)  :: femsp

    call dG_dof_graph_create(1,1,dhand,dof_graph,femsp)

  end subroutine dG_dof_graph_mono_create

  !=============================================================================
  ! This routine identifies the DOF that are coupled in the discontinuous Galerkin method.
  ! In this case two DOF are coupled if they are in the same element or they are involved 
  ! in integration of a common face. In this second case the DOF coupled with the nodes on
  ! the face are coupled with all the DOF of the neighbouring element while the rest of the 
  ! DOF of the element are only coupled with the DOF of the neighbouring element asssociated 
  ! to some node on the face. 
  ! WARNING! There is some places where scond is considered but this will definetely NOT WORK for
  ! the STATIC CONDENSATION CASE
  ! TO DO: Consider the fact that different variables might not be coupled
  ! TO DO: Consider different types of problems.
  subroutine dG_dof_graph_block_create( ibloc, jbloc, dhand, dof_graph,femsp)
    implicit none
    integer(ip)              , intent(in)  :: ibloc,jbloc
    type(dof_handler)        , intent(in)  :: dhand
    type(fem_graph)          , intent(out) :: dof_graph
    type(fem_space)          , intent(in)  :: femsp

    integer(ip)   :: ielem, ivars, gvars, inode, iprob, iface, nnode
    integer(ip)   :: active_nodes_elem(2), active_nodes_face(2)
    integer(ip)   :: jvars, hvars, jnode, jprob, fnode, ft, ndime, ig
    integer(ip)   :: elems(2), nodes(2), nodfs(2), prob, nvars(2), iobje(2), jg(2),i
    integer(ip)   :: g_dof, h_dof, posit(dhand%ndofs(ibloc))
    logical       :: scond

    scond = (femsp%kfl_scond == scond_on)

    ! Initialize the graph
    ndime = femsp%g_mesh%ndime
    dof_graph%type = csr
    dof_graph%nv  = dhand%ndofs(ibloc)
    dof_graph%nv2 = dhand%ndofs(jbloc)
    call memalloc( dof_graph%nv+1, dof_graph%ia, 'dG_dof_graph_block_create::dof_graph%ia' )
    dof_graph%ia = 0

    ! Count DOFs
    
    ! Count the DOF coupled by volumetric integration
    do ielem = 1, femsp%g_mesh%nelem
       iprob = femsp%lelem(ielem)%prob
       
       ! Loop over the variables of the problem
       do ivars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1

          ! Identify the local id of the variable
          gvars = dhand%j_varsxprob(ibloc)%a(ivars)

          ! Identify interpolation for the variable 
          ig    = femsp%lelem(ielem)%iv(gvars)

          ! Compute the number of nodes in the element
          if (scond) then
             nnode = femsp%lelem(ielem)%f_inf(ig)%p%nnode -               &
               femsp%lelem(ielem)%f_inf(ig)%p%nodes_obj(ndime+1)
          else
             nnode = femsp%lelem(ielem)%f_inf(ig)%p%nnode
          end if

          ! Compute the amount of active nodes (nodes with associated DOF for that variable)
          active_nodes_elem(1) = 0
          do inode = 1,nnode
             g_dof = femsp%lelem(ielem)%elem2dof(inode,gvars)
             if (g_dof > 0) then
                active_nodes_elem(1) = active_nodes_elem(1) + 1
             end if
          end do
          
          ! Add this active nodes to the count of the DOF associated to (inode,gvars)
          do inode = 1,nnode
             g_dof = femsp%lelem(ielem)%elem2dof(inode,gvars)
             if (g_dof > 0) then
                dof_graph%ia(g_dof+1) = dof_graph%ia(g_dof+1) + active_nodes_elem(1) * &
                     (dhand%i_varsxprob(jbloc)%a(iprob+1) - dhand%i_varsxprob(jbloc)%a(iprob))
             end if
          end do
       end do
    end do

    ! Count the DOF coupled by face integration
    do iface = 1, femsp%o_mesh%tface
       
       ! Identify neighbouring elements
       elems = femsp%lface(iface)%nei_elem

       ! For the moment we are only considering one type of problem
       assert(femsp%lelem(elems(1))%prob==femsp%lelem(elems(2))%prob)
       prob = femsp%lelem(elems(1))%prob

       ! Loop over the variables of the problem
       do ivars =  dhand%i_varsxprob(ibloc)%a(prob),dhand%i_varsxprob(ibloc)%a(prob+1)-1

          ! Identify the local id of the variable
          gvars = dhand%j_varsxprob(ibloc)%a(ivars)
          
          ! Loop over the neighbours
          do i=1,2
             
             ! Identify interpolation for the variable 
             jg(i) = femsp%lelem(elems(i))%iv(gvars)

             ! Compute the number of nodes in the element i and on the face in that element
             if (scond) then
                nodes(i) = femsp%lelem(elems(i))%f_inf(jg(i))%p%nnode -   &
                           femsp%lelem(elems(i))%f_inf(jg(i))%p%nodes_obj(ndime+1)
             else
                nodes(i) = femsp%lelem(elems(i))%f_inf(jg(i))%p%nnode 
             end if
             iobje(i) = femsp%lelem(elems(i))%f_inf(jg(i))%p%nobje_dim(ndime) +                     &
                  &     femsp%lface(iface)%pos_elem(i)-1
             nodfs(i) = femsp%lelem(elems(i))%f_inf(1)%p%ntxob_i(iobje(i)+1) -                      &
                  &     femsp%lelem(elems(i))%f_inf(1)%p%ntxob_i(iobje(i))
             nvars(i) = dhand%i_varsxprob(jbloc)%a(prob+1) - dhand%i_varsxprob(jbloc)%a(prob)
          end do

          ! Compute the amount of active nodes (nodes with associated DOF for that variable)
          do i = 1,2
             active_nodes_elem(i) = 0
             active_nodes_face(i) = 0
             do inode = 1,nodes(i)
                g_dof = femsp%lelem(elems(i))%elem2dof(inode,gvars)
                if (g_dof > 0) then
                   active_nodes_elem(i) = active_nodes_elem(i) + 1
                end if
             end do
             do inode = femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i)),                    &
                  &     femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i)+1)-1
                g_dof = femsp%lelem(elems(i))%elem2dof(femsp%lelem(elems(i))%f_inf(jg(i))          &
                     &       %p%ntxob_j(inode),gvars)
                if (g_dof > 0) then
                   active_nodes_face(i) = active_nodes_face(i) + 1
                end if
             end do
          end do
          
          ! Add this active nodes to the count of the DOF associated tot (inode,gvars)
          do i = 1,2
             ! DOFS on the face in the neighbouring element
             do inode = 1,nodes(i)
                g_dof = femsp%lelem(elems(i))%elem2dof(inode,gvars)
                if (g_dof > 0) then
                   dof_graph%ia(g_dof+1) = dof_graph%ia(g_dof+1) + active_nodes_face(3-i)*nvars(3-i)
                end if
             end do

             ! DOFS in the neighbouring element but NOT on the face (only for the DOF on the face)
             do inode = femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i)),                    &
                  &     femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i)+1)-1
                g_dof = femsp%lelem(elems(i))%elem2dof(femsp%lelem(elems(i))%f_inf(jg(i))          &
                     &       %p%ntxob_j(inode),gvars)
                if (g_dof > 0) then
                   dof_graph%ia(g_dof+1) = dof_graph%ia(g_dof+1)+(active_nodes_elem(3-i) -         &
                                           active_nodes_face(3-i))*nvars(3-i)
                end if
             end do
          end do
       end do
    end do

    ! Create ia
    dof_graph%ia(1) = 1
    do g_dof = 2,dof_graph%nv+1
       dof_graph%ia(g_dof) = dof_graph%ia(g_dof-1) + dof_graph%ia(g_dof)
    end do

    ! Allocate ja
    call memalloc( dof_graph%ia(dof_graph%nv+1)-1,dof_graph%ja,'dG_dof_graph_create::dof_graph%ja')
    posit = 0; dof_graph%ja = 0

    ! List DOFs in ja
    ! Fill ia for DOFs coupled by the volumetric integration
    do ielem = 1, femsp%g_mesh%nelem
       iprob = femsp%lelem(ielem)%prob

       ! Compute the number of nodes in the element
       do ivars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
          gvars = dhand%j_varsxprob(ibloc)%a(ivars)
          ig = femsp%lelem(ielem)%iv(gvars)
          if (scond) then
             nnode = femsp%lelem(ielem)%f_inf(ig)%p%nnode -                &
                femsp%lelem(ielem)%f_inf(ig)%p%nodes_obj(ndime+1)
          else
             nnode = femsp%lelem(ielem)%f_inf(ig)%p%nnode
          end if

          ! DOFs in the same element
          do inode = 1,nnode
             g_dof = femsp%lelem(ielem)%elem2dof(inode,gvars)
             if (g_dof > 0) then
                do jnode = 1,nnode
                   do jvars = dhand%i_varsxprob(jbloc)%a(iprob),dhand%i_varsxprob(jbloc)%a(iprob+1)-1
                      hvars = dhand%j_varsxprob(jbloc)%a(jvars)
                      h_dof = femsp%lelem(ielem)%elem2dof(jnode,hvars)
                      if (h_dof > 0) then
                         dof_graph%ja(dof_graph%ia(g_dof)+posit(g_dof)) = h_dof
                         posit(g_dof) = posit(g_dof) + 1
                      end if
                   end do
                end do
             end if
          end do
       end do
    end do
   
    ! DOFs coupled by faces
    do iface = 1, femsp%o_mesh%tface
       elems = femsp%lface(iface)%nei_elem
       assert(femsp%lelem(elems(1))%prob==femsp%lelem(elems(2))%prob)
       prob = femsp%lelem(elems(1))%prob

       ! Loop over the variables
       do ivars =  dhand%i_varsxprob(ibloc)%a(prob),dhand%i_varsxprob(ibloc)%a(prob+1)-1
          gvars = dhand%j_varsxprob(ibloc)%a(ivars)
          do i=1,2
             jg(i) = femsp%lelem(elems(i))%iv(gvars)

             ! Compute the number of nodes in the element i and on the face in that element
             if (scond) then
                nodes(i) = femsp%lelem(elems(i))%f_inf(jg(i))%p%nnode -                             &
                           femsp%lelem(elems(i))%f_inf(jg(i))%p%nodes_obj(ndime)
             else
                nodes(i) = femsp%lelem(elems(i))%f_inf(jg(i))%p%nnode
             end if
             iobje(i) = femsp%lelem(elems(i))%f_inf(1)%p%nobje_dim(ndime)  +                        &
                  &     femsp%lface(iface)%pos_elem(i)-1
             nodfs(i) = femsp%lelem(elems(i))%f_inf(1)%p%ntxob_i(iobje(i)+1) -                      &
                  &     femsp%lelem(elems(i))%f_inf(1)%p%ntxob_i(iobje(i))
             nvars(i) = dhand%i_varsxprob(jbloc)%a(prob+1) - dhand%i_varsxprob(jbloc)%a(prob)
          end do
 
          do i = 1,2
             fnode = 1
             ft = femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i))+fnode-1
             do inode = 1,nodes(i)
                g_dof = femsp%lelem(elems(i))%elem2dof(inode,gvars)
                if ( (inode == femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_j(ft)) .and.              &
                     & (fnode <= nodfs(i)) ) then
                   ! DOF on face_idx
                   if (g_dof > 0) then
                      do jvars =  dhand%i_varsxprob(jbloc)%a(prob),                                 & 
                           &      dhand%i_varsxprob(jbloc)%a(prob+1)-1
                         hvars = dhand%j_varsxprob(jbloc)%a(jvars)

                         ! Compute number of nodes in the neighbouring element
                         if (scond) then
                            nnode = femsp%lelem(elems(3-i))%f_inf(femsp%lelem(elems(3-i))%iv(hvars))%p  &
                                 &  %nnode - femsp%lelem(elems(3-i))%f_inf(femsp%lelem(elems(i))%     &
                                 &  iv(gvars))%p%nodes_obj(ndime+1)
                         else
                            nnode = femsp%lelem(elems(3-i))%f_inf(femsp%lelem(elems(i))%iv(hvars))%p  &
                                 &   %nnode
                         end if

                         ! Add all the DOF in the neighbouring element in ia
                         do jnode = 1, nnode
                            h_dof = femsp%lelem(elems(3-i))%elem2dof(jnode,hvars)
                            if (h_dof > 0) then
                               dof_graph%ja(dof_graph%ia(g_dof)+posit(g_dof)) = h_dof
                               posit(g_dof) = posit(g_dof) + 1                               
                            end if
                         end do
                      end do
                   end if
                   
                   ! Increase face node id.
                   fnode = fnode + 1
                   ft = femsp%lelem(elems(i))%f_inf(jg(i))%p%ntxob_i(iobje(i))+fnode-1
                else
                   ! DOF NOT on face_idx
                   if (g_dof > 0) then
                      do jvars =  dhand%i_varsxprob(jbloc)%a(prob),                        &
                           &      dhand%i_varsxprob(jbloc)%a(prob+1)-1
                         hvars = dhand%j_varsxprob(jbloc)%a(jvars)

                         ! Add all DOF in the neighbouring element associated to nodes on the face
                         do jnode = femsp%lelem(elems(3-i))%f_inf(jg(3-i))%p%ntxob_i(iobje(3-i)),   &
                              &     femsp%lelem(elems(3-i))%f_inf(jg(3-i))%p%ntxob_i(iobje(3-i)+1)-1
                            h_dof = femsp%lelem(elems(3-i))%elem2dof(femsp%lelem(elems(3-i))        &
                                 &       %f_inf(jg(3-i))%p%ntxob_j(jnode),hvars)
                            if (h_dof > 0) then
                               dof_graph%ja(dof_graph%ia(g_dof)+posit(g_dof)) = h_dof
                               posit(g_dof) = posit(g_dof) + 1                               
                            end if
                         end do
                      end do
                   end if
                end if
             end do
          end do
       end do
    end do

    ! Sort the entries
    do g_dof = 1, dof_graph%nv
       call psb_hsort(dof_graph%ja( dof_graph%ia(g_dof):dof_graph%ia(g_dof+1)-1 ))
    end do
  end subroutine dG_dof_graph_block_create

  !=============================================================================
  subroutine object2dof_create( ibloc, mesh, dual_mesh, f_dof_handler, femsp, fixities, f_part, new_f_part )
    implicit none
    ! Parameters
    integer(ip)   , intent(in)                    :: ibloc
    type(fem_mesh), intent(in)                    :: mesh
    type(fem_mesh), intent(in)                    :: dual_mesh
    type(dof_handler), intent(inout)              :: f_dof_handler
    type(fem_space), intent(in)                   :: femsp
    type(fem_conditions),intent(in)               :: fixities
    type(fem_partition), optional, intent(in)     :: f_part
    type(fem_partition), optional, intent(inout)  :: new_f_part

    call object2dof_count( ibloc, mesh, dual_mesh, f_dof_handler, femsp, fixities, f_part, new_f_part )

    call object2dof_list ( ibloc, mesh, dual_mesh, f_dof_handler, femsp, fixities )

  end subroutine object2dof_create

  !=============================================================================
  subroutine object2dof_count( ibloc, mesh, dual_mesh, dhand, femsp, fixities, f_part, new_f_part )
    implicit none
    ! Parameters
    integer(ip)   , intent(in)                    :: ibloc
    type(fem_mesh), intent(in)                    :: mesh
    type(fem_mesh), intent(in)                    :: dual_mesh
    type(dof_handler), intent(inout)              :: dhand
    type(fem_space), intent(in)                   :: femsp
    type(fem_conditions),intent(in)               :: fixities
    type(fem_partition), optional, intent(in)     :: f_part
    type(fem_partition), optional, intent(inout)  :: new_f_part

    ! Locals
    integer(ip) :: touch(dhand%nvars,2), iposi, iobje, lelem, gelem, obje_loc, obje_dim
    integer(ip) :: obje_nod, iprob, kvars, lvars, gvars, inode, intcount, l_cing_idx, cing_node, num_cing_dofs
    type(adapt_constraint), pointer      :: constraint_ptr
    logical     :: add_new_DOFs
    ! Allocate ob2dof_p
    call array_create(mesh%npoin+1,dhand%ob2dof_p(ibloc))
    call array_create(dhand%ndofs(ibloc),dhand%dof2ob_l(ibloc))

    iposi = 1
    intcount = 0
    dhand%ob2dof_p(ibloc)%a = 0    
    do iobje = 1,mesh%npoin
       touch = 0    
       do lelem = dual_mesh%pnods(iobje), dual_mesh%pnods(iobje+1)-1
          gelem = dual_mesh%lnods(lelem)
          iprob = femsp%lelem(gelem)%prob
          kvars = 0
          do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
             kvars = kvars+1
             gvars = dhand%j_varsxprob(ibloc)%a(lvars)
             do obje_loc = 1,femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nobje
                if (mesh%lnods(mesh%pnods(gelem)+obje_loc-1) == iobje ) exit
             end do
             do obje_dim = 1,mesh%ndime
                if ( obje_loc < femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nobje_dim(obje_dim+1) ) exit
             end do
             if(obje_dim>mesh%ndime) then
                obje_nod = 0
             else
                obje_nod = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nodes_obj(obje_dim)
             end if
             
             ! constrained nodes
             if (associated(femsp%constraint_list)) then
                if(is_constrained_nodal(femsp%constraint_list, iobje, constraint_ptr) .and.            &
                     & dhand%continuity == cG_cont) then
                   ! for the moment only for scalar equation
                   assert(obje_nod <= 1)
                   obje_nod = 0
                   ! todo: it should be decreased for essential bc in cing nodes
                   do l_cing_idx = 1, constraint_ptr%num_cing_nodes
                      cing_node = constraint_ptr%cing_nodes(l_cing_idx)

                      ! constraining nodes should have smaller id and thus their num of dofs should 
                      ! be allready counted
                      assert(cing_node < iobje)
                      num_cing_dofs = dhand%ob2dof_p(ibloc)%a(cing_node + 1)

                      ! todo: for higher order it will have to be organized slightly differently
                      assert(num_cing_dofs <= 1)
                      obje_nod = obje_nod + num_cing_dofs
                   end do
                end if
             end if
             
             if ( fixities%code(gvars,iobje) /= 1 ) then
                ! TODO: For 3D, on edges, this is not optimal implemented (see elem2dof_create)
                if ( touch(gvars,1) == 0 .or. dhand%continuity == dG_cont ) then
                   add_new_DOFs = .true.
                elseif (obje_dim == 1 .or.femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars)) &
                     &                      %p%order == touch(gvars,2)) then
                   add_new_DOFs = .false.
                else
                   assert(dhand%continuity == cdG_cont )
                   add_new_DOFs = .true.
                end if
                if (add_new_DOFs) then
                   ! BE CAREFUL: Think about boundary conditions (...)
                   ! Create element_dofs(ielem,ivar,inode)
                   dhand%ob2dof_p(ibloc)%a(iobje+1) = dhand%ob2dof_p(ibloc)%a(iobje+1) + obje_nod
                   if (touch(gvars,1) == 0) then
                      touch(gvars,1) = 1
                      touch(gvars,2) = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%order
                   end if
                   iposi = iposi + obje_nod
                end if
             end if
          end do
       end do
       if (present(f_part) .and. present(new_f_part)) then
          if( iobje == f_part%nmap%ni ) then
             intcount = 0
             do gelem = 1, mesh%nelem
                iprob = femsp%lelem(gelem)%prob
                do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
                   gvars = dhand%j_varsxprob(ibloc)%a(lvars)
                   intcount = intcount + femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))     &
                        %p%nodes_obj(mesh%ndime+1)
                end do
             end do
             iposi = iposi + intcount
             new_f_part%nmap%ni = iposi - 1 
          end if
       end if
    end do

    iposi = iposi - 1

    if (present(f_part) .and. present(new_f_part)) then
       new_f_part%nmap%nl = iposi 
       new_f_part%nmap%nb = iposi - new_f_part%nmap%ni
       new_f_part%nmap%ne = 0
    end if
    !dhand%ndofs(ibloc) = iposi

    dhand%ob2dof_p(ibloc)%a(1) = 1
    do iobje = 2,mesh%npoin+1
       dhand%ob2dof_p(ibloc)%a(iobje) = dhand%ob2dof_p(ibloc)%a(iobje-1) + dhand%ob2dof_p(ibloc)%a(iobje)
    end do
    
    ! Allocate ob2dof_l 
    call array_create(iposi-intcount,2,dhand%ob2dof_l(ibloc))

  end subroutine object2dof_count

  !=============================================================================
  subroutine object2dof_list( ibloc, mesh, dual_mesh, dhand, femsp, fixities )
    implicit none
    ! Parameters
    integer(ip)   ,       intent(in)    :: ibloc
    type(fem_mesh),       intent(in)    :: mesh
    type(fem_mesh),       intent(in)    :: dual_mesh
    type(dof_handler),    intent(inout) :: dhand
    type(fem_space),      intent(in)    :: femsp
    type(fem_conditions), intent(in)    :: fixities

    ! Locals
    integer(ip) :: touch(dhand%nvars,2), iposi, iobje, lelem, gelem, obje_loc, obje_dim
    integer(ip) :: obje_nod, obje_sta, obje_end, iprob, kvars, lvars
    integer(ip) :: gvars, inode, jnode, act_elem2dof, l_cing_node, cing_obj_idx, cing_dof, cing_node
    type(adapt_constraint), pointer      :: constraint_ptr
    integer(ip) :: list_idx, list_idx_from, list_idx_next, list_idx_from_cing, list_idx_next_cing
    integer(ip) :: num_cing_dofs
    logical     :: add_new_DOFs

    iposi = 1
    do iobje = 1,mesh%npoin
       touch = 0     
       do lelem = dual_mesh%pnods(iobje), dual_mesh%pnods(iobje+1)-1
          gelem = dual_mesh%lnods(lelem)
          iprob = femsp%lelem(gelem)%prob
          kvars = 0
          do lvars = dhand%i_varsxprob(ibloc)%a(iprob),dhand%i_varsxprob(ibloc)%a(iprob+1)-1
             kvars = kvars+1
             gvars = dhand%j_varsxprob(ibloc)%a(lvars)
             do obje_loc = 1,femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nobje
                if (mesh%lnods(mesh%pnods(gelem)+obje_loc-1) == iobje ) exit
             end do
             do obje_dim = 1,mesh%ndime
                if ( obje_loc < femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nobje_dim(obje_dim+1) ) exit
             end do
             if(obje_dim>mesh%ndime) then
                obje_nod = 0
             else
                obje_nod = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%nodes_obj(obje_dim)
             end if
             obje_sta = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%ndxob_i(obje_loc)
             obje_end = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%ndxob_i(obje_loc+1)-1
             if ( fixities%code(gvars,iobje) /= 1 ) then
                if ( touch(gvars,1) == 0 .or. dhand%continuity == dG_cont ) then
                   add_new_DOFs = .true.
                elseif (obje_dim == 1 .or.femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars)) &
                     &                      %p%order == touch(gvars,2)) then
                   add_new_DOFs = .false.
                else
                   assert(dhand%continuity == cdG_cont )
                   add_new_DOFs = .true.
                end if
                !if ( touch(gvars) == 0 .or. dhand%continuity == dG_cont ) then
                if (add_new_DOFs) then
                   ! BE CAREFUL: Think about boundary conditions (...)
                   ! Create element_dofs(ielem,ivar,inode)
                   do jnode = obje_sta,obje_end
                      inode = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%ndxob_j(jnode)
                      act_elem2dof = femsp%lelem(gelem)%elem2dof(inode,gvars)
                      if(act_elem2dof >= 0) then
                         ! regular node  (should be > 0 ? )
                         dhand%ob2dof_l(ibloc)%a(iposi,1) = act_elem2dof
                         dhand%ob2dof_l(ibloc)%a(iposi,2) = gvars
                         dhand%dof2ob_l(ibloc)%a(dhand%ob2dof_l(ibloc)%a(iposi,1)) = iobje
                         iposi = iposi + 1
                      else
                         !hanging node
                         assert(act_elem2dof<0) 
                         constraint_ptr => femsp%constraint_list%list(-act_elem2dof)
                         !assert(is_constrained(femsp%constraint_list, -act_elem2dof, constraint_ptr))
                         do l_cing_node = 1, constraint_ptr%num_cing_nodes
                            ! we cannot be sure, that dofs of constraining nodes are allready assigned
                            ! we thus put negative object identifiers instead
                            ! they will be replaced in separate loop over dhand%ob2dof_l
                            ! todo: or probably it should allready be assigned...
                            cing_node = constraint_ptr%cing_nodes(l_cing_node)
                            
                            ! is it really done?
                            assert(cing_node < iobje)
                            num_cing_dofs = dhand%ob2dof_p(ibloc)%a(cing_node + 1) - dhand%ob2dof_p(ibloc)%a(cing_node)
                            assert(num_cing_dofs <= 1)
                            ! in the case of linear elements this means vertex on the essential boundary
                            if(num_cing_dofs == 0) then
                               cycle
                            end if
                            !dhand%ob2dof_l(ibloc)%a(iposi,1) = - cing_node
                            dhand%ob2dof_l(ibloc)%a(iposi,1) = dhand%ob2dof_l(ibloc)%a(dhand%ob2dof_p(ibloc)%a(cing_node), 1)
                            dhand%ob2dof_l(ibloc)%a(iposi,2) = gvars
                            
                            ! todo: I do not put anything to a dhand%dof2ob_l array...
                            
                            iposi = iposi + 1
                         end do
                      end if
                   end do
                   if (touch(gvars,1) == 0) then
                      touch(gvars,1) = 1
                      touch(gvars,2) = femsp%lelem(gelem)%f_inf(femsp%lelem(gelem)%iv(gvars))%p%order
                   end if
                end if
             end if
          end do
       end do
    end do

    ! not needed, if constrained nodes are really after constraining
    !
    !! now replace the negative object identifiers by actuall dofs
    !if(fem_space_has_constraints(femsp)) then
    !   ! not ready for anything else yet
    !   assert(ibloc == 1)
    !   do iobje = 1,mesh%npoin
    !      ! find indices range to ob2dof_l
    !      list_idx_from = dhand%ob2dof_p(ibloc)%a(iobje)
    !      list_idx_next = dhand%ob2dof_p(ibloc)%a(iobje+1)
    !      do list_idx = list_idx_from, list_idx_next - 1
    !         cing_obj_idx = -dhand%ob2dof_l(ibloc)%a(list_idx,1)
    !         ! if the value in ob2dof_l was negative, it is a constrained node, which has to be dealt with
    !         if(cing_obj_idx > 0) then
    !            ! indices to ob2dof_l of the constraining node
    !            list_idx_from_cing = dhand%ob2dof_p(ibloc)%a(cing_obj_idx)
    !            list_idx_next_cing = dhand%ob2dof_p(ibloc)%a(cing_obj_idx+1)

    !            ! the constraining node must have exactly one dof
    !            ! for the scalar equation it means, that it is not itself constrained and it should be guarantied by 2:1 rule
    !            ! todo: more complex equations
    !            ! todo: what if the constraining node is on essential BC?
    !            assert( list_idx_next_cing == list_idx_from_cing + 1)
    !            cing_dof = dhand%ob2dof_l(ibloc)%a(list_idx_from_cing, 1)
    !            
    !            ! todo: what for essential BC? 
    !            assert(cing_dof > 0)
    !            dhand%ob2dof_l(ibloc)%a(list_idx, 1) = cing_dof
    !         end if
    !      end do
    !   end do
    !end if

  end subroutine object2dof_list

  !=============================================================================
  subroutine dof_handler_print(lunou, d)
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(dof_handler), intent(in) :: d

    ! Local variables
    integer(ip) :: i,j

    write (lunou, '(a)')     '*** begin dof_handler data structure ***'
    write (lunou, '(a,i10)') 'Number of problems:', d%nprob    
    write (lunou, '(a,i10)') 'Number of physical variables:', d%nvars
    write (lunou, '(a,i10)') 'Number of dofs:', d%ndofs
    write (lunou, '(a,i10)') 'Continuity:', d%continuity

    do j=1,d%blocks%nb
       write (lunou, *) 'Block:',j
       write (lunou, '(a)') '(Problem, prob_code, number of vars, vars associated):'
       do i=1,d%nprob
          write (lunou, '(10i10)') i, d%phys_prob(i), d%blknvarsxprob(j)%a(i), & 
               d%j_varsxprob(j)%a(d%i_varsxprob(j)%a(i): d%i_varsxprob(j)%a(i+1)-1)
       end do

       if ( allocated(d%ob2dof_p(j)%a) ) then
          write(*,*) 'dhand%ob2dof_p',d%ob2dof_p(j)%a
          write(*,*) 'dhand%ob2dof_l',d%ob2dof_l(j)%a
!!$          write (lunou,'(a)') 'Object to dof array (ipoin,1:nvars):'
!!$          do i=1,npoin
!!$             write (lunou, '(i10)', advance='no')  i
!!$             !write (lunou, '(a)' ) 'ob2dof_l: '
!!$             do j=d%ob2dof_p(i), d%ob2dof_p(i+1)-1  !1,d%nvars-1
!!$                !write (*,*) 'j', j
!!$                write (lunou, '(2i10)', advance='no')  d%ob2dof_l(j,1), &
!!$                     d%ob2dof_l(j,2)          
!!$             end do
!!$             write (lunou, *) ' '
!!$             !write (*,*) 'j', j
!!$             !          write (lunou, '(2i10)')  d%ob2dof_l(j+1,1), 'var:', d%ob2dof_l(j+1,2)
!!$          end do
       end if
    end do

    !write (lunou,'(a)') 'Dof to node array (ndofs):'
    !do i=1,d%ndofs
    !   write (lunou, '(2i10)')  i, d%dof2node(i)
    !end do
    !write (lunou, '(a)')     '*** end dof_handler data structure ***'

  end subroutine dof_handler_print

  !=============================================================================
  subroutine dof_handler_free(d)
    implicit none
    type(dof_handler), intent(inout) :: d
    ! Locals
    integer(ip) :: i

    !do i=1,d%nprob
    !   call array_free(d%dof_coupl(i))
    !end do
    call memfree( d%dof_coupl,__FILE__,__LINE__)
    !deallocate(d%dof_coupl)

    d%nprob = 0
    d%nvars = 0
    d%continuity = 0 

    call memfree( d%phys_prob ,__FILE__,__LINE__)!
    call memfree( d%nvarsxprob ,__FILE__,__LINE__)!
    !call memfree( d%i_varsxprob,__FILE__,__LINE__)
    !call memfree( d%j_varsxprob,__FILE__,__LINE__)
    !call memfree( d%continuity,__FILE__,__LINE__)
    !call memfree( d%ob2dof_p,__FILE__,__LINE__)
    !call memfree( d%ob2dof_l,__FILE__,__LINE__)
    !call memfree( d%dof2node,__FILE__,__LINE__)

    do i=1,d%blocks%nb
       call array_free(d%ob2dof_p(i))
       call array_free(d%ob2dof_l(i))
       call array_free(d%dof2ob_l(i))
       call array_free(d%blknvarsxprob(i))!
       call array_free(d%i_varsxprob(i))!
       call array_free(d%j_varsxprob(i))!
    end do
    deallocate(d%ob2dof_p)!
    deallocate(d%ob2dof_l)!
    deallocate(d%dof2ob_l)!
    deallocate(d%blknvarsxprob)!
    deallocate(d%i_varsxprob)!
    deallocate(d%j_varsxprob)!
    call memfree(d%ndofs,__FILE__,__LINE__) !

    call fem_blocks_free(d%blocks)

  end subroutine dof_handler_free

  !=============================================================================
  subroutine dof_handler_create_one_physics ( nvars, code, cont, dhand, blks) 
    implicit none
    integer(ip)      , intent(in)           :: nvars, code, cont
    type(dof_handler), intent(out)          :: dhand
    type(fem_blocks) , optional, intent(in) :: blks

    integer(ip)  ::  i,j,k,ldofs(1),ibloc

    if(present(blks)) then
       call fem_blocks_copy(blks,dhand%blocks)
    else
       ldofs(1) = nvars
       call fem_blocks_alloc(scalar,1,ldofs,dhand%blocks)
    end if

    assert(dhand%blocks%ib(dhand%blocks%nb+1)-1==nvars)

    dhand%nprob = 1
    dhand%nvars = nvars
    dhand%continuity = cont
    call memalloc(dhand%nprob,dhand%phys_prob,__FILE__,__LINE__)
    dhand%phys_prob(1) = code 
    call memalloc(dhand%nprob,dhand%nvarsxprob,__FILE__,__LINE__)
    dhand%nvarsxprob(1) = nvars

    ! Dof coupling
    !call memalloc(dhand%nprob,dhand%dof_coupl,__FILE__,__LINE__)
    !allocate(dhand%dof_coupl(dhand%nprob))
    !call array_create(nvars,nvars,dhand%dof_coupl(1))
    !dhand%dof_coupl(1)%a = 1 ! By default, every dof is coupled with every dof
    call memalloc(nvars, nvars, dhand%dof_coupl,__FILE__,__LINE__)
    dhand%dof_coupl = 1 ! By default, every dof is coupled with every dof
    
    ! Block-problem variables (This only works for one-physics!)
    allocate(dhand%blknvarsxprob(dhand%blocks%nb))
    allocate(dhand%i_varsxprob(dhand%blocks%nb))
    allocate(dhand%j_varsxprob(dhand%blocks%nb))

    do ibloc = 1,dhand%blocks%nb
       call array_create(dhand%nprob,dhand%blknvarsxprob(ibloc))
       call array_create(dhand%nprob+1,dhand%i_varsxprob(ibloc))
       dhand%blknvarsxprob(ibloc)%a(1) = dhand%blocks%ib(ibloc+1)-dhand%blocks%ib(ibloc)
       dhand%i_varsxprob(ibloc)%a(1) = 1
       do i=1,dhand%nprob
          dhand%i_varsxprob(ibloc)%a(i+1) = dhand%i_varsxprob(ibloc)%a(i) + dhand%blknvarsxprob(ibloc)%a(i)
       end do
       call array_create(dhand%i_varsxprob(ibloc)%a(dhand%nprob+1)-1,dhand%j_varsxprob(ibloc))
       k=0
       do i=1,dhand%nprob
          do j=dhand%i_varsxprob(ibloc)%a(i),dhand%i_varsxprob(ibloc)%a(i+1)-1
             k=k+1
             dhand%j_varsxprob(ibloc)%a(k) = dhand%blocks%jb(j+dhand%blocks%ib(ibloc)-1)
          end do
       end do
    end do

    ! Allocate ob2dof variables
    allocate(dhand%ob2dof_p(dhand%blocks%nb))
    allocate(dhand%ob2dof_l(dhand%blocks%nb))
    allocate(dhand%dof2ob_l(dhand%blocks%nb))
    call memalloc(dhand%blocks%nb,dhand%ndofs,__FILE__,__LINE__)

  end subroutine dof_handler_create_one_physics

  !=============================================================================
  subroutine dof_handler_create_multiphysics ( nvars, nprob, prob_code, prob_nunk, prob_varini, &
       &                                       nunkxblock, cont, dhand, blks) 
    implicit none
    integer(ip)         , intent(in)        :: nvars, nprob, prob_code(nprob), prob_nunk(nprob)
    integer(ip)         , intent(in)        :: prob_varini(:,:), cont, nunkxblock(:,:)
    type(dof_handler)   , intent(out)       :: dhand
    type(fem_blocks) , optional, intent(in) :: blks

    integer(ip)  ::  i,j,k,ldofs(1),ibloc,iprob

    if(present(blks)) then
       call fem_blocks_copy(blks,dhand%blocks)
    else
       ldofs(1) = nvars
       call fem_blocks_alloc(scalar,1,ldofs,dhand%blocks)
    end if

    assert(dhand%blocks%ib(dhand%blocks%nb+1)-1==nvars)
    do iprob=1,nprob
       assert(sum(nunkxblock(:,iprob))==prob_nunk(iprob))
    end do

    dhand%nprob = nprob
    dhand%nvars = nvars
    dhand%continuity = cont
    call memalloc(dhand%nprob,dhand%phys_prob,__FILE__,__LINE__)
    call memalloc(dhand%nprob,dhand%nvarsxprob,__FILE__,__LINE__)
    do iprob=1,nprob
       dhand%phys_prob(iprob) = prob_code(iprob) 
       dhand%nvarsxprob(iprob) = prob_nunk(iprob)
    end do
    
    ! Dof coupling
    !do iprob=1,nprob
    !   call array_create(dhand%nvarsxprob(iprob), dhand%nvarsxprob(iprob), dhand%dof_coupl(iprob))
    !   dhand%dof_coupl(iprob)%a = 1 ! By default, every dof is coupled with every dof
    !end do
    call memalloc(nvars, nvars, dhand%dof_coupl,__FILE__,__LINE__)
    dhand%dof_coupl = 1 ! By default, every dof is coupled with every dof

    ! Block-problem variables (This only works for one-physics!)
    allocate(dhand%blknvarsxprob(dhand%blocks%nb))
    allocate(dhand%i_varsxprob(dhand%blocks%nb))
    allocate(dhand%j_varsxprob(dhand%blocks%nb))

    do ibloc = 1,dhand%blocks%nb
       call array_create(dhand%nprob,dhand%blknvarsxprob(ibloc))
       call array_create(dhand%nprob+1,dhand%i_varsxprob(ibloc))
       call array_create(sum(nunkxblock(ibloc,:)),dhand%j_varsxprob(ibloc))
       dhand%i_varsxprob(ibloc)%a(1) = 1
       do i=1,dhand%nprob
          dhand%blknvarsxprob(ibloc)%a(i) = nunkxblock(ibloc,i)
          dhand%i_varsxprob(ibloc)%a(i+1) = dhand%i_varsxprob(ibloc)%a(i) + nunkxblock(ibloc,i)
          do j=1,nunkxblock(ibloc,i)
             dhand%j_varsxprob(ibloc)%a(dhand%i_varsxprob(ibloc)%a(i)-1+j) = prob_varini(ibloc,i) + j-1
          end do
       end do
    end do

    ! Allocate ob2dof variables
    allocate(dhand%ob2dof_p(dhand%blocks%nb))
    allocate(dhand%ob2dof_l(dhand%blocks%nb))
    allocate(dhand%dof2ob_l(dhand%blocks%nb))
    call memalloc(dhand%blocks%nb,dhand%ndofs,__FILE__,__LINE__)

  end subroutine dof_handler_create_multiphysics

end module dof_handler_class

