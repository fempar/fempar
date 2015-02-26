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
  use sort_class
  use fem_mesh_class
  use fem_partition_class
  use fem_graph_class
  use fem_conditions_class
  use fem_materials_class
  use fem_space_class
  use fem_space_types
  use maps_class
  use mesh_graph
  !use fem_blocks_class
  use array_class
  
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

     type(array_ip2), allocatable :: &
          ob2dof_l(:)                         ! (nobje,2) : provides the dofs x object and associated variable

     type(array_ip1), allocatable :: &
          ob2dof_p(:),               &       ! pointer to ob2dof_l
          dof2ob_l(:),               &       ! provides the object associated at each dof
          blknvarsxprob(:),          &       ! Number of vars per problem vars_prob(nprob)
          i_varsxprob(:),            &       ! Physical vars associated to a problem (ia)
          j_varsxprob(:)!,           &       ! (ja)

     !type(fem_blocks)  ::           &
     !     blocks                             ! blocks ordering

     !type(dofh_elem), allocatable :: &
     !     elem2dof(:)                       ! Map from elem to dof

     !logical(lg), allocatable     ::     &
     !     continuity(:)                      ! Continuity of dofs on a given node continuity (npoin)
     
     ! this pointer is at various places: here, in fem_space and in mesh
     ! it is mainly for convenience, since it is needed in multiple places in the code
     ! todo: maybe it should be rewised

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
    
  end subroutine elem2dof_create

 !=============================================================================
  subroutine jvars_create( dhand, femsp, mesh )
    implicit none
    type(dof_handler),   intent(in)    :: dhand
    type(fem_space),     intent(inout) :: femsp
    type(fem_mesh),      intent(in)    :: mesh

    
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


  end subroutine object2dof_list

  !=============================================================================
  subroutine dof_handler_print(lunou, d)
    implicit none
    integer(ip)      , intent(in) :: lunou
    type(dof_handler), intent(in) :: d

  end subroutine dof_handler_print

  !=============================================================================
  subroutine dof_handler_free(d)
    implicit none
    type(dof_handler), intent(inout) :: d

  end subroutine dof_handler_free

  !=============================================================================
  subroutine dof_handler_create_one_physics ( nvars, code, cont, dhand)!, blks) 
    implicit none
    integer(ip)      , intent(in)           :: nvars, code, cont
    type(dof_handler), intent(inout)          :: dhand
    !type(fem_blocks) , optional, intent(in) :: blks

  end subroutine dof_handler_create_one_physics

  !=============================================================================
  subroutine dof_handler_create_multiphysics ( nvars, nprob, prob_code, prob_nunk, prob_varini, &
       &                                       nunkxblock, cont, dhand)!, blks) 
    implicit none
    integer(ip)         , intent(in)        :: nvars, nprob, prob_code(nprob), prob_nunk(nprob)
    integer(ip)         , intent(in)        :: prob_varini(:,:), cont, nunkxblock(:,:)
    type(dof_handler)   , intent(inout)       :: dhand
!    type(fem_blocks) , optional, intent(in) :: blks

  end subroutine dof_handler_create_multiphysics

end module dof_handler_class

