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
module mesh_names
  use types_names
  use memor_names
  use materials_names
  use stdio_names
  use iso_fortran_env, only : output_unit
  !use conditions_names
# include "debug.i90"
  implicit none
  private

  type mesh_t
     integer(ip)                :: &
          nelem,                   &         ! Number of elements
          nnode                              ! Maximum number of nodes per element

     integer(ip), allocatable ::   &
          pnods(:) ,               &         ! pointers to the lnods
          lnods(:)                           ! list of nodes of each element

     integer(ip)                :: &
          ndime,                   &         ! Number of space dimensions
          npoin                              ! Number of nodes

     real(rp), allocatable ::      &
          coord(:,:)                         ! Node coordinates       
  end type mesh_t

  ! Types
  public :: mesh_t

  ! Functions
  public :: mesh_to_dual, mesh_free
  public :: mesh_print

  ! Constants
  public :: count_elements_around_points, list_elements_around_points

contains
  !=============================================================================
  ! Methods defined for mesh are:
  !
  ! - mesh_to_dual(primal_mesh,dual_mesh)
  !   IN: primal_mesh
  !   OUT: dual_mesh
  !
  ! - mesh_dealloc(mesh)
  !   INOUT: mesh
  !
  !-----------------------------------------------------------------------------
  !
  ! These methods call auxiliar routines
  !
  ! -  subroutine count_elemental_graph()
  ! -  subroutine list_elemental_graph()
  ! 
  ! The arguments of these routine are generated in the methods unpacking mesh 
  ! and graph objects.
  !
  !=============================================================================
  subroutine mesh_to_dual(primal_mesh,dual_mesh)
    !-----------------------------------------------------------------------
    ! This routine generates the dual mesh (list of elements around
    ! points) of a given primal mesh. The dual mesh allways has nelty/=1
    ! and, as it is not used, we set it to 0.
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_t), intent(in)     :: primal_mesh
    type(mesh_t), intent(out)    :: dual_mesh
    
    dual_mesh%ndime=primal_mesh%ndime
    dual_mesh%npoin=primal_mesh%nelem
    dual_mesh%nelem=primal_mesh%npoin

    ! Count elements around points and generate pointers to store them (pelpo)
    call memalloc (dual_mesh%nelem+1, dual_mesh%pnods, __FILE__,__LINE__)
    
    call count_elements_around_points( primal_mesh%npoin,primal_mesh%nelem,primal_mesh%pnods, &
         &   primal_mesh%lnods,primal_mesh%nnode,dual_mesh%nnode,dual_mesh%pnods)
    
    ! List elements around points (lelpo)
    call memalloc (dual_mesh%pnods(dual_mesh%nelem+1)-1,dual_mesh%lnods,__FILE__,__LINE__)
    
    call list_elements_around_points( primal_mesh%npoin,primal_mesh%nelem,primal_mesh%pnods, &
         &   primal_mesh%lnods,primal_mesh%nnode,dual_mesh%nnode,dual_mesh%pnods,  &
         &   dual_mesh%lnods)
    
  end subroutine mesh_to_dual

  ! QUESTION (Pavel Kus): Is there a reason, why to two following subrutins, primal_mesh is
  ! not provided directly? It would simplify the interface... 
  
  !==============================================================================
  subroutine count_elements_around_points(npoin,nelem,pnods,lnods, &
       &                                  nnode,nelpo,pelpo)
    !-----------------------------------------------------------------------
    ! This routine counts the number of elements around mesh points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: npoin,nelem
    integer(ip), intent(in)  :: pnods(nelem+1), lnods(pnods(nelem+1)-1)
    integer(ip), intent(in)  :: nnode          ! Max. number of nodes per element 
    integer(ip), intent(out) :: nelpo          ! Max. number of elements around a point
    integer(ip), intent(out) :: pelpo(npoin+1) ! Number of elements around each point

    ! Local variables
    integer(ip)              :: inode, ipoin, size_lnods

    size_lnods = pnods(nelem+1)-1 
    
    ! Compute the number of elements around each point
    pelpo=0
    do inode=1, size_lnods
       ipoin=lnods(inode)
       pelpo(ipoin+1)=pelpo(ipoin+1)+1
    end do
    
    ! Find the maximum number of elements around a point
    nelpo=0
    do ipoin=1,npoin
       nelpo=max(nelpo,pelpo(ipoin+1))
    end do
    
    ! Compute pointers to the starting position of the list
    ! of elements around each point
    pelpo(1)=1
    do ipoin=1,npoin
       pelpo(ipoin+1)=pelpo(ipoin+1)+pelpo(ipoin)
    end do

  end subroutine count_elements_around_points

  !=============================================================================
  subroutine list_elements_around_points(npoin,nelem,pnods,lnods, &
       &                                   nnode,nelpo,pelpo,lelpo)
    !-----------------------------------------------------------------------
    ! This routine lists the number of elements around mesh points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: npoin,nelem
    integer(ip), intent(in)    :: pnods(*),lnods(*) ! Assumed size arrays ... 
    integer(ip), intent(in)    :: nnode           ! Max. number of nodes per element
    integer(ip), intent(in)    :: nelpo           ! Max. number of elements around a point
    integer(ip), intent(inout) :: pelpo(npoin+1)  ! Number of elements around each point
    integer(ip), intent(out)   :: lelpo(pelpo(npoin+1)) ! List of elements around points

    ! Local variables 
    integer(ip)              :: ielem, inode, ipoin

       
    ! Compute the list of elements around each point.
    ! pelpo is used instead of auxiliary work space.
    
    do ielem=1,nelem 
       do inode=pnods(ielem),pnods(ielem+1)-1 
          ipoin=lnods(inode)
          lelpo(pelpo(ipoin))=ielem
          pelpo(ipoin)=pelpo(ipoin)+1
       end do
    end do
    
    
    ! Recover pelpo
    do ipoin=npoin+1, 2, -1
       pelpo(ipoin)=pelpo(ipoin-1)
    end do
    pelpo(1) = 1

    ! write(*,'(a)')      'pelpo='                     ! DBG:
    ! write(*,'(10i10)')  pelpo(1:npoin+1)             ! DBG:

    ! write(*,'(a)')      'lelpo='                     ! DBG:
    ! write(*,'(10i10)')  lelpo(1:(pelpo(npoin+1)-1))  ! DBG:

  end subroutine list_elements_around_points

  !=============================================================================
  subroutine mesh_free (msh)
    !-----------------------------------------------------------------------
    ! This routine generates deallocates a mesh
    !-----------------------------------------------------------------------
    implicit none
    type(mesh_t), intent(inout)  :: msh

    call memfree (msh%pnods,__FILE__,__LINE__)
    call memfree (msh%lnods,__FILE__,__LINE__)

    ! This data structure can hold (depending on the context):
    ! (1) Permanent geometrical mesh (e.g., mesh read from GiD)
    ! (2) Temporary topological mesh (resulting from generate_vefs_mesh_conditions)
    ! (2) Temporary Dual mesh
    ! (3) Coarse-grid mesh in MLBDDC
    ! In the case of (2), (3), and (4) coord is not allocated!!!
    if (allocated(msh%coord)) call memfree (msh%coord,__FILE__,__LINE__)

    msh%ndime=0
    msh%npoin=0
    msh%nelem=0
    msh%nnode=0

  end subroutine mesh_free

  
  subroutine mesh_print(mesh, ounit)
     implicit none
     type(mesh_t), intent(in)          :: mesh
     integer(ip), intent(in), optional   :: ounit
     integer(ip)                         :: idx, num_nodes, ounit_used
     character(len=100)                  :: fmt

     if(present(ounit)) then
        ounit_used = ounit
     else
        ounit_used = output_unit
     end if
     
     do idx = 1, mesh%nelem
        num_nodes = mesh%pnods(idx+1) - mesh%pnods(idx)
        fmt = "(I4, A, "//trim(ch(num_nodes))//"I4)"
        write(ounit_used, fmt) idx, " -> ", mesh%lnods(mesh%pnods(idx) : mesh%pnods(idx+1)-1)
     end do

   end subroutine mesh_print

end module mesh_names
