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
module fem_mesh_names
  use types
  use memor
  use fem_materials_names
  use stdio
  use iso_fortran_env, only : output_unit
  !use fem_conditions_names
# include "debug.i90"
  implicit none
  private

  integer(ip), parameter :: structured = 100
  integer(ip), parameter :: unstructured = 200

  type fem_mesh
     integer(ip)                :: &
          mtype=unstructured,      &         ! Mesh type
          nedir(3),                &         ! Number of elements in each dir (if structured)
          isper(3)                           ! If direction is periodic

     integer(ip)                :: &
          nelty,                   &         ! Number of element types
          nelem,                   &         ! Number of elements
          nboun = 0,               &         ! Number of boundary faces
          nnode,                   &         ! Maximum number of nodes per element
          nnodb                              ! Maximum number of nodes per boundary faces

     integer(ip), allocatable ::   &
          pnods(:) ,               &         ! pointers to the lnods
          lnods(:)                           ! list of nodes of each element

     integer(ip)                :: &
          ndime,                   &         ! Number of space dimensions
          npoin                              ! Number of nodes

     real(rp), allocatable ::      &
          coord(:,:)                         ! Node coordinates

     integer(ip), allocatable ::   &
!          lface_new(:) ,          &         ! list of faces          
          p_face(:),               &         ! start of faces of element (2,nelem)
          pboun(:) ,               &         ! pointers to the boundary faces
          lboun(:) ,               &         ! list of nodes of each boundary face
          lboel(:,:)                         ! information about host volume element

     integer(ip), allocatable ::   &
          pface(:) ,               &         ! pointers to the internal faces
          lface(:)                           ! elements of each face (and their local faces)

     !SHOCK CAPTURING PARAMETERS (may be useful for DG implementation)
     !-----------------------------------------------------------------
     !TODO: Redefine nface,nodfac and permu to take under consideration
     !      meshes with different types of elements. 
     ! bou_tface=nboun but the value is not assigned in nboun because
     ! nboun>0 is used as a criterion to know if pboun and lboun have
     ! been allocated.
     !-----------------------------------------------------------------
     integer(ip)              ::   &
          nelpo = 0                          ! Maximum amount of elements per node

     integer(ip), allocatable ::   &
          pelpo(:) ,               &         ! pointers to the lelpo
          lelpo(:)                           ! elements per node

     integer(ip), allocatable ::   &
          facel(:,:),              &         ! Elements separated by the faces
          permu(:,:)                         ! Permutations of the faces

     integer(ip)                :: &
          nface = 0,               &         ! Number of faces per element
          nodfac= 0,               &         ! Number of nodes per face
          tface = 0,               &         ! Total amount of interior faces
          bou_tface = 0,           &         ! Number of boundary faces     
          shaface = 0                        ! Number of faces in the shared boundary of a subdomain 

     type(fem_materials)      ::  prob       ! Problem to be solved        
     
  end type fem_mesh

  interface fem_mesh_alloc
     module procedure unstructured_mesh_alloc, structured_mesh_alloc
  end interface fem_mesh_alloc

  ! Types
  public :: fem_mesh

  ! Functions
  public :: mesh_to_dual, fem_mesh_alloc, fem_mesh_free
  public :: print_fem_mesh

  ! Constants
  public :: structured, unstructured, count_elements_around_points, list_elements_around_points

contains
  !=============================================================================
  ! Methods defined for fem_mesh are:
  !
  ! - mesh_to_dual(primal_mesh,dual_mesh)
  !   IN: primal_mesh
  !   OUT: dual_mesh
  !
  ! - fem_mesh_dealloc(mesh)
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
    type(fem_mesh), intent(in)     :: primal_mesh
    type(fem_mesh), intent(out)    :: dual_mesh
    
    dual_mesh%nelty=0
    dual_mesh%ndime=primal_mesh%ndime
    dual_mesh%npoin=primal_mesh%nelem
    dual_mesh%nelem=primal_mesh%npoin
    dual_mesh%nelpo=0
    dual_mesh%nboun=0

    ! Count elements around points and generate pointers to store them (pelpo)
    call memalloc (dual_mesh%nelem+1, dual_mesh%pnods, __FILE__,__LINE__)
    
    call count_elements_around_points( &
         &   primal_mesh%nelty,primal_mesh%npoin,primal_mesh%nelem,primal_mesh%pnods, &
         &   primal_mesh%lnods,primal_mesh%nnode,dual_mesh%nnode,dual_mesh%pnods)
    
    ! List elements around points (lelpo)
    call memalloc (dual_mesh%pnods(dual_mesh%nelem+1), dual_mesh%lnods,            __FILE__,__LINE__)
    
    call list_elements_around_points( &
         &   primal_mesh%nelty,primal_mesh%npoin,primal_mesh%nelem,primal_mesh%pnods, &
         &   primal_mesh%lnods,primal_mesh%nnode,dual_mesh%nnode,dual_mesh%pnods,  &
         &   dual_mesh%lnods)
    
    return

  end subroutine mesh_to_dual

  ! QUESTION (Pavel Kus): Is there a reason, why to two following subrutins, primal_mesh is
  ! not provided directly? It would simplify the interface... 
  
  !==============================================================================
  subroutine count_elements_around_points(nelty,npoin,nelem,pnods,lnods, &
       &                                    nnode,nelpo,pelpo)
    !-----------------------------------------------------------------------
    ! This routine counts the number of elements around mesh points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nelty,npoin,nelem
    integer(ip), intent(in)  :: pnods(*), lnods(*)
    integer(ip), intent(in)  :: nnode          ! Max. number of nodes per element 
    integer(ip), intent(out) :: nelpo          ! Max. number of elements around a point
    integer(ip), intent(out) :: pelpo(npoin+1) ! Number of elements around each point

    ! Local variables
    integer(ip)              :: ielem, inode, ipoin, size_lnods, i_cing_node, inode_proper, ipoin_proper
    logical                  :: should_be_added

    
    ! implementation with hanging nodes contains more logic in the inner loops

       ! original implementation without hanging nodes
       if(nelty==1 ) then
          size_lnods = nnode*nelem
       else
          size_lnods = pnods(nelem+1)-1 
       end if

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
  subroutine list_elements_around_points(nelty,npoin,nelem,pnods,lnods, &
       &                                   nnode,nelpo,pelpo,lelpo)
    !-----------------------------------------------------------------------
    ! This routine lists the number of elements around mesh points
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: nelty,npoin,nelem
    integer(ip), intent(in)    :: pnods(*),lnods(*) ! Assumed size arrays ... 
    integer(ip), intent(in)    :: nnode           ! Max. number of nodes per element
    integer(ip), intent(in)    :: nelpo           ! Max. number of elements around a point
    integer(ip), intent(inout) :: pelpo(npoin+1)  ! Number of elements around each point
    integer(ip), intent(out)   :: lelpo(pelpo(npoin+1)) ! List of elements around points

    ! Local variables 
    integer(ip)              :: ielem, inode, ipoin, i_cing_node, inode_proper, ipoin_proper
    logical                  :: should_be_added  
    ! write(*,'(a)')      'pelpo='                     ! DBG:
    ! write(*,'(10i10)')  pelpo(1:npoin+1)             ! DBG:

       ! regular case
       
       ! Compute the list of elements around each point.
       ! pelpo is used instead of auxiliary work space.
       if(nelty==1) then
          do ielem=1,nelem
             do inode=1,nnode
                ipoin=lnods((ielem-1)*nnode+inode)
                lelpo(pelpo(ipoin))=ielem
                pelpo(ipoin)=pelpo(ipoin)+1
             end do
          end do
       else
          do ielem=1,nelem 
             do inode=pnods(ielem),pnods(ielem+1)-1 
                ipoin=lnods(inode)
                lelpo(pelpo(ipoin))=ielem
                pelpo(ipoin)=pelpo(ipoin)+1
             end do
          end do
       end if
       
    
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
  subroutine fem_mesh_free (msh)
    !-----------------------------------------------------------------------
    ! This routine generates deallocates a fem_mesh
    !-----------------------------------------------------------------------
    implicit none
    type(fem_mesh), intent(inout)  :: msh

    ! if(msh%nelty/=1) 
    call memfree (msh%pnods,__FILE__,__LINE__)
    call memfree (msh%lnods,__FILE__,__LINE__)
    if(allocated(msh%coord)) call memfree (msh%coord,__FILE__,__LINE__)
    if(allocated(msh%p_face)) call memfree (msh%p_face,__FILE__,__LINE__)
    
    if (msh%nboun>0) then
!       call memfree (msh%pboun,__FILE__,__LINE__)
!       call memfree (msh%lboun,__FILE__,__LINE__)
!       call memfree (msh%lboel,__FILE__,__LINE__)
    end if

    if (msh%nelpo>0) then
       call memfree (msh%pelpo,__FILE__,__LINE__)
       call memfree (msh%lelpo,__FILE__,__LINE__)
    end if

    msh%ndime=2
    msh%nelty=1
    msh%npoin=0
    msh%nelem=0
    msh%nnode=0
    msh%nelpo=0    
    msh%nboun=0    
    msh%nnodb=0    

    return

  end subroutine fem_mesh_free
  !=============================================================================
  subroutine unstructured_mesh_alloc (ndime,nelty,nnode,npoin,nelem,msh,wcoor)
    !-----------------------------------------------------------------------
    ! This routine allocates a fem_mesh
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)   , intent(in)  :: ndime,nelty,nnode,npoin,nelem
    type(fem_mesh), intent(out) :: msh
    integer(ip)   , optional, intent(in)  :: wcoor

    msh%ndime=ndime
    msh%nelty=nelty
    msh%nnode=nnode
    msh%npoin=npoin
    msh%nelem=nelem
    msh%nboun=0
    msh%nelpo=0

    call memalloc (nelem + 1,msh%pnods, __FILE__,__LINE__)
    call memalloc (nnode*nelem,msh%lnods, __FILE__,__LINE__)
    if(.not.present(wcoor)) call memalloc (ndime,npoin,msh%coord, __FILE__,__LINE__)

    return

  end subroutine unstructured_mesh_alloc

  !=============================================================================
  subroutine structured_mesh_alloc (ndime,nnode,nedir,isper,msh,wcoor)
    !-----------------------------------------------------------------------
    ! This routine allocates a structured mesh
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)   , intent(in)  :: ndime,nnode,nedir(3),isper(3)
    type(fem_mesh), intent(out) :: msh
    integer(ip)   , optional, intent(in)  :: wcoor

    integer(ip) :: npoix,npoiy,npoiz,i

    msh%mtype=structured
    msh%nelty=1
    msh%ndime=ndime
    msh%nnode=nnode
    msh%nedir=nedir
    msh%isper=isper
    msh%nboun=0
    msh%nelpo=0

    if(msh%ndime==2) then

       if(msh%nnode==3.or.msh%nnode==4) then ! P1 or Q1
          npoix=msh%nedir(1)+1
          npoiy=msh%nedir(2)+1
       else if(msh%nnode==6.or.msh%nnode==9) then ! P2 or Q2
          npoix=2*msh%nedir(1)+1
          npoiy=2*msh%nedir(2)+1
       else if(msh%nnode==10.or.msh%nnode==16) then ! P3 or Q3
          npoix=3*msh%nedir(1)+1
          npoiy=3*msh%nedir(2)+1
       end if
       if(isper(1)==1) npoix=npoix-1
       if(isper(2)==1) npoiy=npoiy-1
       msh%npoin=npoix*npoiy

       ! Number of elements
       if(msh%nnode==3.or.msh%nnode==6.or.msh%nnode==10) then ! P1, P2 or P3
          msh%nelem=2*nedir(1)*nedir(2)
       else if(msh%nnode==4.or.msh%nnode==9.or.msh%nnode==16) then ! Q1, Q2 or Q3
          msh%nelem=nedir(1)*nedir(2)
       end if

    else if(msh%ndime==3) then

       if(msh%nnode==8) then ! Q1
          npoix=msh%nedir(1)+1
          npoiy=msh%nedir(2)+1
          npoiz=msh%nedir(3)+1
       end if
       if(isper(1)==1) npoix=npoix-1
       if(isper(2)==1) npoiy=npoiy-1
       if(isper(3)==1) npoiz=npoiz-1
       msh%npoin=npoix*npoiy*npoiz

       ! Number of elements
       if(msh%nnode==8) then ! Q1
          msh%nelem=nedir(1)*nedir(2)*nedir(3)
       end if

    end if

    !call memalloc (1,msh%pnods, __FILE__,__LINE__)
    call memalloc (msh%nelem+1,msh%pnods, __FILE__,__LINE__)

    call memalloc (msh%nnode*msh%nelem,msh%lnods, __FILE__,__LINE__)
    if(.not.present(wcoor))  call memalloc (msh%ndime,msh%npoin,msh%coord, __FILE__,__LINE__)

    !msh%pnods(1) = msh%nnode
    msh%pnods(1) = 1
    do i = 1,msh%nelem
       msh%pnods(i+1) = msh%pnods(i) + msh%nnode
    end do

  end subroutine structured_mesh_alloc
  
  
  subroutine print_fem_mesh(mesh, ounit)
     implicit none
     type(fem_mesh), intent(in)          :: mesh
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

  end subroutine print_fem_mesh

end module fem_mesh_names
