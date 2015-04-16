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
!***********************************************************************
! All allocatable arrays
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc
# include "debug.i90"
!***********************************************************************
module volume_integration_tools_names
  use types
  use memor
#ifdef memcheck
  use iso_c_binding
#endif
  use quadrature_names
  use interpolation_names
  use femap_names
  use fem_space_types
  implicit none
  private

  type volume_integrator
     type(quadrature)         :: quad        ! Quadrature rules for elements
     type(interpolation)      :: uint_ref    ! Unknown interpolation in the reference element domain
     type(interpolation)      :: uint_phy    ! Unknown interpolation in the physical element domain
     type(interpolation)      :: gint_ref    ! Geometry interpolation in the reference element domain
     type(interpolation)      :: gint_phy    ! Geometry interpolation in the physical element domain
     type(femap)              :: femap       ! FE mapping
  end type volume_integrator
  type volume_integrator_pointer
     type(volume_integrator)          , pointer :: p => NULL() 
  end type volume_integrator_pointer

  public :: volume_integrator, volume_integrator_pointer

# define var_type type(volume_integrator_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

  public :: volume_integrator_create, volume_integrator_free, set_integ


contains

  ! =================================================================================================
# include "mem_body.i90"

  !==================================================================================================
  subroutine volume_integrator_create(gtype,utype,ndime,g_ord,u_ord,integ,khie,mnode)
    implicit none
    ! Parameters
    integer(ip)            , intent(in)  :: gtype, utype
    integer(ip)            , intent(in)  :: ndime, g_ord, u_ord
    type(volume_integrator), intent(out) :: integ
    logical(lg),   optional, intent(in)  :: khie
    integer(ip),   optional, intent(in)  :: mnode

    ! Local variables
    integer(ip) :: ngaus,lrule,llapl,gnode,unode,nlocs,i
    real(rp)    :: auxpo(max_nnode)
    logical(lg) :: khier

    ! Store identifiers of element
    !integ%ltype(1) = gtype
    !integ%ltype(2) = utype

    ! Hesssian interpolation
    if(.not.present(khie)) then    
       khier = .false. 
    else
       khier = khie
    end if

    ! Choose integration rules and create them
    call set_integ(utype,ndime,g_ord,u_ord,ngaus,nlocs,lrule,llapl,gnode,unode,mnode)
    call quadrature_create(lrule,ndime,ngaus,integ%quad)

    ! Create unknown interpolations
    call interpolation_create(1,1,    1,ndime,unode,ngaus,integ%uint_ref,khier)
    call interpolation_create(1,1,llapl,ndime,unode,ngaus,integ%uint_phy,khier)

    ! Create geometry interpolations
    call interpolation_create(1,1,    1,ndime,gnode,ngaus,integ%gint_ref,khier)
    call interpolation_create(1,1,llapl,ndime,gnode,ngaus,integ%gint_phy,khier)

    ! Define interpolations in the reference domain
    if (utype == Q_type_id) then
       do i=1,nlocs
          auxpo(i) = integ%quad%pos(1,i)
       end do
       call interpolation_local(auxpo(1:nlocs),integ%uint_ref,ndime,u_ord,nlocs)
       call interpolation_local(auxpo(1:nlocs),integ%gint_ref,ndime,g_ord,nlocs)
    elseif (utype == P_type_id) then
       call interpolation_local(integ%quad%pos,integ%uint_ref)
       call interpolation_local(integ%quad%pos,integ%gint_ref)
    else
       write(*,*) 'integration :: ERROR! Unknown element type.'
       check (1 == 0)
    end if

    ! Create fe map
    call femap_create(llapl,ndime,ngaus,integ%femap)

  end subroutine volume_integrator_create
  ! =================================================================================================
  subroutine volume_integrator_free(integ)
    implicit none
    ! Parameters
    type(volume_integrator), intent(inout) :: integ    
    
    ! Destruct quadratures
    call quadrature_free(integ%quad)

    ! Destruct interpolations
    call interpolation_free(integ%gint_ref)
    call interpolation_free(integ%gint_phy)
    call interpolation_free(integ%uint_ref)
    call interpolation_free(integ%uint_phy)

    ! Destruct fe map
    call femap_free(integ%femap)

    !integ%ltype = 0
  end subroutine volume_integrator_free
  !==================================================================================================
  subroutine set_integ(utype,ndime,g_ord,u_ord,ngaus,nlocs,lrule,llapl,gnode,unode,mnode)
    implicit none
    ! Parameters
    integer(ip),           intent(in)  :: utype,ndime,g_ord,u_ord
    integer(ip),           intent(out) :: ngaus,nlocs,lrule,llapl,gnode,unode
    integer(ip), optional, intent(in)  :: mnode

    integer(ip) ::  ggaus,glocs,grule,glapl

    if (utype == P_type_id) then
       call P_set_integ(ndime,u_ord,unode,ngaus,nlocs,lrule,llapl,mnode)
       call P_set_integ(ndime,g_ord,gnode,ggaus,glocs,grule,glapl,mnode)
    elseif (utype == Q_type_id) then
       call Q_set_integ(ndime,u_ord,unode,ngaus,nlocs,lrule,llapl,mnode)
       call Q_set_integ(ndime,g_ord,gnode,ggaus,glocs,grule,glapl,mnode)
    else
       write(*,*) __FILE__,__LINE__,'ERROR elem type not found!'
       check (1==0)
    end if
    
    ngaus = max(ngaus,ggaus)
    lrule = max(lrule,grule)
    llapl = max(llapl,glapl)
  end subroutine set_integ

end module volume_integration_tools_names

!***********************************************************************
!***********************************************************************
!***********************************************************************

module face_integration_tools_names
  use types
  use memor
#ifdef memcheck
  use iso_c_binding
#endif
  use quadrature_names
  use quadrature_faces
  use interpolation_names
  use face_interpolation_names
  use femap_names
  use bomap_names
  use fem_space_types
  use volume_integration_tools_names
  implicit none
  private

  type face_integrator 
     integer(ip)              :: ltype        ! List of combinations of face types
     type(quadrature)         :: quad         ! Quadrature rules for boundary faces
     type(face_quadrature)    :: fquad        ! Quadrature points on the faces of the element
     type(interpolation)      :: gint_ref     ! Geometry interpolation in the reference face domain
     type(bomap)              :: bomap        ! Boundary face mapping
     type(femap)              :: femap        ! FE mapping
     type(face_interpolation) :: ufint_ref    ! Unknown interpolation in the reference face domain
     type(face_interpolation) :: ufint_phy    ! Unknown interpolation in the physical face domain
     type(face_interpolation) :: gfint_ref    ! Geometry interpolation in the reference face domain
     type(face_interpolation) :: gfint_phy    ! Geometry interpolation in the physical face domain
  end type face_integrator
  type face_integrator_pointer
     type(face_integrator)          , pointer :: p => NULL() 
  end type face_integrator_pointer
  type element_face_integrator
     type(face_integrator_pointer), allocatable :: p(:)  ! Pointer to face integration
  end type element_face_integrator

  public :: face_integrator, face_integrator_pointer, element_face_integrator

# define var_type type(face_integrator_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

  public :: face_integrator_create, face_integrator_free, set_face_type

contains

  ! =================================================================================================
# include "mem_body.i90"

  ! =================================================================================================
  subroutine face_integrator_create(gfinf,ufinf,nd,integ,subface_ids_)
    implicit none
    ! Parameters
    type(fem_fixed_info_pointer), intent(in)    :: gfinf,ufinf
    integer(ip)                 , intent(in)    :: nd
    type(face_integrator)            , intent(inout) :: integ
    integer(ip) ,      optional , intent(in)    :: subface_ids_

    ! Local variables
    integer(ip)              :: gnode, unode, gfnod, ufnod
    integer(ip)              :: nface, utype, ngaux,nlocx,nlocs
    integer(ip)              :: mnode,elnod,i,j, refinement_level
    integer(ip)              :: ng,lrule,llapl,subface_ids
    logical                  :: istat
    real(rp)                 :: auxpo(max_order+1)

       call set_face_type    ( gfinf%p%ftype, gfinf%p%ftype,nd,utype)
       gnode = gfinf%p%nnode
       unode = ufinf%p%nnode
       nface = gfinf%p%nobje_dim(nd+1) -  gfinf%p%nobje_dim(nd)
       call set_integ (utype,nd-1,gfinf%p%order,ufinf%p%order,ng,nlocs,lrule,llapl,           &
            &          gfnod,ufnod)
       ! Store identifiers of element
       integ%ltype = gfinf%p%ftype

       ! Maximum amount of nodes in faces and elements
       mnode = max(gfnod,ufnod)
    
       call quadrature_create(lrule,nd-1,ng,integ%quad) 
       ! Create bf map
       call bomap_create(nd,ng,integ%bomap)
       ! Allocate interpolation
       call interpolation_create(1,1,    1,nd-1,gfnod,ng,integ%gint_ref)
       ! Define interpolations in the reference domain
       call interpolation_local(integ%quad%pos,integ%gint_ref)

       ! Quadrature for the faces of element i
       call face_quadrature_create(nd,nface,integ%quad,gfinf%p,integ%fquad)
       
       ! Allocate interpolation for element i
       call face_interpolation_create (unode,ufnod,nd,ng,nlocs,nface,integ%ufint_ref)
       call face_interpolation_create (unode,ufnod,nd,ng,nlocs,       1,integ%ufint_phy)
       call face_interpolation_create (gnode,gfnod,nd,ng,nlocs,nface,integ%gfint_ref)
       call face_interpolation_create (gnode,gfnod,nd,ng,nlocs,       1,integ%gfint_phy)
       ! Create reference interpolations for element i
       if     (ufinf%p%ftype == Q_type_id) then
          do j=1,nlocs
             auxpo(j) = integ%quad%pos(1,j)
          end do
          call face_interpolation_local_Q (auxpo,integ%ufint_ref,nd,ufinf%p%order,nlocs)
       elseif (ufinf%p%ftype == P_type_id) then
          call face_interpolation_local  (integ%fquad%pos,integ%ufint_ref)
       else
          write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
       end if
       if     (gfinf%p%ftype == Q_type_id) then
          do j=1,nlocs
             auxpo(j) = integ%quad%pos(1,j)
          end do
          call face_interpolation_local_Q (auxpo,integ%gfint_ref,nd,gfinf%p%order,nlocs)
       elseif (gfinf%p%ftype == P_type_id) then
          call face_interpolation_local  (integ%fquad%pos,integ%gfint_ref)
       else
          write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
       end if
       !call face_interpolation_local  (integ%fquad%pos,integ%ufint_ref)
       !call face_interpolation_local  (integ%fquad%pos,integ%gfint_ref)
       ! Create FE map
       call femap_create(llapl,nd,ng,integ%femap)

  end subroutine face_integrator_create
  ! ================================================================================================
  subroutine face_integrator_free(integ)
    implicit none
    ! Parameters
    type(face_integrator), intent(inout) :: integ

    ! Destruct quadratures
    call quadrature_free     (integ%quad)
    call face_quadrature_free(integ%fquad)

    ! Destruct fe map
    call bomap_free(integ%bomap)
    call femap_free(integ%femap)
    call interpolation_free(integ%gint_ref)

    ! Deallocate Face interpolations
    call face_interpolation_free(integ%ufint_ref)
    call face_interpolation_free(integ%gfint_ref)
    call face_interpolation_free(integ%ufint_phy)
    call face_interpolation_free(integ%gfint_phy)

  end subroutine face_integrator_free

 !==================================================================================================
  ! This routine gives the type of the face given the type of its surrounding elements
  subroutine set_face_type( type1, type2, ndime, ftype)
    implicit none
    ! Parameters
    integer(ip) , intent(in)  :: type1, type2
    integer(ip) , intent(in)  :: ndime
    integer(ip) , intent(out) :: ftype

    ! Local variables
    integer(ip)               :: ftype2

    select case (type1)

    case (P_type_id) 
       if (ndime == 1) then
          ftype = Q_type_id
       else
          ftype = P_type_id
       end if

    case (Q_type_id) 
       ftype = Q_type_id

    case default
       write(*,*) 'integration:: ERROR! Element type not recognised.'
       check (1 == 0)
    end select

    select case (type2)

    case (P_type_id) 
       if (ndime == 1) then
          ftype2 = Q_type_id
       else
          ftype2 = P_type_id
       end if

    case (Q_type_id) 
       ftype2 = Q_type_id

    case default
       write(*,*) 'integration:: ERROR! Element type not recognised.'
       check (1 == 0)
    end select

    check (ftype == ftype2)
    
  end subroutine set_face_type

end module face_integration_tools_names

!***********************************************************************
!***********************************************************************
!***********************************************************************

module integration_tools_names
  use volume_integration_tools_names
  use face_integration_tools_names
  implicit none
end module integration_tools_names
