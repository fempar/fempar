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
module volumne_integration_names
  use types
  use memor
#ifdef memcheck
  use iso_c_binding
#endif
  use quadrature_names
  use interpolation_names
  use femap_names
  implicit none
  private

  type vol_integ
     integer(ip)              :: ltype(2)    ! Tags to identify element
     type(quadrature)         :: quad        ! Quadrature rules for elements
     type(interpolation)      :: uint_ref    ! Unknown interpolation in the reference element domain
     type(interpolation)      :: uint_phy    ! Unknown interpolation in the physical element domain
     type(interpolation)      :: gint_ref    ! Geometry interpolation in the reference element domain
     type(interpolation)      :: gint_phy    ! Geometry interpolation in the physical element domain
     type(femap)              :: femap       ! FE mapping
  end type vol_integ
  type vol_integ_pointer
     type(vol_integ)          , pointer :: p => NULL() 
  end type vol_integ_pointer
  public :: vol_integ, vol_integ_pointer

# define var_type type(vol_integ_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module volumne_integration_names

!***********************************************************************
module face_integration_names
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
  implicit none
  private

  type face_integ 
     integer(ip)              :: ltype(2)     ! List of combinations of face types
     type(quadrature)         :: quad         ! Quadrature rules for boundary faces
     type(face_quadrature)    :: fquad(2)     ! Quadrature points on the faces of the element
     type(interpolation)      :: gint_ref     ! Geometry interpolation in the reference face domain
     type(bomap)              :: bomap        ! Boundary face mapping
     type(femap)              :: femap(2)     ! FE mapping
     type(face_interpolation) :: ufint_ref(2) ! Unknown interpolation in the reference face domain
     type(face_interpolation) :: ufint_phy(2) ! Unknown interpolation in the physical face domain
     type(face_interpolation) :: gfint_ref(2) ! Geometry interpolation in the reference face domain
     type(face_interpolation) :: gfint_phy(2) ! Geometry interpolation in the physical face domain
  end type face_integ
  type face_integ_pointer
     type(face_integ)          , pointer :: p => NULL() 
  end type face_integ_pointer
  public :: face_integ, face_integ_pointer

# define var_type type(face_integ_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"

end module face_integration_names
!***********************************************************************

module integration_names
  use types
  use memor
  use array_names
  use fem_mesh_names 
  use quadrature_names
  use quadrature_faces
  use interpolation_names
  use face_interpolation_names
  use femap_names
  use femap_interp
  use bomap_names
  use bomap_interp
  use fem_space_types
  use element_gather_tools
  use volumne_integration_names
  use face_integration_names
  implicit none
  private

  interface integ_create
     module procedure face_integ_create,vol_integ_create
  end interface integ_create
  interface integ_free
     module procedure vol_integ_free, face_integ_free
  end interface integ_free

  ! Functions
  public :: integ_create, integ_element, integ_free, integ_faces 
  public :: vol_integ, vol_integ_pointer
  public :: face_integ, face_integ_pointer
  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

  !==================================================================================================
  subroutine vol_integ_create(gtype,utype,ndime,g_ord,u_ord,integ,khie,mnode)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)    :: gtype, utype
    integer(ip)           , intent(in)    :: ndime, g_ord, u_ord
    type(vol_integ)       , intent(inout) :: integ
    integer(ip),  optional, intent(in)    :: khie,mnode       

    ! Local variables
    integer(ip) :: ngaus,lrule,llapl,khier,gnode,unode,nlocs,i
    real(rp)    :: auxpo(max_nnode)

    ! Store identifiers of element
    integ%ltype(1) = gtype
    integ%ltype(2) = utype

    ! Hesssian interpolation
    if(.not.present(khie)) then    
       khier = 0 
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

  end subroutine vol_integ_create

  !==================================================================================================
  subroutine face_integ_create(gfinf,ufinf,nd,integ,subface_ids_)
    implicit none
    ! Parameters
    type(fem_fixed_info_pointer), intent(in)    :: gfinf(2),ufinf(2)
    integer(ip)                 , intent(in)    :: nd
    type(face_integ)            , intent(inout) :: integ
    integer(ip) ,      optional , intent(in)    :: subface_ids_(2)

    ! Local variables
    integer(ip)              :: gnode(2), unode(2), gfnod(2), ufnod(2)
    integer(ip)              :: nface(2), utype, ngaux(2),nlocx(2),nlocs
    integer(ip)              :: mnode,elnod(2),i,j, refinement_level
    integer(ip)              :: ng,lrule,llapl,subface_ids(2)
    logical                  :: istat
    real(rp)                 :: auxpo(max_order+1)

    ! TO DO: Add refinement level as a variable
    if (present(subface_ids_)) then
       subface_ids=subface_ids_
       refinement_level = 0
       do i =1,2
          if (subface_ids(i) > 0) then
             refinement_level = 1
          end if
       end do
    else
       subface_ids = 0
       refinement_level =0
    end if

    assert (subface_ids(1) == 0)
    ! INTERIOR FACES
    if (associated(gfinf(2)%p)) then

       ! Choose integration rules and create them
       call set_face_type    ( gfinf(1)%p%ftype, gfinf(2)%p%ftype,nd,utype)
       do i = 1, 2
         gnode(i) = gfinf(i)%p%nnode
         unode(i) = ufinf(i)%p%nnode
         nface(i) = gfinf(i)%p%nobje_dim(nd+1) -  gfinf(i)%p%nobje_dim(nd)
         call set_integ (utype,nd-1,gfinf(i)%p%order,ufinf(i)%p%order,ngaux(i),nlocx(i),lrule,   &
              &          llapl,gfnod(i),ufnod(i))
         ! Store identifiers of element
         integ%ltype(i) = gfinf(i)%p%ftype
       end do
       ! Maximum amount of nodes in faces and elements
       mnode = max(maxval(gfnod),maxval(ufnod))
       ng = maxval(ngaux)
       nlocs = maxval(nlocx)
          
       !call set_integration  (nd-1,mnode,ng,lrule,llapl)
       call quadrature_create(lrule,nd-1,ng,integ%quad) 
       call bomap_create     (nd,ng,integ%bomap)
       ! Allocate interpolations
       call interpolation_create (1,1,1,nd-1,maxval(gfnod),ng,integ%gint_ref)
       ! Define interpolations in the reference domain
       call interpolation_local  (integ%quad%pos,integ%gint_ref)

       do i =1, 2
          ! Quadrature for the faces of element i
          call face_quadrature_create(nd,nface(i),integ%quad,gfinf(i)%p,integ%fquad(i),subface_ids(i))
          ! Allocate interpolation for element i
          call face_interpolation_create (unode(i),ufnod(i),nd,ng,nlocs,nface(i),integ%ufint_ref(i))
          call face_interpolation_create (unode(i),ufnod(i),nd,ng,nlocs,       1,integ%ufint_phy(i))
          call face_interpolation_create (gnode(i),gfnod(i),nd,ng,nlocs,nface(i),integ%gfint_ref(i))
          call face_interpolation_create (gnode(i),gfnod(i),nd,ng,nlocs,       1,integ%gfint_phy(i))
          ! Create reference interpolations for element i
          if     (ufinf(i)%p%ftype == Q_type_id) then
                do j=1,nlocs
                   auxpo(j) = integ%quad%pos(1,j)
                end do
             if (subface_ids(i) > 0) then
                !call face_interpolation_local  (integ%fquad(i)%pos,integ%ufint_ref(i))
                
                call subface_interpolation_local_Q (auxpo,integ%ufint_ref(i),nd,ufinf(i)%p%order,   &
                     &                           nlocs,refinement_level,subface_ids(i))
             else
                call face_interpolation_local_Q (auxpo,integ%ufint_ref(i),nd,ufinf(i)%p%order,nlocs)
             end if
          elseif (ufinf(i)%p%ftype == P_type_id) then
             call face_interpolation_local  (integ%fquad(i)%pos,integ%ufint_ref(i))
          else
             write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
          end if
          if     (gfinf(i)%p%ftype == Q_type_id) then
             do j=1,nlocs
                auxpo(j) = integ%quad%pos(1,j)
             end do
             call face_interpolation_local_Q (auxpo,integ%gfint_ref(i),nd,gfinf(i)%p%order,nlocs)
          elseif (gfinf(i)%p%ftype == P_type_id) then
             call face_interpolation_local  (integ%fquad(i)%pos,integ%gfint_ref(i))
          else
             write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
          end if
          ! Create FE map
          call femap_create(llapl,nd,ng,integ%femap(i))
       end do
    else
       i = 1
       call set_face_type    ( gfinf(1)%p%ftype, gfinf(1)%p%ftype,nd,utype)
       gnode(i) = gfinf(i)%p%nnode
       unode(i) = ufinf(i)%p%nnode
       nface(i) = gfinf(i)%p%nobje_dim(nd+1) -  gfinf(i)%p%nobje_dim(nd)
       call set_integ (utype,nd-1,gfinf(i)%p%order,ufinf(i)%p%order,ng,nlocs,lrule,llapl,           &
            &          gfnod(i),ufnod(i))
       ! Store identifiers of element
       integ%ltype(i) = gfinf(i)%p%ftype
       integ%ltype(2) = NULL_type_id

       ! Maximum amount of nodes in faces and elements
       mnode = max(gfnod(1),ufnod(1))
    
       call quadrature_create(lrule,nd-1,ng,integ%quad) 
       ! Create bf map
       call bomap_create(nd,ng,integ%bomap)
       ! Allocate interpolation
       call interpolation_create(1,1,    1,nd-1,gfnod(1),ng,integ%gint_ref)
       ! Define interpolations in the reference domain
       call interpolation_local(integ%quad%pos,integ%gint_ref)

       ! Quadrature for the faces of element i
       call face_quadrature_create(nd,nface(i),integ%quad,gfinf(i)%p,integ%fquad(i))
       
       ! Allocate interpolation for element i
       call face_interpolation_create (unode(i),ufnod(i),nd,ng,nlocs,nface(i),integ%ufint_ref(i))
       call face_interpolation_create (unode(i),ufnod(i),nd,ng,nlocs,       1,integ%ufint_phy(i))
       call face_interpolation_create (gnode(i),gfnod(i),nd,ng,nlocs,nface(i),integ%gfint_ref(i))
       call face_interpolation_create (gnode(i),gfnod(i),nd,ng,nlocs,       1,integ%gfint_phy(i))
       ! Create reference interpolations for element i
       if     (ufinf(i)%p%ftype == Q_type_id) then
          do j=1,nlocs
             auxpo(j) = integ%quad%pos(1,j)
          end do
          call face_interpolation_local_Q (auxpo,integ%ufint_ref(i),nd,ufinf(i)%p%order,nlocs)
       elseif (ufinf(i)%p%ftype == P_type_id) then
          call face_interpolation_local  (integ%fquad(i)%pos,integ%ufint_ref(i))
       else
          write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
       end if
       if     (gfinf(i)%p%ftype == Q_type_id) then
          do j=1,nlocs
             auxpo(j) = integ%quad%pos(1,j)
          end do
          call face_interpolation_local_Q (auxpo,integ%gfint_ref(i),nd,gfinf(i)%p%order,nlocs)
       elseif (gfinf(i)%p%ftype == P_type_id) then
          call face_interpolation_local  (integ%fquad(i)%pos,integ%gfint_ref(i))
       else
          write(*,*) __FILE__,__LINE__, 'ERROR! Unknown elemetn type.'
       end if
       !call face_interpolation_local  (integ%fquad(i)%pos,integ%ufint_ref(i))
       !call face_interpolation_local  (integ%fquad(i)%pos,integ%gfint_ref(i))
       ! Create FE map
       call femap_create(llapl,nd,ng,integ%femap(i))
    end if

  end subroutine face_integ_create

  ! =================================================================================================
  subroutine integ_element(ielem,gmesh,integ)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)    :: ielem      ! Element ID
    type(fem_mesh)        , intent(in)    :: gmesh      ! Geometry interpolation space
    type(vol_integ)       , intent(inout) :: integ      ! Volume integrator of ielem

    integer(ip)              :: nnode
    real(rp), allocatable    :: elcod(:,:)

    nnode = gmesh%pnods(ielem+1) - gmesh%pnods(ielem)
    call memalloc(gmesh%ndime,nnode,elcod, __FILE__,__LINE__ )
    
    call gather (gmesh%ndime,nnode,gmesh%lnods(gmesh%pnods(ielem):gmesh%pnods(ielem+1)-1), &
         & gmesh%coord,elcod)

    ! Define fe map  by interpolation
    call femap_from_interp(integ%gint_ref,elcod,integ%femap)
       
    ! Obtain physical interpolation
    call femap_apply_to_interp(integ%femap,integ%uint_ref,integ%uint_phy)

    call memfree(elcod,__FILE__,__LINE__)

  end subroutine integ_element

  !==================================================================================================
  subroutine integ_faces(elem,face,gmesh,integ,ntxob,nobje)
    implicit none
    ! Parameters
    integer(ip)              , intent(in)    :: elem(2),face(2)
    type(fem_mesh)           , intent(in)    :: gmesh      ! Geometry interpolation space
    type(face_integ)         , intent(inout) :: integ
    type(list)               , intent(in)    :: ntxob
    integer(ip)              , intent(in)    :: nobje

    integer(ip)                      :: i,gfnod,poins(integ%gint_ref%nnode),nelem
    real(rp)                         :: TOL=1e-8
    real(rp)                         :: facod(gmesh%ndime,integ%gint_ref%nnode)
    real(rp)   , allocatable         :: elcod(:,:)
    integer(ip), allocatable         :: elpos(:)
    integer(ip)                      :: nnode,inode,iobje,jobje

    ! Define face map (from reference to physical space) by interpolation
    ! TODO: Assert that the first element is not the one which the subfaces are in
    nnode = 0
    do inode = ntxob%p(nobje+face(1)),ntxob%p(nobje+face(1)+1)-1
       nnode = nnode +1
       poins(nnode) = gmesh%lnods(gmesh%pnods(elem(1))+ntxob%l(inode)-1)
    end do

    call gather (gmesh%ndime,integ%gint_ref%nnode,poins,gmesh%coord,facod)
    call bomap_from_interp(integ%gint_ref,facod,integ%bomap)

    if (elem(2)>0) then
       nelem = 2
    else
       nelem = 1
    end if

    do i = 1, nelem
       ! Define elements map  by interpolation
       call memalloc(gmesh%ndime,integ%gfint_ref(i)%nnode,elcod,__FILE__,__LINE__)
       call memalloc(integ%gfint_ref(i)%nnode,elpos,__FILE__,__LINE__)
       elpos = gmesh%lnods(gmesh%pnods(elem(i)):gmesh%pnods(elem(i)+1)-1)
       call gather (gmesh%ndime,integ%gfint_ref(i)%nnode,elpos,gmesh%coord,elcod)
       call femap_from_face_interp(integ%gfint_ref(i),elcod,face(i),integ%femap(i))
       call femap_face_to_interp(integ%femap(i),integ%ufint_ref(i),face(i),integ%ufint_phy(i))
       call femap_face_to_interp(integ%femap(i),integ%gfint_ref(i),face(i),integ%gfint_phy(i))
       call memfree(elcod,__FILE__,__LINE__)
       call memfree(elpos,__FILE__,__LINE__)
    end do
    
  end subroutine integ_faces

  ! =================================================================================================
  subroutine vol_integ_free(integ)
    implicit none
    ! Parameters
    type(vol_integ), intent(inout) :: integ    
    
    ! Destruct quadratures
    call quadrature_free(integ%quad)

    ! Destruct interpolations
    call interpolation_free(integ%gint_ref)
    call interpolation_free(integ%gint_phy)
    call interpolation_free(integ%uint_ref)
    call interpolation_free(integ%uint_phy)

    ! Destruct fe map
    call femap_free(integ%femap)

    integ%ltype = 0
  end subroutine vol_integ_free

  ! ================================================================================================
  subroutine face_integ_free(integ)
    implicit none
    ! Parameters
    type(face_integ), intent(inout) :: integ

    ! Destruct quadratures
    call quadrature_free     (integ%quad)
    call face_quadrature_free(integ%fquad(1))

    ! Destruct fe map
    call bomap_free(integ%bomap)
    call femap_free(integ%femap(1))
    call interpolation_free(integ%gint_ref)

    ! Deallocate Face interpolations
    call face_interpolation_free(integ%ufint_ref(1))
    call face_interpolation_free(integ%gfint_ref(1))
    call face_interpolation_free(integ%ufint_phy(1))
    call face_interpolation_free(integ%gfint_phy(1))

    if (integ%ltype(2) /= NULL_type_id) then
       call face_quadrature_free(integ%fquad(2))
       call femap_free(integ%femap(2))
       call face_interpolation_free(integ%ufint_ref(2))
       call face_interpolation_free(integ%gfint_ref(2))
       call face_interpolation_free(integ%ufint_phy(2))
       call face_interpolation_free(integ%gfint_phy(2))
    end if
  end subroutine face_integ_free

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
end module integration_names
