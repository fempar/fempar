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
module integrator_names
  use types_names
  use memor_names
#ifdef memcheck
  use iso_c_binding
#endif
  use reference_fe_names
  use SB_quadrature_names
  use SB_interpolation_names
  implicit none
  private

  type fe_map_t

     real(rp), allocatable  :: jacobian(:,:,:)     ! Map Jacobian         (ndime,ndime,nlocs)
     real(rp), allocatable  :: inv_jacobian(:,:,:) ! Map Jacobian inverse (ndime,ndime,nlocs)
     real(rp), allocatable  :: det_jacobian(:)     ! Map Jacobian det     (nlocs)
     real(rp), allocatable  :: d2sdx(:,:,:,:) ! 2nd derivatives (ndime,ndime,ndime,nlocs)
     real(rp), allocatable  :: coordinates_points(:,:)     ! Coordinates of evaluation points (ndime,nlocs)

     contains

       procedure :: get_det_jacobian

  end type fe_map_t

  type SB_volume_integrator_t 

     type(SB_quadrature_t) :: quadrature ! Quadrature rules for elements
     class(reference_fe_t), pointer :: reference_fe
     type(SB_interpolation_t) :: interpolation ! Unknown interpolation_t in the reference element domain
     class(reference_fe_t), pointer :: reference_fe_geometry
     type(SB_interpolation_t) :: interpolation_geometry ! Geometry interpolation_t in the reference element domain

     ! Working arrays
     type(SB_interpolation_t) :: interpolation_o_map ! Unknown interpolation_t in the physical element domain

     ! FE map
     type(fe_map_t) :: fe_map

   contains

     procedure :: create
     !procedure :: free
     procedure :: print
     procedure :: update
     procedure :: set_integration

     procedure :: get_reference_fe
     procedure :: get_quadrature
     procedure :: get_interpolation
     procedure :: get_fe_map

  end type SB_volume_integrator_t

  type SB_p_volume_integrator_t
     type(SB_volume_integrator_t)          , pointer :: p => NULL() 
  end type SB_p_volume_integrator_t

  public :: SB_volume_integrator_t, SB_p_volume_integrator_t, fe_map_t

# define var_type type(SB_p_volume_integrator_t)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

  !public :: volume_integrator_create, volume_integrator_free, volume_integrator_update, set_integ

contains

  ! =================================================================================================
# include "mem_body.i90"

  !==================================================================================================

  subroutine create( this, reference_fe, reference_fe_geometry, max_order )
    implicit none
    ! Parameters
    class(SB_volume_integrator_t), intent(out) :: this
    class(reference_fe_t), target, intent(in) :: reference_fe
    class(reference_fe_t), target, intent(in) :: reference_fe_geometry    
    integer(ip), optional, intent(in)  :: max_order

    integer(ip) :: ndime,ngaus
    ! Create quadrature
    this%reference_fe => reference_fe
    this%reference_fe_geometry => reference_fe_geometry

  end subroutine create

  subroutine set_integration( this, max_order )
    implicit none
    ! Parameters
    class(SB_volume_integrator_t), intent(inout) :: this   
    integer(ip), optional, intent(in)  :: max_order

    integer(ip) :: ndime,ngaus
    ! Create quadrature

    call this%reference_fe%create_quadrature( this%quadrature, max_order )

    call this%reference_fe%create_interpolation( this%quadrature, this%interpolation )

    call this%reference_fe%create_interpolation( this%quadrature, this%interpolation_geometry )

    this%interpolation_o_map = this%interpolation

    ndime = this%reference_fe%get_number_dimensions()
    ngaus = this%interpolation%get_number_evaluation_points()

    call memalloc(ndime,ndime,ngaus,this%fe_map%jacobian,__FILE__,__LINE__)
    call memalloc(ndime,ndime,ngaus,this%fe_map%inv_jacobian,__FILE__,__LINE__)
    call memalloc(ngaus,this%fe_map%det_jacobian,__FILE__,__LINE__)   
    call memalloc(ndime,ngaus,this%fe_map%coordinates_points,__FILE__,__LINE__)


    if( associated( this%interpolation%get_pointer_hessian() ) ) then
       call memalloc(ndime,ndime,ndime,this%interpolation%get_number_evaluation_points(), &
            & this%fe_map%d2sdx,  __FILE__,__LINE__)
    end if

  end subroutine set_integration

  ! =================================================================================================
  subroutine volume_integrator_free(integ)
    implicit none
    type(SB_volume_integrator_t), intent(inout) :: integ    

    !   ! Destruct quadratures
    !   call quadrature_free(integ%quad)

    !   ! Destruct interpolations
    !   call interpolation_free(integ%gint_ref)
    !   call interpolation_free(integ%gint_phy)
    !   call interpolation_free(integ%uint_ref)
    !   call interpolation_free(integ%uint_phy)

    !   ! Destruct fe map
    !   call memfree(map%jacobian,__FILE__,__LINE__)
    !   call memfree(map%inv_jacobian,__FILE__,__LINE__)
    !   call memfree(map%det_jacobian,__FILE__,__LINE__)

  end subroutine volume_integrator_free

  !==================================================================================================
  subroutine update(this, coordinates)
    implicit none
    ! Parameters
    class(SB_volume_integrator_t), intent(inout) :: this
    real(rp), intent(in) :: coordinates(:,:)
    ! Define fe map  by interpolation
    call create_femap( this%interpolation_geometry, coordinates, this%fe_map )
    ! Obtain physical interpolation
    call femap_apply_to_interp( this%fe_map, this%interpolation, this%interpolation_o_map )
  end subroutine update
  !==================================================================================================

  subroutine print( this )
    implicit none
    class(SB_volume_integrator_t), intent(inout) :: this

    write(*,*) 'PRINT VOLUME INTEGRATOR: '
    write(*,*) '%%%%%%%%%%%%%%%% quadrature: %%%%%%%%%%%%%%%%'
    call this%quadrature%print()
    write(*,*) '%%%%%%%%%%%%%%%% reference_fe: %%%%%%%%%%%%%%%%'
    call this%reference_fe%print()
    write(*,*) '%%%%%%%%%%%%%%%% interpolation: %%%%%%%%%%%%%%%%'
    call this%interpolation%print()
    write(*,*) '%%%%%%%%%%%%%%%% reference_fe_geometry: %%%%%%%%%%%%%%%%'
    call this%reference_fe_geometry%print()
    write(*,*) '%%%%%%%%%%%%%%%% interpolation_geometry: %%%%%%%%%%%%%%%%'
    call this%interpolation_geometry%print()
    write(*,*) '%%%%%%%%%%%%%%%% interpolation_o_map: %%%%%%%%%%%%%%%%'
    call this%interpolation_o_map%print()

    write(*,*) '%%%%%%%%%%%%%%%% jacobian: %%%%%%%%%%%%%%%%', this%fe_map%jacobian
    write(*,*) '%%%%%%%%%%%%%%%% inverse jacobian: %%%%%%%%%%%%%%%%',this%fe_map%inv_jacobian
    write(*,*) '%%%%%%%%%%%%%%%% determinant jacobian: %%%%%%%%%%%%%%%%',this%fe_map%det_jacobian
    if ( associated( this%interpolation%get_pointer_hessian() ) ) then 
       write(*,*) '%%%%%%%%%%%%%%%% second derivatives: %%%%%%%%%%%%%%%%',this%fe_map%d2sdx
    end if
    write(*,*) '%%%%%%%%%%%%%%%% coordinates points: %%%%%%%%%%%%%%%%',this%fe_map%coordinates_points

  end subroutine print

  function get_reference_fe ( this )
    implicit none
    class(SB_volume_integrator_t), target, intent(in) :: this
    class(reference_fe_t), pointer :: get_reference_fe
    get_reference_fe => this%reference_fe
  end function get_reference_fe

  function get_quadrature ( this )
    implicit none
    class(SB_volume_integrator_t), target, intent(in) :: this
    type(SB_quadrature_t), pointer :: get_quadrature
    get_quadrature => this%quadrature
  end function get_quadrature

  function get_interpolation ( this )
    implicit none
    class(SB_volume_integrator_t), target, intent(in) :: this
    type(SB_interpolation_t), pointer :: get_interpolation
    get_interpolation => this%interpolation
  end function get_interpolation

  function get_fe_map ( this )
    implicit none
    class(SB_volume_integrator_t), target, intent(in) :: this
    type(fe_map_t), pointer :: get_fe_map
    get_fe_map => this%fe_map
  end function get_fe_map

  function get_det_jacobian ( this, i )
    implicit none
    class(fe_map_t), target, intent(in) :: this
    integer(ip) :: i
    real(rp) :: get_det_jacobian
    get_det_jacobian = this%det_jacobian(i)
  end function get_det_jacobian

  !==============================================================================
  subroutine create_femap(int,elcod,map )
    !-----------------------------------------------------------------------
    ! A map obtained from the (usually isoparametric) interpolation of the geometry
    !-----------------------------------------------------------------------
    implicit none
    type(SB_interpolation_t), intent(in)    :: int
    !    real(rp)           , intent(inout) :: hnatu
    real(rp)           , intent(in)    :: elcod(:,:)
    type(fe_map_t)        , intent(inout) :: map
    ! Locals
    real(rp), allocatable :: wmat1(:,:,:)
    real(rp), allocatable :: wmat2(:,:,:), wvec1(:)
    real(rp)    :: hnatu
    real(rp)    :: enor0,h_tem
    integer(ip) :: ndime,nnode,nlocs,ntens
    integer(ip) :: ilocs,idime,jdime,kdime,ldime,inode,itens
    logical :: khes

    khes = .false.
    if ( associated( int%get_pointer_hessian() ) ) then 
       khes = .true.
    end if

    ! ! Check and get data from int
    ! assert(int%kder==1)
     ndime = int%get_number_dimensions()
     nnode = int%get_number_shape_functions()
     nlocs = int%get_number_evaluation_points()

    ! ! Check elcod
    ! assert(ndime==size(elcod,dim=1))
    ! assert(nnode==size(elcod,dim=2))

    ! ! Check map allocation (jainv, detjm and hleng
    ! ! are assumed to be allocated if jacob is)
    ! assert(ndime==size(map%jacob,dim=1))
    ! assert(nlocs==size(map%jacob,dim=3))

    ! Jacobian and its inverse
    ! write(*,*) 'elcod: ',elcod

    do ilocs=1,nlocs
       ! Matmul is not thread safe
       !map%jacob(:,:,ilocs)=matmul(elcod,transpose(int%deriv(:,:,ilocs)))
       map%jacobian(:,:,ilocs)=0.0_rp
       do inode=1,nnode
          do jdime=1,ndime
             do idime=1,ndime
                map%jacobian(idime,jdime,ilocs) = map%jacobian(idime,jdime,ilocs) &
                     + elcod(idime,inode)*int%get_shape_derivative(jdime,inode,ilocs)
             end do
          end do
       end do
       ! J^(-t)
       call invmtx(map%jacobian(:,:,ilocs),map%inv_jacobian(:,:,ilocs),map%det_jacobian(ilocs),ndime)
    end do

    ! Evaluation (Gauss) point coordinates
    do ilocs=1,nlocs
       map%coordinates_points(:,ilocs)=0.0_rp
       do inode=1,nnode
          do idime=1,ndime
             map%coordinates_points(idime,ilocs) = map%coordinates_points(idime,ilocs) &
                  + elcod(idime,inode)*int%get_shape_function(inode,ilocs)
          end do
       end do
    end do

    ! ! Second derivatives of the map
    if( khes ) then

       ntens=int%get_number_entries_symmetric_tensor()
       ! Check that second derivativesof the map have been allocated.
       assert(ndime==size(map%d2sdx,dim=1))
       assert(nlocs==size(map%d2sdx,dim=4))

       call memalloc(ndime,ndime,nnode,wmat1,__FILE__,__LINE__)
       call memalloc(ndime,ndime,nnode,wmat2,__FILE__,__LINE__)
       call memalloc(ntens,wvec1,__FILE__,__LINE__)


       do ilocs=1,nlocs

          ! Transforms the array HESSI to a symmetric matrix WMAT1
          do inode=1,nnode
             do itens = 1, ntens
                wvec1(itens) = int%get_hessian(itens,inode,ilocs)
             end do
             call vetoma(wvec1,wmat1(1,1,inode),ndime,ntens)
          end do

          ! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for
          ! each node
          do inode=1,nnode
             call btdbma(wmat2(1,1,inode),wmat1(1,1,inode), &
                  &        map%inv_jacobian(:,:,ilocs),ndime,ndime)
          end do

          ! Obtains (d^2 s_k / d x_i d x_j) as the solution of the system
          ! (d x_l / d s_k) (d^2 s_k / d x_i d x_j) 
          !     = - (d^2 x_l / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j), 
          ! for l,i,j = 1,...,NDIME
          do kdime=1,ndime
             do idime=1,ndime
                do jdime=1,ndime
                   map%d2sdx(kdime,idime,jdime,ilocs)=0.0_rp
                   do ldime=1,ndime
                      do inode=1,nnode
                         map%d2sdx(kdime,idime,jdime,ilocs) =    &
                              & map%d2sdx(kdime,idime,jdime,ilocs) &
                              & - map%inv_jacobian(kdime,ldime,ilocs)     &
                              &   * wmat2(idime,jdime,inode) * elcod(ldime,inode)
                      end do
                   end do
                end do
             end do
          end do

       end do

       call memfree(wmat1,__FILE__,__LINE__)
       call memfree(wmat2,__FILE__,__LINE__)

    end if


  end subroutine create_femap


  subroutine invmtx(a,b,deter,nsize)
    !-----------------------------------------------------------------------
    !
    ! This routine inverts a square matrix A -> Mat(nsize,nsize). The
    ! inverse is stored in B. Its determinant is DETER
    !    
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: nsize
    real(rp),    intent(in)  :: a(nsize,nsize)
    real(rp),    intent(out) :: b(nsize,nsize),deter
    integer(ip)              :: isize,jsize
    real(rp)                 :: denom,t1,t2,t3,t4

    select case (nsize)

    case(1)
       deter=a(1,1)
       if(deter==0.0_rp) return
       b(1,1) = 1.0_rp/a(1,1)

    case(2)
       deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
       if(deter/=0.0_rp) then
          denom=1.0_rp/deter
          b(1,1) = a(2,2)*denom
          b(2,2) = a(1,1)*denom
          b(2,1) =-a(2,1)*denom
          b(1,2) =-a(1,2)*denom 
       end if

    case(3)
       t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
       t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
       t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
       deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
       if(deter==0.0_rp) return
       denom = 1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
       b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
       b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
       b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
       b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
       b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

    case(4)
       t1= a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
            + a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
            - a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
       t2=-a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
            - a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
            + a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
       t3=+a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
            + a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
            - a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
       t4=-a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
            - a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
            + a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
       deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
       if(deter==0.0_rp) return
       denom=1.0_rp/deter
       b(1,1) = t1*denom
       b(2,1) = t2*denom
       b(3,1) = t3*denom
       b(4,1) = t4*denom
       b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
            - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
            + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
       b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
            + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
            - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
       b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
            - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
            + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
       b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
            + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
            - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
       b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
            + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
            - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
       b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
            - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
            + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
       b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
            + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
            - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
       b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
            - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
            + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
       b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
            - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
            + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
       b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
            + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
            - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
       b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
            - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
            + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
       b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
            + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
            - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom

    case default
       b=a
       call invert(b,nsize,nsize)

    end select

  end subroutine invmtx

  !-----------------------------------------------------------------------
  subroutine invert(a,nmax,ndm)
    !-----------------------------------------------------------------------
    !
    ! This routine performs the inversion of a ndm*ndm square matrix 
    ! or just part of it (nmax*nmax)
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)    :: ndm,nmax
    real(rp),    intent(inout) :: a(ndm,ndm)
    real(rp)                   :: d
    integer(ip)                :: n,j,i

    do n = 1,nmax
       d = a(n,n)
       do j = 1,nmax
          a(n,j) = -a(n,j)/d
       end do
       do i = 1,nmax
          if(n/=i) then
             do j = 1,nmax
                if(n/=j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
             end do
          end if
          a(i,n) = a(i,n)/d
       end do
       a(n,n) = 1.0_rp/d
    end do

  end subroutine invert


  !-----------------------------------------------------------------------
  subroutine vetoma(vecto,xmatr,ndime,ntens)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine stores a vector VECTO as a symmetric matrix XMATR
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ntens
    real(rp)   , intent(in)  :: vecto(ntens)
    real(rp)   , intent(out) :: xmatr(ndime,ndime)

    if(ndime.eq.2) then
       xmatr(1,1)=vecto(1)
       xmatr(1,2)=vecto(3)
       xmatr(2,1)=vecto(3)
       xmatr(2,2)=vecto(2)
    else
       xmatr(1,1)=vecto(1)
       xmatr(1,2)=vecto(4)
       xmatr(1,3)=vecto(5)
       xmatr(2,1)=vecto(4)
       xmatr(2,2)=vecto(2)
       xmatr(2,3)=vecto(6)
       xmatr(3,1)=vecto(5)
       xmatr(3,2)=vecto(6)
       xmatr(3,3)=vecto(3)
    end if

  end subroutine vetoma

  !-----------------------------------------------------------------------
  subroutine btdbma(aglob,aloca,bmatr,n1,n2)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine computes Ag = Bt Al B  when Ag and Al are stored as full
    ! matrices (Ag := aglob, Al := aloca, B := bmatr). The dimensions are
    ! Al -> Mat(n1,n1), Ag -> Mat(n2,n2), B -> Mat(n2,n1) 
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: n1,n2
    real(rp)   , intent(in)  :: aloca(n1,n1), bmatr(n1,n2)
    real(rp)   , intent(out) :: aglob(n2,n2)
    integer(ip)              :: i,j,k,l

    do i=1,n2
       do j=1,n2
          aglob(i,j)=0.0
          do k=1,n1
             do l=1,n1
                aglob(i,j)=aglob(i,j)+bmatr(k,i)*aloca(k,l)*bmatr(l,j)
             end do
          end do
       end do
    end do

  end subroutine btdbma

  !==============================================================================
  !
  !==============================================================================
  subroutine femap_apply_to_interp(map,ref,phy)
    implicit none
    type(SB_interpolation_t), intent(in)    :: ref
    type(fe_map_t)        , intent(in)    :: map
    type(SB_interpolation_t), intent(inout) :: phy
    real(rp), allocatable :: wmat1(:,:,:)
    real(rp), allocatable :: wmat2(:,:,:), wvec1(:)
    integer(ip) :: ndime,nnode,nlocs,ntens
    integer(ip) :: ilocs,idime,jdime,kdime,ldime,inode,itens
    logical :: khes
    real(rp), pointer :: ref_shape_derivatives(:,:,:), ref_hessian(:,:,:)
    real(rp), pointer :: phy_shape_derivatives(:,:,:), phy_hessian(:,:,:)

    ref_shape_derivatives => ref%get_pointer_shape_derivatives()
    ref_hessian => ref%get_pointer_hessian()

    phy_shape_derivatives => phy%get_pointer_shape_derivatives()
    phy_hessian => phy%get_pointer_hessian()

    khes = .false.
    if ( associated( ref_hessian ) ) then 
       khes = .true.
    end if

    ndime = ref%get_number_dimensions()
    ntens = ref%get_number_entries_symmetric_tensor()
    nlocs = ref%get_number_evaluation_points()
    nnode = ref%get_number_shape_functions()

    ! Shape functions do not change
    ! phy%shape = ref%shape ! SB. not needed

    ! First derivatives do
    !if(phy%kder==1) then
    phy_shape_derivatives=0.0_rp
    do ilocs=1,phy%get_number_evaluation_points()
       do inode=1,phy%get_number_shape_functions()
          do idime=1,ndime
             do jdime=1,ndime
                phy_shape_derivatives(idime,inode,ilocs) = phy_shape_derivatives(idime,inode,ilocs) &
                     + map%inv_jacobian(jdime,idime,ilocs)*ref_shape_derivatives(jdime,inode,ilocs)
             end do
          end do
       end do
    end do
    !end if

    ! Second derivatives are
    !
    !    d^2 N / d x_i d x_j
    !       = (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j)
    !       + (d N / d s_k) (d^2 s_k / d x_i d x_j) 
    !
    if( khes ) then

       call memalloc(ndime,ndime,nnode,wmat1,__FILE__,__LINE__)
       call memalloc(ndime,ndime,nnode,wmat2,__FILE__,__LINE__)

       do ilocs=1,nlocs

          if( khes ) then
             ! Transforms the array HESSI to a symmetric matrix WMAT1
             do inode=1,nnode
                do itens = 1, ntens
                   wvec1(itens) = ref%get_hessian(itens,inode,ilocs)
                end do
                call vetoma(wvec1,wmat1(1,1,inode),ndime,ntens)
             end do
             ! Computes (d^2 N / d s_k d s_l) (d s_k/ d x_i)(d s_l/ d x_j) for each node
             do inode=1,nnode
                call btdbma(wmat2(1,1,inode),wmat1(1,1,inode),map%inv_jacobian(:,:,ilocs), &
                     &        ndime,ndime)
             end do
          end if

          if( khes ) then
             ! Add second cartesian derivatives of the map times 
             ! first derivatives of shape functions
             do inode=1,nnode
                do idime=1,ndime
                   do jdime=1,ndime
                      do kdime=1,ndime
                         wmat2(idime,jdime,inode)=wmat2(idime,jdime,inode) &
                              & + ref_shape_derivatives(kdime,inode,ilocs) &
                              &   * map%d2sdx(kdime,idime,jdime,ilocs)
                      end do
                   end do
                end do
             end do
          end if

          ! Writes the Hessian matrix as an array
          do inode=1,nnode
             do itens = 1, ntens
                wvec1(itens) = phy%get_hessian(itens,inode,ilocs)
             end do
             call matove(wmat2(1,1,inode),wvec1,ndime,ntens)
          end do

       end do

       call memfree(wmat1,__FILE__,__LINE__)
       call memfree(wmat2,__FILE__,__LINE__)

    end if

  end subroutine femap_apply_to_interp


  !-----------------------------------------------------------------------
  subroutine matove(xmatr,vecto,ndime,ntens)
    !-----------------------------------------------------------------------
    !                                      
    ! This routine stores a symmetric matrix XMATR into a vector VECTO
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip), intent(in)  :: ndime,ntens
    real(rp)   , intent(in)  :: xmatr(ndime,ndime)
    real(rp)   , intent(out) :: vecto(ntens)

    if(ndime.eq.2) then
       vecto(1)=xmatr(1,1)
       vecto(3)=xmatr(1,2)
       vecto(2)=xmatr(2,2)
    else
       vecto(1)=xmatr(1,1)
       vecto(4)=xmatr(1,2)
       vecto(2)=xmatr(2,2)
       vecto(5)=xmatr(1,3)
       vecto(6)=xmatr(2,3)
       vecto(3)=xmatr(3,3)
    end if

  end subroutine matove


end module integrator_names
