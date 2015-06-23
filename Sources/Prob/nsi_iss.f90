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
module nsi_cg_iss_names
 use types
 use memor
 use array_names
 use problem_names
 use nsi_names
 use fem_element_names
 use eltrm_gen_names
 use element_fields_names
 use element_tools_names
 use analytical_names
 implicit none
# include "debug.i90"

 private
 
 ! INF-SUP STABLE (iss) NAVIER-STOKES types 
 ! Problem data
 type, extends(discrete_problem) :: nsi_cg_iss_discrete
    integer(ip) ::   & 
         kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
         kfl_lump,   & ! Flag for lumped mass submatrix
         tdimv,      & ! Number of temporal steps stored for velocity
         tdimp         ! Number of temporal steps stored for pressure   
    real(rp) ::         &
         dtinv,         & ! Inverse of time step
         ctime,         & ! Current time
         ktauc,         & ! Constant multiplying stabilization parameter tau_c 
         k1tau,         & ! C1 constant on stabilization parameter tau_m
         k2tau            ! C2 constant on stabilization parameter tau_m
  contains
    procedure :: create  => nsi_create_discrete
 end type nsi_cg_iss_discrete
    
!!$ ! Aproximation
!!$ type, extends(discrete_problem) :: nsi_cg_iss_approximation
!!$    type(nsi_cg_iss_data), pointer :: data
!!$     contains
!!$      procedure :: create  => nsi_create
!!$      procedure :: compute => nsi_matvec ! Default element computation points to matvec
!!$ end type nsi_cg_iss_approximation

 ! Matvec
 type, extends(discrete_integration) :: nsi_cg_iss_matvec
      type(nsi_cg_iss_discrete), pointer :: discret
      type(nsi_problem)        , pointer :: physics
     contains
      procedure :: create  => nsi_matvec_create
      procedure :: compute => nsi_matvec 
      procedure :: free    => nsi_matvec_free
 end type nsi_cg_iss_matvec

 ! Unkno components parameter definition
 integer(ip), parameter :: current   = 1
 integer(ip), parameter :: prev_iter = 2
 integer(ip), parameter :: prev_step = 3

 ! Types
 public :: nsi_cg_iss_matvec, nsi_cg_iss_discrete

 ! Functions
 ! public :: nsi_iss_create, nsi_iss_free, nsi_iss_element_matvec, nsi_iss_element_real, &
 !      &    nsi_iss_result_write, nsi_iss_element_mat, nsi_iss_element_1vec,            &
 !      &    nsi_iss_create_vor, nsi_iss_free_vor, nsi_iss_update_vorti,                 &
 !      &    nsi_iss_create_tder, nsi_iss_free_tder, nsi_iss_create_oss, nsi_iss_face_mat
 
contains

  !=================================================================================================
  subroutine nsi_create_discrete(discret,physics,l2g)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_discrete), intent(out) :: discret
    class(physical_problem)   , intent(in)  :: physics
    integer(ip), optional     , intent(in)  :: l2g(:)
    ! Locals
    integer(ip) :: i

    ! Flags
    discret%kfl_lump = 0 ! Flag for lumped mass submatrix (Off=0, On=1)
    discret%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=1)

    ! Problem variables
    discret%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter tau_m
    discret%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    discret%ktauc  = 0.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Time integration variables
    discret%dtinv  = 1.0_rp ! Inverse of time step
    discret%ctime  = 0.0_rp ! Current time
    discret%tdimv  =  2     ! Number of temporal steps stored for velocity
    discret%tdimp  =  2     ! Number of temporal steps stored for pressure

    discret%nvars = physics%ndime+1
    call memalloc(discret%nvars,discret%l2g_var,__FILE__,__LINE__)
    if ( present(l2g) ) then
       assert ( size(l2g) == discret%nvars )
       discret%l2g_var = l2g
    else
       do i = 1,discret%nvars
          discret%l2g_var(i) = i
       end do
    end if
    
  end subroutine nsi_create_discrete

  !=================================================================================================
  subroutine nsi_matvec_create( approx, physics, discret )
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_matvec)        , intent(inout) :: approx
    ! type(nsi_problem)         , target, intent(in)  :: physics
    ! type(nsi_cg_iss_discrete) , target, intent(in)  :: discret

    ! approx%physics => physics
    ! approx%discret => discret

    class(physical_problem), target , intent(in)   :: physics
    class(discrete_problem), target , intent(in)   :: discret

    select type (physics)
    type is(nsi_problem)
       approx%physics => physics
       class default
       check(.false.)
    end select
    select type (discret)
    type is(nsi_cg_iss_discrete)
       approx%discret => discret
       class default
       check(.false.)
    end select

  end subroutine nsi_matvec_create

  !=================================================================================================
  subroutine nsi_matvec_free(approx)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_matvec), intent(inout) :: approx

    approx%physics => null()
    approx%discret => null()

  end subroutine nsi_matvec_free

  !=================================================================================================
  subroutine nsi_matvec(approx,start,elem)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental matrix-vector integration selection.                !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_matvec), intent(inout) :: approx
    integer(ip)             , intent(in)    :: start(:)
    type(fem_element)       , intent(inout) :: elem
    ! Locals
    real(rp), allocatable :: elmat_vu(:,:,:,:)
    real(rp), allocatable :: elmat_vu_diag(:,:)
    real(rp), allocatable :: elmat_vp(:,:,:,:)
    real(rp), allocatable :: elmat_qu(:,:,:,:)
    real(rp), allocatable :: elvec_u(:,:) 
    integer(ip)           :: igaus,idime,inode,jdime,jnode
    integer(ip)           :: ngaus,ndime,nnodu,nnodp
    real(rp)              :: ctime,dtinv,dvolu,diffu,react
    real(rp)              :: work(4)
    real(rp)              :: agran(elem%integ(1)%p%uint_phy%nnode)
    type(vector)          :: gpvel, gpveln, force

    ! Checks
    !check(elem%f_inf(1)%p%order > elem%f_inf(approx%physics%ndime+1)%p%order)
    check(elem%integ(1)%p%quad%ngaus == elem%integ(approx%physics%ndime+1)%p%quad%ngaus)
    do idime=2,approx%physics%ndime
       check(elem%integ(1)%p%uint_phy%nnode == elem%integ(idime)%p%uint_phy%nnode)
    end do

    ! Initialize matrix and vector
    elem%p_mat%a = 0.0_rp
    elem%p_vec%a = 0.0_rp

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = elem%integ(1)%p%uint_phy%nnode
    nnodp = elem%integ(ndime+1)%p%uint_phy%nnode
    ngaus = elem%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    react = approx%physics%react
    dtinv = approx%discret%dtinv

    ! Allocate auxiliar matrices and vectors
    call memalloc(ndime,ndime,nnodu,nnodu,elmat_vu,__FILE__,__LINE__)
    call memalloc(nnodu,nnodu,elmat_vu_diag,__FILE__,__LINE__)
    call memalloc(ndime,1,nnodu,nnodp,elmat_vp,__FILE__,__LINE__)
    call memalloc(1,ndime,nnodp,nnodu,elmat_qu,__FILE__,__LINE__)
    call memalloc(ndime,nnodu,elvec_u,__FILE__,__LINE__)

    ! Initialize to zero
    elmat_vu      = 0.0_rp
    elmat_vu_diag = 0.0_rp
    elmat_vp      = 0.0_rp
    elmat_qu      = 0.0_rp
    elvec_u       = 0.0_rp

    ! Interpolation operations for velocity
    call create_vector (approx%physics, 1, elem%integ, gpvel)
    call create_vector (approx%physics, 1, elem%integ, gpveln)
    call interpolation (elem%unkno, 1, prev_iter, elem%integ, gpvel)
    gpvel%a=0.0_rp
    if(dtinv == 0.0_rp) then
       call interpolation (elem%unkno, 1, prev_step, elem%integ, gpveln)
    else
       gpveln%a = 0.0_rp
    end if

    ! Set real time
    if(approx%discret%kfl_thet==0) then        ! BE
       ctime = approx%discret%ctime
    elseif(approx%discret%kfl_thet==1) then    ! CN
       ctime = approx%discret%ctime - 1.0_rp/dtinv
    end if

    ! Set force term
    call create_vector(approx%physics,1,elem%integ,force)
    force%a=0.0_rp
    ! Impose analytical solution
    if(approx%physics%case_veloc>0.and.approx%physics%case_press>0) then 
       call nsi_analytical_force(approx,elem,ctime,gpvel,force)
    end if
    
    ! Initializations
    work     = 0.0_rp
    agran    = 0.0_rp 

    ! Loop on Gauss points
    do igaus = 1,elem%integ(1)%p%quad%ngaus
       dvolu = elem%integ(1)%p%quad%weight(igaus)*elem%integ(1)%p%femap%detjm(igaus)

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + &
                     &         gpvel%a(idime,igaus)*elem%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
          end do
       end if

       ! Add external force term
       do idime=1,approx%physics%ndime    
          force%a(idime,igaus) = force%a(idime,igaus) + approx%physics%gravi(idime)
       end do

       ! Computation of elemental terms
       ! ------------------------------
       ! Block U-V
       ! mu * ( grad u, grad v )
       call elmvis_gal(dvolu,diffu,elem%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elmat_vu_diag, &
            &          work)

       ! Add cross terms for symmetric grad
       if(approx%physics%kfl_symg==1) then
          call elmvis_gal_sym(dvolu,diffu,elem%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
               &              elmat_vu,work)
       end if
       if(approx%physics%kfl_skew==0) then
          ! (v, a·grad u) + s*(v,u) + (v, u/dt)
          call elmbuv_gal(dvolu,react,dtinv,elem%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &          elmat_vu_diag,work)
       elseif(approx%physics%kfl_skew==1) then
          ! 1/2(v, a·grad u) - 1/2(u,a·grad v) + s*(v,u) + (v, u/dt)
          call elmbuv_gal_skew1(dvolu,react,dtinv,elem%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &                elmat_vu_diag,work)
       end if

       ! Block P-V
       ! - ( div v, p )
       call elmbpv_gal_div_iss(dvolu,elem%integ(ndime+1)%p%uint_phy%shape(:,igaus),               &
            &                  elem%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_vp,work)

       ! Block U-Q
       ! ( div u, q )
       call elmbuq_gal_div_iss(dvolu,elem%integ(ndime+1)%p%uint_phy%shape(:,igaus),               &
            &                  elem%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_qu,work)

       ! RHS: Block U
       ! ( v, f ) + ( v, u_n/dt )
       call elmrhu_gal(dvolu,dtinv,elem%integ(1)%p%uint_phy%shape(:,igaus),gpveln%a(:,igaus), &
            &          force%a(:,igaus),nnodu,ndime,elvec_u,work)

    end do

    call memfree(gpvel%a,__FILE__,__LINE__)
    call memfree(gpveln%a,__FILE__,__LINE__)
    call memfree(force%a,__FILE__,__LINE__)


    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             do jdime=1,ndime
                ! Block V-U
                elem%p_mat%a(start(idime)+inode-1,start(jdime)+jnode-1) =           &
                     & elem%p_mat%a(start(idime)+inode-1,start(jdime)+jnode-1) +    &
                     & elmat_vu(idime,jdime,inode,jnode)
             end do    
             ! Block V-U (diag)
             elem%p_mat%a(start(idime)+inode-1,start(idime)+jnode-1) =              &
                  & elem%p_mat%a(start(idime)+inode-1,start(idime)+jnode-1) +       &
                  & elmat_vu_diag(inode,jnode)
          end do
          do jnode=1,nnodp
             ! Block V-P
             elem%p_mat%a(start(idime)+inode-1,start(ndime+1)+jnode-1) =            &
                  & elem%p_mat%a(start(idime)+inode-1,start(ndime+1)+jnode-1) +     &
                  & elmat_vp(idime,1,inode,jnode)
             ! Block Q-U
             elem%p_mat%a(start(ndime+1)+jnode-1,start(idime)+inode-1) =            &
                  & elem%p_mat%a(start(ndime+1)+jnode-1,start(idime)+inode-1) +     &
                  & elmat_qu(1,idime,jnode,inode)
          end do
          ! Block U
          elem%p_vec%a(start(idime)+inode-1) = elem%p_vec%a(start(idime)+inode-1) + &
               & elvec_u(idime,inode)
       end do
    end do

    ! Deallocate auxiliar matrices and vectors
    call memfree(elmat_vu,__FILE__,__LINE__)
    call memfree(elmat_vu_diag,__FILE__,__LINE__)
    call memfree(elmat_vp,__FILE__,__LINE__)
    call memfree(elmat_qu,__FILE__,__LINE__)
    call memfree(elvec_u,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(elem) 
    
  end subroutine nsi_matvec

  !==================================================================================================
  subroutine nsi_analytical_force(approx,elem,ctime,gpvel,force)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine computes the elemental force needed to impose an analytical solution        !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_matvec), intent(in)    :: approx
    type(fem_element)       , intent(in)    :: elem
    real(rp)                , intent(in)    :: ctime
    type(vector)            , intent(in)    :: gpvel
    type(vector)            , intent(inout) :: force
    ! Locals
    integer(ip) :: igaus,idime,conve
    real(rp)    :: gpvno
    real(rp)    :: parv(30),parp(10),part(3),part_p(3)

    ! Loop over gauss points
    do igaus=1,elem%integ(1)%p%quad%ngaus
    
       ! Compute velocity norm
       gpvno=0.0_rp
       do idime = 1,approx%physics%ndime
          gpvno = gpvno + gpvel%a(idime,igaus)*gpvel%a(idime,igaus)
       end do
       gpvno = sqrt(gpvno)

       ! Set convection flag
       if(gpvno.lt.1e-4) then
          conve = 0
       else
          conve = approx%physics%kfl_conv
       end if
       
       ! Evaluate unknowns and derivatives
       parv   = 0.0_rp
       parp   = 0.0_rp
       part   = 0.0_rp
       part_p = 0.0_rp
       call analytical_field(approx%physics%case_veloc,approx%physics%ndime, &
            &                elem%integ(1)%p%femap%clocs(:,igaus),ctime,parv)
       call analytical_field(approx%physics%case_press,approx%physics%ndime, &
            &                elem%integ(1)%p%femap%clocs(:,igaus),ctime,parp)
       call analytical_field(approx%physics%case_tempo,approx%physics%ndime, &
            &                elem%integ(1)%p%femap%clocs(:,igaus),ctime,part)
       call analytical_field(approx%physics%case_t_pre,approx%physics%ndime, &
            &                elem%integ(1)%p%femap%clocs(:,igaus),ctime,part_p)

       ! Evaluate force
       call nsi_force(approx%physics%ndime,approx%physics%diffu,approx%physics%react,conve, &
            &         approx%physics%kfl_symg,approx%physics%case_tempo,parv,parp,part,     &
            &         approx%physics%gravi,force%a(:,igaus),part_p)

    end do

  end subroutine nsi_analytical_force

  !=================================================================================================
  subroutine nsi_force(ndime,diffu,react,kfl_conv,kfl_symg,caset,parv,parp,part,gravi,force,part_p)
    !-----------------------------------------------------------------------------------------------!
    !   This subroutine evaluates the elemental force needed to impose an analytical solution       !
    !-----------------------------------------------------------------------------------------------!
    implicit none
    integer(ip)       , intent(in)    :: ndime,kfl_conv,kfl_symg,caset
    real(rp)          , intent(in)    :: diffu,react,parv(:),parp(:),part(:),gravi(3)
    real(rp)          , intent(inout) :: force(ndime)
    real(rp), optional, intent(in)    :: part_p(:)

    real(rp)   :: u,v,w,dpdx,dpdy,dpdz,d2udx,d2udy,d2udz
    real(rp)   :: d2vdx,d2vdy,d2vdz,d2wdx,d2wdy,d2wdz
    real(rp)   :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    real(rp)   :: d2udxy,d2udxz,d2vdxy,d2vdyz,d2wdxz,d2wdyz
    real(rp)   :: dudt,dvdt,dwdt

    !
    if(caset>0) then
       dudt = part(2)*parv(1); dvdt = part(2)*parv(2); dwdt = part(2)*parv(3)
    else
       dudt=0.0_rp; dvdt=0.0_rp; dwdt=0.0_rp
    end if
    !
    u = parv(1)*part(1); v = parv(2)*part(1); w = parv(3)*part(1)
    dudx = parv(4)*part(1); dudy = parv(5)*part(1); dudz = parv(6)*part(1)
    dvdx = parv(7)*part(1); dvdy = parv(8)*part(1); dvdz = parv(9)*part(1)
    dwdx = parv(10)*part(1); dwdy = parv(11)*part(1); dwdz = parv(12)*part(1)
    !
    d2udx = parv(13)*part(1); d2udy = parv(14)*part(1); d2udz = parv(15)*part(1)
    d2vdx = parv(16)*part(1); d2vdy = parv(17)*part(1); d2vdz = parv(18)*part(1)
    d2wdx = parv(19)*part(1); d2wdy = parv(20)*part(1); d2wdz = parv(21)*part(1)
    !
    d2udxy = parv(22)*part(1); d2udxz = parv(23)*part(1)
    d2vdxy = parv(25)*part(1); d2vdyz = parv(27)*part(1)
    d2wdxz = parv(29)*part(1); d2wdyz = parv(30)*part(1)
    !
    if(present(part_p)) then
       dpdx = parp(2)*part_p(1); dpdy = parp(3)*part_p(1); dpdz = parp(4)*part_p(1)
    else
       dpdx = parp(2); dpdy = parp(3); dpdz = parp(4)
    end if
    !

    force(1) = dudt - diffu*(d2udx+d2udy+d2udz) + dpdx + react*u + &
         REAL(kfl_conv)*(u*dudx + v*dudy + w*dudz) - gravi(1) -    &
         REAL(kfl_symg)*diffu*(d2udx+d2vdxy+d2wdxz)
    force(2) = dvdt - diffu*(d2vdx+d2vdy+d2vdz) + dpdy + react*v + &
         REAL(kfl_conv)*(u*dvdx + v*dvdy + w*dvdz) - gravi(2) -    &
         REAL(kfl_symg)*diffu*(d2udxy+d2vdy+d2wdyz)
    if(ndime==3) force(3) = dwdt - diffu*(d2wdx+d2wdy+d2wdz) + dpdz + react*w + &
         REAL(kfl_conv)*(u*dwdx + v*dwdy + w*dwdz) - gravi(3) -    &
         REAL(kfl_symg)*diffu*(d2udxz+d2vdyz+d2wdz)

  end subroutine nsi_force
  
  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat(prob,ielem,elem)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine performs the elemental matrix integration selection.                       !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   integer(ip),           intent(in)    :: ielem
  !   type(fem_element),     intent(inout) :: elem
  !   ! Locals
  !   type(block_elmat)        :: blk_elmat
  !   type(fem_blocks)         :: vars_s
  !   integer(ip)              :: ldofs(2),nblks,nnode(2),nd,nv,mn
  !   integer(ip), allocatable :: nnode_oss(:),ldofs_oss(:)

  !   ! Asserts
  !   assert(elem%nint==2)
  !   assert(elem%iv(prob%ndime+1)==2)
  !   assert(elem%integ(1)%p%quad%ngaus==elem%integ(2)%p%quad%ngaus)
    
  !   ! Allocate blk_elmat & blk_elvec locals
  !   nnode(1) = elem%f_inf(1)%p%nnode
  !   nnode(2) = elem%f_inf(2)%p%nnode
  !   if(prob%kfl_stab==2) then
  !      nblks=3
  !      call memalloc(nblks,nnode_oss,__FILE__,__LINE__)
  !      call memalloc(nblks,ldofs_oss,__FILE__,__LINE__)
  !      ldofs_oss(1) = prob%ndime; ldofs_oss(2)=1; ldofs_oss(3)=prob%ndime;
  !      nnode_oss(1) = elem%f_inf(1)%p%nnode
  !      nnode_oss(2) = elem%f_inf(2)%p%nnode
  !      nnode_oss(3) = elem%f_inf(1)%p%nnode
  !      call fem_blocks_alloc(scalar,nblks,ldofs_oss,vars_s)
  !      call block_elmat_alloc(vars_s,nnode_oss,blk_elmat)
  !   elseif(prob%kfl_stab==3) then
  !      nblks=4
  !      call memalloc(nblks,nnode_oss,__FILE__,__LINE__)
  !      call memalloc(nblks,ldofs_oss,__FILE__,__LINE__)
  !      ldofs_oss(1) = prob%ndime; ldofs_oss(2)=1; ldofs_oss(3)=prob%ndime; ldofs_oss(4)=prob%ndime
  !      nnode_oss(1) = elem%f_inf(1)%p%nnode
  !      nnode_oss(2) = elem%f_inf(2)%p%nnode
  !      nnode_oss(3) = elem%f_inf(1)%p%nnode
  !      nnode_oss(4) = elem%f_inf(1)%p%nnode
  !      call fem_blocks_alloc(scalar,nblks,ldofs_oss,vars_s)
  !      call block_elmat_alloc(vars_s,nnode_oss,blk_elmat)
  !   else
  !      nblks=2
  !      ldofs(1) = prob%ndime; ldofs(2)=1
  !      call fem_blocks_alloc(scalar,nblks,ldofs,vars_s)
  !      call block_elmat_alloc(vars_s,nnode,blk_elmat)
  !   end if

  !   ! Initialize to zero
  !   blk_elmat%blocks(1,1)%data = 0.0_rp
  !   blk_elmat%blocks(1,2)%data = 0.0_rp
  !   blk_elmat%blocks(2,1)%data = 0.0_rp
  !   blk_elmat%blocks(2,2)%data = 0.0_rp
  !   if(prob%kfl_stab==2) then
  !      blk_elmat%blocks(1,3)%data = 0.0_rp
  !      blk_elmat%blocks(2,3)%data = 0.0_rp
  !      blk_elmat%blocks(3,1)%data = 0.0_rp
  !      blk_elmat%blocks(3,2)%data = 0.0_rp
  !      blk_elmat%blocks(3,3)%data = 0.0_rp
  !   elseif(prob%kfl_stab==3) then
  !      blk_elmat%blocks(1,3)%data = 0.0_rp
  !      blk_elmat%blocks(1,4)%data = 0.0_rp
  !      blk_elmat%blocks(2,3)%data = 0.0_rp
  !      blk_elmat%blocks(2,4)%data = 0.0_rp
  !      blk_elmat%blocks(3,1)%data = 0.0_rp
  !      blk_elmat%blocks(3,2)%data = 0.0_rp
  !      blk_elmat%blocks(3,3)%data = 0.0_rp
  !      blk_elmat%blocks(3,4)%data = 0.0_rp
  !      blk_elmat%blocks(4,1)%data = 0.0_rp
  !      blk_elmat%blocks(4,2)%data = 0.0_rp
  !      blk_elmat%blocks(4,3)%data = 0.0_rp
  !      blk_elmat%blocks(4,4)%data = 0.0_rp
  !   end if

  !   ! Select integration subroutine
  !   if(prob%kfl_matr==200) then
  !      call nsi_iss_element_mat_massu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
  !   elseif(prob%kfl_matr==201) then
  !      call nsi_iss_element_mat_diffu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
  !   elseif(prob%kfl_matr==202) then
  !      call nsi_iss_element_mat_conve(prob,elem,nnode(1),blk_elmat%blocks(1,1))
  !   elseif(prob%kfl_matr==203) then
  !      call nsi_iss_element_mat_gradp(prob,elem,nnode,blk_elmat%blocks(1,2))
  !   elseif(prob%kfl_matr==204) then
  !      call nsi_iss_element_mat_react(prob,elem,nnode(1),blk_elmat%blocks(1,1))
  !   elseif(prob%kfl_matr==206) then
  !      call nsi_iss_element_mat_divuq(prob,elem,nnode,blk_elmat%blocks(2,1))
  !   elseif(prob%kfl_matr==207) then
  !      call nsi_iss_element_mat_ossuu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
  !   elseif(prob%kfl_matr==208) then
  !      call nsi_iss_element_mat_ossux(prob,elem,nnode(1),blk_elmat%blocks(1,3))
  !   elseif(prob%kfl_matr==212) then
  !      call nsi_iss_element_mat_massp(prob,elem,nnode(2),blk_elmat%blocks(2,2))
  !   elseif(prob%kfl_matr==214) then
  !      call nsi_iss_element_mat_laplp(prob,elem,nnode(2),blk_elmat%blocks(2,2))
  !   elseif(prob%kfl_matr==216) then
  !      call nsi_iss_element_mat_ossxu(prob,elem,nnode(1),blk_elmat%blocks(3,1))
  !   elseif(prob%kfl_matr==217) then
  !      call nsi_iss_element_mat_ossxx(prob,elem,nnode(1),blk_elmat%blocks(3,3))
  !   end if

  !   ! Permute & reblock to monolithic emat & evec
  !   nd = size(elem%p_nod,1)
  !   nv = size(elem%iv,1)
  !   mn = size(elem%jvars,1)
  !   call block_elmat_permute_reblock(vars_s%ib,vars_s%jb,blk_elmat,elem%nint,nnode,nv,elem%iv, &
  !        &                           nd,elem%p_nod,mn,elem%jvars,elem%p_mat)

  !   ! Deallocate locals
  !   if(prob%kfl_stab>=2) then
  !      call memfree(nnode_oss,__FILE__,__LINE__)
  !      call memfree(ldofs_oss,__FILE__,__LINE__)
  !   end if
  !   call fem_blocks_free(vars_s)
  !   call block_elmat_free(blk_elmat)
    
  ! end subroutine nsi_iss_element_mat

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_massu(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the (u,v) term.                                                   !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus,idime,inode,jnode
  !   real(rp)    :: dvolu
    
  !   ! Initialization
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      ! ( u , v )
  !      if(prob%kfl_lump==0) then
  !         do jnode=1,nnode
  !            do inode=1,nnode
  !               do idime=1,prob%ndime
  !                  emat%data(idime,idime,inode,jnode) = emat%data(idime,idime,inode,jnode) + dvolu * &
  !                       &                              (el%integ(1)%p%uint_phy%shape(inode,igaus) *  &
  !                       &                               el%integ(1)%p%uint_phy%shape(jnode,igaus))
  !               end do
  !            end do
  !         end do
  !      ! Lumped mass matrix
  !      else if(prob%kfl_lump==1) then
  !         do jnode=1,nnode
  !            do inode=1,nnode
  !               do idime=1,prob%ndime
  !                  emat%data(idime,idime,inode,inode) = emat%data(idime,idime,inode,inode) + dvolu * &
  !                       &                              (el%integ(1)%p%uint_phy%shape(inode,igaus) *  &
  !                       &                               el%integ(1)%p%uint_phy%shape(jnode,igaus))
  !               end do
  !            end do
  !         end do
  !      end if
  !   end do
   
  ! end subroutine nsi_iss_element_mat_massu

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_diffu(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the mu*(grad u,grad v) term.                                      !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus
  !   real(rp)    :: dvolu
  !   real(rp)    :: work(4),elmuv(nnode,nnode)
    
  !   ! Initialization
  !   emat%data = 0.0_rp
  !   elmuv=0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      ! mu*(grad u,grad v)
  !      call elmvis_gal(dvolu,prob%diffu,el%integ(1)%p%uint_phy%deriv(:,:,igaus),prob%ndime, &
  !           &          nnode,elmuv,work)
  !      ! Add cross terms for symmetric grad
  !      if(prob%kfl_symg==1) then
  !         call elmvis_gal_sym(dvolu,prob%diffu,el%integ(1)%p%uint_phy%deriv(:,:,igaus),prob%ndime, &
  !              &              nnode,emat%data,work)
  !      end if
  !   end do

  !   ! Assembly elmuv to emat
  !   call elmuv2elmat(prob%ndime,nnode,elmuv,emat%data)
   
  ! end subroutine nsi_iss_element_mat_diffu

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_conve(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the (v, a·grad u) term.                                           !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus,inode,idime,jnode,jdime
  !   real(rp)    :: elvel(prob%ndime,nnode)
  !   real(rp)    :: gpvel(prob%ndime,el%integ(1)%p%quad%ngaus)
  !   real(rp)    :: dtinv,dvolu,react
  !   real(rp)    :: agran(nnode)  ! a.grad N
  !   real(rp)    :: work(4),elmuv(nnode,nnode)

  !   ! Interpolation operations for velocity
  !   call elem2var(prob%ndofn,nnode,1,prob%ndime,el%unkno(:,:,1),elvel)
  !   call interpolate(prob%ndime,nnode,el%integ(1)%p%quad%ngaus,el%integ(1)%p%uint_phy%shape,elvel,gpvel)
  
  !   ! Advection velocity
  !   if(prob%kfl_conv==0) then
  !      gpvel = 0.0_rp
  !   end if

  !   ! Initialization
  !   dtinv = 0.0_rp
  !   react = 0.0_rp
  !   elmuv = 0.0_rp
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      do inode=1,nnode
  !         agran(inode)=0.0_rp
  !         do idime=1,prob%ndime
  !            agran(inode) = agran(inode) + &
  !                 &         gpvel(idime,igaus)*el%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
  !         end do
  !      end do

  !      if(prob%kfl_skew==0) then
  !         ! (v, a·grad u)
  !         call elmbuv_gal(dvolu,react,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),agran,nnode,elmuv,work)
  !      elseif(prob%kfl_skew==1) then
  !          ! 1/2(v, a·grad u) - 1/2(u,a·grad v)
  !         call elmbuv_gal_skew1(dvolu,react,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),agran,nnode, &
  !              &                elmuv,work)
  !      end if
  !   end do

  !   ! Assembly elmuv to emat
  !   call elmuv2elmat(prob%ndime,nnode,elmuv,emat%data)

  ! end subroutine nsi_iss_element_mat_conve

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_gradp(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the - (div v, p) term.                                            !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode(2)
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus
  !   real(rp)    :: dvolu
  !   real(rp)    :: work(4)
    
  !   ! Initialization
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      ! - (div v, p)
  !      call elmbpv_gal_div_iss(dvolu,el%integ(2)%p%uint_phy%shape(:,igaus),                          &
  !           &                  el%integ(1)%p%uint_phy%deriv(:,:,igaus),prob%ndime,nnode(1),nnode(2), &
  !           &                  emat%data,work)
  !   end do
   
  ! end subroutine nsi_iss_element_mat_gradp

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_react(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the s*(v,u) term.                                                 !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus
  !   real(rp)    :: dtinv,dvolu,react
  !   real(rp)    :: agran(nnode)  ! a.grad N
  !   real(rp)    :: work(4),elmuv(nnode,nnode)

  !   ! Initialization
  !   dtinv = 0.0_rp
  !   agran = 0.0_rp
  !   elmuv = 0.0_rp
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu=el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      ! s*(v,u)
  !      call elmbuv_gal(dvolu,react,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),agran,nnode,elmuv,work)
  !   end do

  !   ! Assembly elmuv to emat
  !   call elmuv2elmat(prob%ndime,nnode,elmuv,emat%data)

  ! end subroutine nsi_iss_element_mat_react

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_divuq(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the (div u, q) term.                                              !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode(2)
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus
  !   real(rp)    :: dvolu
  !   real(rp)    :: work(4)
    
  !   ! Initialization
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)
       
  !      ! (div u, q)
  !      call elmbuq_gal_div_iss(dvolu,el%integ(2)%p%uint_phy%shape(:,igaus),                          &
  !           &                  el%integ(1)%p%uint_phy%deriv(:,:,igaus),prob%ndime,nnode(1),nnode(2), &
  !           &                  emat%data,work)
  !   end do

  ! end subroutine nsi_iss_element_mat_divuq
 
  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_massp(prob,el,nnode,emat)
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus,inode,jnode
  !   real(rp)    :: dvolu

  !   ! Initialization
  !   emat%data = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)
  !      if (prob%kfl_lump==0) then
  !         ! ( p , q )
  !         do jnode=1,nnode
  !            do inode=1,nnode
  !               emat%data(1,1,inode,jnode) = emat%data(1,1,inode,jnode) + dvolu * &
  !                    &                       el%integ(2)%p%uint_phy%shape(inode,igaus) *  &
  !                    &                       el%integ(2)%p%uint_phy%shape(jnode,igaus)
  !            end do
  !         end do

  !      else
  !         ! ( p , q )
  !         do jnode=1,nnode
  !            do inode=1,nnode
  !               emat%data(1,1,inode,inode) = emat%data(1,1,inode,inode) + dvolu * &
  !                    &                       el%integ(2)%p%uint_phy%shape(inode,igaus) *  &
  !                    &                       el%integ(2)%p%uint_phy%shape(jnode,igaus)
  !            end do
  !         end do
  !      end if
  !   end do

  ! end subroutine nsi_iss_element_mat_massp

  ! !=================================================================================================
  ! subroutine nsi_iss_element_mat_laplp(prob,el,nnode,emat)
  !   !----------------------------------------------------------------------------------------------!
  !   !   This subroutine computes the (grad p , grad q) term.                                       !
  !   !----------------------------------------------------------------------------------------------!
  !   implicit none
  !   type(nsi_cg_iss_approximation), intent(in)    :: prob
  !   type(fem_element)    , intent(in)    :: el
  !   integer(ip)          , intent(in)    :: nnode
  !   type(elmat)          , intent(inout) :: emat
  !   ! Locals
  !   integer(ip) :: igaus
  !   real(rp)    :: dvolu
  !   real(rp)    :: work(4),elmuv(nnode,nnode)
    
  !   ! Initialization
  !   emat%data = 0.0_rp
  !   elmuv = 0.0_rp

  !   ! Loop on Gauss points
  !   do igaus=1,el%integ(1)%p%quad%ngaus
  !      dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)
       
  !      ! (grad p, grad q)
  !      call elmvis_gal(dvolu,1.0_rp,el%integ(2)%p%uint_phy%deriv(:,:,igaus),prob%ndime,nnode, &
  !           &                  elmuv,work)
  !   end do

  !   ! Assembly elmuv to emat
  !   call elmuv2elmat(1,nnode,elmuv,emat%data)

  ! end subroutine nsi_iss_element_mat_laplp

end module nsi_cg_iss_names
