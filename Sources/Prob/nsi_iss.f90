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
 implicit none
# include "debug.i90"

 private
 
 ! INF-SUP STABLE (iss) NAVIER-STOKES Problem data
 type, extends(discrete_problem) :: nsi_cg_iss_approximation
    type(nsi_problem), pointer :: physics
    integer(ip) ::   & 
         kfl_1vec,   & ! Flag for 1vec integration
         kfl_mtvc,   & ! Flag for matvec integration
         kfl_matr,   & ! Flag for matrix integration
         kfl_real,   & ! Flag for real integration
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
    ! type(array_rp2), allocatable :: &
    !      vorti(:),      & ! Vorticity field
    !      dtunk(:)         ! Time derivative
     contains
      procedure :: create => nsi_create
      procedure :: matvec => nsi_matvec
 end type nsi_cg_iss_approximation

 type block_elmat
    integer(ip)              :: nblocks = 0
    real(rp)   , allocatable :: blocks(:,:)
 end type block_elmat

 ! Types
 public :: nsi_cg_iss_approximation

 ! Functions
 ! public :: nsi_iss_create, nsi_iss_free, nsi_iss_element_matvec, nsi_iss_element_real, &
 !      &    nsi_iss_result_write, nsi_iss_element_mat, nsi_iss_element_1vec,            &
 !      &    nsi_iss_create_vor, nsi_iss_free_vor, nsi_iss_update_vorti,                 &
 !      &    nsi_iss_create_tder, nsi_iss_free_tder, nsi_iss_create_oss, nsi_iss_face_mat
 
contains

  !=================================================================================================
  subroutine nsi_create(approx,prob,l2g)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine contains definitions of the Navier-Stokes problem approximed by a stable   !
    !   finite element formulation with inf-sup stable elemets.                                    !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_approximation)       , intent(out) :: approx
    class(physical_problem), target, intent(in)  :: prob
    integer(ip), intent(in), optional :: l2g(:)
    !type(nsi_problem), target, intent(in)  :: prob
    integer(ip) :: i

    select type (prob)
    type is(nsi_problem)
       approx%physics => prob
       !approx%nsi => prob
    class default
       check(.false.)
    end select

    approx%nvars = prob%ndime+1
    call memalloc(approx%nvars,approx%l2g_var,__FILE__,__LINE__)
    if ( present(l2g) ) then
       assert ( size(l2g) == approx%nvars )
       approx%l2g_var = l2g
    else
       do i = 1,approx%nvars
          approx%l2g_var(i) = i
       end do
    end if

    ! Flags
    approx%kfl_1vec = 1 ! Integrate Full OSS
    approx%kfl_mtvc = 1 ! Integrate nsi_element_matvec
    approx%kfl_matr = 1 ! Integrate nsi_element_mat
    approx%kfl_real = 0 ! Integrate errornorm
    approx%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)

    ! Problem variables
    approx%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter tau_m
    approx%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter tau_m
    approx%ktauc  = 0.0_rp  ! Constant multiplying stabilization parameter tau_c

    ! Time integration variables
    approx%dtinv  = 1.0_rp ! Inverse of time step
    approx%ctime  = 0.0_rp ! Current time
    approx%tdimv  =  2     ! Number of temporal steps stored for velocity
    approx%tdimp  =  2     ! Number of temporal steps stored for pressure
    
  end subroutine nsi_create

  !=================================================================================================
  subroutine nsi_matvec(approx,start,elem)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental matrix-vector integration selection.                !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_approximation), intent(inout) :: approx
    integer(ip)                    , intent(in)    :: start(:)
    type(fem_element)              , intent(inout) :: elem
    ! Locals
    ! type(block_elmat)        :: blk_elmat
    ! type(block_elvec)        :: blk_elvec
    ! type(fem_blocks)         :: vars_s
    ! integer(ip)              :: ldofs(2),nblks,nnode(2),nd,nv,mn
    ! integer(ip), allocatable :: nnode_oss(:),ldofs_oss(:)

    ! Checks
    check(elem%f_inf(1)%p%order > elem%f_inf(approx%physics%ndime+1)%p%order)
    check(elem%integ(1)%p%quad%ngaus == elem%integ(approx%physics%ndime+1)%p%quad%ngaus)
    
    ! ! Allocate blk_elmat & blk_elvec locals
    ! nnode(1) = elem%f_inf(1)%p%nnode
    ! nnode(2) = elem%f_inf(2)%p%nnode
    ! if(prob%kfl_stab==2) then
    !    nblks=3
    !    call memalloc(nblks,nnode_oss,__FILE__,__LINE__)
    !    call memalloc(nblks,ldofs_oss,__FILE__,__LINE__)
    !    ldofs_oss(1) = prob%ndime; ldofs_oss(2)=1; ldofs_oss(3)=prob%ndime;
    !    nnode_oss(1) = elem%f_inf(1)%p%nnode
    !    nnode_oss(2) = elem%f_inf(2)%p%nnode
    !    nnode_oss(3) = elem%f_inf(1)%p%nnode
    !    call fem_blocks_alloc(scalar,nblks,ldofs_oss,vars_s)
    !    call block_elmat_alloc(vars_s,nnode_oss,blk_elmat)
    !    call block_elvec_alloc(vars_s,nnode_oss,blk_elvec)
    ! elseif(prob%kfl_stab==3) then
    !    nblks=4
    !    call memalloc(nblks,nnode_oss,__FILE__,__LINE__)
    !    call memalloc(nblks,ldofs_oss,__FILE__,__LINE__)
    !    ldofs_oss(1) = prob%ndime; ldofs_oss(2)=1; ldofs_oss(3)=prob%ndime; ldofs_oss(4)=prob%ndime
    !    nnode_oss(1) = elem%f_inf(1)%p%nnode
    !    nnode_oss(2) = elem%f_inf(2)%p%nnode
    !    nnode_oss(3) = elem%f_inf(1)%p%nnode
    !    nnode_oss(4) = elem%f_inf(1)%p%nnode
    !    call fem_blocks_alloc(scalar,nblks,ldofs_oss,vars_s)
    !    call block_elmat_alloc(vars_s,nnode_oss,blk_elmat)
    !    call block_elvec_alloc(vars_s,nnode_oss,blk_elvec)
    ! else
    !    nblks=2
    !    ldofs(1) = prob%ndime; ldofs(2)=1
    !    call fem_blocks_alloc(scalar,nblks,ldofs,vars_s)
    !    call block_elmat_alloc(vars_s,nnode,blk_elmat)
    !    call block_elvec_alloc(vars_s,nnode,blk_elvec)
    ! end if

    ! Initialize matrix and vector
    elem%p_mat%a = 0.0_rp
    elem%p_vec%a = 0.0_rp

    ! ! Initialize to zero
    ! blk_elvec%blocks(1)%data   = 0.0_rp
    ! blk_elvec%blocks(2)%data   = 0.0_rp
    ! blk_elmat%blocks(1,1)%data = 0.0_rp
    ! blk_elmat%blocks(1,2)%data = 0.0_rp
    ! blk_elmat%blocks(2,1)%data = 0.0_rp
    ! blk_elmat%blocks(2,2)%data = 0.0_rp
    ! if(prob%kfl_stab==2) then
    !    blk_elvec%blocks(3)%data   = 0.0_rp
    !    blk_elmat%blocks(1,3)%data = 0.0_rp
    !    blk_elmat%blocks(2,3)%data = 0.0_rp
    !    blk_elmat%blocks(3,1)%data = 0.0_rp
    !    blk_elmat%blocks(3,2)%data = 0.0_rp
    !    blk_elmat%blocks(3,3)%data = 0.0_rp
    ! elseif(prob%kfl_stab==3) then
    !    blk_elvec%blocks(3)%data   = 0.0_rp
    !    blk_elvec%blocks(4)%data   = 0.0_rp
    !    blk_elmat%blocks(1,3)%data = 0.0_rp
    !    blk_elmat%blocks(1,4)%data = 0.0_rp
    !    blk_elmat%blocks(2,3)%data = 0.0_rp
    !    blk_elmat%blocks(2,4)%data = 0.0_rp
    !    blk_elmat%blocks(3,1)%data = 0.0_rp
    !    blk_elmat%blocks(3,2)%data = 0.0_rp
    !    blk_elmat%blocks(3,3)%data = 0.0_rp
    !    blk_elmat%blocks(3,4)%data = 0.0_rp
    !    blk_elmat%blocks(4,1)%data = 0.0_rp
    !    blk_elmat%blocks(4,2)%data = 0.0_rp
    !    blk_elmat%blocks(4,3)%data = 0.0_rp
    !    blk_elmat%blocks(4,4)%data = 0.0_rp
    ! end if
    
    ! Select integration subroutine
    if(approx%kfl_mtvc==1) then
       call nsi_iss_element_matvec_stdr(approx,elem,start)
    ! elseif(prob%kfl_mtvc==200) then
    !    call nsi_iss_element_mat_massu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
    ! elseif(prob%kfl_mtvc==201) then
    !    call nsi_iss_element_mat_diffu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
    ! elseif(prob%kfl_mtvc==202) then
    !    call nsi_iss_element_mat_conve(prob,elem,nnode(1),blk_elmat%blocks(1,1))
    ! elseif(prob%kfl_mtvc==203) then
    !    call nsi_iss_element_mat_gradp(prob,elem,nnode,blk_elmat%blocks(1,2))
    ! elseif(prob%kfl_mtvc==204) then
    !    call nsi_iss_element_mat_react(prob,elem,nnode(1),blk_elmat%blocks(1,1))
    ! elseif(prob%kfl_mtvc==206) then
    !    call nsi_iss_element_mat_divuq(prob,elem,nnode,blk_elmat%blocks(2,1))
    ! elseif(prob%kfl_mtvc==207) then
    !    call nsi_iss_element_mat_ossuu(prob,elem,nnode(1),blk_elmat%blocks(1,1))
    ! elseif(prob%kfl_mtvc==208) then
    !    call nsi_iss_element_mat_ossux(prob,elem,nnode(1),blk_elmat%blocks(1,3))
    ! elseif(prob%kfl_mtvc==212) then
    !    call nsi_iss_element_mat_massp(prob,elem,nnode(2),blk_elmat%blocks(2,2))
    ! elseif(prob%kfl_mtvc==214) then
    !    call nsi_iss_element_mat_laplp(prob,elem,nnode(2),blk_elmat%blocks(2,2))
    ! elseif(prob%kfl_mtvc==216) then
    !    call nsi_iss_element_mat_ossxu(prob,elem,nnode(1),blk_elmat%blocks(3,1))
    ! elseif(prob%kfl_mtvc==217) then
    !    call nsi_iss_element_mat_ossxx(prob,elem,nnode(1),blk_elmat%blocks(3,3))
    end if

    ! ! Permute & reblock to monolithic emat & evec
    ! nd = size(elem%p_nod,1)
    ! nv = size(elem%iv,1)
    ! mn = size(elem%jvars,1)
    ! call block_elmat_permute_reblock(vars_s%ib,vars_s%jb,blk_elmat,elem%nint,nnode,nv,elem%iv, &
    !      &                           nd,elem%p_nod,mn,elem%jvars,elem%p_mat)
    ! call block_elvec_permute_reblock(vars_s%ib,vars_s%jb,blk_elvec,elem%nint,nnode,nv,elem%iv, &
    !      &                           nd,elem%p_nod,mn,elem%jvars,elem%p_vec)

    ! ! Deallocate locals
    ! if(prob%kfl_stab>=2) then
    !    call memfree(nnode_oss,__FILE__,__LINE__)
    !    call memfree(ldofs_oss,__FILE__,__LINE__)
    ! end if
    ! call fem_blocks_free(vars_s)
    ! call block_elmat_free(blk_elmat)
    ! call block_elvec_free(blk_elvec)
    
  end subroutine nsi_matvec

  !=================================================================================================
  subroutine nsi_iss_element_matvec_stdr(approx,el,start)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine computes the standard GALERKIN elemental matrix-vector integration.        !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(nsi_cg_iss_approximation), intent(in)    :: approx
    type(fem_element)              , intent(inout) :: el
    integer(ip)                    , intent(in)    :: start(:)
    ! Locals
    integer(ip)  :: igaus,idime,inode,jdime,jnode
    integer(ip)  :: ngaus,ndime,nnodu,nnodp
    real(rp)     :: ctime,dtinv,dvolu,diffu,react
    real(rp)     :: work(4)
    real(rp)     :: force(approx%physics%ndime)
    real(rp)     :: agran(el%integ(1)%p%uint_phy%nnode)
    real(rp)     :: elmuv(el%integ(1)%p%uint_phy%nnode,el%integ(1)%p%uint_phy%nnode)
    real(rp)     :: elmat_uv(approx%physics%ndime,approx%physics%ndime,el%integ(1)%p%uint_phy%nnode, &
         &                   el%integ(1)%p%uint_phy%nnode)
    real(rp)     :: elmat_up(approx%physics%ndime,1,el%integ(1)%p%uint_phy%nnode,el%integ(1)%p%uint_phy%nnode)
    real(rp)     :: elmat_pu(1,approx%physics%ndime,el%integ(1)%p%uint_phy%nnode,el%integ(1)%p%uint_phy%nnode)
    real(rp)     :: elvec_u(approx%physics%ndime,el%integ(1)%p%uint_phy%nnode)
    ! real(rp)    :: parv(30),parp(10),part(3),part_p(3)
    type(vector) :: gpvel, gpveln

    ! Unpack variables
    ndime = approx%physics%ndime
    nnodu = el%integ(1)%p%uint_phy%nnode
    nnodp = el%integ(1)%p%uint_phy%nnode
    ngaus = el%integ(1)%p%quad%ngaus
    diffu = approx%physics%diffu
    react = approx%physics%react
    dtinv = approx%dtinv
    
    ! Interpolation operations for velocity
    call create_vector (approx%physics, 1, el%integ, gpvel)
    call create_vector (approx%physics, 1, el%integ, gpveln)
    call interpolation (el%unkno, 1, 1, el%integ, gpvel)
    if(dtinv == 0.0_rp) then
       call interpolation (el%unkno, 1, 2, el%integ, gpveln)
    else
       gpveln%a = 0.0_rp
    end if

    ! Set real time
    if(approx%kfl_thet==0) then        ! BE
       ctime = approx%ctime
    elseif(approx%kfl_thet==1) then    ! CN
       ctime = approx%ctime - 1.0_rp/dtinv
    end if
    
    ! Initializations
    work     = 0.0_rp
    force    = 0.0_rp
    elmuv    = 0.0_rp
    elmat_uv = 0.0_rp
    elmat_up = 0.0_rp
    elmat_pu = 0.0_rp
    elvec_u  = 0.0_rp
    ! parv     = 0.0_rp
    ! parp     = 0.0_rp
    ! part     = 0.0_rp    

    ! Loop on Gauss points
    do igaus = 1,el%integ(1)%p%quad%ngaus
       dvolu = el%integ(1)%p%quad%weight(igaus)*el%integ(1)%p%femap%detjm(igaus)

  !      ! Analytical solution
  !      if(prob%case_veloc>0.and.prob%case_press>0) then         
  !         ! Velocity norm at gauss point
  !         gpvno=0.0_rp
  !         do idime = 1,prob%ndime
  !            gpvno = gpvno + gpvel(idime,igaus)*gpvel(idime,igaus)
  !         end do
  !         gpvno = sqrt(gpvno)

  !         ! Interpolation of coords
  !         call analytical_field(prob%case_veloc,prob%ndime,el%integ(1)%p%femap%clocs(:,igaus),ctime,parv)
  !         call analytical_field(prob%case_press,prob%ndime,el%integ(1)%p%femap%clocs(:,igaus),ctime,parp)
  !         call analytical_field(prob%case_tempo,prob%ndime,el%integ(1)%p%femap%clocs(:,igaus),ctime,part)
  !         call analytical_field(prob%case_t_pre,prob%ndime,el%integ(1)%p%femap%clocs(:,igaus),ctime,part_p)
  !         if(gpvno.lt.1e-4) then
  !            conve = 0
  !         else
  !            conve = prob%kfl_conv
  !         end if
  !         call nsi_basic_force(prob%ndime,prob%diffu,prob%react,conve,prob%kfl_symg,prob%case_tempo, &
  !              &               parv,parp,part,prob%gravi,force,part_p)
  !      else
  !         force=0.0_rp
  !      end if

       ! Gravity field
       do idime=1,approx%physics%ndime
          force(idime) = force(idime) + approx%physics%gravi(idime)
       end do

       ! Auxiliar variables
       if(approx%physics%kfl_conv.ne.0) then
          do inode = 1,nnodu
             agran(inode) = 0.0_rp
             do idime = 1,ndime
                agran(inode) = agran(inode) + &
                     &         gpvel%a(idime,igaus)*el%integ(1)%p%uint_phy%deriv(idime,inode,igaus)
             end do
          end do
       end if

       ! Computation of elemental terms
       ! ------------------------------
       ! Block U-V
       ! mu * ( grad u, grad v )
       call elmvis_gal(dvolu,diffu,el%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,elmuv,work)

       ! Add cross terms for symmetric grad
       if(approx%physics%kfl_symg==1) then
          call elmvis_gal_sym(dvolu,diffu,el%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu, &
               &              elmat_uv,work)
       end if
       if(approx%physics%kfl_skew==0) then
          ! (v, a·grad u) + s*(v,u) + (v, u/dt)
          call elmbuv_gal(dvolu,react,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu,elmuv,work)
       elseif(approx%physics%kfl_skew==1) then
          ! 1/2(v, a·grad u) - 1/2(u,a·grad v) + s*(v,u) + (v, u/dt)
          call elmbuv_gal_skew1(dvolu,react,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),agran,nnodu, &
               &                elmuv,work)
       end if

       ! Block P-V
       ! - ( div v, p )
       call elmbpv_gal_div_iss(dvolu,el%integ(2)%p%uint_phy%shape(:,igaus),               &
            &                  el%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_up,work)

       ! Block U-Q
       ! ( div u, q )
       call elmbuq_gal_div_iss(dvolu,el%integ(2)%p%uint_phy%shape(:,igaus),               &
            &                  el%integ(1)%p%uint_phy%deriv(:,:,igaus),ndime,nnodu,nnodp, &
            &                  elmat_pu,work)

       ! RHS: Block U
       ! ( v, f ) + ( v, u_n/dt )
       call elmrhu_gal(dvolu,dtinv,el%integ(1)%p%uint_phy%shape(:,igaus),gpveln%a(:,igaus), &
            &          force,nnodu,ndime,elvec_u,work)

    end do

    ! Assembly to elemental p_mat and p_vec
    do inode=1,nnodu
       do idime=1,ndime
          do jnode=1,nnodu
             do jdime=1,ndime
                ! Block V-U
                el%p_mat%a(start(idime)+inode-1,start(jdime)+jnode-1) = &
                     & el%p_mat%a(start(idime)+inode-1,start(jdime)+jnode-1) + elmat_uv(idime,jdime,inode,jnode)
             end do    
             ! Block V-U (diag)
             el%p_mat%a(start(idime)+inode-1,start(idime)+jnode-1) = &
                  & el%p_mat%a(start(idime)+inode-1,start(idime)+jnode-1) + elmuv(inode,jnode)
          end do
          do jnode=1,nnodp
             ! Block V-P
             el%p_mat%a(start(idime)+inode-1,start(ndime+1)+jnode-1) = &
                  & el%p_mat%a(start(idime)+inode-1,start(ndime+1)+jnode-1) + elmat_up(idime,1,inode,jnode)
          end do
          ! Block U
          el%p_vec%a(start(idime)+inode-1) = el%p_vec%a(start(idime)+inode-1) + elvec_u(idime,inode)
       end do
    end do
    do inode=1,nnodp
       do jdime=1,ndime
          do jnode=1,nnodu
             ! Block Q-U
             el%p_mat%a(start(ndime+1)+inode-1,start(jdime)+jnode-1) = &
                  & el%p_mat%a(start(ndime+1)+inode-1,start(jdime)+jnode-1) + elmat_uv(1,jdime,inode,jnode)
          end do
       end do
    end do
    
  end subroutine nsi_iss_element_matvec_stdr

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
