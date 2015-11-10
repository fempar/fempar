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
module cdr_stabilized_continuous_Galerkin_names
  use types_names
  use memor_names
  use problem_names
  use cdr_names
  use element_fields_names
  use element_tools_names
  use finite_element_names
  implicit none
# include "debug.i90"
  private 

  type, extends(discrete_problem_t) :: cdr_discrete_t
    integer(ip) ::   & 
         kfl_thet,   & ! Flag for theta-method (0=BE, 1=CN)
         kfl_lump,   & ! Flag for lumped mass submatrix
         kfl_stab,   & ! Flag for stabilization of the convective term (Off=0; OSS=2)
         kfl_proj,   & ! Flag for Projections weighted with tau's (On=1, Off=0)
         tdimv         ! Number of temporal steps stored
    real(rp) ::      &
         dtinv,      & ! Inverse of time step
         ctime,      & ! Current time 
         k1tau,      & ! C1 constant on stabilization parameter 
         k2tau         ! C2 constant on stabilization parameter 
    contains
      procedure :: create => cdr_create_discrete
   end type cdr_discrete_t

   type, extends(discrete_integration_t) :: cdr_approximation_t
      type(cdr_discrete_t), pointer :: discret
      type(cdr_problem_t) , pointer :: physics
    contains
      procedure :: create  => cdr_matvec_create
      procedure :: compute => cdr_matvec
      procedure :: free    => cdr_matvec_free
   end type cdr_approximation_t   

   type, extends(discrete_integration_t) :: cdr_nonlinear_t
      type(cdr_discrete_t), pointer :: discret
      type(cdr_problem_t) , pointer :: physics
    contains
      procedure :: create  => cdr_nonlinear_create
      procedure :: compute => cdr_nonlinear
      procedure :: free    => cdr_nonlinear_free
   end type cdr_nonlinear_t

  ! Unkno components parameter definition
  integer(ip), parameter :: current   = 1
  integer(ip), parameter :: prev_iter = 2
  integer(ip), parameter :: prev_step = 3


 public :: cdr_approximation_t, cdr_discrete_t, cdr_nonlinear_t

contains

  !=================================================================================================
  subroutine cdr_create_discrete(discret,physics,l2g)
    !---------------------------------------------------------------------------!
    !   This subroutine contains definitions of the cdr problem approximed by   !
    !   a stabilised finite element formulation with inf-sup stable elemets.    !
    !---------------------------------------------------------------------------!
    implicit none
    class(cdr_discrete_t)    , intent(out) :: discret
    class(physical_problem_t), intent(in)  :: physics
    integer(ip), optional  , intent(in)  :: l2g(:)
    ! Locals
    integer(ip) :: i

    ! Flags
    discret%kfl_thet = 0 ! Theta-method time integration (BE=0, CN=0)
    discret%kfl_stab = 0 ! Stabilization of convective term (0: Off, 2: OSS)
    discret%kfl_proj = 0 ! Projections weighted with tau's (On=1, Off=0)

    ! Problem variables
    discret%k1tau  = 4.0_rp  ! C1 constant on stabilization parameter
    discret%k2tau  = 2.0_rp  ! C2 constant on stabilization parameter

    ! Time integration variables
    discret%dtinv  = 1.0_rp ! Inverse of time step
    discret%ctime  = 0.0_rp ! Current time
    discret%tdimv  =  2     ! Number of temporal steps stored

    discret%nvars = 1
    call memalloc(discret%nvars,discret%l2g_var,__FILE__,__LINE__)
    if ( present(l2g) ) then
       assert ( size(l2g) == discret%nvars )
       discret%l2g_var = l2g
    else
       do i = 1,discret%nvars
          discret%l2g_var(i) = i
       end do
    end if
    
  end subroutine cdr_create_discrete

  !=================================================================================================
  subroutine cdr_matvec_create( this, physics, discret )
    implicit none
    class(cdr_approximation_t)       , intent(inout) :: this
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret

    select type (physics)
    type is(cdr_problem_t)
       this%physics => physics
    class default
       check(.false.)
    end select 
    select type (discret)
    type is(cdr_discrete_t)
       this%discret => discret
    class default
       check(.false.)
    end select

    ! Allocate working variables
    call memalloc(1,this%working_vars,__FILE__,__LINE__)
    this%working_vars(1) = 1

    this%domain_dimension = 3 

  end subroutine cdr_matvec_create

  !=================================================================================================
  subroutine cdr_nonlinear_create( this, physics, discret )
    implicit none
    class(cdr_nonlinear_t)       , intent(inout) :: this
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret

    select type (physics)
    type is(cdr_problem_t)
       this%physics => physics
    class default
       check(.false.)
    end select 
    select type (discret)
    type is(cdr_discrete_t)
       this%discret => discret
    class default
       check(.false.)
    end select

    ! Allocate working variables
    call memalloc(1,this%working_vars,__FILE__,__LINE__)
    this%working_vars(1) = 1

    this%domain_dimension = 3 

  end subroutine cdr_nonlinear_create

  !=================================================================================================
  subroutine cdr_matvec_free(this)
    implicit none
    class(cdr_approximation_t)  , intent(inout) :: this

    this%physics => null()
    this%discret => null()

    ! Deallocate working variables
    call memfree(this%working_vars,__FILE__,__LINE__)

  end subroutine cdr_matvec_free

  !=================================================================================================
  subroutine cdr_nonlinear_free(this)
    implicit none
    class(cdr_nonlinear_t)  , intent(inout) :: this

    this%physics => null()
    this%discret => null()

    ! Deallocate working variables
    call memfree(this%working_vars,__FILE__,__LINE__)

  end subroutine cdr_nonlinear_free

  !=================================================================================================
  subroutine cdr_matvec(this,finite_element)
    implicit none
    class(cdr_approximation_t), intent(inout) :: this
    type(finite_element_t)    , intent(inout) :: finite_element

    type(basis_function_t) :: u ! Trial
    type(basis_function_t) :: p 
    type(basis_function_t) :: v ! Test
    type(basis_function_t) :: q 

    type(vector_t) :: a
    type(vector_t) :: u_n
    type(scalar_t) :: h,tau

    integer(ip)  :: ndime
    real(rp)     :: dtinv, c1, c2, mu

    integer(ip)  :: idime,igaus,inode,jnode,ngaus,nnode
    real(rp) :: factor


    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp
    ndime = this%physics%ndime


    nnode = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%uint_phy%nlocs

    do igaus = 1,ngaus
       factor = finite_element%integ(1)%p%femap%detjm(igaus) * finite_element%integ(1)%p%quad%weight(igaus)
       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) + factor * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus) * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
          end do
       end do
    end do

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine cdr_matvec

  !=================================================================================================
  subroutine cdr_nonlinear(this,finite_element)
    implicit none
    class(cdr_nonlinear_t), intent(inout) :: this
    type(finite_element_t)    , intent(inout) :: finite_element
    ! Locals
    type(scalar_t) :: force,gpunk
    integer(ip)    :: idime,igaus,inode,jnode,ngaus,nnode,ndime
    real(rp)       :: factor,work(4)

    work = 0.0_rp
    finite_element%p_mat%a = 0.0_rp
    finite_element%p_vec%a = 0.0_rp
    ndime = this%physics%ndime
    nnode = finite_element%integ(1)%p%uint_phy%nnode
    ngaus = finite_element%integ(1)%p%uint_phy%nlocs

    ! Set force term
    call create_scalar(this%physics,1,finite_element%integ,force)
    force%a=0.0_rp
    ! Impose analytical solution
    if(finite_element%p_analytical_code%a(1,1)>0) then 
       call cdr_analytical_force(this%physics,finite_element,0.0_rp,force)
    end if

    ! Unknown interpolation
    if(this%physics%kfl_react>0) then
       call create_scalar(this%physics,1,finite_element%integ,gpunk)
       call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpunk)
    end if
    
    do igaus = 1,ngaus
       factor = finite_element%integ(1)%p%femap%detjm(igaus) * finite_element%integ(1)%p%quad%weight(igaus)

       ! Compute analytical reaction
       if(this%physics%kfl_react>0) then
          call cdr_analytical_reaction(this%physics%ndime,finite_element%integ(1)%p%femap%clocs(:,igaus), &
               &                       gpunk%a(igaus),this%physics%kfl_react,this%physics%react)
       end if

       ! Compute analytical convection
       if (this%physics%kfl_conv>0) then
          call cdr_analytical_convection(this%physics%ndime,                                      &
               &                         finite_element%integ(1)%p%femap%clocs(:,igaus),            &
               &                         this%physics%kfl_conv,this%discret%ctime,              &
               &                         this%physics%convect)
       end if
       
       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                ! nu (grad u, grad v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                     & factor * this%physics%diffu * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus) * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
                ! (b·grad u, v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                     & factor * this%physics%convect(idime) * &
                     & finite_element%integ(1)%p%uint_phy%shape(inode,igaus) * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
             ! react ( u, v)
             finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                  & factor * this%physics%react * &
                  & finite_element%integ(1)%p%uint_phy%shape(inode,igaus) * &
                  & finite_element%integ(1)%p%uint_phy%shape(jnode,igaus)
          end do
          ! Force term (f,v)
          finite_element%p_vec%a(inode) = finite_element%p_vec%a(inode) +  &
                  & factor * force%a(igaus) * &
                  & finite_element%integ(1)%p%uint_phy%shape(inode,igaus)
       end do
    end do
    call memfree(force%a,__FILE__,__LINE__)
    if(this%physics%kfl_react>0) call memfree(gpunk%a,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine cdr_nonlinear

end module cdr_stabilized_continuous_Galerkin_names
