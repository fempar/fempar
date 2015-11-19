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
  use time_integration_names
  use theta_method_names
  implicit none
# include "debug.i90"
  private 

  type, extends(discrete_problem_t) :: cdr_discrete_t
    integer(ip) ::   & 
         kfl_stab    ! Flag for stabilization of the convective term (Off=0; SUPG = 1)
    !real(rp) ::      &
    !     ctime       ! Current time 

    contains
      procedure :: create        => cdr_create_discrete
      procedure :: get_prev_step => cdr_get_prev_step
      procedure :: get_prev_iter => cdr_get_prev_iter
      procedure :: get_current   => cdr_get_current
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

   type, extends(discrete_integration_t) :: cdr_transient_t
      type(cdr_discrete_t), pointer :: discret
      type(cdr_problem_t) , pointer :: physics
      type(theta_method_t), pointer :: time_integ
    contains
      procedure :: create  => cdr_transient_create
      procedure :: compute => cdr_transient
      procedure :: free    => cdr_transient_free
   end type cdr_transient_t

  ! Unkno components parameter definition
  integer(ip), parameter :: current   = 1
  integer(ip), parameter :: prev_iter = 2
  integer(ip), parameter :: prev_step = 3


 public :: cdr_approximation_t, cdr_discrete_t, cdr_nonlinear_t,  cdr_transient_t

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
    discret%kfl_stab = 0 ! Stabilization of convective term (0: Off, 1: SUPG)

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
  integer function cdr_get_current(discret)
    !---------------------------------------------------------------------------!
    !   Give the value of the storing position of the current solution          !
    !---------------------------------------------------------------------------!
    implicit none
    class(cdr_discrete_t)    , intent(in) :: discret
    
    cdr_get_current = current
  end function cdr_get_current

  !=================================================================================================
  integer function cdr_get_prev_iter(discret)
    !---------------------------------------------------------------------------!
    ! Give the value of the storing position of the previous nonlinear iteration!
    !---------------------------------------------------------------------------!
    implicit none
    class(cdr_discrete_t)    , intent(in) :: discret
    
    cdr_get_prev_iter = prev_iter
  end function cdr_get_prev_iter

  !=================================================================================================
  integer function cdr_get_prev_step(discret)
    !---------------------------------------------------------------------------!
    !   Give the value of the storing position of the previous time step        !
    !---------------------------------------------------------------------------!
    implicit none
    class(cdr_discrete_t)    , intent(in) :: discret
    
    cdr_get_prev_step = prev_step
  end function cdr_get_prev_step

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
  subroutine cdr_transient_create( this, physics, discret)
    implicit none
    class(cdr_transient_t)           , intent(inout) :: this
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

    
  end subroutine cdr_transient_create

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
  subroutine cdr_transient_free(this)
    implicit none
    class(cdr_transient_t)  , intent(inout) :: this

    this%physics => null()
    this%discret => null()

    ! Deallocate working variables
    call memfree(this%working_vars,__FILE__,__LINE__)

  end subroutine cdr_transient_free

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
               &                       gpunk%a(igaus),this%physics%kfl_react,this%physics%reaction)
       end if

       ! Compute analytical convection
       if (this%physics%kfl_conv>0) then
          call cdr_analytical_convection(this%physics%ndime,                                        &
               &                         finite_element%integ(1)%p%femap%clocs(:,igaus),            &
               &                         this%physics%kfl_conv,0.0_rp,this%physics%convection)
       end if
       
       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                ! nu (grad u, grad v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                     & factor * this%physics%diffusion * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus) * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
                ! (b·grad u, v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                     & factor * this%physics%convection(idime) * &
                     & finite_element%integ(1)%p%uint_phy%shape(inode,igaus) * &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
             ! react ( u, v)
             finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +  &
                  & factor * this%physics%reaction * &
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

  !==================================================================================================
  subroutine cdr_transient(this,finite_element)
    implicit none
    class(cdr_transient_t), intent(inout) :: this
    type(finite_element_t), intent(inout) :: finite_element
    ! Locals
    type(scalar_t) :: force,gpunk, prev_unkno
    integer(ip)    :: idime,igaus,inode,jnode,ngaus,nnode,ndime
    real(rp)       :: factor,work(4),tau

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
       call cdr_analytical_force(this%physics,finite_element,this%time_integ%ctime,force)
    end if

    ! Unknown interpolation
    if(this%physics%kfl_react>0) then
       call create_scalar(this%physics,1,finite_element%integ,gpunk)
       call interpolation(finite_element%unkno,1,prev_iter,finite_element%integ,gpunk)
    end if

    ! Previous step unkno interpolation
    call create_scalar(this%physics,1,finite_element%integ,prev_unkno)
    if (this%time_integ%dtinv > 0.0_rp) then
       call interpolation(finite_element%unkno,1,prev_step,finite_element%integ,prev_unkno)
    end if
    
    do igaus = 1,ngaus
       factor = finite_element%integ(1)%p%femap%detjm(igaus) *                                      &
            &   finite_element%integ(1)%p%quad%weight(igaus)

       ! Compute analytical reaction
       if(this%physics%kfl_react>0) then
          call cdr_analytical_reaction(this%physics%ndime,                                          &
               &                       finite_element%integ(1)%p%femap%clocs(:,igaus),              &
               &                       gpunk%a(igaus),this%physics%kfl_react,this%physics%reaction)
       end if

       ! Compute analytical convection
       if (this%physics%kfl_conv>0) then
          call cdr_analytical_convection(this%physics%ndime,                                        &
               &                         finite_element%integ(1)%p%femap%clocs(:,igaus),            &
               &                         this%physics%kfl_conv,this%time_integ%ctime,               &
               &                         this%physics%convection)
       end if

       if (this%discret%kfl_stab == 1) then
          call cdr_cg_tau_supg_compute(this%physics,finite_element,tau,igaus)
       else
          tau = 0.0_rp
       end if
       
       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                ! nu (grad u, grad v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +         &
                     & factor * this%physics%diffusion *                                            &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus) *                &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
                ! (b·grad u, v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +         &
                     & factor * this%physics%convection(idime) *                                    &
                     & finite_element%integ(1)%p%uint_phy%shape(inode,igaus) *                      &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
             ! (react + 1/dt) ( u, v)
             finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +            &
                  & factor * (this%physics%reaction + this%time_integ%dtinv) *                      &
                  & finite_element%integ(1)%p%uint_phy%shape(inode,igaus) *                         &
                  & finite_element%integ(1)%p%uint_phy%shape(jnode,igaus)
             ! SUPG tabilization term
             if (this%discret%kfl_stab == 1) then
                call cdr_cg_supg_matrix_term_compute(finite_element%p_mat%a,inode,jnode,            &
                     & this%physics%convection,this%physics%diffusion,this%physics%reaction,        &
                     & this%time_integ%dtinv,finite_element%integ(1)%p%uint_phy%shape(:,igaus),     &
                     & finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),                         &
                     & finite_element%integ(1)%p%uint_phy%hessi,igaus,tau,factor)
             end if
          end do
          ! Force term (f,v)
          finite_element%p_vec%a(inode) = finite_element%p_vec%a(inode) +                           &
                  & factor * (force%a(igaus) + prev_unkno%a(igaus)*this%time_integ%dtinv) *         &
                  & finite_element%integ(1)%p%uint_phy%shape(inode,igaus)
          ! RHS SUPG tabilization term
          if (this%discret%kfl_stab == 1) then
             call cdr_cg_supg_rhs_term_compute(finite_element%p_vec%a,inode,prev_unkno%a(igaus),    &
                  & force%a(igaus),this%physics%convection,this%time_integ%dtinv,                   &
                     & finite_element%integ(1)%p%uint_phy%deriv(:,:,igaus),tau,factor)
             end if
       end do
    end do

    call memfree(force%a,__FILE__,__LINE__)
    if(this%physics%kfl_react>0) call memfree(gpunk%a,__FILE__,__LINE__)
    call memfree(prev_unkno%a,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine cdr_transient

  !==================================================================================================
  subroutine cdr_cg_supg_matrix_term_compute(matrix,inode,jnode,convection,diffusion,reaction,      &
       &                                     dtinv,shape,deriv,hessi,igaus,tau,factor)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine fills the matrix with the terms corresponding to SUPG stabilization        !
    !----------------------------------------------------------------------------------------------!
    implicit none
    real(rp)             , intent(inout) :: matrix(:,:)
    integer(ip)          , intent(in)    :: inode,jnode,igaus
    real(rp)             , intent(in)    :: convection(:),diffusion,reaction,dtinv
    real(rp)             , intent(in)    :: shape(:), deriv(:,:)
    real(rp), allocatable, intent(in)    :: hessi(:,:,:)
    real(rp)             , intent(in)    :: tau,factor

    real(rp) :: upwinded_gradient,CDR_operator,laplacian

    ! Compute the upwinded gradient
    upwinded_gradient = dot_product(convection,deriv(:,inode))

    ! Compute the laplacian
    if (allocated(hessi)) then
       laplacian = sum(hessi(:,jnode,igaus))
    else
       laplacian = 0.0_rp
    end if

    ! Compute the operator
    CDR_operator = dot_product(convection,deriv(:,jnode)) - diffusion*laplacian +                   &
         &         (reaction+dtinv)*shape(jnode) 

    matrix(inode,jnode) = matrix(inode,jnode) + factor * tau * CDR_operator * upwinded_gradient
    
  end subroutine cdr_cg_supg_matrix_term_compute

  !==================================================================================================
  subroutine cdr_cg_supg_rhs_term_compute(vector,inode,prev_unkno,force,convection,dtinv,deriv,  &
       &                                     tau,factor)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine fills the matrix with the terms corresponding to SUPG stabilization        !
    !----------------------------------------------------------------------------------------------!
    implicit none
    real(rp)   , intent(inout) :: vector(:)
    integer(ip), intent(in)    :: inode
    real(rp)   , intent(in)    :: prev_unkno,force
    real(rp)   , intent(in)    :: convection(:),dtinv
    real(rp)   , intent(in)    :: deriv(:,:)
    real(rp)   , intent(in)    :: tau,factor

    real(rp) :: upwinded_gradient

    upwinded_gradient = dot_product(convection,deriv(:,inode))

    vector(inode) = vector(inode) + factor * tau * (force + dtinv*prev_unkno) * upwinded_gradient
    
  end subroutine cdr_cg_supg_rhs_term_compute

 !=================================================================================================
  subroutine cdr_cg_tau_supg_compute(physics,FE,tau,igaus)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine computes the tau parameter for SUPG.                                       !
    !----------------------------------------------------------------------------------------------!
    implicit none
    type(cdr_problem_t)   , intent(in)    :: physics
    type(finite_element_t), intent(in)    :: FE
    real(rp)              , intent(inout) :: tau
    integer(ip)           , intent(in)    :: igaus

    integer(ip) :: ndime
    real(rp)    :: convection_norm,inverse_tau, h_length
    real(rp)    :: convection(physics%ndime), diffusion, reaction

    ndime    = physics%ndime
    h_length = FE%integ(1)%p%femap%hleng(1,igaus)

    convection = physics%convection
    diffusion  = physics%diffusion
    reaction   = physics%reaction

    ! Compute the norm of the convection
    convection_norm = sqrt(dot_product(convection,convection))

    inverse_tau = 4.0_rp*diffusion/h_length**2.0_rp + 2.0_rp*convection_norm/h_length + reaction

    if (inverse_tau > 1.0e-8) then
       tau = 1.0_rp/inverse_tau      
    else
       tau = 0.0_rp
    end if
  
  end subroutine cdr_cg_tau_supg_compute

end module cdr_stabilized_continuous_Galerkin_names
