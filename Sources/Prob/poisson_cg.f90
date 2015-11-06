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
module poisson_continuous_galerkin_names
  use types_names
  use memor_names
  use problem_names
  use poisson_names
  use element_fields_names
  use element_tools_names
  use finite_element_names
  implicit none
# include "debug.i90"
  private 

  type, extends(discrete_problem_t) :: poisson_discrete_t
     contains
      procedure :: create => poisson_discrete_create
   end type poisson_discrete_t

   type, extends(discrete_integration_t) :: poisson_integration_t
      type(poisson_discrete_t), pointer :: discret
      type(poisson_problem_t) , pointer :: physics
    contains
      procedure :: create  => poisson_integration_create
      procedure :: compute => poisson_integration_compute
      procedure :: free    => poisson_integration_free
   end type poisson_integration_t

 public :: poisson_discrete_t, poisson_integration_t

contains

  !=================================================================================================
  subroutine poisson_discrete_create(discret,physics,l2g)
    !----------------------------------------------------------------------------------------------!
    !This subroutine contains definitions of the poisson problem approximed by continuous Galerkin.!
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(poisson_discrete_t)    , intent(out) :: discret
    class(physical_problem_t), intent(in)  :: physics
    integer(ip), optional  , intent(in)  :: l2g(:)
    ! Locals
    integer(ip) :: i

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
    
  end subroutine poisson_discrete_create

  !=================================================================================================
  subroutine poisson_integration_create( this, physics, discret )
    implicit none
    class(poisson_integration_t)     , intent(inout) :: this
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret

    select type (physics)
    type is(poisson_problem_t)
       this%physics => physics
    class default
       check(.false.)
    end select 
    select type (discret)
    type is(poisson_discrete_t)
       this%discret => discret
    class default
       check(.false.)
    end select

    ! Allocate working variables
    call memalloc(1,this%working_vars,__FILE__,__LINE__)
    this%working_vars(1) = 1

    this%domain_dimension = 3 

  end subroutine poisson_integration_create

  !=================================================================================================
  subroutine poisson_integration_free(this)
    implicit none
    class(poisson_integration_t)  , intent(inout) :: this

    this%physics => null()
    this%discret => null()

    ! Deallocate working variables
    call memfree(this%working_vars,__FILE__,__LINE__)

  end subroutine poisson_integration_free

  !=================================================================================================
  subroutine poisson_integration_compute(this,finite_element)
    implicit none
    class(poisson_integration_t), intent(inout) :: this
    type(finite_element_t)    , intent(inout) :: finite_element
    ! Locals
    type(scalar_t) :: force
    integer(ip)    :: idime,igaus,inode,jnode,ngaus,nnode,ndime
    real(rp)       :: factor

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
       call poisson_analytical_force(this%physics,finite_element,0.0_rp,force)
    end if

    do igaus = 1,ngaus
       factor = finite_element%integ(1)%p%femap%detjm(igaus) *                                      &
            &   finite_element%integ(1)%p%quad%weight(igaus)

       do inode = 1, nnode
          do jnode = 1, nnode
             do idime = 1,ndime
                ! nu (grad u, grad v)
                finite_element%p_mat%a(inode,jnode) = finite_element%p_mat%a(inode,jnode) +         &
                     & factor * this%physics%diffu *                                                &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,inode,igaus) *                &
                     & finite_element%integ(1)%p%uint_phy%deriv(idime,jnode,igaus)
             end do
          end do
          ! Force term (f,v)
          finite_element%p_vec%a(inode) = finite_element%p_vec%a(inode) +                           &
               & factor * force%a(igaus) *                                                          &
               & finite_element%integ(1)%p%uint_phy%shape(inode,igaus)
       end do
    end do

    call memfree(force%a,__FILE__,__LINE__)

    ! Apply boundary conditions
    call impose_strong_dirichlet_data(finite_element) 

  end subroutine poisson_integration_compute

end module poisson_continuous_galerkin_names
