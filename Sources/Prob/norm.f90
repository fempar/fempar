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
module norm_names
  use types_names
  use memor_names
  use problem_names
  use analytical_function_names
  use finite_element_names
  use element_fields_names
  use element_tools_names
  implicit none
# include "debug.i90"
  private
  
  ! Error norm
  type, extends(discrete_integration_t) :: error_norm_t
     class(discrete_problem_t), pointer :: discret
     class(physical_problem_t), pointer :: physics
     integer(ip)                        :: unknown_id
     real(rp)                           :: ctime,alpha,beta,gamma
   contains
     procedure :: create  => error_norm_create
     procedure :: compute => error_norm_compute
     procedure :: free    => error_norm_free
  end type error_norm_t
  
  ! Types
  public :: error_norm_t

contains

  !==================================================================================================
  subroutine error_norm_create(this,physics,discret)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine creates the pointers needed for the discrete integration type              !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(error_norm_t)              , intent(inout) :: this
    class(physical_problem_t), target, intent(in)    :: physics
    class(discrete_problem_t), target, intent(in)    :: discret

    this%physics => physics
    this%discret => discret
    this%unknown_id = 0
    this%ctime = 0.0_rp
    this%alpha = 0.0_rp
    this%beta  = 0.0_rp
    this%gamma = 0.0_rp

    this%domain_dimension = 3
    
  end subroutine error_norm_create

  !=================================================================================================
  subroutine error_norm_free(this)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine deallocates the pointers needed for the discrete integration type          !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(error_norm_t), intent(inout) :: this

    this%physics => null()
    this%discret => null()
    this%unknown_id = 0
    this%ctime = 0.0_rp
    this%alpha = 0.0_rp
    this%beta  = 0.0_rp
    this%gamma = 0.0_rp

  end subroutine error_norm_free

  !=================================================================================================
  subroutine error_norm_compute(this,finite_element)
    !----------------------------------------------------------------------------------------------!
    !   This subroutine performs the elemental error norm integration.                             !
    !----------------------------------------------------------------------------------------------!
    implicit none
    class(error_norm_t)   , intent(inout) :: this
    type(finite_element_t), intent(inout) :: finite_element
    ! Locals
    integer(ip)    :: initial_unknown,final_unknown,iunkn,ivar,jvar,igaus
    real(rp)       :: exact_values(11),dvolu
    type(vector_t) :: gpvector
    type(scalar_t) :: gpscalar

    if(this%unknown_id == 0) then
       initial_unknown = 1
       final_unknown = this%physics%nunks
    else
       initial_unknown = this%unknown_id
       final_unknown   = this%unknown_id
    end if

    ! Initialize scalar
    finite_element%scalar = 0.0_rp

    do iunkn = initial_unknown,final_unknown

       ivar = 1
       if(iunkn>1) ivar = sum(this%physics%vars_of_unk(1:iunkn))
       
       ! Vector unknown
       if(this%physics%vars_of_unk(iunkn)>1) then

          ! Interpolation operations
          call create_vector(this%physics,iunkn,finite_element%integ,gpvector)
          call interpolation(finite_element%unkno, ivar, 1, finite_element%integ, gpvector)

          ! Loop over Gauss points
          do igaus=1,finite_element%integ(1)%p%quad%ngaus
             dvolu=finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

             do ivar=1,this%physics%vars_of_unk(iunkn)
                jvar = ivar + sum(this%physics%vars_of_unk(1:iunkn-1))

                ! Evaluate analytical unknown
                call evaluate_analytical(finite_element%p_analytical_code%a(jvar,1),                      &
                     &                   finite_element%p_analytical_code%a(jvar,2),this%physics%ndime,   &
                     &                   finite_element%integ(1)%p%femap%clocs(:,igaus),                  &
                     &                   this%ctime,exact_values)

                ! Error computation
                finite_element%scalar = finite_element%scalar + &
                     &                  ((gpvector%a(ivar,igaus)-exact_values(1))**2)*dvolu

             end do
          end do
          
          ! Deallocate gpvector
          call memfree(gpvector%a,__FILE__,__LINE__)

       end if

       ! Scalar unknown
       if(this%physics%vars_of_unk(iunkn)==1) then

          ! Interpolation operations
          call create_scalar(this%physics,iunkn,finite_element%integ,gpscalar)
          call interpolation(finite_element%unkno, ivar, 1, finite_element%integ, gpscalar)

          ! Loop over Gauss points
          do igaus=1,finite_element%integ(1)%p%quad%ngaus
             dvolu=finite_element%integ(1)%p%quad%weight(igaus)*finite_element%integ(1)%p%femap%detjm(igaus)

             jvar = 1 + sum(this%physics%vars_of_unk(1:iunkn-1))

             ! Evaluate analytical unknown
             call evaluate_analytical(finite_element%p_analytical_code%a(jvar,1),                      &
                  &                   finite_element%p_analytical_code%a(jvar,2),this%physics%ndime,   &
                  &                   finite_element%integ(1)%p%femap%clocs(:,igaus),                  &
                  &                   this%ctime,exact_values)

             ! Error computation
             finite_element%scalar = finite_element%scalar + &
                  &                  ((gpscalar%a(igaus)-exact_values(1))**2)*dvolu
             
          end do
          
          ! Deallocate gpvector
          call memfree(gpscalar%a,__FILE__,__LINE__)
       
       end if

    end do

  end subroutine error_norm_compute
     
end module norm_names
