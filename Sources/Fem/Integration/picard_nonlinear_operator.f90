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
module picard_nonlinear_operator_names
  use types_names
  use nonlinear_operator_names
  use base_operator_names
  use base_operand_names
  use base_integrable_operator_names
  use abstract_solver_names
  use abstract_environment_names
  use problem_names
  use fe_space_names
  use base_operand_names
  use integration_names
  use update_names
  implicit none
# include "debug.i90"
  private

  type, abstract, extends(nonlinear_operator_t) :: picard_nonlinear_operator_t
     integer(ip)                                :: max_iter = 1          ! Maximum nonlinear iterations
     real(rp)                                   :: nltol    = 1.0e-12_rp ! Nonlinear tolerance
     class(base_operator_t)           , pointer :: A => NULL()           ! Matrix operator pointer
     class(base_operator_t)           , pointer :: M => NULL()           ! Preconditioner operator pointer
     class(base_integrable_operator_t), pointer :: A_int_c => NULL()     ! (par/block) matrix pointer (Constant) 
     class(base_integrable_operator_t), pointer :: A_int_t => NULL()     ! (par/block) matrix pointer (Transient)
     class(base_integrable_operator_t), pointer :: A_int_n => NULL()     ! (par/block) matrix pointer (Nonlinear)
     class(base_operand_t)            , pointer :: b_int_c => NULL()     ! (par/block) vector pointer (Constant) 
     class(base_operand_t)            , pointer :: b_int_t => NULL()     ! (par/block) vector pointer (Transient)
     class(base_operand_t)            , pointer :: b_int_n => NULL()     ! (par/block) vector pointer (Nonlinear)
     class(base_operand_t)            , pointer :: x_sol => NULL()       ! (par/block) vector pointer
   contains
     procedure :: apply          => do_picard_nonlinear_iteration
     procedure :: create         => create_picard_nonlinear_iteration
     procedure :: unassign       => unassign_picard_nonlinear_iteration
     procedure :: fill_constant  => compute_picard_constant_operator
     procedure :: fill_transient => compute_picard_transient_operator
  end type picard_nonlinear_operator_t

  ! Types
  public :: picard_nonlinear_operator_t

contains

  !==================================================================================================
  subroutine create_picard_nonlinear_iteration(nlop,A,M,b,x,x_sol,A_int_c,A_int_t,A_int_n,b_int_c, &
       &                                       b_int_t,b_int_n)
    implicit none
    class(picard_nonlinear_operator_t)                 , intent(inout) :: nlop
    class(base_operator_t)                     , target, intent(in)    :: A
    class(base_operator_t)                     , target, intent(in)    :: M
    class(base_operand_t)                      , target, intent(in)    :: b
    class(base_operand_t)                      , target, intent(in)    :: x
    class(base_operand_t)                      , target, intent(in)    :: x_sol
    class(base_integrable_operator_t), optional, target, intent(in)    :: A_int_c
    class(base_integrable_operator_t), optional, target, intent(in)    :: A_int_t
    class(base_integrable_operator_t), optional, target, intent(in)    :: A_int_n
    class(base_operand_t)            , optional, target, intent(in)    :: b_int_c
    class(base_operand_t)            , optional, target, intent(in)    :: b_int_t
    class(base_operand_t)            , optional, target, intent(in)    :: b_int_n

    ! Assign operators
    nlop%A       => A      
    nlop%M       => M      
    nlop%b       => b      
    nlop%x       => x      
    nlop%x_sol   => x_sol
    if(present(A_int_c)) nlop%A_int_c => A_int_c
    if(present(A_int_t)) nlop%A_int_t => A_int_t
    if(present(A_int_n)) nlop%A_int_n => A_int_n
    if(present(b_int_c)) nlop%b_int_c => b_int_c
    if(present(b_int_t)) nlop%b_int_t => b_int_t
    if(present(b_int_n)) nlop%b_int_n => b_int_n
    
  end subroutine create_picard_nonlinear_iteration
  
  !==================================================================================================
  subroutine unassign_picard_nonlinear_iteration(nlop)
    implicit none
    class(picard_nonlinear_operator_t) , intent(inout) :: nlop

    ! Unassign operators
    nlop%A       => null()  
    nlop%M       => null()  
    nlop%b       => null()  
    nlop%x       => null()  
    nlop%A_int_c => null()
    nlop%A_int_t => null()
    nlop%A_int_n => null()
    nlop%b_int_c => null()
    nlop%b_int_t => null()
    nlop%b_int_n => null()
    nlop%x_sol   => null()

  end subroutine unassign_picard_nonlinear_iteration

  !==================================================================================================
  subroutine compute_picard_constant_operator(nlop,approx,fe_space)
    implicit none
    class(picard_nonlinear_operator_t)  , intent(inout) :: nlop
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    ! Locals
    integer(ip) :: napprox,iapprox

    ! Set approx flags
    napprox = size(approx,1)
    do iapprox = 1,napprox
       approx(iapprox)%p%integration_stage = update_constant
    end do

    ! Initialize
    if(associated(nlop%A_int_c)) call nlop%A_int_c%init()
    if(associated(nlop%b_int_c)) call nlop%b_int_c%init(0.0_rp)

    ! Integrate
    if(associated(nlop%A_int_c).and.associated(nlop%b_int_c)) then
       call volume_integral(approx,fe_space,nlop%A_int_c,nlop%b_int_c)
    elseif(associated(nlop%A_int_c).and.(.not.associated(nlop%b_int_c))) then
       call volume_integral(approx,fe_space,nlop%A_int_c)
    elseif((.not.associated(nlop%A_int_c)).and.associated(nlop%b_int_c)) then
       call volume_integral(approx,fe_space,nlop%b_int_c)
    end if

    ! Fill Preconditioner
    call nlop%M%fill_values(update_constant)

  end subroutine compute_picard_constant_operator

  !==================================================================================================
  subroutine compute_picard_transient_operator(nlop,approx,fe_space)
    implicit none
    class(picard_nonlinear_operator_t)  , intent(inout) :: nlop
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    ! Locals
    integer(ip) :: napprox,iapprox

    ! Set approx flags
    napprox = size(approx,1)
    do iapprox = 1,napprox
       approx(iapprox)%p%integration_stage = update_transient
    end do

    ! Initialize
    if(associated(nlop%A_int_t)) call nlop%A_int_t%init()
    if(associated(nlop%b_int_t)) call nlop%b_int_t%init(0.0_rp)

    ! Integrate
    if(associated(nlop%A_int_t).and.associated(nlop%b_int_t)) then
       call volume_integral(approx,fe_space,nlop%A_int_t,nlop%b_int_t)
    elseif(associated(nlop%A_int_t).and.(.not.associated(nlop%b_int_t))) then
       call volume_integral(approx,fe_space,nlop%A_int_t)
    elseif((.not.associated(nlop%A_int_t)).and.associated(nlop%b_int_t)) then
       call volume_integral(approx,fe_space,nlop%b_int_t)
    end if

    ! Fill Preconditioner
    call nlop%M%fill_values(update_transient)

  end subroutine compute_picard_transient_operator

  !==================================================================================================
  subroutine do_picard_nonlinear_iteration( this, sctrl, env, approx, fe_space )
    implicit none
    class(picard_nonlinear_operator_t), target , intent(inout) :: this
    type(solver_control_t)              , intent(inout) :: sctrl
    class(abstract_environment_t)       , intent(in)    :: env
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    ! Locals
    integer(ip)           :: iiter
    integer(ip)           :: napprox,iapprox
    real(rp)              :: resnorm,ininorm
    class(base_operand_t), allocatable :: y,y2

    ! Checks
    check(associated(this%A))
    check(associated(this%M))
    check(associated(this%b))
    check(associated(this%x))
    check(associated(this%x_sol))

    ! Set approx flags
    napprox = size(approx,1)
    do iapprox = 1,napprox
       approx(iapprox)%p%integration_stage = update_transient
    end do

    ! Allocate y
    allocate(y, mold=this%b); call y%default_initialization()
    allocate(y2, mold=this%b); call y2%default_initialization()
        
    iiter = 0
    do while( iiter < this%max_iter )

       ! Update counter
       iiter = iiter+1

!!$       ! Initialize Matrix and vector
!!$       if(associated(this%A_int_n)) call this%A_int_n%init()
!!$       if(associated(this%b_int_n)) call this%b_int_n%init(0.0_rp)
!!$
!!$       ! Integrate system
!!$       if(associated(this%A_int_n).and.associated(this%b_int_n)) then
!!$          call volume_integral(approx,fe_space,this%A_int_n,this%b_int_n)
!!$       elseif(associated(this%A_int_n).and.(.not.associated(this%b_int_n))) then
!!$          call volume_integral(approx,fe_space,this%A_int_n)
!!$       elseif((.not.associated(this%A_int_n)).and.associated(this%b_int_n)) then
!!$          call volume_integral(approx,fe_space,this%b_int_n)
!!$       end if

       write(*,*) iiter
       ! Check convergence
       !if(iiter==1) ininorm = this%b%nrm2()   
       y = this%A*this%x
       write(*,*) 'xxxx'
       y2 = this%A*this%x
       !y = this%b - this%A*this%x
!!$       resnorm = y%nrm2()
!!$       if( resnorm < this%nltol*ininorm) then
!!$          write(*,*) 'Nonlinear iterations: ', iiter
!!$          write(*,*) 'Nonlinear error norm: ', resnorm
!!$          exit
!!$       end if
!!$       if(iiter==this%max_iter) then
!!$          write(*,*) 'Maximum number of nonlinear iterations reached'
!!$          write(*,*) 'Nonlinear iterations: ',iiter
!!$          write(*,*) 'Nonlinear error norm: ', resnorm
!!$       end if
!!$
!!$       ! Compute Numeric preconditioner
!!$       call this%M%fill_values(update_nonlinear)
!!$
!!$       ! Solve system
!!$       call abstract_solve(this%A,this%M,this%b,this%x,sctrl,env)
!!$       call solver_control_log_conv_his(sctrl)
!!$       call solver_control_free_conv_his(sctrl)
!!$
!!$       ! Free Numeric preconditioner
!!$       call this%M%free_values()
!!$       
!!$       ! Store solution to unkno
!!$       call update_solution(this%x_sol,fe_space)
!!$       
!!$       ! Store nonlinear iteration ( k+1 --> k )
!!$       call update_nonlinear_solution(fe_space)
    call y%free()
       
    end do

    ! Deallocate
    call y%free()

  end subroutine do_picard_nonlinear_iteration

end module picard_nonlinear_operator_names
