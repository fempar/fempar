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
module nonlinear_operator_names
  use types_names
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

  type, abstract :: nonlinear_operator_t
     integer(ip)                                :: max_iter = 1          ! Maximum nonlinear iterations
     real(rp)                                   :: nltol    = 1.0e-12_rp ! Nonlinear tolerance
     class(base_operator_t)           , pointer :: A => NULL()           ! Matrix operator pointer
     class(base_operator_t)           , pointer :: M => NULL()           ! Preconditioner operator pointer
     class(base_operand_t)            , pointer :: b => NULL()           ! RHS operand pointer
     class(base_operand_t)            , pointer :: x => NULL()           ! Solution operand pointer
     class(base_integrable_operator_t), pointer :: A_int => NULL()       ! (par/block) matrix pointer
     class(base_operand_t)            , pointer :: b_int => NULL()       ! (par/block) vector pointer
     class(base_operand_t)            , pointer :: x_sol => NULL()       ! (par/block) vector pointer
   contains
     procedure (set_preconditioner_flags_interface), deferred :: set_preconditioner_flags
     procedure                                                :: do_nonlinear_iteration
  end type nonlinear_operator_t

  ! Abstract interfaces
  abstract interface
     subroutine set_preconditioner_flags_interface(nlop,istep,istge,iiter)
       import :: nonlinear_operator_t, ip
       implicit none
       class(nonlinear_operator_t), intent(inout) :: nlop
       integer(ip)                , intent(in)    :: istep,istge,iiter
     end subroutine set_preconditioner_flags_interface
  end interface

  ! Types
  public :: nonlinear_operator_t

contains

  !==================================================================================================
  subroutine do_nonlinear_iteration( nlop, istep, istge, sctrl, env, approx, fe_space )
    implicit none
    class(nonlinear_operator_t), target , intent(inout) :: nlop
    integer(ip)                         , intent(in)    :: istep,istge
    type(solver_control_t)              , intent(inout) :: sctrl
    class(abstract_environment_t)       , intent(in)    :: env
    type(discrete_integration_pointer_t), intent(inout) :: approx(:)
    type(fe_space_t)                    , intent(inout) :: fe_space
    ! Locals
    ! Locals
    integer(ip)           :: iiter
    real(rp)              :: resnorm,ininorm
    class(base_operand_t), allocatable :: y

    ! Checks
    check(associated(nlop%A))
    check(associated(nlop%M))
    check(associated(nlop%b))
    check(associated(nlop%x))
    check(associated(nlop%A_int))
    check(associated(nlop%b_int))
    check(associated(nlop%x_sol))

    ! Allocate y
    allocate(y, mold=nlop%b); call y%default_initialization()
        
    iiter = 0
    do while( iiter < nlop%max_iter )

       ! Update counter
       iiter = iiter+1

       ! Initialize Matrix and vector
       call nlop%A_int%init()
       call nlop%b_int%init(0.0_rp)

       ! Integrate system
       call volume_integral(approx,fe_space,nlop%A_int,nlop%b_int)

       ! Check convergence
       if(iiter==1) ininorm = nlop%b%nrm2()   
       y = nlop%b - nlop%A*nlop%x
       resnorm = y%nrm2()
       if( resnorm < nlop%nltol*ininorm) then
          write(*,*) 'Nonlinear iterations: ', iiter
          write(*,*) 'Nonlinear error norm: ', resnorm
          exit
       end if

       ! Set preconditioner flags
       call nlop%set_preconditioner_flags(istep,istge,iiter)

       ! Compute Numeric preconditioner
       call nlop%M%fill_values()

       ! Solve system
       call abstract_solve(nlop%A,nlop%M,nlop%b,nlop%x,sctrl,env)
       call solver_control_log_conv_his(sctrl)
       call solver_control_free_conv_his(sctrl)

       ! Free Numeric preconditioner
       call nlop%M%free_values()
       
       ! Store solution to unkno
       call update_solution(nlop%x_sol,fe_space)
       
       ! Store nonlinear iteration ( k+1 --> k )
       call update_nonlinear(fe_space)
       
    end do

    ! Deallocate
    call y%free()

  end subroutine do_nonlinear_iteration

end module nonlinear_operator_names
