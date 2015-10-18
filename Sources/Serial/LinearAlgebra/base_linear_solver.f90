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
module base_linear_solver_names
  use types_names

  ! Abstract modules
  use vector_names
  use operator_names

  implicit none
# include "debug.i90"
  private

  integer(ip), parameter :: start               = 0  ! All parameters set with values and environment set
  integer(ip), parameter :: operators_set       = 1  ! Matrix A and preconditioner M already set
  integer(ip), parameter :: workspace_allocated = 2  ! All workspace required by solve TBP available 
  
  ! List of convergence criteria available for iterative solvers 
  integer(ip), parameter :: res_nrmgiven_rhs_nrmgiven  = 1  ! ||  r(i) ||g <= rtol*||  b    ||g + atol 
  integer(ip), parameter :: res_nrmgiven_res_nrmgiven  = 2  ! ||  r(i) ||g <= rtol*||  r(0) ||g + atol   
  integer(ip), parameter :: delta_rhs                  = 3  ! || dx(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_delta                = 4  ! || dx(i) ||  <= rtol*||dx(1)|| + atol
  integer(ip), parameter :: res_res                    = 5  ! ||  r(i) ||  <= rtol*|| r(0)|| + atol
  integer(ip), parameter :: res_rhs                    = 6  ! ||  r(i) ||  <= rtol*||  b  || + atol
  integer(ip), parameter :: delta_rhs_and_res_res      = 7  ! delta_rhs    AND res_res
  integer(ip), parameter :: delta_rhs_and_res_rhs      = 8  ! delta_rhs    AND res_rhs
  integer(ip), parameter :: delta_delta_and_res_res    = 9  ! delta_delta  AND res_res
  integer(ip), parameter :: delta_delta_and_res_rhs    = 10 ! delta_delta  AND res_rhs 
                                                            ! ||.|| is the 2-norm, dx(i) = x(i) - x(i-1),
                                                            ! r(i) is the residual at the i-th iteration
  
  type, abstract :: base_linear_solver_t
    ! Solver name
    character(len=:), allocatable :: name

    ! State
    integer(ip) :: state = start
  
    ! Matrix and preconditioner
    class(operator_t), pointer :: A 
    class(operator_t), pointer :: M 
  
    ! Parameters
    real(rp)      :: rtol                      ! Relative tolerance
    real(rp)      :: atol                      ! Absolute tolerance
    integer(ip)   :: stopping_criteria         ! Stopping criteria
    integer(ip)   :: output_frequency          ! Message every output_frequency iterations
    integer(ip)   :: max_num_iterations        ! Max. # of iterations
    logical       :: store_convergence_history ! Is the convergence history going to be stored? 
    
    ! Output parameters
    logical               :: converged          ! Converged?
    integer(ip)           :: num_iterations     ! # of iterations to convergence
    real(rp)              :: rhs_convergence_test     
    real(rp)              :: rhs_extra_convergence_test     
    real(rp)              :: error_estimate_convergence_test
    real(rp)              :: error_estimate_extra_convergence_test
    real(rp), allocatable :: error_estimate_history_convergence_test(:)
    real(rp), allocatable :: error_estimate_history_extra_convergence_test(:)
  contains
    ! Concrete TBPs
    procedure          :: set_operators
    procedure          :: did_converge 
    procedure          :: get_num_iterations
    procedure          :: get_error_estimate_convergence_test
    procedure          :: get_error_estimate_extra_convergence_test
    procedure          :: print_convergence_history
    procedure, private :: unset_operators
    procedure, private :: print_convergence_history_header
    procedure, private :: print_convergence_history_footer
    
    ! Deferred TBPs
    procedure (set_parameters_interface), deferred :: set_parameters
    procedure (solve_interface)         , deferred :: solve
    procedure (free_interface)          , deferred :: free 
  end type
  
  abstract interface
    subroutine set_parameters_interface(this) ! Parameter List still missing
     import :: base_linear_solver_t
     implicit none
     class(base_linear_solver_t), intent(in) :: this
    end subroutine set_parameters_interface
    
    subroutine solve_interface(this,x,y)
     import :: base_linear_solver_t, vector_t
     implicit none
     class(base_linear_solver_t), intent(in)    :: this
     class(vector_t)            , intent(in)    :: x
     class(vector_t)            , intent(inout) :: y
    end subroutine solve_interface

    subroutine free_interface(this)
     import :: base_linear_solver_t
     implicit none
     class(base_linear_solver_t), intent(inout) :: this
    end subroutine free_interface
  end interface
  
  ! Data types
  public :: base_linear_solver_t
  
contains

    subroutine set_operators(this,A,M)
     implicit none
     class(base_linear_solver_t)        , intent(inout) :: this
     class(operator_t)          , target, intent(in)    :: A
     class(operator_t)          , target, intent(in)    :: M
     assert(this%state == start .or. this%state == operators_set)
     this%A => A
     this%M => M 
     this%state = operators_set
    end subroutine set_operators
    
    subroutine unset_operators(this)
     implicit none
     class(base_linear_solver_t), intent(inout) :: this
     nullify(this%A)
     nullify(this%M) 
    end subroutine unset_operators
    
    function did_converge(this)
      implicit none
      class(base_linear_solver_t), intent(in) :: this
      logical :: did_converge 
      did_converge = this%converged
    end function did_converge 
    
    function get_num_iterations(this)
      implicit none
      class(base_linear_solver_t), intent(in) :: this
      integer(ip) :: get_num_iterations
      get_num_iterations = this%num_iterations
    end function get_num_iterations
    
    function get_error_estimate_convergence_test(this)
      implicit none
      class(base_linear_solver_t), intent(in) :: this
      real(rp) :: get_error_estimate_convergence_test
      get_error_estimate_convergence_test = this%error_estimate_convergence_test
    end function
    
    function get_error_estimate_extra_convergence_test(this)
      implicit none
      class(base_linear_solver_t), intent(in) :: this
      real(rp) :: get_error_estimate_extra_convergence_test
      get_error_estimate_extra_convergence_test = this%error_estimate_extra_convergence_test
    end function
    
    subroutine print_convergence_history ( this, file_path ) 
      implicit none
      class(base_linear_solver_t), intent(in) :: this
      character(len=*)           , intent(in) :: file_path
      
      
    end subroutine print_convergence_history

    subroutine print_convergence_history_header( this, luout )
      implicit none
      ! Parameters
      class(base_linear_solver_t), intent(in) :: this 
      integer(ip)                , intent(in) :: luout

      ! Local variables
      character(len=*), parameter    :: fmt1='(a18,1x,a4,3(2x,a15))'
      character(len=*), parameter    :: fmt2='(a18,1x,a4,3(2x,a15),3(2x,a15))'
      character(len=:), allocatable  :: outname

      outname = this%name // ':' // '  '
      select case(this%stopping_criteria)
      case ( delta_rhs, delta_delta, res_res, res_rhs, res_nrmgiven_rhs_nrmgiven, &
           & res_nrmgiven_res_nrmgiven )
         write(luout,fmt1) adjustl(outname),'Iteration','Error Estimate','Tolerance'
      case ( delta_rhs_and_res_res, delta_rhs_and_res_rhs,  &
            delta_delta_and_res_res, delta_delta_and_res_rhs )
         write(luout,fmt2) adjustl(outname), 'Iteration', 'Error Estimate', 'Tolerance', &
                                           & 'Error Estimate', 'Tolerance'
      case default
         ! Write an error message and stop ?      
      end select
    end subroutine print_convergence_history_header 

    subroutine print_convergence_history_footer ( this, luout )
      implicit none
      ! Parameters
      class(base_linear_solver_t), intent(in) :: this
      integer(ip)                , intent(in) :: luout 

      character(len=*), parameter  :: fmt11='(a,2x,es16.9,1x,a,1x,i4,1x,a)'
      character(len=*), parameter  :: fmt12='(a,3(2x,es16.9))'

      character(len=*), parameter  :: fmt21='(a,2x,es16.9,1x,es16.9,1x,a,1x,i4,1x,a)'
      character(len=*), parameter  :: fmt22='(a,3(2x,es16.9),3(2x,es16.9))'

      select case( this%stopping_criteria )
      case ( delta_rhs,delta_delta,res_res,res_rhs,&
           & res_nrmgiven_rhs_nrmgiven, res_nrmgiven_res_nrmgiven)
         if ( this%converged ) then
            write(luout,fmt11) this%name //' converged to ', &
                 & this%rhs_convergence_test,' in ',this%num_iterations,' iterations. '
            write(luout,fmt12) 'Last iteration error estimate: ', this%error_estimate_convergence_test
         else
            write(luout,fmt11) this%name //' failed to converge to ', &
                 & this%rhs_convergence_test,' in ',this%num_iterations,' iterations. '
            write(luout,fmt12) 'Last iteration error estimate: ', this%error_estimate_convergence_test
         end if
      case ( delta_rhs_and_res_res  , delta_rhs_and_res_rhs, &
           & delta_delta_and_res_res, delta_delta_and_res_rhs )
         if ( this%converged ) then
            write(luout,fmt21) this%name //' converged to ', &
                & this%rhs_convergence_test, this%rhs_extra_convergence_test, ' in ', this%num_iterations ,' iterations. '
            write(luout,fmt22) 'Last iteration error estimates: ', this%error_estimate_convergence_test, this%error_estimate_extra_convergence_test
         else
            write(luout,fmt21) this%name //' failed to converge to ', &
                & this%rhs_convergence_test, this%rhs_extra_convergence_test, ' in ', this%num_iterations ,' iterations. '
            write(luout,fmt22) 'Last iteration error estimates: ', this%error_estimate_convergence_test, this%error_estimate_extra_convergence_test
         end if
      case default
         ! Write an error message and stop ?      
      end select
    end subroutine print_convergence_history_footer 

end module base_linear_solver_names
