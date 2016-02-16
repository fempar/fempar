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
module rgmres_names
  use types_names
  use stdio_names
  use memor_names
  
  ! Abstract modules
  use vector_names
  use vector_space_names
  use operator_names
  use environment_names
  use base_linear_solver_names  
  use multivector_names

  implicit none
# include "debug.i90"
  private
  
  character(len=*), parameter :: rgmres_name = 'RGMRES'
  character(len=*), parameter :: ls_dkrymax = 'linear_solver_dkrymax'
  
  integer (ip)    , parameter :: default_rgmres_stopping_criteria = res_nrmgiven_res_nrmgiven
  integer (ip)    , parameter :: default_rgmres_dkrymax           = 100
  type, extends(base_linear_solver_t) :: rgmres_t
     ! Parameters
     integer(ip)                    :: dkrymax

     ! Working space data members
     type(multivector_t)            :: bkry
     class(vector_t), allocatable   :: r, z  
     real(rp)       , allocatable   :: hh(:,:), g(:), cs(:,:)
   contains
     procedure          :: allocate_workspace            => rgmres_allocate_workspace
     procedure          :: free_workspace                => rgmres_free_workspace
     procedure          :: set_parameters_from_pl        => rgmres_set_parameters_from_pl
     procedure          :: solve_body                    => rgmres_solve_body
     procedure          :: supports_stopping_criteria    => rgmres_supports_stopping_criteria
     procedure          :: get_default_stopping_criteria => rgmres_get_default_stopping_criteria
  end type rgmres_t
  
  ! Data types
  public :: rgmres_t, create_rgmres
  
contains
  subroutine rgmres_allocate_workspace(this)
    implicit none
    class(rgmres_t), intent(inout) :: this
    type(vector_space_t), pointer :: range
    type(dynamic_state_operator_t), pointer :: A, M
    class(environment_t), pointer :: environment
    A => this%get_A()
    range  => A%get_range_vector_space()
    call range%create_vector(this%r)
    M => this%get_M()
    range  => M%get_range_vector_space()
    call range%create_vector(this%z)
    environment => this%get_environment()
    call this%bkry%create(environment,this%dkrymax+1,mold=this%z)
    call memalloc(this%dkrymax+1,this%dkrymax+1,this%hh,__FILE__,__LINE__)
    call memalloc(this%dkrymax+1,this%g,__FILE__,__LINE__)
    call memalloc(2,this%dkrymax+1,this%cs,__FILE__,__LINE__)
  end subroutine rgmres_allocate_workspace
  
  subroutine rgmres_free_workspace(this)
    implicit none
    class(rgmres_t), intent(inout) :: this
    call this%r%free()
    call this%z%free()
    call this%bkry%free()
    call memfree(this%hh,__FILE__,__LINE__)
    call memfree(this%g,__FILE__,__LINE__)
    call memfree(this%cs,__FILE__,__LINE__)
  end subroutine rgmres_free_workspace

  subroutine rgmres_set_parameters_from_pl(this) 
   implicit none
   class(rgmres_t), intent(inout) :: this
  end subroutine rgmres_set_parameters_from_pl
  
  subroutine rgmres_solve_body(this,x)
    implicit none
    class(rgmres_t), intent(inout) :: this
    class(vector_t)    , intent(inout) :: x 

    ! Locals
    class(vector_t)     , pointer :: initial_solution
    initial_solution          => this%get_initial_solution()

    call x%copy(initial_solution)

!!$    real(rp) :: r_nrm_M      ! |r|_inv(M) 
!!$    real(rp) :: b_nrm_M      ! |b|_inv(M)
!!$    real(rp) :: r_z, Ap_p, alpha, beta
!!$ 
!!$    ! Local variables to store a copy/reference of the corresponding member variables of base class
!!$    class(environment_t), pointer :: environment
!!$    class(operator_t)   , pointer :: A, M 
!!$    class(vector_t)     , pointer :: initial_solution, b
!!$    integer(ip)                   :: stopping_criteria, max_num_iterations, output_frequency, luout
!!$    real(rp)                      :: atol, rtol
!!$    logical                       :: track_convergence_history
!!$
!!$    ! Pointers to freely modify/read private member variables of base class
!!$    integer(ip), pointer :: num_iterations
!!$    logical    , pointer :: did_converge
!!$    real(rp)   , pointer :: rhs_convergence_test, error_estimate_convergence_test
!!$    real(rp)   , pointer :: error_estimate_history_convergence_test(:)
!!$
!!$    environment               => this%get_environment()
!!$    A                         => this%get_A()
!!$    b                         => this%get_rhs()
!!$    M                         => this%get_M()
!!$    initial_solution          => this%get_initial_solution()
!!$    luout                     =  this%get_luout()
!!$    stopping_criteria         =  this%get_stopping_criteria()
!!$    max_num_iterations        =  this%get_max_num_iterations()
!!$    atol                      =  this%get_atol()
!!$    rtol                      =  this%get_rtol()
!!$    output_frequency          =  this%get_output_frequency()
!!$    track_convergence_history =  this%get_track_convergence_history()
!!$     
!!$    num_iterations                          => this%get_pointer_num_iterations()
!!$    did_converge                            => this%get_pointer_did_converge()
!!$    rhs_convergence_test                    => this%get_pointer_rhs_convergence_test()
!!$    error_estimate_convergence_test         => this%get_pointer_error_estimate_convergence_test()
!!$    error_estimate_history_convergence_test => this%get_pointer_error_estimate_history_convergence_test()
!!$
!!$    call x%copy(initial_solution)
!!$
!!$    ! Evaluate |b|_inv(M) if required
!!$    if ( stopping_criteria == res_nrmgiven_rhs_nrmgiven ) then
!!$       ! r = inv(M) b
!!$       call M%apply(b, this%r)
!!$       b_nrm_M = b%dot(this%r)
!!$       if ( environment%am_i_fine_task() ) then ! Am I a fine task ?
!!$          b_nrm_M = sqrt(b_nrm_M)
!!$       end if
!!$    else
!!$       b_nrm_M = 0.0_rp
!!$    endif
!!$
!!$    ! 1) Compute initial residual
!!$    ! 1.a) r=Ax
!!$    call A%apply(x,this%r)
!!$
!!$    ! 1.b) r=b-r
!!$    call this%r%axpby(1.0, b ,-1.0)
!!$
!!$    ! 2) z=inv(M)r
!!$    call M%apply(this%r,this%z)
!!$
!!$    ! 3) <r,z>
!!$    r_z = this%r%dot(this%z)
!!$
!!$    if ( environment%am_i_fine_task() ) then ! Am I a fine task ?
!!$       r_nrm_M = sqrt( r_z )
!!$    end if
!!$
!!$    ! 4) Initializations:
!!$    ! p=z
!!$    call this%p%copy(this%z)
!!$
!!$    ! Init and log convergence 
!!$    call this%init_convergence_data(b, this%r, b_nrm_M, r_nrm_M)
!!$    call this%print_convergence_history_header(luout)
!!$
!!$    ! 5) Iteration
!!$    num_iterations = 0
!!$    loop_rgmres: do
!!$       num_iterations = num_iterations + 1
!!$
!!$       ! Ap = A*p
!!$       call A%apply(this%p,this%Ap)
!!$
!!$       ! <Ap,p>
!!$       Ap_p = this%Ap%dot(this%p)
!!$
!!$       ! Is this correct/appropriate ?
!!$       if (Ap_p /= 0.0_rp) then
!!$          alpha = r_z / Ap_p
!!$       else
!!$          alpha = 0.0_rp
!!$       end if
!!$
!!$       ! x = x + alpha*p
!!$       call x%axpby(alpha, this%p, 1.0_rp)
!!$
!!$       ! r = r - alpha*Ap
!!$       call this%r%axpby (-alpha, this%Ap, 1.0_rp)
!!$
!!$       ! Check and log convergence
!!$       call this%update_convergence_data_and_evaluate_stopping_criteria(this%r, r_nrm_M, alpha, this%p)
!!$       call this%print_convergence_history_new_line(luout)
!!$       
!!$       ! Send did_converge to coarse-grid tasks
!!$       call environment%bcast(did_converge)
!!$
!!$       if(did_converge .or.(num_iterations>=max_num_iterations)) exit loop_rgmres
!!$
!!$       ! z = inv(M) r
!!$       call M%apply(this%r,this%z)
!!$
!!$       if ( environment%am_i_fine_task()) then ! Am I a fine task ?
!!$          beta = 1.0_rp/r_z
!!$       else
!!$          beta = 0.0_rp
!!$       end if
!!$
!!$       r_z=this%r%dot(this%z)
!!$       beta=beta*r_z
!!$
!!$       if (environment%am_i_fine_task()) then ! Am I a fine task ?
!!$          r_nrm_M = sqrt( r_z )
!!$       end if
!!$
!!$       ! p = z + beta*p
!!$       call this%p%axpby(1.0,this%z,beta)
!!$    end do loop_rgmres
!!$
!!$    call this%print_convergence_history_footer(luout)
  end subroutine rgmres_solve_body

  function rgmres_supports_stopping_criteria(this,stopping_criteria)
    implicit none
    class(rgmres_t), intent(in) :: this
    integer(ip), intent(in) :: stopping_criteria
    logical :: rgmres_supports_stopping_criteria
    rgmres_supports_stopping_criteria = ( stopping_criteria == res_nrmgiven_rhs_nrmgiven .or. &
                                          stopping_criteria == res_nrmgiven_res_nrmgiven )
  end function rgmres_supports_stopping_criteria
  
  function rgmres_get_default_stopping_criteria(this)
    implicit none
    class(rgmres_t), intent(in) :: this
    integer(ip) :: rgmres_get_default_stopping_criteria
    rgmres_get_default_stopping_criteria = default_rgmres_stopping_criteria
  end function rgmres_get_default_stopping_criteria
  
  function create_rgmres(environment)
    implicit none
    class(environment_t), intent(in) :: environment
    class(base_linear_solver_t), pointer :: create_rgmres
    type(rgmres_t), pointer :: rgmres
    allocate(rgmres)
    call rgmres%set_environment(environment)
    call rgmres%set_name(rgmres_name)
    call rgmres%set_defaults()
    rgmres%dkrymax = default_rgmres_dkrymax
    call rgmres%set_state(start)
    create_rgmres => rgmres
  end function create_rgmres
  
end module rgmres_names
