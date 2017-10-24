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
module par_nsi_time_integration_names
  use fempar_names
  use par_nsi_nonlinear_operator_names
  use par_nsi_nonlinear_solver_names
  use par_nsi_discrete_integration_names
  implicit none
# include "debug.i90"
  private

!>  Integration of time dependent ADEs of the form 
!> \begin{equation*}
!>  u' + F(u) + G(u) =  0
!> \end{equation*}
!>  where u is a discrete FE function approximating several fields and
!>  F and G are the part of the residual that is evaluated implicitly
!>  or explicitly. The implementation of a general method is described 
!>  [here] (http://dx.doi.org/10.1080/10618562.2011.575368).
!> 
!>  To control de evaluation of F, G or u'+F we rely on functions
!>  provided by the discrete integration. This permits to integrate
!>   
!>  \begin{equation*}
!>  \alfa u^{i} + \beta F(u^{i}) = \alfa \sum_{j=0}^p c_{j} u^{n-j}
!>                                + \sum_{j=1}^{i-1} a_{ij} F(u^{j}) 
!>                                + \sum_{j=1}^{i-1} e_{ij} G(u^{j})
!>  \end{equation*}
!>
!>  This expression can be used for multistage methods (RK) or
!>  linear multistep methods (BDF and AM).

  type, abstract :: time_integration_t
     private
     real(rp)    :: initial_time
     real(rp)    :: current_time
     real(rp)    :: final_time
     real(rp)    :: time_step
     integer(ip) :: current_step
     integer(ip) :: number_time_steps
     type(fe_affine_operator_t)             , pointer     :: residual          => null()
     type(nonlinear_operator_t)                           :: nonlinear_operator
     type(nonlinear_solver_t)               , pointer     :: nonlinear_solver  => null()
     !type(fe_function_t)                                  :: solution
     type(fe_function_t)                                  :: solution_old
     class(vector_t)                        , allocatable :: dof_values
     class(par_nsi_discrete_integration_t), pointer     :: discrete_integration
     class(environment_t)                   , pointer     :: environment 
   contains
     procedure :: time_integration_create
     generic   :: create                   => time_integration_create
     procedure :: free                     => time_integration_free
     procedure :: apply                    => time_integration_apply
     procedure :: get_residual_coefficient => time_integration_get_residual_coefficient
     procedure :: print                    => time_integration_print
     procedure(initialize_time_step_interface), deferred :: initialize_time_step
     procedure(update_solution_interface)     , deferred :: update_solution
     procedure(get_coefficient_interface)     , deferred :: get_mass_coefficient
  end type time_integration_t

  public :: time_integration_t

  abstract interface
     subroutine initialize_time_step_interface(this,fe_space, solution)
       import :: time_integration_t, serial_fe_space_t, fe_function_t
       class(time_integration_t), intent(inout) :: this
       class(serial_fe_space_t) , intent(inout) :: fe_space
       class(fe_function_t)     , intent(inout) :: solution
     end subroutine initialize_time_step_interface
     subroutine update_solution_interface(this,fe_space,solution)
       import :: time_integration_t, serial_fe_space_t, vector_t, fe_function_t
       class(time_integration_t), intent(inout) :: this
       class(serial_fe_space_t) , intent(inout) :: fe_space
       class(fe_function_t)     , intent(inout) :: solution
     end subroutine update_solution_interface
     function get_coefficient_interface(this)
       import :: time_integration_t, rp
       class(time_integration_t), intent(inout) :: this
       real(rp) :: get_coefficient_interface
     end function get_coefficient_interface
  end interface

  type, extends(time_integration_t) :: theta_time_integration_t
     real(rp)    :: theta
   contains
     procedure :: theta_time_integration_create
     generic   :: create                   => theta_time_integration_create
     procedure :: initialize_time_step     => theta_time_integration_initialize_time_step    
     procedure :: update_solution          => theta_time_integration_update_solution         
     procedure :: get_mass_coefficient     => theta_time_integration_get_mass_coefficient    
  end type theta_time_integration_t

  type, extends(time_integration_t) :: bdf2_time_integration_t
     type(fe_function_t) :: solution_n
     type(fe_function_t) :: solution_n_1
   contains
     procedure :: initialize_time_step     => bdf2_time_integration_initialize_time_step    
     procedure :: update_solution          => bdf2_time_integration_update_solution         
     procedure :: get_mass_coefficient     => bdf2_time_integration_get_mass_coefficient    
  end type bdf2_time_integration_t

  public :: theta_time_integration_t, bdf2_time_integration_t

contains

!==============================================================================

  subroutine time_integration_create(this, residual, nonlinear_solver, initial_time, final_time, number_time_steps)
    implicit none
    class(time_integration_t)         , intent(inout) :: this
    type(fe_affine_operator_t), target, intent(in)    :: residual
    type(nonlinear_solver_t)  , target, intent(in)    :: nonlinear_solver
    real(rp)                          , intent(in)    :: initial_time
    real(rp)                          , intent(in)    :: final_time
    integer(ip)                       , intent(in)    :: number_time_steps
    class(discrete_integration_t), pointer :: discrete_integration
    class(serial_fe_space_t)     , pointer :: fe_space

    this%initial_time                = initial_time
    this%final_time                  = final_time
    this%number_time_steps           = number_time_steps
    this%current_step                = 0
    ! Non adaptive so far
    this%time_step                   = ( this%final_time - this%initial_time ) / real(number_time_steps,rp)
    this%current_time                = initial_time ! + this%time_step

    this%residual         => residual
    this%nonlinear_solver => nonlinear_solver
    fe_space              => this%residual%get_fe_space()
    this%environment      => fe_space%get_environment()

    call this%solution_old%create(fe_space)
    call this%residual%create_domain_vector(this%dof_values)

    ! Create a nonlinear fe operator
    call this%nonlinear_operator%create(residual)

    discrete_integration => this%residual%get_discrete_integration()
    select type(discrete_integration)
    class is(par_nsi_discrete_integration_t)
       this%discrete_integration => discrete_integration
    class default
       mcheck(.false.,'nonlinear operator only works with par_nsi_discrete_integration')
    end select

  end subroutine time_integration_create

  subroutine time_integration_free(this)
    implicit none
    class(time_integration_t)         , intent(inout) :: this
    integer(ip) :: istat
    
    this%residual         => null()
    this%nonlinear_solver => null()
    this%environment      => null()
    call this%solution_old%free()
    call this%dof_values%free()
    deallocate(this%dof_values,stat=istat); check(istat==0)
    call this%nonlinear_operator%free()
    this%discrete_integration => null()

  end subroutine time_integration_free
  
  subroutine time_integration_apply(this,solution)
    implicit none
    class(time_integration_t), intent(inout) :: this 
    class(fe_function_t)     , intent(inout) :: solution
    class(serial_fe_space_t) , pointer       :: fe_space
    class(vector_t)          , pointer       :: dof_values
    
    fe_space => this%residual%get_fe_space()
    call this%discrete_integration%set_fe_function(solution) 
    call this%discrete_integration%set_old_fe_function(this%solution_old)
    this%solution_old = solution

    do while(this%current_time < this%final_time )

       call this%initialize_time_step(fe_space, solution)
       call this%print()

       call this%discrete_integration%set_mass_coefficient    (this%get_mass_coefficient())
       call this%discrete_integration%set_residual_coefficient(this%get_residual_coefficient())
       call this%discrete_integration%set_current_time        (this%current_time)

       ! Solve time step
       dof_values => solution%get_free_dof_values() ! initial guess is the previous step
       call this%nonlinear_solver%solve(this%nonlinear_operator,dof_values)

       ! Check if converged
       mcheck( this%nonlinear_solver%has_converged(), 'Nonlinear solver has not converged. Cannot advance to the next step.' )

       call this%update_solution(fe_space,solution)

    end do

  end subroutine time_integration_apply

    ! Needs to be overwritten for some methods
  function time_integration_get_residual_coefficient(this)
    implicit none
    class(time_integration_t)      , intent(inout) :: this 
    real(rp) :: time_integration_get_residual_coefficient
    time_integration_get_residual_coefficient = 1.0_rp 
  end function time_integration_get_residual_coefficient

  subroutine time_integration_print(this)
    implicit none
    class(time_integration_t), intent(in) :: this
    if ( this%environment%am_i_l1_root() ) then
      write(*,*) '========================================================================'
      write(*,'(a10,i6,a20, e10.3,a12,e10.3)') 'Time step ', this%current_step, ': Solving for t=', this%current_time, 'with dt', this%time_step
    end if
  end subroutine time_integration_print

!==============================================================================

  subroutine theta_time_integration_create(this, residual, nonlinear_solver, initial_time, final_time, number_time_steps, theta)
    implicit none
    class(theta_time_integration_t)         , intent(inout) :: this
    type(fe_affine_operator_t), target, intent(in)    :: residual
    type(nonlinear_solver_t)  , target, intent(in)    :: nonlinear_solver
    real(rp)                          , intent(in)    :: theta
    real(rp)                          , intent(in)    :: initial_time
    real(rp)                          , intent(in)    :: final_time
    integer(ip)                       , intent(in)    :: number_time_steps
    massert(theta >  0.0, 'theta method: value of thera should be greater than 0')
    massert(theta <= 1.0, 'theta method: value of thera should be smaller or equal than 1')
    this%theta = theta
    call this%create(residual, nonlinear_solver, initial_time, final_time, number_time_steps)
    this%current_time = this%current_time - (1-this%theta) * this%time_step 
  end subroutine theta_time_integration_create

  subroutine theta_time_integration_initialize_time_step(this,fe_space,solution)
    implicit none
    class(theta_time_integration_t), intent(inout) :: this
    class(serial_fe_space_t)       , intent(inout) :: fe_space
    class(fe_function_t)           , intent(inout) :: solution
    this%current_step = this%current_step + 1
    this%current_time = this%current_time + this%time_step
    ! Update strong Dirichlet Bcs at current time
    call fe_space%interpolate_dirichlet_values (solution, this%current_time)
  end subroutine theta_time_integration_initialize_time_step

  subroutine theta_time_integration_update_solution(this,fe_space,solution)
    implicit none
    class(theta_time_integration_t), intent(inout) :: this
    class(serial_fe_space_t), intent(inout) :: fe_space
    class(fe_function_t)    , intent(inout) :: solution
    class(vector_t)         , pointer       :: dof_values_current
    class(vector_t)         , pointer       :: dof_values_old
    dof_values_current  => solution%get_free_dof_values()
    dof_values_old      => this%solution_old%get_free_dof_values()
    dof_values_current = (1.0_rp/this%theta) * ( dof_values_current - (1.0_rp-this%theta) * dof_values_old )
    ! Update strong Dirichlet Bcs at the end of the step time
    call fe_space%interpolate_dirichlet_values   ( solution, this%current_time + (1.0_rp-this%theta) * this%time_step )
    ! Store solution
    this%solution_old = solution
  end subroutine theta_time_integration_update_solution

  function theta_time_integration_get_mass_coefficient(this)
    implicit none
    class(theta_time_integration_t), intent(inout) :: this
    real(rp) :: theta_time_integration_get_mass_coefficient
    theta_time_integration_get_mass_coefficient = 1.0_rp /(this%theta*this%time_step)
  end function theta_time_integration_get_mass_coefficient

!==============================================================================

  subroutine bdf2_time_integration_initialize_time_step(this, fe_space,solution)
    implicit none
    class(bdf2_time_integration_t), intent(inout) :: this
    class(serial_fe_space_t)      , intent(inout) :: fe_space
    class(fe_function_t)          , intent(inout) :: solution
    this%current_step = this%current_step + 1
    this%current_time = this%current_time + this%time_step
    ! Update strong Dirichlet Bcs
    call fe_space%interpolate_dirichlet_values(solution, this%current_time + this%time_step)
  end subroutine bdf2_time_integration_initialize_time_step

  subroutine bdf2_time_integration_update_solution(this,fe_space,solution)
    implicit none
    class(bdf2_time_integration_t), intent(inout) :: this
    class(serial_fe_space_t), intent(inout) :: fe_space
    class(fe_function_t)    , intent(inout) :: solution
    ! Store solution
    this%solution_n_1 = this%solution_n
    this%solution_n = solution
    this%solution_old = this%solution_n
    call this%solution_old%axpby(4.0_rp/3.0_rp,this%solution_n_1,-1.0_rp/3.0_rp)
  end subroutine bdf2_time_integration_update_solution

  function bdf2_time_integration_get_mass_coefficient(this)
    implicit none
    class(bdf2_time_integration_t), intent(inout) :: this
    real(rp) :: bdf2_time_integration_get_mass_coefficient
    bdf2_time_integration_get_mass_coefficient = 3.0_rp /(2.0_rp*this%time_step)
  end function bdf2_time_integration_get_mass_coefficient

 end module par_nsi_time_integration_names
