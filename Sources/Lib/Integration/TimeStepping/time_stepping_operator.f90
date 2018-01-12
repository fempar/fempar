

!type :: time_stepping_operator_t
!private
!  class(fe_nonlinear_operator_t)              :: fe_op
!  class(mass_discrete_integration_t), pointer :: mass_integration
!  type(time_stepping_scheme_t)                :: scheme
!  class(assembler_t)                , pointer :: assembler => NULL()
!  class(vector_t)                   , pointer :: initial_value
!  integer(ip)                                 :: dt
!  class(assembler_t)                , pointer :: mass => NULL()
!  type(block_vector_t)                        :: dofs_stages(:)
!contains
!  procedure :: create => time_stepping_operator_create
!  procedure :: apply  => time_stepping_operator_apply
!end type time_stepping_operator_t

!!This operator is the global RK operator as follows.
!!In RK methods, for a s-stage method, given u_0, we must compute
!![ v_1, v_s ] st
!!v_0 - u_0 = 0
!!v_1 + A ( t + c_1 dt , v_0 + a_11 v_1 ) = 0
!!...
!!v_i + A ( t + c_i dt, v_0 + \sum_{j=1}^s a_ij v_j ) = 0
!!...
!!to finally get
!!u_1 - u_0 - \sum_{i=1}^s b_i v_i = 0.
!! As a result, the RK problem can be stated in compact form as R(U)=0, where
!! U = [ v_0, v_1, ..., v_s, v_s+1]. The objective is to solve this problem and 
!! extract u_1 = v_s+1.
!! One must compute, given U, its application R(U). So, I want a apply method that, 
!! given the stage i I want to get, it applies the i-stage R_i(U) application.
!! However, in DIRK RK methods, the ones usually used, the coupling is weaker, and 
!! should exploit it. R_i(U) only depends on v_j, for j \leq i.
!! In this case, we want to solve at every stage
!! R_i( v_0, ..., v_i ) = 0, where v_0, ..., v_i-1 are known, i.e.,
!! R*_i(v_i) = 0.
!! In order to solve this potential nonlinear operator, we must provide its tangent too.
!! We need the nonlinear operator: Given i, and U_1, U_i-1, U:i+1,...,U_s, 
!! R_ij(V) = R_i (U_1, ...,U_i-1,V,U_i+1,...,U_s)
!! Next, I compute the apply and tangent of this operator, which is all I need to apply
!! the whole operator...

!subroutine time_stepping_operator_create(...)
!  ! sbadia : get initial condition (to be done)
!  allocate dofs_stages(this%op%scheme%stages) !?

!end subroutine time_stepping_operator_create

!subroutine time_stepping_operator_set_initial_data( this, x0 )
!class(time_stepping_operator_t), intent(inout) :: this
!class(vector_t), intent(in) :: x0
!this%initial_value = x0
!! allocate stages with this ...
!end subroutine time_stepping_operator_set_initial_data

!subroutine time_stepping_operator_apply_row( this, x, y, i )
!class(time_stepping_operator_t), intent(inout) :: this
!integer(ip)                    , intent(in)    :: row
!class(block_vector_t)          , intent(in)    :: x
!class(vector_t)                , intent(inout) :: y_i
!! sbadia : We should set the time in the nonlinear operator too for variable body force / bc's
!call this%set_evaluation_point_row(y,i)
!call this%fe_op%compute_residual()
!y_i = this%fe_op%get_residual() * this%dt
!y_i = y_i + this%mass%get_matrix*dofs_stages(i)
!end subroutine time_stepping_operator_apply_row

!! Compute  \partial R*_i / \partial v_j = M  \delta_ij
!!                 + \partial A / \partial v * a_ij
!subroutine time_stepping_operator_compute_tangent_block( this, x, tangent, i, j )
!class(time_stepping_operator_t), intent(inout) :: this
!integer(ip)                    , intent(in)    :: row
!class(block_vector_t)          , intent(in)    :: x
!class(vector_t)                , intent(inout) :: y_i
!call this%set_evaluation_point_row(y,i)
!this%fe_operator%compute_tangent()
!tangent = this%fe_operator%get_tangent() * ( dt * this%scheme%a(i,j) )
!end subroutine time_stepping_operator_compute_tangent_block

!subroutine time_stepping_operator_set_evaluation_point_row( this, x, y_i, i )
!class(time_stepping_operator_t), intent(inout) :: this
!integer(ip)                    , intent(in)    :: row
!class(block_vector_t)          , intent(in)    :: x
!class(vector_t)                , intent(in)    :: y_i
!class(vector_t)  :: aux
!integer(ip) :: j
!x(i) = y_i 
!aux = this%initial_value
!do j = 1, this%scheme%stages
!   if ( this%scheme%a(i,j) /= 0.0_rp ) aux = aux + dofs_stages(j)*this%scheme%a(i,j)
!end do
!call this%fe_op%set_evaluation_point(aux)
!end subroutine time_stepping_operator_set_evaluation_point_row

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!type, extends(operator_t) :: time_stepping_operator_stage_t
!private
!  type(time_stepping_operator_t) :: op
!  integer(ip) :: i
!  integer(ip) :: j
!contains
!  subroutine time_stepping_operator_block_apply( this, x, y )
!  implicit none
!  class(time_stepping_operator_t) :: op
!  class(vector_t), intent(in) :: x
!  class(vector_t), intent(inout) :: y
!  call op%time_stepping_operator_apply_row(x,y,i)
!  end subroutine time_stepping_operator_block_apply
!  
!  subroutine time_stepping_operator_block_compute_tangent( this )
!  implicit none
!  class(time_stepping_operator_t) :: op
!  call this%op%compute_tangent(i,j)
!  end subroutine time_stepping_operator_compute_tangent
!  
!  subroutine time_stepping_operator_block_get_tangent( this )
!  implicit none
!  class(time_stepping_operator_t) :: op
!  call this%op%get_tangent(i,j)
!  end subroutine time_stepping_operator_get_tangent
!  
!    subroutine time_stepping_operator_block_set_evaluation_point( this, y )
!  implicit none
!  class(time_stepping_operator_t) :: op
!  class(vector_t) :: y
!  call this%op%set_evaluation_point_row( y, i)
!  end subroutine time_stepping_operator_set_evaluation_point
!end type time_stepping_operator_block_stage_t
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!type :: dirk_time_stepping_solver_t
!private
!  type(time_stepping_operator_t) :: op
!  type(time_stepping_operator_block_t) :: op_block
!  type(nonlinear_solver_t) :: nl_solver
!  
!  subroutine dirk_time_stepping_solver_apply( this, x, y )
!  ! Solve the time stepping problem
!  allocate sol(this%op%scheme%stages)
!  call op%set_initial_solution(x)
!  call op_block%create op_block(op)
!  
!  ! loop time
!  ! at every time step
!  y = x
!  do i = 1, this%op%scheme%stages
!    call this%op_block%set_stages(i,i)
!    A => this%op_block%get_operator()                              ! nonlinear operator
!    call nl_solver%apply ( this%op_block, this%op%stages_dofs(i) ) ! nonlinear solver
!    y = x + this%op%stages_dofs(i)*this%op%scheme%b(i)
!  end do
!  end subroutine dirk_time_stepping_solver_apply
!  

!  
!  
!  
!  




















!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!integer(ip), parameter :: forward_euler      = 0
!integer(ip), parameter :: backward_euler     = 1
!integer(ip), parameter :: crank_nicolson     = 2 
!! other integrators here
!type :: time_stepping_scheme_t
!    private
!    integer(ip) :: time_integrator
!    integer(ip) :: order
!    integer(ip) :: stages
!    real(rp)    , allocatable :: a(:,:)
!    real(rp)    , allocatable :: b(:)
!    real(rp)    , allocatable :: c(:)
!  contains
!    procedure :: create => time_stepping_scheme_create
!  public :: time_stepping_scheme_t
!end time_stepping_scheme_t

!contains

!subroutine time_stepping_scheme_create( this, ti_type )
!implicit none
!type(time_stepping_scheme_t), intent(inout) :: this
!integer(ip), intent(in) ti_type

!select case ( ti_type )
!case ( backward_euler )
!  this%order = 1
!  this%stages = 1
!  call this%allocate_arrays()
!  this%a = (/ 1.0_rp /)
!  this%b = (/ 1.0_rp /)
!  this%c = (/ 1.0_rp /)
!case ( trapezoidal_rule )
!  this%order = 2
!  this%stages = 2
!  this%a(1:2,1:2) = (/ (/ 0.0_rp,  0.0_rp /), (/ 0.5_rp, 0.5_rp /) /)
!  this%b(1:2) = (/ 0.5_rp, 0.5_rp /) 
!  ! this%bs = [ 1.0_rp 0.0_rp] for adaptive RK methods
!  this%c(1:2) = (/ 0.0_rp, 1.0_rp /)
!case default
!  check(.false.)
!end select

!end subroutine time_stepping_scheme_create

!subroutine allocate_butcher_tableau ( this )
!implicit none
!type(time_stepping_scheme_t), intent(inout) :: this
!integer(ip) :: aux
!  aux = this%stages
!  call memalloc(aux,aux,this%a,__FILE__,__LINE__)
!  call memalloc(aux,this%b,__FILE__,__LINE__)
!  call memalloc(aux,this%c,__FILE__,__LINE__)
!end subroutine allocate_butcher_tableau

!  
!  
!  
!  

