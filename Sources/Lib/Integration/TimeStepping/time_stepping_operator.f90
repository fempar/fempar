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
module fe_nonlinear_operator_names
  ! sbadia: to check whether all this needed
  use types_names
  use memor_names
  
  use vector_space_names
  use fe_space_names
  use operator_names
  use vector_names
  
  use assembler_names
  use sparse_assembler_names
  use block_sparse_assembler_names  
  use par_sparse_assembler_names
  
  use sparse_matrix_names, only: sparse_matrix_t
  use block_sparse_matrix_names
  use par_sparse_matrix_names
  
  use serial_scalar_array_names
  use serial_block_array_names
  use par_scalar_array_names
  
  use array_names
  use matrix_names
  use discrete_integration_names
  use environment_names
  use direct_solver_names
  use block_layout_names
  
  implicit none
# include "debug.i90"
  
  private
  
  integer(ip), parameter :: forward_euler      = 0
  integer(ip), parameter :: backward_euler     = 1
  integer(ip), parameter :: crank_nicolson     = 2 
  
  ! states to be defined
  integer(ip), parameter :: created             = 0
  integer(ip), parameter :: residual_computed   = 1 
  integer(ip), parameter :: tangent_computed    = 2 
  integer(ip), parameter :: assembler_computed  = 3
  
  
  ! sbadia: to make auto-documentation style
  !
  !This operator is the global RK operator as follows.
  !In RK methods, for a s-stage method, given u_0, we must compute
  ![ v_1, v_s ] st
  !v_1 + A ( t + c_1 dt , v_0 + a_11 v_1 ) = 0
  !...
  !v_i + A ( t + c_i dt, v_0 + \sum_{j=1}^s a_ij v_j ) = 0
  !...
  !to finally get
  !u_1 =  u_0 + \sum_{i=1}^s b_i v_i = 0.
  ! As a result, the RK problem can be stated in compact form as R(U)=0, where
  ! U = [v_1, ..., v_s]. The objective is to solve this problem and 
  ! extract u_1 = u_0 + \sum_{i=1}^s b_i v_i = 0.
  ! One must compute, given U, its application R(U). So, I want a apply method that, 
  ! given the stage i I want to get, it applies the i-stage R_i(U) application.
  ! However, in DIRK RK methods, the ones usually used, the coupling is weaker, and 
  ! should exploit it. R_i(U) only depends on v_j, for j \leq i.
  ! In this case, we want to solve at every stage
  ! R_i( v_0, ..., v_i ) = 0, where v_0, ..., v_i-1 are known, which we can denote as
  ! R*_i(v_i) = 0.
  ! In order to solve this potential nonlinear operator, we must provide its tangent too.
  ! We need the nonlinear operator: Given i, and v_1, ..., v_i-1, v_i+1,...,v_s, 
  ! R_ij(V) = R_i (v_1, ...,v_i-1,V,v_i+1,...,v_s).
  ! Next, I compute the apply and tangent of this operator, which is all I need to apply
  ! the whole operator.

  ! sbadia: I must provide a method to get the solution
  ! u_1 =  u_0 + \sum_{i=1}^s b_i v_i = 0.
  type, extends(operator_t) :: time_stepping_operator_t
     private
     class(fe_nonlinear_operator_t)                     :: fe_op
     class(mass_matrix_discrete_integration_t)          :: mass_integration
     type(time_stepping_scheme_t)                       :: scheme
     class(vector_t)                          , pointer :: initial_value
     real(rp)                                           :: dt
     ! sbadia: I consider the stages vector as a block vector. Good idea?
     type(block_vector_t)                               :: dofs_stages(:)
     ! sbadia: For the mass matrix probably better a fe_affine_operator...
     class(assembler_t)                                 :: mass
     ! sbadia: For the moment, we are not interested in full RK implementations,
     ! even though it would be an easy paper about preconditioning these schemes
     ! So, we don't really need to use the block assembler for the
     ! all-stages operator. We note that the matrix is not needed to be computed
     ! for every stage since it is always the same
     class(assembler_t)                       , pointer :: assembler => NULL()
   contains
     procedure :: create             => time_stepping_operator_create
     ! sbadia: It must be defined since it is an operator, but for the moment
     ! we do not want to use it. Dummy implementation...
     procedure :: apply              => time_stepping_operator_apply
     ! sbadia: to be implemented
	 procedure :: free               => time_stepping_operator_free
     procedure :: set_initial_data   => time_stepping_operator_set_initial_data
     procedure :: set_time_step_size => time_stepping_operator_set_time_step_size
     
     procedure, private :: allocate_dofs_stages      => time_stepping_operator_allocate_dofs_stages
     procedure, private :: apply_row                 => time_stepping_operator_apply_row
     procedure, private :: compute_tangent_block     => time_stepping_operator_compute_tangent_block
     procedure, private :: compute_mass_matrix       => time_stepping_operator_compute_mass_matrix
     procedure, private :: set_evaluation_point_row  => time_stepping_operator_set_evaluation_point_row	 
  end type time_stepping_operator_t
  
  type, extends(operator_t) :: time_stepping_operator_stage_t
private
type(time_stepping_operator_t) :: op
integer(ip) :: i
integer(ip) :: j
class(assembler_t) :: assembler
contains
procedure :: create => time_stepping_operator_block_create
procedure :: apply  => time_stepping_operator_block_apply
procedure :: compute_residual => time_stepping_operator_block_compute_residual
procedure :: compute_tangent => time_stepping_operator_block_compute_tangent( this )
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 call this%op%compute_tangent_block(this%assembler%get_tangent(),i,j)
end subroutine time_stepping_operator_block_compute_tangent

subroutine time_stepping_operator_block_get_residual( this )
 implicit none
 class(time_stepping_operator_stage_t), intent(in) :: this
 call this%op%get_tangent(i,j)
end subroutine time_stepping_operator_block_get_residual

subroutine time_stepping_operator_block_get_residual(this) result(residual) 
 implicit none
 class(fe_nonlinear_operator_t), intent(inout)    :: this
 class(array_t)  , pointer       :: array
 residual = this%assembler%get_array()
end subroutine time_stepping_operator_block_get_residual

function time_stepping_operator_block_get_tangent(this) result(tangent)
 implicit none
 class(time_stepping_operator_stage_t), intent(in) :: this
 type(lvalue_operator_t)          :: tangent
 assert ( associated(this%state)  )
 tangent = this%assembler%get_matrix()
 call tangent%SetTemp()
end function time_stepping_operator_block_get_tangent

subroutine time_stepping_operator_block_set_evaluation_point( this, x )
 implicit none
 class(vector_t) :: x
 call this%op%set_evaluation_point_row( x, i)
end subroutine time_stepping_operator_block_set_evaluation_point
end type time_stepping_operator_block_stage_t
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine time_stepping_operator_create(this, fe_op, scheme_type)
    class(time_stepping_operator_t), intent(inout) :: this
    class(fe_nonlinear_operator_t) , intent(in)    :: fe_op
    integer(ip), intent(in) :: scheme_type
    this%fe_op => fe_op
    call this%scheme%create(scheme_type)
    num_stages = this%scheme%stages
    call this%dofs_stages%create(num_stages)
  end subroutine time_stepping_operator_create
  
  subroutine time_stepping_operator_set_initial_data( this, x0 )
    class(time_stepping_operator_t), intent(inout) :: this
    class(vector_t), intent(in) :: x0
    this%initial_value = x0
    ! check compatible w/ fe_op vector space
    call this%fe_op%abort_if_not_in_domain(x0)
    call this%allocate_dofs_stages()
  end subroutine time_stepping_operator_set_initial_data
  
  subroutine time_stepping_operator_allocate_dofs_stages( this, x0 )
    class(time_stepping_operator_t), intent(inout) :: this
    class(vector_t), intent(in) :: x0
    type(vector_space_t), pointer :: vector_space
    vector_space => this%fe_op%get_domain_vector_space()
    ! check what happens if dofs_stages not allocated
    if ( .not. vector_space%belongs(this%dofs_stages(1)) ) then 
       do i = 1, this%scheme%stages
          call this%blocks(i)%vector%clone(x0)
       end do
    end if
  end subroutine time_stepping_operator_allocate_dofs_stages
  
  subroutine time_stepping_operator_apply_row( this, x, y, i )
    class(time_stepping_operator_t), intent(inout) :: this
    integer(ip)                    , intent(in)    :: row
    class(vector_t)          , intent(in)    :: x
    class(vector_t)                , intent(inout) :: y
    ! sbadia : We should set the time in the nonlinear operator too for variable body force / bc's
    call this%set_evaluation_point_row(x,i)
    call this%fe_op%compute_residual()
    y = this%fe_op%get_residual() * this%dt
    ! Here we should provide also the trial fe space
    ! We should only do it when needed, e.g.,
    call this%compute_mass_matrix()
    y = y + this%mass%get_matrix()*dofs_stages(i)
  end subroutine time_stepping_operator_apply_row
  
  ! Compute  \partial R*_i / \partial v_j = M  \delta_ij
  !                 + \partial A / \partial v * a_ij
  subroutine time_stepping_operator_compute_tangent_block( this, tangent, mass, i, j )
    class(time_stepping_operator_t), intent(inout) :: this
    integer(ip)                    , intent(in)    :: i, j
    class(lvalue_operator_t)       , intent(inout) :: tangent
    class(matrix_t)                , intent(inout) :: mass
    !call this%set_evaluation_point_row(x,i)
    this%fe_operator%compute_tangent()
    ! sbadia : error, get_tangent not a matrix but lvalue operator
    tangent = this%fe_operator%get_tangent() * ( dt * this%scheme%a(i,j) )
    if ( i == j ) then
       ! only when needed probably based on vector spaces compatibility
       call this%compute_mass_matrix()
       tangent = tangent + this%mass%get_matrix()
    end if
  end subroutine time_stepping_operator_compute_tangent_block

  subroutine time_stepping_operator_compute_mass_matrix( this )
    class(time_stepping_operator_t), intent(inout) :: this
    call this%mass_integration%integrate( this%fe_op%get_fe_space(), this%mass )
  end subroutine time_stepping_operator_compute_mass_matrix
  
  subroutine time_stepping_operator_set_evaluation_point_row( this, y_i, i )
    class(time_stepping_operator_t), intent(inout) :: this
    integer(ip)                    , intent(in)    :: row
    class(vector_t)                , intent(in)    :: y_i
    class(vector_t)  :: aux
    integer(ip) :: j
    this%dofs_stages(i) = y_i
    aux = this%initial_value
    do j = 1, this%scheme%stages
       if ( this%scheme%a(i,j) /= 0.0_rp ) aux = aux + this%dofs_stages(j)*this%scheme%a(i,j)
    end do
    call this%fe_op%set_evaluation_point(aux)
  end subroutine time_stepping_operator_set_evaluation_point_row
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type, extends(operator_t) :: time_stepping_operator_stage_t
private
type(time_stepping_operator_t) :: op
integer(ip) :: i
integer(ip) :: j
class(assembler_t) :: assembler
contains
subroutine time_stepping_operator_block_create(this, op)
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 class(time_stepping_operator_t), intent(in)    :: op
 this%op => op
 ! sbadia: allocate assembler here?
end subroutine time_stepping_operator_block_create
subroutine time_stepping_operator_block_apply
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 class(vector_t), intent(in) :: x
 class(vector_t), intent(inout) :: y
 call this%set_evaluation_point(x)
 call this%compute_residual()
 y = this%get_residual()
end subroutine time_stepping_operator_block_apply
subroutine time_stepping_operator_block_apply( this, x, y )
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 class(vector_t), intent(in) :: x
 class(vector_t), intent(inout) :: y
 call this%set_evaluation_point(x)
 call this%compute_residual()
 y = this%get_residual()
end subroutine time_stepping_operator_block_apply

subroutine time_stepping_operator_block_compute_residual( this )
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 call this%op%compute_residual_block(this%assembler%get_array(),i)
end subroutine time_stepping_operator_block_compute_residual

subroutine time_stepping_operator_block_compute_tangent( this )
 implicit none
 class(time_stepping_operator_stage_t), intent(inout) :: this
 call this%op%compute_tangent_block(this%assembler%get_tangent(),i,j)
end subroutine time_stepping_operator_block_compute_tangent

subroutine time_stepping_operator_block_get_residual( this )
 implicit none
 class(time_stepping_operator_stage_t), intent(in) :: this
 call this%op%get_tangent(i,j)
end subroutine time_stepping_operator_block_get_residual

subroutine time_stepping_operator_block_get_residual(this) result(residual) 
 implicit none
 class(fe_nonlinear_operator_t), intent(inout)    :: this
 class(array_t)  , pointer       :: array
 residual = this%assembler%get_array()
end subroutine time_stepping_operator_block_get_residual

function time_stepping_operator_block_get_tangent(this) result(tangent)
 implicit none
 class(time_stepping_operator_stage_t), intent(in) :: this
 type(lvalue_operator_t)          :: tangent
 assert ( associated(this%state)  )
 tangent = this%assembler%get_matrix()
 call tangent%SetTemp()
end function time_stepping_operator_block_get_tangent

subroutine time_stepping_operator_block_set_evaluation_point( this, x )
 implicit none
 class(vector_t) :: x
 call this%op%set_evaluation_point_row( x, i)
end subroutine time_stepping_operator_block_set_evaluation_point
end type time_stepping_operator_block_stage_t
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
type :: dirk_time_stepping_solver_t
private
type(time_stepping_operator_t) :: op
type(time_stepping_operator_block_t) :: op_block
type(nonlinear_solver_t) :: nl_solver

subroutine dirk_time_stepping_solver_apply( this, x, y )
  ! x could be the initial value u_0 and y the final solution u_1
  ! Solve the time stepping problem
 call op%set_initial_solution(x)
 call op_block%create op_block(op)

 ! loop time
 ! at every time step
 y = x
 do i = 1, this%op%scheme%stages
    call this%op_block%set_stages(i,i)
    A => this%op_block%get_operator()                              ! nonlinear operator
    call nl_solver%apply ( this%op_block, this%op%stages_dofs(i) ) ! nonlinear solver
    y = x + this%op%stages_dofs(i)*this%op%scheme%b(i)
 end do
end subroutine dirk_time_stepping_solver_apply


























!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! other integrators here
type :: time_stepping_scheme_t
  private
  integer(ip) :: time_integrator
  integer(ip) :: order
  integer(ip) :: stages
  real(rp)    , allocatable :: a(:,:)
  real(rp)    , allocatable :: b(:)
  real(rp)    , allocatable :: c(:)
contains
  procedure :: create => time_stepping_scheme_create
  public :: time_stepping_scheme_t
end type time_stepping_scheme_t

contains

subroutine time_stepping_scheme_create( this, ti_type )
 implicit none
 type(time_stepping_scheme_t), intent(inout) :: this
 integer(ip), intent(in) ti_type

 select case ( ti_type )
 case ( backward_euler )
    this%order = 1
    this%stages = 1
    call this%allocate_butcher_tableau()
    this%a = (/ 1.0_rp /)
    this%b = (/ 1.0_rp /)
    this%c = (/ 1.0_rp /)
 case ( trapezoidal_rule )
    this%order = 2
    this%stages = 2
    call this%allocate_butcher_tableau()
    this%a(1:2,1:2) = (/ (/ 0.0_rp,  0.0_rp /), (/ 0.5_rp, 0.5_rp /) /)
    this%b(1:2) = (/ 0.5_rp, 0.5_rp /) 
    ! this%bs = [ 1.0_rp 0.0_rp] for adaptive RK methods
    this%c(1:2) = (/ 0.0_rp, 1.0_rp /)
 case default
    check(.false.)
 end select

end subroutine time_stepping_scheme_create

subroutine allocate_butcher_tableau ( this )
 implicit none
 type(time_stepping_scheme_t), intent(inout) :: this
 integer(ip) :: aux
 aux = this%stages
 call memalloc(aux,aux,this%a,__FILE__,__LINE__)
 call memalloc(aux,this%b,__FILE__,__LINE__)
 call memalloc(aux,this%c,__FILE__,__LINE__)
end subroutine allocate_butcher_tableau



module mass_matrix_discrete_integration_names
 use fempar_names

 implicit none
# include "debug.i90"
 private
 type, extends(discrete_integration_t) :: mass_matrix_discrete_integration_t
contains
 procedure :: integrate_galerkin
end type mass_matrix_discrete_integration_t

public :: mass_matrix_discrete_integration_t

contains
  ! sbadia: not ready yet for petrov galerkin 
subroutine integrate_galerkin ( this, fe_space, assembler )
implicit none
class(mass_matrix_discrete_integration_t), intent(in)    :: this
class(serial_fe_space_t)                 ,intent(inout)  :: fe_space
class(assembler_t)                       , intent(inout) :: assembler

! FE space traversal-related data types
class(fe_cell_iterator_t), allocatable :: fe

! FE integration-related data types
type(quadrature_t)       , pointer :: quad
type(point_t)            , pointer :: quad_coords(:)
real(rp)            , allocatable  :: shape_values(:,:)

! FE matrix and vector i.e., A_K + f_K
real(rp), allocatable              :: elmat(:,:)

integer(ip)  :: istat
integer(ip)  :: qpoint, num_quad_points
integer(ip)  :: idof, jdof, num_dofs, max_num_dofs

assert (associated(this%fe_function))   
call fe_space%set_up_cell_integration()
call fe_space%create_fe_cell_iterator(fe)

max_num_dofs = fe_space%get_max_num_dofs_on_a_cell()
call memalloc ( max_num_dofs, max_num_dofs, elmat, __FILE__, __LINE__ )

call fe_space%create_fe_cell_iterator(fe)
do while ( .not. fe%has_finished() )

   ! Update FE-integration related data structures
   ! sbadia: Possible optimization, no grad computation required
   call fe%update_integration()

   quad            => fe%get_quadrature()
   num_quad_points =  quad%get_num_quadrature_points()
   num_dofs        =  fe%get_num_dofs()

   ! Compute element matrix and vector
   elmat = 0.0_rp
   call fe%get_values(shape_values)
   do qpoint = 1, num_quad_points
      factor = fe%get_det_jacobian(qpoint) * quad%get_weight(qpoint)
      do idof = 1, num_dofs
         do jdof = 1, num_dofs
            ! A_K(i,j) = (grad(phi_i),grad(phi_j))
            elmat(idof,jdof) = elmat(idof,jdof) + factor * shape_values(jdof,qpoint) * shape_values(idof,qpoint)
         end do
      end do
   end do
   call fe%assembly( this%fe_function, elmat, assembler )
   !end if
   call fe%next()
end do
call fe_space%free_fe_cell_iterator(fe)

call memfree(shape_values, __FILE__, __LINE__)
call memfree ( elmat, __FILE__, __LINE__ )
end subroutine integrate_galerkin

end module mass_matrix_discrete_integration_names








